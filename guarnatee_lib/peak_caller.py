import os.path
import logging
import numpy as np
import pandas as pd
from scipy import signal, stats
from numba import njit, prange
from numpy.lib.stride_tricks import sliding_window_view
# Assuming this import exists in your environment
from guarnatee_lib.helpers import Helpers

logger = logging.getLogger(__name__)
np.seterr(divide="ignore", invalid="ignore")


class PeakCaller:
    def __init__(
            self,
            raw_signal: np.array,
            min_peak_distance=10,
            min_height=0.0,
            min_step_factor=1.0,
            is_reversed=False,
            prefix="",
            thres_factor=1.0,
            is_coarse=False
    ):
        self.raw_signal = raw_signal
        # Hampel window size (points to look at on either side)
        self.window_size = 30
        self.min_peak_distance = int(min_peak_distance)
        self.min_height = float(min_height)
        self.min_step_factor = float(min_step_factor)
        self.is_reversed = is_reversed
        self.thres_factor = thres_factor
        self.prefix = prefix
        self.is_coarse = is_coarse
        self.peaks_arr = self.call_signal_peaks()

    def call_signal_peaks(self) -> np.array:
        """
        Multi-scale peak detection: runs Savitzky-Golay + Rolling MAD at
        several window sizes and merges results to catch peaks at different scales.
        """
        polyorder = 1
        if self.is_coarse:
            window_sizes = [15, 21, 31, 41, 51]
        else:
            window_sizes = [5, 7, 11, 15, 21]

        all_peaks = []
        all_heights = []

        for smooth_win in window_sizes:
            peaks, heights = self._detect_peaks_single_scale(smooth_win, polyorder)
            if peaks.size > 0:
                all_peaks.append(peaks)
                all_heights.append(heights)

        if not all_peaks:
            return None

        merged_peaks = np.concatenate(all_peaks)
        merged_heights = np.concatenate(all_heights)

        peaks, peak_heights_deriv = self._deduplicate_peaks(
            merged_peaks, merged_heights, self.min_peak_distance
        )

        adjusted_peaks = peaks + 1 if not self.is_reversed else peaks
        peak_indices = adjusted_peaks.astype(int)

        max_idx = len(self.raw_signal)
        valid_mask = (peak_indices >= 3) & (peak_indices < max_idx - 15)
        peak_indices = peak_indices[valid_mask]
        peak_heights_deriv = peak_heights_deriv[valid_mask]

        if len(peak_indices) == 0:
            return None

        raw_heights = self.raw_signal[peak_indices]

        mean_step_before = np.mean([self.raw_signal[i - 3:i] for i in peak_indices], axis=1)
        mean_step_after = np.mean([self.raw_signal[i + 1:i + 4] for i in peak_indices], axis=1)

        step_factor = (mean_step_before / mean_step_after) if self.is_reversed \
            else (mean_step_after / mean_step_before)

        plateau_width = 12
        if self.is_reversed:
            plateau_means = [np.mean(self.raw_signal[i - plateau_width - 1: i]) for i in peak_indices]
        else:
            plateau_means = [np.mean(self.raw_signal[i: i + plateau_width]) for i in peak_indices]

        plateau_means = np.array(plateau_means)

        all_data = np.stack(
            (
                peak_indices,
                np.full_like(peak_indices, self.thres_factor, dtype=float),
                np.round(peak_heights_deriv, 2),
                np.round(raw_heights, 2),
                np.round(mean_step_before, 2),
                np.round(mean_step_after, 2),
                np.round(step_factor, 2),
                np.round(plateau_means, 2),
            ),
            axis=-1,
        )

        mask = (
                (all_data[:, 2] >= self.min_height) &
                (all_data[:, 6] >= self.min_step_factor) &
                (all_data[:, 7] > 0)
        )

        final_peaks = all_data[mask]
        return final_peaks if final_peaks.size > 0 else None

    def _detect_peaks_single_scale(self, smooth_win, polyorder=1, delta=0.1):
        """Run savgol + rolling MAD peak detection at a single window size."""
        if self.is_reversed:
            sig_deriv = np.flipud(
                signal.savgol_filter(
                    np.flipud(self.raw_signal),
                    window_length=smooth_win,
                    polyorder=polyorder,
                    deriv=1,
                    delta=delta,
                )
            )
        else:
            sig_deriv = signal.savgol_filter(
                self.raw_signal,
                window_length=smooth_win,
                polyorder=polyorder,
                deriv=1,
                delta=delta,
            )

        mad_window = smooth_win * 3
        if len(sig_deriv) < mad_window:
            return np.array([]), np.array([])

        results_dict = self.rolling_robust_mad(
            sig_deriv, window_size=mad_window, sigma_cut=3.0, k_factor=1.4826, center=True
        )
        height_threshold = results_dict["upper_bound"]

        peaks, peaks_props = signal.find_peaks(
            sig_deriv,
            height=height_threshold,
        )

        if peaks.size == 0:
            return np.array([]), np.array([])

        return peaks, peaks_props["peak_heights"]

    @staticmethod
    def _deduplicate_peaks(peaks, heights, min_distance):
        """Merge peaks within min_distance, keeping the one with the highest derivative."""
        sort_idx = np.argsort(peaks)
        peaks = peaks[sort_idx]
        heights = heights[sort_idx]

        final_peaks = []
        final_heights = []
        i = 0
        while i < len(peaks):
            group_peaks = [peaks[i]]
            group_heights = [heights[i]]
            j = i + 1
            while j < len(peaks) and peaks[j] - peaks[i] <= min_distance:
                group_peaks.append(peaks[j])
                group_heights.append(heights[j])
                j += 1
            best = int(np.argmax(group_heights))
            final_peaks.append(group_peaks[best])
            final_heights.append(group_heights[best])
            i = j

        return np.array(final_peaks), np.array(final_heights)

    @staticmethod
    @njit(parallel=True, cache=True)
    def _compute_rolling_stats(windows):
        """
        Numba optimized core (unchanged).
        """
        n_windows, window_size = windows.shape
        medians = np.empty(n_windows, dtype=np.float64)
        mads = np.empty(n_windows, dtype=np.float64)

        for i in prange(n_windows):
            row = windows[i]
            med = np.median(row)
            abs_diffs = np.abs(row - med)
            mad = np.median(abs_diffs)
            medians[i] = med
            mads[i] = mad

        return medians, mads

    def rolling_robust_mad(self,
            arr: np.ndarray,
            window_size: int = 21,
            sigma_cut: float = 3.0,
            k_factor: float = 1.4826,
            center: bool = True
    ):
        """
        Calculates Rolling Robust MA with option to Center the window.

        Args:
        - center: If True, window is centered (no lag). If False, window is causal (lagged).
        """
        n = len(arr)
        if n < window_size:
            raise ValueError("Data length must be larger than window size.")

        # --- PADDING LOGIC HERE ---
        if center:
            # Pad both sides to center the window on the current point
            pad_before = (window_size - 1) // 2
            pad_after = window_size // 2
            padded_arr = np.pad(arr, (pad_before, pad_after), mode='edge')
        else:
            # Pad only the beginning (Causal/Real-time mode)
            pad_qty = window_size - 1
            padded_arr = np.pad(arr, (pad_qty, 0), mode='edge')

        # Create views and compute
        windows = sliding_window_view(padded_arr, window_shape=window_size)
        medians, mads = self._compute_rolling_stats(windows)

        # Calculate Bounds
        deviation = mads * k_factor * sigma_cut
        upper_bound = medians + deviation
        lower_bound = medians - deviation

        # Identify Outliers
        outlier_mask = (arr > upper_bound) | (arr < lower_bound)
        outlier_indices = np.where(outlier_mask)[0]

        return {
            "upper_bound": upper_bound,
            "lower_bound": lower_bound,
            "outlier_indices": outlier_indices,
            "medians": medians,
            "mads": mads
        }


    def call_signal_peaks_cwt(self, widths=None) -> np.array:
        """
        Alternative peak detection using Continuous Wavelet Transform.
        Inherently multi-scale — searches across all given widths simultaneously.
        Unused — kept for future experimentation.
        """
        if widths is None:
            widths = np.arange(1, 30)

        if self.is_reversed:
            sig = np.flipud(self.raw_signal)
        else:
            sig = self.raw_signal

        peaks = signal.find_peaks_cwt(sig, widths)

        if peaks.size == 0:
            return None

        if self.is_reversed:
            peaks = len(self.raw_signal) - 1 - peaks

        peak_indices = peaks.astype(int)
        max_idx = len(self.raw_signal)
        valid_mask = (peak_indices >= 3) & (peak_indices < max_idx - 15)
        peak_indices = peak_indices[valid_mask]

        if len(peak_indices) == 0:
            return None

        raw_heights = self.raw_signal[peak_indices]

        mean_step_before = np.mean([self.raw_signal[i - 3:i] for i in peak_indices], axis=1)
        mean_step_after = np.mean([self.raw_signal[i + 1:i + 4] for i in peak_indices], axis=1)

        step_factor = (mean_step_before / mean_step_after) if self.is_reversed \
            else (mean_step_after / mean_step_before)

        plateau_width = 12
        if self.is_reversed:
            plateau_means = [np.mean(self.raw_signal[i - plateau_width - 1: i]) for i in peak_indices]
        else:
            plateau_means = [np.mean(self.raw_signal[i: i + plateau_width]) for i in peak_indices]

        plateau_means = np.array(plateau_means)

        all_data = np.stack(
            (
                peak_indices,
                np.full_like(peak_indices, self.thres_factor, dtype=float),
                np.zeros_like(peak_indices, dtype=float),
                np.round(raw_heights, 2),
                np.round(mean_step_before, 2),
                np.round(mean_step_after, 2),
                np.round(step_factor, 2),
                np.round(plateau_means, 2),
            ),
            axis=-1,
        )

        mask = (
                (all_data[:, 6] >= self.min_step_factor) &
                (all_data[:, 7] > 0)
        )

        final_peaks = all_data[mask]
        return final_peaks if final_peaks.size > 0 else None

    def get_peaks_df(self):
        if self.peaks_arr is None:
            return None
        peaks_df = pd.DataFrame(
            data=self.peaks_arr,
            columns=[
                "peak_index",
                f"{self.prefix}_iqr_factor",  # TODO: Remove
                f"{self.prefix}_diff_height",
                f"{self.prefix}_height",
                f"{self.prefix}_upstream",
                f"{self.prefix}_downstream",
                f"{self.prefix}_step_factor",
                f"{self.prefix}_mean_plateau_height",
            ],
        )
        peaks_df["peak_index"] = peaks_df["peak_index"].astype(int)
        return peaks_df

    def get_bed_str(self, seqid: str):
        gff_df = self.get_peaks_df()
        if gff_df is None:
            return None
        columns = gff_df.columns.tolist()
        gff_df["seqid"] = seqid
        gff_df["start"] = gff_df["peak_index"] + 1
        gff_df["end"] = gff_df["peak_index"] + 1
        gff_df.drop(["peak_index"], inplace=True, axis=1)
        columns.remove("peak_index")
        gff_df["attributes"] = ""
        for indx in gff_df.index:
            gff_df.at[indx, "attributes"] = f"{self.prefix}_id={self.prefix}{indx}"
        for column in columns:
            gff_df["attributes"] += ";" + column + "=" + gff_df[column].astype(str)
        gff_df.drop(columns, inplace=True, axis=1)
        return gff_df.to_csv(index=False, sep="\t", header=False)

    def export_to_gff(
            self, out_path: str, seqid: str, strand: str, anno_source="_", anno_type="_"
    ):
        gff_columns = [
            "seqid",
            "source",
            "type",
            "start",
            "end",
            "score",
            "strand",
            "phase",
            "attributes",
        ]
        gff_df = self.get_peaks_df()
        if gff_df is None:
            logger.info("No peaks to export.")
            return

        non_gff_columns = [x for x in gff_df.columns.tolist() if x not in gff_columns]
        gff_df["peak_index"] += 1
        gff_df.rename({"peak_index": "start"}, inplace=True, axis=1)
        non_gff_columns.remove("peak_index")
        gff_df["end"] = gff_df["start"]
        gff_df["seqid"] = seqid
        gff_df = Helpers.get_gff_df(
            gff_df,
            anno_source=anno_source,
            anno_type=anno_type,
            strand=strand,
            new_id=True,
        )
        gff_df.to_csv(os.path.abspath(out_path), index=False, sep="\t", header=False)
        logger.info("GFF exported")