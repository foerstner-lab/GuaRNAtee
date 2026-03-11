"""
Input validation for GuaRNAtee.
Validates and organizes wiggle file annotations and library modes.
"""
import logging
from typing import List, Dict, Tuple
import pandas as pd

from guarnatee_lib.constants import (
    LibraryModes, WigAnnotationColumns, LibraryTypes
)
from guarnatee_lib.exceptions import LibraryModeError, ValidationError

logger = logging.getLogger(__name__)


class WigfileValidator:
    """
    Validates wiggle file annotations and organizes by library mode.
    """

    @staticmethod
    def parse_and_validate(wig_annotations: List[str]) -> pd.DataFrame:
        """
        Parse wiggle file annotations from command line format.

        Expected format: file_path:condition:replicate:strand:lib_mode

        Args:
            wig_annotations: List of colon-separated annotation strings

        Returns:
            DataFrame with validated wiggle file annotations

        Raises:
            ValidationError: If annotations are malformed
            LibraryModeError: If library modes are invalid
        """
        logger.info(f"Parsing {len(wig_annotations)} wiggle file annotations")

        wigs_df = WigfileValidator._parse_annotations(wig_annotations)
        WigfileValidator._validate_library_modes(wigs_df)

        logger.info(f"Successfully validated wiggle files: "
                   f"{wigs_df['lib_mode'].value_counts().to_dict()}")

        return wigs_df

    @staticmethod
    def _parse_annotations(wig_annotations: List[str]) -> pd.DataFrame:
        """
        Parse annotation strings into DataFrame.

        Args:
            wig_annotations: List of colon-separated strings

        Returns:
            DataFrame with parsed annotations

        Raises:
            ValidationError: If format is incorrect
        """
        parsed_rows = []

        for idx, annotation in enumerate(wig_annotations):
            parts = annotation.split(":")

            if len(parts) != 5:
                raise ValidationError(
                    f"Invalid wiggle annotation format at index {idx}: '{annotation}'\n"
                    f"Expected format: file_path:condition:replicate:strand:lib_mode"
                )

            parsed_rows.append(parts)

        return pd.DataFrame(
            parsed_rows,
            columns=WigAnnotationColumns.ALL
        )

    @staticmethod
    def _validate_library_modes(wigs_df: pd.DataFrame) -> None:
        """
        Validate that all library modes are recognized.

        Args:
            wigs_df: DataFrame with wiggle annotations

        Raises:
            LibraryModeError: If invalid library modes found
        """
        valid_modes = LibraryModes.ALL
        invalid_mask = ~wigs_df[WigAnnotationColumns.LIB_MODE].isin(valid_modes)

        if invalid_mask.any():
            invalid_modes = wigs_df.loc[invalid_mask, WigAnnotationColumns.LIB_MODE].unique()
            raise LibraryModeError(
                f"Invalid library modes found: {list(invalid_modes)}\n"
                f"Allowed modes: {valid_modes}"
            )

    @staticmethod
    def organize_by_library_mode(wigs_df: pd.DataFrame) -> Dict[str, pd.DataFrame]:
        """
        Organize wiggle files by library mode into processing groups.

        Groups:
        - full_length: FL libraries (5' and 3' ends in same file)
        - paired: P1/P2 libraries (separate 5' and 3' end files)
        - dual_lib: 5E/3E libraries (separate enrichment)
        - differential: d5E libraries (differential 5' ends)

        Args:
            wigs_df: DataFrame with wiggle annotations

        Returns:
            Dictionary mapping library type to DataFrame subsets
        """
        library_groups = {}

        # Full-length libraries
        fl_wigs_df = wigs_df[
            wigs_df[WigAnnotationColumns.LIB_MODE] == LibraryModes.FULL_LENGTH
        ].copy()

        if not fl_wigs_df.empty:
            fl_wigs_df = WigfileValidator._prepare_full_length(fl_wigs_df)
            library_groups[LibraryTypes.FULL_LENGTH] = fl_wigs_df
            logger.info(f"Found {len(fl_wigs_df)} full-length library samples")

        # Paired libraries (P1/P2)
        pairs_wigs_df = wigs_df[
            wigs_df[WigAnnotationColumns.LIB_MODE].isin(LibraryModes.PAIRED_MODES)
        ].copy()

        if not pairs_wigs_df.empty:
            pairs_wigs_df = WigfileValidator._prepare_paired(pairs_wigs_df)
            if not pairs_wigs_df.empty:
                library_groups[LibraryTypes.PAIRED] = pairs_wigs_df
                logger.info(f"Found {len(pairs_wigs_df)} paired library samples")

        # End-enriched libraries (5E/3E/d5E)
        ends_wigs_df = wigs_df[
            wigs_df[WigAnnotationColumns.LIB_MODE].isin(LibraryModes.END_MODES)
        ].copy()

        if not ends_wigs_df.empty:
            dual_lib_df, diff_lib_df = WigfileValidator._prepare_ends(ends_wigs_df)

            if not dual_lib_df.empty:
                library_groups[LibraryTypes.DUAL_LIB] = dual_lib_df
                logger.info(f"Found {len(dual_lib_df)} dual-end library samples")

            if not diff_lib_df.empty:
                library_groups[LibraryTypes.DIFFERENTIAL] = diff_lib_df
                logger.info(f"Found {len(diff_lib_df)} differential library samples")

        return library_groups

    @staticmethod
    def _prepare_full_length(fl_wigs_df: pd.DataFrame) -> pd.DataFrame:
        """
        Prepare full-length library DataFrame.

        Full-length libraries use the same file for both 5' and 3' ends.

        Args:
            fl_wigs_df: DataFrame with FL library entries

        Returns:
            Prepared DataFrame with file_path_5e and file_path_3e columns
        """
        fl_wigs_df[WigAnnotationColumns.FILE_PATH_3E] = (
            fl_wigs_df[WigAnnotationColumns.FILE_PATH]
        )
        fl_wigs_df.rename(
            columns={WigAnnotationColumns.FILE_PATH: WigAnnotationColumns.FILE_PATH_5E},
            inplace=True
        )
        fl_wigs_df.drop(columns=[WigAnnotationColumns.LIB_MODE], inplace=True)

        return fl_wigs_df

    @staticmethod
    def _prepare_paired(pairs_wigs_df: pd.DataFrame) -> pd.DataFrame:
        """
        Prepare paired library DataFrame by merging P1 and P2.

        Args:
            pairs_wigs_df: DataFrame with P1/P2 library entries

        Returns:
            Merged DataFrame with file_path_5e and file_path_3e columns
        """
        merge_keys = [
            WigAnnotationColumns.CONDITION,
            WigAnnotationColumns.REPLICATE,
            WigAnnotationColumns.STRAND
        ]

        p1_df = pairs_wigs_df[
            pairs_wigs_df[WigAnnotationColumns.LIB_MODE] == LibraryModes.PAIRED_1
        ].copy()

        p2_df = pairs_wigs_df[
            pairs_wigs_df[WigAnnotationColumns.LIB_MODE] == LibraryModes.PAIRED_2
        ].copy()

        # Remove lib_mode before merge
        p1_df = p1_df[[WigAnnotationColumns.FILE_PATH] + merge_keys]
        p2_df = p2_df[[WigAnnotationColumns.FILE_PATH] + merge_keys]

        merged_df = pd.merge(
            p1_df,
            p2_df,
            on=merge_keys,
            how="inner",
            suffixes=('_5e', '_3e')
        )

        if merged_df.empty:
            logger.warning("No matching P1/P2 pairs found")

        return merged_df

    @staticmethod
    def _prepare_ends(
        ends_wigs_df: pd.DataFrame
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Prepare end-enriched libraries (5E/3E and d5E).

        Returns both dual-lib (5E+3E) and differential (d5E) DataFrames.

        Args:
            ends_wigs_df: DataFrame with 5E/3E/d5E library entries

        Returns:
            Tuple of (dual_lib_df, differential_df)
        """
        merge_keys = [
            WigAnnotationColumns.CONDITION,
            WigAnnotationColumns.REPLICATE,
            WigAnnotationColumns.STRAND
        ]

        # Extract differential 5' end libraries
        diff_ends_wigs_df = ends_wigs_df[
            ends_wigs_df[WigAnnotationColumns.LIB_MODE] == LibraryModes.DIFF_5_PRIME
        ].copy()

        # Remove lib_mode column
        diff_ends_wigs_df = diff_ends_wigs_df[
            [WigAnnotationColumns.FILE_PATH] + merge_keys
        ]

        # Merge 5E and 3E libraries
        e5_df = ends_wigs_df[
            ends_wigs_df[WigAnnotationColumns.LIB_MODE] == LibraryModes.END_5_PRIME
        ][[WigAnnotationColumns.FILE_PATH] + merge_keys]

        e3_df = ends_wigs_df[
            ends_wigs_df[WigAnnotationColumns.LIB_MODE] == LibraryModes.END_3_PRIME
        ][[WigAnnotationColumns.FILE_PATH] + merge_keys]

        merged_ends_df = pd.merge(
            e5_df,
            e3_df,
            on=merge_keys,
            how="inner",
            suffixes=('_5e', '_3e')
        )

        # Determine which samples have differential data
        combined_df = pd.merge(
            diff_ends_wigs_df,
            merged_ends_df,
            on=merge_keys,
            how="outer",
            indicator=True
        )

        # Samples with only 5E+3E (no d5E)
        dual_lib_df = combined_df[
            combined_df["_merge"] == "right_only"
        ].drop(columns=["_merge"]).copy()

        # Samples with d5E data
        diff_lib_df = combined_df[
            combined_df["_merge"] != "right_only"
        ].drop(columns=["_merge"]).copy()

        if not diff_lib_df.empty:
            diff_lib_df.rename(
                columns={WigAnnotationColumns.FILE_PATH: WigAnnotationColumns.FILE_PATH_D5E},
                inplace=True
            )

        return dual_lib_df, diff_lib_df
