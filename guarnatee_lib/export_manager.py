"""
Export management for GuaRNAtee.
Handles all file export operations including GFF, Excel, and statistics.
"""
import os
import logging
from pathlib import Path
from typing import Dict, List
import numpy as np
import pandas as pd

from guarnatee_lib.constants import (
    ExportColumns, AnnotationClasses, FileExtensions,
    StatisticsColumns, GFFColumns
)
from guarnatee_lib.exceptions import ExportError
from guarnatee_lib.helpers import Helpers

logger = logging.getLogger(__name__)


class ExportManager:
    """
    Manages all export operations for GuaRNAtee results.

    Exports GFF files, Excel workbooks, and statistics files
    organized by organism and annotation class.
    """

    def __init__(self, output_dir: str, seqid_groups: Dict[str, List[str]]):
        """
        Initialize ExportManager.

        Args:
            output_dir: Directory for output files (already created by orchestrator)
            seqid_groups: Mapping of organism names to sequence IDs
        """
        self.output_dir = Path(output_dir).resolve()
        self.seqid_groups = seqid_groups

        # Output directory should already exist from orchestrator
        # Ensure it exists in case ExportManager is used standalone
        self.output_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"ExportManager initialized with output directory: {self.output_dir}")

    def export_gff(self, export_df: pd.DataFrame, suffix: str = "candidates") -> None:
        """
        Export GFF files, one per organism.

        Args:
            export_df: DataFrame with GFF-formatted candidates
            suffix: Suffix for output filenames (default: "candidates")

        Raises:
            ExportError: If export fails
        """
        if export_df.empty:
            logger.warning("No data to export to GFF")
            return

        logger.info(f"Exporting GFF files for {len(self.seqid_groups)} organisms")

        try:
            for organism, seqids in self.seqid_groups.items():
                organism_df = export_df[
                    export_df[GFFColumns.SEQID].isin(seqids)
                ].copy()

                if organism_df.empty:
                    logger.warning(f"No candidates for organism: {organism}")
                    continue

                output_path = self.output_dir / f"{organism}_{suffix}{FileExtensions.GFF}"

                organism_df.to_csv(
                    output_path,
                    index=False,
                    sep="\t",
                    header=False
                )

                logger.info(f"Exported {len(organism_df)} candidates to {output_path.name}")

        except Exception as e:
            raise ExportError(f"Failed to export GFF files: {e}")

    def export_excel(
        self,
        export_df: pd.DataFrame,
        suffix: str = "candidates"
    ) -> None:
        """
        Export Excel files with candidates grouped by annotation class.

        Creates one workbook per organism with sheets for:
        - intergenic: Intergenic and ncRNA candidates
        - ORF_int: Candidates internal to ORFs
        - others: Antisense and cross-feature candidates

        Args:
            export_df: DataFrame with candidates
            suffix: Suffix for output filenames (default: "candidates")

        Raises:
            ExportError: If export fails
        """
        if export_df.empty:
            logger.warning("No data to export to Excel")
            return

        logger.info(f"Exporting Excel files for {len(self.seqid_groups)} organisms")

        try:
            # Prepare DataFrame for Excel export
            excel_df = self._prepare_for_excel(export_df)

            # Group candidates by annotation class
            class_groups = self._create_classification_groups(excel_df)

            # Export one workbook per organism
            for organism, seqids in self.seqid_groups.items():
                self._export_organism_workbook(
                    excel_df,
                    class_groups,
                    organism,
                    seqids,
                    suffix
                )

        except Exception as e:
            raise ExportError(f"Failed to export Excel files: {e}")

    def _prepare_for_excel(self, export_df: pd.DataFrame) -> pd.DataFrame:
        """
        Prepare DataFrame for Excel export.

        Expands attributes column and removes unnecessary columns.

        Args:
            export_df: Raw export DataFrame

        Returns:
            Cleaned DataFrame ready for Excel export
        """
        # Expand attributes into separate columns
        excel_df = Helpers.expand_attributes_to_columns(export_df.copy())

        # Drop unnecessary columns
        drop_cols = [
            col for col in ExportColumns.DROP_COLS
            if col in excel_df.columns
        ]
        if drop_cols:
            excel_df.drop(columns=drop_cols, inplace=True)

        # Rename columns (replace underscores with spaces)
        rename_cols = {col: col.replace("_", " ") for col in excel_df.columns}
        excel_df.rename(columns=rename_cols, inplace=True)

        return excel_df

    def _create_classification_groups(
        self,
        export_df: pd.DataFrame
    ) -> Dict[str, List[str]]:
        """
        Group annotation classes into categories.

        Groups:
        - intergenic: ncRNA and intergenic candidates
        - ORF_int: Candidates internal to ORFs
        - others: Antisense and cross-feature candidates

        Args:
            export_df: DataFrame with annotation class column

        Returns:
            Dictionary mapping group names to lists of annotation classes
        """
        annotation_col = "annotation class"  # Space added by Excel prep

        if annotation_col not in export_df.columns:
            logger.warning("No annotation_class column found, using default groups")
            return {
                AnnotationClasses.INTERGENIC: [],
                AnnotationClasses.ORF_INT: [],
                AnnotationClasses.OTHERS: []
            }

        groups = {
            AnnotationClasses.INTERGENIC: [],
            AnnotationClasses.ORF_INT: [],
            AnnotationClasses.OTHERS: []
        }

        for cls in export_df[annotation_col].unique():
            cls_str = str(cls)

            # Intergenic: ncRNA or intergenic
            if "ncRNA" in cls_str or "intergenic" in cls_str:
                # Exclude cross/antisense with ncRNA
                if not (("cross" in cls_str and "ncRNA" in cls_str) or
                        "antisense_to" in cls_str):
                    groups[AnnotationClasses.INTERGENIC].append(cls)
                    continue

            # ORF internal
            if "ORF_int" in cls_str:
                groups[AnnotationClasses.ORF_INT].append(cls)
                continue

            # Everything else
            groups[AnnotationClasses.OTHERS].append(cls)

        return groups

    def _export_organism_workbook(
        self,
        excel_df: pd.DataFrame,
        class_groups: Dict[str, List[str]],
        organism: str,
        seqids: List[str],
        suffix: str
    ) -> None:
        """
        Export Excel workbook for a single organism.

        Args:
            excel_df: Prepared DataFrame
            class_groups: Classification groups
            organism: Organism name
            seqids: List of sequence IDs for this organism
            suffix: Filename suffix
        """
        output_path = self.output_dir / f"{organism}_{suffix}{FileExtensions.XLSX}"

        seqid_col = GFFColumns.SEQID.replace("_", " ")  # Column renamed in prep
        annotation_col = "annotation class"

        organism_df = excel_df[excel_df[seqid_col].isin(seqids)].copy()

        if organism_df.empty:
            logger.warning(f"No candidates for organism {organism}, skipping Excel export")
            return

        with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
            sheets_written = 0

            for group_name, classes in class_groups.items():
                if not classes:
                    continue

                group_df = organism_df[
                    organism_df[annotation_col].isin(classes)
                ].copy()

                if group_df.empty:
                    continue

                # Clean up empty columns
                # Use mask instead of replace to avoid FutureWarning
                group_df = group_df.mask(group_df == "", np.nan)
                group_df = group_df.infer_objects(copy=False)
                group_df.dropna(how='all', axis=1, inplace=True)
                group_df.reset_index(drop=True, inplace=True)

                # Write to sheet
                group_df.to_excel(
                    excel_writer=writer,
                    sheet_name=group_name,
                    index=True,
                    header=True,
                    na_rep="",
                    index_label="index"
                )

                sheets_written += 1

        logger.info(
            f"Exported {len(organism_df)} candidates "
            f"to {output_path.name} ({sheets_written} sheets)"
        )

    def export_statistics(
        self,
        stats_df: pd.DataFrame,
        filename: str = "stats"
    ) -> None:
        """
        Export statistics TSV file.

        Args:
            stats_df: DataFrame with statistics
            filename: Output filename (without extension)

        Raises:
            ExportError: If export fails
        """
        if stats_df.empty:
            logger.warning("No statistics to export")
            return

        try:
            # Add organism labels
            stats_df = self._add_organism_labels(stats_df)

            # Aggregate by organism and library type
            aggregated_stats = stats_df.groupby(
                [
                    StatisticsColumns.ORGANISM,
                    StatisticsColumns.FILE_DESC,
                    StatisticsColumns.TSS_LIB_TYPE
                ],
                as_index=False
            ).agg({
                StatisticsColumns.TSS_LIB_WINDOWS_COUNT: "sum",
                StatisticsColumns.TSS_LIB_PEAKS_COUNT: "sum",
                StatisticsColumns.TTS_LIB_WINDOWS_COUNT: "sum",
                StatisticsColumns.TTS_LIB_PEAKS_COUNT: "sum",
                StatisticsColumns.PEAKS_CONNECTIONS_COUNT: "sum"
            })

            output_path = self.output_dir / f"{filename}{FileExtensions.TSV}"

            aggregated_stats.to_csv(
                output_path,
                sep='\t',
                index=False
            )

            logger.info(f"Exported statistics to {output_path.name}")

        except Exception as e:
            raise ExportError(f"Failed to export statistics: {e}")

    def _add_organism_labels(self, stats_df: pd.DataFrame) -> pd.DataFrame:
        """
        Add organism labels to statistics DataFrame.

        Args:
            stats_df: Statistics DataFrame with seqid column

        Returns:
            DataFrame with Organism column added
        """
        stats_df = stats_df.copy()

        if StatisticsColumns.ORGANISM not in stats_df.columns:
            stats_df[StatisticsColumns.ORGANISM] = ""

        for idx in stats_df.index:
            seqid = stats_df.at[idx, StatisticsColumns.SEQID]

            for organism, seqids in self.seqid_groups.items():
                if seqid in seqids:
                    stats_df.at[idx, StatisticsColumns.ORGANISM] = organism
                    break

        return stats_df


class ORFStatsExporter:
    """
    Exports ORF-internal candidate statistics.

    Separate class for specialized ORF statistics export.
    """

    def __init__(self, output_dir: str, seqid_groups: Dict[str, List[str]]):
        """
        Initialize ORF statistics exporter.

        Args:
            output_dir: Directory for output files
            seqid_groups: Mapping of organism names to sequence IDs
        """
        self.output_dir = Path(output_dir).resolve()
        self.seqid_groups = seqid_groups

    def export(self, export_df: pd.DataFrame) -> None:
        """
        Export ORF-internal statistics to Excel.

        Args:
            export_df: DataFrame with expanded attributes

        Raises:
            ExportError: If export fails
        """
        if export_df.empty or "sub_class" not in export_df.columns:
            logger.warning("No ORF-internal data to export")
            return

        try:
            # Filter to ORF-internal candidates only
            orf_int_df = export_df[export_df["sub_class"] != ""].copy()

            if orf_int_df.empty:
                logger.info("No ORF-internal candidates found")
                return

            # Add organism labels
            orf_int_df["Specie"] = ""
            for organism, seqids in self.seqid_groups.items():
                orf_int_df.loc[
                    orf_int_df["seqid"].isin(seqids),
                    "Specie"
                ] = organism

            # Count by species, library type, condition, and sub-class
            stats_df = orf_int_df.value_counts(
                subset=["Specie", "lib_type", "condition", "sub_class"]
            ).reset_index()

            stats_df.columns = ["Specie", "lib_type", "condition", "sub_class", "count"]

            # Export to Excel with one sheet per organism
            output_path = self.output_dir / f"ORF_int_stats{FileExtensions.XLSX}"

            with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
                for organism in stats_df["Specie"].unique():
                    organism_stats = stats_df[stats_df["Specie"] == organism].copy()

                    organism_stats.to_excel(
                        excel_writer=writer,
                        sheet_name=organism,
                        index=False,
                        header=True,
                        na_rep=""
                    )

            logger.info(f"Exported ORF-internal statistics to {output_path.name}")

        except Exception as e:
            raise ExportError(f"Failed to export ORF statistics: {e}")
