"""
Processor for paired-end and full-length RNA-seq libraries.
Handles FL, P1/P2, and 5E/3E library modes.
"""
import logging
from typing import Dict, Any
import pandas as pd
import pybedtools as pybed

from guarnatee_lib.constants import (
    WigAnnotationColumns, StrandNotation, LibraryTypes,
    AnnotationSource, AnnotationType, StatisticsColumns, GFFColumns
)
from guarnatee_lib.exceptions import ProcessingError
from guarnatee_lib.processing_result import ProcessingResult
from guarnatee_lib.helpers import Helpers
from guarnatee_lib.rna_classifier import RNAClassifier
from guarnatee_lib.transcript_assembler import TranscriptAssembler
from guarnatee_lib.wiggle import Wiggle

logger = logging.getLogger(__name__)


class PairedLibraryProcessor:
    """
    Processes paired-end RNA-seq libraries.

    Handles:
    - Full-length libraries (FL): 5' and 3' ends in same file
    - Paired libraries (P1/P2): Separate 5' and 3' end files
    - Dual-end libraries (5E/3E): Separate end enrichments
    """

    def __init__(
        self,
        gff_obj,
        fastas,
        config_dict: Dict[str, Any],
        threshold: int
    ):
        """
        Initialize paired library processor.

        Args:
            gff_obj: GFF annotation object
            fastas: Fasta sequence object
            config_dict: Configuration dictionary
            threshold: Peak calling threshold factor
        """
        self.gff_obj = gff_obj
        self.fastas = fastas
        self.config_dict = config_dict
        self.threshold = threshold

    def process(
        self,
        wigs_df: pd.DataFrame,
        library_type: str
    ) -> ProcessingResult:
        """
        Process paired-end library samples.

        Args:
            wigs_df: DataFrame with wiggle file annotations
                     Must contain: file_path_5e, file_path_3e, condition, replicate, strand
            library_type: Type of library (full_length, Paired, dual_lib)

        Returns:
            ProcessingResult with candidates and statistics

        Raises:
            ProcessingError: If processing fails
        """
        logger.info(
            f"Processing {len(wigs_df)} {library_type} library samples"
        )

        try:
            export_df = pd.DataFrame()
            stats_df = pd.DataFrame()

            # Add file description for tracking
            wigs_df = self._add_file_descriptions(wigs_df)

            # Process each sample
            for idx in wigs_df.index:
                sample_export_df, sample_stats_df = self._process_single_sample(
                    wigs_df.loc[idx],
                    library_type
                )

                export_df = pd.concat([export_df, sample_export_df], ignore_index=True)
                stats_df = pd.concat([stats_df, sample_stats_df], ignore_index=True)

            # Post-process combined results
            if not export_df.empty:
                export_df = self._cluster_candidates(export_df)
                export_df = self._finalize_candidates(export_df)

            logger.info(
                f"Processed {library_type}: {len(export_df)} candidates "
                f"in {export_df.get('cluster_id', pd.Series()).nunique()} clusters"
            )

            return ProcessingResult(
                candidates_df=export_df,
                stats_df=stats_df,
                library_type=library_type,
                sample_count=len(wigs_df)
            )

        except Exception as e:
            raise ProcessingError(
                f"Failed to process {library_type} libraries: {e}"
            )

    def _add_file_descriptions(self, wigs_df: pd.DataFrame) -> pd.DataFrame:
        """
        Add descriptive labels for each sample.

        Args:
            wigs_df: Wiggle annotations DataFrame

        Returns:
            DataFrame with file_desc column added
        """
        wigs_df = wigs_df.copy()

        for idx in wigs_df.index:
            condition = wigs_df.at[idx, WigAnnotationColumns.CONDITION]
            replicate = wigs_df.at[idx, WigAnnotationColumns.REPLICATE]
            file_desc = f"{condition}_rep_{replicate}"
            wigs_df.at[idx, WigAnnotationColumns.FILE_DESC] = file_desc

        return wigs_df

    def _process_single_sample(
        self,
        sample_row: pd.Series,
        library_type: str
    ) -> tuple:
        """
        Process a single sample's wiggle files.

        Args:
            sample_row: Row from wiggle annotations DataFrame
            library_type: Type of library

        Returns:
            Tuple of (candidates_df, stats_df)
        """
        # Extract sample information
        five_prime_path = sample_row[WigAnnotationColumns.FILE_PATH_5E]
        three_prime_path = sample_row[WigAnnotationColumns.FILE_PATH_3E]
        strand = sample_row[WigAnnotationColumns.STRAND]
        file_desc = sample_row[WigAnnotationColumns.FILE_DESC]

        logger.debug(f"Processing sample: {file_desc}, strand: {strand}")

        # Convert strand notation
        strand_sign = self._convert_strand_notation(strand)

        # Call sRNA candidates
        candidates_df, stats_df = self._call_srnas(
            five_prime_path,
            three_prime_path
        )

        # Annotate candidates
        candidates_df[GFFColumns.STRAND] = strand_sign
        candidates_df["TSS_lib_type"] = LibraryTypes.TSS_PAIR
        candidates_df["condition"] = file_desc

        # Convert to GFF format
        candidates_df = Helpers.get_gff_df(
            candidates_df,
            anno_source=AnnotationSource.GUARNATEE,
            anno_type=AnnotationType.CANDIDATE,
            new_id=True
        )

        # Classify candidates
        candidates_df = Helpers.warp_non_gff_columns(
            RNAClassifier(
                self.gff_obj,
                candidates_df,
                self.fastas,
                self.config_dict
            ).classes
        )

        # Annotate statistics
        stats_df[StatisticsColumns.STRAND] = strand
        stats_df[StatisticsColumns.TSS_LIB_TYPE] = LibraryTypes.TSS_PAIR
        stats_df[StatisticsColumns.FILE_DESC] = file_desc

        return candidates_df, stats_df

    def _convert_strand_notation(self, strand: str) -> str:
        """
        Convert strand notation from f/r to +/-.

        Args:
            strand: Strand notation (f or r)

        Returns:
            Strand sign (+ or -)
        """
        if strand.lower() == StrandNotation.FORWARD:
            return StrandNotation.PLUS
        else:
            return StrandNotation.MINUS

    def _call_srnas(
        self,
        five_end_path: str,
        three_end_path: str
    ) -> tuple:
        """
        Call sRNA candidates from wiggle files.

        Args:
            five_end_path: Path to 5' end wiggle file
            three_end_path: Path to 3' end wiggle file

        Returns:
            Tuple of (candidates_df, stats_df)
        """
        # Load wiggle files
        if five_end_path == three_end_path:
            # Full-length: same file for both ends
            wiggle = Wiggle(five_end_path)
            srnas = TranscriptAssembler(wiggle, wiggle)
        else:
            # Separate files for 5' and 3' ends
            srnas = TranscriptAssembler(
                Wiggle(five_end_path),
                Wiggle(three_end_path)
            )

        # Call window-based sRNA candidates
        srnas.assemble_peaks(self.config_dict, thres_factor=self.threshold)

        return srnas.srna_candidates, srnas.log_df

    def _cluster_candidates(self, export_df: pd.DataFrame) -> pd.DataFrame:
        """
        Cluster overlapping candidates across samples.

        Args:
            export_df: DataFrame with all candidates

        Returns:
            DataFrame with cluster_id column added
        """
        # Wrap attributes before clustering
        export_df = Helpers.warp_non_gff_columns(export_df)

        # Convert to BED format and cluster
        export_bed = pybed.BedTool.from_dataframe(export_df).sort()
        clustered_bed = export_bed.cluster(s=True, d=0)

        # Convert back to DataFrame
        export_df = clustered_bed.to_dataframe(
            names=self.gff_obj.column_names + ["cluster_id"]
        )

        # Unwrap attributes
        export_df = Helpers.warp_non_gff_columns(export_df)

        return export_df

    def _finalize_candidates(self, export_df: pd.DataFrame) -> pd.DataFrame:
        """
        Finalize candidates DataFrame for export.

        Args:
            export_df: Clustered candidates DataFrame

        Returns:
            Sorted and indexed DataFrame ready for export
        """
        # Sort by genomic position
        export_df.sort_values(
            [GFFColumns.SEQID, GFFColumns.START, GFFColumns.END],
            inplace=True
        )
        export_df.reset_index(drop=True, inplace=True)

        # Ensure source is set
        export_df[GFFColumns.SOURCE] = AnnotationSource.GUARNATEE

        return export_df
