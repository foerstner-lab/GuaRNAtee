"""
Pipeline orchestration for GuaRNAtee.
Coordinates loading, validation, processing, and export.
"""
import logging
from typing import List, Dict, Any
from pathlib import Path

from guarnatee_lib.config_manager import ConfigManager
from guarnatee_lib.input_validator import WigfileValidator
from guarnatee_lib.processing_result import ProcessingResult
from guarnatee_lib.processors import PairedLibraryProcessor, DifferentialLibraryProcessor
from guarnatee_lib.export_manager import ExportManager
from guarnatee_lib.gff import GFF
from guarnatee_lib.fasta import Fasta
from guarnatee_lib.constants import LibraryTypes
from guarnatee_lib.exceptions import (
    GuaRNAteeException, InputDataError, ProcessingError
)

logger = logging.getLogger(__name__)


class PipelineOrchestrator:
    """
    Orchestrates the complete GuaRNAtee processing pipeline.

    Coordinates:
    1. Input file loading and validation
    2. Library mode organization
    3. Processing dispatch
    4. Results export
    """

    def __init__(
        self,
        gff_paths: List[str],
        fasta_paths: List[str],
        wig_annotations: List[str],
        config: ConfigManager,
        output_dir: str,
        threshold: int = 1
    ):
        """
        Initialize pipeline orchestrator.

        Args:
            gff_paths: Paths to GFF annotation files
            fasta_paths: Paths to FASTA genome files
            wig_annotations: Wiggle file annotations (colon-separated)
            config: Configuration manager
            output_dir: Output directory for results
            threshold: Peak calling threshold factor

        Raises:
            InputDataError: If inputs are invalid
        """
        self.gff_paths = gff_paths
        self.fasta_paths = fasta_paths
        self.wig_annotations = wig_annotations
        self.config = config
        self.output_dir = Path(output_dir).resolve()
        self.threshold = threshold

        # Will be populated during run()
        self.gff_obj = None
        self.fastas = None
        self.library_groups = {}
        self.results: List[ProcessingResult] = []

        logger.info("PipelineOrchestrator initialized")

    def run(self) -> List[ProcessingResult]:
        """
        Execute the complete processing pipeline.

        Returns:
            List of ProcessingResult objects

        Raises:
            GuaRNAteeException: If pipeline fails at any stage
        """
        logger.info("=" * 80)
        logger.info("Starting GuaRNAtee pipeline")
        logger.info("=" * 80)

        try:
            # Stage 1: Load inputs
            self._load_inputs()

            # Stage 2: Validate and organize
            self._organize_libraries()

            # Stage 3: Process all library groups
            self._process_all_libraries()

            # Stage 4: Export results
            self._export_results()

            logger.info("=" * 80)
            logger.info("Pipeline completed successfully")
            logger.info("=" * 80)

            return self.results

        except GuaRNAteeException:
            # Re-raise GuaRNAtee exceptions as-is
            raise
        except Exception as e:
            # Wrap unexpected exceptions
            raise ProcessingError(f"Unexpected pipeline error: {e}")

    def _load_inputs(self) -> None:
        """
        Load GFF, FASTA, and wiggle file inputs.

        Raises:
            InputDataError: If loading fails
        """
        logger.info("Loading input files...")

        try:
            # Load GFF annotations
            logger.info(f"Loading {len(self.gff_paths)} GFF file(s)")
            self.gff_obj = GFF(gff_paths=self.gff_paths)

            # Load FASTA sequences
            logger.info(f"Loading {len(self.fasta_paths)} FASTA file(s)")
            self.fastas = Fasta(fasta_paths=self.fasta_paths)

            logger.info("Input files loaded successfully")

        except Exception as e:
            raise InputDataError(f"Failed to load input files: {e}")

    def _organize_libraries(self) -> None:
        """
        Parse and organize wiggle files by library mode.

        Raises:
            InputDataError: If validation fails
        """
        logger.info("Organizing library files...")

        # Parse and validate wiggle annotations
        wigs_df = WigfileValidator.parse_and_validate(self.wig_annotations)

        # Organize by library mode
        self.library_groups = WigfileValidator.organize_by_library_mode(wigs_df)

        if not self.library_groups:
            raise InputDataError("No valid library groups found after organization")

        logger.info(f"Organized into {len(self.library_groups)} library group(s)")

    def _process_all_libraries(self) -> None:
        """
        Process all library groups.

        Dispatches to appropriate processor based on library type.
        """
        logger.info("Processing library groups...")

        config_dict = self.config.get_all()

        for library_type, wigs_df in self.library_groups.items():
            logger.info(f"Processing library type: {library_type}")

            try:
                if library_type == LibraryTypes.DIFFERENTIAL:
                    # Use differential processor
                    processor = DifferentialLibraryProcessor(
                        self.gff_obj,
                        self.fastas,
                        config_dict,
                        self.threshold
                    )
                else:
                    # Use paired processor (handles FL, Paired, dual_lib)
                    processor = PairedLibraryProcessor(
                        self.gff_obj,
                        self.fastas,
                        config_dict,
                        self.threshold
                    )

                result = processor.process(wigs_df, library_type)
                self.results.append(result)

                logger.info(
                    f"Completed {library_type}: "
                    f"{result.get_candidate_count()} candidates"
                )

            except Exception as e:
                logger.error(f"Failed to process {library_type}: {e}")
                raise

    def _export_results(self) -> None:
        """
        Export all results to output directory.

        Exports GFF, Excel, and statistics files.
        """
        logger.info("Exporting results...")

        if not self.results:
            logger.warning("No results to export")
            return

        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Initialize export manager
        seqid_groups = self.fastas.organism_seqid_groups
        exporter = ExportManager(str(self.output_dir), seqid_groups)

        # Export each result
        for result in self.results:
            if result.is_empty():
                logger.warning(f"Skipping empty result: {result.library_type}")
                continue

            try:
                # Export GFF files
                exporter.export_gff(result.candidates_df)

                # Export Excel files
                exporter.export_excel(result.candidates_df)

                # Export statistics
                exporter.export_statistics(result.stats_df)

            except Exception as e:
                logger.error(f"Failed to export {result.library_type}: {e}")
                raise

        logger.info(f"All results exported to: {self.output_dir}")

    def get_summary(self) -> Dict[str, Any]:
        """
        Get pipeline execution summary.

        Returns:
            Dictionary with summary statistics
        """
        total_candidates = sum(r.get_candidate_count() for r in self.results)
        total_samples = sum(r.sample_count for r in self.results)

        return {
            "total_candidates": total_candidates,
            "total_samples": total_samples,
            "library_types_processed": [r.library_type for r in self.results],
            "output_directory": str(self.output_dir),
        }
