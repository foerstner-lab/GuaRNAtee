#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
"""
GuaRNAtee - Genome-wide RNA annotation tool for identification of transcripts.

A bioinformatics tool for identifying RNA candidates from RNA-seq data
using peak detection and annotation classification.

Author: Muhammad Elhossary <elhossary@zbmed.de>
        Konrad Förstner <konrad@foerstner.org>
Copyright: 2021 by Muhammad Elhossary <elhossary@zbmed.de>
License: ISC license
"""
import sys
import os
import argparse
import argcomplete
import logging
from pathlib import Path
from colorama import Fore

from guarnatee_lib.config_manager import ConfigManager
from guarnatee_lib.pipeline_orchestrator import PipelineOrchestrator
from guarnatee_lib.exceptions import GuaRNAteeException
from guarnatee_lib.logging_config import setup_logging

__author__ = ("Muhammad Elhossary <elhossary@zbmed.de> "
              "Konrad Förstner <konrad@foerstner.org> ")
__copyright__ = "2021 by Muhammad Elhossary <elhossary@zbmed.de>"
__license__ = "ISC license"
__email__ = "elhossary@zbmed.de"
__version__ = "1.0.0"
__maintainer__ = "Muhammad Elhossary"


def print_welcome_banner() -> None:
    """Print GuaRNAtee welcome banner."""
    banner = r"""
===============================================================================
||    ____                    ____    _   _      _       _                   ||
||   / ___|  _   _    __ _   |  _ \  | \ | |    / \     | |_    ___    ___   ||
||  | |  _  | | | |  / _` |  | |_) | |  \| |   / _ \    | __|  / _ \  / _ \  ||
||  | |_| | | |_| | | (_| |  |  _ <  | |\  |  / ___ \   | |_  |  __/ |  __/  ||
||   \____|  \__,_|  \__,_|  |_| \_\ |_| \_| /_/   \_\   \__|  \___|  \___|  ||
===============================================================================
"""
    print(Fore.RED + banner)
    print(Fore.WHITE + f"Version {__version__}\n")


def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments.

    Returns:
        Parsed arguments namespace
    """
    parser = argparse.ArgumentParser(
        description="GuaRNAtee - Genome-wide RNA annotation tool",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Version
    parser.add_argument(
        "-v", "--version",
        action="version",
        version=f"GuaRNAtee {__version__}"
    )

    # Required inputs
    parser.add_argument(
        "--gffs",
        required=True,
        type=str,
        nargs="+",
        help="GFF annotation files (space separated)"
    )

    parser.add_argument(
        "--fastas",
        required=True,
        type=str,
        nargs="+",
        help="FASTA genome files (space separated)"
    )

    parser.add_argument(
        "--wigs",
        required=True,
        type=str,
        nargs="+",
        help=(
            "Wiggle file annotations in format: "
            "file_path:condition:replicate:strand:lib_mode "
            "(space separated)"
        )
    )

    # Output
    parser.add_argument(
        "--out_dir",
        required=True,
        type=str,
        help="Output directory for results"
    )

    # Optional parameters
    parser.add_argument(
        "--threshold",
        type=int,
        default=1,
        help="Peak calling threshold factor (default: 1)"
    )

    parser.add_argument(
        "--config_file",
        default=None,
        type=str,
        help="Configuration file (default: config.cfg in script directory)"
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )

    # Enable autocompletion
    argcomplete.autocomplete(parser)

    args = parser.parse_args()

    # Set default config file if not provided
    if args.config_file is None:
        script_dir = Path(__file__).parent
        args.config_file = str(script_dir / "config.cfg")

    return args


def validate_inputs(args: argparse.Namespace) -> None:
    """
    Validate that input files exist.

    Args:
        args: Parsed command line arguments

    Raises:
        FileNotFoundError: If required files don't exist
    """
    # Check GFF files
    for gff_path in args.gffs:
        if not Path(gff_path).exists():
            raise FileNotFoundError(f"GFF file not found: {gff_path}")

    # Check FASTA files
    for fasta_path in args.fastas:
        if not Path(fasta_path).exists():
            raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    # Check wiggle files
    for wig_annotation in args.wigs:
        wig_path = wig_annotation.split(":")[0]
        if not Path(wig_path).exists():
            raise FileNotFoundError(f"Wiggle file not found: {wig_path}")

    # Check config file
    if not Path(args.config_file).exists():
        raise FileNotFoundError(f"Config file not found: {args.config_file}")


def main() -> int:
    """
    Main entry point for GuaRNAtee.

    Returns:
        Exit code (0 for success, 1 for failure)
    """
    # Print welcome banner
    print_welcome_banner()

    try:
        # Parse arguments
        args = parse_arguments()

        # Setup logging
        setup_logging(verbose=args.verbose)
        logger = logging.getLogger(__name__)

        # Validate inputs
        logger.info("Validating input files...")
        validate_inputs(args)

        # Load configuration
        logger.info(f"Loading configuration from: {args.config_file}")
        config = ConfigManager(args.config_file)

        # Create and run pipeline
        orchestrator = PipelineOrchestrator(
            gff_paths=args.gffs,
            fasta_paths=args.fastas,
            wig_annotations=args.wigs,
            config=config,
            output_dir=args.out_dir,
            threshold=args.threshold
        )

        # Execute pipeline
        results = orchestrator.run()

        # Print summary
        summary = orchestrator.get_summary()
        logger.info("\n" + "=" * 80)
        logger.info("Pipeline Summary:")
        logger.info(f"  Total candidates: {summary['total_candidates']}")
        logger.info(f"  Total samples: {summary['total_samples']}")
        logger.info(f"  Library types: {', '.join(summary['library_types_processed'])}")
        logger.info(f"  Output directory: {summary['output_directory']}")
        logger.info("=" * 80)

        print(Fore.GREEN + "\n✓ GuaRNAtee completed successfully!\n")
        print(Fore.WHITE + "")

        return 0

    except GuaRNAteeException as e:
        # Handle GuaRNAtee-specific exceptions
        logger = logging.getLogger(__name__)
        logger.error(f"Pipeline failed: {e}")
        print(Fore.RED + f"\n✗ Error: {e}\n")
        print(Fore.WHITE + "")
        return 1

    except KeyboardInterrupt:
        print(Fore.YELLOW + "\n\n✗ Pipeline interrupted by user\n")
        print(Fore.WHITE + "")
        return 1

    except Exception as e:
        # Handle unexpected exceptions
        logger = logging.getLogger(__name__)
        logger.exception("Unexpected error occurred")
        print(Fore.RED + f"\n✗ Unexpected error: {e}\n")
        print(Fore.WHITE + "")
        return 1


if __name__ == "__main__":
    sys.exit(main())
