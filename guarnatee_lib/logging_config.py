"""
Logging configuration for GuaRNAtee.
Handles integration between logging and tqdm progress bars.
"""
import logging
import sys
from typing import Optional
from colorama import Fore, Style


class ColoredFormatter(logging.Formatter):
    """
    Custom formatter that adds colors to log messages.
    """

    COLORS = {
        'DEBUG': Fore.CYAN,
        'INFO': Fore.YELLOW,
        'WARNING': Fore.YELLOW,
        'ERROR': Fore.RED,
        'CRITICAL': Fore.RED + Style.BRIGHT,
    }

    def format(self, record):
        log_color = self.COLORS.get(record.levelname, '')
        record.levelname = f"{log_color}{record.levelname}{Style.RESET_ALL}"
        record.msg = f"{log_color}{record.msg}{Style.RESET_ALL}"
        return super().format(record)


class TqdmLoggingHandler(logging.Handler):
    """
    Custom logging handler that works with tqdm progress bars.

    Ensures log messages don't interfere with progress bar display
    by using tqdm.write() instead of direct print.
    """

    def __init__(self, level=logging.NOTSET):
        super().__init__(level)

    def emit(self, record):
        try:
            msg = self.format(record)
            # Import tqdm here to avoid issues if not installed
            try:
                from tqdm import tqdm
                tqdm.write(msg, file=sys.stderr)
            except ImportError:
                # Fallback to regular print if tqdm not available
                print(msg, file=sys.stderr)
        except Exception:
            self.handleError(record)


def setup_logging(verbose: bool = False, use_tqdm_handler: bool = True) -> None:
    """
    Configure logging for GuaRNAtee with proper tqdm integration.

    Args:
        verbose: Enable debug-level logging
        use_tqdm_handler: Use TqdmLoggingHandler for tqdm compatibility
    """
    log_level = logging.DEBUG if verbose else logging.INFO

    # Create colored formatter
    formatter = ColoredFormatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # Get root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)

    # Remove existing handlers
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    # Add appropriate handler
    if use_tqdm_handler:
        handler = TqdmLoggingHandler()
    else:
        handler = logging.StreamHandler(sys.stderr)

    handler.setLevel(log_level)
    handler.setFormatter(formatter)
    root_logger.addHandler(handler)

    # Suppress verbose third-party libraries
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    logging.getLogger('PIL').setLevel(logging.WARNING)


def get_logger(name: str) -> logging.Logger:
    """
    Get a logger instance for a module.

    Args:
        name: Logger name (typically __name__)

    Returns:
        Logger instance
    """
    return logging.getLogger(name)
