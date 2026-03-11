"""
Configuration management for GuaRNAtee.
Handles loading, validation, and type conversion of configuration files.
"""
import os
import configparser
from pathlib import Path
from typing import Dict, Any, Union
import logging

from guarnatee_lib.constants import ConfigKeys
from guarnatee_lib.exceptions import ConfigurationError

logger = logging.getLogger(__name__)


class ConfigManager:
    """
    Manages configuration file loading and validation.

    Converts string values from config file to appropriate types
    and provides type-safe access to configuration values.
    """

    # Required configuration keys
    REQUIRED_KEYS = [
        ConfigKeys.MIN_HEIGHT,
        ConfigKeys.MIN_LEN,
        ConfigKeys.MAX_LEN,
        ConfigKeys.MIN_STEP_FACTOR,
        ConfigKeys.MIN_DISTANCE,
    ]

    # Type conversion mapping
    TYPE_CONVERTERS = {
        ConfigKeys.MIN_HEIGHT: float,
        ConfigKeys.MIN_LEN: int,
        ConfigKeys.MAX_LEN: int,
        ConfigKeys.READ_LENGTH: int,
        ConfigKeys.MIN_STEP_FACTOR: float,
        ConfigKeys.MIN_ORF_UD_FRAG_RATIO: float,
        ConfigKeys.MIN_DISTANCE: int,
        ConfigKeys.MAX_TSS_LEN: int,
        ConfigKeys.MAX_TTS_LEN: int,
        ConfigKeys.MAX_OUTBOUND_TSS_TOLERANCE: int,
        ConfigKeys.MAX_OUTBOUND_TTS_TOLERANCE: int,
        ConfigKeys.MIN_MFE: float,
        ConfigKeys.DETAILED_OUTPUT: lambda x: x.lower() in ('true', '1', 'yes'),
        ConfigKeys.MERGE_SIMILARITY_RATIO: float,
    }

    def __init__(self, config_file_path: Union[str, Path]):
        """
        Initialize ConfigManager with path to configuration file.

        Args:
            config_file_path: Path to configuration file

        Raises:
            ConfigurationError: If config file doesn't exist
        """
        self.config_path = Path(config_file_path).resolve()
        if not self.config_path.exists():
            raise ConfigurationError(
                f"Configuration file not found: {self.config_path}"
            )

        self.config_dict: Dict[str, Any] = {}
        self._load_config()

    def _load_config(self) -> None:
        """
        Load and parse configuration file.

        Raises:
            ConfigurationError: If config file is malformed
        """
        logger.info(f"Loading configuration from: {self.config_path}")

        try:
            parser = configparser.ConfigParser(strict=True)
            parser.read(self.config_path)

            # Merge all sections into a single dictionary
            for section in parser.sections():
                self.config_dict.update(dict(parser.items(section)))

        except configparser.Error as e:
            raise ConfigurationError(
                f"Failed to parse configuration file: {e}"
            )

        self._validate_config()
        self._convert_types()

        logger.info(f"Configuration loaded successfully with {len(self.config_dict)} parameters")

    def _validate_config(self) -> None:
        """
        Validate that all required keys are present.

        Raises:
            ConfigurationError: If required keys are missing
        """
        missing_keys = [
            key for key in self.REQUIRED_KEYS
            if key not in self.config_dict
        ]

        if missing_keys:
            raise ConfigurationError(
                f"Missing required configuration keys: {', '.join(missing_keys)}"
            )

    def _convert_types(self) -> None:
        """
        Convert configuration values from strings to appropriate types.

        Raises:
            ConfigurationError: If type conversion fails
        """
        for key, value in self.config_dict.items():
            if key in self.TYPE_CONVERTERS:
                converter = self.TYPE_CONVERTERS[key]
                try:
                    self.config_dict[key] = converter(value)
                except (ValueError, TypeError) as e:
                    raise ConfigurationError(
                        f"Failed to convert config key '{key}' with value '{value}': {e}"
                    )

    def get(self, key: str, default: Any = None) -> Any:
        """
        Get configuration value with optional default.

        Args:
            key: Configuration key
            default: Default value if key not found

        Returns:
            Configuration value or default
        """
        return self.config_dict.get(key, default)

    def get_all(self) -> Dict[str, Any]:
        """
        Get all configuration values as dictionary.

        Returns:
            Dictionary of all configuration values
        """
        return self.config_dict.copy()

    def __getitem__(self, key: str) -> Any:
        """
        Dictionary-style access to configuration values.

        Args:
            key: Configuration key

        Returns:
            Configuration value

        Raises:
            KeyError: If key not found
        """
        if key not in self.config_dict:
            raise KeyError(f"Configuration key not found: {key}")
        return self.config_dict[key]

    def __contains__(self, key: str) -> bool:
        """Check if configuration key exists"""
        return key in self.config_dict
