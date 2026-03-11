"""
Custom exceptions for GuaRNAtee.
Provides clear error messages and exception hierarchy.
"""


class GuaRNAteeException(Exception):
    """Base exception for all GuaRNAtee errors"""
    pass


class ConfigurationError(GuaRNAteeException):
    """Raised when configuration file is invalid or missing required keys"""
    pass


class ValidationError(GuaRNAteeException):
    """Raised when input validation fails"""
    pass


class FileValidationError(ValidationError):
    """Raised when file validation fails"""
    pass


class LibraryModeError(ValidationError):
    """Raised when library mode is invalid or incompatible"""
    pass


class ProcessingError(GuaRNAteeException):
    """Raised when processing pipeline encounters an error"""
    pass


class ExportError(GuaRNAteeException):
    """Raised when export operations fail"""
    pass


class InputDataError(GuaRNAteeException):
    """Raised when input data is empty or incompatible"""
    pass
