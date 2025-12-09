class SmudgeplotError(Exception):
    """Base exception for smudgeplot errors."""
    pass
    
class InvalidCoverageDataError(SmudgeplotError):
    """Raised when coverage data is invalid."""
    pass