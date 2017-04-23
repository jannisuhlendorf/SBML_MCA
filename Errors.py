class NoncriticalError(Exception):
    """ error-class for errors that are non-critical (computations can continue) """
    pass


class CriticalError(Exception):
    """ error-class for errors that are critical (computations should be stopped) """
    pass
