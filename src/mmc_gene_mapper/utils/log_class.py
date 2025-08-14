import warnings


class StdoutLog(object):
    """
    Class to mimic log interface that actually just
    writes to stdout
    """
    def info(self, msg, to_stdout=True):
        print(f"===={msg}")

    def warn(self, msg):
        warnings.warn(msg)
