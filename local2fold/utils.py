import logging


class TqdmHandler(logging.StreamHandler):
    """https://stackoverflow.com/a/38895482/3549270"""

    def __init__(self):
        logging.StreamHandler.__init__(self)

    def emit(self, record):
        # We need the native tqdm here
        from tqdm import tqdm

        msg = self.format(record)
        tqdm.write(msg)
