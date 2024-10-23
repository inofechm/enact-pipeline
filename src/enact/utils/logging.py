"""
Created By  : ...
Created Date: DD/MM/YYYY
Description : ...
"""
import os
import logging


def get_logger(app_name, cache_dir):
    """Create and configure logger.

    :return: logger
    :rtype: Logger
    """
    # Clear log
    log_file = os.path.join(cache_dir, f"{app_name}.log")
    
    logger = logging.getLogger(app_name)
    logger.setLevel(logging.DEBUG)

    # Create file handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    
    # Create stream handler
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)
    stream_handler.setFormatter(formatter)

    # Add handlers to the logger
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    return logger
