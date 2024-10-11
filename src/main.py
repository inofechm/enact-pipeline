"""
Created By  : ...
Created Date: DD/MM/YYYY
Description : ...
"""

import argparse
from utils.logging import get_logger


APP_NAME = 'MyProject'
LOGGER = get_logger(APP_NAME)


def dummy(dum):
    """Example function

    :param dum: Text to log.
    :type number: str
    :return: The entry text.
    :rtype: str
    """
    LOGGER.info(f'{dum} in progress')
    return dum

