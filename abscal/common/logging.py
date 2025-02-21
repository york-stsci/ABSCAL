"""Various logging functions for stis_tds_monitor.

Authors
-------

    - Brian York

Use
---

    Import the default logger and use it as follows:

    >>> from stis_tds_monitor.logging import DEFAULT_LOGGER
    >>> DEFAULT_LOGGER.info(f"Message")
"""

import logging.config
import logging
from pathlib import Path
import sys
import yaml

from .utils import get_data_file


def get_default_logger(log_name=__name__):
    """
    Create and configure a logger, then return that logger
    """
    log_config_file = get_data_file("abscal.common", "logging_config.yaml")
    with open(logging_config_file, "rt") as f:
        logging_config = yaml.safe_load(f.read())

    debug_log_file = Path.cwd() / logging_config['handlers']['debugfile']['filename']
    logging_config['handlers']['debugfile']['filename'] = debug_log_file

    history_log_file = Path.cwd() / logging_config['handlers']['processfile']['filename']
    logging_config['handlers']['processfile']['filename'] = history_log_file

    error_log_file = Path.cwd() / logging_config['handlers']['errorfile']['filename']
    logging_config['handlers']['errorfile']['filename'] = error_log_file

    logging.config.dictConfig(logging_config)
    logger = logging.getLogger(log_name)
    logger.debug("Logger has been configured.")

    return logger


DEFAULT_LOGGER = get_default_logger()
