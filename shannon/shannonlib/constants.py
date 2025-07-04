#constants.py

"""Useful constants.

This module contains mathematical constants used throughout the shannonlib package.
"""

import math
import logging
import numpy as np

logger = logging.getLogger(__name__)

try:
    LOG2E = np.log2(math.e)
    logger.debug(f"LOG2E constant initialized: {LOG2E}")
except Exception as e:
    logger.error(f"Failed to initialize LOG2E constant: {e}")
    raise