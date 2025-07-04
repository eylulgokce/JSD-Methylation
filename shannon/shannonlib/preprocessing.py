# -*- coding:utf-8 -*-
# preprocessing.py

"""Functions for preprocessing data.

This module provides utilities for data preprocessing including imputation
and filename generation.
"""

import logging
from typing import Union, List, Tuple, Optional

# Configure module logger
logger = logging.getLogger(__name__)


def impute(data: Optional[any] = None, method: str = 'pseudocount') -> Union[int, float]:
    """Imputes missing values of a data frame.
    
    Args:
        data: Data to impute (currently unused)
        method: Imputation method (default: 'pseudocount')
        
    Returns:
        Imputation value
        
    Raises:
        ValueError: If method is not supported
    """
    logger.debug(f"Imputing missing values using method: {method}")
    
    if method not in ['pseudocount']:
        raise ValueError(f"Unsupported imputation method: {method}")
    
    try:
        if method == 'pseudocount':
            # given by the parameters of a uninformative Dirichlet prior on the
            # probabilities
            impute_value = 1
            logger.debug(f"Using pseudocount imputation value: {impute_value}")
            return impute_value
    except Exception as e:
        logger.error(f"Error in imputation: {e}")
        raise RuntimeError(f"Imputation failed: {e}")


def groupname(by=None, name=None, fname=None):
    """Return filename of the subgroup.
    """

    # ensure that name is a tuple of strings
    name = tuple(str(key) for key in name)

    group = '_and_'.join('_'.join(items) for items in zip(by, name))

    old_suffix = fname.split('.')[-1]

    new_suffix = '.'.join([group, old_suffix])

    return fname.replace(old_suffix, new_suffix)