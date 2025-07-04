# -*- coding:utf-8 -*-
# # estimators.py

"""Estimators of information-theoretic quantities.

This module provides functions for computing Shannon entropy and 
Jensen-Shannon divergence from count matrices.
"""

import logging
from typing import Union, Optional, Any
import numexpr as ne
import numpy as np
import pandas as pd

import shannonlib.constants as constant


logger = logging.getLogger(__name__)


def shannon_entropy(countmatrix: np.ndarray, axis: int = 1, 
                   method: str = 'plug-in') -> np.ndarray:
    """Shannon entropy (in nat) of the feature frequency profile.
    
    Args:
        countmatrix: Count matrix for computing entropy
        axis: Axis along which to compute entropy (default: 1)
        method: Method for entropy estimation (default: 'plug-in')
        
    Returns:
        Array of Shannon entropy values
        
    Raises:
        ValueError: If countmatrix is invalid or method is unsupported
        RuntimeError: If entropy computation fails
    """
    logger.debug(f"Computing Shannon entropy with method: {method}, axis: {axis}")
    
    if countmatrix is None:
        raise ValueError("Count matrix cannot be None")
    
    if method not in ['plug-in']:
        raise ValueError(f"Unsupported method: {method}. Only 'plug-in' is supported.")
    
    try:
        if method == 'plug-in':
            expression = ("sum(where(prob > 0, -prob * log(prob), 0), axis={})"
                         .format(axis))
            count_distribution = countmatrix.sum(axis)[..., np.newaxis]
            
            # Avoid division by zero
            if np.any(count_distribution == 0):
                logger.warning("Zero count distributions detected")
                count_distribution = np.where(count_distribution == 0, 1, count_distribution)
            
            prob = countmatrix / count_distribution
            result = ne.evaluate(expression)
            
            logger.debug(f"Shannon entropy computed successfully, shape: {result.shape}")
            return result
            
    except Exception as e:
        logger.error(f"Error computing Shannon entropy: {e}")
        raise RuntimeError(f"Shannon entropy computation failed: {e}")



def js_divergence(indata, weights=None):
    """
    Compute Jensen-Shannon divergence.
    
    Args:
        indata: Input data frame
        weights: Optional weights for averaging (currently unused)
        
    Returns:
        DataFrame with divergence results
        
    Raises:
        ValueError: If input data is invalid
        RuntimeError: If divergence computation fails

    """

    logger.debug("Starting JSD computation")
    logger.debug(f"Input shape: {indata.shape}")
    logger.debug(f"Input columns: {indata.columns}")

    if not isinstance(indata.columns, pd.MultiIndex):
        logger.warning("Input does not have a MultiIndex — attempting to reconstruct it")
        n_cols = indata.shape[1]
        if n_cols % 2 != 0:
            raise ValueError(f"Expected even number of columns (mC/C pairs), got {n_cols}")

        n_units = n_cols // 2
        units = [f'unit_{i+1}' for i in range(n_units)]
        features = ['mC', 'C']
        indata.columns = pd.MultiIndex.from_product(
            [units, features], names=['sampling_unit', 'feature']
        )
        logger.debug(f"Reconstructed MultiIndex: {indata.columns}")

    # Compute total methylation count per sample
    count_per_unit = indata.sum(axis=1, level='sampling_unit')
    samplesize = count_per_unit.notnull().sum(axis=1)

    # Apply QC filters
    min_samplesize = 2
    min_count = 3
    count_filter = (count_per_unit >= min_count).any(axis=1)
    samplesize_filter = (samplesize >= min_samplesize)
    combined_filter = (count_filter & samplesize_filter)
    data = indata[combined_filter]

    logger.debug(f"Rows before filtering: {len(indata)}")
    logger.debug(f"Rows after filtering: {len(data)}")

    if data.empty:
        logger.warning("No data passed QC filtering — returning empty DataFrame")
        return data

    # Entropy computation
    try:
        data_unit = count_per_unit[combined_filter]
        data_feature = (data.sum(axis=1, level='feature').astype(np.int32))
        counts = data.values.reshape(
            data_unit.shape[0],
            data_unit.shape[1],
            data_feature.shape[1]
        )

        mix_entropy = shannon_entropy(data_feature.values)
        avg_entropy = np.average(
            shannon_entropy(counts, axis=2),
            weights=data_unit.fillna(0),
            axis=1
        )

        div = data_feature.copy()
        div.insert(0, 'JSD_bit_', constant.LOG2E * (mix_entropy - avg_entropy))
        div.insert(1, 'sample size', samplesize[combined_filter])
        div.insert(2, 'HMIX_bit_', constant.LOG2E * mix_entropy)

        logger.debug(f"JSD computation completed for {len(div)} rows")
        return div

    except Exception as e:
        logger.error(f"JSD computation failed: {e}")
        raise RuntimeError(f"JSD divergence computation failed: {e}")
