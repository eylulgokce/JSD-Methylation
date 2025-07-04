
# -*- coding: utf-8 -*-
"""diversity.io
    ~~~~~~~~~~~~

    This module implements functions to load and transform data.
"""

import logging
import subprocess
from typing import Optional, List, Dict, Any, Union

import numpy as np
import pandas as pd

# Configure module logger
logger = logging.getLogger(__name__)


def bedcount_reader(bedcount: str, compression: Optional[str] = None, 
                   chunksize: int = 10000) -> pd.io.parsers.TextFileReader:
    """bedcount_reader returns a dataframe reader of the data.
    
    Args:
        bedcount: Path to bedcount file
        compression: Compression type (optional)
        chunksize: Number of rows per chunk (default: 10000)
        
    Returns:
        DataFrame reader object
        
    Raises:
        FileNotFoundError: If bedcount file doesn't exist
        ValueError: If parameters are invalid
        RuntimeError: If reader creation fails
    """
    logger.debug(f"Creating bedcount reader for file: {bedcount}")
    
    if not bedcount:
        raise ValueError("Bedcount file path cannot be empty")
    
    if chunksize <= 0:
        raise ValueError("Chunksize must be positive")
    
    try:
        reader = pd.read_table(
            bedcount, 
            compression=compression,
            chunksize=chunksize, 
            header=0,
            dtype={'#chrom': str, 'start': np.int64}
        )
        logger.debug(f"Bedcount reader created successfully with chunksize: {chunksize}")
        return reader
        
    except FileNotFoundError:
        logger.error(f"Bedcount file not found: {bedcount}")
        raise
    except Exception as e:
        logger.error(f"Error creating bedcount reader: {e}")
        raise RuntimeError(f"Failed to create bedcount reader: {e}")


def population_filter(metadata: str, subset: Optional[str] = None, 
                     relation: Optional[str] = None) -> Dict[str, Any]:
    """Read metadata and return the population and the quotient set.
    
    Args:
        metadata: Path to metadata file
        subset: Query string for subsetting (optional)
        relation: Column name for grouping (optional)
        
    Returns:
        Dictionary with 'reference' and 'qset' keys
        
    Raises:
        FileNotFoundError: If metadata file doesn't exist
        ValueError: If metadata format is invalid
        RuntimeError: If filtering fails
    """
    logger.debug(f"Filtering population from metadata: {metadata}")
    
    if not metadata:
        raise ValueError("Metadata file path cannot be empty")
    
    try:
        pop = {'reference': None, 'qset': None}
        
        # Read metadata
        meta = pd.read_table(metadata, header=0)
        logger.debug(f"Loaded metadata with {len(meta)} rows")
        
        # Validate metadata format
        if 'sample' not in meta.columns:
            raise ValueError("Metadata must contain 'sample' column")
        
        # Apply subset filter
        if subset is not None:
            logger.debug(f"Applying subset filter: {subset}")
            try:
                filtered_meta = meta.query(subset)
                pop['reference'] = list(filtered_meta['sample'])
                logger.debug(f"Subset filter resulted in {len(pop['reference'])} samples")
            except Exception as e:
                logger.error(f"Error applying subset filter '{subset}': {e}")
                raise ValueError(f"Invalid subset query: {subset}")
        else:
            pop['reference'] = list(meta['sample'])
            logger.debug(f"No subset filter, using all {len(pop['reference'])} samples")

        # Apply relation grouping
        if relation is not None:
            logger.debug(f"Applying relation grouping by: {relation}")
            
            if relation not in meta.columns:
                raise ValueError(f"Relation column '{relation}' not found in metadata")
            
            try:
                reference_meta = meta[meta['sample'].isin(pop['reference'])]
                group = reference_meta.groupby([relation])
                qset = []
                
                for name, df in group:
                    group_samples = list(df['sample'])
                    qset.append(group_samples)
                    logger.debug(f"Group '{name}': {len(group_samples)} samples")
                
                pop['qset'] = qset
                logger.debug(f"Created {len(qset)} groups from relation '{relation}'")
                
            except Exception as e:
                logger.error(f"Error creating relation groups: {e}")
                raise RuntimeError(f"Failed to create relation groups: {e}")

        logger.info(f"Population filter completed successfully")
        return pop
        
    except FileNotFoundError:
        logger.error(f"Metadata file not found: {metadata}")
        raise
    except Exception as e:
        logger.error(f"Error in population filter: {e}")
        raise RuntimeError(f"Population filtering failed: {e}")
