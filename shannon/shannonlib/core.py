# -*- coding:utf-8 -*-
# core.py

"""Core functions.

This module contains the main functions for computing population divergence
using information-theoretic measures.
"""

import os
import logging
from typing import Optional, Dict, List, Any
import pandas as pd

import shannonlib.estimators as est
import shannonlib.gpf_utils as gpf

logger = logging.getLogger(__name__)


def divergence(sample, chrom=None, data_columns=None, outfile=None, chunksize=None):
    """Computes within-group divergence for population.
    
    Args:
        sample: Dictionary containing 'url' and 'label' keys
        chrom: Chromosome identifier (optional)
        data_columns: List of data columns to process (optional)
        outfile: Output file path (optional)
        chunksize: Expected number of sites per chunk (optional)
        
    Returns:
        None
        
    Raises:
        ValueError: If sample dictionary is invalid
        FileNotFoundError: If input files cannot be found
        IOError: If output file cannot be written
    """
    logger.info(f"Starting divergence computation for chromosome: {chrom}")
    
    if isinstance(sample, pd.DataFrame):
        if not {'url', 'label'}.issubset(sample.columns):
            raise ValueError("Sample DataFrame must contain 'url' and 'label' columns")
        sample = {
            "url": list(sample["url"]),
            "label": list(sample["label"])
        }
    elif isinstance(sample, dict):
        if 'url' not in sample or 'label' not in sample:
            raise ValueError("Sample dictionary must contain 'url' and 'label' keys")
    else:
        raise ValueError("Sample must be a dictionary or a DataFrame")
    
    try:
        # Get regions for processing
        logger.debug(f"Getting regions for sample: {sample}")
        regions_result = gpf.get_regions(
            sample['url'], chrom=chrom, exp_numsites=chunksize)
        
        if not regions_result:
            logger.warning("No regions found, skipping divergence computation")
            return None
            
        regions_pct, regions = regions_result
        regions = list(regions)  
        logger.info(f"Found {len(regions)} regions to process")
        
        # Convert regions back to list since it was consumed
        regions = list(regions)
        
        # Get data for the regions
        logger.debug("Retrieving data for regions")
        regions_data = gpf.get_data(
            sample['url'], 
            labels=sample['label'],
            data_columns=data_columns, 
            regions=regions
        )

        processed_count = 0
        skipped_empty = 0
        skipped_quality = 0
        
        for progress, data in zip(regions_pct, regions_data):
            try:
                if data.empty:
                    logger.debug(f"Skipping empty region at {progress}%")
                    print('...{:>5} % (skipped empty region)'.format(progress))
                    skipped_empty += 1
                    continue

                # Compute divergence
                logger.debug(f"Computing JS divergence for region at {progress}%")
                div = est.js_divergence(data)

                if div.empty:
                    logger.debug(f"Skipping low-quality region at {progress}%")
                    print('...{:>5} % (skipped low-quality region)'.format(progress))
                    skipped_quality += 1
                    continue

                # output file
                if outfile:
                    try:
                        if not os.path.isfile(outfile):
                            header = True
                        elif os.stat(outfile).st_size == 0:
                            header = True
                        else:
                            header = False

                        (div
                         .round({'JSD_bit_': 3, 'HMIX_bit_': 3})
                         .to_csv(outfile, header=header, sep='\t', index=True, mode='a'))
                        
                        logger.debug(f"Results written to {outfile}")
                        
                    except IOError as e:
                        logger.error(f"Failed to write to output file {outfile}: {e}")
                        raise

                print('...{:>5} %'.format(progress))
                processed_count += 1
                
            except Exception as e:
                logger.error(f"Error processing region at {progress}%: {e}")
                # continue processing other regions instead of failing completely
                continue

        logger.info(f"Divergence computation completed. Processed: {processed_count}, "
                   f"Skipped empty: {skipped_empty}, Skipped low-quality: {skipped_quality}")
        
    except Exception as e:
        logger.error(f"Fatal error in divergence computation: {e}")
        raise

    return None