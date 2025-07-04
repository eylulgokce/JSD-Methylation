
# -*- coding:utf-8 -*-
# gpf_utils.py

"""Utility functions for handling genome position files.

This module provides utilities for working with tabix-indexed genome position
files, including region extraction and data merging.
"""

import logging
import math
import subprocess
from typing import List, Tuple, Optional, Union, Generator, Any

import numpy as np
import pandas as pd


class InputMismatchError(Exception):
    """Raised when input parameters don't match expected format."""
    pass


class MissingInputError(Exception):
    """Raised when required input parameters are missing."""
    pass


# Configure logging
logging.basicConfig(
    format="=== %(levelname)s === %(asctime)s === %(message)s",
    level=logging.DEBUG, 
    datefmt='%Y-%m-%d %H:%M:%S'
)

logger = logging.getLogger(__name__)


def get_regions(tabixfiles: Union[str, List[str]], chrom: Optional[str] = None, 
                exp_numsites: float = 1e3) -> Union[Tuple[List[float], zip], bool]:
    """Get stepsize and list of regions for tabix-indexed files.
    
    Args:
        tabixfiles: Path(s) to tabix-indexed files
        chrom: Chromosome identifier (optional)
        exp_numsites: Expected number of sites per region (default: 1000)
        
    Returns:
        Tuple of (progress_percentages, regions) or False if no data
        
    Raises:
        ValueError: If input parameters are invalid
        RuntimeError: If region computation fails
    """
    logger.debug(f"Getting regions for chromosome: {chrom}")
    
    if not tabixfiles:
        raise ValueError("Tabix files cannot be empty")
    
    try:
        sup_position = supremum_position(tabixfiles, chrom)

        if sup_position is None:
            logger.info("Skipping because chromosome is missing.")
            return False

        sup_numsites = supremum_numsites(tabixfiles, chrom)

        if sup_numsites is None or sup_numsites == 0:
            logger.info("Skipping because there are no entries.")
            return False

        # Calculate step size
        step = math.ceil(sup_position / sup_numsites * exp_numsites)

        if step < sup_position:
            stepsize = step
        else:
            stepsize = sup_position

        logger.debug(f"Computed stepsize: {stepsize}")

        # Generate position ranges
        pos_start = list(range(0, sup_position, stepsize + 1))
        pos_end = list(range(stepsize, sup_position, stepsize + 1)) + [sup_position]

        progress = [round(100 * pos / sup_position, 1) for pos in pos_end]

        regions = zip([chrom] * len(pos_start), pos_start, pos_end)

        logger.info(f"Generated {len(pos_start)} regions for processing")
        return progress, regions
        
    except Exception as e:
        logger.error(f"Error getting regions: {e}")
        raise RuntimeError(f"Failed to get regions: {e}")


def get_data(files: List[str], labels: Optional[List[str]] = None, 
             data_columns: Optional[List[List[Tuple]]] = None, 
             regions: Optional[List[Tuple]] = None, join: str = 'outer',
             preset: str = 'bed') -> Generator[pd.DataFrame, None, None]:
    """Combines tabix-indexed genome position files.
    
    Args:
        files: List of file paths
        labels: List of labels for files (optional)
        data_columns: List of data columns specifications (optional)
        regions: List of regions to process (optional)
        join: Type of join operation (default: 'outer')
        preset: File format preset (default: 'bed')
        
    Yields:
        DataFrame: Combined data for each region
        
    Raises:
        InputMismatchError: If input parameters don't match
        MissingInputError: If required parameters are missing
        RuntimeError: If data processing fails
    """
    logger.debug(f"Getting data for {len(files)} files with preset: {preset}")
    
    # Validate input arguments
    if not files:
        raise MissingInputError("Files list cannot be empty")
    
    try:
        if labels is None:
            keys = ['unit_{}'.format(pos + 1) for pos, _ in enumerate(files)]
        elif len(labels) == len(files):
            keys = labels
        else:
            raise InputMismatchError('Number of files and labels must match!')

        if data_columns is None:
            raise MissingInputError(
                'The list of data_columns must have at least one entry!')
        elif len(data_columns) == len(files):
            pass
        elif len(data_columns) == 1:
            data_columns = data_columns * len(files)
        else:
            raise InputMismatchError(
                'Either supply a single entry in data_columns or '
                'the number of entries must match the number of files!')

        # Configure preset-specific settings
        if preset == 'bed':
            index = [
                (0, '#chrom', str),
                (1, 'start', np.int64),
                (2, 'end', np.int64)]
            index_col = [i[0] for i in index]
        elif preset == 'gff':
            # TODO: Implement GFF support
            logger.warning("GFF preset not yet implemented")
            raise NotImplementedError("GFF preset not yet implemented")
        elif preset == 'vcf':
            # TODO: Implement VCF support
            logger.warning("VCF preset not yet implemented")
            raise NotImplementedError("VCF preset not yet implemented")
        elif preset == 'sam':
            # TODO: Implement SAM support
            logger.warning("SAM preset not yet implemented")
            raise NotImplementedError("SAM preset not yet implemented")
        else:
            raise ValueError(f"Unsupported preset: {preset}")

        # Output columns
        names = ['sampling_unit', 'feature']
        columns = [index + cols for cols in data_columns]

        if regions is None:
            logger.warning("No regions provided")
            return

        for region in regions:
            try:
                query = '{0}:{1}-{2}'.format(*region)
                logger.debug(f"Processing region: {query}")

                # Create tabix processes
                tabix_processes = []
                for file_ in files:
                    try:
                        process = subprocess.Popen(
                            ['tabix', file_, query],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            universal_newlines=True
                        )
                        tabix_processes.append(process)
                    except FileNotFoundError:
                        logger.error(f"tabix command not found or file not accessible: {file_}")
                        raise RuntimeError("tabix command not found or file not accessible")

                tabix = enumerate(tabix_processes)

                # Create dataframes
                dframes = []
                for i, tbx in tabix:
                    try:
                        df = pd.read_table(
                            tbx.stdout,
                            header=None,
                            index_col=index_col,
                            comment='#',
                            usecols=[f[0] for f in columns[i]],
                            names=[f[1] for f in columns[i]],
                            dtype={f[1]: f[2] for f in columns[i]}
                        )
                        dframes.append(df)
                        
                        # Wait for process to complete and check for errors
                        return_code = tbx.wait()
                        if return_code != 0:
                            stderr_output = tbx.stderr.read()
                            logger.warning(f"tabix process returned code {return_code}: {stderr_output}")
                            
                    except Exception as e:
                        logger.error(f"Error reading data from tabix process: {e}")
                        # Continue with empty dataframe
                        dframes.append(pd.DataFrame())

                # Merge dataframes
                if dframes:
                    merged_dframe = pd.concat(
                        dframes, axis=1, keys=keys, names=names, join=join)
                    logger.debug(f"Merged dataframe shape: {merged_dframe.shape}")
                    yield merged_dframe
                else:
                    logger.warning("No dataframes to merge")
                    yield pd.DataFrame()
                    
            except Exception as e:
                logger.error(f"Error processing region {region}: {e}")
                # Yield empty dataframe for this region
                yield pd.DataFrame()
                
    except Exception as e:
        logger.error(f"Fatal error in get_data: {e}")
        raise RuntimeError(f"Data processing failed: {e}")


def supremum_numsites(tabixfiles: Union[str, List[str]], chrom: str) -> Optional[int]:
    """Return the least upper bound for the number of covered sites.
    
    Args:
        tabixfiles: Path(s) to tabix-indexed files
        chrom: Chromosome identifier
        
    Returns:
        Maximum number of sites or None if no data
        
    Raises:
        RuntimeError: If subprocess execution fails
    """
    logger.debug(f"Computing supremum number of sites for chromosome: {chrom}")
    
    if isinstance(tabixfiles, str):
        tabixfiles = [tabixfiles]
    
    sites = []

    for f in tabixfiles:
        try:
            logger.debug(f"Processing file: {f}")
            
            tabix = subprocess.Popen(
                ["tabix", f, chrom], 
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            wcl = subprocess.Popen(
                ["wc", "-l"], 
                stdin=tabix.stdout, 
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            
            tabix.stdout.close()  # Allow tabix to receive a SIGPIPE if wcl exits.
            
            stdout, stderr = wcl.communicate()
            
            if wcl.returncode != 0:
                logger.warning(f"wc command failed for file {f}: {stderr}")
                continue
            
            try:
                site_count = int(stdout)
                sites.append(site_count)
                logger.debug(f"Site count for {f}: {site_count}")
            except ValueError as e:
                logger.warning(f"Could not parse site count for file {f}: {e}")
                continue
                
        except Exception as e:
            logger.error(f"Error processing file {f}: {e}")
            continue

    try:
        out = np.max(sites) if sites else None
        logger.debug(f"Supremum number of sites: {out}")
        return out
    except ValueError:
        logger.warning("No valid site counts found")
        return None


def supremum_position(tabixfiles: Union[str, List[str]], chrom: str) -> Optional[int]:
    """Return the least upper bound for the chrom end coordinate.
    
    Args:
        tabixfiles: Path(s) to tabix-indexed files
        chrom: Chromosome identifier
        
    Returns:
        Maximum position or None if no data
        
    Raises:
        RuntimeError: If subprocess execution fails
    """
    logger.debug(f"Computing supremum position for chromosome: {chrom}")
    
    if isinstance(tabixfiles, str):
        tabixfiles = [tabixfiles]
    
    end_coordinate = []

    for f in tabixfiles:
        try:
            logger.debug(f"Processing file: {f}")
            
            tabix = subprocess.Popen(
                ["tabix", f, chrom], 
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            tail = subprocess.Popen(
                ["tail", "-1"], 
                stdin=tabix.stdout, 
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            cut = subprocess.Popen(
                ["cut", "-f3"], 
                stdin=tail.stdout, 
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            
            # Allow first process to receive a SIGPIPE if process 2 exits.
            tabix.stdout.close()
            
            stdout, stderr = cut.communicate()
            
            if cut.returncode != 0:
                logger.warning(f"Pipeline failed for file {f}: {stderr}")
                continue
            
            try:
                base_position = int(stdout)
                end_coordinate.append(base_position)
                logger.debug(f"End position for {f}: {base_position}")
            except ValueError as e:
                logger.warning(f"Could not parse end position for file {f}: {e}")
                continue
                
        except Exception as e:
            logger.error(f"Error processing file {f}: {e}")
            continue

    try:
        out = np.max(end_coordinate) if end_coordinate else None
        logger.debug(f"Supremum position: {out}")
        return out
    except ValueError:
        logger.warning("No valid end coordinates found")
        return None


def groupby(by: str, metadata: Optional[pd.DataFrame] = None, 
            data: Optional[pd.DataFrame] = None) -> pd.core.groupby.DataFrameGroupBy:
    """Group merged GPF data frame by levels of a factor.
    
    Args:
        by: Column name to group by
        metadata: Metadata DataFrame (optional)
        data: Data DataFrame (optional)
        
    Returns:
        Grouped DataFrame
        
    Raises:
        ValueError: If required parameters are missing
        KeyError: If grouping column is not found
    """
    logger.debug(f"Grouping data by: {by}")
    
    if metadata is None:
        raise ValueError("Metadata cannot be None")
    if data is None:
        raise ValueError("Data cannot be None")
    if by not in metadata.columns:
        raise KeyError(f"Column '{by}' not found in metadata")
    
    try:
        mapping = metadata.set_index('label')[by].to_dict()
        result = data.groupby(mapping, axis=1, level=0)
        logger.debug(f"Successfully grouped data by {by}")
        return result
    except Exception as e:
        logger.error(f"Error grouping data: {e}")
        raise
