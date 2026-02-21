"""
================================================================================
M15_M92_ATLAS.PY
================================================================================
PURPOSE:
    Generate ATLAS9 model atmospheres for M15 and M92 globular cluster stars
    using BasicATLAS. Creates both nominal models and perturbed models for
    error propagation (±Teff, ±logg, ±[Fe/H], ±vt).

DESCRIPTION:
    This script processes stellar parameters from M15 and M92 photometric and
    spectroscopic catalogs to generate custom ATLAS9 model atmospheres. For
    each star, it:
    1. Loads stellar parameters (Teff, logg, [Fe/H], vt) and abundances
    2. Generates opacity distribution functions (ODF) with DFSYNTHE
    3. Computes ATLAS9 model atmosphere at nominal parameters
    4. Computes perturbed models for error analysis (down-Teff, down-logg,
       down-[Fe/H], down-vt)
    5. Saves all models and settings for later spectral synthesis

INPUTS:
    /raid/caltech/hires/M15_M92/avg_abund.pkl - Pickled catalog containing:
        - cat: Stellar catalog with photometric parameters
        - avg_abunds: Average abundances per cluster and deviations per star
        - star_names: List of star names for M15 and M92

OUTPUTS:
    /raid/caltech/hires/M15_M92/M15_M92_params.pkl - Pickled settings objects
    /raid/atlas/BasicATLAS/DFSYNTHE_ODF_* - ODF directories per star
    /raid/atlas/BasicATLAS/ATLAS_LMHA_* - Model atmosphere directories

DEPENDENCIES:
    - BasicATLAS (custom ATLAS9 wrapper)
    - Intel OneAPI libraries (for ATLAS9 Fortran code)
    - Python: numpy, pandas, multiprocessing, pickle

AUTHOR: E. Kirby
MODIFIED: 2025
================================================================================
"""

import sys
sys.path.append('/raid/BasicATLAS/')
import atlas
import pickle

import os
import numpy as np

# Set Intel OneAPI library path for ATLAS9 Fortran binaries
os.environ['LD_LIBRARY_PATH'] = "/opt/intel/oneapi/2025.2/lib:$LD_LIBRARY_PATH"

from multiprocessing import Pool
import multiprocessing as mp

import shutil
import copy


# =============================================================================
# FUNCTION: setup_atlas
# =============================================================================
def setup_atlas(settings):
    """
    Generate ATLAS9 model atmospheres for a single star.
    
    Creates 5 model atmospheres:
    1. Nominal: using best-fit stellar parameters
    2. Down-Teff: Teff reduced by tefferr (for error propagation)
    3. Down-logg: logg reduced by loggerr
    4. Down-[Fe/H]: [Fe/H] reduced by zscaleerr
    5. Down-vt: vturb reduced by vturberr
    
    Parameters
    ----------
    settings : atlas.Settings
        Settings object containing stellar parameters (teff, logg, zscale, vturb)
        and their uncertainties, plus elemental abundances.
    
    Returns
    -------
    None
        Writes model atmosphere files to disk in /raid/atlas/BasicATLAS/
    
    Notes
    -----
    - Uses DFSYNTHE to generate opacity distribution functions (ODFs)
    - ODFs depend on abundances, so nominal and down-[Fe/H] need separate ODFs
    - Other perturbed models can reuse the nominal ODF
    """
    # -------------------------------------------------------------------------
    # NOMINAL MODEL: Generate ODF and model atmosphere
    # -------------------------------------------------------------------------
    run_dir = os.path.expanduser(f'/raid/atlas/BasicATLAS/DFSYNTHE_ODF_{settings.name}')
    shutil.rmtree(run_dir, ignore_errors=True)
    atlas.dfsynthe(run_dir, settings)

    model_dir = os.path.expanduser(f'/raid/atlas/BasicATLAS/ATLAS_LMHA_{settings.name}')
    shutil.rmtree(model_dir, ignore_errors=True)
    atlas.atlas(model_dir, settings, ODF=run_dir)
    
    # -------------------------------------------------------------------------
    # DOWN-TEFF MODEL: Teff - tefferr (for error analysis)
    # -------------------------------------------------------------------------
    settings_upteff = settings.copy()
    settings_upteff.teff = settings.teff - settings.tefferr
    print(settings.teff, settings.tefferr, settings_upteff.teff)
    model_dir_upteff = os.path.expanduser(f'/raid/atlas/BasicATLAS/ATLAS_LMHA_{settings.name}_downteff')
    shutil.rmtree(model_dir_upteff, ignore_errors=True)
    atlas.atlas(model_dir_upteff, settings_upteff, ODF=run_dir)
    
    # -------------------------------------------------------------------------
    # DOWN-LOGG MODEL: logg - loggerr
    # -------------------------------------------------------------------------
    settings_uplogg = settings.copy()
    settings_uplogg.logg = settings.logg - settings.loggerr
    model_dir_uplogg = os.path.expanduser(f'/raid/atlas/BasicATLAS/ATLAS_LMHA_{settings.name}_downlogg')
    shutil.rmtree(model_dir_uplogg, ignore_errors=True)
    atlas.atlas(model_dir_uplogg, settings_uplogg, ODF=run_dir) 
    
    # -------------------------------------------------------------------------
    # DOWN-[Fe/H] MODEL: [Fe/H] - zscaleerr
    # -------------------------------------------------------------------------
    # Note: Requires new ODF since abundances affect opacity
    settings_upfeh = settings.copy()
    settings_upfeh.zscale = settings.zscale - settings.zscaleerr
    # Scale all elemental abundances by the same amount
    for el in settings_upfeh.abun.keys():
        settings_upfeh.abun[el] = settings_upfeh.abun[el] - settings_upfeh.zscaleerr
    
    run_dir_upfeh = os.path.expanduser(f'/raid/atlas/BasicATLAS/DFSYNTHE_ODF_{settings.name}_downfeh')
    shutil.rmtree(run_dir_upfeh, ignore_errors=True)
    atlas.dfsynthe(run_dir_upfeh, settings_upfeh)
    
    model_dir_upfeh = os.path.expanduser(f'/raid/atlas/BasicATLAS/ATLAS_LMHA_{settings.name}_downfeh')
    shutil.rmtree(model_dir_upfeh, ignore_errors=True)
    atlas.atlas(model_dir_upfeh, settings_upfeh, ODF=run_dir_upfeh)
    
    # -------------------------------------------------------------------------
    # DOWN-VT MODEL: vturb - vturberr
    # -------------------------------------------------------------------------
    settings_upvt = settings.copy()
    settings_upvt.vturb = settings.vturb - settings.vturberr
    model_dir_upvt = os.path.expanduser(f'/raid/atlas/BasicATLAS/ATLAS_LMHA_{settings.name}_downvt')
    shutil.rmtree(model_dir_upvt, ignore_errors=True)
    atlas.atlas(model_dir_upvt, settings_upvt, ODF=run_dir)


# =============================================================================
# FUNCTION: process_single_star
# =============================================================================
def process_single_star(args):
    """
    Process a single star: load parameters and generate ATLAS9 models.
    
    This function is designed to be called in parallel via multiprocessing.
    
    Parameters
    ----------
    args : tuple
        (star, gc_index, star_index, cat, avg_abunds, names) where:
        - star (str): Star name (e.g., 'M15-star-1')
        - gc_index (int): Globular cluster index (0=M15, 1=M92)
        - star_index (int): Index within cluster's star list
        - cat (DataFrame): Catalog with photometric parameters
        - avg_abunds (dict): Average abundances and deviations per cluster
        - names (list): List of all star names in catalog
    
    Returns
    -------
    settings : atlas.Settings or None
        Settings object with all stellar parameters, or None if error
    
    Notes
    -----
    - Stellar parameters (Teff, logg) come from photometry
    - [Fe/H] = -2.41 ± 0.05 (fixed for both M15 and M92)
    - Microturbulence vt = 2 km/s (fixed, logg-dependent relation commented out)
    - Individual elemental abundances from avg_abunds (deviations from cluster mean)
    - O abundance derived from Na via Na-O anticorrelation (Gratton et al. 2019)
    """
    star, gc_index, star_index, cat, avg_abunds, names = args
    
    # -------------------------------------------------------------------------
    # FIND STAR IN CATALOG
    # -------------------------------------------------------------------------
    k = names.index(star)
    if k < 0:
        raise ValueError(f"Star {star} not found in catalog names.")
        return None
    
    # -------------------------------------------------------------------------
    # SETUP STELLAR PARAMETERS
    # -------------------------------------------------------------------------
    settings = atlas.Settings()
    settings.name = star
    
    # Photometric temperature and gravity
    settings.teff = cat.iloc[k]['TEFFPHOT']
    settings.tefferr = cat.iloc[k]['TEFFPHOTERR']
    settings.logg = cat.iloc[k]['LOGGPHOT']
    settings.loggerr = cat.iloc[k]['LOGGPHOTERR']
    
    # Metallicity (fixed for both M15 and M92)
    settings.zscale = -2.41
    settings.zscaleerr = 0.05
    
    # Microturbulence (fixed value; logg-dependent relation commented)
    settings.vturb = 2  # km/s
    # Alternative: 2.13 - 0.23*settings.logg
    settings.vturberr = -1  # Negative = not used for error analysis
    # Alternative: np.sqrt(0.05**2 + (0.03*settings.logg)**2 + (0.23*settings.loggerr)**2)
    
    # -------------------------------------------------------------------------
    # SETUP ELEMENTAL ABUNDANCES
    # -------------------------------------------------------------------------
    # Combine cluster-average abundances with individual deviations
    abunall = {
        el: avgabund + abunddiff for el, avgabund, abunddiff in zip(
            avg_abunds[gc_index]['element'], 
            avg_abunds[gc_index]['avgabund'], 
            [ad[star_index] for ad in avg_abunds[gc_index]['abunddiff']]
        )
    }
    
    # Convert from absolute abundances to [X/H] relative to solar and scaled
    for el, abun in abunall.items():
        if el not in settings.abun_solar().keys():
            raise ValueError(f"Element {el} not found in solar abundances.")
        abunall[el] = abunall[el] - settings.abun_solar()[el] - settings.zscale
    
    # Keep only finite abundances
    settings.abun = {k: v for k, v in abunall.items() if np.isfinite(v)}
    
    # -------------------------------------------------------------------------
    # OXYGEN ABUNDANCE FROM Na-O ANTICORRELATION
    # -------------------------------------------------------------------------
    # Set the O abundance based on Na abundance following the trend for
    # NGC 2808 (Gratton et al. 2019 Fig. 2, based on Carretta 2015)
    # Two regimes: high-Na and low-Na branches
    if settings.abun['Na'] > 0.2:
        settings.abun['O'] = -2.375 * settings.abun['Na'] + 0.625
    else:
        settings.abun['O'] = -2./3. * settings.abun['Na'] + 85./300.
    
    # -------------------------------------------------------------------------
    # GENERATE MODEL ATMOSPHERES
    # -------------------------------------------------------------------------
    setup_atlas(settings)    
    return settings


# =============================================================================
# MAIN FUNCTION
# =============================================================================
def main():
    """
    Main driver: Load stellar catalog and generate ATLAS9 models for all stars.
    
    Workflow:
    1. Load stellar catalog with photometric parameters
    2. Load average abundances and individual deviations
    3. Create task list for all stars in M15 and M92
    4. Process stars in parallel using multiprocessing
    5. Save all settings objects to pickle file
    
    Notes
    -----
    - Uses all available CPU cores for parallel processing
    - Each star gets 5 model atmospheres (nominal + 4 perturbed)
    - Results saved to M15_M92_params.pkl for later use
    """
    # -------------------------------------------------------------------------
    # LOAD INPUT DATA
    # -------------------------------------------------------------------------
    with open('/raid/caltech/hires/M15_M92/avg_abund.pkl', 'rb') as f:
        cat = pickle.load(f)        # Stellar catalog DataFrame
        avg_abunds = pickle.load(f)  # Average abundances per cluster
        star_names = pickle.load(f)  # Star names per cluster
    
    # -------------------------------------------------------------------------
    # PREPARE PROCESSING TASKS
    # -------------------------------------------------------------------------
    gcs = ['M15', 'M92']
    names = [name.decode('utf-8').strip() for name in cat['NAME'].values]
    
    tasks = []
    results = []
    
    # Build task list: one task per star in each globular cluster
    for j, gc in enumerate(gcs):
        for i, star in enumerate(star_names[j]):
            # Option to process serially for debugging:
            # results.append(process_single_star((star, j, i, cat, avg_abunds, names)))
            tasks.append((star, j, i, cat, avg_abunds, names))
    
    # -------------------------------------------------------------------------
    # PARALLEL PROCESSING
    # -------------------------------------------------------------------------
    # Process all stars in parallel using all available CPU cores
    # Each star generates 5 ATLAS9 models (~5-10 min per star)
    with Pool(processes=mp.cpu_count()) as pool:
        results = pool.map(process_single_star, tasks)
    
    # -------------------------------------------------------------------------
    # SAVE RESULTS
    # -------------------------------------------------------------------------
    with open('/raid/caltech/hires/M15_M92/M15_M92_params.pkl', 'wb') as f:
        pickle.dump(results, f)
    
    print(f"Processed {len([r for r in results if r])} stars successfully")
    
    
if __name__ == "__main__":
    main()
