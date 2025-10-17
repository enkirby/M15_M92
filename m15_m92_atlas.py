import sys
sys.path.append('/raid/BasicATLAS/')
import atlas
import pickle

import os
import numpy as np
os.environ['LD_LIBRARY_PATH'] = "/opt/intel/oneapi/2025.2/lib:$LD_LIBRARY_PATH"

from multiprocessing import Pool
import multiprocessing as mp

import shutil
import copy

def setup_atlas(settings):
    run_dir = os.path.expanduser(f'/raid/atlas/BasicATLAS/DFSYNTHE_ODF_{settings.name}')
    shutil.rmtree(run_dir, ignore_errors=True)
    atlas.dfsynthe(run_dir, settings)

    model_dir = os.path.expanduser(f'/raid/atlas/BasicATLAS/ATLAS_LMHA_{settings.name}')
    shutil.rmtree(model_dir, ignore_errors=True)
    atlas.atlas(model_dir, settings, ODF = run_dir)
    
    settings_upteff = settings.copy()
    settings_upteff.teff = settings.teff - settings.tefferr
    print(settings.teff, settings.tefferr, settings_upteff.teff)
    model_dir_upteff = os.path.expanduser(f'/raid/atlas/BasicATLAS/ATLAS_LMHA_{settings.name}_downteff')
    shutil.rmtree(model_dir_upteff, ignore_errors=True)
    atlas.atlas(model_dir_upteff, settings_upteff, ODF = run_dir)
    
    settings_uplogg = settings.copy()
    settings_uplogg.logg = settings.logg - settings.loggerr
    model_dir_uplogg = os.path.expanduser(f'/raid/atlas/BasicATLAS/ATLAS_LMHA_{settings.name}_downlogg')
    shutil.rmtree(model_dir_uplogg, ignore_errors=True)
    atlas.atlas(model_dir_uplogg, settings_uplogg, ODF = run_dir) 
    
    settings_upfeh = settings.copy()
    settings_upfeh.zscale = settings.zscale - settings.zscaleerr
    for el in settings_upfeh.abun.keys():
        settings_upfeh.abun[el] = settings_upfeh.abun[el] - settings_upfeh.zscaleerr
    run_dir_upfeh = os.path.expanduser(f'/raid/atlas/BasicATLAS/DFSYNTHE_ODF_{settings.name}_downfeh')
    shutil.rmtree(run_dir_upfeh, ignore_errors=True)
    atlas.dfsynthe(run_dir_upfeh, settings_upfeh)
    model_dir_upfeh = os.path.expanduser(f'/raid/atlas/BasicATLAS/ATLAS_LMHA_{settings.name}_downfeh')
    shutil.rmtree(model_dir_upfeh, ignore_errors=True)
    atlas.atlas(model_dir_upfeh, settings_upfeh, ODF = run_dir_upfeh)
    
    settings_upvt = settings.copy()
    settings_upvt.vturb = settings.vturb - settings.vturberr
    model_dir_upvt = os.path.expanduser(f'/raid/atlas/BasicATLAS/ATLAS_LMHA_{settings.name}_downvt')
    shutil.rmtree(model_dir_upvt, ignore_errors=True)
    atlas.atlas(model_dir_upvt, settings_upvt, ODF = run_dir)



def process_single_star(args):
    star, gc_index, star_index, cat, avg_abunds, names = args
    
    k = names.index(star)
    if k < 0:
        raise ValueError(f"Star {star} not found in catalog names.")
        return None
    
    settings = atlas.Settings()
    settings.name = star
    settings.teff = cat.iloc[k]['TEFFPHOT']
    settings.tefferr = cat.iloc[k]['TEFFPHOTERR']
    settings.logg = cat.iloc[k]['LOGGPHOT']
    settings.loggerr = cat.iloc[k]['LOGGPHOTERR']
    settings.zscale = -2.41
    settings.zscaleerr = 0.05
    settings.vturb = 2 #2.13 - 0.23*settings.logg
    settings.vturberr = -1 #np.sqrt(0.05**2 + (0.03*settings.logg)**2 + (0.23*settings.loggerr)**2)
    abunall = {el: avgabund+abunddiff for el, avgabund, abunddiff in zip(
        avg_abunds[gc_index]['element'], 
        avg_abunds[gc_index]['avgabund'], 
        [ad[star_index] for ad in avg_abunds[gc_index]['abunddiff']])}
    for el, abun in abunall.items():
        if el not in settings.abun_solar().keys():
            raise ValueError(f"Element {el} not found in solar abundances.")
        abunall[el] = abunall[el] - settings.abun_solar()[el] - settings.zscale
    settings.abun = {k: v for k, v in abunall.items() if np.isfinite(v)}
    
    #Set the O abundance based on Na abundance following the trend for NGC 2808 (Gratton et al. 2019 Fig. 2, which is based on Carretta 2015)
    settings.abun['O'] = -2.375*settings.abun['Na'] + 0.625 if settings.abun['Na'] > 0.2 else -2./3.*settings.abun['Na'] + 85./300.
    
    setup_atlas(settings)    
    return settings

def main():
    with open('/raid/caltech/hires/M15_M92/avg_abund.pkl', 'rb') as f:
        cat = pickle.load(f)
        avg_abunds = pickle.load(f)
        star_names = pickle.load(f)
    
    gcs = ['M15', 'M92']
    names = [name.decode('utf-8').strip() for name in cat['NAME'].values]
    # Prepare tasks
    tasks = []
    results = []
    for j, gc in enumerate(gcs):
        for i, star in enumerate(star_names[j]):
            #results.append(process_single_star((star, j, i, cat, avg_abunds, names)))
            tasks.append((star, j, i, cat, avg_abunds, names))
    
    # Process in parallel
    with Pool(processes=mp.cpu_count()) as pool:
        results = pool.map(process_single_star, tasks)
    
    with open('/raid/caltech/hires/M15_M92/M15_M92_params.pkl', 'wb') as f:
        pickle.dump(results, f)
    
    print(f"Processed {len([r for r in results if r])} stars successfully")
    
    
if __name__ == "__main__":
    main()