; =========================================================================
; MERGE_ALL.PRO
; =========================================================================
; PURPOSE:
; Process and combine multi-order HIRES spectra for M15 and M92 stars.
; This includes:
; - Cross-correlation for radial velocity measurement
; - Continuum normalization using synthetic spectra
; - Merging spectral orders into combined 1D spectra
;
; MAIN PROCEDURES:
; stack_all        - Merge multi-order spectra into single 1D spectrum
; measure_vr       - Measure radial velocities via cross-correlation
; filter_allframes - Filter catalog to entries with existing spectra
;
; DEPENDENCIES:
; - MAKEE reduced spectra (_1.fits, _2.fits, _3.fits + error files)
; - M15_M92_allframes.fits (master catalog from make_allframes.pro)
; - MOOG (synthetic spectra for continuum normalization)
; - IDL astronomy library (x_specrebin, etc.)
;
; CREATED: E. Kirby
; MODIFIED: 2025
; =========================================================================

; -------------------------------------------------------------------------
; FUNCTION: findpix
; -------------------------------------------------------------------------
; PURPOSE:
; Find the pixel index closest to a target wavelength
;
; INPUTS:
; lambda  - Wavelength array
; lambda0 - Target wavelength to find
;
; OUTPUT:
; Pixel index of closest wavelength (zero-indexed, with safeguard)
; -------------------------------------------------------------------------
function findpix, lambda, lambda0
  compile_opt idl2
  diff = abs(lambda - lambda0)
  minval = min(diff, minsub)
  minsub = (minsub - 1) > 0 ; Ensure index >= 0
  return, minsub
end

; -------------------------------------------------------------------------
; FUNCTION: zfind
; -------------------------------------------------------------------------
; PURPOSE:
; Cross-correlate object spectrum with reference to determine redshift
; (Doppler shift). Adapted from SDSS spectroscopic pipeline.
;
; INPUTS:
; lambdaref - Reference wavelength array [Angstroms]
; specref   - Reference flux array
; ivarref   - Reference inverse variance array
; lambda    - Object wavelength array [Angstroms]
; spec      - Object flux array
; ivar      - Object inverse variance array
;
; OPTIONAL KEYWORDS:
; zmin, zmax - Minimum/maximum redshift to search
; zguess     - Initial guess for redshift
; pwidth     - Pixel width for search window around zguess
; nfind      - Number of peaks to find in cross-correlation
; width      - Width of cross-correlation peak
; wvmin, wvmax - Wavelength range to use for cross-correlation
;
; OUTPUT:
; Redshift (z) of best match. For stellar RVs: z = v/c
; zerr - Redshift uncertainty returned via keyword
; -------------------------------------------------------------------------
function zfind, lambdaref, specref, ivarref, lambda, spec, ivar, zmin = zmin, zmax = zmax, $
  zguess = zguess, pwidth = pwidth, nfind = nfind, width = width, $
  wvmin = wvmin, wvmax = wvmax, zerr = zerr
  compile_opt idl2

  ; Common block for template storage (legacy from SDSS pipeline)
  common com_zfind, starflux, starloglam0, stardloglam, $
    nstars, npoints, starflux_corrected, thisfile

  ; -----------------------------------------------------------------
  ; SETUP LOG-LAMBDA GRID FOR REFERENCE SPECTRUM
  ; -----------------------------------------------------------------
  starloglam0 = alog10(lambdaref[0])
  stardloglam = alog10(lambdaref[1]) - starloglam0
  npoints = n_elements(lambdaref)
  starorder = round(3. * npoints / 7000.) ; 3 polynomial terms per 7000 points

  ; -----------------------------------------------------------------
  ; REBIN OBJECT SPECTRUM TO REFERENCE WAVELENGTH GRID
  ; -----------------------------------------------------------------
  ; Only use pixels with valid data in the overlapping region
  w = where(lambda ge min(lambdaref) and lambda le max(lambdaref) and ivar gt 0)
  x_specrebin, lambda[w], spec[w], lambdaref, objflux, $
    var = ivar[w] ^ (-1), nwvar = varrebin, /silent
  objivar = varrebin ^ (-1)

  ; Select valid pixels in overlapping region
  w = where(lambdaref ge min(lambda) and lambdaref le max(lambda) and objivar gt 0, npix)
  objflux = objflux[w]
  objivar = objivar[w]
  objloglam0 = alog10(lambdaref[w[0]])
  objdloglam = stardloglam
  loglam = objloglam0 + dindgen(npix) * objdloglam

  ; -----------------------------------------------------------------
  ; TRIM SPECTRUM TO SPECIFIED WAVELENGTH RANGE
  ; -----------------------------------------------------------------
  ; Check if we need to trim the input spectrum to just a subregion
  if n_elements(wvmin) gt 0 then wvmin = wvmin[0] else wvmin = -1
  if n_elements(wvmax) gt 0 then wvmax = wvmax[0] else wvmax = -1
  if wvmax ge 0 and wvmin ge 0 then begin
    minpix = findpix(loglam, alog10(wvmin))
    maxpix = findpix(loglam, alog10(wvmax))
    loglam = loglam[minpix : maxpix]
    objflux = objflux[minpix : maxpix]
    objivar = objivar[minpix : maxpix]
  endif

  ; -----------------------------------------------------------------
  ; SET REDSHIFT SEARCH LIMITS
  ; -----------------------------------------------------------------
  ; Check if the zmin and zmax arguments were passed. If so, then
  ; convert the redshift values into pixel values pmin and pmax.
  if n_elements(zmin) ne 0 then $
    pmin = floor(alog10(1.0 + zmin) / objdloglam)
  if n_elements(zmax) ne 0 then $
    pmax = ceil(alog10(1.0 + zmax) / objdloglam)

  ; Check if a guess redshift zguess was passed along with a pixel
  ; window pwidth. If so, then reset pmin and pmax according to the
  ; guess value and the window.
  if n_elements(zguess) gt 0 and keyword_set(pwidth) then begin
    if keyword_set(width) then width1 = width $
    else width1 = pwidth
    pmin = floor(alog10(1.0 + zguess) / objdloglam - 0.5 * (pwidth + 1 + width1))
    pmax = floor(alog10(1.0 + zguess) / objdloglam + 0.5 * (pwidth + 1 + width1))
  endif

  ; Limit pmax to ensure sufficient overlap between spectra
  ; Alteration by BJW, 8/21/03
  maxp = long((objloglam0 + objdloglam * 0.99 * n_elements(objflux) - starloglam0) / $
    objdloglam)
  if maxp lt pmax then begin ; Limit upper redshift range to have overlap
    pmax = maxp
    print, 'Resetting pmax to: ', maxp
  endif

  ; Verify wavelength grids match in resolution
  if abs(objdloglam - stardloglam) gt 0.05 * objdloglam then $
    message, 'Template and object lambda resolution do NOT match!'

  ; -----------------------------------------------------------------
  ; COMPUTE CROSS-CORRELATION
  ; -----------------------------------------------------------------
  ; Compute the redshift difference between the first pixel of the object
  ; spectra and the template
  poffset = (objloglam0 - starloglam0) / objdloglam

  ; Call zcompute.pro to compute the redshift(s)
  zans = zcompute(objflux, objivar, specref, poffset = poffset, $
    pmin = pmin, pmax = pmax, nfind = nfind, width = width, $
    plottitle = plottitle, _extra = extra)

  ; -----------------------------------------------------------------
  ; CONVERT PIXEL SHIFTS TO REDSHIFTS
  ; -----------------------------------------------------------------
  ; Convert redshift (and error) from pixels to the conventional dimensionless
  ; value. Do not modify any errors that are less than zero, since those
  ; can be used as just warning flags from the fit.

  indx = where(zans.dof gt 0, npeak)
  if (npeak gt 0) then $
    zans[indx].z = 10. ^ (objdloglam * zans[indx].z) - 1.

  jndx = where(zans.dof gt 0 and zans.z_err ge 0)
  if (jndx[0] ne -1) then $
    zans[jndx].z_err = $
      alog(10d) * objdloglam * zans[jndx].z_err * (1 + zans[jndx].z)

  zerr = zans[indx].z_err
  return, zans[indx].z
end

; ------------------------------------------------------------------------------

; -------------------------------------------------------------------------
; FUNCTION: continuum
; -------------------------------------------------------------------------
; PURPOSE:
; Estimate continuum of a spectrum for normalization. Three methods
; available: spline fitting, polynomial fitting, or synthetic spectrum
; (recommended for high-resolution stellar spectra).
;
; INPUTS:
; lambda - Wavelength array [Angstroms]
; spec   - Flux array
; error  - Error array (1-sigma uncertainties)
;
; OPTIONAL KEYWORDS:
; name - Star name for looking up stellar parameters (required for 'synth')
;
; OUTPUT:
; Continuum flux array (same length as input)
;
; METHOD:
; 'synth' (default): Generate synthetic spectrum using MOOG with stellar
; parameters from M15_M92_allframes.fits, then fit smooth spline
; to observed/synthetic ratio to capture continuum shape
; 'spline': Direct spline fit to observed spectrum
; 'poly': Iteratively sigma-clipped polynomial fit
; -------------------------------------------------------------------------
function continuum, lambda, spec, error, name = name
  compile_opt idl2
  method = 'synth'

  ; Convert errors to inverse variance
  ivar = error ^ (-2)
  w = where(finite(ivar) and finite(spec) and ivar gt 0 and error gt 0, c)

  case method of
    ; -----------------------------------------------------------------
    ; SPLINE METHOD: Direct spline fit to observed spectrum
    ; -----------------------------------------------------------------
    'spline': begin
      bkpt = slatec_splinefit(lambda[w], spec[w], coeff, invvar = ivar[w], $
        everyn = 500 < round(c / 4), upper = 2, lower = 2, /silent)
      cont = slatec_bvalu(lambda, bkpt, coeff)
    end

    ; -----------------------------------------------------------------
    ; POLYNOMIAL METHOD: Iteratively sigma-clipped polynomial fit
    ; -----------------------------------------------------------------
    'poly': begin
      pold = dblarr(4)
      p = pold + 1
      iter = 0
      while ~array_equal(p, pold) do begin
        pold = p
        p = poly_fit(lambda[w], spec[w], 3, measure_errors = error[w])
        cont = poly(lambda, reform(p))
        sigma = stddev(spec / cont)
        ; Reject pixels > 5-sigma from continuum
        w = where(finite(ivar) and finite(spec) and ivar gt 0 and error gt 0 and $
          (spec / cont - 1.0) lt -5.0 * sigma or (spec / cont - 1.0) gt 5.0 * sigma, c)
        iter++
        if c lt 100 or iter ge 5 then break
      endwhile
    end

    ; -----------------------------------------------------------------
    ; SYNTHETIC SPECTRUM METHOD (recommended)
    ; -----------------------------------------------------------------
    'synth': begin
      ; Load stellar parameters from master catalog
      allframes = mrdfits('M15_M92_allframes.fits', 1, /silent)
      allframes = allframes[where(strtrim(allframes.name, 2) eq strtrim(name, 2))]

      ; Prepare line list for synthetic spectrum
      linelist = mrdfits('cont/hires.list.fits.gz', 1, /silent)
      linelist = linelist[where(linelist.lambda ge min(lambda) and $
        linelist.lambda le max(lambda), nlist)]

      ; Write MOOG-format line list
      openw, lun, 'synthcont/hires.list', /get_lun
      printf, lun, 'HIRESLIST'
      for i = 0, nlist - 1 do printf, lun, linelist[i].lambda, linelist[i].species, $
        linelist[i].ep, linelist[i].loggf, $
        format = '(D10.3,D10.1,D10.2,D10.3)'
      close, lun
      free_lun, lun

      ; Create MOOG parameter file for synthesis
      make_par, parfile = 'synthcont/hires.par', linefile = 'synthcont/hires.list', $
        atmfile = 'synthcont/hires.atm', outfile = 'synthcont/hires.out2', $
        driver = 'synth', minlambda = min(lambda), maxlambda = max(lambda), $
        stronglist = 'cont/hires.strong', c12c13 = 10.0

      ; Compute [C/Fe] from log(g) (carbon depletion in giants)
      cfe = 0.2 + ((0.2 + 1.0) / (2.4 - 0.7)) * ((allframes.loggphot < 2.4) - 2.4)

      ; Interpolate model atmosphere at stellar parameters
      interp_atm, allframes.teff_mb20, allframes.loggphot, $
        2.13 - 0.23 * allframes.loggphot, -2.4, 0.3, $
        outfile = 'synthcont/hires.atm', $
        tweakel = [6], tweakabund = [cfe]

      ; Run MOOG to generate synthetic spectrum
      spawn, 'MOOGSILENT synthcont/hires.par'

      ; Read synthetic spectrum and apply Doppler shift
      moog = read_moog_spec('synthcont/hires.out2', /newmoog, /singleabund)
      mooglambda = moog.lambda * (1d + (allframes.vr) / 2.99792d5)
      airtovac, mooglambda ; Convert to vacuum wavelengths
      moogspec = moog.spec

      ; Convolve to HIRES resolution (R ~ 30,000)
      moogspec = smooth_gauss_wrapper(mooglambda, moogspec, lambda, $
        lambda / 30000. / 2.35)

      ; Fit spline to observed/synthetic ratio
      w = where(finite(ivar) and finite(spec) and ivar gt 0 and error gt 0 and $
        finite(moogspec) and moogspec gt 0.01, c)
      if c lt 20 then cont = 1 else begin
        bkpt = slatec_splinefit(lambda[w], spec[w] / moogspec[w], coeff, $
          invvar = ivar[w] * moogspec[w] ^ 2., $
          everyn = 500 < round(c / 2), upper = 1.5, lower = 1.5, /silent)
        cont = slatec_bvalu(lambda, bkpt, coeff)
      endelse

      ; Optional: Diagnostic plot
      if 0 then begin
        splot, lambda, spec
        soplot, lambda, moogspec * median(spec) / median(moogspec), color = fsc_color('red')
        soplot, lambda, spec / moogspec, color = fsc_color('blue')
        soplot, lambda, cont, color = fsc_color('orange')
      endif
    end
  endcase

  ; Optional: Diagnostic plot for specific wavelength region
  if 0 and min(lambda) lt 5160 and max(lambda) gt 5160 then begin
    splot, lambda, spec
    soplot, lambda, cont, color = fsc_color('red')
    stop
  endif

  return, cont
end

; -------------------------------------------------------------------------
; PROCEDURE: stack_all
; -------------------------------------------------------------------------
; PURPOSE:
; Merge multi-order HIRES spectra into single 1D continuum-normalized
; spectra. Processes all stars in M15_M92_allframes.fits.
;
; INPUTS:
; None (reads from M15_M92_allframes.fits)
;
; OUTPUT:
; Individual FITS files for each star in spectra/ directory
; Format: {name, date, exptime, mjd, instrument_config, lambda, spec, ivar}
;
; METHOD:
; For each star and echelle order:
; 1. Read multi-order MAKEE spectra (_1.fits, _2.fits, _3.fits)
; 2. Normalize each order using continuum() function
; 3. Rebin to common log-lambda grid
; 4. Combine using inverse-variance weighting
; 5. Convert to air wavelengths
; 6. Save as gzipped FITS file
; -------------------------------------------------------------------------
pro stack_all
  compile_opt idl2
  allframes = mrdfits('M15_M92_allframes.fits', 1, /silent)
  n = n_elements(allframes)

  ; -----------------------------------------------------------------
  ; LOOP OVER ALL STARS IN CATALOG
  ; -----------------------------------------------------------------
  for i = 0, n - 1 do begin
    file_exist = 0
    ; Skip these problematic stars
    if strtrim(allframes[i].name, 2) ne 'M92-star-5' and $
      strtrim(allframes[i].name, 2) ne 'M92-star-7' then continue

    ; -----------------------------------------------------------------
    ; READ MULTI-ORDER MAKEE SPECTRA
    ; -----------------------------------------------------------------
    ; MAKEE produces 3 files per observation: _1, _2, _3 (different orders)
    for j = 1, 3 do begin
      fitsfile = strtrim(allframes[i].makeefile, 2) + '_' + $
        string(j, format = '(I1)') + '.fits'
      if ~file_test(fitsfile) then begin
        print, fitsfile + ' not found'
        continue
      endif
      file_exist = 1
      efitsfile = strtrim(allframes[i].makeefile, 2) + '_' + $
        string(j, format = '(I1)') + 'e.fits'
      print, fitsfile

      ; Read flux and error arrays
      makeej = mrdfits(fitsfile, 0, hdr, /silent)
      errorj = mrdfits(efitsfile, 0, hdr, /silent)

      ; Concatenate multi-order arrays
      makee = j eq 1 ? makeej : [[makee], [makeej]]
      error = j eq 1 ? errorj : [[error], [errorj]]
      lambda = j eq 1 ? hires_lambda(hdr) : [[lambda], [hires_lambda(hdr)]]
    endfor
    if ~file_exist then continue

    ; -----------------------------------------------------------------
    ; SETUP OUTPUT LOG-LAMBDA GRID
    ; -----------------------------------------------------------------
    npix = (size(lambda))[1]
    norders = (size(lambda))[2]
    lambdarange = minmax(lambda)
    b = (alog10(lambda[1]) - alog10(lambda[0])) ; dlambda/lambda
    a = alog10(lambdarange[0])
    nlambda = ceil((alog10(lambdarange[1]) - a) / b)

    ; Create output structure
    hires = {name: ' ', date: ' ', exptime: 0d, mjd: 0d, xdisp: ' ', decker: ' ', $
      slitwid: 0d, xdangle: 0d, echangle: 0d, $
      lambda: dblarr(nlambda), spec: dblarr(nlambda), ivar: dblarr(nlambda)}
    hires.lambda = 10. ^ (a + b * dindgen(nlambda))

    ; Store metadata from FITS header
    hires.name = allframes[i].name
    hires.date = sxpar(hdr, 'DATE')
    hires.exptime = sxpar(hdr, 'EXPTIME')
    hires.mjd = sxpar(hdr, 'MJD')
    hires.xdisp = sxpar(hdr, 'XDISPERS') ; Cross-disperser
    hires.decker = sxpar(hdr, 'DECKNAME') ; Decker name
    hires.slitwid = sxpar(hdr, 'SLITWID') ; Slit width [arcsec]
    hires.xdangle = sxpar(hdr, 'XDANGL') ; Cross-disperser angle
    hires.echangle = sxpar(hdr, 'ECHANGL') ; Echelle angle

    ; -----------------------------------------------------------------
    ; NORMALIZE AND COMBINE ORDERS
    ; -----------------------------------------------------------------
    for j = 0, norders - 1 do begin
      wave = lambda[*, j]
      flux = makee[*, j]
      fluxerror = error[*, j]

      ; Continuum normalize this order
      cont = continuum(wave, flux, fluxerror, name = hires.name)
      flux /= cont
      fluxerror /= cont

      if n_elements(flux) le 1 then message, 'WTF?'

      ; Rebin to common wavelength grid with inverse-variance weighting
      wj = where(hires.lambda ge min(wave) and hires.lambda le max(wave))
      x_specrebin, wave, flux, hires.lambda[wj], specj, $
        var = (fluxerror) ^ 2d, nwvar = nwvar, /silent
      x_specrebin, wave, dblarr(n_elements(flux)) + 1d, hires.lambda[wj], norm, $
        var = (fluxerror) ^ 2d, /silent
      specj /= norm
      ivarj = norm ^ 2d / nwvar

      ; Set invalid pixels to zero weight
      w0 = where(hires.lambda[wj] lt min(wave) or hires.lambda[wj] gt max(wave) or $
        ~finite(specj) or ~finite(ivarj) or ivarj lt 0d)
      ivarj[w0] = 0d

      ; Combine orders with inverse-variance weighting
      hires.spec[wj] += specj * ivarj
      hires.ivar[wj] += ivarj
    endfor

    ; Normalize by total inverse variance
    w = where(hires.ivar gt 0)
    hires.spec[w] /= hires.ivar[w]

    ; Convert vacuum to air wavelengths for consistency
    lambda_arr = hires.lambda
    vactoair, lambda_arr
    hires.lambda = lambda_arr

    ; Optional: Diagnostic plot
    ; splot, hires.lambda, hires.spec, xrange=[4600, 4800]
    ; soplot, [0, 1d5], [1, 1], color=fsc_color('orange')
    ; stop

    ; -----------------------------------------------------------------
    ; WRITE OUTPUT SPECTRUM
    ; -----------------------------------------------------------------
    specfile = strtrim(allframes[i].specfile, 2)
    mwrfits, hires, specfile, /create
    spawn, 'gzip -f ' + specfile
    allframes[i].specfile = specfile + '.gz'
  endfor
  ; mwrfits, allframes, 'M15_M92_allframes.fits', /create
end

; -------------------------------------------------------------------------
; PROCEDURE: measure_vr
; -------------------------------------------------------------------------
; PURPOSE:
; Measure radial velocities for all stars via cross-correlation with
; a reference star of known RV.
;
; INPUTS:
; None (reads from M15_M92_allframes.fits)
;
; OUTPUT:
; Updates M15_M92_allframes.fits with VR and VRERR columns
;
; METHOD:
; 1. Use M92-star-1 as reference (vr = -39.76 km/s)
; 2. For each star, cross-correlate single echelle order
; 3. Use zfind() to determine velocity shift
; 4. Propagate errors from reference star
;
; NOTES:
; - Uses order 16 by default (good Ca II H/K region)
; - Special handling for problem stars (different orders)
; - Cross-correlation range: 4000-5000 Angstroms
; -------------------------------------------------------------------------
pro measure_vr
  compile_opt idl2
  clight = 2.99792458d5 ; Speed of light [km/s]
  allframes = mrdfits('M15_M92_allframes.fits', 1, /silent)
  n = n_elements(allframes)

  ; -----------------------------------------------------------------
  ; SETUP CROSS-CORRELATION PARAMETERS
  ; -----------------------------------------------------------------
  jorder = 16 ; Default echelle order to use
  wvmin = 4000 ; Minimum wavelength for cross-correlation [Angstroms]
  wvmax = 5000 ; Maximum wavelength for cross-correlation [Angstroms]

  ; Add VR columns to catalog structure
  newstr = {vr: 0d, vrerr: 0d}
  allframes = struct_trimtags(allframes, except_tags = ['vr', 'vrerr'])
  allframes = struct_addtags(allframes, replicate(newstr, n))

  ; -----------------------------------------------------------------
  ; LOAD REFERENCE STAR SPECTRUM
  ; -----------------------------------------------------------------
  vr_ref_star = 'M92-star-1'
  vr_ref = -39.76 ; Reference star RV [km/s] (known from literature)
  vrerr_ref = 0.0 ; Uncertainty in reference RV
  wref = where(strtrim(allframes.name, 2) eq vr_ref_star)

  specfile = getenv('chome') + 'keck/hires/M15_M92/makee/' + vr_ref_star + '_1.fits'
  especfile = getenv('chome') + 'keck/hires/M15_M92/makee/' + vr_ref_star + '_1e.fits'
  specref = mrdfits(specfile, 0, hdr, /silent)
  errorref = mrdfits(especfile, 0, /silent)
  waveref = hires_lambda(hdr)

  ; Rebin reference spectrum to log-lambda grid
  a = alog10(waveref[0, jorder])
  b = alog10(waveref[1, jorder]) - a
  nlambda = ceil((alog10(max(waveref[*, jorder])) - a) / b)
  logwave = 10 ^ (a + b * dindgen(nlambda))

  x_specrebin, waveref[*, jorder], specref[*, jorder], logwave, specref, $
    var = (errorref[*, jorder]) ^ 2, nwvar = varref, /silent
  ivarref = 1. / varref

  ; -----------------------------------------------------------------
  ; MEASURE RV FOR EACH STAR
  ; -----------------------------------------------------------------
  for i = 0, n - 1 do begin
    objname = strtrim(allframes[i].name, 2)

    ; Special handling: use different orders for problem stars
    iorder = jorder
    if objname eq 'M15-star-4' then iorder = 11
    if objname eq 'M92-star-13' then iorder = 18

    ; Read object spectrum
    specfile = getenv('chome') + 'keck/hires/M15_M92/makee/' + objname + '_1.fits'
    especfile = getenv('chome') + 'keck/hires/M15_M92/makee/' + objname + '_1e.fits'
    if ~file_test(specfile) then continue
    spec = mrdfits(specfile, 0, hdr, /silent)
    error = mrdfits(especfile, 0, /silent)
    wave = hires_lambda(hdr)

    ; Rebin object spectrum to same log-lambda grid
    x_specrebin, wave[*, iorder], spec[*, iorder], logwave, spec, $
      var = (error[*, iorder]) ^ 2, nwvar = var, /silent
    ivar = 1. / var

    ; Check for valid pixels in cross-correlation region
    wr = where(logwave ge wvmin and logwave le wvmax and $
      ivarref gt 0 and ivar gt 0 and $
      finite(ivarref) and finite(ivar), c)
    if c lt 20 then begin
      print, objname
      stop
      continue
    endif

    ; -----------------------------------------------------------------
    ; CROSS-CORRELATE WITH REFERENCE STAR
    ; -----------------------------------------------------------------
    if objname ne vr_ref_star then begin
      if c lt 20 then continue

      ; Cross-correlate: search ±100 km/s around reference velocity
      z = zfind(logwave[wr] / (1d + vr_ref / clight), specref[wr], ivarref[wr], $
        logwave[wr], spec[wr], ivar[wr], $
        zmin = -100. / clight, zmax = 100. / clight, zerr = zerr, $
        wvmin = wvmin, wvmax = wvmax)
      vr = clight * z
      vrerr = sqrt((clight * zerr) ^ 2 + vrerr_ref ^ 2)

      ; Optional: Plot comparison of reference and object spectra
      splot, logwave / (1d + vr_ref / clight), specref / median(specref[wr]), $
        xrange = [4440, 4447]
      soplot, logwave / (1d + vr / clight), spec / median(spec[wr]), $
        color = fsc_color('red')
    endif else begin
      ; This is the reference star itself
      vr = vr_ref
      vrerr = vrerr_ref
    endelse

    print, allframes[i].name, vr, vrerr
    ; stop
    allframes[i].vr = vr
    allframes[i].vrerr = vrerr
  endfor

  ; Write updated catalog with radial velocities
  mwrfits, allframes, 'M15_M92_allframes.fits', /create
end

; -------------------------------------------------------------------------
; PROCEDURE: filter_allframes
; -------------------------------------------------------------------------
; PURPOSE:
; Filter M15_M92_allframes.fits to include only stars with existing
; processed spectra files.
;
; INPUTS:
; None (reads M15_M92_allframes.fits)
;
; OUTPUT:
; Overwrites M15_M92_allframes.fits with filtered catalog
;
; METHOD:
; Check for existence of .specfile entries and keep only valid ones
; -------------------------------------------------------------------------
pro filter_allframes
  compile_opt idl2
  allframesfile = 'M15_M92_allframes.fits'
  allframes = mrdfits(allframesfile, 1, /silent)

  ; Ensure .gz extension is present
  allframes.specfile = strtrim(allframes.specfile, 2) + '.gz'

  ; Keep only entries where spectrum file exists
  w = where(file_test(strtrim(allframes.specfile, 2)))
  allframes = allframes[w]
  stop

  ; Write filtered catalog
  mwrfits, allframes, allframesfile, /create
end
