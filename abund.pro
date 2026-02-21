;+
; NAME:
;   abund.pro
;
; PURPOSE:
;   IDL procedures for spectroscopic abundance analysis of M15 and M92 HIRES data.
;   Handles line list preparation, MOOG interface, hyperfine structure, and error analysis.
;
; MAIN PROCEDURES:
;   hyperfine_output     - Output hyperfine structure lines in MOOG format
;   trim_linelist        - Filter and prepare line lists for abundance analysis
;   feh_alphafe          - Calculate [Fe/H] and [alpha/Fe] from abundance structure
;   find_abund           - Fit stellar parameters via MOOG abundance equilibrium
;   find_abund_final     - Like find_abund but using ATLAS model atmospheres
;   calculate_abund      - Run MOOG on multiple elements including hyperfine species
;   error_analysis       - Compute abundance errors from parameter uncertainties
;   error_analysis_final - Like error_analysis but using ATLAS models
;   atlas_to_moog        - Convert ATLAS atmosphere format to MOOG format
;   hires_loop           - Main analysis loop for a single star
;   abund                - Batch process all stars
;   abundmc              - Monte Carlo error propagation for all stars
;
; DEPENDENCIES:
;   - MOOG spectral line analysis code (requires MOOGSILENT executable)
;   - elements() function for atomic data
;   - interp_atm procedure for model atmosphere interpolation
;   - read_moog() function to parse MOOG output
;   - make_par procedure to create MOOG parameter files
;
; USAGE:
;   IDL> abund, ni=0           ; Process first star with fixed [Fe/H]
;   IDL> abund, /final         ; Process all stars with ATLAS models
;   IDL> abundmc, ni=5         ; Run Monte Carlo for star #5
;
; NOTES:
;   - Uses common block 'ews' for equivalent width data
;   - Requires environment variable CALTECH to be set
;   - Output files created in ew/, ew2/, and mc/ subdirectories
;
; MODIFICATION HISTORY:
;   Written: Original analysis for M15/M92 HIRES project
;   Cleaned: 2024 - Improved formatting and documentation
;-

pro hyperfine_output, linelist, hyperfine, atomic, lun
  ; Purpose: Write hyperfine structure lines to MOOG format file
  ; Inputs:
  ; linelist - Structure with spectral line data
  ; hyperfine - Structure with hyperfine component data
  ; atomic - Atomic number of element
  ; lun - Logical unit number for output file

  compile_opt idl2

  ; Get all lines for this element
  w = where(linelist.atomic eq atomic, cw)

  ; Loop over each line
  for i = 0, cw - 1 do begin
    ; Find matching hyperfine components
    hyperfinetemp = hyperfine[where(floor(hyperfine.species * 10.) eq $
      floor(linelist[w[i]].species * 10.))]
    deviation = min(abs(linelist[w[i]].lambda - hyperfinetemp.lambda), wfl)

    ; If no close match (deviation > 1 Angstrom), write as single line
    if deviation gt 1.0 then begin
      printf, lun, $
        linelist[w[i]].lambda, $
        linelist[w[i]].species, $
        linelist[w[i]].ep, $
        linelist[w[i]].loggf, $
        linelist[w[i]].ew, $
        format = '(2X,D8.3,6X,D4.1,1X,D9.2,1X,D9.3,20X,D10.1)'

      ; Otherwise, write with hyperfine components
    endif else begin
      ; Determine range of hyperfine components for this line
      whp = where(hyperfinetemp.lambda gt 0)
      we = where(whp eq wfl)
      wll = wfl eq max(whp) ? n_elements(hyperfinetemp) - 1 : (whp[we + 1])[0] - 1
      hyperfinetemp = hyperfinetemp[wfl : wll]
      nh = wll - wfl + 1

      ; Sort components by wavelength
      hyperfinetemp.lambda = abs(hyperfinetemp.lambda)
      hyperfinetemp = hyperfinetemp[sort(hyperfinetemp.lambda)]

      ; Mark secondary components with negative wavelength (MOOG convention)
      if nh gt 1 then hyperfinetemp[1 : nh - 1].lambda *= -1.0

      ; Write primary component with EW
      printf, lun, $
        abs(hyperfinetemp[0].lambda), $
        hyperfinetemp[0].species, $
        hyperfinetemp[0].ep, $
        hyperfinetemp[0].loggf, $
        linelist[w[i]].ew, $
        format = '(2X,D8.3,3X,D7.4,1X,D9.2,1X,D9.3,20X,D10.1)'

      ; Write additional hyperfine components (if any)
      if nh gt 1 then begin
        for j = 1, nh - 1 do begin
          printf, lun, $
            hyperfinetemp[j].lambda, $
            hyperfinetemp[j].species, $
            hyperfinetemp[j].ep, $
            hyperfinetemp[j].loggf, $
            format = '(1X,D9.3,3X,D7.4,1X,D9.2,1X,D9.3)'
        endfor
      endif
    endelse
  endfor
end

function trim_linelist, star, atomic = atomic, teffphot = teffphot, ti = ti, fe = fe, alpha = alpha, nothyperfine = nothyperfine, upperlimit = upperlimit, dirflag = dirflag, uperr = uperr, downerr = downerr
  ; Purpose: Filter and prepare line list for MOOG abundance analysis
  ; Inputs:
  ; star - Star name
  ; atomic - If set, select only lines of this atomic number
  ; teffphot - Flag for photometric Teff analysis
  ; ti - If set, select only Ti lines
  ; fe - If set, select only Fe lines
  ; alpha - If set, select only alpha element lines
  ; nothyperfine - If set, exclude hyperfine structure elements
  ; upperlimit - If set, select only upper limit lines
  ; dirflag - Directory flag for output files
  ; uperr - If set, add +1sigma EW error
  ; downerr - If set, subtract -1sigma EW error
  ; Returns: Number of lines written to output file

  compile_opt idl2
  common ews, ews

  if ~keyword_set(dirflag) then dirflag = ''

  ; Load list of problematic lines to exclude
  readcol, getenv('CALTECH') + 'hires/badlines.dat', $
    badlambda, badspecies, $
    format = 'X,X,D,D', /silent

  linelist = ews
  newstr = replicate({damping: 0d, keep: 0}, n_elements(linelist))
  linelist = struct_addtags(linelist, newstr)
  linelist = linelist[where(linelist.ew ge 0.1)]

  ; badlines = [-1]
  ; for i=0,n_elements(badlambda)-1 do begin
  ; w = where(round(linelist.species*10.) eq round(badspecies[i]*10.) and abs(linelist.lambda-badlambda[i]) lt 0.2, c)
  ; if c gt 0 then badlines = [badlines, w]
  ; endfor
  ; badlines = badlines[where(badlines ge 0, ckeep)]
  ; wkeep = complement(badlines, n_elements(linelist))
  ; linelist = linelist[wkeep]
  nlines = n_elements(linelist)

  if keyword_set(uperr) then linelist.ew += linelist.ewerr
  if keyword_set(downerr) then linelist.ew -= linelist.ewerr

  linelist = linelist[sort(linelist.lambda)]
  linelist = linelist[sort(linelist.species)]
  if ~keyword_set(upperlimit) then linelist = linelist[where(linelist.upperlimit eq 0)]

  rwcut = -4.5
  keep = bytarr(n_elements(linelist)) + 1
  uniqspecies = uniq(linelist.species, sort(linelist.species))
  for i = 0, n_elements(uniqspecies) - 1 do begin
    if linelist[uniqspecies[i]].species ne 22.0 and linelist[uniqspecies[i]].species ne 22.1 and linelist[uniqspecies[i]].species ne 26.0 and linelist[uniqspecies[i]].species ne 26.1 then continue
    w = where(linelist.species eq linelist[uniqspecies[i]].species, c)
    ww = where(alog10(linelist[w].ew * 1d-3 / linelist[w].lambda) gt rwcut, cbig)
    if cbig gt 0 and c - cbig gt 1 then keep[w[ww]] = 0
  endfor
  linelist = linelist[where(keep)]

  cnothyperfine = 0
  calpha = 0
  chyperfine = 0
  cti = 0
  cfe = 0

  if keyword_set(upperlimit) then begin
    wul = where(linelist.upperlimit eq 1, cul)
    if cul eq 0 then return, 0
    linelist = linelist[wul]
  endif

  if keyword_set(nothyperfine) then wnothyperfine = where(linelist.atomic ne 21 and linelist.atomic ne 23 and linelist.atomic ne 25 and linelist.atomic ne 27 and linelist.atomic ne 29 and linelist.atomic ne 56 and linelist.atomic ne 57 and linelist.atomic ne 60 and linelist.atomic ne 63, cnothyperfine)
  ; DOES NOT INCLUDE OXYGEN
  if keyword_set(alpha) then walpha = where(linelist.atomic eq 10 or linelist.atomic eq 12 or linelist.atomic eq 14 or linelist.atomic eq 20, calpha)
  if keyword_set(atomic) then whyperfine = where(linelist.atomic eq atomic and linelist.ew gt 0.0, chyperfine)
  if keyword_set(ti) then wti = where(linelist.atomic eq 22 and linelist.ew gt 0.0, cti)
  if keyword_set(fe) then wfe = where(linelist.atomic eq 26 and linelist.ew gt 0.0, cfe)
  c = cnothyperfine + calpha + cti + cfe + chyperfine
  if c eq 0 then return, c

  if chyperfine gt 0 then begin
    readcol, getenv('CALTECH') + 'hires/M92_KOA/Ji20_hyperfine.moog', lambda, species, ep, loggf, skipline = 1, format = 'D,D,D,D', /silent
    hyperfine = {lambda: 0d, species: 0d, ep: 0d, loggf: 0d}
    hyperfine = replicate(hyperfine, n_elements(lambda))
    hyperfine.lambda = lambda
    hyperfine.species = species
    hyperfine.ep = ep
    hyperfine.loggf = loggf
    ; hyperfine.damping = damping
  endif

  fout = getenv('CALTECH') + 'hires/M15_M92/' + dirflag + strtrim(star, 2) + (~keyword_set(teffphot) ? '' : '_teffphot') + '_temp.ew'
  openw, lun, fout, /get_lun
  printf, lun, star
  case 1 of
    cnothyperfine gt 0: for i = 0, cnothyperfine - 1 do printf, lun, linelist[wnothyperfine[i]].lambda, linelist[wnothyperfine[i]].species, linelist[wnothyperfine[i]].ep, linelist[wnothyperfine[i]].loggf, linelist[wnothyperfine[i]].damping, linelist[wnothyperfine[i]].ew, format = '(2X,D8.3,6X,D4.1,1X,D9.2,1X,D9.3,1X,G9.4,10X,D10.1)'
    chyperfine gt 0: hyperfine_output, linelist, hyperfine, atomic, lun
    else: begin
      if calpha gt 0 then begin
        for i = 0, calpha - 1 do printf, lun, linelist[walpha[i]].lambda, linelist[walpha[i]].species, linelist[walpha[i]].ep, linelist[walpha[i]].loggf, linelist[walpha[i]].damping, linelist[walpha[i]].ew, format = '(2X,D8.3,6X,D4.1,1X,D9.2,1X,D9.3,1X,G9.4,10X,D10.1)'
      endif
      if cti gt 0 then begin
        for i = 0, cti - 1 do printf, lun, linelist[wti[i]].lambda, linelist[wti[i]].species, linelist[wti[i]].ep, linelist[wti[i]].loggf, linelist[wti[i]].damping, linelist[wti[i]].ew, format = '(2X,D8.3,6X,D4.1,1X,D9.2,1X,D9.3,1X,G9.4,10X,D10.1)'
      endif
      if cfe gt 0 then begin
        for i = 0, cfe - 1 do printf, lun, linelist[wfe[i]].lambda, linelist[wfe[i]].species, linelist[wfe[i]].ep, linelist[wfe[i]].loggf, linelist[wfe[i]].damping, linelist[wfe[i]].ew, format = '(2X,D8.3,6X,D4.1,1X,D9.2,1X,D9.3,1X,G9.4,10X,D10.1)'
      endif
    end
  endcase
  close, lun
  free_lun, lun

  return, c
end

function feh_alphafe, abund, alphafe = alphafe, errfeh = errfeh, erralphafe = erralphafe
  ; Purpose: Calculate [Fe/H] and [alpha/Fe] from abundance structure
  ; Inputs:
  ; abund - Structure with abundances (element, abund, upperlimit fields)
  ; Outputs (via keywords):
  ; alphafe - [alpha/Fe] value
  ; errfeh - Error in [Fe/H]
  ; erralphafe - Error in [alpha/Fe]
  ; Returns: [Fe/H] value

  compile_opt idl2

  ; Solar abundance values (Asplund et al. 2009)
  eps_fe = 7.50

  ; Alpha elements used for [alpha/Fe] calculation
  ; Could also include: ['O', 'Ne', 'Mg', 'Si', 'S', 'Ar', 'Ca', 'Ti']
  ; with eps: [8.93, 8.09, 7.58, 7.55, 7.21, 6.56, 6.36, 4.99]
  alpha = ['Mg', 'Si', 'Ca']
  eps_alpha = [7.58, 7.55, 6.36]
  nalpha = n_elements(eps_alpha)

  ; Calculate [Fe/H] from Fe I lines only (excluding upper limits)
  w = where(abund.element eq 'Fe' and abund.upperlimit eq 0, cfe)
  feh = mean(abund[w].abund) - eps_fe
  errfeh = stddev(abund[w].abund) / sqrt(double(cfe))

  ; Calculate [alpha/Fe] from available alpha elements
  alphastarted = 0
  for i = 0, nalpha - 1 do begin
    w = where(abund.element eq alpha[i] and abund.upperlimit eq 0, c)
    if c gt 0 then begin
      alphaabundi = abund[w].abund - eps_alpha[i]
      alphaabund = alphastarted eq 0 ? alphaabundi : [alphaabund, alphaabundi]
      alphastarted = 1
    endif
  endfor

  ; Compute [alpha/Fe] if any alpha elements are available
  if alphastarted then begin
    alphafe = mean(alphaabund) - feh
    erralphafe = stddev(alphaabund) / sqrt(double(n_elements(alphaabund)))
  endif

  return, feh
end

function find_abund, par, star = star, chimask = chimask, fe1 = fe1, dfe1 = fe1err, fe2 = fe2, dfe2 = fe2err, ti1 = ti1, dti1 = ti1err, ti2 = ti2, dti2 = ti2err, feepslope = feepslope, dfeepslope = feepslopeerr, ferwslope = ferwslope, dferwslope = ferwslopeerr, tiepslope = tiepslope, dtiepslope = tiepslopeerr, tirwslope = tirwslope, dtirwslope = tirwslopeerr, teffphot = teffphot, verbose = verbose, dirflag = dirflag, fixedfeh = fixedfeh, fixedalphafe = fixedalphafe, feh = feh, errfeh = errfeh, alphafe = alphafe, erralphafe = erralphafe
  ; Purpose: Compute chi-squared for spectroscopic equilibrium diagnostics
  ; Used as objective function for MPFIT to find best-fit stellar parameters
  ; Enforces: excitation equilibrium, ionization equilibrium, RW independence
  ; Inputs:
  ; par - [Teff, logg, vt] parameter vector
  ; star - Star name
  ; chimask - Mask for chi-squared terms (1=use, 0=ignore)
  ; teffphot - If set, use photometric Teff (don't fit Teff)
  ; fixedfeh - If set, use this [Fe/H] instead of computing
  ; fixedalphafe - If set, use this [alpha/Fe] instead of computing
  ; Outputs (via keywords):
  ; fe1, fe2 - Mean Fe I and Fe II abundances
  ; fe1err, fe2err - Errors in Fe I and Fe II
  ; feepslope - Slope of Fe I vs. excitation potential (should be ~0)
  ; ferwslope - Slope of Fe I vs. reduced width (should be ~0)
  ; ti1, ti2 - Mean Ti I and Ti II abundances
  ; tiepslope, tirwslope - Ti slopes (similar to Fe)
  ; Returns: Chi-squared vector for equilibrium diagnostics

  compile_opt idl2
  common ews, ews
  if ~keyword_set(dirflag) then dirflag = ''

  ; Extract parameters
  teff = par[0]
  logg = par[1]
  vt = par[2]
  if ~keyword_set(chimask) then chimask = [1, 1, 1, 0, 0, 0]
  teffphot = keyword_set(teffphot) ? 1 : 0

  ; -----------------------------------------------------------------
  ; FIRST ITERATION: Run MOOG with initial metallicity guess
  ; -----------------------------------------------------------------
  interp_atm, teff, logg, vt, -2.41, 0.41, $
    outfile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.atm'
  spawn, 'MOOGSILENT ' + dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.par'
  abund = read_moog(dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.out2')

  ; Compute [Fe/H] and [alpha/Fe] from abundances
  feh = (feh_alphafe(abund, alphafe = alphafe) > (-4.0)) < 0.0
  if keyword_set(fixedfeh) then feh = fixedfeh
  if keyword_set(fixedalphafe) then alphafe = fixedalphafe
  alphafe >= -0.8
  alphafe <= 1.2

  ; -----------------------------------------------------------------
  ; SECOND ITERATION: Run MOOG with updated metallicity
  ; -----------------------------------------------------------------
  interp_atm, teff, logg, vt, feh, alphafe, $
    outfile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.atm'
  spawn, 'MOOGSILENT ' + dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.par'
  abund = read_moog(dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.out2')
  abund = struct_addtags(abund, replicate({ewerr: 0d, rwerr: 0d}, n_elements(abund)))
  for i = 0, n_elements(abund) - 1 do begin
    w = where(abs(ews.lambda - abund[i].lambda) lt (strtrim(abund[i].element, 2) eq 'Li' ? 0.10 : 0.35) and strtrim(strmid(ews.name, 0, 2), 2) eq strtrim(abund[i].element, 2) and abs(ews.ep - abund[i].ep) lt 0.1 and abs(ews.loggf - abund[i].loggf) lt 0.1, c)
    abund[i].ewerr = ews[w].ewerr
    abund[i].rwerr = abund[i].ewerr / abund[i].lambda
  endfor
  if dirflag eq '' then abund = abund[where(abund.ewerr gt 0)]
  ; Add EW errors to abundance structure
  abund = struct_addtags(abund, replicate({ewerr: 0d, rwerr: 0d}, n_elements(abund)))
  for i = 0, n_elements(abund) - 1 do begin
    w = where(abs(ews.lambda - abund[i].lambda) lt (strtrim(abund[i].element, 2) eq 'Li' ? 0.10 : 0.35) and $
      strtrim(strmid(ews.name, 0, 2), 2) eq strtrim(abund[i].element, 2) and $
      abs(ews.ep - abund[i].ep) lt 0.1 and abs(ews.loggf - abund[i].loggf) lt 0.1, c)
    abund[i].ewerr = ews[w].ewerr
    abund[i].rwerr = abund[i].ewerr / abund[i].lambda
  endfor
  if dirflag eq '' then abund = abund[where(abund.ewerr gt 0)]

  ; Recompute [Fe/H] and [alpha/Fe] with proper errors
  feh = feh_alphafe(abund, alphafe = alphafe, errfeh = errfeh, erralphafe = erralphafe)
  if keyword_set(fixedfeh) then feh = fixedfeh
  if keyword_set(fixedalphafe) then alphafe = fixedalphafe

  ; -----------------------------------------------------------------
  ; IRON EQUILIBRIUM DIAGNOSTICS
  ; -----------------------------------------------------------------

  ; Compute slope of Fe I abundance vs. reduced width (RW = log(EW/lambda))
  ; A non-zero slope indicates incorrect microturbulence
  w = where(abund.element eq 'Fe' and abund.ion eq 1 and $
    finite(abund.abund) and abund.upperlimit eq 0, c)
  a = linfit(alog10(abund[w].ew / abund[w].lambda), abund[w].abund, sigma = aerr, /double)
  ferwslope = a[1]

  ; Compute slope of Fe I abundance vs. excitation potential (EP)
  ; A non-zero slope indicates incorrect temperature
  a = linfit(abund[w].ep, abund[w].abund, sigma = aerr, /double)
  feepslope = a[1]

  ; Compute mean Fe I and Fe II abundances
  ; Ionization balance (Fe I - Fe II ~0) indicates correct log g
  w1 = where(abund.element eq 'Fe' and abund.ion eq 1 and $
    finite(abund.abund) and abund.upperlimit eq 0, c1)
  w2 = where(abund.element eq 'Fe' and abund.ion eq 2 and $
    finite(abund.abund) and abund.upperlimit eq 0, c2)
  fe1 = mean(abund[w1].abund)
  fe2 = mean(abund[w2].abund)
  dfe = fe1 - fe2

  ; -----------------------------------------------------------------
  ; JACKKNIFE ERROR ESTIMATION FOR IRON
  ; -----------------------------------------------------------------
  arwjk = dblarr(c1) ; Jackknife RW slopes
  aepjk = dblarr(c1) ; Jackknife EP slopes
  fe1jk = dblarr(c1) ; Jackknife FeI abundances
  fe2jk = dblarr(c2) ; Jackknife FeII abundances

  for i = 0, c1 - 1 do begin
    wjk = lindgen(c1)
    wjk = wjk[where(wjk ne i)] ; Leave out i-th line
    arwjk[i] = (linfit(alog10(abund[w1[wjk]].ew / abund[w1[wjk]].lambda), $
      abund[w1[wjk]].abund, /double))[1]
    aepjk[i] = (linfit(abund[w1[wjk]].ep, abund[w1[wjk]].abund, /double))[1]
    fe1jk[i] = mean(abund[w1[wjk]].abund)
  endfor

  for i = 0, c2 - 1 do begin
    wjk = lindgen(c2)
    wjk = wjk[where(wjk ne i)] ; Leave out i-th line
    fe2jk[i] = mean(abund[w2[wjk]].abund)
  endfor

  ; Calculate jackknife standard errors
  ferwslopeerr = sqrt(double(c1 - 1) / double(c1) * total((arwjk - ferwslope) ^ 2.0))
  feepslopeerr = sqrt(double(c1 - 1) / double(c1) * total((aepjk - feepslope) ^ 2.0))
  fe1err = sqrt(double(c1 - 1) / double(c1) * total((fe1jk - fe1) ^ 2.0))
  fe2err = sqrt(double(c2 - 1) / double(c2) * total((fe2jk - fe2) ^ 2.0))
  dfeerr = sqrt(fe1err ^ 2. + fe2err ^ 2.)

  ; -----------------------------------------------------------------
  ; TITANIUM DIAGNOSTICS
  ; -----------------------------------------------------------------
  ; Compute Ti I excitation potential slope and reduced width slope
  w = where(abund.element eq 'Ti' and abund.ion eq 1 and $
    finite(abund.abund) and abund.upperlimit eq 0)

  ; Reduced width slope (should be ~0 for correct v_t)
  a = linfit(alog10(abund[w].ew / abund[w].lambda), abund[w].abund, sigma = aerr, /double)
  tirwslope = a[1]

  ; EP slope (should be ~0 for correct Teff)
  a = linfit(abund[w].ep, abund[w].abund, sigma = aerr, /double)
  tiepslope = a[1]

  ; Compute mean Ti I and Ti II abundances
  ; Ionization balance (Ti I - Ti II ~0) indicates correct log g
  w1 = where(abund.element eq 'Ti' and abund.ion eq 1 and $
    finite(abund.abund) and abund.upperlimit eq 0, c1)
  w2 = where(abund.element eq 'Ti' and abund.ion eq 2 and $
    finite(abund.abund) and abund.upperlimit eq 0, c2)
  ti1 = mean(abund[w1].abund)
  ti2 = mean(abund[w2].abund)
  dti = ti1 - ti2

  ; -----------------------------------------------------------------
  ; JACKKNIFE ERROR ESTIMATION FOR TITANIUM
  ; -----------------------------------------------------------------
  if c1 gt 0 and c2 gt 0 then begin
    arwjk = dblarr(c1) ; Jackknife RW slopes
    aepjk = dblarr(c1) ; Jackknife EP slopes
    ti1jk = dblarr(c1) ; Jackknife TiI abundances
    ti2jk = dblarr(c2) ; Jackknife TiII abundances

    for i = 0, c1 - 1 do begin
      wjk = lindgen(c1)
      wjk = wjk[where(wjk ne i)] ; Leave out i-th line
      arwjk[i] = (linfit(alog10(abund[w1[wjk]].ew / abund[w1[wjk]].lambda), $
        abund[w1[wjk]].abund, /double))[1]
      aepjk[i] = (linfit(abund[w1[wjk]].ep, abund[w1[wjk]].abund, /double))[1]
      ti1jk[i] = mean(abund[w1[wjk]].abund)
    endfor

    for i = 0, c2 - 1 do begin
      wjk = lindgen(c2)
      wjk = wjk[where(wjk ne i)] ; Leave out i-th line
      ti2jk[i] = mean(abund[w2[wjk]].abund)
    endfor

    ; Calculate jackknife standard errors
    tirwslopeerr = sqrt(double(c1 - 1) / double(c1) * total((arwjk - tirwslope) ^ 2.0))
    tiepslopeerr = sqrt(double(c1 - 1) / double(c1) * total((aepjk - tiepslope) ^ 2.0))
    ti1err = sqrt(double(c1 - 1) / double(c1) * total((ti1jk - ti1) ^ 2.0))
    ti2err = sqrt(double(c2 - 1) / double(c2) * total((ti2jk - ti2) ^ 2.0))
    dtierr = sqrt(ti1err ^ 2. + ti2err ^ 2.)
  endif

  ; Handle cases with too few lines
  if c1 le 2 then begin
    tirwslope = 0.0
    tirwslopeerr = 1.0
    tiepslope = 0.0
    tiepslopeerr = 1.0
  endif
  if c1 le 2 or c2 eq 0 then begin
    dti = 0.0
    dtierr = 1.0
  endif

  ; -----------------------------------------------------------------
  ; VERBOSE OUTPUT (optional)
  ; -----------------------------------------------------------------
  if keyword_set(verbose) then begin
    print, ' dFe/dEP = ' + string(feepslope, format = '(D+5.2)') + ' +/- ' + string(feepslopeerr, format = '(D4.2)')
    print, ' dFe/dRW = ' + string(ferwslope, format = '(D+5.2)') + ' +/- ' + string(ferwslopeerr, format = '(D4.2)')
    print, 'FeI-FeII = ' + string(dfe, format = '(D+5.2)') + ' +/- ' + string(dfeerr, format = '(D4.2)')
    print, ' dTi/dEP = ' + string(tiepslope, format = '(D+5.2)') + ' +/- ' + string(tiepslopeerr, format = '(D4.2)')
    print, ' dTi/dRW = ' + string(tirwslope, format = '(D+5.2)') + ' +/- ' + string(tirwslopeerr, format = '(D4.2)')
    print, 'TiI-TiII = ' + string(dti, format = '(D+5.2)') + ' +/- ' + string(dtierr, format = '(D4.2)')
  endif

  chi = double(chimask) * [feepslope / feepslopeerr, dfe / dfeerr, ferwslope / ferwslopeerr, tiepslope / tiepslopeerr, dti / dtierr, tirwslope / tirwslopeerr]
  return, chi
end

; =================================================================
; FUNCTION: find_abund_final
; =================================================================
; PURPOSE:
; Calculate equilibrium diagnostics for final stellar parameter
; determination. Similar to find_abund but uses ATLAS model
; atmospheres directly (for /final mode).
;
; INPUTS:
; par      - [Teff, log(g), v_t] parameter array
; star     - Star name
; chimask  - Mask for chi-squared terms [dFe/dEP, FeI-FeII, dFe/dRW, dTi/dEP, TiI-TiII, dTi/dRW]
; teffphot - If set, use photometric Teff
; dirflag  - Directory flag for file paths
;
; OUTPUTS:
; chi - Array of normalized equilibrium diagnostics
;
; METHOD:
; Uses ATLAS model atmospheres from /raid/atlas/BasicATLAS/ATLAS_LMHA_{star}/.
; Otherwise identical to find_abund function.
;
; =================================================================
function find_abund_final, par, star = star, chimask = chimask, $
  fe1 = fe1, dfe1 = fe1err, fe2 = fe2, dfe2 = fe2err, $
  ti1 = ti1, dti1 = ti1err, ti2 = ti2, dti2 = ti2err, $
  feepslope = feepslope, dfeepslope = feepslopeerr, $
  ferwslope = ferwslope, dferwslope = ferwslopeerr, $
  tiepslope = tiepslope, dtiepslope = tiepslopeerr, $
  tirwslope = tirwslope, dtirwslope = tirwslopeerr, $
  teffphot = teffphot, dirflag = dirflag
  compile_opt idl2
  common ews, ews
  if ~keyword_set(dirflag) then dirflag = ''

  ; Extract parameters
  teff = par[0]
  logg = par[1]
  vt = par[2]
  if ~keyword_set(chimask) then chimask = [1, 1, 1, 0, 0, 0]
  teffphot = keyword_set(teffphot) ? 1 : 0

  ; -----------------------------------------------------------------
  ; CONVERT ATLAS ATMOSPHERE AND RUN MOOG
  ; -----------------------------------------------------------------
  atlas_to_moog, '/raid/atlas/BasicATLAS/ATLAS_LMHA_' + star + '/output_summary.out', $
    dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.atm', vt = vt
  spawn, 'MOOGSILENT ' + dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.par'
  abund = read_moog(dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.out2')

  ; Add EW error fields
  abund = struct_addtags(abund, replicate({ewerr: 0d, rwerr: 0d}, n_elements(abund)))
  for i = 0, n_elements(abund) - 1 do begin
    w = where(abs(ews.lambda - abund[i].lambda) lt $
      (strtrim(abund[i].element, 2) eq 'Li' ? 0.10 : 0.35) and $
      strtrim(strmid(ews.name, 0, 2), 2) eq strtrim(abund[i].element, 2) and $
      abs(ews.ep - abund[i].ep) lt 0.1 and abs(ews.loggf - abund[i].loggf) lt 0.1, c)
    abund[i].ewerr = ews[w].ewerr
    abund[i].rwerr = abund[i].ewerr / abund[i].lambda
  endfor
  if dirflag eq '' then abund = abund[where(abund.ewerr gt 0)]

  ; -----------------------------------------------------------------
  ; IRON EQUILIBRIUM DIAGNOSTICS
  ; -----------------------------------------------------------------
  ; Select Fe I lines for trend analysis
  w = where(abund.element eq 'Fe' and abund.ion eq 1 and $
    finite(abund.abund) and abund.upperlimit eq 0, c)

  ; Reduced width slope (should be ~0 for correct v_t)
  a = linfit(alog10(abund[w].ew / abund[w].lambda), abund[w].abund, sigma = aerr, /double)
  ferwslope = a[1]

  ; Excitation potential slope (should be ~0 for correct Teff)
  a = linfit(abund[w].ep, abund[w].abund, sigma = aerr, /double)
  feepslope = a[1]

  ; Compute mean Fe I and Fe II abundances
  w1 = where(abund.element eq 'Fe' and abund.ion eq 1 and $
    finite(abund.abund) and abund.upperlimit eq 0, c1)
  w2 = where(abund.element eq 'Fe' and abund.ion eq 2 and $
    finite(abund.abund) and abund.upperlimit eq 0, c2)
  fe1 = mean(abund[w1].abund)
  fe2 = mean(abund[w2].abund)
  dfe = fe1 - fe2 ; Ionization balance (should be ~0 for correct log g)

  ; Jackknife error estimation
  arwjk = dblarr(c1) ; Jackknife RW slopes
  aepjk = dblarr(c1) ; Jackknife EP slopes
  fe1jk = dblarr(c1) ; Jackknife FeI abundances
  fe2jk = dblarr(c2) ; Jackknife FeII abundances

  for i = 0, c1 - 1 do begin
    wjk = lindgen(c1)
    wjk = wjk[where(wjk ne i)] ; Leave out i-th line
    arwjk[i] = (linfit(alog10(abund[w1[wjk]].ew / abund[w1[wjk]].lambda), $
      abund[w1[wjk]].abund, /double))[1]
    aepjk[i] = (linfit(abund[w1[wjk]].ep, abund[w1[wjk]].abund, /double))[1]
    fe1jk[i] = mean(abund[w1[wjk]].abund)
  endfor

  for i = 0, c2 - 1 do begin
    wjk = lindgen(c2)
    wjk = wjk[where(wjk ne i)] ; Leave out i-th line
    fe2jk[i] = mean(abund[w2[wjk]].abund)
  endfor

  ; Calculate jackknife standard errors
  ferwslopeerr = sqrt(double(c1 - 1) / double(c1) * total((arwjk - ferwslope) ^ 2.0))
  feepslopeerr = sqrt(double(c1 - 1) / double(c1) * total((aepjk - feepslope) ^ 2.0))
  fe1err = sqrt(double(c1 - 1) / double(c1) * total((fe1jk - fe1) ^ 2.0))
  fe2err = sqrt(double(c2 - 1) / double(c2) * total((fe2jk - fe2) ^ 2.0))
  dfeerr = sqrt(fe1err ^ 2. + fe2err ^ 2.)

  ; -----------------------------------------------------------------
  ; TITANIUM EQUILIBRIUM DIAGNOSTICS
  ; -----------------------------------------------------------------
  ; Select Ti I lines for trend analysis
  w = where(abund.element eq 'Ti' and abund.ion eq 1 and $
    finite(abund.abund) and abund.upperlimit eq 0)

  ; Reduced width slope
  a = linfit(alog10(abund[w].ew / abund[w].lambda), abund[w].abund, sigma = aerr, /double)
  tirwslope = a[1]

  ; Excitation potential slope
  a = linfit(abund[w].ep, abund[w].abund, sigma = aerr, /double)
  tiepslope = a[1]

  ; Compute mean Ti I and Ti II abundances
  w1 = where(abund.element eq 'Ti' and abund.ion eq 1 and $
    finite(abund.abund) and abund.upperlimit eq 0, c1)
  w2 = where(abund.element eq 'Ti' and abund.ion eq 2 and $
    finite(abund.abund) and abund.upperlimit eq 0, c2)
  ti1 = mean(abund[w1].abund)
  ti2 = mean(abund[w2].abund)
  dti = ti1 - ti2 ; Ionization balance

  ; Jackknife error estimation for titanium
  if c1 gt 0 and c2 gt 0 then begin
    arwjk = dblarr(c1)
    aepjk = dblarr(c1)
    ti1jk = dblarr(c1)
    ti2jk = dblarr(c2)

    for i = 0, c1 - 1 do begin
      wjk = lindgen(c1)
      wjk = wjk[where(wjk ne i)]
      arwjk[i] = (linfit(alog10(abund[w1[wjk]].ew / abund[w1[wjk]].lambda), $
        abund[w1[wjk]].abund, /double))[1]
      aepjk[i] = (linfit(abund[w1[wjk]].ep, abund[w1[wjk]].abund, /double))[1]
      ti1jk[i] = mean(abund[w1[wjk]].abund)
    endfor

    for i = 0, c2 - 1 do begin
      wjk = lindgen(c2)
      wjk = wjk[where(wjk ne i)]
      ti2jk[i] = mean(abund[w2[wjk]].abund)
    endfor

    tirwslopeerr = sqrt(double(c1 - 1) / double(c1) * total((arwjk - tirwslope) ^ 2.0))
    tiepslopeerr = sqrt(double(c1 - 1) / double(c1) * total((aepjk - tiepslope) ^ 2.0))
    ti1err = sqrt(double(c1 - 1) / double(c1) * total((ti1jk - ti1) ^ 2.0))
    ti2err = sqrt(double(c2 - 1) / double(c2) * total((ti2jk - ti2) ^ 2.0))
    dtierr = sqrt(ti1err ^ 2. + ti2err ^ 2.)
  endif

  ; Handle cases with insufficient lines
  if c1 le 2 then begin
    tirwslope = 0.0
    tirwslopeerr = 1.0
    tiepslope = 0.0
    tiepslopeerr = 1.0
  endif
  if c1 le 2 or c2 eq 0 then begin
    dti = 0.0
    dtierr = 1.0
  endif

  ; -----------------------------------------------------------------
  ; COMPUTE CHI-SQUARED VECTOR
  ; -----------------------------------------------------------------
  ; Normalize each diagnostic by its error and apply mask
  chi = double(chimask) * [feepslope / feepslopeerr, $
    dfe / dfeerr, $
    ferwslope / ferwslopeerr, $
    tiepslope / tiepslopeerr, $
    dti / dtierr, $
    tirwslope / tirwslopeerr]
  return, chi
end

; =================================================================
; FUNCTION: calculate_abund
; =================================================================
; PURPOSE:
; Calculate abundances for all elements using MOOG with
; hyperfine structure treatment for selected elements.
;
; INPUTS:
; star     - Star name
; teffphot - If set, use photometric temperature mode
; dirflag  - Directory for EW files and output
; uperr    - If set, use upper EW error (EW + error)
; downerr  - If set, use lower EW error (EW - error)
;
; OUTPUTS:
; Structure array with abundances for each line:
; .lambda, .element, .ion, .ep, .loggf, .abund, .ew
; .ewerr, .rwerr (EW errors)
; .upperlimit (flag for upper limit lines)
;
; METHOD:
; 1. Process normal lines (non-hyperfine) with MOOG 'abfind'
; 2. Process hyperfine elements (Sc,V,Mn,Co,Cu,Ba,La,Nd,Eu) with 'blends'
; 3. Process upper limits separately
; 4. Match results with EW errors from input file
;
; =================================================================
function calculate_abund, star, teffphot = teffphot, dirflag = dirflag, $
  uperr = uperr, downerr = downerr
  compile_opt idl2
  common ews, ews
  dirflag2 = getenv('CALTECH') + 'hires/M15_M92/' + dirflag
  e = elements(/newmoog)

  ; -----------------------------------------------------------------
  ; HYPERFINE STRUCTURE ELEMENTS
  ; -----------------------------------------------------------------
  ; Elements requiring HFS treatment with isotope fractions
  hyperfine_el = [21, 23, 25, 27, 29, 56, 57, 60, 63] ; Atomic numbers
  hyperfine_name = ['sc', 'v', 'mn', 'co', 'cu', 'ba', 'la', 'nd', 'eu']
  nhyper = n_elements(hyperfine_el)
  nlineshyper = intarr(nhyper) ; Count of HFS detections per element
  nlineshyperul = intarr(nhyper) ; Count of HFS upper limits per element

  ; -----------------------------------------------------------------
  ; PROCESS NORMAL LINES (non-hyperfine elements)
  ; -----------------------------------------------------------------
  nlines = trim_linelist(star, /nothyperfine, $
    teffphot = teffphot, dirflag = dirflag, $
    uperr = uperr, downerr = downerr)
  make_par, parfile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.par', $
    atmfile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.atm', $
    driver = 'abfind', $
    linefile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '_temp.ew', $
    outfile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.out2'
  spawn, 'MOOGSILENT ' + dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.par'

  ; -----------------------------------------------------------------
  ; PROCESS HYPERFINE STRUCTURE ELEMENTS (detections)
  ; -----------------------------------------------------------------
  ; Use MOOG 'blends' driver for hyperfine structure
  for i = 0, nhyper - 1 do begin
    nlineshyper[i] = trim_linelist(star, atomic = hyperfine_el[i], $
      teffphot = teffphot, dirflag = dirflag, $
      uperr = uperr, downerr = downerr)
    if nlineshyper[i] gt 0 then begin
      ; Get isotope information for this element
      wel = where(e.atomic eq hyperfine_el[i])
      wiso = where(e[wel].isotope gt 0, niso)

      ; Atomic number with ionization state (add 0.1 for certain elements)
      atomic = hyperfine_el[i] + ((hyperfine_el[i] eq 21 or hyperfine_el[i] gt 30) ? 0.1 : 0.0)
      isotopes = e[wel].isotope[wiso]

      ; Use r-process fractions for heavy elements (Z≥31), solar otherwise
      isofracs = hyperfine_el[i] ge 31 ? e[wel].rfrac[wiso] : e[wel].solarfrac[wiso]

      ; Create parameter file with isotope information (if multiple isotopes)
      if niso gt 1 and hyperfine_el[i] ne 38 then begin
        make_par, parfile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.par', $
          atmfile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.atm', $
          driver = 'blends', $
          linefile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '_temp.ew', $
          outfile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '_' + hyperfine_name[i] + '.out2', $
          atomic = atomic, isotopes = isotopes, isofracs = isofracs
      endif else begin
        make_par, parfile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.par', $
          atmfile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.atm', $
          driver = 'blends', $
          linefile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '_temp.ew', $
          outfile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '_' + hyperfine_name[i] + '.out2', $
          atomic = atomic
      endelse
      spawn, 'MOOGSILENT ' + dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.par'
    endif
  endfor

  ; -----------------------------------------------------------------
  ; PROCESS UPPER LIMITS (normal elements)
  ; -----------------------------------------------------------------
  nlinesul = trim_linelist(star, /nothyperfine, teffphot = teffphot, $
    /upperlimit, dirflag = dirflag, $
    uperr = uperr, downerr = downerr)
  if nlinesul gt 0 then begin
    make_par, parfile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.par', $
      atmfile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.atm', $
      driver = 'abfind', $
      linefile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '_temp.ew', $
      outfile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '_ul.out2'
    spawn, 'MOOGSILENT ' + dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.par'
  endif

  ; -----------------------------------------------------------------
  ; PROCESS UPPER LIMITS (hyperfine elements)
  ; -----------------------------------------------------------------
  for i = 0, nhyper - 1 do begin
    nlineshyperul[i] = trim_linelist(star, atomic = hyperfine_el[i], $
      teffphot = teffphot, /upperlimit, $
      dirflag = dirflag, uperr = uperr, downerr = downerr)
    if nlineshyperul[i] gt 0 then begin
      wel = where(e.atomic eq hyperfine_el[i])
      wiso = where(e[wel].isotope gt 0, niso)
      atomic = hyperfine_el[i] + ((hyperfine_el[i] eq 21 or hyperfine_el[i] gt 30) ? 0.1 : 0.0)
      isotopes = e[wel].isotope[wiso]
      isofracs = hyperfine_el[i] ge 31 ? e[wel].rfrac[wiso] : e[wel].solarfrac[wiso]

      if niso gt 1 and hyperfine_el[i] ne 38 then begin
        make_par, parfile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.par', $
          atmfile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.atm', $
          driver = 'blends', $
          linefile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '_temp.ew', $
          outfile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '_' + hyperfine_name[i] + '_ul.out2', $
          atomic = atomic, isotopes = isotopes, isofracs = isofracs
      endif else begin
        make_par, parfile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.par', $
          atmfile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.atm', $
          driver = 'blends', $
          linefile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '_temp.ew', $
          outfile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '_' + hyperfine_name[i] + '_ul.out2', $
          atomic = atomic
      endelse
      spawn, 'MOOGSILENT ' + dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.par'
    endif
  endfor

  ; -----------------------------------------------------------------
  ; COMBINE ALL MOOG OUTPUTS
  ; -----------------------------------------------------------------
  ; Read normal lines
  abund = read_moog(dirflag2 + star + (teffphot eq 0 ? '' : '_teffphot') + '.out2')

  ; Append hyperfine structure results
  for i = 0, nhyper - 1 do begin
    if nlineshyper[i] gt 0 then begin
      abund = struct_append(abund, $
        read_moog(dirflag2 + star + (teffphot eq 0 ? '' : '_teffphot') + $
          '_' + hyperfine_name[i] + '.out2', /blends))
    endif
  endfor

  ; Append upper limit results (normal lines)
  if nlinesul gt 0 then begin
    abund = struct_append(abund, $
      read_moog(dirflag2 + star + (teffphot eq 0 ? '' : '_teffphot') + '_ul.out2', $
        /upperlimit))
  endif

  ; Append upper limit results (hyperfine lines)
  for i = 0, nhyper - 1 do begin
    if nlineshyperul[i] gt 0 then begin
      abund = struct_append(abund, $
        read_moog(dirflag2 + star + (teffphot eq 0 ? '' : '_teffphot') + $
          '_' + hyperfine_name[i] + '_ul.out2', $
          /upperlimit, /blends))
    endif
  endfor

  ; Filter out invalid abundances (MOOG returns 100 for failed lines)
  abund = abund[where(abund.abund lt 100 and finite(abund.abund))]

  ; -----------------------------------------------------------------
  ; MATCH WITH EW ERRORS
  ; -----------------------------------------------------------------
  ; Add error fields to output structure
  abund = struct_addtags(abund, replicate({ewerr: 0d, rwerr: 0d}, n_elements(abund)))

  ; Match each line with its EW error from the input file
  for i = 0, n_elements(abund) - 1 do begin
    if contains(hyperfine_name, strlowcase(abund[i].element)) then begin
      ; Hyperfine lines: broader wavelength tolerance (1.0 Å)
      w = where(abs(ews.lambda - abund[i].lambda) lt 1.0 and $
        strtrim(strmid(ews.name, 0, 2), 2) eq strtrim(abund[i].element, 2))
      abund[i].lambda = ews[w].lambda
      abund[i].ep = ews[w].ep
      abund[i].loggf = ews[w].loggf
      abund[i].ewerr = ews[w].ewerr
      abund[i].rwerr = abund[i].ewerr / abund[i].lambda
    endif else begin
      ; Normal lines: tighter tolerances (0.35 Å, except Li: 0.10 Å)
      w = where(abs(ews.lambda - abund[i].lambda) lt $
        (strtrim(abund[i].element, 2) eq 'Li' ? 0.10 : 0.35) and $
        strtrim(strmid(ews.name, 0, 2), 2) eq strtrim(abund[i].element, 2) and $
        abs(ews.ep - abund[i].ep) lt 0.1 and abs(ews.loggf - abund[i].loggf) lt 0.1)
      abund[i].ewerr = ews[w].ewerr
      abund[i].rwerr = abund[i].ewerr / abund[i].lambda
    endelse
  endfor

  return, abund
end

function error_analysis, star, abunds, teffphot = teffphot
  compile_opt idl2
  common ews, ews
  teffphot = keyword_set(teffphot) ? 1 : 0
  jlcspecies = [3.0, 8.0, 11.0, 12.0, 13.0, 14.0, 19.0, 20.0, 20.1, 21.1, 22.0, 22.1, 23.0, 23.1, 24.0, 24.1, 25.0, 25.1, 26.0, 26.1, 27.0, 28.0, 29.0, 30.0, 38.0, 38.1, 39.1, 40.0, 40.1, 56.1, 57.1, 58.1, 59.1, 60.1, 62.1, 63.1, 64.1, 66.1, 67.1, 82.0, 6.0, 7.0]
  njlcspecies = n_elements(jlcspecies)
  e = elements(/newmoog)

  hiresall = mrdfits('M15_M92_allframes.fits', 1, /silent)
  hiresall = hiresall[where(hiresall.gmag0 lt 18.5 and strtrim(hiresall.name, 2) ne 'X-20' and strtrim(hiresall.name, 2) ne 'S2303')]

  dirflag = 'ew/'
  dirflag2 = getenv('CALTECH') + 'hires/M15_M92/' + dirflag
  ews = mrdfits(dirflag2 + star + '_Ji20_ew.fits', 1, /silent)

  interp_atm, abunds.teff, abunds.logg, abunds.vt, abunds.feh, abunds.alphafe, outfile = dirflag2 + star + (teffphot eq 0 ? '' : '_teffphot') + '.atm'
  abund = calculate_abund(star, teffphot = teffphot, dirflag = dirflag)
  newstr = {uperr: -999d, downerr: -999d, abunderr: -999d, upteff: -999d, uplogg: -999d, upvt: -999d, upfeh: -999d, upalphafe: -999d, weight: 0d, delta: dblarr(4)}
  na = n_elements(abund)
  abund = struct_addtags(abund, replicate(newstr, na))

  abund_uperr = calculate_abund(star, teffphot = teffphot, dirflag = dirflag, /uperr)
  abund_downerr = calculate_abund(star, teffphot = teffphot, dirflag = dirflag, /downerr)
  match, abund.lambda, abund_uperr.lambda, w1, w2
  abund[w1].uperr = abund_uperr[w2].abund - abund[w1].abund
  match, abund.lambda, abund_downerr.lambda, w1, w2
  abund[w1].downerr = abund_downerr[w2].abund - abund[w1].abund
  w = where(abund.uperr le 0 and abund.downerr lt 0 and abund.downerr gt -10, c)
  if c gt 0 then abund[w].uperr = abund[w].downerr
  w = where(abund.uperr gt 0 and abund.downerr ge 0, c)
  if c gt 0 then abund[w].downerr = abund[w].uperr
  w = where(abund.uperr le 0 or abund.downerr ge 0 or abund.downerr le -10, c)
  if c gt 0 then begin
    abund[w].uperr = 0.1
    abund[w].downerr = -0.1
  endif
  abund.abunderr = (abund.uperr - abund.downerr) / 2.

  interp_atm, abunds.teff + abunds.tefferr, abunds.logg, abunds.vt, abunds.feh, abunds.alphafe, outfile = dirflag2 + star + (teffphot eq 0 ? '' : '_teffphot') + '.atm'
  abund_upteff = calculate_abund(star, teffphot = teffphot, dirflag = dirflag)
  match, abund.lambda, abund_upteff.lambda, w1, w2
  abund[w1].upteff = abund_upteff[w2].abund - abund[w1].abund

  interp_atm, abunds.teff, abunds.logg + abunds.loggerr, abunds.vt, abunds.feh, abunds.alphafe, outfile = dirflag2 + star + (teffphot eq 0 ? '' : '_teffphot') + '.atm'
  abund_uplogg = calculate_abund(star, teffphot = teffphot, dirflag = dirflag)
  match, abund.lambda, abund_uplogg.lambda, w1, w2
  abund[w1].uplogg = abund_uplogg[w2].abund - abund[w1].abund

  interp_atm, abunds.teff, abunds.logg, abunds.vt + abunds.vterr, abunds.feh, abunds.alphafe, outfile = dirflag2 + star + (teffphot eq 0 ? '' : '_teffphot') + '.atm'
  abund_upvt = calculate_abund(star, teffphot = teffphot, dirflag = dirflag)
  match, abund.lambda, abund_upvt.lambda, w1, w2
  abund[w1].upvt = abund_upvt[w2].abund - abund[w1].abund

  interp_atm, abunds.teff, abunds.logg, abunds.vt, abunds.feh + abunds.feherr, abunds.alphafe, outfile = dirflag2 + star + (teffphot eq 0 ? '' : '_teffphot') + '.atm'
  abund_upfeh = calculate_abund(star, teffphot = teffphot, dirflag = dirflag)
  match, abund.lambda, abund_upfeh.lambda, w1, w2
  abund[w1].upfeh = abund_upfeh[w2].abund - abund[w1].abund

  return, abund
end

function error_analysis_final, star, abunds, teffphot = teffphot
  compile_opt idl2
  common ews, ews
  teffphot = keyword_set(teffphot) ? 1 : 0
  jlcspecies = [3.0, 8.0, 11.0, 12.0, 13.0, 14.0, 19.0, 20.0, 20.1, 21.1, 22.0, 22.1, 23.0, 23.1, 24.0, 24.1, 25.0, 25.1, 26.0, 26.1, 27.0, 28.0, 29.0, 30.0, 38.0, 38.1, 39.1, 40.0, 40.1, 56.1, 57.1, 58.1, 59.1, 60.1, 62.1, 63.1, 64.1, 66.1, 67.1, 82.0, 6.0, 7.0]
  njlcspecies = n_elements(jlcspecies)
  e = elements(/newmoog)

  hiresall = mrdfits('M15_M92_allframes.fits', 1, /silent)
  hiresall = hiresall[where(hiresall.gmag0 lt 18.5 and strtrim(hiresall.name, 2) ne 'X-20' and strtrim(hiresall.name, 2) ne 'S2303')]

  dirflag = 'ew2/'
  dirflag2 = getenv('CALTECH') + 'hires/M15_M92/' + dirflag
  ews = mrdfits(dirflag2 + star + '_Ji20_ew.fits', 1, /silent)

  atlas_to_moog, '/raid/atlas/BasicATLAS/ATLAS_LMHA_' + star + '/output_summary.out', dirflag2 + star + (teffphot eq 0 ? '' : '_teffphot') + '.atm', vt = abunds.vt
  abund = calculate_abund(star, teffphot = teffphot, dirflag = dirflag)
  newstr = {uperr: -999d, downerr: -999d, abunderr: -999d, upteff: -999d, uplogg: -999d, upvt: -999d, upfeh: -999d, upalphafe: -999d, weight: 0d, delta: dblarr(4)}
  na = n_elements(abund)
  abund = struct_addtags(abund, replicate(newstr, na))

  abund_uperr = calculate_abund(star, teffphot = teffphot, dirflag = dirflag, /uperr)
  abund_downerr = calculate_abund(star, teffphot = teffphot, dirflag = dirflag, /downerr)
  match, abund.lambda, abund_uperr.lambda, w1, w2
  abund[w1].uperr = abund_uperr[w2].abund - abund[w1].abund
  match, abund.lambda, abund_downerr.lambda, w1, w2
  abund[w1].downerr = abund_downerr[w2].abund - abund[w1].abund
  w = where(abund.uperr le 0 and abund.downerr lt 0 and abund.downerr gt -10, c)
  if c gt 0 then abund[w].uperr = abund[w].downerr
  w = where(abund.uperr gt 0 and abund.downerr ge 0, c)
  if c gt 0 then abund[w].downerr = abund[w].uperr
  w = where(abund.uperr le 0 or abund.downerr ge 0 or abund.downerr le -10, c)
  if c gt 0 then begin
    abund[w].uperr = 0.1
    abund[w].downerr = -0.1
  endif
  abund.abunderr = (abund.uperr - abund.downerr) / 2.

  atlas_to_moog, '/raid/atlas/BasicATLAS/ATLAS_LMHA_' + star + '_upteff/output_summary.out', dirflag2 + star + (teffphot eq 0 ? '' : '_teffphot') + '.atm', vt = abunds.vt
  abund_upteff = calculate_abund(star, teffphot = teffphot, dirflag = dirflag)
  match, abund.lambda, abund_upteff.lambda, w1, w2
  abund[w1].upteff = abund_upteff[w2].abund - abund[w1].abund

  atlas_to_moog, '/raid/atlas/BasicATLAS/ATLAS_LMHA_' + star + '_uplogg/output_summary.out', dirflag2 + star + (teffphot eq 0 ? '' : '_teffphot') + '.atm', vt = abunds.vt
  abund_uplogg = calculate_abund(star, teffphot = teffphot, dirflag = dirflag)
  match, abund.lambda, abund_uplogg.lambda, w1, w2
  abund[w1].uplogg = abund_uplogg[w2].abund - abund[w1].abund

  atlas_to_moog, '/raid/atlas/BasicATLAS/ATLAS_LMHA_' + star + '/output_summary.out', dirflag2 + star + (teffphot eq 0 ? '' : '_teffphot') + '.atm', vt = abunds.vt + abunds.vterr
  abund_upvt = calculate_abund(star, teffphot = teffphot, dirflag = dirflag)
  match, abund.lambda, abund_upvt.lambda, w1, w2
  abund[w1].upvt = abund_upvt[w2].abund - abund[w1].abund

  atlas_to_moog, '/raid/atlas/BasicATLAS/ATLAS_LMHA_' + star + '_upfeh/output_summary.out', dirflag2 + star + (teffphot eq 0 ? '' : '_teffphot') + '.atm', vt = abunds.vt
  abund_upfeh = calculate_abund(star, teffphot = teffphot, dirflag = dirflag)
  match, abund.lambda, abund_upfeh.lambda, w1, w2
  abund[w1].upfeh = abund_upfeh[w2].abund - abund[w1].abund
  return, abund
end

; =================================================================
; PROCEDURE: hires_loop
; =================================================================
; PURPOSE:
; Main analysis loop for determining stellar parameters and abundances
; from HIRES spectra. Optimizes Teff, log(g), and v_t by enforcing:
; 1. Excitation balance (d[Fe/H]/dEP = 0)
; 2. Ionization balance (FeI - FeII = 0)
; 3. No trend with line strength (d[Fe/H]/dRW = 0)
;
; INPUTS:
; star      - Star name (e.g., 'M15-star-10')
; teffphot  - If set, use photometric Teff (fix Teff, optimize log g & v_t)
; fixedmet  - If set, use fixed metallicity instead of Fe lines
; final     - If set, use final EW directory ('ew2/' instead of 'ew/')
;
; OUTPUTS:
; Creates FITS files with abundances and errors:
; - {star}_abund_{mode}.fits        : Line-by-line abundances
; - {star}_abundbyline_{mode}.fits  : Same with all metadata
;
; METHOD:
; Uses MPFIT to minimize chi-squared from equilibrium diagnostics.
; Iterates stellar parameters until Fe I/II slopes and difference converge.
;
; =================================================================
pro hires_loop, star, teffphot = teffphot, fixedmet = fixedmet, final = final
  compile_opt idl2
  common ews, ews

  ; -----------------------------------------------------------------
  ; SETUP DIRECTORIES AND DATA
  ; -----------------------------------------------------------------
  dirflag = keyword_set(final) ? 'ew2/' : 'ew/'
  dirflag2 = getenv('CALTECH') + 'hires/M15_M92/' + dirflag
  dirflag3 = getenv('CALTECH') + 'hires/M15_M92/ew/'
  ews = mrdfits(dirflag3 + star + '_Ji20_ew.fits', 1, /silent)
  allframes = mrdfits(getenv('CALTECH') + 'hires/M15_M92/M15_M92_allframes.fits', 1, /silent)
  w = where(strtrim(allframes.name, 2) eq star)

  ; -----------------------------------------------------------------
  ; INITIAL STELLAR PARAMETERS
  ; -----------------------------------------------------------------
  teff = allframes[w].teffphot ; Photometric temperature
  tefferr = allframes[w].teffphoterr
  logg = allframes[w].loggphot ; Photometric gravity
  loggerr = allframes[w].loggphoterr
  feh = -2.41 ; Cluster metallicity
  alphafe = 0.41 ; Cluster alpha enhancement
  vt = 2.13 - 0.23 * logg ; Initial microturbulence estimate

  teffphot = keyword_set(teffphot) ? 1 : 0

  seed = 785122l
  e = elements(/newmoog)

  ; List of all elements to analyze
  jlcspecies = [3.0, 8.0, 11.0, 12.0, 13.0, 14.0, 19.0, 20.0, 20.1, 21.1, 22.0, 22.1, $
    23.0, 23.1, 24.0, 24.1, 25.0, 25.1, 26.0, 26.1, 27.0, 28.0, 29.0, 30.0, $
    38.0, 38.1, 39.1, 40.0, 40.1, 56.1, 57.1, 58.1, 59.1, 60.1, 62.1, 63.1, $
    64.1, 66.1, 67.1, 82.0, 6.0, 7.0]
  njlcspecies = n_elements(jlcspecies)

  ; Initialize output structure
  abunds = {teff: 0d, tefferr: 0d, logg: 0d, loggerr: 0d, vt: 0d, vterr: 0d, $
    feh: 0d, feherr: 0d, alphafe: 0d, alphafeerr: 0d, $
    abund: dblarr(njlcspecies), abunderr: dblarr(njlcspecies), $
    abundsyserr: dblarr(njlcspecies), nlines: intarr(njlcspecies), $
    upperlimit: dblarr(njlcspecies)}

  ; -----------------------------------------------------------------
  ; METALLICITY CONSTRAINT MODE
  ; -----------------------------------------------------------------
  if keyword_set(fixedmet) then begin
    ; Use cluster mean metallicity instead of optimizing from Fe lines
    feh_mean = -2.39
    feh_mean_err = 0.05
    alphafe_mean = 0.41
    alphafe_mean_err = 0.1
    abunds.teff = teff
    abunds.tefferr = tefferr
    abunds.logg = logg
    abunds.loggerr = loggerr
    abunds.feh = feh_mean
    abunds.feherr = feh_mean_err
    abunds.alphafe = alphafe_mean
    abunds.alphafeerr = alphafe_mean_err
    abunds.vt = 2.13 - 0.23 * abunds.logg
    abunds.vterr = sqrt(0.05 ^ 2. + (0.03 * abunds.logg) ^ 2. + (0.23 * abunds.loggerr) ^ 2.)
  endif else begin
    feh_mean = 0
    alphafe_mean = 0
  endelse
  if ~keyword_set(fixedmet) then abunds.alphafe = 0.41

  ; -----------------------------------------------------------------
  ; PARAMETER OPTIMIZATION SETUP
  ; -----------------------------------------------------------------
  ; Create parameter file for Fe and alpha element lines
  nlines = trim_linelist(star, /fe, /alpha, teffphot = teffphot, dirflag = dirflag)
  make_par, parfile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.par', $
    atmfile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.atm', $
    driver = 'abfind', $
    linefile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '_temp.ew', $
    outfile = dirflag + star + (teffphot eq 0 ? '' : '_teffphot') + '.out2'

  ; MPFIT parameter structure: Teff, log(g), v_t
  pars = {value: 0d, fixed: 0, name: ' ', limited: [1, 1], limits: [0d, 0d], $
    mpprint: 1, mpformat: ' ', step: 0.0}
  pars = replicate(pars, 3)
  pars[0 : 1].fixed = teffphot ; Fix Teff and logg if using photometric values
  pars.name = ['Teff', 'logg', '  vt']
  pars.step = [25, 0.25, 0.2]
  pars.mpformat = ['(I4)', '(D4.2)', '(D4.2)']
  pars.value = [teff, logg, vt]

  ; Parameter limits
  pars.limits[0] = [3600, 0.0001, 0.0] ; Lower limits
  pars.limits[1] = [8000, 4.9, 8.0] ; Upper limits

  pars.value >= pars.limits[0]
  pars.value <= pars.limits[1]

  ; Chi-squared mask: which diagnostics to use
  ; [dFe/dEP, FeI-FeII, dFe/dRW, dTi/dEP, TiI-TiII, dTi/dRW]
  chimask = teffphot ? [0, 0, 1, 0, 0, 0] : [1, 1, 1, 0, 0, 0]

  ; -----------------------------------------------------------------
  ; RUN OPTIMIZATION
  ; -----------------------------------------------------------------
  functargs = {star: star, chimask: chimask, teffphot: teffphot, verbose: 1, $
    dirflag: dirflag2, fixedfeh: feh_mean}
  p = mpfit('find_abund' + (keyword_set(final) ? '_final' : ''), $
    parinfo = pars, functargs = functargs, /nocatch, $
    ftol = 1d-12, xtol = 1d-12, gtol = 1d-12, $
    status = status, errmsg = errmsg, perror = perror)
  chi = find_abund(p, star = star, chimask = chimask, teffphot = teffphot, $
    verbose = verbose, dirflag = dirflag2, fixedfeh = feh_mean, $
    feh = fehmp, errfeh = errfehmp, alphafe = alphafemp, $
    erralphafe = erralphafemp)

  ; -----------------------------------------------------------------
  ; STORE OPTIMIZED PARAMETERS
  ; -----------------------------------------------------------------
  abunds.teff = p[0]
  abunds.tefferr = keyword_set(teffphot) ? tefferr : perror[0]
  abunds.logg = p[1]
  abunds.loggerr = keyword_set(teffphot) ? loggerr : perror[1]
  abunds.vt = p[2]
  abunds.vterr = perror[2]
  if ~keyword_set(fixedmet) then begin
    abunds.feh = fehmp
    abunds.feherr = errfehmp
  endif
  abunds.alphafe = alphafemp
  abunds.alphafeerr = erralphafemp

  ; -----------------------------------------------------------------
  ; CALCULATE FINAL ABUNDANCES AND ERRORS
  ; -----------------------------------------------------------------
  if keyword_set(final) then begin
    abund = error_analysis_final(star, abunds, teffphot = teffphot)
  endif else begin
    abund = error_analysis(star, abunds, teffphot = teffphot)
  endelse
  if dirflag eq '' then abund = abund[where(abund.ewerr gt 0 or abund.upperlimit eq 1)]

  mwrfits, abund, dirflag2 + star + '_abundbyline' + (teffphot ? '_teffphot' : '') + '.fits', /create
  mwrfits, abunds, dirflag2 + star + '_abund' + (teffphot ? '_teffphot' : '') + '.fits', /create
end

; =================================================================
; PROCEDURE: abund
; =================================================================
; PURPOSE:
; Batch driver to run abundance analysis for all stars in catalog.
; Wrapper around hires_loop that processes multiple stars.
;
; INPUTS:
; ni    - If set, process only star number ni (0-indexed)
; If not set, process all stars
; final - If set, use /final mode (optimized parameters, 'ew2/' directory)
; If not set, use /fixedmet mode (fixed metallicity, 'ew/' directory)
;
; OUTPUT:
; Creates FITS files for each star with abundances and parameters
;
; USAGE:
; abund                    ; Process all stars with fixed metallicity
; abund, ni=5              ; Process only star #5
; abund, /final            ; Process all stars with final mode
;
; =================================================================
pro abund, ni = ni, final = final
  compile_opt idl2

  ; Load catalog
  hiresall = mrdfits('M15_M92_allframes.fits', 1, /silent)
  hiresall = hiresall[where(strtrim(hiresall.name, 2) ne 'M92-star-5' and $
    strtrim(hiresall.name, 2) ne 'M92-star-7')]
  hiresall = hiresall[sort(hiresall.name)]
  n = n_elements(hiresall)

  ; Determine which stars to process
  if ~keyword_set(ni) then begin
    istart = 0
    iend = n - 1
  endif else begin
    istart = fix(ni)
    iend = fix(ni)
  endelse

  ; Loop over stars
  for i = istart, iend do begin
    if keyword_set(final) then begin
      hires_loop, strtrim(hiresall[i].name, 2), /teffphot, /final
    endif else begin
      hires_loop, strtrim(hiresall[i].name, 2), /teffphot, /fixedmet
    endelse
  endfor
end

; =================================================================
; PROCEDURE: abundmc
; =================================================================
; PURPOSE:
; Monte Carlo error propagation for stellar abundances.
; Performs 1000 realizations with perturbed stellar parameters,
; equivalent widths, and photometry to estimate full uncertainties.
;
; INPUTS:
; ni - If set, process only star number ni (0-indexed)
; If not set, process all stars
;
; PERTURBATIONS:
; - Teff from photometry (BP-RP color + parallax)
; - log(g) from mass-radius relation
; - v_t from optimized value + error
; - [Fe/H] from optimized value + error
; - [alpha/Fe] from optimized value + error
; - Equivalent widths from measurement errors
;
; OUTPUT:
; Creates FITS files in mc/ directory with MC realizations
;
; USAGE:
; abundmc              ; Process all stars
; abundmc, ni=10       ; Process only star #10
;
; =================================================================
pro abundmc, ni = ni
  compile_opt idl2
  common ews, ews

  ; Load catalog
  hiresall = mrdfits('M15_M92_allframes.fits', 1, /silent)
  hiresall = hiresall[where(strtrim(hiresall.name, 2) ne 'M92-star-5' and $
    strtrim(hiresall.name, 2) ne 'M92-star-7')]
  hiresall = hiresall[sort(hiresall.name)]
  n = n_elements(hiresall)

  ; Determine which stars to process
  if ~keyword_set(ni) then begin
    istart = 0
    iend = n - 1
  endif else begin
    istart = fix(ni)
    iend = fix(ni)
  endelse

  ; Monte Carlo parameters
  nmc = 1000
  dirflag2 = getenv('CALTECH') + 'hires/M15_M92/ew2/'
  dirflag3 = getenv('CALTECH') + 'hires/M15_M92/mc/'

  ; Physical constants
  sigma_SB = 5.6704d-5 ; Stefan-Boltzmann constant (cgs)
  G = 6.674d-8 ; Gravitational constant (cgs)
  Msun = 1.989d33 ; Solar mass (g)
  M = 0.75 * Msun ; RGB star mass
  Merr = 0.1 * Msun ; Mass error

  ; Teff calibration coefficients (BP-RP color relation for giants)
  bprp_giant = [0.5323, 0.4775, -0.0344, -0.0110, -0.0020, -0.0009]

  ; -----------------------------------------------------------------
  ; MONTE CARLO LOOP
  ; -----------------------------------------------------------------
  for i = istart, iend do begin
    star = strtrim(hiresall[i].name, 2)

    ; Load baseline abundances and EWs
    abunds = mrdfits(dirflag2 + star + '_abund_teffphot.fits', 1, /silent)
    ews_orig = mrdfits(dirflag2 + star + '_Ji20_ew.fits', 1, /silent)

    ; Photometric color for Teff calibration
    color = hiresall[i].bpmag0 - hiresall[i].rpmag0
    colorerr = sqrt(hiresall[i].bpmagerr ^ 2. + hiresall[i].rpmagerr ^ 2.)
    vars = [1d, color, (color) ^ 2., abunds.feh, abunds.feh ^ 2., abunds.feh * (color)]
    varserr = [0d, 1d, 2. * (color), 0d, 0d, abunds.feh]
    varserrfeh = [0d, 0d, 0d, 1d, 2. * abunds.feh, color]
    tefferr = sqrt((abunds.teff * total(bprp_giant * varserr * colorerr) / total(bprp_giant * vars)) ^ 2. + (abunds.teff * total(bprp_giant * varserrfeh * abunds.feherr) / total(bprp_giant * vars)) ^ 2.)
    teffmc = abunds.teff + tefferr * randomn(seed, nmc) + 83. * randomn(42, nmc) ; add 83 K systematic error

    teffdiff = abunds.teff - 5772d
    bcG = poly(teffdiff, [6d-2, 6.731d-5, -6.647d-8, 2.859d-11, -7.197d-15])
    logL = (hiresall[i].gmag0 - (hiresall[i].dm - 2.682 * hiresall[i].ebv) + bcG - 4.68) / (-2.5d) + alog10(3.828d33) ; absolute G_sun = 4.68 (Andrae et al. 2018)

    teffdiffmc = teffmc - 5772d
    bcGmc = poly(teffdiffmc, [6d-2, 6.731d-5, -6.647d-8, 2.859d-11, -7.197d-15])
    logLmc = (hiresall[i].gmag0 - (hiresall[i].dm - 2.682 * hiresall[i].ebv) + bcGmc - 4.68) / (-2.5d) + alog10(3.828d33) ; absolute G_sun = 4.68 (Andrae et al. 2018)

    logg = alog10(4 * !dpi * sigma_SB * G) + alog10(M) + 4. * alog10(abunds.teff) - logL
    loggmc = alog10(4 * !dpi * sigma_SB * G) + alog10(M) + 4. * alog10(teffmc) - logLmc
    loggmc += (Merr / (M * alog(10))) * randomn(seed, nmc) + (hiresall[i].gmagerr * randomn(seed, nmc) / 2.5)

    fehmc = abunds.feh + abunds.feherr * randomn(seed, nmc)

    alphafemc = abunds.alphafe + abunds.alphafeerr * randomn(seed, nmc)

    vt = 2.13 - 0.23 * logg
    vtmc = 2.13 - 0.23 * loggmc + (0.03 * loggmc) * randomn(seed, nmc)

    ews = ews_orig
    interp_atm, abunds.teff, logg, vt, abunds.feh, abunds.alphafe, outfile = dirflag3 + star + '_teffphot.atm'
    abunds_mc = calculate_abund(star, /teffphot, dirflag = 'mc/')
    for j = 0, nmc - 1 do begin
      ews = ews_orig
      w = where(ews.ew gt 0, c)
      ews[w].ew = abs(ews[w].ew + ews[w].ewerr * randomn(seed, c))
      interp_atm, teffmc[j], loggmc[j], vtmc[j], fehmc[j], alphafemc[j], outfile = dirflag3 + star + '_teffphot.atm'
      abund = calculate_abund(star, /teffphot, dirflag = 'mc/')
      abunds_mc = [abunds_mc, abund]
    endfor
    mwrfits, abunds_mc, dirflag3 + star + '_abundmc_teffphot.fits', /create
  endfor
end
