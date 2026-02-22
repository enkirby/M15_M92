; possible blends:
; Si I 3905
; Al I 3944
; "bluest" Sc I lines
; all V I and V II lines

pro abundsynth__define
  compile_opt idl2
  abundsynth = {Abundsynth, element: ' ', ion: 0, species: 0.0, lambda: 0., ep: 0., loggf: 0., abund: 0., upperlimit: 0, abunderr: -999d, upteff: -999d, uplogg: -999d, upvt: -999d, upfeh: -999d, upalphafe: -999d, weight: 0d, delta: dblarr(4), vshift: 0d, vshifterr: 0d}
end

function run_synth, x, pars, minlambda = minlambda, maxlambda = maxlambda, star = star, element = element, wave = wave, atmpars = atmpars, r = R, abund_prev = abund_prev, upteff = upteff, uplogg = uplogg, upvt = upvt, upfeh = upfeh
  compile_opt idl2
  clight = 2.99792458d5

  case 1 of
    atmpars.logg le 2.0: c12c13 = 6.0
    atmpars.logg gt 2.7: c12c13 = 50.0
    else: c12c13 = 63.0 * atmpars.logg - 120.0
  endcase

  filename = element + (keyword_set(wave) ? string(wave, format = '(I4)') : '')

  e = elements(/newmoog)
  w = where(e.name eq element)

  if e[w].atomic gt 31 then begin
    wiso = where(e[w].isotope gt 0, niso)
    isotopes = e[w].isotope[wiso]
    isofracs = e[w].atomic ge 31 ? e[w].rfrac[wiso] : e[w].solarfrac[wiso]
    make_par, parfile = 'synth/' + star + '_' + filename + '.par', linefile = 'synth/' + filename + '.list.linemake', atmfile = 'synth/' + star + '_' + filename + '.atm', outfile = 'synth/' + star + '_' + filename + '.out2', driver = 'synth', minlambda = minlambda, maxlambda = maxlambda, c12c13 = c12c13, atomic = e[w].atomic + (e[w].atomic gt 31 ? 0.1 : 0.0), isotopes = isotopes, isofracs = isofracs
  endif else if e[w].atomic eq 3.0 then begin
    isotopes = [6, 7]
    isofracs = [0.0000001, 1.0 - 0.0000001]
    make_par, parfile = 'synth/' + star + '_' + filename + '.par', linefile = 'synth/' + filename + '.list.linemake', atmfile = 'synth/' + star + '_' + filename + '.atm', outfile = 'synth/' + star + '_' + filename + '.out2', driver = 'synth', minlambda = minlambda, maxlambda = maxlambda, c12c13 = c12c13, atomic = 3.0, isotopes = isotopes, isofracs = isofracs
  endif else begin
    make_par, parfile = 'synth/' + star + '_' + filename + '.par', linefile = 'synth/' + filename + '.list.linemake', atmfile = 'synth/' + star + '_' + filename + '.atm', outfile = 'synth/' + star + '_' + filename + '.out2', driver = 'synth', minlambda = minlambda, maxlambda = maxlambda, c12c13 = c12c13
  endelse

  ; jlcspecies = [3.0, 8.0, 11.0, 12.0, 13.0, 14.0, 19.0, 20.0, 20.1, 21.1, 22.0, 22.1, 23.0, 23.1, 24.0, 24.1, 25.0, 25.1, 26.0, 26.1, 27.0, 28.0, 29.0, 30.0, 38.0, 38.1, 39.1, 40.0, 40.1, 56.1, 57.1, 58.1, 59.1, 60.1, 62.1, 63.1, 64.1, 66.1, 67.1, 82.0, 6.0, 7.0]
  ; keep = [0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
  ; wk = where(keep)
  ; wtweak = wk[where(abund_prev.abund[wk] gt -10 and abund_prev.abund[wk] lt 10 and abund_prev.abunderr[wk] gt 0.0 and finite(abund_prev.abunderr[wk]) and floor(jlcspecies[wk]) ne e[w].atomic)]
  ; tweakel = [floor(jlcspecies[wtweak]), e[w].atomic]
  ; match, e.atomic, floor(jlcspecies[wtweak]), ws, wt
  ; ws = ws[sort(wt)]
  ; tweakabund = [abund_prev.abund[wtweak] - atmpars.feh - e[ws].solar, (e[w].atomic le 4 ? pars[0] : (pars[0] - atmpars.feh - e[w].solar))]
  ; s = sort(tweakel)

  flag = ''
  if keyword_set(upteff) then flag += '_upteff'
  if keyword_set(uplogg) then flag += '_uplogg'
  if keyword_set(upvt) then flag += '_upvt'
  if keyword_set(upfeh) then flag += '_upfeh'
  atlas_to_moog, '/raid/atlas/BasicATLAS/ATLAS_LMHA_' + star + flag + '/output_summary.out', 'synth/' + star + '_' + filename + '.atm', vt = atmpars.vt, tweakels = e[w].atomic, tweakabunds = pars[0]

  spawn, 'MOOGSILENT synth/' + star + '_' + filename + '.par'
  moog = read_moog_spec('synth/' + star + '_' + filename + '.out2', /newmoog)
  moogspec = smooth_gauss_wrapper(moog.lambda * (1d + pars[1] / clight), moog.spec, x, x / R / 2.35482)
  return, moogspec
end

function synth_ul, pars, df, _extra = _extra
  compile_opt idl2
  common synth_ul, x, ymp, dymp
  common chisq, chi2

  moogspecul = run_synth(x, pars, _extra = _extra)
  chisq = total(((ymp - moogspecul) / dymp) ^ 2)
  return, abs(chisq - chi2)
end

function calculate_xfe, star, element, abund_prev, wave = wave, perror = perror, r = R, ul = ul, noul = noul, upteff = upteff, uplogg = uplogg, upvt = upvt, upfeh = upfeh, startabund = startabund, vr = vr, vshift = vshift
  compile_opt idl2
  common synth_ul, x, ymp, dymp
  common chisq, chi2

  filename = 'spectra/' + star + '.fits.gz'
  hires = mrdfits(filename, 1, /silent)

  filename = 'ew/' + star + '_abund_teffphot.fits'
  atmpars = mrdfits(filename, 1, /silent)
  if keyword_set(upteff) then atmpars.teff += atmpars.tefferr
  if keyword_set(uplogg) then atmpars.logg += atmpars.loggerr
  if keyword_set(upvt) then atmpars.vt += atmpars.vterr
  if keyword_set(upfeh) then atmpars.feh += atmpars.feherr

  e = elements(/newmoog)
  we = where(e.name eq element)

  if ~keyword_set(startabund) then begin
    maxiter = (element eq 'C' or element eq 'N') ? 50 : 2
    startabund = [element eq 'Li' ? 2.2d : e[we].solar - 2.39, 0.0]
  endif else maxiter = (element eq 'C' or element eq 'N') ? 5 : 1
  pi = {value: 0.0, fixed: 0, limited: [1, 1], limits: [-5.0d, 10.0d], parname: 'A(' + element + ')', mpprint: 1, mpformat: '(D6.3)', step: 0.05d}
  pi = replicate(pi, 2)
  pi.value = startabund
  pi.fixed = [0, keyword_set(vshift)]
  pi.limited = [[1, 1], [1, 1]]
  pi.limits = [[-5.0d, 10.0d], [-3.0, 3.0]]
  pi.parname = ['A(' + element + ')', 'v (km/s)']
  pi.mpprint = [1, 1]
  pi.mpformat = ['(D6.3)', '(D6.3)']
  pi.step = [0.05d, 0.05d]

  linelist = 'synth/' + element + (keyword_set(wave) ? string(wave, format = '(I4)') : '') + '.list.linemake'
  openr, lun, linelist, /get_lun
  skip_lun, lun, 1, /lines
  minlambda = 10000.0
  maxlambda = 0.0
  while ~eof(lun) do begin
    readf, lun, inlambda, format = '(D10)'
    minlambda <= inlambda
    maxlambda >= inlambda
  endwhile
  close, lun
  free_lun, lun

  clight = 2.99792458d5

  lambda = hires.lambda / (1d + vr / clight)
  airtovac, lambda
  w = where(lambda ge minlambda and lambda le maxlambda and finite(hires.ivar) and hires.ivar gt 0.0, c)
  if c lt 10 then begin
    perror = [0.0d, 0.0d]
    ul = 0.0d
    return, [0.0d, 0.0d]
  endif

  x = lambda[w]
  y = hires.spec[w]
  dy = (hires.ivar[w]) ^ (-0.5)
  ymp = y
  dymp = dy

  if ~keyword_set(wave) then wave = 0b

  iter = 0l
  loop = 1b
  oldxfe = 999d
  perrorv = 0.0
  while loop and iter le maxiter do begin
    pars = mpfitfun('run_synth', x, ymp, dymp, parinfo = pi, /nocatch, bestnorm = chisq0, dof = dof, perror = perror, ftol = 1d-10, gtol = 1d-10, xtol = 1d-10, covar = covar, status = status, yfit = ymoog, nprint = 1000, functargs = {minlambda: minlambda, maxlambda: maxlambda, star: star, element: element, wave: wave, atmpars: atmpars, r: R, abund_prev: abund_prev, upteff: keyword_set(upteff), uplogg: keyword_set(uplogg), upvt: keyword_set(upvt), upfeh: keyword_set(upfeh)})
    pi.value = pars
    if pi[1].fixed eq 0 then perrorv = perror[1]
    pi[1].fixed = 1
    perror[1] = perrorv
    ; plot, x, ymp
    ; oplot, x, ymoog, color=fsc_color('red')
    ; stop

    ; ymoog = run_synth(x, pars)
    resid = y / ymoog
    residerr = dy / ymoog
    if 1 or (element eq 'C' or element eq 'N') then begin
      bkpt = slatec_splinefit(x, resid, coeff, invvar = residerr ^ (-2d), bkspace = (element eq 'C' or element eq 'N') ? 10 : 5, upper = 2, lower = 2, /silent)
      cont = slatec_bvalu(x, bkpt, coeff)
    endif else cont = median(resid) ; weightedmean(resid, residerr)

    ympplot = ymp
    ymp = y / cont
    dymp = dy / cont

    xfediff = oldxfe - pars[0]
    loop = abs(xfediff) gt 0.001
    oldxfe = pars[0]
    iter++
  endwhile
  ; plot, x, ympplot, yrange = [0.2, 1.5]
  ; oplot, x, ymoog, color = fsc_color('red')
  ; ; oplot, x, resid, color = fsc_color('blue')
  ; oplot, x, replicate(1.0, c), color = fsc_color('orange')

  if ~keyword_set(noul) then begin
    chi2 = 0.0
    chisq0 = synth_ul(pars, 0.0, minlambda = minlambda, maxlambda = maxlambda, star = star, element = element, wave = wave, atmpars = atmpars, r = R, abund_prev = abund_prev)
    chi2 = chisq0 + 9.0
    pi[0].limits = [pars[0], 5.0]
    pi[0].step = 0.05
    pi[0].value = (pars[0] + pi[0].step)
    result = tnmin('synth_ul', parinfo = pi, quiet = 0, bestmin = bestmin, /autoderivative, functargs = {minlambda: minlambda, maxlambda: maxlambda, star: star, element: element, wave: wave, atmpars: atmpars, r: R, abund_prev: abund_prev})
    ul = result[0]
    ; chisqul3 = bestmin + chi2
  endif

  return, pars
end

pro error_analysis_synth, abundsynth, abund_prev, name = name, r = R, vr = vr
  compile_opt idl2
  nsynths = n_elements(abundsynth)

  for k = 0, nsynths - 1 do begin
    abund = calculate_xfe(name, abundsynth[k].element, abund_prev, perror = abunderr, r = R, ul = abundul, wave = round(abundsynth[k].lambda), vr = vr)
    abundsynth[k].vshift = abund[1]
    abundsynth[k].vshifterr = abunderr[1]
    if abunderr[0] gt 0.15 then begin
      abundsynth[k].abund = abundul
      abundsynth[k].upperlimit = 1
    endif else begin
      abundsynth[k].abund = abund[0]
      abundsynth[k].abunderr = abunderr[0]
      abundsynth[k].upperlimit = 0

      abund_upteff = calculate_xfe(name, abundsynth[k].element, abund_prev, r = R, /noul, wave = round(abundsynth[k].lambda), /upteff, startabund = abund, vr = vr, vshift = abundsynth[k].vshift)
      abund_uplogg = calculate_xfe(name, abundsynth[k].element, abund_prev, r = R, /noul, wave = round(abundsynth[k].lambda), /uplogg, startabund = abund, vr = vr, vshift = abundsynth[k].vshift)
      abund_upvt = calculate_xfe(name, abundsynth[k].element, abund_prev, r = R, /noul, wave = round(abundsynth[k].lambda), /upvt, startabund = abund, vr = vr, vshift = abundsynth[k].vshift)
      abund_upfeh = calculate_xfe(name, abundsynth[k].element, abund_prev, r = R, /noul, wave = round(abundsynth[k].lambda), /upfeh, startabund = abund, vr = vr, vshift = abundsynth[k].vshift)
      abundsynth[k].upteff = abund_upteff[0] - abund[0]
      abundsynth[k].uplogg = abund_uplogg[0] - abund[0]
      abundsynth[k].upvt = abund_upvt[0] - abund[0]
      abundsynth[k].upfeh = abund_upfeh[0] - abund[0]
    endelse
  endfor
end

pro synth, ni = ni
  compile_opt idl2
  e = elements(/newmoog)

  hiresall = mrdfits('M15_M92_allframes.fits', 1, /silent)
  hiresall = hiresall[where(strtrim(hiresall.name, 2) ne 'M92-star-5' and strtrim(hiresall.name, 2) ne 'M92-star-7')]
  n = n_elements(hiresall)
  if ~keyword_set(ni) then begin
    message, 'No star index provided.'
  endif else begin
    istart = fix(ni)
    iend = fix(ni)
    ; filestring = string(fix(ni), format = '(I02)')
    filestring = strtrim(hiresall[ni].name, 2)
  endelse
  synthfilename = 'synth/abundsynth_' + filestring + '.fits'

  cat = mrdfits('M15_M92_catalog.fits', 1, /silent)
  element_names = mrdfits('M15_M92_elements.fits', 1, /silent)
  match, strtrim(e.name), strtrim(element_names.element), w1, w2
  catspecies = e[w1].atomic

  abund1 = cat
  jlcspecies = [3.0, 8.0, 11.0, 12.0, 13.0, 14.0, 19.0, 20.0, 20.1, 21.1, 22.0, 22.1, 23.0, 23.1, 24.0, 24.1, 25.0, 25.1, 26.0, 26.1, 27.0, 28.0, 29.0, 30.0, 38.0, 38.1, 39.1, 40.0, 40.1, 56.1, 57.1, 58.1, 59.1, 60.1, 62.1, 63.1, 64.1, 66.1, 67.1, 82.0, 6.0, 7.0]
  njlcspecies = n_elements(jlcspecies)
  abund1 = struct_trimtags(abund1, except_tags = ['abund'])
  abund1 = struct_addtags(abund1, replicate({abund: dblarr(njlcspecies)}, n_elements(abund1)))
  for i = 0, n_elements(jlcspecies) - 1 do begin
    w = where(floor(catspecies) eq floor(jlcspecies[i]))
    if n_elements(w) eq 1 then abund1.abund[i] = abund1.abunddiffplusavgabund[w[0]] else abund1.abund[i] = -999d
  endfor

  abund_prev = abund1[where(strtrim(abund1.name, 2) eq strtrim(hiresall[ni].name, 2))]

  hiresall = struct_trimtags(hiresall, except_tags = ['R'])
  newstr = replicate({r: 0d}, n)
  hiresall = struct_addtags(hiresall, newstr)
  synthels = ['Li', 'C', 'N', 'Al', 'Sr', 'Y', 'Zr', 'Ba', 'La', 'Eu', 'Dy']
  nsynthels = n_elements(synthels)
  for i = 0, nsynthels - 1 do begin
    hiresall = struct_trimtags(hiresall, except_tags = ['A' + strupcase(synthels[i]), 'A' + strupcase(synthels[i]) + 'ERR', 'A' + strupcase(synthels[i]) + 'UL', 'A' + strupcase(synthels[i]) + 'WEIGHT', 'A' + strupcase(synthels[i]) + '_UPTEFF', 'A' + strupcase(synthels[i]) + '_UPLOGG', 'A' + strupcase(synthels[i]) + '_UPVT', 'A' + strupcase(synthels[i]) + '_UPFEH', 'A' + strupcase(synthels[i]) + 'DELTA'])
  endfor

  readcol, '../M92_KOA/Ji20_linelist.moog', lambda, species, ep, loggf, format = 'D,F,D,D', skipline = 1, /silent
  synthspecies = [13.0, 38.1, 39.1, 40.1, 56.1, 57.1, 63.1, 66.1]
  wlines = [-1]
  nlines = 0
  for j = 0, n_elements(synthspecies) - 1 do begin
    wlinesj = where(round(species * 10.) eq round(synthspecies[j] * 10.), nlinesj)
    if nlinesj eq 0 then message, 'Species not found.'
    wlines = [wlines, wlinesj]
    nlines += nlinesj
  endfor
  wlines = wlines[1 : nlines]
  extrasynths = 3
  abundsynth = replicate({Abundsynth}, nlines + extrasynths)
  if extrasynths gt 0 then begin
    abundsynth[0 : extrasynths - 1].element = ['Li', 'C', 'N']
    abundsynth[0 : extrasynths - 1].ion = [1, 1, 1]
    abundsynth[0 : extrasynths - 1].species = [3.0, 6.0, 7.0]
    abundsynth[0 : extrasynths - 1].lambda = [0.0, 0.0, 0.0]
    abundsynth[0 : extrasynths - 1].ep = [0.0, 0.0, 0.0]
    abundsynth[0 : extrasynths - 1].loggf = [0.0, 0.0, 0.0]
  endif
  for i = 0, nlines - 1 do begin
    we = where(e.atomic eq floor(species[wlines[i]]))
    abundsynth[extrasynths + i].element = e[we].name
    abundsynth[extrasynths + i].ion = round(10. * (species[wlines[i]] - floor(species[wlines[i]]))) + 1
  endfor
  abundsynth[extrasynths : extrasynths + nlines - 1].species = species[wlines]
  abundsynth[extrasynths : extrasynths + nlines - 1].lambda = lambda[wlines]
  abundsynth[extrasynths : extrasynths + nlines - 1].ep = ep[wlines]
  abundsynth[extrasynths : extrasynths + nlines - 1].loggf = loggf[wlines]
  if file_test(synthfilename) then begin
    abundsynth_old = mrdfits(synthfilename, 1, /silent)
    for j = 0, n_elements(abundsynth) - 1 do begin
      w = where(round(abundsynth_old.species * 10.) eq round(abundsynth[j].species * 10.) and round(abundsynth_old.lambda) eq round(abundsynth[j].lambda), c)
      if c eq 1 then begin
        w = w[0]
        abundsynth[j].vshift = abundsynth_old[w].vshift
        abundsynth[j].vshifterr = abundsynth_old[w].vshifterr
        abundsynth[j].abund = abundsynth_old[w].abund
        abundsynth[j].abunderr = abundsynth_old[w].abunderr
        abundsynth[j].upperlimit = abundsynth_old[w].upperlimit
        abundsynth[j].upteff = abundsynth_old[w].upteff
        abundsynth[j].uplogg = abundsynth_old[w].uplogg
        abundsynth[j].upvt = abundsynth_old[w].upvt
        abundsynth[j].upfeh = abundsynth_old[w].upfeh
        abundsynth[j].weight = abundsynth_old[w].weight
      endif
    endfor
  endif
  nsynths = n_elements(abundsynth)

  for i = istart, iend do begin
    name = strtrim(hiresall[i].name, 2)

    ewfile = 'ew/' + name + '_Ji20_ew.fits'
    ew = mrdfits(ewfile, 1, /silent)
    w = where(ew.doppwidth gt 0 and ew.rate lt 0.1)
    hiresall[i].r = 1.1 * median(ew[w].lambda / (ew[w].doppwidth * 1.66511))

    error_analysis_synth, abundsynth, abund_prev, name = name, r = hiresall[i].r, vr = hiresall[i].vr

    mwrfits, abundsynth, synthfilename, /create
  endfor
end

pro plot_synth
  compile_opt idl2
  ; star_indices = [0, 1, 2, 3, 4, 5]
  star_indices = [2, 13, 20]
  ; line_indices = [0]                   ;Li
  ; line_indices = [1]                   ;C
  ; line_indices = [2]                   ;N
  ; line_indices = [3, 4]                ;Al
  ; line_indices = [5, 6]                ;Sr
  ; line_indices = [7, 8, 9, 10, 11, 12, 13, 14] ;Y
  ; line_indices = [15]                  ;Zr
  ; line_indices = [16, 17, 18, 19, 20]  ;Ba
  ; line_indices = [21, 22, 23, 24, 25]  ;La
  ; line_indices = [26, 27, 28, 29, 30]  ;Eu
  line_indices = [26, 27] ; Eu
  ; line_indices = [31, 32, 33]          ;Dy
  star_indices = reverse(star_indices)

  restore, 'ew/abunds.sav'

  setplot
  device, filename = 'M92_synths.eps', xsize = 8, ysize = 8, /inches, /color, /encapsulated
  device, /isolatin

  e = elements(/newmoog)
  nstars = n_elements(star_indices)
  nlines = n_elements(line_indices)

  x0 = 0.10
  x1 = 0.96
  y0 = 0.09
  y1 = 0.99
  xbuf = 0.01
  ybuf = 0.01
  dxp = (x1 - x0 - (nlines - 1) * xbuf) / nlines
  dyp = (y1 - y0 - (nstars - 1) * ybuf) / nstars
  xp = dindgen(nlines) * (dxp + xbuf) + x0
  yp = dindgen(nstars) * (dyp + ybuf) + y0

  xyouts, (x0 + x1) / 2., 0.01, '!7rest wavelength (' + string(197b) + ')', /normal, align = 0.5
  xyouts, 0.03, (y0 + ybuf + y1) / 2., '!7normalized flux', /normal, align = 0.5, orientation = 90

  hiresall = mrdfits('M92_allframes.fits', 1, /silent)
  hiresall = hiresall[sort(hiresall.gmag0 and strtrim(hiresall.name, 2) ne 'X-20' and strtrim(hiresall.name, 2) ne 'S2303')]
  n = n_elements(hiresall)

  for i = 0, nstars - 1 do begin
    ii = star_indices[i]
    filestring = string(fix(ii), format = '(I02)')
    star = strtrim(hiresall[ii].name, 2)
    abund_prev = abund1[i]

    ewfile = 'ew/' + star + '_ew.fits'
    ew = mrdfits(ewfile, 1, /silent)
    w = where(ew.doppwidth gt 0 and ew.rate lt 0.1)
    R = median(ew[w].lambda / (ew[w].doppwidth * 1.66511))

    abundsynth = mrdfits('M92_abundsynth_' + filestring + '.fits', 1, /silent)
    filename = 'spectra/HIRES_' + star + '.fits.gz'
    hires = mrdfits(filename, 1, /silent)

    filename = 'ew/' + star + '_abund_teffphot.fits'
    atmpars = mrdfits(filename, 1, /silent)

    for k = 0, nlines - 1 do begin
      kk = line_indices[k]
      if abundsynth[kk].abunderr le 0.0 then continue

      elk = strtrim(abundsynth[kk].element, 2)
      we = where(strtrim(e.name, 2) eq elk)
      wave = round(abundsynth[kk].lambda)
      linelist = 'synth/' + elk + (keyword_set(wave) ? string(wave, format = '(I4)') : '') + '.list.linemake'
      openr, lun, linelist, /get_lun
      skip_lun, lun, 1, /lines
      minlambda = 10000.0
      maxlambda = 0.0
      while ~eof(lun) do begin
        readf, lun, inlambda, format = '(D10)'
        minlambda <= inlambda
        maxlambda >= inlambda
      endwhile
      close, lun
      free_lun, lun

      clight = 2.99792458d5
      lambda = hires.lambda / (1d + hires.vr / clight)
      w = where(lambda ge minlambda and lambda le maxlambda and finite(hires.ivar) and hires.ivar gt 0.0, c)
      if c lt 10 then continue

      x = lambda[w]
      y = hires.spec[w]
      dy = (hires.ivar[w]) ^ (-0.5)

      ymoog = run_synth(x, [abundsynth[kk].abund, abundsynth[kk].vshift], minlambda = minlambda, maxlambda = maxlambda, star = star, element = elk, wave = wave, atmpars = atmpars, r = R, abund_prev = abund_prev)
      resid = y / ymoog
      residerr = dy / ymoog
      if elk eq 'C' or elk eq 'N' then begin
        bkpt = slatec_splinefit(x, resid, coeff, invvar = residerr ^ (-2d), bkspace = 10, upper = 2, lower = 2, /silent)
        cont = slatec_bvalu(x, bkpt, coeff)
      endif else cont = median(resid) ; weightedmean(resid, residerr)
      y = y / cont
      dy = dy / cont

      case elk of
        'C': xrange = [4279, 4316]
        'N': xrange = [minlambda, maxlambda]
        else: xrange = abundsynth[kk].lambda + [-1, 1]
      endcase

      plot, [0, 1], [0, 1], /nodata, pos = [xp[k], yp[i], xp[k] + dxp, yp[i] + dyp], /normal, xtickname = replicate(i eq 0 ? '' : ' ', 30), ytickname = replicate(k eq 0 ? '' : ' ', 30), xrange = xrange, yrange = [0.0, 1.2], /xstyle, /ystyle, /noerase
      oplot, [abundsynth[kk].lambda, abundsynth[kk].lambda], [-100, 100], linestyle = 1, color = fsc_color('gray')
      oplot, x, y
      oplot, x, ymoog, color = fsc_color('red')
      if k eq 0 then begin
        xyouts, 0.97 * !x.crange[0] + 0.03 * !x.crange[1], 0.95 * !y.crange[0] + 0.05 * !y.crange[1], '!7' + star, charsize = 1.0
        xyouts, 0.65 * !x.crange[0] + 0.35 * !x.crange[1], 0.95 * !y.crange[0] + 0.05 * !y.crange[1], '!19T!7!Deff!N = ' + string(atmpars.teff, format = '(I4)') + ' K', charsize = 1.0
        xyouts, 0.33 * !x.crange[0] + 0.67 * !x.crange[1], 0.95 * !y.crange[0] + 0.05 * !y.crange[1], '!19A!7(' + elk + ') = ' + strtrim(string(abundsynth[kk].abund, format = '(D5.2)'), 2), charsize = 1.0
      endif
    endfor
  endfor

  device, /close
  resetplot
  stop
end
