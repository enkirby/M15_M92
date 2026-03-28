function run_synth, x, pars, minlambda = minlambda, maxlambda = maxlambda, star = star, element = element, wave = wave, atmpars = atmpars, r = R
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

  atlas_to_moog, '/raid/atlas/BasicATLAS/ATLAS_LMHA_' + star + '/output_summary.out', 'synth/' + star + '_' + filename + '.atm', vt = atmpars.vt, tweakels = e[w].atomic, tweakabunds = pars[0]

  spawn, 'MOOGSILENT synth/' + star + '_' + filename + '.par'
  moog = read_moog_spec('synth/' + star + '_' + filename + '.out2', /newmoog)
  moogspec = smooth_gauss_wrapper(moog.lambda * (1d + pars[1] / clight), moog.spec, x, x / R / 2.35482)
  return, moogspec
end

function synth_wrapper, abundsynth, hires, cat, name, r = R
  compile_opt idl2
  minlambda = abundsynth.lambda - 5
  maxlambda = abundsynth.lambda + 5

  element = strtrim(abundsynth.element, 2)
  wave = round(abundsynth.lambda)
  pars = [abundsynth.abund, abundsynth.vshift]

  printwave = ((element ne 'C') and (element ne 'N') and (element ne 'Li'))

  linelist = 'synth/' + element + (printwave ? string(wave, format = '(I4)') : '') + '.list.linemake'
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

  w = where(hires.lambda ge minlambda and hires.lambda le maxlambda and hires.ivar gt 0)
  std = stddev(hires.ivar[w])
  w = w[where(hires.ivar[w] lt (mean(hires.ivar[w]) + 5*std))]
  x = hires.lambda[w]
  y = hires.spec[w]
  dy = (hires.ivar[w]) ^ (-0.5)

  ymoog = run_synth(x, pars, minlambda = minlambda, maxlambda = maxlambda, star = name, element = element, wave = wave, atmpars = cat, r = R)

  resid = y / ymoog
  residerr = dy / ymoog
  bkpt = slatec_splinefit(x, resid, coeff, invvar = residerr ^ (-2d), bkspace = (element eq 'C' or element eq 'N') ? 10 : 5, upper = 2, lower = 2, /silent)
  cont = slatec_bvalu(x, bkpt, coeff)

  y /= cont
  dy /= cont

  return, {lambda: x, spec: y, specerr: dy, moogspec: ymoog}
end

pro synth_gui_load_star, info, star_index
  compile_opt idl2
  clight = 2.99792458d5

  hiresall = *info.hiresall_ptr
  allcat = *info.allcat_ptr

  name = strtrim(hiresall[star_index].name, 2)

  synthfilename = 'synth/abundsynth_' + name + '.fits'
  cat = allcat[where(strtrim(allcat.name, 2) eq name)]
  ewfile = 'ew/' + name + '_Ji20_ew.fits'
  ew = mrdfits(ewfile, 1, /silent)
  w = where(ew.doppwidth gt 0 and ew.rate lt 0.1, cnt)
  r = cnt gt 0 ? 1.1 * median(ew[w].lambda / (ew[w].doppwidth * 1.66511)) : 50000.0

  abundsynth = mrdfits(synthfilename, 1, /silent)
  n_lines = n_elements(abundsynth)

  filename = 'spectra/' + name + '.fits.gz'
  hires = mrdfits(filename, 1, /silent)
  lambda = hires.lambda / (1d + cat.vr / clight)
  airtovac, lambda
  hires.lambda = lambda

  ; Preserve current line index, but clamp to available lines in new star
  preserved_index = info.current_index < (n_lines - 1)

  ; Update pointers with new data
  ptr_free, info.abundsynth_ptr
  ptr_free, info.hires_ptr
  ptr_free, info.cat_ptr
  info.abundsynth_ptr = ptr_new(abundsynth)
  info.hires_ptr = ptr_new(hires)
  info.cat_ptr = ptr_new(cat)
  info.name = name
  info.r = r
  info.n_lines = n_lines
  info.current_index = preserved_index

  ; Update line droplist labels
  droplist_labels = strarr(n_lines)
  for i = 0, n_lines - 1 do begin
    droplist_labels[i] = strtrim(abundsynth[i].element, 2) + ' ' + $
      string(abundsynth[i].lambda, format = '(F8.2)') + ' A'
  endfor
  widget_control, info.droplist, set_value = droplist_labels
  widget_control, info.droplist, set_droplist_select = preserved_index

  ; Update window title
  widget_control, info.tlb, tlb_set_title = name + ' Spectral Line Viewer'
end

pro synth_gui_event, event
  compile_opt idl2

  widget_control, event.top, get_uvalue = info

  case event.id of
    info.star_back_button: begin
      ; Move to previous star
      hiresall = *info.hiresall_ptr
      n_stars = n_elements(hiresall)
      info.current_star_index = (info.current_star_index - 1 + n_stars) mod n_stars
      widget_control, info.star_droplist, set_droplist_select = info.current_star_index
      widget_control, event.top, set_uvalue = info
      synth_gui_load_star, info, info.current_star_index
      widget_control, event.top, set_uvalue = info
      synth_gui_plot, info
    end

    info.star_next_button: begin
      ; Move to next star
      hiresall = *info.hiresall_ptr
      n_stars = n_elements(hiresall)
      info.current_star_index = (info.current_star_index + 1) mod n_stars
      widget_control, info.star_droplist, set_droplist_select = info.current_star_index
      widget_control, event.top, set_uvalue = info
      synth_gui_load_star, info, info.current_star_index
      widget_control, event.top, set_uvalue = info
      synth_gui_plot, info
    end

    info.star_droplist: begin
      ; Update the star
      info.current_star_index = event.index
      widget_control, event.top, set_uvalue = info
      synth_gui_load_star, info, event.index
      widget_control, event.top, set_uvalue = info
      synth_gui_plot, info
    end

    info.droplist: begin
      ; Update the current line index
      info.current_index = event.index
      widget_control, event.top, set_uvalue = info
      synth_gui_plot, info
    end

    info.back_button: begin
      ; Move to previous line
      info.current_index = (info.current_index - 1 + info.n_lines) mod info.n_lines
      widget_control, info.droplist, set_droplist_select = info.current_index
      widget_control, event.top, set_uvalue = info
      synth_gui_plot, info
    end

    info.next_button: begin
      ; Move to next line
      info.current_index = (info.current_index + 1) mod info.n_lines
      widget_control, info.droplist, set_droplist_select = info.current_index
      widget_control, event.top, set_uvalue = info
      synth_gui_plot, info
    end

    info.quit_button: begin
      widget_control, event.top, /destroy
      return
    end

    else:
  endcase
end

pro synth_gui_plot, info
  compile_opt idl2

  ; Set the draw widget as the current plotting window
  wset, info.draw_id

  ; Get the currently selected line
  idx = info.current_index
  abundsynth = (*info.abundsynth_ptr)[idx]

  ; Compute the synthetic spectrum
  spec = synth_wrapper(abundsynth, *info.hires_ptr, *info.cat_ptr, info.name, r = info.r)

  ; Get y-range for the plot
  ymin = min([spec.spec, spec.moogspec], max = ymax)
  yrange = [ymin - 0.05 * (ymax - ymin), ymax + 0.05 * (ymax - ymin)]

  ; Plot the observed spectrum with white background and black axes
  plot, spec.lambda, spec.spec, $
    xtitle = 'Wavelength (Angstroms)', $
    ytitle = 'Normalized Flux', $
    title = info.name + ': ' + strtrim(abundsynth.element, 2) + ' ' + string(abundsynth.lambda, format = '(F8.2)') + ' A', $
    charsize = 1.2, $
    xstyle = 1, $
    ystyle = 1, $
    yrange = yrange, $
    background = fsc_color('white'), $
    color = fsc_color('black')

  ; Add shaded region around target line (+/- 0.15 Angstroms)
  target_lambda = abundsynth.lambda
  shade_lambda_min = target_lambda - 0.15
  shade_lambda_max = target_lambda + 0.15

  ; Create polygon vertices for shading
  shade_x = [shade_lambda_min, shade_lambda_max, shade_lambda_max, shade_lambda_min, shade_lambda_min]
  shade_y = [yrange[0], yrange[0], yrange[1], yrange[1], yrange[0]]

  ; Fill the shaded region with light gray
  polyfill, shade_x, shade_y, color = fsc_color('light gray'), /data

  ; Replot the data on top of the shading
  oplot, spec.lambda, spec.spec, color = fsc_color('black')

  ; Overplot the synthetic spectrum in red
  oplot, spec.lambda, spec.moogspec, color = fsc_color('red'), thick = 2

  ; Update the abundance labels
  abund_text = 'ABUND = ' + string(abundsynth.abund, format = '(F6.2)')
  abunderr_text = 'ABUNDERR = ' + string(abundsynth.abunderr, format = '(F5.2)')
  widget_control, info.abund_label, set_value = abund_text
  widget_control, info.abunderr_label, set_value = abunderr_text

  ; Add a legend
  ; legend, ['Observed', 'Synthetic'], linestyle = [0, 0], color = [fsc_color('black'), fsc_color('red')], $
  ; thick = [1, 2], /top, /right, box = 0, charsize = 1.0
end

pro synth_gui_cleanup, tlb
  compile_opt idl2

  widget_control, tlb, get_uvalue = info

  ; Free the pointers
  ptr_free, info.abundsynth_ptr
  ptr_free, info.hires_ptr
  ptr_free, info.cat_ptr
  ptr_free, info.hiresall_ptr
  ptr_free, info.allcat_ptr
end

pro synth_gui, ni = ni, element = element
  compile_opt idl2
  clight = 2.99792458d5

  ; Load all stars
  hiresall = mrdfits('M15_M92_allframes.fits', 1, /silent)
  allcat = mrdfits('M15_M92_catalog.fits', 1, /silent)

  hiresall = hiresall[where(strtrim(hiresall.name, 2) ne 'M92-star-5' and strtrim(hiresall.name, 2) ne 'M92-star-7')]
  n_stars = n_elements(hiresall)

  ; Set default star index
  if ~keyword_set(ni) then ni = 0
  current_star_index = ni

  ; Create star names for dropdown
  star_names = strtrim(hiresall.name, 2)

  ; Load initial star data
  name = star_names[current_star_index]
  synthfilename = 'synth/abundsynth_' + name + '.fits'
  cat = allcat[where(strtrim(allcat.name, 2) eq name)]
  ewfile = 'ew/' + name + '_Ji20_ew.fits'
  ew = mrdfits(ewfile, 1, /silent)
  w = where(ew.doppwidth gt 0 and ew.rate lt 0.1, cnt)
  r = cnt gt 0 ? 1.1 * median(ew[w].lambda / (ew[w].doppwidth * 1.66511)) : 50000.0

  abundsynth = mrdfits(synthfilename, 1, /silent)
  n_lines = n_elements(abundsynth)

  filename = 'spectra/' + name + '.fits.gz'
  hires = mrdfits(filename, 1, /silent)
  lambda = hires.lambda / (1d + cat.vr / clight)
  airtovac, lambda
  hires.lambda = lambda

  ; Create droplist labels for lines
  droplist_labels = strarr(n_lines)
  for i = 0, n_lines - 1 do begin
    droplist_labels[i] = strtrim(abundsynth[i].element, 2) + ' ' + $
      string(abundsynth[i].lambda, format = '(F8.2)') + ' A'
  endfor

  ; Create the main widget base
  tlb = widget_base(/column, title = name + ' Spectral Line Viewer', $
    mbar = mbar, /tlb_size_events)

  ; Create control panel
  control_base = widget_base(tlb, /row, /align_center)

  ; Create star navigation buttons
  star_back_button = widget_button(control_base, value = '<< Star')
  star_next_button = widget_button(control_base, value = 'Star >>')

  ; Create star selection dropdown
  star_label = widget_label(control_base, value = 'Star: ')
  star_droplist = widget_droplist(control_base, value = star_names, $
    title = 'Stars:', uvalue = 'STAR_DROPLIST')
  widget_control, star_droplist, set_droplist_select = current_star_index

  ; Create navigation buttons for lines
  back_button = widget_button(control_base, value = '< Line')
  next_button = widget_button(control_base, value = 'Line >')

  ; Create droplist for line selection
  label = widget_label(control_base, value = 'Line: ')
  droplist = widget_droplist(control_base, value = droplist_labels, $
    title = 'Lines:', uvalue = 'DROPLIST')

  ; Create abundance labels
  abund_label = widget_label(control_base, value = 'ABUND = ', /align_left, /dynamic_resize, scr_xsize = 150)
  abunderr_label = widget_label(control_base, value = 'ABUNDERR = ', /align_left, /dynamic_resize, scr_xsize = 150)

  ; Create quit button
  quit_button = widget_button(control_base, value = 'Quit')

  ; Create draw widget for plotting - larger size
  draw = widget_draw(tlb, xsize = 1200, ysize = 800, /button_events)

  ; Realize the widget
  widget_control, tlb, /realize

  ; Get the draw window ID
  widget_control, draw, get_value = draw_id

  ; Create info structure to hold state
  info = { $
    tlb: tlb, $
    name: name, $
    star_back_button: star_back_button, $
    star_next_button: star_next_button, $
    star_droplist: star_droplist, $
    droplist: droplist, $
    back_button: back_button, $
    next_button: next_button, $
    quit_button: quit_button, $
    abund_label: abund_label, $
    abunderr_label: abunderr_label, $
    draw_id: draw_id, $
    current_star_index: current_star_index, $
    current_index: 0, $
    n_lines: n_lines, $
    r: r, $
    hiresall_ptr: ptr_new(hiresall), $
    allcat_ptr: ptr_new(allcat), $
    abundsynth_ptr: ptr_new(abundsynth), $
    hires_ptr: ptr_new(hires), $
    cat_ptr: ptr_new(cat) $
    }

  ; Store info in the top-level base
  widget_control, tlb, set_uvalue = info

  ; Plot the first line
  synth_gui_plot, info

  ; Register with the XMANAGER
  xmanager, 'synth_gui', tlb, /no_block, cleanup = 'synth_gui_cleanup'
end
