; =================================================================
; PROCEDURE: atlas_to_moog
; =================================================================
; PURPOSE:
; Convert ATLAS9 model atmosphere format to MOOG format.
; ATLAS files contain stellar atmosphere structure and abundances;
; MOOG requires a different format for input.
;
; INPUTS:
; infile  - Path to ATLAS9 output file
; outfile - Path for output MOOG atmosphere file (default: 'interp.atm')
; vt      - Microturbulence (km/s). If not set, uses relation:
; v_t = 2.13 - 0.23 * log(g)
;
; OUTPUT FILE FORMAT:
; MOOG 'KURUCZ' format atmosphere with:
; - Header: Teff, log(g), [Fe/H], [alpha/Fe], v_t
; - 72 atmosphere layers (T, P, density, etc.)
; - Alpha element abundances
; - Molecular species list
;
; =================================================================
pro atlas_to_moog, infile, outfile, vt = vt, tweakels = tweakels, tweakabunds = tweakabunds
  compile_opt idl2
  if ~keyword_set(outfile) then outfile = 'interp.atm'

  ; -----------------------------------------------------------------
  ; READ ATLAS9 FILE HEADER
  ; -----------------------------------------------------------------
  openr, lun, infile, /get_lun
  readf, lun, teff, logg, format = '(7X,D5,10X,D7)'
  skip_lun, lun, 3, /lines

  ; Parse element abundances from ATLAS file
  el1 = 0
  el2 = 0
  el3 = 0
  el4 = 0
  el5 = 0
  el6 = 0
  readf, lun, Z, el1, abund1, el2, abund2, format = '(18X,D7,17X,2(1X,I1,1X,D7))'
  feh = alog10(Z) ; Metallicity from Z value
  els = [el1, el2]
  abunds = [abund1, abund2]

  ; Read remaining element abundances (16 more lines with 6 elements each)
  for i = 0, 15 do begin
    readf, lun, el1, abund1, el2, abund2, el3, abund3, el4, abund4, el5, abund5, el6, abund6, $
      format = '(17X,6(1X,I2,1X,D6))'
    els = [els, el1, el2, el3, el4, el5, el6]
    abunds = [abunds, abund1, abund2, abund3, abund4, abund5, abund6]
  endfor

  ; Read final element
  readf, lun, el1, abund1, format = '(17X,1X,I2,1X,D6)'
  els = [els, el1]
  abunds = [abunds, abund1]

  ; Convert abundances to absolute scale: log eps = log(N/N_H) + 12.0
  abunds[2 : n_elements(abunds) - 1] += 12.0 + feh
  close, lun
  free_lun, lun

  if keyword_set(tweakels) and keyword_set(tweakabunds) then begin
    for i = 0, n_elements(tweakels) - 1 do begin
      idx = where(els eq tweakels[i], count)
      if count eq 1 then abunds[idx] = tweakabunds[i]
    endfor
  endif

  ; Set microturbulence if not provided
  if ~keyword_set(vt) then vt = 2.13 - 0.23 * logg

  ; -----------------------------------------------------------------
  ; READ ATMOSPHERE STRUCTURE (72 layers)
  ; -----------------------------------------------------------------
  readcol, infile, c1, c2, c3, c4, c5, c6, c7, $
    format = 'D,D,D,D,D,D,D', numline = 72, skipline = 23, /silent

  ; -----------------------------------------------------------------
  ; CALCULATE [ALPHA/FE]
  ; -----------------------------------------------------------------
  e = elements(/newmoog)
  solar = e.solar - 12.0
  alphaels = [8, 10, 12, 14, 16, 18, 20, 22] ; O, Ne, Mg, Si, S, Ar, Ca, Ti
  nalpha = n_elements(alphaels)

  match, els, alphaels, we, wa
  match, e.atomic, alphaels, w1, w2
  alphafe = mean(abunds[we] - e[w1].solar) - feh ; Mean alpha enhancement

  ; -----------------------------------------------------------------
  ; WRITE MOOG FORMAT FILE
  ; -----------------------------------------------------------------
  openw, lun, outfile, /get_lun
  printf, lun, 'KURUCZ'
  printf, lun, teff, logg, feh, alphafe, vt, $
    format = '(D5.0,"/",D4.2,"/",D+5.2,"/",D+5.2,"/",D4.2)'
  printf, lun, 'ntau=      72'

  ; Write atmosphere structure
  for i = 0, n_elements(c1) - 1 do begin
    printf, lun, c1[i], c2[i], c3[i], c4[i], c5[i], c6[i], vt, $
      format = '(1X,E15.9,2X,F8.1,5(1X,E10.4))'
  endfor

  printf, lun, vt, format = '(E13.3)'

  ; Write alpha element abundances
  printels = [3, 6, 7, 8, 11, 12, 13, 14, 19, 20, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 38, 39, 40, 56, 57, 58, 59, 60, 62, 63, 64, 66, 67, 82]
  for i = 0, n_elements(tweakels) - 1 do begin
    if not contains(printels, tweakels[i]) then printels = [printels, tweakels[i]]
  endfor
  printels = printels[sort(printels)]
  nels = n_elements(printels)
  printf, lun, nels, feh, format = '("NATOMS",4X,I2,2X,D8.4)'
  for i = 0, nels - 1 do begin
    printf, lun, printels[i], abunds[where(els eq printels[i])], format = '("      ",I2,"    ",D8.4)'
  endfor

  ; Write molecular species list
  printf, lun, 'NMOL       18'
  printf, lun, '101.0   106.0   107.0   108.0   606.0   607.0   608.0   707.0'
  printf, lun, '708.0   808.0 10108.0 60808.0     6.1     7.1     8.1    22.1'
  printf, lun, ' 23.1   823.0'
  close, lun
  free_lun, lun
end
