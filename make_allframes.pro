pro make_allframes
  compile_opt idl2
  newstr = {dm: 0d, ebv: 0d, feh: 0d, gmag: 0d, bpmag: 0d, rpmag: 0d, gmag0: 0d, bpmag0: 0d, rpmag0: 0d, gmagerr: 0d, bpmagerr: 0d, rpmagerr: 0d, teff_mb20: 0d, teff_mb20_err: 0d, teff_mb20_color: ' ', teffphot: 0d, teffphoterr: 0d, loggphot: 0d, loggphoterr: 0d}

  gaia_m15 = mrdfits('gaia_m15.fits', 1, /silent)
  gaia_m15 = struct_addtags(gaia_m15, replicate(newstr, n_elements(gaia_m15)))
  gaia_m15.dm = 15.42  ;apparent distance modulus (Vandenberg et al. 2016)
  gaia_m15.ebv = 0.10  ;(Vandenberg et al. 2016)
  gaia_m15.feh = -2.37

  gaia_m92 = mrdfits('gaia_m92.fits', 1, /silent)
  gaia_m92 = struct_addtags(gaia_m92, replicate(newstr, n_elements(gaia_m92)))
  gaia_m92.dm = 14.74   ;apparent distance modulus (Vandenberg et al. 2016)
  gaia_m92.ebv = 0.023  ;apparent distance modulus (Vandenberg et al. 2016)
  gaia_m92.feh = -2.31

  gaia = struct_append(gaia_m15, gaia_m92)
  newstr = replicate({name: '', fullname: '', makeefile: '', specfile: ''}, n_elements(gaia))
  gaia = struct_addtags(newstr, gaia)

  readcol, getenv('CALTECH') + 'keck/hires/2022aug13/Kirby.HIRES.2022aug13.starlist', objname, r1, r2, r3, d1, d2, d3, source_id, skip = 4, format = 'A,I,I,D,I,I,D,X,X,X,X,X,L'
  n = n_elements(r1)
  ra = 15d * (double(r1) + double(r2) / 60. + r3 / 3600.)
  dec = double(d1) + double(d2) / 60. + d3 / 3600.
  spherematch, gaia.ra, gaia.dec, ra, dec, 0.1 / 3600., w1, w2
  gaia[w1].fullname = objname[w2]
  gaia = gaia[w1]
  for i = 0, n - 1 do begin
    namearray = strsplit(gaia[i].fullname, '_', /extract)
    gc = namearray[0]
    index = fix(namearray[2])
    gaia[i].name = gc + '-star-' + strtrim(index, 2)
  endfor
  gaia.makeefile = getenv('chome') + 'keck/hires/M15_M92/makee/' + strtrim(gaia.name, 2)
  gaia.specfile = getenv('chome') + 'caltech/hires/M15_M92/spectra/' + strtrim(gaia.name, 2) + '.fits'

  ; f = file_search(getenv('chome')+'keck/hires/M15_M92/makee/*_1.fits', count=c)
  ; for i=0,c-1 do begin
  ; hdr = headfits(f[i])
  ; ras = sxpar(hdr, 'RA')
  ; decs = sxpar(hdr, 'DEC')
  ; get_coords, coords, instring=ras+' '+decs
  ; ra = coords[0]*15.
  ; dec = coords[1]
  ; spherematch, gaia.ra, gaia.dec, ra, dec, 1.0/3600., w1, w2
  ; if w1[0] eq -1 then message, 'I did not find this star.'
  ; names = strsplit(sxpar(hdr, 'OBJECT'), ' ', /extract)
  ; gaia[w1[0]].name = names[0]
  ; endfor

  gaia.gmag = gaia.phot_g_mean_mag
  gaia.bpmag = gaia.phot_bp_mean_mag
  gaia.rpmag = gaia.phot_rp_mean_mag

  gaia.gmagerr = 2.5 / alog(10.) * gaia.phot_g_mean_flux_error / gaia.phot_g_mean_flux
  gaia.bpmagerr = 2.5 / alog(10.) * gaia.phot_bp_mean_flux_error / gaia.phot_bp_mean_flux
  gaia.rpmagerr = 2.5 / alog(10.) * gaia.phot_rp_mean_flux_error / gaia.phot_rp_mean_flux

  ; Gaia Collaboration, Babusiaux et al. 2018, A&A, 616, A10
  kG = [0.9761, -0.1704, 0.0086, 0.0011, -0.0438, 0.0013, 0.0099]
  kBP = [1.1517, -0.0871, -0.0333, 0.0173, -0.0230, 0.0006, 0.0043]
  kRP = [0.6104, -0.0170, -0.0026, -0.0017, -0.0078, 0.00005, 0.0006]

  ; Mucciarelli & Bellazzini 2020, RNAAS
  bprp_dwarf = [0.4929, 0.5092, -0.0353, 0.0192, -0.0020, -0.0395]
  bpk_dwarf = [0.5342, 0.2044, -0.0021, 0.0276, 0.0005, -0.0158]
  bprp_giant = [0.5323, 0.4775, -0.0344, -0.0110, -0.0020, -0.0009]
  bpk_giant = [0.5668, 0.1890, -0.0017, 0.0065, -0.0008, -0.0045]

  for i = 0, n_elements(gaia) - 1 do begin
    ; Gaia Collaboration, Babusiaux et al. 2018, A&A, 616, A10
    vars = [1d, gaia[i].bpmag - gaia[i].rpmag, (gaia[i].bpmag - gaia[i].rpmag) ^ 2., (gaia[i].bpmag - gaia[i].rpmag) ^ 3., 3.1 * gaia[i].ebv, (3.1 * gaia[i].ebv) ^ 2., (gaia[i].bpmag - gaia[i].rpmag) * 3.1 * gaia[i].ebv]
    gaia[i].gmag0 = gaia[i].phot_g_mean_mag - total(kG * vars) * 3.1 * gaia[i].ebv
    gaia[i].bpmag0 = gaia[i].bpmag - total(kBP * vars) * 3.1 * gaia[i].ebv
    gaia[i].rpmag0 = gaia[i].rpmag - total(kRP * vars) * 3.1 * gaia[i].ebv

    ; Mucciarelli & Bellazzini 2020, RNAAS
    if gaia[i].gmag0 gt 17.5 then begin
      color = gaia[i].bpmag0 - gaia[i].rpmag0
      colorerr = sqrt(gaia[i].bpmagerr ^ 2. + gaia[i].rpmagerr ^ 2.)
      vars = [1d, color, (color) ^ 2., gaia[i].feh, gaia[i].feh ^ 2., gaia[i].feh * (color)]
      varserr = [0d, 1d, 2. * (color), 0d, 0d, gaia[i].feh]
      gaia[i].teff_mb20 = 5040d / total(bprp_dwarf * vars)
      gaia[i].teff_mb20_err = sqrt((gaia[i].teff_mb20 * total(bprp_dwarf * varserr * colorerr) / total(bprp_dwarf * vars)) ^ 2. + 61d ^ 2.)
      gaia[i].teff_mb20_color = 'BP-RP(dwarf)'
    endif else begin
      color = gaia[i].bpmag0 - gaia[i].rpmag0
      colorerr = sqrt(gaia[i].bpmagerr ^ 2. + gaia[i].rpmagerr ^ 2.)
      vars = [1d, color, (color) ^ 2., gaia[i].feh, gaia[i].feh ^ 2., gaia[i].feh * (color)]
      varserr = [0d, 1d, 2. * (color), 0d, 0d, gaia[i].feh]
      gaia[i].teff_mb20 = 5040d / total(bprp_giant * vars)
      gaia[i].teff_mb20_err = sqrt((gaia[i].teff_mb20 * total(bprp_giant * varserr * colorerr) / total(bprp_giant * vars)) ^ 2. + 83d ^ 2.)
      gaia[i].teff_mb20_color = 'BP-RP(giant)'
    endelse
  endfor

  gaia.teffphot = gaia.teff_mb20
  gaia.teffphoterr = gaia.teff_mb20_err

  ; Andrae et al. (2018, A&A, 616, A8)
  teffdiff = gaia.teffphot - 5772d
  bcG = poly(teffdiff, [6d-2, 6.731d-5, -6.647d-8, 2.859d-11, -7.197d-15])
  logL = (gaia.gmag0 - (gaia.dm - 2.682 * gaia.ebv) + bcG - 4.68) / (-2.5d) + alog10(3.828d33) ; absolute G_sun = 4.68 (Andrae et al. 2018)

  sigma_SB = 5.6704d-5
  G = 6.674d-8
  Msun = 1.989d33
  M = 0.75 * Msun
  Merr = 0.1
  gaia.loggphot = alog10(4 * !dpi * sigma_SB * G) + alog10(M) + 4. * alog10(gaia.teffphot) - logL
  gaia.loggphoterr = sqrt(Merr ^ 2. + (4. * gaia.teffphoterr / (gaia.teffphot * alog(10.))) ^ 2. + (gaia.gmagerr / 2.5) ^ 2.)

  mwrfits, gaia, 'M15_M92_allframes.fits', /create
  stop
end
