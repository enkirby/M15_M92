; =================================================================
; PROCEDURE: make_allframes
; =================================================================
; PURPOSE:
; Create master catalog combining Gaia photometry/astrometry with
; HIRES spectroscopic observations for M15 and M92 stars.
; Calculates photometric temperatures and surface gravities.
;
; INPUTS:
; Requires these files in current directory:
; - gaia_m15.fits : Gaia data for M15 stars
; - gaia_m92.fits : Gaia data for M92 stars
; - Kirby.HIRES.2022aug13.starlist : HIRES target list
;
; OUTPUT:
; M15_M92_allframes.fits - Combined catalog with:
; - Gaia photometry (G, BP, RP) and errors
; - Dereddened magnitudes
; - Photometric Teff (Mucciarelli & Bellazzini 2020)
; - Photometric log(g) from mass-luminosity relation
; - Cross-matched star names and file paths
;
; REFERENCES:
; - Vandenberg et al. 2016 : Distance moduli, reddening
; - Gaia Collaboration, Babusiaux et al. 2018 : Extinction coefficients
; - Mucciarelli & Bellazzini 2020 : Teff calibrations
; - Andrae et al. 2018 : Bolometric corrections
;
; =================================================================
pro make_allframes
  compile_opt idl2

  ; -----------------------------------------------------------------
  ; INITIALIZE STRUCTURE FOR NEW FIELDS
  ; -----------------------------------------------------------------
  newstr = {dm: 0d, ebv: 0d, feh: 0d, $
    gmag: 0d, bpmag: 0d, rpmag: 0d, $
    gmag0: 0d, bpmag0: 0d, rpmag0: 0d, $
    gmagerr: 0d, bpmagerr: 0d, rpmagerr: 0d, $
    phot_bp_rp_excess_factor: 0d, c_star: 0d, $
    teff_mb20_bprp: 0d, teff_mb20_bprp_err: 0d, teff_mb20_bprp_color: ' ', $
    teff_mb20_bpg: 0d, teff_mb20_bpg_err: 0d, teff_mb20_bpg_color: ' ', $
    teff_mb20_grp: 0d, teff_mb20_grp_err: 0d, teff_mb20_grp_color: ' ', $
    teff_mb20_bpk: 0d, teff_mb20_bpk_err: 0d, teff_mb20_bpk_color: ' ', $
    teff_mb20_rpk: 0d, teff_mb20_rpk_err: 0d, teff_mb20_rpk_color: ' ', $
    teff_mb20_gk: 0d, teff_mb20_gk_err: 0d, teff_mb20_gk_color: ' ', $
    teff_mb20_color: ' ', $
    teffphot: 0d, teffphoterr: 0d, $
    loggphot: 0d, loggphoterr: 0d}

  ; -----------------------------------------------------------------
  ; LOAD M15 DATA
  ; -----------------------------------------------------------------
  gaiadr3_m15 = mrdfits(getenv('CALTECH') + 'gc/n7078/GaiaDR3_M15.fits.gz', 1, /silent)
  gaiadr3_m15 = gaiadr3_m15[where(gaiadr3_m15.phot_g_mean_mag gt 12.5 and gaiadr3_m15.phot_g_mean_mag lt 17.5)]
  gaia_m15 = mrdfits('gaia_m15.fits', 1, /silent)
  spherematch, gaiadr3_m15.ra, gaiadr3_m15.dec, gaia_m15.ra, gaia_m15.dec, 0.5 / 3600., w1, w2

  gaia_m15 = struct_addtags(gaia_m15, replicate(newstr, n_elements(gaia_m15)))
  gaia_m15[w2].phot_bp_rp_excess_factor = gaiadr3_m15[w1].phot_bp_rp_excess_factor
  gaia_m15.dm = 15.42 ; Apparent distance modulus (Vandenberg et al. 2016)
  gaia_m15.ebv = 0.10 ; E(B-V) reddening (Vandenberg et al. 2016)
  gaia_m15.feh = -2.37 ; Cluster metallicity

  ; -----------------------------------------------------------------
  ; LOAD M92 DATA
  ; -----------------------------------------------------------------
  gaiadr3_m92 = mrdfits(getenv('CALTECH') + 'gc/n6341/GaiaDR3_M92.fits.gz', 1, /silent)
  gaiadr3_m92 = gaiadr3_m92[where(gaiadr3_m92.phot_g_mean_mag gt 12.5 and gaiadr3_m92.phot_g_mean_mag lt 17.5)]
  gaia_m92 = mrdfits('gaia_m92.fits', 1, /silent)
  spherematch, gaiadr3_m92.ra, gaiadr3_m92.dec, gaia_m92.ra, gaia_m92.dec, 0.5 / 3600., w1, w2

  gaia_m92 = struct_addtags(gaia_m92, replicate(newstr, n_elements(gaia_m92)))
  gaia_m92[w2].phot_bp_rp_excess_factor = gaiadr3_m92[w1].phot_bp_rp_excess_factor
  gaia_m92.dm = 14.74 ; Apparent distance modulus (Vandenberg et al. 2016)
  gaia_m92.ebv = 0.023 ; E(B-V) reddening (Vandenberg et al. 2016)
  gaia_m92.feh = -2.31 ; Cluster metallicity

  ; -----------------------------------------------------------------
  ; COMBINE CATALOGS AND ADD NAME/FILE FIELDS
  ; -----------------------------------------------------------------
  gaia = struct_append(gaia_m15, gaia_m92)
  newstr = replicate({name: '', fullname: '', makeefile: '', specfile: ''}, $
    n_elements(gaia))
  gaia = struct_addtags(newstr, gaia)

  ; -----------------------------------------------------------------
  ; CROSS-MATCH WITH HIRES TARGET LIST
  ; -----------------------------------------------------------------
  ; Read starlist with RA/Dec and Gaia source IDs
  readcol, getenv('CALTECH') + 'keck/hires/2022aug13/Kirby.HIRES.2022aug13.starlist', $
    objname, r1, r2, r3, d1, d2, d3, source_id, skip = 4, $
    format = 'A,I,I,D,I,I,D,X,X,X,X,X,L'
  n = n_elements(r1)

  ; Convert sexagesimal to decimal degrees
  ra = 15d * (double(r1) + double(r2) / 60. + r3 / 3600.)
  dec = double(d1) + double(d2) / 60. + d3 / 3600.

  ; Match within 0.1 arcsec
  spherematch, gaia.ra, gaia.dec, ra, dec, 0.1 / 3600., w1, w2
  gaia[w1].fullname = objname[w2]
  gaia = gaia[w1]

  ; Parse star names from full names (e.g., "M15_RGB_1" -> "M15-star-1")
  for i = 0, n - 1 do begin
    namearray = strsplit(gaia[i].fullname, '_', /extract)
    gc = namearray[0] ; Cluster name (M15 or M92)
    index = fix(namearray[2]) ; Star index number
    gaia[i].name = gc + '-star-' + strtrim(index, 2)
  endfor

  ; Set file paths
  ; Set file paths
  gaia.makeefile = getenv('chome') + 'keck/hires/M15_M92/makee/' + strtrim(gaia.name, 2)
  gaia.specfile = getenv('chome') + 'caltech/hires/M15_M92_bprp/spectra/' + $
    strtrim(gaia.name, 2) + '.fits'

  ; -----------------------------------------------------------------
  ; LOAD 2MASS CATALOG FOR CROSS-MATCHING
  ; -----------------------------------------------------------------
  newstr = replicate({kmag: 0d, kmag0: 0d, kmagerr: 0d}, $
    n_elements(gaia))
  gaia = struct_addtags(newstr, gaia)
  tmass = mrdfits('M15_M92_catalog_normal_2mass.fits', 1, /silent)
  match, tmass.name, gaia.name, wtmass, wgaia
  gaia[wgaia].kmag = tmass[wtmass].kmag_2Mass
  gaia[wgaia].kmag0 = gaia[wgaia].kmag - 0.310 * gaia[wgaia].ebv
  gaia[wgaia].kmagerr = tmass[wtmass].e_kmag_2Mass

  ; -----------------------------------------------------------------
  ; EXTRACT GAIA PHOTOMETRY
  ; -----------------------------------------------------------------
  gaia.gmag = gaia.phot_g_mean_mag
  gaia.bpmag = gaia.phot_bp_mean_mag
  gaia.rpmag = gaia.phot_rp_mean_mag

  ; Convert flux errors to magnitude errors
  gaia.gmagerr = 2.5 / alog(10.) * gaia.phot_g_mean_flux_error / gaia.phot_g_mean_flux
  gaia.bpmagerr = 2.5 / alog(10.) * gaia.phot_bp_mean_flux_error / gaia.phot_bp_mean_flux
  gaia.rpmagerr = 2.5 / alog(10.) * gaia.phot_rp_mean_flux_error / gaia.phot_rp_mean_flux

  bprp = gaia.bpmag - gaia.rpmag
  gaia.c_star = gaia.phot_bp_rp_excess_factor - (1.162004 + 0.011464 * bprp + 0.049255 * bprp ^ 2. - 0.005879 * bprp ^ 3.)
  ; -----------------------------------------------------------------
  ; EXTINCTION COEFFICIENTS
  ; -----------------------------------------------------------------
  ; Gaia Collaboration, Babusiaux et al. 2018, A&A, 616, A10
  ; Extinction in each Gaia band as function of BP-RP color and E(B-V)
  kG = [0.9761, -0.1704, 0.0086, 0.0011, -0.0438, 0.0013, 0.0099]
  kBP = [1.1517, -0.0871, -0.0333, 0.0173, -0.0230, 0.0006, 0.0043]
  kRP = [0.6104, -0.0170, -0.0026, -0.0017, -0.0078, 0.00005, 0.0006]

  ; -----------------------------------------------------------------
  ; TEMPERATURE CALIBRATION COEFFICIENTS
  ; -----------------------------------------------------------------
  ; Mucciarelli, Bellazzini & Masari 2021, A&A
  ; Photometric Teff relations for metal-poor stars
  bprp_dwarf = [0.4929, 0.5092, -0.0353, 0.0192, -0.0020, -0.0395] ; For dwarfs
  bpg_dwarf = [0.5316, 1.2452, -0.4677, 0.0068, -0.0031, -0.0752]
  grp_dwarf = [0.5050, 0.6532, 0.2284, 0.0260, -0.0011, -0.0726]
  bpk_dwarf = [0.5342, 0.2044, -0.0021, 0.0276, 0.0005, -0.0158]
  rpk_dwarf = [0.5526, 0.3712, -0.0121, 0.0330, 0.0029, -0.0220]
  gk_dwarf = [0.5351, 0.2440, 0.0016, 0.0289, 0.0015, -0.0163]

  bprp_giant = [0.5323, 0.4775, -0.0344, -0.0110, -0.0020, -0.0009] ; For giants
  bpg_giant = [0.5701, 1.1188, -0.3710, -0.0236, -0.0039, 0.0070]
  grp_giant = [0.5472, 0.5914, 0.2347, -0.0019, -0.0012, 0.0060]
  bpk_giant = [0.5668, 0.1890, -0.0017, 0.0065, -0.0008, -0.0045]
  rpk_giant = [0.5774, 0.3637, -0.0226, 0.0346, 0.0007, -0.0221]
  gk_giant = [0.5569, 0.2436, -0.0035, 0.0211, 0.0007, -0.0089]

  ; -----------------------------------------------------------------
  ; CALCULATE DEREDDENED MAGNITUDES AND PHOTOMETRIC TEFF
  ; -----------------------------------------------------------------
  for i = 0, n_elements(gaia) - 1 do begin
    ; Deredenning using Gaia extinction coefficients
    ; Vars = [1, BP-RP, (BP-RP)^2, (BP-RP)^3, A_V, A_V^2, (BP-RP)*A_V]
    vars = [1d, gaia[i].bpmag - gaia[i].rpmag, $
      (gaia[i].bpmag - gaia[i].rpmag) ^ 2., $
      (gaia[i].bpmag - gaia[i].rpmag) ^ 3., $
      3.1 * gaia[i].ebv, $
      (3.1 * gaia[i].ebv) ^ 2., $
      (gaia[i].bpmag - gaia[i].rpmag) * 3.1 * gaia[i].ebv]
    gaia[i].gmag0 = gaia[i].phot_g_mean_mag - total(kG * vars) * 3.1 * gaia[i].ebv
    gaia[i].bpmag0 = gaia[i].bpmag - total(kBP * vars) * 3.1 * gaia[i].ebv
    gaia[i].rpmag0 = gaia[i].rpmag - total(kRP * vars) * 3.1 * gaia[i].ebv

    ; Photometric temperature from color
    ; Use dwarf relation for faint stars (G > 17.5), giant relation otherwise
    if gaia[i].gmag0 gt 17.5 then begin
      ; DWARF calibration
      color = gaia[i].bpmag0 - gaia[i].rpmag0
      colorerr = sqrt(gaia[i].bpmagerr ^ 2. + gaia[i].rpmagerr ^ 2.)
      vars = [1d, color, (color) ^ 2., gaia[i].feh, gaia[i].feh ^ 2., gaia[i].feh * (color)]
      varserr = [0d, 1d, 2. * (color), 0d, 0d, gaia[i].feh]
      gaia[i].teff_mb20_bprp = 5040d / total(bprp_dwarf * vars)
      gaia[i].teff_mb20_bprp_err = sqrt((gaia[i].teff_mb20_bprp * total(bprp_dwarf * varserr * colorerr) / $
        total(bprp_dwarf * vars)) ^ 2. + 61d ^ 2.)
      gaia[i].teff_mb20_bprp_color = 'BP-RP(dwarf)'

      color = gaia[i].bpmag0 - gaia[i].gmag0
      colorerr = sqrt(gaia[i].bpmagerr ^ 2. + gaia[i].gmagerr ^ 2.)
      vars = [1d, color, (color) ^ 2., gaia[i].feh, gaia[i].feh ^ 2., gaia[i].feh * (color)]
      varserr = [0d, 1d, 2. * (color), 0d, 0d, gaia[i].feh]
      gaia[i].teff_mb20_bpg = 5040d / total(bpg_dwarf * vars)
      gaia[i].teff_mb20_bpg_err = sqrt((gaia[i].teff_mb20_bpg * total(bpg_dwarf * varserr * colorerr) / $
        total(bpg_dwarf * vars)) ^ 2. + 58d ^ 2.)
      gaia[i].teff_mb20_bpg_color = 'BP-G(dwarf)'

      color = gaia[i].gmag0 - gaia[i].rpmag0
      colorerr = sqrt(gaia[i].gmagerr ^ 2. + gaia[i].rpmagerr ^ 2.)
      vars = [1d, color, (color) ^ 2., gaia[i].feh, gaia[i].feh ^ 2., gaia[i].feh * (color)]
      varserr = [0d, 1d, 2. * (color), 0d, 0d, gaia[i].feh]
      gaia[i].teff_mb20_grp = 5040d / total(grp_dwarf * vars)
      gaia[i].teff_mb20_grp_err = sqrt((gaia[i].teff_mb20_grp * total(grp_dwarf * varserr * colorerr) / $
        total(grp_dwarf * vars)) ^ 2. + 62d ^ 2.)
      gaia[i].teff_mb20_grp_color = 'G-RP(dwarf)'

      if gaia[i].kmagerr gt 0.0 and finite(gaia[i].kmagerr) then begin
        color = gaia[i].bpmag0 - gaia[i].kmag0
        colorerr = sqrt(gaia[i].bpmagerr ^ 2. + gaia[i].kmagerr ^ 2.)
        vars = [1d, color, (color) ^ 2., gaia[i].feh, gaia[i].feh ^ 2., gaia[i].feh * (color)]
        varserr = [0d, 1d, 2. * (color), 0d, 0d, gaia[i].feh]
        gaia[i].teff_mb20_bpk = 5040d / total(bpk_dwarf * vars)
        gaia[i].teff_mb20_bpk_err = sqrt((gaia[i].teff_mb20_bpk * total(bpk_dwarf * varserr * colorerr) / $
          total(bpk_dwarf * vars)) ^ 2. + 44d ^ 2.)
        gaia[i].teff_mb20_bpk_color = 'BP-Ks(dwarf)'

        color = gaia[i].rpmag0 - gaia[i].kmag0
        colorerr = sqrt(gaia[i].rpmagerr ^ 2. + gaia[i].kmagerr ^ 2.)
        vars = [1d, color, (color) ^ 2., gaia[i].feh, gaia[i].feh ^ 2., gaia[i].feh * (color)]
        varserr = [0d, 1d, 2. * (color), 0d, 0d, gaia[i].feh]
        gaia[i].teff_mb20_rpk = 5040d / total(rpk_dwarf * vars)
        gaia[i].teff_mb20_rpk_err = sqrt((gaia[i].teff_mb20_rpk * total(rpk_dwarf * varserr * colorerr) / $
          total(rpk_dwarf * vars)) ^ 2. + 53d ^ 2.)
        gaia[i].teff_mb20_rpk_color = 'RP-Ks(dwarf)'

        color = gaia[i].gmag0 - gaia[i].kmag0
        colorerr = sqrt(gaia[i].gmagerr ^ 2. + gaia[i].kmagerr ^ 2.)
        vars = [1d, color, (color) ^ 2., gaia[i].feh, gaia[i].feh ^ 2., gaia[i].feh * (color)]
        varserr = [0d, 1d, 2. * (color), 0d, 0d, gaia[i].feh]
        gaia[i].teff_mb20_gk = 5040d / total(gk_dwarf * vars)
        gaia[i].teff_mb20_gk_err = sqrt((gaia[i].teff_mb20_gk * total(gk_dwarf * varserr * colorerr) / $
          total(gk_dwarf * vars)) ^ 2. + 48d ^ 2.)
        gaia[i].teff_mb20_gk_color = 'G-Ks(dwarf)'
      endif else begin
        gaia[i].teff_mb20_bpk = -9999d
        gaia[i].teff_mb20_bpk_err = -9999d

        gaia[i].teff_mb20_rpk = -9999d
        gaia[i].teff_mb20_rpk_err = -9999d

        gaia[i].teff_mb20_gk = -9999d
        gaia[i].teff_mb20_gk_err = -9999d
      endelse
    endif else begin
      ; GIANT calibration
      color = gaia[i].bpmag0 - gaia[i].rpmag0
      colorerr = sqrt(gaia[i].bpmagerr ^ 2. + gaia[i].rpmagerr ^ 2.)
      vars = [1d, color, (color) ^ 2., gaia[i].feh, gaia[i].feh ^ 2., gaia[i].feh * (color)]
      varserr = [0d, 1d, 2. * (color), 0d, 0d, gaia[i].feh]
      gaia[i].teff_mb20_bprp = 5040d / total(bprp_giant * vars)
      gaia[i].teff_mb20_bprp_err = sqrt((gaia[i].teff_mb20_bprp * total(bprp_giant * varserr * colorerr) / $
        total(bprp_giant * vars)) ^ 2. + 83d ^ 2.)
      gaia[i].teff_mb20_bprp_color = 'BP-RP(giant)'

      color = gaia[i].bpmag0 - gaia[i].gmag0
      colorerr = sqrt(gaia[i].bpmagerr ^ 2. + gaia[i].gmagerr ^ 2.)
      vars = [1d, color, (color) ^ 2., gaia[i].feh, gaia[i].feh ^ 2., gaia[i].feh * (color)]
      varserr = [0d, 1d, 2. * (color), 0d, 0d, gaia[i].feh]
      gaia[i].teff_mb20_bpg = 5040d / total(bpg_giant * vars)
      gaia[i].teff_mb20_bpg_err = sqrt((gaia[i].teff_mb20_bpg * total(bpg_giant * varserr * colorerr) / $
        total(bpg_giant * vars)) ^ 2. + 83d ^ 2.)
      gaia[i].teff_mb20_bpg_color = 'BP-G(giant)'

      color = gaia[i].gmag0 - gaia[i].rpmag0
      colorerr = sqrt(gaia[i].gmagerr ^ 2. + gaia[i].rpmagerr ^ 2.)
      vars = [1d, color, (color) ^ 2., gaia[i].feh, gaia[i].feh ^ 2., gaia[i].feh * (color)]
      varserr = [0d, 1d, 2. * (color), 0d, 0d, gaia[i].feh]
      gaia[i].teff_mb20_grp = 5040d / total(grp_giant * vars)
      gaia[i].teff_mb20_grp_err = sqrt((gaia[i].teff_mb20_grp * total(grp_giant * varserr * colorerr) / $
        total(grp_giant * vars)) ^ 2. + 71d ^ 2.)
      gaia[i].teff_mb20_grp_color = 'G-RP(giant)'

      if gaia[i].kmagerr gt 0.0 and finite(gaia[i].kmagerr) then begin
        color = gaia[i].bpmag0 - gaia[i].kmag0
        colorerr = sqrt(gaia[i].bpmagerr ^ 2. + gaia[i].kmagerr ^ 2.)
        vars = [1d, color, (color) ^ 2., gaia[i].feh, gaia[i].feh ^ 2., gaia[i].feh * (color)]
        varserr = [0d, 1d, 2. * (color), 0d, 0d, gaia[i].feh]
        gaia[i].teff_mb20_bpk = 5040d / total(bpk_giant * vars)
        gaia[i].teff_mb20_bpk_err = sqrt((gaia[i].teff_mb20_bpk * total(bpk_giant * varserr * colorerr) / $
          total(bpk_giant * vars)) ^ 2. + 49d ^ 2.)
        gaia[i].teff_mb20_bpk_color = 'BP-Ks(giant)'

        color = gaia[i].rpmag0 - gaia[i].kmag0
        colorerr = sqrt(gaia[i].rpmagerr ^ 2. + gaia[i].kmagerr ^ 2.)
        vars = [1d, color, (color) ^ 2., gaia[i].feh, gaia[i].feh ^ 2., gaia[i].feh * (color)]
        varserr = [0d, 1d, 2. * (color), 0d, 0d, gaia[i].feh]
        gaia[i].teff_mb20_rpk = 5040d / total(rpk_giant * vars)
        gaia[i].teff_mb20_rpk_err = sqrt((gaia[i].teff_mb20_rpk * total(rpk_giant * varserr * colorerr) / $
          total(rpk_giant * vars)) ^ 2. + 61d ^ 2.)
        gaia[i].teff_mb20_rpk_color = 'RP-Ks(giant)'

        color = gaia[i].gmag0 - gaia[i].kmag0
        colorerr = sqrt(gaia[i].gmagerr ^ 2. + gaia[i].kmagerr ^ 2.)
        vars = [1d, color, (color) ^ 2., gaia[i].feh, gaia[i].feh ^ 2., gaia[i].feh * (color)]
        varserr = [0d, 1d, 2. * (color), 0d, 0d, gaia[i].feh]
        gaia[i].teff_mb20_gk = 5040d / total(gk_giant * vars)
        gaia[i].teff_mb20_gk_err = sqrt((gaia[i].teff_mb20_gk * total(gk_giant * varserr * colorerr) / $
          total(gk_giant * vars)) ^ 2. + 46d ^ 2.)
        gaia[i].teff_mb20_gk_color = 'G-Ks(giant)'
      endif else begin
        gaia[i].teff_mb20_bpk = -9999d
        gaia[i].teff_mb20_bpk_err = -9999d

        gaia[i].teff_mb20_rpk = -9999d
        gaia[i].teff_mb20_rpk_err = -9999d

        gaia[i].teff_mb20_gk = -9999d
        gaia[i].teff_mb20_gk_err = -9999d
      endelse
    endelse

    ; if gaia[i].teff_mb20_bpk gt 0 and gaia[i].teff_mb20_bpk_err gt 0 and gaia[i].teff_mb20_bpk_err lt gaia[i].teff_mb20_bprp_err then begin
    ; gaia[i].teffphot = gaia[i].teff_mb20_bpk
    ; gaia[i].teffphoterr = gaia[i].teff_mb20_bpk_err
    ; gaia[i].teff_mb20_color = gaia[i].teff_mb20_bpk_color
    ; endif else begin
    gaia[i].teffphot = gaia[i].teff_mb20_bprp
    gaia[i].teffphoterr = gaia[i].teff_mb20_bprp_err
    gaia[i].teff_mb20_color = gaia[i].teff_mb20_bprp_color
    ; endelse
  endfor

  ; -----------------------------------------------------------------
  ; BOLOMETRIC CORRECTIONS AND SURFACE GRAVITY
  ; -----------------------------------------------------------------
  ; Andrae et al. (2018, A&A, 616, A8)
  ; Bolometric correction as polynomial function of Teff - Tsun
  teffdiff = gaia.teffphot - 5772d
  bcG = poly(teffdiff, [6d-2, 6.731d-5, -6.647d-8, 2.859d-11, -7.197d-15])

  ; Bolometric luminosity from dereddened G magnitude
  ; L = 10^(-0.4*(M_G - M_G,sun)) where M_G,sun = 4.68 (Andrae et al. 2018)
  logL = (gaia.gmag0 - (gaia.dm - 2.682 * gaia.ebv) + bcG - 4.68) / (-2.5d) + $
    alog10(3.828d33) ; L_sun in erg/s

  ; Surface gravity from luminosity-mass-Teff relation
  ; g = G*M/R^2, where L = 4*pi*R^2*sigma*Teff^4
  ; log(g) = log(4*pi*sigma*G) + log(M) + 4*log(Teff) - log(L)
  sigma_SB = 5.6704d-5 ; Stefan-Boltzmann constant (cgs)
  G = 6.674d-8 ; Gravitational constant (cgs)
  Msun = 1.989d33 ; Solar mass (g)

  ; BaSTI isochrones, 12 Gyr, [alpha/Fe] = +0.4, [Fe/H] = -2.448, RGB stars with the appropriate G mag have masses 0.80 Msun at Y=0.247 and 0.76 Msun at Y=0.275
  M = 0.80 * Msun
  Merr = 0.05 ; 0.05 Msun mass uncertainty

  gaia.loggphot = alog10(4 * !dpi * sigma_SB * G) + alog10(M) + $
    4. * alog10(gaia.teffphot) - logL
  gaia.loggphoterr = sqrt(Merr ^ 2. + $
    (4. * gaia.teffphoterr / (gaia.teffphot * alog(10.))) ^ 2. + $
    (gaia.gmagerr / 2.5) ^ 2.)

  ; -----------------------------------------------------------------
  ; WRITE OUTPUT FILE
  ; -----------------------------------------------------------------
  ; Save combined M15+M92 catalog with photometric parameters
  mwrfits, gaia, 'M15_M92_allframes.fits', /create

  ; Halt for inspection
  stop
end
