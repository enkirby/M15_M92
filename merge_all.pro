function findpix, lambda, lambda0
  diff = abs(lambda - lambda0)
  minval = min(diff, minsub)
  minsub = (minsub - 1) > 0
  return, minsub
end

function zfind, lambdaref, specref, ivarref, lambda, spec, ivar, zmin=zmin, zmax=zmax, $
                zguess=zguess, pwidth=pwidth, nfind=nfind, width=width, $
                wvmin=wvmin, wvmax=wvmax, zerr=zerr

    common com_zfind, starflux, starloglam0, stardloglam, $
        nstars, npoints, starflux_corrected, thisfile


    starloglam0 = alog10(lambdaref[0])
    stardloglam = alog10(lambdaref[1]) - starloglam0
    npoints = n_elements(lambdaref)
    starorder =  round(3.*npoints/7000.) ;3 terms/7000 points 
    
    w = where(lambda ge min(lambdaref) and lambda le max(lambdaref) and ivar gt 0)
    x_specrebin, lambda[w], spec[w], lambdaref, objflux, var=ivar[w]^(-1), nwvar=varrebin, /silent
    objivar = varrebin^(-1)
    w = where(lambdaref ge min(lambda) and lambdaref le max(lambda) and objivar gt 0, npix)
    objflux = objflux[w]
    objivar = objivar[w]
    objloglam0 = alog10(lambdaref[w[0]])
    objdloglam = stardloglam
    loglam = objloglam0 + dindgen(npix)*objdloglam
    
; check if we need to trim the input spectrum to just a subregion.
    if n_elements(wvmin) gt 0 then wvmin = wvmin[0] else wvmin = -1
    if n_elements(wvmax) gt 0 then wvmax = wvmax[0] else wvmax = -1
    if wvmax ge 0 and wvmin ge 0 then begin
        minpix = findpix(loglam, alog10(wvmin))
        maxpix = findpix(loglam, alog10(wvmax))
        loglam = loglam[minpix:maxpix]
        objflux = objflux[minpix:maxpix]
        objivar = objivar[minpix:maxpix]
    endif

;;; CHECK IF THE zmin AND zmax ARGUMENTS WERE PASSED. IF SO, THEN
;;; CONVERT THE REDSHIFT VALUES INTO PIXEL VALUES pmin AND pmax.
;;; THIS IS ONLY TRUE IF objloglam0 = temploglam0?
  IF n_elements(zmin) NE 0 THEN $
    pmin = FLOOR( ALOG10(1.0 + zmin) / objdloglam )
  IF n_elements(zmax) NE 0 THEN $
    pmax = CEIL( ALOG10(1.0 + zmax) / objdloglam )


;;; CHECK IF A GUESS REDSHIFT zguess WAS PASSED ALONG WITH A PIXEL
;;; WINDOW pwidth. IF SO, THEN RESET pmin AND pmax ACCORDING TO THE
;;; GUESS VALUE AND THE WINDOW.
  IF N_ELEMENTS(zguess) GT 0 AND KEYWORD_SET(pwidth) THEN BEGIN
    IF KEYWORD_SET(width) THEN width1 = width $
    ELSE width1 = pwidth
      pmin = FLOOR( ALOG10(1.0 + zguess) / objdloglam - 0.5*(pwidth+1+width1))
      pmax = FLOOR( ALOG10(1.0 + zguess) / objdloglam + 0.5*(pwidth+1+width1))
  ENDIF

; if pmax is too large, reset it!
; old version - parentheses error?
;  maxp = fix(objloglam0 + objdloglam*.99*n_elements(objflux) - starloglam0/ $
;          objdloglam)  


;ALTERATION BY BJW, 8/21/03 
 maxp =long((objloglam0 + objdloglam*.99*n_elements(objflux) - starloglam0)/$
            objdloglam)  

  if maxp lt pmax then begin ;limit upper redshift range to have overlap
     pmax = maxp
     print, 'resetting pmax to: ', maxp 
  endif


  IF abs(objdloglam - stardloglam) GT 0.05*objdloglam THEN $
    MESSAGE, 'Template and object lambda resolution do NOT match!'

   ;----------
   ; Compute the redshift difference between the first pixel of the object
   ; spectra and the template.
  poffset = (objloglam0 - starloglam0) / objdloglam

;  print, 'poffset, pmin,pmax :', poffset, pmin, pmax
;;; CALL zcompute.pro TO COMPUTE THE REDSHIFT(S).
   zans = zcompute(objflux, objivar, specref, poffset=poffset, $
                  pmin=pmin, pmax=pmax, nfind=nfind, width=width, $
                  plottitle=plottitle, _EXTRA=EXTRA)
   ;----------
   ; Convert redshift (and error) from pixels to the conventional dimensionless
   ; value.  Do not modify any errors that are less than zero, since those
   ; can be used as just warning flags from the fit.

   indx = where(zans.dof GT 0, npeak)
   if (npeak GT 0) then $
    zans[indx].z = 10.^(objdloglam * zans[indx].z) - 1.

   jndx = where(zans.dof GT 0 and zans.z_err GE 0)
   if (jndx[0] NE -1) then $
    zans[jndx].z_err = $
     alog(10d) * objdloglam * zans[jndx].z_err * (1 + zans[jndx].z)

   zerr = zans[indx].z_err
   return, zans[indx].z
end
;------------------------------------------------------------------------------


function continuum, lambda, spec, error, name=name
    method = 'synth'
    
    ivar = error^(-2)
    w = where(finite(ivar) and finite(spec) and ivar gt 0 and error gt 0, c)

    case method of
        'spline': begin
            bkpt = slatec_splinefit(lambda[w], spec[w], coeff, invvar=ivar[w], everyn=500 < round(c/4), upper=2, lower=2, /silent)
            cont = slatec_bvalu(lambda, bkpt, coeff)
        end
        'poly': begin
            pold = dblarr(4)
            p = pold + 1
            iter = 0
            while ~array_equal(p, pold) do begin
                pold = p
                p = poly_fit(lambda[w], spec[w], 3, measure_errors=error[w])
                cont = poly(lambda, reform(p))
                sigma = stddev(spec/cont)
                w = where(finite(ivar) and finite(spec) and ivar gt 0 and error gt 0 and (spec/cont - 1.0) lt -5.0*sigma or (spec/cont - 1.0) gt 5.0*sigma, c)
                iter++
                if c lt 100 or iter ge 5 then break
            endwhile
        end
        'synth': begin
            allframes = mrdfits('M15_M92_allframes.fits', 1, /silent)
            allframes = allframes[where(strtrim(allframes.name, 2) eq strtrim(name, 2))]
            
            linelist = mrdfits('cont/hires.list.fits.gz', 1, /silent)
            linelist = linelist[where(linelist.lambda ge min(lambda) and linelist.lambda le max(lambda), nlist)]

            openw, lun, 'synthcont/hires.list', /get_lun
            printf, lun, 'HIRESLIST'
            for i=0,nlist-1 do printf, lun, linelist[i].lambda, linelist[i].species, linelist[i].ep, linelist[i].loggf, format='(D10.3,D10.1,D10.2,D10.3)'
            close, lun
            free_lun, lun

            make_par, parfile='synthcont/hires.par', linefile='synthcont/hires.list', atmfile='synthcont/hires.atm', outfile='synthcont/hires.out2', driver='synth', minlambda=min(lambda), maxlambda=max(lambda), stronglist='cont/hires.strong', c12c13=10.0
            cfe = 0.2 + ((0.2 + 1.0)/(2.4 - 0.7))*((allframes.loggphot < 2.4) - 2.4)
            interp_atm, allframes.teff_mb20, allframes.loggphot, 2.13 - 0.23*allframes.loggphot, -2.4, 0.3, outfile='synthcont/hires.atm', tweakel=[6], tweakabund=[cfe]
            spawn, 'MOOGSILENT synthcont/hires.par'

            moog = read_moog_spec('synthcont/hires.out2', /newmoog, /singleabund)
            mooglambda = moog.lambda * (1d + (allframes.vr)/2.99792d5)
            airtovac, mooglambda
            moogspec = moog.spec
            moogspec = smooth_gauss_wrapper(mooglambda, moogspec, lambda, lambda / 30000. / 2.35)

            w = where(finite(ivar) and finite(spec) and ivar gt 0 and error gt 0 and finite(moogspec) and moogspec gt 0.01, c)
            if c lt 20 then cont = 1 else begin
                bkpt = slatec_splinefit(lambda[w], spec[w] / moogspec[w], coeff, invvar=ivar[w] * moogspec[w]^2., everyn=500 < round(c/2), upper=1.5, lower=1.5, /silent)
                cont = slatec_bvalu(lambda, bkpt, coeff)
            endelse

            if 0 then begin
                splot, lambda, spec
                soplot, lambda, moogspec*median(spec)/median(moogspec), color=fsc_color('red')
                soplot, lambda, spec/moogspec, color=fsc_color('blue')
                soplot, lambda, cont, color=fsc_color('orange')
                ;stop
            endif
                        
        end
    endcase
    
    if 0 and min(lambda) lt 5160 and max(lambda) gt 5160 then begin
        splot, lambda, spec
        soplot, lambda, cont, color=fsc_color('red')
        stop
    endif
    
    return, cont
end


pro stack_all        
    allframes = mrdfits('M15_M92_allframes.fits', 1, /silent)
    n = n_elements(allframes)

    for i=0,n-1 do begin
        file_exist = 0
        if strtrim(allframes[i].name, 2) ne 'M92-star-5' and strtrim(allframes[i].name, 2) ne 'M92-star-7' then continue
        for j=1,3 do begin
            fitsfile = strtrim(allframes[i].makeefile, 2)+'_'+string(j, format='(I1)')+'.fits'
            if ~file_test(fitsfile) then begin
                print, fitsfile+' not found'
                continue
            endif
            file_exist = 1
            efitsfile = strtrim(allframes[i].makeefile, 2)+'_'+string(j, format='(I1)')+'e.fits'
            print, fitsfile
            makeej = mrdfits(fitsfile, 0, hdr, /silent)
            errorj = mrdfits(efitsfile, 0, hdr, /silent)
            makee = j eq 1 ? makeej : [[makee], [makeej]]
            error = j eq 1 ? errorj : [[error], [errorj]]
            lambda = j eq 1 ? hires_lambda(hdr) : [[lambda], [hires_lambda(hdr)]]
        endfor
        if ~file_exist then continue

        npix = (size(lambda))[1]
        norders = (size(lambda))[2]
        lambdarange = minmax(lambda)
        b = (alog10(lambda[1]) - alog10(lambda[0]))
        a = alog10(lambdarange[0])
        nlambda = ceil((alog10(lambdarange[1]) - a) / b)

        hires = {name:' ', date:' ', exptime:0d, mjd:0d, xdisp:' ', decker:' ', slitwid:0d, xdangle:0d, echangle:0d, lambda:dblarr(nlambda), spec:dblarr(nlambda), ivar:dblarr(nlambda)}
        hires.lambda = 10.^(a + b*dindgen(nlambda))
        hires.name = allframes[i].name
        hires.date = sxpar(hdr, 'DATE')
        hires.exptime = sxpar(hdr, 'EXPTIME')
        hires.mjd = sxpar(hdr, 'MJD')
        hires.xdisp = sxpar(hdr, 'XDISPERS')
        hires.decker = sxpar(hdr, 'DECKNAME')
        hires.slitwid = sxpar(hdr, 'SLITWID')
        hires.xdangle = sxpar(hdr, 'XDANGL')
        hires.echangle = sxpar(hdr, 'ECHANGL')
        
        for j=0,norders-1 do begin
            wave = lambda[*,j]
            flux = makee[*,j]
            fluxerror = error[*,j]
            cont = continuum(wave, flux, fluxerror, name=hires.name)
            flux /= cont
            fluxerror /= cont

            if n_elements(flux) le 1 then message, 'WTF?'
            
            wj = where(hires.lambda ge min(wave) and hires.lambda le max(wave))
            x_specrebin, wave, flux, hires.lambda[wj], specj, var=(fluxerror)^2d, nwvar=nwvar, /silent
            x_specrebin, wave, dblarr(n_elements(flux))+1d , hires.lambda[wj], norm, var=(fluxerror)^2d, /silent
            specj /= norm
            ivarj = norm^2d / nwvar

            w0 = where(hires.lambda[wj] lt min(wave) or hires.lambda[wj] gt max(wave) or ~finite(specj) or ~finite(ivarj) or ivarj lt 0d)
            ivarj[w0] = 0d
            hires.spec[wj] += specj*ivarj
            hires.ivar[wj] += ivarj
        endfor
        w = where(hires.ivar gt 0)
        hires.spec[w] /= hires.ivar[w]

        lambda_arr = hires.lambda
        vactoair, lambda_arr
        hires.lambda = lambda_arr
        
        ;splot, hires.lambda, hires.spec, xrange=[4600, 4800]
        ;soplot, [0, 1d5], [1, 1], color=fsc_color('orange')
        ;stop
        
        specfile = strtrim(allframes[i].specfile, 2)
        mwrfits, hires, specfile, /create
        spawn, 'gzip -f '+specfile
        allframes[i].specfile = specfile+'.gz'
    endfor
    ;mwrfits, allframes, 'M15_M92_allframes.fits', /create
end


pro measure_vr
    clight = 2.99792458d5
    allframes = mrdfits('M15_M92_allframes.fits', 1, /silent)
    n = n_elements(allframes)

    jorder = 16
    wvmin = 4000
    wvmax = 5000
    
    newstr = {vr:0d, vrerr:0d}
    allframes = struct_trimtags(allframes, except_tags=['vr', 'vrerr'])
    allframes = struct_addtags(allframes, replicate(newstr, n))

    vr_ref_star = 'M92-star-1'
    vr_ref = -39.76
    vrerr_ref = 0.0
    wref = where(strtrim(allframes.name, 2) eq vr_ref_star)
    specfile = getenv('chome')+'keck/hires/M15_M92/makee/'+vr_ref_star+'_1.fits'
    especfile = getenv('chome')+'keck/hires/M15_M92/makee/'+vr_ref_star+'_1e.fits'
    specref = mrdfits(specfile, 0, hdr, /silent)
    errorref = mrdfits(especfile, 0, /silent)
    waveref = hires_lambda(hdr)

    a = alog10(waveref[0,jorder])
    b = alog10(waveref[1,jorder])-a
    nlambda = ceil((alog10(max(waveref[*,jorder])) - a) / b)
    logwave = 10^(a + b*dindgen(nlambda))

    x_specrebin, waveref[*,jorder], specref[*,jorder], logwave, specref, var=(errorref[*,jorder])^2, nwvar=varref, /silent
    ivarref = 1. / varref

    for i=0,n-1 do begin
        objname = strtrim(allframes[i].name, 2)

        iorder = jorder
        if objname eq 'M15-star-4' then iorder = 11
        if objname eq 'M92-star-13' then iorder = 18
        
        specfile = getenv('chome')+'keck/hires/M15_M92/makee/'+objname+'_1.fits'
        especfile = getenv('chome')+'keck/hires/M15_M92/makee/'+objname+'_1e.fits'
        if ~file_test(specfile) then continue
        spec = mrdfits(specfile, 0, hdr, /silent)
        error = mrdfits(especfile, 0, /silent)
        wave = hires_lambda(hdr)

        x_specrebin, wave[*,iorder], spec[*,iorder], logwave, spec, var=(error[*,iorder])^2, nwvar=var, /silent
        ivar = 1. / var
        wr = where(logwave ge wvmin and logwave le wvmax and ivarref gt 0 and ivar gt 0 and finite(ivarref) and finite(ivar), c)
        if c lt 20 then begin
            print, objname
            stop
            continue
        endif

        if objname ne vr_ref_star then begin
            if c lt 20 then continue

            z = zfind(logwave[wr] / (1d + vr_ref/clight), specref[wr], ivarref[wr], logwave[wr], spec[wr], ivar[wr], zmin=-100./clight, zmax=100./clight, zerr=zerr, wvmin=wvmin, wvmax=wvmax)
            vr = clight * z
            vrerr = sqrt((clight * zerr)^2 + vrerr_ref^2)

            splot, logwave / (1d + vr_ref/clight), specref / median(specref[wr]), xrange=[4440, 4447];, xrange=[wvmin, wvmax], yrange=[0, 1.25]
            soplot, logwave / (1d + vr/clight), spec / median(spec[wr]), color=fsc_color('red')
        endif else begin
            vr = vr_ref
            vrerr = vrerr_ref
        endelse

        print, allframes[i].name, vr, vrerr
        ;stop
        allframes[i].vr = vr
        allframes[i].vrerr = vrerr
    endfor
    mwrfits, allframes, 'M15_M92_allframes.fits', /create
end


pro filter_allframes
    allframesfile = 'M15_M92_allframes.fits'
    allframes = mrdfits(allframesfile, 1, /silent)
    allframes.specfile = strtrim(allframes.specfile, 2)+'.gz'
    w = where(file_test(strtrim(allframes.specfile, 2)))
    allframes = allframes[w]
    stop
    mwrfits, allframes, allframesfile, /create
end
