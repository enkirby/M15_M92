function voigt_ul, pars, df
    common voigt_ul, lambda, resid, residerr, chi2

    voigtx = voigtwave(lambda, pars)
    chisq = total(((resid - voigtx) / residerr)^2)
    return, abs(chisq - chi2)
end


; ================= EVENTS =================
pro hiresspec3_event, ev
    widget_control, ev.top, get_uvalue=obj
    widget_control, ev.id, get_uvalue=value
    if (n_elements(value) eq 0) then value = ''
    name = strmid(tag_names(ev, /structure_name), 7, 4)
    
    case (name) of
        'BUTT': obj->handle_button, ev
        'TEXT': obj->handle_text, ev
        'COMB': obj->handle_combobox, ev
        'DRAW': begin
            if ev.type eq 0 then obj->handle_draw_click, ev
            if ev.type eq 5 then obj->handle_draw_key, ev
        end
        'DONE': widget_control, ev.top, /destroy
        else: obj->redraw
    endcase
end


; ================ BUTTONS ================
function line_event, ev
    self->redraw
end


pro hiresspec3::handle_button, ev
    widget_control, ev.top, get_uvalue=obj
    widget_control, ev.id, get_uvalue=uvalue
    
    case (uvalue) of
        'back': self->newspec, increment=-1
        'next': self->newspec, increment=1
        'default_lambda': self->default_lambda
        'measure_all': begin
            if self.change eq 1 then begin
                response = dialog_message('Are you sure you want to abandon your changes to the linelist?', /cancel, /default_cancel, /center, dialog_parent=self.base, title='Overwrite Linelist?')
                if response ne 'OK' then return
            endif
            self->measure_all
        end
        'blue': self->lambdarange, /blue
        'red': self->lambdarange, /red
        'eps': self->make_eps
        'export': self->export_linelist
        'reload': begin
            if self.change eq 1 then begin
                response = dialog_message('Are you sure you want to abandon your changes to the linelist?', /cancel, /default_cancel, /center, dialog_parent=self.base, title='Overwrite Linelist?')
                if response ne 'OK' then return
            endif
            self->getlinelist, /reload
            self->redraw
        end
        'exit': begin
            self->export_linelist
            widget_control, ev.top, /destroy
        end        
    endcase
end


pro hiresspec3::newspec, increment=increment
    n = n_elements(*self.names)
    newi = self.i + increment
    if newi lt 0 then begin
        widget_control, widget_info(self.base, find_by_uname='status'), set_value='This is the first spectrum.'
        return
    endif
    if newi gt n-1 then begin
        widget_control, widget_info(self.base, find_by_uname='status'), set_value='This is the last spectrum.'
        return
    endif
    self->export_linelist
    widget_control, widget_info(self.base, find_by_uname='objlist'), set_combobox_select=newi
    self.i = newi
    self->getscience
    self->redraw
end


pro hiresspec3::export_linelist
    name = (*self.names)[self.i]
    ewfile = self.outdir+name+'_Ji20_ew.fits'
    fname = self.outdir+name+'.ew'
    if self.change eq 0 and file_test(ewfile) and file_test(fname) then begin
        widget_control, widget_info(self.base, find_by_uname='status'), set_value='You have not made any changes to the line list.'
        return
    endif

    linelist = *self.linelist
    linelist = linelist[sort(linelist.lambda)]
    linelist = linelist[sort(linelist.species)]

    mwrfits, linelist, ewfile, /create

    w = where(linelist.ew gt 0.0, c)
    openw, lun, fname, width=10000, /get_lun
    printf, lun, name
    for i=0,c-1 do begin
        case 1 of
            linelist[w[i]].upperlimit eq 1: printf, lun, linelist[w[i]].lambda, linelist[w[i]].species, linelist[w[i]].ep, linelist[w[i]].loggf, linelist[w[i]].ew, 'UPPERLIMIT', format='(2X,D8.3,6X,D4.1,1X,D9.2,1X,D9.3,20X,D10.1,2X,A10)'
            linelist[w[i]].manual eq 2: printf, lun, linelist[w[i]].lambda, linelist[w[i]].species, linelist[w[i]].ep, linelist[w[i]].loggf, linelist[w[i]].ew, 'LOCALCONT ', format='(2X,D8.3,6X,D4.1,1X,D9.2,1X,D9.3,20X,D10.1,2X,A10)'
            else: printf, lun, linelist[w[i]].lambda, linelist[w[i]].species, linelist[w[i]].ep, linelist[w[i]].loggf, linelist[w[i]].ew, 'HIRESSPEC3', format='(2X,D8.3,6X,D4.1,1X,D9.2,1X,D9.3,20X,D10.1,2X,A10)'
        endcase
    endfor
    free_lun, lun
    close, lun
    self.change = 0
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Changes saved.'
end


pro hiresspec3::default_lambda, noredraw=noredraw
    self.ylim = [0, 1.3]
    self.lambdalim = [4500, 4600]
    widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(self.ylim[0], format='(D5.2)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(self.ylim[1], format='(D5.2)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0], format='(D7.1)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1], format='(D7.1)'), /rem)
    if ~keyword_set(noredraw) then self->redraw
end


pro hiresspec3::measure_all
    linelist = *self.linelist
    science = *self.science
    w = where(linelist.lambda gt min(science.lambda) and linelist.lambda lt max(science.lambda), c)
    for i=0,c-1 do begin
        coords = [linelist[w[i]].lambda, 1.0]
        self->measure_ew, coords, /noredraw, /all
    endfor
    self->redraw
end


pro hiresspec3::lambdarange, red=red, blue=blue
    if keyword_set(red)+keyword_set(blue) ne 1 then message, 'You must specify red or blue.'
    lrange = self.lambdalim[1] - self.lambdalim[0]
    if keyword_set(blue) then begin
        if self.lambdalim[0] lt min((*self.science).lambda) then begin
            widget_control, widget_info(self.base, find_by_uname='status'), set_value='The is bluest part of the spectrum.'
            return
        endif
        llownew = self.lambdalim[0] - 0.6*lrange
        lhighnew = self.lambdalim[1] - 0.6*lrange
    endif
    if keyword_set(red) then begin
        if self.lambdalim[1] gt max((*self.science).lambda) then begin
            widget_control, widget_info(self.base, find_by_uname='status'), set_value='The is reddest part of the spectrum.'
            return
        endif
        llownew = self.lambdalim[0] + 0.6*lrange
        lhighnew = self.lambdalim[1] + 0.6*lrange
    endif
    widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(llownew, format='(D7.1)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(lhighnew, format='(D7.1)'), /rem)
    self->lambdalow, llownew, /noredraw
    self->lambdahigh, lhighnew
end


pro hiresspec3::measure_ew, coords, limits=limits, noredraw=noredraw, upperlimit=upperlimit, localcont=localcont, setcont=setcont, all=all
    common voigt_ul, lambda, resid, residerr, chi2

    widget_control, widget_info(self.base, find_by_uname='voigt'), get_value=voigtfit

    linelist = *self.linelist
    science = *self.science

    if ~keyword_set(all) then widget_control, widget_info(self.base, find_by_uname='status'), set_value='Matching cursor position to line list ...'        
    lambda = coords[0]

    if ~keyword_set(all) then begin
        wactive = where(linelist.lambda gt self.lambdalim[0] and linelist.lambda lt self.lambdalim[1], cw)
        if cw eq 0 then begin
            widget_control, widget_info(self.base, find_by_uname='status'), set_value='This window contains no lines in the linelist.'
            return
        endif
        junk = min(abs(linelist[wactive].lambda-lambda), wpress)
        wap = wactive[wpress]
    endif else begin
        junk = min(abs(linelist.lambda-lambda), wap)
    endelse
        
    w = where(linelist.doppwidth gt 0.0 and linelist.doppwidth lt 0.2 and linelist.ew gt 0.0 and linelist.ew lt 1d4, cw)
    if cw gt 2 then doppwidth = median(linelist[w].doppwidth/linelist[w].lambda)*linelist[wap].lambda else doppwidth = 0.05
    if doppwidth gt 0.2 or doppwidth lt 0.01 then doppwidth = 0.05

    wr1 = where(science.lambda gt linelist[wap].lambda-(3.0*doppwidth) and science.lambda lt linelist[wap].lambda+(3.0*doppwidth) and science.ivar gt 0.0 and finite(science.ivar), cr1)
    if cr1 lt 3 then begin
        widget_control, widget_info(self.base, find_by_uname='status'), set_value='No spectrum in range.'
        linelist[wap].ew = -999.0
        linelist[wap].ewerr = -999.0
        linelist[wap].manual = 0
        linelist[wap].upperlimit = 0
        linelist[wap].doppwidth = 0.0
        linelist[wap].doppwidtherr = 0.0
        linelist[wap].rate = 0.0
        linelist[wap].rateerr = 0.0
        goto, doneew
    endif
    if ~keyword_set(limits) then begin
        wrange = where(science.lambda gt linelist[wap].lambda-(1.5*doppwidth) and science.lambda lt linelist[wap].lambda+(1.5*doppwidth) and science.ivar gt 0.0 and finite(science.ivar), cl)
    endif else begin
        wrange = where(science.lambda ge min(limits) and science.lambda le max(limits) and science.ivar gt 0.0 and finite(science.ivar), cl)
    endelse
    lambda = science.lambda[wrange]-linelist[wap].lambda

    s = sort(lambda)
    lambda = lambda[s]
    spec = science.spec[wrange[s]]
    ivar = science.ivar[wrange[s]]
    ;continuum = science.continuum[wrange]

    ;w = where(ivar ge 0 and ivar lt 8.0*stddev(ivar))
    ;lambda = lambda[w]
    ;spec = spec[w]
    ;ivar = ivar[w]
    resid = 1.0 - spec
    
    ewguess = (int_tabulated(lambda, resid)*1d3) > 5.0

    seed = 68484L

    wvalid = where(linelist.ew gt 2 and linelist.upperlimit eq 0, c)
    if c gt 10 then begin
        doppwidthguess = mean(linelist[w].doppwidth)
        rateguess = 0.0 ;mean(linelist[w].rate)
    endif else begin
        doppwidthguess = 0.05
        rateguess = 0.01
    endelse
    if ~voigtfit then rateguess = 0.0
    
    pi = replicate({value:0.0, fixed:0, limited:[1,1], limits:[0.0, 0.0], parname:'EW', mpprint:0, mpformat:'(D6.1)', step:10d, tied:''}, 5)
    pi.parname = ['EW', 'doppwidth', 'rate', 'dlambda', 'localcont']
    pi.value = [ewguess, doppwidthguess, rateguess, 0.0, keyword_set(setcont) ? (1.0 - coords[1]) : 0.0]
    pi.step = [10.0, 0.01, 0.01, 0.05, 0.05]
    pi[0].limited = [1, 0]
    pi[0].limits = [0.0, 10000.0]
    pi[1].limited = [1, 0]
    pi[1].limits = [0.0, 1.0]
    pi[2].limited = [1, 0]
    pi[2].limits = [0.0, 1.0]
    pi[3].limited = [0, 0]
    pi[3].limits = [-100, 100]
    pi[4].limited = [0, 0]
    pi[4].limits = [-1, 1]

    case 1 of
        keyword_set(upperlimit): begin
            if ~keyword_set(all) then widget_control, widget_info(self.base, find_by_uname='status'), set_value='Auto-measuring EW upper limit for '+linelist[wap].name+' ...'
            pi.fixed = [0, 1, 1, 1, ~keyword_set(localcont)]

            residerr = ivar^(-0.5)
            parameterarray = mpfitfun('voigtwave', lambda, resid, residerr, parinfo=pi, /quiet, perror=perror)
            chi2 = 0.0
            chisq0 = voigt_ul(parameterarray, 0.0)
            chi2 = chisq0 + 9.0
            pi[0].limits[0] = parameterarray[0]
            pi[0].value = (parameterarray[0] + pi[0].step)
            parameterarray = tnmin('voigt_ul', parinfo=pi, /quiet, bestmin=bestmin, /autoderivative)
            chisqul3 = bestmin + chi2

            if parameterarray[0] le 0.0 or ~finite(parameterarray[0]) then begin
                linelist[wap].ew = -999.0
                linelist[wap].ewerr = -999.0
                linelist[wap].manual = 0
                linelist[wap].upperlimit = 0
                linelist[wap].upperlimit = 0
                goto, doneew
            endif
            ew = parameterarray[0]
            doppwidth = parameterarray[1]
            rate = parameterarray[2]
            centralw = parameterarray[3]
            dcont = parameterarray[4]
            if ~finite(ew) or doppwidth gt 1.0 then begin
                linelist[wap].ew = -999.0
                linelist[wap].ewerr = -999.0
                linelist[wap].manual = 0
                linelist[wap].upperlimit = 0
                linelist[wap].dlambda = 0
                goto, doneew
            endif
            linelist[wap].upperlimit = 1
            linelist[wap].ew = ew
            linelist[wap].ewerr = ew / 3.
            linelist[wap].doppwidth = doppwidth
            linelist[wap].rate = rate
            linelist[wap].manual = 1
            linelist[wap].dlambda = centralw
            linelist[wap].localcont = dcont
        end
        keyword_set(localcont): begin
            if ~keyword_set(all) then widget_control, widget_info(self.base, find_by_uname='status'), set_value='Auto-measuring EW for '+linelist[wap].name+' ...'
            
            pi.fixed = [0, 0, ~voigtfit, 0, 0]

            array = mpfitfun('voigtwave', lambda, resid, ivar^(-0.5), parinfo=pi, /quiet, perror=perror)
            if (size(perror))[0] eq 0 then begin
                linelist[wap].ew = -999.0
                linelist[wap].ewerr = -999.0
                linelist[wap].manual = 0
                linelist[wap].upperlimit = 0
                linelist[wap].upperlimit = 0
                goto, doneew
            endif
            ew = array[0]
            doppwidth = array[1]
            rate = array[2]
            centralw = array[3]
            dcont = array[4]
            ewerr = perror[0]
            doppwidtherr = perror[1]
            rateerr = perror[2]
            dconterr = perror[4]
            if ~finite(ew) or ew lt 2.*ewerr or doppwidth gt 1.0 or doppwidth lt 0.0001 then begin
                linelist[wap].ew = -999.0
                linelist[wap].ewerr = -999.0
                linelist[wap].manual = 0
                linelist[wap].upperlimit = 0
                linelist[wap].upperlimit = 0
                goto, doneew
            endif

            linelist[wap].upperlimit = 0
            linelist[wap].ew = ew
            linelist[wap].ewerr = ewerr
            linelist[wap].doppwidth = doppwidth
            linelist[wap].doppwidtherr = doppwidtherr
            linelist[wap].rate = rate
            linelist[wap].rateerr = rateerr
            linelist[wap].manual = 2
            linelist[wap].dlambda = centralw
            linelist[wap].localcont = dcont
            linelist[wap].localconterr = dconterr
        end
        else: begin
            if ~keyword_set(all) then widget_control, widget_info(self.base, find_by_uname='status'), set_value='Auto-measuring EW for '+linelist[wap].name+' ...'

            pi.fixed = [0, 0, ~voigtfit, 0, 1]

            perror = replicate(-999d, 5)
            a = mpfitfun('voigtwave', lambda, resid, ivar^(-0.5), parinfo=pi, /quiet, perror=perror)
            if (size(perror))[0] eq 0 then begin
                linelist[wap].ew = -999.0
                linelist[wap].ewerr = -999.0
                linelist[wap].manual = 0
                linelist[wap].upperlimit = 0
                linelist[wap].upperlimit = 0
                goto, doneew
            endif
            ew = a[0]
            doppwidth = a[1]
            rate = a[2]
            centralw = a[3]
            ewerr = perror[0]
            doppwidtherr = perror[1]
            rateerr = perror[2]
            if ~finite(ew) or ew lt 2.*ewerr or doppwidth gt 1.0 or doppwidth lt 0.0001 then begin
                linelist[wap].ew = -999.0
                linelist[wap].ewerr = -999.0
                linelist[wap].manual = 0
                linelist[wap].upperlimit = 0
                linelist[wap].upperlimit = 0
                goto, doneew
            endif
            
            linelist[wap].upperlimit = 0
            linelist[wap].ew = ew[0]
            linelist[wap].ewerr = ewerr
            linelist[wap].doppwidth = doppwidth
            linelist[wap].doppwidtherr = doppwidtherr
            linelist[wap].rate = rate
            linelist[wap].rateerr = rateerr
            linelist[wap].manual = 1
            linelist[wap].dlambda = centralw
            linelist[wap].localcont = keyword_set(setcont) ? (1.0 - coords[1]) : 0.0
        end
    endcase
    self.change = 1

    doneew:
    ptr_free, self.linelist
    self.linelist = ptr_new(linelist)
    if ~keyword_set(noredraw) then self->redraw
    self->line_info, coords, all=all
end

pro hiresspec3::measure_limits, coords, localcont=localcont, upperlimit=upperlimit, setcont=setcont
    case 1 of
        self.measurelim[0] eq -1: begin
            case 1 of
                ~keyword_set(localcont) and ~keyword_set(setcont) and ~keyword_set(upperlimit): keypress = 'x'
                keyword_set(localcont) and ~keyword_set(upperlimit): keypress = 'c'
                ~keyword_set(localcont) and ~keyword_set(setcont) and keyword_set(upperlimit): keypress = 'w'
                keyword_set(localcont) and ~keyword_set(setcont) and keyword_set(upperlimit): keypress = 'y'
                keyword_set(setcont) and ~keyword_set(upperlimit): keypress = 's'
                keyword_set(setcont) and keyword_set(upperlimit): keypress = 'v'
            endcase
            self.measurelim[0] = coords[0]
            widget_control, widget_info(self.base, find_by_uname='status'), set_value='Press '+keypress+' on the other side of the line.'
        end
        self.measurelim[1] eq -1: begin
            self.measurelim[1] = coords[0]
            self->measure_ew, [(self.measurelim[0]+self.measurelim[1])/2., coords[1]], limits=self.measurelim, localcont=localcont, upperlimit=upperlimit, setcont=setcont
            self.measurelim = [-1d, -1d]
        end
        else: begin
            self.measurelim[0] = -1
            self.measurelim[1] = -1
            widget_control, widget_info(self.base, find_by_uname='status'), set_value='Something went wrong.  Line measurement cancelled.'
        end
    endcase
end

; ============== TEXT BOXES =============
pro hiresspec3::handle_text, ev
    widget_control, ev.id, get_uvalue=uvalue
    widget_control, ev.top, get_uvalue=obj
    widget_control, ev.id, get_value=val

    case (uvalue) of
        'ylow': self->ylow, val
        'yhigh': self->yhigh, val
        'lambdalow': self->lambdalow, val
        'lambdahigh': self->lambdahigh, val
        else: 
    end
end


pro hiresspec3::lambdalow, lambdalow, noredraw=noredraw
    if lambdalow gt max((*self.science).lambda) then begin
        lambdalow = max((*self.science).lambda) - (self.lambdalim[1] - lambdalow)
        widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(lambdalow, format='(D7.1)'), /rem)
    endif
    self.lambdalim[0] = lambdalow
    if ~keyword_set(noredraw) then self->redraw
end


pro hiresspec3::lambdahigh, lambdahigh, noredraw=noredraw
    if lambdahigh lt min((*self.science).lambda) then begin
        lambdahigh = min((*self.science).lambda) + (lambdahigh - self.lambdalim[0])
        widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(lambdahigh, format='(D7.1)'), /rem)
    endif
    self.lambdalim[1] = lambdahigh
    if ~keyword_set(noredraw) then self->redraw
end


pro hiresspec3::ylow, ylow, noredraw=noredraw
    self.ylim[0] = ylow
    if ~keyword_set(noredraw) then self->redraw
end


pro hiresspec3::yhigh, yhigh
    self.ylim[1] = yhigh
    self->redraw
end


; ============== COMBOBOX  =============
pro hiresspec3::handle_combobox, ev
    self->export_linelist
    self.i = ev.index
    self->getscience
    self->redraw
end


; ============= DRAW CLICK =============
pro hiresspec3::handle_draw_click, ev
    click_coords = convert_coord(ev.x, ev.y, /device, /to_data)
    case ev.modifiers of
        0: begin
            case ev.press of
                1: lrange = (self.lambdalim[1] - self.lambdalim[0]) / 2.
                4: lrange = (self.lambdalim[1] - self.lambdalim[0]) * 2.
                else: begin
                    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Mouse click not recognized.'
                    return
                end
            endcase
            llownew = click_coords[0] - lrange/2.
            lhighnew = click_coords[0] + lrange/2.
            widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(llownew, format='(D7.1)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(lhighnew, format='(D7.1)'), /rem)
            self->lambdalow, llownew, /noredraw
            self->lambdahigh, lhighnew
        end
        1: begin
            if ev.press ne 1 then begin
                widget_control, widget_info(self.base, find_by_uname='status'), set_value='Mouse click not recognized.'
                return
            endif
            lrange = self.lambdalim[1] - self.lambdalim[0]
            llownew = click_coords[0] - lrange/2.
            lhighnew = click_coords[0] + lrange/2.        
            widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(llownew, format='(D7.1)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(lhighnew, format='(D7.1)'), /rem)
            self->lambdalow, llownew, /noredraw
            self->lambdahigh, lhighnew
        end
        2: begin
            widget_control, widget_info(self.base, find_by_uname='status'), set_value='Mouse click not recognized.'
            return
        end
        3: begin
            case ev.press of
                1: yrange = (self.ylim[1] - self.ylim[0]) / 2.
                4: yrange = (self.ylim[1] - self.ylim[0]) * 2.
                else: begin
                    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Mouse click not recognized.'
                    return
                end
            endcase
            ylownew = click_coords[1] - yrange/2.
            yhighnew = click_coords[1] + yrange/2.
            widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(ylownew, format='(D5.2)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(yhighnew, format='(D5.1)'), /rem)
            self->ylow, ylownew, /noredraw
            self->yhigh, yhighnew
        end
        else: begin
            widget_control, widget_info(self.base, find_by_uname='status'), set_value='Mouse click not recognized.'
            return            
        end
    endcase
end


; ============= DRAW KEYS ==============
pro hiresspec3::handle_draw_key, ev
    if ev.press ne 1 then return
    key = string(ev.ch)
    coords = convert_coord(ev.x, ev.y, /device, /to_data)
    case key of
        'e': self->measure_ew, coords
        'n': self->goto_line, coords, increment=1
        'b': self->goto_line, coords, increment=-1
        'p': self->goto_line, coords, increment=-1
        '1': self->goto_line, coords, /first
        '2': self->goto_line, coords, /last
        'd': self->delete_ew, coords
        'u': self->measure_ew, coords, /upperlimit
        'l': self->measure_ew, coords, /localcont
        'i': self->line_info, coords
        'x': self->measure_limits, coords
        'c': self->measure_limits, coords, /localcont
        'w': self->measure_limits, coords, /upperlimit
        'y': self->measure_limits, coords, /upperlimit, /localcont
        's': self->measure_limits, coords, /setcont
        'v': self->measure_limits, coords, /upperlimit, /setcont
        ' ': self->line_info, coords
        'g': self->goto_lambda
        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
    endcase
end


pro hiresspec3::goto_lambda
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Requesting new wavelength ...'
    lambda = textbox(title='Go to wavelength', group_leader=self.base, label='Wavelength: ', cancel=cancelled, xsize=400, Value='')
    if ~cancelled then begin
        lambda = double(lambda)
        lambdarange = self.lambdalim[1] - self.lambdalim[0]
        lambdalow = lambda - lambdarange/2.0
        lambdahigh = lambda + lambdarange/2.0
        widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(lambdalow, format='(D7.1)'), /rem)
        widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(lambdahigh, format='(D7.1)'), /rem)
        self->lambdalow, lambdalow, /noredraw
        self->lambdahigh, lambdahigh
    endif
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
end


pro hiresspec3::goto_line, coords, increment=increment, first=first, last=last
    lambda = coords[0]
    linelist = *self.linelist
    lines = linelist.lambda
    lines = lines[sort(lines)]
    junk = min(abs(lines-lambda), m)
    lambda = lines[m]
    lines = lines[where(lines ne lambda)]
    if keyword_set(first) then w = 0
    if keyword_set(last) then w = n_elements(lines)-1
    if ~keyword_set(first) and ~keyword_set(last) then begin
        wred = where(lines gt lambda, cred, complement=wblue, ncomplement=cblue)
        if increment gt cred or increment lt -1.0*cblue then begin
            widget_control, widget_info(self.base, find_by_uname='status'), set_value='No line found.'
            return
        endif
        w = increment gt 0 ? wred[increment-1] : wblue[cblue+increment]
    endif
    lambdarange = self.lambdalim[1] - self.lambdalim[0]
    lambdalow = lines[w] - lambdarange/2.0
    lambdahigh = lines[w] + lambdarange/2.0
    widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(lambdalow, format='(D7.1)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(lambdahigh, format='(D7.1)'), /rem)
    self->lambdalow, lambdalow, /noredraw
    self->lambdahigh, lambdahigh
    self->line_info, [lines[w], 1.0]
end


pro hiresspec3::delete_ew, coords
    lambda = coords[0]
    y = coords[1]

    linelist = *self.linelist
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Matching cursor position to measured EWs ...'
    dist = min(abs(lambda-linelist.lambda), wpress)
    if linelist[wpress].lambda gt self.lambdalim[1] or linelist[wpress].lambda lt self.lambdalim[0] then begin
        widget_control, widget_info(self.base, find_by_uname='status'), set_value='No EW measurement found in this window.'
        return
    endif
    linelist[wpress].ew = -999d
    linelist[wpress].ewerr = -999d
    linelist[wpress].manual = 0
    linelist[wpress].upperlimit = 0
    linelist[wpress].doppwidth = 0.0
    linelist[wpress].rate = 0.0
    ptr_free, self.linelist
    self.linelist = ptr_new(linelist)
    self.change = 1
    self->redraw
end


pro hiresspec3::line_info, coords, all=all
    lambda = coords[0]
    linelist = *self.linelist
    if ~keyword_set(all) then begin
        w = where(linelist.lambda gt self.lambdalim[0] and linelist.lambda lt self.lambdalim[1], cw)
        if cw eq 0 then begin
            widget_control, widget_info(self.base, find_by_uname='status'), set_value='No lines satisfy log(gf), excitation potential, and wavelength limits.'
            return
        endif        
    endif
    dist = min(abs(lambda-linelist.lambda), wpress)
    linelist = linelist[wpress]
    outtext = strmid(linelist.name, 0, strlen(linelist.name)-4)+strcompress(string(linelist.lambda, format='(D10.3)'), /rem)+'    EP = '+strcompress(string(linelist.ep, format='(D10.3)'), /rem)+' eV    log(gf) = '+strcompress(string(linelist.loggf, format='(D10.3)'), /rem)
    if linelist.ewerr gt 0 then outtext += '    EW '+(linelist.upperlimit eq 1 ? "<" : "=")+' '+strtrim(string(linelist.ew, format='(D10.2)'), 2)+(linelist.upperlimit eq 1 ? '' : ' +/- '+strtrim(string(linelist.ewerr, format='(D10.2)'), 2))+(linelist.doppwidth le 0 ? '' : '    Dopp = '+strtrim(string(abs(linelist.doppwidth), format='(D10.3)'), 2)+(linelist.doppwidtherr le 0 ? '' : ' +/- '+strtrim(string(linelist.doppwidtherr, format='(D10.3)'), 2)))+(linelist.rate le 0 ? '' : '    Rate = '+strtrim(string(abs(linelist.rate), format='(D10.3)'), 2)+(linelist.rateerr le 0 ? '' : ' +/- '+strtrim(string(linelist.rateerr, format='(D10.3)'), 2)))
    widget_control, widget_info(self.base, find_by_uname='status'), set_value=outtext    
end


; ============== DESTROY ===============
pro hiresspec3_cleanup, ev
    widget_control, ev, get_uvalue=obj
    obj_destroy, obj
end


pro hiresspec3::cleanup
    ptr_free, self.science
    ptr_free, self.linelist
    ptr_free, self.atomic
end


; ============== REDRAW ===============
pro hiresspec3::redraw, plotfile=plotfile
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Redrawing ...'
    science = *self.science
    linelist = *self.linelist
    atomic = *self.atomic
    
    widget_control, widget_info(self.base, find_by_uname='spec'), get_value=index
    if keyword_set(plotfile) then begin
        setplot
        angstrom = '!6!sA!r!u!9 %!6!n'
        device, filename=plotfile, xsize=9, ysize=6, /inches, /encapsulated, /color
        plot, [0, 1], [0, 1], xrange=self.lambdalim, yrange=self.ylim, /nodata, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xtitle='!6rest wavelength (!sA!r!u!9 %!6!n)!3', ytitle='!6flux / continuum!3'
    endif else begin
        wset, index
        plot, [0, 1], [0, 1], xrange=self.lambdalim, yrange=self.ylim, /nodata, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black')
    endelse

    colors = *self.colors
    w = where(linelist.lambda gt self.lambdalim[0] and linelist.lambda lt self.lambdalim[1], cw)

    if cw gt 0 then begin
        for i=0,n_elements(atomic)-1 do begin
            ww = where(linelist[w].atomic eq atomic[i], cww)
            if cww eq 0 then continue            
            usersym, [0, 0], [0, 6]
            oplot, [0.0, linelist[w[ww]].lambda], replicate(1.05, cww+1), color=fsc_color(colors[i]), psym=8
            xyouts, [0.0, linelist[w[ww]].lambda], replicate(1.13, cww+1), [' ', linelist[w[ww]].name], color=fsc_color(colors[i]), orientation=90, charsize=1.2
            
            for j=0,cww-1 do begin
                if linelist[w[ww[j]]].doppwidth ne 0 and linelist[w[ww[j]]].ewerr gt 0 then begin
                    lambdaew = (dindgen(501) - 250)/250.*4.*(linelist[w[ww[j]]].doppwidth > (0.1*linelist[w[ww[j]]].rate))*2.35
                    oplot, linelist[w[ww[j]]].lambda+lambdaew, 1d - voigtwave(lambdaew, [linelist[w[ww[j]]].ew, linelist[w[ww[j]]].doppwidth, linelist[w[ww[j]]].rate, linelist[w[ww[j]]].dlambda, linelist[w[ww[j]]].localcont]), color=fsc_color('red'), thick=3, linestyle=(linelist[w[ww[j]]].upperlimit eq 1 ? 1 : 0)
                endif
            endfor
        endfor
    endif

    w = where(science.lambda ge self.lambdalim[0] and science.lambda le self.lambdalim[1], c)
    ;if 1 then begin
    ;    wsr4215 = where(science.lambda gt 4215.2 and science.lambda lt 4215.85)
    ;    sr4215 = 1d3*int_tabulated(science.lambda[wsr4215], 1.0 - science.spec[wsr4215])
    ;    print, sr4215
    ;endif

    if c lt 2 then begin
        widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
        return
    endif
    lambda = science.lambda[w]
    spec = science.spec[w]
    ivar = science.ivar[w]
    ;continuum = science.continuum[w]
    science = {lambda:lambda, spec:spec, ivar:ivar}

    o2lambda = *self.o2lambda
    wo2 = where(o2lambda gt self.lambdalim[0] and o2lambda lt self.lambdalim[1], co2)
    if co2 gt 0 then begin
        for i=0,co2-1 do oplot, replicate(o2lambda[wo2[i]], 2), [0.5, 1.0], color=fsc_color('blue')
    endif
    
    vsym, 24
    oplot, minmax(science.lambda), [1.0, 1.0], color=fsc_color('orange')
    oplot, science.lambda, science.spec, color=fsc_color('black'), symsize=0.5;, psym=8
    ;oplot, science.lambda, 50*(1. / science.ivar), color=fsc_color('orange')
    
    if keyword_set(plotfile) then begin
        device, /close
        resetplot
    endif

    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
    widget_control, widget_info(self.base, find_by_uname='spec'), /input_focus
end


pro hiresspec3::make_eps
    plotfile = dialog_pickfile(default_extension='eps', path=getenv('CALTECH')+'hiresspec', /overwrite_prompt, dialog_parent=self.base, title='Plot File', /write)
    if plotfile ne '' then self->redraw, plotfile=plotfile
end


pro hiresspec3::getscience
    name = (*self.names)[self.i]
    specfile = (*self.specfiles)[self.i]
    vr = (*self.vr)[self.i]
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Reading spectrum for '+name+' ...'

    print, specfile
    data = mrdfits(specfile, 1, /silent)
    lambda = data.lambda / (1d + vr / 2.99792458d5)
    airtovac, lambda
    data.lambda = lambda
    
    w = where(data.lambda gt 5720 and data.lambda lt 5800, c)
    if c gt 10 then begin
        td = -1.0*ts_diff(data.lambda[w], 1)
        pixscale = median(td[where(td gt 0)])
        print, strtrim(pixscale, 2)+' pix per Ang'
        print, 'S/N = '+strtrim((1d / stddev(data.spec[w]))/sqrt(pixscale), 2)+' per Ang'
        print, minmax(data.lambda)
    endif

    science = {lambda:data.lambda, spec:data.spec, ivar:data.ivar, vr:vr, mjd:data.mjd}
    
    o2lambda = [7697.96, 7698.99]
    ;vactoair, o2lambda
    ;heliocor = helio_deimos(259.28079, 43.135944, 2000., jd=science.mjd+2400000.5)
    heliocor = helio_deimos(322.49304167, 12.167, 2000., jd=science.mjd+2400000.5)
    case name of
        'XI-80': extrav = 3.0
        'XI-19': extrav = 8.0
        'S3108': extrav = 1.0
        'C17333-0832': extrav = 11.5
        'IV-94': extrav = 9.0
        'S61': extrav = 12.5
        'D21': extrav = 12.5
        'S162': extrav = 1.0
        'XII-34': extrav = 1.0
        'XII-8': extrav = -9.0
        'V-45': extrav = -11.0
        else: extrav = 0.0
    endcase
    print, 'vcor: '+string(science.vr + heliocor + extrav)
    o2lambda *= (1d - ((science.vr + heliocor + extrav) / 2.99792458d5))
    self.o2lambda = ptr_new(o2lambda)

    wo2 = where(science.lambda gt o2lambda[0]-0.5 and science.lambda lt o2lambda[1]+0.5, co2)
    if co2 ge 10 and name ne 'AS2302' then begin
        if 1 then begin
            o2i = 1
            o2j = 0
            loggfratio = 0.959
        endif else begin
            o2i = 0
            o2j = 1
            loggfratio = 1 / 0.959
        endelse
        
        w1 = where(science.lambda ge o2lambda[o2i]-0.5 and science.lambda le o2lambda[o2i]+0.5)
        sublambda1 = science.lambda[w1] - o2lambda[o2i]
        subspec1 = 1.0 - science.spec[w1]
        suberr1 = (science.ivar[w1])^(-0.5)
        guess = [0.3, 0.0, 0.05, 0.0]
        fit1 = gaussfit(sublambda1, subspec1, a, measure_errors=suberr1, estimates=guess, nterms=4)
        science.spec[w1] /= (1.0 - fit1)
        science.ivar[w1] *= (1.0 - fit1)^2.
        
        w0 = where(science.lambda ge o2lambda[o2j]-0.5 and science.lambda le o2lambda[o2j]+0.5)
        sublambda0 = science.lambda[w0] - o2lambda[o2j]
        subspec0 = 1.0 - science.spec[w0]
        suberr0 = (science.ivar[w0])^(-0.5)
        fit0 = interpol(fit1, sublambda1, sublambda0) * loggfratio ;0.220 / 0.230
        science.spec[w0] /= (1.0 - fit0)
        science.ivar[w0] *= (1.0 - fit0)^2.
    endif
        
    self.science = ptr_new(science)
    self->getlinelist
end


pro hiresspec3::getlinelist, reload=reload
    e = elements()
    e = [e, e[0]]
    wcn = n_elements(e)-1
    e[wcn].atomic = 607
    e[wcn].name = 'CN'
    ionnames = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X']

    science = *self.science
    name = (*self.names)[self.i]

    linelist_initials = 'APJ'
    
    ewfile = self.outdir+name+'_Ji20_ew.fits'
    if 1 and file_test(ewfile) then begin
        linelist = mrdfits(ewfile, 1, /silent)
    endif else begin
        linelist = {lambda:0d, species:0d, atomic:0l, ion:0l, ep:0d, loggf:0d, ew:-999d, ewerr:-999d, doppwidth:-999d, doppwidtherr:-999d, rate:-999d, rateerr:-999d, name:' ', upperlimit:0, manual:0, dlambda:0d, localcont:0d, localconterr:-999d}
        case linelist_initials of
            'JLC': begin
                readcol, getenv('CALTECH')+'linelist/linelist.jlc', lambdaj, speciesj, epj, loggfj, dampingj, skipline=1, format='D,D,D,D,D', /silent
                nlines = n_elements(lambdaj)
                linelist = replicate(linelist, nlines)
                linelist.lambda = lambdaj
                linelist.species = speciesj
                linelist.atomic = floor(speciesj)
                linelist.ion = long((speciesj - floor(speciesj))*10.)
                linelist.ep = epj
                linelist.loggf = loggfj
                ;linelist.damping = dampingj
                linelist = linelist[sort(linelist.lambda)]
            end
            'AOT': begin
                linelist = replicate(linelist, 10000)
                j = 0L
                openr, lun, getenv('CALTECH')+'hires/aot/GES-Lind-Clean-Fe.list', /get_lun
                skip_lun, lun, 1, /lines
                while ~eof(lun) do begin
                    readf, lun, lambda, species, ep, loggf, format='(D8,1X,D5,1X,D7,1X,D9)'
                    linelist[j].lambda = lambda
                    linelist[j].species = species
                    linelist[j].ep = ep
                    linelist[j].loggf = loggf
                    j++
                endwhile
                linelist = linelist[0:j-1]
                close, lun
                free_lun, lun

                linelist.atomic = floor(linelist.species)
                linelist.ion = long((linelist.species - floor(linelist.species))*10.)
                linelist = linelist[sort(linelist.lambda)]
            end
            'APJ': begin
                readcol, getenv('CALTECH')+'hires/M92_KOA/Ji20_linelist.moog', lambdaj, speciesj, epj, loggfj, skipline=1, format='D,D,D,D', /silent
                nlines = n_elements(lambdaj)
                linelist = replicate(linelist, nlines)
                linelist.lambda = lambdaj
                linelist.species = speciesj
                linelist.atomic = floor(speciesj)
                linelist.ion = long((speciesj - floor(speciesj))*10.)
                linelist.ep = epj
                linelist.loggf = loggfj
                ;linelist.damping = dampingj
                linelist = linelist[sort(linelist.lambda)]
            end            
            'IUR': begin
                readcol, getenv('CALTECH')+'hires/M92_KOA/Roederer2018_linelist.moog', lambdaj, speciesj, epj, loggfj, skipline=1, format='D,D,D,D', /silent
                nlines = n_elements(lambdaj)
                linelist = replicate(linelist, nlines)
                linelist.lambda = lambdaj
                linelist.species = speciesj
                linelist.atomic = floor(speciesj)
                linelist.ion = long((speciesj - floor(speciesj))*10.)
                linelist.ep = epj
                linelist.loggf = loggfj
                ;linelist.damping = dampingj
                linelist = linelist[sort(linelist.lambda)]
            end            
        endcase
        science = *self.science
        ;linelist = linelist[where(linelist.lambda ge min(science.lambda) and linelist.lambda le max(science.lambda))]
        keepe = intarr(n_elements(e.atomic)) + 1
        for i=0,n_elements(e.atomic)-1 do begin
            w = where(linelist.atomic eq e[i].atomic, cw)
            if cw gt 0 then linelist[w].name = e[i].name + ' ' else keepe[i] = 0
        endfor
        e = e[where(keepe)]
        for i=0,n_elements(ionnames)-1 do begin
            w = where(linelist.ion eq i, cw)
            if cw gt 0 then linelist[w].name = linelist[w].name + ionnames[i] + ' '
        endfor
    endelse
    keepe = intarr(n_elements(e.atomic)) + 1
    for i=0,n_elements(e.atomic)-1 do begin
        w = where(linelist.atomic eq e[i].atomic, cw)
        if cw gt 0 then linelist[w].name = e[i].name + ' ' else keepe[i] = 0
    endfor
    e = e[where(keepe)]
    for i=0,n_elements(ionnames)-1 do begin
        w = where(linelist.ion eq i, cw)
        if cw gt 0 then linelist[w].name = linelist[w].name + ionnames[i] + ' '
    endfor

    linelist.name += strcompress(round(linelist.lambda), /rem)
    self.atomic = ptr_new(e.atomic)

    colors = replicate('black', n_elements(e))
    ;alpha
    colors[where(e.name eq 'O' or e.name eq 'Ne' or e.name eq 'Mg' or e.name eq 'Si' or e.name eq 'S' or e.name eq 'Ar' or e.name eq 'Ca' or e.name eq 'Ti')] = 'green'
    ;Fe
    colors[where(e.name eq 'Fe')] = 'red'
    ;Fe peak
    colors[where(e.name eq 'Sc' or e.name eq 'V' or e.name eq 'Cr' or e.name eq 'Mn' or e.name eq 'Co' or e.name eq 'Ni' or e.name eq 'Cu')] = 'blue'
    ;s-process
    colors[where(e.name eq 'Sr' or e.name eq 'Y' or e.name eq 'Zr' or e.name eq 'Ba' or e.name eq 'La' or e.name eq 'Ce' or e.name eq 'Pr' or e.name eq 'Pb' or e.name eq 'Bi')] = 'orange'
    ;r-process
    colors[where(e.name eq 'Eu' or e.name eq 'Se' or e.name eq 'Br' or e.name eq 'Kr' or e.name eq 'Te' or e.name eq 'I' or e.name eq 'Os' or e.name eq 'It' or e.name eq 'Pt')] = 'purple'
    ;molecules
    ;colors[where(e.name eq 'CN')] = 'cyan'
    self.colors = ptr_new(colors)

    ptr_free, self.linelist
    self.linelist = ptr_new(linelist)
end


; =============== INIT ================
function hiresspec3::INIT, catalogfile, indir=indir, outdir=outdir, big=big
    catalog = mrdfits(catalogfile, 1, /silent)
    catalog = catalog[sort(catalog.gmag0)]
    names = strtrim(catalog.name, 2)
    self.names = ptr_new(strtrim(catalog.name, 2))
    self.specfiles = ptr_new(strtrim(catalog.specfile, 2))
    self.vr = ptr_new(catalog.vr)
    
    self.indir = indir
    if keyword_set(outdir) then self.outdir = outdir else self.outdir = self.indir

    base = widget_base(/row, title='hiresspec3', uvalue=self, mbar=menu, tab_mode=0);, /scroll, x_scroll_size=12.5, y_scroll_size=6, units=1)
    file_menu = widget_button(menu, value='File', /menu)
    wexport = widget_button(file_menu, value='Save Changes', uname='export', uvalue='export')
    wreload = widget_button(file_menu, value='Reload Linelist', uname='reload', uvalue='reload')
    wps = widget_button(file_menu, value='Make EPS', uname='eps', uvalue='eps')
    wexit = widget_button(file_menu, value='Exit', uvalue='exit', uname='exit')
    default_menu = widget_button(menu, value='Control', /menu)
    wdefaultlambda = widget_button(default_menu, value='Measure All', uname='measure_all', uvalue='measure_all')
    wdefaultlambda = widget_button(default_menu, value='Default Spectrum Settings', uname='default_lambda', uvalue='default_lambda')
    ;wdefaultboth = widget_button(default_menu, value='Both', uname='default_both', uvalue='default_both')

    wright = widget_base(base, /column, uname='right')
    widget_control, /managed, base

    ; ------ RIGHT -------
    wfile = widget_base(wright, /frame, /row, /align_left, tab_mode=1)
    wback = widget_button(wfile, value='Back', uvalue='back', uname='back', tab_mode=1)
    wobjlist = widget_combobox(wfile, uname='objlist', value='                 ', tab_mode=1, /dynamic_resize)
    wnext = widget_button(wfile, value='Next', uvalue='next', uname='next', tab_mode=1)
    wstatus = widget_text(wfile, xsize=128, value='Initializing ...', uname='status', uvalue='status', tab_mode=0)
    wlines = cw_bgroup(wfile, 'Voigt', /nonexclusive, /return_name, /frame, uname='voigt')

    wspec = widget_base(wright, /frame, /column)
    spawn, 'uname -n', hostname
    if keyword_set(big) or strtrim(hostname, 2) eq 'stravinsky' then begin
        wspecplot = widget_draw(wspec, xsize=3820, ysize=1940, uname='spec', /button_events, keyboard_events=1)
    endif else begin
        wspecplot = widget_draw(wspec, xsize=1900, ysize=870, uname='spec', /button_events, keyboard_events=1)
    endelse

    wspeccontrol = widget_base(wright, /row, /align_center, tab_mode=1)
    wycontrol = widget_base(wspeccontrol, /frame, /row)
    wylow = widget_text(wycontrol, xsize=5, /editable, uname='ylow', uvalue='ylow')
    wylabel = widget_label(wycontrol, value=' < y < ', /align_center, uname='ylabel')
    wyhigh = widget_text(wycontrol, xsize=5, /editable, uname='yhigh', uvalue='yhigh')
    wlambdacontrol = widget_base(wspeccontrol, /frame, /row)
    wblue = widget_button(wlambdacontrol, value='<-', uname='blue', uvalue='blue', /align_center)
    wlambdalow = widget_text(wlambdacontrol, xsize=8, /editable, uname='lambdalow', uvalue='lambdalow')
    wlambdalabel = widget_label(wlambdacontrol, value=' < lambda < ', /align_center, uname='lambdalabel')
    wlambdahigh = widget_text(wlambdacontrol, xsize=8, /editable, uname='lambdahigh', uvalue='lambdahigh')
    wred = widget_button(wlambdacontrol, value='->', uname='red', uvalue='red', /align_center)

    ;widget_control, base, set_uname=self
    widget_control, base, /realize
    self.base = base
    xmanager, 'hiresspec3', self.base, /no_block, cleanup='hiresspec3_cleanup'

    self.measurelim = [-1d, -1d]
    self.i = 0
    widget_control, widget_info(self.base, find_by_uname='objlist'), set_value=*self.names
    self->getscience
    self->default_lambda;, /noredraw
    return, 1
end

pro hiresspec3__define 
    state = {hiresspec3, $
             i:0L, $
             names:ptr_new(), $
             specfiles:ptr_new(), $
             vr:ptr_new(), $
             indir:' ', $
             outdir:' ', $
             base:0L, $
             science:ptr_new(), $
             linelist:ptr_new(), $
             change:0b, $
             atomic:ptr_new(), $
             colors:ptr_new(), $
             lambdalim:[-100d, 100d], $
             ylim:[-100d, 100d], $
             measurelim:[-1d, -1d], $
             o2lambda:ptr_new()}
end


pro hiresspec3_old, big=big
    dir = getenv('CALTECH')+'hires/M92_KOA/'
    indir = dir+'spectra/'
    outdir = dir+'ew/'
    allframesfile = dir+'M92_allframes.fits'
    n = obj_new('hiresspec3', allframesfile, indir=indir, outdir=outdir, big=big)
end


pro hiresspec3, big=big
    dir = getenv('CALTECH')+'hires/M15_M92/'
    indir = getenv('chome')+'keck/hires/M15_M92/'
    outdir = dir+'ew/'
    allframesfile = dir+'M15_M92_allframes.fits'
    n = obj_new('hiresspec3', allframesfile, indir=indir, outdir=outdir, big=big)
end
