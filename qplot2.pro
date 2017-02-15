;+
; NAME:
;   qplot
;
; PURPOSE:
;   Examine reduced hs data
;
; CALLING SEQUENCE:
;
;qplot,'spHect-test.fits',xrange=[3600,9000],/dosnr,$
;        plugmap='/data/cnaw/mmt/2006.0502/plugdir/new_egs_hecto_spec_2',$
;        nsmooth=11,spzbest='spZbest-test.fits',$
;        catalog='/data/cnaw/mmt/2006.0502/plugdir/new_egs_hecto_spec_2_cat'
; INPUTS:
;  filename 6   - filename for for the Hectospec data
;
; OPTIONAL KEYWORDS:
;   nstart - starting aperture
;   nsmooth - number of pixels to smooth over
;   xrange - plotting xrange
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;
;
; EXAMPLES:
;
; BUGS:
;
;;
; REVISION HISTORY:
;   April 2004 - Written by R Cool - UofA
;   June 2005 - Feature added by C Papovich.  If you pass it the
;               spZbest-*.fits file, the code will use that redshift
;               to plot lines automatically.
;
;   June 2006 - Several additions by CNAW. Making it a Frankenstein
;               from RVSAO qplot and DEEP2 zspec
;
;   June 2008 - Cleanup by BJW.  rationalized/made new subroutines,
;               add some em lines, auto-saving, adding comments,
;               try to rationalize some looping/reset plot limits
;               behavior
;
;------------------------------------------------------------------------------
pro  recover_cat, filename, obscat = obscat
;
; read catalogue so that object ids can be tagged to spectra
;
junk       = 'x'
str = strarr(8)
openr, lun, filename,/get_lun
temp  = {obscat}
count = 0
WHILE (NOT EOF(lun)) DO BEGIN
    readf, lun, junk
    str = strsplit(junk,' ',/extract)
    temp.aperture   = fix(str[0])
    temp.ra         = float(str[1])
    temp.dec        = float(str[2])
    temp.mag        = float(str[3])
    temp.fiber      = fix(str[4])
    temp.beam       = fix(str[5])
    temp.id         = str[6]
    if (count eq 0 ) then begin
        obscat = temp
    endif else begin
        obscat  = [obscat, temp]
    endelse
    count = count + 1
ENDWHILE
free_lun, lun
end

;
;=====================================================================
;
pro show_em_lines, z, ymin, diff

; wavelengths are in vacuum.  C IV and Mg II are averages
cyan   = fsc_color("cyan",!D.Table_Size- 5)  
lines = [1033.8, 1215.67, 1240.149, 1399.8, 1549.053, 1640.4,$
         1908.734, 2326.0, 2799.9, 3426.0,$
         3727.1, 3729.9, 3870.2, $
         4341.68, 4862.68, 4960.3, 5008.24, 6549.86,$
         6564.61, 6585.27, 6718.29 ,6732.67]
labels = ['Ly \beta', 'Ly \alpha', 'N V', 'Si IV', 'C IV', $
          'He II', 'C II', 'C II','MgII','NeV',' ','O II', 'NeIII',$
          'H \gamma','H \beta', 'O III', 'O III', 'N II',$
          'H \alpha', 'N II', 'S II', 'S II']
for j = 0, n_elements(lines) -1 do begin
    djs_oplot, (1+Z)* [lines(j), lines(j)], [-1e4,1e5], $
      color=cyan, linestyle=2
    djs_xyouts, (1+Z)*lines(j), ymin+diff*(0.1+0.05*(-1)^j), $
      labels(j), charsize=1, color=cyan
endfor 
end

;
;=====================================================================
;
pro show_abs_lines, z, ymin, diff

; wavelengths are in vacuum.
orange = fsc_color("orange",!D.Table_Size-8)
lines= [3799.0, 3836.5, 3890.2, 3934.8, 3969.6, 4102.9,$
        4305.6, 4341.7, 4862.7, 5176.8, 5270.5, 5894.1,$
        4173.2, 4227.9, 8500.3, 8544.4, 8664.6]
labels = ['H8', 'H7' , 'H6' ,  'K',  'H',  'H \delta',  $
          'G',  'H \gamma',  'H \beta',  'Mg',  'CaFe', $
          'Na',  'CaI',  'CaI','CaI','Ca','Ca','Ca']
for j = 0, n_elements(lines) -1 do begin
    djs_oplot, (1+Z)* [lines(j), lines(j)], [-1e4,1e5], $
      color=orange, linestyle=2
    djs_xyouts, (1+Z)*lines(j), ymin+diff*(0.1+0.05*(-1)^j), $
      labels(j), charsize=1, color=orange
endfor
end

;
;=====================================================================
;
pro show_sky_lines, ymin, diff

; wavelengths should be in vacuum but are not yet
grey   = fsc_color("grey", !D.Table_Size-9)
lines = [3650. , 3663. , 4047. , 4078. , 4358. , 4978. , 5460. ,$
         5539. , 5577. , 5619. , 5683., 5891. , 6300. , 6363.0, 6831.3,$
         6870. , 7242.5, 7276.3, 7600. , 7715.0, 7750.6, 8344.4, 8827.0]
labels = ['Hg','Hg','Hg','Hg','Hg','Na','Hg', $
          'sky','OI','sky','sky','Na','OI','OI','sky', $
          'B','sky','sky','A','sky','sky','sky','sky']
for j = 0, n_elements(lines) -1 do begin
    djs_oplot, [lines(j), lines(j)], [-1e4,1e5], $
      color=grey, linestyle=2
    djs_xyouts, lines(j), ymin+diff*(0.1+0.05*(-1)^j), $
      labels(j), charsize=1, color=grey
endfor 
end

;
;=====================================================================
;
pro save_logfile, logfile, qcat, log_format
  openw, lun, logfile,/get_lun
  for i=0, n_elements(qcat)- 1 do begin
      printf,lun, qcat[i].id, qcat[i].ra, qcat[i].dec, $
             qcat[i].mag, qcat[i].z, qcat[i].z_err, qcat[i].zq, $
             qcat[i].aperture, qcat[i].fiber,$
             qcat[i].zwarning,qcat[i].class,qcat[i].comment,format=log_format
  endfor
  free_lun,lun
end


;
;=====================================================================
;
pro show_commands, showspz
  splog, 'Commands:',/noname
  splog, ' carriage return - next spectrum',/noname
  splog, ' N - redshift from next template',/noname
  splog, ' O - save logfile',/noname
  splog, ' P - previous spectrum',/noname
  splog, ' I - print current index',/noname
  splog, ' J - jump to index J',/noname
  splog, ' Q - quit',/noname
  splog, ' R - restore x wavelength domain',/noname
  splog, ' X - Set X-limits',/noname
  splog, ' Y - Set Y-limits',/noname
  splog, ' W - set X and Y limits clicking mouse ',/noname
  splog, ' A - plot  absorption lines',/noname
  splog, ' Z - plot  emission lines',/noname
  splog, ' M - show pixel mask',/noname
  splog, ' L - clear lines',/noname
  splog, ' S - Smooth',/noname
  splog, ' E - use line to estimate z click mouse on line ',/noname
  splog, ' C - Type a comment',/noname
  splog, ' 4 - z is good ',/noname
  splog, ' 3 - z is almost good',/noname
  splog, ' 2 - z is questionable',/noname
  splog, ' 1 - z cannot be measured',/noname
  splog, ' [ - plot sky lines ',/noname
  if showspz eq 1 then begin
      splog, ' ZA - use Z-best and plot absorption lines',/noname
      splog, ' ZE - use Z-best and plot emission lines',/noname
      splog, ' ]  - use Z-best and plot all galaxy lines',/noname
  endif
  splog, ' ? - show these options',/noname
end

;
;=====================================================================
;
PRO qplot2, filename, nstart=nstart, plugmap=plugmap, ivar=ivar, $
           zmap=zmap, autoline=autoline, smin = smin, synth=synth, $
           nsmooth=nsmooth, psym=psym,title=title, xrange=xrange, $
           clipfact=clipfact, dosnr=dosnr, pseudoflux=pseudoflux, $
           synthmag=synthmag, spzbest=spzbest, $
           doflags=doflags, flags=flags, catalog=catalog, noupdate=noupdate
;
; Number of templates used in the calculation of x
;
ntemplates = 51
n = 0
;
; Make plots colorful
;
white  = fsc_color("white",!D.Table_Size-1)
red    = fsc_color("red", !D.Table_Size-2 )
green  = fsc_color("green",!D.Table_Size- 3)
blue   = fsc_color("blue",!D.Table_Size- 4)  
cyan   = fsc_color("cyan",!D.Table_Size- 5)  
brown  = fsc_color("brown",!D.Table_Size- 6)  
yellow = fsc_color("yellow",!D.Table_Size-7)
orange = fsc_color("orange",!D.Table_Size-8)
grey   = fsc_color("grey", !D.Table_Size-9)
;
; Get a decent window size
;
window,0,xsize=1500,ysize=500
;
; Get results from SDSS z calculation
;
if keyword_set(spzbest) then begin
    print,'reading ', spzbest
    spzbest=spzbest
    spz=mrdfits(spzbest,1)
endif
;
; Read source catalog for this setup
;
if keyword_set(catalog) then begin
    catalog = catalog
    recover_cat, catalog, obscat=obscat
endif

if not keyword_set(xrange) then xrange=[3800,8500]
if not keyword_set(psym)   then psym=0
if not keyword_set(smin)   then smin = 3650
if NOT keyword_set(nstart) then nstart = 0
if NOT keyword_set(clipfact) then clipfact =0.6
if not keyword_set(doflags) then doflags=0

if keyword_set(doflags) then begin
    if (size(flags))[1] ne 300 then $
      flags = intarr(n_elements(spz))
endif

;
; output format

; das_format = '(a25,1x,f13.9,1x,f12.9,1x,f7.2,1x,f9.6,1x,f10.6,1x,i2,1x,i3,1x,i3,1x,i4,1x,a10)'
das_format = '(a25,1x,f13.9,1x,f12.9,1x,f7.2,1x,f9.6,1x,f10.6,1x,i2,1x,i3,1x,i3,1x,i4,1x,a10,1x,a60)'
;
; These variables allow rescaling the plot
;
xmin_org = xrange[0]
xmax_org = xrange[1]
myxmin   = xmin_org
myxmax   = xmax_org

myymin = 1d40
myymax = -1d40
ymin   = -1
ymax   =  1
;
; Read file with reduced spectra
;
print,'reading ', filename
lam = mrdfits(filename, 0,/silent)
object = mrdfits(filename, 1,/silent)
ivar = mrdfits(filename, 2,/silent)
pixmask = mrdfits(filename, 4,/silent) + mrdfits(filename,3,/silent)
mask = mrdfits(filename, 4,/silent)*0.0
plugmap = mrdfits(filename, 5,/silent)


if keyword_set(pseudoflux) then begin
    fcalfile = getenv('HSRED_DIR') + '/etc/average_flux.fits'
    tcalfile = getenv('HSRED_DIR') + '/etc/average_tell.fits'
     fcal = mrdfits(fcalfile, 1,/silent)
    tcal = mrdfits(tcalfile, 1,/silent)
    flux_factor = bspline_valu(alog10(lam), fcal)
    tell_factor = bspline_valu((lam), tcal)
    divideflat, object, flux_factor, invvar=ivar
    divideflat, object, tell_factor, invvar=ivar
endif
;
;**** Open the output file to which results will be written. ****
;  
length   = strlen(filename)
logfile  = strmid(filename,0,length-4)+'zlog'
bkpfile  = logfile+'.bkp'
splog, logfile
; 
; If the output file already exists, read the data so new 
; classifications can be appended.
;
file = findfile(logfile,count=exists)
print,'exists ', exists,' ', logfile
if (exists gt 0 ) then begin
    read_qplot_file, logfile, qcat=qcat
    print,  n_elements(qcat)
    for i=0, n_elements(qcat)-1 do begin
;        if qcat[i].id eq 'noid' then break
        if qcat[i].zq le 0 then break
    endfor
    if i ge 1 then begin
        print, "last object examined was",i-1
        print, qcat[i-1]
    endif
endif else begin
    qcat = {qcat}
    qcat.id='noid'
    qcat.class='none'
    print, qcat
    qcat = replicate(qcat,300)
endelse
;  print, ' qcat ',n_elements(qcat)
openu, lun, bkpfile,/get_lun,/append
;
; Index of commands
;
  test = 'ZA'
  qcomment = ''
  qcommtemp = ''

  if keyword_set(spzbest) then begin
      show_commands, 1
  endif else begin
      show_commands, 0
  endelse
  if  keyword_set(showsky) then begin
      show_lines = ']'
  endif else begin
      show_lines = ''
  endelse

  n          =  0
;  splog, '? - show these options'
  
  i = nstart
  stop = 'go'
  
; calculate s/n if keyword is set
  if keyword_set(dosnr) then hs_snr, lam, object, plugmap, snr=snr

; otherwise ??

  if not keyword_set(dosnr) then snr = fltarr(n_elements(object(0,*)))
  
  if NOT keyword_set(nsmooth) then nsmooth=0
  nsmooth_orig = nsmooth
  splog, '%% Current Smoothing = '+strtrim(string(nsmooth),2)

  plot_pixmask=0
  loop = 0
  newobj = 1
  while stop ne 'stop' do begin
      
      print, 'i ', i
      if i ge 300 then goto, jumpout 

;
; Attempt at not losing an already set value of z. This is
; particularly important when the best z does not come from the first
; template. Needs testing  CNAW 2007-08-03
;
; BJW: there are several things that set loop=1 and should preserve z_choice
; e.g. 'w' 'x' 'y'
; but also things that should reset z, e.g. those that go to a different
; object, 'p' and 'j'
;
;      if(loop ne 2) then begin
      if(loop eq 0 or newobj eq 1) then begin
          z_choice = qcat[i].z
          z_choice_err = qcat[i].z_err
          zqual = qcat[i].zq 
          qcomment = qcat[i].comment
;
; Redshift from Xcor
;
          spzindex = long(i * ntemplates + n)
          print, 'spzindex ',spzindex
          print, ' z = ', spz[spzindex].z
          newobj = 0
      endif
;
; if loop is 0 or if the redshift is not defined get a value from
; the redshift determination.  However, if there is already a
; zquality, preserve the z in qcat, as set above - this happens if you
; go out of order.
;
      if((loop eq 0 or qcat[i].z eq 0) and qcat[i].zq le 0 and zqual le 0) then begin
          z_choice = spz[spzindex].z
          z_choice_err = spz[spzindex].z_err
      endif
      print, ' zplot, z in qcat ', z_choice, qcat[i].z, qcat[i].zq

;Mask pixels were the ivar is lt 0.6 the smoothed value
      temp = ivar[*,i] / smooth(ivar[*,i], 21)
      newmask = temp lt clipfact
      
      junk = min(abs(lam[*,i]-5588),pix)

      mask[(pix-10):(pix+10),i]=mask[(pix-10):(pix+10),i]+$
        newmask[(pix-10):(pix+10)]
      
      junk = min(abs(lam[*,i]-5577), pix)
      
      mask[(pix-10):(pix+10),i] = mask[(pix-10):(pix+10),i] + $
        newmask[(pix-10):(pix+10)]
      junk = min(abs(lam[*,i]-6300), pix)
      
      mask[(pix-10):(pix+10),i] = mask[(pix-10):(pix+10),i] + $
        newmask[(pix-10):(pix+10)]
      
      junk = min(abs(lam[*,i]-6363), pix)
      
      mask[(pix-10):(pix+10),i] = mask[(pix-10):(pix+10),i] + $
        newmask[(pix-10):(pix+10)]
      
      pixels = where(mask[*,i] eq 0)  
      
;
; Rescale x.  BJW: The test for myxmax gt myxmin is in case somebody
; clicks right corner before left corner on 'w'
;   
      if myxmax gt myxmin then begin
          xrange[0] = myxmin
          xrange[1] = myxmax
      endif else begin
          xrange[0] = myxmax
          xrange[1] = myxmin
      endelse
      
      l1 =  min(abs(lam[pixels,i]-xrange(0)),l)
      m1 =  min(abs(lam[pixels,i]-xrange(1)),m)
      
      l = l(0)
      m = m(0)
      
      mypix = ivar[pixels(l):pixels(m),i]
      t = where(mypix gt 0)
      if ( t[0] ge 0 ) then begin 
          mypix=mypix[t]
          ymax2 = max(smooth(1./sqrt(mypix),nsmooth)) *1.1
          ymin2 = min(smooth(1./sqrt(mypix),nsmooth))
          diff2 = ymax2-ymin2    
          t = where(ivar[pixels,i] gt 0)
          ymax = max(smooth(object[pixels(l):pixels(m),i],nsmooth)) *1.1
          ymin = min(smooth(object[pixels(l):pixels(m),i],nsmooth))
      endif else begin
          ymin2 = -0.5 
          ymax2 = 2.5
      endelse
      if myymin lt 1d9 and myymax gt -1d9 then begin
          if myymax gt myymin then begin
              ymax=myymax
              ymin=myymin
          endif else begin
              ymax=myymin
              ymin=myymax
          endelse
          myymin = 1d40
          myymax = -1d40
      endif
      diff = ymax-ymin    
      
      if keyword_set(title) then zetitle1=title(i)
      if not keyword_set(title) then title1=''
;
; Plot box
;
      plot, lam[pixels,i], smooth(object[pixels,i],nsmooth),$
            /xstyle , yrange=[ymin, ymax],  ps=psym, title=title1, $
            xrange=xrange,color=white, /nodata
;
; plot object spectrum
;
      oplot, lam[pixels,i], smooth(object[pixels,i],nsmooth),$
             thick=1,color=red
;
; plot smoothed spectrum ?
;
      if plot_pixmask then begin
          oplot, lam[pixels,i], smooth(pixmask[pixels,i],nsmooth),$
                 thick=1,color=grey,line=1
      endif
;
; Plot inverse variance
;
      if (t[0] ge 0) then begin
          oplot, lam[pixels[t],i], smooth(1./sqrt(ivar[pixels[t],i]),nsmooth),$
                 color=blue,thick=1
      endif
;
; Plot aperture number
;
      djs_xyouts, xrange(0) + (xrange(1)-xrange(0))*0.05,$
                  ymin + (ymax-ymin)*0.95, 'Aperture  ' + strn(fix([i+1])), $
                  color=yellow, thick=2, charsize=2, /isolatin1
      
      if keyword_set(plugmap) then begin
;
; object type 
;
          djs_xyouts, xrange(0) + (xrange(1)-xrange(0))*0.05,$
                      ymin + (ymax-ymin)*0.90  ,  plugmap[i].objtype, $
                      color=yellow, thick=2, charsize=2, /isolatin1
;
; S/N ratio
;
          djs_xyouts,  xrange(0) + (xrange(1)-xrange(0))*0.05,$
                       ymin + (ymax-ymin)*0.85, $
                       'snr = ' + strn(snr(i)), $
                       color=yellow, thick=2, charsize=2, /isolatin1
; Object ID
;
          djs_xyouts,  xrange(0) + (xrange(1)-xrange(0))*0.25,$
                       ymin + (ymax-ymin)*0.95, $
                       'id = '  + strn(obscat[i].id), $
                       color=yellow, thick=2, charsize=2, /isolatin1
;
; apparent magnitude
;
          djs_xyouts,  xrange(0) + (xrange(1)-xrange(0))*0.25,$
                       ymin + (ymax-ymin)*0.9, $
                       'ap mag = ' + strn(plugmap[i].rapmag), $
                       color=yellow, thick=2, charsize=2, /isolatin1
;
; redshift
;
          djs_xyouts,  xrange(0) + (xrange(1)-xrange(0))*0.25,$
                       ymin + (ymax-ymin)*0.85, $
                       'z      = ' + strn(z_choice), $
                       color=yellow, thick=2, charsize=2, /isolatin1
;                       'z      = ' + strn(spz[spzindex].z), $
;
; quality
;
          djs_xyouts,  xrange(0) + (xrange(1)-xrange(0))*0.25,$
                       ymin + (ymax-ymin)*0.80, $
                       'Q     = ' + strn(qcat[i].zq), $
                       color=yellow, thick=2, charsize=2, /isolatin1
;
; identify lines 
;
          if (show_lines eq ']') then begin
;              z = spz[spzindex].z
              z = z_choice
              show_em_lines, z, ymin, diff
              show_abs_lines, z, ymin, diff
              show_lines = ']'
          endif
          
;Here, I am going to do something a bit strange and go ahead and
;integrate the spectrum at the redshift observed.  This only works 
;if the /synthmag is set
        
          if keyword_set(synthmag) then begin
              
              k_load_filters, 'sdss_r0.par', nlam, lam1, trans
              integrand = object[*,i]*1e-17
              intertrans = interpol(trans, lam1, lam[*,i])
              k = where(intertrans lt 0)
              intertrans(k) = 0
              rmag = int_tabulated(lam[*,i], integrand*intertrans)
              djs_xyouts,  xrange(0) + (xrange(1)-xrange(0))*0.45,$
                           ymin + (ymax-ymin)*0.95, $
                           'synth mag = ' + strn(-2.5*alog10(rmag)), $
                           color=7, thick=2, charsize=2, /isolatin1
          endif
      endif
      
;          zqual= 0L
      loop = 0 
      while loop eq 0 do begin
          read, prompt='Command: ', test
          
          if strupcase(test) eq 'I' then begin
              splog, "Current Index : "+strtrim(string(i),2)
          endif
          
          if test eq '?' then begin
              if keyword_set(spzbest) then begin
                  show_commands, 1
              endif else begin
                  show_commands, 0	
              endelse
          endif
          if strupcase(test) eq '1' then begin 
              zqual = 1L
              qcat[i].zq = zqual
          endif
          if strupcase(test) eq '2' then begin
              zqual = 2L
              qcat[i].zq = zqual
          endif
          if strupcase(test) eq '3' then begin
              zqual = 3L
              qcat[i].zq = zqual
          endif
          if strupcase(test) eq '4' then begin
              zqual = 4L
              qcat[i].zq = zqual
          endif          
          
;
; Go on to next object
;
          if strupcase(test) eq '' then begin
;
; however, before doing that, update logfile with redshift quality
;
              aperture = i + 1
              qcat[i].id       = obscat[i].id
              qcat[i].ra       = obscat[i].ra
              qcat[i].dec      = obscat[i].dec
              qcat[i].mag      = obscat[i].mag
;              qcat[i].z        = spz[spzindex].z
;              qcat[i].z_err    = spz[spzindex].z_err
              qcat[i].z        = z_choice
              qcat[i].z_err    = z_choice_err
              qcat[i].zq       = zqual
              qcat[i].zwarning = spz[spzindex].zwarning
              qcat[i].aperture = aperture
              qcat[i].class    = spz[spzindex].class
              qcat[i].comment  = qcomment
              
; this is the backup file
              printf,lun, obscat[i].id, obscat[i].ra, obscat[i].dec,$
                     obscat[i].mag, spz[spzindex].z, spz[spzindex].z_err,$
                     zqual, aperture, obscat[i].fiber,spz[spzindex].zwarning,$
                     spz[spzindex].class,qcomment,format=das_format
; save the actual output file after every object
              if (not keyword_set(noupdate)) then $
                save_logfile, logfile, qcat, das_format
              
              loop = 1
              newobj = 1
              nsmooth=nsmooth_orig
              plot_pixmask=0
              n = 0
              if keyword_set(doflags) then begin
                  read,prompt='Enter Flag value [int]: ', tflag
                  flags[i] = tflag
              endif
; template number 
              i = i+ 1
              if i ge 300 then goto, jumpout 
              spzindex=long(i * ntemplates + n)
              z_choice     = spz[spzindex].z
              z_choice_err = spz[spzindex].z_err
              zqual = 0L
              myxmin = xmin_org
              myxmax = xmax_org
              show_lines = ']'
          endif
;
;         save results without exiting
;

          if strupcase(test) eq 'O' then begin
              aperture = i + 1
              qcat[i].id       = obscat[i].id
              qcat[i].ra       = obscat[i].ra
              qcat[i].dec      = obscat[i].dec
              qcat[i].mag      = obscat[i].mag
;              qcat[i].z        = spz[spzindex].z
;              qcat[i].z_err    = spz[spzindex].z_err
              qcat[i].z        = z_choice
              qcat[i].z_err    = z_choice_err
              qcat[i].zq       = zqual
              qcat[i].zwarning = spz[spzindex].zwarning
              qcat[i].aperture = aperture
              qcat[i].class    = spz[spzindex].class
              qcat[i].comment  = qcomment
              
              printf,lun, obscat[i].id, obscat[i].ra, obscat[i].dec,$
                     obscat[i].mag, spz[spzindex].z, spz[spzindex].z_err,$
                     zqual, aperture, obscat[i].fiber,spz[spzindex].zwarning,$
                     spz[spzindex].class,qcomment,format=das_format
              
;              printf,lun, obscat[i].id, obscat[i].ra, obscat[i].dec, $
;                     obscat[i].mag, spz[spzindex].z, spz[spzindex].z_err, zqual, $
;                     aperture, obscat[i].fiber,$
;                     spz[spzindex].zwarning,spz[spzindex].class,format=das_format
              
              openw, unit, logfile,/get_lun
              for nn=0, n_elements(qcat)- 1 do begin
                  printf,unit, qcat[nn].id, qcat[nn].ra, qcat[nn].dec, $
                         qcat[nn].mag, qcat[nn].z, qcat[nn].z_err, qcat[nn].zq, $
                         qcat[nn].aperture, qcat[nn].fiber,$
                         qcat[nn].zwarning,qcat[nn].class,qcomment,format=das_format
              endfor
              free_lun,unit
              loop = 1
          endif
          
          if strupcase(test) eq 'P' then begin
              i = i - 1
              loop = 1
              newobj = 1
              nsmooth=nsmooth_orig
              plot_pixmask=0
          endif
          
          if strupcase(test) eq 'J' then begin
              read,prompt='Enter aperture: ', j
              i = j-1
              loop = 1
              newobj = 1
              nsmooth=nsmooth_orig
              plot_pixmask=0
          endif
;
; Choose another template
;
          if strupcase(test) eq 'N' then begin
              print,' template, object, z , rchi2'
;              for k=0, ntemplates-1 do begin 
              for k=0, 20 do begin $
                  spzindex = i * ntemplates + k
                  print,k,' ', obscat[i].id,spz[spzindex].z, spz[spzindex].rchi2
              endfor
; catch establishes an error handler to prevent crashes if the 
; user types a letter instead of a number when reading 'n'
; the error is triggered by the read, ndummy statement, 
; not the is_digit statement
              Error_status = 0
              catch, Error_status
              if( Error_status ne 0) then begin
                  print, ' n must be numeric, setting n=0'
                  ndummy = 0
                  catch, /cancel
              endif
              read,prompt='Enter template: ', ndummy
              is_digit = ndummy*0
;              catch, /cancel
              n = ndummy
              myxmin = xmin_org
              myxmax = xmax_org
              loop=2
              spzindex = i * ntemplates + n
              z_choice        = spz[spzindex].z
              z_choice_err    = spz[spzindex].z_err
          endif
;
; Restore wavelength domain to default value
;
          if strupcase(test) eq 'R' then begin
              myxmin = xmin_org
              myxmax = xmax_org
              loop=1
          endif
;
; Blow up in X
;        
          if strupcase(test) eq 'X' then begin
              read,prompt='Enter X limits min, max: ',myxmin,myxmax
              loop=1
          endif
;
; Blow up in Y
;        
          if strupcase(test) eq 'Y' then begin
              read,prompt='Enter Y limits min, max: ',myymin,myymax
              loop=1
          endif
;
; Blow up in X and Y 
;        
          if strupcase(test) eq 'W' then begin
              cursor, myxmin, myymin, /down
              if (!mouse.button eq 1) then begin
                  cursor, myxmax, myymax,/down
                  print, myxmin, myymin,myxmax, myymax
                  loop=1
              endif
              if (!mouse.button eq 4) then begin
                  myxmin = xmin_org
                  myxmax = xmax_org
                  loop=1
              endif
          endif
          
; changed this from C to L
          if strupcase(test) eq 'L' then begin
;              i = i 
              loop = 1
              nsmooth=nsmooth_orig
              plot_pixmask=0
          endif
          
          if strupcase(test) eq 'M' then begin
              loop = 0
              plot_pixmask=1
              djs_oplot, lam[pixels,i], smooth(pixmask[pixels,i],nsmooth),$
                         thick=1,color=green,line=1
          endif

; Add a comment
          if strupcase(test) eq 'C' then begin
             print,'Current comment: ',qcomment
             read, prompt='New comment: ', qcommtemp
; convert spaces to underscores.  This is to make sure it is treated
; as one field when reading/writing qcat
             qcomment = strjoin(strsplit(qcommtemp,/extract),'_')
          endif
          
          if strupcase(test) eq 'STOP' or strupcase(test) eq 'Q' then begin
              aperture = i + 1
              qcat[i].id       = obscat[i].id
              qcat[i].ra       = obscat[i].ra
              qcat[i].dec      = obscat[i].dec
              qcat[i].mag      = obscat[i].mag
;              qcat[i].z        = spz[spzindex].z
;              qcat[i].z_err    = spz[spzindex].z_err
              qcat[i].z        = z_choice
              qcat[i].z_err    = z_choice_err
              qcat[i].zq       = zqual
              qcat[i].zwarning = spz[spzindex].zwarning
              qcat[i].aperture = aperture
              qcat[i].class    = spz[spzindex].class
              qcat[i].comment  = qcomment
              
               printf,lun, obscat[i].id, obscat[i].ra, obscat[i].dec,$
                     obscat[i].mag, spz[spzindex].z, spz[spzindex].z_err,$
                     zqual, aperture, obscat[i].fiber,spz[spzindex].zwarning,$
                     spz[spzindex].class,qcomment,format=das_format

;             printf,lun, obscat[i].id, obscat[i].ra, obscat[i].dec, $
;                     obscat[i].mag, spz[spzindex].z, spz[spzindex].z_err, zqual, $
;                     aperture, obscat[i].fiber,$
;                     spz[spzindex].zwarning,spz[spzindex].class,format=das_format
              loop = 1
              stop = 'stop'
              free_lun,lun
          endif
          
          if strupcase(test) eq 'S' then begin
              splog, '%% Current Smoothing = '+strtrim(string(nsmooth),2)
              read, prompt='Enter new smoothing : ',nsmooth
              loop=1
          endif
;
; Read wavelength and estimate z
;
          if strupcase(test) eq 'E' then begin
              cursor, wl_shifted, yvalue, /down
              read, prompt='Enter rest wavelength : ',wl_rest
              z = (wl_shifted-wl_rest)/wl_rest
              splog,'estimated redshift is '+ string(z)
              show_em_lines, z, ymin, diff
              show_abs_lines, z, ymin, diff
              z_choice        = z
              z_choice_err    = 1.0
              zqual           = 2L
              loop = 2
          endif
              
;
; Kind of find the best z from absorption lines
;        
          if strupcase(test) eq 'A' then begin
              splog, 'Enter the Redshift you would like plotted'
              read, prompt='Redshift : ', z
              show_abs_lines, z, ymin, diff
              z_choice = z
              z_choice_err    = 1.0
              zqual    = 2L
              loop = 2
          endif
;
; Kind of find the best z from emission lines
;        
          
          if strupcase(test) eq 'Z' then begin
              splog, 'Enter the Redshift you would like plotted'
              read, prompt='Redshift : ', z
              show_em_lines, z, ymin, diff
              z_choice = z
              z_choice_err    = 1.0
              zqual    = 2L
              loop = 2
          endif
;
; Plot absorption lines
;        
          if strupcase(test) eq 'ZA' and keyword_set(spzbest) then begin
;           splog, 'Enter the Redshift you would like plotted'
;           read, prompt='Redshift : ', z
;              z = spz[spzindex].z
              z = z_choice
              splog, '%% Z='+strtrim(string(z),2)
              splog, '%% ZERR='+strtrim(string(spz[spzindex].z_err),2)
              splog, '%% ZWARNING='+strtrim(string(spz[spzindex].zwarning),2)
              splog, '%% ZCLASS='+spz[spzindex].class
              show_abs_lines, z, ymin, diff
          endif
;
; Plot emission lines
;        
          if strupcase(test) eq 'ZE' and keyword_set(spzbest) then begin
;           splog, 'Enter the Redshift you would like plotted'
;           read, prompt='Redshift : ', z
;              z = spz[spzindex].z
              z = z_choice
              splog, '%% Using z='+strtrim(string(z),2)
              splog, '%% ZERR='+strtrim(string(spz[spzindex].z_err),2)
              splog, '%% ZWARNING='+strtrim(string(spz[spzindex].zwarning),2)
              splog, '%% ZCLASS='+spz[spzindex].class
              show_em_lines, z, ymin, diff
          endif
;
; plot all galaxy lines
;
          if strupcase(test) eq ']' then begin
;              z = spz[spzindex].z
              z = z_choice
              show_em_lines, z, ymin, diff
              show_abs_lines, z, ymin, diff
          endif
;
; plot sky lines
;
          if strupcase(test) eq '[' then begin
              show_sky_lines, ymin, diff
          endif
          if keyword_set(spzbest) then $
            if i ge n_elements(spz) then stop='stop'
      endwhile
      
  endwhile
  jumpout : print, i
  free_lun,lun
;
; Write out logfile
;
  if (not keyword_set(noupdate)) then $
    save_logfile, logfile, qcat, das_format
;  openw, lun, logfile,/get_lun
;  for i=0, n_elements(qcat)- 1 do begin
;      printf,lun, qcat[i].id, qcat[i].ra, qcat[i].dec, $
;             qcat[i].mag, qcat[i].z, qcat[i].z_err, qcat[i].zq, $
;             qcat[i].aperture, qcat[i].fiber,$
;             qcat[i].zwarning,qcat[i].class,format=das_format
;  endfor
;  free_lun,lun
END
