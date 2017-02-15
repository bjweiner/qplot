pro  recover_cat, filename, obscat = obscat
;
; read catalogue so that object ids can be tagged to spectra
;
junk       = 'x'
openr, lun, filename,/get_lun
temp  = {obscat}
count = 0
WHILE (NOT EOF(lun)) DO BEGIN
    readf, lun, junk
    print, junk
    str = strsplit(junk,' ',/extract)
    temp.aperture   = fix(str[0])
    temp.ra         = float(str[2])
    temp.dec        = float(str[3])
    temp.mag        = float(str[4])
    temp.fiber      = fix(str[6])
    temp.beam       = fix(str[7])
    temp.id         = str[8]
    if (count eq 0 ) then begin
        obscat = temp
    endif else begin
        obscat  = [obscat, temp]
    endelse
    count = count + 1
ENDWHILE
free_lun, lun
end
