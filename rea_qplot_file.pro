;------------------------------------------------------------------------------
pro  read_qplot_file, filename, qcat = qcat
;
; read catalogue so that object ids can be tagged to spectra
;
print,'read_qplot_file: ', filename
junk       = 'x'
str = strarr(12)
openr, lun, filename,/get_lun
temp  = {qcat}
temp.id='noid'
temp.class='none'
print, temp
count = 0
WHILE (count < 300) DO BEGIN
    readf, lun, junk
    print, count,' ', junk
    str = strsplit(junk,' ',/extract)
    temp.id         = str[0]
    temp.ra         = float(str[1])
    temp.dec        = float(str[2])
    temp.mag        = float(str[3])
    temp.z          = float(str[4])
    temp.z_err      = float(str[5])
    temp.zq         = fix(str[6])
    temp.aperture   = fix(str[7])
    temp.fiber      = fix(str[8])
    temp.zwarning   = fix(str[9])
    temp.class      = str[10]
    temp.comment    = str[11]
    if (count eq 0 ) then begin
        qcat = temp
    endif else begin
        qcat  = [qcat, temp]
    endelse
    count = count + 1
ENDWHILE
free_lun, lun
end
