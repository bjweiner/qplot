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
temp.comment=''
temp = replicate(temp,300)
;print, temp
count = 0
WHILE count lt 300 DO BEGIN
    readf, lun, junk
;    print, count,' ', junk
    str = strsplit(junk,' ',/extract)
    temp[count].id         = str[0]
    temp[count].ra         = float(str[1])
    temp[count].dec        = float(str[2])
    temp[count].mag        = float(str[3])
    temp[count].z          = float(str[4])
    temp[count].z_err      = float(str[5])
    temp[count].zq         = fix(str[6])
    temp[count].aperture   = fix(str[7])
    temp[count].fiber      = fix(str[8])
    temp[count].zwarning   = fix(str[9])
    temp[count].class      = str[10]
    if n_elements(str) eq 12 then begin
       temp[count].comment    = str[11]
    endif else begin
       temp[count].comment    = ''
    endelse
    count = count + 1
ENDWHILE
qcat = temp
free_lun, lun
end
;
