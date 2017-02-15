; define catalogue structure for spec1d pipeline.
; function qcat::Init, $
; noid,$
; 0.0,$
; 0.0,$
;  0.0,$
;  0.0,$
;  0.0,$
;  -1, $
;  0,$
;  0,$
;  0,$
;  none
;end
;
pro qcat__define
tmp  = {qcat          , $
          id:     'noid',$
          ra:       0.0,$
          dec:      0.0,$
          mag:      0.0,$
          z:        0.0,$
          z_err:    0.0,$
          zq:       -1, $
          aperture:   0,$
          fiber:      0,$
          zwarning:   0,$
          class:  'none',$
          comment:  ''}
end

