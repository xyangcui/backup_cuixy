
wks = gsn_open_wks("x11", "text")

text = "~F21~-(~F8~w~E~~F8~s~F21~)"
text2 = "~F21~-(v~F34~6~F21~T~F21~/~F34~6~F21~y~F21~)'"
text3 = "~F34~'"
res               = True                      ; text mods desired
;res@txFontHeightF = 0.01                     ; text font height
;res@txJust        = "Center"              ; Default is "CenterCenter".

gsn_text_ndc(wks, text3, 0.5, 0.5, res)

i = 1
do while (i.eq.1)
    i = 1
end do
