;------------------------------------------------------------------------------
; EBM daily temperature series plot for all latitude
;
; last edited: 2019-08-29
;------------------------------------------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin

        type = "oldeps"
        type@wkOrientation = "portrait"
        xwks = gsn_open_wks(type, "plotDailyLatTemp")

        numplots = 1
	nrow = 1
        ncol = 1
	
;------------------------------------------------------------------------------
; Read data from zonal.out
;------------------------------------------------------------------------------
        ;daystart = 2999
        ;dayend   = 4501

	; Zonal temperature data
	filein    = "../model/out/daily.out"
	nfilecols = numAsciiCol( filein )
	filedata  = readAsciiTable( filein, nfilecols, "float", 0 )
	nfilerows = dimsizes( filedata(:,0) ) 

        temp = new( (/ nfilecols-1, nfilerows /), float )
	lat = (/ -85.0, -75.0, -65.0, -55.0, -45.0, -35.0, -25.0, -15.0, -5.0, 5.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0, 85.0 /)            

	time = filedata(:,0)
        do i = 1, nfilecols-1
	  temp(i-1,:) = filedata(:,i)
        end do

        temp!0		= "lat"
	temp!1		= "time"
	temp&lat	= lat
	temp&time	= time

;------------------------------------------------------------------------------
; Create contour plot
;------------------------------------------------------------------------------
        plot = new( numplots, graphic )

        res                                     = True
        res@gsnDraw                             = False
        res@gsnFrame                            = False
        res@gsnRightString                      = ""
        res@tiMainString                        = ""
        res@tiXAxisFontHeightF                  = 0.021
        res@tiYAxisFontHeightF                  = 0.021
        ;res@trXMinF                             = min( latitude )
        ;res@trXMaxF                             = max( latitude )
        res@trYReverse                          = False
        res@tmXTOn                              = True
        res@tmXTMinorOn                         = False
        res@tmYLLabelsOn                        = True
        res@tmYROn                              = True
        res@tmYRMinorOn                         = False
	;res@xyLineThicknesses			= (/ 3, 3, 3 /)
	;res@xyLineColors			= (/ "forestgreen", "dodgerblue", "firebrick" /)
	;res@xyDashPatterns			= (/ 0, 0, 0 /)
        ;res@tmXBMode                            = "Explicit"
        ;res@tmXBValues                          = (/ -90, -60, -30, 0, 30, 60, 90 /)
        ;res@tmXBLabels                          = (/ "90S", "60S", "30S", "0", "30N", "60N", "90N" /)


	;gsn_define_colormap( xwks, "helix" )
	;gsn_define_colormap( xwks, "example" )
	gsn_define_colormap( xwks, "WhiteBlueGreenYellowRed" )

        res@cnFillOn                            = True
        ;res@cnLevelSelectionMode                = "ManualLevels"
        ;res@cnMinLevelValF                      = 0
        ;res@cnMaxLevelValF                      = 3000.0
        ;res@cnLevelSpacingF                     = 20.0
        ;res@cnLineLabelBackgroundColor          = -1
        res@cnLinesOn                           = False
        res@cnInfoLabelOn                       = False
        res@cnLineLabelInterval                 = 0

        res@gsnSpreadColors                     = True
        res@gsnSpreadColorStart                 = 2
        res@gsnSpreadColorEnd                   = 255
	res@lbLabelBarOn                        = False


        botlabel   = (/ True /)
	botstring  = (/ "Day" /)	
	leftstring = (/ "Latitude" /) 

	do j = 0, numplots - 1

	        res@tmXBLabelsOn                        = botlabel(j)
		res@tiXAxisString                       = botstring(j)
        	res@tiYAxisString                       = leftstring(j)

	        ;res@trYMinF   = 270.0
        	;res@trYMaxF   = 310.0
	        res@vpWidthF  = 0.9
        	res@vpHeightF = 0.3

		;res@pmLegendDisplayMode    = "Always"            ; turn on legend
	  	;res@pmLegendSide           = "Top"               ; Change location of 
	  	;res@pmLegendParallelPosF   = .1                  ; move units right
	  	;res@pmLegendOrthogonalPosF = -0.4                ; move units down
	  	;res@pmLegendWidthF         = 0.10                ; Change width and
	  	;res@pmLegendHeightF        = 0.10                ; height of legend.
	  	;res@lgPerimOn              = False               ; turn off box around
	  	;res@lgLabelFontHeightF     = .015                ; label font height
	  	;res@xyExplicitLegendLabels = (/ "T~B~avg", "T~B~min", "T~B~max" /)         ; create explicit labels

		plot(j) = gsn_csm_contour( xwks, temp, res )
		;plot(j) = gsn_csm_contour( xwks, temp(:,daystart:dayend), res )

	end do

;------------------------------------------------------------------------------
; Plot Figures
;------------------------------------------------------------------------------
                             
        resP                                    = True
        resP@txString                           = ""
        resP@lbLabelFontHeightF                 = 0.01
        resP@gsnPanelBottom                     = 0.5
        ;resP@gsnPanelXF                         = (/ -1, -1 /)
        ;resP@gsnPanelYF                         = (/ -1, -1 /)
        resP@gsnPanelLabelBar                   = True
        resP@cnLabelBarEndStyle                 = "ExcludeOuterBoxes"
        resP@lbBoxLinesOn                       = False
        ;resP@lbLabelAlignment                  = "ExternalEdges"
        resP@lbOrientation                      = "Vertical"
        resP@lbTitleOn                          = True
        resP@lbTitleString                      = "K"
        resP@lbTitlePosition                    = "Top"
        resP@lbTitleDirection                   = "Across"
        resP@lbTitleExtentF                     = 0.05
        resP@lbTitleFontHeightF                 = 0.012
        resP@lbLabelFontHeightF                 = 0.012
        ;resP@lbLabelAlignment                   = "BoxCenters"
        resP@lbLabelStride                      = 3
        resP@lbPerimOn                          = False
        resP@pmLabelBarWidthF                   = 0.11
        resP@pmLabelBarHeightF                  = 0.32
        resP@gsnMaximize                        = True
        resP@gsnFrame                           = False
        resP@amJust                             = "TopLeft"
        resP@gsnPanelFigureStringsFontHeightF   = 0.016
        resP@gsnPaperOrientation                = "portrait"       
        
        gsn_panel( xwks, plot, (/ nrow, ncol /), resP)
        frame( xwks )

end
