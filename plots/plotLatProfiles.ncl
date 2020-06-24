;------------------------------------------------------------------------------
; EBM latitudinal temperature profile plot
;
; last edited: 2014-05-23
;------------------------------------------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin

        type = "eps"
        type@wkOrientation = "portrait"
        xwks = gsn_open_wks(type, "plotLatProfiles")

        numplots = 1
	nrow = 1
        ncol = 1
	
;------------------------------------------------------------------------------
; Read data from zonal.out
;------------------------------------------------------------------------------
	temp_f    = 273.16

        filein    = "../model/out/zonal.out"
        nfilecols = numAsciiCol( filein )
        filedata  = readAsciiTable( filein, nfilecols, "float", 0 )
        nfilerows = dimsizes( filedata(:,0) )
        latitude  = filedata(:,0)
        latitude  = filedata(:,0)
	delete( filedata )


	filein    = "../model/data/presentEarth.dat"
	nfilecols = numAsciiCol( filein )
	filedata  = readAsciiTable( filein, nfilecols, "float", 0 )
	nfilerows = dimsizes( filedata(:,0) ) 
	tempPresent = filedata(:,0) - temp_f	

	filein    = "../model/data/bigice.dat"
	nfilecols = numAsciiCol( filein )
	filedata  = readAsciiTable( filein, nfilecols, "float", 0 )
	nfilerows = dimsizes( filedata(:,0) ) 
	tempBigice = filedata(:,0) - temp_f	

	filein    = "../model/data/hothouse.dat"
	nfilecols = numAsciiCol( filein )
	filedata  = readAsciiTable( filein, nfilecols, "float", 0 )
	nfilerows = dimsizes( filedata(:,0) ) 
	tempHothouse = filedata(:,0) - temp_f	


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
        res@trXMinF                             = min( latitude )
        res@trXMaxF                             = max( latitude )
        res@trYReverse                          = False
        res@tmXTOn                              = True
        res@tmXTMinorOn                         = False
        res@tmYLLabelsOn                        = True
        res@tmYROn                              = True
        res@tmYRMinorOn                         = False
	res@xyLineThicknesses			= (/ 3, 3, 3 /)
	res@xyLineColors			= (/ "forestgreen", "dodgerblue", "firebrick" /)
	res@xyDashPatterns			= (/ 0, 0, 0 /)
        res@tmXBMode                            = "Explicit"
        res@tmXBValues                          = (/ -90, -60, -30, 0, 30, 60, 90 /)
        res@tmXBLabels                          = (/ "90S", "60S", "30S", "0", "30N", "60N", "90N" /)


        botlabel   = (/ True /)
	botstring  = (/ "Latitude" /)	
	leftstring = (/ "~F5~D~F~T from freezing [K]" /) 

	do j = 0, numplots - 1

	        res@tmXBLabelsOn                        = botlabel(j)
		res@tiXAxisString                       = botstring(j)
        	res@tiYAxisString                       = leftstring(j)

	        res@trYMinF   = min( tempBigice ) - 5
       		res@trYMaxF   = max( tempHothouse ) + 5
	        res@vpWidthF  = 0.65
       		res@vpHeightF = 0.4

		;res@pmLegendDisplayMode    = "Always"            ; turn on legend
	  	;res@pmLegendSide           = "Top"               ; Change location of 
	  	;res@pmLegendParallelPosF   = .1                  ; move units right
	  	;res@pmLegendOrthogonalPosF = -0.4                ; move units down
	  	;res@pmLegendWidthF         = 0.10                ; Change width and
	  	;res@pmLegendHeightF        = 0.10                ; height of legend.
	  	;res@lgPerimOn              = False               ; turn off box around
	  	;res@lgLabelFontHeightF     = .015                ; label font height
	  	;res@xyExplicitLegendLabels = (/ "T~B~avg", "T~B~min", "T~B~max" /)         ; create explicit labels

		plot(j) = gsn_csm_xy( xwks, latitude, (/ tempPresent, tempBigice, tempHothouse /), res )

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
        resP@gsnPanelLabelBar                   = False
        resP@lbOrientation                      = "Vertical"
        resP@lbTitleOn                          = True
        resP@lbTitleString                      = "K"
        resP@lbTitlePosition                    = "Bottom"
        resP@lbTitleDirection                   = "Across"
        resP@lbTitleExtentF                     = 0.05
        resP@lbTitleFontHeightF                 = 0.012
        resP@lbLabelFontHeightF                 = 0.012
        resP@lbLabelAlignment                   = "BoxCenters"
        resP@lbLabelStride                      = 2
        resP@gsnMaximize                        = True
        resP@gsnFrame                           = False
        resP@amJust                             = "TopLeft"
        resP@gsnPanelFigureStringsFontHeightF   = 0.016
        resP@gsnPaperOrientation                = "portrait"       
        
        gsn_panel( xwks, plot, (/ nrow, ncol /), resP)
        frame( xwks )

end
