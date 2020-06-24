;------------------------------------------------------------------------------
; EBM bistability plot
;
; last edited: 2014-02-14
;------------------------------------------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin

        type = "oldeps"
        type@wkOrientation = "portrait"
        xwks = gsn_open_wks(type, "plotBistability")

        numplots = 1
	nrow = 1
        ncol = 1
	
;------------------------------------------------------------------------------
; Read data from zonal.out
;------------------------------------------------------------------------------

	; present-Earth start
	filein        = "../model/out/bifurcationWarm.out"
	nfilecols     = numAsciiCol( filein )
	filedata      = readAsciiTable( filein, nfilecols, "float", 0 )
	nfilerows     = dimsizes( filedata(:,0) ) 
	solcon        = filedata(:,0)
	icelineSHwarm = filedata(:,1)
	icelineNHwarm = filedata(:,2)


	; present-Earth backward
	fileinB       = "../model/out/bifurcationBack.out"
	nfilecolsB    = numAsciiCol( fileinB )
	filedataB     = readAsciiTable( fileinB, nfilecolsB, "float", 0 )
	nfilerowsB    = dimsizes( filedataB(:,0) ) 
	solconB       = filedataB(:,0)
	icelineSHback = filedataB(:,1)
	icelineNHback = filedataB(:,2)


	; hothouse start
	fileinW      = "../model/out/bifurcationHot.out"
	nfilecolsW   = numAsciiCol( fileinW )
	filedataW    = readAsciiTable( fileinW, nfilecolsW, "float", 0 )
	nfilerowsW   = dimsizes( filedataW(:,0) ) 
	solconW      = filedataW(:,0)
	icelineSHhot = filedataW(:,1)
	icelineNHhot = filedataW(:,2)


	; cold start
	fileinC       = "../model/out/bifurcationCold.out"
	nfilecolsC    = numAsciiCol( fileinC )
	filedataC     = readAsciiTable( fileinC, nfilecolsC, "float", 0 )
	nfilerowsC    = dimsizes( filedataC(:,0) ) 
	solconC       = filedataC(:,0)
	icelineSHcold = filedataC(:,1)
	icelineNHcold = filedataC(:,2)

;------------------------------------------------------------------------------
; Create contour plot
;------------------------------------------------------------------------------
        plot = new( numplots, graphic )

        res                                     = True
        res@gsnDraw                             = False
        res@gsnFrame                            = False
        res@gsnRightString                      = ""
        res@tiMainString                        = ""
        res@tiXAxisFontHeightF                  = 0.029
        res@tiYAxisFontHeightF                  = 0.029
        res@tmXBLabelFontHeightF                = 0.029
        res@tmYLLabelFontHeightF                = 0.029
        ;res@trXMinF                             = min( latitude )
        ;res@trXMaxF                             = max( latitude )
        res@trYReverse                          = False
        res@tmXTOn                              = True
        res@tmXTMinorOn                         = False
        res@tmYLLabelsOn                        = True
        res@tmYROn                              = True
        res@tmYRMinorOn                         = False
        res@tmXBMajorOutwardLengthF             = 0.0
        res@tmYLMajorOutwardLengthF             = 0.0
        res@tmXBMinorOutwardLengthF             = 0.0
        res@tmYLMinorOutwardLengthF             = 0.0
	res@xyLineThicknesses			= (/ 3, 3, 3 /)
	res@xyLineColors			= (/ "deeppink", "firebrick", "dodgerblue" /)
	res@xyDashPatterns			= (/ 16, 16, 16 /)
        ;res@tmXBMode                            = "Explicit"
        ;res@tmXBValues                          = (/ -90, -60, -30, 0, 30, 60, 90 /)
        ;res@tmXBLabels                          = (/ "90S", "60S", "30S", "0", "30N", "60N", "90N" /)
        res@xyMarkLineMode                      = "MarkLines"
        res@xyMarkers                           = (/ 16, 16, 16 /)
        res@xyMarkerColors                      = (/ "deeppink", "firebrick", "dodgerblue" /)

        botlabel   = (/ True /)
	botstring  = (/ "Relative solar flux (S/S~B~0~N~)" /)	
	leftstring = (/ "Ice line latitude (~F35~J~F~)" /) 

	do j = 0, numplots - 1

	        res@tmXBLabelsOn                        = botlabel(j)
		res@tiXAxisString                       = botstring(j)
        	res@tiYAxisString                       = leftstring(j)

	        res@trXMinF   = 0.75
        	res@trXMaxF   = 1.50
	        res@trYMinF   = -9.0
        	res@trYMaxF   = 99.0
	        res@vpWidthF  = 0.6
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

		plot(j) = gsn_csm_xy( xwks, solcon, (/ icelineNHhot, icelineNHwarm, icelineNHcold(::-1) /), res )

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
