;------------------------------------------------------------------------------
; EBM stochastic parameter histogram plot
;
; last edited: 2014-05-22
;------------------------------------------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

        type = "eps"
        type@wkOrientation = "portrait"
        xwks = gsn_open_wks(type, "plotStochHistogram")

        numplots = 1
	nrow = 1
        ncol = 1
	
;------------------------------------------------------------------------------
; Read data from /models/stochastic_ebm/model/out/noise.txt
;------------------------------------------------------------------------------
        numdata  = 100

        params    = new( (/ numdata /), float )
        filein    = "../model/out/noise.txt"
        nfilecols = numAsciiCol( filein )
        filedata  = readAsciiTable( filein, nfilecols, "float", 0 )
	nfilerows = dimsizes( filedata(:,0) )
        params    = filedata(:,0)

	printVarSummary( params )
	print( avg( params ) )
	print( stddev( params ) )
	print( variance( params ) )


;------------------------------------------------------------------------------
; Create plot
;------------------------------------------------------------------------------
        plot = new( numplots, graphic )

        res                                     = True
        res@gsnDraw                             = False
        res@gsnFrame                            = False
        res@tiMainString                        = ""
        res@tiXAxisFontHeightF                  = 0.029
        res@tiYAxisFontHeightF                  = 0.029
	res@tmXBLabelFontHeightF		= 0.029
	res@tmYLLabelFontHeightF		= 0.029
        res@trYReverse                          = False
        res@tmXTOn                              = True
        res@tmXTMinorOn                         = False
        res@tmYLLabelsOn                        = True
        res@tmYROn                              = True
        res@tmYRMinorOn                         = False
	res@tmXBMajorOutwardLengthF 		= 0.0
	res@tmYLMajorOutwardLengthF 		= 0.0
	res@tmXBMinorOutwardLengthF 		= 0.0 
	res@tmYLMinorOutwardLengthF 		= 0.0 

	;res@gsnHistogramNumberOfBins 		= 10
	;res@gsnHistogramClassIntervals 		= (/ 0.10, 0.125, 0.15, 0.175, 0.20, 0.225, 0.25, 0.275, 0.30, 0.325, 0.35, 0.375, 0.40, 0.425, 0.45, 0.475, 0.50, 0.525, 0.55 /)
	res@gsnHistogramClassIntervals 		= (/ 0.10 , 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55 /)
	res@tmXBLabelStride          		= 2
	res@gsnHistogramComputePercentages 	= False
	res@gsnHistogramPercentSign		= False
	res@gsFillColor				= "firebrick"


        botlabel   = (/ True /)
	botstring  = (/ "Noise parameter ~F5~e~F~ (W~S~2~N~ m~S~-4~N~)" /)	
	leftstring = (/ "Frequency" /) 

	do j = 0, numplots - 1

	        res@tmXBLabelsOn                        = botlabel(j)
		res@tiXAxisString                       = botstring(j)
        	;res@tiYAxisString                       = leftstring(j)

	        res@vpWidthF  = 0.6
        	res@vpHeightF = 0.4	

		plot(j) = gsn_histogram( xwks, params, res )

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
