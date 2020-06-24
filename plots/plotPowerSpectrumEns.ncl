;------------------------------------------------------------------------------
; EBM temperature series power spectrum ensemble plot
;
; last edited: 2014-02-21
;------------------------------------------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

        type = "eps"
        type@wkOrientation = "portrait"
        xwks = gsn_open_wks(type, "plotPowerSpectrumEns")

        numplots = 2
	nrow = 1
        ncol = 2
	
;------------------------------------------------------------------------------
; Read data from zonal.out
;------------------------------------------------------------------------------
        numdata  = 101
        numtime  = 2001
        timehalf = 1000
        timequarter = 1500
        timeadj  = 1000000.
        tempadj  = 288.

        temps   = new( (/ numdata, numtime /), float )
        temps1x = new( (/ numdata, numtime /), float )
        temps4x = new( (/ numdata, numtime /), float )
        temps1xdd = new( (/ numdata, numtime /), float )
        temps4xdd = new( (/ numdata, numtime /), float )

        do k = 1, numdata - 1
                ;print( "Reading..." + "tempseries_" + k + ".out" )
                filein     = "../model/out/ens100redo/tempseries_" + k + ".out"
                nfilecols  = numAsciiCol( filein )
                filedata   = readAsciiTable( filein, nfilecols, "float", 0 )
                nfilerows  = dimsizes( filedata(:,0) )
                temps(k,:) = filedata(:,1) - tempadj

                filein1x     = "../model/out/ens1x/tempseries_" + k + ".out"
                nfilecols1x  = numAsciiCol( filein1x )
                filedata1x   = readAsciiTable( filein1x, nfilecols1x, "float", 0 )
                nfilerows1x  = dimsizes( filedata1x(:,0) )
                temps1x(k,:) = filedata1x(:,1) - tempadj

                filein1xdd   = "../model/out/ens1xdd/tempseries_" + k + ".out"
                nfilecols    = numAsciiCol( filein1xdd )
                filedata     = readAsciiTable( filein1xdd, nfilecols, "float", 0 )
                nfilerows    = dimsizes( filedata(:,0) )
                temps1xdd(k,:) = filedata(:,1) - tempadj

                filein4x     = "../model/out/ens4x/tempseries_" + k + ".out"
                nfilecols4x  = numAsciiCol( filein4x )
                filedata4x   = readAsciiTable( filein4x, nfilecols4x, "float", 0 )
                nfilerows4x  = dimsizes( filedata4x(:,0) )
                temps4x(k,:) = filedata4x(:,1) - tempadj

                filein4xdd   = "../model/out/ens4xdd/tempseries_" + k + ".out"
                nfilecols    = numAsciiCol( filein4xdd )
                filedata     = readAsciiTable( filein4xdd, nfilecols, "float", 0 )
                nfilerows    = dimsizes( filedata(:,0) )
                temps4xdd(k,:) = filedata(:,1) - tempadj
        end do


;------------------------------------------------------------------------------
; De-Trend draw-down scenarios
;------------------------------------------------------------------------------

	tave1xdd = dim_avg_n_Wrap( temps1xdd, 0 )
	t1 = dtrend( tave1xdd( timehalf:timequarter-1 ), False )
	t2 = tave1xdd( timequarter:numtime-1 )
	z1 = array_append_record( t1, t2, 0 )
	temps1xddD = array_append_record( tave1xdd( 0:timehalf-1 ), z1, 0 )

	tave4xdd = dim_avg_n_Wrap( temps4xdd, 0 )
	t1 = dtrend( tave4xdd( timehalf:timequarter-1 ), False )
	t2 = tave4xdd( timequarter:numtime-1 )
	z1 = array_append_record( t1, t2, 0 )
	temps4xddD = array_append_record( tave4xdd( 0:timehalf-1 ), z1, 0 )


;------------------------------------------------------------------------------
; Create plot
;------------------------------------------------------------------------------
        plot = new( numplots, graphic )

        res                                     = True
        res@gsnDraw                             = False
        res@gsnFrame                            = False
        res@gsnRightString                      = ""
        res@tiMainString                        = ""
        res@tiXAxisFontHeightF                  = 0.029
        res@tiYAxisFontHeightF                  = 0.029
	res@tmXBLabelFontHeightF		= 0.029
	res@tmYLLabelFontHeightF		= 0.029
        ;res@trXMinF                             = min( latitude )
        ;res@trXMaxF                             = max( latitude )
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
        ;res@tmXBMode                            = "Explicit"
        ;res@tmXBValues                          = (/ -90, -60, -30, 0, 30, 60, 90 /)
        ;res@tmXBLabels                          = (/ "90S", "60S", "30S", "0", "30N", "60N", "90N" /)
	res@trXLog                              = True 
	res@trYLog                              = False


        botlabel   = (/ True, True /)
	botstring  = (/ "Period (yr)", "Period (yr)" /)	
	leftstring = (/ "Power", "" /) 

	do j = 0, numplots - 1

	        res@tmXBLabelsOn                        = botlabel(j)
		res@tiXAxisString                       = botstring(j)
        	res@tiYAxisString                       = leftstring(j)

	        res@trXMinF   = 10000
        	res@trXMaxF   = 1000000
	        res@trYMinF   = 0.0
        	res@trYMaxF   = 1.1
	        res@vpWidthF  = 0.6
        	res@vpHeightF = 0.4	

		;res@pmLegendDisplayMode    = "Always"            ; turn on legend
	  	;res@pmLegendSide           = "Top"               ; Change location of 
	  	;res@pmLegendParallelPosF   = .15                  ; move units right
	  	;res@pmLegendOrthogonalPosF = -0.4                ; move units down
	  	;res@pmLegendWidthF         = 0.10                ; Change width and
	  	;res@pmLegendHeightF        = 0.10                ; height of legend.
	  	;res@lgPerimOn              = False               ; turn off box around
	  	;res@lgLabelFontHeightF     = .015                ; label font height
	  	;res@xyExplicitLegendLabels = (/ "3.7 W/m~S~2", "14.8 W/m~S~2" /)         ; create explicit labels


		;------------------------------------------------------------------------------
		; Calculate power spectrum
		;------------------------------------------------------------------------------
		nfreq      = 1000.

		if ( j .eq. 0 ) then
        		specout    = specx_anal( dim_avg_n_Wrap( temps, 0 ), 0, 3, 0.1 )
	        	specout1x  = specx_anal( dim_avg_n_Wrap( temps( 1:10, : ), 0 ), 0, 3, 0.1 )


                        spectralall = new( (/ numdata, dimsizes( specout@spcx ) /), float )
			do r = 1, numdata - 1
	        		specout        = specx_anal( temps(r,:), 0, 3, 0.1 )
				spectralall(r,:) = specout@spcx               
			end do		
	


	        	specout4x  = specx_anal( dim_avg_n_Wrap( temps4x( :, timehalf:numtime-1 ), 0 ), 0, 3, 0.1 )
        		specout1xdd  = specx_anal( dim_avg_n_Wrap( temps1xdd( :, timehalf:numtime-1 ), 0 ), 0, 3, 0.1 )
	        	specout4xdd  = specx_anal( dim_avg_n_Wrap( temps4xdd( :, 0:timehalf-1 ), 0 ), 0, 3, 0.1 )
		end if
		if ( j .eq. 1 ) then 
		        res@tmYLLabelsOn                        = False

        		specout    = specx_anal( dim_avg_n_Wrap( temps, 0 ), 0, 3, 0.1 )
	        	specout1x  = specx_anal( dim_avg_n_Wrap( temps1x( :, timehalf:numtime-1 ), 0 ), 0, 3, 0.1 )
	        	specout4x  = specx_anal( dim_avg_n_Wrap( temps4x( :, timehalf:numtime-1 ), 0 ), 0, 3, 0.1 )
        		;specout1xdd  = specx_anal( dim_avg_n_Wrap( temps1xddD( :, timequarter:numtime-1 ), 0 ), 0, 3, 0.1 )
	        	;specout4xdd  = specx_anal( dim_avg_n_Wrap( temps4xddD( :, timequarter:numtime-1 ), 0 ), 0, 3, 0.1 )



			; uncomment for detrended temperature series
        		;specout1xdd  = specx_anal( temps1xddD( timehalf:numtime-1 ), 0, 3, 0.1 )
	        	;specout4xdd  = specx_anal( temps4xddD( timehalf:numtime-1 ), 0, 3, 0.1 )
        		specout1xdd  = specx_anal( temps1xddD( timehalf:timequarter-1 ), 0, 3, 0.1 )
	        	specout4xdd  = specx_anal( temps4xddD( timehalf:timequarter-1 ), 0, 3, 0.1 )


		        ;specout  = specx_anal( dim_avg_n_Wrap( temps( :, timequarter:numtime-1 ), 0 ), 0, 3, 0.1 )
        		;specout1x  = specx_anal( dim_avg_n_Wrap( temps1x( :, timequarter:numtime-1 ), 0 ), 0, 3, 0.1 )
        		;specout4x  = specx_anal( dim_avg_n_Wrap( temps4x( :, timequarter:numtime-1 ), 0 ), 0, 3, 0.1 )
        		;specout1xdd  = specx_anal( dim_avg_n_Wrap( temps1xdd( :, timequarter:numtime-1 ), 0 ), 0, 3, 0.1 )
        		;specout4xdd  = specx_anal( dim_avg_n_Wrap( temps4xdd( :, timequarter:numtime-1 ), 0 ), 0, 3, 0.1 )
		end if
	

	        ;specout  = specx_anal( dim_avg_n_Wrap( temps( :, timehalf:timequarter-1 ), 0 ), 0, 3, 0.1 )
	        ;specout1x  = specx_anal( dim_avg_n_Wrap( temps1x( :, timehalf:timequarter-1 ), 0 ), 0, 3, 0.1 )
	        ;specout4x  = specx_anal( dim_avg_n_Wrap( temps4x( :, timehalf:timequarter-1 ), 0 ), 0, 3, 0.1 )
	        ;specout1xdd  = specx_anal( dim_avg_n_Wrap( temps1xdd( :, timehalf:timequarter-1 ), 0 ), 0, 3, 0.1 )
	        ;specout4xdd  = specx_anal( dim_avg_n_Wrap( temps4xdd( :, timehalf:timequarter-1 ), 0 ), 0, 3, 0.1 )


	       	freq       = specout@frq
        	period     = nfreq / freq
        	spectral   = specout@spcx               
        	;spectral1x = specout1x@spcx * 2               
        	;spectral4x = specout4x@spcx * 2           
        	spectral1x = specout1x@spcx
        	spectral4x = specout4x@spcx
        	spectral1xdd = specout1xdd@spcx
        	spectral4xdd = specout4xdd@spcx


		if ( j .eq. 0 ) then

                        tmax = new( (/ dimsizes( spectral ) /), float )
                        tmin = new( (/ dimsizes( spectral ) /), float )
                        do q = 0, dimsizes( spectral ) - 1
                                tmax(q) = max( spectralall(:,q) )
                                tmin(q) = min( spectralall(:,q) )
                        end do


                        res@xyLineThicknesses                   = 1
                        res@xyDashPatterns                      = 0
                        res@gsnXYFillColors                     = "grey90"
                        res@xyLineColors                        = (/ -1, -1 /)
			plot(j) = gsn_csm_xy( xwks, period, (/ tmax / max( tmax ), tmin / max( tmax ) /), res )
                        delete( res@gsnXYFillColors )
			delete( res@xyLineThicknesses )
			delete( res@xyLineColors )
			delete( res@xyDashPatterns )

			res@xyLineThicknesses			= (/ 2, 3 /)
			res@xyLineColors			= (/ "lightskyblue2", "black" /)
			res@xyDashPatterns			= (/ 0, 0 /)
			plotControl = gsn_csm_xy( xwks, period, (/ ( spectral1x * 10 ) / max( tmax ), spectral / max( tmax ) /), res )
			overlay( plot(j), plotControl )

		else
			res@xyLineThicknesses			= (/ 2, 2 /)
			res@xyLineColors			= (/ "firebrick", "darkgoldenrod2" /)
			res@xyDashPatterns			= (/ 0, 0 /)
			plot(j) = gsn_csm_xy( xwks, period, (/ spectral4x / max( spectral4x ), spectral1x / max( spectral1x ) /), res )

			res@xyLineThicknesses			= (/ 2, 2 /)
			res@xyLineColors			= (/ "chocolate4", "darkolivegreen4" /)
			res@xyDashPatterns			= (/ 0, 0 /)
			plotDD = gsn_csm_xy( xwks, period, (/  spectral4xdd / max( spectral4xdd ), spectral1xdd / max( spectral1xdd ) /), res )
			overlay( plot(j), plotDD )
		end if

		delete( res@xyLineThicknesses )
		delete( res@xyLineColors )
		delete( res@xyDashPatterns )
		;delete( res@xyExplicitLegendLabels )

		delete( freq )
		delete( period )
		delete( spectral )
		delete( spectral1x )
		delete( spectral4x )
		delete( spectral1xdd )
		delete( spectral4xdd )
		delete( res@xyLineThicknesses )
		delete( res@xyLineColors )
		delete( res@xyDashPatterns )

	end do

;------------------------------------------------------------------------------
; Plot Figures
;------------------------------------------------------------------------------
                             
        resP                                    = True
        resP@txString                           = ""
        resP@lbLabelFontHeightF                 = 0.01
        resP@gsnPanelBottom                     = 0.5
        ;resP@gsnPanelYWhiteSpacePercent         = 0.1
        ;resP@gsnPanelXWhiteSpacePercent         = 0.1
        resP@gsnPanelXF                         = (/ -1, 0.56 /)
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
	resP@gsnPanelFigureStringsPerimOn 	= False
	resP@gsnPanelFigureStrings              = (/ "", "" /)
        resP@gsnPanelFigureStringsFontHeightF   = 0.016
        resP@gsnPaperOrientation                = "portrait"       
        
        gsn_panel( xwks, plot, (/ nrow, ncol /), resP)
        frame( xwks )

end
