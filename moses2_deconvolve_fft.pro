; Deconvolve using Fourier Translation theorem
; If the measured data is represented as
; f(t) = g(t) + b*g(t + a)
; where g(t) is the single image (the target),
; f(t) is the measured image
; b is a scalar representing the relative intensity of the peaks
; a is a vector representing the spatial distance between the peaks

pro moses2_deconvolve_fft

  ;Load MOSES2 level one dataset into memory
  restore, 'IDLWorkspace/MOSES2_LevelOne_Dataset/mosesLevelOne.sav'

  img = cube_zero[*,*,10]
  shift_img = shift(img, 10)
  double_img = 0.5*img + 0.5*shift_img
  
  i1 = image(bytscl(img))
  i2 = image(bytscl(double_img))
  
  convol_orig = convol_fft(img, temp, /AUTO_CORRELATION)
  convol_double = convol_fft(double_img, temp, /AUTO_CORRELATION)
  
  p1 = plot(convol_orig[*,511])
  p1 = plot(convol_double[*,511], /overplot)
  
  p3 = plot(convol_orig[1023,*])
  p3 = plot(convol_double[1023,*], /overplot)
  
  p2 = plot(deriv(deriv(convol_orig[*,511])))
  p4 = plot(deriv(deriv(convol_double[*,511])))
  
  xtv, deriv(deriv(convol_orig))
  xtv, deriv(deriv(convol_double))
  
  p1.close
  
;  output_dir = 'IDLWorkspace/MOSES2_FFTdeconvol'
;  zero_tiff = bytscl(sqrt(zero_autocor))
;  write_tiff, output_dir+'/autocor.tif', zero_tiff

  ;xtv, zero_autocor,screenwidth=1800,screenheight=900

  ;	sz = size(img)
  ;
  ;	Nx = sz[1]
  ;	Ny = sz[2]
  ;
  ;	;ax = -6
  ;	;ay = -3
  ;	;intensity = 0.8
  ;

  ;
  ;	K = k_arr2d(Nx, Ny, /radians, big_kx=kx, big_ky=ky)
  ;
  ;
  ;	; temp
  ;	;i = 16
  ;	J = complex(0,1) ; imaginary unit
  ;
  ;
  ;
  ;	; Loop through data cube to deconvolve each exposure
  ;	;for i = 0, Ndata-1 do begin
  ;
  ;		;xtv, cube_zero[*,*,i], screenwidth=1800, screenheight=900
  ;
  ;		zero_ft = fft(img)
  ;
  ;
  ;
  ;			zero_ft /= (1 + intensity * exp(-J*(kx*ax+ky*ay)))
  ;
  ;
  ;
  ;		img = real_part(fft(zero_ft, /INVERSE))

  ;xtv, img, screenwidth=1800, screenheight=900

  ;endfor

  ;return, img

end
