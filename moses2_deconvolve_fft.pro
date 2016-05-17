; Deconvolve using Fourier Translation theorem
; If the measured data is represented as
; f(t) = g(t) + b*g(t + a)
; where g(t) is the single image (the target),
; f(t) is the measured image
; b is a scalar representing the relative intensity of the peaks
; a is a vector representing the spatial distance between the peaks

;pro moses2_deconvolve_fft

;Load MOSES2 level one dataset into memory
restore, 'IDLWorkspace/MOSES2_LevelOne_Dataset/mosesDarkSub.sav'

; Define images needed for computation
oimg = cube_zero[*,*,17]
oshift_img = shift(oimg, 20,20)
odouble_img = 0.5*oimg + 0.5*oshift_img
atrous_tr, oimg, 20, decomp=oimgd

; Compute the autocorrelation
oconvol_orig = convol_fft(oimg, temp, /AUTO_CORRELATION)
oconvol_double = convol_fft(odouble_img, temp, /AUTO_CORRELATION)
oconvol_atrous = convol_fft(oimgd[*,*,6], temp, /auto_correlation)

; Compute the Lapalacian of the autocorrelation
ograd_orig = gradient(gradient(oconvol_orig))
ograd_double = gradient(gradient(oconvol_double))
ograd_atrous = gradient(gradient(oconvol_atrous))

ax = 7
ay = 1
intensity = 0.3

J = complex(0,1) ; imaginary unit
K = k_arr2d(Nx, Ny, /radians, big_kx=kx, big_ky=ky)
output_dir = 'IDLWorkspace/MOSES2_LevelOne_Dataset/MOSES2_tiff_images_level1'

; Loop through data cube to deconvolve each exposure
for i = 0, Ndata-1 do begin

  d = data_list[i] ; find index of next data image

  zero_ft = fft(cube_zero[*,*,i])
  zero_ft /= (1 + intensity * exp(-J*(kx*ax+ky*ay)))
  cube_zero[*,*,i] = real_part(fft(zero_ft, /INVERSE))

endfor
;
;ax = 0
;ay = 9
;intensity = 0.10
;
;; Loop again through data cube to deconvolve each exposure
;for i = 0, Ndata-1 do begin
;
;  d = data_list[i] ; find index of next data image
;
;  zero_ft = fft(cube_zero[*,*,i])
;  zero_ft /= (1 + intensity * exp(-J*(kx*ax+ky*ay)))
;  cube_zero[*,*,i] = real_part(fft(zero_ft, /INVERSE))
;
;  zero_tiff = bytscl(sqrt(cube_zero[*,*,i]))
;  write_tiff, output_dir+'/zero/'+strmid(index.filename[d],7,12)+string((index.exptime[d]*1e-6), FORMAT='(f20.2)')+'  corrected.tif', zero_tiff
;
;endfor

; Save images to disk
for i = 0, Ndata-1 do begin
  
    zero_tiff = bytscl(sqrt(cube_zero[*,*,i]))
    write_tiff, output_dir+'/zero/'+strmid(index.filename[d],7,12)+string((index.exptime[d]*1e-6), FORMAT='(f20.2)')+'  corrected.tif', zero_tiff
  
endfor

; Define images needed for computation
img = cube_zero[*,*,17]
shift_img = shift(img, 20,20)
double_img = 0.5*img + 0.5*shift_img
atrous_tr, img, 20, decomp=imgd

; Compute the autocorrelation
convol_orig = convol_fft(img, temp, /AUTO_CORRELATION)
convol_double = convol_fft(double_img, temp, /AUTO_CORRELATION)
convol_atrous = convol_fft(imgd[*,*,6], temp, /auto_correlation)

; Compute the Lapalacian of the autocorrelation
grad_orig = gradient(gradient(convol_orig))
grad_double = gradient(gradient(convol_double))
grad_atrous = gradient(gradient(convol_atrous))


print, systime()+' MOSES FFT deconvolve saving to disk.'
save, filename='IDLWorkspace/MOSES2_LevelOne_Dataset/mosesLevelOne.sav'
print, systime()+' MOSES FFT deconvolve completed.'

end
