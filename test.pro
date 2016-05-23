
restore, 'IDLWorkspace/MOSES2_LevelOne_Dataset/mosesLevelOne.sav'

img = cube_zero[1200:2000,320:800,17]

sz = size(img)

sx = sz[1]
sy = sz[2]

;img -= mean(img)
win = hanning(sx,sy)
img *= win

K = k_arr2d(sx, sy, big_kx=kx, big_ky=ky)

k0 = max(K)/10
width = k0/2

high_pass = (erf((k - k0)/width) + 1) / 2


fimgh = high_pass * fft(img)

imgh = fft(fimgh, /inverse)

ac = real_part(fft(abs(fimgh)^2, /inverse))

ac = shift(ac, sx/2,sy/2)

;atv, ac


end