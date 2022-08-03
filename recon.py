import numpy as np
import cfl
import os
import sys
import time
import scipy.io as sio
from bart import bart
from sys import argv
os.putenv("DEBUG_LEVEL", "5")

# read data
kspace = cfl.readcfl('{}'.format(argv[1]))

# coil compression first, then zero-padding
necho = kspace.shape[-1]
ncoil = 8
nslice = kspace.shape[0]
nrow = kspace.shape[1]
ncol = kspace.shape[2]
kspace_cc = np.zeros((nslice, nrow, ncol, ncoil, necho), dtype=np.complex_)
sens_cc = np.zeros((nslice, nrow, ncol, ncoil, necho), dtype=np.complex_)
for i in range(necho):
    print(i)
    kspace_cc[..., i] = bart(1, 'cc -p 8 -S -r 20', kspace[..., i])
del kspace
print(kspace_cc.shape)

# in ecalib, -c0.8 is the default value of eigenvalues truncation
necho = kspace_cc.shape[-1]
ncoil = 8
nslice = kspace_cc.shape[0]
nrow = kspace_cc.shape[1]
ncol = kspace_cc.shape[2]
sens_cc = np.zeros((nslice, nrow, ncol, ncoil, necho), dtype=np.complex_)
for i in range(necho):
    print(i)
    if i < necho-2:
        sens_cc[..., i] = bart(1, 'ecalib -m1 -c{} -r 20'.format(0.8 - 0.07*i), kspace_cc[..., i])  # sensmap for each echo, 0.8 - 0.07*i
    else:
        sens_cc[..., i] = bart(1, 'ecalib -m1 -c{} -r 20'.format(0.8), kspace_cc[..., i])  # sensmap for each echo, 0.8 - 0.07*i

print(np.mean(abs(kspace_cc[160, :, :, 0, 3]) > 0))
kspace = np.transpose(kspace_cc, (1, 2, 0, 3, 4))
del kspace_cc

# ifft recon
kspace_fft_z = np.fft.fftshift(np.fft.ifft(np.fft.fftshift(kspace, axes=2), axis=2), axes=2)
del kspace
# save data
necho = kspace_fft_z.shape[4]
ncoil = kspace_fft_z.shape[3]
nslice = kspace_fft_z.shape[2]
ncol = kspace_fft_z.shape[1]
nrow = kspace_fft_z.shape[0]

flip_matrix = np.ones((nrow, ncol, 1, ncoil), dtype=np.complex_)
flip_matrix[::2, ...] = -flip_matrix[::2, ...] 
flip_matrix[:, ::2, ...] = -flip_matrix[:, ::2, ...]

iField = np.zeros((nslice, nrow, ncol, necho), dtype=np.complex_)

for i in range(nslice):
    print('Processing slice # {}'.format(i))
    fully_slice = np.zeros((nrow, ncol, necho), dtype=np.complex_)
    cs_slice = np.zeros((nrow, ncol, necho), dtype=np.complex_)
    sensMaps_slice = np.zeros((nrow, ncol, ncoil), dtype=np.complex_)
    kdata_slice = np.zeros((nrow, ncol, ncoil, necho), dtype=np.complex_)
    kspace_slice = np.zeros((nrow, ncol, 1, ncoil), dtype=np.complex_)
    for j in range(necho):
        kspace_slice = kspace_fft_z[:, :, i:i+1, :, j]

        sensMaps_slice = np.concatenate((sens_cc[i, :, ncol//2:, :, j], sens_cc[i, :, :ncol//2, :, j]), axis=1)
        
        sensMaps_slice = sensMaps_slice * np.exp(-1j * np.angle(sensMaps_slice[..., 0:1]))

        recon_full = np.sum(np.fft.fftshift(np.fft.ifft2(kspace_slice, axes=(0, 1)), axes=0) * \
                            flip_matrix * np.conj(sensMaps_slice[:, :, np.newaxis, :]), axis=3)

        fully_slice[:, :, j] = recon_full[..., 0]
        kdata_slice[:, :, :, j] = kspace_fft_z[:, :, i, :, j]
                
    iField[i, ...] = fully_slice

# save iField
adict = {}
adict['iField'] = iField
sio.savemat('{}_recon.mat'.format(argv[1]), adict)