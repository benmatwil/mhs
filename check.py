import matplotlib.pyplot as plt
import numpy as np
import read as rd
# import init

# plt.ioff()

field = rd.field('data/mhs_field_20100601_0081_fft_alpha_1.00_d_0.00-finite.dat')
potfield = rd.field('../msat/data/field_20100601_0081_fft.dat')
synmap = rd.synmap('hmi/synmap_20100601.dat')

minmax = 20
plt.figure(1)
plt.contourf(field[5], field[4], field[0][0,:,:], np.linspace(-minmax,minmax,51), extend='both', cmap=plt.cm.RdBu_r)
plt.colorbar()
plt.title(r'MHS $B_r$ on the base')

minmax = 0.5
plt.figure(2)
plt.contourf(field[5], field[4], field[0][-1,:,:], np.linspace(-minmax,minmax,51), extend='both', cmap=plt.cm.RdBu_r)
plt.colorbar()
plt.contour(field[5], field[4], field[0][-1,:,:], [0])
plt.ylim(np.pi, 0)
plt.title(r'MHS $B_r$ on the top')

plt.figure(3)
plt.contourf(field[5], field[4], field[1][-1,:,:], np.linspace(-minmax,minmax,51), extend='both', cmap=plt.cm.RdBu_r)
plt.colorbar()
plt.contour(field[5], field[4], field[1][-1,:,:], [0])
plt.ylim(np.pi, 0)
plt.title(r'MHS $B_{\theta}$ on the top')

plt.figure(4)
plt.contourf(field[5], field[4], field[2][-1,:,:], np.linspace(-minmax,minmax,51), extend='both', cmap=plt.cm.RdBu_r)
plt.colorbar()
plt.contour(field[5], field[4], field[2][-1,:,:], [0])
plt.ylim(np.pi, 0)
plt.title(r'MHS $B_\phi$ on the top')

minmax = 20
plt.figure(5)
plt.contourf(potfield[5], potfield[4], potfield[0][0,:,:], np.linspace(-minmax,minmax,51), extend='both', cmap=plt.cm.RdBu_r)
plt.colorbar()
plt.title(r'Potential $B_r$ on the base')

minmax = 0.001
plt.figure(6)
plt.contourf(potfield[5], potfield[4], potfield[0][0,:,:] - field[0][0,:,:], np.linspace(-minmax,minmax,51), extend='both', cmap=plt.cm.RdBu_r)
plt.colorbar()
plt.title('Difference between potential and mhs on the solar surface')

print('Max difference is', np.abs(potfield[0][0,:,:] - field[0][0,:,:]).max())

# minmax = 20
# miss = 16
# plt.figure(7)
# plt.contourf(synmap[1][::miss], synmap[2][::miss], synmap[0][::miss,::miss], np.linspace(-minmax,minmax,51), extend='both', cmap=plt.cm.RdBu_r)

br_sin = field[0][:,:,:]
for i, theta in enumerate(field[4]):
    br_sin[:, i, :] = br_sin[:, i, :]*np.sin(theta)

ints = field[3]**2*np.trapz(np.trapz(np.abs(br_sin), x=field[5], axis=-1), x=field[4])

br_sin = potfield[0][:,:,:]
for i, theta in enumerate(potfield[4]):
   br_sin[:, i, :] = br_sin[:, i, :]*np.sin(theta)
potints = potfield[3]**2*np.trapz(np.trapz(np.abs(br_sin), x=field[5], axis=-1), x=field[4])

plt.figure()
plt.plot(potfield[3], potints)
plt.title('Integrals over each radial surface for potential field')
plt.figure()
plt.plot(field[3], ints)
#plt.xlim(1,2.2)
plt.title('Integrals over each radial surface for mhs field')

minmax = 0.05
plt.figure()
plt.contourf(potfield[5], potfield[4], potfield[0][-1,:,:], np.linspace(-minmax,minmax,51), extend='both', cmap=plt.cm.RdBu_r)
plt.colorbar()
plt.contour(potfield[5], potfield[4], potfield[0][-1,:,:], [0])
plt.ylim(np.pi, 0)
plt.title(r'PFSS $B_r$ on the TOP')

# print(ints, 'integrals over surface')
# print(potints, 'integrals over surface')

if np.any(np.isfinite(field[0]) == False) == True: print('Problem in Br')
if np.any(np.isfinite(field[1]) == False) == True: print('Problem in Bt')
if np.any(np.isfinite(field[2]) == False) == True: print('Problem in Bp')

plt.show()
