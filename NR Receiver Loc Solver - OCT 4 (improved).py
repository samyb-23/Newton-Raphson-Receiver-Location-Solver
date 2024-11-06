# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 17:50:36 2024

@author: Samy Boutros

Step-by-step description: 
Goal of this script is to numerically calculate the coordinates (in xyz ECEF 
frame) of a receiver when given three or more satellite pseudorange
measurements. This is useful to give insight into the precision of the 
calculation which can tell us how our ionosphere induces offsets into the
distance measurements due to signal delay and clock differences for the 
receiver and the satellites.

First we begin by accessing the satellite locations
which we get from TLEs downloaded from the following celestrak website:
http://celestrak.org/NORAD/elements/gp.php?GROUP=gnss&FORMAT=tle . We then 
specify which satellites we want to find the coordinates of. I choose these
satellites by looking at the RINEX data files and looking at the line
which says which satellites are in view of the receiver at midnight on a given
day. I then get the satellite coordinates in lat,long,alt at the desired time
of midnight on the desired day, then convert those to xyz coordinates in the 
WGS84 frame (the same as the receiver). 

Next we need the psuedorange measurements which we obtain through RINEX files.
The information we care about from these files is the pseudorange measurements
C1 and C2, which are pseudoranges measured at radio frequencies of 1575.42MHz
and 1227.60 MHz, respectively. Make sure to use RINEX files on days within 
a couple days of the TLE file you download. Once the pseudorange data is 
obtained we can begin numerically solving for the receiver coordinates.

The Newton-Raphson (NR) Solver uses as a convergence check the chi square test. 
If the chi square between consecutive iterations is much less than 1 (in this
case 10^-10) then take these coordinates as the result. We make an initial guess
somewhat near the receiver location (although I checked and it converges even
for a guess of (0,0,0)). The step size change is controlled by the gradient as
well as the curvature of the chi square, following the theory in the least 
squares notes written by my supervisor Prof. Jonathan Sievers. There is also
a variation of the NR method called the Levenberg-Marquart (LM) method, which
is more robust if convergence with NR was not possible, however it could be much
slower. Anyway, it converges with NR so LM was not needed in the end even 
though I kept it in the script (a description of NR and LM is also in the
least squares notes). Finally, when convergence occurs, the script prints 
the calculated coordinates, compares them to the actual receiver coordinates
(which can be found in the RINEX files) and takes the difference between these
lcations as the offset induced by the ionosphere, which ends up decreasing 
with the use of more satellites in the NR solver. 
"""

#%%
import numpy as np
from skyfield.api import load, wgs84
import georinex as gr
import pyproj

#%% GETTING SATELLITE COORDINATES USING TLE DATA
# Specify the desired TLE file name, extension, and location if not in
# the current working directory. Also specify the url that the TLE is
# being downloaded from, I use celestrak 
name = 'gnss 2024-10-04 13-14-08'
ext = '.tle'
loc = 'TLEs'
fname = loc + '/' + name + ext
stations_url='http://celestrak.org/NORAD/elements/gp.php?GROUP=gnss&FORMAT=tle'
    
# Load the TLE objects into a variable called satellites
satellites = load.tle(stations_url, filename = fname)

# Here I specify which SVs are present at midnight on June 25th 
PRNs = ['04', '05', '07', '08', '09', '16', '20', '27', '30'] # USE FOR C1
#PRNs = ['04', '05', '07', '08', '09', '27', '30'] # USE FOR C2


# To use Russian satellites Need to use website: 
# https://www.csno-tarc.cn/en/glonass/constellation to
# convert from conventional PRN R0X representation to the ones used in TLEs.
# Some of them require you to add a 'K' at the end 
# Rs = [R05, R07, R09, R15, R16, R24] (How they appear on RINEX files)
#Rs = ['756', '745', '702K', '757', '761', '760'] (Their codes in TLE files)
Rs = ['745', '702K', '761', '760'] #756 AND 757 MISSING C1 MEASUREMENTS


# The rest of this cell follows the Skyfield Earth Satellites documentation 
#here: https://rhodesmill.org/skyfield/earth-satellites.html 
# USA GPS satellites
SVgs = np.asarray([satellites['PRN' + ' ' + p] for p in PRNs])
# Russian COSMOS satellites
SVrs = np.asarray([satellites[r] for r in Rs])
# Combining satellite arrays
SVs = np.concatenate((SVgs, SVrs))
# Checking epochs to find appropriate TLE
epochs = np.asarray([SV.epoch for SV in SVs])

ts = load.timescale()

# Specify the epoch of the satellite we want the coordinates from 
t = ts.utc(2024, 10, 4, 0, 0, 0)
    
geos = np.asarray([SVs[i].at(t) for i in range(len(SVs))])

# Loads the positions of the satellites (latitude, longitude, altitude)
latlon_sv = np.asarray([wgs84.latlon_of(g) for g in geos])
lat_sv = np.asarray([l[0].degrees for l in latlon_sv])
lon_sv = np.asarray([l[1].degrees for l in latlon_sv])
alt_sv = np.asarray([wgs84.height_of(g).m for g in geos])

# latlongalt to xyz in WGS84/ECEF coordinate frame (same as receiver)
transformer = pyproj.Transformer.from_crs(
    {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
    {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'}
    )
x_sv, y_sv, z_sv = transformer.transform(lon_sv, lat_sv, alt_sv, radians = False)
r_sv = np.asarray([np.asarray([x_sv[i], y_sv[i], z_sv[i]]).reshape((3,1)) 
                    for i in range(len(geos))])

#assert(1==0)
#%% GETTING SATELLITE PSEUDORANGE DATA FROM RINEX FILES
def get_unavco_name(year=2024,day=215,station='drao'):
    year=repr(year)
    assert(len(year)==4)
    day=f'{day:03d}'
    dir='data.unavco.org/archive/gnss/rinex/obs/'+year+'/'+day+'/'
    fname=station+day+'0.'+year[-2:]+'d.Z'
    url='https://'+dir+fname
    return url,fname

station='drao'
year=2024
day=278 # October 3rd = day 277

# List that holds time objects 
stufft=[]
# List that holds C1 frequency objects at midnight
cstuff1 = []
# List that holds C2 frequency objects at midnight
cstuff2 = []

# List that holds C1 frequency objects 718 minutes after midnight 
orbit1 = []
# List that holds C1 frequency objects 1436 minutes after midnight
orbit2 = []

# List that holds C2 frequency objects 718 minutes after midnight
orbit1C2 = []
# List that holds C2 frequency objects 1436 minutes after midnight
orbit2C2 = []

# 32 GPS SVs plus 20 Russian SVs (R06, R10, R13, R23 are missing)
# Specifies which satellite to look at: ind = 0 is G01, 1 is G02 an so on
ind = [3, 4, 6, 7, 8, 15, 19, 26, 29, 37, 39, 44, 51]#, 36, 43]
#ind = [3, 4, 6]#, 7]#, 8]#, 26]#, 29]#, 37]#, 39]#, 44]#, 51]# USE FOR C2 
# The 36th index of SVs which is R05 SV has missing O2 info so we remove that,
# same for the R15 SVwhich is index 43
url,fname=get_unavco_name(year,day,station)
fname=fname[:-2]
dat=gr.load('unavco_files/' + fname) # directory where unavco files are kept
print('read ',fname)
for i in ind:
    cstuff1.append(dat['C1'][:, i][0])
    cstuff2.append(dat['C2'][:, i][0])
    # 1437 for the index where 718 minutes has passed and SV returns to same loc
    orbit1.append(dat['C1'][:, i][1437])
    orbit2.append(dat['C1'][:, i][2873])
    orbit1C2.append(dat['C2'][:, i][1437])
    orbit2C2.append(dat['C2'][:, i][2873])
    tt=dat['time']
    asdf=tt.to_pandas()
    asdf=tt.to_numpy()
    stufft.append(asdf.astype('int64')) 
# Convert unavco file data into numpy arrays for C1 and C2 frequencies
C1=np.hstack([crap.as_numpy() for crap in cstuff1])
C2=np.hstack([crap.as_numpy() for crap in cstuff2])
# O1 is nan so compare C1s with O2s to use in N matrix for Newtons method
O1=np.hstack([crap.as_numpy() for crap in orbit1])
O2=np.hstack([crap.as_numpy() for crap in orbit2])
# OOX same as OX just for C2 frequencies
OO1=np.hstack([crap.as_numpy() for crap in orbit1C2])
OO2=np.hstack([crap.as_numpy() for crap in orbit2C2])
tvec=np.hstack(stufft) # 1e9 times time measurements separated by 30s
t0=tvec[0]/1e9 # Converting time into seconds (initial time)
tt=(tvec-tvec[0])/1e9 # delta t in seconds
#assert(1==0)

#%% DEFINING FUNCTIONS, VECTORS AND MATRICES TO USE IN NEWTON-RAPHSON SOLVER

xr, yr, zr = dat.position

# Newton-Raphson Iterations
iters = 0
maxiter = 100
M = len(ind) # number of equations for number of SVs we want to use in NR method
N = 3 # number of variables/unknowns x, y, z (no model for offset, exclude it)

# receiver location guess
m = x_rc,y_rc,z_rc = np.asarray([-2e6,-3e6,5e6]).reshape((N,1)) 

# define distance calc
def distance(v1, v2):
    return np.sqrt((v1[0]-v2[0])**2+(v1[1]-v2[1])**2+(v1[2]-v2[2])**2)

# define gradients for my A' matrix later on which I call "gradients"
def grad_x(v1, v2):
    return (v1[0]-v2[0])/((v1[0]-v2[0])**2+(v1[1]-v2[1])**2+(v1[2]-v2[2])**2)**0.5
def grad_y(v1, v2):
    return (v1[1]-v2[1])/((v1[0]-v2[0])**2+(v1[1]-v2[1])**2+(v1[2]-v2[2])**2)**0.5
def grad_z(v1, v2):
    return (v1[2]-v2[2])/((v1[0]-v2[0])**2+(v1[1]-v2[1])**2+(v1[2]-v2[2])**2)**0.5

# define curvatures which I use for A'' matrix, which I call "curvatures"
# These are not needed if we use the approximation to the curvature of chi^2
def grad2_x(v1, v2):
    return (1/np.sqrt((v1[0]-v2[0])**2+(v1[1]-v2[1])**2+(v1[2]-v2[2])**2) - 
            (v1[0]-v2[0])**2/(np.sqrt((v1[0]-v2[0])**2+
                                      (v1[1]-v2[1])**2+(v1[2]-v2[2])**2))**3)
def grad2_y(v1, v2):
    return (1/np.sqrt((v1[0]-v2[0])**2+(v1[1]-v2[1])**2+(v1[2]-v2[2])**2) - 
            (v1[1]-v2[1])**2/(np.sqrt((v1[0]-v2[0])**2+
                                      (v1[1]-v2[1])**2+(v1[2]-v2[2])**2))**3)
def grad2_z(v1, v2):
    return (1/np.sqrt((v1[0]-v2[0])**2+(v1[1]-v2[1])**2+(v1[2]-v2[2])**2) - 
            (v1[2]-v2[2])**2/(np.sqrt((v1[0]-v2[0])**2+
                                      (v1[1]-v2[1])**2+(v1[2]-v2[2])**2))**3)

# define chi^2 and dchi^2/dm and d^2chi^2/dm^2 approx and exact equations
def chi_square(d, A, Nmat):
    d = d.reshape((M,1))
    return (d-A).T@np.linalg.inv(Nmat)@(d-A)
def grad_chi(d, A, A_prime, Nmat):
    d = d.reshape((M,1))
    return -2 * A_prime.T @ np.linalg.inv(Nmat) @ (d - A)
def curv_chi_approx(A, A_prime, Nmat):
    return 2*A_prime.T@np.linalg.inv(Nmat)@A_prime
def curv_chi_exact(d, A, A_prime, A_2prime, Nmat):
    d = d.reshape((M,1))
    return -2*A_2prime.T@np.linalg.inv(Nmat)@(d-A)+2*A_prime.T@np.linalg.inv(Nmat)@A_prime

# This will be used for my A vector (model)
dist_vec = np.asarray([distance(m, r_sv[i]) for i in range(M)])
# A' (gradients matrix wrt x, y and z for each satellite)
gradients = np.zeros((M,N))
# A'' (curvatures matrix wrt x, y and z for each satellite)
curvatures = np.zeros((M,N))
for i in range(M):
    gradients[i][0] = grad_x(m, r_sv[i])
    gradients[i][1] = grad_y(m, r_sv[i])
    gradients[i][2] = grad_z(m, r_sv[i])
    curvatures[i][0] = grad2_x(m, r_sv[i])
    curvatures[i][1] = grad2_y(m, r_sv[i])
    curvatures[i][2] = grad2_z(m, r_sv[i])

# Take the variance between the C1 measurements 1436 minutes apart since
# the satellite orbit is half that long (O1 is nan so had to use O2)
N_matrix = np.zeros((M,M))
for i in range(M):
    arr = np.asarray([C1[i], O2[i]])
    var = np.var(arr)
    N_matrix[i][i] = var
    

chi = chi_square(C1, dist_vec, N_matrix)
L = 0 # Lev-Marq Step (variation of NR that is more robust in case NR doesn't 
# converge)

#assert(1==0)
#%% NEWTON-RAPHSON SOLVER
while iters < maxiter:
    iters += 1

    
    # Update newly calculated receiver location and time offset
    m_new = m - (np.linalg.inv(curv_chi_approx(dist_vec, gradients, N_matrix))
                 @ grad_chi(C1, dist_vec, gradients, N_matrix))/(1+L)
    #m_new = m - (np.linalg.inv(curv_chi_exact(C1, dist_vec, gradients, 
      #                                        curvatures, N_matrix))
         #        @ grad_chi(C1, dist_vec, gradients, N_matrix))/(1+L)
    
    # Update dist_vec/radients for new m parameters
    dist_vec_new = np.asarray([distance(m_new, r_sv[i]) for i in range(M)])
        
    gradients_new = np.zeros((M,N))
    curvatures_new = np.zeros((M,N))
    for i in range(M):
        gradients_new[i][0] = grad_x(m_new, r_sv[i])
        gradients_new[i][1] = grad_y(m_new, r_sv[i])
        gradients_new[i][2] = grad_z(m_new, r_sv[i])
        curvatures_new[i][0] = grad2_x(m_new, r_sv[i])
        curvatures_new[i][1] = grad2_y(m_new, r_sv[i])
        curvatures_new[i][2] = grad2_z(m_new, r_sv[i])
    
    # New chi square after updating receiver location
    chi_new = chi_square(C1, dist_vec_new, N_matrix)
    
    # If chi square changes very little, take the new parameters as the 
    # best parameters for the model
    if abs(chi - chi_new) < 1e-10:
        print('-----------')
        print('Iterations = ', iters)
        print('New Receiver Coords: ', np.hstack(m_new))
        print('Actual Receiver Coords = ', np.array([xr, yr, zr]))
        print('Distance btw Calculated and Actual Location (km)= ', 
              float(distance(m_new, dat.position))/1e3)
        print('Chi Square: ', float(chi_new))
        print('Distance from Center of Earth [Actual], [Expected]: ',
              [np.sqrt(xr**2+yr**2+zr**2)], [np.sqrt(m_new[0]**2+
                                                     m_new[1]**2+m_new[2]**2)])
        print('-----------')
        break

    # check if chi square decrease, if not take increase L and restart from the 
    # same spot in parameter space
    if chi_new < chi:
        m = m_new
        chi = chi_new
        dist_vec = dist_vec_new
        gradients = gradients_new
        curvatures = curvatures_new
        if L != 0:
            L -= 0.1
        else:
            L == 0
    else:
        L += 1
        
    #print('New Receiver Coords: ', np.hstack(m_new))
    #print('Actual Receiver Coords = ', np.array([xr, yr, zr]))
    #print('Distance btw Calculated and Actual Location (km)= ', 
    #      float(distance(m_new, dat.position))/1e3)
    #print('Chi Square: ', float(chi_new))
    #print('Delta Chi Square: ', abs(chi_new - chi))
    #print('-----------')

