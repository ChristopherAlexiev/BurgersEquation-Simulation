###IMPORTS
import numpy
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot, cm
from matplotlib.colors import Normalize
'''
import numpy
from matplotlib import pyplot, cm
from matplotlib.colors import Normalize
'''
#NOTES
#currently ignores first index of nozzle

###FUNCTIONS

def equation_of_motion(u,v, dt, nu, dx, dy):
	# generate the next state as a function of the old state
	#snapshots of old u and v lists:
	uXY = u[1:-1,1:-1]
	uXm1Y = u[1:-1,0:-2]
	uXYm1 = u[0:-2, 1:-1]
	uXp1Y = u[1:-1, 2:]
	uXYp1 = u[2:,1:-1]
	vXY = v[1:-1,1:-1]
	vXm1Y = v[1:-1,0:-2]
	vXYm1 = v[0:-2, 1:-1]
	vXp1Y = v[1:-1, 2:]
	vXYp1 = v[2:,1:-1]	

	un = u.copy()
	vn = v.copy()

	u[1:-1,1:-1]=(uXY-dt/dx*uXY*(uXY-uXm1Y)-dt/dy*vXY*(uXY-uXYm1)+nu*dt/dx**2*(uXp1Y-2*uXY+uXm1Y)+nu*dt/dy**2*(uXYp1-2*uXY+uXYm1))
	v[1:-1,1:-1]=(vXY-dt/dx*uXY*(vXY-vXm1Y)-dt/dy*vXY*(vXY-vXYm1)+nu*dt/dx**2*(vXp1Y-2*vXY+vXm1Y)+nu*dt/dy**2*(vXYp1-2*vXY+vXYm1))
	return (u,v)

def boundary(u, v, nozzle_u, nozzle_v, nx, ny, frame):
	u[0, :] = 0
	u[-1, :] = 0
	u[:, 0] = 0
	u[:, -1] = 0

	v[0, :] = 0
	v[-1, :] = 0
	v[:, 0] = 0
	v[:, -1] = 0

	# special nozzle BC
	indexStart=int(ny/2-2)
	indexEnd=int(ny/2+2)
	u[indexStart:indexEnd, 0] = nozzle_u[frame]
	v[indexStart:indexEnd, 0] = nozzle_v[frame]

	return (u, v)

def simulate(u, v, startFrame, endFrame,nx,ny,sigma,nu,nozzle_u,nozzle_v):
	#computer-made vars 
	dx = Lx/float(nx-1)
	dy = Ly/float(ny-1)

	dt = sigma*dx*dy/nu
	for frame in range(startFrame,endFrame):
		(u, v) = equation_of_motion(u, v,  dt, nu, dx, dy)
		#now at frame+1
		(u, v) = boundary(u, v, nozzle_u, nozzle_v, nx, ny, frame+1)
	return (u, v)


###RUN

###VARS
Lx = 2
Ly = 2
nx = 41
ny = 41
sigma = 0.001
nu = 0.01
nt = 2510
initial_u = numpy.zeros((ny, nx))
initial_v = numpy.zeros((ny, nx))
#nozzle located at (0, 1)
nozzle_u = numpy.append(10*numpy.ones(1000), numpy.zeros(nt))
nozzle_v = numpy.append(10*numpy.ones(1000), numpy.zeros(nt))



###RECORD IMAGES
increment = 50 
frame = -1 #shows current frame (-1 is initial without nozzle, nozzle starts at 0), last frame is nt

final_u = initial_u
final_v = initial_v

#print 0 frame
(final_u, final_v) = simulate(final_u, final_v, frame, frame+1,nx,ny,sigma,nu,nozzle_u,nozzle_v)
frame+=1
ax = pyplot.figure()
norm = Normalize()
magnitude = numpy.sqrt(final_u[::2]**2 + final_v[::2]**2)
pyplot.quiver(final_u[::2], final_v[::2], norm(magnitude), scale=60,
cmap=pyplot.cm.jet)
ax.savefig('frame'+str(frame).zfill(5)+'.png',dpi=300)
ax.clear()
while frame+increment<nt:
	(final_u, final_v) = simulate(final_u, final_v, frame, frame+increment,nx,ny,sigma,nu,nozzle_u,nozzle_v)
	frame+=increment
	ax = pyplot.figure()
	norm = Normalize()
	magnitude = numpy.sqrt(final_u[::2]**2 + final_v[::2]**2)
	pyplot.quiver(final_u[::2], final_v[::2], norm(magnitude), scale=60,
	cmap=pyplot.cm.jet)
	ax.savefig('frame'+str(frame).zfill(5)+'.png',dpi=300)
	ax.clear()
