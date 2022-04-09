import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as integrate


def Power(a,b):
	return a**b
def Cos(x):
	return np.cos(x)
def Sin(x):
	return np.sin(x)
def Sec(x):
	return 1/np.cos(x)
def Tan(x):
	return np.tan(x)

def initials_two_metronome(initial_list):
	# Constants
	d = 0       # delta
	b = 0.02  # beta
	mu = 0.02   # dimensionless Van der Pol
	x0 = 0.45   # theta zero

	def two_metronome(X,t):
		x1, v1, x2, v2 = X[0], X[1], X[2], X[3]
		x1dot = v1
		v1dot = -(((-1 + b*Power(Cos(x2),2))*(mu*v1*(-1 + Power(x1,2)/Power(x0,2)) + (1 + d)*Sin(x1) + b*Power(v1,2)*Cos(x1)*Sin(x1) + b*Power(v2,2)*Cos(x1)*Sin(x2)) - b*Cos(x1)*Cos(x2)*(mu*v2*(-1 + Power(x2,2)/Power(x0,2)) + b*Power(v1,2)*Cos(x2)*Sin(x1) + (1 + d)*Sin(x2) + b*Power(v2,2)*Cos(x2)*Sin(x2)))/(-1 + b*Power(Cos(x1),2) + b*Power(Cos(x2),2)))
		x2dot = v2
		v2dot = (Sec(x2)*(-(b*mu*v1*Power(x0,2)*Sec(x1)) + b*mu*v1*Power(x1,2)*Sec(x1) + b*mu*v2*Power(x0,2)*Sec(x2) - b*mu*v2*Power(x2,2)*Sec(x2) - mu*v2*Power(x0,2)*Power(Sec(x1),2)*Sec(x2) + mu*v2*Power(x2,2)*Power(Sec(x1),2)*Sec(x2) + b*Power(v2,2)*Power(x0,2)*Power(Sec(x1),2)*Sin(x2) + b*Power(x0,2)*Tan(x1) + b*d*Power(x0,2)*Tan(x1) + b*Power(v1,2)*Power(x0,2)*Sec(x1)*Tan(x1) - b*Power(x0,2)*Tan(x2) - b*d*Power(x0,2)*Tan(x2) + Power(x0,2)*Power(Sec(x1),2)*Tan(x2) + d*Power(x0,2)*Power(Sec(x1),2)*Tan(x2)))/ (Power(x0,2)*(b*Power(Sec(x1),2) + b*Power(Sec(x2),2) - Power(Sec(x1),2)*Power(Sec(x2),2)))
		return np.array([x1dot, v1dot, x2dot, v2dot])
	
	initial_list = np.array(initial_list)
	color = []
	n, N = 0, len(initial_list)  # to show progress

	for initial in initial_list:
		# Initial condition
		x1_0, v1_0, x2_0, v2_0 = initial		

		# Solve
		tmax, dt = 500, 0.01
		t = np.arange(0,tmax,dt)
		sol = integrate.odeint(two_metronome, initial, t)
		x1_sol = sol[:,0]
		x2_sol = sol[:,2]

		# Detect synchronization
		sync = False if (abs(x1_sol[-10]-x2_sol[-10]) > 0.04 or abs(x1_sol[-100]-x2_sol[-100]) > 0.04) else True
		if sync:
			color.append('r')
		else:
			color.append('b')

		n += 1
		if n%3==0: print(str(100*n/N)+' % Completed')

	x1_0s, x2_0s = initial_list[:,0], initial_list[:,2]
	f = open('variables d{:.2f} b{:.2f} mu{:.2f} x0{:.2f}.py'.format(d,b,mu,x0),'w')
	f.write('x1_0s = '+ repr(x1_0s) + '\n')
	f.write('x2_0s = '+ repr(x2_0s) + '\n')
	f.write('color = '+ repr(color) + '\n')
	f.close()

	# Plot (on time)
	plt.figure(1)
	#plt.legend()
	plt.title('Initial conditions and Sync',fontsize = 12)
	plt.xlabel('Initial angle 1 (x1)'); plt.ylabel('Initial angle 2 (x2)')
	plt.axhline(color='black',lw=1)
	plt.scatter(x1_0s, x2_0s, c=color) # s= 10
	plt.axis('equal')
	plt.grid()
	plt.show()

# initials
# [[x1_0, v1_0, x2_0, v2_0], [x1_0, v1_0, x2_0, v2_0], ...]
v1_0, v2_0 = 0,0
initial_list = [[x1,v1_0,x2,v2_0] for x1 in np.linspace(-1.5,1.5,100) for x2 in np.linspace(-1.5,1.5,100)]

initials_two_metronome(initial_list)






