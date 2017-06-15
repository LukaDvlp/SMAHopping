import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

def CalE():
	if  delta_AE<delta_s+X0_SMA and X_1<delta_s+X0_SMA: 
		(E,X) = (integrate.quad(lambda x:-k*x+k*X0, delta_AE, X_1)[0] \
       -integrate.quad(lambda x:KM*(x-X0_SMA), delta_AE, X_1)[0],X_1)
	   
	elif  delta_AE<delta_s+X0_SMA and delta_s+X0_SMA<X_2 \
	and X_2<delta_f+X0_SMA: 
		(E,X) = (integrate.quad(lambda x:-k*x+k*X0, delta_AE, X_2)[0] \
       -integrate.quad(lambda x:KM*(x-X0_SMA), delta_AE, delta_s+X0_SMA)[0]\
	   -integrate.quad(lambda x:a*x+b1, delta_s+X0_SMA, X_2)[0],X_2)

	elif 	delta_AE<delta_s+X0_SMA and delta_f+X0_SMA<X_3: 
		(E,X) = (integrate.quad(lambda x:-k*x+k*X0, delta_AE, X_3)[0] \
       -integrate.quad(lambda x:KM(x-X0_SMA), delta_AE, delta_s+X0_SMA)[0]\
       -integrate.quad(lambda x:a*x+b1, delta_s+X0_SMA, delta_f+X0_SMA)[0]\
 	   -integrate.quad(lambda x:KM*(x-X0_SMA)+b2,delta_f+X0_SMA,X_3)[0],X_3)

	elif 	delta_s+X0_SMA<delta_AE and delta_AE<delta_f+X0_SMA\
	 and delta_s+X0_SMA<X_2 and X_2<delta_f+X0_SMA: 
		(E,X) = (integrate.quad(lambda x:-k*x+k*X0, delta_AE, X_2)[0] \
       -integrate.quad(lambda x:a*x+b1, delta_AE, X_2)[0],X_2)

	elif delta_s+X0_SMA<delta_AE and delta_AE<delta_f+X0_SMA\
	and delta_f+X0_SMA<X_3:  
		(E,X) = (integrate.quad(lambda x:-k*x+k*X0, delta_AE,X_3)[0] \
       -integrate.quad(lambda x:a*x+b1, delta_AE, delta_f+X0_SMA)[0]\
	   -integrate.quad(lambda x:KM*(x-X0_SMA)+b2,delta_f+X0_SMA, X_3)[0],X_3)

	elif delta_f+X0_SMA<delta_AE and delta_f+X0_SMA<X_3: 
		(E,X) = (integrate.quad(lambda x:-k*x+k*X0, delta_AE,X_3)[0] \
	   -integrate.quad(lambda x:KM*(x-X0_SMA)+b2,delta_AE, X_3)[0],X_3)

	return (E,X)

#UNCHANGABLE PARAMETERS 
g = 1.622
rho = 6.5
c = 440
Af = 80
Ms = 20
Tau_s = 72.4 #MPa
Tau_f = 114.0 #MPa
tau_ca = 550  #Yeilding sheer stress at Austenite [MPa]
tau_cm = 290  #Yeilding sheer stress at Martensite [MPa]
Ga = 22000
Gm = 8000
gamma_l = 0.05
Yeild = 100


#CHANGABLE PARAMETERS 
m = 0.7 
Theta = 30.0
number = 2 
n = 15 
D = 6.0  #Diameter of Coil[mm]
d = 1.0  #Diameter of Wire[mm]
solid_height = d*n ##solid height of SMA spring[mm]
k = 1.03 
X0 = 60 
X0_SMA = solid_height 
#X0_SMA = 45             #Natural Height of SMA[mm]
N_PLOT_MAX = k*X0 
GRANULARITY = 0.01
ELEMENTS = int(N_PLOT_MAX/GRANULARITY)
X_Plot_Max = 5*X0

#CALCULATION of the STATUS of SPRINGS
tau_s = Tau_s*number
tau_f = Tau_f*number
C = np.pi*d**3/8/D      #Conversion from stress to force
A = 8*D**3*n/Gm/d**4/number  #Martensite
B = 8*D**3*n/Ga/d**4/number  #Austenite
H = np.pi*D**2*n*gamma_l/d

mas_sma = number*rho*np.pi**2*(0.5*d)**2*D*n/1000 #Mass of SMA[g]
Q = mas_sma*c*(Af-Ms)/1000 #Energy(heat) Consumption[J]
F_s = C*tau_s  
F_f = C*tau_f
delta_s = A*F_s                            #Unit of delta_s is delta
delta_f = A*F_f + np.pi*D**2*n/d*gamma_l   #Unit of delta_f is delta
#xi = 0.5*np.cos(np.pi/(tau_s-tau_f)*(tau-tau_f)) + 0.5
#F_ = np.pi*d**3/8/D*Gm*gamma_l*xi 
F_ME = (X0-H-X0_SMA)/(A+1/k)
F_AE = (X0-X0_SMA)/(B+1/k)
#delta_ME = -F_ME/k + X0 #Unit of delta_ME is X
delta_AE = -F_AE/k + X0 #Unit of delta_AE is X
strain_a = (delta_AE-X0_SMA)/X0_SMA*100 #Strain of Austenite SMA
#strain_m = (delta_ME-X0_SMA)/X0_SMA*100 #Strain of Martensite SMA
strain_bi = (X0-delta_AE)/X0*100        #Strain of Bias Spring
L_max_a = np.pi*D*n      #Maximum Length of Austenite
L_max_m = np.pi*D*n*1.06 #Maximum Length of Martensite
KM = Gm*d**4/8/D**3/n*number
a = (F_f-F_s)/(delta_f-delta_s)
b1 = (F_s*(delta_f+X0_SMA)-F_f*(delta_s+X0_SMA))/(delta_f-delta_s)
b2 = F_f-KM*delta_f
X_1 = (k*X0+KM*X0_SMA)/(KM+k)
X_2 = (k*X0-b1)/(a+k)
X_3 = (KM*X0_SMA+k*X0-b2)/(KM+k)

#E = integrate.quad(lambda x: -(B+1/k)*x+X0, F_ME, F_AE)
E = 0.0
X = 0.0
E,delta_ME = CalE()
distance_X = 2*E*np.sin(2*Theta/360*2*np.pi)/m/g
distance_Y = E*np.sin(Theta/360*2*np.pi)/m/g
Efficiency = E/Q/1000 * 100
F_ca = tau_ca*C*number
F_cm = tau_cm*C*number
delta_cm = A*F_cm+H+X0_SMA
delta_ca = B*F_ca+X0_SMA

print ("F_s = %f[N], F_f = %f[N]" %(F_s,F_f))
print ("delta_s = %f[mm], delta_f = %f[mm]" %(delta_s,delta_f))
print("F_ME=%f, F_AE=%f" %(F_ME,F_AE))
print("delta_ME=%f, delta_AE=%f, deformation=%f[mm]" %(delta_ME,delta_AE,\
delta_ME-delta_AE))
#print ("Strain_a:%f[per], Strain_m:%f[per]" %(strain_a, strain_m))
print ("L_max_a:%f[mm], L_max_m:%f[mm]" %(L_max_a, L_max_m))
print ("Released energy: %f[mJ], Required energy:%f[J]" %(E,Q))
print ("Distance_X=%f[mm], Distance_Y=%f[mm]" %(distance_X, distance_Y))
print ("Efficiency=%f[per]" %Efficiency)
print ("Critical axial force  at Austenite=%f[N]" %F_ca)
print ("Critical axial force at Martensite=%f[N]" %F_cm)
print ("Critical deformation at Austenite=%f[mm]" %delta_ca)
print ("Critical deformation at Martensite=%f[mm]" %delta_cm)
print ("Mass of SMA=%f[g]" %mas_sma)

#MAKING FORCE-DEFLECTION ARRAY for PLOT
delta = np.zeros(ELEMENTS,np.float64)
F = np.zeros(ELEMENTS,np.float64)
delta_a = np.zeros(ELEMENTS, np.float64)
delta_bi = np.zeros(ELEMENTS, np.float64)
i=0
while i<ELEMENTS:
		F[i] = GRANULARITY*i
		#tau[i] = 8*D/np.pi/d**3 *F[i-1]
		if F[i] < F_s:
			delta[i] = A*F[i]+X0_SMA 
			delta_a[i] = B*F[i]+X0_SMA
			delta_bi[i] = -F[i]/k+X0
			i += 1
		
		elif F[i] < F_f:		 
				delta[i] = A*F[i]+np.pi*D**2*n/d*gamma_l*\
				(0.5*np.cos(np.pi/(tau_s-tau_f)*(8*D/np.pi/d**3*F[i]-tau_f))+0.5)+X0_SMA
				delta_a[i] = B*F[i]+X0_SMA
				delta_bi[i] = -F[i]/k+X0
				i += 1
		elif F[i] < F_cm:
				delta[i] = A*F[i] + H +X0_SMA
				delta_a[i] = B*F[i]+X0_SMA
				delta_bi[i] = -F[i]/k+X0
				i += 1				
		elif F[i] < F_ca:
			delta[i] = Yeild*(F[i]-F_cm) + delta_cm
			delta_a[i] = B*F[i]+X0_SMA
			delta_bi[i] = -F[i]/k+X0
			#if X_Plot_Max<delta[i] or X_Plot_Max<delta_a[i]:
			#		break
			i += 1		
		else:
				delta[i] = Yeild*(F[i]-F_cm) + delta_cm
				delta_a[i] = Yeild*(F[i]-F_ca) + delta_ca
				delta_bi[i] = -F[i]/k+X0
				#if X_Plot_Max<delta[i] or X_Plot_Max<delta_a[i]:
				#		break
				i += 1
print i
Force = F[0:i]
del_m = delta[0:i]
del_a = delta_a[0:i]
del_b = delta_bi[0:i]

plt.plot(del_m,Force, label="Martensite")
plt.plot(del_a,Force, label="Autenite")
plt.plot(del_b,Force, label="Bias Spring")

plt.title('Force-Deflection curve')
plt.xlabel('Length [mm]')
plt.ylabel('Force [N]')

plt.xlim([0,X0])
#plt.plot(F,delta)
#plt.plot(F,delta_a)
#plt.plot(F,delta_bi)
plt.show()
