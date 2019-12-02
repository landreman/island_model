import numpy as np
import math

def fac(x):
	y = math.factorial(x)
	return y	

def alpha(m,l):
	y = (-1.0)**l/(fac(m+l)*fac(l)*2**(2*l+m))
	return y
	
def beta(m,l):
	y = fac(m-l-1)/(fac(l)*2**(2*l-m+1))
	return y

def Dmn(m, n, R, Z):
	y = 0.0
	for k in range(int(n/2) + 1):
		#print(k)
		sumD = 0.0
		for j in range(k+1):
			#print(alpha(m,j), beta(m,k-j), alpha(m, k-j), beta(m,j))
			sumD += ( - (2*k-2*j - m)*alpha(m,j)*beta(m,k-j)*R**(2*j+m) + (2*k-2*j + m)*alpha(m,k-j)*beta(m,j)*R**(2*j-m) )

		y += (Z**(n-2*k)/fac(n-2*k))*sumD

	#print("Dmn is ", y)
	return y

def Nmn(m, n, R, Z):
	y = 0.0
	for k in range(int(n/2) + 1):
		#print(k)
		sumN = 0.0
		for j in range(k+1):
			#print(alpha(m,j), beta(m,k-j), alpha(m, k-j), beta(m,j))
			sumN += ( alpha(m,j)*beta(m,k-j)*R**(2*j+m) - alpha(m,k-j)*beta(m,j)*R**(2*j-m) )

		y += (Z**(n-2*k)/fac(n-2*k))*sumN

	#print("Nmn is ", y)
	return y

def dRDmn(m, n, R, Z):
	y = 0.0
	for k in range(int(n/2) + 1):
		sumD = 0.0
		for j in range(k+1):
			sumD += ( - (2*j+m)*(2*k-2*j - m)*alpha(m,j)*beta(m,k-j)*R**(2*j+m-1) + (2*j-m)*(2*k-2*j + m)*alpha(m,k-j)*beta(m,j)*R**(2*j-m-1) )

		y += (Z**(n-2*k)/fac(n-2*k))*sumD

	return y

def dZDmn(m, n, R, Z):
	y = 0.0
	for k in range(int((n-1)/2) + 1):
		sumD = 0.0
		for j in range(k+1):
			sumD += ( - (2*k-2*j - m)*alpha(m,j)*beta(m,k-j)*R**(2*j+m) + (2*k-2*j + m)*alpha(m,k-j)*beta(m,j)*R**(2*j-m) )

		y += (Z**(n-2*k-1)/fac(n-2*k-1))*sumD

	return y

def dRRDmn(m,n, R, Z):
	y = 0.0
	for k in range(int(n/2) + 1):
		sumD = 0.0
		for j in range(k+1):
			sumD += ( - (2*j+m-1)*(2*j+m)*(2*k-2*j - m)*alpha(m,j)*beta(m,k-j)*R**(2*j+m-2) + (2*j-m-1)*(2*j-m)*(2*k-2*j + m)*alpha(m,k-j)*beta(m,j)*R**(2*j-m-2) )

		y += (Z**(n-2*k)/fac(n-2*k))*sumD

	return y

def dZZDmn(m, n, R, Z):
	y = 0.0
	for k in range(int(n/2)):
		sumD = 0.0
		for j in range(k+1):
			sumD += ( - (2*k-2*j - m)*alpha(m,j)*beta(m,k-j)*R**(2*j+m) + (2*k-2*j + m)*alpha(m,k-j)*beta(m,j)*R**(2*j-m) )

		y += (Z**(n-2*k-2)/fac(n-2*k-2))*sumD

	return y

def dRZDmn(m,n, R, Z):
	y = 0.0
	for k in range(int((n-1)/2) + 1):
		sumD = 0.0
		for j in range(k+1):
			sumD += ( - (2*j+m)*(2*k-2*j - m)*alpha(m,j)*beta(m,k-j)*R**(2*j+m) + (2*j-m)*(2*k-2*j + m)*alpha(m,k-j)*beta(m,j)*R**(2*j-m) )

		y += (Z**(n-2*k-1)/fac(n-2*k-1))*sumD

	return y

def dRNmn(m, n, R, Z):
	y = 0.0
	for k in range(int(n/2) + 1):
		sumN = 0.0
		for j in range(k+1):
			sumN += ( (2*j+m)*alpha(m,j)*beta(m,k-j)*R**(2*j+m-1) - (2*j-m)*alpha(m,k-j)*beta(m,j)*R**(2*j-m-1) )

		y += (Z**(n-2*k)/fac(n-2*k))*sumN

	return y

def dZNmn(m, n, R, Z):
	y = 0.0
	for k in range(int((n-1)/2) + 1):
		sumN = 0.0
		for j in range(k+1):
			sumN += ( alpha(m,j)*beta(m,k-j)*R**(2*j+m) - alpha(m,k-j)*beta(m,j)*R**(2*j-m) )

		y += (Z**(n-2*k-1)/fac(n-2*k-1))*sumN

	return y
		
def dRRNmn(m, n, R, Z):
	y = 0.0
	for k in range(int(n/2) + 1):
		sumN = 0.0
		for j in range(k+1):
			sumN += ( (2*j+m-1)*(2*j+m)*alpha(m,j)*beta(m,k-j)*R**(2*j+m-2) - (2*j-m-1)*(2*j-m)*alpha(m,k-j)*beta(m,j)*R**(2*j-m-2) )

		y += (Z**(n-2*k)/fac(n-2*k))*sumN

	return y

def dZZNmn(m, n, R, Z):
	y = 0.0
	for k in range(int(n/2)):
		sumN = 0.0
		for j in range(k+1):
			sumN += ( alpha(m,j)*beta(m,k-j)*R**(2*j+m) - alpha(m,k-j)*beta(m,j)*R**(2*j-m) )

		y += (Z**(n-2*k-2)/fac(n-2*k-2))*sumN

	return y
		
def dRZNmn(m, n, R, Z):
	y = 0.0
	for k in range(int((n-1)/2) + 1):
		sumN = 0.0
		for j in range(k+1):
			sumN += ( (2*j+m)*alpha(m,j)*beta(m,k-j)*R**(2*j+m-1) - (2*j-m)*alpha(m,k-j)*beta(m,j)*R**(2*j-m-1) )

		y += (Z**(n-2*k-1)/fac(n-2*k-1))*sumN

	return y

def Phi(m, nn, R, Z, phi):
	if (nn%2 == 0):
		a = d = 0
		b = 1
		c = 1
	else:
		a = 1 
		d = -1 
		b = c = 0
	
	y = (a*np.cos(m*phi) + b*np.sin(m*phi))*Dmn(m,nn,R,Z) + (c*np.cos(m*phi) + d*np.sin(m*phi))*Nmn(m,nn-1,R,Z)
	return y

def BR(m, n, R, Z, phi):
	if (n%2 == 0):
		a = d = 0
		b = 1
		c = 1
	else:
		a = 1 
		d = -1 
		b = c = 0

	y = (a*np.cos(m*phi) + b*np.sin(m*phi))*dRDmn(m,n,R,Z) + (c*np.cos(m*phi) + d*np.sin(m*phi))*dRNmn(m,n-1,R,Z)
	return y

def BZ(m, n, R, Z, phi):
	if (n%2 == 0):
		a = d = 0
		b = 1
		c = 1
	else:
		a = 1 
		d = -1 
		b = c = 0

	y = (a*np.cos(m*phi) + b*np.sin(m*phi))*dZDmn(m,n,R,Z) + (c*np.cos(m*phi) + d*np.sin(m*phi))*dZNmn(m,n-1,R,Z)
	return y

def Bphi(m, n, R, Z, phi):
	if (n%2 == 0):
		a = d = 0
		b = 1
		c = 1
	else:
		a = 1 
		d = -1 
		b = c = 0

	y = m*(-a*np.sin(m*phi) + b*np.cos(m*phi))*Dmn(m,n,R,Z)/R + m*(-c*np.sin(m*phi) + d*np.cos(m*phi))*Nmn(m,n-1,R,Z)/R
	return y

def dRBR(m, n, R, Z, phi):
	if (n%2 == 0):
		a = d = 0
		b = 1
		c = 1
	else:
		a = 1 
		d = -1 
		b = c = 0

	y = (a*np.cos(m*phi) + b*np.sin(m*phi))*dRRDmn(m,n,R,Z) + (c*np.cos(m*phi) + d*np.sin(m*phi))*dRRNmn(m,n-1,R,Z)
	return y

def dZBZ(m, n, R, Z, phi):
	if (n%2 == 0):
		a = d = 0
		b = 1
		c = 1
	else:
		a = 1 
		d = -1 
		b = c = 0

	y = (a*np.cos(m*phi) + b*np.sin(m*phi))*dZZDmn(m,n,R,Z) + (c*np.cos(m*phi) + d*np.sin(m*phi))*dZZNmn(m,n-1,R,Z)
	return y

def dRBZ(m, n, R, Z, phi):
	if (n%2 == 0):
		a = d = 0
		b = 1
		c = 1
	else:
		a = 1 
		d = -1 
		b = c = 0

	y = (a*np.cos(m*phi) + b*np.sin(m*phi))*dRZDmn(m,n,R,Z) + (c*np.cos(m*phi) + d*np.sin(m*phi))*dRZNmn(m,n-1,R,Z)
	return y

def dZBR(m, n, R, Z, phi):
	if (n%2 == 0):
		a = d = 0
		b = 1
		c = 1
	else:
		a = 1 
		d = -1 
		b = c = 0

	y = (a*np.cos(m*phi) + b*np.sin(m*phi))*dRZDmn(m,n,R,Z) + (c*np.cos(m*phi) + d*np.sin(m*phi))*dRZNmn(m,n-1,R,Z)
	return y

def dRBphi(m, n, R, Z, phi):
	if (n%2 == 0):
		a = d = 0
		b = 1
		c = 1
	else:
		a = 1 
		d = -1 
		b = c = 0

	y = m*(-a*np.sin(m*phi) + b*np.cos(m*phi))*dRDmn(m,n,R,Z)/R + m*(-c*np.sin(m*phi) + d*np.cos(m*phi))*dRNmn(m,n-1,R,Z)/R - m*(-a*np.sin(m*phi) + b*np.cos(m*phi))*Dmn(m,n,R,Z)/R**2 - m*(-c*np.sin(m*phi) + d*np.cos(m*phi))*Nmn(m,n-1,R,Z)/R**2
	return y

def dZBphi(m, n, R, Z, phi):
	if (n%2 == 0):
		a = d = 0
		b = 1
		c = 1
	else:
		a = 1 
		d = -1 
		b = c = 0

	y = m*(-a*np.sin(m*phi) + b*np.cos(m*phi))*dZDmn(m,n,R,Z)/R + m*(-c*np.sin(m*phi) + d*np.cos(m*phi))*dZNmn(m,n-1,R,Z)/R
	return y
