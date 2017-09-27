#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#--------------------------------------------
# Date: 2017/09/26
# Author: Takumario
# Mail: sakamoto.takuma@ac.jaxa.jp
#--------------------------------------------
#
import numpy as np

class Optimize():
	"""Optimize the evaluate function by brute-force search"""

	def __init__(self):
		self.N = 1
		self.d = 1
		self.D = 1
		self.n = 1
		self.k = 1
		self.L = 1
		self.N_fin = 3 
		self.d_fin = 3  
		self.D_fin = 3 
		self.n_fin = 3 
		self.k_fin = 3 
		self.L_fin = 3 
		self.N_spn = 1
		self.d_spn = 1
		self.D_spn = 1
		self.n_spn = 1
		self.k_spn = 1
		self.L_spn = 1
		self.BL = 500   #Battery Limit
		self.LL = 100   #Latch Limit
		self.Result = np.empty((0, 7), float)
		self.initialise()

	def	BruteSearch(self):
		Value = 0
		while self.N <= self.N_fin:
			while self.d <= self.d_fin:
				while self.D <= self.D_fin:
					while self.n <= self.n_fin:
						while self.k <= self.k_fin:
							while self.L <= self.L_fin:
								self.Evaluate(self.BL, self.LL)
								self.L += 1
							self.L = 0
							self.k += 1
						self.k = 0
						self.n += 1
					self.n = 0
					self.D += 1
				self.D = 0 
				self.d += 1
			self.d = 0
			self.N += 1	
		print Value
		self.get_max(self.Result)
	
	def EF1(self): #I need to modify from here
		blue = (k*L-F/D*(1-A)*d**3)/(A*B*d**4/D**3*n+k)-k*L/(C*d**4/D**3*n + k)
		Value = blue*(-(k/2+A*B/2*d**4/D**3/n)*blue + (k*L-F*(1-A)*d**3/D))
		return Value
	
	def EF2():
		return 1

	def EF3():
		return 1
	
	def Evaluate(self, BL, LL):
		if 0 < self.BL:
				self.Result = np.append(self.Result, np.array([[1,2,3,4,5,6,7]]), axis=0)
	
	def get_max(self, Result):
		index = Result[:,6].argmax()
		print Result[index,:]

	def initialise():
		#Global constants
		global Ga = 18000
		global Gm = 6000
		global tau_s = 200
		global tau_f = 500
		global gamma_l = 0.5
		global A = (tau_f-tau_s)/(tau_f-tau_s+Gm*gamma_l)
		global B = Gm/8
		global C = Ga/8
		global F = np.pi*tau_s/8

if __name__ == '__main__':
	print "hello"
	Body = Optimize()
	Body.BruteSearch()
	print "行列の大きさ:", Body.Result.shape

