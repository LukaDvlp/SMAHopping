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
		self.N_init = 1
		self.d_init = 1
		self.D_init = 1
		self.n_init = 1
		self.k_init = 1
		self.L_init = 1
		self.N = self.N_init
		self.d = self.d_init
		self.D = self.D_init
		self.n = self.n_init
		self.k = self.k_init
		self.L = self.L_init
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
							self.L = self.L_init
							self.k += 1
						self.k = self.k_init
						self.n += 1
					self.n = self.n_init
					self.D += 1
				self.D = self.D_init
				self.d += 1
			self.d = self.d_init
			self.N += 1	
		print Value
		self.get_max(self.Result)
	
	def EF1(self,N,d,D,n,k,L): #Evaluate function 1
		blue = (k*L-self.F/D*(1-self.A)*d**3)/(self.A*self.B*d**4/D**3*n+k)-k*L/(self.C*d**4/D**3*n + k)
		Value = blue*(-(k/2+self.A*self.B/2*d**4/D**3/n)*blue + (k*L-self.F*(1-self.A)*d**3/D))
		return Value
	
	def EF2():
		return 1

	def EF3():
		return 1
	
	def Evaluate(self, BL, LL):
		if 0 < self.BL:
				arr = [[self.N,self.d,self.D,self.n,self.k,self.L,self.EF1(self.N,self.d,self.D,self.n,self.k,self.L)]]
				self.Result = np.append(self.Result, np.array(arr), axis=0)
	
	def get_max(self, Result):
		index = Result[:,6].argmax()
		print Result[index,:]

	def initialise(self):
		#Global constants
		self.Ga = 18000
		self.Gm = 6000
		self.tau_s = 200
		self.tau_f = 500
		self.gamma_l = 0.5
		self.A = (self.tau_f-self.tau_s)/(self.tau_f-self.tau_s+self.Gm*self.gamma_l)
		self.B = self.Gm/8
		self.C = self.Ga/8
		self.F = np.pi*self.tau_s/8

if __name__ == '__main__':
	print "hello"
	Body = Optimize()
	Body.BruteSearch()
	print "行列の大きさ:", Body.Result.shape

