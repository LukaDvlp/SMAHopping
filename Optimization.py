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
		self.d_init = 0.8
		self.D_init = 5
		self.n_init = 3
		self.k_init = 0.8
		self.L_init = 60
		self.N = self.N_init
		self.d = self.d_init
		self.D = self.D_init
		self.n = self.n_init
		self.k = self.k_init
		self.L = self.L_init
		self.N_fin = 10 
		self.d_fin = 1.5 
		self.D_fin = 12
		self.n_fin = 20
		self.k_fin = 6.0 
		self.L_fin = 200 
		self.N_spn = 1
		self.d_spn = 0.1
		self.D_spn = 0.5
		self.n_spn = 1
		self.k_spn = 0.1
		self.L_spn = 2
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
								self.Evaluate()
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
	
	def EF1(self): #Evaluate function 1
		Value = -self.k*(self.dl_A-self.L)*self.dl_M*0.5
		return Value
	
	def EF2(self):
		Fa = -self.k*(self.dl_A-self.L)
		Value = (Fa-self.F1)*(self.dl_M-self.dl_A)*0.5 + self.F1*(self.dl_1-self.dl_A)*0.5
		return Value

	def EF3(self):
		Fa = -self.k*(self.dl_A-self.L)
		Value = (Fa-self.F1)*(self.dl_M-self.dl_A)*0.5 + (self.F1-self.F2)*(self.dl_M-self.dl_2)*0.5+(self.F2-self.F3)*(self.dl_2-self.dl_A)*0.5 + self.F3*(self.dl_1-self.dl_A)*0.5
		return Value
	
	def Evaluate(self):
		self.dl_A = self.k*self.L/(self.N*self.Ga*self.d**4/(8*self.D**3*self.n)+self.k)
		self.dl_Ay = np.pi*self.D**2*self.n/(self.Ga*self.d)*self.tau_ay
		self.dl_1 = self.D**2*self.n/(self.Gm*self.d)*(np.pi*self.tau_s+self.Gm*self.d/(self.D**2*self.n)*self.dl_A)
		self.dl_2 = self.D**2*self.n/(self.Gm*self.d)*(np.pi*self.tau_f+self.Gm*self.d/(self.D**2*self.n)*self.dl_A+np.pi*self.Gm*self.gamma_l)
		self.dl_My = self.D**2*self.n/(self.Gm*self.d)*(np.pi*self.tau_my+self.Gm*self.d/(self.D**2*self.n)*self.dl_A+np.pi*self.Gm*self.gamma_l)
		self.DEL_1 = (self.N*self.Gm*self.d**4/(8*self.D**3*self.n)*self.dl_A+self.k*self.L)/(self.N*self.Gm*self.d**4/(8*self.D**3*self.n)+self.k)
		self.DEL_2 = (self.N*self.Gm*self.d**4/(8*self.D**3*self.n)*self.dl_A+self.k*self.L+self.N*np.pi*self.d**3*self.Gm*self.gamma_l/(8*self.D))/(self.N*self.Gm*self.d**4/(8*self.D**3*self.n)+self.k)
		if self.DEL_1 < self.dl_1:
				self.dl_M = (self.N*self.Gm*self.d**4/(8*self.D**3*self.n)*self.dl_A+self.k*self.L)/(self.N*self.Gm*self.d**4/(8*self.D**3*self.n)+self.k)
				if self.CheckCondition():
						E = self.EF1()
				else:
						E = 0
				arr = [[self.N,self.d,self.D,self.n,self.k,self.L,E]]
				print arr
				self.Result = np.append(self.Result, np.array(arr), axis=0)
				
		elif self.dl_1 < self.DEL_1 and self.DEL_1 < self.dl_2:
				F_s = self.N*np.pi*self.d**3*self.tau_s/(8*self.D)
				F_f = self.N*np.pi*self.d**3*self.tau_f/(8*self.D)
				self.dl_M = ((F_f-F_s)/(self.dl_2-self.dl_1)*self.dl_A -F_s +(F_f-F_s)/(self.dl_2-self.dl_1) + self.k*self.L)/((F_f-F_s)/(self.dl_2-self.dl_1)+self.k)
				self.F1 = (F_f-F_s)/(self.dl_2-self.dl_1)*self.dl_A + F_s - (F_f-F_s)/(self.dl_2-self.dl_1)*self.dl_1
				if self.CheckCondition():
						E = self.EF2()
				else:
						E = 0
				arr = [[self.N,self.d,self.D,self.n,self.k,self.L,E]]
				print arr
				self.Result = np.append(self.Result, np.array(arr), axis=0)

		elif self.dl_2 < self.DEL_2:
				self.dl_M = ((self.N*self.Gm*self.d**4*self.dl_A)/(8*self.D**3*self.n)+(self.N*np.pi*self.d**3*self.Gm*self.gamma_l)/(8*self.D)+self.k*self.L)/((self.N*self.Gm*self.d**4)/(8*self.D**3*self.n)+self.k)
				F_s = self.N*np.pi*self.d**3*self.tau_s/(8*self.D)
				F_f = self.N*np.pi*self.d**3*self.tau_f/(8*self.D)
				self.F1 = -self.k*(self.dl_M-self.L)
				self.F2 = self.N*np.pi*self.d**3*self.tau_f/(8*self.D)
				self.F3 = (F_f-F_s)/(self.dl_2-self.dl_1)*self.dl_A + F_s - (F_f-F_s)/(self.dl_2-self.dl_1)*self.dl_1
				if self.CheckCondition():
						E = self.EF3()
				else:
						E = 0
				arr = [[self.N,self.d,self.D,self.n,self.k,self.L,E]]
				print arr
				self.Result = np.append(self.Result, np.array(arr), axis=0)
	
	def get_max(self, Result):
		index = Result[:,6].argmax()
		print Result[index,:]

	def initialise(self):
		#Global constants
		self.Ga = 16000
		self.Gm = 6000
		self.tau_s = 32.4 
		self.tau_f = 74.0
		self.tau_ay = 550
		self.tau_my = 290
		self.gamma_l = 0.05
		self.rho = 6.5
		self.c = 440
		self.Af = 80
		self.Ms = 20
		self.A = (self.tau_f-self.tau_s)/(self.tau_f-self.tau_s+self.Gm*self.gamma_l)
		self.B = self.Gm/8
		self.C = self.Ga/8
		self.F = np.pi*self.tau_s/8
		#Constraint Conditions
		self.Ltc_max = 200 #Maximum of Latch
		self.Bat_max = 4400   #Maximum of Battery
	
	def CheckCondition(self):
		Latch = -self.k*(self.dl_A-self.L) < self.Ltc_max
		Battery = 6*self.N*self.rho*(np.pi)**2*(0.5*self.d)**2*self.D*self.n/1000*self.c/1000*(self.Af-self.Ms) < self.Bat_max
		SpringIndex = 5 < self.D/self.d < 12
		Yelding = self.dl_A < self.dl_Ay and self.dl_M < self.dl_My
		MaxLength = np.pi*self.D*self.n > self.L
		return Latch and Battery and SpringIndex and Yelding and MaxLength

if __name__ == '__main__':
	print "hello"
	Body = Optimize()
	Body.BruteSearch()
	print "行列の大きさ:", Body.Result.shape

