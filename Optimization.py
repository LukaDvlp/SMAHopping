#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#--------------------------------------------
# Date: 2017/09/26
# Author: Takumario
# Mail: sakamoto.takuma@ac.jaxa.jp
#--------------------------------------------
#

class Optimize():
	"""Optimize the evaluate function by brute-force search"""

	def __init__(self):
		self.N = 1
		self.d = 1
		self.D = 1
		self.n = 1
		self.k = 1
		self.L = 1
		self.N_fin = 10
		self.d_fin = 10
		self.D_fin = 10
		self.n_fin = 10
		self.k_fin = 10
		self.L_fin = 10
		self.N_spn = 1
		self.d_spn = 1
		self.D_spn = 1
		self.n_spn = 1
		self.k_spn = 1
		self.L_spn = 1

	def	BruteSearch(self):
		Value = 0
		while self.N <= 10:
			while self.d <= 10:
				while self.D <=10:
					while self.n <= 10:
						while self.k <= 10:
							while self.L <= 10:
								Value += 1
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

print "hello"
Body = Optimize()
Body.BruteSearch()
