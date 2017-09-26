#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#------------------------------------
# Author: Takumario
# Mail: sakamoto.takuma@ac.jaxa.jp
# Date: 2017/09/26
#------------------------------------
#
import numpy as np

class Optimize:
	"""Optimize the evaluate function by brute-force serch"""
	
	def __init__(self):
		self.N = 1
		self.d = 0.8
		self.D = 5
		self.n = 3
		self.k = 0.5
		self.L = 60
		self.N_fin = 10
		self.d_fin = 1.5 
		self.D_fin = 10 
		self.n_fin = 15 
		self.k_fin = 2.8
		self.L_fin = 120
		self.N_spn = 1
		self.d_spn = 0.1
		self.D_spn = 1 
		self.n_spn = 1 
		self.k_spn = 0.1
		self.L_spn = 5

	def BruteSearch(self)
		while(
