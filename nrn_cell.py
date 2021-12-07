'''
Contains Python Class to instantiate a "Cell" object into NEURON, 
by importing an SWC morphology file, setting biophysical parameters, 
potentially adding synapses, etc. 

NOTE: to use, input your own user_token below, which can be generated
by creating an account on the neuprint website
Code also requires NEURON 7.7.2 and Python 3.7 installed in your environment. 
The requirements.txt provides package versions, but I do 
NOT recommend installing all packages in the txt file with pip/conda, 
as there are many that are not relevant for this project and may introduce
compatibility issues. Instead, use it as a reference in case version
issues arise here with the import statements below. 

Author: Tony X Liu
'''

### USER INPUT: 
user_token = '' 

from neuron import h, gui
from neuron.units import ms, mV # allows us to specify units
import sys
import os
import re
import pdb
import pandas as pd
import numpy as np
import dill as pickle
from datetime import datetime
import seaborn as sns
from matplotlib import cm
from scipy import stats
import matplotlib
from sklearn.neighbors import NearestNeighbors

# set up API connection to neuprint hemibrain server
from neuprint import Client # version 0.4.11
from neuprint import fetch_simple_connections, fetch_synapse_connections, fetch_neurons
from neuprint import SynapseCriteria as SC, NeuronCriteria as NC
try:
	c = Client('neuprint.janelia.org', dataset = 'hemibrain:v1.1', token = user_token)
except:
	print('neuprint client connection failed, possibly no WiFi')

h.load_file('import3d.hoc')

class Cell():
	def __init__(self, fname, gid):
		'''
		fname - string: file path to SWC
		'''
		self._gid = gid
		self.fname = fname
		self.load_morphology(fname)
		self.sec_name = ""
		self.find_sec_name()
		#self.discretize_sections()
		#self.add_biophysics()
		self.tree = pd.DataFrame()
		self.body_id = 0

	def __str__(self):
		return 'Cell[{}]'.format(self._gid)

	def load_morphology(self, nom):
		#print(nom)
		cell = h.Import3d_SWC_read()
		cell.input(nom)
		i3d = h.Import3d_GUI(cell, 0)
		i3d.instantiate(self)

	def find_sec_name(self):
		secnames = []
		for sec in h.allsec():
			name = re.split('\.|\[', sec.name())[2]
			if name not in secnames:
				secnames.append(name)
		if 'dend' in secnames:
			self.sec_name = 'dend'
		elif 'axon' in secnames:
			self.sec_name = 'axon'
		elif 'soma' in secnames:
			self.sec_name = 'soma'

	def discretize_sections(self):
		''' 
			adds at least as many spatial compartments as d_lambda rule
			maximizing segment density also allows better synapse localization 
		'''
		for sec in h.allsec():
			sec.nseg = sec.n3d()

	def add_biophysics(self, ra, cm, gpas, epas):
		# insert passive density mechanism
		mt = h.MechanismType(0)
		mt.select("pas")
		for section in h.allsec():
			# insert distributed mechanism into section
			mt.make(sec=section)	

		change_R_a(ra)
		change_c_m(cm)
		change_g_pas(gpas)
		change_e_pas(epas)

	def trace_tree(self):
		'''
			create table of all specified 3d points (0 to section.n3d()-1), x, y, z coordinates, 
		    (note, 3d point != segment, but arc3d(point i)/section length does give "x position" (0 to 1) of point)
		    and their associated section number (re.split('\[|\]', cell1.axon[192].name())[3] gives 192)
		'''
		tree = [] # tree is initially a list, dict to DataFrame is fastest to create the pandas DataFrame
		for sec in self.axon: # was self.axon
			num_segs = sec.n3d()
			sec_index = re.split('\[|\]', sec.name())[3]
			for i in range(num_segs):
				toAppend = {} 	# each row to add is a dictionary
				loc = sec.arc3d(i) / sec.L
				geodesic_dist = eval("h.distance(self.{}[0](0.5), sec(loc))".format(self.sec_name))
				toAppend.update(sec=sec_index, i3d=i, 
								x=sec.x3d(i), y=sec.y3d(i), z=sec.z3d(i), 
								arc=sec.arc3d(i), gd = geodesic_dist)
				tree.append(toAppend)
		tree = pd.DataFrame(tree)
		return tree

	def add_synapses(self, file_path, syn_strength):
		'''
			add Exp2Syn synapses to model, based on xyz synapse locations
			requires the "tree" DataFrame attribute to be populated
		'''
		### import synaptic locations
		conn = pd.read_csv(file_path)
		num_synapses = conn.shape[0]
		if num_synapses == 0:
			return 0, 0, 0, 0

		### KNN to map each synapse x, y, z (scaled x0.008) to the closest segment
		tree_coords = self.tree.loc[:, 'x':'z']
		syn_coords = conn.loc[:, 'x':'z'] / 125
		nbrs = NearestNeighbors(n_neighbors=1, algorithm='auto').fit(tree_coords)
		distances, indices = nbrs.kneighbors(syn_coords) 
		# indices: index in tree of closest section and point location to a synapse

		### add synapses onto morphology
		syns = h.List()
		j = 0 # index in syns
		for index in indices:
			sec = int(self.tree.loc[index, 'sec'])
			i3d = self.tree.loc[index, 'i3d']	# the 3d point index on the section
			#print("adding synapse " + str(j) + " to section " + str(sec))
			loc = eval("self.{}[sec].arc3d(i3d) / self.{}[sec].L".format(self.sec_name, self.sec_name))
			# 0 to 1, length along section
			#print("about to append")
			syns.append(h.Exp2Syn(self.axon[sec](loc)))

			### synapse parameters from Tobin et al paper: 
			syns.object(j).tau1 = 0.2 #ms
			syns.object(j).tau2 = 1.1 #ms
			syns.object(j).e = -10 #mV, synaptic reversal potential = -10 mV for acetylcholine? 
			# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3125135/
			#syns.object(j).g = 0.0001 #uS default ### seems to have no effect on the result

			h.pop_section() # clear the section stack to avoid overflow (triggered by using ".L" above?)
			j = j + 1

		### use NetStim to activate NetCon
		nc = h.NetStim()
		nc.number = 1
		nc.start = 0
		nc.noise = 0

		ncs = h.List()
		for i in range(len(list(syns))):
			ncs.append(h.NetCon(nc, syns.object(i)))
			ncs.object(i).weight[0] = syn_strength # uS, peak conductance change

		return syns, nc, ncs, num_synapses

	def add_synapses_xyz(self, xyz_locs, syn_strength):
		'''
			add new synapses based on loaded xyz locations
		'''
		num_synapses = xyz_locs.shape[0]
		#print("imported " + str(num_synapses) + " synapses")
		if num_synapses == 0:
			return 0, 0, 0, 0

		### KNN to map each synapse x, y, z (scaled x0.008) to the closest segment
		tree_coords = self.tree.loc[:, 'x':'z']
		syn_coords = xyz_locs.loc[:, 'x_post':'z_post'] / 125
		nbrs = NearestNeighbors(n_neighbors=1, algorithm='auto').fit(tree_coords)
		distances, indices = nbrs.kneighbors(syn_coords) 
		# indices: index in tree of closest section and point location to a synapse

		### add synapses onto morphology
		syns = h.List()
		j = 0 # index in syns
		for index in indices:
			sec = int(self.tree.loc[index, 'sec'])
			i3d = self.tree.loc[index, 'i3d']	# the 3d point index on the section
			#print("adding synapse " + str(j) + " to section " + str(sec))
			# TODO: I think could replace with cell1.axon[`sec`].X (4/2/2021)
			loc = eval("self.{}[sec].arc3d(i3d) / self.{}[sec].L".format(self.sec_name, self.sec_name))
			# 0 to 1, length along section
			#print("about to append")
			syns.append(h.Exp2Syn(self.axon[sec](loc)))

			### synapse parameters from Tobin et al paper: 
			syns.object(j).tau1 = 0.2 #ms
			syns.object(j).tau2 = 1.1 #ms
			syns.object(j).e = -10 #mV, synaptic reversal potential = -10 mV for acetylcholine? 
			# syns.object(j).bodyId_pre = xyz_locs.loc[j, 'bodyId_pre'] unfortunately can't add attrs
			# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3125135/
			#syns.object(j).g = 0.0001 #uS default ### seems to have no effect on the result

			h.pop_section() # clear the section stack to avoid overflow (triggered by using ".L" above?)
			j = j + 1

		### use NetStim to activate NetCon, initially inactive
		nc = h.NetStim()
		nc.number = 0
		nc.start = 0
		nc.noise = 0

		ncs = h.List()
		for i in range(len(list(syns))):
			ncs.append(h.NetCon(nc, syns.object(i)))
			ncs.object(i).weight[0] = syn_strength # uS, peak conductance change

		return syns, nc, ncs, num_synapses

	def add_synapses_subtree(self, sec_for_subtree, syn_count, syn_strength):
		'''
			add <syn_count> synapses to random sections in the subtree of 
			self.axon[<sec_for_subtree>]
		'''

		# get section numbers in the subtree
		subtree_secs = self.axon[sec_for_subtree].subtree()
		subtree_sec_nums_brack = [str(sec).partition('axon')[2] for sec in subtree_secs]
		subtree_sec_nums = [re.findall("\[(.*?)\]", sec)[0] for sec in subtree_sec_nums_brack] # debracket 

		### add synapses onto morphology
		syns = h.List()
		j = 0
		for index in range(syn_count):

			sec = int(random.choice(subtree_sec_nums))
			loc = random.uniform(0, 1)

			syns.append(h.Exp2Syn(self.axon[sec](loc)))

			### synapse parameters from Tobin et al paper: 
			syns.object(j).tau1 = 0.2 #ms
			syns.object(j).tau2 = 1.1 #ms
			syns.object(j).e = -10 #mV, synaptic reversal potential = -10 mV for acetylcholine? 
			# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3125135/
			#syns.object(j).g = 0.0001 #uS default ### seems to have no effect on the result

			h.pop_section() # clear the section stack to avoid overflow 
			j += 1

		### use NetStim to activate NetCon, initially inactive
		nc = h.NetStim()
		nc.number = 0
		nc.start = 0
		nc.noise = 0

		ncs = h.List()
		for i in range(len(list(syns))):
			ncs.append(h.NetCon(nc, syns.object(i)))
			ncs.object(i).weight[0] = syn_strength # uS, peak conductance change

		num_synapses = syn_count

		return syns, nc, ncs, num_synapses

	def total_length(self):
		total_length = 0
		for sec in self.axon: 
			total_length += sec.L
		return total_length

	def surf_area(self):
		total_area = 0
		for sec in h.allsec():
			for seg in sec:
				total_area += seg.area()
		return total_area

###
### change params:
### max synaptic conductance (across all synapses)
### global Ra (axial resistance)
### membrane conductance / membrane resistance
### membrane capacitance
###
param_print = False
def change_syn_strength(strength):
	global syn_strength
	for i in range(len(list(syns))):
		ncs.object(i).weight[0] = strength # uS, peak conductance change
	for i in range(len(list(syns2))):
		ncs2.object(i).weight[0] = strength # uS, peak conductance change
	syn_strength = strength
	if param_print:
		print("new synapse conductance: " + str(strength) + " uS")

def change_R_a(ra):
	global R_a
	for sec in h.allsec():
		sec.Ra = ra
	R_a = ra
	if param_print:
		print("new R_a: " + str(ra) + " ohm-cm")

def change_g_pas(gpas):
	global g_pas, R_m
	for sec in h.allsec():
		for seg in sec.allseg():
			if hasattr(seg, 'pas'):
				seg.pas.g = gpas
	g_pas = gpas
	R_m = 1/gpas
	if param_print:
		print("new g_pas: " + str(gpas) + " S/cm^2, new R_m: " + str(round(R_m/1000,2)) + " kOhm-cm^2")

def change_e_pas(epas):
	global e_pas

	for sec in h.allsec():
		for seg in sec.allseg():
			if hasattr(seg, 'pas'):
				seg.pas.e = epas
	e_pas = epas
	if param_print:
		print("new e_pas: " + str(epas) + " mV")


def change_c_m(cm):
	global c_m
	for sec in h.allsec():
		sec.cm = cm
	c_m = cm
	if param_print:
		print("new c_m: "+ str(cm) + " uF/cm^2")