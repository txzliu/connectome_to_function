'''
Functions contained include: 
- attrs_per_LHNmEPSP: calculate attributes per synapse onto list of LHNs
- visualize_inputs: instantiate a neuron into NEURON GUI environment, 
	potentially with associated synapses
- probe_mEPSPs: compute attributes per synapse onto single LHN
	contains examples of using NEURON to simulate uEPSPs and mEPSPs

Reference Jupyter notebook "generate attrs per LHN synapse" for example
of running the code and examining outputs. 

Author: Tony X Liu
'''

from nrn_cell import *

def attrs_per_LHNmEPSP(nrn_list_path = '21-05-07_LHN_SWCPointNo_to_NEURON.csv',
					   weight_threshold = 0, transf_freq = 0,
					   input_name = 'all_ePN_inputs', 
					   LHLNs_only = False, only_these_body_IDs = [],
					   use_var_time_step = False, fixed_time_step = 0.06, 
					   compt_manips = {'isopot_dendr_arbor': False, 'isopot_ax_arbor': False, 
									   'isopot_prim_ax': False, 'isopot_prim_dendr': False,
									   'isopot_whole_nrn': False, 
									   'isopot_Ra_val': 0.001, 
									   'highpot_prim_ax': False, 
									   'highpot_Ra_val': 10000,
									   'IAC_uniform_diam': False, 
									   'IAC_no_pinch': False, 'IAC_pinch_thresh': 0.14, 
									   'IAC_eucl_dist': False,
									   'IAC_3x_len': False},
					   params = {'R_a': 350, # ohm-cm
							  	'c_m': 0.6, # uF/cm^2
							  	'g_pas': 5.8e-5, # S/cm^2
							  	'e_pas': -55, # mV 
							  	'syn_strength': 5.5e-5},
					   param_str = ''):
	'''	for each neuron in a list (i.e. a list of LHNs), find information about EACH mEPSP of its
		input connections (for all inputs >weight_threshold synapses at all 
		locations on the neuron), such as mEPSP amplitude, impedance measures, synapse counts, etc.

	Args:
	----- 
	nrn_list_path (str): name of file that should be in the same folder, contains n_neurons x n_morphology_attrs, 
		where the morphological attributes include information about axon and dendrite start, extra arbors, etc.
		default value: '21-05-07_LHN_SWCPointNo_to_NEURON.csv', also files from '04-29', '04-18', '03-06', '02-16'
	weight_threshold (int): add synapses from inputs with >threshold total synapses onto the target
		a higher threshold will naturally evaluate fewer mEPSPs
	transf_freq (int): freq at which to evaluate impedance properties
	input_name (str): 'all_inputs' adds all synapses from inputs with threshold > weight_threshold
		'all_ePN_inputs' adds all synapses from ePNs (lPN, adPN)
	only_these_body_IDs (list of int): if list is not empty, will only evaluate mEPSPs for body IDs in the list
	compt_manips (dict): inputs into probe_mEPSPs, morphological manipulations to apply before running mEPSPs,
		i.e. whether to make certain parts of the neuron more 
		or less compartmenalized / change diameters of the IAC / etc. 
		in addition, boolean values here are used to name the output file
	params (dict): to be passed into probe_mEPSPs, to be passed and unpacked into visualize_inputs
	param_str (str): to be added to the save file CSV name, i.e. "Tobin", "GouwensCell1", etc
		when empty refers to default "our" params: 'R_a': 350, # ohm-cm
		'c_m': 0.6, # uF/cm^2
		'g_pas': 5.8e-5, # S/cm^2
		'e_pas': -55, # mV one parameter left same in both, we used -55 in our fits
		'syn_strength': 5.5e-5

	Returns:
	--------
	mEPSPs_allLHNs_df (pd.DataFram): [n_(synapses across all neurons in nrn_list) x n_(attributes)]
	'''
	nrns = pd.read_csv(nrn_list_path)

	# iterate through each target neuron, concatenate relevant file info
	mEPSPs_allLHNs_df = pd.DataFrame()
	for i in range(nrns.shape[0]):
		if LHLNs_only:
			if 'local' not in nrns.iloc[i]['lhn'] and 'LHLN' not in nrns.iloc[i]['lhn']:
				print('----- skipping non-local neurons')
				continue
		if len(only_these_body_IDs) > 0: 
			if nrns.iloc[i]['body_id'] not in only_these_body_IDs: continue

		# 21-04-19 TODO: add filter to only include ePN input synapses
		mEPSPs_oneLHN_df = probe_mEPSPs(target_name = nrns.iloc[i].lhn, target_body_id = nrns.iloc[i].body_id, 
									 input_name = input_name,
									 lhn_morph_path = nrn_list_path,
									 all_input_min_weight = weight_threshold,
									 transf_freq = transf_freq,
									 use_var_time_step = use_var_time_step, fixed_time_step = fixed_time_step,
									 params = params,
									 **compt_manips)

		mEPSPs_allLHNs_df = mEPSPs_allLHNs_df.append(mEPSPs_oneLHN_df)

	date = datetime.today().strftime("%y-%m-%d")
	compt_manip_str = '' 
	for key, val in compt_manips.items():
		if isinstance(val, bool) and val:
			compt_manip_str += '_' + key
	if len(param_str) > 0: compt_manip_str += '_' + param_str

	mEPSPs_allLHNs_df.to_csv(f'{date}_attrs_per_LHNmEPSP{compt_manip_str}.csv')
	return mEPSPs_allLHNs_df

def visualize_inputs(target_name = 'local5', target_body_id = 5813105722, input_name = 'VA6_adPN', 
					get_neuprint_swc = True,
					syn_locs = pd.DataFrame(), all_input_min_weight = 3,
					R_a = 350, # ohm-cm
					c_m = 0.6, # uF/cm^2
					g_pas = 5.8e-5, # S/cm^2
					e_pas = -55, # mV one parameter left same in both, we used -55 in our fits
					syn_strength = 5.5e-5): # uS
	''' given a downstream (target_name) neuron + ID, an upstream neuron, instantiate the
		downstream neuron, and potentially also instantiate synapses
		(potentially pulling from neuprint) for the sake of visualization
		utility function useful for rapidly looking at a neuron and its inputs
		i.e. go to ModelView in the GUI to see the morphology, synapse locations, etc.

		NOTE: preset model parameters (from Nov 2020 parameter fit) can be changed

	Args:
	-----
	target_name (string): name of neuron
	target_body_id (int): body ID of neuron. Function will look for SWC files in the path
		`swc/[target_name]-[target_body_id].swc where the program is executed from. 
		If the file is not present and `get_neuprint_swc` is true, will attempt to pull the SWC
		file from neuprint directly using `fetch_skeleton`
	input_name (string): name/type of upstream pre-synaptic neuron in neuprint, i.e. "VA6_adPN"
		will use neuprint's `fetch_simple_connections` and `fetch_synapse_connections` to 
		find all the body IDs corresponding to the neuron type and pull in their synapse locations
		if "all_inputs", then will read in the synapses locations for all
		upstream neurons to the target neuron
		if "all_ePN_inputs", then will add all synapses from ePN (lPN, adPN) inputs
		if `None`, then will skip adding synapses
	syn_locs (pd.DataFrame): must include x, y, z columns
		allows one to preload the synaptic locations
	all_input_min_weight (int): if `input_name` is 'all_inputs', then will find all input synapses
		onto the target neuron which have total synapse count in the connection >=`all_input_min_weight` 

	Returns:
	--------
	cell1 (Cell class object): custom class defined in `pop_mc_model`
		can query specific segments of the cell via `cell1.axon[section number](segment value)
		where segment value between 0 and 1
	curr_syns (HOC List): list of HOC Exp2Syn objects (synapses w/ double exponential kinetics)
		added to the cell
	netstim (HOC NetStim object): NetStim object, can be changed to specify time when connections in 
		`netcons` are activated
	netcons (HOC List): list of HOC NetCon objects, one per synapse, that connects the Exp2Syn 
		objects to cell1 and specify the conductance/weight (in microSiemens) of synapses
	num (int): number of synapses added
	pre_body_ids (pd.Series): pandas series of the presynaptic (upstream) identity of each synapse
		should be trivially `input_name` if the neuron type is specified, but if the value of 
		`input_name` is 'all_inputs', then this is useful to know which neuron the synapse comes from
		*index within `pre_body_ids` should correspond to index within `curr_syns`
	pre_nrn_info (pd.DataFrame): rows for each unique input body ID returned by neuprint's 
		`fetch_neurons`, columns include body ID, neuron type, presynaptic neuron total input and
		output synapses, etc. 
		(most relevant for `input_name`='all_inputs')
	'''
	global cell1, curr_syns, netstim, netcons, num

	print(f'visualizing target: {target_name} {target_body_id}, input: {input_name}')

	# instantiate target (post-synaptic) cell
	swc_path = f"swc\\{target_name}-{target_body_id}.swc"
	if f"{target_name}-{target_body_id}.swc" not in os.listdir('swc'):
		if get_neuprint_swc:
			print(f'fetching neuprint swc for {target_name} {target_body_id}')
			from neuprint import fetch_skeleton
			swc_str = fetch_skeleton(body = target_body_id, format = 'swc') # import skeleton
			# scale down x, y, z coordinates and radius by 125 times (convert to microns)
			swc_df = pd.DataFrame([s.split(' ') for s in swc_str.split('\n')][3:],
                      columns = ['sample_num', 'structure_id', 'x_pos', 'y_pos', 'z_pos', 'radius', 'parent_sample']).dropna()
			swc_df[['x_pos', 'y_pos', 'z_pos', 'radius']] = swc_df[['x_pos', 'y_pos', 'z_pos', 'radius']].apply(pd.to_numeric) / 125
			swc_df[['x_pos', 'y_pos', 'z_pos', 'radius']] = swc_df[['x_pos', 'y_pos', 'z_pos', 'radius']].astype(str)

			build_str = ''
			# add back first 3 metadata lines:
			for sub_list in [s.split(' ') for s in swc_str.split('\n')][:3]:
			    build_str += '\n' + ' '.join(sub_list)
			# add coordinate rows w/ proper scaling
			for sub_list in swc_df.values:
			    build_str += '\n' + ' '.join(sub_list)

			with open(swc_path, 'w') as new_swc_path:
				new_swc_path.write(build_str[1:]) # subset to remove new line
			print(f'scaled swc fetched from neuprint and saved at: {swc_path}')
		else:
			print('no SWC found, need to upload into folder `swc`')
			return

	# biophysical parameters from our fits
	R_a = R_a # ohm-cm
	c_m = c_m # uF/cm^2
	g_pas = g_pas # S/cm^2
	e_pas = e_pas # mV one parameter left same in both, we used -55 in our fits
	syn_strength = syn_strength # uS

	cell1 = Cell(swc_path, 0) # first argument is name of swc file, second is a gid'
	cell1.discretize_sections()
	cell1.add_biophysics(R_a, c_m, g_pas, e_pas) # ra, cm, gpas, epas
	try:
		cell1.tree = cell1.trace_tree()
	except Exception as e:
		print('could not trace tree bc of ', e)

	# get number of post-synapses on the target neuron
	try:
		target, r = fetch_neurons(target_body_id)
		target_syn_count = target.post[0]
	except:
		print('likely no internet connection')

	# add all input synaptic site locations
	if input_name: # skip if input_name == None
		if not syn_locs.empty: # if synapse locs are pre-loaded, then use those locations
			pass # synapse locs are pre-loaded
		elif input_name == 'all_inputs': # find all incoming synapses
			min_weight = all_input_min_weight
			print(f'fetching all input synapse locations >{min_weight} synapses')
			syn_locs = fetch_synapse_connections(target_criteria = target_body_id, min_total_weight = min_weight)
		elif input_name == 'all_ePN_inputs':
			min_weight = all_input_min_weight
			print(f'fetching all ePN synapse locations >{min_weight} synapses')
			try:
				syn_locs = fetch_synapse_connections(source_criteria = NC(type='.*lPN.*|.*adPN.*', regex = True),
								target_criteria = target_body_id, min_total_weight = min_weight)
			except Exception as e:
				print(e)
				print('target neuron likely has no inputs which fulfill the target criteria')
				syn_locs = pd.DataFrame()
		else:
			conns = fetch_simple_connections(upstream_criteria = input_name, downstream_criteria = target_body_id)
			pre_bodyIds = conns.bodyId_pre
			syn_locs = fetch_synapse_connections(source_criteria = pre_bodyIds, target_criteria = target_body_id)
		curr_syns, netstim, netcons, num = cell1.add_synapses_xyz(xyz_locs = syn_locs, syn_strength = syn_strength)
		# read info about presynaptic neurons
		if 'bodyId_pre' in syn_locs.columns and syn_locs.shape[0] > 0: 
			pre_body_ids = syn_locs.bodyId_pre
			pre_nrn_info, _ = fetch_neurons(pre_body_ids)
		elif any([val==0 for val in [curr_syns, netstim, netcons, num]]):
			print('no synapses were found to be added')
			curr_syns, netstim, netcons, num, pre_body_ids, pre_nrn_info = None, None, None, None, None, None
		else:
			pre_body_ids = [np.nan] * len(curr_syns)
			pre_nrn_info = pd.DataFrame()
			print('need to pass in presynaptic body IDs')
		try:
			print('adding {} synapses from {} to {}; budget = {}'.format(str(num), input_name, target_name, str(num/target_syn_count)))
		except:
			print('potentially due to internet fault, one of previous variables not defined')
		if target_name == 'local5' and target_body_id == 5813105722:
			# dendrite initial section is axon[195] for un-retraced local5
			num_in_dendr = 0
			for syn in curr_syns:
				if str(syn.get_segment()).partition('(')[0] in [str(val) for val in cell1.axon[195].subtree()]:
					num_in_dendr += 1
			print('proportion synapses in dendrite: {}'.format(str(num_in_dendr/num)))
	else:
		print('not adding synapses')
		curr_syns, netstim, netcons, num, pre_body_ids, pre_nrn_info = None, None, None, None, None, None

	# access a random section
	h.load_file('stdrun.hoc')
	x = h.cvode.active(True)
	if cell1.sec_name=='dend':
		print('accessing sec')
		v_siz = h.Vector().record(cell1.dend[0](0.5)._ref_v)
	elif cell1.sec_name=='axon':
		v_siz = h.Vector().record(cell1.axon[0](0.5)._ref_v)
	elif cell1.sec_name=='soma':
		v_siz = h.Vector().record(cell1.soma[0](0.5)._ref_v)

	return cell1, curr_syns, netstim, netcons, num, pre_body_ids, pre_nrn_info

def time_to_percent_peak(t, v, perc):
	'''
		Given a time trace, voltage trace, and a percentage X%, give the closest time in the
		time trace to X% of the voltage's peak amplitude. 
	'''
	base, peak = min(v), max(v)
	peak_loc = np.where(np.array(v) == peak)[0][0] # index of peak location
	perc_of_peak = perc * (peak - base) + base # value of percent of peak amplitude
	if peak_loc == 0: # safety catch for if v trace is all 0, maybe replace with np.nan, then use nansum 
		return 0
	# interpolate t value where v reaches inputted percent of peak amplitude
	time_of_perc = np.interp(perc_of_peak, np.array(v)[0:peak_loc], np.array(t)[0:peak_loc])
	# subsets portion of trace up to the peak
	#ind_of_perc = np.abs(np.array(v)[0:peak_loc] - perc_of_peak).argmin() 
	#time_of_perc = np.array(t)[ind_of_perc]
	return time_of_perc

def probe_mEPSPs(target_name = 'ML9', target_body_id = 542634516, input_name = 'DP1m_adPN',
				siz_sec = 569, siz_seg = 0.01,
				toPlot = False, lhn_morph_path = '21-05-07_LHN_SWCPointNo_to_NEURON.csv',
				all_input_min_weight = 0,
				transf_freq = 0,
				isopot_dendr_arbor = False, isopot_ax_arbor = False, 
				isopot_prim_ax = False, isopot_prim_dendr = False,
				isopot_whole_nrn = False,
				isopot_Ra_val = 0.001,
				highpot_prim_ax = False, highpot_Ra_val = 10000,
				use_var_time_step = False, fixed_time_step = 0.06,
                IAC_uniform_diam = False, IAC_uniform_med_mean = 'median', 
                IAC_no_pinch = False, IAC_pinch_thresh = 0.14, 
                IAC_eucl_dist = False,
                IAC_3x_len = False,
                cable_props_fn = '21-10-20_LHLN_interarbor_cable_lengths.csv',
                params = {'R_a': 350, # ohm-cm
						  'c_m': 0.6, # uF/cm^2
						  'g_pas': 5.8e-5, # S/cm^2
						  'e_pas': -55, # mV one parameter left same in both, we used -55 in our fits
						  'syn_strength': 5.5e-5}):
	''' given a downstream (target_name) neuron + ID, an upstream neuron, instantiate synapses
		(potentially pulling location info from neuprint) and simulate each mEPSP individually, allowing for
		recording of mEPSP amplitude, kinetics, location of synapse in axon vs dendrite 
		(based on lhn_morph_path CSV NEURON coordinates) etc. 
		NOTE: can be re-engineered to be a lot cleaner

	Args:
	-----
	input_name (str): upstream input name, i.e. "DP1m_adPN", or "all_inputs" to simulate
						mEPSPs for all synaptic inputs, or
						'all_ePN_inputs' to simulate all mEPSPs for ePN inputs
	lhn_morph_path (str): path to file with labelled morphology locations 
						contains 'LHN_SWCPointNo_to_NEURON' then is from Jamie's labels
						contains 'LHN_list_siz_dendr_locs.csv' then is my manual labels (slightly older)
	all_input_min_weight (int): if input_name is 'all_inputs', this sets the minimum number of 
						synapses for connections whose mEPSPs are simulated
	transf_freq (int): frequency at which to measure transfer impedance (=0 is transfer resistance)
	isopot_dendr/ax/prim_ax_arbor (bool): whether to set the axial resistance of the dendritic arbor/
						axonal arbor/primary axon to 0.00001 (isopotential) or 1000 (high potential)
    IAC_uniform_diam (bool): if true, make diameter of inter-arbor cable uniform
                        (all sections inclusive from dendr+ax start to dendr+ax first branch)
    IAC_uniform_med_mean (str): 'median' or 'mean', set IAC diameter uniformly to the med/mean original values
    IAC_eucl_dist (bool): if true, shorten the inter-arbor cable to Euclidean dist between 
                        arbors (by dividing each section along main path by necessary factor)
                        Euclidean distance specified by table at `cable_props_fn`
    IAC_no_pinch (bool): if True, change diameter of segments along IAC w/ diameter <pinch_thresh to pinch_thresh
    IAC_pinch_thresh (float): threshold defining a "pinch point" in micrometers
    cable_props_fn (str): CSV file path providing morphological information on IAC median diameter, for example

	Returns:
	--------
	per_synapse_data - pd.DataFrame: each row contains attributes of a particular synapse, incl:
							ax_v_dendr - string: 'dendritic', 'axonal', or 'unknown'
							ATTRIBUTES: (access via df.attr)
							sum_eff_soma - float: uEPSP @ soma / sum of mEPSPs @ soma
							sum_eff_siz - float: uEPSP @ siz / sum of mEPSPs @ siz
							SIZ = first axonal sec for Jamie's labels
	'''

	# read in inter-arbor cable euclidean+geodesic distance, avg diameter
	cable_df = pd.read_csv(cable_props_fn, index_col = 0)
	mean_diam_um = cable_df.loc[target_body_id, 'mean_diam_um']
	median_diam_um = cable_df.loc[target_body_id, 'median_diam_um']
	euclidean_len_um = cable_df.loc[target_body_id, 'euclidean_len_um']
	geodesic_len_um = cable_df.loc[target_body_id, 'geodesic_len_um']

    # instantiate the cell and input synapses of interest
	cell1, curr_syns, netstim, netcons, num, pre_body_ids, pre_nrn_info = visualize_inputs(target_name = target_name, 
															target_body_id = target_body_id, input_name = input_name,
															all_input_min_weight = all_input_min_weight,
															**params)
	if not netstim or not curr_syns or not netcons:
		print(f'no synapses found for {input_name}, returning empty dataframe')
		return pd.DataFrame()

	# read in proximal axon and dendrite section + extra arbor branches: 
	# TODO: implement morphology read out from testTimingDiffs
	if lhn_morph_path:
		lhn_list = pd.read_csv(lhn_morph_path)
		if 'LHN_list_siz_dendr_locs' in lhn_morph_path:
			try:
				morph_locs = lhn_list.loc[(lhn_list.lhn==target_name) & (lhn_list.lhn_id==target_body_id)] \
								[['dendr_branch_out_sec', 'ax_branch_out_sec', 'siz_sec', 'siz_seg']].iloc[0]
				dendr_branch_out_sec, ax_branch_out_sec, siz_sec = [int(m) for m in morph_locs[0:3]]
				siz_seg = float(morph_locs[3])
			except ValueError:
				print('neuron has unclear dendritic branch out point, using prox dendritic point')
				morph_locs = lhn_list.loc[(lhn_list.lhn==target_name) & (lhn_list.lhn_id==target_body_id)] \
								[['prox_dendr_sec', 'prox_ax_sec', 'siz_sec', 'siz_seg']].iloc[0]
				dendr_branch_out_sec, ax_branch_out_sec, siz_sec, siz_seg = [int(m) for m in morph_locs]
			except:
				print('neuron is unlabeled in morphology file')
		elif 'LHN_SWCPointNo_to_NEURON' in lhn_morph_path:
			try:
				morph_locs = lhn_list.loc[(lhn_list.lhn==target_name) & (lhn_list.body_id==target_body_id)] \
								[['N_dendr_first_branch_sec', 'N_ax_first_branch_sec', \
								  'N_ax_start_sec', 'N_ax_start_seg', \
								  'N_dendr_first_branch_seg', 'N_ax_first_branch_seg']].iloc[0]
				dendr_branch_out_sec, ax_branch_out_sec, \
					siz_sec, _, \
					dendr_branch_out_seg, ax_branch_out_seg = [int(m) for m in morph_locs]
				siz_seg = float(morph_locs[3])
				dendr_branch_out_seg = float(morph_locs[4])
				ax_branch_out_seg = float(morph_locs[5])

				# read in dendr/ax start
				dendr_start_sec = int(lhn_list.loc[(lhn_list.body_id==target_body_id)]['N_dendr_start_sec'].iloc[0])
				dendr_start_seg = float(lhn_list.loc[(lhn_list.body_id==target_body_id)]['N_dendr_start_seg'].iloc[0])
				ax_start_sec = int(lhn_list.loc[(lhn_list.body_id==target_body_id)]['N_ax_start_sec'].iloc[0])
				ax_start_seg = float(lhn_list.loc[(lhn_list.body_id==target_body_id)]['N_ax_start_seg'].iloc[0])

				# read in dendr/ax branch out (redundant, clean up later)
				#dendr_first_branch_sec = 

				# read in extra dendr/ax arbors as list of (sec, seg) tuples (length=0 if empty)
				extra_dendr_sec_seg = eval(lhn_list.loc[(lhn_list.body_id==target_body_id)]['N_extra_dendr_sec_seg'].iloc[0])
				extra_ax_sec_seg = eval(lhn_list.loc[(lhn_list.body_id==target_body_id)]['N_extra_ax_sec_seg'].iloc[0])
			except ValueError: # MODIFY FOR THIS CASE
				print('neuron has unclear dendritic branch out point, using prox dendritic point')
				morph_locs = lhn_list.loc[(lhn_list.lhn==target_name) & (lhn_list.lhn_id==target_body_id)] \
								[['N_dendr_first_branch_sec', 'N_ax_first_branch_sec', 'N_ax_start_sec', 'N_ax_start_seg']].iloc[0]
				dendr_branch_out_sec, ax_branch_out_sec, siz_sec, siz_seg = [int(m) for m in morph_locs]
			except:
				print('neuron is unlabeled in morphology file')

	# read in morphology labels (recording locations at proximal axon, dendrite, etc): 
	# 21-06-30: YES this is redundant with above, but don't want to recode things at this hour
	measure_locs = {'dendr_start': [None, 0.5], 'dendr_first_branch': [None, 0.5], 
                	'ax_start': [None, 0.5], 'ax_first_branch': [None, 0.5],
                	'ax_distal': [0, 0.5], 
                	'soma': [0, 0.5]}
	if lhn_morph_path:
	    lhn_list = pd.read_csv(lhn_morph_path)
	    if 'LHN_SWCPointNo_to_NEURON' in lhn_morph_path:
	        try:
	        	for loc, (sec, seg) in measure_locs.items():
	        		try:
	        			measure_locs[loc][0] = int(lhn_list.loc[(lhn_list.lhn==target_name) & (lhn_list.body_id==target_body_id)]['N_'+loc+'_sec'].iloc[0])
	        			measure_locs[loc][1] = float(lhn_list.loc[(lhn_list.lhn==target_name) & (lhn_list.body_id==target_body_id)]['N_'+loc+'_seg'].iloc[0])
	        		except:
	        			print(f'WARNING: program unsuccessfully tried to find morphology label for {loc}, will use (possibly wrong) default')
	        except:
	            print('neuron is unlabeled in morphology file')
	    else:
	        print('unrecognized morphology file type')
	print('measure_locs acquired:', measure_locs)

	# make relevant dendr/ax (plus extra) arbors isopotential if specified
	if isopot_dendr_arbor:
		print(f'setting dendritic arbor axial resistance to {isopot_Ra_val}')
		for sec in [sec for sec in cell1.axon[dendr_branch_out_sec].subtree() 
						if sec != cell1.axon[dendr_branch_out_sec]]:
			sec.Ra = isopot_Ra_val
		if len(extra_dendr_sec_seg) > 0:
			for e_sec, e_seg in extra_dendr_sec_seg:
				for sec in [sec for sec in cell1.axon[e_sec].subtree() \
						if sec != cell1.axon[e_sec]]:
					sec.Ra = isopot_Ra_val
	if isopot_ax_arbor:
		print(f'setting axonal arbor axial resistance to {isopot_Ra_val}')
		for sec in [sec for sec in cell1.axon[ax_branch_out_sec].subtree() 
						if sec != cell1.axon[ax_branch_out_sec]]:
			sec.Ra = isopot_Ra_val
		if len(extra_ax_sec_seg) > 0:
			for e_sec, e_seg in extra_ax_sec_seg:
				for sec in [sec for sec in cell1.axon[e_sec].subtree() \
						if sec != cell1.axon[e_sec]]:
					sec.Ra = isopot_Ra_val
	if isopot_prim_ax:
		print(f'setting primary axon axial resistance to {isopot_Ra_val}')
		prim_ax_secs = [sec for sec in cell1.axon[siz_sec].subtree() \
						if sec not in cell1.axon[ax_branch_out_sec].subtree()[1:]]
		for sec in prim_ax_secs:
			sec.Ra = isopot_Ra_val
	if highpot_prim_ax: # from 'N_ax_start_sec' = siz_sec to 'N_ax_first_branch_sec'
		print(f'setting primary axon axial resistance to {highpot_Ra_val}')
		prim_ax_secs = [sec for sec in cell1.axon[siz_sec].subtree() \
						if sec not in cell1.axon[ax_branch_out_sec].subtree()[1:]]
		for sec in prim_ax_secs:
			sec.Ra = highpot_Ra_val
	if isopot_prim_dendr:
		print(f'setting primary dendrite axial resistance to {isopot_Ra_val}')
		prim_dendr_secs = [sec for sec in cell1.axon[dendr_start_sec].subtree() \
						   if sec not in cell1.axon[dendr_branch_out_sec].subtree()[1:]]
		for sec in prim_dendr_secs:
			sec.Ra = isopot_Ra_val
	if isopot_whole_nrn:
		print(f'setting entire neuron axial resistance to {isopot_Ra_val}')
		for sec in h.allsec():
			sec.Ra = isopot_Ra_val
	if IAC_uniform_diam:
		# find sections along main path, and extraneous twig sections
		for root, distal in [('ax_start', 'ax_first_branch'), \
							 ('dendr_start', 'dendr_first_branch')]:
			print(f'homogenizing diameter of path from {root} to {distal}')
			main_path, extra_twigs = find_path(cell1.axon[measure_locs[root][0]](measure_locs[root][1]), 
												cell1.axon[measure_locs[distal][0]](measure_locs[distal][1]))
			# homogenize diameter along main path
			for sec in main_path:
				assert IAC_uniform_med_mean=='mean' or IAC_uniform_med_mean=='median', 'must set IAC to median or mean'
				if IAC_uniform_med_mean=='median': 
					sec.diam = median_diam_um
				elif IAC_uniform_med_mean=='mean': 
					sec.diam = mean_diam_um
	if IAC_no_pinch: 
		# find sections along main path, and extraneous twig sections
		for root, distal in [('ax_start', 'ax_first_branch'), \
							 ('dendr_start', 'dendr_first_branch')]:
			print(f'homogenizing diameter of path from {root} to {distal}')
			main_path, extra_twigs = find_path(cell1.axon[measure_locs[root][0]](measure_locs[root][1]), 
												cell1.axon[measure_locs[distal][0]](measure_locs[distal][1]))
			# fix pinch point segments along main path
			for sec in main_path:
				for seg in sec: 
					if seg.diam < IAC_pinch_thresh:
						seg.diam = IAC_pinch_thresh 
						# notably, changing diameters of segments (not sections) also affects segs around it
	if IAC_eucl_dist:
		geodesic_len = 0
		for root, distal in [('ax_start', 'ax_first_branch'), \
							 ('dendr_start', 'dendr_first_branch')]:
			print(f'shortening inter-arbor cable to Euclidean distance')
			main_path, extra_twigs = find_path(cell1.axon[measure_locs[root][0]](measure_locs[root][1]), 
												cell1.axon[measure_locs[distal][0]](measure_locs[distal][1]))
			for sec in main_path:
				sec.L = sec.L * euclidean_len_um / geodesic_len_um
	if IAC_3x_len:
		for root, distal in [('ax_start', 'ax_first_branch'), \
							 ('dendr_start', 'dendr_first_branch')]:
			print(f'tripling length of inter-arbor cable')
			main_path, extra_twigs = find_path(cell1.axon[measure_locs[root][0]](measure_locs[root][1]), 
											   cell1.axon[measure_locs[distal][0]](measure_locs[distal][1]))
			for sec in main_path:
				sec.L = sec.L * 3


	# calculate uEPSP size at siz and soma
	# measure uEPSP for connection at sites in measure_locs
	# activate the stim
	netstim.number = 1
	h.load_file('stdrun.hoc')
	x = h.cvode.active(use_var_time_step)
	if not use_var_time_step:
		h.dt = fixed_time_step # shorten step time if using fixed time step
		# https://www.neuron.yale.edu/neuron/static/py_doc/simctrl/programmatic.html#dt
	v_siz = h.Vector().record(cell1.axon[siz_sec](siz_seg)._ref_v)
	v_soma = h.Vector().record(cell1.axon[0](0.5)._ref_v)
	t = h.Vector().record(h._ref_t)                     				# Time stamp vector
	h.finitialize(-55 * mV)
	h.continuerun(30*ms)
	netstim.number = 0
	uEPSP_siz, uEPSP_soma = max(list(v_siz))+55, max(list(v_soma))+55
	t_trace_uEPSP = list(t)
	v_trace_uEPSP_siz, v_trace_uEPSP_soma = list(v_siz), list(v_soma)

	# first set all synapses to weight 0
	for i in range(len(list(curr_syns))):
		netcons.object(i).weight[0] = 0 # uS, peak conductance change

	# sequentially activate all input synapses
	per_synapse_data = []
	# re-run mEPSP simulation for each synapse
	for i in (range(len(list(curr_syns)))):

		# activate a single synapse
		netcons.object(i).weight[0] = params['syn_strength']
		if i % 100 == 1:
			print("probing mEPSP for synapse " + str(i))
		# measure mEPSP for connection at pSIZ 
		# activate the stim
		netstim.number = 1
		h.load_file('stdrun.hoc')
		#if use_var_time_step:
		x = h.cvode.active(use_var_time_step)
		if not use_var_time_step:
			h.dt = fixed_time_step # shorten step time if using fixed time step
			# https://www.neuron.yale.edu/neuron/static/py_doc/simctrl/programmatic.html#dt
		v_siz = h.Vector().record(cell1.axon[siz_sec](siz_seg)._ref_v) # siz_sec = ax_start
		#v_axon = h.Vector().record(cell1.axon[axon_sec](axon_seg)._ref_v)
		v_soma = h.Vector().record(cell1.axon[0](0.5)._ref_v)
		v_dendr_start = h.Vector().record(cell1.axon[dendr_start_sec](dendr_start_seg)._ref_v)
		v_dendr_first_branch = h.Vector().record(cell1.axon[dendr_branch_out_sec](dendr_branch_out_seg)._ref_v)
		v_ax_first_branch = h.Vector().record(cell1.axon[ax_branch_out_sec](ax_branch_out_seg)._ref_v)
		v_synloc = h.Vector().record(curr_syns.object(i).get_segment()._ref_v)
		t = h.Vector().record(h._ref_t)                     				# Time stamp vector
		h.finitialize(-55 * mV)
		h.continuerun(10*ms)	# initiate run
		netstim.number = 0
		# measure rise time of EPSP at pSIZ
		t_10to90_siz = time_to_percent_peak(t, v_siz, 0.90) - time_to_percent_peak(t, v_siz, 0.10)
		t_10to90_soma = time_to_percent_peak(t, v_soma, 0.90) - time_to_percent_peak(t, v_soma, 0.10)
		t_PNspiketo100_siz = time_to_percent_peak(t, v_siz, 0.999)
		t_PNspiketo100_dendr_first_branch = time_to_percent_peak(t, v_dendr_first_branch, 0.999)
		t_PNspiketo100_ax_first_branch = time_to_percent_peak(t, v_ax_first_branch, 0.999)
		# determine if synapse is in axon or dendrite subtree
		# also checks for the other extra dendr/ax arbors
		ax_vs_dendr_loc = 'unknown'
		syn_section = str(curr_syns.object(i).get_segment()).partition('(')[0]
		if syn_section in [str(val) for val in cell1.axon[dendr_start_sec].subtree()]:
			ax_vs_dendr_loc = 'dendritic'
		elif syn_section in [str(val) for val in cell1.axon[ax_start_sec].subtree()]:
			ax_vs_dendr_loc = 'axonal'
		if len(extra_dendr_sec_seg) > 0:
			for e_sec, e_seg in extra_dendr_sec_seg:
				if syn_section in [str(val) for val in cell1.axon[e_sec].subtree()]:
					ax_vs_dendr_loc = 'dendritic'
		if len(extra_ax_sec_seg) > 0:
			for e_sec, e_seg in extra_ax_sec_seg:
				if syn_section in [str(val) for val in cell1.axon[e_sec].subtree()]:
					ax_vs_dendr_loc = 'axonal'

		# if `input_name` is 'all_inputs', use `pre_nrn_info` to get actual identity of presynaptic input
		syn_input_name = input_name
		if input_name=='all_inputs' or input_name=='all_ePN_inputs':
			syn_input_name = pre_nrn_info.loc[pre_nrn_info.bodyId==pre_body_ids[i]]['type'].iloc[0]

		toAppend = {}
		toAppend.update(target_name = target_name, target_body_id = target_body_id,
							input_name = syn_input_name,
							input_body_id = pre_body_ids[i], 
							synapse_number = i, syn_object = curr_syns.object(i),
							syn_loc = curr_syns.object(i).get_segment(),
							ax_v_dendr = ax_vs_dendr_loc,
							syn_count = len(curr_syns),
							local_diam = curr_syns.object(i).get_segment().diam,
							dist_to_siz = h.distance(cell1.axon[siz_sec](siz_seg), curr_syns.object(i).get_segment()),
							dist_to_soma = h.distance(cell1.axon[0](0.5), curr_syns.object(i).get_segment()),
							dist_to_dendr_first_branch = h.distance(cell1.axon[dendr_branch_out_sec](dendr_branch_out_seg), \
																	curr_syns.object(i).get_segment()),
							mEPSP_siz = max(list(v_siz))+55, mEPSP_soma = max(list(v_soma))+55,
							mEPSP_ax_first_branch = max(list(v_ax_first_branch))+55,
							mEPSP_dendr_start = max(list(v_dendr_start))+55, 
							mEPSP_dendr_first_branch = max(list(v_dendr_first_branch))+55,
							mEPSP_synloc = max(list(v_synloc))+55,
							mEPSP_t10to90_siz = t_10to90_siz,
							mEPSP_t10to90_soma = t_10to90_soma,
							mEPSP_tPNspiketo100_siz = t_PNspiketo100_siz,
							mEPSP_tPNspiketo100_dendr_first_branch = t_PNspiketo100_dendr_first_branch,
							mEPSP_tPNspiketo100_ax_first_branch = t_PNspiketo100_ax_first_branch,
							t_trace = list(t), v_trace_soma = list(v_soma), 
							v_trace_siz = list(v_siz), # same as ax_start
							v_trace_ax_first_branch = list(v_ax_first_branch),
							v_trace_dendr_start = list(v_dendr_start),
							v_trace_dendr_first_branch = list(v_dendr_first_branch),
							v_trace_synloc = list(v_synloc),
							t_trace_uEPSP = t_trace_uEPSP,
							v_trace_uEPSP_soma = v_trace_uEPSP_soma, v_trace_uEPSP_siz = v_trace_uEPSP_siz)
		per_synapse_data.append(toAppend)

		# de-activate the synapse
		netcons.object(i).weight[0] = 0
	per_synapse_data = pd.DataFrame(per_synapse_data)

	# reset all synapse strengths before other tests:
	for i in range(len(list(curr_syns))):
		netcons.object(i).weight[0] = params['syn_strength'] # uS, peak conductance change

	# compute summation efficacy and tack them on as ATTRIBUTES of the dataframe (not columns)
	per_synapse_data.sum_eff_soma = uEPSP_soma/per_synapse_data.mEPSP_soma.sum()
	per_synapse_data.sum_eff_siz = uEPSP_siz/per_synapse_data.mEPSP_siz.sum()
	per_synapse_data.uEPSP_siz = uEPSP_siz
	per_synapse_data.uEPSP_soma = uEPSP_soma
	print(f'soma summation efficacy: {per_synapse_data.sum_eff_soma}')
	print(f'SIZ summation efficacy: {per_synapse_data.sum_eff_siz}')

	# add impedance metrics
	imp_measure_sites = ['soma', 'ax_start', 'dendr_first_branch']
	measure_locs = {'soma': (0, 0.5), 'ax_start': (siz_sec, siz_seg), 
					'dendr_first_branch': (dendr_branch_out_sec, 1)} 
					# dendr_first_branch seg set manually, in the spreadsheet they are all = 1 anyways
	for measure_site in imp_measure_sites:
		# set up Impedance measurement class
		imp = h.Impedance()
		if measure_site=='soma': 
			try: meas_sec_seg = cell1.soma[0](0.5)
			except: meas_sec_seg = cell1.axon[0](0.5)
		else: meas_sec_seg = cell1.axon[measure_locs[measure_site][0]](measure_locs[measure_site][1])
		imp.loc(meas_sec_seg)
		imp.compute(transf_freq)	# starts computing transfer impedance @ freq 

		# iterate through each synapse in the connection & measure impedance
		syn_imp_info = []
		for syn in per_synapse_data.syn_object:
			# find Z_c = transfer impedance from synapse to measure_site, # find Z_i = input impedance at synapse
			curr_transf_imp, curr_input_imp = imp.transfer(syn.get_segment()), imp.input(syn.get_segment())
			# find distance from synapse to measure_site
			curr_distance = h.distance(meas_sec_seg, syn.get_segment())
			# find voltage transfer ratio from synapse to measure_site
			curr_transf_ratio = imp.ratio(syn.get_segment())

			# record individual synapse info
			toAppend = {}
			toAppend.update(Zi = curr_input_imp)
			for meas, meas_val in zip(['dist_to_', 'Zc_to_', 'K_to_'], [curr_distance, curr_transf_imp, curr_transf_ratio]):
				toAppend[meas+measure_site] = meas_val
			syn_imp_info.append(toAppend)
		syn_imp_info = pd.DataFrame(syn_imp_info, index = per_synapse_data.index)
		print(f'adding imp. measurements to {measure_site}')
		# concatenate impedance measure columns (at this meas_site) to existing DataFrame
		for col_name, series in syn_imp_info.iteritems():
			per_synapse_data[col_name] = series
		#per_synapse_data = pd.merge(per_synapse_data, syn_imp_info, on = per_synapse_data.index)

	# add dendritic and axonal and total surface area:
	surf_area_dendr_arbor = 0
	for sec in [sec for sec in cell1.axon[dendr_branch_out_sec].subtree() \
						if sec != cell1.axon[dendr_branch_out_sec]]:
		for seg in sec:
			surf_area_dendr_arbor += seg.area()
	surf_area_ax_arbor = 0
	for sec in [sec for sec in cell1.axon[ax_branch_out_sec].subtree() \
						if sec != cell1.axon[ax_branch_out_sec]]:
		for seg in sec:
			surf_area_ax_arbor += seg.area()
	per_synapse_data['target_surf_area_total'] = cell1.surf_area()
	per_synapse_data['target_surf_area_dendr_arbor'] = surf_area_dendr_arbor
	per_synapse_data['target_surf_area_ax_arbor'] = surf_area_ax_arbor

	if toPlot:
		fig, axs = plt.subplots(nrows = 2, ncols = 2)

		axs[0,0].scatter(per_synapse_data.dist_to_siz, per_synapse_data.mEPSP_siz)
		axs[0,0].set_xlabel('distance to SIZ (um)'), axs[0,0].set_ylabel('mEPSP @ SIZ (mV)')
		axs[0,1].scatter(per_synapse_data.local_diam, per_synapse_data.mEPSP_siz)
		axs[0,1].set_xlabel('local diameter (um)'), axs[0,1].set_ylabel('mEPSP @ SIZ (mV)')
		axs[1,0].scatter(per_synapse_data.dist_to_siz, per_synapse_data.mEPSP_t10to90_siz)
		axs[1,0].set_xlabel('distance to SIZ (um)'), axs[1,0].set_ylabel('t 10 to 90% peak @ SIZ (ms)')
		axs[1,1].scatter(per_synapse_data.local_diam, per_synapse_data.mEPSP_t10to90_siz)
		axs[1,1].set_xlabel('local diameter (um)'), axs[1,1].set_ylabel('t 10 to 90% peak @ SIZ (ms)')

		fig.suptitle(f"{input_name} inputs onto {target_name} {target_body_id}")

	return per_synapse_data