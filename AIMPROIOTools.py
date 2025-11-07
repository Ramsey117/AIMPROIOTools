import os
import re
import numpy as np
import math

global ANGSTROMS_PER_BOHR
ANGSTROMS_PER_BOHR = 0.529177210544 # taken from CODATA

def get_output_file_name(target_directory_path):
	for file_name in os.listdir(target_directory_path):
		if "AIM.sh.o" in file_name and file_name.count('.') == 2:
			return file_name # output file name
	return None

def get_final_energy(output_file_name):
	"""
	Gets the final self-consistent energy from an optimisation / energy calculation. Searches for the second instance of the string "etot" by searching backwards from the end of the file.
	
	Args:
		output_file_name (str) : file to find the final energy in.
	Returns:
		final_energy (float) : in atomic units
	"""
	final_energy = 0.0
	with open(output_file_name, 'r') as output_file:
		lines = output_file.readlines()
		etot_count = 0
		for line in reversed(lines):
			if "etot" in line:
				etot_count += 1
				if etot_count == 2:
					final_energy_match = re.search(r"[-+]?\d*\.\d+", line) # searches for optional +/- sign, then digits, then a decimal point, then more digits

					if final_energy_match:
						final_energy = float(final_energy_match.group())
						break
					else:
						print("No final energy found.")
	return final_energy

def get_final_force(file_path):
	with open(file_path, 'r') as file:
		lines = file.readlines()
		for line in reversed(lines):
			if "OPT " in line:
				parts = line.split()
				try:
					final_force = float(parts[-2])
				except ValueError:
					print(f"Error parsing force in {file_path} on line: {line}")
				break
	return final_force

def find_pristine_directory(pwd):
	"""
	Recursively searches backward in the directory structure for a directory named 'Pristine'.

	Args:
		pwd (str): The starting directory to begin the search.
	Returns:
		str or None: The path to the 'Pristine' directory if found, or None if not found.
	"""
	current_directory = os.path.abspath(pwd)

	while True:
		# Check if "Pristine" exists in the current directory
		pristine_path = os.path.join(current_directory, "Pristine")
		if os.path.isdir(pristine_path):
			return pristine_path
		
		# Move up one level
		parent_directory = os.path.dirname(current_directory)
		
		# If we have reached the root directory, stop the search
		if current_directory == parent_directory:
			break
		
		current_directory = parent_directory

	# Return None if "Pristine" string was not found 
	return None

def count_atoms(file_path):
	with open(file_path, 'r') as file_to_be_counted:
		lines = file_to_be_counted.readlines()
	num_atoms = 0
	positions_section_flag = False
	for line in lines:
		if "begin" in line and "{positions}" in line:
			positions_section_flag = True
		elif "end{positions}" in line:
			break
		elif positions_section_flag and line.strip() and not line.strip().startswith("!"): ########## NEED TO ACCOUNT FOR COMMENTED OUT ATOMS
			num_atoms += 1
	return num_atoms

def get_sampling(file_path):
	with open(file_path, 'r') as f:
		lines = f.readlines()
	sampling = [0, 0, 0]
	for line in lines:
		if "sampling{" in line:
			match = re.search(r'grid=\s*(\d+)\s+(\d+)\s+(\d+)', line)
			if match:
				sampling = [int(match.group(1)), int(match.group(2)), int(match.group(3))]
				break
	return sampling

"""
OLD - parses from dat-file section
def get_initial_lat_consts(file_path):
	with open(file_path, 'r') as f:
		lines = f.readlines()
	lat_consts = []
	primitive_lat_consts = []
	for line in lines:
		if "lattice{" in line:
			floats_match = re.search(r'params=([-\d. ]+)', line) # Match 'params=' followed by one or more decimal numbers
			if floats_match:
				float_strings = floats_match.group(1).split()
				primitive_lat_consts = [float(x) for x in float_strings]
				break
	return lat_consts
"""
def get_lattice(file_path, space='real', output='constants'):
	"""
	Retrieves either the lattice constants or lattice vectors (real or reciprocal) from AIMPRO's output.

	Args:
		file_path (str): Path to the AIMPRO output file.
		space (str): 'real' or 'reciprocal' to specify which lattice type to extract.
		output (str): 'constants' to return [a, b, c]; 'vectors' to return 3x3 matrix of lattice vectors.

	Returns:
		list of floats or 3x3 numpy array: Lattice constants [a, b, c] or a 3x3 array of lattice vectors.
	"""
	if space not in ['real', 'reciprocal']:
		raise ValueError("Argument 'space' must be either 'real' or 'reciprocal'.")
	if output not in ['constants', 'vectors']:
		raise ValueError("Argument 'output' must be either 'constants' or 'vectors'.")

	with open(file_path, 'r') as f:
		lines = f.readlines()

	for i, line in enumerate(lines):
		if "unit vectors : real space, reciprocal space" in line:
			index = slice(0, 3) if space == 'real' else slice(-3, None)
			a_vec = np.array(re.findall(r"[-+]?\d*\.\d+", lines[i+1])[index], dtype=float)
			b_vec = np.array(re.findall(r"[-+]?\d*\.\d+", lines[i+2])[index], dtype=float)
			c_vec = np.array(re.findall(r"[-+]?\d*\.\d+", lines[i+3])[index], dtype=float)
			if output == 'constants':
				return [np.linalg.norm(a_vec), np.linalg.norm(b_vec), np.linalg.norm(c_vec)]
			else:  # output == 'vectors'
				return np.vstack([a_vec, b_vec, c_vec]) 
	
	raise RuntimeError("Lattice information not found in file.")

def atom_coords_intp2angstrom_npArray(atom_coords_intp_npArray, file_path):
	lattice_vectors = get_lattice(file_path, space='real', output='vectors')
	return (lattice_vectors.T @ atom_coords_intp_npArray) * ANGSTROMS_PER_BOHR

def atom_coords_angstrom2intp_npArray(atom_coords_angstrom_npArray, file_path):
	lattice_vectors = get_lattice(file_path, space='real', output='vectors')
	return np.linalg.inv(lattice_vectors.T) @ (atom_coords_angstrom_npArray/ANGSTROMS_PER_BOHR) # the reciprocal vectors are not used because they involve a factor of 2pi to convert to momentum space.

def parse_species(file_path):
	"""
	Searches a dat file or AIMPRO output file for the species in the system.
	
	Args:
		file_path (str): path to dat/output file.
	Returns:
		species_list (list of str): each str is the elemental symbol of the associated species.
	"""
	def parse_species_symbol(line):
		"""
		Searches a species line for the elemental symbol for the species with regex.
		
		Args:
			line (str): species line
		Returns:
			species_symbol (str): elemental symbol for the species
		"""
		match = re.search(r'pot=\d+-([A-Z][a-z]*)', line) # Regular expression to capture only the element symbol after "pot="
		if match:
			species_symbol = match.group(1)
			return species_symbol
		else:
			raise ValueError("Input string format is incorrect")
	with open(file_path, 'r') as f:
		lines = f.readlines()
	species_list = []
	species_section_flag = False
	for line in lines:
		if "begin{hghpseudo}" in line:
			species_section_flag = True
		elif "end{hghpseudo}" in line:
			break
		elif species_section_flag and line.strip(): # ignores blank lines 
			species_symbol = parse_species_symbol(line)
			species_list.append(species_symbol)
	return species_list # ordered in the order of the species section of the output file.

class Atom():
	def __init__(self,input_index):
		self.index = input_index
		self.species = "init"
		self.coords_intp = np.array([-10.0, -10.0, -10.0])
		self.coords_angstrom = np.array([-10.0, -10.0, -10.0])
		self.nearest_neighbours = []
		self.RMS_angular_strain_percentage = -10.0 # this is used in the angular strain heatmap
		
	def add_nearest_neighbour(self, neighbour_atom):
		self.nearest_neighbours.append(neighbour_atom)

def parse_atom_data(file_path,species_list):
	"""
	Parses a system's atomic positions. Uses prior knowledge of the species in the species list to assign an elemental system to each atom. BE AWARE - this assumes dat/output file is in hexagonal int-p units, should probably amend to work for both int-p and atomic units.
	
	Args:
		file_path (str): path to the dat/output file being searched.
		species_list (list of str): list of the elemental symbols of the species in the system, in the same order as they are defined in the file.
	Returns:
		system (list of Atom objects): the system represented as a list of Atom objects. This includes the atom's index, species and position in the coords (int-p / atomic) it is represented in the file.
	"""
	with open(file_path, 'r') as f:
		lines = f.readlines()
	system = []
	atoms_flag = False
	intp_flag = False
	for line in lines:
		if "begin" in line and "{positions}" in line:
			atoms_flag = True
			if "[int-p]" in line:
				intp_flag = True
		elif "end{positions}" in line:
			break
		elif atoms_flag:
			parts = line.strip().split()
			if len(parts) < 2 or "!" in line: # this should help ignore blank lines
				continue
			atom = Atom(int(parts[0])) # argument is atom's index
			if intp_flag:
				atom.coords_intp = np.array([float(parts[-3]), float(parts[-2]), float(parts[-1])])
				atom.coords_angstrom = atom_coords_intp2angstrom_npArray(atom.coords_intp, file_path)
			else:
				atom.coords_angstrom = ANGSTROMS_PER_BOHR*np.array([float(parts[-3]), float(parts[-2]), float(parts[-1])])
				atom.coords_intp = angstrom2intp_npArray(atom.coords_angstrom, file_path)
			atom.species = species_list[int(parts[1])-1] # -1 is due to different counting bases
			system.append(atom)
	return system

def get_Rlast_lines(output_file_name):
	"""
	Runs gres -Rlast command on the specified AIMPRO output file. Opens the associated dat.AIM.sh.o.... file and returns its lines.
	
	Args:
        output_file_name (str): The starting directory to begin the search.
	Returns:
        optimised_dat_file_lines (list of str): lines of the dat.AIM.sh.o.... file.
	"""
	if output_file_name is None:
		return []
	else:
		if f"dat.{os.path.basename(output_file_name)}" in os.listdir():
			os.remove(f"dat.{os.path.basename(output_file_name)}") # removes a previous usage's dat.AIM.sh.o
		subprocess.run(["/home/njpg/bin/gres", "-Rlast", output_file_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
		with open(f"dat.{os.path.basename(output_file_name)}", "r") as optimised_dat_file:
			optimised_dat_file_lines = optimised_dat_file.readlines()
		#os.remove(f"dat.{os.path.basename(output_file_name)}")
		return optimised_dat_file_lines


