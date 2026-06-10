import os
import re
import numpy as np
import math
import subprocess
import bz2
from pathlib import Path
from consts import ANG_PER_BOHR, ANG_PER_BOHR, RY_PER_HA, EV_PER_HA

def get_output_file_name(target_directory_path):
	for file_name in os.listdir(target_directory_path):
		if "AIM.sh.o" not in file_name:
			continue
		
		if file_name.endswith(".bz2") and file_name.count('.') == 3:
			return file_name
		if not file_name.endswith(".bz2") and file_name.count('.') == 2:
			return file_name
	return None

def smart_open(file_path, mode="r", **kwargs):
	"""
	Opens a file, automatically handling plain text or .bz2 compressed files.
	
	Args:
		file_path (str or Path) : path to the file to open.
		mode (str, optional)    : file mode, same as built-in open() ('r', 'w', 'a', etc.).
								  Defaults to 'r'.
		**kwargs                : additional keyword arguments passed to open() or bz2.open(
	Returns:
		file-like object        : opened in the requested mode.
	"""
	path = Path(file_path)
	
	if path.suffix == ".bz2":
		if 'b' not in mode and 't' not in mode:
			mode = mode+'t'
		return bz2.open(path, mode, **kwargs)
	else:
		return open(path, mode, **kwargs)

def get_final_energy(output_file_name):
	"""
	Gets the final self-consistent energy from an optimisation / energy calculation. Searches for the second instance of the string "etot" by searching backwards from the end of the file.
	
	Args:
		output_file_name (str) : file to find the final energy in.
	Returns:
		final_energy (float)   : in atomic units
	"""
	final_energy = 0.0
	with smart_open(output_file_name, 'r') as output_file:
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
	with smart_open(file_path, 'r') as output_file:
		lines = output_file.readlines()
	
	for line in reversed(lines):
		if "OPT " in line:
			parts = line.split()
			try:
				final_force = float(parts[-2])
			except ValueError:
				print(f"Error parsing force in {file_path} on line: {line}")
			break
	return final_force

def get_net_charge(file_path):
	with smart_open(file_path, 'r') as file:
		lines = file.readlines()
	for line in lines:
		if "Charged Unit cell, charge :" in line:
			net_charge = float(line.split(":")[1])
			return net_charge
	return 0.0

def find_pristine_directory(pwd):
	"""
	Recursively searches backward in the directory structure for a directory named 'Pristine'.
	Args:
		pwd (str or Path) : The starting directory to begin the search.
	Returns:
		Path or None: The path to the 'Pristine' directory if found, or None if not found.
	"""
	current_directory = Path(pwd).resolve()  # resolves symlinks, absolute path

	while True:
		pristine_path = current_directory / "Pristine"
		if pristine_path.is_dir():
			return pristine_path

		parent_directory = current_directory.parent

		# If we have reached the root directory, stop the search
		if current_directory == parent_directory:
			break

		current_directory = parent_directory
	return None

def find_supercell_dimensions_from_path(pwd):
	"""
	Searches the directory path for a supercell specification of the form '/nxm/'.
	Args:
		pwd (str or Path)       : directory path containing the supercell specification.
	Returns:
		tuple(int, int) or None : (n, m) supercell size if found, otherwise None.
	"""
	path = Path(pwd)

	# search each directory component
	for part in path.parts:
		match = re.fullmatch(r"(\d+)x(\d+)", part)
		if match:
			return int(match.group(1)), int(match.group(2))
	return None

def count_atoms(file_path):
	with smart_open(file_path, 'r') as file_to_be_counted:
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
	with smart_open(file_path, 'r') as f:
		lines = f.readlines()
	sampling = [0, 0, 0]
	for line in lines:
		if "sampling{" in line:
			match = re.search(r'grid=\s*(\d+)\s+(\d+)\s+(\d+)', line)
			if match:
				sampling = [int(match.group(1)), int(match.group(2)), int(match.group(3))]
				break
	return sampling

def get_ecut(file_path,unit=None):
	"""
	Retrieves the cut off energy used in the AIMPRO calculation from the AIMPRO output. Looks for a specific phrase 'Energy cutoff (Hartrees)' in this file.
	Args:
		file_path (str) : path to AIMPRO output file
		unit (str)      : unit of energy that is wished to be returned (Options: 'Ha', 'Ry', 'eV')
	Returns:
		ecut (float)    : cut off energy in the unit specified in the input.
	"""
	if unit == None:
		raise ValueError("Please specify a unit of energy ('Ha', 'Ry', 'eV')")
	if unit not in ["Ha", "Ry", "eV"]:
		raise ValueError("Error: Inappropriate unit label.")
	
	with smart_open(file_path, 'r') as f:
		lines = f.readlines()
	ecut_Ha = 0.0
	for line in lines:
		if "Energy cutoff (Hartrees)" in line:
			ecut_Ha = float(line.split()[0])
			break
	if ecut_Ha == 0.0:
		raise ValueError(f"Ecut not found in {file_path}")
	
	if unit == "Ry":
		return ecut_Ha*RY_PER_HA
	elif unit == "eV":
		return ecut_Ha*EV_PER_HA
	else:
		return ecut_Ha

"""
OLD - parses from dat-file section
def get_initial_lat_consts(file_path):
	with smart_open(file_path, 'r') as f:
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
def get_lattice(*, file_path, space, output, unit, desired_iter, repeats=[0,0,0]): # * indicates that all the arguments are required, no defaults are given, and the user must specify them keyword
	"""
	Retrieves either the lattice constants or lattice vectors (real or reciprocal) from AIMPRO's output.
	Args:
		file_path    (str or Path)   : Path to the AIMPRO output file.
		space        (str)           : 'real' or 'reciprocal' to specify which lattice type to extract. Reciprocal only works with initial lattice parameters.
		output       (str)           : 'constants' to return [a, b, c]; 'vectors' to return 3x3 matrix of lattice vectors.
		unit         (str)           : 'Bohr', 'Ang'.
		desired_iter (str, or int)   : 'initial' or 'final' or integer describing the desired **lattice** optimisation iteration number. If there is no lattice optimisation in the run, the initial parameters will be obtained and the user is notified via the terminal.
		repeats (list of three ints) : the number of desired repeats [na, nb, nc] corresponding to each of the real space lattice directions.
	Returns:
		list of floats : Lattice constants [a, b, c]
		OR
		3x3 numpy array: a 3x3 array of lattice vectors arranged in a vertical stack - i.e. the rows correspond to each lattice vector.
	"""
	if desired_iter == 'initial':
		desired_iter = 0
	elif desired_iter == 'final':
		desired_iter = -1

	if space not in ['real', 'reciprocal']:
		raise ValueError("Argument 'space' must be either 'real' or 'reciprocal'.")
	if output not in ['constants', 'vectors']:
		raise ValueError("Argument 'output' must be either 'constants' or 'vectors'.")
	if unit not in ['Bohr', 'Ang']:
		raise ValueError("Argument 'unit' must be either 'Bohr' or 'Ang'.")
	if not isinstance(desired_iter, int) or desired_iter < -1:
		raise ValueError("Argument 'desired_iter' must be either 'initial' or 'final' or -1, 0, or positive integer that is the optimisation iteration index.")
	if desired_iter == -1 and space != 'real':
		raise ValueError("Reciprocal lattice parameters are only available for desired_iter = 'initial'. This functionality may be added in the future.")
	if isinstance(file_path, (str, Path)):
		if Path(file_path).name == "dat":
			raise ValueError("get_lattice() does not currently work with dat files.")
	else:
		raise TypeError(f"{file_path} needs to be str or Path object type.")

	
	# This dictionary avoids a complicated logical circuit for the return.
	UNIT_FACTORS = {
		("Bohr", "real")       : 1.0,
		("Bohr", "reciprocal") : 1.0,
		("Ang",  "real")       : ANG_PER_BOHR,
		("Ang",  "reciprocal") : 1/ANG_PER_BOHR
	}
	conversion_factor = UNIT_FACTORS[(unit,space)]


	with smart_open(file_path, 'r') as f:
		lines = f.readlines()

	a_vec = None
	b_vec = None
	c_vec = None
	if desired_iter != 0:
		iteration_counter = -1 # initiated as -1 because aucl is present in the initial statement of position data for lattice opt calculations.
		for i, line in enumerate(lines): # this loop will iterate until the end of the file, meaning that these are the final, optimised lattice parameters.
			if line.strip().startswith(f"aucl:"):
				iteration_counter += 1 # after this line the iteration_counter is identical to AIMPRO's description of the iteration number.
				a_values = line.split(":", 1)[1].split()
				a_vec = np.array(a_values[:3], dtype=float)
				b_values = lines[i+1].split(":", 1)[1].split()
				b_vec = np.array(b_values[:3], dtype=float)
				c_values = lines[i+2].split(":", 1)[1].split()
				c_vec = np.array(c_values[:3], dtype=float)
				if iteration_counter == desired_iter:
					break
		if a_vec is None or b_vec is None or c_vec is None: # optimised lattice parameters have not been found. Notify the user that the initials are to be used.
			print(f"No lattice optimisation in {file_path}. Using initial lattice params.")
		if desired_iter > iteration_counter and iteration_counter != -1: # The desired optimisation iteration index is larger than the maximum in the calculation. Raise an error because the user has likely input the wrong number. -1 is excluded to avoid this happening when there is no lattice optimisation i.e. when there is geometry optimisation but not lattice optimisation.
			raise ValueError(f"desired_iter={desired_iter} has not been found in {file_path}. There are only {iteration_counter} lattice optimisation iterations.")

	if a_vec is None or b_vec is None or c_vec is None: # this will happen if desired_iter = 'initial' is selected or if optimised lattice params are not found.
		for i, line in enumerate(lines):
			if "unit vectors : real space, reciprocal space" in line:
				index = slice(0, 3) if space == 'real' else slice(-3, None)
				a_vec = np.array(re.findall(r"[-+]?\d*\.\d+", lines[i+1])[index], dtype=float)
				b_vec = np.array(re.findall(r"[-+]?\d*\.\d+", lines[i+2])[index], dtype=float)
				c_vec = np.array(re.findall(r"[-+]?\d*\.\d+", lines[i+3])[index], dtype=float)
				break

	if a_vec is None or b_vec is None or c_vec is None: # there are no initial lattice parameters. Likely not an AIMPRO output file.
		raise RuntimeError(f"Lattice information not found in {file_path}. Are you sure this is an AIMPRO output file?")

	lattice_vectors = np.vstack([a_vec*conversion_factor, b_vec*conversion_factor, c_vec*conversion_factor])

	if (repeats[0] + repeats[1] + repeats[2]) != 0:
		if space == 'real':
			repeat_matrix = np.diag([repeats[0]+1, repeats[1]+1, repeats[2]+1]) # +1 is to convert between zero-based counting vs. one-based counting.
		elif space == 'reciprocal':
			repeat_matrix = np.diag([1/(repeats[0]+1), 1/(repeats[1]+1), 1/(repeats[2]+1)]) # +1 is to convert between zero-based counting vs. one-based counting.
		lattice_vectors = repeat_matrix @ lattice_vectors

	if output == 'constants':
		return [np.linalg.norm(lattice_vectors[0]), np.linalg.norm(lattice_vectors[1]), np.linalg.norm(lattice_vectors[2])]
	else:  # output == 'vectors'
		return lattice_vectors

def atom_coords_intp2angstrom_npArray(atom_coords_intp_npArray, file_path):
	# depreciated
	lattice_vectors = get_lattice(file_path, space='real', output='vectors')
	return (lattice_vectors.T @ atom_coords_intp_npArray) * ANG_PER_BOHR

def atom_coords_angstrom2intp_npArray(atom_coords_angstrom_npArray, file_path):
	# depreciated
	lattice_vectors = get_lattice(file_path, space='real', output='vectors')
	return np.linalg.inv(lattice_vectors.T) @ (atom_coords_angstrom_npArray/ANG_PER_BOHR) # the reciprocal vectors are not used because they involve a factor of 2pi to convert to momentum space.

def parse_species(file_path):
	"""
	Searches a dat file or AIMPRO output file for the species in the system.
	
	Args:
		file_path (str)           : path to dat/output file.
	Returns:
		species_list (list of str): each str is the elemental symbol of the associated species.
	"""
	def parse_species_symbol(line):
		"""
		Searches a species line for the elemental symbol for the species with regex.
		
		Args:
			line (str)          : species line
		Returns:
			species_symbol (str): elemental symbol for the species
		"""
		match = re.search(r'pot=\d+-([A-Z][a-z]*)', line) # Regular expression to capture only the element symbol after "pot="
		if match:
			species_symbol = match.group(1)
			return species_symbol
		else:
			raise ValueError("Input string format is incorrect")
	with smart_open(file_path, 'r') as f:
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

class Atom:
	def __init__(self,input_index):
		self.index = input_index
		self.species = "init"
		self.coords_intp = np.array([-10.0, -10.0, -10.0])
		self.coords_angstrom = np.array([-10.0, -10.0, -10.0])
		self.nearest_neighbours = []
		self.RMS_angular_strain_percentage = -10.0 # this is used in the angular strain heatmap
		
	def add_nearest_neighbour(self, neighbour_atom):
		self.nearest_neighbours.append(neighbour_atom)

	def print(self):
		print(f"{self.index} {self.species} {self.coords_intp} {self.coords_angstrom}")

def parse_atom_data(*,file_path,species_list,desired_iter,repeats=[0,0,0]):
	"""
	Parses a system's atomic positions. Uses prior knowledge of the species in the species list to assign an elemental system to each atom.

	Args:
		file_path    (str or Path)          : path to the AIMPRO output file being searched.
		species_list (list of str)          : list of the elemental symbols of the species in the system, in the same order as they are defined in the file.
		desired_iter (str or int)           : 'initial' or 'final' or integer describing the index of the desired **geometry** optimisation iteration for which the atom positions will be parsed. Works with co-optimised lattice parameters.
		repeats      (list of three ints)   : the number of desired repeats [na, nb, nc] corresponding to each of the real space lattice directions.
	Returns:
		system       (list of Atom objects) : the system represented as a list of Atom objects. This includes the atom's index, species and position in the coords (int-p / atomic) it is represented in the file.
	"""
	with smart_open(file_path, 'r') as f:
		lines = f.readlines()
	system = []
	atoms_flag = False
	intp_flag = False

	# Convert initial and final to their corresponding iteration integers
	if desired_iter == 'initial':
		desired_iter = 0
	elif desired_iter == 'final':
		desired_iter = -1

	# argument handling
	if not isinstance(desired_iter, int) or desired_iter < -1:
		raise ValueError("Argument 'desired_iter' must be either 'initial' or 'final' or -1, 0, or positive integer that is the optimisation iteration index.")
	if isinstance(file_path, (str, Path)):
		if Path(file_path).name == "dat":
			raise ValueError("get_lattice() does not currently work with dat files.")
	else:
		raise TypeError(f"{file_path} needs to be str or Path object type.")

	if len(repeats) != 3 or not all(isinstance(x, int) for x in repeats):
		raise ValueError("repeats needs to be a 3 item list of integers")

	# Detect whether the job includes optimisation.
	optimisation_flag = False
	for line in lines:
		if "optimise option chosen" in line:
			optimisation_flag = True
			break
		elif "Processing Lattice ..." in line: # this is past the point where the optimisation flag will have been identified.
			break

	# Safeguard against requests for positions of specific iterations for calculation with no optimisation.
	if not optimisation_flag and desired_iter not in [0,-1]:
		raise ValueError(f"No optimisation in calculation. desired_iter={desired_iter} does not exist.")

	lattice_vectors_Ang = get_lattice(file_path=file_path, space='real', output='vectors', unit='Ang', desired_iter=desired_iter)

	# For
	if desired_iter != 0:
		position_start = None
		if desired_iter == -1: # search backwards for the last positions block
			initial_line_index     = len(lines) - 1
			termination_line_index = -1
			direction              = -1
			num_pos_block_to_pass  = 0 # this is the number of atom position blocks it will see before parsing the atom data and stopping.
		else:
			initial_line_index     = 0
			termination_line_index = len(lines)
			direction              = 1
			num_pos_block_to_pass  = desired_iter

		# Search the file to find the position at which the desired atom positions block is.
		num_pos_blocks_passed = 0
		for i in range(initial_line_index, termination_line_index, direction):
			if "positions of atoms :" in lines[i] and num_pos_blocks_passed < num_pos_block_to_pass:
				num_pos_blocks_passed += 1
			elif "positions of atoms :" in lines[i] and num_pos_blocks_passed == num_pos_block_to_pass:
				if lines[i].split('=')[1].strip() == 'intp':
					intp_flag = True
				elif lines[i].split('=')[1].strip() == 'intc':
					raise ValueError("parse_atom_data() is not currently set-up to accept int-c coordinates yet.")
				position_start = i + 1   # first data line
				break
		if num_pos_blocks_passed < num_pos_block_to_pass:
			raise ValueError(f"There are not {desired_iter} optimisation iterations in {file_path}. Only {num_pos_blocks_passed} positions blocks found.")

		if position_start is not None:
			for line in lines[position_start:]:
				line = line.strip()
				if not line: # on first blank line
					break

				parts = line.split()

				atom = Atom(int(parts[0]))
				atom.species = species_list[int(parts[1])-1] # -1 is due to different counting bases
				if intp_flag:
					atom.coords_intp = np.array([float(parts[-3]), float(parts[-2]), float(parts[-1])])
					atom.coords_angstrom = atom.coords_intp @ lattice_vectors_Ang
				else: # atomic units
					atom.coords_angstrom = ANG_PER_BOHR*np.array([float(parts[-3]), float(parts[-2]), float(parts[-1])])
					atom.coords_intp = atom.coords_angstrom @ np.linalg.inv(lattice_vectors_Ang)
				system.append(atom)

	else:
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
					atom.coords_angstrom = atom.coords_intp @ lattice_vectors_Ang
				else:
					atom.coords_angstrom = ANG_PER_BOHR*np.array([float(parts[-3]), float(parts[-2]), float(parts[-1])])
					atom.coords_intp = atom.coords_angstrom @ np.linalg.inv(lattice_vectors_Ang)
				atom.species = species_list[int(parts[1])-1] # -1 is due to different counting bases
				system.append(atom)

	if len(system) == 0:
		raise RuntimeError("No positions blocks found. Are you sure this is an AIMPRO output file?") # unlike the lattice parameters, this *should* find the positions, even in calculations where there is no optimisation.

	if (repeats[0] + repeats[1] + repeats[2]) != 0:
		primitive_system = system.copy()
		for a in range(repeats[0] + 1):
			for b in range(repeats[1] + 1):
				for c in range(repeats[2] + 1):
					if a == 0 and b == 0 and c == 0:
						continue

					print(f"a={a}, b={b}, c={c}")
					for atom in primitive_system:
						new_atom = Atom(len(system)+1)
						new_atom.species = atom.species
						new_atom.coords_angstrom = (
							atom.coords_angstrom + a*lattice_vectors_Ang[0] + b*lattice_vectors_Ang[1] + c*lattice_vectors_Ang[2])
						new_atom.coords_intp = atom.coords_intp.copy()
						new_atom.coords_intp[0] += a
						new_atom.coords_intp[1] += b
						new_atom.coords_intp[2] += c

						system.append(new_atom)
	return system


def get_Rlast_lines(output_file_name):
	"""
	Runs gres -Rlast command on the specified AIMPRO output file. Opens the associated dat.AIM.sh.o.... file and returns its lines.
	
	Args:
        output_file_name (str)                : The starting directory to begin the search.
	Returns:
        optimised_dat_file_lines (list of str): lines of the dat.AIM.sh.o.... file.
	"""
	if output_file_name is None:
		return []
	else:
		if f"dat.{os.path.basename(output_file_name)}" in os.listdir():
			os.remove(f"dat.{os.path.basename(output_file_name)}") # removes a previous usage's dat.AIM.sh.o
		subprocess.run(["/home/njpg/bin/gres", "-Rlast", output_file_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
		with smart_open(f"dat.{os.path.basename(output_file_name)}", "r") as optimised_dat_file:
			optimised_dat_file_lines = optimised_dat_file.readlines()
		#os.remove(f"dat.{os.path.basename(output_file_name)}")
		return optimised_dat_file_lines

def get_bandstructure_eV(output_file_path, what=None):
	if what == None:
		raise TypeError("No bandstructure element defined.")

	initial_pwd = os.getcwd()
	os.chdir(os.path.dirname(output_file_path))
	result = subprocess.run(["perl", "/home/njpg/bin/PlotBYK", output_file_path], capture_output=True, text=True, check=True)

	stdout = result.stdout
	stderr = result.stderr

	lines = result.stderr.splitlines()

	for line in lines:
		if "min empty band value =" in line:
			match = re.search(r"=\s*([-+]?\d*\.\d+)", line)
			if match:
				CBM_energy_eV = float(match.group(1))
			else:
				print("No match found.")
		elif "max occupied band value =" in line:
			match = re.search(r"=\s*([-+]?\d*\.\d+)", line)
			if match:
				VBM_energy_eV = float(match.group(1))
			else:
				print("No match found.")
	os.chdir(initial_pwd)
	if what == "bandgap":
		bandgap_eV = CBM_energy_eV - VBM_energy_eV
		return bandgap_eV
	elif what == "VBM":
		return VBM_energy_eV
	elif what == "CBM":
		return CBM_energy_eV
	else:
		raise ValueError("Unrecognised bandstructure element.")
