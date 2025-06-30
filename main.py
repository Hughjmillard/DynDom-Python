import os
os.environ["OMP_NUM_THREADS"] = '1'
import math
import FileMngr
import numpy as np
import gemmi
from difflib import SequenceMatcher
from Clusterer import Clusterer
from Protein import Protein
from scipy.spatial.transform import Rotation
from arrow_visualization import ArrowGenerator
from clustering_logger import ClusteringLogger

class Engine:
    def __init__(self, input_path, output_path, pdb_1, chain_1, pdb_2, chain_2, k_means_n_init=1, k_means_max_iter=500,
                 window=5, domain=20, ratio=1.0, atoms="backbone"):
        # The atoms to be used for the program
        if atoms == "backbone":
            self.atoms_to_use = ["N", "CA", "C"]
        elif atoms == "ca":
            self.atoms_to_use = ["CA"]
        self.input_path = input_path
        self.output_path = output_path
        self.k_means_n_init = k_means_n_init
        self.k_means_max_iter = k_means_max_iter
        self.window = window
        self.domain = domain
        self.ratio = ratio
        # Initialise proteins
        self.protein_1: Protein = Protein(input_path, pdb_1, chain_1, self.atoms_to_use)
        self.protein_2: Protein = Protein(input_path, pdb_2, chain_2, self.atoms_to_use)
        # Used to move Protein 2 into the coordinate space of Protein 1
        self.fitted_protein_2 = None
        # A object containing superimposition information between the whole protein chains
        self.chain_superimpose_result = None
        # List of arrays of gemmi.SupResult objects containing superimposition information between each residue
        self.slide_window_superimpose_results = []
        # Array of translation vectors
        self.translation_vecs = None
        # Array of (3 x 3) rotation matrices
        self.rotation_mats = None
        # List of rotation vectors created from rotation_mats
        self.rotation_vecs = None
        # List of unit vectors of the rotation vectors
        self.unit_vectors = None
        # Bending residues indices
        self.bending_residues_indices = {}
        # Scaling of rotation vectors
        self.scaling = 0.0
        # Clusterer object
        self.clusterer = None
        self.shared_logger = None

    def run(self):
        """
        Runs the entire program. All operations happen based on Protein 1 conforming to Protein 2.
        :return:
        """
        running = True
        res_ind_1, res_ind_2, res_names_1, res_names_2 = self.get_atoms_indices()
        self.check_sequence_identity(res_names_1, res_names_2)
        self.protein_1.utilised_residues_indices = res_ind_1
        self.protein_2.utilised_residues_indices = res_ind_2

        while running:
            self.get_slide_window_start_end_indices()
            # Superimposes Protein 2 onto Protein 1, so they are in the same "space". The superimposed protein will be
            # called Protein S.
            self.superimpose_chains()
            # self.print_chains_superimposed_result()
            # Slides a window over Protein 1's residues and Protein S's residues and superimposes them in each window.
            # This gets the superimposition results for each sliding window.
            self.sliding_window_superimpose_residues()
            # Obtain the rotation matrices from the superposition results
            self.get_transformations()
            # Convert 3x3 rotation matrices to rotation vectors
            self.convert_rot_mats_to_vecs()
            # Write the rotation vectors to a pdb file
            # FileMngr.write_rotation_vec_to_pdb(self.protein_1.id, self.protein_1.slide_window_residues,
            #                                    self.protein_1.slide_window_residues_indices, self.rotation_vecs)
            FileMngr.write_rotation_vec_to_pdb(self.output_path, self.protein_1, self.protein_2.name, self.protein_2.chain_param, self.rotation_vecs)
            # self.print_unit_vectors()
            # Initialise the Clustering class
            if not self.shared_logger:
                output_base = f"{self.protein_1.name}_{self.protein_1.chain_param}_{self.protein_2.name}_{self.protein_2.chain_param}"
                log_dir = os.path.join(self.output_path, output_base, "clustering_logs")
                protein_name = f"{self.protein_1.name}_{self.protein_1.chain_param}_{self.protein_2.name}_{self.protein_2.chain_param}"
                self.shared_logger = ClusteringLogger(log_dir, protein_name)
            
                # Create clusterer with logging disabled
            self.clusterer: Clusterer = Clusterer(
                self.protein_1, self.protein_2, self.rotation_vecs, self.atoms_to_use,
                self.k_means_n_init, self.k_means_max_iter, self.window, self.domain, self.ratio,
                enable_logging=False  # Disable internal logger
            )
            # Assign shared logger
            self.clusterer.logger = self.shared_logger
            
            self.clusterer.cluster()
            # print(f"Cluster status = {self.clusterer.clusterer_status}")
            if self.clusterer.clusterer_status == -1:
                if self.shared_logger:
                    self.shared_logger.log_window_change(
                        old_window=self.window,
                        new_window=self.window + 2,
                        reason="Clustering failed, trying larger window"
                    )
                
                self.window += 2
                continue
            screws = self.determine_domains_screw_axis()
            for screw in screws:
                print("Screw :")
                print(screw)
            self.determine_bending_residues()

            FileMngr.write_final_output_pdb(
                self.output_path,
                self.protein_1,
                self.fitted_protein_2,
                self.protein_2.name,
                self.protein_2.chain_param,
                self.protein_2.utilised_residues_indices
            )
            FileMngr.write_final_output_pml(self.output_path, self.protein_1, self.protein_2.name,
                                            self.protein_2.chain_param, self.clusterer.domains,
                                            self.clusterer.fixed_domain, self.bending_residues_indices, self.window)
            FileMngr.write_w5_info_file(self.output_path, self.protein_1.name, self.protein_1.chain_param,
                                        self.protein_2.name, self.protein_2.chain_param, self.window, self.domain,
                                        self.ratio, self.atoms_to_use, self.clusterer.domains, self.clusterer.fixed_domain)
            running = False
        return True

    def check_sequence_identity(self, res_names_1, res_names_2):
        """
        Checks the sequence identity of the 2 protein chains.
        :return:
        """
        correct_hit = 0

        for i in range(len(res_names_1)):
            if res_names_1[i] == res_names_2[i]:
                correct_hit += 1
        similarity = correct_hit / len(res_names_1)
        if similarity < 0.4:
            raise ValueError("Sequence Identity less than 40%")

    def get_atoms_indices(self):
        """
        Get the indices of residues containing backbone atoms from the proteins. The residues can only be used if the sequence number of the
        residue the atom is in of Protein 1 and 2 are the same. This is to account for protein chains of different lengths.
        The sequence number of residues in the PDB format usually start at 1, but there will be polymer residues in
        protein chains that start at a number above 1 or even a negative value. Returns 4 lists.
        The indices of the residues with utilised CA atoms from proteins 1 and 2, and the residue names.
        The returned lists are the same size.

        Extra residues: Atoms in residues with sequence numbers 2 and above will be used
        Protein 1 Residue Sequence Numbers [ *  * * * 2 3 ...]
        Protein 2 Residue Sequence Numbers [-2 -1 0 1 2 3 ...]

        Missing residues: Atoms in residues with sequence number 5, 6, 7 will not be used
        Protein 1 Residue Sequence Numbers [ 1 2 3 4 * * * 8 9 ...]
        Protein 2 Residue Sequence Numbers [ 1 2 3 4 5 6 7 8 9 ...]

        :return utilised_res_ind_1: A list of indices of the residues in the
        :return utilised_res_ind_2:
        :return res_names_1:
        :return res_names_2:
        """
        # Get the protein 1 and 2 polymer chains. This excludes residues which are only water.
        protein_1_polymer = self.protein_1.get_polymer()
        protein_2_polymer = self.protein_2.get_polymer()
        # Get the length of the polymers
        protein_1_size = len(protein_1_polymer)
        protein_2_size = len(protein_2_polymer)
        # The indices to iterate the polymers
        index_1 = 0
        index_2 = 0
        # Stores the index of the residues located in the chains
        utilised_res_ind_1 = []
        utilised_res_ind_2 = []
        # Stores the names of the residues
        res_names_1 = []
        res_names_2 = []
        # Stops iterating at the end of the shorter protein chain
        while index_1 < protein_1_size and index_2 < protein_2_size:
            # Get the sequence id number of the residue in the proteins
            res_num_1 = protein_1_polymer[index_1].seqid.num
            res_num_2 = protein_2_polymer[index_2].seqid.num
            # If the residue number of the
            if res_num_1 < res_num_2:
                index_1 += 1
                continue
            if res_num_2 < res_num_1:
                index_2 += 1
                continue
            n_found_1, ca_found_1, c_found_1 = False, False, False
            n_found_2, ca_found_2, c_found_2 = False, False, False
            for a in protein_1_polymer[index_1]:
                if a.name == self.atoms_to_use[0] and not n_found_1:
                    n_found_1 = True
                elif a.name == self.atoms_to_use[1] and not ca_found_1:
                    ca_found_1 = True
                elif a.name == self.atoms_to_use[2] and not c_found_1:
                    c_found_1 = True
                if n_found_1 and ca_found_1 and c_found_1:
                    break
            for a in protein_2_polymer[index_2]:
                if a.name == self.atoms_to_use[0] and not n_found_2:
                    n_found_2 = True
                elif a.name == self.atoms_to_use[1] and not ca_found_2:
                    ca_found_2 = True
                elif a.name == self.atoms_to_use[2] and not c_found_2:
                    c_found_2 = True
                if n_found_2 and ca_found_2 and c_found_2:
                    break
            if n_found_1 and ca_found_1 and c_found_1 and n_found_2 and ca_found_2 and c_found_2:
                utilised_res_ind_1.append(index_1)
                utilised_res_ind_2.append(index_2)
                res_names_1.append(protein_1_polymer[index_1].name)
                res_names_2.append(protein_2_polymer[index_2].name)
            index_1 += 1
            index_2 += 1

        return utilised_res_ind_1, utilised_res_ind_2, res_names_1, res_names_2

    def get_slide_window_start_end_indices(self):
        """
        Gets the start and end indices of the sliding windows using the index of the middle residue in each window.
        If window size is 5 and the length of the utilised residues is 431 (Index 430), the start and end indices are
        2 and 428 respectively.
        :return:
        """
        window_size = self.window
        # The index of the middle residue in the sliding window
        window_mid_index = start_index = (window_size - 1) // 2
        utilised_res_length = len(self.protein_1.utilised_residues_indices)
        final_index = utilised_res_length + 1 - (window_size - window_mid_index)
        self.protein_1.slide_window_residues_indices = (start_index, final_index)
        self.protein_2.slide_window_residues_indices = (start_index, final_index)

    def superimpose_chains(self):
        """
        Superimposes the entire chain of Protein 2 onto Protein 1 using the backbone atoms to get Protein 2 in the same
        coordinate space of Protein 1. Saves the Protein 2 chain after transformation into fitting_protein_2.
        """
        self.fitted_protein_2 = self.protein_2.get_chain()
        fitted_protein_polymer = self.fitted_protein_2.get_polymer()
        ptype = fitted_protein_polymer.check_polymer_type()

        self.chain_superimpose_result: gemmi.SupResult = gemmi.calculate_superposition(self.protein_1.get_polymer(),
                                                                                       fitted_protein_polymer,
                                                                                       ptype, gemmi.SupSelect.MainChain)
        fitted_protein_polymer.transform_pos_and_adp(self.chain_superimpose_result.transform)
        # P
        print(f'rmsd of whole protein best fit: {self.chain_superimpose_result.rmsd:.3f}A')

    def sliding_window_superimpose_residues(self):
        """
        Slides a window of a specified size over both protein chains. The backbone atoms of the residues in the
        sliding windows are superimposed to get the superposition results of the middle residue of each sliding window.
        :return:
        """
        self.slide_window_superimpose_results = []
        window_size = self.window
        fitting_protein_polymer = self.protein_1.get_polymer()  # The protein that will superimpose onto Protein S.
        target_protein_polymer = self.fitted_protein_2.get_polymer()
        # target_protein_polymer = self.protein_1.get_polymer()
        # For each window, get the residues' backbone atoms and get the superposition results.
        for r in range(len(self.protein_1.utilised_residues_indices) - window_size + 1):  # Number of iterations of the window
            # Initialise an empty list that stores the backbone atoms' coordinates of the target protein
            # (The fitted protein 2)
            target_protein_polymer_atoms_pos = []
            # Initialise an empty list that stores the backbone atoms' coordinates of the fitting protein
            # (The protein 1 that will fit onto the fitted protein 2)
            fitting_protein_polymer_atoms_pos = []
            # For each residue in the window,
            for i in range(r, r+window_size):
                fitting_index = self.protein_1.utilised_residues_indices[i]
                target_index = self.protein_2.utilised_residues_indices[i]
                # For each utilised atom in the residue,
                for a in self.atoms_to_use:
                    # Append the atom coordinate to the lists
                    fitting_protein_polymer_atoms_pos.append(fitting_protein_polymer[fitting_index][a][0].pos)
                    target_protein_polymer_atoms_pos.append(target_protein_polymer[target_index][a][0].pos)
                    # target_protein_polymer_atoms_pos.append(target_protein_polymer[i].sole_atom(a).pos)
            # Superimpose the atoms and append the result to list
            self.slide_window_superimpose_results.append(gemmi.superpose_positions(target_protein_polymer_atoms_pos,
                                                                                   fitting_protein_polymer_atoms_pos))

    def get_transformations(self):
        """
        Each window has a gemmi.SupResult object containing information of the superimposition which is representative
        of the middle residue of the window. This function is used to extract the numerical data of the objects
        for KMeans clustering. Specifically, the rotation matrix.
        :return:
        """
        # Initialise the numpy arrays
        # Declare an empty array of rotation matrices
        self.rotation_mats = np.empty(shape=[len(self.slide_window_superimpose_results), 3, 3])
        # Declare an empty array of translation vectors
        self.translation_vecs = np.empty(shape=[len(self.slide_window_superimpose_results), 3])
        # For each SupResult
        for i in range(len(self.slide_window_superimpose_results)):
            # Get the rotation matrix and translation vector of the residue
            rot_mat: gemmi.Mat33 = self.slide_window_superimpose_results[i].transform.mat
            trans_vec: gemmi.Vec3 = self.slide_window_superimpose_results[i].transform.vec
            self.rotation_mats[i] = np.asarray(rot_mat.tolist())
            self.translation_vecs[i] = np.asarray(trans_vec.tolist())

    def convert_rot_mats_to_vecs(self):
        """
        Convert rotation matrices to rotation vectors, angles in degrees, and unit rotation vectors.
        :return:
        """
        self.rotation_vecs = Rotation.from_matrix(self.rotation_mats).as_rotvec(degrees=True)
        self.unit_vectors = self.rotation_vecs / np.linalg.norm(self.rotation_vecs)

    def determine_domains_screw_axis(self):
        """
        Calculates the screw axis of the dynamic domains.
        :return:
        """
        # Get the transformations of the fixed domain of Protein 2 fitting to fixed domain of Protein 1
        fixed_domain_r: gemmi.SupResult = self.get_fixed_domain_transformations()
        # Get the slide chain of Protein 1
        original_protein_1_slide_chain: gemmi.ResidueSpan = self.protein_1.get_slide_window_residues()
        # Get the slide chain of Protein 2 and fit to Protein 1 relative to the fixed domain.
        transformed_protein_2_slide_chain: gemmi.ResidueSpan = self.protein_2.get_slide_window_residues()
        transformed_protein_2_slide_chain.transform_pos_and_adp(fixed_domain_r.transform)

        self.clusterer.domains[self.clusterer.fixed_domain].rmsd = fixed_domain_r.rmsd

        domain_screw_axes = []

        # Go through each dynamic domain
        for domain in self.clusterer.domains:
            # Skip the fixed domain
            if domain.domain_id == self.clusterer.fixed_domain:
                continue
            """
            First thing to do is to fit Protein 1 onto the transformed Protein 2 so that we can find the displacement
            vectors between the initial positions of the dynamic domain atoms of Protein 1 and the final positions of 
            the dynamic domain atoms of Protein 1
            """
            # Prepare an empty Protein chain that will contain the residues of the dynamic domain of Protein 1
            original_protein_1_domain_chain: gemmi.Chain = gemmi.Chain(self.protein_1.chain_param)
            # Prepare an empty Protein chain that will contain the residues of the dynamic domain of Protein 1 fitted on 2
            transformed_protein_1_domain_chain: gemmi.Chain = gemmi.Chain(self.protein_1.chain_param)
            # Prepare an empty Protein chain that will contain the residues of the dynamic domain of Protein 2 fitted on 1
            transformed_protein_2_domain_chain: gemmi.Chain = gemmi.Chain(self.protein_2.chain_param)
            # Go through each segment of the current dynamic domain
            for segment in domain.segments:
                # Add the residues
                for i in range(segment[0], segment[1] + 1):
                    original_protein_1_domain_chain.add_residue(original_protein_1_slide_chain[i])
                    transformed_protein_1_domain_chain.add_residue(original_protein_1_slide_chain[i])
                    transformed_protein_2_domain_chain.add_residue(transformed_protein_2_slide_chain[i])

            # print(f"Ori before : {original_protein_1_slide_chain[4].sole_atom('CA').pos}")
            # print(f"Trans before : {transformed_protein_1_domain_chain[4].sole_atom('CA').pos}")

            # r: gemmi.SupResult = gemmi.superpose_positions(original_atoms, transformed_atoms)
            transformed_protein_1_domain_polymer: gemmi.ResidueSpan = transformed_protein_1_domain_chain.get_polymer()
            transformed_protein_2_domain_polymer: gemmi.ResidueSpan = transformed_protein_2_domain_chain.get_polymer()
            ptype = transformed_protein_1_domain_polymer.check_polymer_type()
            # Fit Protein 1 dynamic domain onto Transformed Protein 2 dynamic domain
            r: gemmi.SupResult = gemmi.calculate_superposition(transformed_protein_2_domain_polymer,
                                                               transformed_protein_1_domain_polymer,
                                                               ptype, gemmi.SupSelect.MainChain)
            domain.rmsd = r.rmsd
            # Transform the domain chain
            transformed_protein_1_domain_polymer.transform_pos_and_adp(r.transform)
            # print(f"Ori after : {original_protein_1_slide_chain[4].sole_atom('CA').pos}")
            # print(f"Trans after : {transformed_protein_1_domain_chain[4].sole_atom('CA').pos}")
            # Get the rotation matrix of the domain transformation and convert to rotation vector
            rot_vec = Rotation.from_matrix(np.asarray(r.transform.mat.tolist())).as_rotvec(degrees=True)
            # Get the unit vector of the rotation vector
            unit_rot_vec = rot_vec / math.sqrt(np.sum(rot_vec**2))
            rot_angle = np.linalg.norm(rot_vec)

            """
            Assuming the dynamic domains move as a rigid body, the translations for all atoms would be the same. No need
            to use all residues. Just 4 atoms/residues is fine. Either backbone atoms from 4 residues or just 4 atoms.
            """
            original_atom_coords = np.asarray(original_protein_1_domain_chain[0][self.atoms_to_use[0]][0].pos.tolist())
            transformed_atom_coords = np.asarray(transformed_protein_1_domain_chain[0][self.atoms_to_use[0]][0].pos.tolist())
            for i in range(1, 4):
                original_atom_coords = np.vstack((original_atom_coords, np.asarray(original_protein_1_domain_chain[i][self.atoms_to_use[0]][0].pos.tolist())))
                transformed_atom_coords = np.vstack((transformed_atom_coords, np.asarray(transformed_protein_1_domain_chain[i][self.atoms_to_use[0]][0].pos.tolist())))
            original_atom_coords = np.mean(original_atom_coords, axis=0)
            transformed_atom_coords = np.mean(transformed_atom_coords, axis=0)

            # Calculate the displacement vector
            disp_vec = transformed_atom_coords - original_atom_coords

            translation_component_value = np.sum(disp_vec * unit_rot_vec)
            parallel_translation = unit_rot_vec * translation_component_value
            # Calculate difference between displacement and parallel translations to get rotational parts
            rotational_part = disp_vec - parallel_translation
            # Calculate the amplitude of rotation
            rotation_amplitude = math.sqrt(np.sum(rotational_part**2))
            # Calculate unit vector in direction of rotational part
            unit_rotational_part = rotational_part/rotation_amplitude

            # Calculate vector in direction from atoms to axis
            cross_prod_axis = np.cross(unit_rot_vec, unit_rotational_part)
            h_tan = 2*math.tan(0.5*rot_angle)
            atoms_to_axis_direction = (rotation_amplitude*cross_prod_axis)/h_tan

            point_on_axis = original_atom_coords + (0.5 * rotational_part) - atoms_to_axis_direction
            domain.rot_angle = rot_angle
            domain.disp_vec = disp_vec
            domain.point_on_axis = point_on_axis
            domain.screw_axis = unit_rot_vec
            domain.translation = translation_component_value
            domain_screw_axes.append((unit_rot_vec, rot_angle, point_on_axis))

        return domain_screw_axes

    def determine_bending_residues(self):
        """
        Determine the bending residues between each fixed-dynamic domain pair.
        The rotation vectors of the fixed and dynamic domains are used to calculate the mean and
        :return:
        """
        fixed_domain = self.clusterer.domains[self.clusterer.fixed_domain]
        fixed_domain_segments = fixed_domain.segments
        print("Fixed segments:", fixed_domain_segments)
        mid_point = (self.window - 1) // 2
        fixed_domain_rot_vecs = self.rotation_vecs[fixed_domain_segments[0][0]+mid_point:fixed_domain_segments[0][1]+mid_point]

        # Get the rotation vectors of the fixed domain
        for i in range(1, fixed_domain_segments.shape[0]):
            rot_vecs = self.rotation_vecs[fixed_domain_segments[i][0]+mid_point:fixed_domain_segments[i][1]+mid_point]
            fixed_domain_rot_vecs = np.append(fixed_domain_rot_vecs, rot_vecs, axis=0)

        # Calculate mean of the fixed domain rotation vectors
        fixed_domain_mean = np.mean(fixed_domain_rot_vecs, axis=0)
        # fixed_domain_std = np.std(fixed_domain_rot_vecs)
        # Center the rotation vectors of the fixed domain
        fixed_domain_centered_vecs = fixed_domain_rot_vecs - fixed_domain_mean
        # Get the covariance and inversed covariance matrices
        fixed_domain_covar = np.cov(fixed_domain_centered_vecs.T)
        fixed_domain_inv_covar = np.linalg.inv(fixed_domain_covar)
        # print(f"Fixed Domain Mean = {fixed_domain_mean}")
        # print(f"Fixed Domain STD = {fixed_domain_std}")
        # print(f"Fixed Domain Var = {fixed_domain_var}")
        # print(f"Fixed Domain Covariance = \n{fixed_domain_covar}")
        # print(f"Fixed Domain Inverse Covariance = \n{fixed_domain_inv_covar}")

        # For each dynamic domain,
        for domain in self.clusterer.domains:
            # Ignore the fixed domain
            if domain.domain_id == self.clusterer.fixed_domain:
                continue
            bend_res_set = set()
            self.bending_residues_indices[domain.domain_id] = []
            print(f"Domain {domain.domain_id}")
            dyn_dom_segments = domain.segments
            print(f"Dyn Segments:", dyn_dom_segments)
            dyn_dom_rot_vecs = self.rotation_vecs[dyn_dom_segments[0][0]+mid_point:dyn_dom_segments[0][1]+mid_point]

            for i in range(1, dyn_dom_segments.shape[0]):
                rot_vecs = self.rotation_vecs[dyn_dom_segments[i][0]+mid_point:dyn_dom_segments[i][1]+mid_point]
                dyn_dom_rot_vecs = np.append(dyn_dom_rot_vecs, rot_vecs, axis=0)

            # Just like the fixed domain, calculate the mean, centered rotation vectors, covariance and inverse
            # covariance matrices
            dyn_dom_mean = np.mean(dyn_dom_rot_vecs, axis=0)
            # dyn_dom_std = np.std(dyn_dom_rot_vecs)
            dyn_dom_centered_vecs = dyn_dom_rot_vecs - dyn_dom_mean
            dyn_dom_covar = np.cov(dyn_dom_centered_vecs.T)
            dyn_dom_inv_covar = np.linalg.inv(dyn_dom_covar)

            # print(f"Dyn Domain Mean = {dyn_dom_mean}")
            # print(f"Dyn Domain STD = {dyn_dom_std}")
            # print(f"Dyn Domain Var = {dyn_dom_var}")
            # print(f"Dyn Domain Covariance = \n{dyn_dom_covar}")
            # print(f"Dyn Domain Inverse Covariance = \n{dyn_dom_inv_covar}")

            # Calculate the indices of the previous and next residues for each segment of the dynamic domain.
            dyn_dom_prev_indices = dyn_dom_segments[:, 0] - 1
            dyn_dom_next_indices = dyn_dom_segments[:, 1] + 1

            # Get the indices of the fixed domain segments that connects the fixed domain to the dynamic domain.
            # 1D Array of booleans where True means next index after fixed domain segment is dyn dom segment.
            fixed_next_is_dyn = np.isin(fixed_domain.segments[:, 1], dyn_dom_prev_indices)
            fixed_next_is_dyn_ind = np.where(fixed_next_is_dyn)[0]
            # 1D Array of booleans where True means previous index before fixed domain segment is dyn dom segment.
            fixed_prev_is_dyn = np.in1d(fixed_domain.segments[:, 0], dyn_dom_next_indices)
            fixed_prev_is_dyn_ind = np.where(fixed_prev_is_dyn)[0]

            # Get the indices of the dynamic domain segments that connects the dynamic domain to the fixed domain.
            # 1D Array of booleans where True means next index after dyn dom segment is fixed domain segment.
            dyn_next_is_fixed = np.in1d(dyn_dom_next_indices, fixed_domain.segments[:, 0])
            dyn_next_is_fixed_ind = np.where(dyn_next_is_fixed)[0]
            # 1D Array of booleans where True means previous index before dyn dom segment is fixed domain segment.
            dyn_prev_is_fixed = np.in1d(dyn_dom_prev_indices, fixed_domain.segments[:, 1])
            dyn_prev_is_fixed_ind = np.where(dyn_prev_is_fixed)[0]

            # print("Fixed domain segment next index is dyn dom segment: ")
            # print("Fixed next is dyn:", fixed_next_is_dyn)
            # print(fixed_next_is_dyn_ind)
            # print("Fixed domain segment prev index is dyn dom segment: ")
            # print("Fixed prev is dyn:", fixed_prev_is_dyn)
            # print(fixed_prev_is_dyn_ind)

            # print("Dyn dom segment next index is fixed domain segment: ")
            # print(dyn_next_is_fixed)
            # print(dyn_next_is_fixed_ind)
            # print("Dyn dom segment prev index is fixed domain segment: ")
            # print(dyn_prev_is_fixed)
            # print(dyn_prev_is_fixed_ind)

            p = 4.6

            # Go backwards through the fixed domain residues of the segments
            # print("Backward fixed")
            for segment_ind in fixed_next_is_dyn_ind:
                segment = fixed_domain.segments[segment_ind]
                for i in range(segment[1], segment[0] - 1, -1):
                    centered_vec = self.rotation_vecs[i+mid_point] - fixed_domain_mean
                    q_value = centered_vec @ fixed_domain_inv_covar @ centered_vec
                    print(q_value)
                    # print("Backward Fixed Q Value =", q_value)
                    # print(i, q_value)
                    if q_value > p:
                        # print("Hit")
                        bend_res_set.add(i)
                    else:
                        break

            # Go forwards through the dyn dom residues of the segments
            # print("Forward Dyn")
            for segment_ind in dyn_prev_is_fixed_ind:
                segment = domain.segments[segment_ind]
                # print(segment)
                for i in range(segment[1], segment[0] + 1):
                    centered_vec = self.rotation_vecs[i+mid_point] - dyn_dom_mean
                    q_value = centered_vec @ dyn_dom_inv_covar @ centered_vec
                    # print("Forward Dyn Q Value =", q_value)
                    # print(i, q_value)
                    if q_value > p:
                        # print("Hit")
                        bend_res_set.add(i)
                    else:
                        break

            # Go forwards through the fixed domain residues of the segments
            # print("Forward Fixed")
            for segment_ind in fixed_prev_is_dyn_ind:
                segment = fixed_domain.segments[segment_ind]
                # print(segment)
                for i in range(segment[0], segment[1] + 1):
                    centered_vec = self.rotation_vecs[i+mid_point] - fixed_domain_mean
                    q_value = centered_vec @ fixed_domain_inv_covar @ centered_vec
                    # print("Forward Fixed Q Value =", q_value)
                    # print(i, q_value)
                    if q_value > p:
                        # print("Hit")
                        bend_res_set.add(i)
                    else:
                        break

            # Go backwards through the dyn dom residues of the segments
            # print("Backward Dyn")
            for segment_ind in dyn_next_is_fixed_ind:
                segment = domain.segments[segment_ind]
                # print(segment)
                for i in range(segment[1], segment[0] - 1, -1):
                    centered_vec = self.rotation_vecs[i+mid_point] - dyn_dom_mean
                    q_value = centered_vec @ dyn_dom_inv_covar @ centered_vec
                    # print("Backward Dyn Q Value =", q_value)
                    # if q_value < p or q_value > 1-p:
                    # print(i, q_value)
                    if q_value > p:
                        # print("Hit")
                        bend_res_set.add(i)
                    else:
                        break

            bend_res_set = list(bend_res_set)
            bend_res_set.sort()
            domain.bend_res = bend_res_set
            self.bending_residues_indices[domain.domain_id] = bend_res_set

    def get_fixed_domain_transformations(self):
        """
        Get the transformation of the fixed domain of protein 2 to protein 1
        :return:
        """
        slide_window_1 = self.protein_1.get_slide_window_residues()
        slide_window_2 = self.protein_2.get_slide_window_residues()
        coords_1 = []
        coords_2 = []
        fixed_domain_id = self.clusterer.fixed_domain
        for s in self.clusterer.domains[fixed_domain_id].segments:
            for i in range(s[0], s[1] + 1):
                for a in self.atoms_to_use:
                    coords_1.append(slide_window_1[i][a][0].pos)
                    coords_2.append(slide_window_2[i][a][0].pos)

        r: gemmi.SupResult = gemmi.superpose_positions(coords_1, coords_2)
        return r

    def get_arrow_coords(self):
        polymer = self.protein_1.get_polymer()
        util_res = self.protein_1.utilised_residues_indices
        domains = self.clusterer.domains
        for d in range(len(domains)):
            if domains[d].domain_id == self.clusterer.fixed_domain:
                continue
            domain = domains[d]
            coords = []
            for s in range(domain.segments.shape[0]):
                for i in range(domain.segments[s][0], domain.segments[s][1]+1):
                    index = util_res[i]
                    res = polymer[index]
                    for a in res:
                        if a.name in self.atoms_to_use:
                            coords.append(a.pos.tolist())

    def print_chains_superimposed_result(self):
        # A gemmi.SupResult object containing superimposition information between 2 chains
        print(f"RMSD =                  {self.chain_superimpose_result.rmsd}")
        print(f"Count =                 {self.chain_superimpose_result.count}")
        print(f"Center 1 =              {self.chain_superimpose_result.center1}")
        print(f"Center 2 =              {self.chain_superimpose_result.center2}")
        print(f"Translation Vector =    {self.chain_superimpose_result.transform.vec}")
        print(f"Rotation Matrix =       {self.chain_superimpose_result.transform.mat}")
        rot_vec = Rotation.from_matrix(self.chain_superimpose_result.transform.mat).as_rotvec(degrees=True)
        print(f"Rotation Vector =       {rot_vec}")

    def print_slide_window_superimpose_results(self, n=None):
        if n is None or n > len(self.slide_window_superimpose_results):
            n = len(self.slide_window_superimpose_results)
        print(f"slide_window_superimpose_result size = {len(self.slide_window_superimpose_results)}")
        for i in range(n):
            item: gemmi.SupResult = self.slide_window_superimpose_results[i]
            print(f"RMSD =                  {item.rmsd}")
            print(f"Count =                 {item.count}")
            print(f"Center 1 =              {item.center1}")
            print(f"Center 2 =              {item.center2}")
            print(f"Translation Vector =    {item.transform.vec}")
            print(f"Rotation Matrix =       {item.transform.mat}")

    def print_slide_window_residue_indices(self):
        print(f"slide_window_residue_indices = {self.protein_1.slide_window_residues_indices}")

    def print_rotation_matrices(self, n=None):
        if n is None or n > self.rotation_mats.shape[0]:
            n = self.rotation_mats.shape[0]
        print(f"rotation_mats shape = {self.rotation_mats.shape}")
        print(f"rotation_mats[0:{n}] = {self.rotation_mats[0:n]}")

    #
    # def print_chains_superimposed_result(self):
    #     print(math.sqrt(sum([r.rmsd for r in self.chain_superimpose_result])))
    #     # for r in self.chain_superimpose_result:
    #     #     print(f"RMSD =                  {r.rmsd}")
    #     #     print(f"Count =                 {r.count}")
    #     #     print(f"Center 1 =              {r.center1}")
    #     #     print(f"Center 2 =              {r.center2}")
    #     #     print(f"Translation Vector =    {r.transform.vec}")
    #     #     print(f"Rotation Matrix =       {r.transform.mat}")


files_dict = FileMngr.read_command_file()
# Read param file to get parameters ( window size, domain size, ratio, etc. )
param_dict = FileMngr.read_param_file()
# Initialise Engine object
engine = Engine(input_path=files_dict["input_path"], output_path=files_dict["output_path"],
                pdb_1=files_dict["filename1"], chain_1=files_dict["chain1id"],
                pdb_2=files_dict["filename2"], chain_2=files_dict["chain2id"],
                k_means_n_init=param_dict["k_means_n_init"], k_means_max_iter=param_dict["k_means_max_iter"],
                window=param_dict["window"], domain=param_dict["domain"],
                ratio=param_dict["ratio"], atoms=param_dict["atoms"])
# Run the Engine
engine.run()


# Generate arrow visualizations if analysis was successful
try:
    # Construct file paths based on your current naming convention
    output_base = f"{files_dict['filename1'].split('.')[0]}_{files_dict['chain1id']}_{files_dict['filename2'].split('.')[0]}_{files_dict['chain2id']}"
    w5_info_file = os.path.join(files_dict["output_path"], output_base, f"{output_base}.w5_info")
    
    # Check if w5_info file exists (indicates successful analysis)
    if os.path.exists(w5_info_file):
        print("\n" + "="*50)
        print("GENERATING ARROW VISUALIZATIONS")
        print("="*50)
        
        # Create arrow generator
        arrow_gen = ArrowGenerator(w5_info_file, os.path.join(files_dict["output_path"], output_base,output_base))
        
        # Generate all arrow files
        results = arrow_gen.generate_all_arrows()
        
        if results:
            print(f"\nArrow files created:")
            print(f"  PDB: {results['pdb']}")
            print(f"  PyMOL: {results['pymol']}")
            print(f"\nTo view: pymol {results['pymol']}")
        else:
            print("No arrows generated (no moving domains found)")
            
    else:
        print("No w5_info file found - analysis may have failed")
        
except Exception as e:
    print(f"Arrow generation failed: {e}")
    # Don't crash the main program if arrow generation fails
    import traceback
    traceback.print_exc()