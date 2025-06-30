import numpy as np
import gemmi
from scipy.spatial.transform import Rotation

class SingleResidueBackbonePlusCBCalculator:
    def __init__(self, protein_1, protein_2, window_size=5):
        """
        Calculate rotation for a single residue using N, CA, C, CB atoms (like original DynDom)
        
        Args:
            protein_1: First protein structure
            protein_2: Second protein structure (should be fitted_protein_2)
            window_size: Size of sliding window (default 5)
        """
        self.protein_1 = protein_1
        self.protein_2 = protein_2
        self.window_size = window_size
        self.atoms_to_use = ["N", "CA", "C", "CB"]  # Original DynDom atom set
        
    def find_residue_without_nearby_glycine(self, utilised_residues_indices, min_distance=3):
        """
        Find a residue that doesn't have glycine within the sliding window
        
        Args:
            utilised_residues_indices: List of residue indices to check
            min_distance: Minimum distance from glycine (default 3)
            
        Returns:
            suitable_residue_index: Index of a suitable residue, or None if not found
        """
        polymer = self.protein_1.get_polymer()
        half_window = self.window_size // 2
        
        # Get positions that can be analyzed with sliding window
        start_analysis = half_window
        end_analysis = len(utilised_residues_indices) - half_window
        
        for i in range(start_analysis, end_analysis):
            residue_idx = utilised_residues_indices[i]
            residue = polymer[residue_idx]
            
            # Check if this residue and its window neighbors are not glycine
            window_start = max(0, i - half_window)
            window_end = min(len(utilised_residues_indices), i + half_window + 1)
            
            has_glycine_in_window = False
            for j in range(window_start, window_end):
                check_residue_idx = utilised_residues_indices[j]
                check_residue = polymer[check_residue_idx]
                if check_residue.name == "GLY":
                    has_glycine_in_window = True
                    break
            
            if not has_glycine_in_window:
                print(f"Found suitable residue at index {i} (residue {residue.seqid.num}, {residue.name})")
                print(f"Window covers residues {window_start} to {window_end-1}")
                return i, residue_idx
                
        return None, None
    
    def calculate_backbone_plus_cb_rotation(self, target_residue_position, utilised_residues_indices):
        """
        Calculate rotation vector for a specific residue using N, CA, C, CB atoms in sliding window
        
        Args:
            target_residue_position: Position in utilised_residues_indices list
            utilised_residues_indices: List of residue indices
            
        Returns:
            rotation_magnitude: Magnitude of rotation vector
            rotation_vector: 3D rotation vector in degrees
            window_info: Dictionary with details about the calculation
        """
        half_window = self.window_size // 2
        
        # Calculate window bounds
        window_start = target_residue_position - half_window
        window_end = target_residue_position + half_window + 1
        
        # Get polymers
        polymer_1 = self.protein_1.get_polymer()
        polymer_2 = self.protein_2.get_polymer()  # This should be fitted_protein_2
        
        # Collect coordinates for all atoms in the window
        coords_1 = []
        coords_2 = []
        residue_info = []
        atoms_found = []
        
        for i in range(window_start, window_end):
            residue_idx_1 = utilised_residues_indices[i]
            residue_idx_2 = utilised_residues_indices[i]  # Same index for both proteins
            
            residue_1 = polymer_1[residue_idx_1]
            residue_2 = polymer_2[residue_idx_2]
            
            residue_atoms_1 = []
            residue_atoms_2 = []
            residue_atom_names = []
            
            # Check for each required atom (N, CA, C, CB)
            for atom_name in self.atoms_to_use:
                try:
                    atom_1 = residue_1[atom_name][0]
                    atom_2 = residue_2[atom_name][0]
                    
                    coords_1.append(atom_1.pos)
                    coords_2.append(atom_2.pos)
                    residue_atoms_1.append(atom_1.pos)
                    residue_atoms_2.append(atom_2.pos)
                    residue_atom_names.append(atom_name)
                    
                except KeyError:
                    if atom_name == "CB" and residue_1.name == "GLY":
                        # This should not happen since we filter out glycines
                        raise ValueError(f"Found glycine {residue_1.name} at position {i} (seq {residue_1.seqid.num}) - this should have been filtered out")
                    else:
                        raise ValueError(f"Residue {residue_1.name} at position {i} (seq {residue_1.seqid.num}) lacks {atom_name} atom")
            
            residue_info.append({
                'position': i,
                'residue_num': residue_1.seqid.num,
                'residue_name': residue_1.name,
                'is_target': i == target_residue_position,
                'atoms_found': residue_atom_names.copy()
            })
            
            atoms_found.extend(residue_atom_names)
        
        expected_atoms = len(self.atoms_to_use) * self.window_size
        if len(coords_1) != expected_atoms:
            raise ValueError(f"Expected {expected_atoms} atoms ({self.atoms_to_use} × {self.window_size}), found {len(coords_1)}")
        
        # Perform superposition
        sup_result = gemmi.superpose_positions(coords_1, coords_2)
        
        # Extract rotation matrix and convert to rotation vector
        rotation_matrix = np.array(sup_result.transform.mat.tolist())
        rotation_vector = Rotation.from_matrix(rotation_matrix).as_rotvec(degrees=True)
        rotation_magnitude = np.linalg.norm(rotation_vector)
        
        window_info = {
            'residue_info': residue_info,
            'target_residue': residue_info[half_window],  # Middle residue is the target
            'rmsd': sup_result.rmsd,
            'window_size': self.window_size,
            'atoms_per_residue': len(self.atoms_to_use),
            'total_atoms_used': len(coords_1),
            'atom_types': self.atoms_to_use
        }
        
        return rotation_magnitude, rotation_vector, window_info
    
    def find_and_calculate_suitable_residue(self, utilised_residues_indices):
        """
        Find a suitable residue (no nearby glycines) and calculate its rotation using N,CA,C,CB
        
        Returns:
            result: Dictionary with calculation results
        """
        # Find suitable residue
        position, residue_idx = self.find_residue_without_nearby_glycine(utilised_residues_indices)
        
        if position is None:
            return None
        
        # Calculate rotation for this residue
        magnitude, rot_vector, window_info = self.calculate_backbone_plus_cb_rotation(
            position, utilised_residues_indices
        )
        
        result = {
            'target_position': position,
            'target_residue_index': residue_idx,
            'target_residue_info': window_info['target_residue'],
            'rotation_magnitude': magnitude,
            'rotation_vector': rot_vector,
            'window_info': window_info
        }
        
        return result
    
    def compare_with_backbone_only(self, target_residue_position, utilised_residues_indices):
        """
        Compare rotation calculation using N,CA,C,CB vs N,CA,C only
        This helps verify the CB contribution
        """
        # Calculate with backbone + CB
        mag_full, vec_full, info_full = self.calculate_backbone_plus_cb_rotation(
            target_residue_position, utilised_residues_indices
        )
        
        # Calculate with backbone only (temporarily change atoms)
        original_atoms = self.atoms_to_use.copy()
        self.atoms_to_use = ["N", "CA", "C"]  # Backbone only
        
        try:
            mag_backbone, vec_backbone, info_backbone = self.calculate_backbone_plus_cb_rotation(
                target_residue_position, utilised_residues_indices
            )
        finally:
            self.atoms_to_use = original_atoms  # Restore original
        
        return {
            'backbone_plus_cb': {'magnitude': mag_full, 'vector': vec_full, 'info': info_full},
            'backbone_only': {'magnitude': mag_backbone, 'vector': vec_backbone, 'info': info_backbone},
            'cb_contribution': {
                'magnitude_diff': mag_full - mag_backbone,
                'vector_diff': vec_full - vec_backbone
            }
        }

# Integration with your Engine class
def add_dyndom_comparison_to_engine(engine):
    """
    Example of how to use this with your existing Engine for DynDom comparison
    """
    # Make sure you've run superimpose_chains first
    if engine.fitted_protein_2 is None:
        raise ValueError("Must run superimpose_chains() first")
    
    # Create calculator
    calc = SingleResidueBackbonePlusCBCalculator(engine.protein_1, engine.fitted_protein_2, 
                                                 window_size=engine.window)
    
    # Find and analyze a suitable residue
    result = calc.find_and_calculate_suitable_residue(engine.protein_1.utilised_residues_indices)
    
    if result:
        print("=== DynDom Comparison: N,CA,C,CB Rotation Analysis ===")
        print(f"Target residue: {result['target_residue_info']['residue_name']} {result['target_residue_info']['residue_num']}")
        print(f"Position in analysis: {result['target_position']}")
        print(f"Window size: {result['window_info']['window_size']}")
        print(f"Atoms per residue: {result['window_info']['atoms_per_residue']} ({result['window_info']['atom_types']})")
        print(f"Total atoms used: {result['window_info']['total_atoms_used']}")
        print()
        print("RESULTS FOR FORTRAN DYNDOM COMPARISON:")
        print(f"Rotation magnitude: {result['rotation_magnitude']:.6f} degrees")
        print(f"Rotation vector: [{result['rotation_vector'][0]:.6f}, {result['rotation_vector'][1]:.6f}, {result['rotation_vector'][2]:.6f}]")
        print(f"Window RMSD: {result['window_info']['rmsd']:.6f}Å")
        print()
        print("FOR QUICK COMPARISON:")
        print(f"Magnitude: {result['rotation_magnitude']:.6f}")
        print(f"Vector:    [{result['rotation_vector'][0]:.6f}, {result['rotation_vector'][1]:.6f}, {result['rotation_vector'][2]:.6f}]")
        print()
        print(f"Window residues:")
        for res_info in result['window_info']['residue_info']:
            marker = " <-- TARGET" if res_info['is_target'] else ""
            atoms_str = ",".join(res_info['atoms_found'])
            print(f"  {res_info['residue_name']} {res_info['residue_num']} ({atoms_str}){marker}")
        
        return result
    else:
        print("No suitable residue found (all have glycine in window)")
        return None

# Method to add to your Engine class
def get_dyndom_rotation_for_comparison(self, target_residue_position=None):
    """
    Method to add to Engine class for DynDom-style rotation calculation
    
    Args:
        target_residue_position: Specific position to analyze, or None to auto-find
    """
    calc = SingleResidueBackbonePlusCBCalculator(self.protein_1, self.fitted_protein_2, 
                                                 window_size=self.window)
    
    if target_residue_position is None:
        # Auto-find suitable residue
        return calc.find_and_calculate_suitable_residue(self.protein_1.utilised_residues_indices)
    else:
        # Calculate for specific residue
        magnitude, rot_vector, window_info = calc.calculate_backbone_plus_cb_rotation(
            target_residue_position, self.protein_1.utilised_residues_indices
        )
        return {
            'target_position': target_residue_position,
            'rotation_magnitude': magnitude,
            'rotation_vector': rot_vector,
            'window_info': window_info
        }