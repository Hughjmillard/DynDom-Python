import numpy as np
from typing import List, Tuple
import os

class ArrowGenerator:
    """Simple arrow generator for DynDom screw axes"""
    
    def __init__(self, w5_info_file: str, output_base: str):
        """
        Initialize with w5 info file and output base name
        
        Args:
            w5_info_file: Path to the .w5_info file
            output_base: Base name for output files (e.g., "4ake_A_2eck_B")
        """
        self.w5_info_file = w5_info_file
        self.output_base = output_base
        self.domain_data = []
        self.fixed_domain_id = None
        self.protein_coordinates = []
        
    def parse_w5_info(self):
        """Parse the w5_info file to extract domain and screw axis data"""
        with open(self.w5_info_file, 'r') as f:
            lines = f.readlines()
        
        current_domain = None
        in_moving_domain = False
        
        for i, line in enumerate(lines):
            line = line.strip()
            
            # Debug print to see what we're parsing
            if any(keyword in line for keyword in ["FIXED DOMAIN", "MOVING DOMAIN", "DOMAIN NUMBER:", "SCREW AXIS:", "POINT ON AXIS:"]):
                print(f"Parsing line: {line}")
            
            # Check for section headers - ORDER MATTERS!
            if "MOVING DOMAIN" in line:
                # Save previous domain if it exists and is complete
                if (current_domain is not None and 
                    "screw_axis" in current_domain and 
                    "point_on_axis" in current_domain and 
                    "domain_id" in current_domain):
                    print(f"Completed domain {current_domain['domain_id']} before new moving domain")
                    self.domain_data.append(current_domain)
                
                in_moving_domain = True
                current_domain = {"is_fixed": False}
                print("Entered MOVING DOMAIN section")
                continue
                
            elif "FIXED DOMAIN" in line and "MOVING DOMAIN" not in line:
                in_moving_domain = False
                print("Entered FIXED DOMAIN section")
                continue
                
            # Parse domain numbers
            elif "DOMAIN NUMBER:" in line:
                domain_num = int(line.split()[2])
                
                if not in_moving_domain:
                    # This is the fixed domain
                    self.fixed_domain_id = domain_num
                    print(f"Found fixed domain: {domain_num}")
                else:
                    # This is a moving domain
                    if current_domain is None:
                        current_domain = {"is_fixed": False}
                    current_domain["domain_id"] = domain_num
                    print(f"Found moving domain: {domain_num}")
                
            # Parse screw axis data (only for moving domains)
            elif "SCREW AXIS:" in line and in_moving_domain and current_domain is not None:
                parts = line.split()
                screw_axis = [float(parts[2]), float(parts[3]), float(parts[4])]
                current_domain["screw_axis"] = np.array(screw_axis)
                print(f"Found screw axis: {screw_axis}")
                
            elif "POINT ON AXIS:" in line and in_moving_domain and current_domain is not None:
                parts = line.split()
                point_on_axis = [float(parts[3]), float(parts[4]), float(parts[5])]
                current_domain["point_on_axis"] = np.array(point_on_axis)
                print(f"Found point on axis: {point_on_axis}")
                
            elif "ANGLE OF ROTATION:" in line and in_moving_domain and current_domain is not None:
                angle = float(line.split()[3])
                current_domain["rotation_angle"] = angle
                print(f"Found rotation angle: {angle}")
                
            elif "TRANSLATION ALONG AXIS:" in line and in_moving_domain and current_domain is not None:
                translation = float(line.split()[3])
                current_domain["translation"] = translation
                print(f"Found translation: {translation}")
                
        # Don't forget the last domain if file ends
        if (current_domain is not None and 
            "screw_axis" in current_domain and 
            "point_on_axis" in current_domain and 
            "domain_id" in current_domain):
            print(f"Completed final domain {current_domain['domain_id']}")
            self.domain_data.append(current_domain)
            
        print(f"Total domains found: {len(self.domain_data)}")
        print(f"Fixed domain ID: {self.fixed_domain_id}")
        for i, domain in enumerate(self.domain_data):
            print(f"Domain {i+1}: ID={domain['domain_id']}, has_screw_axis={('screw_axis' in domain)}")
            
        # DEBUG: Print the actual color mappings
        domain_colors = ["blue", "red", "yellow", "green", "orange", "purple", "pink"]
        print(f"DEBUG: Fixed domain {self.fixed_domain_id} -> color index {(self.fixed_domain_id - 1) % len(domain_colors)} -> {domain_colors[(self.fixed_domain_id - 1) % len(domain_colors)]}")
        for domain in self.domain_data:
            print(f"DEBUG: Moving domain {domain['domain_id']} -> color index {(domain['domain_id'] - 1) % len(domain_colors)} -> {domain_colors[(domain['domain_id'] - 1) % len(domain_colors)]}")
            

    def read_protein_coordinates(self):
        """Read protein coordinates from PDB file to calculate proper arrow scaling"""
        pdb_file = f"{self.output_base}.pdb"
        if not os.path.exists(pdb_file):
            print(f"Warning: PDB file {pdb_file} not found, using default arrow scaling")
            return
            
        coordinates = []
        try:
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith('ATOM') and line[12:16].strip() in ['CA', 'N', 'C']:
                        x = float(line[30:38])
                        y = float(line[38:46]) 
                        z = float(line[46:54])
                        coordinates.append([x, y, z])
            
            if coordinates:
                self.protein_coordinates = np.array(coordinates)
                print(f"Read {len(coordinates)} protein atoms for scaling")
            else:
                print("No suitable atoms found in PDB file")
                
        except Exception as e:
            print(f"Error reading PDB file: {e}")

    def calculate_arrow_length(self, screw_axis, point_on_axis):
        """Calculate appropriate arrow length like axleng.f does"""
        if len(self.protein_coordinates) == 0:
            print("Using default arrow length (no protein coordinates)")
            return 40.0  # Larger default shaft length
            
        # Normalize screw axis
        unit_axis = screw_axis / np.linalg.norm(screw_axis)
        
        # Calculate projections of all protein atoms onto the screw axis
        # This is equivalent to the dot product calculation in axleng.f
        relative_coords = self.protein_coordinates - point_on_axis
        projections = np.dot(relative_coords, unit_axis)
        
        # Find min and max projections
        tmin = np.min(projections)
        tmax = np.max(projections)
        
        # Add 10% padding on each end (like axleng.f does)
        padding = (tmax - tmin) / 10.0
        tmin_padded = tmin - padding
        tmax_padded = tmax + padding
        
        # Total length along the axis - this should be the FULL span
        total_length = tmax_padded - tmin_padded
        
        # Don't cap it too aggressively - let it be the protein size + 20%
        total_length = max(total_length, 30.0)  # Minimum 30Å
        # Remove the maximum cap to let it scale properly with large proteins
        
        print(f"Protein projection range: {tmin:.1f} to {tmax:.1f}")
        print(f"With padding: {tmin_padded:.1f} to {tmax_padded:.1f}")
        print(f"Calculated arrow shaft length: {total_length:.1f}Å")
        
        return total_length

    def create_arrow_atoms(self, domain_data: dict, arrow_index: int, start_atom_id: int = 1000) -> List[str]:
        """Create PDB atom lines for an arrow visualization with DynDom-style arrow head"""
        screw_axis = domain_data["screw_axis"]
        point_on_axis = domain_data["point_on_axis"]
        
        # Normalize the screw axis
        screw_axis = screw_axis / np.linalg.norm(screw_axis)
        
        # Calculate appropriate arrow length based on protein size
        base_shaft_length = self.calculate_arrow_length(screw_axis, point_on_axis)
        head_length = 8.0  # Slightly larger head length
        
        # Extend shaft backwards by 10% to balance the arrow head extension
        shaft_extension = base_shaft_length * 0.1
        total_shaft_length = base_shaft_length + shaft_extension
        
        shaft_spacing = 1.0  # Even denser spacing for smoother shaft
        
        pdb_lines = []
        atom_id = start_atom_id
        
        # Use widely separated residue IDs for each arrow to prevent cross-connections
        shaft_res_id = 100 + arrow_index * 50
        head_res_id = shaft_res_id + 20
        
        # Use different chain IDs for each arrow
        chain_id = chr(ord('A') + arrow_index)
        
        # === CREATE SHAFT (linear, extended backwards) ===
        n_shaft_atoms = int(total_shaft_length / shaft_spacing)
        for i in range(n_shaft_atoms):
            # Start further back: -base_length/2 - extension, end at +base_length/2
            t = -(base_shaft_length/2 + shaft_extension) + i * shaft_spacing
            pos = point_on_axis + t * screw_axis
            
            pdb_line = f"ATOM  {atom_id:5d}  CA  SHF {chain_id}{shaft_res_id:4d}    {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}  1.00 30.00           C"
            pdb_lines.append(pdb_line)
            atom_id += 1
            
        pdb_lines.append("TER")
        
        # === CREATE ARROW HEAD (overlapping with shaft end) ===
        # Create two perpendicular vectors to the screw axis for the cone base
        if abs(screw_axis[0]) < 0.9:
            temp_vec = np.array([1.0, 0.0, 0.0])
        else:
            temp_vec = np.array([0.0, 1.0, 0.0])
        
        # Create perpendicular vectors using cross product
        perp1 = np.cross(screw_axis, temp_vec)
        perp1 = perp1 / np.linalg.norm(perp1)
        perp2 = np.cross(screw_axis, perp1)
        perp2 = perp2 / np.linalg.norm(perp2)
        
        # Arrow head starts INSIDE the shaft, overlapping the final portion
        head_overlap = 3.0  # How far back into the shaft the head starts
        head_start_pos = point_on_axis + (base_shaft_length/2 - head_overlap) * screw_axis
        
        # Create smaller, more compact cone structure
        n_layers = 5  # Fewer layers for smaller head
        n_points_per_layer = 6  # Fewer points per layer
        
        for layer in range(n_layers + 1):
            # Position along the cone axis (0 = base, 1 = tip)
            layer_progress = layer / n_layers
            layer_pos = head_start_pos + layer_progress * head_length * screw_axis
            
            # Smaller radius that decreases more aggressively
            max_radius = 1.8  # Much smaller maximum radius at base
            layer_radius = max_radius * (1.0 - layer_progress) ** 1.5  # Faster taper
            
            if layer == n_layers:
                # Tip of arrow - single point
                pdb_line = f"ATOM  {atom_id:5d}  CA  ARH {chain_id}{head_res_id:4d}    {layer_pos[0]:8.3f}{layer_pos[1]:8.3f}{layer_pos[2]:8.3f}  1.00 50.00           C"
                pdb_lines.append(pdb_line)
                atom_id += 1
            else:
                # Create circular layer
                for point in range(n_points_per_layer):
                    angle = 2 * np.pi * point / n_points_per_layer
                    
                    # Calculate position on the circle
                    circle_pos = (layer_pos + 
                                 layer_radius * np.cos(angle) * perp1 + 
                                 layer_radius * np.sin(angle) * perp2)
                    
                    pdb_line = f"ATOM  {atom_id:5d}  CA  ARH {chain_id}{head_res_id:4d}    {circle_pos[0]:8.3f}{circle_pos[1]:8.3f}{circle_pos[2]:8.3f}  1.00 40.00           C"
                    pdb_lines.append(pdb_line)
                    atom_id += 1
        
        pdb_lines.append("TER")
            
        return pdb_lines
    
    def create_arrows_pdb(self) -> str:
        """Create PDB file with arrow visualizations"""
        if not self.domain_data:
            self.parse_w5_info()
            
        # Read protein coordinates for proper scaling
        self.read_protein_coordinates()
            
        pdb_filename = f"{self.output_base}_arrows.pdb"
        
        pdb_lines = [
            "REMARK DynDom Arrow Visualization",
            "REMARK Shaft atoms = SHF residues",
            "REMARK Head atoms = ARH residues", 
            "REMARK Each arrow uses separate chain ID (A, B, C...)",
            "REMARK Arrow length calculated from protein dimensions",
            "REMARK Generated from: " + self.w5_info_file,
        ]
        
        atom_id = 1000
        for i, domain in enumerate(self.domain_data):
            pdb_lines.append(f"REMARK Arrow {i+1}: Domain {domain['domain_id']} moving relative to fixed domain {self.fixed_domain_id}")
            pdb_lines.append(f"REMARK   Chain {chr(ord('A') + i)}: Rotation angle: {domain['rotation_angle']:.1f} degrees")
            
            arrow_atoms = self.create_arrow_atoms(domain, i, atom_id)
            pdb_lines.extend(arrow_atoms)
            atom_id += 200  # Large gap between arrows to prevent any connections
            
        pdb_lines.append("END")
            
        # Write to file
        with open(pdb_filename, 'w') as f:
            f.write('\n'.join(pdb_lines))
            
        return pdb_filename
    
    def create_pymol_script(self, pdb_filename: str) -> str:
        """Create PyMOL script for arrow visualization"""
        if not self.domain_data:
            self.parse_w5_info()
            
        pymol_filename = f"{self.output_base}_arrows.pml"
        
        # Domain colors (matching your current system)
        domain_colors = ["blue", "red", "yellow", "green", "orange", "purple", "pink"]
        
        script_lines = [
            "# DynDom Arrow Visualization Script",
            f"# Load main protein structure",
            f"load {self.output_base}.pdb",
            f"",
            f"# Load arrows",
            f"load {os.path.basename(pdb_filename)}",
            f"",
            f"# Basic protein display",
            f"hide everything, {self.output_base}",
            f"show cartoon, {self.output_base}",
            f"color gray80, {self.output_base}",
            f"",
            f"# Hide arrow atoms initially",
            f"hide everything, {os.path.splitext(os.path.basename(pdb_filename))[0]}",
            f"",
        ]
        
        # Add arrow visualization for each domain
        for i, domain in enumerate(self.domain_data):
            # CORRECTED COLOR MAPPING FOR PYTHON DYNDOM:
            # From your .w5_info: Fixed domain = 1, Moving domain = 0
            # We want: Shaft = blue (fixed), Head = red (moving)
            
            domain_id = domain['domain_id']
            fixed_id = self.fixed_domain_id
            
            # Direct mapping without special cases:
            # Domain 0 -> index 0 -> blue
            # Domain 1 -> index 1 -> red  
            # Domain 2 -> index 2 -> yellow, etc.
            moving_color = domain_colors[domain_id % len(domain_colors)]
            fixed_color = domain_colors[fixed_id % len(domain_colors)]
            
            # Use chain-specific selections to completely separate arrows
            chain_id = chr(ord('A') + i)
            shaft_res_id = 100 + i * 50
            head_res_id = shaft_res_id + 20
            
            script_lines.extend([
                f"# Arrow {i+1}: Domain {domain['domain_id']} (moving) -> Domain {self.fixed_domain_id} (fixed)",
                f"# Shaft color: {moving_color} (moving domain), Head color: {fixed_color} (fixed domain)",
                f"# Rotation: {domain.get('rotation_angle', 'unknown'):.1f}°, Translation: {domain.get('translation', 'unknown'):.1f}Å",
                f"",
                f"# Select shaft and head atoms by chain and residue",
                f"select shaft_{i+1}, chain {chain_id} and resn SHF and resi {shaft_res_id}",
                f"select head_{i+1}, chain {chain_id} and resn ARH and resi {head_res_id}",
                f"",
                f"# Display shaft as thick licorice stick (MOVING domain color: {moving_color})",
                f"show sticks, shaft_{i+1}",
                f"color {moving_color}, shaft_{i+1}",
                f"set stick_radius, 0.3, shaft_{i+1}",  # Thicker for licorice style
                f"",
                f"# Display arrow head as clean cone (FIXED domain color: {fixed_color})",
                f"show sticks, head_{i+1}",
                f"color {fixed_color}, head_{i+1}",
                f"set stick_radius, 0.25, head_{i+1}",
                f"",
                f"# Connect atoms ONLY within each section",
                f"bond shaft_{i+1}, shaft_{i+1}",
                f"bond head_{i+1}, head_{i+1}",
                f"",
            ])
            
        script_lines.extend([
            "# Disable automatic bonding between different chains",
            "set auto_bond, 0",
            "",
            "# Make arrows more prominent",
            "set stick_transparency, 0.0",
            "set stick_quality, 15",
            "set sphere_quality, 3",
            "set surface_quality, 2",
            "",
            "# Final settings",
            "bg_color white",
            "set depth_cue, 0",
            "set ray_shadows, 1",
            "set ray_shadow_decay_factor, 0.1",
            "",
            "# Better lighting for 3D arrow heads",
            "set ambient, 0.2",
            "set direct, 0.8",
            "set reflect, 0.5",
            "set shininess, 10",
            "",
            "# Center view",
            "zoom all",
            "orient",
            "",
            "# Clean up selections", 
            "delete shaft_*",
            "delete head_*",
            "",
            "print 'DynDom arrows with 3D heads loaded successfully!'",
            f"print 'Fixed domain: {self.fixed_domain_id} ({domain_colors[self.fixed_domain_id % len(domain_colors)]})'",
        ])
        
        # Add info about each moving domain
        for i, domain in enumerate(self.domain_data):
            domain_id = domain['domain_id']
            fixed_id = self.fixed_domain_id
            
            # Apply same direct mapping as above
            moving_color = domain_colors[domain_id % len(domain_colors)]
            fixed_color = domain_colors[fixed_id % len(domain_colors)]
                
            chain_id = chr(ord('A') + i)
            script_lines.append(f"print 'Moving domain {domain['domain_id']}: Chain {chain_id}, {fixed_color} shaft with {moving_color} head, {domain.get('rotation_angle', 0):.1f}° rotation'")
        
        # Write script
        with open(pymol_filename, 'w') as f:
            f.write('\n'.join(script_lines))
            
        return pymol_filename
    
    def generate_all_arrows(self) -> dict:
        """Generate complete arrow visualization set"""
        print(f"Generating arrows from: {self.w5_info_file}")
        
        # Parse the w5 info
        self.parse_w5_info()
        
        if not self.domain_data:
            print("No moving domains found - no arrows to generate")
            return {}
            
        print(f"Found {len(self.domain_data)} moving domains")
        
        # Create PDB with arrows
        pdb_file = self.create_arrows_pdb()
        print(f"Created arrow PDB: {pdb_file}")
        
        # Create PyMOL script
        pymol_file = self.create_pymol_script(pdb_file)
        print(f"Created PyMOL script: {pymol_file}")
        
        return {
            'pdb': pdb_file,
            'pymol': pymol_file,
            'domains': len(self.domain_data)
        }

