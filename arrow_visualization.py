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
                    
    def create_arrow_atoms(self, domain_data: dict, start_atom_id: int = 1000) -> List[str]:
        """Create PDB atom lines for an arrow visualization"""
        screw_axis = domain_data["screw_axis"]
        point_on_axis = domain_data["point_on_axis"]
        
        # Normalize the screw axis
        screw_axis = screw_axis / np.linalg.norm(screw_axis)
        
        # Length of arrow components currently hardcoded - will need to change
        arrow_length = 50.0  
        shaft_length = 35.0  
        head_length = 15.0   
        spacing = 3
        
        pdb_lines = []
        atom_id = start_atom_id
        res_id = start_atom_id // 100  # Use different residue IDs
        
        # Create shaft atoms (every 2.5Å for smoother line)
        n_shaft_atoms = int(shaft_length / spacing)
        for i in range(n_shaft_atoms):
            t = -shaft_length/2 + i * spacing
            pos = point_on_axis + t * screw_axis
            
            # Use different residue numbers for shaft vs head
            pdb_line = f"ATOM  {atom_id:5d}  CA  SHF A{res_id:4d}    {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}  1.00 30.00           C"
            pdb_lines.append(pdb_line)
            atom_id += 1
            
        # Create arrowhead atoms (larger and more visible)
        n_head_atoms = int(head_length / spacing)
        for i in range(n_head_atoms):
            t = shaft_length/2 + i * spacing
            pos = point_on_axis + t * screw_axis
            
            # Different residue ID for head
            pdb_line = f"ATOM  {atom_id:5d}  CA  ARH A{res_id+1:4d}    {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}  1.00 40.00           C"
            pdb_lines.append(pdb_line)
            atom_id += 1
            
        return pdb_lines
    
    def create_arrows_pdb(self) -> str:
        """Create PDB file with arrow visualizations"""
        if not self.domain_data:
            self.parse_w5_info()
            
        pdb_filename = f"{self.output_base}_arrows.pdb"
        
        pdb_lines = [
            "REMARK DynDom Arrow Visualization",
            "REMARK Shaft atoms = SHF residues (fixed domain)",
            "REMARK Head atoms = ARH residues (moving domain)",
            "REMARK Generated from: " + self.w5_info_file,
        ]
        
        atom_id = 1000
        for i, domain in enumerate(self.domain_data):
            pdb_lines.append(f"REMARK Arrow {i+1}: Domain {domain['domain_id']} moving relative to fixed domain {self.fixed_domain_id}")
            pdb_lines.append(f"REMARK   Rotation angle: {domain['rotation_angle']:.1f} degrees")
            
            arrow_atoms = self.create_arrow_atoms(domain, atom_id)
            pdb_lines.extend(arrow_atoms)
            pdb_lines.append("TER")
            atom_id += 100  # Space out atom IDs
            
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
            fixed_color = domain_colors[self.fixed_domain_id % len(domain_colors)]
            moving_color = domain_colors[domain['domain_id'] % len(domain_colors)]
            
            res_id_base = (1000 + i*100) // 100  # Calculate residue ID base
            
            script_lines.extend([
                f"# Arrow {i+1}: Domain {domain['domain_id']} -> Domain {self.fixed_domain_id}",
                f"# Rotation: {domain.get('rotation_angle', 'unknown'):.1f}°, Translation: {domain.get('translation', 'unknown'):.1f}Å",
                f"",
                f"# Select shaft and head atoms",
                f"select shaft_{i+1}, resn SHF and resi {res_id_base}",
                f"select head_{i+1}, resn ARH and resi {res_id_base+1}",
                f"",
                f"# Display shaft (fixed domain color: {fixed_color})",
                f"show spheres, shaft_{i+1}",
                f"show sticks, shaft_{i+1}",
                f"color {fixed_color}, shaft_{i+1}",
                f"set sphere_scale, 0.8, shaft_{i+1}",  # Larger spheres
                f"set stick_radius, 0.3, shaft_{i+1}",  # Thicker sticks
                f"",
                f"# Display head (moving domain color: {moving_color})",
                f"show spheres, head_{i+1}",
                f"show sticks, head_{i+1}",
                f"color {moving_color}, head_{i+1}",
                f"set sphere_scale, 1.2, head_{i+1}",   # Even larger for head
                f"set stick_radius, 0.5, head_{i+1}",   # Much thicker for head
                f"",
                f"# Connect atoms within each section",
                f"bond shaft_{i+1}, shaft_{i+1}",
                f"bond head_{i+1}, head_{i+1}",
                f"# Connect shaft to head",
                f"bond (resi {res_id_base} and name CA), (resi {res_id_base+1} and name CA)",
                f"",
            ])
            
        script_lines.extend([
            "# Make arrows more prominent",
            "set stick_transparency, 0.0",
            "set sphere_transparency, 0.0",
            "set stick_quality, 15",
            "set sphere_quality, 3",
            "",
            "# Final settings",
            "bg_color white",
            "set depth_cue, 0",
            "set ray_shadows, 0",
            "",
            "# Show domain assignments on protein",
            f"color {domain_colors[self.fixed_domain_id % len(domain_colors)]}, {self.output_base} and not (resn SHF or resn ARH)",
            "",
            "# Center view",
            "zoom all",
            "orient",
            "",
            "# Clean up selections", 
            "delete shaft_*",
            "delete head_*",
            "",
            "print 'DynDom arrows loaded successfully!'",
            f"print 'Fixed domain: {self.fixed_domain_id} ({domain_colors[self.fixed_domain_id % len(domain_colors)]})'",
        ])
        
        # Add info about each moving domain
        for i, domain in enumerate(self.domain_data):
            moving_color = domain_colors[domain['domain_id'] % len(domain_colors)]
            script_lines.append(f"print 'Moving domain {domain['domain_id']}: {moving_color} arrow, {domain.get('rotation_angle', 0):.1f}° rotation'")
        
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

