import traceback
import urllib.request
from itertools import groupby
from operator import itemgetter
from pathlib import Path

input_command_file_path = "data"


def read_command_file():
    temp_dict = {}
    try:
        fr = open(f"{input_command_file_path}/command.txt", "r")
        lines = fr.readlines()
        for line in lines:
            if not ("#" in line):
                line = line.replace("\n", "")
                line = line.replace(" ", "")
                tokens = line.split("=")
                param_name = tokens[0]
                param_val = tokens[1]
                if "chain" in param_name:
                    param_val = param_val.upper()
                elif "filename" in param_name:
                    param_val = param_val.lower()
                temp_dict[param_name] = param_val
        fr.close()
        check_pdb_exists(temp_dict['input_path'], temp_dict['filename1'])
        check_pdb_exists(temp_dict['input_path'], temp_dict['filename2'])
    except Exception as e:
        print(e)
    return temp_dict


def check_pdb_exists(input_path, file_name):
    """
    Uses HTTPS request to make sure that the pdb file is in the directory
    :param input_path: The directory to the file
    :param file_name:
    :return:
    """
    file_path = f"{input_path}/{file_name}.pdb"
    path = Path(file_path)
    if not path.exists():
        try:
            urllib.request.urlretrieve(
                f"https://files.rcsb.org/download/{file_name}.pdb",
                file_path
            )
        except Exception as e:
            traceback.print_exc()
            print(e)


def read_param_file():
    temp_dict = {}
    try:
        fr = open(f"{input_command_file_path}/param.txt", "r")
        lines = fr.readlines()
        for line in lines:
            if not ("#" in line):
                line = line.replace("\n", "")
                line = line.replace(" ", "")
                tokens = line.split("=")
                param_name = tokens[0]
                param_val = tokens[1]
                if param_name == "window":
                    param_val = int(tokens[1])
                elif param_name == "domain":
                    param_val = int(tokens[1])
                elif param_name == "ratio":
                    param_val = float(tokens[1])
                elif param_name == "k_means_n_init":
                    param_val = int(tokens[1])
                elif param_name == "k_means_max_iter":
                    param_val = int(tokens[1])
                temp_dict[param_name] = param_val
        fr.close()
    except Exception as e:
        print(e)
    return temp_dict


def write_rotation_vec_to_pdb(output_path, protein_1, protein_2_name, chain_2, rotation_vectors):
    """
    Writes the rotation vectors of each residue of the slide window into a pdb file
    :param output_path:
    :param protein_1: The Protein object of protein 1
    :param protein_2_name:
    :param chain_2:
    :param rotation_vectors:
    :return:
    """
    try:
        dir_path_str = f"{output_path}/{protein_1.name}_{protein_1.chain_param}_{protein_2_name}_{chain_2}"
        dir_path = Path(dir_path_str)
        if not dir_path.exists():
            dir_path.mkdir(parents=True)
        fw = open(f"{dir_path_str}/{protein_1.name}_rot_vecs.pdb", "w")
        slide_window_indices = protein_1.slide_window_residues_indices
        protein_polymer = protein_1.get_polymer()
        for i in range(rotation_vectors.shape[0]):
            index = protein_1.utilised_residues_indices[i+slide_window_indices[0]]
            residue_name = protein_polymer[index].name
            residue_num = protein_polymer[index].seqid.num
            x = str(round(rotation_vectors[i][0], 3)).rjust(8, " ")
            y = str(round(rotation_vectors[i][1], 3)).rjust(8, " ")
            z = str(round(rotation_vectors[i][2], 3)).rjust(8, " ")
            row = f"ATOM         CA  {residue_name} A {residue_num}    {x}{y}{z}\n"
            fw.write(row)
        fw.close()
    except Exception as e:
        print(e)
        return False
    return True


def write_final_output_pdb(output_path, protein_1, fitted_protein_2, fitted_protein_2_name, fitted_protein_2_chain, fitted_protein_2_res_ind):
    """
    :param output_path:
    :param protein_1: The Protein object of protein 1
    :param fitted_protein_2: The Chain object of protein 2 fitted to protein 1
    :param fitted_protein_2_name: The name of the protein_2
    :param fitted_protein_2_chain: The chain of protein 2 fitted to protein 1
    :param fitted_protein_2_res_ind: The indices of residues used in the clustering
    :return:
    """
    try:
        folder_name = f"{protein_1.name}_{protein_1.chain_param}_{fitted_protein_2_name}_{fitted_protein_2_chain}"
        dir_path_str = f"{output_path}/{folder_name}"
        dir_path = Path(dir_path_str)
        if not dir_path.exists():
            dir_path.mkdir(parents=True)
        fw = open(f"{dir_path_str}/{folder_name}.pdb", "w")
        fw.write(f"MODEL{'1'.rjust(9, ' ')}\n")
        protein_1_res_ind = protein_1.utilised_residues_indices
        protein_1_residues = protein_1.get_polymer()
        atom_count = 1
        subchain = protein_1.chain_param
        for i in protein_1_res_ind:
            r = protein_1_residues[i]
            res_name = r.name.rjust(3, " ")
            res_num = str(r.seqid.num).rjust(4, " ")
            for a in r:
                atom_num = str(atom_count).rjust(5, " ")
                atom_name = a.name.ljust(4, " ")
                x = str(round(a.pos.x, 3)).rjust(8, " ")
                y = str(round(a.pos.y, 3)).rjust(8, " ")
                z = str(round(a.pos.z, 3)).rjust(8, " ")
                row = f"ATOM  {atom_num} {atom_name} {res_name} {subchain}{res_num}    {x}{y}{z}\n"
                fw.write(row)
                atom_count += 1
        fw.write("ENDMDL\n")

        fw.write(f"MODEL{'2'.rjust(9, ' ')}\n")
        atom_count = 1
        protein_2_residues = fitted_protein_2.get_polymer()
        for i in fitted_protein_2_res_ind:
            r = protein_2_residues[i]
            res_name = r.name.rjust(3, " ")
            res_num = str(r.seqid.num).rjust(4, " ")
            for a in r:
                atom_num = str(atom_count).rjust(5, " ")
                atom_name = a.name.ljust(4, " ")
                x = str(round(a.pos.x, 3)).rjust(8, " ")
                y = str(round(a.pos.y, 3)).rjust(8, " ")
                z = str(round(a.pos.z, 3)).rjust(8, " ")
                row = f"ATOM  {atom_num} {atom_name} {res_name} {fitted_protein_2_chain}{res_num}    {x}{y}{z}\n"
                fw.write(row)
                atom_count += 1
        fw.write("ENDMDL\n")
    except Exception as e:
        traceback.print_exc()
        print(e)


def write_final_output_pml(output_path, protein_1, protein_2_name, protein_2_chain, domains, fixed_domain_id, bending_residues, window_size):
    bend_res_colour = "[0  ,255,0  ]"
    dom_colours = ["[0  ,0  ,255]", "[255,0  ,0  ]", "[255,255,0  ]", "[255,100,255]", "[0  ,255,255]"]

    try:
        folder_name = f"{protein_1.name}_{protein_1.chain_param}_{protein_2_name}_{protein_2_chain}"
        dir_path_str = f"{output_path}/{folder_name}"
        dir_path = Path(dir_path_str)
        if not dir_path.exists():
            dir_path.mkdir(parents=True)
        fw = open(f"{dir_path}/{folder_name}.pml", "w")
        fw.write("reinitialize\n")
        fw.write(f"load {folder_name}.pdb\n")
        fw.write(f"bg_color white\n")
        fw.write("color grey\n")

        mid_point = (window_size - 1) // 2
        fixed_dom_segments = domains[fixed_domain_id].segments
        util_res = protein_1.utilised_residues_indices
        polymer = protein_1.get_polymer()

        fixed_dom_res_reg = []
        all_bend_res_indices = []
        for b in bending_residues.values():
            for i in b:
                bb = i+mid_point
                index = polymer[bb].seqid.num
                all_bend_res_indices.append(index)

        # Colour the fixed domains blue
        for s in range(fixed_dom_segments.shape[0]):
            reg = []
            for i in range(fixed_dom_segments[s][0], fixed_dom_segments[s][1]+1):
                j = i + mid_point
                index = util_res[j]
                res_num = polymer[index].seqid.num
                if res_num not in all_bend_res_indices:
                    reg.append(res_num)
            fixed_dom_res_reg.extend(group_continuous_regions(reg))

        for s in range(len(fixed_dom_res_reg)):
            print(s)
            if s == 0:
                sel_reg_str = f"select region0, resi {fixed_dom_res_reg[s][0]}-{fixed_dom_res_reg[s][1]}\n"
            else:
                print(fixed_dom_res_reg[s])
                sel_reg_str = f"select region0, region0 + resi {fixed_dom_res_reg[s][0]}-{fixed_dom_res_reg[s][1]}\n"
            fw.write(sel_reg_str)
        fw.write(f"set_color colour0 = {dom_colours[0]}\n")
        fw.write("color colour0, region0\n")

        # Colour the dynamic domains
        region_count = 1
        for domain in domains:
            dyn_dom_res_reg = []
            if domain.domain_id == fixed_domain_id:
                continue
            segments = domain.segments
            dom_bend_res = bending_residues[domain.domain_id]
            for s in range(segments.shape[0]):
                reg = []
                for i in range(segments[s][0], segments[s][1]+1):
                    j = i + mid_point
                    index = util_res[j]
                    res_num = polymer[index].seqid.num
                    if res_num not in dom_bend_res:
                        reg.append(res_num)

                dyn_dom_res_reg.extend(group_continuous_regions(reg))

            for s in range(len(dyn_dom_res_reg)):
                if s == 0:
                    sel_reg_str = f"select region{region_count}, resi {dyn_dom_res_reg[s][0]}-{dyn_dom_res_reg[s][1]}\n"
                else:
                    sel_reg_str = f"select region{region_count}, region{region_count} + resi {dyn_dom_res_reg[s][0]}-{dyn_dom_res_reg[s][1]}\n"
                fw.write(sel_reg_str)
            fw.write(f"set_color colour{region_count} = {dom_colours[region_count]}\n")
            fw.write(f"color colour{region_count}, region{region_count}\n")

            region_count += 1

        # Colour the bending residues
        bend_res_groups = group_continuous_regions(all_bend_res_indices)
        for g in bend_res_groups:
            fw.write(f"select region{region_count}, resi {g[0]}-{g[1]}\n")
            fw.write(f"set_color colour{region_count} = {bend_res_colour}\n")
            fw.write(f"color colour{region_count}, region{region_count}\n")
            region_count += 1

        fw.write("set dash_gap, 0\n")
        fw.write("set dash_radius, 0.2\n")

    except Exception as e:
        traceback.print_exc()
        print(e)


def write_w5_info_file(output_path, protein_1_name: str, chain_1, protein_2_name: str, chain_2, window, domain_size, ratio, atoms,
                       domains: list, fixed_domain_id: int):
    try:
        protein_folder = f"{protein_1_name}_{chain_1}_{protein_2_name}_{chain_2}"
        fw = open(f"{output_path}/{protein_folder}/{protein_folder}.w5_info", "w")
        fw.write("DynDom Python Version 1.0\n")
        fw.write(f"{protein_1_name}{chain_1}_{protein_2_name}{chain_2}.w5\n")
        fw.write(f"file name of conformer 1: {protein_1_name}.pdb\n")
        fw.write(f"chain id: {chain_1}\n")
        fw.write(f"file name of conformer 2: {protein_2_name}.pdb\n")
        fw.write(f"chain id: {chain_2}\n")
        fw.write(f"window length: {window}\n")
        fw.write(f"minimum ratio of external to internal motion: {ratio}\n")
        fw.write(f"minimum domain size: {domain_size}\n")
        fw.write(f"atoms to use: {atoms}\n")
        fw.write(f"THERE ARE {len(domains)} DOMAINS\n")
        fw.write("================================================================================\n")
        domain_colours = ["blue", "red", "yellow", "pink", "cyan"]
        fixed_domain = domains[fixed_domain_id]
        fw.write("FIXED DOMAIN\n")
        fw.write(f"DOMAIN NUMBER: \t {fixed_domain_id} (coloured {domain_colours[0]} for rasmol)\n")
        residue_str = ""
        for s in range(fixed_domain.segments.shape[0]):
            if s == 0:
                residue_str = f"{str(fixed_domain.segments[s][0])} - {str(fixed_domain.segments[s][1])}"
            else:
                residue_str = residue_str + f", {str(fixed_domain.segments[s][0])} - {str(fixed_domain.segments[s][1])}"
        fw.write(f"RESIDUE NUMBERS: \t{residue_str}\n")
        fw.write(f"SIZE: \t{fixed_domain.num_residues}\n")
        fw.write(f"BACKBONE RMSD ON THIS DOMAIN: \t{round(fixed_domain.rmsd, 3)}A\n")

        domain_count = 1
        for domain in domains:
            if domain.domain_id != fixed_domain_id:
                fw.write("------------------------------------------------------------------------------\n")
                fw.write(f"MOVING DOMAIN (RELATIVE TO FIXED DOMAIN),  PAIR {domain_count}\n")
                fw.write(f"DOMAIN NUMBER: \t {domain.domain_id} (coloured {domain_colours[domain_count]} for rasmol)\n")
                residue_str = ""
                for s in range(domain.segments.shape[0]):
                    if s == 0:
                        residue_str = f"{str(domain.segments[s][0])} - {str(domain.segments[s][1])}"
                    else:
                        residue_str = residue_str + f", {str(domain.segments[s][0])} - {str(domain.segments[s][1])}"
                fw.write(f"RESIDUE NUMBERS: \t{residue_str}\n")
                fw.write(f"SIZE: \t{domain.num_residues}\n")
                fw.write(f"BACKBONE RMSD ON THIS DOMAIN: \t{round(domain.rmsd, 3)}A\n")
                fw.write(f"RATIO OF INTERDOMAIN TO INTRADOMAIN DISPLACEMENT: \t{round(domain.ratio, 3)}\n")
                fw.write(f"ANGLE OF ROTATION: \t{round(domain.rot_angle, 3)} DEGREES\n")
                fw.write(f"TRANSLATION ALONG AXIS:\t{round(domain.translation, 3)} A\n")
                fw.write(f"SCREW AXIS: \t{round(domain.screw_axis[0], 3)} \t{round(domain.screw_axis[1], 3)} \t{round(domain.screw_axis[2], 3)}\n")
                fw.write(f"POINT ON AXIS: \t{round(domain.point_on_axis[0], 3)} \t{round(domain.point_on_axis[1], 3)} \t{round(domain.point_on_axis[2], 3)}\n")
                groups = group_continuous_regions(domain.bend_res)
                for group in groups:
                    fw.write(f"BENDING RESIDUES: \t{group[0]} - {group[-1]}\n")
                domain_count += 1

    except Exception as e:
        print(e)
        return False
    return True


def group_continuous_regions(data: list):
    groups = []
    for k, g in groupby(enumerate(data), lambda ix: ix[0] - ix[1]):
        temp = list(map(itemgetter(1), g))
        groups.append([temp[0], temp[-1]+1])
    return groups


