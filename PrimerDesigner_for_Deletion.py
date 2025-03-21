import pandas as pd
import argparse
import os
import subprocess
import time
import logging
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
import snapgene_reader
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import nt_search
from Bio.SeqUtils import GC
import sys

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
excel_log = []
MAX_LINE_WIDTH = 80
QUALIFIER_START = 21
FEATURE_WIDTH = 15

def read_genbank_table(file_path):
    logging.info(f"Reading GenBank table from {file_path}")
    start_time = time.time()
    if file_path.endswith('.xlsx'):
        df = pd.read_excel(file_path)
    else:
        df = pd.read_csv(file_path, sep='\t')
    logging.info(f"GenBank table loaded in {time.time() - start_time:.2f} seconds")
    return df

def read_genbank_file(file_path):
    logging.info(f"Reading GenBank file from {file_path}")
    start_time = time.time()
    records = []
    with open(file_path, 'r') as f:
        for record in SeqIO.parse(f, "genbank"):
            records.append(record)
    if not records:
        raise ValueError(f"No records found in GenBank file: {file_path}")
    logging.info(f"Found {len(records)} records in GenBank file, loaded in {time.time() - start_time:.2f} seconds")
    return records

def get_deletion_arms(genbank_records, genbank_df, locus_tag, upstream_bp, downstream_bp):
    logging.info(f"Fetching deletion arms for locus_tag: {locus_tag} with upstream {upstream_bp} bp and downstream {downstream_bp} bp")
    start_time = time.time()
    
    target_record = genbank_records[0]
    seq_length = len(target_record.seq)
    
    if '-' in locus_tag:
        start_tag, end_tag = locus_tag.split('-', 1)
        start_row = genbank_df[genbank_df['locus_tag'] == start_tag]
        end_row = genbank_df[genbank_df['locus_tag'] == end_tag]
        
        if start_row.empty or end_row.empty:
            raise ValueError(f"One or both locus tags ({start_tag}, {end_tag}) not found in the Excel table.")
        
        target_start = int(start_row['start'].iloc[0])
        target_end = int(end_row['end'].iloc[0])
        
        if target_start >= target_end:
            raise ValueError(f"Start locus_tag {start_tag} must be before end locus_tag {end_tag} in the genome.")
    else:
        target_row = genbank_df[genbank_df['locus_tag'] == locus_tag]
        if target_row.empty:
            raise ValueError(f"Locus tag {locus_tag} not found in the Excel table.")
        
        target_start = int(target_row['start'].iloc[0])
        target_end = int(target_row['end'].iloc[0])
    
    upstream_start = max(1, target_start - upstream_bp)
    upstream_end = target_start - 1
    downstream_start = target_end + 1
    downstream_end = min(seq_length, target_end + downstream_bp)
    
    upstream_seq = str(target_record.seq[upstream_start - 1:upstream_end])
    downstream_seq = str(target_record.seq[downstream_start - 1:downstream_end])
    
    upstream_length = upstream_end - (upstream_start - 1)
    downstream_length = downstream_end - (downstream_start - 1)
    if upstream_length != upstream_bp:
        logging.warning(f"Upstream length {upstream_length} does not match requested {upstream_bp} bp")
    if downstream_length != downstream_bp:
        logging.warning(f"Downstream length {downstream_length} does not match requested {downstream_bp} bp")
    
    combined_seq = upstream_seq + downstream_seq
    cds_positions = {
        f"{locus_tag}_upstream": (0, len(upstream_seq)),
        f"{locus_tag}_downstream": (len(upstream_seq), len(upstream_seq) + len(downstream_seq))
    }
    
    if '-' in locus_tag:
        product = f"deletion arm for region {start_tag} to {end_tag}"
        translation = ""
        protein_id = f"{start_tag}_{end_tag}_deletion_arm"
    else:
        product = target_row['product'].iloc[0] if 'product' in target_row.columns else f"deletion arm for {locus_tag}"
        translation = target_row['translation'].iloc[0] if 'translation' in target_row.columns else ""
        protein_id = target_row['protein_id'].iloc[0] if 'protein_id' in target_row.columns else f"{locus_tag}_deletion_arm"
    
    direction = '+'
    
    logging.debug(f"Target range: {target_start}..{target_end}")
    logging.debug(f"Upstream: {upstream_start}..{upstream_end}, Length: {len(upstream_seq)} bp")
    logging.debug(f"Downstream: {downstream_start}..{downstream_end}, Length: {len(downstream_seq)} bp")
    logging.info(f"Deletion arms fetched in {time.time() - start_time:.2f} seconds")
    return combined_seq, direction, product, translation, protein_id, cds_positions

def read_snapgene_dna(file_path):
    logging.info(f"Reading SnapGene file from {file_path}")
    start_time = time.time()
    snap_dict = snapgene_reader.snapgene_file_to_dict(file_path)
    seq = Seq.Seq(snap_dict['seq'])
    record = SeqRecord(seq, id=snap_dict.get('name', 'vector'), description='')
    logging.info(f"SnapGene file loaded in {time.time() - start_time:.2f} seconds")
    logging.info("Assuming vector promoter direction is 5'->3'")
    return record

def convert_dna_to_genbank(dna_file_path, locus_tag):
    logging.info(f"Converting {dna_file_path} to GenBank format for {locus_tag}")
    snapgene_path = "/Applications/SnapGene.app/Contents/MacOS/SnapGene"
    desktop_path = os.path.expanduser("~/Desktop")
    temp_gb_file = os.path.abspath(os.path.join(desktop_path, f"temp_genbank_{locus_tag}.gbk"))
    
    if os.path.exists(temp_gb_file):
        logging.info(f"{temp_gb_file} already exists, skipping conversion")
        return temp_gb_file
    
    cmd = f'{snapgene_path} --convert "GenBank - SnapGene" --input "{dna_file_path}" --output "{temp_gb_file}"'
    logging.debug(f"Executing command: {cmd}")
    print(f"Executing SnapGene CLI command: {cmd}")
    
    try:
        exit_code = os.system(cmd)
        if exit_code != 0:
            raise RuntimeError(f"SnapGene CLI failed with exit code {exit_code}")
        if not os.path.exists(temp_gb_file):
            raise FileNotFoundError(f"Output file {temp_gb_file} not generated")
        if os.path.getsize(temp_gb_file) == 0:
            raise ValueError(f"Output file {temp_gb_file} is empty")
        logging.info(f"Converted {dna_file_path} to {temp_gb_file}")
    except Exception as e:
        logging.error(f"Error during conversion: {str(e)}")
        if os.path.exists(temp_gb_file):
            os.remove(temp_gb_file)
        raise
    return temp_gb_file

def replace_sequence(vector_seq, start, end, insert_seq):
    start = int(start) - 1
    end = int(end)
    logging.info(f"Replacing sequence from {start + 1} to {end} with target sequence of length {len(insert_seq)}")
    start_time = time.time()
    if start < 0 or end > len(vector_seq) or start >= end:
        raise ValueError("Invalid start or end position.")
    new_seq = vector_seq[:start] + insert_seq + vector_seq[end:]
    logging.info(f"New sequence length: {len(new_seq)}. Sequence replaced in {time.time() - start_time:.2f} seconds")
    return new_seq

def check_self_dimerization(primer_seq, threshold=8):
    primer_seq = primer_seq.upper()
    rev_comp = str(Seq.Seq(primer_seq).reverse_complement())
    matches = nt_search(primer_seq, rev_comp[1:])
    for pos in matches:
        if isinstance(pos, int) and pos > 0 and len(primer_seq) - pos >= threshold:
            logging.warning(f"Self-dimerization detected in {primer_seq} at position {pos}")
            return True
    return False

def design_deletion_primers(vector_seq, insert_seq, start, end, upstream_length, downstream_length, tm_target=60, overlap_length=16, min_target_length=18, max_iterations=100):
    logging.info("Designing Deletion Gibson Assembly primers")
    start_time = time.time()
    start = int(start) - 1
    end = int(end)
    
    vector_seq = vector_seq.lower()
    insert_seq = insert_seq.upper()
    
    upstream_seq = insert_seq[:upstream_length]
    downstream_seq = insert_seq[upstream_length:]
    
    upstream_vector_overlap = str(vector_seq[start - overlap_length:start])
    downstream_vector_overlap = str(vector_seq[end:end + overlap_length])
    
    upstream_left = upstream_seq[:min_target_length]
    upstream_forward_primer = upstream_vector_overlap + upstream_left
    upstream_forward_target = upstream_left
    upstream_forward_tm_target = mt.Tm_NN(upstream_forward_target)
    upstream_forward_overlap_tm = mt.Tm_NN(upstream_vector_overlap)
    upstream_forward_tm_full = mt.Tm_NN(upstream_forward_primer)
    best_upstream_forward = {"primer": upstream_forward_primer, "tm_target": upstream_forward_tm_target, "tm_full": upstream_forward_tm_full, "diff": abs(upstream_forward_tm_target - tm_target)}
    
    iteration = 0
    while upstream_forward_tm_target < 55 and iteration < max_iterations:
        if len(upstream_left) < 25:
            upstream_forward_primer += upstream_seq[len(upstream_left):len(upstream_left) + 1]
            upstream_left += upstream_seq[len(upstream_left):len(upstream_left) + 1]
            upstream_forward_target = upstream_left
            upstream_forward_tm_target = mt.Tm_NN(upstream_forward_target)
            upstream_forward_tm_full = mt.Tm_NN(upstream_forward_primer)
            diff = abs(upstream_forward_tm_target - tm_target)
            if upstream_forward_tm_target > 60:
                upstream_forward_primer = upstream_vector_overlap + upstream_left[:-1]
                upstream_left = upstream_left[:-1]
                upstream_forward_target = upstream_left
                upstream_forward_tm_target = mt.Tm_NN(upstream_forward_target)
                upstream_forward_tm_full = mt.Tm_NN(upstream_forward_primer)
                break
            if not check_self_dimerization(upstream_forward_primer):
                best_upstream_forward = {"primer": upstream_forward_primer, "tm_target": upstream_forward_tm_target, "tm_full": upstream_forward_tm_full, "diff": diff}
        iteration += 1
    
    if upstream_forward_overlap_tm < 55 and len(upstream_vector_overlap) < 20:
        while upstream_forward_overlap_tm < 55 and len(upstream_vector_overlap) < 20:
            upstream_vector_overlap = str(vector_seq[start - (len(upstream_vector_overlap) + 1):start])
            upstream_forward_primer = upstream_vector_overlap + upstream_left
            upstream_forward_overlap_tm = mt.Tm_NN(upstream_vector_overlap)
            upstream_forward_tm_full = mt.Tm_NN(upstream_forward_primer)
            if not check_self_dimerization(upstream_forward_primer):
                best_upstream_forward = {"primer": upstream_forward_primer, "tm_target": upstream_forward_tm_target, "tm_full": upstream_forward_tm_full, "diff": abs(upstream_forward_tm_target - tm_target)}
    
    downstream_right = downstream_seq[-min_target_length:]
    downstream_reverse_primer = str(Seq.Seq(downstream_right + downstream_vector_overlap).reverse_complement())
    downstream_reverse_target = str(Seq.Seq(downstream_right).reverse_complement())
    downstream_reverse_tm_target = mt.Tm_NN(downstream_reverse_target)
    downstream_reverse_overlap_tm = mt.Tm_NN(downstream_vector_overlap)
    downstream_reverse_tm_full = mt.Tm_NN(downstream_reverse_primer)
    best_downstream_reverse = {"primer": downstream_reverse_primer, "tm_target": downstream_reverse_tm_target, "tm_full": downstream_reverse_tm_full, "diff": abs(downstream_reverse_tm_target - tm_target)}
    
    iteration = 0
    while downstream_reverse_tm_target < 55 and iteration < max_iterations:
        if len(downstream_right) < 25:
            downstream_right = downstream_seq[-len(downstream_right) - 1:-len(downstream_right)] + downstream_right
            downstream_reverse_primer = str(Seq.Seq(downstream_right + downstream_vector_overlap).reverse_complement())
            downstream_reverse_target = str(Seq.Seq(downstream_right).reverse_complement())
            downstream_reverse_tm_target = mt.Tm_NN(downstream_reverse_target)
            downstream_reverse_tm_full = mt.Tm_NN(downstream_reverse_primer)
            diff = abs(downstream_reverse_tm_target - tm_target)
            if downstream_reverse_tm_target > 60:
                downstream_right = downstream_right[1:]
                downstream_reverse_primer = str(Seq.Seq(downstream_right + downstream_vector_overlap).reverse_complement())
                downstream_reverse_target = str(Seq.Seq(downstream_right).reverse_complement())
                downstream_reverse_tm_target = mt.Tm_NN(downstream_reverse_target)
                downstream_reverse_tm_full = mt.Tm_NN(downstream_reverse_primer)
                break
            if not check_self_dimerization(downstream_reverse_primer):
                best_downstream_reverse = {"primer": downstream_reverse_primer, "tm_target": downstream_reverse_tm_target, "tm_full": downstream_reverse_tm_full, "diff": diff}
        iteration += 1
    
    if downstream_reverse_overlap_tm < 55 and len(downstream_vector_overlap) < 20:
        while downstream_reverse_overlap_tm < 55 and len(downstream_vector_overlap) < 20:
            downstream_vector_overlap = str(vector_seq[end:end + len(downstream_vector_overlap) + 1])
            downstream_reverse_primer = str(Seq.Seq(downstream_right + downstream_vector_overlap).reverse_complement())
            downstream_reverse_overlap_tm = mt.Tm_NN(downstream_vector_overlap)
            downstream_reverse_tm_full = mt.Tm_NN(downstream_reverse_primer)
            if not check_self_dimerization(downstream_reverse_primer):
                best_downstream_reverse = {"primer": downstream_reverse_primer, "tm_target": downstream_reverse_tm_target, "tm_full": downstream_reverse_tm_full, "diff": abs(downstream_reverse_tm_target - tm_target)}
    
    min_target_length = 18
    max_target_length = 50
    max_primer_length = 70
    target_tm_min = 55
    target_tm_max = 65
    target_tm_optimal = 60
    
    upstream_right = upstream_seq[-min_target_length:]
    upstream_reverse_target = str(Seq.Seq(upstream_right).reverse_complement())
    upstream_reverse_tm_target = mt.Tm_NN(upstream_reverse_target)
    iteration = 0
    best_upstream_reverse_target = upstream_right
    best_upstream_reverse_tm_target = upstream_reverse_tm_target
    
    while (upstream_reverse_tm_target < target_tm_min or abs(upstream_reverse_tm_target - target_tm_optimal) > 1) and iteration < max_iterations:
        if len(upstream_right) < max_target_length:
            if upstream_reverse_tm_target < target_tm_min:
                upstream_right = upstream_seq[-(len(upstream_right) + 1):] + upstream_right
            elif upstream_reverse_tm_target > target_tm_max:
                upstream_right = upstream_right[1:]
            else:
                if upstream_reverse_tm_target < target_tm_optimal and len(upstream_right) < max_target_length:
                    upstream_right = upstream_seq[-(len(upstream_right) + 1):] + upstream_right
                elif upstream_reverse_tm_target > target_tm_optimal and len(upstream_right) > min_target_length:
                    upstream_right = upstream_right[1:]
                else:
                    break
            upstream_reverse_target = str(Seq.Seq(upstream_right).reverse_complement())
            upstream_reverse_tm_target = mt.Tm_NN(upstream_reverse_target)
            if target_tm_min <= upstream_reverse_tm_target <= target_tm_max:
                best_upstream_reverse_target = upstream_right
                best_upstream_reverse_tm_target = upstream_reverse_tm_target
            logging.debug(f"Upstream_Reverse iteration {iteration}: Length: {len(upstream_right)} bp, Tm: {upstream_reverse_tm_target:.2f}°C")
        else:
            break
        iteration += 1
    
    downstream_left = downstream_seq[:min_target_length]
    downstream_forward_tm_target = mt.Tm_NN(downstream_left)
    iteration = 0
    best_downstream_forward_target = downstream_left
    best_downstream_forward_tm_target = downstream_forward_tm_target
    
    while (downstream_forward_tm_target < target_tm_min or abs(downstream_forward_tm_target - target_tm_optimal) > 1) and iteration < max_iterations:
        if len(downstream_left) < max_target_length:
            if downstream_forward_tm_target < target_tm_min:
                downstream_left += downstream_seq[len(downstream_left)]
            elif downstream_forward_tm_target > target_tm_max:
                downstream_left = downstream_left[:-1]
            else:
                if downstream_forward_tm_target < target_tm_optimal and len(downstream_left) < max_target_length:
                    downstream_left += downstream_seq[len(downstream_left)]
                elif downstream_forward_tm_target > target_tm_optimal and len(downstream_left) > min_target_length:
                    downstream_left = downstream_left[:-1]
                else:
                    break
            downstream_forward_tm_target = mt.Tm_NN(downstream_left)
            if target_tm_min <= downstream_forward_tm_target <= target_tm_max:
                best_downstream_forward_target = downstream_left
                best_downstream_forward_tm_target = downstream_forward_tm_target
            logging.debug(f"Downstream_Forward iteration {iteration}: Length: {len(downstream_left)} bp, Tm: {downstream_forward_tm_target:.2f}°C")
        else:
            break
        iteration += 1
    
    max_overlap_length = min(20, max_primer_length - len(best_upstream_reverse_target), max_primer_length - len(best_downstream_forward_target))
    tail_upstream_reverse_length = min(10, max_overlap_length)
    tail_downstream_forward_length = min(10, max_overlap_length)
    tail_upstream_reverse = downstream_seq[:tail_upstream_reverse_length]
    tail_downstream_forward = upstream_seq[-tail_downstream_forward_length:]
    
    upstream_reverse_primer = str(Seq.Seq(best_upstream_reverse_target + tail_upstream_reverse).reverse_complement())
    downstream_forward_primer = tail_downstream_forward + best_downstream_forward_target
    
    upstream_reverse_tm_full = mt.Tm_NN(upstream_reverse_primer)
    downstream_forward_tm_full = mt.Tm_NN(downstream_forward_primer)
    overlap_combined = tail_upstream_reverse + tail_downstream_forward
    overlap_tm = mt.Tm_NN(overlap_combined)
    overlap_length_sum = len(tail_upstream_reverse) + len(tail_downstream_forward)
    
    best_upstream_reverse = {"primer": upstream_reverse_primer, "tm_target": best_upstream_reverse_tm_target, "tm_full": upstream_reverse_tm_full, "overlap_tm": overlap_tm}
    best_downstream_forward = {"primer": downstream_forward_primer, "tm_target": best_downstream_forward_tm_target, "tm_full": downstream_forward_tm_full, "overlap_tm": overlap_tm}
    
    iteration = 0
    while (overlap_length_sum < 20 or overlap_tm < 55) and iteration < max_iterations:
        if len(upstream_reverse_primer) < max_primer_length and len(downstream_forward_primer) < max_primer_length:
            if iteration % 2 == 0 and len(tail_upstream_reverse) < max_overlap_length:
                tail_upstream_reverse_length += 2
                tail_upstream_reverse = downstream_seq[:tail_upstream_reverse_length]
            elif len(tail_downstream_forward) < max_overlap_length:
                tail_downstream_forward_length += 2
                tail_downstream_forward = upstream_seq[-tail_downstream_forward_length:]
            
            upstream_reverse_primer = str(Seq.Seq(best_upstream_reverse_target + tail_upstream_reverse).reverse_complement())
            downstream_forward_primer = tail_downstream_forward + best_downstream_forward_target
            
            upstream_reverse_tm_full = mt.Tm_NN(upstream_reverse_primer)
            downstream_forward_tm_full = mt.Tm_NN(downstream_forward_primer)
            overlap_combined = tail_upstream_reverse + tail_downstream_forward
            overlap_tm = mt.Tm_NN(overlap_combined)
            overlap_length_sum = len(tail_upstream_reverse) + len(tail_downstream_forward)
            
            if len(upstream_reverse_primer) <= max_primer_length and len(downstream_forward_primer) <= max_primer_length and not check_self_dimerization(upstream_reverse_primer) and not check_self_dimerization(downstream_forward_primer):
                best_upstream_reverse = {"primer": upstream_reverse_primer, "tm_target": best_upstream_reverse_tm_target, "tm_full": upstream_reverse_tm_full, "overlap_tm": overlap_tm}
                best_downstream_forward = {"primer": downstream_forward_primer, "tm_target": best_downstream_forward_tm_target, "tm_full": downstream_forward_tm_full, "overlap_tm": overlap_tm}
                if overlap_tm >= 55 and overlap_length_sum >= 20:
                    break
        iteration += 1
    
    logging.debug(f"Upstream_Reverse target: {best_upstream_reverse_target}, Length: {len(best_upstream_reverse_target)} bp, Tm: {best_upstream_reverse_tm_target:.2f}°C, GC: {GC(best_upstream_reverse_target):.2f}%")
    logging.debug(f"Upstream_Reverse full: {upstream_reverse_primer}, Length: {len(upstream_reverse_primer)} bp, Tm: {upstream_reverse_tm_full:.2f}°C")
    logging.debug(f"Downstream_Forward target: {best_downstream_forward_target}, Length: {len(best_downstream_forward_target)} bp, Tm: {best_downstream_forward_tm_target:.2f}°C, GC: {GC(best_downstream_forward_target):.2f}%")
    logging.debug(f"Downstream_Forward full: {downstream_forward_primer}, Length: {len(downstream_forward_primer)} bp, Tm: {downstream_forward_tm_full:.2f}°C")
    logging.debug(f"Final Overlap sequence (5'→3'): {overlap_combined}, Length: {overlap_length_sum} bp, Tm: {overlap_tm:.2f}°C, GC: {GC(overlap_combined):.2f}%")
    
    primers = {
        "upstream_forward_primer": best_upstream_forward["primer"],
        "upstream_forward_tm_target": best_upstream_forward["tm_target"],
        "upstream_forward_tm_full": best_upstream_forward["tm_full"],
        "upstream_reverse_primer": best_upstream_reverse["primer"],
        "upstream_reverse_tm_target": best_upstream_reverse["tm_target"],
        "upstream_reverse_tm_full": best_upstream_reverse["tm_full"],
        "downstream_forward_primer": best_downstream_forward["primer"],
        "downstream_forward_tm_target": best_downstream_forward["tm_target"],
        "downstream_forward_tm_full": best_downstream_forward["tm_full"],
        "downstream_reverse_primer": best_downstream_reverse["primer"],
        "downstream_reverse_tm_target": best_downstream_reverse["tm_target"],
        "downstream_reverse_tm_full": best_downstream_reverse["tm_full"],
        "upstream_forward_start": start - len(upstream_vector_overlap) + 1,
        "upstream_forward_end": start + len(upstream_left),
        "upstream_reverse_start": start + upstream_length - len(best_upstream_reverse_target) + 1,
        "upstream_reverse_end": start + upstream_length + len(tail_upstream_reverse),
        "downstream_forward_start": start + upstream_length - len(tail_downstream_forward) + 1,
        "downstream_forward_end": start + upstream_length + len(best_downstream_forward_target),
        "downstream_reverse_start": start + len(insert_seq) - len(downstream_right) + 1,
        "downstream_reverse_end": start + len(insert_seq) + len(downstream_vector_overlap)
    }
    
    logging.info(f"Primers designed in {time.time() - start_time:.2f} seconds")
    return primers

def adjust_feature_location(line, start, end, length_diff):
    if not line.strip() or line.strip().startswith("/"):
        return line
    
    parts = line.strip().split(maxsplit=1)
    if len(parts) < 2:
        return line
    
    feature_type = parts[0]
    location = parts[1].split()[0]
    logging.debug(f"Original feature: {feature_type} at {location}")
    
    start_1based = start - 1
    end_1based = end
    
    if "join(" in location:
        if "complement(join(" in location:
            loc_parts = location.split("complement(join(")[1].rstrip("))").split(",")
            adjusted_locs = []
            for loc in loc_parts:
                loc_start, loc_end = map(int, loc.split(".."))
                if loc_end <= start_1based:
                    adjusted_locs.append(f"{loc_start}..{loc_end}")
                elif loc_start >= end_1based:
                    loc_start += length_diff
                    loc_end += length_diff
                    adjusted_locs.append(f"{loc_start}..{loc_end}")
                else:
                    logging.debug(f"Removing feature {feature_type} at {loc} overlapping insertion {start_1based + 1}..{end_1based}")
                    return None
            if adjusted_locs:
                new_location = f"complement(join({','.join(adjusted_locs)}))"
                logging.debug(f"Adjusted feature: {feature_type} at {new_location}")
                return f"     {feature_type:<{FEATURE_WIDTH}} {new_location}\n"
            return None
        else:
            loc_parts = location.split("join(")[1].rstrip(")").split(",")
            adjusted_locs = []
            for loc in loc_parts:
                loc_start, loc_end = map(int, loc.split(".."))
                if loc_end <= start_1based:
                    adjusted_locs.append(f"{loc_start}..{loc_end}")
                elif loc_start >= end_1based:
                    loc_start += length_diff
                    loc_end += length_diff
                    adjusted_locs.append(f"{loc_start}..{loc_end}")
                else:
                    logging.debug(f"Removing feature {feature_type} at {loc} overlapping insertion {start_1based + 1}..{end_1based}")
                    return None
            if adjusted_locs:
                new_location = f"join({','.join(adjusted_locs)})"
                logging.debug(f"Adjusted feature: {feature_type} at {new_location}")
                return f"     {feature_type:<{FEATURE_WIDTH}} {new_location}\n"
            return None
    
    elif "complement(" in location:
        loc_start, loc_end = map(int, location.split("complement(")[1].rstrip(")").split(".."))
        if loc_end <= start_1based:
            new_location = f"complement({loc_start}..{loc_end})"
        elif loc_start >= end_1based:
            loc_start += length_diff
            loc_end += length_diff
            new_location = f"complement({loc_start}..{loc_end})"
        else:
            logging.debug(f"Removing feature {feature_type} at {location} overlapping insertion {start_1based + 1}..{end_1based}")
            return None
        logging.debug(f"Adjusted feature: {feature_type} at {new_location}")
        return f"     {feature_type:<{FEATURE_WIDTH}} {new_location}\n"
    
    elif ".." in location:
        loc_start, loc_end = map(int, location.split(".."))
        if loc_end <= start_1based:
            new_location = f"{loc_start}..{loc_end}"
        elif loc_start >= end_1based:
            loc_start += length_diff
            loc_end += length_diff
            new_location = f"{loc_start}..{loc_end}"
        else:
            logging.debug(f"Removing feature {feature_type} at {location} overlapping insertion {start_1based + 1}..{end_1based}")
            return None
        logging.debug(f"Adjusted feature: {feature_type} at {new_location}")
        return f"     {feature_type:<{FEATURE_WIDTH}} {new_location}\n"
    
    elif "^" in location:
        loc_start, loc_end = map(int, location.split("^"))
        if loc_end <= start_1based:
            new_location = f"{loc_start}^{loc_end}"
        elif loc_start >= end_1based:
            loc_start += length_diff
            loc_end += length_diff
            new_location = f"{loc_start}^{loc_end}"
        else:
            logging.debug(f"Removing feature {feature_type} at {location} overlapping insertion {start_1based + 1}..{end_1based}")
            return None
        logging.debug(f"Adjusted feature: {feature_type} at {new_location}")
        return f"     {feature_type:<{FEATURE_WIDTH}} {new_location}\n"
    
    elif location.isdigit():
        loc = int(location)
        if loc <= start_1based:
            new_location = str(loc)
        elif loc >= end_1based:
            loc += length_diff
            new_location = str(loc)
        else:
            logging.debug(f"Removing feature {feature_type} at {location} overlapping insertion {start_1based + 1}..{end_1based}")
            return None
        logging.debug(f"Adjusted feature: {feature_type} at {new_location}")
        return f"     {feature_type:<{FEATURE_WIDTH}} {new_location}\n"
    
    return line

def format_long_qualifier(value, first_line_prefix, subsequent_prefix, use_quotes=False, is_target_translation=False):
    lines = []
    
    if is_target_translation and use_quotes:
        value = value.replace('\n', '')
        remaining_value = value
        first_line_done = False
        
        while remaining_value:
            if not first_line_done:
                prefix = first_line_prefix
                max_width = MAX_LINE_WIDTH - len(first_line_prefix) - 1
                split_value = remaining_value[:max_width]
                lines.append(f"{prefix}\"{split_value}\n")
                remaining_value = remaining_value[max_width:]
                first_line_done = True
            else:
                prefix = subsequent_prefix
                max_width = MAX_LINE_WIDTH - QUALIFIER_START
                if len(remaining_value) > max_width:
                    split_value = remaining_value[:max_width]
                    lines.append(f"{prefix}{split_value}\n")
                    remaining_value = remaining_value[max_width:]
                else:
                    lines.append(f"{prefix}{remaining_value}\"\n")
                    break
    else:
        original_lines = value.split('\n')
        for i, line in enumerate(original_lines):
            if i == 0:
                prefix = first_line_prefix
                max_width = MAX_LINE_WIDTH - len(first_line_prefix)
            else:
                prefix = subsequent_prefix
                max_width = MAX_LINE_WIDTH - QUALIFIER_START
            
            remaining_value = line.strip()
            while remaining_value:
                if len(remaining_value) > max_width:
                    split_value = remaining_value[:max_width]
                    lines.append(f"{prefix}{split_value}\n")
                    remaining_value = remaining_value[max_width:]
                else:
                    lines.append(f"{prefix}{remaining_value}\n")
                    break
                prefix = subsequent_prefix
    
    return lines

def modify_genbank(temp_gb_file, new_seq, start, end, locus_tag, target_seq, product, translation, protein_id, primers, output_path, cds_positions, direction):
    logging.info(f"Modifying GenBank file {temp_gb_file} to {output_path}")
    start_time = time.time()
    
    if not os.path.exists(temp_gb_file):
        raise FileNotFoundError(f"File {temp_gb_file} does not exist")
    with open(temp_gb_file, 'r') as f:
        lines = f.readlines()
    
    header_lines = []
    feature_lines = []
    origin_lines = []
    in_features = False
    in_origin = False
    
    for line in lines:
        if line.startswith("FEATURES"):
            in_features = True
            header_lines.append(line)
            continue
        elif line.startswith("ORIGIN"):
            in_features = False
            in_origin = True
            origin_lines.append(line)
            continue
        if in_features:
            feature_lines.append(line)
        elif in_origin:
            origin_lines.append(line)
        else:
            header_lines.append(line)
    
    new_length = len(new_seq)
    for i, line in enumerate(header_lines):
        if line.startswith("LOCUS"):
            parts = line.split()
            parts[1] = f"{parts[1]}_{locus_tag}_deletion"
            parts[2] = f"{new_length} bp"
            header_lines[i] = " ".join(parts) + "\n"
        elif line.strip().startswith("REFERENCE") and "(bases" in line:
            parts = line.split("(bases")
            parts[1] = f" (bases 1 to {new_length})\n"
            header_lines[i] = "".join(parts)
    
    length_diff = len(target_seq) - (end - (start - 1))
    logging.debug(f"Length difference for feature adjustment: {length_diff}")
    adjusted_features = []
    current_feature_lines = []
    
    for line in feature_lines:
        if line.strip() and not line.strip().startswith("/"):
            if current_feature_lines:
                adjusted_line = adjust_feature_location(current_feature_lines[0], start, end, length_diff)
                if adjusted_line:
                    adjusted_features.append(adjusted_line)
                    qualifiers = []
                    for qual_line in current_feature_lines[1:]:
                        qual_content = qual_line.strip()
                        if qual_content:
                            key, value = qual_content.split("=", 1)
                            qualifiers.append((key, value))
                    for i, (key, value) in enumerate(qualifiers):
                        prefix = f"{' ' * QUALIFIER_START}{key}="
                        subsequent_prefix = " " * QUALIFIER_START
                        use_quotes = False
                        formatted_lines = format_long_qualifier(value, prefix, subsequent_prefix, use_quotes)
                        adjusted_features.extend(formatted_lines)
            current_feature_lines = [line]
        else:
            current_feature_lines.append(line)
    
    if current_feature_lines:
        adjusted_line = adjust_feature_location(current_feature_lines[0], start, end, length_diff)
        if adjusted_line:
            adjusted_features.append(adjusted_line)
            qualifiers = []
            for qual_line in current_feature_lines[1:]:
                qual_content = qual_line.strip()
                if qual_content:
                    key, value = qual_content.split("=", 1)
                    qualifiers.append((key, value))
            for i, (key, value) in enumerate(qualifiers):
                prefix = f"{' ' * QUALIFIER_START}{key}="
                subsequent_prefix = " " * QUALIFIER_START
                use_quotes = False
                formatted_lines = format_long_qualifier(value, prefix, subsequent_prefix, use_quotes)
                adjusted_features.extend(formatted_lines)
    
    locus_tags = [f"{locus_tag}_upstream", f"{locus_tag}_downstream"]
    products = [f"upstream arm for {locus_tag} deletion", f"downstream arm for {locus_tag} deletion"]
    translations = ["", ""]
    protein_ids = [f"{locus_tag}_upstream_arm", f"{locus_tag}_downstream_arm"]
    
    for i, tag in enumerate(locus_tags):
        rel_start, rel_end = cds_positions[tag]
        abs_start = start + rel_start
        abs_end = start + rel_end - 1
        adjusted_features.append(f"     misc_feature    {abs_start}..{abs_end}\n")
        
        qualifiers = [
            ("/label", tag),
            ("/note", products[i]),
            ("/protein_id", protein_ids[i]),
            ("/translation", translations[i])
        ]
        for key, value in qualifiers:
            prefix = f"{' ' * QUALIFIER_START}{key}="
            subsequent_prefix = " " * QUALIFIER_START
            use_quotes = (key == "/translation")
            is_target_translation = (key == "/translation")
            formatted_lines = format_long_qualifier(value, prefix, subsequent_prefix, use_quotes, is_target_translation)
            adjusted_features.extend(formatted_lines)
    
    adjusted_features.append(f"     primer          {primers['upstream_forward_start']}..{primers['upstream_forward_end']}\n")
    qualifiers = [
        ("/label", "Upstream_Forward"),
        ("/note", f"Tm (Target): {primers['upstream_forward_tm_target']:.2f}°C, Tm (Full): {primers['upstream_forward_tm_full']:.2f}°C")
    ]
    for key, value in qualifiers:
        prefix = f"{' ' * QUALIFIER_START}{key}="
        subsequent_prefix = " " * QUALIFIER_START
        use_quotes = False
        formatted_lines = format_long_qualifier(value, prefix, subsequent_prefix, use_quotes)
        adjusted_features.extend(formatted_lines)
    
    adjusted_features.append(f"     primer          complement({primers['upstream_reverse_start']}..{primers['upstream_reverse_end']})\n")
    qualifiers = [
        ("/label", "Upstream_Reverse"),
        ("/note", f"Tm (Target): {primers['upstream_reverse_tm_target']:.2f}°C, Tm (Full): {primers['upstream_reverse_tm_full']:.2f}°C")
    ]
    for key, value in qualifiers:
        prefix = f"{' ' * QUALIFIER_START}{key}="
        subsequent_prefix = " " * QUALIFIER_START
        use_quotes = False
        formatted_lines = format_long_qualifier(value, prefix, subsequent_prefix, use_quotes)
        adjusted_features.extend(formatted_lines)
    
    adjusted_features.append(f"     primer          {primers['downstream_forward_start']}..{primers['downstream_forward_end']}\n")
    qualifiers = [
        ("/label", "Downstream_Forward"),
        ("/note", f"Tm (Target): {primers['downstream_forward_tm_target']:.2f}°C, Tm (Full): {primers['downstream_forward_tm_full']:.2f}°C")
    ]
    for key, value in qualifiers:
        prefix = f"{' ' * QUALIFIER_START}{key}="
        subsequent_prefix = " " * QUALIFIER_START
        use_quotes = False
        formatted_lines = format_long_qualifier(value, prefix, subsequent_prefix, use_quotes)
        adjusted_features.extend(formatted_lines)
    
    adjusted_features.append(f"     primer          complement({primers['downstream_reverse_start']}..{primers['downstream_reverse_end']})\n")
    qualifiers = [
        ("/label", "Downstream_Reverse"),
        ("/note", f"Tm (Target): {primers['downstream_reverse_tm_target']:.2f}°C, Tm (Full): {primers['downstream_reverse_tm_full']:.2f}°C")
    ]
    for key, value in qualifiers:
        prefix = f"{' ' * QUALIFIER_START}{key}="
        subsequent_prefix = " " * QUALIFIER_START
        use_quotes = False
        formatted_lines = format_long_qualifier(value, prefix, subsequent_prefix, use_quotes)
        adjusted_features.extend(formatted_lines)
    
    origin_lines = ["ORIGIN\n"]
    seq_str = str(new_seq).lower()
    for i in range(0, len(seq_str), 60):
        chunk = seq_str[i:i+60]
        formatted_chunk = " ".join(chunk[j:j+10] for j in range(0, len(chunk), 10))
        pos = i + 1
        origin_lines.append(f"{' ' * (9 - len(str(pos)))}{pos} {formatted_chunk}\n")
    origin_lines.append("//\n")
    
    new_lines = header_lines + adjusted_features + origin_lines
    
    output_path = os.path.abspath(output_path)
    with open(output_path, 'w') as f:
        f.writelines(new_lines)
    
    logging.info(f"GenBank file modified in {time.time() - start_time:.2f} seconds")
    logging.debug(f"Final adjusted features: {adjusted_features}")

def main():
    parser = argparse.ArgumentParser(description="Design Gibson Assembly primers for deletion vector construction for multiple targets.")
    parser.add_argument("--genbank_file", required=True, help="Path to the GenBank file containing the target genome sequence")
    parser.add_argument("--genbank_table", required=True, help="Path to the table file (TSV or XLSX) for locus_tag position and metadata")
    parser.add_argument("--locus_tag", required=True, help="Target locus_tag(s) to delete (e.g., 'A, B' or 'A-B' or single 'A')")
    parser.add_argument("--upstream_bp", type=int, required=True, help="Number of base pairs upstream of the target locus_tag")
    parser.add_argument("--downstream_bp", type=int, required=True, help="Number of base pairs downstream of the target locus_tag")
    parser.add_argument("--vector_file", required=True, help="Path to the SnapGene DNA vector file (.dna)")
    parser.add_argument("--start", required=True, type=int, help="Start position in vector for insertion (1-based)")
    parser.add_argument("--end", required=True, type=int, help="End position in vector for insertion (1-based)")
    parser.add_argument("--output_dir", default="deletion_results", help="Directory to save output files")
    parser.add_argument("--tm", type=float, default=60, help="Target Tm for primers (default: 60°C, min 55°C)")
    
    args = parser.parse_args()
    
    args.genbank_file = os.path.abspath(args.genbank_file)
    args.genbank_table = os.path.abspath(args.genbank_table)
    args.vector_file = os.path.abspath(args.vector_file)
    args.output_dir = os.path.abspath(args.output_dir)
    
    # locus_tag을 쉼표로 분리
    locus_tags = [tag.strip() for tag in args.locus_tag.split(',')]
    
    genbank_records = read_genbank_file(args.genbank_file)
    genbank_df = read_genbank_table(args.genbank_table)
    vector_record = read_snapgene_dna(args.vector_file)
    vector_seq = vector_record.seq
    
    os.makedirs(args.output_dir, exist_ok=True)
    excel_log.clear()  # 여러 타겟을 위해 로그 초기화
    
    for locus_tag in locus_tags:
        print(f"\nProcessing locus_tag: {locus_tag}")
        
        target_seq, direction, product, translation, protein_id, cds_positions = get_deletion_arms(
            genbank_records, genbank_df, locus_tag, args.upstream_bp, args.downstream_bp
        )
        
        new_seq = replace_sequence(vector_seq, args.start, args.end, target_seq)
        
        primers = design_deletion_primers(
            vector_seq, target_seq, args.start, args.end, args.upstream_bp, args.downstream_bp, tm_target=args.tm
        )
        
        print(f"Results for deletion of locus_tag: {locus_tag}")
        print("Upstream Forward Primer (Vector 5' + Upstream 5'):", primers["upstream_forward_primer"])
        print("Upstream Forward Tm (Target):", primers["upstream_forward_tm_target"])
        print("Upstream Forward Tm (Full):", primers["upstream_forward_tm_full"])
        print("Upstream Reverse Primer (Upstream 3' rev_comp + Downstream 5' rev_comp):", primers["upstream_reverse_primer"])
        print("Upstream Reverse Tm (Target):", primers["upstream_reverse_tm_target"])
        print("Upstream Reverse Tm (Full):", primers["upstream_reverse_tm_full"])
        print("Downstream Forward Primer (Upstream 3' + Downstream 5'):", primers["downstream_forward_primer"])
        print("Downstream Forward Tm (Target):", primers["downstream_forward_tm_target"])
        print("Downstream Forward Tm (Full):", primers["downstream_forward_tm_full"])
        print("Downstream Reverse Primer (Downstream 3' rev_comp + Vector 3' rev_comp):", primers["downstream_reverse_primer"])
        print("Downstream Reverse Tm (Target):", primers["downstream_reverse_tm_target"])
        print("Downstream Reverse Tm (Full):", primers["downstream_reverse_tm_full"])
        
        excel_log.append({
            "Locus Tag": locus_tag,
            "Upstream Forward Primer": primers["upstream_forward_primer"],
            "Upstream Forward Tm (Target)": primers["upstream_forward_tm_target"],
            "Upstream Forward Tm (Full)": primers["upstream_forward_tm_full"],
            "Upstream Reverse Primer": primers["upstream_reverse_primer"],
            "Upstream Reverse Tm (Target)": primers["upstream_reverse_tm_target"],
            "Upstream Reverse Tm (Full)": primers["upstream_reverse_tm_full"],
            "Downstream Forward Primer": primers["downstream_forward_primer"],
            "Downstream Forward Tm (Target)": primers["downstream_forward_tm_target"],
            "Downstream Forward Tm (Full)": primers["downstream_forward_tm_full"],
            "Downstream Reverse Primer": primers["downstream_reverse_primer"],
            "Downstream Reverse Tm (Target)": primers["downstream_reverse_tm_target"],
            "Downstream Reverse Tm (Full)": primers["downstream_reverse_tm_full"]
        })
        
        temp_gb_file = convert_dna_to_genbank(args.vector_file, f"{locus_tag}_deletion")
        output_path = os.path.join(args.output_dir, f"{locus_tag}_deletion_vector.gbk")
        modify_genbank(
            temp_gb_file, new_seq, args.start, args.end, locus_tag, target_seq, product, translation, 
            protein_id, primers, output_path, cds_positions, direction
        )
        
        print(f"Updated GenBank file saved to: {output_path}")
    
    # 모든 타겟 처리 후 로그 저장
    log_df = pd.DataFrame(excel_log)
    log_file = os.path.join(args.output_dir, "log_deletions.xlsx")
    log_df.to_excel(log_file, index=False)
    print(f"Combined log saved to: {log_file}")

if __name__ == "__main__":
    main()