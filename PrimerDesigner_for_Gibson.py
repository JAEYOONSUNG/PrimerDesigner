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

# Set up console logging (including DEBUG level)
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

# List to store Excel log data
excel_log = []

# GenBank format settings
MAX_LINE_WIDTH = 80  # Maximum line width for GenBank format
QUALIFIER_START = 21  # Starting column for qualifiers (21 spaces)
FEATURE_WIDTH = 15    # Maximum width for feature names

def read_genbank_table(file_path):
    # Log the start of reading the GenBank table
    logging.info(f"Reading GenBank table from {file_path}")
    start_time = time.time()
    if file_path.endswith('.xlsx'):
        df = pd.read_excel(file_path)
    else:
        df = pd.read_csv(file_path, sep='\t')
    # Log the completion of loading the table
    logging.info(f"GenBank table loaded in {time.time() - start_time:.2f} seconds")
    return df

def get_target_sequence(genbank_df, locus_tag):
    # Log the start of fetching sequence data
    logging.info(f"Fetching sequence and details for locus_tag: {locus_tag}")
    start_time = time.time()
    target_row = genbank_df[genbank_df['locus_tag'] == locus_tag]
    if target_row.empty:
        raise ValueError(f"Locus tag {locus_tag} not found in GenBank table.")
    seq = target_row['rearranged_nt_seq'].iloc[0]
    direction = target_row['direction'].iloc[0]
    product = target_row['product'].iloc[0] if 'product' in target_row.columns else "unknown protein"
    translation = target_row['translation'].iloc[0] if 'translation' in target_row.columns else ""
    protein_id = target_row['protein_id'].iloc[0] if 'protein_id' in target_row.columns else f"{locus_tag}_protein"
    # Log the completion of sequence fetching
    logging.info(f"Sequence fetched in {time.time() - start_time:.2f} seconds")
    return seq, direction, product, translation, protein_id

def read_snapgene_dna(file_path):
    # Log the start of reading the SnapGene file
    logging.info(f"Reading SnapGene file from {file_path}")
    start_time = time.time()
    snap_dict = snapgene_reader.snapgene_file_to_dict(file_path)
    seq = Seq.Seq(snap_dict['seq'])
    record = SeqRecord(seq, id=snap_dict.get('name', 'vector'), description='')
    # Log the completion of loading the SnapGene file
    logging.info(f"SnapGene file loaded in {time.time() - start_time:.2f} seconds")
    logging.info("Assuming vector promoter direction is 5'->3'")
    return record

def convert_dna_to_genbank(dna_file_path, locus_tag):
    # Log the start of converting the DNA file to GenBank format
    logging.info(f"Converting {dna_file_path} to GenBank format")
    snapgene_path = "/Applications/SnapGene.app/Contents/MacOS/SnapGene"
    desktop_path = os.path.expanduser("~/Desktop")
    temp_gb_file = os.path.abspath(os.path.join(desktop_path, f"temp_genbank_{locus_tag}.gbk"))
    
    if os.path.exists(temp_gb_file):
        # Skip conversion if the file already exists
        logging.info(f"{temp_gb_file} already exists, skipping conversion")
        return temp_gb_file
    
    cmd = f'{snapgene_path} --convert "GenBank - SnapGene" --input "{dna_file_path}" --output "{temp_gb_file}"'
    # Log the command being executed
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
        # Log successful conversion
        logging.info(f"Converted {dna_file_path} to {temp_gb_file}")
    except Exception as e:
        # Log any errors during conversion
        logging.error(f"Error during conversion: {str(e)}")
        if os.path.exists(temp_gb_file):
            os.remove(temp_gb_file)
        raise
    return temp_gb_file

def replace_sequence(vector_seq, start, end, insert_seq):
    start, end = int(start) - 1, int(end)
    # Log the start of sequence replacement
    logging.info(f"Replacing sequence from {start + 1} to {end} with target sequence of length {len(insert_seq)}")
    start_time = time.time()
    if start < 0 or end > len(vector_seq) or start >= end:
        raise ValueError("Invalid start or end position.")
    new_seq = vector_seq[:start] + insert_seq + vector_seq[end:]
    # Log the completion of sequence replacement
    logging.info(f"New sequence length: {len(new_seq)}. Sequence replaced in {time.time() - start_time:.2f} seconds")
    return new_seq

def check_self_dimerization(primer_seq, threshold=8):
    primer_seq = primer_seq.upper()
    rev_comp = str(Seq.Seq(primer_seq).reverse_complement())
    matches = nt_search(primer_seq, rev_comp[1:])
    for pos in matches:
        if isinstance(pos, int) and pos > 0 and len(primer_seq) - pos >= threshold:
            # Warn if self-dimerization is detected
            logging.warning(f"Self-dimerization detected in {primer_seq} at position {pos}")
            return True
    return False

def design_gibson_primers(vector_seq, insert_seq, start, tm_target=60, overlap_length=20, min_target_length=16, max_target_length=50, max_iterations=100):
    # Log the start of primer design
    logging.info("Designing Gibson Assembly primers")
    start_time = time.time()
    start = int(start) - 1
    
    vector_seq = vector_seq.lower()
    insert_seq = insert_seq.upper()
    
    upstream_overlap = str(vector_seq[start - overlap_length:start])
    downstream_overlap = str(vector_seq[start:start + overlap_length])
    
    insert_left = insert_seq[:max(min_target_length, overlap_length)]
    insert_right = insert_seq[-max(min_target_length, overlap_length):]
    
    forward_primer = upstream_overlap + insert_left
    reverse_primer = str(Seq.Seq(insert_right + downstream_overlap).reverse_complement())
    
    forward_target_part = insert_left
    reverse_target_part = str(Seq.Seq(insert_right).reverse_complement())
    forward_tm_target = mt.Tm_NN(forward_target_part)
    reverse_tm_target = mt.Tm_NN(reverse_target_part)
    forward_tm_full = mt.Tm_NN(forward_primer)
    reverse_tm_full = mt.Tm_NN(reverse_primer)
    
    best_forward = {"primer": forward_primer, "tm_target": forward_tm_target, "tm_full": forward_tm_full, "diff": abs(forward_tm_target - tm_target)}
    best_reverse = {"primer": reverse_primer, "tm_target": reverse_tm_target, "tm_full": reverse_tm_full, "diff": abs(reverse_tm_target - tm_target)}
    
    iteration = 0
    while (forward_tm_target < 55 or abs(forward_tm_target - tm_target) > 2) and iteration < max_iterations:
        if forward_tm_target < 55 and len(insert_left) < max_target_length:
            forward_primer += insert_seq[len(insert_left):len(insert_left) + 1]
            insert_left += insert_seq[len(insert_left):len(insert_left) + 1]
            forward_target_part = insert_left
            forward_tm_target = mt.Tm_NN(forward_target_part)
            forward_tm_full = mt.Tm_NN(forward_primer)
            diff = abs(forward_tm_target - tm_target)
            if diff < best_forward["diff"] and not check_self_dimerization(forward_primer):
                best_forward = {"primer": forward_primer, "tm_target": forward_tm_target, "tm_full": forward_tm_full, "diff": diff}
        elif forward_tm_target > tm_target + 2 and len(insert_left) > min_target_length:
            forward_primer = upstream_overlap + insert_left[:-1]
            insert_left = insert_left[:-1]
            forward_target_part = insert_left
            forward_tm_target = mt.Tm_NN(forward_target_part)
            forward_tm_full = mt.Tm_NN(forward_primer)
            diff = abs(forward_tm_target - tm_target)
            if diff < best_forward["diff"] and not check_self_dimerization(forward_primer):
                best_forward = {"primer": forward_primer, "tm_target": forward_tm_target, "tm_full": forward_tm_full, "diff": diff}
        iteration += 1
    
    iteration = 0
    while (reverse_tm_target < 55 or abs(reverse_tm_target - tm_target) > 2) and iteration < max_iterations:
        if reverse_tm_target < 55 and len(insert_right) < max_target_length:
            insert_right = insert_seq[-len(insert_right) - 1:-len(insert_right)] + insert_right
            reverse_primer = str(Seq.Seq(insert_right + downstream_overlap).reverse_complement())
            reverse_target_part = str(Seq.Seq(insert_right).reverse_complement())
            reverse_tm_target = mt.Tm_NN(reverse_target_part)
            reverse_tm_full = mt.Tm_NN(reverse_primer)
            diff = abs(reverse_tm_target - tm_target)
            if diff < best_reverse["diff"] and not check_self_dimerization(reverse_primer):
                best_reverse = {"primer": reverse_primer, "tm_target": reverse_tm_target, "tm_full": reverse_tm_full, "diff": diff}
        elif reverse_tm_target > tm_target + 2 and len(insert_right) > min_target_length:
            insert_right = insert_right[1:]
            reverse_primer = str(Seq.Seq(insert_right + downstream_overlap).reverse_complement())
            reverse_target_part = str(Seq.Seq(insert_right).reverse_complement())
            reverse_tm_target = mt.Tm_NN(reverse_target_part)
            reverse_tm_full = mt.Tm_NN(reverse_primer)
            diff = abs(reverse_tm_target - tm_target)
            if diff < best_reverse["diff"] and not check_self_dimerization(reverse_primer):
                best_reverse = {"primer": reverse_primer, "tm_target": reverse_tm_target, "tm_full": reverse_tm_full, "diff": diff}
        iteration += 1
    
    if iteration >= max_iterations:
        # Warn if maximum iterations are reached
        logging.warning(f"Max iterations reached. Using best approximations: Forward Tm (Target)={best_forward['tm_target']:.2f}, Reverse Tm (Target)={best_reverse['tm_target']:.2f}")
    forward_primer = best_forward["primer"]
    forward_tm_target = best_forward["tm_target"]
    forward_tm_full = best_forward["tm_full"]
    reverse_primer = best_reverse["primer"]
    reverse_tm_target = best_reverse["tm_target"]
    reverse_tm_full = best_reverse["tm_full"]
    
    # Log the completion of primer design
    logging.info(f"Primers designed in {time.time() - start_time:.2f} seconds")
    logging.info(f"Forward Primer: {forward_primer} (Target Length: {len(insert_left)})")
    logging.info(f"Reverse Primer: {reverse_primer} (Target Length: {len(insert_right)})")
    return {
        "forward_primer": forward_primer,
        "forward_tm_target": forward_tm_target,
        "forward_tm_full": forward_tm_full,
        "reverse_primer": reverse_primer,
        "reverse_tm_target": reverse_tm_target,
        "reverse_tm_full": reverse_tm_full,
        "forward_start": start - overlap_length + 1,  # 1-based
        "forward_end": start + len(insert_left),      # 1-based
        "reverse_start": start + len(insert_seq) - len(insert_right) + 1,  # 1-based
        "reverse_end": start + len(insert_seq) + overlap_length  # 1-based
    }

def adjust_feature_location(line, start, end, length_diff):
    """Adjust feature locations and maintain strict formatting"""
    if not line.strip() or line.strip().startswith("/"):
        return line
    
    parts = line.strip().split(maxsplit=1)
    if len(parts) < 2:
        return line
    
    feature_type = parts[0]
    location = parts[1].split()[0]  # Extract only the location part
    # Log the original feature details
    logging.debug(f"Original feature: {feature_type} at {location}")
    
    start_1based = start
    end_1based = end
    
    if "complement(join(" in location:
        loc_parts = location.split("complement(join(")[1].rstrip("))").split(",")
        adjusted_locs = []
        for loc in loc_parts:
            loc_start, loc_end = map(int, loc.split(".."))
            if loc_end < start_1based:
                adjusted_locs.append(f"{loc_start}..{loc_end}")
            elif loc_start >= end_1based:
                loc_start += length_diff
                loc_end += length_diff
                adjusted_locs.append(f"{loc_start}..{loc_end}")
            else:
                # Log removal of overlapping feature
                logging.debug(f"Removing feature {feature_type} at {loc} overlapping insertion {start_1based}..{end_1based}")
                return None
        if adjusted_locs:
            new_location = f"complement(join({','.join(adjusted_locs)}))"
            return f"     {feature_type:<{FEATURE_WIDTH}} {new_location}\n"
        return None
    
    elif "complement(" in location:
        loc_start, loc_end = map(int, location.split("complement(")[1].rstrip(")").split(".."))
        if loc_end < start_1based:
            pass
        elif loc_start >= end_1based:
            loc_start += length_diff
            loc_end += length_diff
            new_location = f"complement({loc_start}..{loc_end})"
            return f"     {feature_type:<{FEATURE_WIDTH}} {new_location}\n"
        else:
            # Log removal of overlapping feature
            logging.debug(f"Removing feature {feature_type} at {location} overlapping insertion {start_1based}..{end_1based}")
            return None
    
    elif ".." in location:
        loc_start, loc_end = map(int, location.split(".."))
        if loc_end < start_1based:
            pass
        elif loc_start >= end_1based:
            loc_start += length_diff
            loc_end += length_diff
            new_location = f"{loc_start}..{loc_end}"
            return f"     {feature_type:<{FEATURE_WIDTH}} {new_location}\n"
        else:
            # Log removal of overlapping feature
            logging.debug(f"Removing feature {feature_type} at {location} overlapping insertion {start_1based}..{end_1based}")
            return None
    
    elif "^" in location:
        loc_start, loc_end = map(int, location.split("^"))
        if loc_end < start_1based:
            pass
        elif loc_start >= end_1based:
            loc_start += length_diff
            loc_end += length_diff
            new_location = f"{loc_start}^{loc_end}"
            return f"     {feature_type:<{FEATURE_WIDTH}} {new_location}\n"
        else:
            # Log removal of overlapping feature
            logging.debug(f"Removing feature {feature_type} at {location} overlapping insertion {start_1based}..{end_1based}")
            return None
    
    elif location.isdigit():
        loc = int(location)
        if loc < start_1based:
            pass
        elif loc >= end_1based:
            loc += length_diff
            new_location = str(loc)
            return f"     {feature_type:<{FEATURE_WIDTH}} {new_location}\n"
        else:
            # Log removal of overlapping feature
            logging.debug(f"Removing feature {feature_type} at {location} overlapping insertion {start_1based}..{end_1based}")
            return None
    
    return line

def format_long_qualifier(value, first_line_prefix, subsequent_prefix, use_quotes=False, is_target_translation=False):
    """Format qualifier values, preserving newlines and splitting long values"""
    lines = []
    
    if is_target_translation and use_quotes:
        # For target gene's /translation: remove newlines and split
        value = value.replace('\n', '')  # Remove newlines
        remaining_value = value
        first_line_done = False
        
        while remaining_value:
            if not first_line_done:
                prefix = first_line_prefix
                max_width = MAX_LINE_WIDTH - len(first_line_prefix) - 1  # Space for opening quote
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
        # For other qualifiers: preserve newlines
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

def modify_genbank(temp_gb_file, new_seq, start, end, locus_tag, target_seq, product, translation, protein_id, primers, output_path):
    # Log the start of GenBank file modification
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
            parts[1] = f"{parts[1]}_{locus_tag}"
            parts[2] = f"{new_length} bp"
            header_lines[i] = " ".join(parts) + "\n"
        elif line.strip().startswith("REFERENCE") and "(bases" in line:
            parts = line.split("(bases")
            parts[1] = f" (bases 1 to {new_length})\n"
            header_lines[i] = "".join(parts)
    
    length_diff = len(target_seq) - (end - start)
    # Log the length difference for feature adjustment
    logging.debug(f"Length difference for feature adjustment: {length_diff}")
    adjusted_features = []
    current_feature_lines = []
    
    # Adjust existing feature locations and maintain formatting
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
                        use_quotes = False  # No quotes for existing features
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
                use_quotes = False  # No quotes for existing features
                formatted_lines = format_long_qualifier(value, prefix, subsequent_prefix, use_quotes)
                adjusted_features.extend(formatted_lines)
    
    # Add target gene CDS
    target_end = start + len(target_seq) - 1
    adjusted_features.append(f"     CDS             {start}..{target_end + 1}\n")
    qualifiers = [
        ("/locus_tag", locus_tag),
        ("/product", product),
        ("/protein_id", protein_id),
        ("/translation", translation)
    ]
    for i, (key, value) in enumerate(qualifiers):
        prefix = f"{' ' * QUALIFIER_START}{key}="
        subsequent_prefix = " " * QUALIFIER_START
        use_quotes = (key == "/translation")  # Use quotes only for /translation
        is_target_translation = (key == "/translation")  # Indicate target gene's /translation
        formatted_lines = format_long_qualifier(value, prefix, subsequent_prefix, use_quotes, is_target_translation)
        adjusted_features.extend(formatted_lines)
    
    # Add primers (using 'primer' feature) - exclude /primer qualifier
    adjusted_features.append(f"     primer          {primers['forward_start']}..{primers['forward_end']}\n")
    qualifiers = [
        ("/label", "Gibson_Forward"),
        ("/note", f"Tm (Target): {primers['forward_tm_target']:.2f}°C, Tm (Full): {primers['forward_tm_full']:.2f}°C")
    ]
    for i, (key, value) in enumerate(qualifiers):
        prefix = f"{' ' * QUALIFIER_START}{key}="
        subsequent_prefix = " " * QUALIFIER_START
        use_quotes = False
        formatted_lines = format_long_qualifier(value, prefix, subsequent_prefix, use_quotes)
        adjusted_features.extend(formatted_lines)
    
    adjusted_features.append(f"     primer          complement({primers['reverse_start']}..{primers['reverse_end']})\n")
    qualifiers = [
        ("/label", "Gibson_Reverse"),
        ("/note", f"Tm (Target): {primers['reverse_tm_target']:.2f}°C, Tm (Full): {primers['reverse_tm_full']:.2f}°C")
    ]
    for i, (key, value) in enumerate(qualifiers):
        prefix = f"{' ' * QUALIFIER_START}{key}="
        subsequent_prefix = " " * QUALIFIER_START
        use_quotes = False
        formatted_lines = format_long_qualifier(value, prefix, subsequent_prefix, use_quotes)
        adjusted_features.extend(formatted_lines)
    
    # Update ORIGIN section
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
    
    # Log the completion of GenBank file modification
    logging.info(f"GenBank file modified in {time.time() - start_time:.2f} seconds")
    logging.debug(f"Final adjusted features: {adjusted_features}")

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Design Gibson Assembly primers and update DNA sequence in GenBank format with primer features.")
    parser.add_argument("--genbank_table", required=True, help="Path to GenBank table file (TSV or XLSX format)")
    parser.add_argument("--locus_tag", required=True, help="Target locus_tag to clone")
    parser.add_argument("--vector_file", required=True, help="Path to SnapGene DNA vector file (.dna)")
    parser.add_argument("--start", required=True, type=int, help="Start position in vector for insertion (1-based)")
    parser.add_argument("--end", required=True, type=int, help="End position in vector for insertion (1-based)")
    parser.add_argument("--output_dir", default="cloning_results", help="Directory to save output files")
    parser.add_argument("--tm", type=float, default=60, help="Target Tm for primers (default: 60°C, min 55°C)")
    
    args = parser.parse_args()
    
    args.genbank_table = os.path.abspath(args.genbank_table)
    args.vector_file = os.path.abspath(args.vector_file)
    args.output_dir = os.path.abspath(args.output_dir)
    
    # Read GenBank table
    genbank_df = read_genbank_table(args.genbank_table)
    
    # Fetch target sequence and metadata
    target_seq, direction, product, translation, protein_id = get_target_sequence(genbank_df, args.locus_tag)
    if direction == '-':
        target_seq = str(Seq.Seq(target_seq).reverse_complement())
        # Log sequence reversal
        logging.info(f"Target sequence reversed to match 5'->3' direction")
    
    # Read SnapGene vector file
    vector_record = read_snapgene_dna(args.vector_file)
    vector_seq = vector_record.seq
    
    # Replace sequence in vector
    new_seq = replace_sequence(vector_seq, args.start, args.end, target_seq)
    
    # Design Gibson Assembly primers
    primers = design_gibson_primers(vector_seq, target_seq, args.start, tm_target=args.tm)
    
    # Print results
    print("Forward Primer (Backbone 5' + Target 5'):", primers["forward_primer"])
    print("Forward Tm (Target part only):", primers["forward_tm_target"])
    print("Forward Tm (Full):", primers["forward_tm_full"])
    print("Reverse Primer (Target 3' rev_comp + Backbone 3' rev_comp):", primers["reverse_primer"])
    print("Reverse Tm (Target part only):", primers["reverse_tm_target"])
    print("Reverse Tm (Full):", primers["reverse_tm_full"])
    
    # Append to Excel log
    excel_log.append({
        "Forward Primer (Backbone 5' + Target 5')": primers["forward_primer"],
        "Forward Tm (Target part only)": primers["forward_tm_target"],
        "Forward Tm (Full)": primers["forward_tm_full"],
        "Reverse Primer (Target 3' rev_comp + Backbone 3' rev_comp)": primers["reverse_primer"],
        "Reverse Tm (Target part only)": primers["reverse_tm_target"],
        "Reverse Tm (Full)": primers["reverse_tm_full"]
    })
    
    # Convert and modify GenBank file
    temp_gb_file = convert_dna_to_genbank(args.vector_file, args.locus_tag)
    os.makedirs(args.output_dir, exist_ok=True)
    output_path = os.path.join(args.output_dir, f"{args.locus_tag}_vector.gbk")
    modify_genbank(temp_gb_file, new_seq, args.start, args.end, args.locus_tag, target_seq, product, translation, protein_id, primers, output_path)
    
    # Save log to Excel
    log_df = pd.DataFrame(excel_log)
    log_file = os.path.join(args.output_dir, f"log_{args.locus_tag}.xlsx")
    log_df.to_excel(log_file, index=False)
    
    print(f"Updated GenBank file saved to: {output_path}")
    print(f"Log saved to: {log_file}")

if __name__ == "__main__":
    main()