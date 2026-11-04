import re

def read_fasta(path):
    """Reads a FASTA file and returns a single continuous DNA string."""
    with open(path) as f:
        return "".join(l.strip().upper() for l in f if not l.startswith(">"))

def analyze_transposons(dna_seq, tir_left="ATGCA", tir_right="TACGT", min_len=10):
    """Detects Normal, Nested, and Overlapping TEs, filtering by minimum length."""
    starts = [m.start() for m in re.finditer(f"(?={tir_left})", dna_seq)]
    ends = [m.start() + len(tir_right) - 1 for m in re.finditer(f"(?={tir_right})", dna_seq)]
    
    events = [(s, 'start') for s in starts] + [(e, 'end') for e in ends]
    events.sort(key=lambda x: x[0])
    
    normal_tes = []
    complex_events = []
    
    i = 0
    while i < len(events):
        if i < len(events) - 1 and events[i][1] == 'start' and events[i+1][1] == 'end':
            s_pos = events[i][0]
            e_pos = events[i+1][0]
            te_len = e_pos - s_pos + 1
            if te_len >= min_len:
                normal_tes.append({
                    "start": s_pos, "end": e_pos, 
                    "length": te_len, "seq": dna_seq[s_pos:e_pos+1]
                })
            i += 2 
        else:
            complex_events.append(events[i])
            i += 1

    stack = []
    nested_tes = []
    for pos, etype in complex_events:
        if etype == 'start':
            stack.append(pos)
        elif etype == 'end' and stack:
            s_pos = stack.pop()
            te_len = pos - s_pos + 1
            if te_len >= min_len:
                nested_tes.append({
                    "start": s_pos, "end": pos, 
                    "length": te_len, "seq": dna_seq[s_pos:pos+1]
                })

    queue = []
    overlapping_tes = []
    for pos, etype in complex_events:
        if etype == 'start':
            queue.append(pos)
        elif etype == 'end' and queue:
            s_pos = queue.pop(0)
            te_len = pos - s_pos + 1
            if te_len >= min_len:
                overlapping_tes.append({
                    "start": s_pos, "end": pos, 
                    "length": te_len, "seq": dna_seq[s_pos:pos+1]
                })

    return normal_tes, nested_tes, overlapping_tes


fasta_file = 'GCF_000857045.1_ViralProj15142_genomic.fna'

try:
    genome = read_fasta(fasta_file)
    print(f"Loaded genome with {len(genome):,} base pairs.")
    
    normal, nested, overlapping = analyze_transposons(genome, min_len=10)

    print(f"NORMAL CASES")
    print(f"Found {len(normal)} normal transposons.")
    for i, te in enumerate(normal[:5], 1): # Printing first 5 for brevity
        print(f"TE {i} [{te['start']}:{te['end']}] (Len: {te['length']}): {te['seq']}")

    print(f"\nNESTED INTERPRETATION")
    print(f"Found {len(nested)} nested transposons.")
    for i, te in enumerate(nested[:5], 1):
        print(f"Nested TE {i} [{te['start']}:{te['end']}] (Len: {te['length']}): {te['seq']}")

    print(f"\nOVERLAPPING INTERPRETATION")
    print(f"Found {len(overlapping)} overlapping transposons.")
    for i, te in enumerate(overlapping[:5], 1):
        print(f"Overlapping TE {i} [{te['start']}:{te['end']}] (Len: {te['length']}): {te['seq']}")

except FileNotFoundError:
    print(f"Error: The file {fasta_file} was not found in the current directory.")