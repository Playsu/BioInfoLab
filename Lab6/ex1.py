import re

def analyze_transposons(dna_seq, tir_left="ATGCA", tir_right="TACGT"):
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
            normal_tes.append({"start": s_pos, "end": e_pos, "seq": dna_seq[s_pos:e_pos+1]})
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
            nested_tes.append({"start": s_pos, "end": pos, "seq": dna_seq[s_pos:pos+1]})

    queue = []
    overlapping_tes = []
    for pos, etype in complex_events:
        if etype == 'start':
            queue.append(pos)
        elif etype == 'end' and queue:
            s_pos = queue.pop(0)
            overlapping_tes.append({"start": s_pos, "end": pos, "seq": dna_seq[s_pos:pos+1]})

    return normal_tes, nested_tes, overlapping_tes


dna_sequence = (
    "GGATCCGTAAGCGCGTA"                           
    "ATGCAGGGGGGGGGGGGGGGGTACGT"                  
    "CTACTACTACTACTACTACTACTACTACTA"
    "ATGCATTTTTTACGT"
    "AGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG"
    "ATGCACTGACTGACTGACTGACTGACTGACTGACTACGT"
    "CCCGGGAAATTT"
    "ATCGATGTCATTCATGCATATCGTTACGGTTTACGT"
    "ATGCACCCATGCAGGGGTACGTTTTTACGT"
) 

normal, nested, overlapping = analyze_transposons(dna_sequence)

print("NORMAL CASES")
for i, te in enumerate(normal, 1):
    print(f"TE {i} [{te['start']}:{te['end']}]: {te['seq']}")

print("\nNESTED INTERPRETATION")
for i, te in enumerate(nested, 1):
    print(f"Nested TE {i} [{te['start']}:{te['end']}]: {te['seq']}")

print("\nOVERLAPPING INTERPRETATION")
for i, te in enumerate(overlapping, 1):
    print(f"Overlapping TE {i} [{te['start']}:{te['end']}]: {te['seq']}")