import random

fasta_file = 'GCF_000006945.2_ASM694v2_genomic.fna'
N_READS = 10
MIN_LEN = 100
MAX_LEN = 3000

def read_fasta(path):
    with open(path) as f:
        return "".join(l.strip().upper() for l in f if not l.startswith(">"))
    
genome = read_fasta(fasta_file)

seq_len = random.randint(1000, 3000)
seq_start = random.randint(0, len(genome) - seq_len)
sequence = genome[seq_start:seq_start + seq_len]
G = len(sequence)

counts = {}
for k in range(3, 7):
    for i in range(len(sequence) - k + 1):
        piece = sequence[i:i + k]
        counts[piece] = counts.get(piece, 0) + 1
        
for piece, n in counts.items():
    if 6 <= n <= 10:
        print(piece, n)