import random
import matplotlib.pyplot as plt
import numpy as np

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

reads = []
for _ in range(N_READS):
    L = random.randint(MIN_LEN, G)
    start = random.randint(0, G - L)
    reads.append(sequence[start:start + L])
    
for i in range(len(reads)):
    print(len(reads[i]))
    
sizes = [len(r) for r in reads]

def distance(size):
    return 1 / np.log10(size)
    
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_facecolor("black")

for s in sizes:
    d = distance(s)
    ax.hlines(d, 0.3, 0.7, colors="green", linewidth=4)
    ax.text(0.75, d, f"{s} bp", color="white", va="center", fontsize=8)

ax.set_xlim(0, 1.2)
ax.invert_yaxis()
ax.set_xticks([])
ax.set_ylabel("Migration distance")
ax.set_title("Simulated gel")
plt.show()