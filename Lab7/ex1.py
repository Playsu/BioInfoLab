import random
import matplotlib.pyplot as plt

fasta_file = 'GCF_000857045.1_ViralProj15142_genomic.fna'

def read_fasta(path):
    with open(path) as f:
        return "".join(l.strip().upper() for l in f if not l.startswith(">"))
    
genome = read_fasta(fasta_file)

seq_len = random.randint(1000, 3000)
seq_start = random.randint(0, len(genome) - seq_len)
sequence = genome[seq_start:seq_start + seq_len]


enzymes = {
    'EcoRI': 'GAATTC', 
    'BamHI': 'GGATCC', 
    'HindIII': 'AAGCTT'
}

def get_fragment_lengths(seq, site):
    cuts = []
    idx = seq.find(site)
    while idx != -1:
        cuts.append(idx + 1)
        idx = seq.find(site, idx + 1)
        
    fragments = []
    prev = 0
    for cut in cuts:
        fragments.append(cut - prev)
        prev = cut
    fragments.append(len(seq) - prev)
    
    return fragments

lanes = []
lane_names = ['Ladder']
for name, site in enzymes.items():
    lanes.append(get_fragment_lengths(sequence, site))
    lane_names.append(name)

ladder = [100, 200, 300, 400, 500, 750, 1000, 1500, 2000, 3000]

plt.figure(figsize=(6, 8), facecolor='#111111')
ax = plt.gca()
ax.set_facecolor('#111111')

plt.scatter([1] * len(ladder), ladder, marker='_', color='#00FFcc', s=2000, linewidth=3)

for i, frags in enumerate(lanes):
    plt.scatter([i + 2] * len(frags), frags, marker='_', color='#00FFcc', s=2000, linewidth=3)

plt.yscale('log')
plt.ylim(50, 4000)
plt.gca().invert_yaxis()

plt.xticks([1, 2, 3, 4], lane_names, color='white', fontsize=12)
plt.yticks(color='white')
plt.ylabel("Fragment Size (bp)", color='white')
plt.title("Simulated Gel Electrophoresis", color='white', pad=15)

for spine in ax.spines.values():
    spine.set_visible(False)

plt.show()