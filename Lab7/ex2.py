import matplotlib.pyplot as plt
import random

def read_fasta(path):
    with open(path) as f:
        return "".join(l.strip().upper() for l in f if not l.startswith(">"))

h1_genome = read_fasta('Helicobacter1.fna')
h2_genome = read_fasta('Helicobacter2.fna')

seq_len = random.randint(1000, 3000)
h1_seq = h1_genome[seq_len : seq_len + 3000]
h2_seq = h2_genome[seq_len : seq_len + 3000]

# h1_seq = h1_genome[seq_len : seq_len + 3000]
# h2_seq = h2_genome[seq_len : seq_len + 3000]

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
    lanes.append(get_fragment_lengths(h1_seq, site))
    lane_names.append(f"H1\n{name}")

for name, site in enzymes.items():
    lanes.append(get_fragment_lengths(h2_seq, site))
    lane_names.append(f"H2\n{name}")

ladder = [100, 200, 300, 400, 500, 750, 1000, 1500, 2000, 3000]

plt.figure(figsize=(9, 7), facecolor='#111111')
ax = plt.gca()
ax.set_facecolor('#111111')

plt.scatter([1] * len(ladder), ladder, marker='_', color='#00FFcc', s=1200, linewidth=3)

for i, frags in enumerate(lanes):
    plt.scatter([i + 2] * len(frags), frags, marker='_', color='#00FFcc', s=1200, linewidth=3)

plt.yscale('log')
plt.ylim(50, 4000)
plt.gca().invert_yaxis()

plt.xticks(range(1, 8), lane_names, color='white', fontsize=11)
plt.yticks(color='white')
plt.ylabel("Fragment Size (bp)", color='white')
plt.title("Comparative Gel: Helicobacter1 vs. Helicobacter2", color='white', pad=20)

for spine in ax.spines.values():
    spine.set_visible(False)

plt.tight_layout()
plt.show()