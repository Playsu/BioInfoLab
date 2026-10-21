import random
from collections import defaultdict

fasta_file = 'GCF_000857045.1_ViralProj15142_genomic.fna'
N_READS, MIN_LEN, MAX_LEN = 2000, 100, 120
MIN_OVERLAP, K = 20, 12
random.seed(42)


def read_fasta(path):
    with open(path) as f:
        return "".join(l.strip().upper() for l in f if not l.startswith(">"))


def wrap(seq, width=80):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def assemble(reads, min_ov=MIN_OVERLAP):
    reads = sorted(set(reads), key=len, reverse=True)
    reads = [r for i, r in enumerate(reads)
             if not any(r in o for o in reads[:i])]
    n = len(reads)

    prefix_idx = defaultdict(list)
    for j, r in enumerate(reads):
        prefix_idx[r[:K]].append(j)

    edges = []
    for i, a in enumerate(reads):
        for pos in range(len(a) - min_ov + 1):
            for j in prefix_idx.get(a[pos:pos + K], []):
                if j != i and a[pos:] == reads[j][:len(a) - pos]:
                    edges.append((len(a) - pos, i, j))
    edges.sort(reverse=True)

    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    nxt, prv = {}, set()
    for ov, i, j in edges:
        if i in nxt or j in prv or find(i) == find(j):
            continue
        nxt[i] = (j, ov)
        prv.add(j)
        parent[find(i)] = find(j)

    contigs = []
    for s in range(n):
        if s in prv:
            continue
        contig, cur = reads[s], s
        while cur in nxt:
            cur, ov = nxt[cur]
            contig += reads[cur][ov:]
        contigs.append(contig)
    return sorted(contigs, key=len, reverse=True)


genome = read_fasta(fasta_file)
G = len(genome)

reads = []
for _ in range(N_READS):
    L = random.randint(MIN_LEN, MAX_LEN)
    start = random.randint(0, G - L)
    reads.append(genome[start:start + L])

contigs = assemble(reads)

with open("rebuilt_contigs.fasta", "w") as out:
    for idx, c in enumerate(contigs, 1):
        record = f">contig_{idx} length={len(c)}\n{wrap(c)}\n"
        out.write(record)
        print(record)

if len(contigs) == 1 and contigs[0] == genome:
    print("Perfect reconstruction.")
else:
    print(f"Partial reconstruction: {len(contigs)} contig(s), "
          f"{sum(map(len, contigs)):,} of {G:,} bp.")