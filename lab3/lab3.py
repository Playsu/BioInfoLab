import os
import matplotlib.pyplot as plt
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq

from dotenv import load_dotenv
from openai import OpenAI

def codon_to_amino_acid(codon):
    genetic_code = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'Stop', 'TAG':'Stop', 'TGC':'C', 'TGT':'C', 'TGA':'Stop', 'TGG':'W',
    }
    return genetic_code.get(codon.upper(), "Unknown")

def load_local_genome(filepath):
    if not os.path.exists(filepath):
        print(f"[!] Error: File not found at '{filepath}'")
        return None
    print(f"Loading sequence from: {filepath}...")
    try:
        record = SeqIO.read(filepath, "fasta")
        print("Genome loaded successfully!")
        return record
    except Exception as e:
        print(f"[!] Failed to read FASTA file: {e}")
        return None

def analyze_codon_frequencies(sequence):
    seq_str = str(sequence).upper()
    rev_seq_str = str(sequence.reverse_complement()).upper()
    
    frames = [
        seq_str[0:], seq_str[1:], seq_str[2:],
        rev_seq_str[0:], rev_seq_str[1:], rev_seq_str[2:]
    ]
    
    all_codons = []
    
    for frame in frames:
        codons = [frame[i:i+3] for i in range(0, len(frame) - len(frame)%3, 3)]
        valid_codons = [c for c in codons if len(c) == 3 and 'N' not in c]
        all_codons.extend(valid_codons)
        
    codon_counts = Counter(all_codons)
    return codon_counts.most_common()

def plot_top_10_codons(sorted_codons):
    top_10 = sorted_codons[:10]
    codons = [item[0] for item in top_10]
    counts = [item[1] for item in top_10]
    
    plt.figure(figsize=(10, 6))
    plt.bar(codons, counts, color='teal', edgecolor='black')
    plt.title('Top 10 Most Frequent Codons (All 6 Reading Frames)', fontsize=14)
    plt.xlabel('Codons', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    
    for i, count in enumerate(counts):
        plt.text(i, count + (max(counts) * 0.01), str(count), ha='center')
        
    plt.tight_layout()
    plt.show()
    
    return top_10

def generate_llm_prompt(top_codons):
    amino_acids = set()
    for codon, _ in top_codons:
        aa = codon_to_amino_acid(codon)
        if aa not in ['Stop', 'Unknown']:
            amino_acids.add(aa)
            
    aa_list_str = ", ".join(sorted(list(amino_acids)))
    
    prompt = (
        f"I have analyzed the genome of Helicobacter pylori strain CHC155. "
        f"The most frequent codons code for the following amino acids (1-letter codes): {aa_list_str}. "
        f"Please provide a comprehensive list of foods that contain the LOWEST amounts of these specific amino acids. "
        f"Format the output as a clear table."
    )
    return prompt

def ask_openai(prompt_text):
    """Securely calls the OpenAI API using the hidden .env file."""
    load_dotenv()
    api_key = os.getenv("OPENAI_API_KEY")
    
    if not api_key:
        print("\n[!] ERROR: API key not found.")
        print("Make sure you created a .env file with OPENAI_API_KEY=your_key_here")
        return
        
    client = OpenAI(api_key=api_key)
    
    print("\nSending prompt to OpenAI API...")
    try:
        response = client.chat.completions.create(
            model="gpt-4o-mini",
            messages=[
                {"role": "system", "content": "You are a helpful bioinformatics and nutrition assistant."},
                {"role": "user", "content": prompt_text}
            ],
            temperature=0.7
        )
        
        print("\n=== AI RESPONSE ===")
        print(response.choices[0].message.content)
        
    except Exception as e:
        print(f"\n[!] API Call Failed: {e}")

if __name__ == "__main__":
    
    current_dir = os.path.dirname(os.path.abspath(__file__))
    local_file_path = os.path.join(current_dir, "data", "GCF_025998455.1_ASM2599845v1_genomic.fna")    
    
    genome_record = load_local_genome(local_file_path)
    
    if genome_record:
        print(f"\nAnalyzing Sequence: {genome_record.description}")
        print(f"Total Sequence Length: {len(genome_record.seq)} bp")
        
        sorted_frequencies = analyze_codon_frequencies(genome_record.seq)
        
        print("\nGenerating bar chart...")
        top_10_codons = plot_top_10_codons(sorted_frequencies)
        
        final_prompt = generate_llm_prompt(top_10_codons)
        
        ask_openai(final_prompt)