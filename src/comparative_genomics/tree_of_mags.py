import argparse
import os
import subprocess
import time
from pathlib import Path
import importlib.resources as pkg_resources
from multiprocessing import cpu_count

from comparative_genomics.fasta import FastaParser, write_fasta
from comparative_genomics.blast import TabularBlastParser
from comparative_genomics.orthologues import compute_orthologues, write_orthologues_to_fasta
import comparative_genomics.database


VERSION = "0.1"
START_TIME = time.monotonic()
LOG_FILE = Path('log.txt')


def log(log_message, values=()):
    if len(values):
        final_msg = f'{format_runtime()} {log_message.format(*values)}'
    else:
        final_msg = f'{format_runtime()} {log_message}'
    print(final_msg)
    try:
        with open(LOG_FILE, 'a') as log_handle:
            log_handle.write(final_msg)
            log_handle.write('\n')
    except FileNotFoundError:
        pass


def format_runtime():
    runtime = time.monotonic() - START_TIME
    return f'[{int(runtime / 3600):02d}h:{int((runtime % 3600) / 60):02d}m:{int(runtime % 60):02d}s]'


def run_external(exec, stdin=None, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, log_cmd=True):
    if log_cmd:
        log(exec)
    result = subprocess.run(exec.split(), stdout=stdout, stdin=stdin, stderr=stderr)
    if result.returncode != 0:
        raise Exception(f'Error while trying to run "{exec}"')


def parse_arguments():
    parser = argparse.ArgumentParser(description='comparative_genomics.py. (C) Marc Strous, 2022')
    parser.add_argument('--dir', help='Folder with aminoacid fasta files, one for each genome or mag.')
    parser.add_argument('--cpus', default=cpu_count(), help='How many cpus/threads to use (default: all = 0).')
    parser.add_argument('--file_extension', default='.faa', help='extension of aminoacid fasta files (default ".faa")')
    parser.add_argument('--delimiter', default='|', help='character to separate filenames and orfnames during '
                                                         'orthologue calling')
    parser.add_argument('--min_frequency', default=0, help='The minimum fraction of genes a taxon should have to be'
                                                             ' included in the final multiple sequence alignment '
                                                             '(default 0)')
    return parser.parse_args()


def prep_hmms(hmm_dir):
    log('extracting hmms included in python package...')
    hmm_file = hmm_dir / 'conserved_genes.hmm'
    with open(hmm_file, 'w') as handle:
        for hmms in ('gtdb-pfam.hmm', 'gtdb-tigr.hmm', 'ribosomal.pfam.hmm', 'ribosomal.tigr.hmm', 'rpoABC.hmm'):
            handle.write(pkg_resources.read_text(comparative_genomics.database, hmms))
    run_external(f'hmmpress -f {hmm_file}')
    return hmm_file


def collect_seqs(hmm_file: Path, fasta_dir: Path, genes_dir: Path, file_extension, delimiter: str = ''):
    log('Now collecting target sequences from fasta files using hmmscan...')
    all_orfs = {}
    hmm_lengths = {}
    with open(hmm_file) as handle:
        for line in handle:
            if line.startswith('NAME'):
                name = line.split()[1]
            elif line.startswith('LENG'):
                hmm_lengths[name] = int(line[4:].strip())
    log(f'HMM database has {len(hmm_lengths)} HMM profiles.')

    for file in sorted(fasta_dir.glob(f'*{file_extension}')):
        fasta_file = fasta_dir / file.name
        hmm_results_file = genes_dir / f'{file.name}.hmm-results'
        gene_file = genes_dir / file.name
        unique_orf_ids = set()
        with FastaParser(fasta_file) as fasta_reader:
            for orf in fasta_reader:
                all_orfs[orf.id] = orf
                if delimiter and delimiter in orf['id']:
                    raise Exception('Delimiter should not occur in ids of any fasta seq.')
                if orf['id'] in unique_orf_ids:
                    raise Exception(f'Fasta seq ids should be unique within each file. Duplicate: {orf["id"]}')
                unique_orf_ids.add(orf['id'])
        run_external(f'hmmscan -E 1e-25 --domtblout {hmm_results_file} {hmm_file} {fasta_file}')
        orfs_already_done = set()
        with TabularBlastParser(hmm_results_file, 'HMMSCAN') as handle:
            count = 0
            with open(gene_file, 'w') as writer:
                for blast_result in handle:
                    for hit in blast_result:
                        if hit.aligned_length / hmm_lengths[hit] < 0.7:
                            continue
                        if hit.query in orfs_already_done:
                            continue
                        orfs_already_done.add(hit.query)
                        orf = all_orfs[hit.query]
                        write_fasta(writer, orf)
                        count += 1
                        break
        log(f'{file}: Wrote {count} seqs to fasta.')


def align_seqs(src_dir: Path, dest_dir: Path, file_extension):
    log(f'Now aligning fasta files in {src_dir}...')
    for file in sorted(src_dir.glob(f'*{file_extension}')):
        raw_result_file = dest_dir / f'{file.name}.raw'
        final_result_file = dest_dir / file.name
        run_external('clustalo ...')
        run_external('prune alignments ...')


def concatenate_alignments(src_dir: Path, file_extension: str, delimiter: str, min_frequency = 0):
    log(f'Now aligning fasta files in {src_dir}...')
    alignments = []
    lengths = []
    unique_taxa = set()
    for file in sorted(src_dir.glob(f'*{file_extension}')):
        file = src_dir / file
        alignment = {}
        length = 0
        with FastaParser(file) as fasta_reader:
            for orf in fasta_reader:
                taxon, orf_id = orf['id'].split(delimiter)
                alignment[taxon] = orf
                unique_taxa.add(taxon)
                if length:
                    if len(orf['seq']) != length:
                        raise Exception(f'Unequal alignment length in file {file}')
                else:
                    length = len(orf['seq'])
            alignments.append(alignment)
            lengths.append(length)
    occurences = {}
    sorted_taxa = sorted([unique_taxa])
    concatenated_alignment = {}
    for taxon in sorted_taxa:
        concatenated_fasta_seq = {'id': taxon, 'seq': ''}
        for i in range(len(alignments)):
            alignment = alignments[i]
            if taxon in alignment.keys():
                concatenated_fasta_seq['seq'] += alignment[taxon]['seq']
                prev_occurence = occurences.get(taxon, 0)
                occurences[taxon] = prev_occurence + 1
            else:
                concatenated_fasta_seq['seq'] += '-' * lengths[i]
        concatenated_alignment[taxon] = concatenated_fasta_seq
    with open(src_dir / 'concatenated_alignment', 'w') as writer:
        for taxon in sorted_taxa:
            frequency = occurences[taxon]/len(alignments)
            if occurences[taxon] >= 1 and frequency >= min_frequency:
                write_fasta(writer, concatenated_alignment[taxon])
                log(f'[taxon] Written, frequency {frequency:.2f}')
            else:
                log(f'[taxon] Skipped, frequency {frequency:.2f} < minimum ({min_frequency:.2f})')


def main():
    print(f'This is comparative_genomics.py {VERSION}')
    args = parse_arguments()
    fasta_dir = Path(args.dir)
    cpus = int(args.cpus)
    file_extension = args.file_extension
    delimiter = args.delimiter
    min_frequency = args.min_frequency

    working_dir = Path(os.getcwd())
    hmm_dir = working_dir/ 'hmm'
    genes_dir = working_dir/ 'genes'
    orthologue_dir = working_dir/ 'orthologues'
    alignments_dir = working_dir/ 'alignments'
    temp_dir = working_dir/ 'temp'

    for d in (hmm_dir, genes_dir, orthologue_dir, alignments_dir, temp_dir):
        d.mkdir(exist_ok=True)

    hmm_file = prep_hmms(hmm_dir)
    collect_seqs(hmm_file, fasta_dir, genes_dir, file_extension, delimiter)
    merged_and_coded_fasta_file, taxa_by_orf_id, unique_blast_results, orthologues, orthologues_by_orf_id = \
        compute_orthologues(genes_dir, cpus, file_extension, delimiter)
    write_orthologues_to_fasta(merged_and_coded_fasta_file, orthologues_by_orf_id, taxa_by_orf_id, orthologue_dir,
                               include_paralogues=False)
    align_seqs(orthologue_dir, alignments_dir, file_extension)
    concatenate_alignments(alignments_dir, file_extension, delimiter, min_frequency)

if __name__ == "__main__":
    main()