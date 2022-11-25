import argparse
import os
import subprocess
import time
from pathlib import Path
from multiprocessing import cpu_count
from concurrent import futures

from comparative_genomics.fasta import FastaParser, write_fasta
from comparative_genomics.blast import TabularBlastParser
from comparative_genomics.orthologues import compute_orthologues, write_orthologues_to_fasta, SetOfOrthologues


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
    parser.add_argument('--dir', required=True, help='Folder with aminoacid fasta files, one for each genome or mag.')
    parser.add_argument('--cpus', default=cpu_count(), help='How many cpus/threads to use (default: all = 0).')
    parser.add_argument('--file_extension', default='.faa', help='extension of aminoacid fasta files (default ".faa")')
    parser.add_argument('--delimiter', default='|', help='character to separate filenames and orfnames during '
                                                         'orthologue calling')
    parser.add_argument('--predict_orfs_to_dir', default='', help='predict orfs and store in dir')
    parser.add_argument('--min_frequency', default=0.6, help='The minimum fraction of genes a taxon should have to be'
                                                             ' included in the final multiple sequence alignment '
                                                             '(default 0)')
    return parser.parse_args()


def prep_hmms(hmm_dir):
    log('extracting hmms included in python package...')
    hmm_file = hmm_dir / 'conserved_genes.hmm'
    with open(hmm_file, 'w') as writer:
        for hmms in (Path(__file__).parent / 'database').glob('*.hmm'):
            if hmms.name.startswith('ribosomal') or hmms.name.startswith('rpo'):
                continue
            print(hmms)
            run_external(f'hmmconvert {hmms}', stdout=writer)
    run_external(f'hmmpress -f {hmm_file}')
    return hmm_file


def predict_orfs(nt_dir: Path, aa_dir: Path, file_extension: str, cpus: int):
    aa_dir.mkdir(exist_ok=True, parents=True)
    future_obj = []
    with futures.ProcessPoolExecutor(max_workers=cpus) as executor:
        for nt_file in nt_dir.glob('*' + file_extension):
            if not (aa_dir / nt_file.name).exists():
                future_obj.append(
                    executor.submit(run_external, f'prodigal -m -f gff -q -i {nt_file} -a {aa_dir / nt_file.name}'))
        for f in future_obj:
            f.result()


def collect_seqs(hmm_file: Path, fasta_dir: Path, genes_dir: Path, file_extension, delimiter: str = '', cpus: int = 1):
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

    with futures.ProcessPoolExecutor(max_workers=cpus) as executor:
        future_obj = {}
        for file in sorted(fasta_dir.glob(f'*{file_extension}')):
            fasta_file = fasta_dir / file.name
            hmm_results_file = genes_dir / f'{file.name}.hmm-results'
            gene_file = genes_dir / file.name
            unique_orf_ids = set()
            with FastaParser(fasta_file) as fasta_reader:
                for orf in fasta_reader:
                    all_orfs[orf['id']] = orf
                    if delimiter and delimiter in orf['id']:
                        raise Exception('Delimiter should not occur in ids of any fasta seq.')
                    if orf['id'] in unique_orf_ids:
                        raise Exception(f'Fasta seq ids should be unique within each file. Duplicate: {orf["id"]}')
                    unique_orf_ids.add(orf['id'])
                if hmm_results_file.exists() and hmm_results_file.stat().st_size:
                    future_obj[hmm_results_file] = 1
                else:
                    future_obj[hmm_results_file] = executor.submit(run_external,
                                        f'hmmscan -E 1e-25 --domtblout {hmm_results_file} {hmm_file} {fasta_file}')

        for f in future_obj.keys():
            gene_file = genes_dir / (f.name + file_extension)
            if f.exists() and f.stat().st_size:
                if not gene_file.exists() or not gene_file.stat().st_size:
                    orfs_already_done = set()
                    with TabularBlastParser(f, 'HMMSCAN_DOM_TABLE') as handle:
                        count = 0
                        with open(gene_file, 'w') as writer:
                            for blast_result in handle:
                                for hit in blast_result:
                                    if hit.aligned_length / hmm_lengths[hit.hit] < 0.7:
                                        continue
                                    if hit.query in orfs_already_done:
                                        continue
                                    orfs_already_done.add(hit.query)
                                    orf = all_orfs[hit.query]
                                    write_fasta(writer, orf)
                                    count += 1
                                    break
                        log(f'{gene_file}: Wrote {count} seqs to fasta.')


def filter_orthologues(taxa_by_orf_id: list, orthologues: list[SetOfOrthologues], orthologues_by_orf_id: dict,
                       min_frequency: float, fasta_dir: Path, file_extension: str):
    unique_taxa = set(taxa_by_orf_id)
    log(f'Now filtering {len(orthologues)} orthologues for representation of and among {len(unique_taxa)} taxa.')
    remaining_orthologues = {}
    removed_orthologues = set()
    for i in range(len(orthologues)):
        o = orthologues[i]
        if len(o.orthologues) >= min_frequency * len(unique_taxa):
            remaining_orthologues[i] = o
        else:
            removed_orthologues.add(i)
    log(f'Keeping {len(remaining_orthologues)}/{len(orthologues)} orthologues with >= {int(min_frequency * len(unique_taxa))} taxa.')

    remaining_taxa = []
    removed_taxa = set()
    taxa_names = sorted(fasta_dir.glob(f'*{file_extension}'))
    for taxon in unique_taxa:
        freq = 0
        for o in remaining_orthologues.values():
            for orf_id in o.orthologues:
                if taxa_by_orf_id[orf_id] == taxon:
                    freq += 1
                    break
        if freq >= min_frequency * len(remaining_orthologues):
            remaining_taxa.append(taxon)
        else:
            removed_taxa.add(taxon)
            log(f'{taxa_names[taxon]} just got removed from analysis.')
    log(f'Keeping {len(remaining_taxa)}/{len(unique_taxa)} taxa with >= {int(min_frequency * len(remaining_orthologues))} remaining orthologues.')

    # Now sort through alignment source files and weed out orthologues and taxa with poor representation
    prev_orf_count = len(orthologues_by_orf_id)
    for orf_id in range(len(taxa_by_orf_id)):
        if orthologue_id := orthologues_by_orf_id.get(orf_id, 0):
            if int(orthologue_id[1:]) in removed_orthologues:
                orthologues_by_orf_id.pop(orf_id)
            elif taxa_by_orf_id[orf_id] in removed_taxa:
                orthologues_by_orf_id.pop(orf_id)
    log(f'{len(orthologues_by_orf_id)}/{prev_orf_count} total proteins remaining in analysis.')



def run_align_programs(src_file: Path, raw_result_file: Path, final_result_file: Path, cpus: int):
    run_external(f'clustalo -i {src_file} -o {raw_result_file} -t Protein --threads={cpus}')
    run_external(f'java -jar /bio/bin/BMGE/src/BMGE.jar -i {raw_result_file} -t AA -o {final_result_file}')


def align_seqs(src_dir: Path, dest_dir: Path, file_extension, cpus: int):
    log(f'Now aligning fasta files in {src_dir}...')
    with futures.ProcessPoolExecutor(max_workers=cpus) as executor:
        future_obj = {}
        for src_file in sorted(src_dir.glob(f'*{file_extension}')):
            raw_result_file = dest_dir / f'{src_file.name}.raw'
            final_result_file = dest_dir / src_file.name
            if not final_result_file.exists() or not final_result_file.stat().st_size:
                run_align_programs(src_file, raw_result_file, final_result_file, cpus)
                # future_obj[src_file] = executor.submit(run_align_programs, src_file, raw_result_file, final_result_file)
        for ff in future_obj.values():
            ff.result()


def concatenate_alignments(src_dir: Path, file_extension: str, delimiter: str, min_frequency = 0):
    log(f'Now concatenating alignments in {src_dir}...')
    alignments = []
    unique_taxa = set()
    file_extension = '.faa.raw'
    for file in sorted(src_dir.glob(f'*{file_extension}')):
        file = src_dir / file
        alignment = {}
        length = 0
        with FastaParser(file, cleanup_seq=False) as fasta_reader:
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
    taxa_removed = set()
    log(f'Parsed {len(alignments)} alignments from fasta files with {len(unique_taxa)} taxa.')

    concatenated_alignment = {}
    concatenated_alignment_length = 0
    for taxon in unique_taxa:
        concatenated_fasta_seq = {'id': taxon, 'seq': ''}
        count = 0
        for alignment in alignments:
            length = 0
            for t, s in alignment.items():
                length = len(s['seq'])
                break
            try:
                concatenated_fasta_seq['seq'] += alignment[taxon]['seq']
                count += 1
            except KeyError:
                concatenated_fasta_seq['seq'] += '-' * length
        if not concatenated_alignment:
            concatenated_alignment_length = len(concatenated_fasta_seq["seq"])
            log(f'Concatenated alignment has {concatenated_alignment_length} positions.')
        concatenated_alignment[taxon] = concatenated_fasta_seq
        #log(f'{count:4d} {taxon}')
    success_count = 0
    for pos in range(concatenated_alignment_length):
        count = 0
        for taxon in unique_taxa:
            seq = concatenated_alignment[taxon]['seq']
            if seq[pos].isalnum():
                count += 1
        if count >= min_frequency * len(unique_taxa):
            success_count += 1
    log(f'{success_count}/{concatenated_alignment_length} positions ({success_count/concatenated_alignment_length:.1%}) '
        f'have < {1-min_frequency:.1%} gaps.')

    seqs_done = set()
    with open(src_dir / 'concatenated_alignment', 'w') as writer:
        for taxon in unique_taxa:
            a = concatenated_alignment[taxon]
            if a['seq'] in seqs_done:
                log(f'skipping {taxon} - duplicated sequence.')
            else:
                seqs_done.add(a['seq'])
            write_fasta(writer, a)


def main():
    print(f'This is comparative_genomics.py {VERSION}')
    args = parse_arguments()
    cpus = int(args.cpus)
    fasta_dir = Path(args.dir)
    fasta_aa_dir = Path(args.predict_orfs_to_dir)
    file_extension = args.file_extension
    delimiter = args.delimiter
    min_frequency = args.min_frequency
    os.environ["PATH"] += ':/bio/bin:/bio/bin/hmmer3/bin'

    if args.predict_orfs_to_dir:
        predict_orfs(fasta_dir, fasta_aa_dir, file_extension, cpus)
        fasta_dir = fasta_aa_dir
    working_dir = Path(os.getcwd())
    hmm_dir = working_dir/ 'hmm'
    genes_dir = working_dir/ 'genes'
    orthologue_dir = working_dir/ 'orthologues'
    alignments_dir = working_dir/ 'alignments'
    temp_dir = working_dir/ 'temp'

    for d in (hmm_dir, genes_dir, orthologue_dir, alignments_dir, temp_dir):
        d.mkdir(exist_ok=True)

    hmm_file = prep_hmms(hmm_dir)
    collect_seqs(hmm_file, fasta_dir, genes_dir, file_extension, delimiter, cpus)
    merged_and_coded_fasta_file, taxa_by_orf_id, unique_blast_results, orthologues, orthologues_by_orf_id = \
        compute_orthologues(genes_dir, cpus, file_extension, delimiter)
    filter_orthologues(taxa_by_orf_id, orthologues, orthologues_by_orf_id, min_frequency, genes_dir, file_extension)
    write_orthologues_to_fasta(merged_and_coded_fasta_file, orthologues_by_orf_id, taxa_by_orf_id, orthologue_dir,
                               include_paralogues=False)
    align_seqs(orthologue_dir, alignments_dir, '.faa', cpus)
    # number_of_taxa = set(taxa_by_orf_id.values())
    concatenate_alignments(alignments_dir, '.faa', delimiter, min_frequency)

if __name__ == "__main__":
    main()