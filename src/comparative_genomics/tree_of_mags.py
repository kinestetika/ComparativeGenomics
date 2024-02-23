import argparse
import os
import subprocess
import time
import shutil
from pathlib import Path

from comparative_genomics.fasta import FastaParser, write_fasta
from comparative_genomics.blast import TabularBlastParser
from comparative_genomics.orthologues import merge_and_code_fasta_input, cluster, align, concatenate_alignments


VERSION = "0.17"
START_TIME = time.monotonic()
LOG_FILE = Path('log.txt')


def parse_arguments():
    parser = argparse.ArgumentParser(description='comparative_genomics.py. (C) Marc Strous, 2024')
    parser.add_argument('--mag_fna_dir', default='', help='Folder with contig nucleotide fasta files, one for each '
                                                          'genome or mag. ORFs will be predicted from these files '
                                                          'with faa files saved in --mag_faa_dir. If you already '
                                                          'have .faa files, you can omit this argument.')
    parser.add_argument('--mag_faa_dir', required=True, help='Folder with aminoacid fasta files, one for each genome or '
                                                             'mag.')
    parser.add_argument('--mag_file_extension', default='.faa', help='extension of input aminoacid or nucleotide fasta '
                                                                     'files (default ".faa")')
    parser.add_argument('--delimiter', default='~', help='character to separate filenames and orfnames during '
                                                         'orthologue calling, default ~')
    parser.add_argument('--min_fraction_of_taxa_represented', default=0.1, help='Minimum fraction of total taxa that should be '
                                                                                'represented in a cluster (default: 0.1)')
    parser.add_argument('--min_taxa_represented', default=3, help='Minimum # of taxa that should be represented in a cluster '
                                                                    '(default: 3)')
    parser.add_argument('--min_fraction_of_genes_per_taxon', default=0.0, help='Reject taxa that are not represented in a '
                                                                              'minimum fraction of proteins.')
    parser.add_argument('--keep_identical', default=False, action='store_true', help='Whether identical sequences are'
                                                                                     'included in the final concatenated '
                                                                                     'alignment (default: False)')
    return parser.parse_args()


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


def prep_hmms(hmm_dir):
    log('extracting hmms included in python package...')
    hmm_file = hmm_dir / 'conserved_genes.hmm'
    print(Path(__file__))
    with open(hmm_file, 'w') as writer:
        for hmms in (Path(__file__).parent / 'database').glob('*.hmm'):

            if hmms.name.startswith('ribosomal') or hmms.name.startswith('rpo'):
                continue
            print(hmms)
            run_external(f'hmmconvert {hmms}', stdout=writer)
    run_external(f'hmmpress -f {hmm_file}')
    return hmm_file


def predict_orfs(nt_dir: Path, aa_dir: Path, file_extension: str):
    aa_dir.mkdir(exist_ok=True, parents=True)
    for nt_file in nt_dir.glob('*' + file_extension):
        if " " in nt_file.name:
            new_name = nt_file.name.replace(' ', '_')
            shutil.move(nt_file, nt_file.parent / new_name)
            log(f'renamed "{nt_file.name}" to "{new_name}"')
            nt_file = nt_file.parent / new_name
        result_file  =aa_dir / (nt_file.stem + ".faa")
        if not result_file.exists():
            run_external(f'prodigal -m -f gff -q -i {nt_file} -a {result_file}')


def collect_seqs(hmm_file: Path, fasta_dir: Path, genes_dir: Path, file_extension, delimiter: str = ''):
    log(f'Now collecting target sequences from fasta files ending in {file_extension} in {fasta_dir} using hmmscan...')
    all_orfs = {}
    hmm_lengths = {}
    with open(hmm_file) as handle:
        for line in handle:
            if line.startswith('NAME'):
                name = line.split()[1]
            elif line.startswith('LENG'):
                hmm_lengths[name] = int(line[4:].strip())
    log(f'HMM database has {len(hmm_lengths)} HMM profiles.')


    if file_number := len([f for f in fasta_dir.glob(f'*{file_extension}')]):
        log(f'Processing {file_number} files...')
    else:
        log(f'No files with extension {file_extension} in dir {fasta_dir}, aborting!')
        exit(1)

    for file in sorted(fasta_dir.glob(f'*{file_extension}')):
        fasta_file = fasta_dir / file.name
        hmm_results_file = genes_dir / f'{file.name}.hmm-results'
        gene_file = genes_dir / (file.name + file_extension)
        unique_orf_ids = set()
        with FastaParser(fasta_file) as fasta_reader:
            for orf in fasta_reader:
                all_orfs[orf['id']] = orf
                if delimiter and delimiter in orf['id']:
                    raise Exception('Delimiter should not occur in ids of any fasta seq.')
                if orf['id'] in unique_orf_ids:
                    raise Exception(f'Fasta seq ids should be unique within each file. Duplicate: {orf["id"]}')
                unique_orf_ids.add(orf['id'])
        if not hmm_results_file.exists() or not hmm_results_file.stat().st_size:
            run_external(f'hmmscan -E 1e-25 --domtblout {hmm_results_file} {hmm_file} {fasta_file}')
        if not gene_file.exists() or not gene_file.stat().st_size:
            orfs_already_done = set()
            with TabularBlastParser(hmm_results_file, 'HMMSCAN_DOM_TABLE') as handle:
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


def main():
    print(f'This is comparative_genomics.py {VERSION}')
    args = parse_arguments()

    working_dir = Path(os.getcwd())
    genes_dir = working_dir/ 'genes'
    orthologue_dir = working_dir/ 'orthologues'
    alignments_dir = working_dir/ 'alignments'
    tmp_dir = working_dir/ 'tmp'
    hmm_dir = working_dir/ 'hmm'
    for d in (hmm_dir, genes_dir, orthologue_dir, alignments_dir, tmp_dir):
        d.mkdir(exist_ok=True)

    hmm_file = prep_hmms(hmm_dir)

    fasta_aa_dir = Path(args.mag_faa_dir)
    input_file_extension = args.mag_file_extension
    if args.mag_fna_dir:
        fasta_nt_dir = Path(args.mag_fna_dir)
        predict_orfs(fasta_nt_dir, fasta_aa_dir, input_file_extension)
        input_file_extension = '.faa'

    delimiter = args.delimiter

    collect_seqs(hmm_file, fasta_aa_dir, genes_dir, input_file_extension, delimiter)

    cluster_dir = working_dir / 'hmm_clusters'
    cluster_dir.mkdir(exist_ok=True)
    for file in cluster_dir.glob('*'):
        if not file.is_dir():
            file.unlink()
    merged_fasta_file = cluster_dir / 'merged_and_coded.fasta'
    taxa_by_orf_id = []
    merge_and_code_fasta_input(mag_faa_dir = Path(genes_dir),
                               mag_faa_file_extension = input_file_extension,
                               delimiter = args.delimiter,
                               taxa_by_orf_id = taxa_by_orf_id,
                               merged_fasta_file = merged_fasta_file
                               )
    cluster(fraction_id = 0.6,
            fraction_overlap = 0.8,
            min_fraction_orthologues = 0.8,
            min_fraction_of_taxa_represented = float(args.min_fraction_of_taxa_represented),
            min_taxa_represented = int(args.min_taxa_represented),
            min_fraction_of_genes_per_taxon = float(args.min_fraction_of_genes_per_taxon),
            include_paralogues_in_fasta_output = False,
            taxa_by_orf_id = taxa_by_orf_id,
            input_fasta_file = merged_fasta_file,
            tmp_dir = tmp_dir,
            fasta_output_dir=orthologue_dir,
            )

    align(orthologue_dir, alignments_dir)
    concatenate_alignments(alignments_dir, args.delimiter)


    log('Done!')
    log('After running this script, you can for example use:')
    log('raxmlHPC-PTHREADS -s alignments/concatenated_alignment -n tree -m PROTGAMMALG -f a -p 13 -x 123 -# 100 -T 20')
    log('iqtree2 -T 20 -m TEST -s alignments/concatenated_alignment -B 1000 -alrt 1000 --keep-ident')

if __name__ == "__main__":
    main()