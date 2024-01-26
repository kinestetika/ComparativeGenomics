import argparse
import subprocess
import time
from pathlib import Path

from comparative_genomics.fasta import FastaParser, write_fasta
from comparative_genomics.blast import TabularBlastParser

VERSION = "0.14"
START_TIME = time.monotonic()
LOG_FILE = Path('log.txt')

def parse_arguments():
    parser = argparse.ArgumentParser(description='orthologues.py. (C) Marc Strous, 2024')
    parser.add_argument('--mag_faa_dir', help='Dir with aminoacid fasta files, one for each genome or mag.')
    parser.add_argument('--mag_faa_file_extension', default='.faa', help='Extension of genome aminoacid fasta files (default ".faa")')
    parser.add_argument('--delimiter', default='~', help='This script will create a fasta file for each cluster of homologuous'
                                                         'proteins. In those files, sequence ids will consist of a file name and'
                                                         'the previous sequence id, separated with the delimiter (default(~)')
    parser.add_argument('--cluster_dir', default='clustering', help='This dir will be used to store intermediate files that '
                                                                    'could potentially be useful (default: clustering).')
    parser.add_argument('--tmp_dir', default='tmp', help='This dir will be used to store intermediate files (default: tmp)')
    parser.add_argument('--fasta_output_dir', help='This script will write a fasta file for each cluster of homologous proteins '
                                                   'to this dir.')
    parser.add_argument('--min_fraction_id', default=0.5, help='Minimum fraction id during clustering by mmseqs (default: 0.5)')
    parser.add_argument('--min_fraction_overlap', default=0.8, help='Minimum coverage during clustering by mmseqs (default: 0.8)')
    parser.add_argument('--min_fraction_orthologues', default=0.5, help='Minimum fraction of proteins in a cluster that should '
                                                                        'be orthologues (default: 0.5)')
    parser.add_argument('--min_fraction_of_taxa_represented', default=0.5, help='Minimum fraction of total taxa that should be '
                                                                                'represented in a cluster (default: 0.1)')
    parser.add_argument('--min_taxa_represented', default=0.5, help='Minimum # of taxa that should be represented in a cluster '
                                                                    '(default: 3)')
    parser.add_argument('--include_paralogues_in_fasta_output', default=True,  action="store_false",
                        help='Whether to include paralogues in the output cluster fasta files (default: True).')
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


def merge_and_code_fasta_input(mag_faa_dir: Path, mag_faa_file_extension: str, delimiter: str, taxa_by_orf_id: list,
                               merged_fasta_file: Path):
    log(f'Now merging and coding fasta input as {merged_fasta_file}..')
    unique_ids = set()
    with open(merged_fasta_file, 'w') as writer:
        file_count = 0
        orf_count = 0
        for file in sorted(mag_faa_dir.glob(f'*{mag_faa_file_extension}')):
            fasta_file = mag_faa_dir / file.name
            with FastaParser(fasta_file) as fasta_reader:
                for orf in fasta_reader:
                    if orf['id'] in unique_ids:
                        log(f'warning: duplicate id for {orf["id"]} in {file.name}: skipping.')
                        continue
                    unique_ids.add(orf['id'])
                    if orf['seq'][-1] == '*':
                        orf['seq'] = orf['seq'][:-1]
                    recoded_orf = {'id': f'{orf_count}',
                                   'seq': orf['seq'],
                                   'descr': f'{file.stem}{delimiter}{orf["id"]} {orf["descr"]}'.strip()
                                  }
                    write_fasta(writer, recoded_orf)
                    taxa_by_orf_id.append(file_count)
                    orf_count += 1
            file_count += 1
    log(f'Merging and coding fasta input complete; recoded and wrote {len(taxa_by_orf_id)} proteins to fasta, '
        f'{file_count} unique taxa.')


class Cluster:
    def __init__(self, seq_ids: list, taxa_by_orf_id: list, blast_scores: dict):
        self.id = seq_ids[0]
        try:
            self.seq_ids = sorted(seq_ids, key=lambda k: blast_scores[self.id][k][0], reverse=True)
            self.fraction_id = sum([blast_scores[self.id][id2][1] for id2 in self.seq_ids[1:]]) / (len(self.seq_ids)-1)
            self.error = False
        except KeyError:
            log('Warning: no alignment information for cluster - could not compute percentage id and orthologues '
                  'could not be called for this cluster.')
            self.error = True
            self.seq_ids = seq_ids
            self.fraction_id = 0
        self.paralogues = {}
        self.taxa = set()
        for seq_id in self.seq_ids:
            self.paralogues[seq_id] = taxa_by_orf_id[seq_id] in self.taxa
            self.taxa.add(taxa_by_orf_id[seq_id])
        self.fraction_orthologues = len(self.taxa) / len(self.seq_ids)
        # print(f'cluster at {self.fraction_id:.1%}')
        # for i in range(len(self.seq_ids)):
        #    print(self.seq_ids[i], taxa_by_orf_id[self.seq_ids[i]], blast_scores[self.id][self.seq_ids[i]], self.paralogues[i])


class MMSeqsClusterParser:
    # This class parses a *_cluster.tsv file created by "mmseqs easy-cluster ..."
    def __init__(self, path, taxa_by_orf_id: list, blast_scores: dict):
        self.path = path
        self.handle = None
        self.taxa_by_orf_id = taxa_by_orf_id
        self.blast_scores = blast_scores

    def __enter__(self):
        self.handle = open(self.path)
        return self

    def __exit__(self, exception_type, exception_value, exception_traceback):
        self.handle.close()

    def __iter__(self):
        current_cluster_seq_id_list = []
        while line := self.handle.readline():
            id1, id2 = line.strip().split('\t')
            if id1 == id2:
                if len(current_cluster_seq_id_list) > 1:
                    yield Cluster(current_cluster_seq_id_list, self.taxa_by_orf_id, self.blast_scores)
                current_cluster_seq_id_list = [int(id1)]
            elif len(current_cluster_seq_id_list):
                current_cluster_seq_id_list.append(int(id2))
        if len(current_cluster_seq_id_list) > 1:
            yield Cluster(current_cluster_seq_id_list, self.taxa_by_orf_id, self.blast_scores)


def cluster(taxa_by_orf_id: list, input_fasta_file:Path, tmp_dir: Path, fasta_output_dir: Path,
            fraction_id:float=0.5, fraction_overlap:float=0.8,
            min_fraction_orthologues:float=0.5, min_fraction_of_taxa_represented:float=0.1,
            min_taxa_represented:float=3, include_paralogues_in_fasta_output=True):
    unique_taxa = {taxon for taxon in taxa_by_orf_id}
    log(f'Detected {len(unique_taxa)} unique taxa.')
    # prep files
    tmp_dir.mkdir(exist_ok=True)
    cluster_file_base = input_fasta_file.parent / f'{input_fasta_file.stem}.clustering'
    cluster_file_tsv = input_fasta_file.parent / f'{input_fasta_file.stem}.clustering.tsv'
    cluster_file_align = input_fasta_file.parent / f'{input_fasta_file.stem}.clustering.align'
    cluster_file_blast = input_fasta_file.parent / f'{input_fasta_file.stem}.clustering.blast'
    input_fasta_db = input_fasta_file.parent / (input_fasta_file.stem + '.mmseqdb')
    mmseqs_path = '' # '/bio/bin/mmseqs/bin/'
    # run programs
    run_external(f'{mmseqs_path}mmseqs createdb {input_fasta_file} {input_fasta_db}')
    run_external(f'{mmseqs_path}mmseqs cluster -c {fraction_overlap} --cov-mode 0 --min-seq-id {fraction_id} '
                 f'{input_fasta_db} {cluster_file_base} {tmp_dir}')
    run_external(f'{mmseqs_path}mmseqs createtsv {input_fasta_db} {input_fasta_db} {cluster_file_base} {cluster_file_tsv}')
    run_external(f'{mmseqs_path}mmseqs align {input_fasta_db} {input_fasta_db} {cluster_file_base} '
                 f'{cluster_file_align} -a')
    run_external(f'{mmseqs_path}mmseqs convertalis {input_fasta_db} {input_fasta_db} {cluster_file_align} '
                 f'{cluster_file_blast}')
    # parse results
    blast_scores = {}
    with TabularBlastParser(cluster_file_blast, 'BLAST') as reader:
        for r in reader:
            if len(r.hits) > 1:
                blast_scores[int(r.hits[0].query)] = {int(h.hit): (h.score, h.percent_id) for h in r.hits}
    seq_id2cluster = {}
    seq_id2cluster_rejected = {}
    rejected_cluster_count = 0
    error_count = 0
    cluster_count = 0
    percent_id = 0
    with MMSeqsClusterParser(cluster_file_tsv, taxa_by_orf_id, blast_scores) as reader:
        for cluster in reader:
            if cluster.fraction_orthologues < min_fraction_orthologues \
                    or len(cluster.taxa) < max(min_taxa_represented, min_fraction_of_taxa_represented * len(unique_taxa))\
                    or cluster.error:
                rejected_cluster_count += 1
                error_count += cluster.error
                for seq_id in cluster.seq_ids:
                    seq_id2cluster_rejected[seq_id] = cluster
                cluster.id = f'rejected_{rejected_cluster_count}'
            else:
                cluster.id = cluster_count
                cluster_count += 1
                percent_id += cluster.fraction_id
                for seq_id in cluster.seq_ids:
                    seq_id2cluster[seq_id] = cluster
    log(f'{len(seq_id2cluster)}/{len(taxa_by_orf_id)} seqs clustered ({len(seq_id2cluster)/len(taxa_by_orf_id):.1%}).')
    log(f'Accepted {cluster_count} clusters, rejected {rejected_cluster_count}, {error_count} due to errors.')
    log(f'Estimated average percent id: {percent_id/cluster_count:.1%}')
    # write output
    fasta_output_dir.mkdir(exist_ok=True)
    rejected_fasta_output_dir = fasta_output_dir / 'rejected'
    rejected_fasta_output_dir.mkdir(exist_ok=True)
    for file in rejected_fasta_output_dir.glob('*'):
        if not file.is_dir():
            file.unlink()
    for file in fasta_output_dir.glob('*'):
        if not file.is_dir():
            file.unlink()
    with FastaParser(input_fasta_file) as fasta_reader:
        for orf in fasta_reader:
            if int(orf['id']) in seq_id2cluster.keys():
                cluster = seq_id2cluster[int(orf['id'])]
                target_fasta_file = fasta_output_dir / f'{cluster.id}.faa'
            elif int(orf['id']) in seq_id2cluster_rejected.keys():
                cluster = seq_id2cluster_rejected[int(orf['id'])]
                target_fasta_file = rejected_fasta_output_dir / f'{cluster.id}.faa'
            else:
                continue
            if include_paralogues_in_fasta_output or not cluster.paralogues[int(orf['id'])]:
                paralogue_text = '(PARALOGUE)' if cluster.paralogues[int(orf['id'])] else ''
                with open(target_fasta_file, 'a') as writer:
                    try:
                        descr_space_index = orf['descr'].index(' ')
                        orf['id'] = orf['descr'][:descr_space_index]
                        orf['descr'] = paralogue_text + orf['descr'][descr_space_index+1:]
                    except ValueError:
                        orf['id'] = orf['descr']
                        orf['descr'] = paralogue_text
                    write_fasta(writer, orf)


def main():
    print(f'This is orthologues_v2.py {VERSION}')
    args = parse_arguments()
    cluster_dir = Path(args.cluster_dir)
    merged_fasta_file = cluster_dir / 'merged_and_coded.fasta'
    taxa_by_orf_id = []
    merge_and_code_fasta_input(mag_faa_dir = Path(args.mag_faa_dir),
                               mag_faa_file_extension = args.mag_faa_file_extension,
                               delimiter = args.delimiter,
                               taxa_by_orf_id = taxa_by_orf_id,
                               merged_fasta_file = merged_fasta_file
                               )
    cluster(fraction_id = float(args.min_fraction_id),
            fraction_overlap = float(args.min_fraction_overlap),
            min_fraction_orthologues = float(args.min_fraction_orthologues),
            min_fraction_of_taxa_represented = float(args.min_fraction_of_taxa_represented),
            min_taxa_represented = int(args.min_taxa_represented),
            include_paralogues_in_fasta_output = bool(args.include_paralogues_in_fasta_output),
            taxa_by_orf_id = taxa_by_orf_id,
            input_fasta_file = merged_fasta_file,
            tmp_dir = Path(args.tmp_dir),
            fasta_output_dir=Path(args.fasta_output_dir),
            )


if __name__ == "__main__":
    main()