import argparse
import subprocess
import time
from tqdm import tqdm
from pathlib import Path
from math import log10
from multiprocessing import cpu_count

from fasta import FastaParser, write_fasta
from blast import TabularBlastParser

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
    parser = argparse.ArgumentParser(description='orthologues.py. (C) Marc Strous, 2022')
    parser.add_argument('--input_dir', help='Dir with aminoacid fasta files, one for each genome or mag.')
    parser.add_argument('--output_dir', help='Dir to write output to')
    parser.add_argument('--cpus', default=cpu_count(), help='How many cpus/threads to use (default: all).')
    parser.add_argument('--file_extension', default='.faa', help='extension of aminoacid fasta files (default ".faa")')
    parser.add_argument('--delimiter', default='|', help='character to separate filenames and orfnames in output')

    return parser.parse_args()


class Cluster:
    def __init__(self, parent = None):
        self.parent = parent
        self.children = []
        self.members = []
        self.is_done = False


class SetOfOrthologues:
    def __init__(self):
        self.orthologues = []
        self.paralogues = []

def make_blast_key(id1: int, id2: int) -> str:
    if id1 > id2:
        return f'{id2} {id1}'
    else:
        return f'{id1} {id2}'


def merge_and_code_fasta_input(fasta_dir: Path, file_extension: str, delimiter: str, taxa_by_orf_id: list) -> Path:
    log('Now merging and coding fasta input..')
    unique_ids = set()
    merged_fasta_file = fasta_dir / 'merged_and_coded_fasta'
    with open(merged_fasta_file, 'w') as writer:
        file_count = 0
        orf_count = 0
        for file in sorted(fasta_dir.glob(f'*{file_extension}')):
            fasta_file = fasta_dir / file.name
            with FastaParser(fasta_file) as fasta_reader:
                for orf in fasta_reader:
                    if orf['id'] in unique_ids:
                        log(f'warning: duplicate id for {orf["id"]} in {file.name}: skipping.')
                        continue
                    unique_ids.add(orf['id'])
                    recoded_orf = {'id': f'{orf_count}',
                                   'seq': orf['seq'],
                                   'descr': f' {file.name}{delimiter}{orf["id"]} {orf["descr"]}'.strip()
                                  }
                    write_fasta(writer, recoded_orf)
                    taxa_by_orf_id.append(file_count)
                    orf_count += 1
            file_count += 1
    log(f'Merging and coding fasta input complete; recoded and wrote {len(taxa_by_orf_id)} proteins to fasta.')
    return merged_fasta_file


def run_diamond(fasta_file: Path, cpus: int, unique_blast_results: dict) -> Path:
    log('Now computing homology between genes with diamond..')
    diamond_result_file = fasta_file.parent / 'diamond_results'
    run_external(f'diamond makedb --in {fasta_file} --db {fasta_file}')
    run_external(f'diamond blastp -d {fasta_file} -q {fasta_file} -o {diamond_result_file} -f 6 '
                         f'--threads {cpus} --fast --max-target-seqs 10')
    mcl_cluster_input_file = fasta_file.parent / 'mcl_input'
    with TabularBlastParser(diamond_result_file, 'BLAST') as handle:
        with open(mcl_cluster_input_file, 'w') as writer:
            for blast_result in handle:
                for hit in blast_result:
                    hit.hit = int(hit.hit)
                    hit.query = int(hit.query)
                    if hit.hit == hit.query:
                        continue
                    unique_blast_results[make_blast_key(hit.hit, hit.query)] = hit.score
                    writer.write(f'{hit.hit}\t{hit.query}\t{min(200, round(-log10(hit.evalue)))}')
    log(f'Homology computation complete; wrote {len(unique_blast_results)} pairwise evalues to file for clustering..')
    return mcl_cluster_input_file


def run_mcl(graph_input_file: Path, output_file, fineness, cpus):
    run_external(f'mcxload -abc {graph_input_file} --write-binary --stream-mirror -write-tab {graph_input_file}.tab'
                 f'-o {graph_input_file}.mci')
    run_external(f'mcl {graph_input_file}.mci-I {fineness:.1f} -t {cpus} -o {graph_input_file}.clusters')
    run_external(f'mcxdump -icl {graph_input_file}.clusters -tabr {graph_input_file}.tab  -o {output_file}')


def mcl_cluster(mcl_cluster_input_file: Path, cpus, fineness_steps=None) -> list[Cluster]:
    log('Now clustering homologous genes with mcl..')
    clustering_data_by_orf = {}
    all_clusters = []
    if not fineness_steps:
        fineness_steps = [1.2, 1.4, 2, 4, 6]
    # use mcl to cluster data de novo at each fineness step
    for fineness in range(len(fineness_steps)):
        mcl_output_file = mcl_cluster_input_file.parent / f'mcl.{fineness_steps[fineness]}'
        run_mcl(mcl_cluster_input_file, mcl_output_file, fineness_steps[fineness], cpus)
        with open(mcl_output_file) as reader:
            for line in reader:
                new_cluster = Cluster()
                for id in line.split('\t'):
                    id = int(id)
                    new_cluster.members.append(id)
                    try:
                        prev_clustering_data_for_orf = clustering_data_by_orf[id]
                    except KeyError:
                        prev_clustering_data_for_orf = {}
                        clustering_data_by_orf[id] = prev_clustering_data_for_orf
                    prev_clustering_data_for_orf[fineness] = len(all_clusters)
                all_clusters.append(new_cluster)
    log('now creating cluster hierarchy...')

    for orf_id in clustering_data_by_orf:
        parent_cluster = None
        for cluster_id in clustering_data_by_orf[orf_id].values():
            cluster_at_fineness = all_clusters[cluster_id]
            if parent_cluster:
                parent_cluster.children.append(cluster_at_fineness)
                cluster_at_fineness.parent = parent_cluster
            parent_cluster = cluster_at_fineness
    log(f'Cluster generation complete; created {sum((1 for c in all_clusters if c.parent == None))} root clusters')
    return all_clusters


def compute_cluster_score(cluster, taxa_by_orf_id) -> float:
    taxa_collected = {taxa_by_orf_id[orf_id] for orf_id in cluster['members']}
    return len(taxa_collected) - (len(cluster['members'] - len(taxa_collected)))


def make_orthologue_from_cluster(cluster:Cluster, taxa_by_orf_id: list, unique_blast_results: dict)\
        -> SetOfOrthologues:
    members_scores = {orf_id: sum((unique_blast_results.get(make_blast_key(orf_id, orf_id_2), 0)
                                   for orf_id_2 in cluster.members)) for orf_id in cluster.members}
    cluster.members.sort(key=lambda orf_id: members_scores[orf_id], reverse=True)
    taxa_collected = set()
    my_set = SetOfOrthologues()
    for orf_id in cluster.members:
        taxon = taxa_by_orf_id[orf_id]
        if taxon in taxa_collected:
            my_set.paralogues.append(orf_id)
        else:
            my_set.orthologues.append(orf_id)
            taxa_collected.add(taxon)
    return my_set


def make_orthologues_from_cluster_family(cluster: Cluster, taxa_by_orf_id: list, unique_blast_results: dict) \
        -> list[SetOfOrthologues]:
    log('Now calling orthologues from cluster families...')
    orthologues = []
    while cluster.children:
        if cluster.is_done:
            return []
        cluster.is_done = True
        score = compute_cluster_score(cluster, taxa_by_orf_id)
        best_child_score = 0
        best_child_cluster = None
        for child_cluster in cluster.children:
            child_score = compute_cluster_score(child_cluster, taxa_by_orf_id)
            if child_score > best_child_score:
                best_child_score = child_score
                best_child_cluster = child_cluster
        if best_child_score > score:
            cluster = best_child_cluster
        else:
            break
    orthologues.append(make_orthologue_from_cluster(cluster, taxa_by_orf_id, unique_blast_results))
    if cluster.parent:
        for sibling in cluster.parent.children:
            if sibling is cluster:
                continue
            orthologues.extend(make_orthologues_from_cluster_family(sibling, taxa_by_orf_id, unique_blast_results))
    log(f'Orthologue calling complete; created {len(orthologues)} sets of orthologous genes')
    return orthologues


def write_orthologues_to_fasta(merged_and_coded_fasta_file:Path, orthologues_by_orf_id: dict[int, str],
                               taxa_by_orf_id: dict, output_dir: Path, include_paralogues: bool = True):
    log('Now writing a fasta file for each set of orthologous genes...')
    if include_paralogues:
        H = {'O': '[Orthologue]', 'P' : '[Paralogue]'}
    else:
        H = {'O': '', 'P' : ''}
    with FastaParser(merged_and_coded_fasta_file) as fasta_reader:
        for orf in tqdm(fasta_reader, len(taxa_by_orf_id)):
            if result := orthologues_by_orf_id.get(int(orf['id']), 0):
                if 'P' == result[0:1] and not include_paralogues:
                    continue
                output_file = output_dir / f'{result[1:]}.faa'
                space_index = orf['descr'].index[' ']
                recoded_orf = {'id': orf['descr'][0, space_index],
                               'seq': orf['seq'],
                               'descr':  "{} {}".format(H[result[0:1]], orf['descr'][space_index + 1]).strip()
                               }
                with open(output_file, 'a') as writer:
                    write_fasta(writer, recoded_orf)


def compute_orthologues(fasta_dir: Path, cpus: int, file_extension: str = '.faa', delimiter: str = '|') -> tuple:
    taxa_by_orf_id = []
    unique_blast_results = {}
    merged_and_coded_fasta_file = merge_and_code_fasta_input(fasta_dir, file_extension, delimiter, taxa_by_orf_id)
    mcl_cluster_input_file = run_diamond(merged_and_coded_fasta_file, cpus, unique_blast_results)
    clusters = mcl_cluster(mcl_cluster_input_file, cpus)
    orthologues = []
    for cluster in clusters:
        if cluster.parent:
            continue  # only take the "root" clusters, the kids will be processed by
        orthologues.extend(make_orthologues_from_cluster_family(cluster, taxa_by_orf_id, unique_blast_results))
    orthologues_by_orf_id = {}
    for i in range(len(orthologues)):
        for orf_id in orthologues[i].orthologues:
            orthologues_by_orf_id[orf_id] = f'O{i}'
        for orf_id in orthologues[i].paralogues:
            orthologues_by_orf_id[orf_id] = f'P{i}'
    return merged_and_coded_fasta_file, taxa_by_orf_id, unique_blast_results, orthologues, orthologues_by_orf_id


def main():
    print(f'This is orthologues.py {VERSION}')
    args = parse_arguments()
    merged_and_coded_fasta_file, taxa_by_orf_id, unique_blast_results, orthologues, orthologues_by_orf_id = \
        compute_orthologues(args.input_dir, int(args.cpus), args.file_extension)
    write_orthologues_to_fasta(merged_and_coded_fasta_file, orthologues_by_orf_id, taxa_by_orf_id,
                               args.output_dir)


if __name__ == "__main__":
    main()