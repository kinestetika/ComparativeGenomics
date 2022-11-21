import gzip
import re
from pathlib import Path

NON_IUPAC_RE_NT = re.compile(r'[^ACTGN]')
NON_IUPAC_RE_AA = re.compile(r'[^RHKDESTNQCUGPAVILMFYW]')
ALL_MASK_TARGETS = set('CDS rRNA tRNA tmRNA ncRNA repeat crispr_repeat retrotransposon'.split())
COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',}


def reverse_complement(seq:str) -> str:
    return ''.join(COMPLEMENT.get(b, 'N') for b in reversed(seq))


def _parse_fasta_header(line:str) -> dict:
    si = line.find(' ')
    if si > 0:
        return {'id': line[1:si], 'seq': '', 'descr': line[si + 1:].strip()}
    else:
        return {'id': line[1:].strip(), 'seq': '', 'descr': ''}


class FastaParser:
    def __init__(self, path, cleanup_seq=False):
        self.path = path
        self.handle = None
        self.cleanup_seq = cleanup_seq
        self.alphabet = None
        self.unknown_char = ''

    def __enter__(self):
        if str(self.path).endswith('.gz'):
            self.handle = gzip.open(self.path, 'rt')
        else:
            self.handle = open(self.path)
        return self

    def __exit__(self, exception_type, exception_value, exception_traceback):
        self.handle.close()

    def __iter__(self):
        seq_rec = None
        seq = []
        while line := self.handle.readline():
            line = line.strip()
            if line.startswith('>'):
                if len(seq) and seq_rec is not None:
                    seq_rec['seq'] = self._cleanup(''.join(seq))
                    seq = []
                    yield seq_rec
                seq_rec = _parse_fasta_header(line)
            elif seq_rec is not None:
                seq.append(line)
                if len(line) > 10 and not self.alphabet and self.cleanup_seq:
                    seq_nt_errors = sum(1 for match in NON_IUPAC_RE_NT.finditer(line))
                    seq_aa_errors = sum(1 for match in NON_IUPAC_RE_AA.finditer(line))
                    if seq_nt_errors <= seq_aa_errors:
                        self.alphabet = NON_IUPAC_RE_NT
                        self.unknown_char = 'N'
                    else:
                        self.alphabet = NON_IUPAC_RE_AA
                        self.unknown_char = 'X'
        if seq_rec is not None:
            seq_rec['seq'] = self._cleanup(''.join(seq))
            yield seq_rec

    def _cleanup(self, seq) -> str:
        if self.cleanup_seq:
            seq = seq.upper()
            seq = ''.join(seq.split())
            if self.alphabet:
                seq = self.alphabet.sub(self.unknown_char, seq)
        return seq


def write_fasta(handle, fasta, line_length=80):
    if not fasta:
        raise Exception('Attempt to write zero-length sequence to fasta.')
    handle.write(f'>{fasta["id"]}')
    try:
        if fasta["descr"] and isinstance(fasta["descr"], str):
            handle.write(f' {fasta["descr"]}')
    except KeyError:
        pass
    try:
        if fasta["subsystems"]:
            handle.write(f' ({fasta["subsystems"]})')
    except KeyError:
        pass
    try:
        if fasta["taxon"]:
            handle.write(f' ({fasta["taxon"]})')
    except KeyError:
        pass
    handle.write('\n')
    for i in range(0, len(fasta['seq']), line_length):
        handle.write(fasta['seq'][i:i+line_length])
        handle.write('\n')
