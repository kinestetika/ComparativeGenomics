import gzip


class BlastHit:
    FIELDS = ('query', 'hit', 'percent_id', 'aligned_length', 'mismatches', 'gaps', 'query_start',
              'query_end', 'hit_start', 'hit_end', 'evalue', 'score')

    def __init__(self, query: str, hit: str, percent_id: float, aligned_length: int,
                 query_start: int, query_end: int, hit_start: int, hit_end: int, evalue: float, score: float,
                 mismatches: int = 0, gaps: int = 0):
        self.query = query
        self.hit = hit
        self.percent_id = percent_id
        self.aligned_length = aligned_length
        self.mismatches = mismatches
        self.gaps = gaps
        self.query_start = query_start
        self.query_end = query_end
        self.hit_start = hit_start
        self.hit_end = hit_end
        self.evalue = evalue
        self.score = score

    def __iter__(self):
        return ((k, v) for k, v in zip(BlastHit.FIELDS, (self.query, self.hit, self.percent_id, self.aligned_length,
                                                         self.mismatches, self.gaps, self.query_start, self.query_end,
                                                         self.hit_start, self.hit_end, self.evalue, self.score)))

    def __repr__(self):
        return '{}({})'.format(type(self).__name__, ', '.join(f'{k}={v!r}' for (k, v) in self))

    def __len__(self):
        return self.aligned_length


class BlastResult:
    def __init__(self, hits: tuple[BlastHit]):
        self.hits = hits
        if not len(hits):
            raise Exception('Attempt to create empty blast result.')

    def __iter__(self):
        return self.hits.__iter__()

    def __len__(self):
        return len(self.hits)

    def __repr__(self):
        return '{}(({},))'.format(type(self).__name__, ',\n'.join(f'{h!r}' for h in self))

    def query(self):
        return self.hits[0].query


class TabularBlastParser:
    def __init__(self, filename, mode):
        self.filename = filename
        self.mode = mode
        self.handle = None

    def __enter__(self):
        if str(self.filename).endswith('.gz'):
            self.handle = gzip.open(self.filename, 'rt')
        else:
            self.handle = open(self.filename)
        return self

    def __exit__(self, exception_type, exception_value, exception_traceback):
        self.handle.close()

    def __iter__(self):
        all_hits: list[BlastHit] = []
        while next_hit := self.load_next_hit_from_file():
            if not len(all_hits) or all_hits[-1].query == next_hit.query:
                all_hits.append(next_hit)
            else:
                yield BlastResult(tuple(all_hits))
                all_hits = [next_hit]
        if len(all_hits):
            yield BlastResult(tuple(all_hits))

    def load_next_hit_from_file(self) -> BlastHit | None:
        while line := self.handle.readline():
            words = line.strip().split()  # \t works for blast, but hmmer and cmscan use spaces
            match words:
                case [word, *_] if word.startswith('#'):
                    continue
                case [query, hit, percent_id, aligned_length, mismatches, gaps, query_start, query_end, hit_start,
                      hit_end, evalue, score] if 'BLAST' == self.mode:
                    # print(hit_db_entry)
                    b = BlastHit(query, hit, float(percent_id), int(aligned_length), int(query_start),
                                 int(query_end), int(hit_start), int(hit_end), float(evalue), float(score),
                                 int(mismatches), int(gaps))
                    return(b)
                case [hit, _, query, _, evalue, score, _, _, _, _, _, _, _, _, _, _, _, _, *_] if 'HMMSCAN' == self.mode:
                    return BlastHit(query, hit, 0, 0, 0, 0, 0, 0, float(evalue), float(score), 0, 0)
                case [query, _, hit, _, evalue, score, _, _, _, _, _, _, _, _, _, _, _, _, *_] if 'HMMSEARCH' == self.mode:
                    return BlastHit(query, hit, 0, 0, 0, 0, 0, 0, float(evalue), float(score), 0, 0)
                case [hit, _, _, query, _, _, evalue, score, _, _, _, _, _, _, _, hit_start, hit_end,
                      query_start, query_end, *_] if 'HMMSCAN_DOM_TABLE' == self.mode:
                    hit_start = int(hit_start)
                    hit_end = int(hit_end)
                    return BlastHit(query, hit, 0, hit_end - hit_start, int(query_start), int(query_end),
                                    hit_start, hit_end, float(evalue), float(score))
                case [*_]:
                    continue
        return None
