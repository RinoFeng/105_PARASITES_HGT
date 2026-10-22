#!/usr/bin/env python3
from __future__ import annotations

import os
import re
import csv
import math
import shutil
import subprocess
import gzip
from collections import defaultdict, Counter
from typing import Dict, List, Tuple, Iterable, Optional

AA_VALID = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*-.")

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}
RC = str.maketrans('ACGTUacgtuNn', 'TGCAAtgcaaNn')


def safe_id(s: str, max_len: int = 180) -> str:
    s = str(s)
    s = re.sub(r'\s+', '_', s)
    s = re.sub(r'[^A-Za-z0-9_.:|=,+\-]', '_', s)
    return s[:max_len]


def read_taxonomy_file(path: str) -> Dict[str, Dict[str, str]]:
    tax = {}
    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            row = line.split('\t')
            if len(row) < 3:
                # Some files may be whitespace separated in previews, but should be tab separated.
                row = re.split(r'\s+', line, maxsplit=2)
            if len(row) < 3:
                continue
            taxid = row[0]
            detail_field = row[1]
            detail = detail_field.split(' ', 1)[1] if ' ' in detail_field else detail_field
            short = row[2]
            tax[taxid] = {'detail': detail, 'short': short}
    return tax


def read_genome_taxid_map(path: str) -> Dict[str, Dict[str, str]]:
    """Read user's 105_name2accfil_withTaxid.txt.

    Expected columns, tab separated:
        species_name    genome_file    taxid

    Returns both genome->record and species->record under two nested keys.
    The sampling script uses genome_file first, then species_name as fallback.
    """
    by_genome: Dict[str, Dict[str, str]] = {}
    by_species: Dict[str, Dict[str, str]] = {}
    if not path:
        return {'by_genome': by_genome, 'by_species': by_species}
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            row = line.split('\t')
            if len(row) < 3:
                row = re.split(r'\s+', line, maxsplit=2)
            if len(row) < 3:
                continue
            species, genome, taxid = row[0].strip(), row[1].strip(), row[2].strip()
            rec = {'species': species, 'genome': genome, 'taxid': taxid}
            by_genome[genome] = rec
            by_species[species] = rec
    return {'by_genome': by_genome, 'by_species': by_species}


def get_self_taxid_for_target(hgt_row, genome_taxid_map: Dict[str, Dict[str, Dict[str, str]]]) -> Tuple[str, str]:
    """Return recipient/self taxid using the user's genome-taxid table.

    Priority:
      1. Genome column match to 105_name2accfil_withTaxid.txt
      2. Specie column match to 105_name2accfil_withTaxid.txt
      3. Existing self_taxid/TaxID columns if present
      4. Empty string
    """
    genome = str(hgt_row.get('Genome', '')).strip()
    specie = str(hgt_row.get('Specie', '')).strip()
    by_genome = genome_taxid_map.get('by_genome', {}) if genome_taxid_map else {}
    by_species = genome_taxid_map.get('by_species', {}) if genome_taxid_map else {}
    if genome in by_genome:
        return by_genome[genome]['taxid'], 'GENOME_TAXID_MAP:Genome'
    if specie in by_species:
        return by_species[specie]['taxid'], 'GENOME_TAXID_MAP:Specie'
    for col in ['self_taxid', 'TaxID', 'genome_taxid']:
        try:
            val = str(hgt_row.get(col, '')).strip()
            if val and val.lower() != 'nan':
                return val, f'target_column:{col}'
        except Exception:
            pass
    return '', 'missing'


def get_host_levels(self_detail: str) -> List[str]:
    if ';Lepeophtheirus salmonis' in self_detail:
        return [';Actinopterygii']
    if ';Varroa destructor' in self_detail:
        return [';Insecta;']
    if ';Aphanomyces astaci' in self_detail:
        return [';Malacostraca;']
    if ';Perkinsus marinus' in self_detail:
        return [';Bivalvia;']
    if ';Leishmania' in self_detail:
        return [';Diptera;', ';Mammalia;']
    if any(x in self_detail for x in [';Angiostrongylus', ';Fasciola', ';Schistosoma']):
        return [';Gastropoda;', ';Mammalia;']
    if any(x in self_detail for x in [';Clonorchis', ';Opisthorchis', ';Anisakis']):
        return [';Gastropoda;', ';Mammalia;', ';Actinopterygii;']
    if ';Babesia' in self_detail:
        return [';Acari;', ';Mammalia;']
    if any(x in self_detail for x in [';Plasmodium', ';Brugia', ';Wuchereria']):
        return [';Culicidae;', ';Mammalia;']
    if ';Cryptosporidium' in self_detail:
        return [';Aves;', ';Mammalia;']
    return [';Mammalia;']


def get_self_level(taxonomy_info: Dict[str, Dict[str, str]], self_taxid: str) -> str:
    info = taxonomy_info.get(str(self_taxid), {'detail': '', 'short': ''})
    detail = info.get('detail', '')
    short = info.get('short', '')
    try:
        if ';Metazoa;' in detail:
            return short.split(';')[1]
        return detail.split(';')[2]
    except Exception:
        return ''


def classify_hit(taxonomy_info: Dict[str, Dict[str, str]], self_taxid: str, target_taxid: str) -> Tuple[str, str, str]:
    """Return group, taxonomy_detail, taxonomy_short.

    Groups match HGT_verify_parasite.py exactly in spirit:
    - self: hit lineage contains recipient self_level. For Metazoa this is short[1];
      otherwise this is detail.split(';')[2].
    - out_host: known host lineages hard-coded by parasite group.
    - out_other: Bacteria/Viruses/Fungi/Archaea/Viridiplantae.
    - close: other cellular organisms not self and not out groups.
    - unknown: no usable taxonomy or not cellular organisms.

    IMPORTANT: self_taxid must be the recipient parasite taxid, supplied from
    105_name2accfil_withTaxid.txt, not the qtaxid column inferred from hits.
    """
    self_taxid = str(self_taxid)
    target_taxid = str(target_taxid)
    self_info = taxonomy_info.get(self_taxid, {'detail': '', 'short': ''})
    target_info = taxonomy_info.get(target_taxid, {'detail': '', 'short': ''})
    taxonomy = target_info.get('detail', '')
    short = target_info.get('short', '')
    if not taxonomy or not short:
        return 'unknown', taxonomy, short

    self_level = get_self_level(taxonomy_info, self_taxid)
    distant_levels = ['Bacteria;', 'Viruses;', ';Fungi;', ';Archaea;', ';Viridiplantae;']
    host_levels = get_host_levels(self_info.get('detail', ''))

    if self_level and f';{self_level};' in taxonomy:
        return 'self', taxonomy, short
    if 'organisms;' not in taxonomy:
        return 'unknown', taxonomy, short
    if any(h in taxonomy for h in host_levels):
        return 'out_host', taxonomy, short
    if any(d in taxonomy for d in distant_levels):
        return 'out_other', taxonomy, short
    return 'close', taxonomy, short


def taxonomy_rank(short: str, idx: int, default: str = 'NA') -> str:
    parts = str(short).split(';')
    if len(parts) > idx and parts[idx]:
        return parts[idx]
    return default


def species_name(short: str) -> str:
    parts = str(short).split(';')
    return parts[-1] if parts and parts[-1] else 'unknown_species'


def genus_name(short: str) -> str:
    parts = str(short).split(';')
    if len(parts) >= 2:
        return parts[-2] if parts[-2] else 'unknown_genus'
    return 'unknown_genus'


def open_text(path: str):
    return open(path, 'r', encoding='utf-8', errors='replace')


def read_fasta(path: str) -> Dict[str, str]:
    seqs = {}
    name = None
    chunks = []
    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            seqs[name] = ''.join(chunks)
    return seqs


def write_fasta(records: Iterable[Tuple[str, str]], path: str) -> None:
    with open(path, 'w') as out:
        for name, seq in records:
            seq = re.sub(r'\s+', '', str(seq))
            if not seq:
                continue
            out.write(f'>{safe_id(name)}\n')
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + '\n')


def crop_aa(seq: str, start: Optional[int], end: Optional[int], flank: int = 0) -> str:
    seq = re.sub(r'\s+', '', str(seq)).upper()
    if not start or not end:
        return seq
    a, b = sorted([int(start), int(end)])
    a = max(1, a - flank)
    b = min(len(seq), b + flank)
    if a > b or a > len(seq):
        return seq
    return seq[a-1:b]


def revcomp(nt: str) -> str:
    return nt.translate(RC)[::-1].replace('U', 'T').replace('u', 't')


def translate_nt(nt: str) -> str:
    nt = re.sub(r'[^ACGTUacgtu]', 'N', nt).upper().replace('U', 'T')
    aa = []
    for i in range(0, len(nt) - 2, 3):
        codon = nt[i:i+3]
        aa.append(CODON_TABLE.get(codon, 'X'))
    return ''.join(aa)


def translate_query_region(nt_seq: str, qstart: int, qend: int, flank_nt: int = 0) -> str:
    qstart, qend = int(qstart), int(qend)
    if qstart <= qend:
        a = max(1, qstart - flank_nt)
        b = min(len(nt_seq), qend + flank_nt)
        frag = nt_seq[a-1:b]
    else:
        a = max(1, qend - flank_nt)
        b = min(len(nt_seq), qstart + flank_nt)
        frag = revcomp(nt_seq[a-1:b])
    # Trim to the first complete codon according to blastx qstart frame.
    return translate_nt(frag)


def run_cmd(cmd: List[str], log_path: Optional[str] = None) -> None:
    if log_path:
        with open(log_path, 'w') as log:
            subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
    else:
        subprocess.run(cmd, check=True)


def _open_maybe_gzip(path: str):
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8', errors='replace')
    return open(path, 'r', encoding='utf-8', errors='replace')


def _read_id_set(id_file: str) -> set:
    ids = set()
    with open(id_file, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            x = line.strip()
            if x:
                ids.add(x)
    return ids


def extract_fasta_by_ids_python(fasta_path: str, id_file: str, out_fasta: str) -> None:
    """Stream a FASTA/FASTA.GZ file and keep records whose first header token is in id_file.

    This is a fallback when seqkit is unavailable. It is slower than seqkit on very large
    databases, but it avoids relying on DIAMOND getseq accession-list support, which is
    absent in many DIAMOND versions.
    """
    wanted = _read_id_set(id_file)
    if not wanted:
        open(out_fasta, 'w').close()
        return
    found = 0
    keep = False
    with _open_maybe_gzip(fasta_path) as fh, open(out_fasta, 'w') as out:
        for line in fh:
            if line.startswith('>'):
                header_id = line[1:].split()[0]
                keep = header_id in wanted
                if keep:
                    found += 1
                    out.write(line)
            else:
                if keep:
                    out.write(line)
    missing = len(wanted) - found
    if missing > 0:
        with open(out_fasta + '.missing_count.txt', 'w') as m:
            m.write(f'wanted\t{len(wanted)}\nfound\t{found}\nmissing\t{missing}\n')


def extract_fasta_by_ids(fasta_path: str, id_file: str, out_fasta: str, extractor: str = 'auto', seqkit_threads: int = 4, seqkit_match_name: int = 0) -> None:
    """Extract FASTA records by accession/sequence ID from the original FASTA.

    extractor: auto | seqkit | python
    - auto: use seqkit if available, otherwise Python streaming.
    - seqkit: use seqkit grep -f id_file. seqkit_threads controls -j. If seqkit_match_name=1, add -n.
    - python: pure Python streaming fallback.
    """
    if os.path.getsize(id_file) == 0:
        open(out_fasta, 'w').close()
        return
    if not fasta_path or not os.path.exists(fasta_path):
        raise FileNotFoundError(
            f'Original FASTA not found: {fasta_path}. DIAMOND getseq in many versions cannot extract by accession list; '
            f'please set PROTEIN_DB_FASTA/NT_DB_FASTA to the FASTA used to build the .dmnd database.'
        )
    use_seqkit = extractor in ('auto', 'seqkit') and shutil.which('seqkit') is not None
    if extractor == 'seqkit' and not use_seqkit:
        raise RuntimeError('EXTRACTOR=seqkit was requested, but seqkit is not in PATH.')
    if use_seqkit:
        cmd = ['seqkit', 'grep', '-j', str(seqkit_threads)]
        if int(seqkit_match_name) == 1:
            cmd.append('-n')
        cmd += ['-f', id_file, fasta_path, '-o', out_fasta]
        run_cmd(cmd)
    else:
        extract_fasta_by_ids_python(fasta_path, id_file, out_fasta)


# Backward-compatible name. This no longer calls diamond getseq by default because
# many DIAMOND versions only support `getseq --seq <internal sequence numbers>`, not
# accession-list extraction.
def diamond_getseq(db: str, id_file: str, out_fasta: str, diamond_bin: str = 'diamond') -> None:
    raise RuntimeError(
        'diamond getseq by accession list is not supported by this script version. '
        'Use extract_fasta_by_ids() with the original FASTA database instead.'
    )
