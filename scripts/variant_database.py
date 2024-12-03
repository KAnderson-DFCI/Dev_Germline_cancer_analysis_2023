from gzip import open as gzopen
from multiprocessing import Pool, Lock, Queue
import sys, os, re, sqlite3
from os.path import basename, exists
from functools import partial

from collections import defaultdict as ddict, Counter
import random # used for randomly naming child processes in log output.

# name
_ADD_CONTIG_ = """
INSERT INTO contig (name, loci) VALUES (?, 0)
"""
_ADD_COUNT_ = """
UPDATE contig
SET loci = contig.loci + ?
WHERE cid = ?
"""
_GET_CONTIG_ID_ = """
SELECT cid
FROM contig
WHERE name = ?
"""

# data, sample column, locus line
_ADD_VARIANT_ = """
INSERT INTO variant VALUES (?, ?, ?, ?)
"""
# line, contig id, position, format
_ADD_LOCUS_ = """
INSERT INTO locus VALUES (?, ?, ?, ?)
"""
# column, name
_ADD_SAMPLE_ = """
INSERT INTO sample VALUES (?, ?)
"""

dblock = None
nameq = None

name = None
vcf = None
vcfopen = None
data_top = None
data_len = None
nchunks = None

dbn = None
dbcon = None
dbcrs = None
transact = None

names1 = ['Abby', 'Beni', 'Cari', 'Davy', 'Ezri', 'Fuji', 'Gary',
          'Hidi', 'Iggy', 'Jeri', 'Kody', 'Lemi', 'Mary',
          'Niki', 'Olli', 'Puny', 'Quty', 'Remy', 'Sadi', 'Toby',
          'Urie', 'Viki', 'Wily', 'Xixi', 'Yogi', 'Zoey']
names2 = ['Arvo', 'Brio', 'Cujo', 'Dido', 'Echo', 'Fafo', 'Geko',
          'Hero', 'Iroh', 'Jelo', 'Kado', 'Lobo', 'Medo',
          'Nero', 'Oreo', 'Polo', 'Quno', 'Reno', 'Silo', 'Taro',
          'Urso', 'Volo', 'Wako', 'Xylo', 'Yoyo', 'Zero']

def main():
  vcf = sys.argv[1]
  ncpu = int(sys.argv[2]) if len(sys.argv) > 2 else 1
  chunks = max(ncpu, int(sys.argv[3])) if len(sys.argv) > 3 else ncpu

  gz = vcf.endswith('.gz')
  base = re.sub(r'\.vcf(\.gz)?', "", basename(vcf))
  vcfopen = partial(gzopen, mode='rt') if gz else open

  contigs = []
  with vcfopen(vcf) as inp:
    print('Parsing Header...')
    while True:
      line = inp.readline()
      if line.startswith('##contig'):
        contigs.append(re.search('ID=([^,>]+),?', line)[1])
      if line.startswith('##'): continue
      data = line.strip().split()
      samples = data[9:]
      break
    data_top = inp.tell()
    print(f'Found top of data at byte {data_top}')
    data_end = inp.seek(0, 2)
    print(f'Found end of data at byte {data_end}')

  data_len = data_end - data_top

  print('Initializing database...')
  dbn = base + '.issues.db'
  create_database(dbn)
  con = sqlite3.connect(dbn)
  con.executemany(_ADD_CONTIG_, [(c,) for c in contigs])
  con.executemany(_ADD_SAMPLE_, [(i, s) for i,s in enumerate(samples)])
  con.commit()
  con.close()

  print('Initializing Worker(s)...')
  global dblock, nameq
  nameq = Queue()
  for s in random.sample(range(26**2), ncpu):
    nameq.put_nowait(s)

  conf = (vcf, dbn, gz, data_top, data_len, chunks if ncpu > 1 else 1)

  if ncpu == 1:
    print('Executing in single process.')
    init_worker(*conf)
    process_vcf()
  else:
    print(f'Executing in process pool({ncpu}).')
    dblock = Lock()
    
    p = Pool(ncpu, initializer=init_worker, initargs=conf)
    for i in range(chunks):
      p.apply_async(process_vcf, (i,))

    p.close()
    p.join()


def build_name(a: int) -> str:
  return names1[a%26] + '-' + names2[a//26]

def init_worker(v, d, gz, dt, dl, nc):
  global name, vcf, data_top, data_len, nchunks, \
    dbn, dbcon, dbcrs, vcfopen, transact
  name = build_name(nameq.get())
  nchunks = nc

  vcf = v
  vcfopen = partial(gzopen, mode='rt') if gz else open
  data_top = dt
  data_len = dl

  dbn = d
  dbcon = sqlite3.connect(dbn)
  dbcrs = dbcon.cursor()
  transact = safe_execute if nchunks > 1 else execute


def merge_count_maps():
  pass


def safe_execute(stmt, args, many=False):
  print(name, 'executing')
  try:
    dblock.acquire()
    print(name, 'obtained lock')
    execute(stmt, args, many)
  except:
    print('Bonk')
    print(sys.exc_info())
  finally:
    try:
      dblock.release()
      print(name, 'released lock')
    except:
      print('Bonked')

def execute(stmt, args, many=False):
  if many:
    dbcrs.executemany(stmt, args)
  else:
    dbcrs.execute(stmt, args)
  dbcon.commit()


def get_chunk_b(funk):
  return int(data_top + data_len * (funk/nchunks))

def process_vcf(funk=0):
  if nchunks == 1:
    sb, eb = data_len, data_top
  else:
    sb = get_chunk_b(funk)
    se = get_chunk_b(funk+1)

  loci = 0
  lxrm = None
  cid = None
  with vcfopen(vcf) as inp:
    inp.seek(sb-1)
    if inp.read(1) != '\n': inp.readline()

    while inp.tell() < se:
      line = inp.readline()
      xrm, pos, ref, info, vrts = process_line(line)

      # contig transition
      if xrm != lxrm:
        if lxrm:
          transact(_ADD_COUNT_, (loci, cid))
          loci = 0
        transact(_GET_CONTIG_ID_, (xrm,))
        cid = int(dbcrs.fetchone())
        lxrm = xrm

      loci += 1
      transact(_ADD_LOCUS_, (cid, pos, ref, info))
      lid = dbcrs.lastrowid
      transact(_ADD_VARIANT_, [(*acs, lid) for acs in vrts], many=True)

    if lxrm:
      transact(_ADD_COUNT_, (loci, cid))
  
  dbcon.close()


def tabulate_rowform(form, sample_data):
  fields = form.split(':')
  return [{f: v for f,v in zip(fields, s.split(':'))} for s in sample_data]

def process_line(line):
  data = line.strip().split()
  xrm, pos = data[0:2]
  ref = data[3]
  alts = data[4].split(',')
  alcs = [vrt_type(ref, a) for a in alts]
  info = data[7]
  sample_forms = tabulate_rowform(data[8], data[9:])

  vrts = []
  for i, sdm in enumerate(sample_forms):
    gtf = sdm['GT']
    ali = int(gtf[0]), int(gtf[2])
    vals = [a-1 for a in ali if a > 0]
    if vals:
      vrts += [(alts[a], alcs[a], i) for a in vals]

  return xrm, pos, ref, info, vrts

tshape = {'a': 'ag', 'g': 'ag', 'c': 'ct', 't': 'ct'}

def vrt_type(ref, alt):
  if ref == '.':
    return 'nr'
  
  lr, la = len(ref), len(alt)
  if la < lr: return 'de'
  if la > lr: return 'in'
  if la > 1: return 'mp'
  if alt in tshape[ref]: return 'ts'
  else: return 'tv'


def check_line(line):
  data = line.strip().split()
  xrm, pos = data[0:2]
  form = data[8]
  samples = data[9:]

  fields = form.split(':')
  n = len(fields)
  
  missing = []
  for i, s in enumerate(samples):
    sdata = s.split(':')
    if len(sdata) == n: continue
    missing.append((i, s))
  
  return xrm, pos, form, missing


def create_database(dbn):
  if exists(dbn):
    os.remove(dbn)

  con = sqlite3.connect(dbn)
  cur = con.cursor()

  cur.execute("PRAGMA foreign_keys = ON")
  cur.execute("""
              CREATE TABLE contig(
                name VARCHAR(20) UNIQUE,
                loci INT,
                cid INTEGER PRIMARY KEY
              )""")
  
  cur.execute("""
              CREATE TABLE locus(
                cid REFERENCES contig,
                pos INT,
                ref VARCHAR(20),
                info VARCHAR(100),
                line INTEGER PRIMARY KEY
              )""")
  
  cur.execute("""
              CREATE TABLE sample(
                column INTEGER PRIMARY KEY,
                name VARCHAR(20)
              )""")
  
  cur.execute("""
              CREATE TABLE variant(
                alt VARCHAR(20),
                class CHAR(2),
                column INT REFERENCES sample,
                line INT REFERENCES locus
              )""")
  
  con.commit()
  con.close()


if __name__ == '__main__':
  main()

# python variant_database.py "..\\misc_data\\1KGQ_common_pop_phased.vcf.gz" 6 100