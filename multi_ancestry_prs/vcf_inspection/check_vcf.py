from gzip import open as gzopen
from multiprocessing import Pool, Lock
import sys, os, re, sqlite3
from os.path import basename, exists
from functools import partial

from collections import defaultdict as ddict, Counter

# name
_ADD_CONTIG_ = """
INSERT INTO contig (name) VALUES (?)
"""
_SET_COUNT_ = """
UPDATE contig
SET loci = ?
WHERE cid = ?
"""

# data, sample column, locus line
_ADD_ISSUE_ = """
INSERT INTO issue VALUES (?, ?, ?)
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
dbcrs = None
dbcon = None
name = ""

def main():
  global dblock
  f = sys.argv[1]
  threads = int(sys.argv[2]) if len(sys.argv) > 2 else 1

  gz = f.endswith('.gz')
  base = re.sub('\.vcf(\.gz)?', "", basename(f))
  dbn = base + '.issues.db'

  create_database(dbn)
  vcfopen = partial(gzopen, mode='rt') if gz else open
  contigs = []
  with vcfopen(f) as inp:
    for line in inp:
      if line.startswith('##contig'):
        contigs.append(re.search('ID=([^,]+),', line)[1])
      if line.startswith('##'): continue
      data = line.strip().split()
      samples = data[9:]
      break

  con = sqlite3.connect(dbn)
  con.executemany(_ADD_SAMPLE_, [(i, s) for i,s in enumerate(samples)])
  con.commit()
  con.close()

  if threads == 1:
    process_vcf(f, dbn, gz=gz)
  else:
    dblock = Lock()
    p = Pool(threads)
    for c in contigs:
      p.apply_async(process_vcf, (f, dbn, c, gz, True))

    p.close()
    p.join()



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


def process_vcf(f, dbn, contig=None, gz=True, lock=False):
  global name, dbcon, dbcrs
  transact = safe_execute if lock else execute
  vcfopen = partial(gzopen, mode='rt') if gz else open
  if contig:
    print("Processing", contig)
    name = contig
  
  dbcon = sqlite3.connect(dbn)
  dbcrs = dbcon.cursor()

  loci = 0
  lxrm = None
  cid = None
  with vcfopen(f) as inp:

    switch = False
    for j, line in enumerate(inp):
      if not line: continue
      if line.startswith('#'): continue
      if contig and not line.startswith(contig): # there's probably a better way to do this using the index
        if switch: break
        continue
      switch = True

      xrm, pos, form, issues = check_line(line)

      if xrm != lxrm:
        if lxrm:
          transact(_SET_COUNT_, (loci, cid))
          loci = 0
        transact(_ADD_CONTIG_, (xrm,))
        cid = dbcrs.lastrowid
        lxrm = xrm

      loci += 1

      if not issues: continue
      transact(_ADD_LOCUS_, (j, cid, pos, form))
      transact(_ADD_ISSUE_, [(d, c, j) for c, d in issues], many=True)

    if lxrm:
      transact(_SET_COUNT_, (loci, cid))
  
  dbcon.close()


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
                line INTEGER PRIMARY KEY,
                cid REFERENCES contig,
                pos INT,
                format VARCHAR(100)
              )""")
  
  cur.execute("""
              CREATE TABLE sample(
                column INTEGER PRIMARY KEY,
                name VARCHAR(20)
              )""")
  
  cur.execute("""
              CREATE TABLE issue(
                entry VARCHAR(100),
                column INT REFERENCES sample,
                line INT REFERENCES locus
              )""")
  
  con.commit()
  con.close()


if __name__ == '__main__':
  main()
  