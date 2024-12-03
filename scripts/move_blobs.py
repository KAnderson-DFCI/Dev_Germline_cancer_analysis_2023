import kterra as kt
from google.cloud import storage
from multiprocessing import Pool

from posixpath import join, basename

dest_bucket_name = "aldubayan-qatari-wgs-2024"
ws = kt.nbWorkspace

def export_file(eid, srce_file):
  """Moves a blob from one bucket to another with a new name."""

  srce_bucket = ws.folder.gbucket
  client = srce_bucket.client

  # generic
  srce_blob_name = srce_file.removeprefix(ws.folder.gcloud_path())
  if srce_blob_name == srce_file:
    print(srce_file, 'not in srce_bucket')
    return
  
  srce_blob = srce_bucket.blob(srce_blob_name)
  if not srce_blob.exists():
    print(srce_blob_name, 'not found')
    return

  # case
  dest_blob_name = join(eid, basename(srce_file)) 
  
  # generic
  dest_bucket = client.bucket(dest_bucket_name)
  dest_blob = dest_bucket.blob(dest_blob_name)
  if dest_blob.exists():
    print(dest_blob_name, 'already exists')
    return

  blob_copy = srce_bucket.copy_blob(
    srce_blob, dest_bucket, dest_blob_name,
  )
  if not blob_copy.exists():
    print('Failed to copy', srce_blob_name, 'to destination')
    return

  srce_bucket.delete_blob(srce_blob_name)

if __name__ == "__main__":
  entities = ws.get_table('sample_qat')

  with Pool(6) as p:
    crams = entities['cram']
    p.starmap_async(export_file, crams.items())

    crais = entities['crai']
    p.starmap_async(export_file, crais.items())

    p.join()



