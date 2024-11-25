import kterra as kt
from google.cloud import storage

from multiprocessing import Pool
from posixpath import join, basename
import argparse


c_workspace: kt.Workspace = None
c_wbucket: storage.Bucket = None

c_gclient: storage.Client = None
c_archive: storage.Bucket = None


def setup_connections(c_ws_n, c_a_n):
  global c_workspace, c_gclient, c_archive

  c_workspace = kt.Workspace(c_ws_n)
  c_wbucket = c_workspace.folder.gbucket

  c_gclient = c_wbucket.client
  c_archive = c_gclient.bucket(c_a_n)


def stubborn_update(entity_type, entity, _tries=3, **data):
  for i in range(_tries):
    try:
      c_workspace.update_entity(entity_type, entity, **data)
      if c_workspace.check_request(): break
    except:
      pass


def extern_model_file(entity_type, entity, column, datum, ext_form):
  
  srce_blob_name = datum.removeprefix(f"gs://{c_wbucket.name}")
  srce_blob = c_wbucket.blob(srce_blob_name)
  srce_exists = srce_blob.exists()
  
  dest_blob_name = ext_form.format(
    entity_type = entity_type,
    entity = entity,
    column = column,
    file_name = basename(datum),
    workspace = c_workspace.name
  )
  dest_blob = c_archive.blob(dest_blob_name)
  dest_exists = dest_blob.exists()
  new_datum = join(f"gs://{c_archive.name}", dest_blob_name)

  if srce_exists:
    if not dest_exists:
      bucket_transfer_blob(srce_blob, dest_blob)
    stubborn_update(entity_type, entity, **{column: new_datum})
    c_wbucket.delete_blob(srce_blob_name)
  else:
    if dest_exists:
      stubborn_update(entity_type, entity, **{column: new_datum})
    else:
      print("Could not locate", datum, "for transfer!")


def bucket_transfer_blob(srce_blob: storage.Blob, dest_blob: storage.Blob):
  token, _, _ = dest_blob.rewrite(srce_blob)
  while token != None:
    token, _, _ = dest_blob.rewrite(srce_blob, token=token)

  if not dest_blob.exists():
    print('Failed to copy', srce_blob.name, 'to', dest_blob.name)
    raise RuntimeError(f"Failed to copy {srce_blob.name} from {srce_blob.bucket.name}" + \
                       f" to {dest_blob.name} in {dest_blob.bucket.name}")


def parallel_extern(fromWorkspace, toArchive, entity_type, column, ext_form, N):
  setup_connections(fromWorkspace, toArchive)

  entities = c_workspace.get_table(entity_type)
  data = entities[column]
  arglist = [(entity_type, entity, column, datum, ext_form) for entity, datum in data.items()]

  with Pool(N, 
            initializer=setup_connections,
            initargs=(fromWorkspace, toArchive)) as p:
    res = p.starmap_async(extern_model_file, arglist)


if __name__ == "__main__":
  parcel = argparse.ArgumentParser()

  parcel.add_argument('workspace',
    help="The name of the workspace storing the files")
  parcel.add_argument('entity_type',
    help="The relevant entity table")
  parcel.add_argument('column',
    help="The target column of files")
  parcel.add_argument('archive',
    help="The bucket name of the archive")
  parcel.add_argument('ext_form',
    help="A template string for the destination file names. " +
      "E.g. \"{entity}/{file_name}\"; available components are " +
      "workspace, entity_type, entity, column, and file_name")
  parcel.add_argument('-p', '--processes', default=2,
    help="Number of processes to use. For large files (>2GB), " +
        "significant time may be spent waiting for server responses. " +
        "In this case it is safe to exceed available cpus.")
  
  args = parcel.parse_args()
  parallel_extern(args.workspace,
                  args.archive,
                  args.entity_type,
                  args.column,
                  args.ext_form,
                  args.processes)

