###
# This module should be entirely functional from within terra's notebook environment
# Instructions are provided below for supporting local operation.
#
# **Windows Instructions**
# 1. Install gcloud sdk
# 2. open gcloud commnd line<sup>*</sup> to initialize defaults and account info
# 3. run `gcloud auth application-default login` to set up application default credentials ("ADC")
# 4. open python environment and 
#    i. `pip install google-cloud-storage`
#   ii. `pip install firecloud`

import os, sys
from posixpath import join, relpath, basename, normpath
from os.path import join as osjoin
from collections import defaultdict as ddict
from collections.abc import Mapping

import pandas as pd
import numpy as np
import json

from google.cloud.storage import Client, transfer_manager as gctm, Bucket
from firecloud.fiss import fapi as fcl


### google storage clients
class ClientCache(dict):
    def __missing__(self, key):
        client = Client(key)
        self[key] = client
        return client


### firecloud workspaces
def list_workspaces(refresh=False, fields=None, fmap=None):
    global last_request, fc_workspaces
    if not refresh and not fields\
        and fc_workspaces is not None: return fc_workspaces

    preserve = False
    if not fields and not fmap:
        fields = [
            'accessLevel',
            'workspace.name',
            'workspace.namespace',
            'workspace.workspaceId',
            'workspace.bucketName',
            'workspace.createdBy',
            'workspace.billingAccount',
            'workspace.googleProject'
        ]
        fmap = key_ep
        preserve = True

    last_request = fcl.list_workspaces(",".join(fields))
    if last_request.status_code != 200:
        print('Bad Request', last_request.status_code, last_request.reason, file=sys.stderr)
        return None
    
    table = tabulate(last_request.json(), fields, fmap)
    table.set_index('name', inplace=True)
    if preserve:
        fc_workspaces = table
    return table


def getWorkspaceBucket(name):
    return fc_workspaces.loc[name].bucketName

def getBucketWorkspace(bucket):
    bucket = bucket.removeprefix('gs://')
    return fc_workspaces[fc_workspaces.bucketName == bucket].index[0]


class BucketFolder():
    def __init__(self, gbucket, path):
        if not type(gbucket) == Bucket:
            raise ValueError('Must provide a google Bucket instance to BucketFolder.')
        self.gbucket = gbucket
        if path in ['', '.', './']: self.path = ''
        else: self.path = path + ('/' if not path.endswith('/') else '')
        

    def cloud_path(self):
        return join('gs://', self.gbucket.name, self.path)
    
    def getdir(self, dest):
        _, folders = self.list_files(dest.strip('/'))
        if not len(folders) == 1:
            print(f"Could not locate folder {dest} in {self.cloud_path()}")
            return
        return BucketFolder(self.gbucket, self.join(dest))

    def get_blob(self, path):
        path = self.join(path)
        return self.gbucket.get_blob(path)

    def join(self, suf):
        normed = normpath(join(self.path, suf))
        if normed == '.': return ''
        return normed + ('/' if suf.endswith('/') or suf == '' else '')
    
    def list_files(self, prefix='', delimiter='/', glob=''):
        prefix = self.join(prefix)
        blobs = self.gbucket.list_blobs(prefix=prefix, delimiter=delimiter, match_glob=glob)
        files = set(relpath(b.name, prefix) for b in blobs)
        folders = set(relpath(b, prefix) for b in blobs.prefixes)
        return files, folders
    
    def download_files(self, srce_files, srce_pref='', dest_pref='.'):
        """
          This method downloads the specified files,
          including the folder tree relative to `srce_pref`,
          rebasing the folders at `dest_pref`.
        """

        srce_pref = self.join(srce_pref)
        gctm.download_many_to_path(self.gbucket, srce_files, dest_pref, srce_pref)

    def download_glob(self, srce_glob, srce_pref='', dest_pref='.'):
        """
          This method searches beneath `srce_pref` for `srce_glob`.
          Any matched files, including their folder tree relative to 'src_pref',
          are downloaded to `dest_pref`.
        """
        # both of these methods adjust the prefix, so don't here
        files, _ = self.list_files(srce_pref, glob=srce_glob)
        self.download_files(files, srce_pref, dest_pref)


class Workspace:

    def __init__(self, name, gproject='same'):
        self.name = name
        try:
            wsdata = fc_workspaces.loc[name]
        except KeyError as e:
            print(f"Workspace {name} not found!", file=sys.stderr)
            raise KeyError(f"Workspace {name} not found!") from e

        self.project = wsdata.namespace
        self.bucket = wsdata.bucketName

        self.load_atts()
        self.folder = None
        if gproject:
            if gproject == 'same': gproject = self.project
            gclient = clients[gproject]
            self.folder = BucketFolder(gclient.get_bucket(self.bucket), '')
    
    def load_atts(self):
        self.last_request = fcl.get_workspace(self.project, self.name, fields=('workspace, workspace.attributes'))
        if self.check_request():
            return
        
        data = self.last_request.json()
        self.attr_table = tabulate_fcattrs([data['workspace']])
        keys = self.attr_table.columns
        self.reference_keys = set(k for k in keys if k.startswith("referenceData_"))
        self.system_keys = set(k for k in keys if k.startswith("system:"))
        self.data_keys = (set(keys) - self.reference_keys) - self.system_keys

        data['workspace'].pop('attributes')
        data['workspace'].pop('name')
        self.metadata = data
    
    def setGProject(self, gproject):
        gclient = clients[gproject]
        self.folder = BucketFolder(gclient.get_bucket(self.bucket), '')
        
    def check_request(self):
        return check_request(self.last_request)

    # Firecloud
    ## Tables and Entities

    def list_tables(self):
        self.last_request = fcl.list_entity_types(self.project, self.name)
        if self.check_request(): return None
        tables = self.last_request.json()
        return list(tables.keys())

    def get_table(self, tab):
        self.last_request = fcl.get_entities(self.project, self.name, tab)
        if self.check_request(): return None
        entities = self.last_request.json()
        
        table = tabulate_fcattrs(entities)
        table.set_index('name', inplace=True)
        table.index.name = 'id'
        return table
    
    def update_entity(self, tab, ent, **attr_val):
        req = []
        for attr, val in attr_val.items():
            if type(val) == ReferenceList:
                val = reflist(val.entity_list, val.entity_type)

            if type(val) == list:
                val = attlist(val)

            if type(val) == JSONEntry:
                val = val.value
        
            req.append({
                "op": "AddUpdateAttribute",
                "attributeName": attr,
                "addUpdateAttribute": val
            })
        self.last_request = fcl.update_entity(self.project, self.name, tab, ent, req)
        return self.check_request()


    def get_workspace_data(self):
        return self.attr_table[list(self.data_keys)]
    
    def get_system_data(self):
        return self.attr_table[list(self.system_keys)]
    
    def get_reference_data(self):
        return self.attr_table[list(self.reference_keys)]


    ## Submissions
    def list_submissions(self, fields='brief', fmap=None):
        self.last_request = fcl.list_submissions(self.project, self.name)
        if self.check_request():
            return

        if fields == 'brief':
            fields = [
                'submissionId',
                'methodConfigurationName',
                'submissionDate',
                'status',
                'submissionRoot',
                'submissionEntity.entityName',
                'submissionEntity.entityType']
            if not fmap:
                fmap = lambda x: key_ep(x).removeprefix('submission').lower()
        
        table = tabulate(self.last_request.json(), fields)
        table.set_index('submissionId', inplace=True)
        if 'submissionDate' in table:
            table['submissionDate'] = pd.to_datetime(table['submissionDate'])
            table.sort_values(by="submissionDate", ascending=False, inplace=True)
        if 'submissionRoot' in table:
            pre = join('gs://', self.bucket)
            table['submissionRoot'].apply(lambda d: relpath(d, pre) + '/')
        if fmap:
            table.rename(columns=fmap, inplace=True)
        return table
    
    def get_submission(self, sid, fields='brief', workflows=True):
        return Submission(self, sid)


class Submission():

    def __init__(self, ws, id):
        self.ws = ws
        self.id = id

        req = fcl.get_submission(ws.project, ws.name, self.id)
        if check_request(req):
            raise KeyError(f'Submission {id} not found in workspace.')
        
        raw_data = req.json()
        raw_data.pop('submissionId')

        self.config = dict()
        for k in ['memoryRetryMultiplier',
                  'ignoreEmptyOutputs',
                  'useCallCache',
                  'deleteIntermediateOutputFiles',
                  'useReferenceDisks']:
            self.config[k] = raw_data.pop(k)
        
        wfd = raw_data.pop('workflows')
        
        self.data = dict()
        self.data['date'] = pd.to_datetime(raw_data.pop('submissionDate'))
        self.data['path'] = relpath(raw_data.pop('submissionRoot'), ws.folder.cloud_path())
        self.data |= {key_ep(k): navkey(raw_data, k) for k in flatten(raw_data)}
          
        wft = tabulate(wfd,
                       fields = ['workflowId', 'status', 'cost',
                                 'workflowEntity.entityName',
                                 'workflowEntity.entityType', 'messages'],
                       fmap = key_ep)
        wft.set_index('workflowId', inplace=True)

        self.workflows = wft

        self.folder = None
        if ws.folder is not None:
            self.folder = ws.folder.getdir(self.data['path'])

    def get_workflow(self, wid):
        return Workflow(self.ws, self, wid)
    
    def __getattr__(self, name):
        if name in self.data:
            return self.data[name]
        raise AttributeError(f'{name} not in submission details')
    
class Workflow():

    def __init__(self, ws, sub, id):
        self.ws = ws
        self.sub = sub
        self.id = id

        raw_data = self.get_metadata(exclude_key=['calls', 'workflowProcessingEvents',
                                             'inputs', 'outputs', 'labels',
                                             'submittedFiles'])
        
        raw_data.pop('calls', None) # because it still likes to include an empty dict
        raw_data.pop('id')

        self.data = dict()
        self.data['queue'] = pd.to_datetime(raw_data.pop('submission'))
        if 'start' in raw_data: self.data['start'] = pd.to_datetime(raw_data.pop('start'))
        if 'end' in raw_data: self.data['end'] = pd.to_datetime(raw_data.pop('end'))
        self.data['path'] = relpath(raw_data.pop('workflowRoot'), sub.folder.cloud_path())
        self.data |= raw_data

        self.folder = None
        if sub.folder is not None:
            self.folder = sub.folder.getdir(self.data['path'])

    def get_metadata(self, include_key=None, exclude_key=None):
        if type(include_key) == str:
            include_key = [include_key]
        
        self.last_request = fcl.get_workflow_metadata(
            self.ws.project, self.ws.name,
            self.sub.id, self.id,
            include_key, exclude_key)
        
        if check_request(self.last_request):
            raise KeyError(f'Workflow {id} not found in submission.')
        
        return self.last_request.json()
    
    def __getattr__(self, name):
        if name in self.data:
            return self.data[name]
        raise AttributeError(f'{name} not in submission details')


#######################################################
### Methods for analyzing and converting firecloud json

def check_request(req):
    if req.status_code != 200:
        print('Bad Request:',
                req.status_code,
                req.reason,
                req.content)
        return True
    return False

class JSONEntry():
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return f'JSONEntry({self.value})'
    
    def __repr__(self):
        return f'JSONEntry({self.value})'
    
class ReferenceList():
    def __init__(self, et, el):
        self.entity_type = et
        self.entity_list = el

    def __str__(self):
        return f'ReferenceList({self.entity_type}, {self.entity_list})'
    
    def __repr__(self):
        return f'ReferenceList({self.entity_type}, {self.entity_list})'

def tabulate_fcattrs(dlist: list[dict], fields=None, fmap=None, delim='.'):
    '''
      First eliminates itemsType nesting for lists and entity references
      while wrapping json entries in a benign object.
      then tabulates normally.
    '''
    flatter = []
    for i, listing in enumerate(dlist):
        new_listing = {}
        new_listing['name'] = listing.get('name', i)

        atts = listing['attributes']
        for k in atts:
            v = atts[k]
            if type(v) == dict:
                if 'items' in v: # proper lists and entity references
                    v = v['items'] # now v is a list
                    if len(v) == 0:
                        v = None
                    elif type(v[0]) == dict: # entity references
                        v = ReferenceList(v[0]['entityType'], [d['entityName'] for d in v])
                else: # a json object
                    v = JSONEntry(v)
            elif type(v) == list: # also json object
                v = JSONEntry(v)
            new_listing[k] = v
        flatter.append(new_listing)
    return tabulate(flatter, fields, fmap, delim, 0)

# list of dicts to table
def tabulate(dlist: list[dict], fields=None, fmap=None, delim='.', max_depth=None):
    '''
      Converts a list of nested dictionaries into a table.
        Each top-level dictionary in the list becomes a row.
        The nesting is flattened by concatenating child keys recursively.
        All keys and values are retained by default,
        with np.nan for missing values in rows.
      
      fields: list of keys to retain,
      fmap: dict to change field names in post
      delim: key concat delimiter
      max_depth: levels to unnest
    '''

    if fields is None:
        fields = agg_keys(dlist, delim=delim, max_depth=max_depth)
    
    by_field = ddict(list)
    for d in dlist:
        for f in fields:
            by_field[f].append(navkey(d, f, delim=delim))

    df = pd.DataFrame(by_field)
    if fmap: df.rename(columns=fmap, inplace=True)
    return df

def key_ep(k, delim='.'):
    return k.split(delim)[-1]

def navkey(d, k, delim='.'):
    a, _, b = k.partition(delim)
    if a not in d: return np.nan
    if not b: return d[a]
    return navkey(d[a], b, delim=delim)

def agg_keys(dlist, delim='.', max_depth=None):
    keys = set()
    for d in dlist:
        keys |= flatten(d, delim=delim, max_depth=max_depth)
    return keys
    
def flatten(d, prefix="", delim='.', _depth=0, max_depth=None):
    if (max_depth is not None) and _depth >= max_depth:
        return set(prefix+k for k in d)
    
    keys = set()
    for k in d:
        if isinstance(d[k], Mapping):
            keys |= flatten(d[k], prefix+k+delim, delim, _depth+1, max_depth)
        else:
            keys.add(prefix + k)
    return keys

def get_types(ndict):
    bytype = ddict(list)
    for k in flatten(ndict):
        v = navkey(ndict, k)
        w = type(v) if v != np.nan else None
        bytype[w].append(k)

    return bytype

# value/reference list to DataModel json
def attlist(l):
    return {'itemsType': 'AttributeValue',
            'items': l}

def reflist(l, t='sample'):
    return {'itemsType': 'EntityReference',
            'items': [{'entityType': t,
                       'entityName': i} for i in l] }


### Initialization
clients = ClientCache()
last_request = None
fc_workspaces = None
list_workspaces()


# environment variables set in terra notebook environments
#   if these exist, hook up to the notebook for convenience
nbWorkspace = None
if 'WORKSPACE_NAMESPACE' in os.environ and \
   'WORKSPACE_NAME' in os.environ:
    nbWorkspace = Workspace(os.environ['WORKSPACE_NAME'])
    print('kterra:')
    print('  Discovered Environemnt Workspace:', nbWorkspace.name)
    print("  Exposed in module variable 'nbWorkspace'")

