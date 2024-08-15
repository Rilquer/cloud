## Creating ancestral trees with msprime
def eDem_anctree(params):
  import msprime
  ts = msprime.sim_ancestry(samples=params['samples'],recombination_rate=params['rec_rate'],sequence_length=params['seq_length'],population_size=params['pop_size'])
  f = open(params['path'] + '.tre', 'w')
  ts.dump(f,position_transform='legacy')
  f.close()

## Make ancestral pop by merging independent trees into one tree
def eDem_ancPop(seq_len):
  import tskit
  import numpy as np
  anc = tskit.TableCollection(sequence_length=seq_len)
  edge_table = tables.edges
edge_table.set_columns(
    left=np.array([0.0, 0.0, 0.0, 0.0]),
    right=np.array([1e3, 1e3, 1e3, 1e3]),
    parent=np.array([2, 2, 4, 4], dtype=np.int32),  # References IDs in the node table
    child=np.array([0, 1, 2, 3], dtype=np.int32),  # References IDs in the node table
)
edge_table


def demBurninRep(params):
  import msprime
  import concurrent
  with concurrent.futures.ProcessPoolExecutor() as executor:
    executor.map(demBurnin, params['ms_path'],params['samples'],params['rec_rate'],params['seq_length'],params['pop_size'],params['mut_rate'])




####  Recapitation - reticulate   ####
def recap(x):
  import msprime, tskit, pyslim
  words = [x['outpath'],x['outfile'],'.trees']
  treefile = "".join(words)
  ts = tskit.load(treefile)
  # Recapitate
  if x['ev_type']=='none':
    rts = pyslim.recapitate(ts, ancestral_Ne=x['Ne_recap'], recombination_rate=0)
    
  if x['ev_type']=='split':
    demography = msprime.Demography.from_tree_sequence(ts)
    for pop in demography.populations:
    # must set their effective population sizes
      pop.initial_size = x['Ne_recap']
    demography.add_migration_rate_change(time=ts.metadata['SLiM']['tick'],rate=x['mig_rate'], source="p1", dest="p2")
    demography.add_migration_rate_change(time=ts.metadata['SLiM']['tick'],rate=x['mig_rate'], source="p2", dest="p1")
    rts = pyslim.recapitate(ts, demography=demography,recombination_rate=0)
  
  # Adding mutations with default JukesCantor matrix
  rts = msprime.sim_mutations(rts,rate=x['mut_rate'])
  return rts

def piCalc(rts):
  return rts.diversity()

def tajDCalc(rts):
  return rts.Tajimas_D()

def segCalc(rts):
  return rts.segregating_sites()

def afsCalc(rts):
  return rts.allele_frequency_spectrum()

# Function below from:
# https://github.com/tskit-dev/tskit/issues/504
def alleCount(rts, sample_sets=None):
  if sample_sets is None:
    sample_sets = [rts.samples()]
  def f(x):
    return x
  return rts.sample_count_stat(sample_sets, f, len(sample_sets), windows='sites', polarised=True, mode='site', strict=False, span_normalise=False)

