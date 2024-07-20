def demBurnin(params):
  import msprime
  ts = msprime.sim_ancestry(samples=params['samples'],recombination_rate=params['rec_rate'],sequence_length=params['seq_length'],population_size=params['pop_size'])
  ts = msprime.sim_mutations(ts, rate=params['mut_rate'])
  f = open(params['ms_path'] + '.vcf', 'w')
  ts.write_vcf(f,position_transform='legacy')
  f.close()

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
    ts = pyslim.recapitate(ts, ancestral_Ne=x['Ne_recap'], recombination_rate=0)
    
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
  return ts

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

