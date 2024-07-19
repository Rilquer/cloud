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
