# Function below from:
# https://github.com/tskit-dev/tskit/issues/504
def alleCount(rts, sample_sets=None):
    if sample_sets is None:
       sample_sets = [rts.samples()]
    def f(x):
       return x
    return rts.sample_count_stat(sample_sets, f, len(sample_sets), windows='sites', polarised=True, mode='site', strict=False, span_normalise=False)
