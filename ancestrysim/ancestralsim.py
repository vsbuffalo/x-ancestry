from __future__ import print_function
from numpy import random as rand
import numpy as np
from scipy.stats import poisson, binom
from collections import namedtuple, deque, Counter, defaultdict
import operator
import sys
import json
import click
import random


DEBUG = True
EPS = 1.e-6

width = operator.attrgetter('width')


def centiMorgans2Morgans(x):
    return [l/100.0 for l in x]


XGENLEN = 1.96
AUTO_GENLENS = centiMorgans2Morgans([179.0, 160.0, 173.0, 127.0, 117.0, 131.0,
                                     135.0, 129.0, 120.0, 108.0, 278.0, 108.0,
                                     62.0, 74.0, 263.0, 225.0, 213.0, 204.0,
                                     193.0, 187.0, 170.0, 168.0])
AUTO_GENLENS_DICT = dict((i, l) for i, l in enumerate(AUTO_GENLENS))

AUTO_GENLEN = sum(AUTO_GENLENS)

MAX_CHROM_LEN = max(AUTO_GENLENS[:] + [XGENLEN])

Recombination = namedtuple("Recombination",
                           ('segments', 'states', 'nbreaks', 'breaks'))

SegmentPair = namedtuple('SegmentPair', ('mum', 'dad'))


def genlen(mode='auto'):
    return {'auto': AUTO_GENLEN, 'x': XGENLEN}[mode]


def _default(self, obj):
    # an extremely elegant monkey-patch
    # http://bit.ly/1UbGiKS
    return getattr(obj.__class__, "to_json", _default.default)(obj)


_default.default = json.JSONEncoder().default
json.JSONEncoder.default = _default


def fib(k):
    phi = (1 + np.sqrt(5))/2.0
    psi = (1 - np.sqrt(5))/2.0
    return int((phi**k - psi**k)/np.sqrt(5))


def prob_xancestor_female(k, sex=0):
    """
    The formula here onle works for a present-day female individual that's
    female!
    """
    assert(sex == 0)
    pf = (fib(k+1)/float(fib(k+2)))
    pm = 1 - pf
    return pf**2/(pf**2 + pm**2)


def prob_xancestor_male(k, sex=0):
    return 1 - prob_xancestor_female(k, sex)


def probs_ancestor_sex(gen, mode='auto'):
    """Returns probability of female, male ancestor.
    """
    if mode == 'auto':
        return (0.5, 0.5)
    elif mode == 'x':
        return (prob_xancestor_female(gen),
                prob_xancestor_male(gen))
    else:
        raise NotImplementedError("mode must be 'auto' or 'x'")


def find_overlaps(seglist1, seglist2):
    """Find all overlaps of segments

    >>> a = [Segment(0, 4), Segment(5, 9)]
    >>> b = [Segment(3, 5), Segment(8, 13)]
    >>> find_overlaps(a, b)
    [(3.0, 4.0], (8.0, 9.0]]
    """
    overlaps = list()
    for seg in seglist1:
        overlaps.extend(seg.intersection(seglist2))
    return overlaps


class Individual(object):
    __slots__ = ['id', 'gen', 'sex', 'child', 'nrec', 'segments', 'lineage',
                 'xsegments', 'xnrec', 'is_X']

    def __init__(self, segments, id=0, gen=0, sex=0, child=None,
                 nrec=0, lineage=None, xnrec=None, xsegments=None, is_X=None):
        self.id = id
        self.gen = gen
        self.sex = sex
        self.child = child
        self.nrec = nrec
        self.segments = segments
        self.xsegments = xsegments
        self.lineage = lineage
        self.xnrec = xnrec
        self.is_X = is_X

    def __repr__(self):
        slotstrs = ["%s=%s" % (s, getattr(self, s)) for s in self.__slots__]
        return "Individual(%s)" % ", ".join(slotstrs)

    def json_prob_anc(self):
        """return the probability this individual is an X ancestor"""
        breakpoints = np.arange(1e3)  # how many breakpoints there are
        prob = 1 - np.sum(poisson.pmf(breakpoints, (self.xnrec-1)*XGENLEN) *
                          binom.pmf(0, n=breakpoints+1, p=1/2.0**(self.xnrec-1)))
        if np.isnan(prob):
            prob = 1
        all_slots = dict([(s, getattr(self, s)) for s in self.__slots__])
        return dict(prob_anc=prob, **all_slots)

    def to_json(self):
        segs = dict(mum=map(lambda x: x.to_json(), self.segments.mum),
                    dad=map(lambda x: x.to_json(), self.segments.dad))
        all_slots = dict([(s, getattr(self, s)) for s in self.__slots__
                          if s != 'segments'])
        return dict(segments=segs, **all_slots)

    def ibd_segments(self, other_ind):
        """Return the segments that are overlapping across
        two individuals.

        Note: we search both maternal/paternal segments for IBD overlap.
        In the case of a male individual and searching for X chromosome
        segments in Individual.segments.dad, we are assuming that the
        simulation code is either only tracing X *or* autosome ancestry
        through Individual.segments. If both need to be traced, both
        Individual.segements and Individual.xsegments need to be used.

        In other words, if X chrom segments are in Individual.segments,
        the guarantee that there are no X segments in a male's
        Individual.segements.dad must come from simulation, and is not
        ensured here.
        """
        mum = find_overlaps(self.segments.mum, other_ind.segments.mum)
        dad = find_overlaps(self.segments.dad, other_ind.segments.dad)
        return list(mum) + list(dad)


def block_counts_nrec2pd(x):
    """
    """
    import pandas as pd
    max_blocks = max([max(v.keys()) for v in x.values()])
    block_rng = range(0, max_blocks+1)
    dt = defaultdict(list)
    for (gen, nrec1, nrec2), blocks in x.iteritems():
        dt['gen'].append(gen)
        dt['nrec1'].append(nrec1)
        dt['nrec2'].append(nrec2)
        for nblocks in block_rng:
            dt[nblocks].append(blocks.get(nblocks, 0))
    cols = ['gen', 'nrec1', 'nrec2'] + block_rng
    return pd.DataFrame(dt, columns=cols)


def auto_haplotypes():
    """Create a set of segments corresponding to autosomes.
    """
    segs = list()
    last_chrom_len = 0
    for chrom_len in AUTO_GENLENS:
        end = last_chrom_len + chrom_len
        segs.append(Segment(last_chrom_len, end))
        last_chrom_len = end
    return segs


def auto_segments():
    """Create segments for mother and father autosome haplotypes.
    """
    return SegmentPair(auto_haplotypes(), auto_haplotypes())


def single_auto_segments(genlen):
    mum, dad = ([Segment(0, genlen)], [Segment(0, genlen)])
    return SegmentPair(mum, dad)


def seglengths(ind):
    mum = [s.width for s in ind.segments.mum]
    dad = [s.width for s in ind.segments.dad]
    return sum(mum) + sum(dad)


def x_segments(sex=0):
    mum, dad = ([Segment(0, XGENLEN)], [Segment(0, XGENLEN)])
    if sex:
        # males only have maternal X
        dad = []
    return SegmentPair(mum, dad)


class Segment(object):
    def __init__(self, start, end):
        if (end < start):
            raise ValueError('error: end < start; start <= end')
        self.start = float(start)
        self.end = float(end)

    @property
    def width(self):
        return self.end - self.start

    def contains(self, position):
        """Returns whether the segment contains the position
        >>> a = Segment(3, 4)
        >>> a.contains(3)
        True
        >>> a.contains(4)
        False
        >>> b = Segment(3, 10)
        >>> b.contains(5)
        True
        """
        return position >= self.start and position < self.end

    def __repr__(self):
        return "(%s, %s]" % (round(self.start, 4), round(self.end, 4))

    def overlaps(self, other):
        """Returns whether `other` overlaps this segment.
        >>> a = Segment(3, 4)
        >>> b = Segment(4, 10)
        >>> c = Segment(0, 3)
        >>> d = Segment(6, 12)
        >>> b.overlaps(a)
        False
        >>> a.overlaps(b)
        False
        >>> c.overlaps(a)
        False
        >>> d.overlaps(b)
        True
        >>> b.overlaps(d)
        True
        """
        overlaps = other.start < self.end and self.start < other.end
        return overlaps

    def intersection(self, ranges):
        """Returns the overlapping region as a Segment().
        This is defined as the region of maximal overlap.

        >>> a = Segment(0, 10)
        >>> b = Segment(3, 4)
        >>> c = Segment(0, 1)
        >>> a.intersection(b)
        (3.0, 4.0]
        >>> a.intersection([b, c])
        [(3.0, 4.0], (0.0, 1.0]]
        """
        single_range = isinstance(ranges, Segment)
        if single_range:
            ranges = [ranges]
        overlaps = []
        for range in ranges:
            if not self.overlaps(range):
                if single_range:
                    return None
                else:
                    continue
            seg = Segment(max(self.start, range.start),
                          min(self.end, range.end))
            overlaps.append(seg)
        if single_range:
            return overlaps[0]
        return overlaps

    def to_json(self):
        return dict(start=self.start, end=self.end)


def recombine(segments, total_genlen):
    """Given a genetic length and some ancestral segments, simulate recombination.

    Input is a list (ordered by genetic position) of segments and a total
    genetic length of the chromosome.
    """
    if not len(segments):
        msg = "segments must be a list with more one or more elements"
        raise ValueError(msg)

    curr_state = int(random.choice((0, 1)))
    nbreaks = rand.poisson(total_genlen)
    if nbreaks == 0:
        # no recombination
        return Recombination(segments, [curr_state]*len(segments), 0, 0)
    rand_breakpoints = sorted(rand.uniform(0., total_genlen, nbreaks))

    breakpoints = deque(rand_breakpoints)
    new_segments = list()
    new_states = list()
    process_segments = deque(segments)  # segments to process

    # guaranteed at least one breakpoint and segment
    segment = process_segments.popleft()
    breakpoint = breakpoints.popleft()

    need_breakpoint = False
    while True:
        if need_breakpoint or breakpoint < segment.start:
            # eat break points, but switch state accordingly. Note that on the
            # ends, this sequence is non-ancestral so we could get away without
            # switching state.
            try:
                curr_state = int(not curr_state)
                breakpoint = breakpoints.popleft()
                need_breakpoint = False
            except IndexError:
                # no more breakpoints, so push all segments and the current
                # state to the results and exit
                new_segments.append(segment)
                new_states.append(curr_state)
                while len(process_segments):
                    new_segments.append(process_segments.popleft())
                    new_states.append(curr_state)
                break
        elif segment.contains(breakpoint):
            left = Segment(segment.start, breakpoint)
            # right segment is next segment
            segment = Segment(breakpoint, segment.end)
            new_segments.append(left)
            new_states.append(curr_state)
            need_breakpoint = True
        elif segment.end < breakpoint:
            # get a new segment and save old segment
            new_segments.append(segment)
            new_states.append(curr_state)
            try:
                segment = process_segments.popleft()
            except IndexError:
                while len(breakpoints):
                    breakpoints.popleft()
                # no more segments, but maybe some breakpoints
                # but these are all in non-ancestral sequence
                # so we can ignore.
                break
        else:
            raise Exception
    # a whole load of assert statements to check that everything looks good
    assert(len(process_segments) == len(breakpoints) == 0)
    assert(len(new_segments) == len(new_states))
    assert(len(new_segments) >= len(segments))  # all segments survive
    # now, make sure no segments are lost
    assert(abs(sum([s.width for s in new_segments]) -
           sum(s.width for s in segments)) < EPS)
    # rounded = map(lambda x: round(x, 4), rand_breakpoints)
    return Recombination(new_segments, new_states, nbreaks, rand_breakpoints)


def segregate(recomb):
    """Take a Recombination named tuple and return a SegmentPair() named tuple of
    the segments from the mother and father.
    """
    seg_states = zip(recomb.segments, recomb.states)
    mum = [sg for sg, st in seg_states if st == 0]
    dad = [sg for sg, st in seg_states if st == 1]
    assert(len(mum) + len(dad) == len(recomb.segments) ==
           len(recomb.states))
    return SegmentPair(mum, dad)


def meiosis(segments, genlen):
    """Simulate meiosis (in reverse)"""
    if not len(segments):
        return SegmentPair([], [])
    return segregate(recombine(segments, genlen))


def segregate2(recomb):
    """Take a Recombination named tuple and return a SegmentPair() named tuple of
    the segments from the mother and father.

    Unlike segregate(), this segregae2() function makes a further
    approximation; all blocks survive segregation independent of
    the others (the "coin-flipping" approximation).
    """
    seg_states = zip(recomb.segments, recomb.states)
    bern_states = rand.choice((0, 1), size=len(seg_states), replace=True)
    mum = [sg for i, (sg, st) in enumerate(seg_states) if not bern_states[i]]
    dad = [sg for i, (sg, st) in enumerate(seg_states) if bern_states[i]]
    assert(len(mum) + len(dad) == len(recomb.segments) ==
           len(bern_states))
    return SegmentPair(mum, dad)


def meiosis2(segments, genlen):
    """Simulate meiosis (in reverse)"""
    if not len(segments):
        return SegmentPair([], [])
    return segregate2(recombine(segments, genlen))


# Some general notes on the recursive ancestry tree
# creation functions below.
#
# All state must be passed through arguments
def x_recursion(seed, k, cid, tree, id, ngens, meiosfun):
    """The ancestral recursion function for the X chromosome
    """
    if k == 2:  # parent generation, insert lineage marker
        seed.lineage = seed.sex
    if k > ngens:
        return
    cid = id[0]
    id[0] += 1
    # in female, run X segments through meiosis
    # both_segments = seed.segments.mum + seed.segments.dad
    mum_segments = meiosfun(seed.segments.mum, XGENLEN)
    mum_nrec = seed.nrec+1
    mother = Individual(id=id[0], gen=k, sex=0, child=cid,
                        nrec=mum_nrec, segments=mum_segments,
                        lineage=seed.lineage, xnrec=mum_nrec, is_X=True)
    tree[k].append(mother)
    # both females and males have X inherited from mother
    x_recursion(mother, k+1, cid, tree, id, ngens, meiosfun)
    if seed.sex == 0:  # female
        # if female, also trace through father
        id[0] += 1
        # don't throw X segments in male through meiosis
        dad_segments = SegmentPair(list(seed.segments.dad), [])
        dad_nrec = seed.nrec
        father = Individual(id=id[0], gen=k, sex=1, child=cid,
                            nrec=dad_nrec, segments=dad_segments,
                            lineage=seed.lineage, xnrec=dad_nrec, is_X=True)
        tree[k].append(father)
        x_recursion(father, k+1, cid, tree, id, ngens, meiosfun)


def parent_is_X(seed, parent_sex):
    """
    Return True/False whether this seed is an X ancestor.
    """
    ind_sex, ind_is_X = seed.sex, seed.is_X
    if not ind_is_X:
        return False
    if ind_sex == 1 and parent_sex == 1:
        return False
    if ind_sex == 1 and parent_sex == 0:
        return True
    if ind_sex == 0:
        return True
    raise AssertionError("execution should not reach this point")


def parent_xnrec(seed, parent_sex):
    if seed.xnrec is None:
        return None
    if not parent_is_X(seed, parent_sex):
        return None
    # everything below is X ancestor
    if parent_sex == 0:
        return seed.xnrec+1
    else:
        return seed.xnrec
    raise AssertionError("execution should not reach this point")


def auto_recursion(seed, k, cid, tree, id, ngens, meiosfun):
    """The main recursion function for the autosomes"""
    if k == 2:  # parents, insert lineage marker
        seed.lineage = seed.sex
    if k > ngens:
        return
    cid = id[0]
    # mother
    id[0] += 1
    # 'mum' spelled properly for GC
    mum_segments = meiosfun(seed.segments.mum, AUTO_GENLEN)
    mother = Individual(id=id[0], gen=k, sex=0, child=cid,
                        nrec=seed.nrec+1, segments=mum_segments,
                        lineage=seed.lineage, is_X=parent_is_X(seed, 0),
                        xnrec=parent_xnrec(seed, 0))
    tree[k].append(mother)
    # both females and males have X inherited from mother
    auto_recursion(mother, k+1, cid, tree, id, ngens, meiosfun)
    # father
    id[0] += 1
    dad_segments = meiosfun(seed.segments.dad, AUTO_GENLEN)
    father = Individual(id=id[0], gen=k, sex=1, child=cid,
                        nrec=seed.nrec+1, segments=dad_segments,
                        lineage=seed.lineage, is_X=parent_is_X(seed, 1),
                        xnrec=parent_xnrec(seed, 1))
    tree[k].append(father)
    auto_recursion(father, k+1, cid, tree, id, ngens, meiosfun)


def both_recursion(seed, k, cid, tree, id, ngens, meiosfun):
    """The main recursion function for both the autosomes and X"""
    if k == 2:  # parents, insert lineage marker
        seed.lineage = seed.sex
    if k > ngens:
        return
    cid = id[0]
    # mother
    id[0] += 1
    mum_auto_segments = meiosfun(seed.segments.mum, AUTO_GENLEN)
    is_mum_x = parent_is_X(seed, 0)
    if is_mum_x:
        mum_x_segments = meiosfun(seed.xsegments.mum, XGENLEN)
    else:
        mum_x_segments = None
    mum_xnrec = is_mum_x if is_mum_x else None
    mother = Individual(id=id[0], gen=k, sex=0, child=cid,
                        nrec=seed.nrec+1, segments=mum_auto_segments,
                        lineage=seed.lineage, is_X=is_mum_x,
                        xsegments=mum_x_segments,
                        xnrec=mum_xnrec)
    tree[k].append(mother)
    # both females and males have X inherited from mother
    both_recursion(mother, k+1, cid, tree, id, ngens, meiosfun)
    # father
    id[0] += 1
    dad_auto_segments = meiosfun(seed.segments.dad, AUTO_GENLEN)
    is_dad_x = parent_is_X(seed, 1)
    if is_dad_x:
        dad_x_segments = SegmentPair(list(seed.xsegments.dad), [])
    else:
        dad_x_segments = None
    dad_xnrec = is_dad_x if is_dad_x else None
    father = Individual(id=id[0], gen=k, sex=1, child=cid,
                        nrec=seed.nrec+1, segments=dad_auto_segments,
                        lineage=seed.lineage, is_X=is_dad_x,
                        xsegments=dad_x_segments,
                        xnrec=dad_xnrec)
    tree[k].append(father)
    both_recursion(father, k+1, cid, tree, id, ngens, meiosfun)


def init_segments(mode, sex, genlen=None):
    """Initialize segments.
    genlen is for single, which is a single autosomal chromosome
    genelen Morgans long
    """
    initializers = {'auto': auto_segments,
                    'x': x_segments,
                    'single': single_auto_segments}
    if mode == 'single':  # optional genetic length allowed
        assert(genlen is not None)
        return initializers[mode](genlen)
    elif mode == 'x':
        return initializers[mode](sex)
    return initializers[mode]()


def sim_ancestry(ngens, mode='auto', genlen=None, seed=None,
                 meiosfun=meiosis):
    """Take a female seed individual and reconstruct their X ancestry ngens
    generations back.

    Focal individual's sex not a parameter, since this doesn't affect autosomal
    ancestry.
    """
    if (mode == 'single'):
        assert(genlen is not None)
    if mode != 'both':
        initial_segs = init_segments(mode, 0, genlen)
        if seed is None:
            seed = Individual(segments=initial_segs, xnrec=0, is_X=True)
    else:
        if seed is None:
            x_segs = init_segments('x', 0)
            auto_segs = init_segments('auto', 0)
            seed = Individual(segments=auto_segs, xsegments=x_segs, xnrec=0,
                              is_X=True)
    tree = defaultdict(list)
    # pop seed into tree
    tree[0] = [seed]
    id = [seed.id]

    # select the appropriate ancestry recurrence function
    anc_recur = {'auto': auto_recursion,
                 'x': x_recursion,
                 'both': both_recursion,
                 'single': auto_recursion}[mode]

    # not that we send k=1 generation since this k is set to parents gen
    # (parents gen back = 1)
    anc_recur(seed, 1, seed.child, tree, id, ngens, meiosfun)

    # check that we don't lose any material
    last_gen_seq = sum([seglengths(ind) for ind in tree[ngens-1]])
    first_gen_seq = sum([seglengths(ind) for ind in tree[0]])
    if abs(last_gen_seq - first_gen_seq) > EPS:
        msg = "ancestral sequence lost (1st gen: %0.2f, last gen: %0.2f)"
        raise AssertionError(msg % (first_gen_seq, last_gen_seq))

    return tree


def pair_x_with_auto_sims(xsim, autosim):
    """Given an X simulation and an autosome simulation, pair such that
    RMs are same across ancestors. This allows for comparison between X
    blocks and autosome blocks. This uses mutability of Individual().

    This treats X ancestors exchangeable once we condition on R.
    """
    max_gen = max(xsim.keys())
    assert(max_gen == max(autosim.keys()))
    rms = filter(lambda x: x is not None,
                 map(operator.attrgetter('xnrec'), autosim[max_gen]))
    assert(Counter(rms) == Counter(map(operator.attrgetter('nrec'),
                                       xsim[max_gen])))
    for rm in set(rms):
        # map xsim individuals for this value of RM to autosome
        # ancestors
        auto_inds = filter(lambda x: x.xnrec == rm, autosim[max_gen])
        x_inds = filter(lambda x: x.xnrec == rm, xsim[max_gen])
        x_inds = rand.permutation(x_inds)
        for i, ind in enumerate(auto_inds):
            # add X segments to autosomes
            ind.xsegments = x_inds[i].segments


def divide_segments(breakpoints):
    """
    Divided a set of segments into smaller segments based on
    breakpoints over entire length.
    """
    new_segments = list()
    breakpoints = deque(breakpoints)
    last_bpoint = breakpoints.popleft()
    while len(breakpoints):
        bpoint = breakpoints.popleft()
        new_segments.append(Segment(last_bpoint, bpoint))
        last_bpoint = bpoint
    return new_segments


def get_blocks(ind, mode='auto'):
    """Get blocks in Individual.segments attribute

    Note that for x genealogies, we only count the maternal haplotype in men,
    to avoid listing no blocks in a male's Individual.segments.dad (Y chrom)
    as a zero.
    """
    # gets number of blocks
    if mode == 'auto':
        return [len(ind.segments.dad), len(ind.segments.mum)]
    elif mode == 'x':
        # condition on sex so we don't count zero segments in
        # Individual.segments.dad (where Y chrom is) as zero X segments
        if ind.sex == 0:  # female
            if ind.segments is None:
                return [0, 0]
            return [len(ind.segments.dad), len(ind.segments.mum)]
        else:  # male
            # TODO BUG HERE
            if ind.segments is None:
                return [0]
            try:
                assert(len(ind.segments.dad) == 0)
            except AssertionError:
                import pdb; pdb.set_trace()
            return [len(ind.segments.mum)]
    raise AssertionError("execution should not reach this point")


def get_xblocks(ind, mode=None):
    """Get blocks in Individual.xsegments attribute.

    Note that for x genealogies, we only count the maternal haplotype in men,
    to avoid listing no blocks in a male's Individual.segments.dad (Y chrom)
    as a zero.

    mode is ignored, TODO in a refactor, should merge get_blocks and
    get_xblocks
    """
    # TODO: should be renamed to get_numblocks.
    # get blocks, not counting male's blocks
    if ind.sex == 0:
        if ind.xsegments is None:
            return [0, 0]
        return [len(ind.xsegments.dad), len(ind.xsegments.mum)]
    else:
        if ind.xsegments is None:
            return [0]
        return [len(ind.xsegments.mum)]


def get_blocks_across_gen(anc, gen, mode='auto', seg_mode='auto',
                          only_X=False, offset=True):
    """For an a defaultdict(list) containing the genealogy for certain number of
    generations, count the number of segment blocks back in a certain
    generation.

    Very important note: if we ask for the number of blocks in generation gen,
    this returns the number of blocks in generation gen - 1. This is because by
    convention we say a parent and offspring share 22 blocks, even though these
    are on different haplotypes in the parent. This is identical to offsetting
    by 1. We call this `offset`.
    """
    # segment mode, whether to access ind.segments or ind.xsegments
    blockfun = {'auto': get_blocks, 'x': get_xblocks}[seg_mode]
    # Below, note the minus one -- seet note above
    inds = anc[gen-int(offset)]
    if mode == 'both' and only_X:
        inds = filter(operator.attrgetter('is_X'), inds)
    nblocks = list()
    # below is some code for sampling an individual rather than
    # averaging across all individuals in a generation
    # i = int(rand.choice(range(len(inds))))
    # nblocks.extend(blockfun(inds[i]))
    # return nblocks
    for ind in inds:
        # use the blockfun to extend blocks
        nblocks.extend(blockfun(ind, mode))
    return nblocks


def runsims(gens, nsims, mode, genlen=None, store_sims=True,
            meiosfun=meiosis, seg_mode='auto', only_X=False,
            verbose=True):
    sims = defaultdict(list)
    block_counts = defaultdict(Counter)
    for gen in gens:
        # generation are 1-indexed!
        for i in range(nsims):
            anc = sim_ancestry(gen, mode, genlen=genlen, meiosfun=meiosfun)
            if store_sims:
                sims[gen].append(anc)
            block_counts[gen] += Counter(get_blocks_across_gen(anc, gen,
                                                               mode,
                                                               seg_mode,
                                                               only_X))
        if verbose:
            sys.stderr.write("completed generation %d\n" % gen)
    return sims, block_counts


def ancestry2dot(tree, html=True, fill=True):
    """
    Convert an ancestry tree to dot format.
    """
    # get all nodes
    raw_nodes = list()
    for gen_nodes in tree.values():
        raw_nodes.extend(gen_nodes)
    node_dict = dict([(x.id, x) for x in raw_nodes])
    nodes = list()
    edges = list()
    ranks = defaultdict(list)
    for node in raw_nodes:
        sex = 'M' if node.sex else 'F'
    nid = '    %s_%d_%d_%d' % (sex, node.id, node.gen, node.nrec)
    style = "shape=square" if node.sex else 'shape=circle'
    if fill and not node.sex:
        style += ', style=filled, fillcolor="#B2C2F0"'
    label = ""
    if html and fill:
        label = ('fontname="Latin Modern", fontsize=10, label=' +
                 '<<i>%s<sub>%d, %d</sub></i>>,' % (sex, node.gen, node.nrec))
    if not fill:
        label = 'label=""'
    vals = (sex, node.gen, node.nrec, label, style)
    nstr = nid + (' [texlbl="$%s_{%d, %d}$", %s fixedsize=true, penwidth=1, ' +
                  '%s, color="#888888"]' % vals)
    nodes.append(nstr)
    if node.child is not None:
        child = node_dict[node.child]
        csex = 'M' if child.sex else 'F'
        evals = (sex, node.id, node.gen, node.nrec, csex, child.id,
                 child.gen, child.nrec)
        estr = '    %s_%d_%d_%d -- %s_%d_%d_%d' % evals
        edges.append(estr)
    ranks[node.gen].append(nid)
    nodestr = "\n".join(nodes)
    edgestr = "\n".join(edges)
    # format ranks based on generation
    rank_nodes = list()
    for rnodes in ranks.values():
        rank_nodes.append("{rank=same; %s}" % "; ".join(rnodes))
    rankstr = "\n".join(rank_nodes)
    dotstr = "graph G {\n%s\n%s\n%s\n}" % (nodestr, edgestr, rankstr)
    return dotstr


def tree2json(tree):
    """Serialize an ancestry tree from one simulation into a JSON"""
    return json.dumps(tree)


def sims2json(sims, type, genlen, chromlens=None, prob_anc=False):
    """Serialize the sims from runsims()

    The data is reformatted to make it easier to visualize with javascript.
    """
    sims_out = {}
    max_gen = max(sims.keys())
    for sim in range(len(sims[max_gen])):
        tree = []
        for gen_inds in sims[max_gen][sim].values():
            gid = 0
            for ind in gen_inds:
                if prob_anc:
                    ind_obj = ind.json_prob_anc()
                else:
                    ind_obj = ind.to_json()
                ind_obj['gid'] = gid
                # gid += 1 + ind.sex
                gid += 1
                tree.append(ind_obj)
        sims_out[sim] = tree
    out = dict(sims=sims_out, type=type, genlen=genlen,
               chromlens=chromlens)
    return json.dumps(out)


def tab_delim_results(block_counts):
    """Output the simulation results as a tab-delimited file.
    """
    # don't output anything if interactive mode
    print("\t".join(["gen", "nblocks", "count"]))
    for gen, block_counts in block_counts.iteritems():
        for blocks, counts in block_counts.iteritems():
            print("\t".join(map(str, [gen, blocks, counts])))


def sample_ind(tree, sex=None, gen=None, nrec=None, size=None):
    """Sample an individual from a gen (default last gen) and
    optionally a specific sex"""
    if gen is None:
        gen = max(tree.keys())  # last gen
    inds = tree[gen]
    if sex is not None:
        inds = filter(lambda x: x.sex == sex, inds)
    if nrec is not None:
        inds = filter(lambda x: x.nrec == nrec, inds)
    try:
        sample = rand.choice(inds, size, replace=False)
    except ValueError:
        msg = ("no individuals to sample -- perhaps your configuration " +
               "is impossible (e.g. two male cousins one "
               + "generation apart cannot be X fullsibs).")
        raise ValueError(msg)
    return sample


def nrec(ind):
    return ind.nrec


def get_fullsib_overlaps(mums, dads):
    """Return the segments that overlap in mom and dad.

    Two backwards genealogies give us two instatiations of the ancestral
    genealogical process. In the fullsib case, this gives us two moms and two
    dads, that we compute the overlaps for. The number of shared blocks is the
    sum of these two groups of shared blocks.

    Two individuals are IBD iff they share an overlap on the
    haplotype, for either paternal haplotype.
    """
    mum_shared = mums[0].ibd_segments(mums[1])
    dad_shared = dads[0].ibd_segments(dads[1])
    return list(mum_shared) + list(dad_shared)


def get_halfsib_overlaps(shared_parents):
    """Return the segments that overlap in mom and dad.
    """
    segs = shared_parents[0].ibd_segments(shared_parents[1])
    return list(segs)


def runcousinsims(gens, nsims, sib_type, mode, cousins_sex=None,
                  shared=0, include_rms=False,
                  meiosfun=meiosis, verbose=True):
    """Like runsims() but runs two seperate ancestral genealogies,
    and chooses two (fullsib) or three (halfsib) individuals to
    be parents. Number of shared blocks are overlapping blocks.

    Returns a dict of block lengths, and a summary of block counts.

    TODO (low priority): refactor so just returns block add counting elsewhere.

    Beware: this code does not follow specific haplotypes, so may misbave
    with <2 generations.

    include_rms: whether to include how many recombinational meioses occur down
    each lineage.
    """
    if mode == 'single':
        msg = 'cousins with single chroms is not implemented'
        raise NotImplementedError(msg)

    if cousins_sex is not None and tuple(cousins_sex) != (0, 0):
        raise NotImplementedError('cousins_sex None or (0, 0) supported')
    # pre-allocate dict for dataframe
    block_lengths = defaultdict(list)
    block_counts = defaultdict(Counter)
    block_counts_nrec = None
    if include_rms:
        block_counts_nrec = defaultdict(Counter)
    female, male = 0, 1
    for gen in gens:
        for i in range(nsims):
            # Comparing relatedness of two individuals, we simulate two
            # genealogies.
            # The shared ancestor is chosen at random from the final gen.
            # Overlaps are shared blocks.
            # if gen > 3:
            #     import pdb; pdb.set_trace()  # XXX BREAKPOINT
            sims = list()
            seeds = (None, None)  # use default seed of female (None)
            if cousins_sex is not None:
                # We need to explicitly seed the sim ancestry routine
                # with individuals of the specified sex.
                # We also pass the sex to the initial segment generator,
                # as males and females start out with different initial
                # segments
                sex0, sex1 = cousins_sex
                seeds = (Individual(segments=init_segments(mode, sex0),
                                    sex=sex0),
                         Individual(segments=init_segments(mode, sex1),
                                    sex=sex1))
            for i in range(2):
                sims.append(sim_ancestry(gen, mode, seed=seeds[i],
                                         meiosfun=meiosfun))
            if sib_type == 'hs':
                ssex = shared
                if shared is None:
                    # random shared parent, conditioning on X ancestry
                    # TODO this *DOES NOT* handle different sexes for
                    # present day cousins
                    probs = probs_ancestor_sex(gen, mode)
                    ssex = rand.choice([female, male], size=1, p=probs)
                # the following assertion checks for a fence post error
                # damn fence post errors
                assert(max(sims[0].keys()) == gen)
                # sample two random parents (half-sib)
                shared_parent = (sample_ind(sims[0], sex=ssex, gen=gen),
                                 sample_ind(sims[1], sex=ssex, gen=gen))
                overlaps = get_halfsib_overlaps(shared_parent)
                block_counts[gen][len(overlaps)] += 1
                # also hatch number of recombinational meioses
                if include_rms:
                    nrec1, nrec2 = map(nrec, shared_parent)
                    block_counts_nrec[(gen, nrec1, nrec2)][len(overlaps)] += 1
                # incorporate all block lengths and their sim number
                block_lengths[gen].extend(map(width, overlaps))
            elif sib_type == 'fs':
                # not need to specify full-sib shared ancestory
                # due to symmetry of full-sibs
                # TODO:  FIXME, we need to pin on RMs too?
                mums = (sample_ind(sims[0], female),
                        sample_ind(sims[1], female))
                dads = (sample_ind(sims[0], male),
                        sample_ind(sims[1], male))
                overlaps = get_fullsib_overlaps(mums, dads)
                block_counts[gen][len(overlaps)] += 1
                if include_rms:
                    # also hatch number of recombinational meioses
                    # TODO: all we need to report is len(overlaps) += 1?
                    # TODO THIS IS WRONG FIXME
                    mum_nrec1, mum_nrec2 = map(nrec, mums)
                    dad_nrec1, dad_nrec2 = map(nrec, dads)
                    key = (gen, mum_nrec1, mum_nrec2, dad_nrec1, dad_nrec2)
                    block_counts_nrec[key][len(overlaps)] += 1
                # incorporate all block lengths and their sim number
                block_lengths[gen].extend(map(width, overlaps))
            else:
                msg = "sib_type not valid (either 'hs' or 'fs')"
                raise ValueError(msg)
        if verbose:
            sys.stderr.write("completed generation %d\n" % gen)

    return block_lengths, block_counts, block_counts_nrec


def parse_cousin_sex(arg):
    sexes = arg.split(',')
    if not all([s in 'mf' for s in sexes]):
        raise click.BadParameter("sex option must be in format 'm,f'")
    return tuple([dict(m=1, f=0)[s] for s in sexes])


def parse_gens(arg):
    if ':' in arg:
        start, end = map(int, arg.split(':'))
        return range(start, end+1)
    elif ',' in arg:
        return map(int, arg.split(','))
    else:
        return [int(arg)]


def sim_hscousins_coinflip(ngens, nsims, mode='auto', shared=None):
    """
    Simulation procedure for half-sib cousins following the
    "coin-flipping" approximation.
    """
    import pandas as pd
    # create one X genealogy that can be used to gather number RMs
    sims = dict(gen=[], nblocks=[], count=[])
    total_genlen = genlen(mode)
    for gen in ngens:
        # shared ancestor sex probabilities, which depend on generation
        nblocks = Counter()
        for nsim in range(nsims):
            # sample a sex of shared ancestor
            probs = probs_ancestor_sex(gen, mode)
            if shared is None:
                shared = int(rand.choice((0, 1), size=1, p=probs))
            # simulate a tree, only for nrec of shared ancestors
            # (not to use segments!)
            sim = sim_ancestry(gen, mode=mode)
            ind1, ind2 = (sample_ind(sim, sex=shared),
                          sample_ind(sim, sex=shared))
            # get number of recombinations in these two individuals
            r1, r2 = ind1.nrec, ind2.nrec
            nbreaks = rand.poisson(total_genlen*(r1+r2))
            # add additional breakpoints for chromosome ends
            rand_breakpoints = list(rand.uniform(0, total_genlen, nbreaks))
            rand_breakpoints.extend({'auto': AUTO_GENLENS + [0, AUTO_GENLEN],
                                     'x': (0, total_genlen)}[mode])
            breaks = sorted(rand_breakpoints)
            segments = divide_segments(breaks)
            # segregation surival probability, with
            # male ancestors have one haplotype:
            pseg = 1/2.0**(r1 + r2 - int(not shared))
            psegs = (1 - pseg, pseg)
            blocks = rand.choice((0, 1), size=len(segments),
                                 p=psegs, replace=True)
            nblocks[sum(blocks)] += 1

        for nblock, counts in nblocks.iteritems():
            sims['gen'].append(gen)
            sims['nblocks'].append(nblock)
            sims['count'].append(counts)
    df = pd.DataFrame(sims)
    df['prob'] = df.groupby(('gen'))['count'].transform(lambda x: x/sum(x))
    return df


@click.group()
def cli():
    pass


@click.command(help="run cousin simulations")
@click.option('--full-sibs', 'sib_type', flag_value='fs',
              default=True, help="full sib cousin")
@click.option('--half-sibs', 'sib_type', flag_value='hs',
              default=False, help="half sib cousin")
@click.option('--shared-mother', 'shared', flag_value=0, is_flag=True,
              help="for half-sibs, whether the shared parent is a mother")
@click.option('--shared-father', 'shared', flag_value=1, is_flag=True,
              help="for half-sibs, whether the shared parent is a father")
@click.option('--shared-random', 'shared', flag_value=None, is_flag=True,
              help=("for half-sibs, whether the shared parent is " +
                    "randomly drawn"))
@click.option('--auto', 'mode', flag_value='auto',
              default=True,
              help="cousins through autosomal chromosome genealogy")
@click.option('--x', 'mode', flag_value='x',
              default=False, help="cousins through X chromosome genealogy")
@click.option('--seed', default=None, type=click.INT,
              help="the seed to use for the PRNG")
@click.option('--ngens', default='1:9', type=str,
              help="number of generations to run simulation for (default: 9)")
@click.option('--sexes', default='f,f', type=str,
              help="sex of cousins in format 'm,f' (default: f,f)")
@click.option('--nsims', default=800,
              help="number of simulations to run (default: 800)")
@click.option('--coin-flip', '-c', is_flag=True, default=False,
              help='"enable coin-flipping" approximation to recombination')
@click.option('--verbose', '-v', is_flag=True, default=True,
              help="print status messages during simulations")
def cousin(sib_type, ngens, seed, nsims, shared, sexes, mode,
           coin_flip, verbose):
    """Simulate X chromosome cousins.
    """
    shared = None  # FIXME
    meiosfun = (meiosis, meiosis2)[int(coin_flip)]
    if seed is not None:
        rand.seed(seed)
    gens = parse_gens(ngens)
    if verbose:
        msg = (("genealogy: %s, cousin type: %s, shared parent: %s, " +
               "cousin sexes: %s\n") %
               (dict(auto='autosomal', x='X')[mode],
                dict(hs='half-sib', fs='full-sib')[sib_type],
                {0: 'mother', 1: 'father', None: 'random'}[shared],
                sexes))
        sys.stderr.write(msg)
    assert(shared in (0, 1, None))
    # don't store block_counts_nrec
    block_lengths, block_counts, _ = runcousinsims(gens, nsims, sib_type, mode,
                                                   parse_cousin_sex(sexes),
                                                   shared, include_rms=False,
                                                   meiosfun=meiosfun,
                                                   verbose=verbose)
    tab_delim_results(block_counts)


@click.command(help="run X chromosome simulations")
@click.option('--ngens', default='1:12',
              help="number of generations to run simulation for (default: 12)")
@click.option('--nsims', default=2000,
              help="number of simulations to run (default: 2000)")
@click.option('--seed', default=None, type=click.INT,
              help="the seed to use for the PRNG")
@click.option('--coin-flip', '-c', is_flag=True, default=False,
              help='"enable coin-flipping" approximation to recombination')
@click.option('--verbose', '-v', is_flag=True, default=True,
              help="print status messages during simulations")
@click.option('--json', is_flag=True, default=False,
              help="output simulations as JSON (default: block counts as tab" +
              " delimited")
@click.option('--json-prob-anc', is_flag=True, default=False,
              help="output the probability of ancestry")
def x(ngens, nsims, seed, coin_flip, verbose, json, json_prob_anc):
    if seed is not None:
        rand.seed(seed)

    meiosfun = (meiosis, meiosis2)[int(coin_flip)]
    if json_prob_anc:
        # run one sim, and return the probability of ancestry, which depends
        # only on RMs.
        sims, block_counts = runsims(gens=parse_gens(ngens),
                                     nsims=1, mode='x',
                                     store_sims=True,
                                     meiosfun=meiosfun,
                                     verbose=verbose)
        print(sims2json(sims, 'x', XGENLEN, dict(x=XGENLEN), prob_anc=True))
        return


    sims, block_counts = runsims(gens=parse_gens(ngens),
                                 nsims=nsims, mode='x',
                                 store_sims=json,
                                 meiosfun=meiosfun,
                                 verbose=verbose)
    if not json:
        tab_delim_results(block_counts)
    else:
        print(sims2json(sims, 'x', XGENLEN, dict(x=XGENLEN)))


@click.command(help="run autosomal chromosome simulations")
@click.option('--ngens', default='1:12',
              help="number of generations to run simulation for")
@click.option('--nsims', default=2000,
              help="number of simulations to run")
@click.option('--seed', default=None, type=click.INT,
              help="the seed to use for the PRNG")
@click.option('--coin-flip', '-c', is_flag=True, default=False,
              help='"enable coin-flipping" approximation to recombination')
@click.option('--verbose', '-v', is_flag=True, default=True,
              help="print status messages during simulations")
@click.option('--json', is_flag=True, default=False,
              help="output simulations as JSON (default: block counts as tab" +
              " delimited")
def auto(ngens, nsims, seed, coin_flip, verbose, json):
    if seed is not None:
        rand.seed(seed)
    sims, block_counts = runsims(gens=parse_gens(ngens),
                                 nsims=nsims,
                                 mode='auto', store_sims=json,
                                 verbose=verbose)
    if not json:
        tab_delim_results(block_counts)
    else:
        print(sims2json(sims, 'auto', AUTO_GENLEN, AUTO_GENLENS_DICT))


@click.command(help="run single autosomal chromosome simulations")
@click.option('--genlen', default=7.0, type=click.FLOAT,
              help="genetic length in Morgans")
@click.option('--ngens', default='1:12',
              help="number of generations to run simulation for")
@click.option('--nsims', default=2000,
              help="number of simulations to run")
@click.option('--verbose', '-v', is_flag=True, default=True,
              help="print status messages during simulations")
@click.option('--json', is_flag=True, default=False,
              help="output simulations as JSON (default: block counts as tab" +
              " delimited")
@click.option('--seed', default=None, type=click.INT,
              help="the seed to use for the PRNG")
def single(genlen, seed, ngens, nsims, verbose, json):
    if seed is not None:
        rand.seed(seed)
    sims, block_counts = runsims(gens=parse_gens(ngens),
                                 nsims=nsims,
                                 mode='single',
                                 genlen=genlen,
                                 store_sims=json,
                                 verbose=verbose)
    if not json:
        tab_delim_results(block_counts)
    else:
        print(sims2json(sims, 'single', genlen, 1))


@click.command(help="run the doctests")
@click.option('--verbose', '-v', is_flag=True, default=False,
              help="verbose doctesting")
def test(verbose):
    import doctest
    doctest.testmod(verbose=verbose)


cli.add_command(x)
cli.add_command(auto)
cli.add_command(single)
cli.add_command(cousin)
cli.add_command(test)


if __name__ == "__main__":
    cli()
