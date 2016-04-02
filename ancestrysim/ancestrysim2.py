import sys
import random
from collections import namedtuple, defaultdict, Counter
from operator import attrgetter
import click
import numpy as np
from numpy import random as rand
import itertools

import _ancestrysim2 as _anc

SegmentPair = namedtuple('SegmentPair', ('mum', 'dad'))
Sim = namedtuple('Sim', ('ngens', 'only_x', 'tree'))

SEGTYPES_NBLOCKS = {'auto': 'auto_nblocks', 'x': 'x_nblocks'}
SEGTYPES_LENS = {'auto': 'auto_blocklens', 'x': 'x_blocklens'}
SEGTYPES_HASGEN = {'auto': 'has_asegments', 'x': 'has_xsegments'}
LENS = {'x': 1.96, 'auto': 35.239999999999995}


def fib(k):
    phi = (1 + np.sqrt(5))/2.0
    psi = (1 - np.sqrt(5))/2.0
    return int((phi**k - psi**k)/np.sqrt(5))


def prob_xancestor_female(k, sex=0):
    """
    The formula here onle works for a present-day female individual that's
    female!
    """
    pf = (fib(k+int(not sex))/float(fib(k+2)))
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
        probs = (prob_xancestor_female(gen), prob_xancestor_male(gen))
        return probs
    else:
        raise NotImplementedError("mode must be 'auto' or 'x'")


def print_inds(inds):
    for ind in inds:
        print(ind)


def width(seg):
    return seg[1] - seg[0]


def segwidths(segs):
    return [(s[1] - s[0]) for s in segs]


def get_segs(ind, type='auto'):
    segfuns = {'x': attrgetter('xsegments'), 'auto': attrgetter('asegments')}
    segs = segfuns[type]
    return segs(ind)


def gen_material(ind, type='auto'):
    def gm(ind):
        return (sum(segwidths(get_segs(ind, type)[0])) +
                sum(segwidths(get_segs(ind, type)[1])))
    if isinstance(ind, list):
        return sum([gm(i) for i in ind])
    elif isinstance(ind, Individual):
        return gm(ind)
    else:
        raise ValueError("ind must be list or Individual")


def nrec(ind, type='auto'):
    if type == 'auto':
        return ind.anrec
    elif type == 'x':
        return ind.xnrec
    else:
        raise ValueError("type must be 'auto' or 'x'")


def find_overlaps(seglist1, seglist2):
    """Find all overlaps of segments

    >>> a = [(0, 4), (5, 9)]
    >>> b = [(3, 5), (8, 13)]
    >>> find_overlaps(a, b)
    [(3, 4), (8, 9)]
    """
    overlaps = list()
    for seg in seglist1:
        overlaps.extend(intersection(seg, seglist2))
    return overlaps


def overlaps(seg, other):
    """Returns whether `other` overlaps this segment.
    >>> a = (3, 4)
    >>> b = (4, 10)
    >>> c = (0, 3)
    >>> d = (6, 12)
    >>> e = (0, 1)
    >>> f = (0, 10)
    >>> overlaps(f, e)
    True
    >>> overlaps(b, a)
    False
    >>> overlaps(a, b)
    False
    >>> overlaps(c, a)
    False
    >>> overlaps(d, b)
    True
    >>> overlaps(b, d)
    True
    """
    olaps = other[0] < seg[1] and seg[0] < other[1]
    return olaps


def intersection(seg, segs):
    """Returns the overlapping region as a Segment().
    This is defined as the region of maximal overlap.

    >>> a = (0, 10)
    >>> b = (3, 4)
    >>> c = (0, 1)
    >>> intersection(a, b)
    (3, 4)
    >>> intersection(a, c)
    (0, 1)
    >>> intersection(a, [b, c])
    [(3, 4), (0, 1)]
    """
    single_range = not isinstance(segs, list)
    if single_range:
        segs = [segs]
    olaps = []
    for range in segs:
        if not overlaps(seg, range):
            if single_range:
                return None
            else:
                continue
        sg = (max(seg[0], range[0]), min(seg[1], range[1]))
        olaps.append(sg)
    if single_range:
        return olaps[0]
    return olaps


def ancestry2dot(tree):
    """
    Convert an ancestry tree to dot format.
    """
    raw_nodes = list()
    for gen_nodes in tree.values():
        raw_nodes.extend(gen_nodes)
    node_dict = dict([(x.id, x) for x in raw_nodes])
    nodes = list()
    edges = list()
    ranks = defaultdict(list)
    for node in raw_nodes:
        sex = 'M' if node.sex else 'F'
        nid = '    %s_%d_%d_%d' % (sex, node.id, node.gen, node.xnrec)
        style = "shape=square" if node.sex else 'shape=circle'
        if node.is_x:
            style += ', style=filled, fillcolor="#DDDDDD"'
        if node.is_x and node.child != -1:
            label = 'label="%d"' % node.xnrec
        else:
            label = 'label=""'
        vals = (label, style)
        nstr = nid + ((' [%s fontname="Latin Modern", fontsize=20, ' +
                       'penwidth=1, %s, color="#888888"]') % vals)
        nodes.append(nstr)
        if node.child != -1:
            child = node_dict[node.child]
            csex = 'M' if child.sex else 'F'
            evals = (sex, node.id, node.gen, node.xnrec, csex, child.id,
                     child.gen, child.xnrec)
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


def get_fullsib_overlaps(mums, dads, type='auto'):
    """Return the segments that overlap in mom and dad.

    Two backwards genealogies give us two instatiations of the ancestral
    genealogical process. In the fullsib case, this gives us two moms and two
    dads, that we compute the overlaps for. The number of shared blocks is the
    sum of these two groups of shared blocks.

    Two individuals are IBD iff they share an overlap on the
    haplotype, for either paternal haplotype.
    """
    mum_shared = mums[0].ibd_segments(mums[1], type)
    dad_shared = dads[0].ibd_segments(dads[1], type)
    return list(mum_shared) + list(dad_shared)


def get_halfsib_overlaps(shared_parents, type='auto'):
    """Return the segments that overlap in mom and dad.
    """
    segs = shared_parents[0].ibd_segments(shared_parents[1], type)
    return list(segs)


class Individual(object):
    __slots__ = ['id', 'sex', 'csex', 'gen', 'child', 'anrec', 'asegments',
                 'xsegments', 'xnrec', 'is_x']

    def __init__(self, id=0, sex=0, csex=None, gen=0, child=None, anrec=0,
                 asegments=None, xsegments=None, xnrec=None, is_x=False):
        # Note: it's assumed that the order of these arguments is exact same
        # as from _anc._simtree()
        self.id = id
        self.gen = gen
        self.sex = sex
        self.csex = csex
        self.child = child
        self.anrec = anrec
        # TODO this is slow, could be done in C.
        # self.asegments = ([Segment(*s) for s in asegments[0]],
        #                   [Segment(*s) for s in asegments[1]])
        # self.xsegments = ([Segment(*s) for s in xsegments[0]],
        #                   [Segment(*s) for s in xsegments[1]])
        self.asegments = asegments
        self.xsegments = xsegments
        self.xnrec = xnrec
        self.is_x = is_x

    def __repr__(self):
        slotstrs = ["%s=%s" % (s, getattr(self, s)) for s in self.__slots__]
        return "Individual(%s)" % ", ".join(slotstrs)

    def to_json(self, type='auto'):
        all_slots = dict([(s, getattr(self, s)) for s in self.__slots__])
        return dict(**all_slots)

    @property
    def auto_blocks(self):
        mum, dad = self.asegments
        return mum + dad

    @property
    def x_blocks(self):
        mum, dad = self.xsegments
        if self.sex == 1:
            return mum
        else:
            return mum + dad

    @property
    def auto_nblocks(self):
        mum, dad = self.asegments
        return [len(mum), len(dad)]

    @property
    def auto_blocklens(self):
        mum, dad = self.asegments
        return segwidths(mum) + segwidths(dad)

    @property
    def x_blocklens(self):
        mum, dad = self.xsegments
        if self.sex == 1:
            return segwidths(mum)
        else:
            return segwidths(mum) + segwidths(dad)

    @property
    def x_nblocks(self):
        mum, dad = self.xsegments
        if self.sex == 1:
            return [len(mum)]
        else:
            return [len(mum), len(dad)]

    @property
    def has_xsegments(self):
        return sum(self.x_nblocks) > 0

    @property
    def has_asegments(self):
        return sum(self.auto_nblocks) > 0

    def ibd_segments(self, other_ind, type='auto'):
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
        mum = find_overlaps(get_segs(self, type)[0],
                            get_segs(other_ind, type)[0])
        dad = find_overlaps(get_segs(self, type)[1],
                            get_segs(other_ind, type)[1])
        return list(mum) + list(dad)


def set_seed(seed):
    _anc._set_seed(int(seed))


def simtree(ngens, xchrom=True, autochroms=False, sex=0, only_x=True,
            random_seed=False):
    """
    Call C routine to generate genealogical trees.

    random_seed=True should only be used for testing.
    """
    if random_seed:
        seed = int(100000*random.random())
        set_seed(seed)
    ng, ox, sim = _anc._simtree(int(ngens), int(xchrom),
                                int(autochroms), int(sex), int(only_x))
    assert(ng == ngens)
    assert(ox == only_x)
    nsim = defaultdict(list)
    for gen, inds in enumerate(sim):
        for j, ind in enumerate(inds):
            nsim[gen].append(Individual(*sim[gen][j]))
    return nsim


def results2dataframe(d, keycol, include_rms=False):
    """Transform a Counter() object to a Pandas DataFrame
    """
    import pandas as pd
    if keycol in ('gen_anc', 'nblocks'):
        pd_dict = {'gen': list(), keycol: list(), 'count': list()}
        if include_rms:
            pd_dict['rm'] = list()
        max_nblocks = max([max(v.keys()) for v in d.values()])
        for genkey, counts in d.iteritems():
            for key in range(max_nblocks+1):
                if include_rms:
                    gen, rm = genkey
                else:
                    gen = genkey
                count = counts[key]
                pd_dict['gen'].append(gen)
                if include_rms:
                    pd_dict['rm'].append(rm)
                pd_dict[keycol].append(key)
                pd_dict['count'].append(count)
        return pd.DataFrame(pd_dict)
    if keycol == 'lens':
        pd_dict = {'gen': list(), keycol: list()}
        for genkey, lengths in d.iteritems():
            if include_rms:
                gen, rm = genkey
            else:
                gen = genkey
            pd_dict['gen'].extend([gen]*len(lengths))
            if include_rms:
                pd_dict['rm'].extend([rm]*len(lengths))
            pd_dict[keycol].extend(lengths)
        return pd.DataFrame(pd_dict)


def runsims(gens, nsims, autochroms=False, xchrom=True, only_x=True,
            type='auto', store_sims=False, gen_ancestors=False, verbose=True,
            asserts=True):
    sims = None
    if store_sims:
        sims = defaultdict(list)
    lens = defaultdict(list)
    nblocks = defaultdict(Counter)
    genetic_ancestors = defaultdict(Counter)
    segtype_nblocks = SEGTYPES_NBLOCKS[type]
    segtype_lens = SEGTYPES_LENS[type]
    segtype_hasgen = SEGTYPES_HASGEN[type]
    for gen in gens:
        for i in xrange(nsims):
            sim = simtree(gen, autochroms=autochroms, xchrom=xchrom,
                          only_x=only_x)
            if asserts:
                assert(sim[gen][0].gen == gen)
            if store_sims:
                sims[gen].append(sim)
            blocks = itertools.chain(*[getattr(i, segtype_nblocks) for i in
                                       sim[gen]])
            nblocks[gen] += Counter(blocks)
            blocklens = list()
            for ind in sim[gen]:
                blocklens.extend(getattr(ind, segtype_lens))
                if gen_ancestors:
                    has_gen = getattr(ind, segtype_hasgen)
                    genetic_ancestors[gen][int(has_gen)] += 1
            lens[gen].extend(blocklens)
            if asserts and gens[-1] == gen:
                # check we don't lose any material
                assert(abs(sum(blocklens) - LENS[type]*2) < 0.001)

        if verbose:
            sys.stderr.write("completed generation %d\n" % gen)

    out = dict(nblocks=nblocks, lens=lens)
    if gen_ancestors:
        out['gen_anc'] = genetic_ancestors
    if store_sims:
        out['sims'] = sims
    return out


def sample_ind(tree, sex=None, child_sex=None, gen=None, nrec=None, size=None):
    """Sample an individual from a gen (default last gen) and
    optionally a specific sex"""
    if gen is None:
        gen = max(tree.keys())  # last gen
    inds = tree[gen]
    if sex is not None:
        inds = filter(lambda x: x.sex == sex, inds)
    if child_sex is not None:
        inds = filter(lambda x: x.csex == child_sex, inds)
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


def runcousinsims(gens, nsims, sib_type, mode='auto',
                  cousins_sex=(0, 0), shared=None,
                  include_rms=False, verbose=True):
    """
    """
    autochroms, xchrom = {'auto': (True, False), 'x': (False, True)}[mode]
    cousins_sex = (int(cousins_sex[0]), int(cousins_sex[1]))
    # results from summary functions
    nblocks = defaultdict(Counter)
    lens = defaultdict(list)
    for gen in gens:
        for i in range(nsims):
            sims = list()
            for i in range(2):
                sims.append(simtree(gen, xchrom=xchrom, autochroms=autochroms,
                                    sex=cousins_sex[i], only_x=(mode == 'x'),
                                    # seed set elsewhere
                                    random_seed=False))
            # process the two trees depending on type of sib relationship
            if sib_type == 'hs':
                # run the two sims:
                ssex = shared  # shared sex
                if shared is None:
                    # if shared sex is None, sample random sex
                    # for this simulation
                    probs = probs_ancestor_sex(gen, mode=mode)
                    ssex = int(rand.choice((0, 1), size=1, p=probs))
                    # sample the shared ancestor from sim individuals
                    shared_parents = (sample_ind(sims[0], sex=ssex, gen=gen),
                                      sample_ind(sims[1], sex=ssex, gen=gen))
                    overlaps = get_halfsib_overlaps(shared_parents, mode)
                    if include_rms:
                        rms = sum(nrec(i, mode) for i in shared_parents)
                        nblocks[(gen, rms)][len(overlaps)] += 1
                    else:
                        nblocks[gen][len(overlaps)] += 1
                    lens[gen].extend(map(width, overlaps))
            elif sib_type == 'fs':
                # sample two mothers and fathers from each tree
                # with full sibs, the siblings must be female (male shared
                # ancestor), and thus we need to sample mums with only female
                # children.
                mums = (sample_ind(sims[0], sex=0, child_sex=0, gen=gen),
                        sample_ind(sims[1], sex=0, child_sex=0, gen=gen))
                dads = (sample_ind(sims[0], sex=1, child_sex=0, gen=gen),
                        sample_ind(sims[1], sex=1, child_sex=0, gen=gen))
                overlaps = get_fullsib_overlaps(mums, dads, mode)
                nblocks[gen][len(overlaps)] += 1
                lens[gen].extend(map(width, overlaps))
            else:
                raise ValueError("sib_type must be 'hs' or 'fs'")
        if verbose:
            sys.stderr.write("completed generation %d\n" % gen)

    return dict(nblocks=nblocks, lens=lens)


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


def write_rows(blocks, file=sys.stdout):
    blocks.to_csv(file, sep='\t', encoding='utf-8', index=False)


@click.group()
def cli():
    pass


@click.command(help="run ancestry simulations")
@click.option('--ngens', default='1:9',
              help="number of generations to run simulation for (default: 9)")
@click.option('--nsims', default=1000,
              help="number of simulations to run (default: 1000)")
@click.option('--auto', 'mode', flag_value='auto',
              default=True,
              help="cousins through autosomal chromosome genealogy")
@click.option('--x', 'mode', flag_value='x',
              default=False, help="cousins through X chromosome genealogy")
@click.option('--gen-anc', is_flag=True,
              default=False, help="count genetic ancestors")
@click.option('--seed', default=None, type=click.INT,
              help="the seed to use for the PRNG")
@click.option('--len', 'out', flag_value='lens',
              default=False,
              help="output segment lengths")
@click.option('--num', 'out', flag_value='nblocks',
              default=True, help="output segment numbers")
@click.option('--all-ancestors', is_flag=True,
              default=False,
              help="whether to include all genealogical individuals")
@click.option('--verbose', '-v', is_flag=True, default=True,
              help="print status messages during simulations")
def ancestor(ngens, nsims, seed, mode, out, gen_anc, all_ancestors, verbose):
    if seed is None:
        seed = int(100000*random.random())
    set_seed(seed)
    rand.seed(seed)
    only_x = not all_ancestors and mode == 'x'
    res = runsims(gens=parse_gens(ngens),
                  nsims=nsims, autochroms=(mode == 'auto'),
                  xchrom=(mode == 'x'), only_x=only_x,
                  type=mode, gen_ancestors=gen_anc,
                  verbose=verbose)
    if gen_anc:
        out = 'gen_anc'
    write_rows(results2dataframe(res[out], out))


@click.command(help="run cousin simulations")
@click.option('--full-sibs', 'sib_type', flag_value='fs',
              default=False, help="full sib cousin")
@click.option('--half-sibs', 'sib_type', flag_value='hs',
              default=True, help="half sib cousin")
@click.option('--shared-random', 'shared', is_flag=True, flag_value='random',
              default=True,
              help=("for half-sibs, whether the shared parent is " +
                    "randomly drawn"))
@click.option('--shared-mother', 'shared', is_flag=True,
              flag_value='mother', default=False,
              help="for half-sibs, whether the shared parent is a mother")
@click.option('--shared-father', 'shared', is_flag=True,
              flag_value='father', default=False,
              help="for half-sibs, whether the shared parent is a father")
@click.option('--auto', 'mode', flag_value='auto',
              default=True,
              help="cousins through autosomal chromosome genealogy")
@click.option('--rms', 'include_rms', is_flag=True,
              default=False,
              help="include total recombination meioses")
@click.option('--x', 'mode', flag_value='x',
              default=False, help="cousins through X chromosome genealogy")
@click.option('--len', 'out', flag_value='lens',
              default=False,
              help="output segment lengths")
@click.option('--num', 'out', flag_value='nblocks',
              default=True, help="output segment numbers")
@click.option('--seed', default=None, type=click.INT,
              help="the seed to use for the PRNG")
@click.option('--ngens', default='1:9', type=str,
              help="number of generations to run simulation for (default: 9)")
@click.option('--sexes', default='f,f', type=str,
              help="sex of cousins in format 'm,f' (default: f,f)")
@click.option('--nsims', default=1000,
              help="number of simulations to run (default: 800)")
@click.option('--verbose', '-v', is_flag=True, default=True,
              help="print status messages during simulations")
def cousin(sib_type, shared, mode, include_rms, ngens, out, seed, sexes,
           nsims, verbose):
    """Simulate X chromosome cousins.
    """
    if seed is None:
        seed = int(100000*random.random())
    set_seed(seed)
    rand.seed(seed)
    gens = parse_gens(ngens)
    shared = dict(random=None, mother=0, father=1)[shared]
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
    sims = runcousinsims(gens=gens, nsims=nsims, sib_type=sib_type,
                         mode=mode, cousins_sex=parse_cousin_sex(sexes),
                         shared=shared, include_rms=include_rms,
                         verbose=verbose)
    write_rows(results2dataframe(sims[out], out, include_rms=include_rms))


@click.command(help="run the doctests")
@click.option('--verbose', '-v', is_flag=True, default=False,
              help="verbose doctesting")
def test(verbose):
    import doctest
    doctest.testmod(verbose=verbose)


@click.command(help="generate a DOT file for genealogy")
@click.option('--ngens', default=5, type=int,
              help="number of generations back")
def dot(ngens):
    seed = int(100000*random.random())
    set_seed(seed)
    tree = simtree(ngens, True, True, 0, False, False)
    print(ancestry2dot(tree))


cli.add_command(test)
cli.add_command(cousin)
cli.add_command(ancestor)
cli.add_command(dot)

if __name__ == "__main__":
    cli()
    # import ancestralsim as anc
    # nsims, ngens = 500, 8,
    # if sys.argv[1] == 'old':
    #     for _ in range(nsims):
    #         s = anc.sim_ancestry(ngens, 'auto')
    # elif sys.argv[1] == 'new':
    #     for _ in range(nsims):
    #         s = simtree(ngens, autochroms=True, only_x=False)
    #         # print(gen_material(s[ngens-1]))
