# Faster exact minimum cut

This release speeds up VieCut's exact minimum cut (`mincut_parallel <g> exact`
and `mincut <g> noi`) without changing a single computed cut. Every change below
was accepted only after confirming it leaves the minimum cut of every test
instance unchanged.

## Results

Over 53 instances (33 generated families from the 1997 Chekuri-Goldberg-Karger-
Levine-Stein study, 20 real-world k-cores of large web and social graphs), on a
64-core AMD EPYC:

| configuration     | before  | after   | speedup |
|-------------------|---------|---------|---------|
| exact, 1 thread   | 25.0 s  | 2.80 s  | 8.9x    |
| exact, 32 threads | 5.76 s  | 0.65 s  | 8.8x    |

Geometric mean of the per-instance running time. Against independent
reimplementations of the two classical competitors (Nagamochi-Ibaraki and
Hao-Orlin, as shipped with the 1997 study), the tuned exact algorithm is about
five times faster at 32 threads over the instances those codes can represent,
and faster in the mean even on a single thread.

## What changed, and why

**The priority queue is sized by the cut, not the graph.** The parallel exact
path always used a bucket queue with one bucket per attainable key, a
`std::deque` sized by the value of the cut. On a graph with a large minimum cut
that is hundreds of gigabytes of empty buckets per round. The code now selects
between that queue, a heap, and an intrusive-list queue by the cut value, as the
sequential `noi_minimum_cut` already did.

**The graph is contracted in place, not rebuilt per round.** Both exact loops
rebuilt the whole graph each round and retained every intermediate. On families
that contract one vertex per round this was quadratic in time and memory; it is
now a single copy, overwritten each round.

**The contraction is parallel where the round collapses the graph, and dense.**
`contractGraphSparse` partitioned by original vertex, which forces a concurrent
hash map and a serialized build of the result. A second path partitions by
contracted vertex, needs neither, and sums parallel edges in a dense array that
fits in cache when the contraction leaves few vertices.

**The pass count is the thread count.** The parallel capforest runs independent
traversals, one per thread. A feedback rule that reduced the count drove the
quantity it measured and collapsed to one pass; it is removed. Whether the
passes share their visited array is decided from the graph size before the first
round.

**The pre-run's reductions are made to cascade.** The Padberg-Rinaldi tests only
fire once contraction has raised the degrees they compare against. The batch
versions measure a frozen graph and so find almost nothing; a new `pr_cascade`
contracts on the spot and re-reads degrees, which on the cycle families removes
the whole graph in a fraction of a second where the capforest needed a round per
vertex. It runs on sparse graphs only, retires when unproductive, and all its
merges ride the round's existing contraction.

**Smaller edge record, lighter inner loop, cheaper contraction primitive.** The
edge record drops two fields the exact path never reads, under `COMPACT_EDGES`,
for the parallel binaries that cannot reach the flow algorithms. The traversal
scans the adjacency as one range, uses byte flags, and resolves a vertex's queue
position once per update. `contractEdge` and `contractVertexSet` find shared
neighbourhoods through stamp arrays rather than a hash map built per merge.

**The pre-run is skipped where the minimum degree is already the cut**, on small
sparse graphs, where it was pure overhead.

## New

`graph2bin` converts a METIS text graph to a binary form that loads without
parsing; `graph_io` reads it transparently.

## Compatibility

No public interface changed. `exact`, `noi`, and the output format are as
before; every change is internal and every computed cut is identical.
