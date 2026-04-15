VieCut 1.00
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![C++](https://img.shields.io/badge/C++-17-blue.svg)](https://isocpp.org/)
[![CMake](https://img.shields.io/badge/CMake-3.9+-064F8C.svg)](https://cmake.org/)
[![Linux](https://img.shields.io/badge/Linux-supported-success.svg)](https://github.com/KaHIP/VieCut)
[![GitHub Stars](https://img.shields.io/github/stars/KaHIP/VieCut)](https://github.com/KaHIP/VieCut/stargazers)
[![GitHub Issues](https://img.shields.io/github/issues/KaHIP/VieCut)](https://github.com/KaHIP/VieCut/issues)
[![Last Commit](https://img.shields.io/github/last-commit/KaHIP/VieCut)](https://github.com/KaHIP/VieCut/commits)
[![JEA 2018](https://img.shields.io/badge/JEA'18-10.1145/3274662-blue)](https://doi.org/10.1145/3274662)
[![IPDPS 2019](https://img.shields.io/badge/IPDPS'19-10.1109/IPDPS.2019.00013-blue)](https://doi.org/10.1109/IPDPS.2019.00013)
[![ESA 2020](https://img.shields.io/badge/ESA'20-10.4230/LIPIcs.ESA.2020.59-blue)](https://doi.org/10.4230/LIPIcs.ESA.2020.59)
[![arXiv](https://img.shields.io/badge/arXiv-1708.06127-b31b1b.svg)](https://arxiv.org/abs/1708.06127)
[![arXiv](https://img.shields.io/badge/arXiv-1808.05458-b31b1b.svg)](https://arxiv.org/abs/1808.05458)
[![arXiv](https://img.shields.io/badge/arXiv-2002.06948-b31b1b.svg)](https://arxiv.org/abs/2002.06948)
[![Homebrew](https://img.shields.io/badge/homebrew-available-orange)](https://github.com/KaHIP/homebrew-kahip)
[![Heidelberg University](https://img.shields.io/badge/Heidelberg-University-c1002a)](https://www.uni-heidelberg.de)
=====

<p align="center">
  <img src="https://raw.githubusercontent.com/KaHIP/VieCut/master/logo/viecut-banner.png" alt="VieCut Logo" width="900"/>
</p>

**VieCut** is a library of shared-memory parallel algorithms for the minimum cut problem on undirected edge-weighted graphs. Part of the [KaHIP](https://github.com/KaHIP) organization.

> **Python Interface:** An easy-to-use Python interface for this software is available in [CHSZLabLib](https://github.com/CHSZLab/CHSZLabLib).

| | |
|:--|:--|
| **What it solves** | Find the minimum edge cut that separates a graph into two (or more) components |
| **Algorithms** | Inexact heuristic (VieCut), exact parallel (NOI-based), cactus representation, multiterminal cut |
| **Parallelism** | Shared-memory via OpenMP; sequential versions included |
| **Interfaces** | CLI |
| **Requires** | C++17 compiler, CMake 3.9+, OpenMP, MPI (multiterminal cut only) |

## Quick Start

### Install via Homebrew (Linux only)

```bash
brew install KaHIP/kahip/viecut
```

### Or build from source

```bash
git clone --recursive https://github.com/KaHIP/VieCut.git && cd VieCut
mkdir build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DUSE_TCMALLOC=OFF && make
```

### Run

```bash
# Heuristic minimum cut (fast, near-optimal)
./build/mincut network.graph vc

# Exact minimum cut (parallel)
./build/mincut_parallel network.graph exact

# Find all minimum cuts (cactus representation)
./build/mincut_parallel network.graph cactus

# Multiterminal cut (NP-hard, branch-and-reduce)
./build/multiterminal_cut network.graph -k 4
```

When installed via Homebrew, use `viecut_mincut`, `viecut_mincut_parallel`, etc. directly.

---

## Executables

| Binary | Parallel | Description |
|:-------|:---------|:------------|
| `mincut` | `mincut_parallel` | Minimum cut with choice of algorithm |
| `multiterminal_cut` | | Branch-and-reduce multiterminal cut (MPI) |
| `kcore` | `kcore_parallel` | k-core decomposition |
| `mincut_contract` | `mincut_contract_parallel` | Minimum cut with random edge contraction |
| `mincut_recursive` | `mincut_recursive_parallel` | Recursive minimum cut on largest SCC |

---

## Command Line Usage

```
./build/mincut [options] <graph-file> <algorithm>
./build/mincut_parallel [options] <graph-file> <algorithm>
```

### Algorithms

**Sequential:**

| Algorithm | Flag | Description |
|:----------|:-----|:------------|
| VieCut | `vc` | Fast heuristic, near-optimal in practice |
| Nagamochi-Ono-Ibaraki | `noi` | Classic exact algorithm |
| Karger-Stein | `ks` | Randomized contraction |
| Matula | `matula` | (2+e)-approximation |
| Padberg-Rinaldi | `pr` | Contraction rules |
| Cactus | `cactus` | Find all minimum cuts |

**Parallel:**

| Algorithm | Flag | Description |
|:----------|:-----|:------------|
| Inexact | `inexact` | Parallel VieCut heuristic |
| Exact | `exact` | Parallel exact minimum cut |
| Cactus | `cactus` | Parallel cactus of all minimum cuts |

### Options

| Option | Description | Default |
|:-------|:-----------|:--------|
| `-p <int>` | Number of processors (parallel only) | all |
| `-i <int>` | Number of iterations | `1` |
| `-q <type>` | Priority queue: `bqueue`, `bstack`, `heap` | `bqueue` |
| `-s` | Compute and save the cut | off |
| `-o <file>` | Output file for the cut (requires `-s`) | |
| `-b` | Find most balanced minimum cut (cactus only, requires `-s`) | off |
| `-l` | Disable priority queue limiting | off |

### Examples

```bash
# Parallel exact minimum cut with 8 threads, 5 iterations
./build/mincut_parallel -p 8 -i 5 network.graph exact

# Save cut to file
./build/mincut -s -o cut.txt network.graph noi

# Most balanced minimum cut
./build/mincut_parallel -s -b network.graph cactus
```

---

## Multiterminal Cut

The multiterminal cut separates a set of terminal vertices from each other with minimum total cut weight. For |T|=2 this equals the minimum s-t-cut; for |T|>2 it is NP-hard. VieCut uses a shared-memory parallel branch-and-reduce approach.

```
./build/multiterminal_cut <graph-file> [options]
```

| Option | Description |
|:-------|:-----------|
| `-f <file>` | Partition file assigning terminals |
| `-t <int>` | Add vertex as terminal (repeatable) |
| `-k <int>` | Use `k` highest-degree vertices as terminals |
| `-r <int>` | Use `r` random vertices as terminals |
| `-b <int>` | BFS expansion around terminals |
| `-p <int>` | Number of threads |

---

## Graph Format

VieCut uses the [METIS graph format](http://people.sc.fsu.edu/~jburkardt/data/metis_graph/metis_graph.html). Graphs can be unweighted or edge-weighted (format flag `1`).

```
<num_nodes> <num_edges> [format]
<neighbors of node 1> [weights]
<neighbors of node 2> [weights]
...
```

---

## Building from Source

### Requirements

- C++17 compiler (GCC 7+ or Clang 11+)
- CMake 3.9+
- OpenMP
- MPI (for multiterminal cut)
- TCMalloc (optional, for performance)

### Build

```bash
git clone --recursive https://github.com/KaHIP/VieCut.git && cd VieCut
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```

To build without TCMalloc:
```bash
cmake .. -DCMAKE_BUILD_TYPE=Release -DUSE_TCMALLOC=OFF
```

Sequential and parallel executables are placed in `build/`.

---

## Related Projects

| Project | Description |
|:--------|:------------|
| [KaHIP](https://github.com/KaHIP/KaHIP) | Karlsruhe High Quality Graph Partitioning |
| [fpt-max-cut](https://github.com/KaHIP/fpt-max-cut) | FPT-based maximum cut solvers |

---

## Licence

VieCut is free software provided under the MIT License.
If you publish results using our algorithms, please cite the applicable papers:

```bibtex
@article{DBLP:journals/jea/HenzingerNSS18,
  author    = {Monika Henzinger and Alexander Noe and Christian Schulz and Darren Strash},
  title     = {Practical Minimum Cut Algorithms},
  journal   = {{ACM} J. Exp. Algorithmics},
  volume    = {23},
  year      = {2018},
  doi       = {10.1145/3274662}
}

@inproceedings{DBLP:conf/ipps/HenzingerN019,
  author    = {Monika Henzinger and Alexander Noe and Christian Schulz},
  title     = {Shared-Memory Exact Minimum Cuts},
  booktitle = {2019 {IEEE} International Parallel and Distributed Processing Symposium ({IPDPS})},
  pages     = {13--22},
  publisher = {{IEEE}},
  year      = {2019},
  doi       = {10.1109/IPDPS.2019.00013}
}

@inproceedings{DBLP:conf/esa/HenzingerN0S20,
  author    = {Monika Henzinger and Alexander Noe and Christian Schulz and Darren Strash},
  title     = {Finding All Global Minimum Cuts in Practice},
  booktitle = {28th Annual European Symposium on Algorithms ({ESA} 2020)},
  series    = {LIPIcs},
  volume    = {173},
  pages     = {59:1--59:20},
  publisher = {Schloss Dagstuhl - Leibniz-Zentrum f{\"{u}}r Informatik},
  year      = {2020},
  doi       = {10.4230/LIPIcs.ESA.2020.59}
}

@inproceedings{DBLP:conf/alenex/HenzingerN022,
  author    = {Monika Henzinger and Alexander Noe and Christian Schulz},
  title     = {Practical Fully Dynamic Minimum Cut Algorithms},
  booktitle = {Proceedings of the 24th Symposium on Algorithm Engineering and Experiments ({ALENEX} 2022)},
  pages     = {13--26},
  publisher = {{SIAM}},
  year      = {2022},
  doi       = {10.1137/1.9781611977042.2}
}

@inproceedings{DBLP:conf/alenex/HenzingerN020,
  author    = {Monika Henzinger and Alexander Noe and Christian Schulz},
  title     = {Shared-Memory Branch-and-Reduce for Multiterminal Cuts},
  booktitle = {Proceedings of the 22nd Symposium on Algorithm Engineering and Experiments ({ALENEX} 2020)},
  pages     = {42--55},
  publisher = {{SIAM}},
  year      = {2020},
  doi       = {10.1137/1.9781611976007.4}
}

@inproceedings{DBLP:conf/acda/HenzingerN021,
  author    = {Monika Henzinger and Alexander Noe and Christian Schulz},
  title     = {Faster Parallel Multiterminal Cuts},
  booktitle = {Proceedings of the 2021 {SIAM} Conference on Applied and Computational Discrete Algorithms ({ACDA} 2021)},
  pages     = {100--110},
  publisher = {{SIAM}},
  year      = {2021},
  doi       = {10.1137/1.9781611976830.10}
}
```
