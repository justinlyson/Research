# Research

Code and data for Justin Lyson's work on Hamiltonian packings with prescribed symmetry in de Bruijn graphs.

Paper: [Hamiltonian packings with prescribed symmetry in de Bruijn graphs](https://zenodo.org/records/23717476).

## Contents

| Folder | What it holds |
|---|---|
| [`kasa-debruijn/`](kasa-debruijn/) | Constructor and independent verifier. For an alphabet of size `q` and vertex-word length `k`, it exports `q-1` arc-disjoint Hamiltonian cycles related by the prescribed alphabet permutation (0 stays fixed; 1 through q-1 cycle). Pure Python 3.10+, no dependencies. See its [README](kasa-debruijn/README.md) for full usage. |
| [`output_q31_k5/`](output_q31_k5/) | Exported cycles for `q = 31`, `k = 5`: 30 files `cycle_0000.txt` to `cycle_0029.txt`, each a circular sequence of `31^5` symbols. |

## Quick start

```sh
cd kasa-debruijn
python generate.py --q 4 --k 5 --out outputs/q4_k5
python verify.py outputs/q4_k5 --q 4 --k 5 --report outputs/q4_k5_check.json
```

To check the included q31 data:

```sh
cd kasa-debruijn
python verify.py ../output_q31_k5 --q 31 --k 5 --workers 4 --report outputs/q31_k5_check.json
```

Verifying q31, k5 needs roughly 3.4 GB of temporary disk space for arc bitsets; see the [resource notes](kasa-debruijn/README.md#workers-and-resource-settings).

## License

The code in `kasa-debruijn/` is released under the Apache License 2.0; see [`kasa-debruijn/LICENSE`](kasa-debruijn/LICENSE).
