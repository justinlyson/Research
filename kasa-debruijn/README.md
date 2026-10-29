# Hamiltonian packings with prescribed symmetry

A constructor and independent verifier for Justin Lyson's revised de Bruijn graph construction. Given an alphabet of size `q` and vertex-word length `k`, the constructor exports `q-1` arc-disjoint Hamiltonian cycles related by the prescribed alphabet permutation: **0 stays fixed; 1 through q-1 cycle**.

Paper: [Hamiltonian packings with prescribed symmetry in de Bruijn graphs](https://zenodo.org/records/23717476).

## Setup

Use Python **3.10 or later**, preferably 64-bit for large jobs. Everything uses the standard library: no pip dependencies, compiler or solver. Extract the complete folder, open a terminal there, and keep its Python modules together. Windows multiprocessing uses spawned processes; the supplied commands handle this automatically.

## Quick start

```sh
python generate.py --q 4 --k 5 --out outputs/q4_k5
python verify.py outputs/q4_k5 --q 4 --k 5 --report outputs/q4_k5_check.json
```

Choose a new output directory for each generation run. Existing directories are protected. Generation stages files privately and publishes the complete output by an atomic directory rename on Windows. Verification can also read an existing compatible export; it does not import the constructor or trust its metadata.

Larger examples, using four workers:

```sh
python generate.py --q 30 --k 5 --out outputs/q30_k5 --workers 4
python verify.py outputs/q30_k5 --q 30 --k 5 --workers 4 --report outputs/q30_k5_check.json

python generate.py --q 31 --k 5 --out outputs/q31_k5 --workers 4
python verify.py outputs/q31_k5 --q 31 --k 5 --workers 4 --report outputs/q31_k5_check.json
```

The defaults permit these **k=5** examples by preflight arithmetic. These particular large runs are not claimed as completed benchmarks. Plan disk space before starting: q30,k5 uses about 1.88 GB for exported text and 2.82 GB for verifier arc bitsets; q31,k5 uses about 2.30 GB and 3.44 GB respectively. Temporary references and reports need additional space. Allow generous headroom, such as 8 GiB of free space when outputs and temporary data share a drive. Actual memory and runtime depend on the machine and graph.

## Live progress and timing

The parent process displays stage progress on **stderr**, including measured counts, percentages, elapsed time and rate/estimated remaining time when meaningful. Construction identifies its memory level, raw lift, actual-cycle labelling and repair. Export shows completed phase/chunk counts. Verification identifies per-cycle checks, the **GLOBAL collision merge**, and residual analysis. Worker processes do not print competing progress lines. Updates are throttled to limit overhead. A stage reaches 100% only when it finishes; verification announces PASS only after the global checks and residual analysis finish.

The final generation summary reports total wall time, construction time, export preparation/writing time and overhead. `generation.json` also contains measured construction sub-stages, worker settings and file hashes. Total wall time includes worker startup, worker shutdown, temporary cleanup and output commit; the final console/receipt formatting is outside that interval. Verification records its file-check, global-merge, residual and total times in the optional report.

Use `--quiet` to suppress live progress. Use `--json` for a complete machine-readable stdout result; progress remains on stderr. For example:

```sh
python generate.py --q 4 --k 5 --out outputs/quiet_run --quiet --json
python verify.py outputs/quiet_run --q 4 --k 5 --quiet --json
```

## Workers and resource settings

`--workers 0` selects automatic mode: serial for small jobs and up to four processes for large jobs. `--workers 1` is the deterministic serial fallback. Explicit worker counts are limited by available processors, useful tasks and the memory estimate. More processes can make small jobs slower.

| Option | Default | Applies to |
|---|---:|---|
| `--workers` | `0` (automatic) | Both |
| `--memory-limit-mib` | `8192` (8 GiB) | Both |
| `--max-states` | `100000000` | Both |
| `--max-arcs` | `10000000000` | Verifier |
| `--max-temp-mib` | `16384` (16 GiB) | Verifier arc bitsets |
| `--temp-dir` | System temporary directory | Both |
| `--quiet` | Off | Both |

These are **caps, not preallocated RAM or reserved disk**. Memory uses a conservative sum of process working-set peaks; shared mapped pages may be counted more than once. The verifier's temporary budget covers arc bitsets; its reference, final exports and reports require additional disk space. Error messages state required/allowed sizes and the flag to raise. You can raise resource caps explicitly. The genuine representation limits remain `q <= 255` and `q^k <= 4294967295`.

To use another drive, first create a directory there and pass it to both commands with `--temp-dir`. Generation uses it for raw-lift mappings; export staging stays beside `--out` to preserve the atomic commit. The verifier uses it for arc bitsets and its reference. Overrides can be combined, for example:

```sh
python verify.py outputs/q30_k5 --q 30 --k 5 --workers 4 --max-states 100000000 --max-arcs 10000000000 --memory-limit-mib 8192 --max-temp-mib 16384
```

The induction between memories, actual cycle labels and repair-tree choices stay sequential and deterministic. Large raw lifts and phase exports use processes and shared file mappings. Verification workers check whole files; the parent explicitly checks their global arc union. Large arrays are not serialized into every worker. Worker failures and Ctrl+C stop the operation and join owned workers. An interrupted operation emits no new PASS result; any already committed export is retained and can be verified separately.

## Output and verification

Each `cycle_NNN.txt` contains exactly `q^k` decimal integer symbols separated by whitespace. Line breaks are for readability. The sequence is **circular**: after its last symbol, continue at its first. No extra closing symbols are appended. All phases start at the all-zero length-k vertex. File phase j is the prescribed permutation applied j times to phase zero.

The verifier independently checks every file's alphabet and exact length, every cyclic length-k word exactly once, closing arcs, every cross-family length-(k+1) arc, and prescribed symmetry separately. `--no-symmetry` disables only the prescribed-action test; all Hamiltonicity and global arc-disjointness checks remain active. The same checks apply with one or several workers.

The residual graph is derived from actual omitted arcs. It has one incoming and one outgoing arc at each vertex, including the q forced constant-word loops. Its other components need not form one cycle. A full verification report includes their lengths and a lexicographically smallest representative word for every component.

## Mathematical scope

The construction implements the revised prescribed-symmetry theorem for `q > 2`, `k >= 2`, subject to representation and resource limits. Here `k` is the length of a graph vertex word. It includes the separate ternary construction, the finite small-alphabet initializers and the current general even/odd inductions. It does not use generic Hamilton-cycle search or an older construction as a fallback.

The `k=1` boundary is not implemented. In particular, three arc-disjoint Hamilton cycles at q4,k1 are impossible. An ordinary nonsymmetric packing at q20,k1 is a separate valid problem; it is not supplied or claimed as passed here. The binary ordinary-packing case is also outside this runtime's scope.

The mathematical paper and its separately distributed formal companion have their own verification scope. This Python runtime is not claimed to have been Lean-verified.

## License

Copyright 2026 Justin Lyson. Original runtime software is licensed under the Apache License, Version 2.0; see [LICENSE](LICENSE) and [NOTICE](NOTICE). The manuscript is separately licensed under CC BY 4.0 and is not included here.
