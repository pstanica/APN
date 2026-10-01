# Trace modifications of bivariate APN functions

Supplementary material for the paper by Pantelimon Stănică.
Uploaded this entire `trace_bivariate_apn` directory to
https://github.com/pstanica/APN. 
## Reproduction

Requires Python 3, numpy, and g++ supporting C++17. The symbolic audit
additionally requires python-flint. Run from this directory:

```sh
python compare_all.py
python audit_invariants.py
python audit_switching.py
python audit_cubic_extension.py
python audit_identities.py
```

All finite-family comparisons and the cubic audit use the `BinaryField`
implementation and construction formulas in `compare_all.py`.
The APN polar and ortho-derivative comparisons use `invn_canonical.cpp`.
The switching equations use `switchn.cpp` with the same field builders.
The cubic audit directly enumerates derivative fibers because its output is
nonquadratic; this verifies the actual APN condition rather than applying the
quadratic polar test outside its scope.

`comparison_manifest.json` records all 1,395 dimension-12 function tables:
exact parameters, formula source, field modulus, binary encoding, table SHA256,
algebraic degree, APN result, and complete invariant histogram.
There are 21 target tables, 1,334 APN comparator tables with different
histograms, and 40 non-APN trivariate tables. Counts are tables, not classes.
Equal histograms are never used as proofs of equivalence.
The smaller audit JSON files reproduce the permutation exclusions, switching
quotients, and uniform cubic recipe. Recorded logs are also included.

The comparison coverage is exactly the manuscript's table. Still pending:
other cyclic quadrivariate types, the Gold-switching classes, and BBMM
parameters beyond the fixed primitive element. No full novelty classification
is claimed. Additional exact tables can be checked through the shared driver:

```sh
python compare_all.py --external-manifest extra.json
```

The external manifest is a JSON list whose entries have `family`, `parameters`,
`source`, `field_degree`, `field_modulus`, and `raw_table` (relative to its file).
Each table has 4096 little-endian uint32 entries, in input integer order,
with values 0..4095. Every table passes the same degree, APN, and invariant
checks. A deterministic subset can be selected with `--families` and an
independent output prefix with `--output-prefix`.

Before journal submission, upload this directory and cite a fixed repository
commit or archived release in addition to the repository URL.
