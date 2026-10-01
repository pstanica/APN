The tools share one field implementation (BinaryField in compare_all.py) and
one invariant engine (invn_canonical.cpp). Requires Python 3, numpy, and g++.
From tools/ run:
    python compare_all.py
    python audit_invariants.py
    python audit_switching.py

compare_all.py produces comparison_manifest.json and comparison_summary.json.
Each record includes exact parameters, polynomial field modulus, coordinate
encoding, SHA256 of the table, ANF degree, APN rank check, and full histogram.
Equal invariants are never treated as proofs of equivalence.

Optional deterministic batch selection:
    python compare_all.py --families Gold Taniguchi_normalized --output-prefix batch

Additional raw tables use --external-manifest path.json. That file is a JSON
list; each object has family, parameters, source, field_degree, field_modulus,
and raw_table (a path relative to the manifest). Tables must contain exactly
4096 little-endian uint32 entries in input integer order, values 0..4095.
They pass the same degree, APN, and invariant checks as built-in families.

Coverage is explicit in the manuscript and manifest. Remaining reported
cyclic types and Gold-switching classes need exact builders or raw tables.
BBMM currently fixes one primitive u; no reduction over all primitive u is
asserted. Thus the data do not establish novelty against every known family.

verify_identities.py (python-flint) and verify_trace_apn.py (numpy) retain the
symbolic and cubic-formula checks from the previous package. They are not
separate engines for the new novelty comparisons.
