# Corpus body: a plate with a square through-hole

A 60 × 40 × 4 mm aluminium plate with a 12 × 12 mm through-hole at its
centre, dissipating 2 W and cooled on every face (hole walls included) by a
declared h = 10 W/m²/K to 293.15 K. It is the first genus-1 body in the
corpus (bead q61wp.53): the through-opening gives each face an inner
boundary loop, which is what the mesher's facet recovery has to get right.

```bash
python3 examples/plate-hole/generate_plate_hole_stl.py examples/plate-hole/plate-hole.stl
cargo run -p fs-cli --bin frankensim -- --json validate examples/plate-hole/plate-hole.fsim
```

The generator cuts every face on one global tensor grid into diagonal-split
rectangles (no T-junctions, no fans) and asserts that the shell is closed and
consistently oriented, with the analytic volume (9.024 cm³), exterior area
(5504 mm², hole walls included) and Euler characteristic 0.

## No air network

The project declares no fan, fan system, vent or airflow leakage. The
flow-network stage records that absence as its receipt (`"status":
"not-declared"`, authority `declared-absence-of-an-air-network`) and passes
conduction no operating point. Every coefficient is the declared one. A
project that declares only part of an air network still refuses by name.

## The exact check

At steady state, all 2 W leave by convection, so the area-weighted mean
surface temperature is exactly `293.15 + 2 / (10 · A)` = 329.4872 K. It has
to lie between the solved minimum and maximum. If the hole were lost
(A = 5600 mm²), the mean would be 328.864 K. MEASURED 2026-09-29: T_min
329.4643 K, T_max 329.4992 K on 180 tets, energy closure 5e-11. The
freshness lane (`scripts/ci/examples_freshness_e2e.sh`, section 11) checks
the bracket on every run.
