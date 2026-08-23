# Worked example: the minimal project, section by section

This is the smallest project the validator accepts, and its purpose is to
teach the schema: all seventeen sections the validator requires, each
present with the least content that means something. Run it:

```bash
cargo run -p fs-cli --bin frankensim -- --json validate examples/heated-plate/heated-plate.fsim
```

Expected: exit 0, `"status":"ok"`, `"finding_count":0`, plus the canonical
`project_hash` of exactly these bytes. `scripts/ci/examples_freshness_e2e.sh`
runs this continuously and compares the hash against a frozen expectation,
so this walkthrough cannot silently rot.

## What each section declares

| Section | Here | Why it exists |
|---------|------|---------------|
| `metadata` | name, context of use, intended decision, gate, consequence | What question this project may be used to answer |
| `versions` | schema 2, constellation/workspace placeholders | The Five Explicits: which contract versions this assumes |
| `seeds` | root `0x7` | Deterministic randomness is keyed, never ambient |
| `budgets` | 60 s, 64 MiB, 1% relative accuracy | Compute and accuracy are declared, not discovered |
| `capabilities` | `thermal.conduction-solve` | Which registered capability this project targets |
| `units` | SI storage, engineering display | One storage doctrine; display never redefines it |
| `geometry` | one tracked STL artifact by content hash | Geometry is admitted by bytes, not by filename trust |
| `assignments` | half-space selector → region `core` | Mesh-index-free region assignment |
| `assembly` | part `plate`, region `core` | The part/region tree the assignment names |
| `materials` | one binding to the tracked AA6061 card pack | Materials bind by content hash over a validity range |
| `interface-cards` | empty | Declared fact: no TIM/contact cards (absence is not a fact; an empty list is) |
| `power` | 2 W at duty 1.0 into `core` | The heat source, with a duty cycle in [0, 1] |
| `cooling` | no fans, no vents, zero leakage | Declared facts: `(fans)` and `(vents)` with no rows mean "none" |
| `envelope` | ambient 293.15–313.15 K at 1 atm | The operating envelope the requirement is read against |
| `requirements` | one temperature limit with sourced margin | A limit is only meaningful with its authority trail |
| `solver` | auto fidelity, 1e-6 relative tolerance | Solver intent, declared up front |
| `outputs` | one scalar `temperature-max` QoI | What may be asked for, nothing more |

Two deliberate teaching points hide in this minimality:

1. **Every section is mandatory.** Omit any of them and the validator
   refuses with `project-<section>-missing` and the fix. FrankenSim never
   defaults physics.
2. **Empty means empty.** `(interface-cards)` with no rows, `(fans)` with
   no fans, and `(vents)` with no vents are declarations, not omissions —
   the solve prefix consumes these facts instead of inventing defaults.

## What this example deliberately does not do

It does not solve. The solve prefix runs import-verify, assign,
material-resolve, and the flow-network operating point; with no fan system
declared, the airflow stage refuses rather than inventing a fan curve, and
the conduction stage refuses as a typed gap owned by an open bead. See
`examples/cooling-enclosure/README.md` for the full prefix walkthrough and
`docs/QUICKSTART.md` for the end-to-end path.
