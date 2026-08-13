# Worked example: the refusal/fix loop

Refusals are a feature. FrankenSim's validator is deliberately fail-closed:
a project that declares something the physics stack cannot honor is refused
at the front door, with a named code and a concrete fix, instead of
contaminating a solve. This example normalizes that loop.

`broken.fsim` is the tracked reference project with exactly one token
changed: the power dissipation's duty cycle is declared as `:duty 2.0` —
a 200% duty cycle, which is not a thing. Everything else, including the
canonical spelling the parser enforces, is untouched (the delta is proven
one-token-wide by the e2e lane).

## Step 1 — watch it refuse

```bash
cargo run -p fs-cli --bin frankensim -- --json validate examples/refusal-loop/broken.fsim
```

Expected: a nonzero exit (`4`, the schema/semantic-refusal class) and a
finding on stderr carrying:

- the code `project-duty-range`,
- what happened (`duty for \`air\` is 2`), and
- the fix (`duty must lie in 0.0..=1.0`).

That triple — code, what, fix — is the shape of every FrankenSim refusal,
from validation findings through solve-stage refusals to package checks.
Agents branch on the code; humans read the fix.

## Step 2 — fix it and re-validate

Change `:duty 2.0` back to `:duty 1.0` (or any value in `[0, 1]`) and run
the same command against your edited copy: `"status":"ok"`,
`"finding_count":0`. The fixed file is byte-for-byte the tracked
`data/reference-project/cooling-reference.fsim`, which the full worked
example (`examples/cooling-enclosure/`) then takes through import and
solve.

## Why the example is this small

One token is the point. A refusal loop teaches nothing if the broken
fixture is exotic; the lesson is that ONE bad declared value in an
otherwise perfect project is caught, named, and mapped to its repair
before any stage runs. The e2e lane validates both directions
continuously: `broken.fsim` must keep refusing with exactly this code, and
its one-token repair must keep validating clean.
