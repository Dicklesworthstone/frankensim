# CONTRACT: fs-couple

> **Layout (restored 2026-11-09).** This contract is an index plus part files. A file larger than
> 1 MiB (1,048,576 bytes) comes back with empty content from the GitHub contents API; a writer using
> that API then saw an empty file and replaced the whole contract with its one new section.
> This file was cut to a single section 18 times between 2026-10-30 and 2026-11-09; the text written before 2026-11-05 19:45 (2.3 MB) was absent from then until this restore.
>
> - All earlier contract text is in `contracts/part-NN.md`, unchanged. A higher part number is newer;
>   inside a part, sections are newest first, as they were in the single file.
> - Add a new section under "New sections" below, newest first.
> - When this file passes 600 KB, move the sections under "New sections" into a new
>   highest-numbered part and add its row to the table. No file here may exceed 900 KB.
> - The ten required sections are in the parts named under their headings below.

## New sections (newest first)

### Adaptive laminar enthalpy and external hydraulic continuation

`LaminarEnthalpy::advance_windkessel_adaptive` reuses the original coupled
BDF1 coarse/two-fine controller. Coarse and first-half solves borrow the same
accepted flow/compliance history; the second half borrows the actual first-half
pair. Rejected comparisons advance no external state. The existing ordinary,
Boussinesq and driven entries retain their original semantics and leave the
new optional hydraulic observations empty. No second time controller or
extrapolated physical field is introduced.

An explicit hydraulic norm joins the original caloric, velocity and work
criteria. It checks each port's compliance and added pressure, actual integrated
native Q, separate C/L storage, and signed plus resolved-magnitude resistor,
reservoir and port work. Signed and magnitude Q*dt*depth throughput have their
own cubic-metre tolerance, so equal endpoint flow cannot conceal different
transported volume. Pressure scaling removes the held reservoir offset and
uses the start compliance/Ohmic state; it invents no initial acceleration.
Storage changes and errors use difference of squares before dimensional
multiplication. All spans and local corrections come from the original
first-order controller, including exact target clipping and one-ULP cap tails.

Only accepted fine receipts contribute compensated per-port work and signed/
absolute transported volume totals. Initial Q is independently obtained from
the original accepted BDM edge transport; initial/final C*x²/2 and L*Q²/2
changes use the paired physical states and declared depth. Each port and their
aggregate must satisfy port work minus reservoir work, both resistor losses,
both storage changes and the separate native BDF defect within the explicit
window allowance. This external-circuit identity does not add heat, replace
the independent native caloric balance, or certify whole-system passivity.

All coarse/rejected work consumes the shared native-step and thermal-endpoint
ceilings. Source, material, momentum, stream and thermal refusals remain
terminal. Cancellation or a failed final circuit account rejects the complete
window, including any internally accepted prefix. The final step's retained
compliance coordinates continue with the returned original velocity and clock.
Circuit parameters/references and physical boundary data remain held; this
entry adds no scheduled circuit events or Boussinesq/load combination.

Native regression definitions cover symmetric compliance discharge with zero
fluid/thermal change, hydraulic-only refinement and an independent accepted BE
recurrence; signed equilibrium all-port plug transport and exact work/storage;
nonzero resistor losses excluded from heat; decimal target completion; and
late source, work-budget and cancellation rollback after an accepted prefix.
All seven changed Rust files passed the pinned parser, formatter and
`rustfmt --check`. Bounded independent source review found no paired-history
or accounting blocker. Full native FEM/thermal compilation and these native
integration test executions remain unverified in the current partial restored
environment; no runtime pass is inferred from those checks.

## Parts

| Part | Bytes | `##` sections | First (newest) section | Last (oldest) section |
|---|---|---|---|---|
| [`contracts/part-08.md`](contracts/part-08.md) | 449466 | 145 | Driven full-curve Hertz friction with original thermal sources | CST/DKT hereditary material memory on the original modal basis |
| [`contracts/part-07.md`](contracts/part-07.md) | 442691 | 141 | Thermal porous follower pressure through enthalpy feedback | Electrical contact heat in thermally live actuator supports |
| [`contracts/part-06.md`](contracts/part-06.md) | 446940 | 130 | Complete material and acoustic memory in exponential impact | Shared material cavitation, films and bubble acoustic history |
| [`contracts/part-05.md`](contracts/part-05.md) | 420894 | 27 | Finite solvent heat during gas absorption and exsolution | Conformance tests |
| [`contracts/part-04.md`](contracts/part-04.md) | 449410 | 65 | No-claim boundaries | Live gas changes through the shared acoustic port |
| [`contracts/part-03.md`](contracts/part-03.md) | 444523 | 95 | Finite thermochemical feedback through the complete acoustic load | Nonlinear spatial current with complete material heat exchange |
| [`contracts/part-02.md`](contracts/part-02.md) | 425077 | 113 | Source-resolved Hertz islands on shared finite thermal cells — v0 | Complete sliding-film layer heat on the impact clock |
| [`contracts/part-01.md`](contracts/part-01.md) | 285661 | 57 | Generalized-Newtonian film resistance in impact mechanics | Native joint friction with finite-cell thermal feedback |

## Purpose and layer

Full text: [`contracts/part-05.md`](contracts/part-05.md), under the same heading.

## Public types and semantics

Full text: [`contracts/part-05.md`](contracts/part-05.md), under the same heading.

## Invariants

Full text: [`contracts/part-05.md`](contracts/part-05.md), under the same heading.

## Error model

Full text: [`contracts/part-05.md`](contracts/part-05.md), under the same heading.

## Determinism class

Full text: [`contracts/part-05.md`](contracts/part-05.md), under the same heading.

## Cancellation behavior

Full text: [`contracts/part-05.md`](contracts/part-05.md), under the same heading.

## Unsafe boundary

Full text: [`contracts/part-05.md`](contracts/part-05.md), under the same heading.

## Feature flags

Full text: [`contracts/part-05.md`](contracts/part-05.md), under the same heading.

## Conformance tests

Full text: [`contracts/part-05.md`](contracts/part-05.md), under the same heading.

## No-claim boundaries

Full text: [`contracts/part-04.md`](contracts/part-04.md), under the same heading.
