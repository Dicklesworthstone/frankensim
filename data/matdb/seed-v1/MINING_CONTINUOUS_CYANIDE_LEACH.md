# Continuous gold cyanidation with shared aqueous chemistry

`fs_material::cyanide_leach::flow::CyanidationCstr` models one full,
isothermal, ideally mixed leach tank with equal incoming and outgoing volume
flow. It predicts current tank inventories and the amounts exported during
each interval while the existing chemistry determines pH, free cyanide and
optional competing-metal complexes.

This is the physical connection needed to distinguish gold dissolved in a
tank from dissolved gold actually leaving it. Feed composition and throughput
can change between calls; the tank's inventory, clock and cumulative accounts
continue across those changes. A stopped flow can retain active chemistry.

## Physical and chemical bases

For each analytical amount `n`, the ideal balance is
`dn/dt = Q*C_in - Q*n/V + reaction contribution`. The volume-flow basis follows
the standard stirred-tank mole balance in [Fogler, chapter 1,
equations 1-4/1-8](https://www.pearson.com/content/dam/one-dot-com/one-dot-com/hei-samples/science/9780135337554.pdf).
All entrained solids and solution components share that residence in this
model. `Q`, `V` and feed concentrations use the same declared **solution-volume
basis**; solid displacement, slurry density and ore tonnage are separate inputs
not inferred here.

The independently declared local kinetic law is
`r_Au = k*n_Au*c_free_CN*c_O2`, in mol Au/s. The original conserved dissolved-O2
Elsner reaction supplies its contributions:

| Analytical component | Reaction contribution |
|---|---:|
| Accessible solid Au | `-r_Au` |
| Cyanide ligand outside the gold complex | `-2*r_Au` |
| Dissolved Au(CN)2− | `+r_Au` |
| Dissolved O2 | `-r_Au/4` |
| Other families, spectator ions and metal totals | `0` |

The cyanide-ligand inventory includes HCN, free CN and cyanide bound to the
declared competing metals. It excludes the ligands already counted in the
separate gold-complex inventory. The original acid/base and metal-equilibrium
owners partition the current analytical totals jointly. Transporting a
constant initial free-CN fraction would give a different, unsupported model.

## Inputs and source authority

Construct `LocalGoldCyanidation::from_card(chemistry, cyanide_family, source)`
with the distinct version-one law `local-speciated-gold-cyanidation` and one
`rate-constant` in **m6/(mol2 s)**. The supplied card and its evidence remain
unchanged. Required typed validity coordinates cover actual pH, free-CN and
DO concentrations, accessible-Au inventory, temperature and solution volume.
The original owned total-CN, gold-complex and activity-basis coordinates can
also be constrained. Extra source coordinates remain explicit caller inputs.

Batch-age, elapsed-time and residence-time cards are refused as local rate
authority. Neither equilibrium data nor the chemical equation provides `k`.
The physical source must support the local model being invoked.

Initialize `CyanidationCstr::new` with a `CyanidationCharge`, source conditions
and bounded speciation controls, or use `new_with_competing_metals` for the
original source-backed metal families and their finite initial totals.
The caller's original initial charge is retained. Internal inventories have
a canonical component order.

`CyanidationFeed` supplies `AmountConcentration` values in **mol/m3** for
accessible gold, cyanide ligand, dissolved oxygen and dissolved gold complex.
It also supplies every initialized background acid-family ID, spectator
species/charge and competing-metal ID exactly once, including zero feeds.
Missing, duplicate, extra or charge-mismatched feed components refuse.

Advance with `advance(duration, volumetric_flow, &feed, numerics, cancelled)`.
The common activity-reference amount is independent of volume. For the NBS
unit-molality convention it corresponds to solvent mass times 1 mol/kg; a
physical liquid volume is not itself an activity basis. Temperature, solution
volume and this solvent basis remain fixed throughout a run.

## Results and numerical meaning

Each accepted result retains the current analytical inventory, actual shared
speciation, original kinetic conditions, optional metal-equilibrium report,
and original source owners. It reports interval and cumulative inlet/outlet
amounts, gold dissolution and original-equation amounts separately.

For every component the independently recomputed account is
`final + outlet - initial - inlet - reaction contribution`. Au and total
cyanide accounts combine the appropriate components; the actual ligand
solver's signed residual remains visible. Sodium and dissolved oxygen have
their own balances. Dissolved-gold export is solution inventory, not carbon
recovery, doré production or refined metal.

The numerical method uses exact positive transport on each side of a bounded
RK4 reaction step. Full-step/two-half-step differences divided by three estimate
local error for the second-order split in an L1 amount norm over stocks, independently
accumulated exports and reaction progress. `FlowNumerics` provides an absolute
amount allowance, relative allowance, maximum step and bounded work/attempts.
The relative scale includes all admitted component amounts; select an
appropriate absolute allowance when trace-metal accuracy matters.

Source support is checked at all evaluated split endpoints and reaction
stages. This sampling is not an enclosure of the continuous unsplit path.
Numerical estimates do not include equilibrium tolerance, source uncertainty
or model error. Negative or unresolved stocks, invalid source points, exhausted
work and cancellation leave the entire previously accepted tank unchanged.

## Execute the actual-source example

From the repository root:

```bash
continuous_leach_pack_dir=$(mktemp -d)
for pack in hydrogen-cyanide-aqueous-nbs1989 hydrogen-ion-aqueous-nbs1989 cyanide-ion-aqueous-nbs1989 hydroxide-ion-aqueous-nbs1989 water-liquid-nbs1989-thermochemistry; do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$continuous_leach_pack_dir/$pack.fsmatpk"
done

cargo run -p fs-material --example continuous_gold_leach -- \
  --materials "$continuous_leach_pack_dir"

cargo test -p fs-material cyanide_leach::flow::tests
```

The example uses actual NBS water/HCN/CN/OH equilibrium sources and prints their
original receipts. The local `k=0.005 m6/(mol2 s)`, inventories, flows and validity
boxes are explicitly authored numerical inputs, not measured ore calibration.
It increases throughput, stops flow while dissolution continues, and restarts
with a changed dissolved NaCN feed. It prints retained gold, dissolved and
solid gold export, pH/free CN and Au/CN/O2/Na balance residuals.

The focused regressions cover analytic nonreactive wash-in/wash-out, zero-flow
agreement with equivalent batch chemistry, reaction-plus-transport accounts,
finite copper competition, feed/flow changes, source ownership and atomic
refusal. These are implementation checks, not plant validation.

The base `CyanidationCstr` receives oxygen through the initial liquid and feed.
The aeration extension below adds a separately declared gas boundary. Carbon
adsorption/CIL, solids slip, liberation, competing-metal dissolution/redox,
mineral precipitation, variable volume, thermal feedback and finite solvent
H/O balance remain separate models. The single tank is not a fitted tank train
or an RTD deconvolution result. Finite hydrated-lime preparation is available
in [MINING_LIME_EQUILIBRIUM.md](MINING_LIME_EQUILIBRIUM.md); carrying that solid
phase inside a flowing leach requires its own phase-transport balance.

## Add prescribed-pressure oxygen exchange

`flow::aeration::AeratedCyanidationCstr::from_tank(tank, transfer)` attaches
the existing `cyanide_leach::aeration::OxygenTransfer` to this continuous tank.
It preserves the original local kinetic card, charge, elapsed time and all
previous liquid/reaction accounts, including an already advanced tank.
The inner tank remains read-only through `tank()`.

Advance with
`advance(duration, flow, &feed, oxygen_partial_pressure, &transfer_conditions, numerics, cancelled)`.
The maintained gas reservoir contributes to the dissolved-oxygen balance:

```text
dn_O2/dt = Q*C_O2,in - Q*n_O2/V + kLa*(n_sat - n_O2) - r_Au/4
n_sat   = n_ref*K*p_O2/p_ref
```

The original NBS gas/liquid formation-Gibbs sources supply K. The original,
independent transfer card supplies kLa and its applicability; local leach
kinetics retain their own source. See
[MINING_OXYGEN_TRANSFER.md](MINING_OXYGEN_TRANSFER.md) for the source/activity
conventions and the transfer-law reference. The activity-reference amount
remains distinct from physical liquid volume.

Feed, flow and pressure can change between calls. A stopped liquid flow leaves
gas exchange and reaction active. Setting pO2 to zero declares a zero-oxygen
reservoir that can strip liquid oxygen; it does not disable gas contact.
The reservoir has prescribed pressure, with no inferred finite gas supply or
aerator-power calculation.

`AeratedFlowAdvance.flow` retains the liquid/reaction report.
`oxygen_in` and `oxygen_out` separately report gross gas exchange during the
interval, with corresponding cumulative fields. These amounts are integrated
independently of liquid inlet/outlet and Au consumption. All participate in
the reported interval and cumulative oxygen/component balances; the original
charge and liquid streams are never altered to absorb a gas imbalance.
The result retains the original transfer source and final admitted
`transfer_conditions`.

Symmetric exact liquid transport and exact gas relaxation surround the
existing reaction integration. Raw full/two-half differences control stocks,
liquid exports, Au progress and both gross gas amounts without a smooth-order
divisor at flux reversals. The L1 estimate and work limits use
`FlowNumerics`; its relative scale includes the component and gas-exchange
amounts. Choose explicit tolerances for the trace inventory of interest.
Source checks and cancellation cover evaluated stages, and refused advances
publish no partial history. The numerical estimate remains separate from
source uncertainty, model error and a continuous-path validity enclosure.

### Run the coupled actual-source example

```bash
continuous_oxygen_pack_dir=$(mktemp -d)
for pack in hydrogen-cyanide-aqueous-nbs1989 hydrogen-ion-aqueous-nbs1989 cyanide-ion-aqueous-nbs1989 hydroxide-ion-aqueous-nbs1989 water-liquid-nbs1989-thermochemistry oxygen-gas-nbs1989-thermochemistry oxygen-aqueous-nbs1989-thermochemistry; do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$continuous_oxygen_pack_dir/$pack.fsmatpk"
done

cargo run -p fs-material --example continuous_aerated_gold_leach -- \
  --materials "$continuous_oxygen_pack_dir"

cargo test -p fs-material cyanide_leach::flow::aeration::tests
```

The example starts with zero initial and feed dissolved oxygen; explicit gas
uptake supports gold dissolution. It then changes throughput and gas pressure,
stops liquid flow while stripping oxygen, and restarts with changed cyanide
feed. Its local k=0.005 m6/(mol2 s), kLa=0.02/s, source validity boxes and
operating schedule are authored numerical inputs, not measured calibration.
It prints the seven original pack receipts, liquid/gas oxygen accounts,
dissolved-gold inventory and export, pH/free CN and material-balance residuals.

Tests compare nonreacting flow/gas behavior with an independently solved
unsplit ODE, including uptake-to-stripping reversal within one interval.
They also check zero-transfer and zero-flow limits, gas-enabled reaction,
finite copper competition, preattachment history, changed schedules and
atomic source/input/work/cancellation refusal. The extension retains the
fixed-temperature, single-tank scope and does not add carbon adsorption,
finite-headspace dynamics, phase precipitation or a tank-train model.
