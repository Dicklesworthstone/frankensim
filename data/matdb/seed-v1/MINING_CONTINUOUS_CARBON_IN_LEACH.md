# Continuous aerated CIL with explicitly retained carbon

`cyanide_leach::flow::cil::AeratedCilCstr` combines the existing continuous
cyanidation tank, original oxygen-transfer source, and independently declared
carbon populations. Slurry feed and liquid exports, gold dissolution, oxygen
uptake/stripping, and reversible carbon loading share one evolving solution.
Free cyanide, pH, and optional Cu complexes are re-solved from actual inventories.

Every carbon population in this owner is explicitly **retained in the tank**.
A `PregRobbing` role classifies competing capture; it does not imply that
carbonaceous ore enters or washes out with the slurry. Product-carbon loading
remains in-circuit inventory until a separate recovery operation removes it.

## Run the actual-source example

The [example](../../../crates/fs-material/examples/continuous_aerated_carbon_in_leach.rs)
uses eleven original NBS profiles. From the repository root:

```bash
CIL_PACK_DIR=target/mining-continuous-cil-packs
mkdir -p "$CIL_PACK_DIR"
for pack in \
  hydrogen-cyanide-aqueous-nbs1989 \
  hydrogen-ion-aqueous-nbs1989 \
  cyanide-ion-aqueous-nbs1989 \
  hydroxide-ion-aqueous-nbs1989 \
  water-liquid-nbs1989-thermochemistry \
  copper-i-ion-aqueous-nbs1989 \
  dicyanocuprate-i-ion-aqueous-nbs1989 \
  tricyanocuprate-i-ion-aqueous-nbs1989 \
  tetracyanocuprate-i-ion-aqueous-nbs1989 \
  oxygen-aqueous-nbs1989-thermochemistry \
  oxygen-gas-nbs1989-thermochemistry
do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$CIL_PACK_DIR/$pack.fsmatpk"
done
cargo run -p fs-material --example continuous_aerated_carbon_in_leach -- \
  --materials "$CIL_PACK_DIR"
```

Canonical packs already compiled for the earlier examples can be reused.
The loader verifies their IDs, content hashes, and canonical serialization.

The example begins with zero dissolved oxygen in both tank and feed. It runs
two otherwise identical scenarios: retained competing carbon alone, and that
same population alongside retained product carbon.

| Interval | Slurry flow | Feed or gas change |
| --- | --- | --- |
| 0–30 s | 0.0001 m³/s | Oxygen supply at 21 kPa O₂ partial pressure starts leaching. |
| 30–70 s | 0.0001 m³/s | Delivered cyanide and sodium concentrations increase together. |
| 70–90 s | Zero | Slurry exports stop; a maintained zero-O₂ boundary strips dissolved oxygen. |
| 90–120 s | 0.0002 m³/s | Original feed and 21 kPa O₂ boundary resume at higher throughput. |

Each phase uses one-second continuation windows with its original accuracy
allowances apportioned across the calls. Component accounts are checked after
every window; output distinguishes the final real one-second interval from
cumulative totals.

Zero oxygen partial pressure represents an oxygen-free maintained reservoir;
it is not a command to turn gas contact off. Zero slurry flow stops liquid
transport while reaction, carbon exchange, and gas transfer remain active.

## Explicit inputs and component accounts

The incoming `CyanidationFeed` declares accessible Au, cyanide ligand outside
gold complexes, aqueous Au, dissolved O₂, every background acid family, every
spectator, and every competing-metal total on the same liquid-volume basis.
Changing delivered cyanide or strong-ion concentrations changes calculated pH;
the caller does not prescribe a pH result. Entrained solids share the declared
liquid residence; solids displacement and slip are not inferred.

The local leach card explicitly supplies `r_Au = k*n_Au*c_free_CN*c_DO`.
Each retained population supplies its own Langmuir capacity/affinity, dry mass,
preloaded amount, and local loading-relaxation rate. The gas model supplies its
original O₂ solubility and a separate kLa source. No coefficient is inferred
from an ore or carbon name.

The same runtime accepts `OxygenTransfer::from_henry_source` for an explicitly
selected measured Hcp source. That route uses `n_sat = V * Hcp * pO2`, requires
the original pure-water/infinite-dilution qualification flags, and retains its
pinned query receipt in `henry_answer`. It does not synthesize Gibbs data from
solubility; the NBS example above deliberately exercises the thermodynamic route.

`IonPairTransfer::sodium()` transfers one Au(CN)2− and one Na+ into a neutral
sorbed pair. Both cyanide ligands remain bound; desorption returns the complete
pair. Counterion inventory and feed concentration must cover the corresponding
accessible-plus-aqueous gold pools, preserving the declared model without
inventing a counterion-limited adsorption rate.

Reports distinguish retained carbon, current mobile stocks, integrated liquid
inlets/outlets, cumulative gold dissolution, and gross oxygen input/output.
Whole-history accounts include Au, cyanide, the selected counterion, Cu,
dissolved oxygen plus stoichiometric consumption, and occupied/vacant sites.
The example independently recomputes these mobile/carbon/gas component balances.

Constructors require an explicit initial oxygen partial pressure and additional
transfer conditions, so the initial carbon and gas source points are admitted.
Attachment to an existing aerated tank retains its accepted liquid and gas
history. Continued advances preserve inventories and numerical history;
source refusal, cancellation, or exhausted work budgets publish no partial state.

## Interpretation

The NBS records supply equilibrium thermochemistry. The example's local leach,
kLa, and carbon coefficients are authored numerical inputs, not measured ore,
slurry-aerator, or carbon calibration. Its calculated capture comparison is not
a plant-recovery forecast. Numerical local-error estimates do not quantify
source or model uncertainty.

This example performs no carbon harvest, carbon replacement, entrained-carbon
transport, inter-tank carbon transfer, tank cascade, variable-volume operation,
gas-supply-capacity calculation, nonideal activity correction, cyanide
destruction, or heat balance. Related background:
[finite simultaneous CIL](MINING_CARBON_IN_LEACH.md) and
[copper–cyanide sources](MINING_COPPER_CYANIDE.md).
