# Simultaneous gold leaching and carbon capture

The `cyanide_leach::cil::CarbonInLeachBatch` runtime advances accessible-gold
dissolution, dissolved-gold concentration, and separately declared sorbent
populations on one timeline. Every rate stage uses the original HCN/water
equilibrium owner and optional dissolved-metal competition to resolve free
cyanide, pH, and Cu complexes. Dissolved oxygen remains a finite liquid stock.

The [runnable example](../../../crates/fs-material/examples/simultaneous_gold_carbon_in_leach.rs)
compares the same natural-carbon population with and without added product
carbon. It uses nine actual NBS thermochemical profiles and explicit synthetic
leach, isotherm, and local-transfer coefficients. It continues each bath through
10, 30, and 80 second intervals and checks whole-bath Au, cyanide, sodium, Cu,
oxygen, and charge accounts.

## Run the example

From the repository root, compile the original profiles:

```bash
CIL_PACK_DIR=target/mining-cil-packs
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
  tetracyanocuprate-i-ion-aqueous-nbs1989
do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$CIL_PACK_DIR/$pack.fsmatpk"
done

cargo run -p fs-material --example simultaneous_gold_carbon_in_leach -- \
  --materials "$CIL_PACK_DIR"
```

The executable's `--help` also lists these packs. Existing canonical compiled
packs can be reused. The loader verifies pack IDs, content hashes, and canonical
serialization before constructing the source-backed equilibria.

## What is coupled

The leach law remains the caller's explicit
`r_Au = k * n_accessible_Au * c_free_CN * c_DO`. Each retained sorbent has
its own original Langmuir equilibrium card and a separate
`local-aqueous-langmuir-linear-driving-force` card:

```text
q_eq(c) = q_max * c / (c_half + c)
dA_i/dt = k_i * (M_i * q_eq(n_aqueous_Au / V) - A_i)
```

Dry mass, capacity, affinity, relaxation rate, and any preloaded gold are explicit
inputs. Population roles classify the signed transfer accounts:

| Role | Meaning |
| --- | --- |
| `ProductCarbon` | Gold captured on the explicitly supplied recoverable carbon population. |
| `PregRobbing` | Gold captured by the separately supplied natural or other competing sorbent. |

Roles supply no isotherm, kinetic coefficient, or ore calibration. Both kinds
share the actual finite dissolved-gold stock and can desorb. Preloaded sorbents
contribute their original Au, cyanide, and counterion inventories.

Adsorption chemistry is an explicit modeling choice. `IonPairTransfer::sodium()`
declares transfer of one Au(CN)2− ion and one Na+ ion into a neutral sorbed pair.
Desorption returns both ions and both cyanide ligands. A custom pair requires
an exactly conserved reaction with admitted species formulas and charges; labels
alone do not establish chemical identity.

The aqueous Au and counterion amounts change together. The remaining available
HCN/CN/metal-ligand total changes only through leaching, by two cyanide ligands
per Au dissolved. The neutral pair therefore preserves electroneutrality and
free-ligand equilibrium at fixed leach progress. Because the stated dissolution
law has no product inhibition, capture changes the destination of released gold
without increasing that law's dissolution rate.

## Stocks, sources, and continuation

The complete gold account includes remaining accessible metal, dissolved
Au(CN)2−, and every sorbent. The cyanide account also includes uncomplexed HCN/CN,
Cu-bound ligands, and two ligands per aqueous or sorbed Au. The counterion account
includes both dissolved ions and retained pairs. Oxygen use is one quarter mole
O2 per mole Au dissolved. Every sorbent retains occupied and vacant capacity.

Attach to a fresh gold bath with zero elapsed time and zero leach progress.
The chosen monovalent counterion must be an explicit original spectator.
Construction rejects insufficient counterion support for the full admitted
gold inventory; the model does not invent cation-limited adsorption kinetics.
Subsequent advances preserve all accepted stocks and cumulative error estimates.

Every kinetic stage rechecks the original source domains using actual aqueous
gold concentration, pH, free cyanide, counterion concentration, dissolved oxygen,
loading, temperature and liquid volume. The remaining known leach coordinates
also come from the actual common solution; caller inputs cannot freeze them.
Unsupported
sources, exhausted work budgets, or cancellation publish no partial state.
Temporal error estimates concern the numerical integration; they do not quantify
uncertainty in the chosen chemistry, isotherms, or rate coefficients.

## Scope and interpretation

This is a simultaneous finite-bath CIL model. The older
`gold_leach_to_carbon` example runs leaching first and then hands analytical
dissolved Au to a separate carbon contact with externally specified pH.

The new example's product capture and preg-robbing comparison are calculated
outcomes of authored coefficients. The NBS records supply equilibrium
thermochemistry; they supply no measured ore dissolution rate, carbon affinity,
natural-carbon uptake, particle liberation, or plant recovery forecast.

This owner does not supply continuous solution or slurry flow, a CSTR train,
inter-tank carbon movement, gas-phase oxygen replenishment, simultaneous reagent
dosing, cyanide destruction, copper-complex adsorption, nonideal activities,
or thermal feedback. Oxygen exhaustion remains physical: carbon does not
replenish dissolved O2.

Related source details: [aqueous cyanide](MINING_AQUEOUS_CYANIDE.md),
[copper–cyanide competition](MINING_COPPER_CYANIDE.md), and
[finite gold leaching](MINING_CYANIDE_LEACH_KINETICS.md).
