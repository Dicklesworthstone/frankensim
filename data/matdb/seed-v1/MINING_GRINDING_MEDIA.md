# Grinding media: discrete product data and charge inventory

Two manufacturer profiles supply eight claims through the existing material
compiler and catalog. They support explicit media inventories and geometry
inputs to the existing milling models. They do not infer ore grindability,
impact breakage, liner life or media consumption from hardness.

## Forged steel balls

`gerdau-forged-grinding-balls-high-hardness` retains seven discrete nominal
size rows from page 2 of Gerdau's [Grinding Balls brochure](https://gsn.gerdau.com/sites/gsn_gerdau/files/downloadable_files/GrindingBalls_Brochure_v1.pdf).
For each row it provides nominal ball mass, average-diameter specification
bounds, out-of-round tolerance, separate surface and average volumetric
Rockwell C nominal hardness, and the stated hardness tolerance.

Nominal millimetre labels select product rows. They are not interchangeable
with converted inch labels: the 50 mm row has nominal mass 0.536 kg and
average-diameter bounds 51.08–52.02 mm. The 76 mm row distinguishes nominal
surface hardness 62 HRC from nominal average volumetric hardness 61 HRC.
The ±4 HRC specification tolerance is a separate offset, never a statistical
uncertainty. These dimensionless specification records cannot satisfy the
measured-hardness API, which requires original apparatus and loading history.
No interpolation, intrinsic steel density or wear coefficient
is supplied. Temperature and lot qualification remain unknown.

```bash
cargo run -p xtask -- matdb-catalog \
  --pack gerdau-forged-grinding-balls-high-hardness \
  --property grinding-ball-nominal-mass \
  --at source-nominal-ball-size-label=0.05 \
  --at source-gerdau-high-hardness-forged-steel=1 \
  --at source-product-specification=1 \
  --at source-test-temperature-known=0 --at source-lot-qualified=0
```

An explicitly authored inventory of 2,000 balls in that row therefore has
nominal mass 1,072 kg. An actual charge still requires its counts, sizes,
wear state and inventory measurements. For an ideal spherical critical-speed
calculation, `fs_material::comminution::ball_mill_critical_speed` takes the
free diameter inside the mill liner, an explicitly selected ball diameter,
and gravity. Product average-diameter limits do not certify the surface of
an out-of-round ball, and the ideal threshold is not a mill operating setting.

## Apparent charge density

`molycop-sag-charge-packing-reference` retains the **theoretical apparent
packing density** 4,520 kg/m³ from [Molycop's SAG ball specification](https://www.molycop.com/products-and-services/consumables/grinding-media/molycop/sag-balls).
Its property name explicitly identifies that basis. A declared 30 m³ bulk
ball-charge volume, including voids, corresponds to 135,600 kg under this
reference assumption. This is not intrinsic steel density, an ore/slurry
density or a site-measured fill law. Querying generic `density` refuses.

Required selectors are `source-molycop-sag-ball-charge=1`,
`source-theoretical-packing-reference=1`, `source-site-packing-qualified=0`
and `source-test-temperature-known=0`.

## Verification

`cargo test -p xtask --test mining_grinding_media_sources` compiles the original
manifests, round-trips and verifies their packs and receipts, independently
checks every source row, exercises the catalog, and rejects unsupported sizes,
unit/scale substitutions and invented ore properties. G0/G3/G5 cover this
source-to-query behavior; no plant-calibration or wear-performance claim follows.
