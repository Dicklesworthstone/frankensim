# Supply hammer materials to the physical piano

`--hammers materials.fsh` replaces the per-key felt and relaxation cards in the
actual hammer/string contact solve. It composes with `--scale`, either board
input, the Steinway-D reconstruction, and performance scores. Without this
option the existing material selection is unchanged.

```sh
cargo run --release -p fs-couple --example grand_piano -- \
  --scale strings.csv --hammers materials.fsh --render piano.wav
```

A file covering a scale containing only key 69 could be:

```text
frankensim-hammer-materials-v1
# Illustrative parameters, NOT a measured specimen or a fitted coupon.
felt,69,400000,0.2,2.5,3.2,0.25,0.8,2500000
branch,69,2000000,0.0002
branch,69,500000,0.004
```

The `felt` columns after the tag are key, reference stress [Pa], reference strain,
loading exponent, unloading exponent, crush fraction, densification strain, and
Prony equilibrium modulus [Pa]. The `branch` columns are key, branch modulus
[Pa], and relaxation time [s]. Zero to eight branches may follow each felt card;
no branches means no creep memory, not selection of default relaxation.
`#` introduces a comment. All scalars must be finite and the existing material
owners enforce their physical domains. The creep-solid equilibrium modulus must
be positive. Unknown records, duplicate cards, branches before their felt card,
and missing or extra keys are errors.

Supply exactly one card for **every key in the scale**, including unplayed keys.
For the default or preset 88-key scale that means 88 cards. Card order need not
match scale order. Every unison member retains its own contact/crush/creep state;
sharing a material does not share its history.

These are the existing `WoolFelt` and `GeneralizedMaxwell` inputs. The existing
creep realization replaces the Prony solid's instantaneous spring with the
unilateral felt envelope and couples the resulting deformation to the hammer
solve. This is not a parallel stress addition, an output filter, a sample bank,
or a newly invented material law. Changing the cards can change attack, rebound,
loss and radiated response while preserving string tensions and geometry.

Hammer mass, felt patch area/thickness, and unison geometry still come from the
string scale. The preset still supplies its geometry-derived shank and jack port;
material import does not turn it into the point-mass comparison image.

The file is a route for declared or externally identified material parameters,
not a fitting tool or a measurement certificate. Dynamic coupon agreement,
frequency-band validity, instrument calibration and real-time deadlines still
need their own evidence. Tests exercise identity/order, refusal, actual contact
response, energy closure and preset-shank composition; native results must be
obtained by running `cargo test --release -p fs-couple --example grand_piano`.
