# Evaluated radioactive-decay fixtures

These two unmodified ENDF/B-VIII.0 files exercise the radioactive-decay reader
and the existing physical inventory with published data, including a stable
daughter and the radiation records following the decay-mode list.

Source: the [Sandia National Laboratories cross-section data archive](https://github.com/sandialabs/xsec/tree/aaafe23d5f80a5eefbf49421a8e7ff591676c13a/ENDFB-VIII-0/ENDF-B-VIII.0_decay),
pinned to commit `aaafe23d5f80a5eefbf49421a8e7ff591676c13a`.
The original files identify their NNDC/ENSDF or Nuclear Wallet Cards evaluations
and authors. The official library is described by the
[NNDC ENDF/B-VIII.0 release](https://www.nndc.bnl.gov/endf-releases/?version=B-VIII.0).

| File | MAT | Original Git blob |
| --- | ---: | --- |
| `dec-001_H_003.endf` | 4 | `306f60015644d4ba758217d6b4e66ac96986998f` |
| `dec-002_He_003.endf` | 9 | `6941a4d431ca4a9b48a77d8fad807a695c84f957` |

The tritium record supplies a half-life of `3.887896e8` seconds, uncertainty
`6.311520e5` seconds, and average electron-related energy `5690` eV per decay.
Its branch Q value is `18590` eV. These are separate quantities: assigning Q to
deposited heat would incorrectly include energy outside the supplied aggregate
decay-heat categories. The three-category record does not supply a numerical
neutrino-energy component, so the example does not infer one from Q.

The helium-3 record explicitly marks the daughter stable. The model bridge
requires that record; it does not infer stability from missing daughter data.

These fixtures test import and model interoperability with one identified
evaluation. They are not a complete isotope catalogue or independent empirical
validation. Source uncertainties are retained as metadata; the deterministic
example propagates the reported central values.

Run the data-to-inventory example with:

```sh
cargo run -p fs-nuclear --example evaluated_decay
```
