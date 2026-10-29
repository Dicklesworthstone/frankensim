# Violin varnish: historical layers and measured modern films

Seven source packs retain **52 claims**: six modern DMA specimens (48 claims)
and one historical stratigraphy pack (four tabulated claims). They compile
through the ordinary `xtask matdb-pack` path and are discoverable through
`matdb-catalog`. Historical identification and a modern mechanical surrogate
have different source identities.

## Historical instrument evidence

The [Echard et al. study](https://doi.org/10.1002/anie.200905131) identified a
penetrating drying-oil treatment beneath an oil/Pinaceae-resin surface film.
Its sampled instruments have these approximate observed ranges:

| `specimen_code` | Stradivari specimen | External film, µm | Oil penetration into wood, µm |
| --- | --- | ---: | ---: |
| 1 | Long Pattern, circa 1692 | 10–30 | 30–100, maple |
| 2 | Davidoff, 1708 | 10–30 | 10–30, spruce |
| 3 | Provigny, 1716 | 10–30 | 10–30, spruce |
| 4 | Viola d'amore head, circa 1720 | 10–30 | 30–100, maple |
| 5 | Sarasate, 1724 | 10–30 | Not supplied |

The source reports no pigment detected in Long Pattern; probable cochineal
lake in Provigny; lake plus iron oxides in Davidoff and the viola head; and
previously identified vermilion in Sarasate. These are identifications, not
pigment fractions or evidence for every Stradivari instrument. Mineral ground,
protein sealer and shellac are not established by these samples.

`varnish-stradivari-echard2010-stratigraphy` retains range endpoints separately
and refuses interpolation between specimen codes. Its penetration depth must
not be added as an external film or counted twice as wood plus solid varnish.
Density, modulus and damping cannot be queried from this historical pack.

## Modern mechanical source

[Lämmlein et al.](https://doi.org/10.1016/j.rinma.2020.100137), Tables 1–2,
characterized three modern systems, with two independent DMA fits apiece:

| Pack family, ending `-lammlein2021-s1` or `-s2` | Tested composition | Density mean / SD, kg/m³ | Density sample count |
| --- | --- | ---: | ---: |
| `varnish-oil-copal` | Hammerl golden-brown standard oil varnish | 1097 / 26 | 5 |
| `varnish-spirit-shellac` | Orange shellac, elemi and spike-lavender oil in ethanol; 1704 type | 984 / 30 | 7 |
| `varnish-mineral-ground` | OldWood classical amber oil varnish with 29 wt% pumice | 1155 / 30 | 4 |

Films dried at 20 °C/65% RH, then conditioned at 20 °C/35% RH for at least
three months. DMA used 0.03% strain and 0.1–30 Hz sweeps at multiple
temperatures; time-temperature superposition supplied the 20 °C master curve.
Table 2's E0 values came from separate static tests at 65% RH. The packs retain
the published combined fit, not an independently measured static law at 35% RH.

The executable fit band is deliberately **20–20,000 Hz**. This is a conservative
subset of the source's fit coverage, not direct DMA measurements at each audio
frequency. A different temperature or conditioning state refuses. Ageing,
fresh-film curing and changing humidity require additional laws and data.

For cyclic frequency `f` in Hz, `x = 2π f τ`, the published complex law is

`E*(f) = E0 + (Einf − E0) (i x)^alpha / (1 + (i x)^beta)`.

`Einf` is a fit scale, **not a finite infinite-frequency modulus** when
`alpha > beta`. Storage is `Re E*`; loss modulus is `Im E*`; loss tangent is
their ratio. None of the five coefficients is a Kelvin–Voigt viscosity.

## Query and source semantics

Every modern pack exposes `varnish_e0` and `varnish_einf` [Pa],
`varnish_alpha` and `varnish_beta` [1], and `varnish_tau` [s]. Query these with
typed coordinates `temperature=293.15` K, `conditioning_relative_humidity=0.35`,
`source_fitted_response=1` and cyclic `frequency` within the admitted band.
For example:

```sh
cargo run -p xtask -- matdb-catalog \
  --pack varnish-oil-copal-lammlein2021-s1 --property varnish_einf \
  --at temperature=293.15 --at conditioning_relative_humidity=0.35 \
  --at source_fitted_response=1 --at frequency=1000
```

Density additionally requires `source_density_sample_matched=0`: Table 1's
group mean is not paired with an individual DMA specimen. The separate
`density_standard_deviation` and `density_sample_count` preserve the reported
statistics. All uncertainty models remain `Unstated`; the SD is not silently
converted into a confidence interval, parameter covariance or prediction bound.

## Applying these data to a violin

Layer geometry belongs to the specimen. A 20 µm coat is a declared scenario,
not the exact thickness of any listed historical instrument. A mechanical use
of the modern oil-varnish fit for historical oil/pine resin is explicitly a
surrogate. No source supplies Poisson ratio or a complete isotropic tensor;
the consumer must name such a closure and its assumptions.

The shared laminate owner must integrate each layer's mass and A/B/D stiffness
at its actual offset. A one-sided coat has membrane–bending coupling and a
shifted mass centroid. Loss follows layer strain energy. Both stiffness and
added inertia influence each eigenfrequency; no universal upward/downward
pitch shift or acoustic superiority follows merely from the word “varnish”.

[Moisture experiments](https://doi.org/10.1038/s41598-019-54991-5) separately
show that complete modern finishes and worn pretreatments protect wood
differently. Their observed 40–70 µm external films do not replace the
historical thickness ranges or provide a vapor-permeability coefficient.
Open f-holes and uncoated internal wood also prevent treating a violin as
hermetically sealed.

Missing capabilities remain explicit: historical bulk viscoelasticity,
thermal capacity/conductivity/expansion, frequency-dependent Poisson ratio,
moisture diffusion/permeability, cure and ageing kinetics, adhesion/fracture,
and wavelength-resolved optical constants. Pigment and resin identification
does not manufacture any of these numerical laws.

`cargo test -p xtask --test varnish_materials` checks source compilation,
typed queries, exact coefficients, categorical specimen selection and refusals.
It establishes data transport; audible or historical validation belongs to
independent physical comparisons.
