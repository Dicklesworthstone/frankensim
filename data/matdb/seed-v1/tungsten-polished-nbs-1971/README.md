# Polished tungsten: NBS pulse-heating measurements, 1971

This pack supplies existing electrical and caloric consumers with one characterized
hot tungsten specimen, alongside separately named radiation measurements. It retains five claims
from [Cezairliyan and McClure, NBS 75A, 283–290 (1971)](https://doi.org/10.6028/jres.075A.027),
using the [official primary PDF](https://nvlpubs.nist.gov/nistpubs/jres/75A/jresv75An4p283_A1b.pdf).

| Property key | Support | Source representation |
| --- | --- | --- |
| `specific-heat-capacity-cp` | 2000–3600 K | All 17 Table 2 molar ordinates divided by the paper's 0.18385 kg/mol; explicit linear interpolation between printed rows |
| `electrical-resistivity-reference-geometry-298k` | 2000–3600 K | Original affine Eq. 2, in Ohm m; **298 K reference geometry** |
| `hemispherical-total-emissivity-reference-area-298k` | 2300–3400 K | Original affine Eq. 3; **298 K reference emitting area** |
| `normal-spectral-emissivity` | 2000–3000 K | Original affine Eq. 4; only the reported 650 nm effective wavelength and normal direction |
| `density` | 293 K only | The specimen's room-temperature measurement, 19230 kg/m³ |

The temperature coordinate is named **`T_IPTS68`** deliberately. The source used
the 1968 International Practical Temperature Scale; this pack makes no conversion
to ITS-90 or an exact thermodynamic-temperature scale. Hot-property queries also
require **`p=0.013332236842105263` Pa**, the SI conversion of the reported nominal
vacuum of approximately (10^{-4}) torr. The conversion digits do not add
measurement accuracy, and the single pressure point establishes no pressure law.

## Specimen and approximation boundaries

The specimen was a tungsten tube machined from rod by electroerosion, polished
externally and conditioned by about 30 heating pulses to 3200 K. The paper
reports post-experiment impurities between 450 and 740 ppm by mass and an initial
residual-resistivity ratio of 41. These data identify that specimen and treatment;
they do not qualify a commercial wire grade, coating or arbitrary surface finish.

The Cp ordinates are published evaluations of a fit to measurements. Their
piecewise-linear interpolation is an engineering representation of the rounded
table, and the native caloric owner integrates that representation exactly. It
does not silently substitute a constant Cp or claim exact evaluation of the
source's cubic fit. The resistivity and emittance fits are affine, so two endpoint
knots retain their published temperature dependence without table-rounding loss.

The authors estimate inaccuracies of 2–3% for Cp, 1% for resistivity and 3% for
both emittances. Their separately stated fit standard deviations have a different
meaning. Neither supplies a confidence level or covariance; the normalized
claims retain `Unstated` uncertainty and preserve those estimates in observations.

Total emittance below 2300 K is marked extrapolated in Table 2 and is excluded.
The 650 nm measurement used a 10 nm instrument bandwidth; it does not establish
a constant spectrum across that bandwidth or the visible band. Its additional
query coordinates are `effective-wavelength=6.5e-7` m and `direction-cosine=1`.
The paper does not establish a vacuum-wavelength conversion. Total and normal
spectral emittance remain separate quantities; no diffuse or gray optical model
is inferred from either measurement.

There is no cold-start trajectory, high-temperature density, thermal conductivity,
expansion, melting law or broad spectral-color calibration in this pack.
Reference resistivity must use the original reference L/A convention; total
emittance must use the reference emitting area. Applying either to independently
expanded geometry would change the source's physical meaning.

## Existing consumers

The source-specific electrical key must be selected explicitly through
`resolve_conductor_state_with_property`, then retained by `MaterialOhmicNode`
with fixed reference dimensions on the explicit `T_IPTS68` axis. The coefficient
must not be paired with independently changing geometry. The generic
`electrical-resistivity` lookup has no claim in this pack.

`IsobaricHeatCapacityLaw` integrates and inverts the Cp source without a
density-temperature law or conductivity. The reference-area emittance remains
available through the original typed claim query. Generic
`SurfaceEmissivity::from_condition` refuses: no intrinsic hot-area
`hemispherical-total-emissivity` is supplied. A consumer must explicitly retain
the reference emitting area before using the measured coefficient. These are
independent source uses; a complete thermal or optical simulation must provide
the additional physics and compatible conditions it actually consumes.

```bash
cargo run -p xtask -- matdb-catalog \
  --pack tungsten-polished-nbs-1971 \
  --property electrical-resistivity-reference-geometry-298k \
  --property specific-heat-capacity-cp \
  --property hemispherical-total-emissivity-reference-area-298k \
  --range T_IPTS68=2300:3400 \
  --at p=0.013332236842105263

cargo test -p xtask --test material_electrical_sources tungsten -- --nocapture
```

The focused tests use the actual compiler, native source owners and independent
source values. They check signed Joule power, Cp integration/inversion, retained
conditions, and unsupported temperature, pressure, wavelength and direction.

The NBS authors' work is a U.S. Government publication. NIST's
[Journal of Research reuse statement](https://www.nist.gov/nist-research-library/journal-research-nist/about-journal)
allows reuse with attribution; the source manifest retains that decision and the
original citation.
