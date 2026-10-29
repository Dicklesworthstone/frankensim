# Soot radiation with wavelength-dependent absorption

This runnable example compares the existing Planck-mean gray calculation with
wavelength-dependent Rayleigh-soot transfer through the same two gas layers.
Every wall and layer emits its own Planck spectrum. The spectral calculation
uses the original reciprocal optical paths, diffuse-gray surface reflections,
prepared linear solve, and signed gas/wall heat accounts.

## Run

```sh
cargo run -p fs-airflow --features compartment-fire \
  --example soot_layer_radiation -- \
  1.75 1.03 1800 0.001 200 2500 \
  0.00018 0.0018 600 1000 600 1000 \
  1 0.25 0.75 300 400 0.8 0.6 authored-demo \
  spectral 1e-7 1e-3 1e-5 1e-5 8192 1e-4
```

The first 20 arguments retain the original gray invocation. `--help` lists
them and their units. The optional seven-argument suffix supplies:

| Argument | Meaning | Example |
|---|---|---:|
| `spectral` | Enable the comparison | — |
| `lambda_min_m` | Lower vacuum wavelength | `1e-7` m |
| `lambda_max_m` | Upper vacuum wavelength | `1e-3` m |
| `heat_tol_W` | Absolute per-owner heat quadrature tolerance | `1e-5` W |
| `radiosity_tol_W_m2` | Absolute per-surface radiosity quadrature tolerance | `1e-5` W/m² |
| `max_evaluations` | Wavelength work budget | `8192` |
| `max_omitted_fraction` | Largest admitted omitted blackbody fraction | `1e-4` |

The two supplied spectral-weighting temperatures affect only the gray
comparison. Actual gas and wall temperatures determine every spectral emitter.
Set both soot concentrations to zero for the transparent comparison.

## What changes physically

A gray Planck mean compresses wavelength-dependent absorption into one
coefficient. Its weighting temperature describes one assumed incident or
emitted spectrum. When surfaces and layers have different temperatures, and
absorption filters those spectra along the path, one coefficient cannot
preserve every exchange. The new calculation evaluates the original
`RayleighSoot::spectral_coefficient` at each integration wavelength.

The reusable API also accepts explicit molecular/background absorption in
contiguous wavelength bands. This example declares that contribution zero;
it does not invent molecular opacity from gas names. For coupled room/wall
calculations, use `evaluate_room_spectral`, add its gas powers once through
`add_gas_heat`, and apply the matching signed surface powers to their wall
owners in the same accepted network interval.

## Retained run

Actual native results for the command above are retained in
[`spectral_soot_example.json`](../../../../data/reference/compartment_fire/spectral_soot_example.json).
Positive power enters the named thermal owner.

| Owner | Spectral heat, W | Gray heat, W |
|---|---:|---:|
| Lower gas layer | +2,191.714109488 | +1,172.641529627 |
| Upper gas layer | −59,957.24708007 | −64,059.11712089 |
| Wall 0 | +33,499.61123510 | +36,682.77795559 |
| Wall 1 | +24,265.92173548 | +26,203.69763567 |

The spectral run used 288 wavelength evaluations/factorizations and 11 accepted
panels. Its enclosure energy residual was −7.28e−12 W, the integrated surface
equation residual was 1.59e−12 W/m², and the largest reported heat quadrature
estimate was 9.45e−7 W. The largest omitted blackbody fraction was 5.57e−6.
These are numerical checks, not validation against a fire experiment. The
focused native radiation suite passed 22 tests, including six spectral groups;
this is not a full-workspace or DSR result.

## Assumptions and limits

All inputs in this invocation are authored demonstration values. In particular,
`n=1.75`, absorptive `k=1.03`, and particle density 1800 kg/m³ are not promoted
to universal or validated soot constants. A real application must supply an
optical source supporting its indices, loading, wavelength and temperature
range, small-particle limit and negligible scattering.

The geometry uses an infinite-parallel-plate view factor with one reciprocal
normal chord. This conservative angular approximation is not the exact
hemispherical slab solution. Surface emissivities remain explicitly gray;
particles and gas share each layer's temperature. No soot production, transport,
deposition, gas-pressure contribution or molecular spectroscopy is predicted.

The positive 8/16-point quadrature differences are error estimates, not certified
bounds. Omitted blackbody emission is reported separately and is not a bound
on omitted net heat. Spectral band/data and angular errors remain separate
from the reported numerical residuals. Constant band coefficients can reuse
one factorization per band; this varying-soot example needs one per wavelength.

Planck radiation follows [NIST SP 250-1, §2](https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nbsspecialpublication250-1.pdf),
with the original conduction module's Stefan--Boltzmann normalization.
The small-argument tail expansion uses [DLMF §24.2.1](https://dlmf.nist.gov/24.2.E1).
