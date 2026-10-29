# Steckler compartment-fire reference cases

The adjacent `steckler_1982.json` retains six open, steady-state experiments
with a centered methane burner: tests 10, 11, 12, 13, 14 and 16. Every case has
the original 19-height temperature/velocity profile, the published reduced
layer CSV, and the complete authored CFAST input deck in `steckler_1982_raw/`.
The JSON records the upstream repository, commit, path and Git blob SHA for
each copied file. No model was fitted to these observations. Dalmarnock Test
One is excluded from these data and from calibration.

## Sources and interpretation

- Experiment: K. D. Steckler, J. G. Quintiere and W. J. Rinkinen,
  *Flow Induced by Fire in a Compartment*, NBSIR 82-2520 (1982),
  [original report](https://doi.org/10.6028/NBS.IR.82-2520).
  Appendix B contains the measured profiles; the selected cases are
  Tables B-3, B-4, B-5, B-7, B-8 and B-12.
- Numerical experiment data: [firemodels/exp, commit
  704e0d91f15088ae41c6a30a662a7524bf5ff4e8](https://github.com/firemodels/exp/tree/704e0d91f15088ae41c6a30a662a7524bf5ff4e8/Steckler_Compartment).
- Authored comparison inputs: [firemodels/cfast, commit
  1b5ad2cb3cd9db885382d262a27324a4495c96b2](https://github.com/firemodels/cfast/tree/1b5ad2cb3cd9db885382d262a27324a4495c96b2/Validation/Steckler_Compartment).
- The same revision's [experimental description](https://github.com/firemodels/cfast/blob/1b5ad2cb3cd9db885382d262a27324a4495c96b2/Manuals/CFAST_Validation_Guide/Experiment_Chapter.tex)
  corrects the original report's 2.18 m room-height misprint to 2.13 m.

All six decks prescribe 62.9 kW from 0 to 1800 s in a 2.8 m by 2.8 m by
2.13 m room, with a floor-level door 1.83 m high. The original burner is
described as 0.30 m in diameter. CFAST instead authors a 0.09 m² fire area,
a 0.02 m fire-base height, a 0.20 radiative fraction, and 50 MJ/kg heat of
combustion. The fixture preserves these modeling choices. The per-case
ambient temperatures are also taken from the individual CFAST decks.

| Test | Door width (m) | CFAST ambient (°C) | Reduced upper T (°C) | Reduced lower T (°C) | Interface above floor (m) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 10 | 0.24 | 26 | 180.1 | 65.5 | 0.48 |
| 11 | 0.36 | 28 | 153.7 | 56.4 | 0.58 |
| 12 | 0.49 | 22 | 125.9 | 43.8 | 0.67 |
| 13 | 0.62 | 23 | 118.4 | 40.3 | 0.71 |
| 14 | 0.74 | 29 | 115.3 | 44.2 | 0.78 |
| 16 | 0.86 | 26 | 107.4 | 37.8 | 0.82 |

The temperature and interface values above are copied from the upstream
layer files. Celsius-to-kelvin conversion uses +273.15; hot-layer depth is
2.13 m minus the interface height. No profile was redigitized or reduced
again. Original `NaN` and empty fields remain unchanged in the raw CSVs.
The layer CSV's `Time` entries at 1500 and 1800 s are used by the
[CFAST plotting configuration](https://github.com/firemodels/cfast/blob/1b5ad2cb3cd9db885382d262a27324a4495c96b2/Utilities/Python/CFAST_validation_dataplot_inputs.csv)
to draw a constant experimental interface over a comparison interval; they
are not two transient temperature measurements. That configuration uses
the maximum raw room-temperature profile for its temperature scatterplot.
The retained reduced upper-layer mean is a different, explicitly named QoI.

## Wall approximation and units

Every retained CFAST deck uses the same single-layer material on all room
surfaces: conductivity 1 W/(m K), density 200 kg/m³, heat capacity
1000 J/(kg K), thickness 0.013 m, and emissivity 0.9. These are authored
model properties. The experimental construction described in the CFAST
guide has 12.7 mm ceramic fiber insulation over 0.66 mm aluminum on walls
and ceiling, and 19 mm calcium silicate over 12.7 mm plywood on the floor.
The available sources do not supply a complete measured property set for
that assembly. Do not label the deck approximation as measured material
data or silently change it to improve agreement.

For SI conversion, the pinned
[CFAST namelist reader](https://github.com/firemodels/cfast/blob/1b5ad2cb3cd9db885382d262a27324a4495c96b2/Source/CFAST/input_namelist.f90)
assigns `conductivity` directly to the material conductivity and multiplies
`specific_heat` by 1000. Its conductivity declaration is W/(m K).
This resolves a conflicting kW/(m K) label in that revision's user guide.
The deck's heat of combustion is in kJ/kg and HRR is in kW.

## Measurement uncertainty and model comparisons

Appendix B, printed page 27 (PDF page 39), of NBSIR 82-2520 estimates
aspirated room-thermocouple error limits of +1/−3 °C, including radiation.
The ±1 °C wire/recording limit alone excludes installation errors. Neither
is a standard deviation, a confidence level, or a complete uncertainty for
the later upper/lower-layer reduction. The reduced QoI uncertainty remains
`null` in the fixture. The report's 8–50% range for visually inferred plume
entrainment heights is not assigned to these later reduced interfaces.

The [FDS validation uncertainty derivation](https://github.com/firemodels/fds/blob/f3523f3e09b565a5a16fdf36fc40c6ebaeec9789/Manuals/FDS_Validation_Guide/Error_Chapter.tex)
uses positive measured/predicted pairs E and M, the sample mean of
r = ln(M/E), and its sample variance s² (denominator n−1). With a supported
experimental relative standard uncertainty σE, it estimates
σM² = s² − σE² and δ = exp(mean(r) + (σM² − σE²)/2), with the guide's
additional constraint σM ≥ σE. The geometric mean ratio exp(mean(r)) and
combined log scatter can be reported without assigning σE. They must not
be mislabeled as separated model uncertainty or as the adjusted FDS bias
δ. Temperature-rise and hot-layer-depth comparisons must state their
baseline; ratios of Celsius temperatures are not a valid substitute.

`firemodels/out` describes itself as an FDS output repository. Its pinned
Test 10 output at commit `02f6cc3fbcdac0566104abbd7b7691b3a72bc139`
identifies FDS `6.11.1-77-gc444053-master` and starts at 22 °C, whereas the
retained CFAST Test 10 deck specifies 26 °C. It is therefore not included
as an equivalent CFAST reference run. No CFAST simulation output or model
validation pass is claimed by this data fixture.

The intact upstream NIST notice accompanies the raw files. Extracted
metadata and SI conversions were added on 2026-10-29; the source data and
input decks were copied without modification.
