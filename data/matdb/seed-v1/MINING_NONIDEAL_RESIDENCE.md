# Nonideal plant residence and retained inventory

`fs_material::recovery::residence::nonideal` propagates one conserved passive
constituent through unequal, full tanks with an explicit bypass and a stagnant
region in each tank. It extends the existing equal, ideal-tank hydraulic owner
without turning a tracer response into a leach-recovery law.

This addresses a concrete precursor to FrankenMining plan §17.3 / BD-313:
the metal leaving during a reporting interval can have entered in an earlier
interval, and some metal can remain in poorly exchanged regions. Pairing an
outlet assay with only the current head assay can therefore misstate recovery.

## Physical model and sources

The model families are described in [Fogler, chapter 18, sections 18.7–18.8](https://public.websites.umich.edu/~elements/5e/18chap/Fogler_Web_Ch18_final.pdf).
[Dagadu et al. (2012), *Radiotracer investigation in gold leaching tanks*](https://doi.org/10.1016/j.apradiso.2011.09.003)
reports a tanks-in-series model with exchange between active and stagnant
regions for an aqueous-phase gold-plant tracer investigation. These citations
support the model family. No numerical plant parameters are transcribed or
inferred from them here.

For bypass fraction `a`, the declared total feed splits as:

```text
Q_bypass = a * Q_total
Q_train  = (1 - a) * Q_total
```

Each tank has a positive active volume `Va`, a nonnegative stagnant volume
`Vs`, and a nonnegative exchange flow `qex` in each direction. The two retained
constituent masses obey:

```text
Ca = Ma / Va
Cs = Ms / Vs                 (only when Vs > 0)
dMa/dt = Q_train * (C_in - Ca) - qex * (Ca - Cs)
dMs/dt = qex * (Ca - Cs)
```

The preceding tank's active concentration supplies `C_in` for downstream
tanks. The first tank receives the feed concentration. The bypass carries the
feed directly to the combined outlet with zero modeled delay.

**The flow basis matters.** Active residence time is `Va / Q_train`. With a
bypass, substituting total feed flow for train flow changes the residence
prediction. A source or downstream interface that calls both quantities `Q`
must disambiguate that basis before using the result.

An isolated dead region has `qex = 0` and retains its original mass exactly.
An exchanging stagnant region can release a long tail after a feed pulse.
Exchange can continue when external flow stops, if the caller explicitly
supplies a nonzero exchange flow. When all flows stop, the complete inventory
is unchanged. A zero stagnant volume requires zero stagnant mass and exchange.

## Inputs, outputs and use

`TankGeometry` retains physical volumes. `NonidealTankCascade` retains the
constituent identity, original geometry/initial-state provenance, and separate
active/stagnant inventories. `ResidenceFlows` explicitly declares the total
flow, bypass fraction and one exchange flow per tank for the current interval.

`advance` integrates constant conditions over one interval; conditions can
change between calls. `advance_with_cancel` adds bounded cancellation checks.
The interval reports the inlet, bypass and train volumes, independent inlet
and outlet masses, both directional exchange totals, both inventories and
their signed conservation residuals. The combined outlet mean is absent when
there is no discharged volume; it is not replaced with zero concentration.

```bash
cargo run -p fs-material --example nonideal_residence
```

The runnable example is explicitly authored: three unequal tanks, a feed
pulse, a flow stop with continuing mixing, restart and clean-water flushing.
It reports the retained material during each interval. It contains no
empirical rate calibration or named-site operating prediction.

## Numerical and scientific boundaries

Positive uniformization evaluates the compartment equations and independently
integrates the final active outlet and both exchange fluxes. It does not solve
for the outlet as a balancing remainder or normalize inventories afterward.
Disconnected regions use their exact identity solution. Inputs, physical
quantities, positive reachable paths and conservation are checked before any
state is published.

The evaluator admits 1–64 tanks, with at most 4096 of the largest transition-rate
exposure per call. Excess work, numerical range loss or cancellation leaves
the original state unchanged. The bounded binary64 result is not an interval
enclosure. The direct G0/G1/G3 tests compare unequal-tank and closed-exchange
solutions, independent outlet integrals, the original ideal owner, and
partitioned intervals; they also exercise bypass, storage tails and refusals.

Equal cumulative throughput gives the same history only when all relevant
internal flows scale with the external flow. An exchange rate prescribed in
clock time must not silently be changed to a throughput-volume clock.

The model has no settling, phase slip, changing fill level, chemistry, adsorption
or inferred metal recovery. Fitting hydraulic parameters to tracer data,
deconvolution, and calibration against the named mining-plan upset fixtures
remain separate work. A forward hydraulic example does not close BD-313.
