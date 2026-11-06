//! Spontaneous radiative recombination with escaping photons in a rigid parcel.
//!
//! The reciprocal L1 source supplies both alpha=<v sigma_RR> and the captured
//! kinetic-energy moment L. For positive capture extent xi [mol], electron
//! sensible energy loses NA*(L/alpha)*d_xi and chemical storage loses chi*d_xi.
//! Their sum leaves as photons; chemical release is not deposited as heat.
//! Preferential capture of slow electrons can RAISE the surviving electron
//! temperature even while electron sensible energy decreases.
//!
//! Each substep predicts its midpoint temperature, freezes both source moments
//! there, and solves binary reactant depletion exactly. This midpoint method
//! is second order for smooth source coefficients; full/two-half-step comparison
//! controls its local error. Initial rates are shared and source moments are
//! refreshed at actual/predicted stage temperatures. No trial allocates a species
//! vector. Publication retains the original bank and unchanged heavy temperature.
//!
//! This is spontaneous capture into the one bound state represented by the
//! photon source, with immediate electron thermalization and escaping photons.
//! It does not supply stimulated emission, photon trapping, total Case A/B
//! recombination, excited levels, recoil, or a closed-gas entropy inequality.
//! Local integration estimates and conditional source-tail bounds are separate;
//! neither is a global solution-error or empirical model-error certificate.

use fs_math::det::{exp, expm1, ln};
use fs_qty::{Amount, Energy, Temperature, Time};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::electron_collision::{AVOGADRO_PER_MOL as NA, MaxwellianRateControlsV1};
use fs_thermochem::equilibrium::composition::{
    FrozenIonizedGasStateErrorV1, FrozenIonizedGasStateV1,
};
use fs_thermochem::kinetics::atomic_ionization::radiative_recombination::{
    AtomicRadiativeRecombinationV1, RadiativeRecombinationErrorV1,
    RadiativeRecombinationEvaluationV1,
};
use fs_thermochem::kinetics::{KineticsErrorV1, admit_caloric_interval};

/// Electron source interval, independent quadrature and time-step allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RadiativeRecombinationSolve {
    /// Supported electron temperature lower bound [K].
    pub minimum_electron_temperature_k: f64,
    /// Supported electron temperature upper bound [K].
    pub maximum_electron_temperature_k: f64,
    /// Local relative amount/temperature error allowance, in (0,0.1].
    pub relative_tolerance: f64,
    /// Local absolute amount error scale [mol].
    pub absolute_amount_tolerance_mol: f64,
    /// Local absolute electron temperature error scale [K].
    pub temperature_tolerance_k: f64,
    /// Final original-source gas plus escaped-radiation energy allowance [J].
    pub energy_tolerance_j: f64,
    /// Charge residual scaled by actual carriers, in (0,1).
    pub relative_charge_tolerance: f64,
    /// First proposed step [s].
    pub initial_step_s: f64,
    /// Smallest retry [s]; a shorter final remainder is allowed.
    pub minimum_step_s: f64,
    /// Largest proposed step [s].
    pub maximum_step_s: f64,
    /// Bound on accepted plus rejected attempts.
    pub maximum_steps: usize,
    /// Rate and energy-moment quadrature/tail qualification at every temperature.
    pub quadrature: MaxwellianRateControlsV1,
}
impl RadiativeRecombinationSolve {
    /// Source bounds are explicit; numerical defaults are not physical accuracy.
    pub fn new(minimum: Temperature, maximum: Temperature) -> Self {
        Self {
            minimum_electron_temperature_k: minimum.value(),
            maximum_electron_temperature_k: maximum.value(),
            relative_tolerance: 1e-5,
            absolute_amount_tolerance_mol: 1e-14,
            temperature_tolerance_k: 1e-6,
            energy_tolerance_j: 1e-7,
            relative_charge_tolerance: 1e-10,
            initial_step_s: 1e-8,
            minimum_step_s: 1e-24,
            maximum_step_s: 1e-3,
            maximum_steps: 100_000,
            quadrature: MaxwellianRateControlsV1 {
                relative_tolerance: 1e-10,
                maximum_subdivisions: 4096,
            },
        }
    }
    fn validate(self) -> Result<(), RadiativeRecombinationError> {
        if [
            self.minimum_electron_temperature_k,
            self.maximum_electron_temperature_k,
            self.relative_tolerance,
            self.absolute_amount_tolerance_mol,
            self.temperature_tolerance_k,
            self.energy_tolerance_j,
            self.relative_charge_tolerance,
            self.initial_step_s,
            self.minimum_step_s,
            self.maximum_step_s,
        ]
        .iter()
        .any(|x| !positive(*x))
            || self.minimum_electron_temperature_k >= self.maximum_electron_temperature_k
            || self.relative_tolerance > 0.1
            || self.relative_charge_tolerance >= 1.0
            || self.minimum_step_s > self.initial_step_s
            || self.initial_step_s > self.maximum_step_s
            || !(1..=1_000_000).contains(&self.maximum_steps)
        {
            return Err(RadiativeRecombinationError::Invalid(
                "radiative recombination controls",
            ));
        }
        Ok(())
    }
}

/// Source or evolution refusal; incoming state and source remain immutable.
#[derive(Clone, Debug, PartialEq)]
pub enum RadiativeRecombinationError {
    /// Invalid controls, exposure or original species binding.
    Invalid(&'static str),
    /// Original reciprocal source, numerical moment or tail qualification refused.
    Source(RadiativeRecombinationErrorV1),
    /// Original species cannot support the authored caloric path.
    Admission(KineticsErrorV1),
    /// Positive electron temperature left the supported interval.
    OutsideDomain,
    /// Bounded time-step work or a representable forward time increment exhausted.
    StepBudget,
    /// Required arithmetic is not representable.
    Arithmetic(&'static str),
    /// Actual sensible-energy change fails the escaped-energy account.
    EnergyBalance { residual_j: f64 },
    /// Original common charged-state publication refused.
    State(FrozenIonizedGasStateErrorV1),
    /// Cancellation observed before a complete endpoint.
    Cancelled,
}
impl core::fmt::Display for RadiativeRecombinationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "radiative recombination gas evolution: {self:?}")
    }
}
impl std::error::Error for RadiativeRecombinationError {}

/// One consistent population, gas-energy and escaped-photon account.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RadiativeRecombinationBalance {
    /// One escaping photon per captured ion/electron pair.
    pub emitted_photons: f64,
    /// Original formation-energy change, negative for recombination [J].
    pub chemical_storage_change_j: f64,
    /// Actual electron sensible-energy change [J].
    pub electron_sensible_energy_change_j: f64,
    /// Integrated kinetic energy of captured electrons, nonnegative [J].
    pub captured_electron_energy_j: f64,
    /// Released chemical plus captured kinetic energy in escaping photons [J].
    pub emitted_radiation_energy_j: f64,
    /// Actual gas-energy change plus escaping radiation [J].
    pub energy_residual_j: f64,
    /// Largest accepted full/two-half-step local error divided by its allowance.
    pub maximum_normalized_local_error: f64,
    /// Largest relative quadrature estimate encountered in evaluated coefficients.
    pub maximum_relative_quadrature_error: f64,
    /// Largest conditional omitted-tail bound relative to an evaluated moment.
    pub maximum_relative_tail_bound: f64,
    /// Complete accepted steps, each retaining the two-half-step endpoint.
    pub accepted_steps: usize,
    /// Rejected attempts.
    pub rejected_steps: usize,
}

/// Original-source endpoint reusable by existing collision and material owners.
#[derive(Clone, Debug, PartialEq)]
pub struct RadiativeRecombinationResult {
    initial: FrozenIonizedGasStateV1,
    state: FrozenIonizedGasStateV1,
    source: AtomicRadiativeRecombinationV1,
    controls: RadiativeRecombinationSolve,
    elapsed: Time,
    extent: Amount,
    balance: RadiativeRecombinationBalance,
}
impl RadiativeRecombinationResult {
    /// Accepted original-bank gas state, without an LTE projection.
    pub const fn state(&self) -> &FrozenIonizedGasStateV1 {
        &self.state
    }
    /// Original input state.
    pub const fn initial_state(&self) -> &FrozenIonizedGasStateV1 {
        &self.initial
    }
    /// Original photon source and declared cross-section envelope.
    pub const fn source(&self) -> &AtomicRadiativeRecombinationV1 {
        &self.source
    }
    /// Applied source, quadrature and time-step allowances.
    pub const fn controls(&self) -> RadiativeRecombinationSolve {
        self.controls
    }
    /// Completed physical duration [s].
    pub const fn elapsed_time(&self) -> Time {
        self.elapsed
    }
    /// Nonnegative recombination extent [mol], one ion and electron removed per event.
    pub const fn recombination_extent(&self) -> Amount {
        self.extent
    }
    /// Population, energy and numerical work accounts.
    pub const fn balance(&self) -> RadiativeRecombinationBalance {
        self.balance
    }
    /// Escaping photon energy [J].
    pub const fn emitted_radiation_energy(&self) -> Energy {
        Energy::new(self.balance.emitted_radiation_energy_j)
    }
}

#[derive(Clone, Copy)]
struct Point {
    minor: f64,
    extent: f64,
    kinetic: f64,
    te: f64,
}

/// Advance spontaneous radiative recombination at fixed volume.
///
/// The source resolves one final bound state; no total recombination coefficient
/// is inferred. All emitted photons escape. Source intervals and both L1 moment
/// tail bounds are rechecked as electron temperature changes. Exact absence is
/// an identity, but original-source and numerical-control admission still apply.
pub fn advance_radiative_recombination(
    initial: &FrozenIonizedGasStateV1,
    source: &AtomicRadiativeRecombinationV1,
    duration: Time,
    controls: RadiativeRecombinationSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<RadiativeRecombinationResult, RadiativeRecombinationError> {
    poll(&mut cancelled)?;
    controls.validate()?;
    let channel = source.source().channel();
    let end = duration.value();
    if !nonnegative(end)
        || initial.thermodynamics() != channel.thermodynamics()
        || initial.electron_index() != Some(channel.electron_index())
        || initial.relative_charge_residual() > controls.relative_charge_tolerance
    {
        return Err(RadiativeRecombinationError::Invalid(
            "duration or original atomic bank",
        ));
    }
    let [a, i, e] = [
        channel.neutral_index(),
        channel.ion_index(),
        channel.electron_index(),
    ];
    let initial_te = initial.electron_temperature().value();
    if initial_te < controls.minimum_electron_temperature_k
        || initial_te > controls.maximum_electron_temperature_k
    {
        return Err(RadiativeRecombinationError::OutsideDomain);
    }
    for index in [a, i, e] {
        let (lo, hi) = if index == e {
            (
                Temperature::new(controls.minimum_electron_temperature_k),
                Temperature::new(controls.maximum_electron_temperature_k),
            )
        } else {
            (initial.heavy_temperature(), initial.heavy_temperature())
        };
        admit_caloric_interval(
            &channel.thermodynamics().models()[index],
            lo,
            hi,
            &mut cancelled,
        )
        .map_err(|error| match error {
            KineticsErrorV1::Cancelled => RadiativeRecombinationError::Cancelled,
            error => RadiativeRecombinationError::Admission(error),
        })?;
    }
    let ne = initial.amounts()[e].value();
    let ni = initial.amounts()[i].value();
    let difference = ne - ni;
    let electron_excess = difference.max(0.0);
    let volume = initial.volume().value();
    let mut point = Point {
        minor: ne.min(ni),
        extent: 0.0,
        kinetic: 0.0,
        te: initial_te,
    };
    let mut balance = RadiativeRecombinationBalance::default();
    let mut rates = evaluate(source, point.te, controls, &mut balance, &mut cancelled)?;
    let mut time = 0.0;
    let mut step = controls.initial_step_s;
    while time < end && point.minor > 0.0 && rates.coefficient_m3_s > 0.0 {
        poll(&mut cancelled)?;
        if balance.accepted_steps + balance.rejected_steps >= controls.maximum_steps {
            return Err(RadiativeRecombinationError::StepBudget);
        }
        let dt = step.min(end - time);
        if time + dt <= time {
            return Err(RadiativeRecombinationError::StepBudget);
        }
        let trial = (|| {
            let full = midpoint_step(
                point,
                rates,
                source,
                difference,
                volume,
                dt,
                controls,
                &mut balance,
                &mut cancelled,
            )?;
            let half = midpoint_step(
                point,
                rates,
                source,
                difference,
                volume,
                0.5 * dt,
                controls,
                &mut balance,
                &mut cancelled,
            )?;
            let midpoint = evaluate(source, half.te, controls, &mut balance, &mut cancelled)?;
            let next = midpoint_step(
                half,
                midpoint,
                source,
                difference,
                volume,
                0.5 * dt,
                controls,
                &mut balance,
                &mut cancelled,
            )?;
            let amount_scale = controls.absolute_amount_tolerance_mol
                + controls.relative_tolerance * point.minor.max(next.minor);
            let temperature_scale = controls.temperature_tolerance_k
                + controls.relative_tolerance * full.te.abs().max(next.te.abs());
            let error = ((full.minor - next.minor).abs() / amount_scale)
                .max((full.extent - next.extent).abs() / amount_scale)
                .max((full.te - next.te).abs() / temperature_scale);
            if !nonnegative(error) {
                return Err(RadiativeRecombinationError::Arithmetic(
                    "local step estimate",
                ));
            }
            Ok((next, error))
        })();
        match trial {
            Ok((next, error)) if error <= 1.0 => {
                point = next;
                time += dt;
                balance.accepted_steps += 1;
                balance.maximum_normalized_local_error =
                    balance.maximum_normalized_local_error.max(error);
                // Also qualify the final endpoint's moments, even at end time.
                rates = evaluate(source, point.te, controls, &mut balance, &mut cancelled)?;
                let factor = if error == 0.0 {
                    2.0
                } else {
                    (0.85 * exp(-ln(error) / 3.0)).clamp(0.25, 2.0)
                };
                step = (dt * factor)
                    .min(controls.maximum_step_s)
                    .max(controls.minimum_step_s);
            }
            Ok((_, error)) => {
                balance.rejected_steps += 1;
                step = dt * (0.8 * exp(-ln(error) / 3.0)).clamp(0.1, 0.5);
                if step < controls.minimum_step_s {
                    return Err(RadiativeRecombinationError::StepBudget);
                }
            }
            Err(RadiativeRecombinationError::OutsideDomain) => {
                balance.rejected_steps += 1;
                step = 0.5 * dt;
                if step < controls.minimum_step_s {
                    return Err(RadiativeRecombinationError::OutsideDomain);
                }
            }
            Err(error) => return Err(error),
        }
    }
    let remaining_e = point.minor + electron_excess;
    let sensible = 1.5 * R * (ne * (point.te - initial_te) + (remaining_e - ne) * point.te);
    // Conserved elemental/charge reference shifts cancel before subtraction.
    let chemical = -channel.ionization_energy_j_mol() * point.extent;
    let emitted = -chemical + point.kinetic;
    let residual = sensible + point.kinetic;
    if !nonnegative(emitted)
        || !nonnegative(point.kinetic)
        || !sensible.is_finite()
        || !residual.is_finite()
    {
        return Err(RadiativeRecombinationError::Arithmetic(
            "escaping radiation energy",
        ));
    }
    if residual.abs() > controls.energy_tolerance_j {
        return Err(RadiativeRecombinationError::EnergyBalance {
            residual_j: residual,
        });
    }
    let state = if point.extent == 0.0 {
        initial.clone()
    } else {
        let mut amounts = initial.amounts().to_vec();
        amounts[a] = Amount::new(amounts[a].value() + point.extent);
        amounts[i] = Amount::new(point.minor + (-difference).max(0.0));
        amounts[e] = Amount::new(remaining_e);
        FrozenIonizedGasStateV1::new(
            channel.thermodynamics(),
            &amounts,
            initial.volume(),
            initial.heavy_temperature(),
            Temperature::new(point.te),
            controls.relative_charge_tolerance,
            &mut cancelled,
        )
        .map_err(|error| match error {
            FrozenIonizedGasStateErrorV1::Cancelled => RadiativeRecombinationError::Cancelled,
            error => RadiativeRecombinationError::State(error),
        })?
    };
    balance.emitted_photons = point.extent * NA;
    balance.chemical_storage_change_j = chemical;
    balance.electron_sensible_energy_change_j = sensible;
    balance.captured_electron_energy_j = point.kinetic;
    balance.emitted_radiation_energy_j = emitted;
    balance.energy_residual_j = residual;
    poll(&mut cancelled)?;
    Ok(RadiativeRecombinationResult {
        initial: initial.clone(),
        state,
        source: source.clone(),
        controls,
        elapsed: duration,
        extent: Amount::new(point.extent),
        balance,
    })
}

fn evaluate(
    source: &AtomicRadiativeRecombinationV1,
    te: f64,
    controls: RadiativeRecombinationSolve,
    balance: &mut RadiativeRecombinationBalance,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<RadiativeRecombinationEvaluationV1, RadiativeRecombinationError> {
    let result = source
        .evaluate(Temperature::new(te), controls.quadrature, cancelled)
        .map_err(|error| match error {
            RadiativeRecombinationErrorV1::Cancelled => RadiativeRecombinationError::Cancelled,
            error => RadiativeRecombinationError::Source(error),
        })?;
    balance.maximum_relative_quadrature_error = balance
        .maximum_relative_quadrature_error
        .max(result.relative_quadrature_error_estimate);
    balance.maximum_relative_tail_bound = balance
        .maximum_relative_tail_bound
        .max(result.relative_tail_bound);
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
fn midpoint_step(
    point: Point,
    initial_rates: RadiativeRecombinationEvaluationV1,
    source: &AtomicRadiativeRecombinationV1,
    difference: f64,
    volume: f64,
    dt: f64,
    controls: RadiativeRecombinationSolve,
    balance: &mut RadiativeRecombinationBalance,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Point, RadiativeRecombinationError> {
    let predictor = frozen_step(point, initial_rates, difference, volume, 0.5 * dt, controls)?;
    let midpoint = evaluate(source, predictor.te, controls, balance, cancelled)?;
    frozen_step(point, midpoint, difference, volume, dt, controls)
}

fn frozen_step(
    point: Point,
    rates: RadiativeRecombinationEvaluationV1,
    difference: f64,
    volume: f64,
    dt: f64,
    controls: RadiativeRecombinationSolve,
) -> Result<Point, RadiativeRecombinationError> {
    if point.minor == 0.0 || rates.coefficient_m3_s == 0.0 || dt == 0.0 {
        return Ok(point);
    }
    let mean = rates
        .mean_captured_energy_j()
        .ok_or(RadiativeRecombinationError::Arithmetic(
            "captured electron energy",
        ))?;
    let maximum = point.minor + difference.abs();
    let dose = ((rates.coefficient_m3_s * NA) * maximum / volume) * dt;
    if !nonnegative(dose) {
        return Err(RadiativeRecombinationError::OutsideDomain);
    }
    if dose == 0.0 {
        return Err(RadiativeRecombinationError::Arithmetic(
            "positive recombination dose underflow",
        ));
    }
    let (captured, minor) = binary_depletion(point.minor, maximum, dose)?;
    let remaining_e = minor + difference.max(0.0);
    if !positive(remaining_e) {
        return Err(RadiativeRecombinationError::OutsideDomain);
    }
    let equivalent_temperature = (mean * NA) / (1.5 * R);
    let te = point.te + (captured / remaining_e) * (point.te - equivalent_temperature);
    if !te.is_finite()
        || te < controls.minimum_electron_temperature_k
        || te > controls.maximum_electron_temperature_k
    {
        return Err(RadiativeRecombinationError::OutsideDomain);
    }
    let next = Point {
        minor,
        extent: point.extent + captured,
        kinetic: point.kinetic + (captured * NA) * mean,
        te,
    };
    if !nonnegative(next.extent) || !nonnegative(next.kinetic) {
        return Err(RadiativeRecombinationError::Arithmetic(
            "integrated capture",
        ));
    }
    Ok(next)
}

// Exact fixed-coefficient d(minor)/dt=-K*minor*(minor+constant_difference).
// dose=K*maximum_initial*dt. The equal-reactant limit is continuous. Compute
// the smaller side directly to preserve tiny captures and surviving tails.
fn binary_depletion(
    minor: f64,
    maximum: f64,
    dose: f64,
) -> Result<(f64, f64), RadiativeRecombinationError> {
    let ratio = (maximum - minor) / maximum;
    let z = dose * ratio;
    let phi = if z == 0.0 { 1.0 } else { -expm1(-z) / z };
    let a = dose * phi;
    let denominator = 1.0 + (minor / maximum) * a;
    let fraction = a / denominator;
    let (captured, remaining) = if fraction <= 0.5 {
        let captured = minor * fraction;
        (captured, minor - captured)
    } else {
        let remaining = minor * (exp(-z) / denominator);
        (minor - remaining, remaining)
    };
    if !nonnegative(captured) || captured > minor || !nonnegative(remaining) || remaining > minor {
        return Err(RadiativeRecombinationError::Arithmetic(
            "binary reactant depletion",
        ));
    }
    Ok((captured, remaining))
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn nonnegative(x: f64) -> bool {
    x.is_finite() && x >= 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), RadiativeRecombinationError> {
    if cancelled() {
        Err(RadiativeRecombinationError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g1_binary_capture_preserves_trace_limits_and_unequal_reactants() {
        for (minor, maximum, dose, capture, remaining) in [
            (1.0, 1.0, 1.0, 0.5, 0.5),
            (1.0, 1.0, 1e-12, 9.99999999999e-13, 0.999999999999),
            (1.0, 1.0, 1e200, 1.0, 1e-200),
            (1.0, 2.0, 2.0, 0.7746003264394359, 0.2253996735605641),
        ] {
            let (a, b) = binary_depletion(minor, maximum, dose).unwrap();
            assert!((a - capture).abs() <= 3e-14 * capture);
            assert!((b - remaining).abs() <= 3e-14 * remaining);
            assert!(a >= 0.0 && b >= 0.0 && a <= minor && b <= minor);
        }
        let (_, first) = binary_depletion(1.0, 2.0, 0.8).unwrap();
        let (_, split) = binary_depletion(first, first + 1.0, 0.6 * (first + 1.0)).unwrap();
        let (_, whole) = binary_depletion(1.0, 2.0, 2.0).unwrap();
        assert!((split - whole).abs() < 3e-15);
    }
}
