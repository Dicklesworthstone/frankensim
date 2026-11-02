//! Nonlinear native enthalpy evolution under prescribed interference power.
//!
//! The original body's phase curve, mass, transport and convection/radiation
//! helpers remain authoritative. Each accepted backward-Euler step consumes
//! the complete source's analytic interval energy, without replacing its
//! spectrum by a sampled or carrier-averaged power.
//!
//! A separate analytic primitive-ripple bound limits how far the true
//! intra-step input can depart from its interval-mean replacement. This
//! prevents full carrier cycles from fooling endpoint or step-doubling tests.
//! It is a source-energy bound, not a certified nonlinear temperature error.
//! Maximum thermal step and native enthalpy/root-residual tolerances remain
//! independent numerical controls.
//!
//! The spectrum is a prescribed applied-thermal-power declaration. Signed
//! extraction is allowed; reversible field storage is not inferred to be
//! heat. An ensemble-mean source evolved through nonlinear thermal laws is a
//! trajectory driven by that mean source, NOT an ensemble-mean temperature.
//! Temperature-dependent electromagnetic feedback requires an explicit
//! rebuild by the original field owner; this march never invents it.

use fs_fft::harmonic_window::{PreparedHarmonicProductSpectrum, ProductSpectrumError, WindowError};
use fs_material::phase::EquilibriumPhaseState;
use fs_math::dd::Dd;

use super::super::{
    BiotGate, LumpedEnthalpyBody, LumpedEnthalpyMarchConfig, MAX_LUMPED_ENTHALPY_STEPS,
    bisect_enthalpy_step, enthalpy_power, lumped_error, phase_state, validate_enthalpy_config,
};
use crate::ConductionError;

/// Hard bound on ripple-admission trials plus implicit thermal steps.
pub const MAX_ENTHALPY_SPECTRUM_ATTEMPTS: usize = 10_000_000;
/// Hard bound on tone/DC visits across all ripple and causal-integral calls.
pub const MAX_ENTHALPY_SPECTRUM_WORK: usize = 268_435_456;

/// Explicit application of one prescribed power spectrum to an original body.
///
/// A W spectrum uses scale 1. A spatially uniform W/m³ spectrum may use an
/// explicitly justified volume. Labels and conversion are caller declarations;
/// no absorption, spatial homogenization, source phase or source statistics
/// are inferred. A Biot gate alone does not justify homogenizing rapid local
/// deposition or establish an initially isothermal body.
#[derive(Clone, Copy, Debug)]
pub struct LumpedEnthalpyPowerBinding<'a> {
    pub body: &'a LumpedEnthalpyBody<'a>,
    pub spectrum: &'a PreparedHarmonicProductSpectrum,
    pub input_unit: &'a str,
    pub watts_per_spectrum_unit: f64,
    pub source_clock: &'a str,
}

/// Source admission, numerical residual and resource controls.
///
/// The native config supplies the maximum THERMAL step, initial enthalpy,
/// reservoirs, maximum accepted steps and bisection tolerance. Its constant
/// internal_power_w must be zero: include all DC power in the bound spectrum.
#[derive(Clone, Copy, Debug)]
pub struct LumpedEnthalpySpectrumControls {
    /// Absolute origin in the declared source clock. Elapsed time is kept in
    /// DD when added to this origin, preserving late short substeps.
    pub start_s: f64,
    /// Refinement below this positive step refuses. A final remainder can be
    /// shorter; it inherits the already-admitted larger step's ripple bound.
    pub minimum_step_s: f64,
    /// Maximum analytic intra-step input primitive ripple [J], nonnegative.
    /// This is separate from nonlinear integration and root-solve errors.
    pub source_ripple_tolerance_j: f64,
    /// Maximum absolute discrete energy closure error of each accepted step [J].
    pub maximum_step_energy_residual_j: f64,
    pub maximum_attempts: usize,
    pub maximum_spectral_work: usize,
    /// Includes the initial sample. Accepted states are never silently skipped.
    pub maximum_output_samples: usize,
    /// Logical retained sample-vector payload; borrowed bodies/spectra,
    /// source-owned data, stack values and allocator overhead are excluded.
    pub maximum_workspace_bytes: usize,
}

/// Signed discrete energy accounting; positive boundary exchange is INWARD.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LumpedEnthalpySpectrumEnergy {
    source_input_j: Dd,
    convection_into_body_j: Dd,
    radiation_into_body_j: Dd,
    stored_energy_change_j: Dd,
    residual_j: Dd,
}

impl LumpedEnthalpySpectrumEnergy {
    const ZERO: Self = Self {
        source_input_j: Dd::ZERO,
        convection_into_body_j: Dd::ZERO,
        radiation_into_body_j: Dd::ZERO,
        stored_energy_change_j: Dd::ZERO,
        residual_j: Dd::ZERO,
    };

    /// Analytic integral of the actual prescribed spectrum times its W/unit scale.
    pub const fn source_input_j_dd(self) -> Dd {
        self.source_input_j
    }

    /// Native endpoint convection times the accepted duration (backward Euler).
    pub const fn convection_into_body_j_dd(self) -> Dd {
        self.convection_into_body_j
    }

    /// Native endpoint radiation times the accepted duration (backward Euler).
    pub const fn radiation_into_body_j_dd(self) -> Dd {
        self.radiation_into_body_j
    }

    /// Original mass times the actual accepted change in specific enthalpy.
    pub const fn stored_energy_change_j_dd(self) -> Dd {
        self.stored_energy_change_j
    }

    /// Storage minus source, convection and radiation; not a temporal-error bound.
    pub const fn residual_j_dd(self) -> Dd {
        self.residual_j
    }
}

/// One actual initial/accepted state; no trial state appears in the trajectory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LumpedEnthalpySpectrumSample {
    elapsed_s: Dd,
    step_duration_s: f64,
    phase_state: EquilibriumPhaseState,
    convection_into_body_w: f64,
    radiation_into_body_w: f64,
    source_ripple_bound_j: f64,
    energy: LumpedEnthalpySpectrumEnergy,
}

impl LumpedEnthalpySpectrumSample {
    pub const fn elapsed_s_dd(self) -> Dd {
        self.elapsed_s
    }

    pub const fn step_duration_s(self) -> f64 {
        self.step_duration_s
    }

    pub const fn phase_state(self) -> EquilibriumPhaseState {
        self.phase_state
    }

    pub const fn convection_into_body_w(self) -> f64 {
        self.convection_into_body_w
    }

    pub const fn radiation_into_body_w(self) -> f64 {
        self.radiation_into_body_w
    }

    /// Conservative analytic input-ripple bound for this step. The final
    /// shorter step may retain the admitted maximum step's looser bound.
    pub const fn source_ripple_bound_j(self) -> f64 {
        self.source_ripple_bound_j
    }

    pub const fn energy(self) -> LumpedEnthalpySpectrumEnergy {
        self.energy
    }
}

/// Complete atomic trajectory on the original body and phase/transport owners.
#[derive(Debug)]
pub struct LumpedEnthalpySpectrumMarch<'a> {
    binding: LumpedEnthalpyPowerBinding<'a>,
    gate: BiotGate,
    config: LumpedEnthalpyMarchConfig,
    controls: LumpedEnthalpySpectrumControls,
    source_start_s: Dd,
    maximum_biot: f64,
    admitted_maximum_step_s: f64,
    attempts: usize,
    spectral_work: usize,
    samples: Vec<LumpedEnthalpySpectrumSample>,
    energy: LumpedEnthalpySpectrumEnergy,
    cumulative_absolute_energy_residual_j: Dd,
}

impl<'a> LumpedEnthalpySpectrumMarch<'a> {
    pub const fn binding(&self) -> LumpedEnthalpyPowerBinding<'a> {
        self.binding
    }

    pub const fn gate(&self) -> BiotGate {
        self.gate
    }

    pub const fn config(&self) -> LumpedEnthalpyMarchConfig {
        self.config
    }

    pub const fn controls(&self) -> LumpedEnthalpySpectrumControls {
        self.controls
    }

    /// Canonical absolute start in the declared source clock, including its
    /// low component. Use this plus the last sample's elapsed DD time when
    /// resuming or refreshing a prescribed source.
    pub const fn source_start_s_dd(&self) -> Dd {
        self.source_start_s
    }

    pub const fn maximum_biot(&self) -> f64 {
        self.maximum_biot
    }

    pub const fn admitted_maximum_step_s(&self) -> f64 {
        self.admitted_maximum_step_s
    }

    pub const fn attempts(&self) -> usize {
        self.attempts
    }

    pub const fn spectral_work(&self) -> usize {
        self.spectral_work
    }

    pub fn samples(&self) -> &[LumpedEnthalpySpectrumSample] {
        &self.samples
    }

    /// Move the accepted history out of a march whose source may be temporary.
    /// No samples are cloned or discarded. Capture any desired energy/work
    /// report before consuming the march; the source binding is not retained.
    pub fn into_samples(self) -> Vec<LumpedEnthalpySpectrumSample> {
        self.samples
    }

    pub const fn energy(&self) -> LumpedEnthalpySpectrumEnergy {
        self.energy
    }

    pub const fn cumulative_absolute_energy_residual_j_dd(&self) -> Dd {
        self.cumulative_absolute_energy_residual_j
    }
}

pub fn solve_lumped_enthalpy_spectrum<'a>(
    binding: LumpedEnthalpyPowerBinding<'a>,
    gate: BiotGate,
    config: LumpedEnthalpyMarchConfig,
    controls: LumpedEnthalpySpectrumControls,
) -> Result<LumpedEnthalpySpectrumMarch<'a>, ConductionError> {
    solve_lumped_enthalpy_spectrum_cancellable(binding, gate, config, controls, || false)
}

/// Solve with the original native phase, transport and bounded implicit solver.
///
/// Ripple admission refines one maximum step, independent of source phase.
/// Its analytic bound is monotone in duration, so all smaller accepted steps
/// inherit it without repeating a full spectral bound traversal. Thermal
/// backward-Euler error remains controlled separately by maximum_step_s;
/// an energy residual alone does not certify nonlinear temporal accuracy.
///
/// Refusal/cancellation drops every partial sample; neither the original body,
/// material sources, spectrum nor any external accepted state is mutated.
pub fn solve_lumped_enthalpy_spectrum_cancellable<'a>(
    binding: LumpedEnthalpyPowerBinding<'a>,
    gate: BiotGate,
    config: LumpedEnthalpyMarchConfig,
    controls: LumpedEnthalpySpectrumControls,
    cancelled: impl FnMut() -> bool,
) -> Result<LumpedEnthalpySpectrumMarch<'a>, ConductionError> {
    solve_lumped_enthalpy_spectrum_from_clock_dd_cancellable(
        binding,
        gate,
        config,
        controls,
        Dd::from_f64(controls.start_s),
        cancelled,
    )
}

/// Native enthalpy evolution from a DD absolute source-clock origin.
///
/// `source_start_s` is authoritative. Its finite components are normalized
/// without an intermediate f64 time conversion. The returned controls retain
/// its rounded f64 view in `start_s`; the supplied controls' `start_s` is
/// replaced. All other controls and original material/thermal admission apply.
pub fn solve_lumped_enthalpy_spectrum_from_clock_dd<'a>(
    binding: LumpedEnthalpyPowerBinding<'a>,
    gate: BiotGate,
    config: LumpedEnthalpyMarchConfig,
    controls: LumpedEnthalpySpectrumControls,
    source_start_s: Dd,
) -> Result<LumpedEnthalpySpectrumMarch<'a>, ConductionError> {
    solve_lumped_enthalpy_spectrum_from_clock_dd_cancellable(
        binding,
        gate,
        config,
        controls,
        source_start_s,
        || false,
    )
}

/// Cancellable DD-clock entrypoint for phase-continuous feedback intervals.
///
/// The absolute origin remains in DD for every analytic source integral.
/// A caller can refresh the original field/material owner between intervals
/// and continue at `march.source_start_s_dd() + last.elapsed_s_dd()` without
/// dropping a short elapsed time at a large epoch. This only preserves the
/// declared source clock: refreshing constitutive coefficients still requires
/// its own physical time-level and feedback-accuracy contract.
///
/// Like the f64 entrypoint, this returns the complete accepted trajectory
/// atomically. Nonfinite or unrepresentable origins refuse before source work.
pub fn solve_lumped_enthalpy_spectrum_from_clock_dd_cancellable<'a>(
    binding: LumpedEnthalpyPowerBinding<'a>,
    gate: BiotGate,
    config: LumpedEnthalpyMarchConfig,
    mut controls: LumpedEnthalpySpectrumControls,
    source_start_s: Dd,
    mut cancelled: impl FnMut() -> bool,
) -> Result<LumpedEnthalpySpectrumMarch<'a>, ConductionError> {
    let body = binding.body;
    poll(&mut cancelled, 0)?;
    if !source_start_s.hi.is_finite() || !source_start_s.lo.is_finite() {
        return Err(lumped_error(
            &body.name,
            "absolute source clock must have finite DD components".to_owned(),
        ));
    }
    let clock_origin = finite(
        body,
        Dd::from_f64(source_start_s.hi) + Dd::from_f64(source_start_s.lo),
        "absolute source-clock origin",
    )?;
    controls.start_s = clock_origin.to_f64();
    validate_enthalpy_config(body, config)?;
    admit_controls(binding, config, controls)?;
    let maximum_biot = body.maximum_biot(config.radiation_temperature_k);
    if !maximum_biot.is_finite() || maximum_biot < 0.0 || maximum_biot > gate.ceiling() {
        return Err(lumped_error(
            &body.name,
            format!(
                "convection-plus-radiation Biot {maximum_biot} is not admitted by ceiling {}; \
                 escalate to spatial transient enthalpy transport",
                gate.ceiling(),
            ),
        ));
    }
    let initial = phase_state(body, config.initial_specific_enthalpy_j_kg)?;
    let (initial_convection, initial_radiation) =
        boundary_power(body, config, initial.temperature_k())?;
    let terms = binding
        .spectrum
        .tones()
        .len()
        .checked_add(1)
        .ok_or_else(|| lumped_error(&body.name, "spectral work size overflow".to_owned()))?;
    let mut attempts = 0usize;
    let mut spectral_work = 0usize;
    let mut admitted_step = config.maximum_step_s.min(config.duration_s);
    let mut ripple_j = 0.0;
    if config.duration_s > 0.0 {
        loop {
            poll(&mut cancelled, attempts)?;
            charge_attempt(body, &mut attempts, controls)?;
            charge_spectrum(body, &mut spectral_work, terms, controls)?;
            let raw = match binding
                .spectrum
                .input_ripple_bound_cancellable(admitted_step, &mut cancelled)
            {
                Ok(bound) => bound,
                Err(ProductSpectrumError::Window(WindowError::Unrepresentable)) => f64::INFINITY,
                Err(error) => return Err(source_error(body, error, attempts)),
            };
            // One outward f64 rounding protects the explicit positive scale.
            // Exact darkness stays zero; overflow requests refinement.
            ripple_j = if raw == 0.0 {
                0.0
            } else {
                (raw * binding.watts_per_spectrum_unit).next_up()
            };
            if ripple_j <= controls.source_ripple_tolerance_j {
                break;
            }
            let smaller = 0.5 * admitted_step;
            if smaller < controls.minimum_step_s || smaller <= 0.0 || smaller == admitted_step {
                return Err(lumped_error(
                    &body.name,
                    "source ripple cannot meet its joule tolerance above the minimum step"
                        .to_owned(),
                ));
            }
            admitted_step = smaller;
        }
    }
    // Reserve only the estimated trajectory, with two spare boundary samples
    // for exact DD horizon partitioning. Runtime caps remain authoritative.
    let estimated = if config.duration_s == 0.0 {
        0.0
    } else {
        (config.duration_s / admitted_step).ceil()
    };
    if !estimated.is_finite() || estimated > MAX_LUMPED_ENTHALPY_STEPS as f64 + 1.0 {
        return Err(lumped_error(
            &body.name,
            "ripple-admitted thermal trajectory exceeds the hard step bound".to_owned(),
        ));
    }
    let step_limit = config.maximum_steps.min(MAX_LUMPED_ENTHALPY_STEPS);
    let capacity = (estimated as usize)
        .saturating_add(3)
        .min(step_limit.saturating_add(1))
        .min(controls.maximum_output_samples);
    admit_output_bytes(body, capacity, controls.maximum_workspace_bytes)?;
    let mut samples = Vec::new();
    samples.try_reserve_exact(capacity).map_err(|_| {
        lumped_error(
            &body.name,
            "accepted trajectory allocation failed".to_owned(),
        )
    })?;
    admit_output_bytes(body, samples.capacity(), controls.maximum_workspace_bytes)?;
    samples.push(LumpedEnthalpySpectrumSample {
        elapsed_s: Dd::ZERO,
        step_duration_s: 0.0,
        phase_state: initial,
        convection_into_body_w: initial_convection,
        radiation_into_body_w: initial_radiation,
        source_ripple_bound_j: 0.0,
        energy: LumpedEnthalpySpectrumEnergy::ZERO,
    });
    let total_duration = Dd::from_f64(config.duration_s);
    let mut elapsed = Dd::ZERO;
    let mut previous_h = config.initial_specific_enthalpy_j_kg;
    let mut totals = LumpedEnthalpySpectrumEnergy::ZERO;
    let mut cumulative_absolute_residual = Dd::ZERO;
    while positive(total_duration - elapsed) {
        poll(&mut cancelled, attempts)?;
        let accepted_steps = samples.len() - 1;
        if accepted_steps >= step_limit || samples.len() >= capacity {
            return Err(lumped_error(
                &body.name,
                "accepted-step or output-sample budget exhausted".to_owned(),
            ));
        }
        charge_attempt(body, &mut attempts, controls)?;
        let remaining = finite(body, total_duration - elapsed, "remaining duration")?;
        let mut dt_s = admitted_step.min(remaining.to_f64());
        if positive(Dd::from_f64(dt_s) - remaining) {
            dt_s = dt_s.next_down();
        }
        if !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(lumped_error(
                &body.name,
                "the remaining source-clock interval is not representable".to_owned(),
            ));
        }
        let next_elapsed = finite(body, elapsed + Dd::from_f64(dt_s), "accepted time")?;
        if !positive(next_elapsed - elapsed) || positive(next_elapsed - total_duration) {
            return Err(lumped_error(
                &body.name,
                "accepted source-clock interval failed to make bounded progress".to_owned(),
            ));
        }
        let clock_start = finite(body, clock_origin + elapsed, "absolute source clock")?;
        charge_spectrum(body, &mut spectral_work, terms, controls)?;
        let source = binding
            .spectrum
            .causal_integral_from_clock_dd_scaled_cancellable(
                clock_start,
                dt_s,
                0.0,
                binding.watts_per_spectrum_unit,
                &mut cancelled,
            )
            .map_err(|error| source_error(body, error, attempts))?;
        let source_j = finite(body, source.integral_dd(), "integrated source energy")?;
        let mut residual = |specific_enthalpy: f64| -> Result<f64, ConductionError> {
            poll(&mut cancelled, attempts)?;
            let state = phase_state(body, specific_enthalpy)?;
            let energy = step_energy(body, config, previous_h, state, dt_s, source_j)?;
            Ok(energy.residual_j.to_f64())
        };
        // Preserve exact dark/equilibrium states instead of injecting a root
        // bisection interval's rounding into an already exact zero balance.
        let next_h = if residual(previous_h)? == 0.0 {
            previous_h
        } else {
            bisect_enthalpy_step(body, config.enthalpy_tolerance_j_kg, &mut residual)?
        };
        let state = phase_state(body, next_h)?;
        let (convection, radiation) = boundary_power(body, config, state.temperature_k())?;
        let energy = step_energy(body, config, previous_h, state, dt_s, source_j)?;
        let absolute_residual = abs_dd(energy.residual_j);
        if positive(absolute_residual - Dd::from_f64(controls.maximum_step_energy_residual_j)) {
            return Err(lumped_error(
                &body.name,
                format!(
                    "accepted candidate energy residual {} J exceeds {}; tighten native enthalpy \
                     tolerance or revise the numerical budget",
                    absolute_residual.to_f64(),
                    controls.maximum_step_energy_residual_j,
                ),
            ));
        }
        let next_totals = sum_energy(body, totals, energy)?;
        let next_absolute_residual = finite(
            body,
            cumulative_absolute_residual + absolute_residual,
            "cumulative absolute energy residual",
        )?;
        poll(&mut cancelled, attempts)?;
        samples.push(LumpedEnthalpySpectrumSample {
            elapsed_s: next_elapsed,
            step_duration_s: dt_s,
            phase_state: state,
            convection_into_body_w: convection,
            radiation_into_body_w: radiation,
            source_ripple_bound_j: ripple_j,
            energy,
        });
        elapsed = next_elapsed;
        previous_h = next_h;
        totals = next_totals;
        cumulative_absolute_residual = next_absolute_residual;
    }
    // Total storage is also evaluated directly from the actual original and
    // final phase coordinates, rather than assuming every step closed exactly.
    totals.stored_energy_change_j = finite(
        body,
        Dd::from_f64(body.mass_kg)
            * (Dd::from_f64(previous_h) - Dd::from_f64(config.initial_specific_enthalpy_j_kg)),
        "total stored energy change",
    )?;
    totals.residual_j = balance(body, totals)?;
    poll(&mut cancelled, attempts)?;
    Ok(LumpedEnthalpySpectrumMarch {
        binding,
        gate,
        config,
        controls,
        source_start_s: clock_origin,
        maximum_biot,
        admitted_maximum_step_s: admitted_step,
        attempts,
        spectral_work,
        samples,
        energy: totals,
        cumulative_absolute_energy_residual_j: cumulative_absolute_residual,
    })
}

fn admit_controls(
    binding: LumpedEnthalpyPowerBinding<'_>,
    config: LumpedEnthalpyMarchConfig,
    controls: LumpedEnthalpySpectrumControls,
) -> Result<(), ConductionError> {
    let body = binding.body;
    if binding.input_unit.trim().is_empty() || binding.source_clock.trim().is_empty() {
        return Err(lumped_error(
            &body.name,
            "spectrum units and common source-clock label must be explicit".to_owned(),
        ));
    }
    if !binding.watts_per_spectrum_unit.is_finite() || binding.watts_per_spectrum_unit <= 0.0 {
        return Err(lumped_error(
            &body.name,
            "watts-per-spectrum-unit scale must be finite and positive".to_owned(),
        ));
    }
    if config.internal_power_w != 0.0 {
        return Err(lumped_error(
            &body.name,
            "put all constant/internal power into the spectrum; native internal_power_w must be zero"
                .to_owned(),
        ));
    }
    if !controls.start_s.is_finite()
        || !controls.minimum_step_s.is_finite()
        || controls.minimum_step_s <= 0.0
        || controls.minimum_step_s > config.maximum_step_s
        || !controls.source_ripple_tolerance_j.is_finite()
        || controls.source_ripple_tolerance_j < 0.0
        || !controls.maximum_step_energy_residual_j.is_finite()
        || controls.maximum_step_energy_residual_j < 0.0
    {
        return Err(lumped_error(
            &body.name,
            "source clock, minimum step, ripple tolerance or energy residual budget is invalid"
                .to_owned(),
        ));
    }
    if controls.maximum_attempts == 0
        || controls.maximum_spectral_work == 0
        || controls.maximum_output_samples == 0
    {
        return Err(lumped_error(
            &body.name,
            "attempt, spectral-work and output-sample bounds must be positive".to_owned(),
        ));
    }
    Ok(())
}

fn step_energy(
    body: &LumpedEnthalpyBody<'_>,
    config: LumpedEnthalpyMarchConfig,
    previous_h: f64,
    state: EquilibriumPhaseState,
    dt_s: f64,
    source_j: Dd,
) -> Result<LumpedEnthalpySpectrumEnergy, ConductionError> {
    let (convection, radiation) = boundary_power(body, config, state.temperature_k())?;
    let duration = Dd::from_f64(dt_s);
    let mut energy = LumpedEnthalpySpectrumEnergy {
        source_input_j: source_j,
        convection_into_body_j: finite(
            body,
            duration * Dd::from_f64(convection),
            "convective energy",
        )?,
        radiation_into_body_j: finite(
            body,
            duration * Dd::from_f64(radiation),
            "radiative energy",
        )?,
        stored_energy_change_j: finite(
            body,
            Dd::from_f64(body.mass_kg)
                * (Dd::from_f64(state.specific_enthalpy_j_kg()) - Dd::from_f64(previous_h)),
            "stored energy change",
        )?,
        residual_j: Dd::ZERO,
    };
    energy.residual_j = balance(body, energy)?;
    Ok(energy)
}

fn boundary_power(
    body: &LumpedEnthalpyBody<'_>,
    config: LumpedEnthalpyMarchConfig,
    temperature_k: f64,
) -> Result<(f64, f64), ConductionError> {
    let (convection, radiation, net) = enthalpy_power(body, config, temperature_k)?;
    if [convection, radiation, net]
        .iter()
        .all(|value| value.is_finite())
    {
        Ok((convection, radiation))
    } else {
        Err(lumped_error(
            &body.name,
            "original convection/radiation power is not finite".to_owned(),
        ))
    }
}

fn balance(
    body: &LumpedEnthalpyBody<'_>,
    energy: LumpedEnthalpySpectrumEnergy,
) -> Result<Dd, ConductionError> {
    finite(
        body,
        energy.stored_energy_change_j
            - energy.source_input_j
            - energy.convection_into_body_j
            - energy.radiation_into_body_j,
        "energy balance residual",
    )
}

fn sum_energy(
    body: &LumpedEnthalpyBody<'_>,
    a: LumpedEnthalpySpectrumEnergy,
    b: LumpedEnthalpySpectrumEnergy,
) -> Result<LumpedEnthalpySpectrumEnergy, ConductionError> {
    Ok(LumpedEnthalpySpectrumEnergy {
        source_input_j: finite(
            body,
            a.source_input_j + b.source_input_j,
            "total source input",
        )?,
        convection_into_body_j: finite(
            body,
            a.convection_into_body_j + b.convection_into_body_j,
            "total convective energy",
        )?,
        radiation_into_body_j: finite(
            body,
            a.radiation_into_body_j + b.radiation_into_body_j,
            "total radiative energy",
        )?,
        stored_energy_change_j: finite(
            body,
            a.stored_energy_change_j + b.stored_energy_change_j,
            "summed storage change",
        )?,
        residual_j: finite(body, a.residual_j + b.residual_j, "summed residual")?,
    })
}

fn charge_attempt(
    body: &LumpedEnthalpyBody<'_>,
    attempts: &mut usize,
    controls: LumpedEnthalpySpectrumControls,
) -> Result<(), ConductionError> {
    if *attempts
        >= controls
            .maximum_attempts
            .min(MAX_ENTHALPY_SPECTRUM_ATTEMPTS)
    {
        return Err(lumped_error(
            &body.name,
            "attempt budget exhausted".to_owned(),
        ));
    }
    *attempts += 1;
    Ok(())
}

fn charge_spectrum(
    body: &LumpedEnthalpyBody<'_>,
    work: &mut usize,
    terms: usize,
    controls: LumpedEnthalpySpectrumControls,
) -> Result<(), ConductionError> {
    *work = work
        .checked_add(terms)
        .ok_or_else(|| lumped_error(&body.name, "spectral work count overflow".to_owned()))?;
    if *work
        > controls
            .maximum_spectral_work
            .min(MAX_ENTHALPY_SPECTRUM_WORK)
    {
        return Err(lumped_error(
            &body.name,
            "spectral work budget exhausted".to_owned(),
        ));
    }
    Ok(())
}

fn admit_output_bytes(
    body: &LumpedEnthalpyBody<'_>,
    capacity: usize,
    maximum: usize,
) -> Result<(), ConductionError> {
    if capacity
        .checked_mul(core::mem::size_of::<LumpedEnthalpySpectrumSample>())
        .is_none_or(|bytes| bytes > maximum)
    {
        Err(lumped_error(
            &body.name,
            "sample storage budget exceeded".to_owned(),
        ))
    } else {
        Ok(())
    }
}

fn finite(body: &LumpedEnthalpyBody<'_>, value: Dd, quantity: &str) -> Result<Dd, ConductionError> {
    if value.hi.is_finite() && value.lo.is_finite() && value.to_f64().is_finite() {
        Ok(value)
    } else {
        Err(lumped_error(
            &body.name,
            format!("{quantity} is not representable"),
        ))
    }
}

fn positive(value: Dd) -> bool {
    value.hi > 0.0 || (value.hi == 0.0 && value.lo > 0.0)
}

fn abs_dd(value: Dd) -> Dd {
    if positive(-value) { -value } else { value }
}

fn poll(cancelled: &mut impl FnMut() -> bool, at: usize) -> Result<(), ConductionError> {
    if cancelled() {
        Err(ConductionError::Cancelled {
            stage: "lumped enthalpy interference",
            at,
        })
    } else {
        Ok(())
    }
}

fn source_error(
    body: &LumpedEnthalpyBody<'_>,
    error: ProductSpectrumError,
    at: usize,
) -> ConductionError {
    match error {
        ProductSpectrumError::Cancelled => ConductionError::Cancelled {
            stage: "lumped enthalpy interference spectrum",
            at,
        },
        error => lumped_error(
            &body.name,
            format!("thermal power spectrum refused: {error}"),
        ),
    }
}

#[cfg(test)]
mod tests;
