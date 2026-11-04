//! Finite thermal contact driving the original two-fluid channel.
//!
//! A mass-blocked endpoint exchanges reference-state heat Q=T0*Sdot with one
//! original thermal cell. Its temperature and the channel's midpoint endpoint
//! temperature solve Q=G*(T_channel-T_cell) together. The finite cell's actual
//! caloric chart supplies its response; no replacement heat capacity is fitted.
//! Native sources, contacts and events march first, then this midpoint contact
//! exchanges heat. Coupling to that network is a first-order operator split.
//!
//! This is a SMALL-SIGNAL thermal connection about the channel's frozen T0.
//! Reference heat/entropy transport and quadratic wave availability have separate
//! accounts. In particular, boundary wave work is -T'*Q*dt/T0, not -Q*dt.
//! Native caloric joules balance the transported reference heat. This does not
//! supply a nonlinear helium internal-energy law or a total-energy theorem for
//! finite temperature excursions. Conductance and its valid temperature range
//! must be supplied for the actual interface; Kapitza conductance is not guessed.

use std::collections::BTreeMap;

use fs_duct::superfluid::{
    BodyForce, ChannelEnd, ChannelError, ChannelStep, HeatStepError, MutualFriction,
    TwoFluidChannel,
};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionSimulation, EnergyAudit, StepConfig,
};

/// Heat-only contact; the opposite channel end is closed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CounterflowThermalContact {
    pub end: ChannelEnd,
    pub cell: CellId,
    /// Total interface conductance [W/K], including its physical area. Zero
    /// opens the thermal contact without changing either material's properties.
    pub conductance_w_k: f64,
}

/// Physical small-signal range, bounded scalar solve and native heat admission.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CounterflowHeatControls {
    pub thermal: StepConfig,
    /// Maximum scalar caloric probes; 1..=128, without repeated wave solves.
    pub max_iterations: usize,
    pub max_abs_heat_power_w: f64,
    /// Allowed thermal-cell departure from T0, at most 10% of T0. This is a
    /// caller admission limit, not a material phase or linearization certificate.
    pub max_temperature_departure_k: f64,
    pub contact_power_tolerance_w: f64,
    pub cell_energy_tolerance_j: f64,
    pub transfer_energy_tolerance_j: f64,
}

/// Original owner reports plus the heat, entropy and availability port accounts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CounterflowHeatStep {
    pub channel: ChannelStep,
    pub thermal_march: EnergyAudit,
    pub heat_intake: EnergyAudit,
    /// Positive from the channel into the finite thermal cell.
    pub heat_outward_w: f64,
    pub heat_to_thermal_j: f64,
    pub channel_midpoint_temperature_k: f64,
    pub thermal_midpoint_temperature_k: f64,
    /// Q-G*(T_channel_mid-T_cell_mid), W.
    pub contact_power_residual_w: f64,
    /// Native caloric heat intake plus channel reference heat import, J.
    pub reference_heat_transfer_residual_j: f64,
    /// Channel boundary entropy import plus Q*dt/T0, J/K.
    pub reference_entropy_transfer_residual_j_k: f64,
    /// Original boundary wave work plus T'_mid*Q*dt/T0, J.
    pub wave_port_residual_j: f64,
    /// Linearized interface availability loss G*(delta T_mid)^2*dt/T0.
    /// This is NOT extra deposited heat and is never added to heat_intake.
    pub contact_availability_loss_j: f64,
    pub iterations: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum CounterflowHeatError {
    Invalid(&'static str),
    Channel(ChannelError),
    Thermal(DepositionError),
    PowerLimit,
    Unresolved { residual_w: f64, iterations: usize },
    Cancelled,
}
impl core::fmt::Display for CounterflowHeatError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "superfluid thermal contact refused: {self:?}")
    }
}
impl std::error::Error for CounterflowHeatError {}
impl From<ChannelError> for CounterflowHeatError {
    fn from(value: ChannelError) -> Self {
        match value {
            ChannelError::Cancelled => Self::Cancelled,
            other => Self::Channel(other),
        }
    }
}
impl From<DepositionError> for CounterflowHeatError {
    fn from(value: DepositionError) -> Self {
        match value {
            DepositionError::Cancelled => Self::Cancelled,
            other => Self::Thermal(other),
        }
    }
}
impl From<HeatStepError<CounterflowHeatError>> for CounterflowHeatError {
    fn from(value: HeatStepError<CounterflowHeatError>) -> Self {
        match value {
            HeatStepError::Channel(error) => error.into(),
            HeatStepError::Heat(error) => error,
        }
    }
}

/// Advance a heat-driven channel and its actual finite thermal source or sink.
///
/// The two clocks must match. Native thermal sources, links and scheduled events
/// advance once, and the reciprocal contact uses the mean of the original and
/// final cell temperatures. The original immutable caloric chart is probed at
/// candidate energies, including any supported latent plateau. A scalar bracket
/// follows from its required monotonicity and the positive wave impedance.
/// All root iterations are O(1) caloric queries, not channel or network solves.
///
/// The other channel end is closed and both endpoints have zero total mass
/// current. Optional mutual-friction losses remain exported by the original
/// channel report and per-cell heat field; they are NOT silently deposited into
/// this boundary cell. Geometry and finite thermal storage do not relinearize
/// helium properties. Both histories publish together only after all original
/// source, caloric, amplitude, transfer and cancellation checks succeed.
#[allow(clippy::too_many_arguments)]
pub fn advance_counterflow_heat(
    channel: &mut TwoFluidChannel,
    thermal: &mut DepositionSimulation<'_>,
    contact: CounterflowThermalContact,
    dt_s: f64,
    friction: MutualFriction,
    controls: CounterflowHeatControls,
    cancelled: impl Fn() -> bool,
) -> Result<CounterflowHeatStep, CounterflowHeatError> {
    check(&cancelled)?;
    let t0 = channel.material().temperature_k();
    if channel.time_s() != thermal.time_s()
        || !nonnegative(contact.conductance_w_k)
        || !(1..=128).contains(&controls.max_iterations)
        || !positive(controls.max_abs_heat_power_w)
        || !positive(controls.max_temperature_departure_k)
        || controls.max_temperature_departure_k > 0.1 * t0
        || ![
            controls.contact_power_tolerance_w,
            controls.cell_energy_tolerance_j,
            controls.transfer_energy_tolerance_j,
        ]
        .into_iter()
        .all(nonnegative)
    {
        return Err(CounterflowHeatError::Invalid(
            "matching clocks and finite admitted thermal-contact controls required",
        ));
    }
    let temperature_before = thermal.temperature_k(contact.cell)?;
    admit_temperature(temperature_before, t0, controls)?;
    let slot = match contact.end {
        ChannelEnd::Left => 0,
        ChannelEnd::Right => 1,
    };
    let response =
        channel
            .counterflow_ends_response(dt_s, [0.0; 2], BodyForce::None, friction, || cancelled())?;
    let (marched, thermal_march) =
        thermal.advance_to(channel.time_s() + dt_s, controls.thermal, || {
            thermal_check(&cancelled)
        })?;
    let cell = marched
        .cells()
        .get(&contact.cell)
        .ok_or(CounterflowHeatError::Invalid("inactive contact cell"))?;
    let initial = thermal
        .cells()
        .get(&contact.cell)
        .ok_or(CounterflowHeatError::Invalid(
            "contact cell must be active at both endpoints",
        ))?;
    if cell.mass_kg != initial.mass_kg || cell.material_id != initial.material_id {
        return Err(CounterflowHeatError::Invalid(
            "thermal contact requires unchanged cell mass and caloric law",
        ));
    }
    let impedance = response.temperature_heat_impedance_k_w[slot][slot];
    let theta_reference = response.reference_temperature_perturbation_k[slot];
    let evaluate = |q: f64| -> Result<ContactSolution, CounterflowHeatError> {
        check(&cancelled)?;
        let heat_j = q * dt_s;
        let proposed_h = cell.specific_h_j_kg + heat_j / cell.mass_kg;
        let temperature_after =
            marched.temperature_at_specific_enthalpy_k(contact.cell, proposed_h)?;
        let thermal_theta = 0.5 * (temperature_before - t0) + 0.5 * (temperature_after - t0);
        let channel_theta = theta_reference - impedance * q;
        let difference = channel_theta - thermal_theta;
        let residual = q - contact.conductance_w_k * difference;
        if ![heat_j, thermal_theta, channel_theta, residual]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(CounterflowHeatError::Invalid(
                "thermal contact arithmetic is not representable",
            ));
        }
        Ok(ContactSolution {
            q,
            channel_theta,
            thermal_theta,
            temperature_after,
            residual,
            iterations: 0,
        })
    };
    let solution = solve_contact(evaluate, impedance, contact.conductance_w_k, controls)?;
    admit_temperature(solution.temperature_after, t0, controls)?;
    let heat_j = solution.q * dt_s;
    let (next_thermal, heat_intake) = marched.apply_heat_j(
        &BTreeMap::from([(contact.cell, heat_j)]),
        controls.cell_energy_tolerance_j,
        controls.transfer_energy_tolerance_j,
        || thermal_check(&cancelled),
    )?;
    let mut q_out = [0.0; 2];
    q_out[slot] = solution.q;
    let boundary = channel.counterflow_boundary_velocities(q_out)?;
    let (_, (next, report)) = channel.step_with_heat(
        dt_s,
        boundary,
        friction,
        || cancelled(),
        |wave, _| {
            let heat_residual = heat_intake.external_heat_j + wave.boundary_counterflow_heat_j;
            let entropy_residual = wave.boundary_entropy_j_k + heat_j / t0;
            let wave_residual = wave.boundary_wave_work_j + solution.channel_theta / t0 * heat_j;
            let delta = solution.channel_theta - solution.thermal_theta;
            let availability_loss = dt_s * contact.conductance_w_k * delta * delta / t0;
            if ![
                heat_residual,
                entropy_residual,
                wave_residual,
                availability_loss,
            ]
            .into_iter()
            .all(f64::is_finite)
                || heat_residual.abs() > controls.transfer_energy_tolerance_j
                || (t0 * entropy_residual).abs() > controls.transfer_energy_tolerance_j
                || wave_residual.abs() > controls.transfer_energy_tolerance_j
            {
                return Err(CounterflowHeatError::Invalid(
                    "reference heat, entropy or wave port transfer did not close",
                ));
            }
            check(&cancelled)?;
            Ok((
                next_thermal,
                CounterflowHeatStep {
                    channel: *wave,
                    thermal_march,
                    heat_intake,
                    heat_outward_w: solution.q,
                    heat_to_thermal_j: heat_j,
                    channel_midpoint_temperature_k: t0 + solution.channel_theta,
                    thermal_midpoint_temperature_k: t0 + solution.thermal_theta,
                    contact_power_residual_w: solution.residual,
                    reference_heat_transfer_residual_j: heat_residual,
                    reference_entropy_transfer_residual_j_k: entropy_residual,
                    wave_port_residual_j: wave_residual,
                    contact_availability_loss_j: availability_loss,
                    iterations: solution.iterations,
                },
            ))
        },
    )?;
    *thermal = next;
    Ok(report)
}

#[derive(Clone, Copy)]
struct ContactSolution {
    q: f64,
    channel_theta: f64,
    thermal_theta: f64,
    temperature_after: f64,
    residual: f64,
    iterations: usize,
}

fn solve_contact(
    evaluate: impl Fn(f64) -> Result<ContactSolution, CounterflowHeatError>,
    impedance: f64,
    conductance: f64,
    controls: CounterflowHeatControls,
) -> Result<ContactSolution, CounterflowHeatError> {
    let zero = evaluate(0.0)?;
    if zero.residual.abs() <= controls.contact_power_tolerance_w {
        return Ok(ContactSolution {
            iterations: 1,
            ..zero
        });
    }
    let difference = zero.channel_theta - zero.thermal_theta;
    let sign = difference.signum();
    // Holding the finite cell at its zero-exchange midpoint gives the maximum
    // possible heat magnitude. A monotone caloric law can only reduce it.
    let resistance = 1.0 / conductance + impedance;
    let bound = difference.abs() / resistance;
    if !positive(resistance) || !positive(bound) {
        return Err(CounterflowHeatError::Invalid(
            "thermal contact bracket is unresolved",
        ));
    }
    let mut lower = 0.0;
    let mut upper = bound.min(controls.max_abs_heat_power_w);
    let mut x = upper;
    let mut last_residual = zero.residual;
    let mut domain_refused = false;
    for iteration in 1..controls.max_iterations {
        match evaluate(sign * x) {
            Ok(mut trial) => {
                trial.iterations = iteration + 1;
                last_residual = trial.residual;
                if trial.residual.abs() <= controls.contact_power_tolerance_w {
                    return Ok(trial);
                }
                if sign * trial.residual < 0.0 {
                    if iteration == 1 && upper < bound {
                        return Err(CounterflowHeatError::PowerLimit);
                    }
                    lower = x;
                } else {
                    upper = x;
                }
            }
            Err(CounterflowHeatError::Thermal(DepositionError::MaterialDomain)) => {
                // The root can remain inside the original chart even when the
                // ideal-reservoir upper bound lies outside it. Never extrapolate.
                domain_refused = true;
                upper = x;
            }
            Err(error) => return Err(error),
        }
        let next = 0.5 * lower + 0.5 * upper;
        if next == lower || next == upper {
            break;
        }
        x = next;
    }
    if domain_refused {
        Err(DepositionError::MaterialDomain.into())
    } else {
        Err(CounterflowHeatError::Unresolved {
            residual_w: last_residual,
            iterations: controls.max_iterations,
        })
    }
}

fn admit_temperature(
    temperature: f64,
    reference: f64,
    controls: CounterflowHeatControls,
) -> Result<(), CounterflowHeatError> {
    if !positive(temperature)
        || (temperature - reference).abs() > controls.max_temperature_departure_k
    {
        Err(CounterflowHeatError::Invalid(
            "thermal cell left the declared small-signal range",
        ))
    } else {
        Ok(())
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn check(cancelled: &impl Fn() -> bool) -> Result<(), CounterflowHeatError> {
    if cancelled() {
        Err(CounterflowHeatError::Cancelled)
    } else {
        Ok(())
    }
}
fn thermal_check(cancelled: &impl Fn() -> bool) -> Result<(), DepositionError> {
    if cancelled() {
        Err(DepositionError::Cancelled)
    } else {
        Ok(())
    }
}
