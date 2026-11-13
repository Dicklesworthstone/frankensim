//! Two finite thermal contacts on the original reciprocal two-fluid channel.
//!
//! One native thermal cell is connected to each mass-blocked channel end. The
//! two heat powers are solved together with the complete 2x2 endpoint response,
//! so waves communicating between the ends affect both contacts. Each cell's
//! original monotone caloric chart supplies its finite storage, including latent
//! plateaux; no replacement heat capacity or duplicate wave solver is used.
//!
//! Native thermal sources, links and events advance once before the contact
//! exchange, a first-order operator split. Midpoint contact temperatures use
//! the accepted initial and final cell temperatures. The two clocks, both
//! caloric updates and the physical wave step publish atomically.
//!
//! As with the single contact, this transports reference heat Q=T0*Sdot in the
//! frozen small-signal channel. Quadratic wave availability, reference heat and
//! entropy have separate accounts. It is not a finite-excursion helium EOS or
//! total internal-energy theorem. Wall/mutual-friction heat stays in the original
//! channel's spatial heat field; it is not assigned to either endpoint cell.
//! The optional [`dissipation`] composition routes that field to explicitly
//! mapped distinct native caloric receivers within the same atomic update.

use std::collections::BTreeMap;

use fs_duct::superfluid::{BodyForce, ChannelEnd, ChannelStep, MutualFriction, TwoFluidChannel};
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation, EnergyAudit};

use super::counterflow::{
    CounterflowContactLaw, CounterflowHeatControls, CounterflowHeatError,
    CounterflowThermalContact, admit_temperature, check, nonnegative, positive, thermal_check,
};

/// Explicit spatial dissipation destinations on the same finite thermal timeline.
pub mod dissipation;

/// Original owner reports and both finite-contact accounts, ordered [left,right].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CounterflowPairHeatStep {
    pub channel: ChannelStep,
    pub thermal_march: EnergyAudit,
    pub heat_intake: EnergyAudit,
    /// Positive from each channel end into its finite thermal cell [W].
    pub heat_outward_w: [f64; 2],
    pub heat_to_thermal_j: [f64; 2],
    pub channel_midpoint_temperature_k: [f64; 2],
    pub thermal_midpoint_temperature_k: [f64; 2],
    /// Q_i minus the original interface law at both midpoint temperatures, W.
    pub contact_power_residual_w: [f64; 2],
    pub reference_heat_transfer_residual_j: f64,
    pub reference_entropy_transfer_residual_j_k: f64,
    pub wave_port_residual_j: f64,
    /// Each interface's `law(T_channel_mid_i,T_cell_mid_i)*delta T_mid_i*dt/T0`. These availability
    /// losses are not extra caloric heat and are never deposited twice.
    pub contact_availability_loss_j: [f64; 2],
    /// Caloric probes per endpoint, each bounded by controls.max_iterations.
    pub iterations: [usize; 2],
}

/// Advance a finite heater-to-sink channel with two original caloric owners.
///
/// Supply two contacts in either order, one per end, attached to DISTINCT cells
/// of the same thermal simulation. Both cells must remain active with unchanged
/// mass and material law over the interval. Each end independently selects a
/// linear or Kapitza interface. Zero linear conductance closes that heat port;
/// one open contact reduces to the single-contact equations.
///
/// The original wave owner supplies one reciprocal response. The bounded
/// nonlinear iteration makes only scalar caloric probes, followed by one native
/// combined heat intake and one physical wave stage. Its projected coordinate
/// solves use the positive diagonal response and monotone caloric laws. A
/// temporary coordinate at a power/chart edge is allowed to move inward as the
/// opposite contact converges; it does not trigger premature refusal. Only a
/// solution satisfying BOTH physical contact equations is published.
#[allow(clippy::too_many_arguments)]
pub fn advance_counterflow_heat_pair(
    channel: &mut TwoFluidChannel,
    thermal: &mut DepositionSimulation<'_>,
    contacts: [CounterflowThermalContact; 2],
    dt_s: f64,
    friction: MutualFriction,
    controls: CounterflowHeatControls,
    cancelled: impl Fn() -> bool,
) -> Result<CounterflowPairHeatStep, CounterflowHeatError> {
    advance_pair_prepared(
        channel,
        thermal,
        contacts,
        dt_s,
        friction,
        controls,
        cancelled,
        |next, report, _| Ok((next, report)),
    )
}

#[allow(clippy::too_many_arguments)]
fn advance_pair_prepared<'a, T>(
    channel: &mut TwoFluidChannel,
    thermal: &mut DepositionSimulation<'a>,
    mut contacts: [CounterflowThermalContact; 2],
    dt_s: f64,
    friction: MutualFriction,
    controls: CounterflowHeatControls,
    cancelled: impl Fn() -> bool,
    prepare: impl FnOnce(
        DepositionSimulation<'a>,
        CounterflowPairHeatStep,
        &[f64],
    ) -> Result<(DepositionSimulation<'a>, T), CounterflowHeatError>,
) -> Result<T, CounterflowHeatError> {
    check(&cancelled)?;
    let t0 = channel.material().temperature_k();
    if contacts[0].end == ChannelEnd::Right {
        contacts.swap(0, 1);
    }
    for contact in contacts {
        contact.law.validate()?;
    }
    if contacts[0].end != ChannelEnd::Left
        || contacts[1].end != ChannelEnd::Right
        || contacts[0].cell == contacts[1].cell
        || channel.time_s() != thermal.time_s()
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
            "opposite ends, distinct thermal cells, matching clocks and admitted contact controls required",
        ));
    }
    let temperature_before = [
        thermal.temperature_k(contacts[0].cell)?,
        thermal.temperature_k(contacts[1].cell)?,
    ];
    for temperature in temperature_before {
        admit_temperature(temperature, t0, controls)?;
    }
    let channel_theta_before = [
        channel.cell(0).unwrap().temperature_perturbation_k,
        channel
            .cell(channel.cells() - 1)
            .unwrap()
            .temperature_perturbation_k,
    ];
    for i in 0..2 {
        contacts[i]
            .law
            .heat_w(t0, channel_theta_before[i], temperature_before[i] - t0)?;
    }
    let response =
        channel
            .counterflow_ends_response(dt_s, [0.0; 2], BodyForce::None, friction, || cancelled())?;
    let (marched, thermal_march) =
        thermal.advance_to(channel.time_s() + dt_s, controls.thermal, || {
            thermal_check(&cancelled)
        })?;
    let mut mass = [0.0; 2];
    let mut marched_h = [0.0; 2];
    for i in 0..2 {
        let initial =
            thermal
                .cells()
                .get(&contacts[i].cell)
                .ok_or(CounterflowHeatError::Invalid(
                    "contact cells must be active at both endpoints",
                ))?;
        let cell = marched
            .cells()
            .get(&contacts[i].cell)
            .ok_or(CounterflowHeatError::Invalid("inactive contact cell"))?;
        if cell.mass_kg != initial.mass_kg || cell.material_id != initial.material_id {
            return Err(CounterflowHeatError::Invalid(
                "thermal contacts require unchanged cell masses and caloric laws",
            ));
        }
        mass[i] = cell.mass_kg;
        marched_h[i] = cell.specific_h_j_kg;
    }
    let evaluate = |i: usize, q: f64| -> Result<CaloricPoint, CounterflowHeatError> {
        check(&cancelled)?;
        let heat = q * dt_s;
        let h = marched_h[i] + heat / mass[i];
        let temperature_after = marched.temperature_at_specific_enthalpy_k(contacts[i].cell, h)?;
        let thermal_theta = 0.5 * (temperature_before[i] - t0) + 0.5 * (temperature_after - t0);
        if ![q, heat, thermal_theta, temperature_after]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(CounterflowHeatError::Invalid(
                "finite-contact caloric arithmetic is unrepresentable",
            ));
        }
        Ok(CaloricPoint {
            q,
            thermal_theta,
            temperature_after,
        })
    };
    let laws = contacts.map(|c| c.law);
    let mut solution = solve_pair(
        evaluate,
        response.reference_temperature_perturbation_k,
        response.temperature_heat_impedance_k_w,
        laws,
        t0,
        controls,
    )?;
    let mut interface_power = [0.0; 2];
    for (i, point) in solution.point.into_iter().enumerate() {
        admit_temperature(point.temperature_after, t0, controls)?;
        interface_power[i] = laws[i].heat_w(t0, solution.channel_theta[i], point.thermal_theta)?;
        solution.residual[i] = point.q - interface_power[i];
        if solution.residual[i].abs() > controls.contact_power_tolerance_w {
            return Err(CounterflowHeatError::Unresolved {
                residual_w: solution.residual[i],
                iterations: solution.probes[i],
            });
        }
        // The channel temperature is unchanged by either friction stage;
        // recover and admit the original wave endpoint from its midpoint.
        laws[i].heat_w(
            t0,
            2.0 * solution.channel_theta[i] - channel_theta_before[i],
            point.temperature_after - t0,
        )?;
    }
    let q = solution.point.map(|point| point.q);
    let heat_j = q.map(|power| power * dt_s);
    let (next_thermal, heat_intake) = marched.apply_heat_j(
        &BTreeMap::from([(contacts[0].cell, heat_j[0]), (contacts[1].cell, heat_j[1])]),
        controls.cell_energy_tolerance_j,
        controls.transfer_energy_tolerance_j,
        || thermal_check(&cancelled),
    )?;
    let boundary = channel.counterflow_boundary_velocities(q)?;
    let (_, (next, report)) = channel.step_with_heat(
        dt_s,
        boundary,
        friction,
        || cancelled(),
        |wave, per_cell_heat| {
            let total_heat = heat_j[0] + heat_j[1];
            let heat_residual = heat_intake.external_heat_j + wave.boundary_counterflow_heat_j;
            let entropy_residual = wave.boundary_entropy_j_k + total_heat / t0;
            let wave_residual = wave.boundary_wave_work_j
                + solution.channel_theta[0] / t0 * heat_j[0]
                + solution.channel_theta[1] / t0 * heat_j[1];
            let availability_loss = std::array::from_fn(|i| {
                let delta = solution.channel_theta[i] - solution.point[i].thermal_theta;
                dt_s * interface_power[i] * delta / t0
            });
            if ![
                heat_residual,
                entropy_residual,
                wave_residual,
                availability_loss[0],
                availability_loss[1],
            ]
            .into_iter()
            .all(f64::is_finite)
                || heat_residual.abs() > controls.transfer_energy_tolerance_j
                || (t0 * entropy_residual).abs() > controls.transfer_energy_tolerance_j
                || wave_residual.abs() > controls.transfer_energy_tolerance_j
            {
                return Err(CounterflowHeatError::Invalid(
                    "paired reference heat, entropy or wave port transfer did not close",
                ));
            }
            check(&cancelled)?;
            prepare(
                next_thermal,
                CounterflowPairHeatStep {
                    channel: *wave,
                    thermal_march,
                    heat_intake,
                    heat_outward_w: q,
                    heat_to_thermal_j: heat_j,
                    channel_midpoint_temperature_k: solution.channel_theta.map(|theta| t0 + theta),
                    thermal_midpoint_temperature_k: solution.point.map(|p| t0 + p.thermal_theta),
                    contact_power_residual_w: solution.residual,
                    reference_heat_transfer_residual_j: heat_residual,
                    reference_entropy_transfer_residual_j_k: entropy_residual,
                    wave_port_residual_j: wave_residual,
                    contact_availability_loss_j: availability_loss,
                    iterations: solution.probes,
                },
                per_cell_heat,
            )
        },
    )?;
    *thermal = next;
    Ok(report)
}

#[derive(Clone, Copy)]
struct CaloricPoint {
    q: f64,
    thermal_theta: f64,
    temperature_after: f64,
}

struct PairSolution {
    point: [CaloricPoint; 2],
    channel_theta: [f64; 2],
    residual: [f64; 2],
    probes: [usize; 2],
}

struct ContactEquations {
    reference: [f64; 2],
    resistance: [[f64; 2]; 2],
    laws: [CounterflowContactLaw; 2],
    reference_temperature_k: f64,
    minimum_slope: [f64; 2],
}
impl ContactEquations {
    fn channel_theta(&self, q: [f64; 2]) -> [f64; 2] {
        std::array::from_fn(|i| {
            self.reference[i] - self.resistance[i][0] * q[0] - self.resistance[i][1] * q[1]
        })
    }

    fn residual(
        &self,
        i: usize,
        point: CaloricPoint,
        other_q: f64,
    ) -> Result<f64, CounterflowHeatError> {
        let theta = self.reference[i]
            - self.resistance[i][i] * point.q
            - self.resistance[i][1 - i] * other_q;
        let residual = point.q
            - self.laws[i].probe_heat_w(
                self.reference_temperature_k,
                theta,
                point.thermal_theta,
            )?;
        if residual.is_finite() && theta.is_finite() {
            Ok(residual)
        } else {
            Err(CounterflowHeatError::Invalid(
                "paired contact residual is unrepresentable",
            ))
        }
    }
}

fn solve_pair(
    evaluate: impl Fn(usize, f64) -> Result<CaloricPoint, CounterflowHeatError>,
    reference: [f64; 2],
    resistance: [[f64; 2]; 2],
    laws: [CounterflowContactLaw; 2],
    reference_temperature_k: f64,
    controls: CounterflowHeatControls,
) -> Result<PairSolution, CounterflowHeatError> {
    let minimum_slope = std::array::from_fn(|i| laws[i].minimum_probe_slope(resistance[i][i]));
    for i in 0..2 {
        let maximum_tangent = laws[i].maximum_tangent_w_k()?;
        let diagonal_slope = 1.0 + maximum_tangent * resistance[i][i];
        let cross_slope = maximum_tangent * resistance[i][1 - i];
        // For any actual channel-side tangent g in [0,g_max], the scalar
        // response map's derivative is bounded by g*R_ij/(1+g*R_ii).
        // The finite caloric response only increases its denominator. This
        // bound also covers numerical probes with a flat interface potential.
        if !positive(minimum_slope[i])
            || !positive(diagonal_slope)
            || !nonnegative(cross_slope)
            || diagonal_slope <= cross_slope
        {
            return Err(CounterflowHeatError::Invalid(
                "paired thermal response must be finite and strictly diagonally dominant",
            ));
        }
    }
    let equations = ContactEquations {
        reference,
        resistance,
        laws,
        reference_temperature_k,
        minimum_slope,
    };
    let mut point = [evaluate(0, 0.0)?, evaluate(1, 0.0)?];
    let mut probes = [1; 2];
    let mut domain_refused = false;
    for _ in 0..controls.max_iterations {
        let last_residual = [
            equations.residual(0, point[0], point[1].q)?,
            equations.residual(1, point[1], point[0].q)?,
        ];
        let largest = last_residual[0].abs().max(last_residual[1].abs());
        if largest <= controls.contact_power_tolerance_w {
            return Ok(PairSolution {
                point,
                channel_theta: equations.channel_theta(point.map(|p| p.q)),
                residual: last_residual,
                probes,
            });
        }
        // Avoid solving the first coordinate to final precision while the
        // other contact is still far away. Monotone diagonal dominance makes
        // the exact coordinate maps contractive; the local tolerance tightens
        // with their joint residual, with no fitted caloric derivative.
        let local_tolerance = (0.05 * largest).max(0.25 * controls.contact_power_tolerance_w);
        let before_probes = probes;
        for i in 0..2 {
            point[i] = solve_coordinate(
                i,
                point[i],
                point[1 - i].q,
                &equations,
                &evaluate,
                controls,
                local_tolerance,
                &mut probes[i],
                &mut domain_refused,
            )?;
        }
        if probes == before_probes {
            break;
        }
    }
    // Test the final pair before classifying any bounded-solve refusal. An
    // intermediate coordinate outside a power/chart limit proves nothing
    // about the eventual jointly coupled root.
    let last_residual = [
        equations.residual(0, point[0], point[1].q)?,
        equations.residual(1, point[1], point[0].q)?,
    ];
    let largest = last_residual[0].abs().max(last_residual[1].abs());
    if largest <= controls.contact_power_tolerance_w {
        return Ok(PairSolution {
            point,
            channel_theta: equations.channel_theta(point.map(|p| p.q)),
            residual: last_residual,
            probes,
        });
    }
    for i in 0..2 {
        if point[i].q.abs() == controls.max_abs_heat_power_w
            && point[i].q * last_residual[i] < 0.0
            && last_residual[i].abs() > controls.contact_power_tolerance_w
        {
            return Err(CounterflowHeatError::PowerLimit);
        }
    }
    if domain_refused {
        Err(DepositionError::MaterialDomain.into())
    } else {
        Err(CounterflowHeatError::Unresolved {
            residual_w: largest,
            iterations: probes[0].max(probes[1]),
        })
    }
}

#[allow(clippy::too_many_arguments)]
fn solve_coordinate(
    i: usize,
    initial: CaloricPoint,
    other_q: f64,
    equations: &ContactEquations,
    evaluate: &impl Fn(usize, f64) -> Result<CaloricPoint, CounterflowHeatError>,
    controls: CounterflowHeatControls,
    tolerance: f64,
    probes: &mut usize,
    domain_refused: &mut bool,
) -> Result<CaloricPoint, CounterflowHeatError> {
    let initial_residual = equations.residual(i, initial, other_q)?;
    if initial_residual.abs() <= tolerance || *probes == controls.max_iterations {
        return Ok(initial);
    }
    // The law's global probe slope and the original caloric monotonicity bound
    // the scalar root from any admitted cached coordinate. For a nonlinear
    // interface this includes flat numerical continuation at its domain edges.
    let bound = (initial.q - initial_residual / equations.minimum_slope[i]).clamp(
        -controls.max_abs_heat_power_w,
        controls.max_abs_heat_power_w,
    );
    if bound == initial.q {
        return Ok(initial);
    }
    let mut near = initial;
    let mut near_residual = initial_residual;
    let mut far_q = bound;
    let mut far_residual: Option<f64> = None;
    let mut best = initial;
    let mut best_residual = initial_residual.abs();
    // Bounded inexact coordinates keep one temporarily chart-limited end from
    // consuming the full budget before the other end can move its target.
    for iteration in 0..8 {
        if *probes == controls.max_iterations {
            break;
        }
        let q = if iteration == 0 {
            far_q
        } else if let Some(residual) = far_residual.filter(|_| iteration % 4 != 0) {
            let a = near_residual.abs();
            let b = residual.abs();
            let fraction = if a <= b {
                (a / b) / (1.0 + a / b)
            } else {
                1.0 / (1.0 + b / a)
            };
            let secant = (1.0 - fraction) * near.q + fraction * far_q;
            if secant != near.q && secant != far_q {
                secant
            } else {
                0.5 * near.q + 0.5 * far_q
            }
        } else {
            0.5 * near.q + 0.5 * far_q
        };
        if q == near.q || (iteration != 0 && q == far_q) {
            break;
        }
        *probes += 1;
        let trial = match evaluate(i, q) {
            Ok(trial) => trial,
            Err(CounterflowHeatError::Thermal(DepositionError::MaterialDomain)) => {
                *domain_refused = true;
                far_q = q;
                far_residual = None;
                continue;
            }
            Err(error) => return Err(error),
        };
        let residual = equations.residual(i, trial, other_q)?;
        if residual.abs() < best_residual {
            best = trial;
            best_residual = residual.abs();
        }
        if residual.abs() <= tolerance {
            return Ok(trial);
        }
        if residual.signum() == initial_residual.signum() {
            near = trial;
            near_residual = residual;
            if q == bound {
                // A projected coordinate at the power limit remains eligible
                // for an inward correction after the opposite end advances.
                return Ok(best);
            }
        } else {
            far_q = q;
            far_residual = Some(residual);
        }
    }
    Ok(best)
}
