//! Native finite-body rolling and its original sphere's finite caloric state.
//!
//! Viscoelastic rolling generates bulk heat in the deforming sphere, not at a
//! zero-capacity interface. This owner deposits the original moment work into
//! one explicitly bound whole-sphere cell and retains ordinary network links.
//! The original rolling card has constant material parameters with a declared
//! temperature/frequency envelope. Actual endpoint temperature admits those
//! same parameters; no viscosity curve, material refresh or heat split is added.

use fs_exec::CancelGate;
use fs_mbd::tangential_friction::RigidContactBinding;
use fs_mbd::tangential_friction::rolling::{
    RigidRollingContact, RigidRollingControls, RigidRollingError, RigidRollingHistory,
    RigidRollingReport, apply_rigid_rolling,
};
use fs_mbd::{MassProperties, RigidBodyState};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::{BackwardEulerStep, StepCoefficients};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionSimulation, EnergyAudit, MaterialId,
};
use fs_tribo::rolling_loss::viscoelastic::DeformingSphereSurface;
use fs_tribo::rolling_loss::{RollingKinematics, RollingLossLaw};
use std::cell::Cell;

/// Explicit correspondence between the native sphere and its original cell.
/// The cell represents the entire sphere at one bulk temperature and must have
/// the native sphere's mass. Existing thermal body/material IDs are checked;
/// the caller declares that the rolling card describes this material, since
/// rigid `MassProperties` has no material identity to independently compare.
/// No other active cell may carry this thermal body ID. Unborn scheduled cells
/// hold no accepted energy yet; the original thermal owner handles activation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RigidRollingHeatBinding {
    pub body_index: usize,
    pub cell_id: CellId,
    pub thermal_body_id: u64,
    pub material_id: MaterialId,
}

/// Native mechanics, finite-cell nonlinear work and independent joule gates.
#[derive(Debug, Clone, Copy)]
pub struct RigidRollingHeatControls {
    pub mechanical: RigidRollingControls,
    pub maximum_cells: usize,
    pub maximum_trials: usize,
    pub maximum_newton_iterations: usize,
    pub difference_step_j: f64,
    pub cell_energy_tolerance_j: f64,
    pub thermal_energy_tolerance_j: f64,
    pub combined_energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
}

#[derive(Debug)]
pub enum RigidRollingHeatError<E> {
    Invalid(&'static str),
    Mechanical(RigidRollingError),
    Thermal(DepositionError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
}
impl<E> From<RigidRollingError> for RigidRollingHeatError<E> {
    fn from(error: RigidRollingError) -> Self {
        Self::Mechanical(error)
    }
}
impl<E> From<DepositionError> for RigidRollingHeatError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E> From<NewtonError> for RigidRollingHeatError<E> {
    fn from(error: NewtonError) -> Self {
        Self::Newton(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for RigidRollingHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native rolling/bulk heat: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for RigidRollingHeatError<E> {}

/// Same-endpoint physical rolling work and original caloric energy account.
pub struct RigidRollingHeatReport<P> {
    pub mechanical: RigidRollingReport,
    /// Original inner thermal audit. Its source input includes rolling heat;
    /// use `external_heat_j` below for the complete mechanics/thermal system.
    pub thermal: EnergyAudit,
    pub binding: RigidRollingHeatBinding,
    pub sphere_temperature_before_k: f64,
    pub sphere_temperature_after_k: f64,
    /// Authored boundaries and caller's other heat, measured before adding D.
    pub external_heat_j: f64,
    /// Actual caloric change minus external heat minus native rolling D.
    pub heat_transfer_residual_j: f64,
    /// Actual contacted-pair kinetic change plus caloric change minus other heat.
    pub combined_energy_residual_j: f64,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub trials: usize,
    /// Actual shared interval, from the original thermal clock and endpoint.
    pub duration_s: f64,
    pub target_time_s: f64,
}

/// Complete private candidate. Publish both original indexed bodies, native
/// rolling history and thermal state together in the surrounding integrator.
/// Unrelated bodies are not copied, evaluated or changed by this owner.
pub struct RigidRollingHeatCandidate<'m, P> {
    body_indices: [usize; 2],
    bodies: [RigidBodyState; 2],
    history: RigidRollingHistory,
    thermal: DepositionSimulation<'m>,
    report: RigidRollingHeatReport<P>,
}
impl<'m, P> RigidRollingHeatCandidate<'m, P> {
    pub const fn body_indices(&self) -> [usize; 2] {
        self.body_indices
    }
    /// Candidate native states in the original contact's A/B order.
    pub const fn bodies(&self) -> &[RigidBodyState; 2] {
        &self.bodies
    }
    pub const fn history(&self) -> &RigidRollingHistory {
        &self.history
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    pub const fn report(&self) -> &RigidRollingHeatReport<P> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        [usize; 2],
        [RigidBodyState; 2],
        RigidRollingHistory,
        DepositionSimulation<'m>,
        RigidRollingHeatReport<P>,
    ) {
        (
            self.body_indices,
            self.bodies,
            self.history,
            self.thermal,
            self.report,
        )
    }
}

struct Endpoint<P> {
    temperature_k: f64,
    external_heat_j: f64,
    source: P,
}

/// Propose physical rolling recoil and whole-sphere bulk heating together.
///
/// Every accepted owner is borrowed. `contact.temperature_k` must equal the
/// actual initial sphere-cell temperature. The original constant card supplies
/// one cached native mechanical solve; every thermal trial admits its actual
/// sphere temperature through that same card. After the thermal root converges,
/// original mechanics and history are restaged at its endpoint temperature and
/// exact physical agreement with the cached solution is required.
///
/// `resolve` supplies only other signed heat and original conductance updates,
/// without publishing source receipts or adding rolling heat again. Authored
/// links, sources and boundaries remain owned by the original thermal network.
/// Energy guesses follow original sorted CellId order. The actual interval is
/// `(accepted_time_s + dt_s) - accepted_time_s`, shared by both original owners.
///
/// This is a one-cell uniform-temperature approximation of sphere bulk loss,
/// with fixed card coefficients, normal patch and geometry. It does not resolve
/// spatial loss density, flash temperature, thermoelastic storage, expansion,
/// material changes or temperature-dependent viscosity. Endpoint heat transport
/// is first order and requires time/spatial refinement independently of Newton.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn propose_rigid_rolling_heat<'m, P, E>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    history: &RigidRollingHistory,
    contact: RigidRollingContact<'_>,
    thermal: &DepositionSimulation<'m>,
    accepted_time_s: f64,
    binding: RigidRollingHeatBinding,
    dt_s: f64,
    controls: RigidRollingHeatControls,
    energy_guess_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>) -> Result<StepCoefficients<P>, E>,
    gate: &CancelGate,
) -> Result<RigidRollingHeatCandidate<'m, P>, RigidRollingHeatError<E>> {
    let checkpoint = || {
        if gate.is_requested() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    checkpoint()?;
    controls.newton.validate()?;
    let cells = thermal.cells().len();
    let target_time_s = accepted_time_s + dt_s;
    if !accepted_time_s.is_finite()
        || accepted_time_s != thermal.time_s()
        || !dt_s.is_finite()
        || dt_s <= 0.0
        || !target_time_s.is_finite()
        || target_time_s <= accepted_time_s
        || bodies.len() != properties.len()
        || controls.mechanical.maximum_bodies == 0
        || bodies.len() > controls.mechanical.maximum_bodies
        || contact.binding.body_a >= bodies.len()
        || contact.binding.body_b >= bodies.len()
        || contact.binding.body_a == contact.binding.body_b
        || cells == 0
        || cells > controls.maximum_cells
        || controls.maximum_trials == 0
        || controls.maximum_newton_iterations == 0
        || !controls.difference_step_j.is_finite()
        || controls.difference_step_j <= 0.0
        || [
            controls.cell_energy_tolerance_j,
            controls.thermal_energy_tolerance_j,
            controls.combined_energy_tolerance_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
        || energy_guess_j.len() != cells
        || energy_guess_j.iter().any(|x| !x.is_finite())
    {
        return Err(RigidRollingHeatError::Invalid(
            "bounded original bodies/cells and one synchronized finite interval required",
        ));
    }
    let RollingLossLaw::ViscoelasticSphere(card) = contact.law else {
        return Err(RigidRollingHeatError::Invalid(
            "bulk sphere heat requires the original viscoelastic sphere card",
        ));
    };
    let sphere_body = match card.parameters().sphere_surface {
        DeformingSphereSurface::SurfaceA => contact.binding.body_a,
        DeformingSphereSurface::SurfaceB => contact.binding.body_b,
    };
    let cell = thermal
        .cells()
        .get(&binding.cell_id)
        .ok_or(RigidRollingHeatError::Invalid("inactive sphere bulk cell"))?;
    if binding.body_index != sphere_body
        || cell.body_id != binding.thermal_body_id
        || cell.material_id != binding.material_id
        || cell.mass_kg != properties[sphere_body].mass()
    {
        return Err(RigidRollingHeatError::Invalid(
            "original whole-sphere cell mass and declared body/material identities required",
        ));
    }
    // The admitted cell already carries the sphere's entire mass and capacity.
    // Inspect only active inventory, within the previously admitted cell bound;
    // the original backward-Euler owner separately enforces event boundaries.
    for (&id, other) in thermal.cells() {
        checkpoint()?;
        if id != binding.cell_id && other.body_id == binding.thermal_body_id {
            return Err(RigidRollingHeatError::Invalid(
                "whole-sphere thermal body must have exactly one active caloric cell",
            ));
        }
    }
    let temperature_before = thermal.temperature_k(binding.cell_id)?;
    if contact.temperature_k != temperature_before {
        return Err(RigidRollingHeatError::Invalid(
            "rolling temperature must belong to the actual accepted sphere cell",
        ));
    }
    let thermal_step = BackwardEulerStep::new(thermal, target_time_s, checkpoint)?;
    let duration_s = thermal_step.duration_s();
    let body_indices = [contact.binding.body_a, contact.binding.body_b];
    let accepted_pair = body_indices.map(|i| bodies[i]);
    let pair_properties = body_indices.map(|i| properties[i]);
    let pair_contact = RigidRollingContact {
        binding: RigidContactBinding {
            body_a: 0,
            body_b: 1,
            ..contact.binding
        },
        ..contact
    };
    let mut cached_bodies = accepted_pair;
    let mut cached_history = history.clone();
    let mut cached = apply_rigid_rolling(
        &mut cached_bodies,
        &pair_properties,
        &mut cached_history,
        pair_contact,
        duration_s,
        controls.mechanical,
        gate,
    )?;
    cached.deforming_body = sphere_body;
    let heat_j = cached.dissipated_heat_j;
    let heat_w = finite(heat_j / duration_s)?;
    if heat_j > 0.0 && heat_w == 0.0 {
        return Err(DepositionError::Arithmetic.into());
    }
    let initial_rate =
        cached.rolling_rate_before_rad_s[0].hypot(cached.rolling_rate_before_rad_s[1]);
    let trials = Cell::new(0usize);
    let evaluate = |x: &[f64]| {
        checkpoint()?;
        if trials.get() >= controls.maximum_trials {
            return Err(RigidRollingHeatError::TrialBudget);
        }
        trials.set(trials.get() + 1);
        thermal_step.evaluate(
            x,
            |endpoint| {
                let temperature_k = endpoint.temperature_k(binding.cell_id)?;
                // Constants are immutable. This original response enforces the
                // actual trial temperature/frequency and full Hertz envelope;
                // native substeps already checked all rolling rates once.
                let kinematics = RollingKinematics::new(
                    card.parameters().radius_m * initial_rate,
                    initial_rate,
                    0.0,
                    temperature_k,
                    contact.excitation_frequency_hz,
                    duration_s,
                )
                .map_err(RigidRollingError::from)?;
                card.response(contact.patch, contact.interface, kinematics)
                    .map_err(RigidRollingError::from)?;
                let StepCoefficients {
                    mut sources_w,
                    conductances_w_k,
                    payload,
                } = resolve(endpoint).map_err(RigidRollingHeatError::Source)?;
                checkpoint()?;
                let external_heat_j = finite(
                    duration_s
                        * endpoint.external_step_power_w(
                            &sources_w,
                            &conductances_w_k,
                            checkpoint,
                        )?,
                )?;
                if heat_w != 0.0 {
                    let power = sources_w.entry(binding.cell_id).or_default();
                    let combined = finite(*power + heat_w)?;
                    if combined == *power {
                        return Err(DepositionError::Arithmetic.into());
                    }
                    *power = combined;
                }
                Ok(StepCoefficients {
                    sources_w,
                    conductances_w_k,
                    payload: Endpoint {
                        temperature_k,
                        external_heat_j,
                        source: payload,
                    },
                })
            },
            checkpoint,
        )
    };
    let problem =
        FallibleFiniteDifference::new(cells, controls.difference_step_j, |x: &[f64]| {
            evaluate(x).map(|trial| trial.residual_j().to_vec())
        })
        .map_err(numerical)?;
    let construction = NewtonKrylovState::new(&problem, energy_guess_j.to_vec(), controls.newton);
    if let Some(error) = problem.take_failure() {
        return Err(numerical(error));
    }
    let mut solver = construction?;
    let result = solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
        gate.is_requested() || problem.has_failed()
    });
    if let Some(error) = problem.take_failure() {
        return Err(numerical(error));
    }
    checkpoint()?;
    if result.progress != SolverRunProgress::Complete {
        return Err(RigidRollingHeatError::Interrupted(result.progress));
    }
    if !result.report.converged {
        return Err(RigidRollingHeatError::NotConverged(result.report));
    }
    // Retain the exact source evaluation at the final root, never a probe.
    let trial = evaluate(&solver.x)?;
    let mut next_bodies = accepted_pair;
    let mut next_history = history.clone();
    let mut mechanical = apply_rigid_rolling(
        &mut next_bodies,
        &pair_properties,
        &mut next_history,
        RigidRollingContact {
            temperature_k: trial.payload().temperature_k,
            ..pair_contact
        },
        duration_s,
        controls.mechanical,
        gate,
    )?;
    mechanical.deforming_body = sphere_body;
    if next_bodies != cached_bodies
        || next_history.state() != cached_history.state()
        || !same_physics(&mechanical, &cached)
    {
        return Err(RigidRollingHeatError::Invalid(
            "original fixed rolling card changed its cached physical response",
        ));
    }
    let heat_transfer_residual_j =
        finite(trial.audit().stored_energy_change_j - trial.payload().external_heat_j - heat_j)?;
    let combined_energy_residual_j = finite(
        (mechanical.kinetic_after_j - mechanical.kinetic_before_j)
            + trial.audit().stored_energy_change_j
            - trial.payload().external_heat_j,
    )?;
    if heat_transfer_residual_j.abs() > controls.thermal_energy_tolerance_j
        || combined_energy_residual_j.abs() > controls.combined_energy_tolerance_j
    {
        return Err(RigidRollingHeatError::Invalid(
            "native rolling work and actual bulk caloric change do not balance",
        ));
    }
    let (thermal, endpoint, thermal_audit) = trial.accept(
        controls.cell_energy_tolerance_j,
        controls.thermal_energy_tolerance_j,
        checkpoint,
    )?;
    checkpoint()?;
    Ok(RigidRollingHeatCandidate {
        body_indices,
        bodies: next_bodies,
        history: next_history,
        thermal,
        report: RigidRollingHeatReport {
            mechanical,
            thermal: thermal_audit,
            binding,
            sphere_temperature_before_k: temperature_before,
            sphere_temperature_after_k: endpoint.temperature_k,
            external_heat_j: endpoint.external_heat_j,
            heat_transfer_residual_j,
            combined_energy_residual_j,
            source: endpoint.source,
            nonlinear: result.report,
            trials: trials.get(),
            duration_s,
            target_time_s,
        },
    })
}

/// Compare physical response, retaining the final source/checkpoint receipts.
fn same_physics(a: &RigidRollingReport, b: &RigidRollingReport) -> bool {
    a.deforming_body == b.deforming_body
        && a.mobility == b.mobility
        && a.rolling_mobility == b.rolling_mobility
        && a.rolling_rate_before_rad_s == b.rolling_rate_before_rad_s
        && a.rolling_rate_after_rad_s == b.rolling_rate_after_rad_s
        && a.traction_impulse_world_ns == b.traction_impulse_world_ns
        && a.rolling_angular_impulse_world_nms == b.rolling_angular_impulse_world_nms
        && a.substeps == b.substeps
        && a.dissipated_heat_j == b.dissipated_heat_j
        && a.kinetic_before_j == b.kinetic_before_j
        && a.kinetic_after_j == b.kinetic_after_j
        && a.midpoint_impulse_work_j == b.midpoint_impulse_work_j
        && a.static_reaction_work_j == b.static_reaction_work_j
        && a.linear_momentum_change_ns == b.linear_momentum_change_ns
        && a.angular_momentum_change_nms == b.angular_momentum_change_nms
        && a.final_step.wrench == b.final_step.wrench
        && a.final_step.stored_energy_j == b.final_step.stored_energy_j
        && a.final_step.storage_change_j == b.final_step.storage_change_j
        && a.final_step.dissipation == b.final_step.dissipation
}
fn finite(value: f64) -> Result<f64, DepositionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic)
    }
}
fn numerical<E>(
    error: FiniteDifferenceError<RigidRollingHeatError<E>>,
) -> RigidRollingHeatError<E> {
    match error {
        FiniteDifferenceError::Source(source) => source,
        FiniteDifferenceError::Invalid(message) => RigidRollingHeatError::Invalid(message),
    }
}
