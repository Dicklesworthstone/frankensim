//! Adhesive sliding/spin work and finite thermal transport on one body island.
//!
//! Each normal owner supplies its actual attached disk. Each thermal trial
//! refreshes that disk's original selected junction-strength claims at its own
//! solved interface temperature. The complete prepared body response is reused;
//! simultaneous impulses and all finite-cell temperatures publish together.
//! Normal coefficients and geometry are held for this split phase. Their
//! thermal evolution and reversible surface/normal energy remain with the normal
//! owner, and are not silently inferred from the tangential heat calculation.

/// Curved finite-time slip/spin, arrest and original finite-cell heat together.
pub mod transient;
/// Original normal snap release and explicitly declared finite-cell heat.
pub mod snap;

use super::FrictionHeatRoute;
use super::network::{
    FrictionHeatNetworkInput, FrictionHeatNetworkReport, FrictionHeatNetworkStep,
    InterfaceHeatInput,
};
use fs_contact::interface_binding::adhesive_shear::spinning::temperature::{
    SourceThermalAdhesivePatch, SourceThermalAdhesivePatchError,
};
use fs_exec::CancelGate;
use fs_math::det;
use fs_mbd::tangential_friction::adhesive::RigidAdhesiveContact;
use fs_mbd::tangential_friction::adhesive::multiple::{
    PreparedRigidAdhesiveIsland, RigidAdhesiveIslandCandidate, RigidAdhesiveIslandContact,
    RigidAdhesiveIslandError, RigidAdhesiveIslandThermalLaw,
};
use fs_mbd::{MassProperties, RigidBodyState};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{NewtonError, NewtonKrylovState, NewtonReport, SolverRunProgress};
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use std::cell::Cell;

/// Endpoint conductivities keyed by the original thermal contact link.
pub use super::hertz::multiple::RigidHertzIslandHeatConductivity as RigidAdhesiveIslandHeatConductivity;
/// Shared dimensional island/root/energy controls contain no normal-load law.
pub use super::hertz::multiple::RigidHertzIslandHeatControls as RigidAdhesiveIslandHeatControls;
/// Original bulk transport and other heat; adhesive work has its own owner.
pub use super::hertz::multiple::RigidHertzIslandHeatSource as RigidAdhesiveIslandHeatSource;

/// One original adhesive disk and its coincident physical thermal face.
#[derive(Clone, Copy)]
pub struct RigidAdhesiveIslandHeatContact<'a> {
    /// Unique mechanical identity, independent of caller enumeration.
    pub patch_id: &'a str,
    /// Original source, body pair, arms and tangent frame.
    pub contact: RigidAdhesiveContact<'a>,
    /// Original admitted absolute-temperature axis of the junction strengths.
    pub temperature_axis: &'a str,
    /// The same disk's existing thermal link, in the body's A/B order.
    pub route: &'a FrictionHeatRoute,
    /// Positive solver seed; the actual interface temperature is an unknown.
    pub temperature_guess_k: f64,
}

#[derive(Debug)]
pub enum RigidAdhesiveIslandHeatError<E> {
    Invalid(&'static str),
    Mechanical(RigidAdhesiveIslandError),
    FrictionSource(SourceThermalAdhesivePatchError),
    Thermal(DepositionError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
}
impl<E> From<RigidAdhesiveIslandError> for RigidAdhesiveIslandHeatError<E> {
    fn from(error: RigidAdhesiveIslandError) -> Self {
        Self::Mechanical(error)
    }
}
impl<E> From<SourceThermalAdhesivePatchError> for RigidAdhesiveIslandHeatError<E> {
    fn from(error: SourceThermalAdhesivePatchError) -> Self {
        Self::FrictionSource(error)
    }
}
impl<E> From<DepositionError> for RigidAdhesiveIslandHeatError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E> From<NewtonError> for RigidAdhesiveIslandHeatError<E> {
    fn from(error: NewtonError) -> Self {
        Self::Newton(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for RigidAdhesiveIslandHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adhesive island friction/heat: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for RigidAdhesiveIslandHeatError<E> {}

/// Original disk identity and the exact endpoint junction-property usage.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidAdhesiveIslandHeatReceipt<'a> {
    pub patch_id: String,
    pub friction_source: SourceThermalAdhesivePatch<'a>,
}

/// All body impulses, source receipts and the caloric endpoint accept together.
pub struct HeatedRigidAdhesiveIslandCandidate<'a, 'm, P> {
    pub mechanical: RigidAdhesiveIslandCandidate,
    pub thermal: DepositionSimulation<'m>,
    /// Ascending original link order, including signed heat into each body.
    pub heat: FrictionHeatNetworkReport,
    /// Same thermal order, with explicit IDs preserving mechanical mapping.
    pub friction_sources: Vec<RigidAdhesiveIslandHeatReceipt<'a>>,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub trials: usize,
    pub target_time_s: f64,
    /// Actual kinetic plus caloric change, minus external heat.
    pub combined_energy_residual_j: f64,
    /// Independent sum of mechanical and thermal audit residuals.
    pub summed_energy_residual_j: f64,
}

struct Staged<'a, P> {
    mechanical: RigidAdhesiveIslandCandidate,
    friction_sources: Vec<RigidAdhesiveIslandHeatReceipt<'a>>,
    source: P,
}

/// Solve attached adhesive patches and their shared finite heat cells together.
///
/// Adhesion at zero or tensile net normal load remains active whenever the
/// original normal response has a finite solid-contact area. A detached disk
/// cannot own a finite conducting face; normal/topology owners must remove that
/// route before this phase. The attractive noncontact annulus is never heated
/// as a solid junction. Each route retains the original body's A/B assignment.
///
/// Interface temperatures are solved from actual sliding PLUS spin work and
/// the original per-side k*A/d conductances. Body temperatures may differ and
/// heat into one body may be negative. No heat fraction, added heat store,
/// Coulomb normal load or artificial stopping loss is introduced. The callback
/// supplies endpoint bulk conductivities and OTHER heat only, staging receipts
/// until the complete candidate is accepted.
///
/// This is the frozen-normal/geometry split of a thermoadhesive calculation:
/// only junction strengths refresh here. Normal surface energy, changing area,
/// pose evolution and stopping-event subdivision have their own owners. The
/// native event refusals remain intact. The orthogonal finite-cell, perfect-face
/// and backward-Euler thermal model requires spatial/time refinement; nonlinear
/// convergence does not establish uniqueness or experimental validation.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn propose_rigid_adhesive_island_heat<'a, 'm, P, E>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[RigidAdhesiveIslandHeatContact<'a>],
    thermal: &DepositionSimulation<'m>,
    accepted_time_s: f64,
    dt_s: f64,
    controls: RigidAdhesiveIslandHeatControls,
    energy_guess_j: &[f64],
    resolve: impl Fn(
        &DepositionSimulation<'m>,
        &[u64],
        &[f64],
    ) -> Result<RigidAdhesiveIslandHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedRigidAdhesiveIslandCandidate<'a, 'm, P>, RigidAdhesiveIslandHeatError<E>> {
    let checkpoint = || {
        if gate.is_requested() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    checkpoint()?;
    controls.newton.validate()?;
    controls
        .mechanical
        .validate()
        .map_err(RigidAdhesiveIslandError::from)?;
    let count = contacts.len();
    let cells = thermal.cells().len();
    let dimension = cells
        .checked_add(count)
        .ok_or(RigidAdhesiveIslandHeatError::Invalid(
            "thermal dimension overflow",
        ))?;
    let target = accepted_time_s + dt_s;
    if !accepted_time_s.is_finite()
        || accepted_time_s != thermal.time_s()
        || !dt_s.is_finite()
        || dt_s <= 0.0
        || !target.is_finite()
        || target <= accepted_time_s
        || count == 0
        || count > controls.mechanical.maximum_contacts
        || controls.maximum_cells == 0
        || cells > controls.maximum_cells
        || dimension > controls.maximum_unknowns
        || controls.maximum_trials == 0
        || controls.maximum_newton_iterations == 0
        || !controls.difference_step_j.is_finite()
        || controls.difference_step_j <= 0.0
        || !controls.temperature_scale_j_per_k.is_finite()
        || controls.temperature_scale_j_per_k <= 0.0
        || energy_guess_j.len() != cells
        || energy_guess_j.iter().any(|v| !v.is_finite())
        || contacts.iter().any(|c| {
            !c.temperature_guess_k.is_finite()
                || c.temperature_guess_k <= 0.0
                || c.contact.source.contact_area_m2() <= 0.0
        })
        || [
            controls.heat.cell_j,
            controls.heat.thermal_j,
            controls.heat.interface_j,
            controls.heat.total_j,
            controls.combined_energy_tolerance_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
    {
        return Err(RigidAdhesiveIslandHeatError::Invalid(
            "bounded attached disks and synchronized finite thermal coordinates required",
        ));
    }
    let duration_s = target - accepted_time_s;
    let mut ordered: Vec<_> = contacts.iter().collect();
    ordered.sort_by_key(|c| c.route.link_id());
    if ordered
        .windows(2)
        .any(|p| p[0].route.link_id() == p[1].route.link_id())
    {
        return Err(RigidAdhesiveIslandHeatError::Invalid(
            "each original thermal link requires one mechanical heat owner",
        ));
    }
    let native_contacts: Vec<_> = ordered
        .iter()
        .map(|c| RigidAdhesiveIslandContact {
            patch_id: c.patch_id,
            contact: c.contact,
        })
        .collect();
    let prepared = PreparedRigidAdhesiveIsland::new(
        bodies,
        properties,
        &native_contacts,
        duration_s,
        controls.mechanical,
        gate,
    )?;
    let mut mechanical_indices = Vec::with_capacity(count);
    for c in &ordered {
        checkpoint()?;
        let index = prepared
            .contacts()
            .binary_search_by(|p| p.patch_id.as_str().cmp(c.patch_id))
            .map_err(|_| {
                RigidAdhesiveIslandHeatError::Invalid("missing original attached patch")
            })?;
        mechanical_indices.push(index);
        let point = prepared.contacts()[index].contact_points_world[0];
        validate_geometry(
            [point.x, point.y, point.z],
            c.contact,
            c.route,
            controls.mechanical.mechanical.geometry_tolerance_m,
        )
        .map_err(RigidAdhesiveIslandHeatError::Invalid)?;
    }
    let routes: Vec<_> = ordered.iter().map(|c| c.route).collect();
    let step = FrictionHeatNetworkStep::new(
        thermal,
        target,
        &routes,
        controls.mechanical.maximum_contacts,
        checkpoint,
    )?;
    let trials = Cell::new(0usize);
    let evaluate = |x: &[f64], entropy: bool| {
        checkpoint()?;
        if trials.get() >= controls.maximum_trials {
            return Err(RigidAdhesiveIslandHeatError::TrialBudget);
        }
        trials.set(trials.get() + 1);
        if x.len() != dimension || x.iter().any(|v| !v.is_finite()) {
            return Err(RigidAdhesiveIslandHeatError::Invalid(
                "complete finite thermal iterate required",
            ));
        }
        let temperatures: Vec<_> = ordered
            .iter()
            .enumerate()
            .map(|(i, c)| c.temperature_guess_k + x[cells + i] / controls.temperature_scale_j_per_k)
            .collect();
        let endpoint = |endpoint: &DepositionSimulation<'m>, temperatures: &[f64]| {
            let mut receipts = Vec::with_capacity(count);
            for (c, &temperature) in ordered.iter().zip(temperatures) {
                checkpoint()?;
                receipts.push(RigidAdhesiveIslandHeatReceipt {
                    patch_id: c.patch_id.to_owned(),
                    friction_source: c
                        .contact
                        .source
                        .at_temperature(c.temperature_axis, temperature)?,
                });
            }
            let sources: Vec<_> = receipts
                .iter()
                .map(|r| RigidAdhesiveIslandThermalLaw {
                    patch_id: &r.patch_id,
                    source: &r.friction_source,
                })
                .collect();
            let mechanical = prepared.advance_thermal(&sources, gate)?;
            let RigidAdhesiveIslandHeatSource {
                interfaces,
                transport,
            } = resolve(endpoint, step.link_ids(), temperatures)
                .map_err(RigidAdhesiveIslandHeatError::Source)?;
            if interfaces.len() != count
                || interfaces
                    .iter()
                    .zip(step.link_ids())
                    .any(|(source, id)| source.link_id != *id)
            {
                return Err(RigidAdhesiveIslandHeatError::Invalid(
                    "conductivities must retain every canonical original link identity",
                ));
            }
            let mut inputs = Vec::with_capacity(count);
            for ((c, &index), source) in ordered.iter().zip(&mechanical_indices).zip(interfaces) {
                checkpoint()?;
                let native = &mechanical.report.mechanical.contacts[index];
                if native.patch_id != c.patch_id {
                    return Err(RigidAdhesiveIslandHeatError::Invalid(
                        "changed original adhesive heat owner",
                    ));
                }
                inputs.push(InterfaceHeatInput {
                    dissipated_j: native.dissipated_energy_j,
                    conductivity_w_m_k: source.conductivity_w_m_k,
                });
            }
            let StepCoefficients {
                sources_w,
                conductances_w_k,
                payload,
            } = transport;
            Ok(FrictionHeatNetworkInput {
                interfaces: inputs,
                transport: StepCoefficients {
                    sources_w,
                    conductances_w_k,
                    payload: Staged {
                        mechanical,
                        friction_sources: receipts,
                        source: payload,
                    },
                },
            })
        };
        if entropy {
            step.evaluate_with_entropy(&x[..cells], &temperatures, endpoint, checkpoint)
        } else {
            step.evaluate(&x[..cells], &temperatures, endpoint, checkpoint)
        }
    };
    let problem =
        FallibleFiniteDifference::new(dimension, controls.difference_step_j, |x: &[f64]| {
            evaluate(x, false).map(|trial| trial.residual_j().to_vec())
        })
        .map_err(numerical)?;
    let mut initial = energy_guess_j.to_vec();
    initial.resize(dimension, 0.0);
    let construction = NewtonKrylovState::new(&problem, initial, controls.newton);
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
        return Err(RigidAdhesiveIslandHeatError::Interrupted(result.progress));
    }
    if !result.report.converged {
        return Err(RigidAdhesiveIslandHeatError::NotConverged(result.report));
    }
    let trial = evaluate(&solver.x, controls.audit_entropy)?;
    let native = &trial.payload().mechanical.report.mechanical;
    let heat = trial.report();
    let combined = native.kinetic_change_j
        + heat.thermal.stored_energy_change_j
        - heat.external_heat_j;
    let summed = native.total_energy_residual_j + heat.total_energy_residual_j;
    if [combined, summed]
        .iter()
        .any(|v| !v.is_finite() || v.abs() > controls.combined_energy_tolerance_j)
    {
        return Err(RigidAdhesiveIslandHeatError::Invalid(
            "original adhesive island kinetic and finite-cell energy does not balance",
        ));
    }
    let (thermal, staged, heat) = trial.accept(controls.heat, checkpoint)?;
    checkpoint()?;
    Ok(HeatedRigidAdhesiveIslandCandidate {
        mechanical: staged.mechanical,
        thermal,
        heat,
        friction_sources: staged.friction_sources,
        source: staged.source,
        nonlinear: result.report,
        trials: trials.get(),
        target_time_s: target,
        combined_energy_residual_j: combined,
        summed_energy_residual_j: summed,
    })
}

fn validate_geometry(
    point: [f64; 3],
    contact: RigidAdhesiveContact<'_>,
    route: &FrictionHeatRoute,
    tolerance: f64,
) -> Result<(), &'static str> {
    let geometry = route.geometry;
    let centre = [geometry[0].start_m, geometry[0].end_m]
        .into_iter()
        .find(|p| *p == geometry[1].start_m || *p == geometry[1].end_m)
        .ok_or("thermal route has no shared face")?;
    let area = contact.source.contact_area_m2();
    let thermal_area = route.contact_area_m2();
    if norm(std::array::from_fn(|i| point[i] - centre[i])) > tolerance
        || (area - thermal_area).abs() > 64.0 * f64::EPSILON * area.max(thermal_area)
    {
        return Err(
            "original adhesive disk and thermal face must have the same centre and solid area",
        );
    }
    let n = contact.frame.normal();
    for g in geometry {
        let a = std::array::from_fn::<_, 3, _>(|i| g.end_m[i] - g.start_m[i]);
        if norm([
            a[1] * n[2] - a[2] * n[1],
            a[2] * n[0] - a[0] * n[2],
            a[0] * n[1] - a[1] * n[0],
        ]) > tolerance
        {
            return Err("adhesive normal and finite-cell centrelines differ");
        }
    }
    Ok(())
}
fn norm(v: [f64; 3]) -> f64 {
    det::hypot(det::hypot(v[0], v[1]), v[2])
}
fn numerical<E>(
    error: FiniteDifferenceError<RigidAdhesiveIslandHeatError<E>>,
) -> RigidAdhesiveIslandHeatError<E> {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(field) => RigidAdhesiveIslandHeatError::Invalid(field),
    }
}
