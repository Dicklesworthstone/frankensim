//! Native rigid cohesive motion and finite-body heating accept together.
//!
//! Each thermal face is explicitly attached to the same original material
//! quadrature patch. Geometry and conduction use that fixed reference pairing;
//! the native mechanics still supplies objective finite-rotation work. This is
//! not an Eulerian moving-footprint transfer or a remap of sliding thermal cells.
//! Fracture creates surfaces and stays outside heat. Only original bond viscosity
//! and damaged-area sliding work feed the existing finite enthalpy network.

use super::FrictionHeatRoute;
use super::cohesive::{CohesiveHeatError, solve_frozen_cohesive_heat};
use super::network::{FrictionHeatNetworkReport, FrictionHeatNetworkStep};
use crate::cohesive_contact::rigid::{
    RigidCohesiveError, RigidCohesiveFrame, RigidCohesiveLoad, RigidCohesiveSystem,
};
use fs_exec::CancelGate;
use fs_matdb::MaterialStateId;
use fs_solver::NewtonReport;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation, MaterialId};
use std::collections::BTreeSet;

pub use super::cohesive::{CohesiveHeatConductivity, CohesiveHeatControls, CohesiveHeatSource};

/// Explicit map from one original native material to its thermal cells.
/// The two entries passed to `apply_rigid_cohesive_heat` follow native A/B order.
#[derive(Clone, Copy)]
pub struct RigidCohesiveHeatBody<'a> {
    pub body_id: u64,
    pub thermal_material: MaterialId,
    pub mechanical_material: &'a MaterialStateId,
}

/// One fixed material attachment and its existing reference thermal face.
/// Route A/B follows the original native body A/B, with exactly the same area.
#[derive(Clone, Copy)]
pub struct RigidCohesiveHeatContact<'a> {
    pub site_index: usize,
    pub route: &'a FrictionHeatRoute,
    pub temperature_guess_k: f64,
}

#[derive(Debug)]
pub enum RigidCohesiveHeatError<E> {
    Invalid(&'static str),
    Mechanical(RigidCohesiveError),
    Heat(CohesiveHeatError<E>),
}
impl<E> From<RigidCohesiveError> for RigidCohesiveHeatError<E> {
    fn from(error: RigidCohesiveError) -> Self {
        Self::Mechanical(error)
    }
}
impl<E> From<CohesiveHeatError<E>> for RigidCohesiveHeatError<E> {
    fn from(error: CohesiveHeatError<E>) -> Self {
        Self::Heat(error)
    }
}
impl<E> From<DepositionError> for RigidCohesiveHeatError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Heat(error.into())
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for RigidCohesiveHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native rigid cohesive heat: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for RigidCohesiveHeatError<E> {}

/// Observations from the same accepted native and thermal increment.
pub struct HeatedRigidCohesiveReport<P> {
    pub mechanical: RigidCohesiveFrame,
    pub heat: FrictionHeatNetworkReport,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub thermal_trials: usize,
    /// Original site indices in canonical ascending thermal-link order.
    pub site_indices: Vec<usize>,
    pub fracture_work_j: f64,
    /// Delta(native kinetic/gravity/bond + H) + fracture - load work - other heat.
    pub combined_energy_residual_j: f64,
    pub summed_energy_residual_j: f64,
}

/// Advance both original owners through one common native timestep.
///
/// All sites require distinct thermal routes. Signed heat into A/B comes from
/// solved common face temperatures and original k*A/d conductances, including
/// ordinary conduction and shared finite cells. No heat fraction or new capacity
/// is supplied. Source, native and caloric refusal leaves both accepted owners
/// unchanged. The callback must return staged source payloads without publishing.
///
/// Conductive geometry is explicitly a fixed material-reference attachment.
/// Normal jump must stay closed along the original straight constitutive path.
/// A moving contact footprint, changing coverage/area, finite bond resistance,
/// or Eulerian thermal geometry requires another model. Rigid rotations alone
/// do not transform this reference network into current world coordinates.
/// Mechanical coefficients remain frozen within the increment; explicit source
/// temperature continuation and timestep refinement own feedback accuracy.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn apply_rigid_cohesive_heat<'a, 't, P, E>(
    mechanical: &mut RigidCohesiveSystem<'a>,
    loads: [RigidCohesiveLoad; 2],
    thermal: &mut DepositionSimulation<'t>,
    bodies: [RigidCohesiveHeatBody<'_>; 2],
    contacts: &[RigidCohesiveHeatContact<'_>],
    controls: CohesiveHeatControls,
    energy_guess_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'t>, &[u64], &[f64]) -> Result<CohesiveHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedRigidCohesiveReport<P>, RigidCohesiveHeatError<E>> {
    let check = || {
        if gate.is_requested() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    check()?;
    if mechanical.time_s() != thermal.time_s()
        || contacts.is_empty()
        || contacts.len() != mechanical.sites().len()
        || contacts.len() > controls.maximum_interfaces
        || bodies[0].body_id == bodies[1].body_id
        || bodies
            .iter()
            .zip(mechanical.bodies())
            .any(|(a, b)| a.mechanical_material != b.material())
    {
        return Err(RigidCohesiveHeatError::Invalid(
            "complete synchronized original body/site assignment required",
        ));
    }
    let mut ordered: Vec<_> = contacts.iter().collect();
    ordered.sort_by_key(|c| c.route.link_id());
    let mut seen = BTreeSet::new();
    for contact in &ordered {
        check()?;
        let site =
            mechanical
                .sites()
                .get(contact.site_index)
                .ok_or(RigidCohesiveHeatError::Invalid(
                    "unknown original cohesive site",
                ))?;
        let original = site.point();
        let point = original.numerical();
        let area = point.area_m2();
        let other = contact.route.contact_area_m2();
        if !seen.insert(contact.site_index)
            || point.material().jump_m()[0] > 0.0
            || (area - other).abs() > 64.0 * f64::EPSILON * area.max(other)
            || !contact.temperature_guess_k.is_finite()
            || contact.temperature_guess_k <= 0.0
        {
            return Err(RigidCohesiveHeatError::Invalid(
                "unique closed reference patches with matching face areas required",
            ));
        }
        for (id, body) in contact.route.cells().into_iter().zip(&bodies) {
            let cell = thermal
                .cells()
                .get(&id)
                .ok_or(RigidCohesiveHeatError::Invalid(
                    "missing original thermal receiver",
                ))?;
            if cell.body_id != body.body_id || cell.material_id != body.thermal_material {
                return Err(RigidCohesiveHeatError::Invalid(
                    "thermal face changed native A/B body or material assignment",
                ));
            }
        }
    }
    let old_energy = mechanical.total_energy_j()?;
    let next_step = mechanical
        .steps()
        .checked_add(1)
        .ok_or(RigidCohesiveHeatError::Invalid("native clock overflow"))?;
    let target = next_step as f64 * mechanical.step_duration_s();
    let routes: Vec<_> = ordered.iter().map(|c| c.route).collect();
    let guesses: Vec<_> = ordered.iter().map(|c| c.temperature_guess_k).collect();
    let network =
        FrictionHeatNetworkStep::new(thermal, target, &routes, controls.maximum_interfaces, check)?;
    let (solved, report) = mechanical.step_with_acceptance(loads, Some(gate), |frame, sites| {
        if frame.time_s != target || frame.parameter_work_j != 0.0 {
            return Err(RigidCohesiveHeatError::Invalid(
                "frozen source native increment required",
            ));
        }
        let mut losses = Vec::with_capacity(ordered.len());
        let mut total_heat = 0.0;
        for contact in &ordered {
            check()?;
            let site =
                sites[contact.site_index]
                    .as_ref()
                    .ok_or(RigidCohesiveHeatError::Invalid(
                        "incomplete native site observation",
                    ))?;
            let friction = site.friction.map_or(0.0, |f| f.friction_dissipation_j);
            if site.jump_m[0] > 0.0
                || [
                    friction,
                    site.viscous_dissipation_j,
                    site.fracture_dissipation_j,
                ]
                .iter()
                .any(|v| !v.is_finite() || *v < 0.0)
            {
                return Err(RigidCohesiveHeatError::Invalid(
                    "closed nonnegative physical native site losses required",
                ));
            }
            let loss = friction + site.viscous_dissipation_j;
            total_heat += loss;
            losses.push(loss);
        }
        if !total_heat.is_finite()
            || (total_heat - frame.friction_dissipation_j - frame.viscous_dissipation_j).abs()
                > controls.combined_energy_tolerance_j
        {
            return Err(RigidCohesiveHeatError::Invalid(
                "site heat differs from original native loss",
            ));
        }
        let solved = solve_frozen_cohesive_heat(
            thermal,
            &network,
            &guesses,
            &losses,
            controls,
            energy_guess_j,
            resolve,
            gate,
        )?;
        let combined = frame.mechanical_energy_j + frame.cohesive_stored_energy_j - old_energy
            + solved.heat.thermal.stored_energy_change_j
            + frame.fracture_dissipation_j
            - frame.external_work_j
            - solved.heat.external_heat_j;
        let summed = frame.energy_residual_j + solved.heat.total_energy_residual_j;
        if [combined, summed]
            .iter()
            .any(|r| !r.is_finite() || r.abs() > controls.combined_energy_tolerance_j)
        {
            return Err(RigidCohesiveHeatError::Invalid(
                "original native mechanics and finite heat do not balance",
            ));
        }
        let report = (*frame, combined, summed);
        Ok((solved, report))
    })?;
    // Native publication has passed its final cancellation gate. This move is
    // infallible and publishes the corresponding accepted caloric endpoint.
    *thermal = solved.thermal;
    Ok(HeatedRigidCohesiveReport {
        mechanical: report.0,
        fracture_work_j: report.0.fracture_dissipation_j,
        heat: solved.heat,
        source: solved.source,
        nonlinear: solved.nonlinear,
        thermal_trials: solved.trials,
        site_indices: ordered.iter().map(|c| c.site_index).collect(),
        combined_energy_residual_j: report.1,
        summed_energy_residual_j: report.2,
    })
}
