//! Shared-body Hertz/rolling contact islands on the original finite heat network.
//!
//! Each original interface receives only its own surface-friction work. Bulk
//! rolling work from every patch on a sphere is accumulated into that sphere's
//! single original caloric inventory. Interface temperature and sphere bulk
//! temperature remain distinct constitutive coordinates in every trial.

use super::RigidHertzRollingHeatContact;
use crate::friction_heat::hertz::multiple::{
    RigidHertzIslandHeatReceipt, RigidHertzIslandHeatSource,
};
use crate::friction_heat::hertz::{RigidHertzHeatContact, validate_geometry};
use crate::friction_heat::network::{
    FrictionHeatNetworkInput, FrictionHeatNetworkReport, FrictionHeatNetworkStep,
    InterfaceHeatInput,
};
use crate::friction_heat::rolling::RigidRollingHeatBinding;
use crate::friction_heat::{FrictionHeatRoute, FrictionHeatTolerances, finite};
use fs_contact::interface_binding::thermal_friction::hertz::SourceHertzWrenchError;
use fs_exec::CancelGate;
use fs_math::det;
use fs_mbd::tangential_friction::hertz_rolling::RigidHertzRollingContact;
use fs_mbd::tangential_friction::hertz_rolling::multiple::{
    PreparedRigidHertzRollingIsland, RigidHertzRollingIslandCandidate,
    RigidHertzRollingIslandContact, RigidHertzRollingIslandControls, RigidHertzRollingIslandError,
    RigidHertzRollingIslandLaw,
};
use fs_mbd::tangential_friction::rolling::RigidRollingHistory;
use fs_mbd::{MassProperties, RigidBodyState};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use std::cell::Cell;
use std::collections::BTreeMap;

/// One original patch, surface heat route and independently owned rolling history.
/// Multiple patches may share the exact same sphere bulk binding. A different
/// binding for that same native body, or a shared cell for different bodies,
/// refuses rather than duplicating or conflating caloric inventory.
#[derive(Clone, Copy)]
pub struct RigidHertzRollingIslandHeatContact<'a> {
    pub patch_id: &'a str,
    pub contact: RigidHertzRollingHeatContact<'a>,
    pub history: &'a RigidRollingHistory,
    pub route: &'a FrictionHeatRoute,
    pub bulk_binding: RigidRollingHeatBinding,
    /// Positive initial guess only, never a prescribed physical temperature.
    pub temperature_guess_k: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct RigidHertzRollingIslandHeatControls {
    pub mechanical: RigidHertzRollingIslandControls,
    pub maximum_cells: usize,
    /// Original cell increments plus distinct interface temperatures.
    pub maximum_unknowns: usize,
    pub maximum_trials: usize,
    pub maximum_newton_iterations: usize,
    pub difference_step_j: f64,
    /// Numerical scaling [J/K], not an additional heat capacity.
    pub temperature_scale_j_per_k: f64,
    pub heat: FrictionHeatTolerances,
    pub combined_energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
}

#[derive(Debug)]
pub enum RigidHertzRollingIslandHeatError<E> {
    Invalid(&'static str),
    Mechanical(RigidHertzRollingIslandError),
    FrictionSource(SourceHertzWrenchError),
    Thermal(DepositionError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
}
impl<E> From<RigidHertzRollingIslandError> for RigidHertzRollingIslandHeatError<E> {
    fn from(error: RigidHertzRollingIslandError) -> Self {
        Self::Mechanical(error)
    }
}
impl<E> From<SourceHertzWrenchError> for RigidHertzRollingIslandHeatError<E> {
    fn from(error: SourceHertzWrenchError) -> Self {
        Self::FrictionSource(error)
    }
}
impl<E> From<DepositionError> for RigidHertzRollingIslandHeatError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E> From<NewtonError> for RigidHertzRollingIslandHeatError<E> {
    fn from(error: NewtonError) -> Self {
        Self::Newton(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for RigidHertzRollingIslandHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Hertz/rolling island finite heat: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for RigidHertzRollingIslandHeatError<E> {}

/// One whole-sphere inventory, regardless of the number of its contact patches.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzRollingIslandBulkHeat {
    pub binding: RigidRollingHeatBinding,
    pub temperature_before_k: f64,
    pub temperature_after_k: f64,
    pub dissipated_energy_j: f64,
}

/// Publish all original body, history and thermal endpoints together.
pub struct HeatedRigidHertzRollingIslandCandidate<'m, P> {
    pub mechanical: RigidHertzRollingIslandCandidate,
    pub thermal: DepositionSimulation<'m>,
    /// Original surface-network audit. Its source heat includes injected bulk
    /// rolling work; actual outside heat is retained separately below.
    /// No bulk-source entropy classification is claimed by this adapter.
    pub heat: FrictionHeatNetworkReport,
    /// Ascending original thermal link order; explicit IDs preserve native mapping.
    pub friction_sources: Vec<RigidHertzIslandHeatReceipt>,
    /// Ascending native body index, each physical sphere occurring once.
    pub bulk_heat: Vec<RigidHertzRollingIslandBulkHeat>,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub trials: usize,
    pub target_time_s: f64,
    /// Authored and callback outside heat measured before bulk injection.
    pub external_heat_j: f64,
    /// Actual caloric change minus outside heat and both physical loss channels.
    pub heat_transfer_residual_j: f64,
    /// Actual unique-body kinetic plus original caloric change minus outside heat.
    pub combined_energy_residual_j: f64,
    /// Independent sum of the original native and heat-network audit residuals.
    pub summed_energy_residual_j: f64,
}

struct Staged<P> {
    mechanical: RigidHertzRollingIslandCandidate,
    friction_sources: Vec<RigidHertzIslandHeatReceipt>,
    bulk_heat: Vec<RigidHertzRollingIslandBulkHeat>,
    external_heat_j: f64,
    source: P,
}

/// Solve original surface friction, material rolling and finite heat together.
///
/// One prepared full five-coordinate response per patch retains all shared-body
/// recoil. Each trial queries the pinned friction source at its own common face
/// temperature and explicit pressure, and admits each original rolling card at
/// the actual endpoint temperature of its original sphere. The callback supplies
/// the original endpoint conductivities and other heat only, in ascending link
/// order as in the existing Hertz-island heat source API.
///
/// Every patch needs a distinct original thermal link and rolling work/history
/// owner. Shared sphere bindings must agree exactly, carry that native sphere's
/// entire mass, and be its only active thermal cell. That original cell receives
/// the sum of its patches' bulk work once. No mean face temperature, duplicate
/// capacity, prescribed heat fraction or extra material loss is introduced.
///
/// The native normal/shear decoupling approximation, frozen normal load/geometry
/// and unresolved endpoint-event limits remain in force. Each sphere is one
/// uniform bulk temperature; the original orthogonal finite-cell/common-face
/// transport is backward Euler. No spatial thermal validity, variable viscosity,
/// partial-slip memory, entropy production or root uniqueness claim is added.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn propose_rigid_hertz_rolling_island_heat<'m, P, E>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[RigidHertzRollingIslandHeatContact<'_>],
    thermal: &DepositionSimulation<'m>,
    accepted_time_s: f64,
    dt_s: f64,
    controls: RigidHertzRollingIslandHeatControls,
    energy_guess_j: &[f64],
    resolve: impl Fn(
        &DepositionSimulation<'m>,
        &[u64],
        &[f64],
    ) -> Result<RigidHertzIslandHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedRigidHertzRollingIslandCandidate<'m, P>, RigidHertzRollingIslandHeatError<E>> {
    let checkpoint = || {
        if gate.is_requested() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    checkpoint()?;
    controls.newton.validate()?;
    controls.mechanical.validate()?;
    let count = contacts.len();
    let cells = thermal.cells().len();
    let dimension = cells
        .checked_add(count)
        .ok_or(RigidHertzRollingIslandHeatError::Invalid(
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
        || cells == 0
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
        || contacts
            .iter()
            .any(|c| !c.temperature_guess_k.is_finite() || c.temperature_guess_k <= 0.0)
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
        return Err(RigidHertzRollingIslandHeatError::Invalid(
            "synchronized bounded native, caloric and interface coordinates required",
        ));
    }
    let duration_s = target - accepted_time_s;
    let mut ordered: Vec<_> = contacts.iter().collect();
    ordered.sort_by_key(|c| c.route.link_id());
    if ordered
        .windows(2)
        .any(|p| p[0].route.link_id() == p[1].route.link_id())
    {
        return Err(RigidHertzRollingIslandHeatError::Invalid(
            "each original thermal link requires exactly one surface-work owner",
        ));
    }
    let mut initial_sources = Vec::with_capacity(count);
    let mut initial_bulk_temperatures = Vec::with_capacity(count);
    for c in &ordered {
        checkpoint()?;
        if c.contact.interface != c.contact.source.interface().interface() {
            return Err(RigidHertzRollingIslandHeatError::Invalid(
                "rolling and surface friction must retain the same original interface",
            ));
        }
        let radius = det::sqrt(c.contact.patch.contact_area_m2() / core::f64::consts::PI);
        initial_sources.push(c.contact.source.hertz_at(
            c.temperature_guess_k,
            c.contact.nominal_pressure_pa,
            c.contact.patch.normal_load_n(),
            radius,
        )?);
        initial_bulk_temperatures.push(thermal.temperature_k(c.bulk_binding.cell_id)?);
    }
    let native_contacts: Vec<_> = ordered
        .iter()
        .zip(&initial_sources)
        .zip(&initial_bulk_temperatures)
        .map(
            |((c, source), &temperature_k)| RigidHertzRollingIslandContact {
                patch_id: c.patch_id,
                contact: RigidHertzRollingContact {
                    binding: c.contact.binding,
                    frame: c.contact.frame,
                    hertz: source.law(),
                    patch: c.contact.patch,
                    interface: c.contact.interface,
                    rolling: c.contact.rolling,
                    temperature_k,
                    excitation_frequency_hz: c.contact.excitation_frequency_hz,
                    work: c.contact.work,
                },
                history: c.history,
            },
        )
        .collect();
    let prepared = PreparedRigidHertzRollingIsland::new(
        bodies,
        properties,
        &native_contacts,
        duration_s,
        controls.mechanical,
        gate,
    )?;
    let tolerance = controls
        .mechanical
        .mechanical
        .mechanical
        .geometry_tolerance_m;
    let mut mechanical_indices = Vec::with_capacity(count);
    let mut spheres = BTreeMap::new();
    let mut cell_owners = BTreeMap::new();
    let mut thermal_body_owners = BTreeMap::new();
    for c in &ordered {
        checkpoint()?;
        let index = prepared
            .contacts()
            .binary_search_by(|p| p.patch_id.as_str().cmp(c.patch_id))
            .map_err(|_| {
                RigidHertzRollingIslandHeatError::Invalid("missing native patch identity")
            })?;
        mechanical_indices.push(index);
        let patch = &prepared.contacts()[index];
        let sphere = patch.deforming_body;
        let binding = c.bulk_binding;
        let side = usize::from(sphere == c.contact.binding.body_b);
        let expected_cell = [c.route.patch.link.a, c.route.patch.link.b][side];
        let cell = thermal.cells().get(&binding.cell_id).ok_or(
            RigidHertzRollingIslandHeatError::Invalid("inactive sphere bulk cell"),
        )?;
        if binding.body_index != sphere
            || binding.cell_id != expected_cell
            || cell.body_id != binding.thermal_body_id
            || cell.material_id != binding.material_id
            || cell.mass_kg != properties[sphere].mass()
        {
            return Err(RigidHertzRollingIslandHeatError::Invalid(
                "whole-sphere mass and declared body/material must occupy the original route side",
            ));
        }
        if let Some((previous, _)) = spheres.get(&sphere) {
            if *previous != binding {
                return Err(RigidHertzRollingIslandHeatError::Invalid(
                    "all patches of one sphere must retain exactly one original bulk binding",
                ));
            }
        } else {
            if cell_owners.insert(binding.cell_id, sphere).is_some()
                || thermal_body_owners
                    .insert(binding.thermal_body_id, sphere)
                    .is_some()
            {
                return Err(RigidHertzRollingIslandHeatError::Invalid(
                    "different native spheres cannot share one original caloric inventory",
                ));
            }
            for (&id, other) in thermal.cells() {
                checkpoint()?;
                if id != binding.cell_id && other.body_id == binding.thermal_body_id {
                    return Err(RigidHertzRollingIslandHeatError::Invalid(
                        "whole-sphere thermal body must have exactly one active caloric cell",
                    ));
                }
            }
            spheres.insert(sphere, (binding, thermal.temperature_k(binding.cell_id)?));
        }
        let p = patch.contact_points_world[0];
        let point = [p.x, p.y, p.z];
        validate_geometry(
            point,
            RigidHertzHeatContact {
                binding: c.contact.binding,
                frame: c.contact.frame,
                source: c.contact.source,
                normal_force_n: patch.normal_force_n,
                contact_radius_m: patch.contact_radius_m,
                nominal_pressure_pa: c.contact.nominal_pressure_pa,
            },
            c.route,
            tolerance,
        )
        .map_err(RigidHertzRollingIslandHeatError::Invalid)?;
        let sign = if side == 0 { 1.0 } else { -1.0 };
        let normal = c.contact.frame.normal();
        let delta = std::array::from_fn::<_, 3, _>(|i| {
            0.5 * cell.geometry.start_m[i] + 0.5 * cell.geometry.end_m[i]
                - (point[i] + sign * patch.sphere_radius_m * normal[i])
        });
        if delta.iter().any(|v| !v.is_finite())
            || det::hypot(det::hypot(delta[0], delta[1]), delta[2]) > tolerance
        {
            return Err(RigidHertzRollingIslandHeatError::Invalid(
                "whole-sphere caloric centre must coincide with its original native COM",
            ));
        }
    }
    let sphere_bindings: Vec<_> = spheres.values().copied().collect();
    let sphere_indices: Vec<_> = ordered
        .iter()
        .map(|c| {
            sphere_bindings
                .binary_search_by_key(&c.bulk_binding.body_index, |(b, _)| b.body_index)
                .map_err(|_| {
                    RigidHertzRollingIslandHeatError::Invalid("missing original sphere binding")
                })
        })
        .collect::<Result<_, _>>()?;
    let routes: Vec<_> = ordered.iter().map(|c| c.route).collect();
    let step = FrictionHeatNetworkStep::new(
        thermal,
        target,
        &routes,
        controls.mechanical.maximum_contacts,
        checkpoint,
    )?;
    let trials = Cell::new(0usize);
    let evaluate = |x: &[f64]| {
        checkpoint()?;
        if trials.get() >= controls.maximum_trials {
            return Err(RigidHertzRollingIslandHeatError::TrialBudget);
        }
        trials.set(trials.get() + 1);
        if x.len() != dimension || x.iter().any(|v| !v.is_finite()) {
            return Err(RigidHertzRollingIslandHeatError::Invalid(
                "complete finite thermal iterate required",
            ));
        }
        let temperatures: Vec<_> = ordered
            .iter()
            .enumerate()
            .map(|(i, c)| c.temperature_guess_k + x[cells + i] / controls.temperature_scale_j_per_k)
            .collect();
        step.evaluate(
            &x[..cells],
            &temperatures,
            |endpoint, temperatures| {
                let mut receipts = Vec::with_capacity(count);
                for (c, &temperature) in ordered.iter().zip(temperatures) {
                    checkpoint()?;
                    let native = &prepared.contacts()[mechanical_indices[receipts.len()]];
                    receipts.push(RigidHertzIslandHeatReceipt {
                        patch_id: c.patch_id.to_owned(),
                        friction_source: c.contact.source.hertz_at(
                            temperature,
                            c.contact.nominal_pressure_pa,
                            native.normal_force_n,
                            native.contact_radius_m,
                        )?,
                    });
                }
                let mut bulk_heat = Vec::with_capacity(sphere_bindings.len());
                for &(binding, temperature_before_k) in &sphere_bindings {
                    checkpoint()?;
                    bulk_heat.push(RigidHertzRollingIslandBulkHeat {
                        binding,
                        temperature_before_k,
                        temperature_after_k: endpoint.temperature_k(binding.cell_id)?,
                        dissipated_energy_j: 0.0,
                    });
                }
                let laws: Vec<_> = receipts
                    .iter()
                    .zip(&sphere_indices)
                    .map(|(r, &i)| RigidHertzRollingIslandLaw {
                        patch_id: &r.patch_id,
                        law: r.friction_source.law(),
                        temperature_k: bulk_heat[i].temperature_after_k,
                    })
                    .collect();
                let mechanical = prepared.advance(&laws, gate)?;
                let RigidHertzIslandHeatSource {
                    interfaces,
                    transport,
                } = resolve(endpoint, step.link_ids(), temperatures)
                    .map_err(RigidHertzRollingIslandHeatError::Source)?;
                if interfaces.len() != count
                    || interfaces
                        .iter()
                        .zip(step.link_ids())
                        .any(|(source, id)| source.link_id != *id)
                {
                    return Err(RigidHertzRollingIslandHeatError::Invalid(
                        "conductivities must retain every canonical original link identity",
                    ));
                }
                let mut inputs = Vec::with_capacity(count);
                for (((c, &i), &sphere_index), source) in ordered
                    .iter()
                    .zip(&mechanical_indices)
                    .zip(&sphere_indices)
                    .zip(interfaces)
                {
                    checkpoint()?;
                    let native = &mechanical.report.contacts[i];
                    if native.patch_id != c.patch_id
                        || native.deforming_body != c.bulk_binding.body_index
                    {
                        return Err(RigidHertzRollingIslandHeatError::Invalid(
                            "changed original surface or bulk work owner",
                        ));
                    }
                    inputs.push(InterfaceHeatInput {
                        dissipated_j: native.surface_heat_j,
                        conductivity_w_m_k: source.conductivity_w_m_k,
                    });
                    let original = bulk_heat[sphere_index].dissipated_energy_j;
                    let total = finite(original + native.bulk_heat_j)?;
                    if native.bulk_heat_j > 0.0 && total == original {
                        return Err(DepositionError::Arithmetic.into());
                    }
                    bulk_heat[sphere_index].dissipated_energy_j = total;
                }
                let StepCoefficients {
                    mut sources_w,
                    conductances_w_k,
                    payload,
                } = transport;
                checkpoint()?;
                // Outside heat is measured before adding internal mechanical loss.
                // Subtracting the bulk contribution from a larger sum can erase a
                // small signed outside supply and changes the physical boundary.
                let external_heat_j = finite(
                    duration_s
                        * endpoint.external_step_power_w(
                            &sources_w,
                            &conductances_w_k,
                            checkpoint,
                        )?,
                )?;
                for bulk in &bulk_heat {
                    checkpoint()?;
                    let power = finite(bulk.dissipated_energy_j / duration_s)?;
                    if bulk.dissipated_energy_j > 0.0 && power == 0.0 {
                        return Err(DepositionError::Arithmetic.into());
                    }
                    if power != 0.0 {
                        let outside = sources_w.entry(bulk.binding.cell_id).or_default();
                        let sum = finite(*outside + power)?;
                        if sum == *outside {
                            return Err(DepositionError::Arithmetic.into());
                        }
                        *outside = sum;
                    }
                }
                Ok(FrictionHeatNetworkInput {
                    interfaces: inputs,
                    transport: StepCoefficients {
                        sources_w,
                        conductances_w_k,
                        payload: Staged {
                            mechanical,
                            friction_sources: receipts,
                            bulk_heat,
                            external_heat_j,
                            source: payload,
                        },
                    },
                })
            },
            checkpoint,
        )
    };
    let problem =
        FallibleFiniteDifference::new(dimension, controls.difference_step_j, |x: &[f64]| {
            evaluate(x).map(|trial| trial.residual_j().to_vec())
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
        return Err(RigidHertzRollingIslandHeatError::Interrupted(
            result.progress,
        ));
    }
    if !result.report.converged {
        return Err(RigidHertzRollingIslandHeatError::NotConverged(
            result.report,
        ));
    }
    let trial = evaluate(&solver.x)?;
    let staged = trial.payload();
    let native = &staged.mechanical.report;
    let caloric_change = trial.report().thermal.stored_energy_change_j;
    let transfer = finite(
        caloric_change - staged.external_heat_j - native.surface_heat_j - native.bulk_heat_j,
    )?;
    let combined = finite(
        native.kinetic_after_j - native.kinetic_before_j + caloric_change - staged.external_heat_j,
    )?;
    let summed = finite(native.total_energy_residual_j + trial.report().total_energy_residual_j)?;
    if transfer.abs() > controls.heat.total_j
        || combined.abs() > controls.combined_energy_tolerance_j
        || summed.abs() > controls.combined_energy_tolerance_j
    {
        return Err(RigidHertzRollingIslandHeatError::Invalid(
            "original island surface/bulk work and finite caloric inventory do not balance",
        ));
    }
    let (thermal, staged, heat) = trial.accept(controls.heat, checkpoint)?;
    checkpoint()?;
    Ok(HeatedRigidHertzRollingIslandCandidate {
        mechanical: staged.mechanical,
        thermal,
        heat,
        friction_sources: staged.friction_sources,
        bulk_heat: staged.bulk_heat,
        source: staged.source,
        nonlinear: result.report,
        trials: trials.get(),
        target_time_s: target,
        external_heat_j: staged.external_heat_j,
        heat_transfer_residual_j: transfer,
        combined_energy_residual_j: combined,
        summed_energy_residual_j: summed,
    })
}

fn numerical<E>(
    error: FiniteDifferenceError<RigidHertzRollingIslandHeatError<E>>,
) -> RigidHertzRollingIslandHeatError<E> {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(field) => RigidHertzRollingIslandHeatError::Invalid(field),
    }
}
