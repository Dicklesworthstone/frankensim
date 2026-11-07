//! Sustained native finite-solid friction and the original finite enthalpy network.
//!
//! Each endpoint thermal trial restarts the SAME force/contact/drift trajectory
//! and original source histories. Only its unique, loaded, nonopening midpoint
//! tangential impulse work supplies heat. Normal/contact remainder, numerical
//! integration defects and declared exported work remain separate accounts.
//! Original body/common/declared-mean temperatures feed the source; native speed,
//! distance, material directors, simultaneous impulses and path refinement retain
//! their existing owners. One balanced candidate returns all accepted states.
//!
//! This is endpoint temperature feedback on fixed reference thermal patches and
//! frozen bulk hyperelastic coefficients. It adds no bulk expansion, moving heat
//! footprint, pressure closure or impact heating. Pressure-dependent friction
//! coordinates refuse: an instantaneous normal impulse is not a pressure source.
//! Refinement of the existing first-order trajectory and thermal interval is
//! required; nonlinear convergence alone does not certify time accuracy.

mod account;

use super::directional::{DirectionalHeatError, DirectionalHeatMaterial};
use super::network::{
    FrictionHeatNetworkInput, FrictionHeatNetworkReport, FrictionHeatNetworkStep,
    InterfaceHeatInput,
};
use super::{FrictionHeatRoute, FrictionHeatTolerances};
use crate::finite_solid::contact::scene::{
    FiniteSolidSphereScene, FiniteSolidSphereSceneError, FiniteSolidSphereSceneKey,
};
use crate::finite_solid::contact::step::FiniteSolidSphereTrajectoryControls;
use crate::finite_solid::contact::step::sustained::source::*;
use crate::finite_solid::contact::surface::FiniteSolidContactSurface;
use crate::finite_solid::{FiniteSolidBody, FiniteSolidError};
use fs_contact::interface_binding::thermal_friction::elliptic::SourceFrictionConditions;
use fs_exec::Cx;
use fs_matdb::QueryPoint;
use fs_mbd::{DynamicsError, Gravity, MassProperties, RigidBodyState, Wrench};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::{StepCoefficients, resolvable_increment_j};
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use std::cell::{Cell, RefCell};
use std::convert::Infallible;

pub use super::hertz::multiple::{
    RigidHertzIslandHeatConductivity as FiniteSolidHeatConductivity,
    RigidHertzIslandHeatSource as FiniteSolidHeatSource,
};

/// Borrowed original mechanics and held loads; no second mechanical state owner.
pub struct FiniteSolidHeatMotion<'s, 'a> {
    pub solid: &'s FiniteSolidBody<'a>,
    pub surface: &'s FiniteSolidContactSurface,
    pub scene: &'s FiniteSolidSphereScene,
    pub bodies: &'s [RigidBodyState],
    pub properties: &'s [MassProperties],
    pub nodal_forces_n: &'s [(usize, f64)],
    pub body_force_n_m3: [f64; 3],
    pub wrenches: &'s [Wrench],
    pub gravity: Gravity,
    pub origin_m: [f64; 3],
}

/// One original scene feature, source law and ordered thermal face.
/// Route A belongs to the first rigid witness; B belongs to its counterpart.
#[derive(Clone, Copy, Debug)]
pub struct FiniteSolidHeatContact<'a> {
    pub contact_id: u64,
    pub key: FiniteSolidSphereSceneKey,
    pub binding: FiniteSolidSourceContactBinding<'a>,
    pub route: &'a FrictionHeatRoute,
    pub material: DirectionalHeatMaterial<'a>,
    /// Declared physical conversion; the remaining tangential work is exported.
    pub heat_fraction: f64,
    pub temperature_guess_k: f64,
    /// Caller-declared circular neighborhood INSIDE the original thermal face.
    /// Its area must fit the face; its centre/normal come from that route.
    /// Witness migration outside this fixed neighborhood refuses, never remaps.
    pub patch_radius_m: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct FiniteSolidHeatControls {
    pub mechanical: FiniteSolidSphereTrajectoryControls,
    pub history: FiniteSolidSourceTrajectoryControls,
    pub maximum_cells: usize,
    pub maximum_interfaces: usize,
    pub maximum_unknowns: usize,
    pub maximum_trials: usize,
    pub maximum_newton_iterations: usize,
    pub difference_step_j: f64,
    /// Numerical coordinate scale only; not an interface heat capacity.
    pub temperature_scale_j_per_k: f64,
    pub geometry_tolerance_m: f64,
    /// Norm of the difference between unit contact and original face normals.
    pub normal_tolerance: f64,
    pub heat: FrictionHeatTolerances,
    pub combined_energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
    pub audit_entropy: bool,
}

#[derive(Debug)]
pub enum FiniteSolidHeatError<E> {
    Invalid(&'static str),
    Mechanical(FiniteSolidSustainedSourceTrajectoryError<DepositionError>),
    Geometry(FiniteSolidSphereSceneError),
    Solid(FiniteSolidError),
    Rigid(DynamicsError),
    Material(DirectionalHeatError<Infallible>),
    Thermal(DepositionError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
}
impl<E> From<DepositionError> for FiniteSolidHeatError<E> {
    fn from(e: DepositionError) -> Self {
        Self::Thermal(e)
    }
}
impl<E> From<NewtonError> for FiniteSolidHeatError<E> {
    fn from(e: NewtonError) -> Self {
        Self::Newton(e)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for FiniteSolidHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-solid sustained friction/heat: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for FiniteSolidHeatError<E> {}

pub struct HeatedFiniteSolidCandidate<'a, 't, P> {
    pub solid: FiniteSolidBody<'a>,
    pub bodies: Vec<RigidBodyState>,
    pub history: FiniteSolidSourceTrajectoryHistory<'a>,
    pub mechanical: FiniteSolidSustainedSourceTrajectoryReport,
    pub thermal: DepositionSimulation<'t>,
    pub heat: FrictionHeatNetworkReport,
    pub exported_tangential_work_j: f64,
    /// Signed original contact loss minus ALL eligible tangential work.
    /// Contains the remaining normal/impact channel; never deposited as heat.
    pub unconverted_contact_energy_j: f64,
    /// Delta(K+U+H)-W_source-Q_outside+exported+unconverted-integration_defect.
    pub combined_energy_residual_j: f64,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub thermal_trials: usize,
    pub mechanical_solves: usize,
}

#[derive(Clone)]
struct Mechanics<'a> {
    solid: FiniteSolidBody<'a>,
    bodies: Vec<RigidBodyState>,
    history: FiniteSolidSourceTrajectoryHistory<'a>,
    report: FiniteSolidSustainedSourceTrajectoryReport,
    heat_j: Vec<f64>,
    exported_j: f64,
    unconverted_j: f64,
}
struct Staged<'a, P> {
    mechanical: Mechanics<'a>,
    source: P,
}

/// Close native sustained motion and finite heat on one accepted physical clock.
///
/// All scene contacts must retain their original feature, thermal neighborhood
/// and source binding at every native source evaluation and final endpoint.
/// These patch checks are not a continuous swept thermal-footprint certificate.
/// Unloaded,
/// opening, indeterminate or nonpassive events refuse this fixed-contact lane.
/// The two caloric objects and ordered material cards are checked by the existing
/// thermal interface owner. No source temperature may silently stay frozen.
/// The callback supplies only original endpoint conductivities and OTHER heat,
/// in ascending link order, and must stage its payload until candidate acceptance.
/// No accepted body, source history or thermal state is mutated on refusal.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn propose_finite_solid_contact_heat<'a, 't, P, E>(
    motion: &FiniteSolidHeatMotion<'_, 'a>,
    history: &FiniteSolidSourceTrajectoryHistory<'a>,
    thermal: &DepositionSimulation<'t>,
    duration_s: f64,
    contacts: &[FiniteSolidHeatContact<'_>],
    controls: FiniteSolidHeatControls,
    energy_guess_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'t>, &[u64], &[f64]) -> Result<FiniteSolidHeatSource<P>, E>,
    cx: &Cx<'_>,
) -> Result<HeatedFiniteSolidCandidate<'a, 't, P>, FiniteSolidHeatError<E>> {
    let check = || cx.checkpoint().map_err(|_| DepositionError::Cancelled);
    check()?;
    controls.newton.validate()?;
    let cells = thermal.cells().len();
    let count = contacts.len();
    let dimension = cells
        .checked_add(count)
        .ok_or(FiniteSolidHeatError::Invalid("thermal dimension overflow"))?;
    let target = thermal.time_s() + duration_s;
    if motion.solid.native_state().t != thermal.time_s()
        || history.time_s() != thermal.time_s()
        || !duration_s.is_finite()
        || duration_s <= 0.0
        || !target.is_finite()
        || target <= thermal.time_s()
        || cells == 0
        || cells > controls.maximum_cells
        || count == 0
        || count > controls.maximum_interfaces
        || dimension > controls.maximum_unknowns
        || controls.maximum_trials == 0
        || controls.maximum_newton_iterations == 0
        || energy_guess_j.len() != cells
        || energy_guess_j.iter().any(|x| !x.is_finite())
        || [
            controls.difference_step_j,
            controls.temperature_scale_j_per_k,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        || [
            controls.geometry_tolerance_m,
            controls.normal_tolerance,
            controls.combined_energy_tolerance_j,
            controls.heat.cell_j,
            controls.heat.thermal_j,
            controls.heat.interface_j,
            controls.heat.total_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
        || controls.normal_tolerance >= 1.0
    {
        return Err(FiniteSolidHeatError::Invalid(
            "bounded synchronous native/thermal inputs required",
        ));
    }
    let mut ordered: Vec<_> = contacts.iter().collect();
    ordered.sort_by_key(|c| c.route.link_id());
    for (index, contact) in ordered.iter().enumerate() {
        check()?;
        if ordered[..index]
            .iter()
            .any(|c| c.contact_id == contact.contact_id || c.key == contact.key)
            || !contact.heat_fraction.is_finite()
            || !(0.0..=1.0).contains(&contact.heat_fraction)
            || !contact.temperature_guess_k.is_finite()
            || contact.temperature_guess_k <= 0.0
            || !contact.patch_radius_m.is_finite()
            || contact.patch_radius_m < 0.0
            || core::f64::consts::PI * contact.patch_radius_m.powi(2)
                > contact.route.contact_area_m2()
        {
            return Err(FiniteSolidHeatError::Invalid(
                "unique fixed source patches and physical conversion fractions required",
            ));
        }
        contact
            .material
            .validate(contact.binding.source.law)
            .map_err(FiniteSolidHeatError::Material)?;
        account::source(contact)?;
    }
    account::scene(motion, motion.solid, motion.bodies, &ordered, controls, cx)?;
    let routes: Vec<_> = ordered.iter().map(|c| c.route).collect();
    let step =
        FrictionHeatNetworkStep::new(thermal, target, &routes, controls.maximum_interfaces, check)?;
    let trials = Cell::new(0usize);
    let mechanical_solves = Cell::new(0usize);
    // One bounded native result. Complete sampled queries include every body,
    // interface and mean T; transport receipts refresh even on a cache hit.
    let cache = RefCell::new(None::<(Vec<QueryPoint>, Mechanics<'a>)>);
    let evaluate = |x: &[f64], entropy: bool| {
        check()?;
        if trials.get() >= controls.maximum_trials {
            return Err(FiniteSolidHeatError::TrialBudget);
        }
        trials.set(trials.get() + 1);
        if x.len() != dimension || x.iter().any(|x| !x.is_finite()) {
            return Err(FiniteSolidHeatError::Invalid(
                "complete finite thermal iterate required",
            ));
        }
        let temperatures: Vec<_> = ordered
            .iter()
            .enumerate()
            .map(|(i, c)| c.temperature_guess_k + x[cells + i] / controls.temperature_scale_j_per_k)
            .collect();
        let increments: Vec<_> = thermal
            .cells()
            .values()
            .zip(&x[..cells])
            .map(|(c, &q)| resolvable_increment_j(c.specific_h_j_kg, q, c.mass_kg))
            .collect();
        let endpoint = |state: &DepositionSimulation<'t>, temperatures: &[f64]| {
            let queries: Vec<_> = ordered
                .iter()
                .zip(temperatures)
                .map(|(c, &t)| {
                    check()?;
                    c.material
                        .sample(state, c.route, c.binding.source.law, t, 0.0)
                        .map(|(_, q)| q)
                        .map_err(FiniteSolidHeatError::Material)
                })
                .collect::<Result<_, _>>()?;
            let cached = cache
                .borrow()
                .as_ref()
                .filter(|(q, _)| *q == queries)
                .map(|(_, m)| m.clone());
            let mechanical = if let Some(m) = cached {
                m
            } else {
                let find = |id| {
                    ordered.iter().position(|c| c.contact_id == id).ok_or(
                        DepositionError::InvalidInput("unmapped native contact cannot supply heat"),
                    )
                };
                mechanical_solves.set(mechanical_solves.get() + 1);
                let (solid, bodies, history, report) = motion
                    .solid
                    .sustained_source_elliptic_sphere_trajectory_candidate(
                        motion.surface,
                        motion.scene,
                        motion.bodies,
                        motion.properties,
                        motion.nodal_forces_n,
                        motion.body_force_n_m3,
                        motion.wrenches,
                        motion.gravity,
                        step.duration_s(),
                        motion.origin_m,
                        controls.mechanical,
                        history,
                        controls.history,
                        |key, c| {
                            let selected = ordered[find(c.id)?];
                            if selected.key != *key {
                                return Err(DepositionError::InvalidInput(
                                    "native contact changed thermal feature",
                                ));
                            }
                            Ok(selected.binding)
                        },
                        |context| {
                            check()?;
                            let i = find(context.kinematics.contact_id)?;
                            account::patch(
                                ordered[i],
                                context.kinematics.axes_world[0],
                                [
                                    context.kinematics.first_point_world,
                                    context.kinematics.second_point_world,
                                ],
                                controls,
                            )?;
                            let primary = &ordered[i].binding.source.law.principal_sources()[0]
                                .axes()
                                .temperature;
                            let t = *queries[i].axes().get(primary).ok_or(
                                DepositionError::InvalidInput("primary source temperature missing"),
                            )?;
                            Ok(SourceFrictionConditions {
                                query: queries[i].clone(),
                                contact_temperature_k: Some(t),
                                nominal_pressure_pa: None,
                            })
                        },
                        cx,
                    )
                    .map_err(FiniteSolidHeatError::Mechanical)?;
                if solid.native_state().t != target || history.time_s() != target {
                    return Err(FiniteSolidHeatError::Invalid(
                        "native source and finite heat endpoints must share the exact clock",
                    ));
                }
                account::scene(motion, &solid, &bodies, &ordered, controls, cx)?;
                let (heat_j, exported_j, unconverted_j) = account::work(&report, &ordered, cx)?;
                let m = Mechanics {
                    solid,
                    bodies,
                    history,
                    report,
                    heat_j,
                    exported_j,
                    unconverted_j,
                };
                *cache.borrow_mut() = Some((queries, m.clone()));
                m
            };
            let source = resolve(state, step.link_ids(), temperatures)
                .map_err(FiniteSolidHeatError::Source)?;
            if source.interfaces.len() != count
                || source
                    .interfaces
                    .iter()
                    .zip(step.link_ids())
                    .any(|(c, id)| c.link_id != *id)
            {
                return Err(FiniteSolidHeatError::Invalid(
                    "original conductivities must retain complete canonical link order",
                ));
            }
            let interfaces = source
                .interfaces
                .iter()
                .zip(&mechanical.heat_j)
                .map(|(c, &heat)| InterfaceHeatInput {
                    dissipated_j: heat,
                    conductivity_w_m_k: c.conductivity_w_m_k,
                })
                .collect();
            Ok(FrictionHeatNetworkInput {
                interfaces,
                transport: StepCoefficients {
                    sources_w: source.transport.sources_w,
                    conductances_w_k: source.transport.conductances_w_k,
                    payload: Staged {
                        mechanical,
                        source: source.transport.payload,
                    },
                },
            })
        };
        if entropy {
            step.evaluate_with_entropy(&increments, &temperatures, endpoint, check)
        } else {
            step.evaluate(&increments, &temperatures, endpoint, check)
        }
    };
    let numerical = |e| match e {
        FiniteDifferenceError::Source(e) => e,
        FiniteDifferenceError::Invalid(s) => FiniteSolidHeatError::Invalid(s),
    };
    let problem =
        FallibleFiniteDifference::new(dimension, controls.difference_step_j, |x: &[f64]| {
            evaluate(x, false).map(|trial| trial.residual_j().to_vec())
        })
        .map_err(numerical)?;
    let mut initial = energy_guess_j.to_vec();
    initial.resize(dimension, 0.0);
    let construction = NewtonKrylovState::new(&problem, initial, controls.newton);
    if let Some(e) = problem.take_failure() {
        return Err(numerical(e));
    }
    let mut solver = construction?;
    let result = solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
        cx.checkpoint().is_err() || problem.has_failed()
    });
    if let Some(e) = problem.take_failure() {
        return Err(numerical(e));
    }
    check()?;
    if result.progress != SolverRunProgress::Complete {
        return Err(FiniteSolidHeatError::Interrupted(result.progress));
    }
    if !result.report.converged {
        return Err(FiniteSolidHeatError::NotConverged(result.report));
    }
    let trial = evaluate(&solver.x, controls.audit_entropy)?;
    let staged = &trial.payload().mechanical;
    let m = &staged.report;
    let heat = trial.report();
    let combined = m.kinetic_change_j + m.strain_change_j + heat.thermal.stored_energy_change_j
        - m.source_work_j
        - heat.external_heat_j
        + staged.exported_j
        + staged.unconverted_j
        - m.integration_defect_j;
    if !combined.is_finite() || combined.abs() > controls.combined_energy_tolerance_j {
        return Err(FiniteSolidHeatError::Invalid(
            "native contact and finite heat accounts do not balance",
        ));
    }
    let (thermal, staged, heat) = trial.accept(controls.heat, check)?;
    check()?;
    let m = staged.mechanical;
    Ok(HeatedFiniteSolidCandidate {
        solid: m.solid,
        bodies: m.bodies,
        history: m.history,
        mechanical: m.report,
        thermal,
        heat,
        exported_tangential_work_j: m.exported_j,
        unconverted_contact_energy_j: m.unconverted_j,
        combined_energy_residual_j: combined,
        source: staged.source,
        nonlinear: result.report,
        thermal_trials: trials.get(),
        mechanical_solves: mechanical_solves.get(),
    })
}
