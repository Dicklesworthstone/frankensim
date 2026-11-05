//! Native quasistatic solid contact and the original finite enthalpy network.
//!
//! Every thermal trial starts from the same accepted displacement and sliding
//! history. Original interface claims see the solved common face temperature;
//! native contact resolves the corresponding normal load, slip and friction.
//! Only a complete mechanical/material/thermal candidate leaves this owner.
//!
//! This is a fixed-footprint plane-strain reduction with an explicit physical
//! depth and declared conversion of friction loss to heat. The existing thermal
//! owner solves the two-body heat partition; no side fraction or interface heat
//! capacity is invented. The remaining work leaves this thermal control volume.
//! Mechanical Huber creep and endpoint work/path quadrature remain unchanged.
//! Contact activation, moving thermal footprints, thermal strain/softening of
//! the elastic bulk, finite contact resistance and flash temperatures are outside
//! this lane. Its port balance is not a whole-solid elastic energy integral.

use super::network::{
    FrictionHeatNetworkInput, FrictionHeatNetworkReport, FrictionHeatNetworkStep,
    InterfaceHeatInput,
};
use super::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_contact::interface_binding::thermal_friction::FrictionQueryAxisKind;
use fs_exec::CancelGate;
use fs_qty::semantic::QuantityKind;
use fs_qty::{QuantitySpec, Temperature};
use fs_solid::SolidError;
use fs_solid::contact::{ContactFrictionBinding, ContactSolution, ContactTrajectory};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::resolvable_increment_j;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use std::cell::{Cell, RefCell};
use std::collections::BTreeSet;

mod work;

/// Original link-keyed conductivities and other thermal sources.
pub use super::hertz::multiple::{
    RigidHertzIslandHeatConductivity as SolidContactHeatConductivity,
    RigidHertzIslandHeatSource as SolidContactHeatSource,
};

/// One native solid node and its explicitly corresponding original thermal face.
/// The caller owns face coincidence and ordered source/body material assignment.
#[derive(Clone)]
pub struct SolidContactHeatContact<'a> {
    pub node: usize,
    pub friction: ContactFrictionBinding<'a>,
    pub route: &'a FrictionHeatRoute,
    /// Fraction converted to interface heat; the rest is reported outside heat.
    /// This does not prescribe the heat split between the two finite bodies.
    pub heat_fraction: f64,
    /// Positive numerical seed only; the common face temperature is solved.
    pub temperature_guess_k: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct SolidContactHeatControls {
    pub maximum_cells: usize,
    pub maximum_interfaces: usize,
    pub maximum_unknowns: usize,
    pub maximum_trials: usize,
    pub maximum_newton_iterations: usize,
    pub difference_step_j: f64,
    /// Numerical coordinate scaling only, never physical interface capacity.
    pub temperature_scale_j_per_k: f64,
    pub heat: FrictionHeatTolerances,
    pub port_energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
    pub audit_entropy: bool,
}

#[derive(Debug)]
pub enum SolidContactHeatError<E> {
    Invalid(&'static str),
    Mechanical(SolidError),
    Thermal(DepositionError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
}
impl<E> From<DepositionError> for SolidContactHeatError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E> From<NewtonError> for SolidContactHeatError<E> {
    fn from(error: NewtonError) -> Self {
        Self::Newton(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for SolidContactHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native solid friction/heat refused: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for SolidContactHeatError<E> {}

/// Complete new owners; the borrowed accepted trajectory and thermal state never
/// mutate during this solve, including on a late source or balance refusal.
pub struct HeatedSolidContactCandidate<'p, 'a, 't, P> {
    pub trajectory: ContactTrajectory<'p, 'a>,
    pub mechanical: ContactSolution,
    pub thermal: DepositionSimulation<'t>,
    pub heat: FrictionHeatNetworkReport,
    /// Same canonical link order as `heat.interfaces`.
    pub nodes: Vec<usize>,
    pub exported_nonthermal_work_j: f64,
    pub source: P,
    pub nonlinear: NewtonReport,
    /// Thermal residual evaluations, including repeated mechanical temperatures.
    pub trials: usize,
    /// Actual native equilibrium solves; exact repeated temperatures reuse one
    /// bounded candidate, while conductivity and thermal sources still refresh.
    pub mechanical_solves: usize,
    /// W_body_friction + delta(H) + D_other - W_support - Q_external.
    /// This closes the friction/thermal ports, not bulk elastic/load work.
    pub port_energy_residual_j: f64,
}

struct Staged<'p, 'a, P> {
    trajectory: ContactTrajectory<'p, 'a>,
    mechanical: ContactSolution,
    exported_work_j: f64,
    source: P,
}

struct CachedMechanical<'p, 'a> {
    temperatures_k: Vec<f64>,
    trajectory: ContactTrajectory<'p, 'a>,
    solution: ContactSolution,
}

/// Solve temperature-dependent native contact and finite heat simultaneously.
///
/// Every active node needs exactly one route; a mapped node losing contact also
/// refuses. Each thermal face is the declared nodal tributary patch. Where the
/// source has a nominal-pressure axis, width * depth must equal that face area.
/// The primary source temperature is the solved common face temperature.
/// Additional temperature axes refuse until their body/side ownership is supplied.
/// Pressure, speed and accumulated distance retain their native mechanical owner.
///
/// The callback supplies original endpoint conductivities and OTHER heat only,
/// in ascending original link-ID order. It must stage external history until the
/// complete candidate is accepted. Original mechanical iteration bounds remain
/// in force; cancellation is checked between complete mechanical trials. Refine
/// time for work/path accuracy and resolve changes to the fixed footprint first.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn propose_solid_contact_heat<'p, 'a, 't, P, E>(
    accepted: &ContactTrajectory<'p, 'a>,
    thermal: &DepositionSimulation<'t>,
    load: f64,
    surface_displacement: &[[f64; 2]],
    dt_s: f64,
    depth_m: f64,
    contacts: &[SolidContactHeatContact<'_>],
    controls: SolidContactHeatControls,
    energy_guess_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'t>, &[u64], &[f64]) -> Result<SolidContactHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedSolidContactCandidate<'p, 'a, 't, P>, SolidContactHeatError<E>> {
    let check = || {
        if gate.is_requested() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    check()?;
    controls.newton.validate()?;
    let nodes = accepted.displacement().len() / 2;
    let cells = thermal.cells().len();
    let count = contacts.len();
    let dimension = cells
        .checked_add(count)
        .ok_or(SolidContactHeatError::Invalid("thermal dimension overflow"))?;
    let target_s = accepted.time_s() + dt_s;
    if accepted.time_s() != thermal.time_s()
        || !dt_s.is_finite()
        || dt_s <= 0.0
        || !target_s.is_finite()
        || target_s <= accepted.time_s()
        || !depth_m.is_finite()
        || depth_m <= 0.0
        || cells == 0
        || cells > controls.maximum_cells
        || count == 0
        || count > controls.maximum_interfaces
        || dimension > controls.maximum_unknowns
        || controls.maximum_trials == 0
        || controls.maximum_newton_iterations == 0
        || energy_guess_j.len() != cells
        || energy_guess_j.iter().any(|q| !q.is_finite())
        || !controls.difference_step_j.is_finite()
        || controls.difference_step_j <= 0.0
        || !controls.temperature_scale_j_per_k.is_finite()
        || controls.temperature_scale_j_per_k <= 0.0
        || [
            controls.heat.cell_j,
            controls.heat.thermal_j,
            controls.heat.interface_j,
            controls.heat.total_j,
            controls.port_energy_tolerance_j,
        ]
        .iter()
        .any(|t| !t.is_finite() || *t < 0.0)
    {
        return Err(SolidContactHeatError::Invalid(
            "bounded synchronous solid/thermal inputs required",
        ));
    }
    let mut ordered: Vec<_> = contacts.iter().collect();
    ordered.sort_by_key(|c| c.route.link_id());
    let mut used = BTreeSet::new();
    let mut bindings = vec![None; nodes];
    for contact in &ordered {
        check()?;
        if contact.node >= nodes
            || !used.insert(contact.node)
            || !contact.temperature_guess_k.is_finite()
            || contact.temperature_guess_k <= 0.0
        {
            return Err(SolidContactHeatError::Invalid(
                "unique native nodes and positive temperature guesses required",
            ));
        }
        if !accepted
            .is_contact_active(contact.node)
            .map_err(SolidContactHeatError::Mechanical)?
        {
            return Err(SolidContactHeatError::Invalid(
                "fixed heat footprint is not initially in contact",
            ));
        }
        work::partition(0.0, depth_m, contact.heat_fraction)
            .map_err(SolidContactHeatError::Invalid)?;
        // A source may admit semantic temperatures as fixed external conditions
        // without declaring a query axis. They still need a thermal owner here.
        let primary_temperature = &contact.friction.law.axes().temperature;
        let conditions = &contact.friction.conditions;
        let extra_temperature = conditions.axes().keys().any(|axis| {
            if axis == primary_temperature {
                return false;
            }
            let quantity = conditions.axis_quantities().get(axis).copied();
            quantity
                .and_then(QuantitySpec::semantic_type)
                .is_some_and(|q| q.kind() == QuantityKind::AbsoluteTemperature)
                || (axis == "T"
                    && (quantity.is_none()
                        || quantity == Some(QuantitySpec::dimensional(Temperature::DIMS))))
        });
        if contact
            .friction
            .law
            .query_axes()
            .iter()
            .any(|axis| axis.kind == FrictionQueryAxisKind::Temperature)
            || extra_temperature
        {
            return Err(SolidContactHeatError::Invalid(
                "additional source temperatures require explicit thermal ownership",
            ));
        }
        if let Some(width) = contact.friction.nominal_contact_width_m {
            let area = width * depth_m;
            let face = contact.route.contact_area_m2();
            if !area.is_finite()
                || area <= 0.0
                || (area - face).abs() > 64.0 * f64::EPSILON * area.max(face)
            {
                return Err(SolidContactHeatError::Invalid(
                    "source pressure width and physical heat face disagree",
                ));
            }
        }
        bindings[contact.node] = Some(contact.friction.clone());
    }
    let routes: Vec<_> = ordered.iter().map(|c| c.route).collect();
    let step = FrictionHeatNetworkStep::new(
        thermal,
        target_s,
        &routes,
        controls.maximum_interfaces,
        check,
    )?;
    let trials = Cell::new(0usize);
    let mechanical_solves = Cell::new(0usize);
    // The accepted state, load, surface motion, bulk law, depth and source
    // bindings are immutable in this call. Only interface T affects mechanics.
    // Keep one O(nodes) result, never the mechanical Hessian or thermal payload.
    let cache = RefCell::new(None::<CachedMechanical<'p, 'a>>);
    let evaluate = |x: &[f64], entropy: bool| {
        check()?;
        if trials.get() >= controls.maximum_trials {
            return Err(SolidContactHeatError::TrialBudget);
        }
        trials.set(trials.get() + 1);
        if x.len() != dimension || x.iter().any(|v| !v.is_finite()) {
            return Err(SolidContactHeatError::Invalid(
                "complete finite coupled iterate required",
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
            .map(|(cell, &q)| resolvable_increment_j(cell.specific_h_j_kg, q, cell.mass_kg))
            .collect();
        let endpoint = |state: &DepositionSimulation<'t>, temperatures: &[f64]| {
            check()?;
            let cached = cache
                .borrow()
                .as_ref()
                .filter(|entry| entry.temperatures_k.as_slice() == temperatures)
                .map(|entry| (entry.trajectory.clone(), entry.solution.clone()));
            let (trajectory, mechanical) = if let Some(cached) = cached {
                cached
            } else {
                let mut local = bindings.clone();
                for (contact, &temperature) in ordered.iter().zip(temperatures) {
                    let binding = local[contact.node].as_mut().expect("validated bound node");
                    binding.conditions = binding
                        .conditions
                        .clone()
                        .with(&binding.law.axes().temperature, temperature)
                        .map_err(|error| {
                            SolidContactHeatError::Mechanical(SolidError::InvalidInput {
                                what: error.to_string(),
                            })
                        })?;
                }
                let mut trajectory = accepted.clone();
                mechanical_solves.set(mechanical_solves.get() + 1);
                let mechanical = trajectory
                    .advance_with_sources(step.duration_s(), load, surface_displacement, &local)
                    .map_err(SolidContactHeatError::Mechanical)?;
                check()?;
                *cache.borrow_mut() = Some(CachedMechanical {
                    temperatures_k: temperatures.to_vec(),
                    trajectory: trajectory.clone(),
                    solution: mechanical.clone(),
                });
                (trajectory, mechanical)
            };
            check()?;
            if trajectory.time_s() != target_s
                || mechanical.friction_work_per_node.len() != nodes
                || mechanical.friction_sources.len() != nodes
            {
                return Err(SolidContactHeatError::Invalid(
                    "incomplete native nodal friction output",
                ));
            }
            let source = resolve(state, step.link_ids(), temperatures)
                .map_err(SolidContactHeatError::Source)?;
            if source.interfaces.len() != count
                || source
                    .interfaces
                    .iter()
                    .zip(step.link_ids())
                    .any(|(s, id)| s.link_id != *id)
            {
                return Err(SolidContactHeatError::Invalid(
                    "original conductivities must retain canonical link order",
                ));
            }
            let mut interfaces = Vec::with_capacity(count);
            let (mut exported_work_j, mut accounted_j) = (0.0, 0.0);
            for (contact, conductivity) in ordered.iter().zip(&source.interfaces) {
                if mechanical.friction_sources[contact.node].is_none() {
                    return Err(SolidContactHeatError::Invalid(
                        "fixed heat footprint lost native contact",
                    ));
                }
                let (heat, other) = work::partition(
                    mechanical.friction_work_per_node[contact.node],
                    depth_m,
                    contact.heat_fraction,
                )
                .map_err(SolidContactHeatError::Invalid)?;
                accounted_j += heat + other;
                exported_work_j += other;
                interfaces.push(InterfaceHeatInput {
                    dissipated_j: heat,
                    conductivity_w_m_k: conductivity.conductivity_w_m_k,
                });
            }
            let original_j = depth_m * mechanical.friction_work;
            if !original_j.is_finite()
                || !accounted_j.is_finite()
                || !exported_work_j.is_finite()
                || (original_j - accounted_j).abs() > controls.port_energy_tolerance_j
            {
                return Err(SolidContactHeatError::Invalid(
                    "nodal heat and outside work do not close native friction loss",
                ));
            }
            let transport = fs_thermal::enthalpy::backward_euler::StepCoefficients {
                sources_w: source.transport.sources_w,
                conductances_w_k: source.transport.conductances_w_k,
                payload: Staged {
                    trajectory,
                    mechanical,
                    exported_work_j,
                    source: source.transport.payload,
                },
            };
            Ok(FrictionHeatNetworkInput {
                interfaces,
                transport,
            })
        };
        if entropy {
            step.evaluate_with_entropy(&increments, &temperatures, endpoint, check)
        } else {
            step.evaluate(&increments, &temperatures, endpoint, check)
        }
    };
    let numerical = |error| match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(field) => SolidContactHeatError::Invalid(field),
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
    let solved = solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
        gate.is_requested() || problem.has_failed()
    });
    if let Some(error) = problem.take_failure() {
        return Err(numerical(error));
    }
    check()?;
    if solved.progress != SolverRunProgress::Complete {
        return Err(SolidContactHeatError::Interrupted(solved.progress));
    }
    if !solved.report.converged {
        return Err(SolidContactHeatError::NotConverged(solved.report));
    }
    let trial = evaluate(&solver.x, controls.audit_entropy)?;
    let staged = trial.payload();
    let heat = trial.report();
    let body_work_j = depth_m
        * staged
            .mechanical
            .friction_reactions
            .iter()
            .enumerate()
            .map(|(i, force)| {
                force[0] * (staged.mechanical.u[2 * i] - accepted.displacement()[2 * i])
                    + force[1]
                        * (staged.mechanical.u[2 * i + 1] - accepted.displacement()[2 * i + 1])
            })
            .sum::<f64>();
    let port_energy_residual_j =
        body_work_j + heat.thermal.stored_energy_change_j + staged.exported_work_j
            - depth_m * staged.mechanical.friction_support_work
            - heat.external_heat_j;
    if !port_energy_residual_j.is_finite()
        || port_energy_residual_j.abs() > controls.port_energy_tolerance_j
    {
        return Err(SolidContactHeatError::Invalid(
            "native solid friction and finite heat ports do not balance",
        ));
    }
    let (thermal, staged, heat) = trial.accept(controls.heat, check)?;
    check()?;
    Ok(HeatedSolidContactCandidate {
        trajectory: staged.trajectory,
        mechanical: staged.mechanical,
        thermal,
        heat,
        nodes: ordered.iter().map(|c| c.node).collect(),
        exported_nonthermal_work_j: staged.exported_work_j,
        source: staged.source,
        nonlinear: solved.report,
        trials: trials.get(),
        mechanical_solves: mechanical_solves.get(),
        port_energy_residual_j,
    })
}
