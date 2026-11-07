//! Sustained plastic-shell contact and simultaneous finite-body heat transfer.
//!
//! Original shell layers, contact forces, material histories, cell enthalpies
//! and common interface temperatures belong to the same tentative endpoint.
//! Coulomb sliding heat enters the existing massless-interface network once;
//! each body's conductivity and finite caloric inventory determine its share.
//! Bulk physical/reversible heat uses the original section audit. Numerical
//! drilling, inertia and normal approach work never heat either body.
//!
//! This is a fixed-reference skin footprint against prescribed rigid obstacle
//! motion, with finite obstacle thermal cells. Each mapped interface must remain
//! closed throughout its declared affine gap path. Opening, changing footprints,
//! contact resistance, wear, flash-temperature models and unknown-force-dependent
//! friction need their own coupled owners. The source callback resolves explicit
//! endpoint Coulomb coefficients at the solved interface temperature; it must
//! not freeze an unevaluated pressure/rate/history dependence into that number.

mod geometry;

use super::feedback::Mechanical;
use super::*;
use crate::friction_heat::network::{
    FrictionHeatNetworkInput, FrictionHeatNetworkReport, FrictionHeatNetworkStep,
    InterfaceHeatInput,
};
use crate::friction_heat::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_solid::plastic_shell::contact::{
    PlasticShellContacts, ShellContactControls, ShellContactReport, ShellContactStep,
    ShellContactWitness,
};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::{StepCoefficients, resolvable_increment_j};
use std::cell::{Cell, RefCell};
use std::collections::BTreeMap;

/// One physical skin witness and its original finite thermal contact patch.
pub struct ShellContactHeatBinding<'a> {
    pub contact_id: u64,
    pub route: &'a FrictionHeatRoute,
    /// Index in `route.cells()`: the shell's exposed outer layer is this side.
    pub shell_side: usize,
    /// Positive nonlinear seed, never a prescribed interface temperature.
    pub temperature_guess_k: f64,
}

/// Explicit endpoint law and original-body conductivities for one thermal link.
#[derive(Clone, Copy, Debug)]
pub struct ShellContactInterfaceSource {
    pub link_id: u64,
    pub friction_coefficient: f64,
    pub conductivity_w_m_k: [f64; 2],
}

pub struct ShellContactHeatSource<P> {
    /// Exactly the callback's ascending original link order.
    pub interfaces: Vec<ShellContactInterfaceSource>,
    /// Other heat and other links only. Do not include bulk plastic or friction
    /// heat, or override a bound interface conductance a second time.
    pub transport: StepCoefficients<P>,
}

/// Independent budgets and physical acceptance gates for one coupled interval.
#[derive(Clone, Copy, Debug)]
pub struct ShellContactHeatControls {
    pub mechanical: ShellContactControls,
    pub maximum_cells: usize,
    pub maximum_unknowns: usize,
    pub maximum_trials: usize,
    pub maximum_newton_iterations: usize,
    pub difference_step_j: f64,
    /// Numerical coordinate scale only; no interface heat capacity is added.
    pub temperature_scale_j_per_k: f64,
    pub heat: FrictionHeatTolerances,
    pub material_energy_tolerance_j: f64,
    pub yield_tolerance_pa: f64,
    pub combined_energy_tolerance_j: f64,
    pub geometry_tolerance_m: f64,
    pub newton: NewtonKrylovConfig,
    /// Exact caloric entropy is audited on the final trial only, when available.
    pub audit_entropy: bool,
}

#[derive(Debug)]
pub enum ShellContactHeatError<E> {
    Invalid(&'static str),
    Shell(ShellThermoplasticError),
    Thermal(DepositionError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
}
impl<E> From<DepositionError> for ShellContactHeatError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E> From<ShellThermoplasticError> for ShellContactHeatError<E> {
    fn from(error: ShellThermoplasticError) -> Self {
        Self::Shell(error)
    }
}
impl<E> From<NewtonError> for ShellContactHeatError<E> {
    fn from(error: NewtonError) -> Self {
        Self::Newton(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for ShellContactHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "shell contact/thermal coupling refused: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for ShellContactHeatError<E> {}

/// Complete next owners and independent material, interface and system accounts.
pub struct ShellContactHeatCandidate<
    'm,
    't,
    P,
    S: ThermalSectionSourceFamily<'m> = ThermoplasticSource<'m>,
> {
    pub shell: ThermoplasticShell<'m, S>,
    pub thermal: DepositionSimulation<'t>,
    pub mechanical: ShellThermoplasticEnergy,
    pub cell_heat_j: BTreeMap<CellId, ShellLayerHeat>,
    pub motion: Option<ShellThermoplasticMotion>,
    /// Original mechanical witness order, independent of canonical thermal IDs.
    pub contacts: Vec<ShellContactReport>,
    pub contact_solve: fs_contact::sustained::solve::SolveReport,
    pub heat: FrictionHeatNetworkReport,
    pub interfaces: Vec<ShellContactInterfaceSource>,
    pub source: P,
    /// Outside supply, excluding internally transferred bulk and friction heat.
    pub external_heat_j: f64,
    /// Full material + caloric + kinetic + numerical accounts, less actual
    /// dead-load/support/obstacle work and outside heat.
    pub combined_residual_j: f64,
    pub nonlinear: NewtonReport,
    pub trials: usize,
    pub mechanical_solves: usize,
}

struct Staged<S, P> {
    mechanical: Mechanical<S>,
    interfaces: Vec<ShellContactInterfaceSource>,
    source: P,
}

struct Cached<S> {
    enthalpies: Vec<Vec<f64>>,
    friction: Vec<f64>,
    mechanical: Mechanical<S>,
}

impl<'a, S: ThermalSectionSourceFamily<'a>> ThermoplasticShell<'a, S> {
    /// Prepare immutable physical point/force rows once, independently of later
    /// temperature-dependent section laws. The native owner validates witnesses,
    /// budgets and the original mesh; no mechanical state is advanced here.
    pub fn prepare_contacts(
        &self,
        witnesses: &[ShellContactWitness],
        maximum_contacts: usize,
        cx: &Cx<'_>,
    ) -> Result<PlasticShellContacts<'a>, ShellThermoplasticError> {
        poll(cx)?;
        let d = self.discretization;
        let mut sections = Vec::with_capacity(self.facets.len());
        for facet in self.facets.iter() {
            poll(cx)?;
            sections.push(ThermoplasticSectionStep::new(
                &facet.layers,
                &facet.reference_h_j_kg,
                &facet.reference_h_j_kg,
                1.0,
                d.quadrature,
                d.reduction,
                d.section_budget,
            )?);
        }
        let references: Vec<_> = sections
            .iter()
            .map(ThermoplasticSectionStep::section)
            .collect();
        let fixed: Vec<_> = (0..6 * self.mesh.nodes.len()).map(|d| (d, 0.0)).collect();
        let problem = SmallStrainShellProblem {
            mesh: self.mesh,
            sections: PlasticShellSections::PerElementBorrowed(&references),
            coordinate_length_m: d.coordinate_length_m,
            rotation_limit_rad: d.rotation_limit_rad,
            drilling_alpha: d.drilling_alpha,
            prescribed: &fixed,
            nodal_loads: &[],
            load_path: &[1.0],
            budget: d.shell_budget,
            settings: IncrementSettings::default(),
        };
        PlasticShellContacts::prepare(self.mesh, &problem, witnesses, maximum_contacts, cx)
            .map_err(ShellThermoplasticError::Contact)
    }

    /// Solve native sustained contact, sourced bulk plasticity and original
    /// enthalpy cells at one shared endpoint. `motion` selects the existing
    /// consistent-mass backward-Euler path without changing reference mass.
    /// Every trial starts from the same accepted histories and velocities.
    ///
    /// All contact witnesses need one original thermal route. Each route must
    /// match its skin position, outward face, area and layer inventory and remain
    /// closed; this method does not keep a thermal bridge across an open gap.
    /// `resolve` sees canonical link IDs and their solved common temperatures;
    /// its source receipt stays tentative until all owners accept. Fixed original
    /// thermal topology is required; split the interval at scheduled events.
    /// Energy guesses follow the thermal owner's ascending CellId order.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance_contact_heat<'t, P, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        prescribed: &[(usize, f64)],
        nodal_loads: &[(usize, f64)],
        contacts: &PlasticShellContacts<'_>,
        contact_steps: &[ShellContactStep],
        bindings: &[ShellContactHeatBinding<'_>],
        motion: Option<(&[[f64; 6]], &PlasticShellReferenceInertia<'_>)>,
        settings: IncrementSettings,
        controls: ShellContactHeatControls,
        energy_guess_j: &[f64],
        resolve: impl Fn(
            &DepositionSimulation<'t>,
            &[u64],
            &[f64],
        ) -> Result<ShellContactHeatSource<P>, E>,
        cx: &Cx<'_>,
    ) -> Result<ShellContactHeatCandidate<'a, 't, P, S>, ShellContactHeatError<E>> {
        let check = || cx.checkpoint().map_err(|_| DepositionError::Cancelled);
        check()?;
        controls.newton.validate()?;
        let dt = target_s - self.time_s;
        let cells = thermal.cells().len();
        let dimension = cells
            .checked_add(bindings.len())
            .ok_or(ShellContactHeatError::Invalid("coupled dimension overflow"))?;
        if thermal.time_s() != self.time_s
            || !target_s.is_finite()
            || !positive(dt)
            || cells == 0
            || cells > controls.maximum_cells
            || bindings.is_empty()
            || bindings.len() != contacts.witnesses().len()
            || contact_steps.len() != bindings.len()
            || bindings.len() > controls.mechanical.solve.maximum_contacts
            || dimension > controls.maximum_unknowns
            || controls.maximum_trials == 0
            || controls.maximum_newton_iterations == 0
            || energy_guess_j.len() != cells
            || energy_guess_j.iter().any(|q| !q.is_finite())
            || !positive(controls.difference_step_j)
            || !positive(controls.temperature_scale_j_per_k)
            || !positive(controls.material_energy_tolerance_j)
            || !positive(controls.yield_tolerance_pa)
            || !positive(controls.geometry_tolerance_m)
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
            return Err(ShellContactHeatError::Invalid(
                "synchronized bounded shell/contact/thermal interval and finite controls required",
            ));
        }
        self.enthalpies::<Infallible>(thermal, true, cx)?;
        if let Some((_, reference)) = motion {
            for (e, facet) in self.facets.iter().enumerate() {
                check()?;
                let density = reference.reference_layer_corner_densities(e).ok_or(
                    ShellContactHeatError::Invalid("missing shell reference mass inventory"),
                )?;
                if density.len() != facet.layers.len()
                    || density
                        .zip(&facet.layers)
                        .any(|(rho, layer)| rho != [layer.law.reference_density_kg_m3(); 3])
                {
                    return Err(ShellContactHeatError::Invalid(
                        "contact motion must retain every original sourced layer density",
                    ));
                }
            }
        }
        let mut ordered: Vec<_> = bindings.iter().collect();
        ordered.sort_by_key(|binding| binding.route.link_id());
        let indices = geometry::validate_bindings(
            self,
            thermal,
            contacts,
            &ordered,
            controls.geometry_tolerance_m,
            cx,
        )?;
        let routes: Vec<_> = ordered.iter().map(|b| b.route).collect();
        let network = FrictionHeatNetworkStep::new(
            thermal,
            target_s,
            &routes,
            controls.mechanical.solve.maximum_contacts,
            check,
        )?;
        let trials = Cell::new(0usize);
        let mechanical_solves = Cell::new(0usize);
        // One bounded endpoint is retained. Exact repeated thermal states and
        // Coulomb coefficients reuse its native response; transport/source
        // receipts still refresh on every thermal trial.
        let cache = RefCell::new(None::<Cached<History<'a, S>>>);
        let material_controls = ThermoplasticConfig {
            maximum_iterations: 1,
            temperature_tolerance_k: 1.0,
            energy_tolerance_j: controls.material_energy_tolerance_j,
            yield_tolerance_pa: controls.yield_tolerance_pa,
            enthalpy_relaxation: 1.0,
        };
        let evaluate = |x: &[f64], entropy: bool| {
            check()?;
            if trials.get() >= controls.maximum_trials {
                return Err(ShellContactHeatError::TrialBudget);
            }
            trials.set(trials.get() + 1);
            if x.len() != dimension || x.iter().any(|v| !v.is_finite()) {
                return Err(ShellContactHeatError::Invalid(
                    "finite complete coupled iterate required",
                ));
            }
            let temperatures: Vec<_> = ordered
                .iter()
                .enumerate()
                .map(|(i, b)| {
                    b.temperature_guess_k + x[cells + i] / controls.temperature_scale_j_per_k
                })
                .collect();
            let increments: Vec<_> = thermal
                .cells()
                .values()
                .zip(&x[..cells])
                .map(|(cell, &q)| resolvable_increment_j(cell.specific_h_j_kg, q, cell.mass_kg))
                .collect();
            let endpoint = |state: &DepositionSimulation<'t>, temperatures: &[f64]| {
                check()?;
                let enthalpies = self.enthalpies::<Infallible>(state, false, cx)?;
                let source = resolve(state, network.link_ids(), temperatures)
                    .map_err(ShellContactHeatError::Source)?;
                if source.interfaces.len() != ordered.len()
                    || source
                        .interfaces
                        .iter()
                        .zip(network.link_ids())
                        .any(|(s, id)| {
                            s.link_id != *id
                                || !s.friction_coefficient.is_finite()
                                || s.friction_coefficient < 0.0
                        })
                {
                    return Err(ShellContactHeatError::Invalid(
                        "complete canonical link-resolved Coulomb coefficients required",
                    ));
                }
                let friction: Vec<_> = source
                    .interfaces
                    .iter()
                    .map(|s| s.friction_coefficient)
                    .collect();
                let cached = cache
                    .borrow()
                    .as_ref()
                    .filter(|old| old.enthalpies == enthalpies && old.friction == friction)
                    .map(|old| old.mechanical.clone());
                let mechanical = if let Some(cached) = cached {
                    cached
                } else {
                    let mut steps = contact_steps.to_vec();
                    for (&i, source) in indices.iter().zip(&source.interfaces) {
                        steps[i].friction_coefficient = source.friction_coefficient;
                    }
                    mechanical_solves.set(mechanical_solves.get() + 1);
                    let mechanical = self.mechanics::<Infallible>(
                        &enthalpies,
                        prescribed,
                        nodal_loads,
                        settings,
                        material_controls,
                        dt,
                        motion,
                        Some((contacts, &steps, controls.mechanical)),
                        cx,
                    )?;
                    *cache.borrow_mut() = Some(Cached {
                        enthalpies,
                        friction,
                        mechanical: mechanical.clone(),
                    });
                    mechanical
                };
                let reports = &mechanical
                    .contacts
                    .as_ref()
                    .expect("native contact response retains contact reports")
                    .0;
                let mut interfaces = Vec::with_capacity(ordered.len());
                for ((&i, binding), source) in indices.iter().zip(&ordered).zip(&source.interfaces)
                {
                    check()?;
                    let report = &reports[i];
                    if report.id != binding.contact_id
                        || report.initial_gap_m.abs() > controls.geometry_tolerance_m
                        || report.final_gap_m.abs() > controls.geometry_tolerance_m
                    {
                        return Err(ShellContactHeatError::Invalid(
                            "fixed heat footprint must remain closed on its original affine path",
                        ));
                    }
                    interfaces.push(InterfaceHeatInput {
                        dissipated_j: report.friction_loss_j,
                        conductivity_w_m_k: source.conductivity_w_m_k,
                    });
                }
                let mut transport = source.transport;
                for (&cell, heat) in &mechanical.heat {
                    let power = transport.sources_w.entry(cell).or_insert(0.0);
                    *power += heat.total_j() / dt;
                    if !power.is_finite() {
                        return Err(ShellContactHeatError::Invalid(
                            "unrepresentable bulk heat power",
                        ));
                    }
                }
                Ok(FrictionHeatNetworkInput {
                    interfaces,
                    transport: StepCoefficients {
                        sources_w: transport.sources_w,
                        conductances_w_k: transport.conductances_w_k,
                        payload: Staged {
                            mechanical,
                            interfaces: source.interfaces,
                            source: transport.payload,
                        },
                    },
                })
            };
            if entropy {
                network.evaluate_with_entropy(&increments, &temperatures, endpoint, check)
            } else {
                network.evaluate(&increments, &temperatures, endpoint, check)
            }
        };
        let numerical = |error| match error {
            FiniteDifferenceError::Source(error) => error,
            FiniteDifferenceError::Invalid(field) => ShellContactHeatError::Invalid(field),
        };
        let residual =
            FallibleFiniteDifference::new(dimension, controls.difference_step_j, |x: &[f64]| {
                evaluate(x, false).map(|trial| trial.residual_j().to_vec())
            })
            .map_err(numerical)?;
        let mut initial = energy_guess_j.to_vec();
        initial.resize(dimension, 0.0);
        let construction = NewtonKrylovState::new(&residual, initial, controls.newton);
        if let Some(error) = residual.take_failure() {
            return Err(numerical(error));
        }
        let mut solver = construction?;
        let solved =
            solver.run_controlled(&residual, controls.maximum_newton_iterations, &mut || {
                cx.checkpoint().is_err() || residual.has_failed()
            });
        if let Some(error) = residual.take_failure() {
            return Err(numerical(error));
        }
        check()?;
        if solved.progress != SolverRunProgress::Complete {
            return Err(ShellContactHeatError::Interrupted(solved.progress));
        }
        if !solved.report.converged {
            return Err(ShellContactHeatError::NotConverged(solved.report));
        }
        let trial = evaluate(&solver.x, controls.audit_entropy)?;
        let staged = trial.payload();
        let a = &staged.mechanical.energy;
        let heat = trial.report();
        let bulk_heat_j: f64 = staged.mechanical.heat.values().map(|q| q.total_j()).sum();
        let external_heat_j = heat.external_heat_j - bulk_heat_j;
        let contact_reports = &staged
            .mechanical
            .contacts
            .as_ref()
            .expect("native contact response retains contact reports")
            .0;
        let obstacle_work_j: f64 = contact_reports.iter().map(|r| r.obstacle_work_j).sum();
        let approach_loss_j: f64 = contact_reports
            .iter()
            .map(|r| r.normal_numerical_loss_j)
            .sum();
        let combined_residual_j = a.stored_change_j
            + a.thermoelastic_storage_change_j
            + heat.thermal.stored_energy_change_j
            + a.material_numerical_loss_j
            + a.drilling.final_j
            - a.drilling.initial_j
            + a.drilling.numerical_loss_j
            + staged.mechanical.motion.as_ref().map_or(0.0, |m| {
                m.inertia.total.kinetic_change_j + m.inertia.total.numerical_loss_j
            })
            + approach_loss_j
            - a.dead_load_work_j
            - a.support_work_j
            - obstacle_work_j
            - a.transverse_residual_work_j
            - external_heat_j
            - heat.thermal.incoming_energy_j;
        if !external_heat_j.is_finite()
            || !combined_residual_j.is_finite()
            || combined_residual_j.abs() > controls.combined_energy_tolerance_j
        {
            return Err(ShellContactHeatError::Invalid(
                "combined shell/contact/thermal first law failed",
            ));
        }
        let (next_thermal, staged, heat) = trial.accept(controls.heat, check)?;
        self.enthalpies::<Infallible>(&next_thermal, false, cx)?;
        let accepted = self
            .accepted
            .iter()
            .map(|cells| {
                cells
                    .iter()
                    .map(|cell| next_thermal.cells()[&cell.id].clone())
                    .collect()
            })
            .collect();
        let mechanical = staged.mechanical;
        let shell = Self {
            mesh: self.mesh,
            facets: self.facets.clone(),
            accepted,
            reference_drilling: self.reference_drilling.clone(),
            solution: mechanical.solution,
            discretization: self.discretization,
            area_m2: self.area_m2,
            volume_m3: self.volume_m3,
            time_s: target_s,
        };
        let (contacts, contact_solve) = mechanical
            .contacts
            .expect("accepted contact response retains contact reports");
        check()?;
        Ok(ShellContactHeatCandidate {
            shell,
            thermal: next_thermal,
            mechanical: mechanical.energy,
            cell_heat_j: mechanical.heat,
            motion: mechanical.motion,
            contacts,
            contact_solve,
            heat,
            interfaces: staged.interfaces,
            source: staged.source,
            external_heat_j,
            combined_residual_j,
            nonlinear: solved.report,
            trials: trials.get(),
            mechanical_solves: mechanical_solves.get(),
        })
    }
}
