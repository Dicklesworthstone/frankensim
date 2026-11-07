//! One frozen-history shell/enthalpy interval, using the original owners.

use super::*;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::thermal::MAX_FEEDBACK_ITERATIONS;
use fs_solid::plastic_shell::contact::{
    PlasticShellContacts, ShellContactControls, ShellContactReport, ShellContactStep,
};
use fs_thermal::enthalpy::{DepositionObserver, Event, StepConfig};
use fs_thermal::observers::{CoupledObservers, ObserverCouplingError};
use std::collections::BTreeMap;

/// Actual heat supplied to one original layer cell, in joules.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ShellLayerHeat {
    pub physical_j: f64,
    /// Creep work is a subaccount of `physical_j`, not additional heat.
    pub creep_j: f64,
    pub reversible_j: f64,
}
impl ShellLayerHeat {
    pub fn total_j(self) -> f64 {
        self.physical_j + self.reversible_j
    }
}

/// Independent thermal and whole-system work accounts for one accepted step.
#[derive(Clone, Debug)]
pub struct ShellThermoplasticReport {
    pub iterations: usize,
    pub thermal_audit: EnergyAudit,
    pub mechanical: ShellThermoplasticEnergy,
    /// Native endpoint heat; never an inferred residual or numerical loss.
    pub cell_heat_j: BTreeMap<CellId, ShellLayerHeat>,
    pub applied_heat_j: f64,
    pub external_heat_j: f64,
    pub thermal_residual_j: f64,
    pub combined_residual_j: f64,
    pub temperature_defect_k: f64,
    pub heat_defect_j: f64,
}

/// Complete accepted physical velocities and the original inertia partition.
#[derive(Clone, Debug)]
pub struct ShellThermoplasticMotion {
    /// Translation [m/s] followed by axial rotation [rad/s].
    pub velocity: Vec<[f64; 6]>,
    pub inertia: ShellInertialEnergy,
}

/// Mechanical, caloric and observer states are accepted together.
pub struct ShellThermoplasticCandidate<
    'm,
    't,
    O,
    S: ThermalSectionSourceFamily<'m> = ThermoplasticSource<'m>,
> {
    shell: ThermoplasticShell<'m, S>,
    thermal: DepositionSimulation<'t>,
    observer: O,
    report: ShellThermoplasticReport,
    motion: Option<ShellThermoplasticMotion>,
}
impl<'m, 't, O, S: ThermalSectionSourceFamily<'m>> ShellThermoplasticCandidate<'m, 't, O, S> {
    pub const fn shell(&self) -> &ThermoplasticShell<'m, S> {
        &self.shell
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn observer(&self) -> &O {
        &self.observer
    }
    pub const fn report(&self) -> &ShellThermoplasticReport {
        &self.report
    }
    /// Present only for an actual backward-Euler momentum step.
    pub fn motion(&self) -> Option<&ShellThermoplasticMotion> {
        self.motion.as_ref()
    }
    pub fn into_parts(
        self,
    ) -> (
        ThermoplasticShell<'m, S>,
        DepositionSimulation<'t>,
        O,
        ShellThermoplasticReport,
        Option<ShellThermoplasticMotion>,
    ) {
        (
            self.shell,
            self.thermal,
            self.observer,
            self.report,
            self.motion,
        )
    }
}

#[derive(Clone)]
pub(super) struct Mechanical<S> {
    pub(super) solution: SmallStrainShellSolution<S>,
    pub(super) temperatures_k: Vec<f64>,
    pub(super) heat: BTreeMap<CellId, ShellLayerHeat>,
    pub(super) energy: ShellThermoplasticEnergy,
    pub(super) motion: Option<ShellThermoplasticMotion>,
    pub(super) contacts: Option<(
        Vec<ShellContactReport>,
        fs_contact::sustained::solve::SolveReport,
    )>,
}

impl<'a, S: ThermalSectionSourceFamily<'a>> ThermoplasticShell<'a, S> {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn mechanics<E>(
        &self,
        enthalpies: &[Vec<f64>],
        prescribed: &[(usize, f64)],
        nodal_loads: &[(usize, f64)],
        settings: IncrementSettings,
        config: ThermoplasticConfig,
        duration_s: f64,
        motion: Option<(&[[f64; 6]], &PlasticShellReferenceInertia<'_>)>,
        contact: Option<(
            &PlasticShellContacts<'_>,
            &[ShellContactStep],
            ShellContactControls,
        )>,
        cx: &Cx<'_>,
    ) -> Result<Mechanical<History<'a, S>>, ShellThermoplasticError<E>> {
        let d = self.discretization;
        let old_h: Vec<Vec<_>> = self
            .accepted
            .iter()
            .map(|cells| cells.iter().map(|cell| cell.specific_h_j_kg).collect())
            .collect();
        let mut steps = Vec::with_capacity(self.facets.len());
        for (e, facet) in self.facets.iter().enumerate() {
            poll(cx)?;
            steps.push(ThermoplasticSectionStep::new(
                &facet.layers,
                &old_h[e],
                &enthalpies[e],
                duration_s,
                d.quadrature,
                d.reduction,
                d.section_budget,
            )?);
        }
        let references: Vec<_> = steps
            .iter()
            .map(ThermoplasticSectionStep::section)
            .collect();
        let problem = SmallStrainShellProblem {
            mesh: self.mesh,
            sections: PlasticShellSections::PerElementBorrowed(&references),
            coordinate_length_m: d.coordinate_length_m,
            rotation_limit_rad: d.rotation_limit_rad,
            drilling_alpha: d.drilling_alpha,
            prescribed,
            nodal_loads,
            load_path: &[1.0],
            budget: d.shell_budget,
            settings,
        };
        let (solution, drilling, motion, contacts) =
            if let Some((contacts, steps, controls)) = contact {
                let solved = if let Some((velocity, reference)) = motion {
                    problem.advance_backward_euler_with_contacts(
                        Some(&self.solution),
                        velocity,
                        reference,
                        &self.reference_drilling,
                        contacts,
                        steps,
                        duration_s,
                        controls,
                        cx,
                    )
                } else {
                    problem.solve_from_with_contacts(
                        &self.solution.coordinates,
                        &self.solution.section_states,
                        &self.reference_drilling,
                        contacts,
                        steps,
                        controls,
                        cx,
                    )
                }
                .map_err(ShellThermoplasticError::Contact)?;
                (
                    solved.solution,
                    solved.drilling,
                    solved.motion.map(|m| ShellThermoplasticMotion {
                        velocity: m.velocity,
                        inertia: m.inertia,
                    }),
                    Some((solved.contacts, solved.contact_solve)),
                )
            } else if let Some((velocity, reference)) = motion {
                let solved = problem
                    .advance_backward_euler_with_reference_drilling(
                        Some(&self.solution),
                        velocity,
                        reference,
                        &self.reference_drilling,
                        duration_s,
                        cx,
                    )
                    .map_err(ShellThermoplasticError::Solid)?;
                (
                    solved.solution,
                    solved.drilling,
                    Some(ShellThermoplasticMotion {
                        velocity: solved.velocity,
                        inertia: solved.inertia,
                    }),
                    None,
                )
            } else {
                let solved = problem
                    .solve_from_with_reference_drilling(
                        &self.solution.coordinates,
                        &self.solution.section_states,
                        &self.reference_drilling,
                        cx,
                    )
                    .map_err(ShellThermoplasticError::Solid)?;
                (solved.solution, solved.drilling, None, None)
            };
        let mut energy = ShellThermoplasticEnergy {
            drilling,
            contact_work_j: finite(contacts.as_ref().map_or(0.0, |(reports, _)| {
                reports.iter().map(|report| report.shell_work_j).sum()
            }))?,
            inertial_work_j: motion
                .as_ref()
                .map_or(0.0, |m| m.inertia.total.inertial_work_j),
            ..Default::default()
        };
        let mut heat = BTreeMap::new();
        for (e, (step, facet)) in steps.iter().zip(self.facets.iter()).enumerate() {
            let area = facet.area_m2 / 3.0;
            for q in 0..3 {
                poll(cx)?;
                let audit = step.audit_energy(
                    &self.solution.section_states[3 * e + q],
                    &solution.section_states[3 * e + q],
                    EnergyTolerances {
                        energy_j_m3: config.energy_tolerance_j / self.volume_m3,
                        yield_pa: config.yield_tolerance_pa,
                    },
                    config.energy_tolerance_j / self.area_m2,
                )?;
                if audit.total.resultants != solution.section_resultants[3 * e + q] {
                    return Err(ShellThermoplasticError::Invalid(
                        "thermal work and shell equilibrium must use the same native resultant",
                    ));
                }
                for (layer, cell) in audit.layers.iter().zip(&self.accepted[e]) {
                    let entry = heat.entry(cell.id).or_insert(ShellLayerHeat::default());
                    entry.physical_j = finite(entry.physical_j + area * layer.physical_heat_j_m2)?;
                    entry.creep_j = finite(entry.creep_j + area * layer.creep_heat_j_m2)?;
                    entry.reversible_j =
                        finite(entry.reversible_j + area * layer.reversible_heat_j_m2)?;
                }
                let a = audit.total;
                energy.physical_heat_j =
                    finite(energy.physical_heat_j + area * a.physical_heat_j_m2)?;
                energy.creep_heat_j = finite(energy.creep_heat_j + area * a.creep_heat_j_m2)?;
                energy.reversible_heat_j =
                    finite(energy.reversible_heat_j + area * a.reversible_heat_j_m2)?;
                energy.stored_change_j =
                    finite(energy.stored_change_j + area * a.stored_change_j_m2)?;
                energy.coefficient_storage_change_j = finite(
                    energy.coefficient_storage_change_j + area * a.coefficient_storage_change_j_m2,
                )?;
                energy.thermoelastic_storage_change_j = finite(
                    energy.thermoelastic_storage_change_j
                        + area * a.thermoelastic_storage_change_j_m2,
                )?;
                energy.material_numerical_loss_j =
                    finite(energy.material_numerical_loss_j + area * a.numerical_loss_j_m2)?;
                energy.material_work_j =
                    finite(energy.material_work_j + area * a.mechanical_work_j_m2)?;
                energy.generalized_work_j =
                    finite(energy.generalized_work_j + area * a.generalized_work_j_m2)?;
                energy.transverse_residual_work_j = finite(
                    energy.transverse_residual_work_j + area * a.transverse_residual_work_j_m2,
                )?;
            }
        }
        let delta =
            |d: usize| solution.coordinates[d / 6][d % 6] - self.solution.coordinates[d / 6][d % 6];
        for &(dof, load) in nodal_loads {
            energy.dead_load_work_j = finite(energy.dead_load_work_j + load * delta(dof))?;
        }
        for (&(dof, _), &reaction) in prescribed.iter().zip(&solution.increments[0].reactions) {
            energy.support_work_j = finite(energy.support_work_j + reaction * delta(dof))?;
        }
        energy.equilibrium_residual_j = finite(
            energy.external_work_j()
                - energy.generalized_work_j
                - energy.drilling.work_j
                - energy.inertial_work_j,
        )?;
        if energy.equilibrium_residual_j.abs() > config.energy_tolerance_j {
            return Err(ShellThermoplasticError::Invalid(
                "shell endpoint misses the independent force/support work gate",
            ));
        }
        Ok(Mechanical {
            solution,
            temperatures_k: steps
                .iter()
                .flat_map(|s| s.temperatures_k().iter().copied())
                .collect(),
            heat,
            energy,
            motion,
            contacts,
        })
    }

    /// Solve endpoint loads/support motion and reciprocal heating together.
    /// Original physical and reversible heat is distributed over this interval;
    /// it is not re-applied for each equilibrium iteration or thermal substep.
    /// Each trial clones the same original observer. Its own conductances,
    /// events, reservoirs and interval hooks run on the shared thermal timeline.
    #[allow(clippy::too_many_arguments)]
    pub fn advance<'t, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        observer: &O,
        target_s: f64,
        prescribed: &[(usize, f64)],
        nodal_loads: &[(usize, f64)],
        settings: IncrementSettings,
        thermal_config: StepConfig,
        config: ThermoplasticConfig,
        cx: &Cx<'_>,
    ) -> Result<ShellThermoplasticCandidate<'a, 't, O, S>, ShellThermoplasticError<O::Error>> {
        self.advance_with_motion(
            thermal,
            observer,
            target_s,
            prescribed,
            nodal_loads,
            settings,
            thermal_config,
            config,
            None,
            cx,
        )
    }

    /// Evolve the same layer heat/history feedback with native consistent-mass
    /// momentum. The actual physical interval owns both thermal rate laws and
    /// backward Euler. Physical director kinetics and numerical drilling/inertial
    /// loss stay in their original separate accounts; neither heats a layer.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_backward_euler<'t, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        observer: &O,
        target_s: f64,
        prescribed: &[(usize, f64)],
        nodal_loads: &[(usize, f64)],
        velocity: &[[f64; 6]],
        reference: &PlasticShellReferenceInertia<'_>,
        settings: IncrementSettings,
        thermal_config: StepConfig,
        config: ThermoplasticConfig,
        cx: &Cx<'_>,
    ) -> Result<ShellThermoplasticCandidate<'a, 't, O, S>, ShellThermoplasticError<O::Error>> {
        self.advance_with_motion(
            thermal,
            observer,
            target_s,
            prescribed,
            nodal_loads,
            settings,
            thermal_config,
            config,
            Some((velocity, reference)),
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_with_motion<'t, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        observer: &O,
        target_s: f64,
        prescribed: &[(usize, f64)],
        nodal_loads: &[(usize, f64)],
        settings: IncrementSettings,
        thermal_config: StepConfig,
        config: ThermoplasticConfig,
        motion: Option<(&[[f64; 6]], &PlasticShellReferenceInertia<'_>)>,
        cx: &Cx<'_>,
    ) -> Result<ShellThermoplasticCandidate<'a, 't, O, S>, ShellThermoplasticError<O::Error>> {
        poll(cx)?;
        let dt = target_s - self.time_s;
        if thermal.time_s() != self.time_s
            || !target_s.is_finite()
            || !positive(dt)
            || config.maximum_iterations == 0
            || config.maximum_iterations > MAX_FEEDBACK_ITERATIONS
            || !positive(config.temperature_tolerance_k)
            || !positive(config.energy_tolerance_j)
            || !positive(config.yield_tolerance_pa)
            || !positive(config.enthalpy_relaxation)
            || config.enthalpy_relaxation > 1.0
            || !thermal_config.energy_tolerance_j.is_finite()
            || thermal_config.energy_tolerance_j < 0.0
            || thermal_config.energy_tolerance_j > config.energy_tolerance_j
        {
            return Err(ShellThermoplasticError::Invalid(
                "matching clocks, positive physical interval and bounded thermal/work controls required",
            ));
        }
        let mut enthalpies = self.enthalpies(thermal, true, cx)?;
        if let Some((_, reference)) = motion {
            // Equal total mass is insufficient: swapping layer densities also
            // changes eccentric translation/director inertia. Retain every
            // admitted source density, in the original physical layer order.
            for (e, facet) in self.facets.iter().enumerate() {
                poll(cx)?;
                let density = reference.reference_layer_corner_densities(e).ok_or(
                    ShellThermoplasticError::Invalid("missing original shell layer mass inventory"),
                )?;
                if density.len() != facet.layers.len()
                    || density
                        .zip(&facet.layers)
                        .any(|(rho, layer)| rho != [layer.law.reference_density_kg_m3(); 3])
                {
                    return Err(ShellThermoplasticError::Invalid(
                        "shell inertia must retain the original sourced density of every layer",
                    ));
                }
            }
        }
        let mut defects = [f64::INFINITY; 3];
        for iteration in 1..=config.maximum_iterations {
            poll(cx)?;
            let trial = self.mechanics(
                &enthalpies,
                prescribed,
                nodal_loads,
                settings,
                config,
                dt,
                motion,
                None,
                cx,
            )?;
            let rates = trial
                .heat
                .iter()
                .map(|(&cell, heat)| Ok((cell, finite(heat.total_j() / dt)?)))
                .collect::<Result<_, ShellThermoplasticError<O::Error>>>()?;
            let heat = CoupledObservers::new(
                LayerHeat {
                    shell: self,
                    rates,
                    cx,
                },
                observer.clone(),
            );
            let (next_thermal, heat, thermal_audit) = thermal
                .advance_observed(&heat, target_s, thermal_config, || {
                    cx.checkpoint().map_err(|_| DepositionError::Cancelled)
                })
                .map_err(coupled_error)?;
            let final_h = self.enthalpies(&next_thermal, false, cx)?;
            let actual = self.mechanics(
                &final_h,
                prescribed,
                nodal_loads,
                settings,
                config,
                dt,
                motion,
                None,
                cx,
            )?;
            let applied_heat_j = finite(trial.heat.values().map(|q| q.total_j()).sum())?;
            let caloric_heat_j = finite(actual.heat.values().map(|q| q.total_j()).sum())?;
            let external_heat_j = finite(thermal_audit.external_heat_j - applied_heat_j)?;
            let thermal_residual_j = finite(
                thermal_audit.stored_energy_change_j
                    - thermal_audit.incoming_energy_j
                    - external_heat_j
                    - caloric_heat_j,
            )?;
            let a = &actual.energy;
            // The full tensor work contains the signed finite plane-stress
            // residual. It is an explicit constraint-work term, never heat.
            let combined_residual_j = finite(
                a.stored_change_j
                    + a.thermoelastic_storage_change_j
                    + thermal_audit.stored_energy_change_j
                    + a.material_numerical_loss_j
                    + a.drilling.final_j
                    - a.drilling.initial_j
                    + a.drilling.numerical_loss_j
                    + actual.motion.as_ref().map_or(0.0, |m| {
                        m.inertia.total.kinetic_change_j + m.inertia.total.numerical_loss_j
                    })
                    - a.external_work_j()
                    - a.transverse_residual_work_j
                    - external_heat_j
                    - thermal_audit.incoming_energy_j,
            )?;
            defects = [
                trial
                    .temperatures_k
                    .iter()
                    .zip(&actual.temperatures_k)
                    .map(|(a, b)| (a - b).abs())
                    .fold(0.0, f64::max),
                finite(
                    trial
                        .heat
                        .iter()
                        .map(|(cell, q)| (q.total_j() - actual.heat[cell].total_j()).abs())
                        .sum(),
                )?,
                combined_residual_j,
            ];
            if defects[0] <= config.temperature_tolerance_k
                && defects[1] <= config.energy_tolerance_j
                && thermal_residual_j.abs() <= config.energy_tolerance_j
                && combined_residual_j.abs() <= config.energy_tolerance_j
            {
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
                let shell = Self {
                    mesh: self.mesh,
                    facets: self.facets.clone(),
                    accepted,
                    reference_drilling: self.reference_drilling.clone(),
                    solution: actual.solution,
                    discretization: self.discretization,
                    area_m2: self.area_m2,
                    volume_m3: self.volume_m3,
                    time_s: target_s,
                };
                let report = ShellThermoplasticReport {
                    iterations: iteration,
                    thermal_audit,
                    mechanical: actual.energy,
                    cell_heat_j: actual.heat,
                    applied_heat_j,
                    external_heat_j,
                    thermal_residual_j,
                    combined_residual_j,
                    temperature_defect_k: defects[0],
                    heat_defect_j: defects[1],
                };
                poll(cx)?;
                return Ok(ShellThermoplasticCandidate {
                    shell,
                    thermal: next_thermal,
                    observer: heat.into_parts().1,
                    report,
                    motion: actual.motion,
                });
            }
            for (h, actual) in enthalpies
                .iter_mut()
                .flatten()
                .zip(final_h.into_iter().flatten())
            {
                *h = finite(
                    (1.0 - config.enthalpy_relaxation) * *h + config.enthalpy_relaxation * actual,
                )?;
            }
        }
        poll(cx)?;
        Err(ShellThermoplasticError::NotConverged {
            iterations: config.maximum_iterations,
            temperature_defect_k: defects[0],
            heat_defect_j: defects[1],
            combined_residual_j: defects[2],
        })
    }
}

#[derive(Clone)]
struct LayerHeat<'s, 'a, 'cx, S: ThermalSectionSourceFamily<'a>> {
    shell: &'s ThermoplasticShell<'a, S>,
    rates: BTreeMap<CellId, f64>,
    cx: &'s Cx<'cx>,
}
impl<'a, S: ThermalSectionSourceFamily<'a>> DepositionObserver for LayerHeat<'_, 'a, '_, S> {
    type Error = ShellThermoplasticError;
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        check()?;
        self.shell.enthalpies(thermal, false, self.cx)?;
        Ok(self.rates.clone())
    }
    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.shell.enthalpies(thermal, false, self.cx).map(|_| ())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        check()?;
        self.shell.enthalpies(thermal, false, self.cx).map(|_| ())
    }
}

fn coupled_error<E>(
    error: ObserverCouplingError<ShellThermoplasticError, E>,
) -> ShellThermoplasticError<E> {
    match error {
        ObserverCouplingError::Thermal(e) => e.into(),
        ObserverCouplingError::Second(e) => ShellThermoplasticError::Observer(e),
        ObserverCouplingError::First(e) => match e {
            ShellThermoplasticError::Invalid(e) => ShellThermoplasticError::Invalid(e),
            ShellThermoplasticError::Material(e) => ShellThermoplasticError::Material(e),
            ShellThermoplasticError::Section(e) => ShellThermoplasticError::Section(e),
            ShellThermoplasticError::Solid(e) => ShellThermoplasticError::Solid(e),
            ShellThermoplasticError::Contact(e) => ShellThermoplasticError::Contact(e),
            ShellThermoplasticError::Thermal(e) => e.into(),
            ShellThermoplasticError::Observer(e) => match e {},
            ShellThermoplasticError::Cancelled => ShellThermoplasticError::Cancelled,
            ShellThermoplasticError::NotConverged {
                iterations,
                temperature_defect_k,
                heat_defect_j,
                combined_residual_j,
            } => ShellThermoplasticError::NotConverged {
                iterations,
                temperature_defect_k,
                heat_defect_j,
                combined_residual_j,
            },
        },
        ObserverCouplingError::ConflictingConductance(_)
        | ObserverCouplingError::ConflictingReservoirConductance(_)
        | ObserverCouplingError::ConflictingAmbientConductance(_) => {
            ShellThermoplasticError::Invalid("conflicting thermal conductance owners")
        }
    }
}
