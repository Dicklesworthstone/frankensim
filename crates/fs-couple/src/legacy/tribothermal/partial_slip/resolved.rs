//! Conduction-derived heat partition for the original loaded partial-slip owner.
//!
//! Only irreversible microslip work enters the common, massless temperature
//! node. Reversible -Ti*dsi exchange stays with each original caloric body,
//! using the SAME two-temperature affine contact free energy. Ordinary signed
//! body-to-body conduction is distinguished from generated heat. Numerical
//! return-map loss remains numerical loss and never heats a cell.
//!
//! Original material friction queries retain their declared body-temperature,
//! load and actual translational-speed semantics. A common interface temperature
//! is not silently substituted for a source's separately declared surface axes.
//! Conductivity and other sources see the common node and the whole endpoint.
//! Fixed orthogonal finite-volume geometry and negligible interface resistance
//! have the same scope as the original `FrictionHeatRoute`.

use super::*;
use crate::friction_heat::{
    FrictionHeatInput, FrictionHeatRoute, FrictionHeatStep, FrictionHeatTolerances,
};

/// Same-endpoint body conductivities and independent thermal transport.
pub struct ResolvedCaloricContactSource<S> {
    pub conductivity_w_m_k: [f64; 2],
    /// OTHER cell sources only; neither friction nor reversible contact heat
    /// belongs here. The route's direct link has exclusive transport ownership.
    pub transport: StepCoefficients<S>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ResolvedCaloricContactSolve {
    pub contact: CaloricContactSolve,
    /// Numerical interface-temperature coordinate scale, J/K; never storage.
    pub temperature_scale_j_per_k: f64,
    pub interface_energy_tolerance_j: f64,
}

pub struct ResolvedCaloricContactReport<S> {
    /// Original contact/internal-energy account. Its physical heat allocation
    /// and any material sample's partition are derived from this same network.
    pub caloric: CaloricContactReport<S>,
    pub interface_temperature_k: f64,
    pub conductances_w_k: [f64; 2],
    /// Net node heat into each body, including ordinary unequal-T exchange.
    pub body_heat_j: [f64; 2],
    /// Ordinary signed body-to-body transfer, excluding generated microslip heat.
    pub conduction_heat_j: [f64; 2],
    pub interface_energy_residual_j: f64,
}

pub struct ResolvedCaloricContactCandidate<'a, 't, S> {
    thermal: DepositionSimulation<'t>,
    contact: CaloricPartialSlip<'a>,
    report: ResolvedCaloricContactReport<S>,
}
impl<'a, 't, S> ResolvedCaloricContactCandidate<'a, 't, S> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn contact(&self) -> &CaloricPartialSlip<'a> {
        &self.contact
    }
    pub const fn report(&self) -> &ResolvedCaloricContactReport<S> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        CaloricPartialSlip<'a>,
        ResolvedCaloricContactReport<S>,
    ) {
        (self.thermal, self.contact, self.report)
    }
}

impl<'a> CaloricPartialSlip<'a> {
    /// Bind an original loaded checkpoint and caloric objects to an existing
    /// physical interface. There is no prescribed fraction in this constructor.
    /// `curves` follow the existing thermal link's A/B order, which must be the
    /// original ordered mechanical surfaces. Reference conductivities admit
    /// the link's existing half-cell series conductance before any motion.
    pub fn new_with_resolved_interface(
        thermal: &DepositionSimulation<'_>,
        model: AffineContactStiffness,
        checkpoint: PartialSlipCheckpoint,
        curves: [&'a EquilibriumEnthalpyPhaseCurve; 2],
        link_id: u64,
        reference_conductivity_w_m_k: [f64; 2],
    ) -> Result<Self, CaloricContactError> {
        let route = FrictionHeatRoute::from_link(thermal, link_id, reference_conductivity_w_m_k)?;
        let cells = route.cells();
        Self::from_checkpoint(
            thermal,
            model,
            checkpoint,
            curves,
            cells,
            CaloricHeatRouting::Resolved(route),
        )
    }

    /// Switch only heat routing on an accepted original contact. Loaded strain,
    /// branches, clock, free energy and pinned material claims remain intact.
    /// This also upgrades existing material-bound contacts without reconstructing
    /// their history. Subsequent stepping must use `step_resolved_candidate`.
    pub fn with_resolved_interface(
        mut self,
        thermal: &DepositionSimulation<'_>,
        link_id: u64,
        reference_conductivity_w_m_k: [f64; 2],
    ) -> Result<Self, CaloricContactError> {
        let route = FrictionHeatRoute::from_link(thermal, link_id, reference_conductivity_w_m_k)?;
        if route.cells() != self.cells
            || thermal.time_s() != self.time_s
            || temperatures(thermal, self.curves, self.cells)? != self.temperature_k
        {
            return Err(CaloricContactError::Invalid(
                "resolved route must retain accepted ordered cells, temperatures and clock",
            ));
        }
        self.routing = CaloricHeatRouting::Resolved(route);
        Ok(self)
    }

    pub fn resolved_heat_route(&self) -> Option<&FrictionHeatRoute> {
        match &self.routing {
            CaloricHeatRouting::Resolved(route) => Some(route),
            CaloricHeatRouting::Prescribed(_) => None,
        }
    }

    /// Joint native contact/enthalpy solve with one common interface-temperature
    /// coordinate. Every trial starts from the accepted mechanical checkpoint.
    /// The returned contact and thermal owners must be published together.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_resolved_candidate<'t, S, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        frame: TangentFrame,
        kinematics: PartialSlipKinematics,
        ownership: &GeneralizedWorkOwnership,
        controls: ResolvedCaloricContactSolve,
        guess_j: &[f64],
        temperature_guess_k: f64,
        resolve: impl Fn(&DepositionSimulation<'t>, f64) -> Result<ResolvedCaloricContactSource<S>, E>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ResolvedCaloricContactCandidate<'a, 't, S>, CaloricContactError<E>> {
        let CaloricHeatRouting::Resolved(route) = &self.routing else {
            return Err(CaloricContactError::Invalid(
                "contact has no admitted resolved heat route",
            ));
        };
        let cancel = RefCell::new(cancelled);
        let poll = || {
            if (*cancel.borrow_mut())() {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        };
        poll()?;
        let contact_controls = controls.contact;
        contact_controls
            .newton
            .validate()
            .map_err(CaloricContactError::Newton)?;
        if thermal.time_s() != self.time_s
            || temperatures(thermal, self.curves, self.cells)? != self.temperature_k
            || contact_controls.maximum_cells == 0
            || thermal.cells().len() > contact_controls.maximum_cells
            || guess_j.len() != thermal.cells().len()
            || guess_j.iter().any(|v| !v.is_finite())
            || contact_controls.maximum_contact_trials == 0
            || contact_controls.maximum_newton_iterations == 0
            || !controls.temperature_scale_j_per_k.is_finite()
            || controls.temperature_scale_j_per_k <= 0.0
            || !temperature_guess_k.is_finite()
            || temperature_guess_k <= 0.0
            || [
                contact_controls.cell_energy_tolerance_j,
                contact_controls.thermal_energy_tolerance_j,
                contact_controls.total_energy_tolerance_j,
                controls.interface_energy_tolerance_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(CaloricContactError::Invalid(
                "incompatible accepted contact, thermal coordinates or resolved-interface controls",
            ));
        }
        let step = FrictionHeatStep::new(thermal, target_s, route, poll)?;
        let dt = step.duration_s();
        if kinematics.dt_s != dt {
            return Err(CaloricContactError::Invalid(
                "contact and caloric physical intervals must match",
            ));
        }
        if self.material_friction.is_some() && kinematics.torsional_spin_rad_per_s != 0.0 {
            return Err(CaloricContactError::Invalid(
                "material friction has no torsional-speed mapping",
            ));
        }
        let count = step.cell_ids().len();
        let trials = Cell::new(0usize);
        let evaluate = |x: &[f64]| {
            poll()?;
            if x.len() != count + 1 {
                return Err(CaloricContactError::Invalid(
                    "incomplete resolved-interface iterate",
                ));
            }
            if trials.get() >= contact_controls.maximum_contact_trials {
                return Err(CaloricContactError::TrialBudget);
            }
            trials.set(trials.get() + 1);
            let interface_temperature =
                temperature_guess_k + x[count] / controls.temperature_scale_j_per_k;
            step.evaluate(
                &x[..count],
                interface_temperature,
                |state, interface_temperature| {
                    let mechanical = self.mechanical_endpoint(
                        state,
                        frame,
                        kinematics,
                        ownership,
                        contact_controls,
                    )?;
                    let ResolvedCaloricContactSource {
                        conductivity_w_m_k,
                        transport,
                    } = resolve(state, interface_temperature)
                        .map_err(CaloricContactError::Source)?;
                    let StepCoefficients {
                        mut sources_w,
                        conductances_w_k,
                        payload,
                    } = transport;
                    // Evaluate outside supply BEFORE adding internal reversible heat.
                    let external_heat_j = finite(
                        dt * state.external_step_power_w(&sources_w, &conductances_w_k, poll)?,
                    )?;
                    for (id, heat) in self.cells.into_iter().zip(mechanical.reversible_heat_j) {
                        let power = finite(heat / dt)?;
                        if heat != 0.0 && power == 0.0 {
                            return Err(CaloricContactError::Invalid(
                                "reversible contact heat vanished in power conversion",
                            ));
                        }
                        let total = sources_w.entry(id).or_insert(0.0);
                        *total = finite(*total + power)?;
                    }
                    Ok(FrictionHeatInput {
                        dissipated_j: mechanical.contact.dissipation.heat_j,
                        conductivity_w_m_k,
                        transport: StepCoefficients {
                            sources_w,
                            conductances_w_k,
                            payload: (mechanical, external_heat_j, payload),
                        },
                    })
                },
                poll,
            )
        };
        let problem = FallibleFiniteDifference::new(
            count + 1,
            contact_controls.difference_step_j,
            |x: &[f64]| evaluate(x).map(|t| t.residual_j().to_vec()),
        )
        .map_err(difference_error)?;
        let mut initial = guess_j.to_vec();
        initial.push(0.0);
        let construction = NewtonKrylovState::new(&problem, initial, contact_controls.newton);
        if let Some(e) = problem.take_failure() {
            return Err(difference_error(e));
        }
        let mut solver = construction.map_err(CaloricContactError::Newton)?;
        let result = solver.run_controlled(
            &problem,
            contact_controls.maximum_newton_iterations,
            &mut || (*cancel.borrow_mut())() || problem.has_failed(),
        );
        if let Some(e) = problem.take_failure() {
            return Err(difference_error(e));
        }
        poll()?;
        if result.progress != SolverRunProgress::Complete {
            return Err(CaloricContactError::Interrupted(result.progress));
        }
        if !result.report.converged {
            return Err(CaloricContactError::NotConverged(result.report));
        }
        let trial = evaluate(&solver.x)?;
        let maximum_cell_residual_j = trial.residual_j()[..count]
            .iter()
            .copied()
            .map(f64::abs)
            .fold(0.0, f64::max);
        let (next, (mechanical, external_heat_j, source), heat) = trial.accept(
            FrictionHeatTolerances {
                cell_j: contact_controls.cell_energy_tolerance_j,
                thermal_j: contact_controls.thermal_energy_tolerance_j,
                interface_j: controls.interface_energy_tolerance_j,
                total_j: contact_controls.total_energy_tolerance_j,
            },
            poll,
        )?;
        let MechanicalEndpoint {
            contact,
            mut material,
            thermodynamics,
            temperature_k,
            internal_change,
            reversible_heat_j,
        } = mechanical;
        let conductance_sum = finite(heat.conductances_w_k[0] + heat.conductances_w_k[1])?;
        let fraction_a = finite(heat.conductances_w_k[0] / conductance_sum)?;
        let derived_partition =
            HeatPartition::new(fraction_a, 1.0 - fraction_a, 0.0).map_err(|_| {
                CaloricContactError::Invalid("unrepresentable conduction-derived heat share")
            })?;
        if let Some(sample) = &mut material {
            sample.partition = derived_partition;
        }
        let physical_a = finite(contact.dissipation.heat_j * fraction_a)?;
        let physical_heat_j = [physical_a, finite(contact.dissipation.heat_j - physical_a)?];
        let deposited_heat_j = [
            finite(physical_heat_j[0] + reversible_heat_j[0])?,
            finite(physical_heat_j[1] + reversible_heat_j[1])?,
        ];
        let conduction_heat_j = [
            finite(heat.body_heat_j[0] - physical_heat_j[0])?,
            finite(heat.body_heat_j[1] - physical_heat_j[1])?,
        ];
        let total_energy_residual_j = finite(
            heat.thermal.stored_energy_change_j
                + internal_change
                + contact.generalized_work.numerical_dissipation_j
                - contact.generalized_work.work_into_interface_j
                - external_heat_j,
        )?;
        if total_energy_residual_j.abs() > contact_controls.total_energy_tolerance_j {
            return Err(CaloricContactError::Invalid(
                "resolved contact internal energy and finite caloric storage do not balance",
            ));
        }
        let mut next_contact = self.clone();
        next_contact.checkpoint = contact.checkpoint.clone();
        next_contact.temperature_k = temperature_k;
        next_contact.thermodynamics = thermodynamics;
        next_contact.time_s = target_s;
        if let (Some(binding), Some(sample)) = (&mut next_contact.material_friction, &material) {
            binding.distance_m = sample.accumulated_distance_m;
            binding.time_s = target_s;
        }
        poll()?;
        Ok(ResolvedCaloricContactCandidate {
            thermal: next,
            contact: next_contact,
            report: ResolvedCaloricContactReport {
                caloric: CaloricContactReport {
                    contact,
                    material,
                    thermodynamics,
                    contact_internal_energy_change_j: internal_change,
                    reversible_heat_j,
                    physical_heat_j,
                    deposited_heat_j,
                    thermal: heat.thermal,
                    external_heat_j,
                    source,
                    nonlinear: result.report,
                    contact_trials: trials.get(),
                    maximum_cell_residual_j,
                    total_energy_residual_j,
                },
                interface_temperature_k: heat.interface_temperature_k,
                conductances_w_k: heat.conductances_w_k,
                body_heat_j: heat.body_heat_j,
                conduction_heat_j,
                interface_energy_residual_j: heat.interface_residual_j,
            },
        })
    }
}
