//! Source-bound nonlinear meridional current and finite thermal feedback.
//!
//! The native electrical response closes its actual E/J coordinates at each
//! thermal trial. Exact swept bulk and contact losses enter the existing
//! enthalpy residual once; constitutive iteration error is never heat.
//! Geometry, source pins and caloric charts remain with their original owners.

mod heat;
pub use heat::{
    AxisymmetricConductionResponseDrive, step_axisymmetric_conduction_response_candidate,
};

use super::super::{
    ConductionDrive, ConductionThermalError, ConductionThermalNetwork, Error, check, contacts,
    invalid, map_conduction_error, material,
    response::{ConductionResponseThermalStep, terminals::Drive},
    same_material,
};
use crate::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, CircuitThermalCoefficients,
    DrivenCoefficients, DriverAccount, ThermalDriver, cell_driver_candidate, finite,
    material::{MaterialNetwork, MaterialNetworkSample},
};
use fs_em::conduction::{
    ConductionControls,
    axisymmetric::response::AxisymmetricConductionResponse,
    response::{
        ConductionResponseControls, ConductionResponseError, MAX_RESPONSE_ITERATIONS,
        terminals::ConductionTerminalDrive,
    },
};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_thermal::enthalpy::DepositionSimulation;
use std::cell::RefCell;

struct ResponseCandidate {
    problem: AxisymmetricConductionResponse,
    step: ConductionResponseThermalStep,
}
impl ThermalDriver for ResponseCandidate {
    type Model = AxisymmetricConductionResponse;
    type Step = ConductionResponseThermalStep;

    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.step.heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: 0.0,
            supplied_j: self.step.supplied_work_j,
        }
    }

    fn into_parts(self) -> (Self::Model, Self::Step) {
        (self.problem, self.step)
    }
}

/// Jointly accepted finite enthalpy and nonlinear meridional electrical state.
/// Consume both together to continue; rejected trials publish neither owner.
pub struct AxisymmetricConductionResponseEnthalpyCandidate<'t, S = MaterialNetworkSample> {
    thermal: DepositionSimulation<'t>,
    problem: AxisymmetricConductionResponse,
    report: CircuitEnthalpyReport<S, ConductionResponseThermalStep>,
}
impl<'t, S> AxisymmetricConductionResponseEnthalpyCandidate<'t, S> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }

    pub fn problem(&self) -> &AxisymmetricConductionResponse {
        &self.problem
    }

    pub fn report(&self) -> &CircuitEnthalpyReport<S, ConductionResponseThermalStep> {
        &self.report
    }

    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        AxisymmetricConductionResponse,
        CircuitEnthalpyReport<S, ConductionResponseThermalStep>,
    ) {
        (self.thermal, self.problem, self.report)
    }
}

impl ConductionThermalNetwork<'_> {
    /// Advance scalar sigma/rho(T, phase, |E|, |J|), native swept current and
    /// finite enthalpy to one backward-Euler endpoint. The existing voltage,
    /// total-current, regulated-power and distributed-current-density drives
    /// are solved at their actual amplitudes; nonlinear fields are not scaled.
    ///
    /// Original complete swept-cell routes, fixed geometry and constant-density
    /// caloric charts remain required. Source selections, fixed coordinates and
    /// phase axes stay pinned. E/J-dependent caloric or Fourier laws refuse:
    /// this memoryless conduction adapter adds no reversible field storage.
    ///
    /// Every native field solve, including constitutive and power iterations
    /// and final fresh thermal admission, consumes maximum_circuit_trials.
    /// Failure or cancellation leaves both supplied owners reusable. Contacts
    /// require the explicit contact entry and its receiving heat partitions.
    #[allow(clippy::too_many_arguments)]
    pub fn step_axisymmetric_response_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &AxisymmetricConductionResponse,
        drive: ConductionDrive<'_>,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        response: ConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<AxisymmetricConductionResponseEnthalpyCandidate<'t>, Error> {
        self.step_axisymmetric_response_contact_candidate(
            thermal,
            target_s,
            problem,
            drive,
            &[],
            assembly,
            field_controls,
            response,
            controls,
            initial_energy_increments_j,
            cancelled,
        )
    }

    /// Close native nonlinear bulk/interface response and both receiving bodies'
    /// heat through the same residual. Each contact retains its original card,
    /// ordered material pair, fixed pressure/history and declared phase axes.
    /// The native owner alone closes any explicitly owned electrical interface
    /// coordinates against its solved traces; no interface law is inferred.
    ///
    /// Exact swept-surface loss is additional to bulk cell heat. The existing
    /// authored two-sided fractions route it once to the original adjacent
    /// receivers. Original material Fourier links remain active independently.
    /// Uniqueness, nonlinear convergence and continuum accuracy are not implied
    /// by acceptance; time/mesh refinement remains the caller's responsibility.
    #[allow(clippy::too_many_arguments)]
    pub fn step_axisymmetric_response_contact_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &AxisymmetricConductionResponse,
        drive: ConductionDrive<'_>,
        contact_bindings: &[contacts::ConductionContactThermal<'_>],
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        response: ConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<AxisymmetricConductionResponseEnthalpyCandidate<'t>, Error> {
        self.step_axisymmetric_response_driven_candidate(
            thermal,
            target_s,
            problem,
            Drive::Field(drive),
            contact_bindings,
            assembly,
            field_controls,
            response,
            controls,
            initial_energy_increments_j,
            cancelled,
        )
    }

    /// Advance ideal voltage/current/mixed electrodes, nonlinear bulk and
    /// declared interface laws, and finite enthalpy to one common endpoint.
    /// A floating equipotential electrode can transfer current between faces
    /// while its total current remains zero. Original source pins, explicit
    /// contact heat partitions, Fourier links and phase ownership remain bound.
    ///
    /// Every actual grouped-electrode field solve consumes the original
    /// thermal circuit budget, including discarded nonlinear iterates and
    /// final fresh admission. The accepted report retains the ordered terminal
    /// operating point and independent source work. Failure publishes neither
    /// owner; consume both accepted owners together to continue.
    #[allow(clippy::too_many_arguments)]
    pub fn step_axisymmetric_response_terminal_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &AxisymmetricConductionResponse,
        drive: ConductionTerminalDrive<'_>,
        contact_bindings: &[contacts::ConductionContactThermal<'_>],
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        response: ConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<AxisymmetricConductionResponseEnthalpyCandidate<'t>, Error> {
        self.step_axisymmetric_response_driven_candidate(
            thermal,
            target_s,
            problem,
            Drive::Terminals(drive),
            contact_bindings,
            assembly,
            field_controls,
            response,
            controls,
            initial_energy_increments_j,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_axisymmetric_response_driven_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &AxisymmetricConductionResponse,
        drive: Drive<'_>,
        contact_bindings: &[contacts::ConductionContactThermal<'_>],
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        response: ConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<AxisymmetricConductionResponseEnthalpyCandidate<'t>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        };
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        self.step_axisymmetric_response_resolved_candidate(
            thermal,
            target_s,
            problem,
            drive,
            contact_bindings,
            assembly,
            field_controls,
            response,
            controls,
            initial_energy_increments_j,
            |state| {
                prepared
                    .resolve(state, || (*cancel.borrow_mut())())
                    .map_err(material)
            },
            |sample| sample,
            &cancel,
        )
    }

    // The original numerical driver accepts a complete thermal resolver; all
    // mechanisms see one proposed enthalpy state and retain their source type.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::circuit_enthalpy::conduction) fn step_axisymmetric_response_resolved_candidate<
        't,
        S,
        E: From<ConductionThermalError>,
    >(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &AxisymmetricConductionResponse,
        drive: Drive<'_>,
        contact_bindings: &[contacts::ConductionContactThermal<'_>],
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        response: ConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        resolve_sources: impl Fn(
            &DepositionSimulation<'t>,
        )
            -> Result<CircuitThermalCoefficients<S>, CircuitEnthalpyError<E>>,
        sample: impl Fn(&S) -> &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<AxisymmetricConductionResponseEnthalpyCandidate<'t, S>, CircuitEnthalpyError<E>>
    {
        check(cancel).map_err(map_conduction_error)?;
        // Admit the authored cap before restricting it by the outer budget.
        if !(1..=MAX_RESPONSE_ITERATIONS).contains(&response.maximum_iterations) {
            return Err(map_conduction_error(invalid(
                "invalid bounded axisymmetric response iteration limit",
            )));
        }
        let dt = self
            .admit_mesh_step(
                thermal,
                target_s,
                problem.mesh(),
                assembly,
                controls,
                initial_energy_increments_j,
            )
            .map_err(map_conduction_error)?;
        let routes = self
            .admit_mesh_routes(thermal, problem.mesh(), cancel)
            .map_err(map_conduction_error)?;
        self.admit_response_caloric_sources(problem.sources(), cancel)
            .map_err(map_conduction_error)?;
        let resolve = |state: &DepositionSimulation<'t>| {
            let coefficients = resolve_sources(state)?;
            if !coefficients.resistors.is_empty() {
                return Err(map_conduction_error(invalid(
                    "axisymmetric response cannot omit graph electrical work",
                )));
            }
            self.admit_response_heat(sample(&coefficients.source), cancel)
                .map_err(map_conduction_error)?;
            Ok::<_, CircuitEnthalpyError<E>>(coefficients)
        };
        let accepted = resolve(thermal)?;
        let contact_axes = contacts::response_axes(problem.contact_response_sources());
        let contacts = contacts::PreparedContacts::new_for_response_contacts(
            self,
            problem.mesh(),
            problem.contacts(),
            problem.contact_liquid_mass_fraction_axes(),
            contact_bindings,
            &routes,
            sample(&accepted.source),
            &contact_axes,
            cancel,
        )
        .map_err(map_conduction_error)?;
        let materials_at = |sample: &MaterialNetworkSample| {
            self.response_materials_from(
                problem.conductors(),
                problem.sources(),
                sample,
                &routes,
                cancel,
            )
            .map_err(map_conduction_error)
        };
        let initial = materials_at(sample(&accepted.source))?;
        if initial.iter().zip(problem.conductors()).any(|(a, b)| {
            !a.has_same_electrical_source(b) || !same_material(a.material(), b.material())
        }) {
            return Err(map_conduction_error(invalid(
                "accepted axisymmetric response sources must match the actual thermal state",
            )));
        }
        let (thermal, problem, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |state, remaining| {
                check(cancel).map_err(map_conduction_error)?;
                if remaining == 0 {
                    return Err(CircuitEnthalpyError::CircuitBudget);
                }
                let coefficients = resolve(state)?;
                let sources = materials_at(sample(&coefficients.source))?;
                let contact_sources = contacts
                    .refresh_contacts(problem.contacts(), sample(&coefficients.source), cancel)
                    .map_err(map_conduction_error)?;
                let bounded = ConductionResponseControls {
                    maximum_iterations: response.maximum_iterations.min(remaining),
                    ..response
                };
                let (candidate, terminals) = match drive {
                    Drive::Field(field_drive) => (match field_drive {
                        ConductionDrive::Voltage(voltage) => problem
                            .solve_with_materials_and_contacts_candidate(
                                sources,
                                contact_sources,
                                voltage,
                                &[],
                                assembly,
                                field_controls,
                                bounded,
                                &mut || (*cancel.borrow_mut())(),
                            ),
                        ConductionDrive::CurrentDensity { voltage, currents } => problem
                            .solve_with_materials_and_contacts_candidate(
                                sources,
                                contact_sources,
                                voltage,
                                currents,
                                assembly,
                                field_controls,
                                bounded,
                                &mut || (*cancel.borrow_mut())(),
                            ),
                        ConductionDrive::Current {
                            positive,
                            negative,
                            current_a,
                        } => problem.solve_current_with_materials_and_contacts_candidate(
                            sources,
                            contact_sources,
                            positive,
                            negative,
                            current_a,
                            assembly,
                            field_controls,
                            bounded,
                            &mut || (*cancel.borrow_mut())(),
                        ),
                        ConductionDrive::Power {
                            positive,
                            negative,
                            source,
                        } => problem.solve_power_with_materials_and_contacts_candidate(
                            sources,
                            contact_sources,
                            positive,
                            negative,
                            source,
                            assembly,
                            field_controls,
                            bounded,
                            &mut || (*cancel.borrow_mut())(),
                        ),
                    })
                    .map(|candidate| (candidate, None)),
                    Drive::Terminals(terminal_drive) => problem
                        .solve_terminals_with_materials_and_contacts_candidate(
                            sources,
                            contact_sources,
                            terminal_drive,
                            assembly,
                            field_controls,
                            bounded,
                            remaining,
                            &mut || (*cancel.borrow_mut())(),
                        )
                        .map(|(candidate, terminals)| (candidate, Some(terminals))),
                }
                .map_err(|error| match error {
                    ConductionResponseError::FieldBudget { .. } => {
                        CircuitEnthalpyError::CircuitBudget
                    }
                    ConductionResponseError::NotConverged { .. }
                        if bounded.maximum_iterations < response.maximum_iterations =>
                    {
                        CircuitEnthalpyError::CircuitBudget
                    }
                    other => map_conduction_error(Error::from(other)),
                })?;
                let (next, field) = candidate.into_parts();
                let heat_j = finite(dt * field.power_w)?;
                let supplied_work_j = finite(
                    dt * drive
                        .supplied_power_w(&field, terminals.as_ref())
                        .map_err(map_conduction_error)?,
                )?;
                if (field.power_w > 0.0 && heat_j == 0.0)
                    || (heat_j - supplied_work_j).abs() > controls.total_energy_tolerance_j
                {
                    return Err(map_conduction_error(invalid(
                        "axisymmetric nonlinear Joule heat and terminal source work do not close",
                    )));
                }
                let circuit_trials = field.field_solves;
                let contact_heat = contacts
                    .heat(&field.contact_response)
                    .map_err(map_conduction_error)?;
                Ok(DrivenCoefficients {
                    circuit: ResponseCandidate {
                        problem: next,
                        step: ConductionResponseThermalStep {
                            start_s: thermal.time_s(),
                            end_s: target_s,
                            field,
                            terminals,
                            contacts: contact_heat,
                            heat_j,
                            supplied_work_j,
                        },
                    },
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: coefficients.source,
                    circuit_trials,
                })
            },
            |candidate, _dt, checkpoint| {
                let mut powers = self
                    .deposit_power(&candidate.step.field.cell_power_w, &routes, checkpoint)
                    .map_err(map_conduction_error)?;
                contacts::deposit(&candidate.step.contacts, &mut powers, checkpoint)
                    .map_err(map_conduction_error)?;
                Ok(powers)
            },
            || (*cancel.borrow_mut())(),
        )?;
        check(cancel).map_err(map_conduction_error)?;
        Ok(AxisymmetricConductionResponseEnthalpyCandidate {
            thermal,
            problem,
            report,
        })
    }
}
