//! Source-bound axisymmetric current and the original finite thermal history.
//!
//! Every trial refreshes the original conductivity at the same temperature and
//! declared liquid fraction as its caloric receiver. The native meridional field
//! integrates E.J over the actual swept cells, including the angular span.
//! Those cell powers enter the existing backward-Euler enthalpy residual once.
//! Supplied electrical work and thermal energy close before joint publication.
//!
//! Thermal links retain their authored areas and material half-paths; this
//! bridge does not infer a radial thermal mesh or add a second heat integrator.
//! Geometry and material frames remain fixed. Electrical contacts require an
//! explicit source-refresh and receiving-surface partition and are refused here.

use super::super::{
    CircuitEnthalpyReport, CircuitEnthalpySolve, DrivenCoefficients, DriverAccount, ThermalDriver,
    cell_driver_candidate, finite,
    material::{MaterialNetwork, MaterialNetworkSample},
};
use super::{
    ConductionDrive, ConductionThermalNetwork, ConductionThermalStep, Error, check, field_error,
    invalid, material, same_material,
};
use fs_em::conduction::{
    ConductionControls,
    axisymmetric::{AxisymmetricConduction, AxisymmetricConductionField},
};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_material::conductor::{ResolvedConductorState, planar::PlanarConductorMaterial};
use fs_thermal::enthalpy::DepositionSimulation;
use std::cell::RefCell;

impl ConductionDrive<'_> {
    fn solve_axisymmetric<'p, C: PlanarConductorMaterial>(
        self,
        problem: &'p AxisymmetricConduction<C>,
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricConductionField<'p, C>, Error> {
        match self {
            Self::Voltage(boundary) => problem
                .solve(boundary, controls, cancelled)
                .map_err(field_error),
            Self::CurrentDensity { voltage, currents } => problem
                .solve_with_current_traces(voltage, currents, controls, cancelled)
                .map_err(field_error),
            Self::Current {
                positive,
                negative,
                current_a,
            } => problem
                .resistance(positive, negative, controls, cancelled)
                .and_then(|response| {
                    response.field_at_current(
                        current_a,
                        controls.current_tolerance_a,
                        controls.power_tolerance_w,
                        cancelled,
                    )
                })
                .map_err(field_error),
            Self::Power {
                positive,
                negative,
                source,
            } => problem
                .resistance(positive, negative, controls, cancelled)
                .and_then(|response| {
                    response.field_at_power(
                        source,
                        controls.current_tolerance_a,
                        controls.power_tolerance_w,
                        cancelled,
                    )
                })
                .map_err(field_error),
        }
    }
}

struct FieldCandidate<C> {
    problem: AxisymmetricConduction<C>,
    step: ConductionThermalStep,
}
impl<C: PlanarConductorMaterial> ThermalDriver for FieldCandidate<C> {
    type Model = AxisymmetricConduction<C>;
    type Step = ConductionThermalStep;

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

/// Jointly accepted thermal history and the current endpoint material operator.
/// Electrical observations use local radial and axial components; watts
/// and joules refer to the original angular span. Consume both owners together
/// to continue; a failed trial changes neither supplied input.
pub struct AxisymmetricConductionEnthalpyCandidate<'t, C = ResolvedConductorState> {
    thermal: DepositionSimulation<'t>,
    problem: AxisymmetricConduction<C>,
    report: CircuitEnthalpyReport<MaterialNetworkSample, ConductionThermalStep>,
}
impl<'t, C: PlanarConductorMaterial> AxisymmetricConductionEnthalpyCandidate<'t, C> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }

    pub fn problem(&self) -> &AxisymmetricConduction<C> {
        &self.problem
    }

    pub fn report(&self) -> &CircuitEnthalpyReport<MaterialNetworkSample, ConductionThermalStep> {
        &self.report
    }

    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        AxisymmetricConduction<C>,
        CircuitEnthalpyReport<MaterialNetworkSample, ConductionThermalStep>,
    ) {
        (self.thermal, self.problem, self.report)
    }
}

impl ConductionThermalNetwork<'_> {
    /// Advance finite thermal cells and a source-bound axisymmetric DC field
    /// to one common endpoint under voltage, total-current, power or distributed
    /// inward-current-density drive. All unlisted electrical surfaces insulate.
    ///
    /// Original ConductionRegion triangle indices route the actual swept cell
    /// volumes to complete thermal receivers. Several triangles may share one
    /// receiver; passive thermal cells and source-resolved Fourier links remain
    /// supported. No planar thickness or extra 2*pi factor is inserted.
    ///
    /// Complete cards, source selections, fixed coordinates and meridional
    /// frames stay pinned. Temperature and an explicitly owned liquid fraction
    /// evolve from the SAME original caloric chart. Density must match the
    /// receiving volume throughout that chart; no expansion or phase mixing
    /// rule is inferred. Contact-bearing electrical owners refuse because this
    /// entry point has no receiving-surface heat partition or contact refresh.
    ///
    /// Every residual evaluation rebuilds the native electrical operator and
    /// solves the actual drive. maximum_circuit_trials bounds these field
    /// evaluations, including the final fresh endpoint; controls.circuit is
    /// unused because no circuit storage is introduced. Original material,
    /// enthalpy, assembly, electrical balance and cancellation limits apply.
    /// The returned interval is first order in time; spatial/temporal accuracy
    /// remains the caller's refinement decision.
    #[allow(clippy::too_many_arguments)]
    pub fn step_axisymmetric_candidate<'t, C: PlanarConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &AxisymmetricConduction<C>,
        drive: ConductionDrive<'_>,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<AxisymmetricConductionEnthalpyCandidate<'t, C>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        if !problem.contacts().is_empty() {
            return Err(invalid(
                "axisymmetric contact heat requires explicit surface partitions and source refresh",
            ));
        }
        let dt = self.admit_mesh_step(
            thermal,
            target_s,
            problem.mesh(),
            assembly,
            controls,
            initial_energy_increments_j,
        )?;
        let routes = self.admit_mesh_routes(thermal, problem.mesh(), &cancel)?;
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        };
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        let accepted = prepared
            .resolve(thermal, || (*cancel.borrow_mut())())
            .map_err(material)?;
        let initial =
            self.materials_from_sources(problem.conductors(), &accepted.source, &routes, &cancel)?;
        if initial.iter().zip(problem.conductors()).any(|(a, b)| {
            !a.has_same_electrical_source(b)
                || a.frame() != b.frame()
                || !same_material(a.material(), b.material())
        }) {
            return Err(invalid(
                "accepted axisymmetric conductor sources must match the actual thermal state",
            ));
        }
        let (thermal, problem, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |state, _remaining| {
                check(&cancel)?;
                let coefficients = prepared
                    .resolve(state, || (*cancel.borrow_mut())())
                    .map_err(material)?;
                let sources = self.materials_from_sources(
                    problem.conductors(),
                    &coefficients.source,
                    &routes,
                    &cancel,
                )?;
                let next = AxisymmetricConduction::new(
                    problem.mesh().clone(),
                    sources,
                    assembly,
                    &mut || (*cancel.borrow_mut())(),
                )
                .map_err(field_error)?;
                let field = drive
                    .solve_axisymmetric(&next, field_controls, &mut || (*cancel.borrow_mut())())?;
                let heat_j = finite(dt * field.power_w())?;
                let supplied_work_j = finite(
                    dt * drive.supplied_power_w(
                        field.voltage_v(),
                        field.power_w(),
                        field.power_balance_residual_w(),
                    )?,
                )?;
                if (field.power_w() > 0.0 && heat_j == 0.0)
                    || (heat_j - supplied_work_j).abs() > controls.total_energy_tolerance_j
                {
                    return Err(invalid(
                        "axisymmetric Joule heat and terminal source work do not close",
                    ));
                }
                let step = ConductionThermalStep {
                    start_s: thermal.time_s(),
                    end_s: target_s,
                    voltage_v: field.voltage_v().to_vec(),
                    electric_field_v_m: field.electric_field_v_m().to_vec(),
                    current_density_a_m2: field.current_density_a_m2().to_vec(),
                    cell_power_w: field.cell_power_w().to_vec(),
                    contacts: Vec::new(),
                    nodal_boundary_current_a: field.nodal_boundary_current_a().to_vec(),
                    nodal_applied_current_a: field.nodal_applied_current_a().to_vec(),
                    boundary_current_a: field.boundary_current_a(),
                    maximum_free_current_residual_a: field.maximum_free_current_residual_a(),
                    power_balance_residual_w: field.power_balance_residual_w(),
                    heat_j,
                    supplied_work_j,
                    linear: field.report().cloned(),
                };
                Ok(DrivenCoefficients {
                    circuit: FieldCandidate {
                        problem: next,
                        step,
                    },
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: coefficients.source,
                    circuit_trials: 1,
                })
            },
            |candidate, _dt, checkpoint| {
                self.deposit_power(&candidate.step.cell_power_w, &routes, checkpoint)
            },
            || (*cancel.borrow_mut())(),
        )?;
        check(&cancel)?;
        Ok(AxisymmetricConductionEnthalpyCandidate {
            thermal,
            problem,
            report,
        })
    }
}
