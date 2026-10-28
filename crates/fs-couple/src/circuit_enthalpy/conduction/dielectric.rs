//! Native dielectric storage/current and finite enthalpy share one endpoint.
//!
//! The existing conduction material/volume binding updates every common source
//! bundle at actual (T, liquid mass fraction). The retained electroquasistatic
//! owner advances charge and deposits its ORIGINAL local midpoint E.J. It
//! admits only validity-wide constant permittivity, so all changing electrical
//! storage is represented by the native potential history, never hidden
//! parameter work. No independent charge or heat integrator is introduced.

use super::*;
use fs_em::electrostatic::relaxation::{
    DielectricRelaxationCandidate, DielectricRelaxationControls, DielectricRelaxationError,
    DielectricRelaxationStep, PlanarLeakyDielectric,
};

/// Source-backed reversible dielectric thermal work on the same charge owner.
pub mod thermal;

/// One jointly accepted charge/potential history, material state and enthalpy.
pub struct DielectricEnthalpyCandidate<'t, C = ResolvedConductorState> {
    thermal: DepositionSimulation<'t>,
    field: PlanarLeakyDielectric<C>,
    report: CircuitEnthalpyReport<MaterialNetworkSample, DielectricRelaxationStep>,
}
impl<'t, C: PlanarConductorMaterial> DielectricEnthalpyCandidate<'t, C> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn field(&self) -> &PlanarLeakyDielectric<C> {
        &self.field
    }
    pub fn report(
        &self,
    ) -> &CircuitEnthalpyReport<MaterialNetworkSample, DielectricRelaxationStep> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        PlanarLeakyDielectric<C>,
        CircuitEnthalpyReport<MaterialNetworkSample, DielectricRelaxationStep>,
    ) {
        (self.thermal, self.field, self.report)
    }
}

impl<C: PlanarConductorMaterial> ThermalDriver for DielectricRelaxationCandidate<C> {
    type Model = PlanarLeakyDielectric<C>;
    type Step = DielectricRelaxationStep;
    fn account(&self) -> DriverAccount<'_> {
        let step = self.step();
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: step.heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: step.stored_after_j - step.stored_before_j,
            supplied_j: step.supplied_work_j + step.parameter_work_j,
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        DielectricRelaxationCandidate::into_parts(self)
    }
}

fn dielectric(error: DielectricRelaxationError) -> Error {
    match error {
        DielectricRelaxationError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(ConductionThermalError::Dielectric(other)),
    }
}

fn deposit_heat<C>(
    proposal: &DielectricRelaxationCandidate<C>,
    routes: &[CellId],
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<BTreeMap<CellId, f64>, Error> {
    let mut powers = BTreeMap::new();
    if proposal.step().cell_power_w.len() != routes.len() {
        return Err(invalid("native dielectric heat lost its original cell correspondence"));
    }
    for (&power, &cell) in proposal.step().cell_power_w.iter().zip(routes) {
        checkpoint()?;
        let previous = *powers.get(&cell).unwrap_or(&0.0);
        let combined = finite(previous + power)?;
        if power != 0.0 && combined == previous {
            return Err(invalid("local dielectric conduction heat disappeared in its receiving cell"));
        }
        powers.insert(cell, combined);
    }
    Ok(powers)
}

impl ConductionThermalNetwork<'_> {
    /// Continue actual spatial charge relaxation, finite heat and source-backed
    /// conductivity together. Each candidate restarts from the SAME old charge
    /// and potential. Endpoint T/phase controls scalar or oriented sigma; the
    /// common bundle must also retain its constant absolute epsilon source.
    /// Both original Hodge operators and their cell integrals remain native.
    ///
    /// `endpoint_voltage_v` ramps linearly from the retained initial trace.
    /// Other boundaries isolate total current. Joule heat is deposited once
    /// into each original triangle's whole receiving cell. The same fixed-volume
    /// admission allows equal-density confined phase change and rejects motion
    /// or phase-density change. Physical heat contains no numerical damping.
    ///
    /// Thermal rates use backward Euler and conductivity uses the thermal
    /// endpoint; overall coupled accuracy is first order in general. Material
    /// dispersion, variable permittivity, induction, breakdown, charge-carrier
    /// chemistry and electrodes with changing topology require their own laws.
    #[allow(clippy::too_many_arguments)]
    pub fn step_dielectric_candidate<'t, C: PlanarConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &PlanarLeakyDielectric<C>,
        endpoint_voltage_v: &[(usize, f64)],
        assembly: WeightedAssemblyLimits,
        field_controls: DielectricRelaxationControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<DielectricEnthalpyCandidate<'t, C>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        super::super::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        if field.thermal_dielectric_laws().is_some() {
            return Err(invalid("temperature-dependent dielectric storage requires step_dielectric_thermal_candidate"));
        }
        if field.time_s().to_bits() != thermal.time_s().to_bits()
            || self.cells.len() != thermal.cells().len()
        {
            return Err(invalid(
                "matching dielectric/thermal clocks and complete material bindings required",
            ));
        }
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        };
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        let problem = field.conduction();
        let routes = self.admit_routes(thermal, problem, &cancel)?;
        let initial = prepared
            .resolve(thermal, || (*cancel.borrow_mut())())
            .map_err(material)?;
        let sources = self.materials_at(problem, &initial.source, &routes, &cancel)?;
        if sources.iter().zip(problem.conductors()).any(|(a, b)| {
            !a.has_same_electrical_source(b)
                || a.frame() != b.frame()
                || !same_material(a.material(), b.material())
        }) {
            return Err(invalid(
                "accepted dielectric source bundle must match the actual thermal state",
            ));
        }
        let (thermal, field, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |state, _| {
                check(&cancel)?;
                let coefficients = prepared
                    .resolve(state, || (*cancel.borrow_mut())())
                    .map_err(material)?;
                let sources = self.materials_at(problem, &coefficients.source, &routes, &cancel)?;
                let proposal = field
                    .step_to_candidate(
                        target_s,
                        sources,
                        endpoint_voltage_v,
                        assembly,
                        field_controls,
                        &mut || (*cancel.borrow_mut())(),
                    )
                    .map_err(dielectric)?;
                Ok(DrivenCoefficients {
                    circuit: proposal,
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: coefficients.source,
                    circuit_trials: 1,
                })
            },
            |proposal: &DielectricRelaxationCandidate<C>, _dt, checkpoint| {
                deposit_heat(proposal, &routes, checkpoint)
            },
            || (*cancel.borrow_mut())(),
        )?;
        check(&cancel)?;
        Ok(DielectricEnthalpyCandidate {
            thermal,
            field,
            report,
        })
    }
}
