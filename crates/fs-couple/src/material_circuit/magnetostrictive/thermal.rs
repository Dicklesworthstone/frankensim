//! Actual biased winding heat and reciprocal rod motion on the native enthalpy owner.
//!
//! The magnetostrictive material remains at its prepared, externally held core
//! temperature and DC-biased reference. Its fixed quadratic tangent describes
//! perturbations, while this series winding is driven by TOTAL terminal voltage.
//! For I=I_bias+delta_i, the original incremental circuit receives
//! delta_V=V_total-R(T) I_bias. The finite coil receives the original conductor's
//! R(T) I^2, including DC and cross terms. Physical magnetic/mechanical storage
//! CHANGE is delta_E_incremental+I_bias*delta_lambda; the unknown absolute energy
//! datum at the prepared bias is neither needed nor invented. The magnetic
//! work convention is dW=i d(lambda), as in MIT 6.685 Class Notes 1:
//! <https://web.mit.edu/6.685/www/chapter1.pdf>.
//!
//! One stationary winding, its source-pinned resistance and finite thermal
//! network share the existing circuit descriptor and cell enthalpy/Newton owner.
//! All temperature trials start from the same accepted linkage and rod motion.
//! Endpoint R(T) plus thermal backward Euler is first order in general. The
//! core's magnetic/mechanical coefficients, reference geometry and bias remain
//! fixed; this does not infer magnetic heating, thermal strain, hysteresis,
//! skin/eddy loss, thermal expansion or a moving coil's mass.

use super::{MagnetostrictiveCircuit, MagnetostrictiveCircuitError, MagnetostrictiveCircuitStep};
use crate::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, DrivenCoefficients,
    DriverAccount, ThermalDriver, cell_driver_candidate,
    material::{
        MaterialNetwork, MaterialNetworkError, MaterialNetworkSample, PreparedMaterialNetwork,
    },
};
use crate::magnetostrictive::{MagnetostrictiveCouplingError, MagnetostrictiveRod};
use crate::material_circuit::MaterialCircuitError;
use fs_material::conductor::{ElectricalPropertyBasis, ResolvedConductor};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::{cell::RefCell, collections::BTreeMap, sync::Arc};

#[derive(Debug)]
pub enum MagnetostrictiveThermalError {
    Material(MaterialNetworkError),
    Magnetic(MagnetostrictiveCircuitError),
}
impl core::fmt::Display for MagnetostrictiveThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "biased magnetostrictive circuit/enthalpy: {self:?}")
    }
}
impl std::error::Error for MagnetostrictiveThermalError {}
type Error = CircuitEnthalpyError<MagnetostrictiveThermalError>;
fn invalid(message: &'static str) -> Error {
    CircuitEnthalpyError::Invalid(message)
}
fn material(error: MaterialNetworkError) -> Error {
    match error {
        MaterialNetworkError::Thermal(DepositionError::Cancelled) => {
            CircuitEnthalpyError::Cancelled
        }
        other => CircuitEnthalpyError::Source(MagnetostrictiveThermalError::Material(other)),
    }
}
fn magnetic(error: MagnetostrictiveCircuitError) -> Error {
    match error {
        MagnetostrictiveCircuitError::Circuit(MaterialCircuitError::Cancelled)
        | MagnetostrictiveCircuitError::Magnetostrictive(
            MagnetostrictiveCouplingError::Cancelled,
        ) => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(MagnetostrictiveThermalError::Magnetic(other)),
    }
}

/// Original incremental solution together with its actual DC-biased work account.
pub struct MagnetostrictiveThermalStep {
    /// Original native circuit/rod result. Its resistance losses are incremental.
    pub incremental: MagnetostrictiveCircuitStep,
    /// Source-bound winding at the same accepted thermal endpoint.
    pub coil: ResolvedConductor,
    /// Actual signed midpoint winding current, I_bias+delta_i [A].
    pub midpoint_current_a: f64,
    /// Total held source voltage, including the DC operating voltage [V].
    pub total_voltage_v: f64,
    /// Actual positive R(T)*(I_bias+delta_i)^2*dt, already deposited once [J].
    pub coil_heat_j: f64,
    /// Linear term I_bias*delta_lambda absent from the incremental Hamiltonian [J].
    pub bias_storage_change_j: f64,
    /// Full local quadratic-model magnetic/mechanical storage change [J].
    /// Its signed value does not define an absolute biased-device energy datum.
    pub stored_change_j: f64,
    /// Actual source V_total*I_mid*dt [J]; can be negative on regeneration.
    pub electrical_supply_j: f64,
    pub mechanical_supply_j: f64,
    /// Total storage + actual coil heat - electrical/mechanical supply [J].
    pub energy_residual_j: f64,
}

/// One source-bound stationary winding, original moving rod and thermal clock.
/// The fixed three-branch series topology makes the DC operating state explicit;
/// no unknown bias currents of other RLC branches are inferred.
pub struct ThermalMagnetostrictiveCircuit<'n, 'a> {
    circuit: MagnetostrictiveCircuit,
    prepared: Arc<PreparedMaterialNetwork<'n, 'a>>,
    coil_cell: CellId,
    coil: ResolvedConductor,
    time_s: f64,
}
impl<'n, 'a> ThermalMagnetostrictiveCircuit<'n, 'a> {
    /// Prepare the original descriptor with winding R at branch 0, its ideal
    /// magnetostrictive coil at branch 1, and total-voltage source at branch 2.
    ///
    /// The material network must contain exactly one CellResistor at branch 0.
    /// It describes the whole stationary winding cell, with a wholly solid,
    /// constant-density caloric chart. Other cells and reservoirs are passive
    /// thermal recipients. The rod remains at its independently held prepared
    /// temperature; none of these cells implicitly refreshes its coefficients.
    /// Initial linkage, extension and velocity are increments about that rod.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        rod: MagnetostrictiveRod,
        linkage_increment_wb: f64,
        extension_m: f64,
        velocity_m_s: f64,
        total_voltage_v: f64,
        thermal: &DepositionSimulation<'_>,
        network: &'n MaterialNetwork<'a>,
        basis: ElectricalPropertyBasis,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        if !total_voltage_v.is_finite()
            || network.resistors.len() != 1
            || network.resistors[0].destination.branch != 0
        {
            return Err(invalid(
                "finite total voltage and one winding binding at branch 0 required",
            ));
        }
        let prepared = network
            .prepare_with_bases(thermal, &[(0, basis)], &mut cancelled)
            .map_err(material)?;
        let coil_cell = network.resistors[0].destination.cell;
        let cell = &thermal.cells()[&coil_cell];
        let binding = network
            .cells
            .iter()
            .find(|b| b.cell == coil_cell)
            .ok_or_else(|| invalid("winding lacks its original caloric binding"))?;
        let density = cell.mass_kg / cell.geometry.volume_m3;
        if !density.is_finite()
            || density <= 0.0
            || binding.curve.knots().iter().any(|k| {
                k.liquid_mass_fraction != 0.0
                    || (k.bulk_density_kg_m3 - density).abs() > 1e-10 * density
            })
        {
            return Err(invalid(
                "fixed winding geometry requires a wholly solid constant-density chart",
            ));
        }
        let mut initial = prepared
            .resolve(thermal, &mut cancelled)
            .map_err(material)?;
        let resistor = initial.resistors.remove(0);
        let coil = resistor.conductor.clone();
        let increment_v = total_voltage_v - coil.resistance_ohm() * rod.bias_current_a();
        if !increment_v.is_finite() {
            return Err(invalid("biased winding voltage is unrepresentable"));
        }
        let circuit = MagnetostrictiveCircuit::new(
            CircuitGraph {
                node_count: 3,
                branches: vec![
                    (
                        1,
                        2,
                        Branch::Resistor {
                            ohms: coil.resistance_ohm(),
                        },
                    ),
                    (
                        2,
                        0,
                        Branch::Inductor {
                            henries: rod.clamped_inductance_h(),
                        },
                    ),
                    (1, 0, Branch::VoltageSource { port: 0 }),
                ],
                transformers: vec![],
            },
            vec![resistor],
            vec![linkage_increment_wb],
            vec![],
            1,
            rod,
            extension_m,
            velocity_m_s,
            &[increment_v],
            &mut cancelled,
        )
        .map_err(magnetic)?;
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        Ok(Self {
            circuit,
            prepared: Arc::new(prepared),
            coil_cell,
            coil,
            time_s: thermal.time_s(),
        })
    }
    pub const fn circuit(&self) -> &MagnetostrictiveCircuit {
        &self.circuit
    }
    pub const fn coil(&self) -> &ResolvedConductor {
        &self.coil
    }
    pub const fn coil_cell(&self) -> CellId {
        self.coil_cell
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    /// Actual endpoint current in the held-bias local material model.
    pub fn current_a(&self) -> Result<f64, Error> {
        let current = self.circuit.rod().bias_current_a()
            + self.circuit.current_increment_a().map_err(magnetic)?;
        if current.is_finite() {
            Ok(current)
        } else {
            Err(invalid("total winding current is unrepresentable"))
        }
    }

    /// Joint voltage, reciprocal motion, actual DC-biased heat and cooling step.
    ///
    /// Total voltage and force are held over this interval. Force retains the
    /// original rod's stress-free biased-reference convention. Every endpoint
    /// thermal proposal re-evaluates the same pinned conductor source and reruns
    /// the original circuit from identical accepted storage/motion. Thermal
    /// feedback therefore changes physical current and motion on this same step.
    /// Cut source/event discontinuities and refine the step for trajectory accuracy.
    /// Both histories remain unchanged on source, budget, work or cancellation
    /// failure. The returned thermal history already contains the full coil heat.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        total_voltage_v: f64,
        force_increment_n: f64,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<MagnetostrictiveEnthalpyCandidate<'t, 'n, 'a>, Error> {
        let cancel = RefCell::new(cancelled);
        let check = || {
            if (*cancel.borrow_mut())() {
                Err(CircuitEnthalpyError::Cancelled)
            } else {
                Ok(())
            }
        };
        check()?;
        if self.time_s.to_bits() != thermal.time_s().to_bits()
            || !target_s.is_finite()
            || target_s <= self.time_s
            || !total_voltage_v.is_finite()
            || !force_increment_n.is_finite()
            || [
                controls.circuit.energy_tolerance_j,
                controls.circuit.power_tolerance_w,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(invalid(
                "matching clocks, positive interval and finite voltage/force/work controls required",
            ));
        }
        let initial = self
            .prepared
            .resolve(thermal, || (*cancel.borrow_mut())())
            .map_err(material)?;
        if initial.resistors[0].conductor != self.coil {
            return Err(invalid(
                "accepted winding source must match the actual thermal endpoint",
            ));
        }
        let dt = target_s - self.time_s;
        let bias = self.circuit.rod().bias_current_a();
        let initial_linkage = self.circuit.flux_linkage_increment_wb();
        let (thermal, driver, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |state, _remaining| {
                check()?;
                let coefficients = self
                    .prepared
                    .resolve(state, || (*cancel.borrow_mut())())
                    .map_err(material)?;
                let coil = coefficients.resistors[0].conductor.clone();
                let resistance = coil.resistance_ohm();
                let increment_v = total_voltage_v - resistance * bias;
                if !increment_v.is_finite() {
                    return Err(invalid("biased source voltage is unrepresentable"));
                }
                let proposal = self
                    .circuit
                    .step_candidate(
                        coefficients.resistors,
                        &[increment_v],
                        force_increment_n,
                        dt,
                        controls.circuit,
                        || (*cancel.borrow_mut())(),
                    )
                    .map_err(magnetic)?;
                let (circuit, incremental) = proposal.into_parts();
                if incremental.incremental_resistors.len() != 1
                    || incremental.incremental_resistors[0].material.branch != 0
                {
                    return Err(invalid("original winding branch loss is incomplete"));
                }
                let delta_current = incremental.incremental_resistors[0].current_a;
                let midpoint_current_a = bias + delta_current;
                let coil_heat_j = dt
                    * coil
                        .joule_power_w(midpoint_current_a)
                        .map_err(|e| material(MaterialNetworkError::Conductor(e.into())))?;
                let delta_linkage = circuit.flux_linkage_increment_wb() - initial_linkage;
                let bias_storage_change_j = bias * delta_linkage;
                let stored_change_j = incremental.stored_after_j - incremental.stored_before_j
                    + bias_storage_change_j;
                let electrical_supply_j = dt * total_voltage_v * midpoint_current_a;
                let mechanical_supply_j = incremental.incremental_mechanical_supply_j;
                let energy_residual_j =
                    stored_change_j + coil_heat_j - electrical_supply_j - mechanical_supply_j;
                let bias_kvl_residual_w = bias
                    * ((delta_linkage / dt) + resistance * midpoint_current_a - total_voltage_v);
                if [
                    midpoint_current_a,
                    coil_heat_j,
                    bias_storage_change_j,
                    stored_change_j,
                    electrical_supply_j,
                    mechanical_supply_j,
                    energy_residual_j,
                    bias_kvl_residual_w,
                ]
                .iter()
                .any(|v| !v.is_finite())
                    || coil_heat_j < 0.0
                    || energy_residual_j.abs() > controls.circuit.energy_tolerance_j
                    || bias_kvl_residual_w.abs() > controls.circuit.power_tolerance_w
                {
                    return Err(invalid(
                        "actual biased winding work, stored field/motion and heat do not close",
                    ));
                }
                let step = MagnetostrictiveThermalStep {
                    incremental,
                    coil: coil.clone(),
                    midpoint_current_a,
                    total_voltage_v,
                    coil_heat_j,
                    bias_storage_change_j,
                    stored_change_j,
                    electrical_supply_j,
                    mechanical_supply_j,
                    energy_residual_j,
                };
                let model = Self {
                    circuit,
                    prepared: Arc::clone(&self.prepared),
                    coil_cell: self.coil_cell,
                    coil,
                    time_s: target_s,
                };
                check()?;
                Ok(DrivenCoefficients {
                    circuit: ThermalCandidate { model, step },
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: coefficients.source,
                    circuit_trials: 1,
                })
            },
            |candidate, dt, checkpoint| {
                checkpoint()?;
                let power = candidate.step.coil_heat_j / dt;
                if !power.is_finite() || candidate.step.coil_heat_j != 0.0 && power == 0.0 {
                    return Err(invalid("actual winding heat power is unrepresentable"));
                }
                Ok(BTreeMap::from([(self.coil_cell, power)]))
            },
            || (*cancel.borrow_mut())(),
        )?;
        check()?;
        Ok(MagnetostrictiveEnthalpyCandidate {
            thermal,
            driver,
            report,
        })
    }
}
struct ThermalCandidate<'n, 'a> {
    model: ThermalMagnetostrictiveCircuit<'n, 'a>,
    step: MagnetostrictiveThermalStep,
}
impl<'n, 'a> ThermalDriver for ThermalCandidate<'n, 'a> {
    type Model = ThermalMagnetostrictiveCircuit<'n, 'a>;
    type Step = MagnetostrictiveThermalStep;
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.step.coil_heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: self.step.stored_change_j,
            supplied_j: self.step.electrical_supply_j + self.step.mechanical_supply_j,
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        (self.model, self.step)
    }
}

/// Jointly accepted rod/circuit and finite winding thermal histories.
pub struct MagnetostrictiveEnthalpyCandidate<'t, 'n, 'a> {
    thermal: DepositionSimulation<'t>,
    driver: ThermalMagnetostrictiveCircuit<'n, 'a>,
    report: CircuitEnthalpyReport<MaterialNetworkSample, MagnetostrictiveThermalStep>,
}
impl<'t, 'n, 'a> MagnetostrictiveEnthalpyCandidate<'t, 'n, 'a> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn driver(&self) -> &ThermalMagnetostrictiveCircuit<'n, 'a> {
        &self.driver
    }
    pub const fn report(
        &self,
    ) -> &CircuitEnthalpyReport<MaterialNetworkSample, MagnetostrictiveThermalStep> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        ThermalMagnetostrictiveCircuit<'n, 'a>,
        CircuitEnthalpyReport<MaterialNetworkSample, MagnetostrictiveThermalStep>,
    ) {
        (self.thermal, self.driver, self.report)
    }
}

#[cfg(test)]
mod tests;
