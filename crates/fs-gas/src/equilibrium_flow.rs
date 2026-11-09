//! Steady adiabatic pressure changes with reacting ideal-gas equilibrium.
//!
//! The original element/charge-constrained chemistry owns both the reversible
//! SP reference and the actual HP outlet. Shaft work into the fluid is Delta H;
//! flow work and formation energy already belong to H. No extra pV or reaction
//! heat is added. A throttle imposes zero work, while a compressor/expander uses
//! the authored isentropic efficiency to select work before solving the outlet.
//!
//! The supplied feed declares an elemental inventory. Its equilibrium at the
//! inlet T,p defines the inlet, not the enthalpy of an unequilibrated feed. This
//! model neglects changes in kinetic/potential energy and does not predict mass
//! flow, geometry, choking, finite-rate reaction, phase separation or transport.
//! The optional material adapter retains original frozen-product transport;
//! the separately reported response and entropy remain chemical equilibrium.

use fs_qty::{Amount, Energy, MassFlowRate, Power, Pressure, Temperature};
use fs_thermochem::SpeciesId;
use fs_thermochem::equilibrium::composition::{
    CompositionEquilibriumErrorV1, EquilibriumEnthalpySolveV1, EquilibriumEntropySolveV1,
    EquilibriumExpansionErrorV1, EquilibriumThermodynamicStateV1, IdealGasEquilibriumV1,
    enthalpy::{EquilibriumEnthalpyPressurePathV1, EquilibriumEnthalpyPressureSampleV1},
    isentrope::EquilibriumIsentropeV1,
};

/// Existing source-bearing material parcels over this same process balance.
#[cfg(feature = "material")]
pub mod material;

/// Pressure direction and convention for work into the fluid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum EquilibriumPressureProcess {
    /// Passive pressure drop with zero heat and shaft work.
    Throttle,
    /// Pressure rises; actual work is reversible work / efficiency.
    Compressor { isentropic_efficiency: f64 },
    /// Pressure falls; actual signed work is reversible work * efficiency.
    Expander { isentropic_efficiency: f64 },
}

/// Original source domains, work bounds and extensive physical allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumFlowSolve {
    /// Actual HP outlet; all original Gibbs/temperature limits remain active.
    pub enthalpy: EquilibriumEnthalpySolveV1,
    /// Reversible SP outlet for machinery. Its Gibbs policy must equal HP's.
    /// A throttle and an exactly unchanged pressure do not run an SP solve.
    pub isentropic: EquilibriumEntropySolveV1,
    /// Maximum absolute actual Delta H minus declared shaft work [J].
    /// Never scaled by absolute source formation enthalpies.
    pub energy_tolerance_j: f64,
    /// Allowed numerical shortfall from nonnegative chemical Delta S [J/K].
    /// The measured value is retained without clipping.
    pub entropy_production_tolerance_j_k: f64,
}

/// Original stage responsible for a refusal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EquilibriumFlowStage {
    /// Inlet equilibrium and its original source reference.
    Inlet,
    /// Reversible pressure change at the inlet chemical entropy.
    Isentropic,
    /// Actual outlet at the imposed work and pressure.
    Outlet,
}

/// No failed or cancelled process publishes a partial outlet.
#[derive(Clone, Debug, PartialEq)]
pub enum EquilibriumFlowError {
    /// Invalid input, incompatible references or unrepresentable work.
    Invalid(&'static str),
    /// Original source, conservation, inverse, derivative or work-budget refusal.
    Chemistry {
        stage: EquilibriumFlowStage,
        source: EquilibriumExpansionErrorV1,
    },
    /// Original material source/transport/freezing refusal.
    #[cfg(feature = "material")]
    Material {
        stage: EquilibriumFlowStage,
        source: fs_material::thermochemical_gas::energy::equilibrium::pressure::GasEquilibriumEnthalpyError,
    },
    /// Actual source energy change failed the independent process allowance.
    Energy { residual_j: f64, tolerance_j: f64 },
    /// Chemical entropy production failed the adiabatic second-law allowance.
    Entropy { change_j_k: f64, tolerance_j_k: f64 },
    /// Cancellation before complete publication.
    Cancelled,
}
impl core::fmt::Display for EquilibriumFlowError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "reacting equilibrium pressure change: {self:?}")
    }
}
impl std::error::Error for EquilibriumFlowError {}

/// Extensive accounts for the original conserved representative feed mass.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumFlowBalance {
    /// Original feed mass [kg], independent of numerical product residuals.
    pub mass_kg: f64,
    /// Authored process; equal inlet/outlet pressure always means zero work.
    pub process: EquilibriumPressureProcess,
    /// Source-relative reversible Delta H [J]; absent for throttling.
    pub isentropic_enthalpy_change_j: Option<f64>,
    /// Declared shaft work INTO the fluid [J]; negative means power output.
    /// Exactly zero for throttles, including when their numerical residual is not.
    pub shaft_work_j: f64,
    /// Actual original-source or material Delta H [J].
    pub enthalpy_change_j: f64,
    /// Actual Delta H minus declared work [J], not an extra heat/work source.
    pub energy_residual_j: f64,
    /// Chemical-equilibrium outlet minus inlet entropy [J/K].
    /// With no heat transfer this is entropy production for the modeled stream.
    pub entropy_production_j_k: f64,
    /// Native chemical product-minus-feed mass / feed mass.
    pub relative_mass_residual: f64,
    /// Actual Gibbs iterations across inlet, SP preparation/query and HP query.
    pub composition_iterations: usize,
    /// Actual temperature iterations across SP and HP queries (not brackets).
    pub temperature_iterations: usize,
}

/// Steady rates obtained from the accepted representative feed accounts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumFlowRates {
    /// Positive prescribed mass flow; this endpoint model does not solve it.
    pub mass_flow: MassFlowRate,
    /// Shaft power into the fluid, negative for an expander [W].
    pub shaft_power: Power,
    /// Actual stream enthalpy power difference [W].
    pub enthalpy_power_change: Power,
    /// Numerical enthalpy-power closure defect [W].
    pub energy_residual: Power,
    /// Chemical entropy production rate [W/K].
    pub entropy_production_w_k: f64,
}
impl EquilibriumFlowBalance {
    /// Scale this same conserved inventory to an explicit positive mass flow.
    pub fn rates(self, mass_flow: MassFlowRate) -> Result<EquilibriumFlowRates, EquilibriumFlowError> {
        if !positive(mass_flow.value()) || !positive(self.mass_kg) {
            return Err(EquilibriumFlowError::Invalid("positive finite feed mass and mass flow required"));
        }
        let factor = mass_flow.value() / self.mass_kg;
        if !positive(factor) {
            return Err(EquilibriumFlowError::Invalid("unrepresentable stream inventory rate"));
        }
        let scale = |value: f64| {
            let result = value * factor;
            if !result.is_finite() || (value != 0.0 && result == 0.0) {
                Err(EquilibriumFlowError::Invalid("unrepresentable stream power or entropy rate"))
            } else {
                Ok(result)
            }
        };
        Ok(EquilibriumFlowRates {
            mass_flow,
            shaft_power: Power::new(scale(self.shaft_work_j)?),
            enthalpy_power_change: Power::new(scale(self.enthalpy_change_j)?),
            energy_residual: Power::new(scale(self.energy_residual_j)?),
            entropy_production_w_k: scale(self.entropy_production_j_k)?,
        })
    }
}

/// Complete original-source inlet/outlet and the admitted process balance.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumPressureChange<S> {
    /// Immutable equilibrium inlet, including the original source receipts.
    pub inlet: S,
    /// Actual equilibrium outlet, with original material products when selected.
    pub outlet: S,
    /// Reversible native chemical endpoint; absent for throttles and zero work.
    pub isentropic_outlet: Option<EquilibriumThermodynamicStateV1>,
    /// Energy, mass, entropy and bounded-work accounts for the same feed.
    pub balance: EquilibriumFlowBalance,
}

/// Solve a source-bound adiabatic pressure change, including chemical reactions.
///
/// There are at most two original reference TP solves and one SP plus one HP
/// query; their explicit native Gibbs, backtracking and temperature limits bound
/// all work. A throttle needs only the HP reference and query. Borrowed source
/// models/feed are unchanged, and cancellation is polled before publication.
pub fn equilibrium_pressure_change(
    chemistry: &IdealGasEquilibriumV1,
    inlet_temperature: Temperature,
    inlet_pressure: Pressure,
    feed: &[(SpeciesId, Amount)],
    outlet_pressure: Pressure,
    process: EquilibriumPressureProcess,
    controls: EquilibriumFlowSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<EquilibriumPressureChange<EquilibriumThermodynamicStateV1>, EquilibriumFlowError> {
    poll(&mut cancelled)?;
    admit(controls, process, inlet_pressure, outlet_pressure)?;
    let mut path = EquilibriumEnthalpyPressurePathV1::new(
        chemistry, inlet_temperature, inlet_pressure, feed, controls.enthalpy, &mut cancelled,
    ).map_err(|source| chemistry_error(EquilibriumFlowStage::Inlet, source))?;
    resolve(&mut path, chemistry, feed, outlet_pressure, process, controls, &mut cancelled)
}

struct Endpoint<S> {
    state: S,
    enthalpy_change_j: f64,
    sample: EquilibriumEnthalpyPressureSampleV1,
}

// This adapter shares process physics, not thermodynamic inversions. Each
// implementation delegates all chemistry/material publication to its owner.
trait EnthalpyOwner {
    type State: Clone;
    fn reference(&self) -> &Self::State;
    fn chemical_reference(&self) -> &EquilibriumThermodynamicStateV1;
    fn mass_kg(&self) -> f64;
    fn query(&mut self, pressure: Pressure, work: Energy, cancelled: &mut impl FnMut() -> bool)
        -> Result<Endpoint<Self::State>, EquilibriumFlowError>;
}
impl EnthalpyOwner for EquilibriumEnthalpyPressurePathV1<'_> {
    type State = EquilibriumThermodynamicStateV1;
    fn reference(&self) -> &Self::State { self.reference_state() }
    fn chemical_reference(&self) -> &EquilibriumThermodynamicStateV1 { self.reference_state() }
    fn mass_kg(&self) -> f64 { self.feed_mass_kg() }
    fn query(&mut self, pressure: Pressure, work: Energy, cancelled: &mut impl FnMut() -> bool)
        -> Result<Endpoint<Self::State>, EquilibriumFlowError> {
        let sample = self.sample(pressure, work, &mut *cancelled)
            .map_err(|source| chemistry_error(EquilibriumFlowStage::Outlet, source))?;
        let state = self.publish_current(cancelled)
            .map_err(|source| chemistry_error(EquilibriumFlowStage::Outlet, source))?;
        Ok(Endpoint { state, enthalpy_change_j: sample.enthalpy_change_j, sample })
    }
}

fn resolve<O: EnthalpyOwner>(
    path: &mut O, chemistry: &IdealGasEquilibriumV1, feed: &[(SpeciesId, Amount)],
    outlet_pressure: Pressure, process: EquilibriumPressureProcess,
    controls: EquilibriumFlowSolve, cancelled: &mut impl FnMut() -> bool,
) -> Result<EquilibriumPressureChange<O::State>, EquilibriumFlowError> {
    poll(cancelled)?;
    let reference = path.chemical_reference();
    let inlet = path.reference().clone();
    let initial = reference.equilibrium();
    let mass = path.mass_kg();
    admit(controls, process, initial.pressure(), outlet_pressure)?;
    if !positive(mass) || !initial.is_bound_to(chemistry) {
        return Err(EquilibriumFlowError::Invalid("process inlet must retain the original source and feed mass"));
    }
    let mut balance = EquilibriumFlowBalance {
        mass_kg: mass, process, isentropic_enthalpy_change_j: None,
        shaft_work_j: 0.0, enthalpy_change_j: 0.0, energy_residual_j: 0.0,
        entropy_production_j_k: 0.0, relative_mass_residual: initial.relative_mass_residual(),
        composition_iterations: initial.composition_iterations(), temperature_iterations: 0,
    };
    if initial.pressure() == outlet_pressure {
        if !matches!(process, EquilibriumPressureProcess::Throttle) {
            balance.isentropic_enthalpy_change_j = Some(0.0);
        }
        poll(cancelled)?;
        return Ok(EquilibriumPressureChange {
            outlet: inlet.clone(), inlet, isentropic_outlet: None, balance,
        });
    }
    let isentropic_outlet = if matches!(process, EquilibriumPressureProcess::Throttle) {
        None
    } else {
        let mut reversible = EquilibriumIsentropeV1::new(
            chemistry, initial.temperature(), initial.pressure(), feed,
            controls.isentropic, &mut *cancelled,
        ).map_err(|source| chemistry_error(EquilibriumFlowStage::Isentropic, source))?;
        // Identical original Gibbs policy must establish the same physical and
        // source reference, including trace species; no fitted reference offset.
        if reversible.reference_state() != reference || reversible.feed_mass_kg() != mass {
            return Err(EquilibriumFlowError::Invalid("HP and SP must share the identical original equilibrium inlet"));
        }
        balance.composition_iterations += reversible.reference_state().equilibrium().composition_iterations();
        let sample = reversible.sample(outlet_pressure, &mut *cancelled)
            .map_err(|source| chemistry_error(EquilibriumFlowStage::Isentropic, source))?;
        balance.composition_iterations += sample.composition_iterations;
        balance.temperature_iterations += sample.temperature_iterations;
        let ideal = sample.specific_enthalpy_change_j_kg * mass;
        admit_work(process, ideal)?;
        balance.isentropic_enthalpy_change_j = Some(ideal);
        balance.shaft_work_j = match process {
            EquilibriumPressureProcess::Compressor { isentropic_efficiency } => ideal / isentropic_efficiency,
            EquilibriumPressureProcess::Expander { isentropic_efficiency } => ideal * isentropic_efficiency,
            EquilibriumPressureProcess::Throttle => unreachable!(),
        };
        admit_work(process, balance.shaft_work_j)?;
        Some(reversible.publish_current(&mut *cancelled)
            .map_err(|source| chemistry_error(EquilibriumFlowStage::Isentropic, source))?)
    };
    let endpoint = path.query(outlet_pressure, Energy::new(balance.shaft_work_j), cancelled)?;
    balance.enthalpy_change_j = endpoint.enthalpy_change_j;
    balance.energy_residual_j = endpoint.enthalpy_change_j - balance.shaft_work_j;
    balance.entropy_production_j_k = endpoint.sample.entropy_change_j_k;
    balance.relative_mass_residual = endpoint.sample.relative_mass_residual;
    balance.composition_iterations += endpoint.sample.composition_iterations;
    balance.temperature_iterations += endpoint.sample.temperature_iterations;
    if !balance.energy_residual_j.is_finite() || balance.energy_residual_j.abs() > controls.energy_tolerance_j {
        return Err(EquilibriumFlowError::Energy {
            residual_j: balance.energy_residual_j, tolerance_j: controls.energy_tolerance_j,
        });
    }
    if !balance.entropy_production_j_k.is_finite()
        || balance.entropy_production_j_k < -controls.entropy_production_tolerance_j_k {
        return Err(EquilibriumFlowError::Entropy {
            change_j_k: balance.entropy_production_j_k,
            tolerance_j_k: controls.entropy_production_tolerance_j_k,
        });
    }
    if !matches!(process, EquilibriumPressureProcess::Throttle) {
        admit_work(process, balance.enthalpy_change_j)?;
    }
    poll(cancelled)?;
    Ok(EquilibriumPressureChange { inlet, outlet: endpoint.state, isentropic_outlet, balance })
}

fn admit(controls: EquilibriumFlowSolve, process: EquilibriumPressureProcess,
    inlet: Pressure, outlet: Pressure) -> Result<(), EquilibriumFlowError> {
    if !positive(inlet.value()) || !positive(outlet.value())
        || !positive(controls.energy_tolerance_j)
        || !positive(controls.entropy_production_tolerance_j_k) {
        return Err(EquilibriumFlowError::Invalid("positive finite pressures and physical residual allowances required"));
    }
    let efficiency = match process {
        EquilibriumPressureProcess::Throttle => None,
        EquilibriumPressureProcess::Compressor { isentropic_efficiency }
        | EquilibriumPressureProcess::Expander { isentropic_efficiency } => Some(isentropic_efficiency),
    };
    if let Some(efficiency) = efficiency {
        if !positive(efficiency) || efficiency > 1.0 {
            return Err(EquilibriumFlowError::Invalid("isentropic efficiency must be finite and in (0,1]"));
        }
        if controls.enthalpy.equilibrium != controls.isentropic.equilibrium {
            return Err(EquilibriumFlowError::Invalid("HP and SP require identical original Gibbs controls"));
        }
    }
    let direction_ok = match process {
        EquilibriumPressureProcess::Compressor { .. } => outlet.value() >= inlet.value(),
        EquilibriumPressureProcess::Throttle | EquilibriumPressureProcess::Expander { .. } => outlet.value() <= inlet.value(),
    };
    if !direction_ok {
        return Err(EquilibriumFlowError::Invalid("outlet pressure contradicts the declared process"));
    }
    Ok(())
}
fn admit_work(process: EquilibriumPressureProcess, work: f64) -> Result<(), EquilibriumFlowError> {
    let sign_ok = match process {
        EquilibriumPressureProcess::Compressor { .. } => work > 0.0,
        EquilibriumPressureProcess::Expander { .. } => work < 0.0,
        EquilibriumPressureProcess::Throttle => work == 0.0,
    };
    if !work.is_finite() || !sign_ok {
        Err(EquilibriumFlowError::Invalid("nonzero pressure change requires resolved correctly signed source work"))
    } else {
        Ok(())
    }
}
fn positive(value: f64) -> bool { value.is_finite() && value > 0.0 }
fn chemistry_error(stage: EquilibriumFlowStage, source: EquilibriumExpansionErrorV1) -> EquilibriumFlowError {
    if matches!(source, EquilibriumExpansionErrorV1::Equilibrium(CompositionEquilibriumErrorV1::Cancelled)) {
        EquilibriumFlowError::Cancelled
    } else {
        EquilibriumFlowError::Chemistry { stage, source }
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), EquilibriumFlowError> {
    if cancelled() { Err(EquilibriumFlowError::Cancelled) } else { Ok(()) }
}

#[cfg(test)]
mod tests;
