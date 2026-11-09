//! Original material-product publication for reacting steady pressure changes.
//!
//! The common process owner uses native chemical SP work and entropy. Actual
//! HP products, source cards, transport choices, fixed feed mass and independent
//! material enthalpy/freezing residuals remain owned by the existing material
//! pressure/enthalpy path. Frozen material Cp/acoustics are not reacting ones.

use super::*;
use fs_matdb::QueryPoint;
use fs_material::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::equilibrium::pressure::{
        EquilibriumGasEnthalpyPath, GasEquilibriumEnthalpyError, GasEquilibriumEnthalpySolve,
        GasEquilibriumEnthalpyState,
    },
};

/// Physical flow allowances plus the original material publication bounds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialEquilibriumFlowSolve {
    /// Shared pressure-change physics and original HP/SP controls.
    pub flow: EquilibriumFlowSolve,
    /// Independent product/feed and density-volume mass allowance [kg].
    pub mass_tolerance_kg: f64,
    /// Independent material HP and composition-freezing allowance [J].
    pub material_enthalpy_tolerance_j: f64,
}

/// Resolve reacting machinery/throttling directly into original material parcels.
///
/// The context supplies inlet T,p through the existing material coordinate
/// owner; feed elements equilibrate there before the process starts. The bank
/// must retain every admitted product and original source/transport identity.
pub fn equilibrium_material_pressure_change(
    bank: &ThermochemicalGasModel,
    chemistry: &IdealGasEquilibriumV1,
    context: &QueryPoint,
    feed: &[(SpeciesId, Amount)],
    outlet_pressure: Pressure,
    process: EquilibriumPressureProcess,
    controls: MaterialEquilibriumFlowSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<EquilibriumPressureChange<GasEquilibriumEnthalpyState>, EquilibriumFlowError> {
    poll(&mut cancelled)?;
    let mut path = EquilibriumGasEnthalpyPath::new(
        bank, chemistry, context, feed,
        GasEquilibriumEnthalpySolve {
            equilibrium: controls.flow.enthalpy,
            mass_tolerance_kg: controls.mass_tolerance_kg,
            enthalpy_tolerance_j: controls.material_enthalpy_tolerance_j,
        },
        &mut cancelled,
    ).map_err(|source| material_error(EquilibriumFlowStage::Inlet, source))?;
    resolve(&mut path, chemistry, feed, outlet_pressure, process, controls.flow, &mut cancelled)
}

impl EnthalpyOwner for EquilibriumGasEnthalpyPath<'_> {
    type State = GasEquilibriumEnthalpyState;
    fn reference(&self) -> &Self::State { self.reference_state() }
    fn chemical_reference(&self) -> &EquilibriumThermodynamicStateV1 {
        self.reference_state().equilibrium()
    }
    fn mass_kg(&self) -> f64 { self.reference_state().balance().mass_kg }
    fn query(&mut self, pressure: Pressure, work: Energy, cancelled: &mut impl FnMut() -> bool)
        -> Result<Endpoint<Self::State>, EquilibriumFlowError> {
        let state = self.sample(pressure, work, cancelled)
            .map_err(|source| material_error(EquilibriumFlowStage::Outlet, source))?;
        Ok(Endpoint {
            enthalpy_change_j: state.balance().enthalpy_change_j,
            sample: state.native_sample(),
            state: state.clone(),
        })
    }
}
fn material_error(stage: EquilibriumFlowStage, source: GasEquilibriumEnthalpyError) -> EquilibriumFlowError {
    if matches!(source, GasEquilibriumEnthalpyError::Cancelled) {
        EquilibriumFlowError::Cancelled
    } else {
        EquilibriumFlowError::Material { stage, source }
    }
}
