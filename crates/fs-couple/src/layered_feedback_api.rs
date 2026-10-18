//! Public bindings and numerical controls for layered thermal feedback.
use std::collections::BTreeMap;
use fs_material::phase::ResolvedHeatCapacityCurve;
use fs_material::thermovisco::source::ResolvedThermoMaxwell;
use fs_solid::deposition;
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::CellId;

#[derive(Clone,Copy,Debug)]
pub struct LayerThermalBinding<'a>{
    pub layer:u64,
    pub cell:CellId,
    pub material:&'a ResolvedThermoMaxwell,
    pub caloric:&'a ResolvedHeatCapacityCurve,
}

#[derive(Debug)]
pub struct ThermalCoefficients<S>{
    pub sources_w:BTreeMap<CellId,f64>,
    pub conductances_w_k:BTreeMap<u64,f64>,
    pub source:S,
}

#[derive(Clone,Copy,Debug)]
pub struct StripEnthalpyConfig{
    pub maximum_cells:usize,
    pub maximum_structural_trials:usize,
    pub maximum_newton_iterations:usize,
    pub difference_step_j:f64,
    pub cell_energy_tolerance_j:f64,
    pub thermal_energy_tolerance_j:f64,
    pub total_energy_tolerance_j:f64,
    pub mechanics:deposition::Config,
    pub newton:NewtonKrylovConfig,
}
