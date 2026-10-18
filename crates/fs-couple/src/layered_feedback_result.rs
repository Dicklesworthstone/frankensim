//! Candidate histories and final reports for layered thermal feedback.
use std::collections::BTreeMap;
use fs_material::thermovisco::Energy;
use fs_solid::deposition::{self,LayeredStrip};
use fs_solver::NewtonReport;
use fs_thermal::enthalpy::{DepositionSimulation,EnergyAudit};

pub struct StripEnthalpyReport<S>{
    pub mechanical:deposition::Report,
    pub layer_energy_j:BTreeMap<u64,Energy>,
    pub thermal:EnergyAudit,
    pub source:S,
    pub nonlinear:NewtonReport,
    pub structural_trials:usize,
    pub maximum_cell_residual_j:f64,
    pub deposited_heat_j:f64,
    pub non_mechanical_heat_j:f64,
    pub external_mechanical_work_j:f64,
    pub total_energy_residual_j:f64,
}

pub struct StripEnthalpyCandidate<'t,'m,S>{
    pub(in crate::strip_enthalpy) thermal:DepositionSimulation<'t>,
    pub(in crate::strip_enthalpy) strip:LayeredStrip<'m>,
    pub(in crate::strip_enthalpy) report:StripEnthalpyReport<S>,
}
impl<'t,'m,S> StripEnthalpyCandidate<'t,'m,S>{
    pub const fn thermal(&self)->&DepositionSimulation<'t>{&self.thermal}
    pub const fn strip(&self)->&LayeredStrip<'m>{&self.strip}
    pub const fn report(&self)->&StripEnthalpyReport<S>{&self.report}
    pub fn into_parts(self)->(DepositionSimulation<'t>,LayeredStrip<'m>,StripEnthalpyReport<S>){
        (self.thermal,self.strip,self.report)
    }
}
