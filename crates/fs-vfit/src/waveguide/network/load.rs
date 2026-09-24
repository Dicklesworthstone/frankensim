//! Two admitted terminal realizations sharing one simultaneous network step.
use super::NetworkNode;
use super::admittance::{AdmittanceFrame, RelaxationAdmittance};
use crate::impedance::{ImpedanceFrame, ImpedanceState, SeriesImpedance};
use crate::relaxation::{RelaxationFrame, RelaxationImpedance};
use crate::waveguide::WaveguideError;

// Inline bounded histories are included by the parent's memory admission.
#[derive(Clone, Copy)]
#[allow(clippy::large_enum_variant)] // bounded inline storage is admitted in the node payload
pub(super) enum TerminalLoad {
    Series(SeriesImpedance),
    Relaxation(RelaxationImpedance),
    Admittance(RelaxationAdmittance),
}
#[derive(Clone, Copy)]
pub(super) enum TerminalFrame {
    Series(ImpedanceFrame),
    Relaxation(RelaxationFrame),
    Admittance(AdmittanceFrame),
}
impl Default for TerminalFrame {
    fn default() -> Self { Self::Series(ImpedanceFrame::default()) }
}
impl TerminalFrame {
    pub(super) fn port(self) -> ImpedanceFrame {
        match self { Self::Series(f) => f, Self::Relaxation(f) => f.port,
            Self::Admittance(f) => ImpedanceFrame {
                // Generic private port observations carry no base R-L-C state.
                // Typed pressure histories are exposed separately by the owner.
                state: ImpedanceState::default(), reflected_pressure_pa:f.reflected_pressure_pa,
                pressure_pa:f.pressure_pa, flow_m3_s:f.flow_m3_s, stored_energy_j:f.stored_energy_j,
                storage_change_j:f.storage_change_j,dissipated_energy_j:f.dissipated_energy_j,
                supplied_work_j:f.supplied_work_j,
            }
        }
    }
}
impl TerminalLoad {
    pub(super) fn new(node: NetworkNode, z: f64, dt: f64) -> Result<Option<Self>, WaveguideError> {
        Ok(match node {
            NetworkNode::Impedance { load } => Some(Self::Series(SeriesImpedance::new(load, z, dt)?)),
            NetworkNode::Relaxation { load } | NetworkNode::Shunt { load } | NetworkNode::Series { load } => Some(Self::Relaxation(RelaxationImpedance::new(load, z, dt)?)),
            NetworkNode::ShuntAdmittance {load} => Some(Self::Admittance(RelaxationAdmittance::new(load,z,dt)?)),
            _ => None,
        })
    }
    pub(super) fn stored_energy_j(&self) -> f64 {
        match self { Self::Series(l) => l.stored_energy_j(), Self::Relaxation(l) => l.stored_energy_j(), Self::Admittance(l) => l.stored_energy_j() }
    }
    pub(super) fn state(&self) -> Option<ImpedanceState> {
        match self { Self::Series(l) => Some(l.state()), Self::Relaxation(l) => Some(l.base_state()), Self::Admittance(_) => None }
    }
    pub(super) fn relaxation_flows(&self) -> Option<&[f64]> {
        match self { Self::Series(_) | Self::Admittance(_) => None, Self::Relaxation(l) => Some(l.branch_flows_m3_s()) }
    }
    pub(super) fn relaxation_pressures(&self) -> Option<&[f64]> {
        match self {Self::Admittance(l)=>Some(l.branch_pressures_pa()), _=>None}
    }
    pub(super) fn preview_step(&self, wave: f64) -> Result<TerminalFrame, WaveguideError> {
        match self {
            Self::Series(l) => l.preview_step(wave).map(TerminalFrame::Series),
            Self::Relaxation(l) => l.preview_step(wave).map(TerminalFrame::Relaxation),
            Self::Admittance(l) => l.preview_step(wave).map(TerminalFrame::Admittance),
        }
    }
    pub(super) fn accept_frame(&mut self, frame: TerminalFrame) {
        match (self, frame) {
            (Self::Series(l), TerminalFrame::Series(f)) => l.accept_frame(f),
            (Self::Relaxation(l), TerminalFrame::Relaxation(f)) => l.accept_frame(f),
            (Self::Admittance(l), TerminalFrame::Admittance(f)) => l.accept_frame(f),
            // Node kinds and their loads cannot change after admission. Every
            // occupied slot was previewed in the same call before publication.
            _ => unreachable!("terminal preview kind must match its fixed admitted load"),
        }
    }
}
