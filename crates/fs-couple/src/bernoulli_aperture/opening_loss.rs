//! Separated opening flow in series with the original passive outlet load.
//!
//! For actual circular area A, the declared contraction gives
//! `dp = rho Q|Q| / (2 (Cc A)^2)`. The same complete local gas supplies density
//! and the solved-flow Mach limit. The existing RLC/relaxation solver advances
//! its original history at that flow. The loss is internal dissipation, separate
//! from radiation; it is never appended to already rendered audio.
//!
//! This explicit quasi-steady, symmetric closure is not a resolved vortex model
//! or measured tone-hole coefficient. No instrument label chooses Cc. Mutual
//! exterior terminals require a coupled nonlinear solve and are refused here.

use super::{contraction::ContractionSpec, network::{NetworkNode, TubeNetworkSpec}};
use crate::acoustic_realize::AcousticRealizeError;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_vfit::{quadratic_loss::QuadraticLossSpec, waveguide::network::WaveguideNetwork};

/// Original endpoint and an explicitly sourced dimensionless flow closure.
#[derive(Clone, Debug, PartialEq)]
pub struct OpeningLossInput {
    /// Original degree-one graph node; subdivision preserves this address.
    pub node: usize,
    /// Contraction, flow limit and source; no fitted instrument preset.
    pub contraction: ContractionSpec,
}

/// Immutable source binding retained alongside the physical runtime.
#[derive(Clone, Debug, PartialEq)]
pub struct BoundOpeningLoss {
    input: OpeningLossInput,
    radius_m: f64,
    gas: GasState,
    loss: QuadraticLossSpec,
}
impl BoundOpeningLoss {
    /// Original node, closure and its provenance.
    pub const fn input(&self) -> &OpeningLossInput { &self.input }
    /// Actual circular outlet radius [m].
    pub const fn radius_m(&self) -> f64 { self.radius_m }
    /// Complete local source gas, including viscosity and thermochemical state.
    pub const fn gas(&self) -> &GasState { &self.gas }
    /// Actual dimensional law installed in the retained wave network.
    pub const fn loss(&self) -> QuadraticLossSpec { self.loss }
}

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}
fn wave(error: fs_vfit::waveguide::WaveguideError) -> AcousticRealizeError {
    AcousticRealizeError::Nonlinear(error.to_string())
}

/// Bind all requested openings atomically without replacing their radiation.
/// Graph and gases must be the matched retained runtime image (including any
/// loss-cell subdivision); original endpoint indices remain stable.
pub fn bind_opening_losses(graph: &TubeNetworkSpec, gases: &[GasState],
    network: &mut WaveguideNetwork, inputs: Vec<OpeningLossInput>, gate: &CancelGate)
    -> Result<Vec<BoundOpeningLoss>, AcousticRealizeError> {
    graph.validate_section_gases(gases)?;
    if inputs.len() > 64 { return Err(invalid("at most 64 explicit opening losses are admitted")); }
    let mut bindings = Vec::with_capacity(inputs.len());
    for input in inputs {
        if gate.is_requested() { return Err(AcousticRealizeError::Cancelled); }
        if !matches!(graph.nodes.get(input.node), Some(NetworkNode::Impedance { .. } | NetworkNode::Relaxation { .. }))
            || network.node(input.node).as_ref() != graph.nodes.get(input.node) {
            return Err(invalid("opening loss requires the original local passive terminal, including its radiation load"));
        }
        let mut adjacent = graph.sections.iter().enumerate()
            .filter(|(_, section)| section.nodes.contains(&input.node));
        let (index, section) = adjacent.next().ok_or_else(|| invalid("opening has no physical section"))?;
        if adjacent.next().is_some() { return Err(invalid("opening loss requires a degree-one terminal")); }
        if !section.radius_m.is_finite() || section.radius_m <= 0.0 {
            return Err(invalid("opening loss requires a positive finite local radius"));
        }
        let gas = gases[index];
        let area = core::f64::consts::PI*section.radius_m*section.radius_m;
        let loss = super::contraction::quadratic_opening_loss(area, &gas, &input.contraction)?;
        bindings.push(BoundOpeningLoss { input, radius_m: section.radius_m, gas, loss });
    }
    if gate.is_requested() { return Err(AcousticRealizeError::Cancelled); }
    let assignments: Vec<_> = bindings.iter().map(|b| (b.input.node, b.loss)).collect();
    network.install_terminal_quadratic_losses(&assignments).map_err(wave)?;
    Ok(bindings)
}

/// Shared file-record tail: NODE AREA_RATIO MAXIMUM_MACH SOURCE_TOKEN.
/// The enclosing file owns line numbers and row delimiters.
pub(crate) fn parse_record(mut fields: std::str::SplitAsciiWhitespace<'_>)
    -> Result<OpeningLossInput, &'static str> {
    let node = fields.next().and_then(|v| v.parse().ok())
        .ok_or("opening loss needs a node index")?;
    let area_ratio = fields.next().and_then(|v| v.parse::<f64>().ok())
        .filter(|v| v.is_finite() && *v > 0.0 && *v <= 1.0)
        .ok_or("opening contraction area ratio must be in (0,1]")?;
    let maximum_contracted_mach = fields.next().and_then(|v| v.parse::<f64>().ok())
        .filter(|v| v.is_finite() && *v > 0.0 && *v <= 0.3)
        .ok_or("opening loss Mach limit must be in (0,0.3]")?;
    let provenance = fields.next().filter(|v| !v.is_empty() && v.len() <= 1024)
        .ok_or("opening loss needs a source token of at most 1024 bytes")?.to_string();
    if fields.next().is_some() { return Err("unexpected opening loss fields"); }
    Ok(OpeningLossInput { node,
        contraction: ContractionSpec { area_ratio, maximum_contracted_mach, provenance } })
}
