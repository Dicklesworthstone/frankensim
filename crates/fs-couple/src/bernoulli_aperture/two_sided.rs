//! One material valve loaded by BOTH upstream and downstream acoustic networks.
//!
//! Each network's inlet points away from the valve. For through-flow Q,
//! b_up = a_up - Z_up Q and b_down = a_down + Z_down Q. An explicit series
//! pressure source P_s supplies P_s Q, so the pressure difference across the
//! moving faces is P_s + 2(a_up-a_down) - (Z_up+Z_down)Q. The existing valve
//! therefore solves the complete load at Z_up+Z_down; no extra nonlinear law,
//! iteration over lagged pressures, or sample of feedback delay is introduced.
//!
//! The series-impedance model is the reed/tract approximation in Fritz & Wolfe,
//! JASA 118 (2005), arXiv:physics/0505195. This adapter does not model lungs,
//! mean-flow convection, compressible valve flow, or changing tract geometry.
//! The two inlet gas states must match: one small-signal volume flow cannot
//! conserve mass across arbitrarily different densities. Remote sections may
//! retain their own frozen, common-pressure gas states and physical loads.

use super::dynamic::{ApertureDrive, ApertureFrame, ApertureProgress, ApertureTerminal,
    ApertureTrial, DynamicAperture};
use super::network::{NetworkNode, NodeFrame, SectionRealization, TubeNetworkSpec};
use super::paired::{PairedAperture, PairedTrial};
use super::viscothermal::{ViscothermalSection, ViscothermalSectionSpec,
    with_regional_viscothermal_sections};
use crate::acoustic_realize::AcousticRealizeError;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_vfit::waveguide::network::{NetworkFrame, WaveguideNetwork};

/// Physical side of the valve; each side's network inlet points AWAY from it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ValveSide {
    /// Positive through-flow leaves this network.
    Upstream,
    /// Positive through-flow enters this network.
    Downstream,
}
impl ValveSide {
    fn index(self) -> usize { match self { Self::Upstream => 0, Self::Downstream => 1 } }
}

/// Original physical inputs for one acoustic side, before loss subdivision.
#[derive(Clone, Debug)]
pub struct ApertureDuctSpec {
    /// Geometry, passive loads and original wave-payload budget.
    pub graph: TubeNetworkSpec,
    /// One complete coherent gas state per original section.
    pub gases: Vec<GasState>,
    /// Existing distributed-loss selections; empty explicitly means none.
    pub viscothermal: Vec<ViscothermalSectionSpec>,
}

/// Existing mechanical owners. Source geometry, contacts and material histories
/// remain inside these objects; this enum adds no alternative reed model.
pub enum LoadedAperture {
    /// One opening coordinate, including the existing symmetric-pair reduction.
    Single(Box<DynamicAperture>),
    /// Two independently moving blades and their independent material memories.
    Paired(Box<PairedAperture>),
}
impl From<DynamicAperture> for LoadedAperture {
    fn from(value: DynamicAperture) -> Self { Self::Single(Box::new(value)) }
}
impl From<PairedAperture> for LoadedAperture {
    fn from(value: PairedAperture) -> Self { Self::Paired(Box::new(value)) }
}
impl LoadedAperture {
    /// Original valve clock [s]. Its load declaration is the COMBINED impedance.
    #[must_use]
    pub fn time_step_s(&self) -> f64 {
        match self { Self::Single(a) => a.spec().time_step_s, Self::Paired(a) => a.spec().time_step_s }
    }
    /// Accepted samples; a coupled refusal never changes this value.
    #[must_use]
    pub fn accepted_steps(&self) -> u64 {
        match self { Self::Single(a) => a.accepted_steps(), Self::Paired(a) => a.accepted_steps() }
    }
    /// Original total lifetime sample budget.
    #[must_use]
    pub fn max_steps(&self) -> u64 {
        match self { Self::Single(a) => a.spec().max_steps, Self::Paired(a) => a.spec().max_steps }
    }
    fn load(&self) -> (f64, f64) {
        match self {
            Self::Single(a) => (a.spec().density_kg_m3, a.spec().impedance_pa_s_m3),
            Self::Paired(a) => (a.spec().density_kg_m3, a.spec().impedance_pa_s_m3),
        }
    }
    fn energy(&self) -> f64 {
        match self { Self::Single(a) => a.stored_energy_j(), Self::Paired(a) => a.stored_energy_j() }
    }
    fn preview(&self, drive: ApertureDrive, forces: [f64; 2]) -> Result<ValveTrial, AcousticRealizeError> {
        match self {
            Self::Single(a) => {
                if forces[1] != 0.0 { return Err(invalid("single-coordinate valve has no second force port")); }
                a.preview_step_with_force(drive, forces[0]).map(ValveTrial::Single)
            }
            Self::Paired(a) => a.preview(drive, forces).map(ValveTrial::Paired),
        }
    }
    fn accept(&mut self, trial: ValveTrial) {
        match (self, trial) {
            (Self::Single(a), ValveTrial::Single(t)) => a.accept_frame(t),
            (Self::Paired(a), ValveTrial::Paired(t)) => a.accept(t),
            _ => unreachable!("candidate belongs to the exclusively held valve"),
        }
    }
}
enum ValveTrial { Single(ApertureTrial), Paired(PairedTrial) }
impl ValveTrial {
    fn frame(&self) -> ApertureFrame {
        match self { Self::Single(t) => t.frame, Self::Paired(t) => t.frame.aperture }
    }
}

/// Explicit whole-system limits, additional to every participant's original gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TwoSidedLimits {
    /// Maximum actual valve + both wave/load storage [J].
    pub maximum_energy_j: f64,
    /// Absolute full-step energy-balance allowance [J].
    pub energy_absolute_tolerance_j: f64,
    /// Relative full-step allowance, in (0, 1e-3].
    pub energy_relative_tolerance: f64,
    /// Absolute volume-flow mismatch allowance [m³/s].
    pub flow_absolute_tolerance_m3_s: f64,
    /// Relative volume-flow mismatch allowance, in (0, 1e-3].
    pub flow_relative_tolerance: f64,
}
impl TwoSidedLimits {
    fn validate(self) -> Result<(), AcousticRealizeError> {
        if !self.maximum_energy_j.is_finite() || self.maximum_energy_j <= 0.0
            || [self.energy_absolute_tolerance_j, self.flow_absolute_tolerance_m3_s]
                .iter().any(|x| !x.is_finite() || *x < 0.0)
            || [self.energy_relative_tolerance, self.flow_relative_tolerance]
                .iter().any(|x| !x.is_finite() || *x <= 0.0 || *x > 1e-3)
        { return Err(invalid("two-sided valve needs explicit finite storage, work and flow limits")); }
        Ok(())
    }
}

/// Held physical supplies, not pitch or output-amplitude controls.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct TwoSidedDrive {
    /// Pressure rise [Pa] across the ideal SERIES source on the upstream side.
    /// The actual pressure on the upstream blade face also includes tract feedback.
    pub source_pressure_pa: f64,
    /// Independent opening-positive forces [N]; second must be zero for Single.
    pub blade_forces_n: [f64; 2],
}

/// Complete accepted observation. No acoustic load is counted as a source.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct TwoSidedFrame {
    /// Original mechanics/contact/material accounts. Pressure and outgoing wave
    /// have been mapped back to the PHYSICAL downstream port, not the eliminated
    /// combined load. The local pressure-work term uses the actual face difference.
    pub aperture: ApertureFrame,
    /// Actual upstream face pressure including the series source [Pa].
    pub upstream_pressure_pa: f64,
    /// Actual upstream-minus-downstream face pressure [Pa].
    pub pressure_drop_pa: f64,
    /// Upstream network port/work/storage; its inlet flow is minus through-flow.
    pub upstream: NetworkFrame,
    /// Downstream network port/work/storage; its inlet flow is through-flow.
    pub downstream: NetworkFrame,
    /// Actual valve + both networks' wave and reactive-load storage [J].
    pub stored_energy_j: f64,
    /// Independently evaluated change of total storage [J].
    pub storage_change_j: f64,
    /// All valve, contact, material and two-network load losses, counted once [J].
    pub dissipated_energy_j: f64,
    /// Explicit series-source work P_s Q dt [J], which may be negative.
    pub source_work_j: f64,
}
impl TwoSidedFrame {
    /// No junction or nonlinear-solve defect is subtracted from this balance.
    #[must_use]
    pub fn balance_residual_j(&self) -> f64 {
        self.storage_change_j + self.dissipated_energy_j
            - self.source_work_j - self.aperture.mechanical_work_j
    }
}

struct DuctSide {
    source: ApertureDuctSpec,
    graph: TubeNetworkSpec,
    gases: Vec<GasState>,
    represented: Vec<SectionRealization>,
    losses: Vec<ViscothermalSection>,
    network: WaveguideNetwork,
}
fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}
fn wave(error: fs_vfit::waveguide::WaveguideError) -> AcousticRealizeError {
    AcousticRealizeError::Nonlinear(error.to_string())
}
fn poll(gate: &CancelGate) -> Result<(), AcousticRealizeError> {
    if gate.is_requested() { Err(AcousticRealizeError::Cancelled) } else { Ok(()) }
}

/// Resolve the one inlet fluid and SUM of characteristic impedances used by
/// the original mechanical valve constructor. No already-running valve is retuned.
/// Only the two inlet states must agree; remote regional states remain explicit.
pub fn combined_inlet_state(sides: &[ApertureDuctSpec; 2])
    -> Result<(GasState, f64), AcousticRealizeError>
{
    let mut inlets = [None; 2];
    let mut z = 0.0;
    for (i, side) in sides.iter().enumerate() {
        z += side.graph.inlet_impedance_with_gases(&side.gases)?;
        inlets[i] = Some(side.gases[side.graph.inlet_section_index()?]);
    }
    let gas = inlets[0].expect("validated nonempty inlet");
    if inlets[1] != Some(gas) || !z.is_finite() || z <= 0.0 {
        return Err(invalid("two-sided valve needs matching coherent inlet gases and finite summed impedance"));
    }
    Ok((gas, z))
}

/// One valve and two retained acoustic loads. Neither network is an output filter.
/// The caller's source graphs/gases/loss selections remain available unchanged.
/// No additional per-step allocation beyond the existing valve owners is added.
pub struct TwoSidedAperture {
    valve: LoadedAperture,
    sides: [DuctSide; 2],
    limits: TwoSidedLimits,
}
impl TwoSidedAperture {
    /// Attach two quiescent networks to an unadvanced valve, preserving initial
    /// blade/contact/material energy. Construct its load using combined_inlet_state.
    /// A source already stepped, a stale gas, or a one-sided impedance refuses.
    pub fn new(valve: LoadedAperture, sources: [ApertureDuctSpec; 2],
        limits: TwoSidedLimits, gate: &CancelGate) -> Result<Self, AcousticRealizeError>
    {
        poll(gate)?;
        limits.validate()?;
        let (gas, z) = combined_inlet_state(&sources)?;
        let (density, declared_z) = valve.load();
        if valve.accepted_steps() != 0 || density.to_bits() != gas.density.to_bits()
            || declared_z.to_bits() != z.to_bits()
            || matches!(&valve, LoadedAperture::Paired(a) if a.gas().is_some_and(|g| *g != gas))
        { return Err(invalid("two-sided coupling requires an unadvanced valve with the actual inlet gas and summed load")); }
        let dt = valve.time_step_s();
        let prepare = |source: ApertureDuctSpec| -> Result<DuctSide, AcousticRealizeError> {
            poll(gate)?;
            let (graph, gases, losses) = with_regional_viscothermal_sections(
                source.graph.clone(), source.gases.clone(), dt, &source.viscothermal, gate)?;
            let (network, represented) = graph.realize_waveguide(&gases, dt, gate)?;
            Ok(DuctSide { source, graph, gases, represented, losses, network })
        };
        let [up, down] = sources;
        let model = Self { valve, sides: [prepare(up)?, prepare(down)?], limits };
        // A lowering must preserve the exact boundary used in the mechanical solve.
        let actual_z = model.sides.iter().map(|s| s.network.inlet_impedance_pa_s_m3()).sum::<f64>();
        if actual_z.to_bits() != z.to_bits() || !model.stored_energy_j().is_finite()
            || model.stored_energy_j() > limits.maximum_energy_j {
            return Err(invalid("two-sided lowering changed the valve load or exceeded initial storage"));
        }
        poll(gate)?;
        Ok(model)
    }
    /// Read-only valve including its original material/geometry and accepted state.
    #[must_use]
    pub const fn valve(&self) -> &LoadedAperture { &self.valve }
    /// Exact source data before distributed-loss subdivision or terminal controls.
    #[must_use]
    pub fn source(&self, side: ValveSide) -> &ApertureDuctSpec { &self.sides[side.index()].source }
    /// Actual graph, including lowered loss sites and current terminal controls.
    #[must_use]
    pub fn graph(&self, side: ValveSide) -> &TubeNetworkSpec { &self.sides[side.index()].graph }
    /// Original-owner represented lengths, delays and characteristic impedances.
    #[must_use]
    pub fn represented_sections(&self, side: ValveSide) -> &[SectionRealization] {
        &self.sides[side.index()].represented
    }
    /// Gas states in realized section order, not a substitute average atmosphere.
    #[must_use]
    pub fn section_gases(&self, side: ValveSide) -> &[GasState] { &self.sides[side.index()].gases }
    /// Existing material-gas loss spectra and their approximation diagnostics.
    #[must_use]
    pub fn viscothermal_sections(&self, side: ValveSide) -> &[ViscothermalSection] {
        &self.sides[side.index()].losses
    }
    /// Last accepted physical node pressure and flow; zero before the first step.
    #[must_use]
    pub fn node_frame(&self, side: ValveSide, node: usize) -> Option<&NodeFrame> {
        self.sides[side.index()].network.node_frame(node)
    }
    /// Current stored energy, evaluated from the original state owners.
    #[must_use]
    pub fn stored_energy_j(&self) -> f64 {
        self.valve.energy() + self.sides[0].network.stored_energy_j()
            + self.sides[1].network.stored_energy_j()
    }
    /// Change a passive memoryless load without clearing either acoustic history.
    /// Reactive loads refuse; this is not moving-pad mechanics.
    pub fn set_terminal_reflection(&mut self, side: ValveSide, node: usize, reflection: f64)
        -> Result<(), AcousticRealizeError>
    {
        let side = &mut self.sides[side.index()];
        side.network.set_terminal_reflection(node, reflection).map_err(wave)?;
        side.graph.nodes[node] = NetworkNode::Termination { reflection };
        Ok(())
    }
    /// Extend the shared lifetime horizon without clearing valve or wave history.
    pub fn extend_step_budget(&mut self, total: u64) -> Result<(), AcousticRealizeError> {
        match &mut self.valve {
            LoadedAperture::Single(a) => a.extend_step_budget(total),
            LoadedAperture::Paired(a) => a.extend_step_budget(total),
        }
    }

    /// Solve both acoustic loads and the original moving valve at the same step.
    /// Invalid drives, budget/physics failures and cancellation leave ALL accepted
    /// states/observations unchanged. The network owner's paired commit performs
    /// every fallible wave check before it publishes either side.
    pub fn step_under_gate(&mut self, drive: TwoSidedDrive, gate: &CancelGate)
        -> Result<TwoSidedFrame, AcousticRealizeError>
    {
        poll(gate)?;
        if !drive.source_pressure_pa.is_finite()
            || drive.blade_forces_n.iter().any(|f| !f.is_finite()) {
            return Err(invalid("two-sided pressure and force supplies must be finite"));
        }
        let before = self.stored_energy_j();
        let a = self.sides.each_ref().map(|s| s.network.incoming_pressure_pa());
        let z = self.sides.each_ref().map(|s| s.network.inlet_impedance_pa_s_m3());
        let effective_source = drive.source_pressure_pa + 2.0 * a[0];
        let trial = self.valve.preview(ApertureDrive {
            upstream_pressure_pa: effective_source, incoming_pressure_pa: a[1], body_flow_m3_s: 0.0,
        }, drive.blade_forces_n)?;
        let mut aperture = trial.frame();
        let q = aperture.bore_flow_m3_s;
        let outgoing = [a[0] - z[0] * q, a[1] + z[1] * q];
        let up = self.sides[0].network.preview_step(outgoing[0]).map_err(wave)?;
        let down = self.sides[1].network.preview_step(outgoing[1]).map_err(wave)?;
        let upstream_pressure = drive.source_pressure_pa + up.inlet_pressure_pa;
        let dp = upstream_pressure - down.inlet_pressure_pa;
        // Only observation coordinates change. The immutable candidate still owns
        // the original mechanics/material histories until all participants pass.
        aperture.outgoing_pressure_pa = outgoing[1];
        aperture.bore_pressure_pa = down.inlet_pressure_pa;
        let dt = self.valve.time_step_s();
        let frame = TwoSidedFrame {
            aperture, upstream_pressure_pa: upstream_pressure, pressure_drop_pa: dp,
            upstream: up, downstream: down,
            stored_energy_j: aperture.stored_energy_j + up.stored_energy_j + down.stored_energy_j,
            storage_change_j: aperture.storage_change_j + up.storage_change_j + down.storage_change_j,
            dissipated_energy_j: aperture.dissipated_energy_j + up.terminal_loss_j + up.interior_loss_j
                + down.terminal_loss_j + down.interior_loss_j,
            source_work_j: drive.source_pressure_pa * (q * dt),
        };
        let flow_scale = q.abs() + up.inlet_flow_m3_s.abs() + down.inlet_flow_m3_s.abs()
            + aperture.jet_flow_m3_s.abs() + aperture.swept_flow_m3_s.abs();
        let flow_allowance = self.limits.flow_absolute_tolerance_m3_s
            + self.limits.flow_relative_tolerance * flow_scale;
        let scale = before + frame.stored_energy_j + frame.dissipated_energy_j.abs()
            + frame.source_work_j.abs() + aperture.mechanical_work_j.abs();
        let allowance = self.limits.energy_absolute_tolerance_j + self.limits.energy_relative_tolerance * scale;
        if ![upstream_pressure, dp, frame.stored_energy_j, frame.storage_change_j,
            frame.dissipated_energy_j, frame.source_work_j, frame.balance_residual_j(),
            flow_scale, flow_allowance, scale, allowance].iter().all(|x| x.is_finite())
            || frame.stored_energy_j > self.limits.maximum_energy_j
            || frame.balance_residual_j().abs() > allowance
            || (up.inlet_flow_m3_s + q).abs() > flow_allowance
            || (down.inlet_flow_m3_s - q).abs() > flow_allowance
            || aperture.flow_residual_m3_s.abs() > flow_allowance
        { return Err(invalid("two-sided valve refused the complete storage, work or flow balance")); }
        poll(gate)?;
        let [upstream, downstream] = &mut self.sides;
        upstream.network.step_pair(&mut downstream.network, outgoing).map_err(wave)?;
        self.valve.accept(trial);
        Ok(frame)
    }
    /// Caller-sized blocks with exact accepted prefixes on cancellation/exhaustion.
    /// A physical refusal preserves the failing sample; prior samples remain accepted.
    pub fn advance_block(&mut self, drives: &[TwoSidedDrive], out: &mut [TwoSidedFrame], gate: &CancelGate)
        -> Result<ApertureProgress, AcousticRealizeError>
    {
        if drives.len() != out.len() { return Err(invalid("two-sided input/output lengths must match")); }
        for (completed, (drive, slot)) in drives.iter().zip(out).enumerate() {
            if gate.is_requested() {
                return Ok(ApertureProgress { completed, terminal: ApertureTerminal::Cancelled });
            }
            if self.valve.accepted_steps() >= self.valve.max_steps() {
                return Ok(ApertureProgress { completed, terminal: ApertureTerminal::BudgetExhausted });
            }
            match self.step_under_gate(*drive, gate) {
                Ok(frame) => *slot = frame,
                Err(AcousticRealizeError::Cancelled) => {
                    return Ok(ApertureProgress { completed, terminal: ApertureTerminal::Cancelled });
                }
                Err(error) => return Err(error),
            }
        }
        Ok(ApertureProgress { completed: drives.len(), terminal: ApertureTerminal::Complete })
    }
}
