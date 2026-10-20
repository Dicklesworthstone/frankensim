//! Independent blade mechanics on the existing regional acoustic graph.
//!
//! TubeNetworkSpec and fs-vfit retain geometry lowering, simultaneous junctions,
//! reactive loads and traveling-wave storage. Viscothermal subdivision reuses
//! the same gas-derived loss owner as the single-reed and flue paths. Only the
//! mechanical inlet changes. Both blades and every wave/load commit together.
//! No oscillator pitch, artificial pressure output, gas-name switch, or second
//! bore integrator is introduced. Piecewise gas states are frozen at a common
//! static pressure; this is not gas mixing, advection or thermal evolution.

use super::dynamic::{ApertureDrive, ApertureProgress, ApertureTerminal};
use super::network::{ApertureNetworkFrame, NetworkNode, NodeFrame, SectionRealization, TubeNetworkSpec};
use super::network::boundaries::{BoundaryAperture, BoundaryApertures, BoundaryDrive};
use super::network::keyed;
use super::paired::{PairedAperture, PairedApertureFrame};
use super::tube::TubeDrive;
use super::viscothermal::{ViscothermalSection, ViscothermalSectionSpec, with_regional_viscothermal_sections};
use crate::acoustic_realize::AcousticRealizeError;
use crate::pcm_wav::baffled::{BaffledPressure, CircularOutletReceiver, RayleighMedium};
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_vfit::waveguide::network::WaveguideNetwork;

/// One jointly accepted sample, including both blade momentum diagnostics.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct PairedNetworkFrame {
    /// Individual blade motion and local pressure/contact work.
    pub mechanics: PairedApertureFrame,
    /// Existing combined mechanical/wave/load accounting, in joules.
    pub coupled: ApertureNetworkFrame,
}

/// The original source declarations and the retained physical realization.
pub struct PairedApertureNetwork {
    aperture: PairedAperture,
    source_graph: TubeNetworkSpec,
    source_gases: Vec<GasState>,
    graph: TubeNetworkSpec,
    gases: Vec<GasState>,
    represented: Vec<SectionRealization>,
    losses: Vec<ViscothermalSection>,
    network: WaveguideNetwork,
    boundaries: Option<BoundaryApertures>,
}
fn invalid(what: &'static str) -> AcousticRealizeError { AcousticRealizeError::InvalidDescription { what } }
fn wave(e: fs_vfit::waveguide::WaveguideError) -> AcousticRealizeError { AcousticRealizeError::Nonlinear(e.to_string()) }
fn poll(gate: Option<&CancelGate>) -> Result<(), AcousticRealizeError> {
    if gate.is_some_and(|g| g.is_requested()) { Err(AcousticRealizeError::Cancelled) } else { Ok(()) }
}

impl PairedApertureNetwork {
    /// Bind quiescent waves to a not-yet-advanced paired source. Initial blade
    /// vibration is preserved. There are no distributed wall losses unless
    /// explicitly selected through `with_viscothermal`.
    pub fn new(aperture: PairedAperture, graph: TubeNetworkSpec, gases: Vec<GasState>, gate: &CancelGate)
        -> Result<Self, AcousticRealizeError>
    {
        Self::with_viscothermal(aperture,graph,gases,&[],gate)
    }

    /// Use the original geometry/gas loss lowering. Original node indices survive
    /// section subdivision, so tone holes, reactive loads and outlet observations
    /// remain tied to the same physical nodes. Unsupported loss bandwidths refuse.
    /// Source input copies are retained separately; they are not included in the
    /// wave payload budget, which has the original owner's documented scope.
    pub fn with_viscothermal(aperture: PairedAperture, graph: TubeNetworkSpec,
        gases: Vec<GasState>, selections: &[ViscothermalSectionSpec], gate: &CancelGate)
        -> Result<Self, AcousticRealizeError>
    {
        poll(Some(gate))?;
        if aperture.accepted_steps()!=0 { return Err(invalid("a running paired source cannot acquire empty tube waves")); }
        graph.validate_section_gases(&gases)?;
        let inlet=graph.inlet_section_index()?;
        let a=aperture.spec();
        if gases[inlet].density.to_bits()!=a.density_kg_m3.to_bits()
            || graph.inlet_impedance_with_gases(&gases)?.to_bits()!=a.impedance_pa_s_m3.to_bits()
            || aperture.gas().is_some_and(|g| *g!=gases[inlet]) {
            return Err(invalid("paired pressure faces and physical inlet must use the same gas and impedance"));
        }
        let (lowered,actual_gases,losses)=with_regional_viscothermal_sections(
            graph.clone(),gases.clone(),a.time_step_s,selections,gate)?;
        let (network,represented)=lowered.realize_waveguide(&actual_gases,a.time_step_s,gate)?;
        poll(Some(gate))?;
        Ok(Self {aperture,source_graph:graph,source_gases:gases,graph:lowered,gases:actual_gases,
            represented,losses,network,boundaries:None})
    }
    /// Bind the same material/contact key owner used by single reeds and flues.
    /// Keys use the actual gas and geometry after viscothermal subdivision, with
    /// original node addresses and the paired source's complete sample horizon.
    pub fn with_boundary_apertures(mut self,inputs:Vec<BoundaryAperture>,max_apertures:usize)
        -> Result<Self,AcousticRealizeError>
    {
        if self.aperture.accepted_steps()!=0 || self.boundaries.is_some() {
            return Err(invalid("paired moving keys must bind once before any accepted sample"));
        }
        keyed::admit_horizon(&inputs,max_apertures,self.aperture.spec().max_steps)?;
        self.boundaries=Some(BoundaryApertures::new(&self.graph,&self.gases,&self.network,inputs,max_apertures)?);
        Ok(self)
    }
    /// Original material/contact state at each explicitly bound key.
    #[must_use]
    pub const fn boundary_apertures(&self) -> Option<&BoundaryApertures> { self.boundaries.as_ref() }
    /// Original independent solid specimens and accepted blade states.
    #[must_use]
    pub const fn aperture(&self) -> &PairedAperture { &self.aperture }
    /// Original graph, before selected wall-loss subdivision or terminal controls.
    #[must_use]
    pub const fn source_graph(&self) -> &TubeNetworkSpec { &self.source_graph }
    /// Complete original regional gas states, in source-section order.
    #[must_use]
    pub fn source_gases(&self) -> &[GasState] { &self.source_gases }
    /// Actual physical graph and current memoryless terminal declarations.
    #[must_use]
    pub const fn graph(&self) -> &TubeNetworkSpec { &self.graph }
    /// Actual gas for every realized section, including loss-cell subdivisions.
    #[must_use]
    pub fn section_gases(&self) -> &[GasState] { &self.gases }
    /// Represented lengths, integer transits and physical impedances.
    #[must_use]
    pub fn represented_sections(&self) -> &[SectionRealization] { &self.represented }
    /// Source-derived loss spectra and original approximation diagnostics.
    #[must_use]
    pub fn viscothermal_sections(&self) -> &[ViscothermalSection] { &self.losses }
    /// Last completed node pressure/flow, initially zero. These are not external microphones.
    #[must_use]
    pub fn node_frame(&self,node:usize) -> Option<&NodeFrame> { self.network.node_frame(node) }
    /// Actual combined mechanical/contact/wave/reactive-load storage [J].
    #[must_use]
    pub fn stored_energy_j(&self) -> f64 {
        let original=self.aperture.stored_energy_j()+self.network.stored_energy_j();
        match &self.boundaries {
            Some(boundaries)=>original+boundaries.stored_energy_j(),
            None=>original,
        }
    }
    /// Increase the complete source horizon without resetting mechanics or waves.
    pub fn extend_step_budget(&mut self,total:u64) -> Result<(),AcousticRealizeError> {
        keyed::validate_step_budget(self.boundaries.as_ref(),total)?;
        self.aperture.extend_step_budget(total)?;
        keyed::extend_step_budget(self.boundaries.as_mut(),total)
    }
    /// Change a memoryless passive terminal without erasing in-flight sound.
    /// Reactive loads refuse rather than lose their stored states.
    pub fn set_terminal_reflection(&mut self,node:usize,reflection:f64) -> Result<(),AcousticRealizeError> {
        self.network.set_terminal_reflection(node,reflection).map_err(wave)?;
        self.graph.nodes[node]=NetworkNode::Termination {reflection};
        Ok(())
    }

    /// Existing one-way, retarded baffled receiver of a physical terminal FLOW.
    /// Radius and exterior density/speed come from that outlet's actual section,
    /// never an unrelated inlet gas or a fitted audio gain. The receiver does not
    /// add radiation back-loading; select a physical load explicitly in the graph.
    /// Caller feeds accepted `net_flow_into_node_m3_s` once per complete sample.
    /// If receiver stepping refuses, discard the observed run rather than replay
    /// an already-advanced mechanical sample into the receiver.
    pub fn baffled_receiver(&self,node:usize,rate:u32,receiver:CircularOutletReceiver)
        -> Result<BaffledPressure,AcousticRealizeError>
    {
        if self.aperture.accepted_steps()!=0 || rate==0
            || (1.0/f64::from(rate)).to_bits()!=self.aperture.spec().time_step_s.to_bits() {
            return Err(invalid("new paired receiver requires quiescent flow history and the exact mechanical clock"));
        }
        if !matches!(self.graph.nodes.get(node),Some(NetworkNode::Termination {..}
            | NetworkNode::Impedance {..} | NetworkNode::Relaxation {..})) {
            return Err(invalid("paired exterior receiver requires a passive degree-one outlet"));
        }
        let mut connected=self.graph.sections.iter().enumerate().filter(|(_,s)|s.nodes.contains(&node));
        let (index,section)=connected.next().ok_or_else(||invalid("paired outlet has no section"))?;
        if connected.next().is_some() { return Err(invalid("paired exterior outlet must have exactly one section")); }
        let gas=self.gases[index];
        BaffledPressure::circular_outlet(section.radius_m,rate,receiver,
            RayleighMedium {density:gas.density,sound_speed:gas.sound_speed})
            .map_err(AcousticRealizeError::Nonlinear)
    }

    /// Advance the actual reflected-wave feedback with no additional blade forces.
    pub fn step(&mut self,drive:TubeDrive) -> Result<PairedNetworkFrame,AcousticRealizeError> {
        self.step_inner(drive,[0.0;2],&[],None)
    }
    /// Pressure plus independently prescribed opening-positive blade forces [N].
    pub fn step_with_forces(&mut self,drive:TubeDrive,forces:[f64;2])
        -> Result<PairedNetworkFrame,AcousticRealizeError> { self.step_inner(drive,forces,&[],None) }
    /// Same transaction with cancellation before trials and before publication.
    pub fn step_under_gate(&mut self,drive:TubeDrive,forces:[f64;2],gate:&CancelGate)
        -> Result<PairedNetworkFrame,AcousticRealizeError> { self.step_inner(drive,forces,&[],Some(gate)) }
    /// Advance both independent reed blades and all moving keys in one physical
    /// transaction. Every key requires an explicit ordered reservoir/force drive.
    pub fn step_with_boundaries(&mut self,drive:TubeDrive,forces:[f64;2],drives:&[BoundaryDrive])
        -> Result<PairedNetworkFrame,AcousticRealizeError> { self.step_inner(drive,forces,drives,None) }
    /// The complete keyed transaction with cancellation before publication.
    pub fn step_with_boundaries_under_gate(&mut self,drive:TubeDrive,forces:[f64;2],
        drives:&[BoundaryDrive],gate:&CancelGate) -> Result<PairedNetworkFrame,AcousticRealizeError>
    { self.step_inner(drive,forces,drives,Some(gate)) }
    fn step_inner(&mut self,drive:TubeDrive,forces:[f64;2],drives:&[BoundaryDrive],gate:Option<&CancelGate>)
        -> Result<PairedNetworkFrame,AcousticRealizeError>
    {
        poll(gate)?;
        let before=self.stored_energy_j();
        let trial=self.aperture.preview(ApertureDrive {
            upstream_pressure_pa:drive.upstream_pressure_pa,
            incoming_pressure_pa:self.network.incoming_pressure_pa(),body_flow_m3_s:drive.body_flow_m3_s,
        },forces)?;
        let mechanics=trial.frame;
        poll(gate)?;
        let aperture=mechanics.aperture;
        let network_trial=keyed::preview(aperture,drive,self.aperture.spec().time_step_s,
            &mut self.network,self.boundaries.as_ref(),drives)?;
        let coupled=network_trial.frame;
        let mut scale=before+coupled.stored_energy_j+coupled.dissipated_energy_j.abs()
            +coupled.upstream_work_j.abs()+coupled.body_work_j.abs()+aperture.mechanical_work_j.abs();
        if let Some(boundaries)=coupled.boundaries {
            scale+=boundaries.exterior_work_j.abs()+boundaries.mechanical_work_j.abs()
                +coupled.network.coupled_terminal_work_j.abs();
        }
        let allowance=256.0*self.aperture.spec().limits.relative_tolerance.max(f64::EPSILON)
            *scale.max(f64::MIN_POSITIVE);
        if ![coupled.stored_energy_j,coupled.storage_change_j,coupled.dissipated_energy_j,
            coupled.upstream_work_j,coupled.body_work_j,coupled.balance_residual_j(),scale,allowance]
            .iter().all(|x|x.is_finite()) || coupled.balance_residual_j().abs()>allowance {
            return Err(invalid("paired network refused its actual combined pressure/work balance"));
        }
        poll(gate)?;
        // Original wave owner revalidates and commits; all blade checks are
        // already complete. Its failure leaves waves and both blades unchanged.
        keyed::accept(&mut self.network,self.boundaries.as_mut(),network_trial)?;
        self.aperture.accept(trial);
        Ok(PairedNetworkFrame {mechanics,coupled})
    }
    /// Bounded caller-sized blocks with exact accepted-prefix cancellation and
    /// budget exhaustion. All output suffixes and participant clocks are retained.
    pub fn advance_block(&mut self,inputs:&[TubeDrive],out:&mut[PairedNetworkFrame],gate:&CancelGate)
        -> Result<ApertureProgress,AcousticRealizeError>
    {
        if inputs.len()!=out.len() { return Err(invalid("paired network input/output lengths must match")); }
        for (completed,(drive,slot)) in inputs.iter().zip(out.iter_mut()).enumerate() {
            if gate.is_requested() { return Ok(ApertureProgress {completed,terminal:ApertureTerminal::Cancelled}); }
            if self.aperture.accepted_steps()>=self.aperture.spec().max_steps {
                return Ok(ApertureProgress {completed,terminal:ApertureTerminal::BudgetExhausted});
            }
            *slot=self.step(*drive)?;
        }
        Ok(ApertureProgress {completed:inputs.len(),terminal:ApertureTerminal::Complete})
    }
}

#[cfg(test)]
mod tests;


/// Canonical pressure phrases and independent blade forces on the common audio pipeline.
pub mod performance;
