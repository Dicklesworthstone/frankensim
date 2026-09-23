//! Moving valve coupled to a geometry-defined, branched acoustic tube network.
//!
//! Each cylindrical section derives its load and integer transit from physical
//! geometry and the aperture's declared fluid. fs-vfit owns simultaneous
//! scattering and traveling-wave storage; the existing dynamic aperture owns
//! nonlinear mechanics/contact. A sample commits both or neither. Ideal branch
//! nodes have zero volume and no end corrections. Memoryless terminal controls
//! are not moving-pad mechanics, radiation loads or measured instrument models.

use super::dynamic::{ApertureDrive, ApertureFrame, ApertureProgress, ApertureTerminal, DynamicAperture};
use super::tube::{TubeDrive, UniformTubeSpec};
use crate::acoustic_realize::AcousticRealizeError;
use fs_exec::CancelGate;
use fs_vfit::waveguide::network::{NetworkFrame, NetworkSegment, WaveguideNetwork};
pub use fs_vfit::waveguide::network::{NetworkNode, NodeFrame};

/// Physical cylindrical section between two network nodes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TubeSection {
    /// Zero-based endpoint indices. Positive-length cycles are supported.
    pub nodes: [usize; 2],
    /// Requested centerline length [m].
    pub length_m: f64,
    /// Internal radius [m].
    pub radius_m: f64,
    /// Explicit maximum length change from integer transit [m].
    pub max_length_error_m: f64,
}

impl TubeSection {
    fn uniform(&self, speed: f64) -> UniformTubeSpec {
        UniformTubeSpec {
            length_m: self.length_m, radius_m: self.radius_m,
            sound_speed_m_s: speed, terminal_reflection: 0.0,
            max_length_error_m: self.max_length_error_m, max_wave_memory_bytes: 0,
        }
    }

    fn realize(&self, speed: f64, density: f64, dt: f64) -> Result<SectionRealization, AcousticRealizeError> {
        let z = self.uniform(speed).characteristic_impedance(density)?;
        if !self.max_length_error_m.is_finite() || self.max_length_error_m < 0.0 {
            return Err(invalid("each network section needs a finite nonnegative length-error allowance"));
        }
        let cell = speed * dt;
        let samples = (self.length_m / cell).round();
        if !cell.is_finite() || cell <= 0.0 || !samples.is_finite()
            || samples < 1.0 || samples >= usize::MAX as f64
        {
            return Err(invalid("network section transit must fit a positive integer delay"));
        }
        let represented_length_m = samples * cell;
        if !represented_length_m.is_finite()
            || (represented_length_m - self.length_m).abs() > self.max_length_error_m
        {
            return Err(invalid("network section transit exceeds its length-error allowance"));
        }
        Ok(SectionRealization { represented_length_m, one_way_samples: samples as usize,
            impedance_pa_s_m3: z })
    }
}

/// Requested topology, geometry, medium and whole-network payload budget.
#[derive(Debug, Clone, PartialEq)]
pub struct TubeNetworkSpec {
    /// One degree-one inlet, degree >=2 junctions, and degree-one passive loads.
    pub nodes: Vec<NetworkNode>,
    /// Fixed topology and requested cylindrical geometry, in reduction order.
    pub sections: Vec<TubeSection>,
    /// Common sound speed in the aperture's fluid [m/s].
    pub sound_speed_m_s: f64,
    /// Network payload and geometry-lowering scratch allowance. Excludes caller
    /// input vectors, aperture allocations, allocator overhead and process RSS.
    pub max_wave_memory_bytes: usize,
}

impl TubeNetworkSpec {
    /// Derive the unique inlet section's impedance for constructing the aperture.
    /// No average or substitute impedance is selected at a branched inlet.
    ///
    /// # Errors
    /// Missing/repeated inlet, inlet degree other than one or invalid geometry.
    pub fn inlet_impedance(&self, density: f64) -> Result<f64, AcousticRealizeError> {
        let mut inlet = None;
        for (n, kind) in self.nodes.iter().enumerate() {
            if matches!(kind, NetworkNode::Inlet) {
                if inlet.replace(n).is_some() { return Err(invalid("tube network requires one inlet")); }
            }
        }
        let inlet = inlet.ok_or_else(|| invalid("tube network requires one inlet"))?;
        let mut selected = None;
        for section in &self.sections {
            for &node in &section.nodes {
                if node == inlet && selected.replace(section).is_some() {
                    return Err(invalid("tube network inlet must meet exactly one section"));
                }
            }
        }
        selected.ok_or_else(|| invalid("tube network inlet is disconnected"))?
            .uniform(self.sound_speed_m_s).characteristic_impedance(density)
    }
}

/// Explicit mapping from requested geometry to the discrete section.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SectionRealization {
    /// Length actually represented [m]; compare to the same-index request.
    pub represented_length_m: f64,
    /// Positive one-way sample delay.
    pub one_way_samples: usize,
    /// Derived rho c / (pi r^2) [Pa s/m^3].
    pub impedance_pa_s_m3: f64,
}

/// One accepted valve/network sample with both external power supplies.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ApertureNetworkFrame {
    /// Midpoint valve observations and end mechanical/contact state.
    pub aperture: ApertureFrame,
    /// Same-step network ports and all-section wave storage.
    pub network: NetworkFrame,
    /// Mechanical/contact plus all traveling-wave energy [J].
    pub stored_energy_j: f64,
    /// Change in independently evaluated total storage [J].
    pub storage_change_j: f64,
    /// Valve dissipation plus all terminal and interior load dissipation [J].
    pub dissipated_energy_j: f64,
    /// P_upstream (U_inlet - U_body) dt [J].
    pub upstream_work_j: f64,
    /// P_inlet U_body dt [J], kept separate from the upstream source.
    pub body_work_j: f64,
}
impl ApertureNetworkFrame {
    /// Total balance diagnostic; neither junction nor solver error is removed.
    #[must_use]
    pub fn balance_residual_j(&self) -> f64 {
        self.storage_change_j + self.dissipated_energy_j - self.upstream_work_j - self.body_work_j
    }
}

/// Stateful branch/step feedback into the existing nonlinear moving aperture.
/// The graph is fixed; passive terminal values may change between samples.
/// Wave propagation allocates no per-step storage; aperture trials still do.
pub struct ApertureNetwork {
    aperture: DynamicAperture,
    network: WaveguideNetwork,
    spec: TubeNetworkSpec,
    represented: Vec<SectionRealization>,
}

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

impl ApertureNetwork {
    /// Bind a not-yet-advanced valve to a quiescent network. Initial mechanical
    /// motion is allowed; a running valve cannot silently acquire empty waves.
    ///
    /// # Errors
    /// Clock/load mismatch, geometry/tolerance, topology, memory or allocation.
    pub fn new(aperture: DynamicAperture, spec: TubeNetworkSpec) -> Result<Self, AcousticRealizeError> {
        if aperture.accepted_steps() != 0 {
            return Err(invalid("quiescent network coupling requires an unadvanced aperture"));
        }
        let a = aperture.spec();
        let z = spec.inlet_impedance(a.density_kg_m3)?;
        if z.to_bits() != a.impedance_pa_s_m3.to_bits() {
            return Err(invalid("aperture load must equal the physical network inlet impedance"));
        }
        let lower_bytes = spec.sections.len().checked_mul(
            core::mem::size_of::<NetworkSegment>() + core::mem::size_of::<SectionRealization>())
            .ok_or_else(|| invalid("network geometry allocation size overflow"))?;
        let network_budget = spec.max_wave_memory_bytes.checked_sub(lower_bytes)
            .ok_or_else(|| invalid("network geometry payload exceeds memory allowance"))?;
        let mut represented = Vec::new();
        let mut segments = Vec::new();
        represented.try_reserve_exact(spec.sections.len())
            .map_err(|_| invalid("network geometry allocation failed"))?;
        segments.try_reserve_exact(spec.sections.len())
            .map_err(|_| invalid("network geometry allocation failed"))?;
        for section in &spec.sections {
            let realized = section.realize(spec.sound_speed_m_s, a.density_kg_m3, a.time_step_s)?;
            represented.push(realized);
            segments.push(NetworkSegment { nodes: section.nodes,
                one_way_samples: realized.one_way_samples,
                impedance_pa_s_m3: realized.impedance_pa_s_m3 });
        }
        let network = WaveguideNetwork::new(&spec.nodes, &segments, a.time_step_s, network_budget)
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        Ok(Self { aperture, network, spec, represented })
    }

    /// Immutable accepted mechanical state, contact law and sample count.
    #[must_use]
    pub const fn aperture(&self) -> &DynamicAperture { &self.aperture }

    /// Requested geometry and current memoryless terminal values.
    #[must_use]
    pub const fn spec(&self) -> &TubeNetworkSpec { &self.spec }

    /// Actual represented section lengths, transit times and loads.
    #[must_use]
    pub fn represented_sections(&self) -> &[SectionRealization] { &self.represented }

    /// Last accepted physical node pressure and flow, initially zero.
    #[must_use]
    pub fn node_frame(&self, node: usize) -> Option<&NodeFrame> { self.network.node_frame(node) }

    /// Accepted physical R-L-C coordinates at an endpoint or interior load.
    /// No mutable access is exposed: waves and load clocks remain synchronized.
    #[must_use]
    pub fn load_state(&self, node: usize) -> Option<fs_vfit::impedance::ImpedanceState> {
        self.network.terminal_state(node)
    }

    /// Actual mechanical/contact/wave storage [J].
    #[must_use]
    pub fn stored_energy_j(&self) -> f64 {
        self.aperture.stored_energy_j() + self.network.stored_energy_j()
    }

    /// Change a passive memoryless endpoint without clearing sound already in
    /// flight. This does not model pad displacement, junction volume or radiation.
    ///
    /// # Errors
    /// Nonterminal index or nonfinite/active reflection; nothing changes.
    pub fn set_terminal_reflection(&mut self, node: usize, reflection: f64) -> Result<(), AcousticRealizeError> {
        self.network.set_terminal_reflection(node, reflection)
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        self.spec.nodes[node] = NetworkNode::Termination { reflection };
        Ok(())
    }

    /// Extend the accepted sample budget without resetting either participant.
    ///
    /// # Errors
    /// The aperture's finite, increasing-horizon admission.
    pub fn extend_step_budget(&mut self, total: u64) -> Result<(), AcousticRealizeError> {
        self.aperture.extend_step_budget(total)
    }

    /// Solve the valve against the retained incoming wave, scatter the graph,
    /// and validate the combined observation before committing either state.
    ///
    /// # Errors
    /// Input, budget, nonlinear solve, propagation or finite-set refusal. Private
    /// network scratch may change, but accepted waves/state/clock do not.
    pub fn step(&mut self, drive: TubeDrive) -> Result<ApertureNetworkFrame, AcousticRealizeError> {
        let trial = self.aperture.preview_step(ApertureDrive {
            upstream_pressure_pa: drive.upstream_pressure_pa,
            incoming_pressure_pa: self.network.incoming_pressure_pa(),
            body_flow_m3_s: drive.body_flow_m3_s,
        })?;
        let aperture = trial.frame;
        let network = self.network.preview_step(aperture.outgoing_pressure_pa)
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        let dt = self.aperture.spec().time_step_s;
        let frame = ApertureNetworkFrame {
            aperture, network,
            stored_energy_j: aperture.stored_energy_j + network.stored_energy_j,
            storage_change_j: aperture.storage_change_j + network.storage_change_j,
            dissipated_energy_j: aperture.dissipated_energy_j + network.terminal_loss_j + network.interior_loss_j,
            upstream_work_j: drive.upstream_pressure_pa
                * ((aperture.bore_flow_m3_s - drive.body_flow_m3_s) * dt),
            body_work_j: aperture.bore_pressure_pa * (drive.body_flow_m3_s * dt),
        };
        if ![frame.stored_energy_j, frame.storage_change_j, frame.dissipated_energy_j,
            frame.upstream_work_j, frame.body_work_j, frame.balance_residual_j()]
            .iter().all(|v| v.is_finite())
        {
            return Err(invalid("coupled network observation left the finite set"));
        }
        self.network.step(aperture.outgoing_pressure_pa)
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        self.aperture.accept_frame(trial);
        Ok(frame)
    }

    /// Advance a caller-sized block; cancellation or exhaustion returns an
    /// accepted prefix with untouched output suffix. Resume with that suffix and
    /// a fresh gate; explicitly extend the total budget after exhaustion.
    ///
    /// # Errors
    /// Shape or sample refusal. The accepted-step count identifies completed
    /// work, and every section remains synchronized with the valve on failure.
    pub fn advance_block(
        &mut self, inputs: &[TubeDrive], out: &mut [ApertureNetworkFrame], gate: &CancelGate,
    ) -> Result<ApertureProgress, AcousticRealizeError> {
        if inputs.len() != out.len() { return Err(invalid("network input/output block lengths must match")); }
        for (completed, (drive, slot)) in inputs.iter().zip(out.iter_mut()).enumerate() {
            if gate.is_requested() {
                return Ok(ApertureProgress { completed, terminal: ApertureTerminal::Cancelled });
            }
            if self.aperture.accepted_steps() >= self.aperture.spec().max_steps {
                return Ok(ApertureProgress { completed, terminal: ApertureTerminal::BudgetExhausted });
            }
            *slot = self.step(*drive)?;
        }
        Ok(ApertureProgress { completed: inputs.len(), terminal: ApertureTerminal::Complete })
    }
}
