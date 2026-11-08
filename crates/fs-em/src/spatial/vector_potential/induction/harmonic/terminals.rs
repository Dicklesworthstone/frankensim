//! Native AC circuits on resolved solid-conductor electrodes.
//!
//! Each non-reference electrode supplies one current coordinate, returning
//! through its own conducting component's original reference electrode. The
//! complete complex impedance includes magnetic coupling between every port;
//! prescribed motion may make it nonreciprocal or active. Only this bounded
//! interface is dense. Every unit response and actual drive uses the original
//! sparse A-phi owner, source materials and physical acceptance checks.
//!
//! The shared circuit descriptor retains RLC storage, transformers, sources and
//! Kirchhoff constraints. Its numerical reference resistors cancel before any
//! physical heat is reported. A physical electrode has exactly one graph node,
//! including references shared by several ports of a multi-electrode body.
//!
//! This is a homogeneous settled AC law V = Z I: independent impressed-current
//! excitation is absent. Material state, geometry and steady Eulerian velocity
//! are fixed. Neither transient history nor thermal/mechanical state advances.

use super::super::terminals::current::SpatialCurrentInduction;
use super::{SpatialHarmonicInductionControls, SpatialHarmonicInductionResponse};
use crate::planar::eddy::harmonic::circuit::{
    HarmonicCircuitField, HarmonicCircuitResponse, MAX_HARMONIC_CIRCUIT_WINDINGS, circuit_response,
};
use crate::spatial::SpatialMagneticError;
use core::mem::size_of;
use fs_material::conductor::{ResolvedConductorState, spatial::SpatialConductorMaterial};
use fs_math::c64::C64;
use fs_phs::circuit::{Branch, CircuitGraph, MutualInductance};
use fs_solver::SolveReport;

pub use super::circuit::{HarmonicCircuitControls, HarmonicCircuitError};

/// Native field allowances and additional cached port-response storage.
#[derive(Clone, Copy, Debug)]
pub struct SpatialHarmonicTerminalControls {
    pub field: SpatialHarmonicInductionControls,
    /// Number of independent electrode-current coordinates, from 1 through 16.
    pub maximum_ports: usize,
    /// Additional cached cochains, velocity, zero source and interface vectors
    /// in real scalar slots. Native field and circuit budgets apply separately.
    pub maximum_response_scalars: usize,
    /// Positive real current used for each physical unit-response trial [A].
    pub current_probe_a: f64,
}

#[derive(Debug)]
pub enum SpatialHarmonicTerminalError {
    Invalid(&'static str),
    Field(SpatialMagneticError),
    Circuit(HarmonicCircuitError),
    Cancelled,
}
impl core::fmt::Display for SpatialHarmonicTerminalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial harmonic electrodes: {self:?}")
    }
}
impl std::error::Error for SpatialHarmonicTerminalError {}
impl From<SpatialMagneticError> for SpatialHarmonicTerminalError {
    fn from(error: SpatialMagneticError) -> Self {
        match error {
            SpatialMagneticError::Cancelled => Self::Cancelled,
            other => Self::Field(other),
        }
    }
}
impl From<HarmonicCircuitError> for SpatialHarmonicTerminalError {
    fn from(error: HarmonicCircuitError) -> Self {
        match error {
            HarmonicCircuitError::Cancelled => Self::Cancelled,
            other => Self::Circuit(other),
        }
    }
}
type Result<T> = core::result::Result<T, SpatialHarmonicTerminalError>;
fn invalid(message: &'static str) -> SpatialHarmonicTerminalError {
    SpatialHarmonicTerminalError::Invalid(message)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(SpatialHarmonicTerminalError::Cancelled)
    } else {
        Ok(())
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("unrepresentable harmonic electrode arithmetic"))
    }
}
fn phasor(value: C64) -> Result<C64> {
    finite(value.re)?;
    finite(value.im)?;
    Ok(value)
}

/// Positive current enters terminal and returns through reference_terminal.
/// Voltage is the physical first-minus-reference difference [V].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SpatialHarmonicTerminalPort {
    pub terminal: usize,
    pub reference_terminal: usize,
}

/// Immutable frequency/material/motion-specific port map and original cochains.
pub struct SpatialHarmonicTerminalFrequency<
    's,
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    source: &'s SpatialCurrentInduction<'m, C>,
    omega_rad_s: f64,
    velocity_m_s: Vec<[f64; 3]>,
    controls: SpatialHarmonicTerminalControls,
    ports: Vec<SpatialHarmonicTerminalPort>,
    impedance_ohm: Vec<C64>,
    solutions: Vec<Vec<f64>>,
    reports: Vec<SolveReport>,
    zero_impressed: Vec<[C64; 3]>,
}

/// Accepted original field with independently recovered physical port values.
pub struct SpatialHarmonicTerminalResponse<
    's,
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    pub field: SpatialHarmonicInductionResponse<'s, 'm, C>,
    /// Recovered currents, not a copy of the requested circuit amplitudes [A].
    pub currents_a: Vec<C64>,
    pub terminal_voltages_v: Vec<C64>,
}

pub type SpatialHarmonicTerminalCircuitResponse<'s, 'm, C = ResolvedConductorState> =
    HarmonicCircuitResponse<SpatialHarmonicTerminalResponse<'s, 'm, C>>;

impl<C: SpatialConductorMaterial> HarmonicCircuitField
    for SpatialHarmonicTerminalResponse<'_, '_, C>
{
    fn currents_a(&self) -> &[C64] {
        &self.currents_a
    }
    fn terminal_voltages_v(&self) -> &[C64] {
        &self.terminal_voltages_v
    }
    fn mean_dissipated_power_w(&self) -> f64 {
        self.field.balance().mean_ohmic_power_w
    }
    fn mean_mechanical_power_w(&self) -> f64 {
        self.field.balance().mean_mechanical_power_w
    }
    fn absorbed_power_va(&self) -> C64 {
        let balance = self.field.balance();
        C64::new(
            balance.mean_terminal_power_w,
            balance.reactive_terminal_power_var,
        )
    }
}

impl<'s, 'm, C: SpatialConductorMaterial> SpatialHarmonicTerminalFrequency<'s, 'm, C> {
    /// Prepare one native field trial per non-reference electrode, in original
    /// electrode order. Multi-electrode components share their original return;
    /// magnetically coupled separate components retain the complete mutual Z.
    /// No symmetry or positive-real assumption is imposed under motion.
    pub fn new(
        source: &'s SpatialCurrentInduction<'m, C>,
        omega_rad_s: f64,
        velocity_m_s: &[[f64; 3]],
        controls: SpatialHarmonicTerminalControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        let mesh = source.problem().mesh();
        let references = source.electrode_references();
        let count = references
            .iter()
            .enumerate()
            .filter(|&(i, r)| i != *r)
            .count();
        if !positive(omega_rad_s)
            || !positive(controls.current_probe_a)
            || controls.maximum_ports == 0
            || controls.maximum_ports > MAX_HARMONIC_CIRCUIT_WINDINGS
            || count == 0
            || count > controls.maximum_ports
            || velocity_m_s.len() != mesh.vertex_count
            || velocity_m_s.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "bounded electrode ports, positive frequency/probe and complete motion required",
            ));
        }
        let coordinates = source.harmonic_coordinate_count();
        let scalars = coordinates
            .checked_mul(2)
            .and_then(|v| v.checked_mul(count + 1))
            .and_then(|v| {
                count
                    .checked_mul(count)
                    .and_then(|n| n.checked_mul(2))
                    .and_then(|n| v.checked_add(n))
            })
            .and_then(|v| {
                mesh.vertex_count
                    .checked_mul(3)
                    .and_then(|n| v.checked_add(n))
            })
            .and_then(|v| {
                mesh.tets
                    .len()
                    .checked_mul(6)
                    .and_then(|n| v.checked_add(n))
            })
            .and_then(|v| {
                references
                    .len()
                    .checked_mul(12)
                    .and_then(|n| v.checked_add(n))
            })
            .and_then(|v| {
                count
                    .checked_mul(32 + size_of::<SolveReport>().div_ceil(size_of::<f64>()))
                    .and_then(|n| v.checked_add(n))
            })
            .and_then(|v| {
                // Every cached native report retains its restart history.
                controls
                    .field
                    .maximum_cycles
                    .checked_mul(8)
                    .and_then(|n| n.checked_mul(count))
                    .and_then(|n| v.checked_add(n))
            });
        if scalars.is_none_or(|n| n > controls.maximum_response_scalars) {
            return Err(invalid(
                "harmonic electrode response storage allowance exceeded",
            ));
        }
        let ports: Vec<_> = references
            .iter()
            .enumerate()
            .filter_map(|(terminal, &reference_terminal)| {
                (terminal != reference_terminal).then_some(SpatialHarmonicTerminalPort {
                    terminal,
                    reference_terminal,
                })
            })
            .collect();
        let zero_impressed = vec![[C64::ZERO; 3]; mesh.tets.len()];
        let mut impedance_ohm = vec![C64::ZERO; count * count];
        let mut solutions = Vec::with_capacity(count);
        let mut reports = Vec::with_capacity(count);
        let probe = controls.current_probe_a;
        for (column, port) in ports.iter().enumerate() {
            poll(cancelled)?;
            let mut currents = vec![C64::ZERO; references.len()];
            currents[port.terminal] = C64::from_re(probe);
            currents[port.reference_terminal] = C64::from_re(-probe);
            let field = source.solve_harmonic_current(
                &currents,
                &zero_impressed,
                velocity_m_s,
                omega_rad_s,
                controls.field,
                cancelled,
            )?;
            for (row, receiver) in ports.iter().enumerate() {
                poll(cancelled)?;
                let voltage = field.terminals()[receiver.terminal].voltage_to_reference_v;
                impedance_ohm[row * count + column] =
                    phasor(C64::new(voltage.re / probe, voltage.im / probe))?;
            }
            solutions.push(field.solution);
            reports.push(field.report);
        }
        poll(cancelled)?;
        Ok(Self {
            source,
            omega_rad_s,
            velocity_m_s: velocity_m_s.to_vec(),
            controls,
            ports,
            impedance_ohm,
            solutions,
            reports,
            zero_impressed,
        })
    }

    pub fn source(&self) -> &'s SpatialCurrentInduction<'m, C> {
        self.source
    }
    pub fn omega_rad_s(&self) -> f64 {
        self.omega_rad_s
    }
    pub fn ports(&self) -> &[SpatialHarmonicTerminalPort] {
        &self.ports
    }
    /// Row-major coupled physical port impedance [ohm], using exp(-i omega t).
    pub fn impedance_ohm(&self) -> &[C64] {
        &self.impedance_ohm
    }
    pub fn unit_reports(&self) -> &[SolveReport] {
        &self.reports
    }

    /// Recover one actual-amplitude native field. A linear combination of
    /// cached cochains seeds its original solver, which can skip Krylov work
    /// when that seed already meets the requested residual. Full charge,
    /// Ampere, skin/motion and real/reactive-power checks always run.
    pub fn current_response(
        &self,
        currents_a: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHarmonicTerminalResponse<'s, 'm, C>> {
        poll(cancelled)?;
        if currents_a.len() != self.ports.len()
            || currents_a
                .iter()
                .any(|v| !v.re.is_finite() || !v.im.is_finite())
        {
            return Err(invalid(
                "one finite peak current per physical electrode port required",
            ));
        }
        let dimension = self.solutions[0].len() / 2;
        let mut seed = vec![0.0; 2 * dimension];
        let mut currents = vec![C64::ZERO; self.source.electrode_faces().len()];
        for (column, (port, &current)) in self.ports.iter().zip(currents_a).enumerate() {
            poll(cancelled)?;
            currents[port.terminal] = current;
            currents[port.reference_terminal] =
                phasor(currents[port.reference_terminal] - current)?;
            let amplitude = phasor(C64::new(
                current.re / self.controls.current_probe_a,
                current.im / self.controls.current_probe_a,
            ))?;
            for coordinate in 0..dimension {
                poll(cancelled)?;
                let contribution = phasor(
                    C64::new(
                        self.solutions[column][coordinate],
                        self.solutions[column][dimension + coordinate],
                    ) * amplitude,
                )?;
                seed[coordinate] = finite(seed[coordinate] + contribution.re)?;
                seed[dimension + coordinate] =
                    finite(seed[dimension + coordinate] + contribution.im)?;
            }
        }
        let field = self.source.solve_harmonic_current_seeded(
            &currents,
            &self.zero_impressed,
            &self.velocity_m_s,
            self.omega_rad_s,
            self.controls.field,
            Some(&seed),
            cancelled,
        )?;
        let mut recovered_currents = Vec::with_capacity(self.ports.len());
        let mut terminal_voltages_v = Vec::with_capacity(self.ports.len());
        for port in &self.ports {
            poll(cancelled)?;
            let transfer = &field.terminals()[port.terminal];
            recovered_currents.push(transfer.incoming_current_a);
            terminal_voltages_v.push(transfer.voltage_to_reference_v);
        }
        poll(cancelled)?;
        Ok(SpatialHarmonicTerminalResponse {
            field,
            currents_a: recovered_currents,
            terminal_voltages_v,
        })
    }

    /// Close the original RLC/transformer/source graph on every physical port.
    /// `terminal_branches` names distinct inductor placeholders in ports() order.
    /// Each branch runs from its non-reference electrode to that component's
    /// reference. All occurrences of one physical electrode must share one
    /// graph node; references cannot silently split across circuit nodes.
    ///
    /// Original descriptor assembly, terminal closure and final circuit replay
    /// are shared with native harmonic windings. Only physical graph resistance
    /// and recovered conductor heat dissipate power; prescribed motion remains
    /// signed, and numerical reference terminations contribute no physical heat.
    #[allow(clippy::too_many_arguments)]
    pub fn circuit_response(
        &self,
        graph: &CircuitGraph,
        terminal_branches: &[usize],
        external_mutual: &[MutualInductance],
        source_phasors: &[C64],
        controls: HarmonicCircuitControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHarmonicTerminalCircuitResponse<'s, 'm, C>> {
        poll(cancelled)?;
        if terminal_branches.len() != self.ports.len() {
            return Err(invalid(
                "one original circuit branch per physical electrode port required",
            ));
        }
        let mut nodes = vec![None; self.source.electrode_faces().len()];
        for (&branch, port) in terminal_branches.iter().zip(&self.ports) {
            poll(cancelled)?;
            let &(a, b, element) = graph
                .branches
                .get(branch)
                .ok_or_else(|| invalid("electrode branch is outside the original graph"))?;
            if !matches!(element, Branch::Inductor { .. }) {
                return Err(invalid(
                    "physical electrode ports require native inductor placeholders",
                ));
            }
            for (electrode, node) in [(port.terminal, a), (port.reference_terminal, b)] {
                if nodes[electrode].is_some_and(|previous| previous != node) {
                    return Err(invalid(
                        "one physical electrode cannot represent different circuit nodes",
                    ));
                }
                nodes[electrode] = Some(node);
            }
        }
        circuit_response(
            self.omega_rad_s,
            &self.impedance_ohm,
            graph,
            terminal_branches,
            external_mutual,
            source_phasors,
            controls,
            cancelled,
            |currents, cancel| self.current_response(currents, cancel),
        )
    }
}
