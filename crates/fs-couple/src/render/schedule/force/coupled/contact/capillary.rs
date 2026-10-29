//! Source-backed capillary forces driving the original coupled modal bodies.
//!
//! A complete work column maps mechanical motion into increasing plate gap;
//! its transpose applies the SAME work-conjugate mean bridge force. The
//! original fixed-volume liquid owner supplies topology, work and unresolved
//! free-energy release. Native mechanical propagation retains all modes,
//! connections and original energy/force limits. Only an accepted joint
//! mechanics/surface account publishes the bridge and body successors.
//!
//! The bridge remains the planar quasistatic liquid model, with its own
//! capillary/Weber/geometry screens. This integrates body inertia; it does not
//! add resolved liquid inertia, three-dimensional ends or physical pinch-off.

use super::{
    CoupledModalSystem, ModalCouplingError, dot, finite, invalid, limit, network_response, poll,
};
use crate::modal_acoustic_time::ModalAcousticTimeModel;
use fs_contact::capillary_bridge::load_line::{
    BridgeLoadLineConfig, BridgeLoadLineError, PlanarBridgeLoadLine,
};
use fs_contact::capillary_bridge::{
    BridgeError, BridgeState, BridgeStep, SourcePlanarCapillaryBridge,
};
use fs_exec::CancelGate;
use fs_matdb::MaterialStateId;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CapillaryModalConfig {
    /// Ceiling includes one final bridge proposal on actual staged motion.
    pub contact: BridgeLoadLineConfig,
    pub force_absolute_tolerance_n: f64,
    pub force_relative_tolerance: f64,
}

#[derive(Debug)]
pub enum CapillaryModalError {
    Network(ModalCouplingError),
    Contact(BridgeLoadLineError),
    Bridge(BridgeError),
    GapMismatch { residual_m: f64, tolerance_m: f64 },
    ForceMismatch { residual_n: f64, tolerance_n: f64 },
}
impl core::fmt::Display for CapillaryModalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native capillary trial refused: {self:?}")
    }
}
impl std::error::Error for CapillaryModalError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Network(e) => Some(e),
            Self::Contact(e) => Some(e),
            Self::Bridge(e) => Some(e),
            _ => None,
        }
    }
}
impl From<ModalCouplingError> for CapillaryModalError {
    fn from(value: ModalCouplingError) -> Self {
        Self::Network(value)
    }
}
impl From<BridgeLoadLineError> for CapillaryModalError {
    fn from(value: BridgeLoadLineError) -> Self {
        Self::Contact(value)
    }
}
impl From<BridgeError> for CapillaryModalError {
    fn from(value: BridgeError) -> Self {
        Self::Bridge(value)
    }
}

/// Complete accepted sample. Topology release has no automatic heat recipient.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CapillaryModalFrame {
    pub sample: u64,
    pub observer_pressure_pa: f64,
    pub network_energy_j: f64,
    /// Signed original surface energy relative to dry substrates.
    pub surface_energy_j: f64,
    pub external_work_j: f64,
    /// -F_mean times the caller-prescribed opening increment.
    pub prescribed_work_j: f64,
    pub network_dissipation_j: f64,
    pub topology_release_j: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    pub applied_mean_force_n: f64,
    /// Independent original constitutive proposal on actual staged body motion.
    pub bridge: BridgeStep,
    pub compatibility_residual_m: f64,
    pub compatibility_tolerance_m: f64,
    pub constitutive_residual_n: f64,
    pub force_tolerance_n: f64,
    pub iterations: usize,
    pub evaluations: usize,
}

/// One fixed-direction liquid bridge and the complete unchanged native network.
pub struct CapillaryModalSystem<'a> {
    network: CoupledModalSystem,
    load_line: PlanarBridgeLoadLine<'a>,
    state: BridgeState,
    column: Vec<f64>,
    forces: Vec<f64>,
    config: CapillaryModalConfig,
    last: Option<CapillaryModalFrame>,
}
impl<'a> CapillaryModalSystem<'a> {
    /// Attach an explicit bridge state to the current sample-zero configuration.
    /// Its gap is the current physical separation even if initial modal q!=0.
    /// Ordered substrate identities must match the two original wetting cards.
    /// The complete column includes every translational/elastic/rotational
    /// participation in increasing gap; its reciprocal force is +B*F_bridge.
    pub fn new(
        network: CoupledModalSystem,
        source: &'a SourcePlanarCapillaryBridge,
        substrate_materials: [&MaterialStateId; 2],
        state: BridgeState,
        column: Vec<f64>,
        config: CapillaryModalConfig,
        gate: &CancelGate,
    ) -> Result<Self, CapillaryModalError> {
        poll(Some(gate))?;
        if network.samples_rendered() != 0
            || column.len() != network.mode_count()
            || column.iter().any(|v| !v.is_finite())
            || config.contact.max_evaluations < 2
            || !config.force_absolute_tolerance_n.is_finite()
            || config.force_absolute_tolerance_n <= 0.0
            || !config.force_relative_tolerance.is_finite()
            || !(0.0 < config.force_relative_tolerance && config.force_relative_tolerance < 1.0)
        {
            return Err(invalid("capillary binding needs a sample-zero network, complete finite column and positive solve/force budgets").into());
        }
        for i in 0..2 {
            if &source.wetting()[i].surfaces()[0].material != substrate_materials[i] {
                return Err(invalid("capillary substrate assignment differs from the original ordered wetting sources").into());
            }
        }
        source.evaluate(state)?;
        let response = network_response(&network, &column, gate)?;
        let compliance = dot(&column, &response)?;
        let mut solve = config.contact;
        solve.max_evaluations -= 1;
        let load_line = PlanarBridgeLoadLine::new(source, compliance, solve)?;
        let owner = Self {
            forces: vec![0.0; network.mode_count()],
            network,
            load_line,
            state,
            column,
            config,
            last: None,
        };
        limit(
            "initial native/capillary energy magnitude",
            owner.total_energy_j()?.abs(),
            owner.network.config.maximum_total_energy_j,
        )?;
        poll(Some(gate))?;
        Ok(owner)
    }
    #[must_use]
    pub fn source(&self) -> &'a SourcePlanarCapillaryBridge {
        self.load_line.source()
    }
    #[must_use]
    pub const fn state(&self) -> BridgeState {
        self.state
    }
    #[must_use]
    pub fn components(&self) -> &[ModalAcousticTimeModel] {
        self.network.components()
    }
    #[must_use]
    pub fn mode_count(&self) -> usize {
        self.network.mode_count()
    }
    #[must_use]
    pub fn sample_period_s(&self) -> f64 {
        self.network.sample_period_s()
    }
    #[must_use]
    pub fn samples_rendered(&self) -> u64 {
        self.network.samples_rendered()
    }
    #[must_use]
    pub fn compliance_m_per_n(&self) -> f64 {
        self.load_line.compliance_m_per_n()
    }
    #[must_use]
    pub fn work_column(&self) -> &[f64] {
        &self.column
    }
    #[must_use]
    pub const fn config(&self) -> CapillaryModalConfig {
        self.config
    }
    #[must_use]
    pub const fn last_frame(&self) -> Option<&CapillaryModalFrame> {
        self.last.as_ref()
    }
    /// Signed combined storage. Surface energy retains the original dry reference.
    pub fn total_energy_j(&self) -> Result<f64, CapillaryModalError> {
        Ok(finite(
            self.network.total_energy_j()?
                + self.source().evaluate(self.state)?.surface_free_energy_j,
        )?)
    }

    /// Hold external modal forces and prescribed extra opening velocity [m/s].
    /// Set velocity to zero for wholly internal adhesion between moving bodies.
    pub fn step(
        &mut self,
        external: &[f64],
        prescribed_gap_rate_m_s: f64,
    ) -> Result<&CapillaryModalFrame, CapillaryModalError> {
        let gate = CancelGate::new_clock_free();
        self.step_under_gate(external, prescribed_gap_rate_m_s, &gate)
    }
    pub fn step_under_gate(
        &mut self,
        external: &[f64],
        prescribed_gap_rate_m_s: f64,
        gate: &CancelGate,
    ) -> Result<&CapillaryModalFrame, CapillaryModalError> {
        poll(Some(gate))?;
        let source = self.load_line.source();
        let before = self.total_energy_j()?;
        let dt = self.network.dt;
        let prescribed_delta = finite(prescribed_gap_rate_m_s * dt)?;
        self.network.prepare_prediction(external, Some(gate))?;
        let free_gap = finite(
            self.state.gap_m() + prescribed_delta + dot(&self.column, &self.network.free_delta)?,
        )?;
        let solved = self.load_line.trial(self.state, free_gap, dt, gate)?;
        let force = solved.bridge.mean_normal_force_n;
        self.forces.copy_from_slice(external);
        for (f, &b) in self.forces.iter_mut().zip(&self.column) {
            *f = finite(*f + b * force)?;
        }
        self.network.stage_inner(&self.forces, Some(gate))?;
        let actual_gap = finite(
            self.state.gap_m() + prescribed_delta + dot(&self.column, &self.network.free_delta)?,
        )?;
        let gap_residual = finite(actual_gap - solved.bridge.state().gap_m())?;
        if gap_residual.abs() > solved.compatibility_tolerance_m {
            return Err(CapillaryModalError::GapMismatch {
                residual_m: gap_residual,
                tolerance_m: solved.compatibility_tolerance_m,
            });
        }
        limit(
            "actual capillary gap increment",
            (actual_gap - self.state.gap_m()).abs(),
            self.config.contact.maximum_gap_increment_m,
        )?;
        let actual = source.advance(self.state, actual_gap, dt, gate)?;
        let force_residual = finite(force - actual.mean_normal_force_n)?;
        let force_tolerance = finite(
            self.config.force_absolute_tolerance_n
                + self.config.force_relative_tolerance
                    * force.abs().max(actual.mean_normal_force_n.abs()),
        )?;
        if force_residual.abs() > force_tolerance {
            return Err(CapillaryModalError::ForceMismatch {
                residual_n: force_residual,
                tolerance_n: force_tolerance,
            });
        }
        let frame = &self.network.staged_frame;
        let network_energy = finite(frame.total_energy_j())?;
        let surface_energy = actual.profile.surface_free_energy_j;
        let after = finite(network_energy + surface_energy)?;
        limit(
            "native/capillary energy magnitude",
            after.abs(),
            self.network.config.maximum_total_energy_j,
        )?;
        let external_work = dot(external, &self.network.free_delta)?;
        let prescribed_work = finite(-force * prescribed_delta)?;
        let network_loss = finite(frame.total_dissipation_j())?;
        let total_loss = finite(network_loss + actual.irreversible_release_j)?;
        let residual = finite((after - before) + total_loss - external_work - prescribed_work)?;
        let scale = before
            .abs()
            .max(after.abs())
            .max(external_work.abs())
            .max(prescribed_work.abs())
            .max(total_loss.abs());
        let tolerance = finite(
            self.network.config.energy_absolute_tolerance_j
                + self.network.config.energy_relative_tolerance * scale,
        )?;
        if residual.abs() > tolerance {
            return Err(ModalCouplingError::EnergyBalance {
                residual_j: residual,
                tolerance_j: tolerance,
            }
            .into());
        }
        let accepted = CapillaryModalFrame {
            sample: frame.sample,
            observer_pressure_pa: frame.observer_pressure_pa,
            network_energy_j: network_energy,
            surface_energy_j: surface_energy,
            external_work_j: external_work,
            prescribed_work_j: prescribed_work,
            network_dissipation_j: network_loss,
            topology_release_j: actual.irreversible_release_j,
            energy_residual_j: residual,
            energy_tolerance_j: tolerance,
            applied_mean_force_n: force,
            bridge: actual,
            compatibility_residual_m: gap_residual,
            compatibility_tolerance_m: solved.compatibility_tolerance_m,
            constitutive_residual_n: force_residual,
            force_tolerance_n: force_tolerance,
            iterations: solved.iterations,
            evaluations: solved.evaluations + 1,
        };
        poll(Some(gate))?;
        self.network.publish_staged();
        self.state = actual.state();
        Ok(self.last.insert(accepted))
    }
}
