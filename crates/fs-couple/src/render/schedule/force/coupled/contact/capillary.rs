//! Source-backed capillary forces driving the original coupled modal bodies.
//!
//! A complete work column maps mechanical motion into increasing plate gap;
//! its transpose applies the SAME work-conjugate mean bridge force. The
//! original fixed-volume liquid owner supplies topology, work and unresolved
//! free-energy release. Native mechanical propagation retains all modes,
//! connections and original energy/force limits. Only an accepted joint
//! mechanics/surface account publishes the bridge and body successors.
//!
//! The planar, axisymmetric and dynamic contact-line models share this SAME
//! mechanical transaction. Contact-line friction has a distinct loss port.
//! Each original liquid owner retains its geometric and numerical scope.
//! Body inertia does not add liquid inertia or dynamic contact-line physics;
//! the axisymmetric branch requires caller-established quasistatic applicability
//! and does not manufacture rupture at a branch refusal.

use super::{
    CoupledModalSystem, ModalCouplingError, dot, finite, invalid, limit, network_response, poll,
};
use crate::modal_acoustic_time::ModalAcousticTimeModel;
use fs_contact::capillary_bridge::asymmetric_contact_line::{
    AsymmetricContactLineLoadLine, AsymmetricContactLineState, AsymmetricContactLineStep,
    SourceAsymmetricContactLineBridge,
};
use fs_contact::capillary_bridge::axisymmetric::asymmetric_contact_line::{
    AsymmetricCircularContactLineLoadLine, AsymmetricCircularContactLineState,
    AsymmetricCircularContactLineStep, SourceAsymmetricCircularContactLineBridge,
};
use fs_contact::capillary_bridge::axisymmetric::contact_line::{
    AxisymmetricContactLineLoadLine, AxisymmetricContactLineState, AxisymmetricContactLineStep,
    SourceAxisymmetricContactLineBridge,
};
use fs_contact::capillary_bridge::axisymmetric::load_line::{
    AxisymmetricBridgeLoadLine, AxisymmetricLoadLineConfig, AxisymmetricLoadLineError,
};
use fs_contact::capillary_bridge::axisymmetric::{
    AxisymmetricError, AxisymmetricState, AxisymmetricStep, SourceAxisymmetricCapillaryBridge,
};
use fs_contact::capillary_bridge::contact_line::load_line::{
    ContactLineBridgeLoadLine, ContactLineLoadLineConfig, ContactLineLoadLineError,
};
use fs_contact::capillary_bridge::contact_line::{
    ContactLineError, ContactLineState, ContactLineStep, SourceContactLineBridge,
};
use fs_contact::capillary_bridge::load_line::{
    BridgeLoadLineConfig, BridgeLoadLineError, PlanarBridgeLoadLine,
};
use fs_contact::capillary_bridge::{
    BridgeError, BridgeState, BridgeStep, PlanarBridgeSourceRefresh, SourcePlanarCapillaryBridge,
};
use fs_exec::CancelGate;
use fs_matdb::MaterialStateId;

mod model;
pub use model::NativeCapillaryModel;

/// Simultaneous original circular wetting bridges on one complete native network.
pub mod multiple;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CapillaryModalConfig<C = BridgeLoadLineConfig> {
    /// Ceiling includes one final bridge proposal on actual staged motion.
    pub contact: C,
    pub force_absolute_tolerance_n: f64,
    pub force_relative_tolerance: f64,
}

#[derive(Debug)]
pub enum CapillaryModalError {
    Network(ModalCouplingError),
    Contact(BridgeLoadLineError),
    Bridge(BridgeError),
    Axisymmetric(AxisymmetricError),
    AxisymmetricContact(AxisymmetricLoadLineError),
    ContactLine(ContactLineError),
    ContactLineContact(ContactLineLoadLineError),
    GapMismatch {
        residual_m: f64,
        tolerance_m: f64,
    },
    ForceMismatch {
        residual_n: f64,
        tolerance_n: f64,
    },
    GapAccuracy {
        residual_m: f64,
        estimated_error_m: f64,
        tolerance_m: f64,
    },
    ForceAccuracy {
        residual_n: f64,
        estimated_error_n: f64,
        tolerance_n: f64,
    },
    EnergyAccuracy {
        residual_j: f64,
        estimated_error_j: f64,
        tolerance_j: f64,
    },
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
            Self::Axisymmetric(e) => Some(e),
            Self::AxisymmetricContact(e) => Some(e),
            Self::ContactLine(e) => Some(e),
            Self::ContactLineContact(e) => Some(e),
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
impl From<AxisymmetricError> for CapillaryModalError {
    fn from(value: AxisymmetricError) -> Self {
        Self::Axisymmetric(value)
    }
}
impl From<AxisymmetricLoadLineError> for CapillaryModalError {
    fn from(value: AxisymmetricLoadLineError) -> Self {
        Self::AxisymmetricContact(value)
    }
}
impl From<ContactLineError> for CapillaryModalError {
    fn from(value: ContactLineError) -> Self {
        Self::ContactLine(value)
    }
}
impl From<ContactLineLoadLineError> for CapillaryModalError {
    fn from(value: ContactLineLoadLineError) -> Self {
        Self::ContactLineContact(value)
    }
}

/// Complete accepted sample. Topology release has no automatic heat recipient.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CapillaryModalFrame<S = BridgeStep> {
    pub sample: u64,
    pub observer_pressure_pa: f64,
    pub network_energy_j: f64,
    /// Signed original surface energy relative to dry substrates.
    pub surface_energy_j: f64,
    pub external_work_j: f64,
    /// -F_mean times the caller-prescribed opening increment.
    pub prescribed_work_j: f64,
    pub network_dissipation_j: f64,
    /// Physical contact-line friction loss, separate from topology release.
    /// A thermal recipient is not inferred by the mechanical owner.
    pub contact_line_dissipation_j: f64,
    pub topology_release_j: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    /// Original numerical estimate admitted inside the combined energy budget.
    /// Zero for the analytic equilibrium planar law; never an interval claim.
    pub estimated_constitutive_energy_error_j: f64,
    pub applied_mean_force_n: f64,
    /// Independent original constitutive proposal on actual staged body motion.
    pub bridge: S,
    pub compatibility_residual_m: f64,
    pub compatibility_tolerance_m: f64,
    /// Estimated quadrature/force contribution to gap error; not an interval bound.
    pub estimated_compatibility_error_m: f64,
    /// Actual mechanical gap minus the final constitutive equilibrium gap.
    pub constitutive_gap_residual_m: f64,
    pub constitutive_residual_n: f64,
    pub force_tolerance_n: f64,
    /// Sum of both original constitutive force-error estimates, admitted within
    /// the force budget. Zero for the closed-form planar path.
    pub estimated_constitutive_force_error_n: f64,
    pub iterations: usize,
    pub evaluations: usize,
}

/// Common proposal metadata retaining the complete original bridge step.
#[doc(hidden)]
pub struct CapillaryProposal<S> {
    bridge: S,
    compatibility_tolerance_m: f64,
    estimated_compatibility_error_m: f64,
    iterations: usize,
    evaluations: usize,
}

/// Work ports read directly from the original constitutive step.
#[doc(hidden)]
pub struct CapillaryStepView<S> {
    state: S,
    gap_m: f64,
    surface_energy_j: f64,
    mean_force_n: f64,
    contact_line_dissipation_j: f64,
    topology_release_j: f64,
    estimated_energy_error_j: f64,
}

/// Native 3-D capillary mechanics with the same network and work audit as the
/// planar owner. Construct using `new_axisymmetric`.
pub type AxisymmetricCapillaryModalSystem<'a> =
    CapillaryModalSystem<'a, AxisymmetricBridgeLoadLine<'a>>;
pub type AxisymmetricCapillaryModalConfig = CapillaryModalConfig<AxisymmetricLoadLineConfig>;
pub type AxisymmetricCapillaryModalFrame = CapillaryModalFrame<AxisymmetricStep>;

/// Native finite-rate wetting with actual contact-line history and loss.
pub type ContactLineCapillaryModalSystem<'a> =
    CapillaryModalSystem<'a, ContactLineBridgeLoadLine<'a>>;
pub type ContactLineCapillaryModalConfig = CapillaryModalConfig<ContactLineLoadLineConfig>;
pub type ContactLineCapillaryModalFrame = CapillaryModalFrame<ContactLineStep>;

/// Native finite-rate wetting between two independently sourced substrates.
/// Each original footprint, pinning regime and loss remains in the frame.
pub type AsymmetricContactLineCapillaryModalSystem<'a> =
    CapillaryModalSystem<'a, AsymmetricContactLineLoadLine<'a>>;
pub type AsymmetricContactLineCapillaryModalConfig =
    CapillaryModalConfig<ContactLineLoadLineConfig>;
pub type AsymmetricContactLineCapillaryModalFrame = CapillaryModalFrame<AsymmetricContactLineStep>;

/// Native three-dimensional liquid adhesion with evolving circular footprints.
/// Both actual contact lines share the original symmetric Young--Laplace shape.
pub type AxisymmetricContactLineCapillaryModalSystem<'a> =
    CapillaryModalSystem<'a, AxisymmetricContactLineLoadLine<'a>>;
pub type AxisymmetricContactLineCapillaryModalConfig =
    CapillaryModalConfig<ContactLineLoadLineConfig>;
pub type AxisymmetricContactLineCapillaryModalFrame =
    CapillaryModalFrame<AxisymmetricContactLineStep>;

/// Three-dimensional wetting with two independently moving circular footprints.
/// The ordered original surfaces retain their own pinning and line resistance.
pub type AsymmetricCircularContactLineCapillaryModalSystem<'a> =
    CapillaryModalSystem<'a, AsymmetricCircularContactLineLoadLine<'a>>;
pub type AsymmetricCircularContactLineCapillaryModalConfig =
    CapillaryModalConfig<ContactLineLoadLineConfig>;
pub type AsymmetricCircularContactLineCapillaryModalFrame =
    CapillaryModalFrame<AsymmetricCircularContactLineStep>;

/// One fixed-direction liquid bridge and the complete unchanged native network.
pub struct CapillaryModalSystem<'a, L: NativeCapillaryModel<'a> = PlanarBridgeLoadLine<'a>> {
    network: CoupledModalSystem,
    source: &'a L::Source,
    load_line: L,
    state: L::State,
    // Preserve the mechanical coordinate independently of numerical profile
    // root residuals, so those residuals cannot accumulate as geometric drift.
    gap_m: f64,
    column: Vec<f64>,
    forces: Vec<f64>,
    config: CapillaryModalConfig<L::Config>,
    last: Option<CapillaryModalFrame<L::Step>>,
}
impl<'a> CapillaryModalSystem<'a> {
    /// Adopt an explicit material endpoint without resetting body or liquid history.
    ///
    /// The candidate must start at this exact accepted bridge state and original
    /// source. All native coordinates, connections, clocks, gap and work columns
    /// remain unchanged. The last rendered frame remains a historical sample;
    /// this parameter event does not emit a new acoustic observation.
    ///
    /// Returns signed surface free-energy parameter work [J]. The caller must
    /// retain it separately from subsequent mechanical work, dissipation and
    /// heat; there is no implicit thermal recipient or cumulative heat ledger.
    /// Source, energy or cancellation refusal leaves the complete owner unchanged.
    pub fn apply_source_refresh(
        &mut self,
        refresh: &'a PlanarBridgeSourceRefresh,
        gate: &CancelGate,
    ) -> Result<f64, CapillaryModalError> {
        poll(Some(gate))?;
        if refresh.previous_state() != self.state {
            return Err(BridgeError::StateMismatch.into());
        }
        self.source.evaluate(refresh.previous_state())?;
        let source = refresh.source();
        let state = refresh.state();
        let profile = source.evaluate(state)?;
        if state.gap_m() != self.gap_m {
            return Err(BridgeError::StateMismatch.into());
        }
        let load_line = PlanarBridgeLoadLine::new(
            source,
            self.load_line.compliance_m_per_n(),
            self.load_line.config(),
        )?;
        limit(
            "native/capillary energy magnitude after source refresh",
            finite(self.network.total_energy_j()? + profile.surface_free_energy_j)?.abs(),
            self.network.config.maximum_total_energy_j,
        )?;
        poll(Some(gate))?;
        self.source = source;
        self.state = state;
        self.load_line = load_line;
        Ok(refresh.parameter_work_j())
    }

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
        Self::from_source(
            network,
            source,
            substrate_materials,
            state,
            column,
            config,
            gate,
        )
    }
}

impl<'a> CapillaryModalSystem<'a, AxisymmetricBridgeLoadLine<'a>> {
    /// Couple the original fixed-volume 3-D bridge to actual native body motion.
    /// Its admitted branch, source accuracy and substrate identities remain in
    /// force. The caller must establish quasistatic liquid applicability.
    pub fn new_axisymmetric(
        network: CoupledModalSystem,
        source: &'a SourceAxisymmetricCapillaryBridge,
        substrate_materials: [&MaterialStateId; 2],
        state: AxisymmetricState,
        column: Vec<f64>,
        config: AxisymmetricCapillaryModalConfig,
        gate: &CancelGate,
    ) -> Result<Self, CapillaryModalError> {
        Self::from_source(
            network,
            source,
            substrate_materials,
            state,
            column,
            config,
            gate,
        )
    }
}

impl<'a> CapillaryModalSystem<'a, ContactLineBridgeLoadLine<'a>> {
    /// Bind sourced finite-rate contact-line friction/pinning to actual bodies.
    /// Original connected planar geometry and source bounds remain active.
    pub fn new_contact_line(
        network: CoupledModalSystem,
        source: &'a SourceContactLineBridge,
        substrate_materials: [&MaterialStateId; 2],
        state: ContactLineState,
        column: Vec<f64>,
        config: ContactLineCapillaryModalConfig,
        gate: &CancelGate,
    ) -> Result<Self, CapillaryModalError> {
        Self::from_source(
            network,
            source,
            substrate_materials,
            state,
            column,
            config,
            gate,
        )
    }
}

impl<'a> CapillaryModalSystem<'a, AsymmetricContactLineLoadLine<'a>> {
    /// Bind both original wetting histories to the actual native gap motion.
    /// The original connected planar geometry, independent line resistances
    /// and source bounds remain active under the shared mechanical transaction.
    pub fn new_asymmetric_contact_line(
        network: CoupledModalSystem,
        source: &'a SourceAsymmetricContactLineBridge,
        substrate_materials: [&MaterialStateId; 2],
        state: AsymmetricContactLineState,
        column: Vec<f64>,
        config: AsymmetricContactLineCapillaryModalConfig,
        gate: &CancelGate,
    ) -> Result<Self, CapillaryModalError> {
        Self::from_source(
            network,
            source,
            substrate_materials,
            state,
            column,
            config,
            gate,
        )
    }
}

impl<'a> CapillaryModalSystem<'a, AxisymmetricContactLineLoadLine<'a>> {
    /// Bind the original 3-D shape and circular-line pinning/friction to the
    /// actual mechanical gap. Both source identities and the original shape's
    /// explicit symmetric branch and accuracy limits remain active.
    pub fn new_axisymmetric_contact_line(
        network: CoupledModalSystem,
        source: &'a SourceAxisymmetricContactLineBridge,
        substrate_materials: [&MaterialStateId; 2],
        state: AxisymmetricContactLineState,
        column: Vec<f64>,
        config: AxisymmetricContactLineCapillaryModalConfig,
        gate: &CancelGate,
    ) -> Result<Self, CapillaryModalError> {
        Self::from_source(
            network,
            source,
            substrate_materials,
            state,
            column,
            config,
            gate,
        )
    }
}

impl<'a> CapillaryModalSystem<'a, AsymmetricCircularContactLineLoadLine<'a>> {
    /// Bind both original circular contact-line histories to the actual native
    /// gap. Source replay and the independent body/work checks use the shared
    /// mechanical transaction, including the original liquid shape uncertainty.
    pub fn new_asymmetric_circular_contact_line(
        network: CoupledModalSystem,
        source: &'a SourceAsymmetricCircularContactLineBridge,
        substrate_materials: [&MaterialStateId; 2],
        state: AsymmetricCircularContactLineState,
        column: Vec<f64>,
        config: AsymmetricCircularContactLineCapillaryModalConfig,
        gate: &CancelGate,
    ) -> Result<Self, CapillaryModalError> {
        Self::from_source(
            network,
            source,
            substrate_materials,
            state,
            column,
            config,
            gate,
        )
    }
}

impl<'a, L: NativeCapillaryModel<'a>> CapillaryModalSystem<'a, L> {
    fn from_source(
        network: CoupledModalSystem,
        source: &'a L::Source,
        substrate_materials: [&MaterialStateId; 2],
        state: L::State,
        column: Vec<f64>,
        config: CapillaryModalConfig<L::Config>,
        gate: &CancelGate,
    ) -> Result<Self, CapillaryModalError> {
        poll(Some(gate))?;
        if network.samples_rendered() != 0
            || column.len() != network.mode_count()
            || column.iter().any(|v| !v.is_finite())
            || L::maximum_evaluations(config.contact) < 2
            || !config.force_absolute_tolerance_n.is_finite()
            || config.force_absolute_tolerance_n <= 0.0
            || !config.force_relative_tolerance.is_finite()
            || !(0.0 < config.force_relative_tolerance && config.force_relative_tolerance < 1.0)
        {
            return Err(invalid("capillary binding needs a sample-zero network, complete finite column and positive solve/force budgets").into());
        }
        let original_materials = L::substrate_materials(source);
        for i in 0..2 {
            if original_materials[i] != substrate_materials[i] {
                return Err(invalid("capillary substrate assignment differs from the original ordered wetting sources").into());
            }
        }
        let response = network_response(&network, &column, gate)?;
        let compliance = dot(&column, &response)?;
        let load_line = L::create(source, compliance, config.contact)?;
        let (gap_m, _) = load_line.gap_energy(state)?;
        let owner = Self {
            forces: vec![0.0; network.mode_count()],
            network,
            source,
            load_line,
            state,
            gap_m,
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
    pub fn source(&self) -> &'a L::Source {
        self.source
    }
    #[must_use]
    pub const fn state(&self) -> L::State {
        self.state
    }
    /// Actual accepted mechanical separation, distinct from the source's
    /// numerically converged profile gap and its reported root residual.
    #[must_use]
    pub const fn gap_m(&self) -> f64 {
        self.gap_m
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
        self.load_line.compliance()
    }
    #[must_use]
    pub fn work_column(&self) -> &[f64] {
        &self.column
    }
    #[must_use]
    pub const fn config(&self) -> CapillaryModalConfig<L::Config> {
        self.config
    }
    #[must_use]
    pub const fn last_frame(&self) -> Option<&CapillaryModalFrame<L::Step>> {
        self.last.as_ref()
    }
    /// Signed combined storage. Surface energy retains the original dry reference.
    pub fn total_energy_j(&self) -> Result<f64, CapillaryModalError> {
        Ok(finite(
            self.network.total_energy_j()? + self.load_line.gap_energy(self.state)?.1,
        )?)
    }

    /// Hold external modal forces and prescribed extra opening velocity [m/s].
    /// Set velocity to zero for wholly internal adhesion between moving bodies.
    pub fn step(
        &mut self,
        external: &[f64],
        prescribed_gap_rate_m_s: f64,
    ) -> Result<&CapillaryModalFrame<L::Step>, CapillaryModalError> {
        let gate = CancelGate::new_clock_free();
        self.step_under_gate(external, prescribed_gap_rate_m_s, &gate)
    }
    pub fn step_under_gate(
        &mut self,
        external: &[f64],
        prescribed_gap_rate_m_s: f64,
        gate: &CancelGate,
    ) -> Result<&CapillaryModalFrame<L::Step>, CapillaryModalError> {
        poll(Some(gate))?;
        let before = self.total_energy_j()?;
        let old_gap = self.gap_m;
        let dt = self.network.dt;
        let prescribed_delta = finite(prescribed_gap_rate_m_s * dt)?;
        self.network.prepare_prediction(external, Some(gate))?;
        let free_gap =
            finite(old_gap + prescribed_delta + dot(&self.column, &self.network.free_delta)?)?;
        let solved = self.load_line.propose(self.state, free_gap, dt, gate)?;
        let solved_view = L::view(&solved.bridge);
        let force = solved_view.mean_force_n;
        self.forces.copy_from_slice(external);
        for (f, &b) in self.forces.iter_mut().zip(&self.column) {
            *f = finite(*f + b * force)?;
        }
        self.network.stage_inner(&self.forces, Some(gate))?;
        let actual_gap =
            finite(old_gap + prescribed_delta + dot(&self.column, &self.network.free_delta)?)?;
        let gap_residual = finite(actual_gap - solved_view.gap_m)?;
        if gap_residual.abs() > solved.compatibility_tolerance_m {
            return Err(CapillaryModalError::GapMismatch {
                residual_m: gap_residual,
                tolerance_m: solved.compatibility_tolerance_m,
            });
        }
        limit(
            "actual capillary gap increment",
            (actual_gap - old_gap).abs(),
            self.load_line.maximum_gap_increment(),
        )?;
        if finite(gap_residual.abs() + solved.estimated_compatibility_error_m)?
            > solved.compatibility_tolerance_m
        {
            return Err(CapillaryModalError::GapAccuracy {
                residual_m: gap_residual,
                estimated_error_m: solved.estimated_compatibility_error_m,
                tolerance_m: solved.compatibility_tolerance_m,
            });
        }
        let actual = self
            .load_line
            .advance_actual(self.state, actual_gap, dt, gate)?;
        let actual_view = L::view(&actual);
        let constitutive_gap_residual = finite(actual_gap - actual_view.gap_m)?;
        if constitutive_gap_residual.abs() > solved.compatibility_tolerance_m {
            return Err(CapillaryModalError::GapMismatch {
                residual_m: constitutive_gap_residual,
                tolerance_m: solved.compatibility_tolerance_m,
            });
        }
        let force_residual = finite(force - actual_view.mean_force_n)?;
        let force_tolerance = finite(
            self.config.force_absolute_tolerance_n
                + self.config.force_relative_tolerance
                    * force.abs().max(actual_view.mean_force_n.abs()),
        )?;
        if force_residual.abs() > force_tolerance {
            return Err(CapillaryModalError::ForceMismatch {
                residual_n: force_residual,
                tolerance_n: force_tolerance,
            });
        }
        let force_error = finite(
            self.load_line.force_error(self.state, &solved.bridge)?
                + self.load_line.force_error(self.state, &actual)?,
        )?;
        if finite(force_residual.abs() + force_error)? > force_tolerance {
            return Err(CapillaryModalError::ForceAccuracy {
                residual_n: force_residual,
                estimated_error_n: force_error,
                tolerance_n: force_tolerance,
            });
        }
        let frame = &self.network.staged_frame;
        let network_energy = finite(frame.total_energy_j())?;
        let surface_energy = actual_view.surface_energy_j;
        let after = finite(network_energy + surface_energy)?;
        limit(
            "native/capillary energy magnitude",
            after.abs(),
            self.network.config.maximum_total_energy_j,
        )?;
        let external_work = dot(external, &self.network.free_delta)?;
        let prescribed_work = finite(-force * prescribed_delta)?;
        let network_loss = finite(frame.total_dissipation_j())?;
        let contact_line_loss = finite(actual_view.contact_line_dissipation_j)?;
        if contact_line_loss < 0.0 {
            return Err(invalid("capillary contact-line loss must be nonnegative").into());
        }
        let total_loss = finite(network_loss + contact_line_loss + actual_view.topology_release_j)?;
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
        let energy_error = finite(actual_view.estimated_energy_error_j)?;
        if energy_error < 0.0 {
            return Err(invalid("capillary energy error estimate must be nonnegative").into());
        }
        if finite(residual.abs() + energy_error)? > tolerance {
            return Err(CapillaryModalError::EnergyAccuracy {
                residual_j: residual,
                estimated_error_j: energy_error,
                tolerance_j: tolerance,
            });
        }
        let accepted = CapillaryModalFrame {
            sample: frame.sample,
            observer_pressure_pa: frame.observer_pressure_pa,
            network_energy_j: network_energy,
            surface_energy_j: surface_energy,
            external_work_j: external_work,
            prescribed_work_j: prescribed_work,
            network_dissipation_j: network_loss,
            contact_line_dissipation_j: contact_line_loss,
            topology_release_j: actual_view.topology_release_j,
            energy_residual_j: residual,
            energy_tolerance_j: tolerance,
            estimated_constitutive_energy_error_j: energy_error,
            applied_mean_force_n: force,
            bridge: actual,
            compatibility_residual_m: gap_residual,
            compatibility_tolerance_m: solved.compatibility_tolerance_m,
            estimated_compatibility_error_m: solved.estimated_compatibility_error_m,
            constitutive_gap_residual_m: constitutive_gap_residual,
            constitutive_residual_n: force_residual,
            force_tolerance_n: force_tolerance,
            estimated_constitutive_force_error_n: force_error,
            iterations: solved.iterations,
            evaluations: solved.evaluations + 1,
        };
        self.network.validate_staged_material_points(Some(gate))?;
        poll(Some(gate))?;
        self.network.publish_staged();
        self.state = actual_view.state;
        self.gap_m = actual_gap;
        Ok(self.last.insert(accepted))
    }
}
