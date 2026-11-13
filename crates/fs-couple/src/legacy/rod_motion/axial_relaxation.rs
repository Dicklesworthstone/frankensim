//! Source-bound axial Maxwell memory on the original Cosserat trajectory.
//!
//! Each material arm stores its elastic axial strain q. Its resultant is
//! delta(EA)*q, its storage is L*delta(EA)*q²/2, and q_dot=gamma_dot-q/tau.
//! The native geodesic-midpoint axial strain supplies both gamma and its full
//! translation/rotation pullback; rigid motion cannot create branch strain.
//!
//! Symmetric exact reciprocal momentum/memory flows surround the original
//! elastic mechanical step. Each frozen-pose arm uses the existing pHS exact
//! quadratic propagator and its independently positive heat integral. Fast
//! relaxation therefore tends to its physical viscous limit, without a separate
//! explicit branch-force kick. Native elastic-wave and geometry splitting
//! accuracy still require refinement; numerical defects never become heat.

use std::sync::Arc;

use fs_exec::CancelGate;
use fs_material::section::relaxation::axial::ConcentricPronyAxialSection;
use fs_material::visco::PronyRamp;
use fs_math::det;
use fs_phs::linear_exact::{ExactZohWorkspace, LinearSplitBudget, PreparedExactZoh};
use fs_phs::{PortHamiltonian, QuadraticStorage};
use fs_solid::rod::MAX_MATERIAL_ROD_SEGMENTS;

use super::axial_loss::{
    DisplacementMode, apply_displacement_flow, apply_displacement_impulse, displacement_mode,
    displacement_rate,
};
use super::torsional_loss::same_elastic_source;
use super::{
    KickAccounts, PreparedRodMotion, RodMotion, RodMotionError, RodMotionInertia, checkpoint,
    finite,
};

#[derive(Clone, Copy, Debug)]
struct Branch {
    segment: usize,
    rigidity_n: f64,
    storage_scale_j: f64,
    root_storage: f64,
    force_gain: f64,
    relaxation_time_s: f64,
}

/// Immutable original axial spectra, in segment/layer/declaration order.
/// The native rod's equilibrium EA is retained; arms add their own delta(EA).
#[derive(Clone, Debug)]
pub struct PreparedRodAxialRelaxation {
    sections: Vec<ConcentricPronyAxialSection>,
    lengths_m: Vec<f64>,
    branches: Vec<Branch>,
}

impl PreparedRodAxialRelaxation {
    /// Bind the original scalar Young spectra to the actual reference segments.
    /// Geometry, elastic selectors, source identities and complete held queries
    /// must agree. Transient use requires explicit [0, f64::MAX] rad/s support;
    /// finite-band samples do not become a broadband hereditary law.
    pub fn new<I>(
        model: &PreparedRodMotion<I>,
        sections: Vec<ConcentricPronyAxialSection>,
        gate: &CancelGate,
    ) -> Result<Self, RodMotionError> {
        checkpoint(gate)?;
        if sections.is_empty()
            || sections.len() > MAX_MATERIAL_ROD_SEGMENTS
            || sections.len() != model.source().sections().len()
        {
            return Err(RodMotionError::Invalid(
                "axial relaxation needs matching bounded source sections",
            ));
        }
        let lengths_m = model.source().reference_mass().segment_lengths_m().to_vec();
        let mut branches = Vec::new();
        for (segment, (section, &length)) in sections.iter().zip(&lengths_m).enumerate() {
            checkpoint(gate)?;
            if section.omega_band_rad_s() != (0.0, f64::MAX) {
                return Err(RodMotionError::Invalid(
                    "transient axial relaxation requires explicit broadband source support",
                ));
            }
            for branch in section.branches() {
                let storage_scale_j = branch.axial_stiffness_n * length;
                if !storage_scale_j.is_finite()
                    || storage_scale_j < 0.0
                    || (branch.axial_stiffness_n > 0.0 && storage_scale_j == 0.0)
                {
                    return Err(RodMotionError::Invalid(
                        "axial branch storage is unrepresentable",
                    ));
                }
                let root_storage = det::sqrt(storage_scale_j);
                let force_gain = if branch.axial_stiffness_n == 0.0 {
                    0.0
                } else {
                    branch.axial_stiffness_n / root_storage
                };
                if !force_gain.is_finite() || (branch.axial_stiffness_n > 0.0 && force_gain == 0.0)
                {
                    return Err(RodMotionError::Invalid(
                        "axial branch force normalization is unrepresentable",
                    ));
                }
                branches.push(Branch {
                    segment,
                    rigidity_n: branch.axial_stiffness_n,
                    storage_scale_j,
                    root_storage,
                    force_gain,
                    relaxation_time_s: branch.relaxation_time_s,
                });
            }
        }
        let law = Self {
            sections,
            lengths_m,
            branches,
        };
        law.admit(model, gate)?;
        Ok(law)
    }

    /// Source sections, retaining every layer, spectrum and material receipt.
    pub fn sections(&self) -> &[ConcentricPronyAxialSection] {
        &self.sections
    }
    /// Number of history coordinates in segment/layer/declaration order.
    pub fn branch_count(&self) -> usize {
        self.branches.len()
    }

    fn admit<I>(
        &self,
        model: &PreparedRodMotion<I>,
        gate: &CancelGate,
    ) -> Result<(), RodMotionError> {
        if self.sections.len() != model.source().sections().len()
            || self.lengths_m != model.source().reference_mass().segment_lengths_m()
        {
            return Err(RodMotionError::Invalid(
                "axial relaxation must match the reference rod segments",
            ));
        }
        for (section, actual) in self.sections.iter().zip(model.source().sections()) {
            checkpoint(gate)?;
            if !same_elastic_source(actual.elastic(), section.elastic()) {
                return Err(RodMotionError::Invalid(
                    "axial relaxation source differs from the actual rod material state",
                ));
            }
        }
        Ok(())
    }
}

/// Independent physical relaxation heat over one accepted native step.
#[derive(Clone, Debug, PartialEq)]
pub struct RodAxialRelaxationHeat {
    pub total_heat_j: f64,
    /// Original source-segment order [J].
    pub segment_heat_j: Vec<f64>,
    /// Segment/layer/declaration order, matching the retained source arms [J].
    pub branch_heat_j: Vec<f64>,
}

#[derive(Clone)]
enum PreparedArmFlow {
    Held {
        duration_s: f64,
        ramp: PronyRamp,
    },
    Coupled {
        duration_s: f64,
        mobility: f64,
        flow: Arc<PreparedExactZoh<'static>>,
    },
}

/// Accepted axial memory. Only native motion can change this history.
/// Clones retain its strains, physical heat and the exact original sources.
#[derive(Clone)]
pub struct RodAxialRelaxationState {
    law: Arc<PreparedRodAxialRelaxation>,
    elastic_strains: Vec<f64>,
    stored_energy_j: f64,
    accumulated_heat_j: f64,
    prepared: Vec<Option<PreparedArmFlow>>,
    workspace: Option<ExactZohWorkspace>,
}

impl core::fmt::Debug for RodAxialRelaxationState {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("RodAxialRelaxationState")
            .field("law", &self.law)
            .field("elastic_strains", &self.elastic_strains)
            .field("stored_energy_j", &self.stored_energy_j)
            .field("accumulated_heat_j", &self.accumulated_heat_j)
            .finish_non_exhaustive()
    }
}

impl RodAxialRelaxationState {
    pub fn law(&self) -> &PreparedRodAxialRelaxation {
        &self.law
    }
    /// Dimensionless elastic branch strains, in source arm order.
    pub fn elastic_strains(&self) -> &[f64] {
        &self.elastic_strains
    }
    pub const fn stored_energy_j(&self) -> f64 {
        self.stored_energy_j
    }
    pub const fn accumulated_heat_j(&self) -> f64 {
        self.accumulated_heat_j
    }

    pub(super) fn heat_record(&self) -> RodAxialRelaxationHeat {
        RodAxialRelaxationHeat {
            total_heat_j: 0.0,
            segment_heat_j: vec![0.0; self.law.sections.len()],
            branch_heat_j: vec![0.0; self.law.branches.len()],
        }
    }

    fn held_ramp(&mut self, index: usize, duration_s: f64) -> PronyRamp {
        if let Some(PreparedArmFlow::Held {
            duration_s: previous,
            ramp,
        }) = &self.prepared[index]
            && *previous == duration_s
        {
            return *ramp;
        }
        let ramp = PronyRamp::new(duration_s / self.law.branches[index].relaxation_time_s);
        self.prepared[index] = Some(PreparedArmFlow::Held { duration_s, ramp });
        ramp
    }

    fn coupled_flow(
        &mut self,
        index: usize,
        duration_s: f64,
        mobility: f64,
    ) -> Result<Arc<PreparedExactZoh<'static>>, RodMotionError> {
        if let Some(PreparedArmFlow::Coupled {
            duration_s: previous,
            mobility: old,
            flow,
        }) = &self.prepared[index]
            && *previous == duration_s
            && *old == mobility
        {
            return Ok(flow.clone());
        }
        let branch = self.law.branches[index];
        // Unit-clock preparation represents the physical interval directly.
        // It avoids separately forming a potentially overflowing 1/tau.
        let coupling = (duration_s * det::sqrt(mobility)) * branch.force_gain;
        let rate = duration_s / branch.relaxation_time_s;
        finite(&[coupling, rate])?;
        if coupling <= 0.0 || rate <= 0.0 {
            return Err(RodMotionError::Invalid(
                "axial memory rate/coupling is unrepresentable",
            ));
        }
        let storage = QuadraticStorage::new(vec![1.0, 0.0, 0.0, 1.0], 2)
            .map_err(RodMotionError::AxialRelaxation)?;
        let system = PortHamiltonian::new(
            2,
            1,
            vec![0.0, -coupling, coupling, 0.0],
            vec![0.0, 0.0, 0.0, rate],
            vec![0.0, branch.force_gain],
            Box::new(storage),
        )
        .map_err(RodMotionError::AxialRelaxation)?;
        let flow = Arc::new(
            PreparedExactZoh::new_owned(
                system,
                1.0,
                LinearSplitBudget {
                    max_modes: 2,
                    max_omega_dt: 1.0,
                },
            )
            .map_err(RodMotionError::AxialRelaxation)?,
        );
        if self.workspace.is_none() {
            self.workspace = Some(flow.workspace());
        }
        self.prepared[index] = Some(PreparedArmFlow::Coupled {
            duration_s,
            mobility,
            flow: flow.clone(),
        });
        Ok(flow)
    }
}

impl<I: RodMotionInertia> RodMotion<I> {
    /// Install fully relaxed axial memory at this initial configuration.
    /// Existing equilibrium prestress is retained; no extra arm stress or stored
    /// energy is silently inferred from the current deformation.
    pub fn with_relaxed_axial_relaxation(
        self,
        law: Arc<PreparedRodAxialRelaxation>,
        gate: &CancelGate,
    ) -> Result<Self, RodMotionError> {
        self.initialize_axial_relaxation(law, None, gate)
    }

    /// Install caller-declared elastic branch strains at initialization.
    /// Values are dimensionless and ordered by segment, layer and declaration,
    /// exactly as the source sections. They describe actual initial material
    /// memory, including its storage and stress. A continued or already bound
    /// trajectory refuses instead of resetting its history or injecting energy.
    pub fn with_axial_relaxation(
        self,
        law: Arc<PreparedRodAxialRelaxation>,
        elastic_strains: &[f64],
        gate: &CancelGate,
    ) -> Result<Self, RodMotionError> {
        self.initialize_axial_relaxation(law, Some(elastic_strains), gate)
    }

    pub fn axial_relaxation_state(&self) -> Option<&RodAxialRelaxationState> {
        self.axial_relaxation.as_ref()
    }

    fn initialize_axial_relaxation(
        mut self,
        law: Arc<PreparedRodAxialRelaxation>,
        strains: Option<&[f64]>,
        gate: &CancelGate,
    ) -> Result<Self, RodMotionError> {
        checkpoint(gate)?;
        if self.steps != 0 || self.axial_relaxation.is_some() {
            return Err(RodMotionError::Invalid(
                "axial memory may only be installed once at initialization",
            ));
        }
        law.admit(&self.model, gate)?;
        if strains
            .is_some_and(|s| s.len() != law.branches.len() || s.iter().any(|v| !v.is_finite()))
        {
            return Err(RodMotionError::Invalid(
                "initial axial memory needs one finite elastic strain per source branch",
            ));
        }
        let elastic_strains =
            strains.map_or_else(|| vec![0.0; law.branch_count()], <[f64]>::to_vec);
        let mut stored_energy_j = 0.0;
        for (branch, &strain) in law.branches.iter().zip(&elastic_strains) {
            stored_energy_j += (0.5 * branch.storage_scale_j * strain) * strain;
        }
        finite(&[stored_energy_j])?;
        self.axial_relaxation = Some(RodAxialRelaxationState {
            prepared: vec![None; law.branch_count()],
            law,
            elastic_strains,
            stored_energy_j,
            accumulated_heat_j: 0.0,
            workspace: None,
        });
        self.restore_axial_relaxation_response(gate)?;
        checkpoint(gate)?;
        Ok(self)
    }

    // A fresh native response contains equilibrium forces only. B already owns
    // ALL reciprocal branch impulses and history; adding branch forces to A's
    // kicks or advecting its strain during A's drift would double-count them.
    pub(super) fn restore_axial_relaxation_response(
        &mut self,
        gate: &CancelGate,
    ) -> Result<(), RodMotionError> {
        let Some(memory) = &self.axial_relaxation else {
            return Ok(());
        };
        memory.law.admit(&self.model, gate)?;
        if memory.law.branches.is_empty() {
            return Ok(());
        }
        self.response.potential_energy_j += memory.stored_energy_j;
        finite(&[self.response.potential_energy_j])
    }

    // Each frozen-pose branch has one exact reciprocal momentum/memory flow.
    // Symmetric original-arm sweeps surround the native elastic mechanics.
    // Shared centroid/spin covectors retain eccentric inertia and supports.
    pub(super) fn flow_axial_relaxation(
        &mut self,
        dt: f64,
        origin: [f64; 3],
        accounts: &mut KickAccounts,
        heat: &mut RodAxialRelaxationHeat,
        gate: &CancelGate,
    ) -> Result<(), RodMotionError> {
        let Some(mut memory) = self.axial_relaxation.take() else {
            return Ok(());
        };
        let result = (|| {
            memory.law.admit(&self.model, gate)?;
            let n = memory.law.branches.len();
            if n == 0 {
                return Ok(());
            }
            let count = memory.law.sections.len();
            let mut gradients = [[[0.0; 6]; 2]; MAX_MATERIAL_ROD_SEGMENTS];
            self.model
                .axial_extension_gradients_into(&self.pose, &mut gradients[..count], gate)?;
            let mut modes = [DisplacementMode::ZERO; MAX_MATERIAL_ROD_SEGMENTS];
            for segment in 0..count {
                checkpoint(gate)?;
                modes[segment] = displacement_mode(self, segment, gradients[segment])?;
            }
            let mut storage_change = 0.0;
            let mut physical_heat = 0.0;
            let mut advance = |index: usize, h: f64| -> Result<(), RodMotionError> {
                checkpoint(gate)?;
                let branch = memory.law.branches[index];
                let mode = modes[branch.segment];
                let rate = displacement_rate(self, branch.segment, mode)?;
                let old = memory.elastic_strains[index];
                if old == 0.0 && rate == 0.0 {
                    return Ok(());
                }
                let (next, impulse, released, endpoint) =
                    if mode.mobility() == 0.0 || branch.rigidity_n == 0.0 {
                        let ramp = memory.held_ramp(index, h);
                        let increment = (h * rate) / memory.law.lengths_m[branch.segment];
                        finite(&[increment])?;
                        let (next, mean, integral) = ramp.apply(old, increment);
                        (
                            next,
                            (h * branch.rigidity_n) * mean,
                            branch.storage_scale_j * integral,
                            None,
                        )
                    } else {
                        let root_mobility = det::sqrt(mode.mobility());
                        let initial = [rate / root_mobility, branch.root_storage * old];
                        finite(&initial)?;
                        let flow = memory.coupled_flow(index, h, mode.mobility())?;
                        let mut next = [0.0; 2];
                        let mut output = [0.0];
                        let record = flow
                            .step_with_workspace_into(
                                memory.workspace.as_mut().expect("prepared above"),
                                &initial,
                                &[0.0],
                                &mut next,
                                &mut output,
                            )
                            .map_err(RodMotionError::AxialRelaxation)?;
                        // The original pHS port integral supplies force impulse;
                        // subtracting almost equal endpoint rates would lose it.
                        (
                            next[1] / branch.root_storage,
                            h * output[0],
                            record.dissipated,
                            Some(root_mobility * next[0]),
                        )
                    };
                finite(&[next, impulse, released])?;
                if released < 0.0 {
                    return Err(RodMotionError::Invalid("negative axial relaxation heat"));
                }
                if let Some(endpoint_rate) = endpoint {
                    apply_displacement_flow(
                        self,
                        branch.segment,
                        mode,
                        impulse,
                        endpoint_rate,
                        origin,
                        accounts,
                        gate,
                    )?;
                } else {
                    apply_displacement_impulse(
                        self,
                        branch.segment,
                        mode,
                        impulse,
                        origin,
                        accounts,
                        gate,
                    )?;
                }
                let change = (0.5 * branch.storage_scale_j * (old + next)) * (next - old);
                storage_change += change;
                physical_heat += released;
                heat.branch_heat_j[index] += released;
                heat.segment_heat_j[branch.segment] += released;
                finite(&[
                    storage_change,
                    physical_heat,
                    heat.branch_heat_j[index],
                    heat.segment_heat_j[branch.segment],
                ])?;
                memory.elastic_strains[index] = next;
                Ok(())
            };
            for index in 0..n - 1 {
                advance(index, 0.5 * dt)?;
            }
            advance(n - 1, dt)?;
            for index in (0..n - 1).rev() {
                advance(index, 0.5 * dt)?;
            }
            let mut energy = 0.0;
            for (branch, &strain) in memory.law.branches.iter().zip(&memory.elastic_strains) {
                energy += (0.5 * branch.storage_scale_j * strain) * strain;
            }
            memory.stored_energy_j = energy;
            memory.accumulated_heat_j += physical_heat;
            heat.total_heat_j += physical_heat;
            self.response.potential_energy_j += storage_change;
            finite(&[
                energy,
                memory.accumulated_heat_j,
                heat.total_heat_j,
                self.response.potential_energy_j,
            ])?;
            checkpoint(gate)
        })();
        self.axial_relaxation = Some(memory);
        result
    }
}
