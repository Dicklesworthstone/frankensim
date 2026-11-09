//! Native bulk-history hardness mapped explicitly to surface contact and wear.
//!
//! An accepted field is a frozen property transfer from an independently evolved
//! or preconditioned bulk history. The caller supplies the original geometry and
//! a cell-to-material-point map; no nearest-point interpolation, residual-stress
//! transfer, or conversion from indentation depth to accumulated strain occurs.
//! Each mapped value must represent its material column through the declared
//! recession bound. A resolved bulk indentation solve and this reduced surface
//! indentation law must not own the same deformation or plastic work.
//!
//! The homogeneous elastic half-space kernel remains independent of hardness.
//! Local pressure caps use min(H_A,H_B), while each body's Archard removal uses
//! its own H. External field refresh retains the original native indentation,
//! removal inventories and signed parameter-work account at the chosen boundary.

mod binding;
mod wear;
/// Adaptive pressure feedback along prescribed local sliding paths.
pub mod integration;

use std::sync::Arc;

use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_matdb::PropertyKey;
use fs_material::hardness::IndentationHardnessStatePoint;
use fs_material::hardness::history_field::AcceptedHardnessField;
use fs_material::state_point::ResolvedInterfaceStatePoint;
use fs_tribo::{TriboError, WearAdvance, WearState};

use crate::interface_binding::BoundNormalInterfaceState;
use crate::normal_patch::{
    FiniteGapContactError, FiniteGapGrid, FiniteGapPlasticConfig, FiniteGapPlasticContinuation,
    FiniteGapPlasticContinuationReceipt, FiniteGapPlasticDrive, FiniteGapPlasticDriveError,
    FiniteGapPlasticDriveReceipt, FiniteGapPlasticError, FiniteGapPlasticLinearSolver,
    FiniteGapPlasticLoadControls, FiniteGapPlasticLoadError, FiniteGapPlasticLoadReceipt,
    FiniteGapPlasticModel, FiniteGapPlasticReceipt, FiniteGapPlasticState,
};

use super::{
    PlasticContactControls, PlasticContactEnvelope, SourcePlasticContactError,
    admit_native_candidate, within,
};

/// Original fixed hardness or a complete accepted native-history field.
#[derive(Clone, Debug, PartialEq)]
pub enum ContactHardnessField {
    /// One original indentation test at fixed external conditions.
    Fixed(IndentationHardnessStatePoint),
    /// Ordered accepted material points from one original constitutive source.
    Accepted(Arc<AcceptedHardnessField>),
}

impl ContactHardnessField {
    /// Number of explicitly addressable original material points (one if fixed).
    #[must_use]
    pub fn len(&self) -> usize {
        match self {
            Self::Fixed(_) => 1,
            Self::Accepted(field) => field.len(),
        }
    }

    /// Whether this source contains no material points.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Original immutable source descriptor, independent of accepted alpha.
    #[must_use]
    pub fn source_identity(&self) -> ContentHash {
        match self {
            Self::Fixed(source) => source.identity(),
            Self::Accepted(field) => field.source_identity(),
        }
    }

    /// Original indentation source, retaining its test, claim and conditions.
    #[must_use]
    pub fn hardness_source(&self) -> &IndentationHardnessStatePoint {
        match self {
            Self::Fixed(source) => source,
            Self::Accepted(field) => field.hardness_source(),
        }
    }

    /// Hardness at an explicit original material-point index [Pa].
    #[must_use]
    pub fn hardness_pa(&self, point: usize) -> Option<f64> {
        match self {
            Self::Fixed(source) => (point == 0).then(|| source.hardness_pa()),
            Self::Accepted(field) => field.hardness_pa(point),
        }
    }
}

/// One native normal boundary. Wear and field refresh hold this boundary fixed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum HistoryHardnessControl {
    /// Prescribed approach [m]; a wear/refresh call must retain current approach.
    Approach(f64),
    /// Maintained compressive resultant with its original native search bounds.
    Load {
        force_n: f64,
        controls: FiniteGapPlasticLoadControls,
    },
    /// Fixed finite loading spring and drive position.
    Drive {
        drive: FiniteGapPlasticDrive,
        controls: FiniteGapPlasticLoadControls,
    },
}

/// Original native boundary receipt attached to an accepted physical step.
#[derive(Clone, Debug, PartialEq)]
pub enum HistoryHardnessControlReceipt {
    /// Prescribed approach supplied no additional load-search receipt.
    Approach,
    /// Actual native maintained-force equilibrium.
    Load(FiniteGapPlasticLoadReceipt),
    /// Actual native finite-spring equilibrium and storage.
    Drive(FiniteGapPlasticDriveReceipt),
}

/// Source, mapping, applicability, native solve or inventory refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum HistoryHardnessContactError {
    /// Original source, field ordering, arithmetic or boundary was incompatible.
    Invalid(&'static str),
    /// Original native normal-contact or physical-domain refusal.
    Contact(SourcePlasticContactError),
    /// Native prescribed-force refusal.
    Load(FiniteGapPlasticLoadError),
    /// Native finite-loading-spring refusal.
    Drive(FiniteGapPlasticDriveError),
    /// Original local-hardness Archard inventory refusal.
    Tribology(TriboError),
}

impl From<SourcePlasticContactError> for HistoryHardnessContactError {
    fn from(error: SourcePlasticContactError) -> Self {
        Self::Contact(error)
    }
}
impl From<FiniteGapPlasticError> for HistoryHardnessContactError {
    fn from(error: FiniteGapPlasticError) -> Self {
        Self::Contact(error.into())
    }
}
impl From<FiniteGapPlasticLoadError> for HistoryHardnessContactError {
    fn from(error: FiniteGapPlasticLoadError) -> Self {
        Self::Load(error)
    }
}
impl From<FiniteGapPlasticDriveError> for HistoryHardnessContactError {
    fn from(error: FiniteGapPlasticDriveError) -> Self {
        Self::Drive(error)
    }
}
impl From<TriboError> for HistoryHardnessContactError {
    fn from(error: TriboError) -> Self {
        Self::Tribology(error)
    }
}
impl core::fmt::Display for HistoryHardnessContactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "accepted-history hardness contact/wear: {self:?}")
    }
}
impl std::error::Error for HistoryHardnessContactError {}

#[derive(Clone, Debug, PartialEq)]
struct Sources {
    normal: BoundNormalInterfaceState,
    coefficients: ResolvedInterfaceStatePoint,
    coefficient_keys: [PropertyKey; 2],
    k: [f64; 2],
    cell_to_material_point: Vec<[usize; 2]>,
    envelope: PlasticContactEnvelope,
    maximum_recession_m: [f64; 2],
    identity: ContentHash,
}

/// Complete immutable contact, mapped hardness and separate body wear state.
/// Source data and compliance kernels are shared by accepted successors.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceHistoryHardnessWear {
    sources: Arc<Sources>,
    hardness: [ContactHardnessField; 2],
    model: FiniteGapPlasticModel,
    contact: FiniteGapPlasticState,
    volumes: [WearState; 2],
    recession_m: [Vec<f64>; 2],
}

/// Ordinary normal loading without changing material history or body inventories.
#[derive(Clone, Debug, PartialEq)]
pub struct HistoryHardnessNormalAdvance {
    /// Complete immutable successor.
    pub candidate: SourceHistoryHardnessWear,
    /// Original native contact, elastic storage and indentation work.
    pub receipt: FiniteGapPlasticReceipt,
    /// Actual native boundary result.
    pub control: HistoryHardnessControlReceipt,
}

/// Re-equilibration after one externally supplied accepted bulk-history update.
#[derive(Clone, Debug, PartialEq)]
pub struct HistoryHardnessRefresh {
    /// Complete successor with the same map, body inventories and original sources.
    pub candidate: SourceHistoryHardnessWear,
    /// Original native signed parameter work and retained indentation history.
    pub receipt: FiniteGapPlasticContinuationReceipt,
    /// Actual native boundary result.
    pub control: HistoryHardnessControlReceipt,
}

/// Local two-body material removal and the native normal re-equilibration.
#[derive(Clone, Debug, PartialEq)]
pub struct HistoryHardnessWearAdvance {
    /// Complete successor; rejection leaves the predecessor untouched.
    pub candidate: SourceHistoryHardnessWear,
    /// Separate original Archard inventory receipts in A/B order.
    pub wear: [WearAdvance; 2],
    /// Actual retained removal at each original body/grid cell [m].
    pub recession_increment_m: [Vec<f64>; 2],
    /// Sum of beginning pressure times actual unsigned slip times cell area [J].
    /// This is a sliding measure and is not automatically heat.
    pub pressure_weighted_sliding_j: f64,
    /// Original native normal storage, work and irreversible indentation receipt.
    pub receipt: FiniteGapPlasticContinuationReceipt,
    /// Actual native boundary result.
    pub control: HistoryHardnessControlReceipt,
}

impl SourceHistoryHardnessWear {
    /// Bind original fixed elasticity/interface coefficients and an explicit
    /// original-point map. Both bodies retain their independent hardness source.
    /// Mapped columns must remain representative through each recession bound.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        normal: BoundNormalInterfaceState,
        coefficients: ResolvedInterfaceStatePoint,
        coefficient_keys: [PropertyKey; 2],
        hardness: [ContactHardnessField; 2],
        cell_to_material_point: Vec<[usize; 2]>,
        grid: FiniteGapGrid,
        undeformed_gap_m: Vec<f64>,
        controls: PlasticContactControls,
        envelope: PlasticContactEnvelope,
        maximum_recession_m: [f64; 2],
        linear_solver: FiniteGapPlasticLinearSolver,
        cx: &Cx<'_>,
    ) -> Result<Self, HistoryHardnessContactError> {
        checkpoint(cx)?;
        let (k, selected) = binding::admit_sources(
            &normal,
            &coefficients,
            &coefficient_keys,
            &hardness,
            envelope,
            maximum_recession_m,
        )?;
        if cell_to_material_point.len() != undeformed_gap_m.len() {
            return Err(invalid(
                "one explicit original material-point pair per grid cell required",
            ));
        }
        let caps = binding::mapped_caps(&hardness, &cell_to_material_point, cx)?;
        let count = caps.len();
        let model = FiniteGapPlasticModel::new_with_linear_solver(
            FiniteGapPlasticConfig {
                grid,
                undeformed_gap_m,
                reduced_modulus_pa: normal.elastic().reduced_modulus_pa(),
                hardness_pa: caps,
                maximum_active_set_iterations: controls.maximum_active_set_iterations,
                complementarity_tolerance_m: controls.complementarity_tolerance_m,
                boundary_clearance_cells: controls.boundary_clearance_cells,
            },
            linear_solver,
            cx,
        )?;
        let sources = Sources {
            identity: binding::identity(
                &normal,
                &coefficients,
                &selected,
                &hardness,
                &cell_to_material_point,
                &model,
                envelope,
                maximum_recession_m,
                cx,
            )?,
            normal,
            coefficients,
            coefficient_keys,
            k,
            cell_to_material_point,
            envelope,
            maximum_recession_m,
        };
        checkpoint(cx)?;
        Ok(Self {
            sources: Arc::new(sources),
            hardness,
            contact: model.initial_state(),
            model,
            volumes: [WearState::default(); 2],
            recession_m: [vec![0.0; count], vec![0.0; count]],
        })
    }

    /// Complete original source/configuration identity, preserved by continuations.
    #[must_use]
    pub fn source_identity(&self) -> ContentHash {
        self.sources.identity
    }
    /// Original fixed elastic and dry-interface sources.
    #[must_use]
    pub fn normal_source(&self) -> &BoundNormalInterfaceState {
        &self.sources.normal
    }
    /// Original interface bundle supplying both independent wear coefficients.
    #[must_use]
    pub fn coefficient_source(&self) -> &ResolvedInterfaceStatePoint {
        &self.sources.coefficients
    }
    /// Exact original A/B coefficient keys.
    #[must_use]
    pub fn coefficient_keys(&self) -> &[PropertyKey; 2] {
        &self.sources.coefficient_keys
    }
    /// Current immutable hardness fields in original body order.
    #[must_use]
    pub fn hardness_fields(&self) -> &[ContactHardnessField; 2] {
        &self.hardness
    }
    /// Explicit retained mapping; each row contains original A/B point indices.
    #[must_use]
    pub fn cell_to_material_point(&self) -> &[[usize; 2]] {
        &self.sources.cell_to_material_point
    }
    /// Native model with its current local pressure caps and original compliance.
    #[must_use]
    pub const fn model(&self) -> &FiniteGapPlasticModel {
        &self.model
    }
    /// Complete native accepted indentation and pressure history.
    #[must_use]
    pub const fn contact_state(&self) -> &FiniteGapPlasticState {
        &self.contact
    }
    /// Actual cumulative per-cell material removal in original A/B order [m].
    #[must_use]
    pub fn recession_m(&self) -> [&[f64]; 2] {
        [&self.recession_m[0], &self.recession_m[1]]
    }
    /// Cumulative original body-volume inventories [m³].
    #[must_use]
    pub fn volumes_m3(&self) -> [f64; 2] {
        self.volumes.map(|state| state.volume_m3())
    }

    /// Native normal loading at the fixed accepted hardness snapshot.
    pub fn advance_normal(
        &self,
        control: HistoryHardnessControl,
        cx: &Cx<'_>,
    ) -> Result<HistoryHardnessNormalAdvance, HistoryHardnessContactError> {
        self.admit_control(control, false, cx)?;
        let (step, control) = match control {
            HistoryHardnessControl::Approach(approach) => (
                self.model.advance(&self.contact, approach, cx)?,
                HistoryHardnessControlReceipt::Approach,
            ),
            HistoryHardnessControl::Load { force_n, controls } => {
                let step = self
                    .model
                    .advance_at_load(&self.contact, force_n, controls, cx)?;
                (step.contact, HistoryHardnessControlReceipt::Load(step.load))
            }
            HistoryHardnessControl::Drive { drive, controls } => {
                let step = self
                    .model
                    .advance_with_drive(&self.contact, drive, controls, cx)?;
                (
                    step.contact,
                    HistoryHardnessControlReceipt::Drive(step.drive),
                )
            }
        };
        self.admit(&self.model, &step.candidate, &self.recession_m, cx)?;
        let candidate = Self {
            sources: Arc::clone(&self.sources),
            hardness: self.hardness.clone(),
            model: self.model.clone(),
            contact: step.candidate,
            volumes: self.volumes,
            recession_m: self.recession_m.clone(),
        };
        Ok(HistoryHardnessNormalAdvance {
            candidate,
            receipt: step.receipt,
            control,
        })
    }

    /// Adopt an independently accepted successor field without changing its
    /// original source, native point order/count, cell map or fixed environment.
    /// Every native alpha must be nondecreasing; H may harden or soften according
    /// to its original test. Native cap refresh owns all ensuing indentation and
    /// signed parameter work; no bulk work, residual stress or heat is inferred.
    pub fn refresh_hardness(
        &self,
        hardness: [ContactHardnessField; 2],
        control: HistoryHardnessControl,
        cx: &Cx<'_>,
    ) -> Result<HistoryHardnessRefresh, HistoryHardnessContactError> {
        self.admit_control(control, true, cx)?;
        binding::admit_successor(&self.hardness, &hardness, cx)?;
        let caps = binding::mapped_caps(&hardness, self.cell_to_material_point(), cx)?;
        let model = self.model.with_material_properties(
            self.model.config().reduced_modulus_pa,
            caps,
            cx,
        )?;
        let (step, control) = match control {
            HistoryHardnessControl::Approach(_) => (
                model.refresh_from(&self.model, &self.contact, cx)?,
                HistoryHardnessControlReceipt::Approach,
            ),
            HistoryHardnessControl::Load { force_n, controls } => {
                let step = model.refresh_from_at_load(
                    &self.model,
                    &self.contact,
                    force_n,
                    controls,
                    cx,
                )?;
                (
                    step.continuation,
                    HistoryHardnessControlReceipt::Load(step.load),
                )
            }
            HistoryHardnessControl::Drive { drive, controls } => {
                let step = model.refresh_from_with_drive(
                    &self.model,
                    &self.contact,
                    drive,
                    controls,
                    cx,
                )?;
                (
                    step.continuation,
                    HistoryHardnessControlReceipt::Drive(step.drive),
                )
            }
        };
        self.admit(&model, &step.candidate, &self.recession_m, cx)?;
        let candidate = Self {
            sources: Arc::clone(&self.sources),
            hardness,
            model,
            contact: step.candidate,
            volumes: self.volumes,
            recession_m: self.recession_m.clone(),
        };
        Ok(HistoryHardnessRefresh {
            candidate,
            receipt: step.receipt,
            control,
        })
    }

    fn admit_control(
        &self,
        control: HistoryHardnessControl,
        held: bool,
        cx: &Cx<'_>,
    ) -> Result<(), HistoryHardnessContactError> {
        checkpoint(cx)?;
        match control {
            HistoryHardnessControl::Approach(approach) => {
                within(
                    "approach [m]",
                    approach,
                    self.sources.envelope.maximum_approach_m,
                )?;
                if held && approach != self.contact.approach_m() {
                    return Err(invalid(
                        "wear and hardness refresh must retain the accepted approach",
                    ));
                }
            }
            HistoryHardnessControl::Load { controls, .. }
            | HistoryHardnessControl::Drive { controls, .. } => {
                within(
                    "load-control maximum approach [m]",
                    controls.maximum_approach_m,
                    self.sources.envelope.maximum_approach_m,
                )?;
            }
        }
        Ok(())
    }

    fn admit(
        &self,
        model: &FiniteGapPlasticModel,
        state: &FiniteGapPlasticState,
        recession: &[Vec<f64>; 2],
        cx: &Cx<'_>,
    ) -> Result<(), HistoryHardnessContactError> {
        let radius = admit_native_candidate(model, self.sources.envelope, state, cx)?;
        for side in 0..2 {
            let maximum = recession[side].iter().copied().fold(0.0_f64, f64::max);
            let remaining = self.sources.envelope.half_space_depth_m[side] - maximum;
            if remaining <= 0.0 {
                return Err(invalid("wear exhausted original homogeneous elastic depth"));
            }
            within(
                "patch radius / remaining homogeneous elastic depth",
                radius / remaining,
                self.sources.envelope.maximum_patch_to_depth[side],
            )?;
        }
        checkpoint(cx)
    }

    fn recede(
        &self,
        increment: &[f64],
        control: HistoryHardnessControl,
        cx: &Cx<'_>,
    ) -> Result<
        (FiniteGapPlasticContinuation, HistoryHardnessControlReceipt),
        HistoryHardnessContactError,
    > {
        Ok(match control {
            HistoryHardnessControl::Approach(_) => (
                self.model.recede(&self.contact, increment, cx)?,
                HistoryHardnessControlReceipt::Approach,
            ),
            HistoryHardnessControl::Load { force_n, controls } => {
                let step =
                    self.model
                        .recede_at_load(&self.contact, increment, force_n, controls, cx)?;
                (
                    step.continuation,
                    HistoryHardnessControlReceipt::Load(step.load),
                )
            }
            HistoryHardnessControl::Drive { drive, controls } => {
                let step =
                    self.model
                        .recede_with_drive(&self.contact, increment, drive, controls, cx)?;
                (
                    step.continuation,
                    HistoryHardnessControlReceipt::Drive(step.drive),
                )
            }
        })
    }
}

fn invalid(reason: &'static str) -> HistoryHardnessContactError {
    HistoryHardnessContactError::Invalid(reason)
}
fn checkpoint(cx: &Cx<'_>) -> Result<(), HistoryHardnessContactError> {
    cx.checkpoint()
        .map_err(|_| FiniteGapPlasticError::Contact(FiniteGapContactError::Cancelled).into())
}
