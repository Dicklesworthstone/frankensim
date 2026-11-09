//! Actual native Coulomb slip driving original depth or bulk-history hardness.
//!
//! These concrete adapters share the existing cell closure with ordinary wear.
//! Only the unsigned length of actual irreversible slip removes material;
//! released body motion and net signed travel are not wear distances. Original
//! spatial hardness and normal continuation remain with their existing owners.
//! This is a beginning-pressure split with fixed axes and diagonal external
//! body compliance, not cross-cell tangential elasticity or a heat model.

pub mod temperature;
/// Full cross-cell mechanical recoil driving actual slip and own-hardness wear.
pub mod coupled;

use std::sync::Arc;

use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_material::state_point::ResolvedInterfaceStatePoint;
use fs_tribo::{
    DryTangentialState, FrictionLaw, TangentialLaw, TangentialSolveConfig, WearAdvance,
};

use crate::interface_binding::{BoundNormalInterfaceState, thermal_friction::SourceThermalCoulomb};
use crate::normal_patch::plastic_source::depth_wear::{
    DepthHardnessWearError, SourceDepthHardnessWear,
};
use crate::normal_patch::plastic_source::history_field::{
    ContactHardnessField, HistoryHardnessContactError, SourceHistoryHardnessWear,
};
pub use crate::normal_patch::plastic_source::history_field::{
    HistoryHardnessControl as HardnessSlidingControl,
    HistoryHardnessControlReceipt as HardnessSlidingBoundaryReceipt,
};
use crate::normal_patch::{
    FiniteGapPlasticContinuationReceipt, FiniteGapPlasticModel, FiniteGapPlasticReceipt,
    FiniteGapPlasticState,
};

use super::{
    PlasticCellSliding, PlasticSlidingInput, SourcePlasticSlidingError, cells, checkpoint,
};

/// Refusal leaves both the original wear owner and shear history unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum HardnessSlidingError {
    /// Original Coulomb source, local body response or native shear refusal.
    Sliding(SourcePlasticSlidingError),
    /// Original depth-profile, wear inventory or native normal refusal.
    Depth(DepthHardnessWearError),
    /// Original bulk-history field, wear inventory or native normal refusal.
    History(HistoryHardnessContactError),
}
impl From<SourcePlasticSlidingError> for HardnessSlidingError {
    fn from(error: SourcePlasticSlidingError) -> Self {
        Self::Sliding(error)
    }
}
impl From<DepthHardnessWearError> for HardnessSlidingError {
    fn from(error: DepthHardnessWearError) -> Self {
        Self::Depth(error)
    }
}
impl From<HistoryHardnessContactError> for HardnessSlidingError {
    fn from(error: HistoryHardnessContactError) -> Self {
        Self::History(error)
    }
}
impl core::fmt::Display for HardnessSlidingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "original spatial-hardness sliding/wear: {self:?}")
    }
}
impl std::error::Error for HardnessSlidingError {}

/// Original shear responses and source samples for one actual sliding interval.
#[derive(Clone, Debug, PartialEq)]
pub struct HardnessSlidingReceipt {
    /// Native responses at beginning pressure, including unloaded cells.
    pub cells: Vec<PlasticCellSliding>,
    /// Original pinned friction samples, shared across equal loaded queries.
    pub friction_sources: Vec<ResolvedInterfaceStatePoint>,
    /// Actual unsigned local distances consumed by the wear owner [m].
    pub unsigned_slip_m: Vec<f64>,
    /// Original physical duration used in each native cell solve [s].
    pub duration_s: f64,
    /// Original finite per-cell solver and effort budgets.
    pub friction_controls: TangentialSolveConfig,
    /// Compensated native friction dissipation, without heat allocation [J].
    pub frictional_work_j: f64,
    /// Beginning-pressure integral of actual unsigned slip over cell area [J].
    pub pressure_weighted_sliding_j: f64,
}

/// One original normal continuation paired with its actual boundary receipt.
#[derive(Clone, Debug, PartialEq)]
pub struct HardnessSlidingContinuation {
    /// Native normal storage, irreversible indentation and parameter work.
    pub normal: FiniteGapPlasticContinuationReceipt,
    /// Held approach, maintained force or finite loading-spring result.
    pub boundary: HardnessSlidingBoundaryReceipt,
}

/// Complete immutable successor plus the distinct shear, wear and normal accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct HardnessSlidingAdvance<C> {
    /// Original normal, hardness, wear and friction registration.
    pub source_identity: ContentHash,
    /// Accept this entire original-source candidate, including its shear history.
    pub candidate: C,
    /// Native actual slip, traction, source samples and friction work.
    pub sliding: HardnessSlidingReceipt,
    /// Original A/B Archard inventory receipts, each using its own hardness.
    pub wear: [WearAdvance; 2],
    /// Actual retained A/B material recession consumed by native contact [m].
    pub recession_increment_m: [Vec<f64>; 2],
    /// One recession for frozen fields; recession then cap exposure for depth H.
    pub normal: Vec<HardnessSlidingContinuation>,
}

/// Ordinary normal loading while actual shear travel and body wear remain intact.
#[derive(Clone, Debug, PartialEq)]
pub struct HardnessSlidingNormalAdvance<C> {
    /// Complete coupled successor.
    pub candidate: C,
    /// Original native normal response; no unresolved actuator path is inferred.
    pub receipt: FiniteGapPlasticReceipt,
    /// Actual native normal boundary result.
    pub boundary: HardnessSlidingBoundaryReceipt,
}

/// External bulk-history refresh with retained actual shear and wear histories.
#[derive(Clone, Debug, PartialEq)]
pub struct HardnessSlidingRefresh<C> {
    /// Complete coupled successor.
    pub candidate: C,
    /// Native material continuation and its separate signed parameter work.
    pub normal: HardnessSlidingContinuation,
}

#[derive(Clone, Debug, PartialEq)]
struct SlidingState {
    friction: Arc<SourceThermalCoulomb>,
    initial_law: FrictionLaw,
    temperature_k: f64,
    tangential: Vec<DryTangentialState>,
}
impl SlidingState {
    fn new(
        normal: &BoundNormalInterfaceState,
        count: usize,
        friction: SourceThermalCoulomb,
    ) -> Result<Self, HardnessSlidingError> {
        let query = friction.initial_source().query();
        cells::admit_query(normal, &friction, &query)?;
        let temperature_k = query.axes()[&friction.axes().temperature];
        let (_, initial_law) = friction
            .coefficients_at_query(&query)
            .map_err(SourcePlasticSlidingError::Source)?;
        let tangential = vec![initial_law.initial_state(); count];
        Ok(Self {
            friction: Arc::new(friction),
            initial_law,
            temperature_k,
            tangential,
        })
    }
    fn advance(
        &self,
        model: &FiniteGapPlasticModel,
        contact: &FiniteGapPlasticState,
        motion: &[PlasticSlidingInput],
        duration: f64,
        controls: TangentialSolveConfig,
        cx: &Cx<'_>,
    ) -> Result<cells::CellSlidingAdvance, HardnessSlidingError> {
        let grid = model.config().grid;
        Ok(cells::advance(
            &self.friction,
            &self.initial_law,
            self.temperature_k,
            &self.tangential,
            contact.pressure_pa(),
            grid.cell_width_m * grid.cell_depth_m,
            motion,
            duration,
            controls,
            cx,
        )?)
    }
    fn finish(
        &self,
        step: cells::CellSlidingAdvance,
        area: f64,
        duration_s: f64,
        friction_controls: TangentialSolveConfig,
        cx: &Cx<'_>,
    ) -> Result<(Self, HardnessSlidingReceipt), HardnessSlidingError> {
        let mut measure = 0.0;
        let mut correction = 0.0;
        for (index, (cell, distance)) in step.cells.iter().zip(&step.unsigned_slip_m).enumerate() {
            if index % 64 == 0 {
                checkpoint(cx).map_err(SourcePlasticSlidingError::from)?;
            }
            let contribution = (cell.pressure_pa * area) * distance;
            if !contribution.is_finite() {
                return Err(SourcePlasticSlidingError::Invalid(
                    "finite actual pressure-sliding measure",
                )
                .into());
            }
            super::add_compensated(&mut measure, &mut correction, contribution);
        }
        if !measure.is_finite() {
            return Err(SourcePlasticSlidingError::Invalid(
                "finite total actual pressure-sliding measure",
            )
            .into());
        }
        checkpoint(cx).map_err(SourcePlasticSlidingError::from)?;
        Ok((
            Self {
                friction: Arc::clone(&self.friction),
                initial_law: self.initial_law.clone(),
                temperature_k: self.temperature_k,
                tangential: step.tangential,
            },
            HardnessSlidingReceipt {
                cells: step.cells,
                friction_sources: step.friction_sources,
                unsigned_slip_m: step.unsigned_slip_m,
                duration_s,
                friction_controls,
                frictional_work_j: step.frictional_work_j,
                pressure_weighted_sliding_j: measure,
            },
        ))
    }
    fn identity(&self, wear: ContentHash) -> ContentHash {
        let mut hash = DomainHasher::new("org.frankensim.fs-contact.spatial-hardness-sliding.v1");
        hash.update(wear.as_bytes());
        hash.update(self.friction.identity().as_bytes());
        hash.update(&self.temperature_k.to_bits().to_le_bytes());
        hash.finalize()
    }
}

/// Complete original depth-hardness wear owner and its native actual shear history.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceDepthSlidingWear {
    wear: SourceDepthHardnessWear,
    sliding: SlidingState,
}
impl SourceDepthSlidingWear {
    /// Begin zero signed shear travel on the supplied accepted wear checkpoint.
    /// Existing dents, body recession and inventories remain original; this does
    /// not import a shear history from another registration.
    pub fn new(
        wear: SourceDepthHardnessWear,
        friction: SourceThermalCoulomb,
    ) -> Result<Self, HardnessSlidingError> {
        let sliding = SlidingState::new(
            wear.normal_source(),
            wear.contact_state().pressure_pa().len(),
            friction,
        )?;
        Ok(Self { wear, sliding })
    }
    /// Original complete depth-hardness/contact/wear checkpoint.
    pub const fn wear_state(&self) -> &SourceDepthHardnessWear {
        &self.wear
    }
    /// Native signed actual shear histories in unchanged cell/tangent order.
    pub fn tangential_states(&self) -> &[DryTangentialState] {
        &self.sliding.tangential
    }
    /// Original pinned friction source.
    pub fn friction_source(&self) -> &SourceThermalCoulomb {
        &self.sliding.friction
    }
    /// Original coupled source/configuration identity.
    pub fn source_identity(&self) -> ContentHash {
        self.sliding.identity(self.wear.source_identity())
    }

    /// Change only the normal boundary, retaining actual shear and body histories.
    pub fn advance_normal(
        &self,
        control: HardnessSlidingControl,
        cx: &Cx<'_>,
    ) -> Result<HardnessSlidingNormalAdvance<Self>, HardnessSlidingError> {
        let (normal, boundary) = match control {
            HardnessSlidingControl::Approach(value) => (
                self.wear.advance_normal(value, cx)?,
                HardnessSlidingBoundaryReceipt::Approach,
            ),
            HardnessSlidingControl::Load { force_n, controls } => {
                let next = self.wear.advance_normal_at_load(force_n, controls, cx)?;
                (next.normal, HardnessSlidingBoundaryReceipt::Load(next.load))
            }
            HardnessSlidingControl::Drive { drive, controls } => {
                let next = self.wear.advance_normal_with_drive(drive, controls, cx)?;
                (
                    next.normal,
                    HardnessSlidingBoundaryReceipt::Drive(next.drive),
                )
            }
        };
        Ok(HardnessSlidingNormalAdvance {
            candidate: Self {
                wear: normal.candidate,
                sliding: self.sliding.clone(),
            },
            receipt: normal.receipt,
            boundary,
        })
    }

    /// Resolve native actual slip, integrate each original H(z) through removal,
    /// and apply the requested held normal boundary at both native stages.
    /// Refine motion when pressure changes; no global coupled-energy theorem or
    /// wear from released/sticking motion is inferred.
    pub fn advance(
        &self,
        motion: &[PlasticSlidingInput],
        duration_s: f64,
        controls: TangentialSolveConfig,
        control: HardnessSlidingControl,
        cx: &Cx<'_>,
    ) -> Result<HardnessSlidingAdvance<Self>, HardnessSlidingError> {
        if let HardnessSlidingControl::Approach(value) = control
            && value != self.wear.contact_state().approach_m()
        {
            return Err(DepthHardnessWearError::Invalid(
                "sliding wear must retain the accepted approach",
            )
            .into());
        }
        let shear = self.sliding.advance(
            self.wear.model(),
            self.wear.contact_state(),
            motion,
            duration_s,
            controls,
            cx,
        )?;
        let (wear, boundaries) = match control {
            HardnessSlidingControl::Approach(_) => (
                self.wear
                    .advance_wear_integrated(&shear.unsigned_slip_m, cx)?,
                vec![
                    HardnessSlidingBoundaryReceipt::Approach,
                    HardnessSlidingBoundaryReceipt::Approach,
                ],
            ),
            HardnessSlidingControl::Load { force_n, controls } => {
                let next = self.wear.advance_wear_integrated_at_load(
                    &shear.unsigned_slip_m,
                    force_n,
                    controls,
                    cx,
                )?;
                (
                    next.wear,
                    vec![
                        HardnessSlidingBoundaryReceipt::Load(next.recession_load),
                        HardnessSlidingBoundaryReceipt::Load(next.hardness_load),
                    ],
                )
            }
            HardnessSlidingControl::Drive { drive, controls } => {
                let next = self.wear.advance_wear_integrated_with_drive(
                    &shear.unsigned_slip_m,
                    drive,
                    controls,
                    cx,
                )?;
                (
                    next.wear,
                    vec![
                        HardnessSlidingBoundaryReceipt::Drive(next.recession_drive),
                        HardnessSlidingBoundaryReceipt::Drive(next.hardness_drive),
                    ],
                )
            }
        };
        let grid = self.wear.model().config().grid;
        let (sliding, receipt) = self.sliding.finish(
            shear,
            grid.cell_width_m * grid.cell_depth_m,
            duration_s,
            controls,
            cx,
        )?;
        Ok(HardnessSlidingAdvance {
            source_identity: self.source_identity(),
            candidate: Self {
                wear: wear.candidate,
                sliding,
            },
            sliding: receipt,
            wear: wear.wear,
            recession_increment_m: wear.recession_increment_m,
            normal: [wear.recession, wear.hardness]
                .into_iter()
                .zip(boundaries)
                .map(|(normal, boundary)| HardnessSlidingContinuation { normal, boundary })
                .collect(),
        })
    }
}

/// Complete frozen native bulk-history hardness owner and actual shear history.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceHistorySlidingWear {
    wear: SourceHistoryHardnessWear,
    sliding: SlidingState,
}
impl SourceHistorySlidingWear {
    /// Begin zero signed shear travel without changing the accepted material map,
    /// native indentation, body recession or volume inventories.
    pub fn new(
        wear: SourceHistoryHardnessWear,
        friction: SourceThermalCoulomb,
    ) -> Result<Self, HardnessSlidingError> {
        let sliding = SlidingState::new(
            wear.normal_source(),
            wear.contact_state().pressure_pa().len(),
            friction,
        )?;
        Ok(Self { wear, sliding })
    }
    /// Original complete mapped-hardness/contact/wear checkpoint.
    pub const fn wear_state(&self) -> &SourceHistoryHardnessWear {
        &self.wear
    }
    /// Native signed actual shear histories in unchanged cell/tangent order.
    pub fn tangential_states(&self) -> &[DryTangentialState] {
        &self.sliding.tangential
    }
    /// Original pinned friction source.
    pub fn friction_source(&self) -> &SourceThermalCoulomb {
        &self.sliding.friction
    }
    /// Original coupled source/configuration identity.
    pub fn source_identity(&self) -> ContentHash {
        self.sliding.identity(self.wear.source_identity())
    }

    /// Change only normal loading, preserving the original shear and wear histories.
    pub fn advance_normal(
        &self,
        control: HardnessSlidingControl,
        cx: &Cx<'_>,
    ) -> Result<HardnessSlidingNormalAdvance<Self>, HardnessSlidingError> {
        let next = self.wear.advance_normal(control, cx)?;
        Ok(HardnessSlidingNormalAdvance {
            candidate: Self {
                wear: next.candidate,
                sliding: self.sliding.clone(),
            },
            receipt: next.receipt,
            boundary: next.control,
        })
    }

    /// Adopt an original accepted field through the native history owner while
    /// retaining all actual signed travel; no bulk deformation or work is added.
    pub fn refresh_hardness(
        &self,
        fields: [ContactHardnessField; 2],
        control: HardnessSlidingControl,
        cx: &Cx<'_>,
    ) -> Result<HardnessSlidingRefresh<Self>, HardnessSlidingError> {
        let next = self.wear.refresh_hardness(fields, control, cx)?;
        Ok(HardnessSlidingRefresh {
            candidate: Self {
                wear: next.candidate,
                sliding: self.sliding.clone(),
            },
            normal: HardnessSlidingContinuation {
                normal: next.receipt,
                boundary: next.control,
            },
        })
    }

    /// Resolve native actual slip and wear each mapped body at its own frozen H.
    /// Recession uses the original held approach, force or finite spring; fields
    /// remain frozen until a separately accepted original-history refresh.
    pub fn advance(
        &self,
        motion: &[PlasticSlidingInput],
        duration_s: f64,
        controls: TangentialSolveConfig,
        control: HardnessSlidingControl,
        cx: &Cx<'_>,
    ) -> Result<HardnessSlidingAdvance<Self>, HardnessSlidingError> {
        let shear = self.sliding.advance(
            self.wear.model(),
            self.wear.contact_state(),
            motion,
            duration_s,
            controls,
            cx,
        )?;
        let wear = self
            .wear
            .advance_wear(&shear.unsigned_slip_m, control, cx)?;
        let grid = self.wear.model().config().grid;
        let (sliding, receipt) = self.sliding.finish(
            shear,
            grid.cell_width_m * grid.cell_depth_m,
            duration_s,
            controls,
            cx,
        )?;
        Ok(HardnessSlidingAdvance {
            source_identity: self.source_identity(),
            candidate: Self {
                wear: wear.candidate,
                sliding,
            },
            sliding: receipt,
            wear: wear.wear,
            recession_increment_m: wear.recession_increment_m,
            normal: vec![HardnessSlidingContinuation {
                normal: wear.receipt,
                boundary: wear.control,
            }],
        })
    }
}
