//! Native bulk-profile refresh with exactly retained actual sliding history.
//!
//! The depth owner admits the original source/history continuation and solves
//! the requested normal boundary. Every successor loaded pressure must then
//! admit the original fixed-temperature friction source before publication.
//! No shear interval, removal, frictional heat or synthetic travel is introduced.

use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_material::hardness::depth::ResolvedDepthHardness;
use fs_material::state_point::ResolvedInterfaceStatePoint;

use super::{
    DepthHardnessWearError, HardnessSlidingBoundaryReceipt, HardnessSlidingContinuation,
    HardnessSlidingControl, HardnessSlidingError, SourceDepthSlidingWear,
    SourcePlasticSlidingError, cells, checkpoint,
};

/// Complete profile/contact successor and unchanged genuine sliding memory.
#[derive(Clone, Debug, PartialEq)]
pub struct DepthSlidingHistoryRefresh {
    pub original_source_identity: ContentHash,
    pub successor_source_identity: ContentHash,
    pub candidate: SourceDepthSlidingWear,
    /// Native normal storage, indentation and signed parameter work, with boundary.
    pub normal: HardnessSlidingContinuation,
    /// Original friction claims at every actual successor loaded pressure.
    pub friction_sources: Vec<ResolvedInterfaceStatePoint>,
    /// Shared query index per cell; detached cells make no zero-pressure claim.
    pub cell_source_indices: Vec<Option<usize>>,
}

impl SourceDepthSlidingWear {
    /// Adopt independently accepted depth/history profiles without a shear step.
    ///
    /// The native depth owner retains both removal inventories and admits exact
    /// material/test/family pins, nondecreasing alpha and continuous depth support.
    /// Held approach must equal the accepted approach; force and finite-spring
    /// control delegate directly to their existing normal continuation owners.
    /// Complete signed shear memory and its last genuine temperature are unchanged.
    /// Normal parameter work remains distinct from frictional heat.
    pub fn refresh_hardness_profiles(
        &self,
        hardness: [ResolvedDepthHardness; 2],
        control: HardnessSlidingControl,
        cx: &Cx<'_>,
    ) -> Result<DepthSlidingHistoryRefresh, HardnessSlidingError> {
        checkpoint(cx).map_err(SourcePlasticSlidingError::from)?;
        let (history, boundary) = match control {
            HardnessSlidingControl::Approach(value) => {
                if value != self.wear.contact_state().approach_m() {
                    return Err(DepthHardnessWearError::Invalid(
                        "sliding profile refresh must retain the accepted approach",
                    )
                    .into());
                }
                (
                    self.wear.refresh_hardness_profiles(hardness, cx)?,
                    HardnessSlidingBoundaryReceipt::Approach,
                )
            }
            HardnessSlidingControl::Load { force_n, controls } => {
                let next = self
                    .wear
                    .refresh_hardness_profiles_at_load(hardness, force_n, controls, cx)?;
                (
                    next.history,
                    HardnessSlidingBoundaryReceipt::Load(next.load),
                )
            }
            HardnessSlidingControl::Drive { drive, controls } => {
                let next = self
                    .wear
                    .refresh_hardness_profiles_with_drive(hardness, drive, controls, cx)?;
                (
                    next.history,
                    HardnessSlidingBoundaryReceipt::Drive(next.drive),
                )
            }
        };
        let wear = history.candidate;
        let query = cells::query_at_temperature(&self.sliding.friction, self.sliding.temperature_k)
            .map_err(SourcePlasticSlidingError::Source)?;
        cells::admit_query(wear.normal_source(), &self.sliding.friction, &query)?;
        let grid = wear.model().config().grid;
        let samples = cells::sample_pressures(
            &self.sliding.friction,
            self.sliding.temperature_k,
            wear.contact_state().pressure_pa(),
            grid.cell_width_m * grid.cell_depth_m,
            cx,
        )?;
        let candidate = Self {
            wear,
            sliding: self.sliding.clone(),
        };
        checkpoint(cx).map_err(SourcePlasticSlidingError::from)?;
        Ok(DepthSlidingHistoryRefresh {
            original_source_identity: self.source_identity(),
            successor_source_identity: candidate.source_identity(),
            candidate,
            normal: HardnessSlidingContinuation {
                normal: history.receipt,
                boundary,
            },
            friction_sources: samples.sources,
            cell_source_indices: samples.cell_source_indices,
        })
    }
}
