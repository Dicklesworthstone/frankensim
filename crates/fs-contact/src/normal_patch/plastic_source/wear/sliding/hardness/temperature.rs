//! Original depth, elasticity, wear and friction sources at prescribed temperatures.
//!
//! A native normal thermal continuation is provisional until every loaded cell
//! admits its original friction query at the actual successor pressure. No shear
//! step is taken, so signed travel, wear inventory and last genuine shear-state
//! temperature remain unchanged. Normal parameter work is not frictional heat.

use std::sync::Arc;

use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_matdb::MaterialCard;
use fs_material::state_point::ResolvedInterfaceStatePoint;

use crate::normal_patch::plastic_source::depth_wear::temperature::DepthHardnessTemperatureError;

use super::{
    HardnessSlidingBoundaryReceipt, HardnessSlidingContinuation, HardnessSlidingControl,
    SlidingState, SourceDepthSlidingWear, SourcePlasticSlidingError, cells, checkpoint,
};

/// An original thermal or friction refusal leaves the complete predecessor intact.
#[derive(Clone, Debug, PartialEq)]
pub enum DepthSlidingTemperatureError {
    /// Original body/interface source or native normal boundary refusal.
    Depth(DepthHardnessTemperatureError),
    /// Original friction source, actual loaded pressure or cancellation refusal.
    Sliding(SourcePlasticSlidingError),
}
impl From<DepthHardnessTemperatureError> for DepthSlidingTemperatureError {
    fn from(error: DepthHardnessTemperatureError) -> Self {
        Self::Depth(error)
    }
}
impl From<SourcePlasticSlidingError> for DepthSlidingTemperatureError {
    fn from(error: SourcePlasticSlidingError) -> Self {
        Self::Sliding(error)
    }
}
impl core::fmt::Display for DepthSlidingTemperatureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "original depth-hardness sliding temperature continuation: {self:?}"
        )
    }
}
impl std::error::Error for DepthSlidingTemperatureError {}

/// Complete thermal endpoint, exact friction applicability and native normal work.
#[derive(Clone, Debug, PartialEq)]
pub struct DepthSlidingThermalRefresh {
    /// Coupled original source/configuration before the prescribed change.
    pub original_source_identity: ContentHash,
    /// Complete updated owner with exactly retained genuine shear and wear history.
    pub candidate: SourceDepthSlidingWear,
    /// Native normal storage, indentation and signed parameter work, with boundary.
    pub normal: HardnessSlidingContinuation,
    /// Original friction claims at actual successor loaded T/p queries.
    pub friction_sources: Vec<ResolvedInterfaceStatePoint>,
    /// Per-cell shared query index; detached cells make no zero-pressure claim.
    pub cell_source_indices: Vec<Option<usize>>,
}

impl SourceDepthSlidingWear {
    /// Prescribed interface temperature held through each subsequent shear interval.
    /// Native shear states still record the temperature of their last genuine step.
    pub const fn interface_temperature_k(&self) -> f64 {
        self.sliding.temperature_k
    }

    /// Continue original E(T), nu(T), H(T,z), k(T) and mu(T,p) without a shear step.
    ///
    /// Body temperatures independently select the original hardness depth profiles
    /// and elastic sources; interface temperature selects original wear/friction
    /// claims. The requested held approach, maintained force or fixed loading
    /// spring is applied by the native thermal owner directly. All loaded friction
    /// queries must accept before this complete candidate can be adopted.
    ///
    /// No reference-pressure-only admission, synthetic slip, material removal or
    /// frictional work is introduced. The caller owns temperature evolution and
    /// heat allocation; this endpoint update preserves the original tangent axes.
    pub fn refresh_at_temperatures(
        &self,
        material_cards: [&MaterialCard; 2],
        bulk_temperatures_k: [f64; 2],
        interface_temperature_k: f64,
        control: HardnessSlidingControl,
        cx: &Cx<'_>,
    ) -> Result<DepthSlidingThermalRefresh, DepthSlidingTemperatureError> {
        checkpoint(cx).map_err(SourcePlasticSlidingError::from)?;
        let query = cells::query_at_temperature(&self.sliding.friction, interface_temperature_k)
            .map_err(SourcePlasticSlidingError::Source)?;
        let (thermal, boundary) = match control {
            HardnessSlidingControl::Approach(value) => {
                if value != self.wear.contact_state().approach_m() {
                    return Err(DepthHardnessTemperatureError::Invalid(
                        "sliding temperature continuation must retain the accepted approach",
                    )
                    .into());
                }
                (
                    self.wear.refresh_at_temperatures(
                        material_cards,
                        bulk_temperatures_k,
                        interface_temperature_k,
                        cx,
                    )?,
                    HardnessSlidingBoundaryReceipt::Approach,
                )
            }
            HardnessSlidingControl::Load { force_n, controls } => {
                let next = self.wear.refresh_at_temperatures_at_load(
                    material_cards,
                    bulk_temperatures_k,
                    interface_temperature_k,
                    force_n,
                    controls,
                    cx,
                )?;
                (
                    next.thermal,
                    HardnessSlidingBoundaryReceipt::Load(next.load),
                )
            }
            HardnessSlidingControl::Drive { drive, controls } => {
                let next = self.wear.refresh_at_temperatures_with_drive(
                    material_cards,
                    bulk_temperatures_k,
                    interface_temperature_k,
                    drive,
                    controls,
                    cx,
                )?;
                (
                    next.thermal,
                    HardnessSlidingBoundaryReceipt::Drive(next.drive),
                )
            }
        };
        let wear = thermal.candidate;
        cells::admit_query(wear.normal_source(), &self.sliding.friction, &query)?;
        let grid = wear.model().config().grid;
        let samples = cells::sample_pressures(
            &self.sliding.friction,
            interface_temperature_k,
            wear.contact_state().pressure_pa(),
            grid.cell_width_m * grid.cell_depth_m,
            cx,
        )?;
        let candidate = Self {
            wear,
            sliding: SlidingState {
                friction: Arc::clone(&self.sliding.friction),
                // Only detached cells use this original fallback; N=0 carries
                // no source claim. Actual loaded motion resolves its own query.
                initial_law: self.sliding.initial_law.clone(),
                temperature_k: interface_temperature_k,
                tangential: self.sliding.tangential.clone(),
            },
        };
        checkpoint(cx).map_err(SourcePlasticSlidingError::from)?;
        Ok(DepthSlidingThermalRefresh {
            original_source_identity: self.source_identity(),
            candidate,
            normal: HardnessSlidingContinuation {
                normal: thermal.receipt,
                boundary,
            },
            friction_sources: samples.sources,
            cell_source_indices: samples.cell_source_indices,
        })
    }
}
