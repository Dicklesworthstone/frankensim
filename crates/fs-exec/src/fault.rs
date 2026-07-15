//! Deterministic G4 tile-fault plans.
//!
//! A plan identifies exactly one logical tile and one bounded touch within
//! that tile. Selection depends only on the declared seed and dimensions, so
//! a failure receipt can be replayed without preserving worker or arrival
//! order.

use crate::TileFailure;

/// Schema version stamped beside retained tile-fault evidence.
pub const TILE_FAULT_PLAN_VERSION: u32 = 1;

/// A deterministic request to fail one logical tile at one numbered touch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TileFaultPlan {
    seed: u64,
    tile: u64,
    touch: u32,
}

impl TileFaultPlan {
    /// Derive a plan from a seed and non-empty tile/touch domains.
    ///
    /// # Errors
    /// Returns [`FaultPlanError::ZeroTiles`] when `tiles == 0` and
    /// [`FaultPlanError::ZeroTouches`] when `touches_per_tile == 0`.
    pub fn seeded(seed: u64, tiles: u64, touches_per_tile: u32) -> Result<Self, FaultPlanError> {
        if tiles == 0 {
            return Err(FaultPlanError::ZeroTiles);
        }
        if touches_per_tile == 0 {
            return Err(FaultPlanError::ZeroTouches);
        }

        let tile_word = mix(seed ^ 0x4653_2d54_494c_4521);
        let touch_word = mix(seed ^ 0x4653_2d54_4f55_4348);
        Ok(Self {
            seed,
            tile: tile_word % tiles,
            touch: u32::try_from(touch_word % u64::from(touches_per_tile))
                .expect("bounded touch index fits u32")
                + 1,
        })
    }

    /// Seed carried into the structured failure receipt.
    #[must_use]
    pub const fn seed(self) -> u64 {
        self.seed
    }

    /// Logical tile selected by this plan.
    #[must_use]
    pub const fn tile(self) -> u64 {
        self.tile
    }

    /// One-based touch selected within the logical tile.
    #[must_use]
    pub const fn touch(self) -> u32 {
        self.touch
    }

    /// Return the typed failure exactly at the selected tile/touch pair.
    #[must_use]
    pub fn failure_at(self, tile: u64, touch: u32) -> Option<TileFailure> {
        (tile == self.tile && touch == self.touch).then_some(TileFailure::InjectedFault {
            plan_seed: self.seed,
            touch,
        })
    }
}

/// Invalid dimensions for a deterministic tile-fault plan.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FaultPlanError {
    /// A plan cannot select from an empty tile domain.
    ZeroTiles,
    /// A plan cannot select from an empty per-tile touch domain.
    ZeroTouches,
}

impl core::fmt::Display for FaultPlanError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::ZeroTiles => write!(f, "tile fault plan requires at least one tile"),
            Self::ZeroTouches => {
                write!(f, "tile fault plan requires at least one touch per tile")
            }
        }
    }
}

impl core::error::Error for FaultPlanError {}

fn mix(mut value: u64) -> u64 {
    value = value.wrapping_add(0x9e37_79b9_7f4a_7c15);
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}

#[cfg(test)]
mod tests {
    use super::{FaultPlanError, TileFaultPlan};

    #[test]
    fn plan_is_replayable_and_fail_closed_on_empty_domains() {
        let plan = TileFaultPlan::seeded(0xF404, 17, 3).expect("valid plan");
        assert_eq!(plan, TileFaultPlan::seeded(0xF404, 17, 3).unwrap());
        assert!(plan.tile() < 17);
        assert!((1..=3).contains(&plan.touch()));
        assert_eq!(
            TileFaultPlan::seeded(0xF404, 0, 3),
            Err(FaultPlanError::ZeroTiles)
        );
        assert_eq!(
            TileFaultPlan::seeded(0xF404, 17, 0),
            Err(FaultPlanError::ZeroTouches)
        );
    }
}
