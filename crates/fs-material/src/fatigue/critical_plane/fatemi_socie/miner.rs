//! Optional fixed-plane Palmgren-Miner sums over declared closed tensor cycles.
//!
//! This is a maximum-accumulated-damage FS/Miner approximation, distinct from
//! selecting the maximum-shear plane of one representative cycle. Every cycle
//! contributes to the SAME material planes; per-cycle maxima are never added.
//! The caller resolves and declares one complete cycle per finalization. No
//! multiaxial rainflow decomposition, missed repetitions, sequence interaction,
//! stiffness degradation, or experimental life validation is inferred.

use super::{FatemiSocieError, FatigueError, SearchPlane};

/// Accumulated empirical life consumption on one fixed material plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FatemiSociePlaneDamage {
    pub normal_in_frame: [f64; 3],
    pub miner_damage: f64,
    compensation: f64,
}

/// Completed, explicitly declared cycles, including the current window only
/// after successful finalization. Open or discarded windows add no damage.
/// Damage >= 1 is the empirical rule's exhaustion threshold, not a certificate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FatemiSocieMinerSummary {
    pub completed_cycles: u64,
    pub maximum_cycles: u64,
    pub evaluated_planes: usize,
    /// Greatest TOTAL damage; exact ties select the first declared plane.
    /// Before the first completed cycle all planes have zero damage, selecting 0.
    pub selected_plane: usize,
    pub normal_in_frame: [f64; 3],
    pub miner_damage: f64,
}

#[derive(Debug, Clone)]
pub(super) struct DamageState {
    completed_cycles: u64,
    planes: Vec<FatemiSociePlaneDamage>,
}
impl DamageState {
    fn try_clone(&self) -> Result<Self, FatemiSocieError> {
        let mut planes = Vec::new();
        planes
            .try_reserve_exact(self.planes.len())
            .map_err(|_| FatemiSocieError::Allocation)?;
        planes.extend_from_slice(&self.planes);
        Ok(Self {
            completed_cycles: self.completed_cycles,
            planes,
        })
    }

    pub(super) fn add_life(
        &mut self,
        plane: usize,
        life: Option<f64>,
    ) -> Result<(), FatemiSocieError> {
        let Some(life) = life else {
            return Ok(());
        };
        // The calibration already admits a finite life >= 0.5 full cycles.
        // Keep the same compensated positive accumulation as scalar Miner.
        let state = &mut self.planes[plane];
        let corrected = 1.0 / life - state.compensation;
        let sum = state.miner_damage + corrected;
        let compensation = (sum - state.miner_damage) - corrected;
        if !sum.is_finite() || !compensation.is_finite() {
            return Err(FatigueError::Arithmetic.into());
        }
        state.miner_damage = sum;
        state.compensation = compensation;
        Ok(())
    }
}

/// A carried prefix and a separately reversible current-window contribution.
/// This preserves the native solvers' existing append-and-discard rollback:
/// finalization publishes current, fresh promotes it, discard removes it.
#[derive(Debug, Clone)]
pub(super) struct MinerHistory {
    maximum_cycles: u64,
    prefix: DamageState,
    current: Option<DamageState>,
}
impl MinerHistory {
    pub(super) fn new(
        planes: &[SearchPlane],
        maximum_cycles: u64,
    ) -> Result<Self, FatemiSocieError> {
        if maximum_cycles == 0 {
            return Err(FatemiSocieError::Invalid(
                "declare a positive complete-cycle limit",
            ));
        }
        let mut damage = Vec::new();
        damage
            .try_reserve_exact(planes.len())
            .map_err(|_| FatemiSocieError::Allocation)?;
        for plane in planes {
            damage.push(FatemiSociePlaneDamage {
                normal_in_frame: plane.normal,
                miner_damage: 0.0,
                compensation: 0.0,
            });
        }
        Ok(Self {
            maximum_cycles,
            prefix: DamageState {
                completed_cycles: 0,
                planes: damage,
            },
            current: None,
        })
    }

    fn latest(&self) -> &DamageState {
        self.current.as_ref().unwrap_or(&self.prefix)
    }

    pub(super) fn plane_damage(&self) -> &[FatemiSociePlaneDamage] {
        &self.latest().planes
    }

    pub(super) fn summary(&self) -> FatemiSocieMinerSummary {
        let state = self.latest();
        let mut selected_plane = 0;
        for (index, plane) in state.planes.iter().enumerate().skip(1) {
            if plane.miner_damage > state.planes[selected_plane].miner_damage {
                selected_plane = index;
            }
        }
        let plane = state.planes[selected_plane];
        FatemiSocieMinerSummary {
            completed_cycles: state.completed_cycles,
            maximum_cycles: self.maximum_cycles,
            evaluated_planes: state.planes.len(),
            selected_plane,
            normal_in_frame: plane.normal_in_frame,
            miner_damage: plane.miner_damage,
        }
    }

    pub(super) fn begin_cycle(&self) -> Result<DamageState, FatemiSocieError> {
        if self.prefix.completed_cycles >= self.maximum_cycles {
            return Err(FatemiSocieError::Budget);
        }
        let mut candidate = self.prefix.try_clone()?;
        candidate.completed_cycles += 1;
        Ok(candidate)
    }

    pub(super) fn publish_cycle(&mut self, candidate: DamageState) {
        self.current = Some(candidate);
    }

    pub(super) fn fresh(&self) -> Result<Self, FatemiSocieError> {
        Ok(Self {
            maximum_cycles: self.maximum_cycles,
            prefix: self.latest().try_clone()?,
            current: None,
        })
    }

    pub(super) fn discard_current(&mut self) {
        self.current = None;
    }
}
