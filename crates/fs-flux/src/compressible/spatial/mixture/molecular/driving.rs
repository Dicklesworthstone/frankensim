//! Prescribed tangential acceleration on the original fixed molecular walls.

use super::{Boundaries3d, BoundaryCondition, BoundarySide, DissipationError, MolecularWall3d};
use super::{physical_condition, poll, side_index};

/// Constant Cartesian acceleration of one existing no-slip wall [m/s²].
///
/// The current accepted wall velocity is the initial value. On successful
/// advancement it becomes `u_wall + acceleration * elapsed`; a refused call
/// preserves that velocity. Each call authors only its accepted interval, so
/// callers can compose piecewise-constant acceleration histories. Unlisted
/// walls retain their velocity. Geometry and thermal conditions remain fixed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallAcceleration3d {
    /// A reflecting physical side already configured as no-slip.
    pub side: BoundarySide,
    /// Finite tangential acceleration; the normal component must be zero.
    pub acceleration: [f64; 3],
}

/// Bounded immutable history shared by both halves and every shortening retry.
pub(in crate::compressible::spatial::mixture) struct WallDrive3d {
    initial: [MolecularWall3d; 6],
    rates: [[f64; 3]; 6],
}

impl WallDrive3d {
    pub(super) fn new(
        initial: [MolecularWall3d; 6],
        boundaries: Boundaries3d,
        accelerations: &[WallAcceleration3d],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<Self>, DissipationError> {
        if accelerations.len() > 6 {
            return Err(DissipationError::Invalid(
                "at most six distinct molecular wall accelerations",
            ));
        }
        let mut seen = [false; 6];
        let mut rates = [[0.0; 3]; 6];
        let mut active = false;
        for drive in accelerations {
            poll(cancelled)?;
            let side = drive.side;
            let index = side_index(side);
            if seen[index] {
                return Err(DissipationError::Invalid(
                    "duplicate molecular wall acceleration",
                ));
            }
            seen[index] = true;
            if physical_condition(boundaries, side) != Some(BoundaryCondition::Reflecting)
                || !matches!(initial[index], MolecularWall3d::NoSlip { .. })
            {
                return Err(DissipationError::Invalid(
                    "wall acceleration requires an existing reflecting no-slip wall",
                ));
            }
            if drive.acceleration.iter().any(|a| !a.is_finite())
                || drive.acceleration[side.axis()] != 0.0
            {
                return Err(DissipationError::Invalid(
                    "molecular wall acceleration must be finite and tangential",
                ));
            }
            rates[index] = drive.acceleration;
            active |= drive.acceleration.iter().any(|a| *a != 0.0);
        }
        // A valid all-zero history takes the original held-wall path bit for bit.
        Ok(active.then_some(Self { initial, rates }))
    }

    pub(in crate::compressible::spatial::mixture) fn at(
        &self,
        elapsed: f64,
    ) -> Result<[MolecularWall3d; 6], DissipationError> {
        if !elapsed.is_finite() || elapsed < 0.0 {
            return Err(DissipationError::Invalid("molecular wall stage time"));
        }
        let mut walls = self.initial;
        if elapsed == 0.0 {
            return Ok(walls);
        }
        for (wall, rate) in walls.iter_mut().zip(self.rates) {
            if let MolecularWall3d::NoSlip { velocity, .. } = wall {
                for (component, acceleration) in velocity.iter_mut().zip(rate) {
                    if acceleration != 0.0 {
                        *component = acceleration.mul_add(elapsed, *component);
                        if !component.is_finite() {
                            return Err(DissipationError::Invalid(
                                "accelerated molecular wall velocity overflow",
                            ));
                        }
                    }
                }
            }
        }
        Ok(walls)
    }
}
