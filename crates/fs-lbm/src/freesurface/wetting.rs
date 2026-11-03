//! Geometric planar wetting for the existing D2Q9 curvature operator.
//!
//! With n_w directed from solid into liquid and n=-grad(phi)/|grad(phi)|,
//! n.n_w=cos(theta) gives dphi/dn_w=-|grad_t(phi)| cot(theta). Extend the
//! frozen smoothed fill through the planar wall; do not clip ghost values,
//! which would change the prescribed derivative. No tracked mass is added.
//! This is the dense counterpart of the existing planar D3Q19 closure, not
//! a height-function boundary scheme. The material owner supplies the optional
//! source-pinned advancing/receding field and its explicit rate-law limits.

use super::{Cell, ContactModel, FreeSurface};
use std::collections::BTreeMap;
use fs_exec::CancelGate;

mod dynamic;
pub(crate) use dynamic::WallFill2;

/// A complete, stationary outer wall row/column in the dense D2Q9 lattice.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PlanarWall2 {
    /// Wall cells at x=0; inward normal +x.
    XMin,
    /// Wall cells at x=nx-1; inward normal -x.
    XMax,
    /// Wall cells at y=0; inward normal +y.
    YMin,
    /// Wall cells at y=ny-1; inward normal -y.
    YMax,
}

impl PlanarWall2 {
    const fn axis(self) -> usize {
        match self {
            Self::XMin | Self::XMax => 0,
            Self::YMin | Self::YMax => 1,
        }
    }
    const fn upper(self) -> bool {
        matches!(self, Self::XMax | Self::YMax)
    }
}

/// Refusal before changing a planar wetting boundary.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ContactAngleError2 {
    /// Angle must be finite, strictly in (0,pi), with a representable stencil.
    InvalidAngle,
    /// The selected face is periodic, incomplete, or lacks its inward stencil.
    InvalidPlanarWall,
    /// A dynamic field must cover exactly the admitted wall-adjacent row.
    InvalidField,
    /// Explicit cancellation before publishing a candidate field.
    Cancelled,
}
impl core::fmt::Display for ContactAngleError2 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "D2Q9 planar wetting refusal: {self:?}")
    }
}
impl std::error::Error for ContactAngleError2 {}

#[derive(Clone, Copy)]
pub(super) struct PlanarContactAngle2 {
    pub(super) wall: PlanarWall2,
    angle_radians: f64,
    cotangent: f64,
}

// Direct tangential indexing keeps the hot curvature/smoothing stencil O(1).
// Corner entries are absent; the other intersecting wall retains its closure.
#[derive(Clone)]
pub(super) struct PlanarContactField2 {
    pub(super) wall: PlanarWall2,
    angles: Vec<Option<PlanarContactAngle2>>,
}

impl FreeSurface {
    /// Set one equilibrium contact angle, measured through the liquid.
    /// The selected face must consist of wall cells and be nonperiodic along
    /// its normal. Other walls retain the neutral closure. Exactly pi/2 uses
    /// the original neutral arithmetic, including its spatial approximation.
    ///
    /// This changes capillary pressure reconstruction, never liquid mass.
    /// It does not prescribe contact-line speed, hysteresis or wall slip.
    pub fn set_planar_contact_angle(
        &mut self,
        wall: PlanarWall2,
        angle_radians: f64,
    ) -> Result<(), ContactAngleError2> {
        if !angle_radians.is_finite()
            || angle_radians <= 0.0
            || angle_radians >= core::f64::consts::PI
        {
            return Err(ContactAngleError2::InvalidAngle);
        }
        let cotangent = if angle_radians == core::f64::consts::FRAC_PI_2 {
            0.0
        } else {
            let (sin, cos) = fs_math::det::sin_cos(angle_radians);
            cos / sin
        };
        let bound = 4.0 * cotangent;
        if !bound.is_finite() || !(bound * bound).is_finite() {
            return Err(ContactAngleError2::InvalidAngle);
        }
        self.admit_planar_wall(wall)?;
        self.planar_contact = Some(PlanarContactAngle2 {
            wall,
            angle_radians,
            cotangent,
        });
        self.planar_contact_field = None;
        self.contact = ContactModel::Neutral;
        Ok(())
    }

    /// Explicitly assigned wall and angle in radians.
    #[must_use]
    pub fn planar_contact_angle(&self) -> Option<(PlanarWall2, f64)> {
        self.planar_contact.map(|c| (c.wall, c.angle_radians))
    }

    pub(super) fn admit_planar_wall(&self, wall: PlanarWall2) -> Result<(), ContactAngleError2> {
        let dims = [self.grid.nx, self.grid.ny];
        let periodic = [self.grid.periodic_x, self.grid.periodic_y];
        let axis = wall.axis();
        let tangent = 1 - axis;
        if periodic[axis] || dims[axis] < 4 || dims[tangent] < 3 {
            return Err(ContactAngleError2::InvalidPlanarWall);
        }
        let edge = if wall.upper() { dims[axis] - 1 } else { 0 };
        for t in 0..dims[tangent] {
            let mut p = [0; 2];
            p[axis] = edge;
            p[tangent] = t;
            if self.grid.flags[self.grid.idx(p[0], p[1])] != Cell::Wall {
                return Err(ContactAngleError2::InvalidPlanarWall);
            }
            // Tangential end walls meet at corners, whose compound geometry
            // is not assigned a new contact angle by this planar closure.
            if periodic[tangent] || (t > 0 && t + 1 < dims[tangent]) {
                p[axis] = if wall.upper() { edge - 1 } else { 1 };
                if self.grid.flags[self.grid.idx(p[0], p[1])] == Cell::Wall {
                    return Err(ContactAngleError2::InvalidPlanarWall);
                }
            }
        }
        Ok(())
    }

    pub(super) fn planar_angle_ghost(&self, position: [i64; 2]) -> Option<f64> {
        self.planar_angle_ghost_in(position, &self.fill_smooth)
    }

    // The same geometric extension must precede smoothing. Dropping wall
    // samples from that average would bend an otherwise planar interface
    // before its curvature stencil even sees the prescribed angle.
    pub(super) fn planar_angle_ghost_in(&self, position: [i64; 2], fill: &[f64]) -> Option<f64> {
        let contact = self.planar_contact_at(position)?;
        if contact.cotangent == 0.0 {
            return None;
        }
        let axis = contact.wall.axis();
        let tangent = 1 - axis;
        let dims = [self.grid.nx, self.grid.ny];
        let periodic = [self.grid.periodic_x, self.grid.periodic_y];
        let interior = if contact.wall.upper() {
            dims[axis] as i64 - 2
        } else {
            1
        };
        let distance = if contact.wall.upper() {
            position[axis] - interior
        } else {
            interior - position[axis]
        };
        if distance <= 0 {
            return None;
        }
        let sample = |t: i64| -> Option<f64> {
            let t = if periodic[tangent] {
                t.rem_euclid(dims[tangent] as i64)
            } else {
                t
            };
            if t < 0 || t >= dims[tangent] as i64 {
                return None;
            }
            let mut p = [0; 2];
            p[axis] = interior as usize;
            p[tangent] = t as usize;
            let i = self.grid.idx(p[0], p[1]);
            (self.grid.flags[i] != Cell::Wall).then_some(fill[i])
        };
        let t = position[tangent];
        let value = sample(t)?;
        let derivative = match (sample(t + 1), sample(t - 1)) {
            (Some(p), Some(m)) => 0.5 * (p - m),
            (Some(p), None) => p - value,
            (None, Some(m)) => value - m,
            (None, None) => 0.0,
        };
        Some(value + distance as f64 * derivative.abs() * contact.cotangent)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core2::Grid;
    use core::f64::consts::{FRAC_PI_2, PI};

    fn box_surface() -> FreeSurface {
        let mut grid = Grid::uniform(12, 12, 0.8);
        grid.periodic_x = false;
        grid.periodic_y = false;
        for y in 0..12 {
            for x in 0..12 {
                let i = grid.idx(x, y);
                grid.flags[i] = if x == 0 || y == 0 || x == 11 || y == 11 {
                    Cell::Wall
                } else {
                    Cell::Interface
                };
            }
        }
        FreeSurface::new(grid, 0.01, ContactModel::Neutral)
    }

    #[test]
    fn g1_planar_normal_and_zero_curvature_at_all_four_walls() {
        for wall in [
            PlanarWall2::XMin,
            PlanarWall2::XMax,
            PlanarWall2::YMin,
            PlanarWall2::YMax,
        ] {
            for theta in [PI / 3.0, 2.0 * PI / 3.0] {
                for slope in [-0.02_f64, 0.02] {
                    let mut flow = box_surface();
                    flow.set_planar_contact_angle(wall, theta).unwrap();
                    let axis = wall.axis();
                    let tangent = 1 - axis;
                    let normal_slope = -slope.abs() * theta.cos() / theta.sin();
                    // An independent smooth plane, rather than a curvature
                    // field generated by the operator under test.
                    for y in 1..11 {
                        for x in 1..11 {
                            let p = [x, y];
                            let distance = if wall.upper() {
                                10.0 - p[axis] as f64
                            } else {
                                p[axis] as f64 - 1.0
                            };
                            let i = flow.grid.idx(x, y);
                            let fill =
                                0.5 + slope * (p[tangent] as f64 - 6.0) + distance * normal_slope;
                            flow.mass[i] = fill * flow.grid.f[i].iter().sum::<f64>();
                        }
                    }
                    let mass = flow.ledger_mass();
                    let populations = flow.grid.f.clone();
                    // Exercise the actual frozen-fill/smoothing path too;
                    // neutral averaging would destroy this planar invariant.
                    flow.refresh_fill();
                    let mut p = [6; 2];
                    p[axis] = if wall.upper() { 10 } else { 1 };
                    for depth in [1, 2] {
                        let mut ghost = [p[0] as i64, p[1] as i64];
                        ghost[axis] += if wall.upper() { depth } else { -depth };
                        let value = flow.planar_angle_ghost(ghost).unwrap();
                        assert!((value - (0.5 - depth as f64 * normal_slope)).abs() < 2e-15);
                        let measured = (0.5 - value) / depth as f64;
                        let cosine = -measured / measured.hypot(slope);
                        assert!((cosine - theta.cos()).abs() < 2e-14);
                    }
                    assert!(flow.curvature(p[0], p[1]).abs() < 2e-14);
                    assert_eq!(flow.grid.f, populations);
                    assert_eq!(flow.ledger_mass().to_bits(), mass.to_bits());
                }
            }
        }
    }

    #[test]
    fn dynamic_geometry_tracks_fill_and_not_the_imposed_angle_field() {
        let gate = CancelGate::new_clock_free();
        for wall in [PlanarWall2::XMin, PlanarWall2::XMax, PlanarWall2::YMin, PlanarWall2::YMax] {
            let mut flow = box_surface();
            let axis = wall.axis();
            let tangent = 1 - axis;
            let theta = PI / 3.0;
            let slope: f64 = 0.02;
            let normal_slope = -slope * theta.cos() / theta.sin();
            for y in 1..11 {
                for x in 1..11 {
                    let p = [x, y];
                    let distance = if wall.upper() { 10 - p[axis] } else { p[axis] - 1 };
                    let fill = 0.5 + slope * (p[tangent] as f64 - 6.0)
                        + distance as f64 * normal_slope;
                    let i = flow.grid.idx(x, y);
                    flow.mass[i] = fill * flow.grid.f[i].iter().sum::<f64>();
                }
            }
            let original = flow.planar_wall_fill(wall, 10, &gate).unwrap();
            let mut p = [6; 2];
            p[axis] = if wall.upper() { 10 } else { 1 };
            assert!((original[&p].angle_radians - theta).abs() < 2e-14);
            assert!((original[&p].tangent_gradient - slope).abs() < 2e-15);
            for applied in [PI / 6.0, 5.0 * PI / 6.0] {
                flow.set_planar_contact_field(wall,
                    original.keys().map(|&p| (p, applied)).collect(), &gate).unwrap();
                flow.refresh_fill(); // Actual ghost-dependent curvature cache.
                let actual = flow.planar_wall_fill(wall, 10, &gate).unwrap();
                for (p, before) in &original {
                    let after = actual[p];
                    assert_eq!(before.fill.to_bits(), after.fill.to_bits());
                    assert_eq!(before.tangent_gradient.to_bits(), after.tangent_gradient.to_bits());
                    assert_eq!(before.angle_radians.to_bits(), after.angle_radians.to_bits());
                }
            }
        }
    }

    #[test]
    fn invalid_angle_or_geometry_preserves_the_previous_assignment() {
        let mut flow = box_surface();
        flow.set_planar_contact_angle(PlanarWall2::YMin, FRAC_PI_2)
            .unwrap();
        let before = flow.planar_contact_angle();
        for angle in [f64::NAN, f64::INFINITY, -1.0, 0.0, PI, f64::MIN_POSITIVE] {
            assert_eq!(
                flow.set_planar_contact_angle(PlanarWall2::XMin, angle),
                Err(ContactAngleError2::InvalidAngle)
            );
            assert_eq!(flow.planar_contact_angle(), before);
        }
        flow.grid.periodic_x = true;
        assert_eq!(
            flow.set_planar_contact_angle(PlanarWall2::XMin, 1.0),
            Err(ContactAngleError2::InvalidPlanarWall)
        );
        assert_eq!(flow.planar_contact_angle(), before);
        flow.grid.periodic_x = false;
        let gap = flow.grid.idx(0, 6);
        flow.grid.flags[gap] = Cell::Gas;
        assert_eq!(
            flow.set_planar_contact_angle(PlanarWall2::XMin, 1.0),
            Err(ContactAngleError2::InvalidPlanarWall)
        );
        assert_eq!(flow.planar_contact_angle(), before);
    }
}
