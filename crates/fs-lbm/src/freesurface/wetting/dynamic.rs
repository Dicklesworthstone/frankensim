//! Complete wall history and per-cell angle ghosts for dense planar wetting.

use super::*;

/// Tangentially smoothed first-cell geometry; gradients are per lattice cell.
#[derive(Clone, Copy, Debug)]
pub(crate) struct WallFill2 {
    pub fill: f64,
    pub tangent_gradient: f64,
    pub angle_radians: f64,
}

fn check(gate: &CancelGate) -> Result<(), ContactAngleError2> {
    if gate.is_requested() {
        Err(ContactAngleError2::Cancelled)
    } else {
        Ok(())
    }
}

impl FreeSurface {
    pub(crate) fn admits_dynamic_planar_wetting(&self) -> bool {
        self.contact == ContactModel::Neutral
            && self.planar_contact.is_none()
            && self.planar_contact_field.is_none()
    }

    /// Preserve gas coordinates as well as liquid/interface cells so advancing
    /// and receding fronts use the same accepted geometric history.
    pub(crate) fn planar_wall_fill(
        &self,
        wall: PlanarWall2,
        maximum_cells: usize,
        gate: &CancelGate,
    ) -> Result<BTreeMap<[usize; 2], WallFill2>, ContactAngleError2> {
        check(gate)?;
        self.admit_planar_wall(wall)?;
        let axis = wall.axis();
        let tangent = 1 - axis;
        let dims = [self.grid.nx, self.grid.ny];
        let periodic = [self.grid.periodic_x, self.grid.periodic_y];
        let interior = if wall.upper() { dims[axis] - 2 } else { 1 };
        let needed = dims[tangent] - if periodic[tangent] { 0 } else { 2 };
        if needed > maximum_cells {
            return Err(ContactAngleError2::InvalidField);
        }
        // Measure physical fill, never the angle-dependent curvature ghosts.
        // Otherwise changing a boundary angle would manufacture a contact-line
        // speed even with identical populations and mass. Tangential smoothing
        // preserves affine planes and needs only O(wall cells) work, independent
        // of the volume occupied by this dense lattice.
        let sample = |p: [usize; 2]| {
            let i = self.grid.idx(p[0], p[1]);
            if self.grid.flags[i] == Cell::Wall {
                return None;
            }
            let mut sum = self.fill(i);
            let mut count = 1.0;
            for offset in [-1, 1] {
                if let Some(t) =
                    super::super::offset_coord(p[tangent], offset, dims[tangent], periodic[tangent])
                {
                    let mut q = p;
                    q[tangent] = t;
                    let j = self.grid.idx(q[0], q[1]);
                    if self.grid.flags[j] != Cell::Wall {
                        sum += self.fill(j);
                        count += 1.0;
                    }
                }
            }
            Some(sum / count)
        };
        let mut result = BTreeMap::new();
        for t in 0..dims[tangent] {
            check(gate)?;
            let mut p = [0; 2];
            p[axis] = interior;
            p[tangent] = t;
            let Some(fill) = sample(p) else { continue };
            let at = |offset| {
                let mut q = p;
                q[tangent] =
                    super::super::offset_coord(t, offset, dims[tangent], periodic[tangent])?;
                sample(q)
            };
            let derivative = match (at(1), at(-1)) {
                (Some(a), Some(b)) => (a - b) / 2.0,
                (Some(a), None) => a - fill,
                (None, Some(b)) => fill - b,
                (None, None) => 0.0,
            };
            let tangent_gradient = derivative.abs();
            let mut inward = p;
            inward[axis] = if wall.upper() {
                interior - 1
            } else {
                interior + 1
            };
            let normal_gradient =
                sample(inward).ok_or(ContactAngleError2::InvalidPlanarWall)? - fill;
            let angle_radians = if tangent_gradient == 0.0 {
                core::f64::consts::FRAC_PI_2
            } else {
                fs_math::det::atan2(tangent_gradient, -normal_gradient)
            };
            if !fill.is_finite() || !tangent_gradient.is_finite() || !angle_radians.is_finite() {
                return Err(ContactAngleError2::InvalidField);
            }
            if result.len() == maximum_cells {
                return Err(ContactAngleError2::InvalidField);
            }
            result.insert(
                p,
                WallFill2 {
                    fill,
                    tangent_gradient,
                    angle_radians,
                },
            );
        }
        Ok(result)
    }

    pub(crate) fn set_planar_contact_field(
        &mut self,
        wall: PlanarWall2,
        angles: BTreeMap<[usize; 2], f64>,
        gate: &CancelGate,
    ) -> Result<(), ContactAngleError2> {
        check(gate)?;
        if self.contact != ContactModel::Neutral
            || self.planar_contact.is_some()
            || self
                .planar_contact_field
                .as_ref()
                .is_some_and(|f| f.wall != wall)
        {
            return Err(ContactAngleError2::InvalidField);
        }
        self.admit_planar_wall(wall)?;
        let axis = wall.axis();
        let tangent = 1 - axis;
        let dims = [self.grid.nx, self.grid.ny];
        let interior = if wall.upper() { dims[axis] - 2 } else { 1 };
        let mut admitted = vec![None; dims[tangent]];
        for (p, angle_radians) in angles {
            check(gate)?;
            if p[axis] != interior
                || p[tangent] >= dims[tangent]
                || self.grid.flags[self.grid.idx(p[0], p[1])] == Cell::Wall
                || !angle_radians.is_finite()
                || angle_radians <= 0.0
                || angle_radians >= core::f64::consts::PI
            {
                return Err(ContactAngleError2::InvalidField);
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
            admitted[p[tangent]] = Some(PlanarContactAngle2 {
                wall,
                angle_radians,
                cotangent,
            });
        }
        for (t, value) in admitted.iter().enumerate() {
            check(gate)?;
            let mut p = [0; 2];
            p[axis] = interior;
            p[tangent] = t;
            if value.is_some() != (self.grid.flags[self.grid.idx(p[0], p[1])] != Cell::Wall) {
                return Err(ContactAngleError2::InvalidField);
            }
        }
        check(gate)?;
        self.planar_contact_field = Some(PlanarContactField2 {
            wall,
            angles: admitted,
        });
        Ok(())
    }

    pub(super) fn planar_contact_at(&self, position: [i64; 2]) -> Option<PlanarContactAngle2> {
        let Some(field) = &self.planar_contact_field else {
            return self.planar_contact;
        };
        let tangent = 1 - field.wall.axis();
        let dims = [self.grid.nx, self.grid.ny];
        let periodic = [self.grid.periodic_x, self.grid.periodic_y];
        let t = if periodic[tangent] {
            position[tangent].rem_euclid(dims[tangent] as i64)
        } else {
            position[tangent]
        };
        field
            .angles
            .get(usize::try_from(t).ok()?)
            .copied()
            .flatten()
    }
}
