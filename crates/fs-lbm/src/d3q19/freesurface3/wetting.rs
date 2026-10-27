//! Planar contact geometry and spatial angle ghosts for the native VOF owner.

use super::*;

#[derive(Clone)]
pub(super) struct PlanarContactField3 {
    wall: PlanarWall3,
    angles: BTreeMap<[usize; 3], PlanarContactAngle3>,
}

/// Smoothed, wall-adjacent VOF geometry. Gradients use lattice cell coordinates.
/// This is a first-cell contact-line estimate, not a subcell wall intersection.
#[derive(Clone, Copy, Debug)]
pub(crate) struct WallFill3 {
    pub fill: f64,
    pub tangent_gradient: f64,
    pub angle_radians: f64,
}

impl FreeSurface3 {
    pub(crate) fn admits_dynamic_planar_wetting(&self) -> bool {
        self.contact == ContactModel3::Neutral && self.planar_contact.is_none()
    }

    /// The complete wall-adjacent fill history includes gas cells. A later
    /// activated tile was gas, so a missing history coordinate means zero fill.
    pub(crate) fn planar_wall_fill(
        &mut self,
        wall: PlanarWall3,
        maximum_cells: usize,
        gate: &CancelGate,
    ) -> Result<BTreeMap<[usize; 3], WallFill3>, FreeSurfaceError3> {
        self.refresh_fill(gate)?;
        self.planar_wall_fill_from_cache(wall, maximum_cells, gate)
    }

    fn planar_wall_fill_from_cache(
        &self,
        wall: PlanarWall3,
        maximum_cells: usize,
        gate: &CancelGate,
    ) -> Result<BTreeMap<[usize; 3], WallFill3>, FreeSurfaceError3> {
        let axis = wall.axis();
        let dims = self.grid.tile_dims();
        let extents = [dims.0 * TILE, dims.1 * TILE, dims.2 * TILE];
        let edge = if wall.upper() { extents[axis] - 1 } else { 0 };
        let mut result = BTreeMap::new();
        let sample = |position: [i64; 3]| -> Option<f64> {
            let (slot, lane) = self
                .grid
                .resolve_source(position[0], position[1], position[2])?;
            Some(self.fill_smooth[&self.grid.active_keys()[slot]][lane])
        };
        for &key in self.tiles.keys() {
            check_cancel(gate)?;
            for lane in 0..TILE_CELLS {
                let (x, y, z) = Self::global_of(key, lane);
                let p = [x, y, z];
                if p[axis] != edge as i64 {
                    continue;
                }
                if result.len() == maximum_cells {
                    return Err(FreeSurfaceError3::Transport(
                        "planar wetting exceeds its wall-cell budget",
                    ));
                }
                let fill = self.fill_smooth[&key][lane];
                let mut tangent_squared = 0.0;
                for tangent in 0..3 {
                    if tangent == axis {
                        continue;
                    }
                    let mut plus = p;
                    let mut minus = p;
                    plus[tangent] += 1;
                    minus[tangent] -= 1;
                    let derivative = match (sample(plus), sample(minus)) {
                        (Some(a), Some(b)) => (a - b) / 2.0,
                        (Some(a), None) => a - fill,
                        (None, Some(b)) => fill - b,
                        (None, None) => 0.0,
                    };
                    tangent_squared += derivative * derivative;
                }
                let tangent_gradient = tangent_squared.sqrt();
                let mut inward = p;
                inward[axis] += if wall.upper() { -1 } else { 1 };
                let normal_gradient = sample(inward).ok_or(FreeSurfaceError3::Transport(
                    "planar wetting requires its inward first-cell stencil",
                ))? - fill;
                // With n = -grad(phi)/|grad(phi)| and n_w into the domain,
                // cos(theta) = -d(phi)/d(n_w) / |grad(phi)|. atan2 avoids
                // division by a vanishing total gradient. A zero tangent has
                // no resolved contact line. Its stationary source point is
                // retained so a collapsed bracket reproduces a uniform wall.
                let angle_radians = if tangent_gradient == 0.0 {
                    core::f64::consts::FRAC_PI_2
                } else {
                    fs_math::det::atan2(tangent_gradient, -normal_gradient)
                };
                if !fill.is_finite() || !tangent_gradient.is_finite() || !angle_radians.is_finite()
                {
                    return Err(FreeSurfaceError3::Transport(
                        "nonfinite planar contact geometry",
                    ));
                }
                result.insert(
                    p.map(|v| v as usize),
                    WallFill3 {
                        fill,
                        tangent_gradient,
                        angle_radians,
                    },
                );
            }
        }
        Ok(result)
    }

    /// Install only a complete validated angle field. The existing ghost,
    /// curvature, pressure reconstruction, collision and mass routes consume it.
    pub(crate) fn set_planar_contact_field(
        &mut self,
        wall: PlanarWall3,
        angles: BTreeMap<[usize; 3], f64>,
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceError3> {
        if !self.admits_dynamic_planar_wetting() {
            return Err(FreeSurfaceError3::InvalidContactAngle);
        }
        let axis = wall.axis();
        let dims = self.grid.tile_dims();
        let extents = [dims.0 * TILE, dims.1 * TILE, dims.2 * TILE];
        let edge = if wall.upper() { extents[axis] - 1 } else { 0 };
        let mut admitted = BTreeMap::new();
        for (p, angle_radians) in angles {
            check_cancel(gate)?;
            if p[axis] != edge
                || (0..3).any(|d| p[d] >= extents[d])
                || !angle_radians.is_finite()
                || angle_radians <= 0.0
                || angle_radians >= core::f64::consts::PI
            {
                return Err(FreeSurfaceError3::InvalidContactAngle);
            }
            let cotangent = if angle_radians == core::f64::consts::FRAC_PI_2 {
                0.0
            } else {
                let (sin, cos) = fs_math::det::sin_cos(angle_radians);
                cos / sin
            };
            let bound = 4.0 * cotangent;
            if !bound.is_finite() || !(bound * bound).is_finite() {
                return Err(FreeSurfaceError3::InvalidContactAngle);
            }
            admitted.insert(
                p,
                PlanarContactAngle3 {
                    wall,
                    angle_radians,
                    cotangent,
                },
            );
        }
        check_cancel(gate)?;
        self.planar_contact_field = Some(PlanarContactField3 {
            wall,
            angles: admitted,
        });
        Ok(())
    }

    pub(super) fn planar_contact_at(&self, mut position: [i64; 3]) -> Option<PlanarContactAngle3> {
        let Some(field) = &self.planar_contact_field else {
            return self.planar_contact;
        };
        let axis = field.wall.axis();
        let dims = self.grid.tile_dims();
        let extents = [dims.0 * TILE, dims.1 * TILE, dims.2 * TILE];
        let edge = if field.wall.upper() {
            extents[axis] as i64 - 1
        } else {
            0
        };
        if field.wall.upper() && position[axis] <= edge
            || !field.wall.upper() && position[axis] >= 0
        {
            return None;
        }
        position[axis] = edge;
        let p = [
            usize::try_from(position[0]).ok()?,
            usize::try_from(position[1]).ok()?,
            usize::try_from(position[2]).ok()?,
        ];
        field.angles.get(&p).copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g1_wall_fill_angle_and_per_cell_ghost_follow_all_six_normal_conventions() {
        for wall in [
            PlanarWall3::XMin,
            PlanarWall3::XMax,
            PlanarWall3::YMin,
            PlanarWall3::YMax,
            PlanarWall3::ZMin,
            PlanarWall3::ZMax,
        ] {
            for angle in [
                core::f64::consts::PI / 3.0,
                2.0 * core::f64::consts::PI / 3.0,
            ] {
                for sign in [-1.0, 1.0] {
                    let mut surface = FreeSurface3::new(
                        8,
                        8,
                        8,
                        0.8,
                        [0.0; 3],
                        0.001,
                        ContactModel3::Neutral,
                        |_, _, _| true,
                    )
                    .unwrap();
                    let axis = wall.axis();
                    let tangent = (axis + 1) % 3;
                    let slope: f64 = 0.03;
                    for &key in surface.grid.active_keys() {
                        surface.fill_smooth.insert(
                            key,
                            core::array::from_fn(|lane| {
                                let (x, y, z) = FreeSurface3::global_of(key, lane);
                                let p = [x, y, z];
                                let inward = if wall.upper() { 7 - p[axis] } else { p[axis] };
                                0.5 - sign * slope * p[tangent] as f64
                                    - slope * (angle.cos() / angle.sin()) * inward as f64
                            }),
                        );
                    }
                    let gate = CancelGate::new_clock_free();
                    let geometry = surface
                        .planar_wall_fill_from_cache(wall, 64, &gate)
                        .unwrap();
                    let mut p = [3; 3];
                    p[axis] = if wall.upper() { 7 } else { 0 };
                    let actual = geometry[&p];
                    assert!((actual.tangent_gradient - slope).abs() < 1e-14);
                    assert!(
                        (actual.angle_radians - angle).abs() < 2e-12,
                        "{wall:?} {sign}"
                    );
                    surface
                        .set_planar_contact_field(wall, BTreeMap::from([(p, angle)]), &gate)
                        .unwrap();
                    let mut ghost = p.map(|v| v as i64);
                    ghost[axis] += if wall.upper() { 1 } else { -1 };
                    let value = surface.planar_angle_ghost(ghost).unwrap();
                    assert!(
                        (value - actual.fill - slope * angle.cos() / angle.sin()).abs() < 2e-12
                    );
                }
            }
        }
    }
}
