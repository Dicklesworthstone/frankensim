//! Shared exact cell crossings for accepted, piecewise-constant physical grids.

use super::{VolumeGrid, VolumeTransportError, checkpoint};
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};

/// Finite forward support in the native grid's original face convention.
/// Return entry distance and the rebased length, avoiding a second box law.
pub(crate) fn grid_forward_interval(
    grid: &VolumeGrid<'_>,
    origin: Point3,
    direction: Vec3,
) -> Result<Option<(f64, f64)>, VolumeTransportError> {
    validate_ray(origin, direction, f64::INFINITY)?;
    let Some(ray) = CellRay::clip(grid, origin, direction, f64::INFINITY, 1)? else {
        return Ok(None);
    };
    if !ray.end.is_finite() || ray.end <= 0.0 || ray.entry_m + ray.end <= ray.entry_m {
        return Err(VolumeTransportError::Unrepresentable);
    }
    Ok(Some((ray.entry_m, ray.end)))
}

/// Integrate the original piecewise-constant mass density to infinity.
/// The positive coefficient is m²/kg; every visited cell consumes its native
/// crossing budget. Zero rate and geometric misses require no cell visits.
pub(crate) fn grid_scattering_optical_depth(
    cx: &Cx<'_>,
    grid: &VolumeGrid<'_>,
    origin: Point3,
    direction: Vec3,
    mass_scattering_m2_kg: f64,
    max_cells: u32,
) -> Result<(f64, u32), VolumeTransportError> {
    grid_column_with_checkpoint(
        grid, origin, direction, mass_scattering_m2_kg, max_cells,
        &mut || checkpoint(cx),
    )
}

fn grid_column_with_checkpoint(
    grid: &VolumeGrid<'_>,
    origin: Point3,
    direction: Vec3,
    mass_scattering_m2_kg: f64,
    max_cells: u32,
    check: &mut impl FnMut() -> Result<(), VolumeTransportError>,
) -> Result<(f64, u32), VolumeTransportError> {
    check()?;
    validate_ray(origin, direction, f64::INFINITY)?;
    if !mass_scattering_m2_kg.is_finite() || mass_scattering_m2_kg < 0.0 {
        return Err(VolumeTransportError::ScatteringSpectrum);
    }
    if mass_scattering_m2_kg == 0.0 {
        return Ok((0.0, 0));
    }
    let Some(mut ray) = CellRay::clip(grid, origin, direction, f64::INFINITY, max_cells)? else {
        return Ok((0.0, 0));
    };
    if !ray.end.is_finite() || ray.end <= 0.0 {
        return Err(VolumeTransportError::Unrepresentable);
    }
    let (mut depth, mut correction, mut visits) = (0.0_f64, 0.0_f64, 0_u32);
    while let Some(segment) = ray.next_with_checkpoint(check)? {
        visits = visits.checked_add(1).ok_or(VolumeTransportError::Budget)?;
        let density = grid.data[segment.cell];
        if !density.is_finite() || density < 0.0 {
            return Err(VolumeTransportError::Density { index: segment.cell });
        }
        if density == 0.0 || segment.length_m == 0.0 {
            continue;
        }
        let increment = super::positive_product(density, mass_scattering_m2_kg, segment.length_m);
        if increment.is_nan() || increment <= 0.0 {
            return Err(VolumeTransportError::Unrepresentable);
        }
        // Keep traversing after a range failure so cancellation, invalid
        // input and work budgets remain exact. A finite column overflow must
        // not impersonate the separately declared infinite uniform tail.
        if depth == f64::INFINITY {
            continue;
        }
        let next = depth + increment;
        if next == f64::INFINITY {
            depth = next;
            correction = 0.0;
        } else {
            correction += if depth >= increment {
                (depth - next) + increment
            } else {
                (increment - next) + depth
            };
            depth = next;
        }
    }
    check()?;
    let total = depth + correction;
    if !total.is_finite() || total < 0.0 {
        return Err(VolumeTransportError::Unrepresentable);
    }
    Ok((total, visits))
}

pub(in crate::volumes) fn validate_ray(
    origin: Point3,
    direction: Vec3,
    distance: f64,
) -> Result<(), VolumeTransportError> {
    let o = [origin.x, origin.y, origin.z];
    let d = [direction.x, direction.y, direction.z];
    if o.iter().chain(d.iter()).any(|v| !v.is_finite())
        || (direction.dot(direction) - 1.0).abs() > 2.0e-10
        || distance.is_nan()
        || distance < 0.0
    {
        Err(VolumeTransportError::Ray)
    } else {
        Ok(())
    }
}

pub(in crate::volumes) struct CellSegment {
    pub cell: usize,
    /// Parameter from the rebased entry point, avoiding distant-camera loss.
    pub start_m: f64,
    pub length_m: f64,
}

pub(in crate::volumes) struct CellRay<'g, 'a> {
    grid: &'g VolumeGrid<'a>,
    start: [f64; 3],
    direction: [f64; 3],
    entry_m: f64,
    end: f64,
    index: [usize; 3],
    next: [f64; 3],
    current: f64,
    steps_left: u32,
    done: bool,
}

impl<'g, 'a> CellRay<'g, 'a> {
    /// The caller has validated the ray and immutable grid once already.
    pub fn clip(
        grid: &'g VolumeGrid<'a>,
        origin: Point3,
        direction: Vec3,
        distance: f64,
        max_steps: u32,
    ) -> Result<Option<Self>, VolumeTransportError> {
        if distance == 0.0 {
            return Ok(None);
        }
        let o = [origin.x, origin.y, origin.z];
        let d = [direction.x, direction.y, direction.z];
        let (lo, hi) = grid.bounds();
        let (mut enter, mut exit) = (0.0_f64, distance);
        let mut collapsed_interval = false;
        for a in 0..3 {
            if d[a] == 0.0 {
                if o[a] < lo[a] || o[a] >= hi[a] {
                    return Ok(None);
                }
            } else {
                let t1 = (lo[a] - o[a]) / d[a];
                let t2 = (hi[a] - o[a]) / d[a];
                collapsed_interval |= t1 == t2;
                enter = enter.max(t1.min(t2));
                exit = exit.min(t1.max(t2));
            }
        }
        if exit < enter {
            return Ok(None);
        }
        if exit == enter {
            return if collapsed_interval {
                Err(VolumeTransportError::Unrepresentable)
            } else {
                Ok(None)
            };
        }
        if !enter.is_finite() || !exit.is_finite() {
            return Err(VolumeTransportError::Unrepresentable);
        }
        let start =
            std::array::from_fn::<_, 3, _>(|a| d[a].mul_add(enter, o[a]).clamp(lo[a], hi[a]));
        let mut end = exit - enter;
        for a in 0..3 {
            if d[a] != 0.0 {
                let far = if d[a] > 0.0 { hi[a] } else { lo[a] };
                end = end.min((far - start[a]) / d[a]);
            }
        }
        let mut index = [0usize; 3];
        for a in 0..3 {
            let q = (start[a] - lo[a]) / grid.cell[a];
            let floor = q.floor();
            let cell = if d[a] < 0.0 && q == floor {
                floor - 1.0
            } else {
                floor
            };
            index[a] = cell.clamp(0.0, (grid.dims[a] - 1) as f64) as usize;
            // Division can round a point on a stored decimal cell face just
            // below (or above) its integer lattice coordinate. Reconcile the
            // fast estimate with the SAME face coordinates used by crossing.
            // Parallel/positive rays own the upper cell; negative rays own
            // the lower one. A rare mismatch uses a bounded binary search.
            let past_face = |boundary: usize| {
                let coordinate = (boundary as f64).mul_add(grid.cell[a], grid.origin[a]);
                if d[a] < 0.0 {
                    start[a] > coordinate
                } else {
                    start[a] >= coordinate
                }
            };
            if (index[a] > 0 && !past_face(index[a]))
                || (index[a] + 1 < grid.dims[a] && past_face(index[a] + 1))
            {
                let (mut lower, mut upper) = (0, grid.dims[a] - 1);
                while lower < upper {
                    let middle = lower + (upper - lower).div_ceil(2);
                    if past_face(middle) {
                        lower = middle;
                    } else {
                        upper = middle - 1;
                    }
                }
                index[a] = lower;
            }
        }
        let mut result = Self {
            grid,
            start,
            direction: d,
            entry_m: enter,
            end,
            index,
            next: [0.0; 3],
            current: 0.0,
            steps_left: max_steps,
            done: false,
        };
        result.next = std::array::from_fn(|a| result.crossing(a));
        Ok(Some(result))
    }

    fn crossing(&self, axis: usize) -> f64 {
        if self.direction[axis] == 0.0 {
            return f64::INFINITY;
        }
        let boundary = if self.direction[axis] > 0.0 {
            self.index[axis] + 1
        } else {
            self.index[axis]
        };
        let coordinate = (boundary as f64).mul_add(self.grid.cell[axis], self.grid.origin[axis]);
        (coordinate - self.start[axis]) / self.direction[axis]
    }

    pub fn next(&mut self, cx: &Cx<'_>) -> Result<Option<CellSegment>, VolumeTransportError> {
        self.next_with_checkpoint(&mut || checkpoint(cx))
    }

    pub fn next_with_checkpoint(
        &mut self,
        check: &mut (impl FnMut() -> Result<(), VolumeTransportError> + ?Sized),
    ) -> Result<Option<CellSegment>, VolumeTransportError> {
        if self.done {
            return Ok(None);
        }
        if self.steps_left == 0 {
            return Err(VolumeTransportError::Budget);
        }
        check()?;
        self.steps_left -= 1;
        let stop = self
            .next
            .into_iter()
            .fold(self.end, f64::min)
            .max(self.current);
        let result = CellSegment {
            cell: self.index[0]
                + self.grid.dims[0] * (self.index[1] + self.grid.dims[1] * self.index[2]),
            start_m: self.current,
            length_m: stop - self.current,
        };
        self.current = stop;
        if self.current >= self.end {
            self.done = true;
            return Ok(Some(result));
        }
        for a in 0..3 {
            if self.next[a] <= stop {
                if (self.direction[a] < 0.0 && self.index[a] == 0)
                    || (self.direction[a] > 0.0 && self.index[a] + 1 == self.grid.dims[a])
                {
                    return Err(VolumeTransportError::Unrepresentable);
                }
                if self.direction[a] > 0.0 {
                    self.index[a] += 1;
                } else {
                    self.index[a] -= 1;
                }
                self.next[a] = self.crossing(a);
            }
        }
        Ok(Some(result))
    }

    pub fn position(&self, rebased_distance_m: f64) -> Point3 {
        Point3::new(
            self.direction[0].mul_add(rebased_distance_m, self.start[0]),
            self.direction[1].mul_add(rebased_distance_m, self.start[1]),
            self.direction[2].mul_add(rebased_distance_m, self.start[2]),
        )
    }

    pub fn distance_from_origin(&self, rebased_distance_m: f64) -> f64 {
        self.entry_m + rebased_distance_m
    }
}

#[cfg(test)]
mod finite_support_tests;
