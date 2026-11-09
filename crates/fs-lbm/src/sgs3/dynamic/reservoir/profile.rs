//! Spatial primitive ghosts for the original transverse Germano filter.
//! A physical normal derivative uses its adjacent ghost. A test-level normal
//! derivative instead uses the SAME nine-cell density-weighted transverse
//! filter on the exterior plane; filtering interior values alone is incorrect.

use super::*;

/// Primitive state at a full-spacing exterior ghost centre. Units and scalar
/// chart must match the interior fields supplied to the combined preparation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReservoirPrimitive {
    pub density: f64,
    pub velocity: [f64; 3],
    pub scalar: f64,
}

pub(in crate::sgs3::dynamic) struct ProfileBoundary<'a> {
    shape: [usize; 3],
    normal: usize,
    sample: &'a dyn Fn(usize, usize) -> ReservoirPrimitive,
    velocity_reference: [f64; 3],
    scalar_reference: f64,
    pub velocity_extent: f64,
    pub scalar_extent: f64,
}

impl<'a> ProfileBoundary<'a> {
    fn new(
        shape: [usize; 3],
        normal: usize,
        sample: &'a dyn Fn(usize, usize) -> ReservoirPrimitive,
        velocity_reference: [f64; 3],
        scalar_reference: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, DynamicSgsError> {
        let mut result = Self {
            shape,
            normal,
            sample,
            velocity_reference,
            scalar_reference,
            velocity_extent: 0.0,
            scalar_extent: 0.0,
        };
        for side in 0..2 {
            for tile in 0..shape.iter().product::<usize>() / shape[normal] {
                poll(cancelled, tile)?;
                let cell = face_cell(shape, normal, side, tile);
                let ghost = sample(side, cell);
                if !(ghost.density.is_finite() && ghost.density > 0.0) {
                    return Err(DynamicSgsError::Invalid(
                        "positive finite reservoir density required",
                    ));
                }
                for a in 0..3 {
                    let difference = ghost.velocity[a] - velocity_reference[a];
                    if !difference.is_finite() {
                        return Err(DynamicSgsError::Invalid(
                            "finite representable reservoir velocity differences required",
                        ));
                    }
                    result.velocity_extent = result.velocity_extent.max(difference.abs());
                }
                let difference = ghost.scalar - scalar_reference;
                if !difference.is_finite() {
                    return Err(DynamicSgsError::Invalid(
                        "finite representable reservoir scalar differences required",
                    ));
                }
                result.scalar_extent = result.scalar_extent.max(difference.abs());
            }
        }
        Ok(result)
    }

    fn side(&self, cell: usize) -> Option<usize> {
        let p = coordinates(self.shape, cell)[self.normal];
        if p == 0 {
            Some(0)
        } else if p + 1 == self.shape[self.normal] {
            Some(1)
        } else {
            None
        }
    }

    fn value<const N: usize>(
        &self,
        cell: usize,
        side: usize,
        filtered: bool,
        primitive: impl Fn(ReservoirPrimitive) -> [f64; N],
    ) -> Result<[f64; N], DynamicSgsError> {
        if !filtered {
            return Ok(primitive((self.sample)(side, cell)));
        }
        let mut donors = [ReservoirPrimitive {
            density: 0.0,
            velocity: [0.0; 3],
            scalar: 0.0,
        }; 9];
        let mut count = 0;
        let mut density_scale = 0.0_f64;
        FilterGeometry::Channel(self.normal).each_box(self.shape, cell, |j| {
            let ghost = (self.sample)(side, j);
            density_scale = density_scale.max(ghost.density);
            donors[count] = ghost;
            count += 1;
        });
        let mut mass = 0.0;
        let mut mean = [0.0; N];
        for ghost in donors {
            let weight = ghost.density / density_scale / 9.0;
            if !(weight.is_finite() && weight > 0.0) {
                return Err(DynamicSgsError::Cell {
                    cell,
                    reason: "unrepresentable reservoir Favre weight",
                });
            }
            mass += weight;
            for (sum, value) in mean.iter_mut().zip(primitive(ghost)) {
                *sum += weight * value;
            }
        }
        for value in &mut mean {
            *value /= mass;
        }
        if mean.iter().any(|value| !value.is_finite()) {
            return Err(DynamicSgsError::Cell {
                cell,
                reason: "unrepresentable filtered reservoir primitive",
            });
        }
        Ok(mean)
    }

    pub fn velocity(
        &self,
        cell: usize,
        filtered: bool,
        scale: f64,
    ) -> Result<[[f64; 3]; 2], DynamicSgsError> {
        let mut values = [[0.0; 3]; 2];
        if let Some(side) = self.side(cell) {
            values[side] = self.value(cell, side, filtered, |ghost| {
                core::array::from_fn(|a| (ghost.velocity[a] - self.velocity_reference[a]) / scale)
            })?;
        }
        Ok(values)
    }

    pub fn scalar(
        &self,
        cell: usize,
        filtered: bool,
        scale: f64,
    ) -> Result<[ScalarWall; 2], DynamicSgsError> {
        let mut values = [ScalarWall::Value(0.0); 2];
        if let Some(side) = self.side(cell) {
            values[side] = ScalarWall::Value(
                self.value(cell, side, filtered, |ghost| {
                    [(ghost.scalar - self.scalar_reference) / scale]
                })?[0],
            );
        }
        Ok(values)
    }
}

fn face_cell(shape: [usize; 3], normal: usize, side: usize, tile: usize) -> usize {
    let axes = match normal {
        0 => [1, 2],
        1 => [0, 2],
        _ => [0, 1],
    };
    let mut p = [0; 3];
    p[normal] = side * (shape[normal] - 1);
    p[axes[0]] = tile % shape[axes[0]];
    p[axes[1]] = tile / shape[axes[0]];
    p[0] + shape[0] * (p[1] + shape[1] * p[2])
}

impl ReservoirDynamicSgs {
    /// Fit the current interior field and spatially varying exterior states.
    /// `sample(side, cell)` supplies the low/high ghost beside the original
    /// flattened boundary fluid cell. It must be repeatable and read-only.
    /// Every exterior tile is admitted before fitting, including exact-zero
    /// interior fields. Physical gradients use individual ghosts; test-level
    /// gradients use their nine-cell Favre averages. No halo arrays are stored.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_with_profile(
        &mut self,
        momentum: DynamicSgsConfig,
        heat: DynamicScalarConfig,
        density: &[f64],
        velocity: &[[f64; 3]],
        scalar: &[f64],
        sample: impl Fn(usize, usize) -> ReservoirPrimitive,
        scalar_workspace: &mut ReservoirDynamicScalar,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(DynamicSgsReport, DynamicScalarReport), DynamicSgsError> {
        poll(&mut cancelled, 0)?;
        momentum.validate()?;
        heat.validate()?;
        let n = self.inner.velocity.len();
        if density.len() != n
            || velocity.len() != n
            || scalar.len() != n
            || scalar_workspace.normal_axis != self.normal_axis
            || scalar_workspace.inner.shape != self.inner.shape
        {
            return Err(DynamicSgsError::Invalid(
                "reservoir profile and workspace geometry differ",
            ));
        }
        let boundary = ProfileBoundary::new(
            self.shape(),
            self.normal_axis,
            &sample,
            velocity[0],
            scalar[0],
            &mut cancelled,
        )?;
        let representative = core::array::from_fn::<_, 2, _>(|side| {
            sample(side, face_cell(self.shape(), self.normal_axis, side, 0))
        });
        let scalar_scale = scalar_workspace.inner.normalize_boundary(
            scalar,
            Some(representative.map(|g| ScalarWall::Value(g.scalar))),
            boundary.scalar_extent,
            &mut cancelled,
        )?;
        let (report, scales) = self.inner.prepare_scaled_with_profile(
            momentum,
            density,
            velocity,
            representative.map(|g| g.velocity),
            Some(&boundary),
            &mut cancelled,
        )?;
        let heat_report = scalar_workspace.inner.prepare_from_profile(
            &self.inner,
            heat,
            density,
            scales,
            scalar_scale,
            Some(&boundary),
            &mut cancelled,
        )?;
        Ok((report, heat_report))
    }
}

#[cfg(test)]
#[path = "profile/tests.rs"]
mod tests;
