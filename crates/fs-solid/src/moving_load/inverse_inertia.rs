//! Coupled point inverse inertia from the original constrained shell mass.

use super::{MovingLoadError, MovingShellMotion, ShellLoadPoint, finite, invalid, linear, poll};

/// Immutable `B M_ff^-1 B^T` for one ordered set of reference material points.
///
/// Each six-component row/column block follows the corresponding entry of
/// `points()`: Cartesian x/y/z translation followed by x/y/z rotation.
/// Columns receive linear impulse (N s) then free angular impulse about the
/// point (N m s); rows return material velocity change (m/s) then angular
/// velocity change (rad/s). Blocks therefore have mixed physical units.
///
/// Cross-point and translation/rotation couplings are retained. Constrained
/// and redundant directions may be singular; no regularization, diagonal
/// replacement or artificial stiffness is added. Reciprocity holds to the
/// rounding of the original native mass solve, without symmetrizing entries.
///
/// This is the linear response BEFORE adding increments to an accepted
/// velocity. `apply_impulses` still checks representability of that addition.
/// The operator can be reused while the model and ordered material points
/// stay fixed, including through prescribed pin motion. It retains the
/// original reference-frame, small-rotation shell scope.
#[derive(Clone, Debug, PartialEq)]
pub struct ShellPointInverseInertia {
    points: Vec<ShellLoadPoint>,
    values: Vec<f64>,
}

impl ShellPointInverseInertia {
    /// Original ordered material locations, including their reference offsets.
    pub fn points(&self) -> &[ShellLoadPoint] {
        &self.points
    }

    /// Matrix dimension: six times the admitted point count.
    pub fn dimension(&self) -> usize {
        6 * self.points.len()
    }

    /// Complete row-major response, including both off-diagonal block halves.
    pub fn values(&self) -> &[f64] {
        &self.values
    }

    /// Apply simultaneous point impulses in this operator's six-component
    /// ordering, returning every point's coupled linear/angular response.
    pub fn apply(&self, impulses: &[[f64; 6]]) -> Result<Vec<[f64; 6]>, MovingLoadError> {
        if impulses.len() != self.points.len() {
            return Err(invalid("point inverse inertia impulse dimension mismatch"));
        }
        for impulse in impulses {
            finite(impulse)?;
        }
        let dimension = self.dimension();
        let mut response = vec![[0.0; 6]; self.points.len()];
        for (row, value) in response.iter_mut().flatten().enumerate() {
            *value = self.values[row * dimension..(row + 1) * dimension]
                .iter()
                .zip(impulses.iter().flatten())
                .map(|(a, j)| a * j)
                .sum();
        }
        for value in &response {
            finite(value)?;
        }
        Ok(response)
    }
}

impl MovingShellMotion<'_> {
    /// Prepare the coupled point inverse inertia using the retained mass LU.
    ///
    /// Every unit point impulse uses the original force/impulse transpose;
    /// its native free-mass response uses the same material sampler as
    /// observation. Payload inertia and all cross-point blocks are retained.
    /// Prescribed velocities have zero event increment and supply no affine
    /// term here. Their physical reaction/work belongs to `apply_impulses`.
    ///
    /// Point count is bounded by the model's original load-point budget.
    /// Checked matrix dimensions/bytes and the original shell matrix-entry
    /// budget bound dense storage before allocation. Preparation costs six
    /// native mass solves per point, with no new factorization; callers can
    /// prepare once and reuse the immutable operator.
    ///
    /// Only solver scratch is borrowed mutably. Accepted state, endpoint
    /// forces, supports and reactions stay unchanged, including on refusal
    /// or cancellation. No inverse of this possibly singular response is
    /// inferred, and no collision or actuator law is selected.
    pub fn point_inverse_inertia(
        &mut self,
        points: &[ShellLoadPoint],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ShellPointInverseInertia, MovingLoadError> {
        poll(&mut cancelled)?;
        let dimension = points
            .len()
            .checked_mul(6)
            .ok_or_else(|| invalid("point inverse inertia dimension overflow"))?;
        let entries = dimension
            .checked_mul(dimension)
            .ok_or_else(|| invalid("point inverse inertia matrix size overflow"))?;
        let bytes = entries
            .checked_mul(core::mem::size_of::<f64>())
            .ok_or_else(|| invalid("point inverse inertia byte size overflow"))?;
        if points.len() > self.model.max_load_points
            || entries > self.model.plate.budget.max_matrix_entries
            || bytes > isize::MAX as usize
        {
            return Err(invalid(
                "point inverse inertia exceeds original shell budgets",
            ));
        }
        let mut nodes = Vec::with_capacity(points.len());
        for &point in points {
            poll(&mut cancelled)?;
            nodes.push(self.model.point(point)?.0);
        }
        let mut values = Vec::new();
        values
            .try_reserve_exact(entries)
            .map_err(|_| invalid("point inverse inertia allocation refused"))?;
        values.resize(entries, 0.0);
        let mut velocity = vec![0.0; self.state.velocity.len()];
        for column in 0..dimension {
            poll(&mut cancelled)?;
            let mut impulse = [0.0; 6];
            impulse[column % 6] = 1.0;
            let load = self.model.project_components(
                1,
                core::iter::once((
                    points[column / 6],
                    [impulse[0], impulse[1], impulse[2]],
                    [impulse[3], impulse[4], impulse[5]],
                )),
                &mut cancelled,
            )?;
            self.mass_solver
                .solve_factored_into(&load.reduced, &mut velocity)
                .map_err(linear)?;
            finite(&velocity)?;
            for (index, (&point, &nodes)) in points.iter().zip(&nodes).enumerate() {
                let response = self
                    .model
                    .material_components(point, nodes, &velocity, None);
                finite(&response)?;
                for axis in 0..6 {
                    values[(6 * index + axis) * dimension + column] = response[axis];
                }
            }
        }
        poll(&mut cancelled)?;
        Ok(ShellPointInverseInertia {
            points: points.to_vec(),
            values,
        })
    }
}
