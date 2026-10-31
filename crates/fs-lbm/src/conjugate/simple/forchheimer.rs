//! A material-frame nonlinear potential on the original staggered faces.
//!
//! For the eight choices of u=(u_x±,u_y±,u_z±), use
//! Phi = V/8 sum_choices sum_i c_i |r_i.u|^3/3, c_i=rho beta_i.
//! Its derivative supplies each face force; its homogeneous cubic power is
//! 3 Phi. This positive quadrature retains face coercivity and reduces to
//! the original V/2 separable face law for signed grid-axis permutations.
//! It is NOT the rotationally invariant isotropic vector-speed closure.

use super::*;

/// Separable principal-direction quadratic resistance with its actual frame.
/// Coefficients are rho*beta [kg/m4], with no factor one half. Rows of the
/// frame are orthonormal material axes in world coordinates. Equal principal
/// coefficients do not erase that frame for this nonlinear constitutive law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ForchheimerResistance {
    principal_kg_m4: [f64; 3],
    axes: [[f64; 3]; 3],
    grid: Option<[f64; 3]>,
}

fn invalid(reason: &str) -> ChtError {
    ChtError::InvalidInput {
        field: "simple.forchheimer_resistance",
        reason: reason.into(),
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn tuple(faces: [[f64; 2]; 3], choice: usize) -> [f64; 3] {
    std::array::from_fn(|a| faces[a][(choice >> a) & 1])
}

impl ForchheimerResistance {
    pub fn try_new(
        principal_kg_m4: [f64; 3],
        material_axes_in_world: [[f64; 3]; 3],
    ) -> Result<Self, ChtError> {
        let axes = material_axes_in_world;
        if principal_kg_m4.iter().any(|c| !c.is_finite() || *c < 0.0)
            || axes
                .iter()
                .any(|a| a.iter().any(|v| !v.is_finite()) || (dot(*a, *a) - 1.0).abs() > 2e-12)
            || dot(axes[0], axes[1]).abs() > 2e-12
            || dot(axes[0], axes[2]).abs() > 2e-12
            || dot(axes[1], axes[2]).abs() > 2e-12
        {
            return Err(invalid(
                "finite nonnegative principal rho*beta and orthonormal material axes required",
            ));
        }
        let mut grid = Some([0.0; 3]);
        for (i, axis) in axes.iter().enumerate() {
            if let Some(a) = axis.iter().position(|v| v.abs() == 1.0)
                && axis.iter().enumerate().all(|(b, v)| a == b || *v == 0.0)
            {
                grid.as_mut().expect("grid candidate")[a] = principal_kg_m4[i];
            } else {
                grid = None;
                break;
            }
        }
        Ok(Self {
            principal_kg_m4,
            axes,
            grid,
        })
    }

    pub const fn principal_kg_m4(&self) -> &[f64; 3] {
        &self.principal_kg_m4
    }
    pub const fn material_axes_in_world(&self) -> &[[f64; 3]; 3] {
        &self.axes
    }
    pub(super) const fn grid_coefficients(&self) -> Option<[f64; 3]> {
        self.grid
    }
    pub(super) fn is_zero(&self) -> bool {
        self.principal_kg_m4 == [0.0; 3]
    }

    /// Maximum of every principal velocity in this SAME constitutive cubature.
    pub(in crate::conjugate) fn maximum_component_speed(
        &self,
        faces: [[f64; 2]; 3],
    ) -> Result<f64, ChtError> {
        if faces.iter().flatten().any(|v| !v.is_finite()) {
            return Err(invalid("nonfinite native face velocity"));
        }
        if self.grid.is_some() {
            return Ok(faces.iter().flatten().map(|v| v.abs()).fold(0.0, f64::max));
        }
        let mut maximum: f64 = 0.0;
        for choice in 0..8 {
            let u = tuple(faces, choice);
            for axis in self.axes {
                let q = dot(axis, u);
                finite("simple.forchheimer_principal_speed", q)?;
                maximum = maximum.max(q.abs());
            }
        }
        Ok(maximum)
    }

    // Per unit physical cell volume: [implicit diagonal, explicit cross force]
    // for this component's minus and plus faces. Freeze q at the current
    // iterate; only a nonnegative diagonal enters the scalar native matrix.
    fn row(&self, faces: [[f64; 2]; 3], a: usize) -> Result<[[f64; 2]; 2], ChtError> {
        let mut row = [[0.0; 2]; 2];
        for choice in 0..8 {
            let u = tuple(faces, choice);
            let side = (choice >> a) & 1;
            for (c, r) in self.principal_kg_m4.into_iter().zip(self.axes) {
                if c == 0.0 {
                    continue;
                }
                let q = dot(r, u);
                let coefficient = c * q.abs();
                let cross: f64 = (0..3).filter(|b| *b != a).map(|b| r[b] * u[b]).sum();
                row[side][0] += 0.125 * coefficient * r[a] * r[a];
                row[side][1] += 0.125 * coefficient * r[a] * cross;
            }
        }
        for value in row.iter().flatten() {
            finite("simple.forchheimer_row", *value)?;
        }
        Ok(row)
    }

    fn power_density(&self, faces: [[f64; 2]; 3]) -> Result<f64, ChtError> {
        let mut power = 0.0;
        for choice in 0..8 {
            let u = tuple(faces, choice);
            for (c, r) in self.principal_kg_m4.into_iter().zip(self.axes) {
                if c == 0.0 {
                    continue;
                }
                let q = dot(r, u);
                power += 0.125 * c * q.abs() * q * q;
            }
        }
        finite("simple.forchheimer_power_density", power)?;
        Ok(power)
    }

    // Independent derivative traversal for the momentum-work heat audit.
    fn face_force_density(
        &self,
        faces: [[f64; 2]; 3],
        a: usize,
        side: usize,
    ) -> Result<f64, ChtError> {
        let mut force = 0.0;
        for choice in 0..8 {
            if (choice >> a) & 1 != side {
                continue;
            }
            let u = tuple(faces, choice);
            for (c, r) in self.principal_kg_m4.into_iter().zip(self.axes) {
                if c == 0.0 {
                    continue;
                }
                let q = dot(r, u);
                force += 0.125 * c * q.abs() * q * r[a];
            }
        }
        finite("simple.forchheimer_face_force_density", force)?;
        Ok(force)
    }
}

impl Solver<'_> {
    fn forchheimer_faces(&self, cell: usize) -> [[f64; 2]; 3] {
        let at = self.domain.coords(cell);
        std::array::from_fn(|a| {
            std::array::from_fn(|side| self.vel[a][self.cell_face(a, at, side == 1)])
        })
    }

    pub(super) fn refresh_forchheimer_rows(
        &mut self,
        a: usize,
        gate: &CancelGate,
    ) -> Result<(), ChtError> {
        for cell in 0..self.porous_framed.len() {
            if cell % 256 == 0 {
                poll(gate)?;
            }
            let Some(law) = self.porous_framed[cell] else {
                continue;
            };
            let mut row = law.row(self.forchheimer_faces(cell), a)?;
            for value in row.iter_mut().flatten() {
                *value *= self.domain.volume(cell);
                finite("simple.forchheimer_face_row", *value)?;
            }
            self.porous_framed_rows[cell] = row;
        }
        Ok(())
    }

    pub(super) fn forchheimer_face_row(&self, minus: [usize; 3], plus: [usize; 3]) -> [f64; 2] {
        let m = self.porous_framed_rows[cell_of(self.n, minus)][1];
        let p = self.porous_framed_rows[cell_of(self.n, plus)][0];
        [m[0] + p[0], m[1] + p[1]]
    }

    pub(super) fn framed_forchheimer_power(&self, cell: usize) -> Result<f64, ChtError> {
        let Some(Some(law)) = self.porous_framed.get(cell) else {
            return Ok(0.0);
        };
        let power = self.domain.volume(cell) * law.power_density(self.forchheimer_faces(cell))?;
        finite("simple.forchheimer_cell_power", power)?;
        Ok(power)
    }

    pub(super) fn framed_forchheimer_force(
        &self,
        a: usize,
        minus: [usize; 3],
        plus: [usize; 3],
    ) -> Result<f64, ChtError> {
        let mut force = 0.0;
        for (at, side) in [(minus, 1), (plus, 0)] {
            let cell = cell_of(self.n, at);
            if let Some(Some(law)) = self.porous_framed.get(cell) {
                force += self.domain.volume(cell)
                    * law.face_force_density(self.forchheimer_faces(cell), a, side)?;
            }
        }
        finite("simple.forchheimer_face_force", force)?;
        Ok(force)
    }

    pub(super) fn admit_forchheimer_opening(&self, cell: usize) -> Result<(), ChtError> {
        if (0..6).any(|side| {
            self.domain.neighbor(cell, side).is_none()
                && matches!(
                    self.faces[side],
                    FvBoundary::Inlet { .. } | FvBoundary::Outlet | FvBoundary::Pressure { .. }
                )
        }) {
            return Err(invalid(
                "rotated Forchheimer cells need closed normal boundaries or a clear-fluid layer before pressure/inlet openings",
            ));
        }
        Ok(())
    }
}
