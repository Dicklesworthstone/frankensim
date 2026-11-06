//! Explicit nonlinear porous potentials on the original staggered faces.
//!
//! For the eight choices of u=(u_x±,u_y±,u_z±), the separable law uses
//! Phi = V/8 sum_choices sum_i c_i |r_i.u|^3/3, c_i=rho beta_i;
//! the isotropic vector law uses Phi = V/8 sum_choices c |u|^3/3.
//! Its derivative supplies each face force; its homogeneous cubic power is
//! 3 Phi. This positive quadrature retains face coercivity and reduces to
//! the original V/2 separable face law for signed grid-axis permutations of
//! the separable model. Equal principal coefficients do not select vector drag.

use super::*;

/// Explicit separable or isotropic vector-speed quadratic resistance.
/// Coefficients are rho*beta [kg/m4], with no factor one half. A separable
/// law retains its orthonormal material frame, even for equal coefficients;
/// the isotropic vector law depends on the Euclidean speed instead.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ForchheimerResistance {
    law: InertialLaw,
    grid: Option<[f64; 3]>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum InertialLaw {
    Principal {
        coefficients: [f64; 3],
        axes: [[f64; 3]; 3],
    },
    Isotropic(f64),
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

fn speed(u: [f64; 3]) -> Result<f64, ChtError> {
    let scale = u.into_iter().map(f64::abs).fold(0.0_f64, f64::max);
    if scale == 0.0 {
        return Ok(0.0);
    }
    let normalized = u.map(|value| value / scale);
    let speed = scale * fs_math::det::sqrt(dot(normalized, normalized));
    finite("simple.forchheimer_vector_speed", speed)?;
    Ok(speed)
}

impl ForchheimerResistance {
    /// Calibrated separable principal-direction coefficients and material frame.
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
            law: InertialLaw::Principal {
                coefficients: principal_kg_m4,
                axes,
            },
            grid,
        })
    }

    /// Isotropic `-c |u| u`, where c is rho*beta and |u| is Euclidean speed.
    /// Linear permeability may remain anisotropic. This law has no principal
    /// frame and never uses the separable grid-axis fast path.
    pub fn try_isotropic(coefficient_kg_m4: f64) -> Result<Self, ChtError> {
        if !coefficient_kg_m4.is_finite() || coefficient_kg_m4 < 0.0 {
            return Err(invalid("finite nonnegative isotropic rho*beta required"));
        }
        Ok(Self {
            law: InertialLaw::Isotropic(coefficient_kg_m4),
            grid: None,
        })
    }

    /// Coefficients only for the separable principal-direction convention.
    pub const fn principal_kg_m4(&self) -> Option<&[f64; 3]> {
        match &self.law {
            InertialLaw::Principal { coefficients, .. } => Some(coefficients),
            InertialLaw::Isotropic(_) => None,
        }
    }
    /// Material frame only for the separable principal-direction convention.
    pub const fn material_axes_in_world(&self) -> Option<&[[f64; 3]; 3]> {
        match &self.law {
            InertialLaw::Principal { axes, .. } => Some(axes),
            InertialLaw::Isotropic(_) => None,
        }
    }
    /// Scalar coefficient only for the Euclidean vector-speed convention.
    pub const fn isotropic_kg_m4(&self) -> Option<f64> {
        match self.law {
            InertialLaw::Isotropic(coefficient) => Some(coefficient),
            InertialLaw::Principal { .. } => None,
        }
    }
    pub(super) const fn grid_coefficients(&self) -> Option<[f64; 3]> {
        self.grid
    }
    pub(super) fn is_zero(&self) -> bool {
        match self.law {
            InertialLaw::Principal { coefficients, .. } => coefficients == [0.0; 3],
            InertialLaw::Isotropic(coefficient) => coefficient == 0.0,
        }
    }

    /// Maximum speed in this SAME cubature, using the law's actual convention:
    /// a principal component for separable drag, Euclidean speed for vector drag.
    pub(in crate::conjugate) fn maximum_calibrated_speed(
        &self,
        faces: [[f64; 2]; 3],
    ) -> Result<f64, ChtError> {
        if faces.iter().flatten().any(|v| !v.is_finite()) {
            return Err(invalid("nonfinite native face velocity"));
        }
        let InertialLaw::Principal { axes, .. } = self.law else {
            // Every coordinate can choose either face independently. Thus the
            // largest of eight Euclidean norms needs only this one norm.
            return speed(faces.map(|pair| pair[0].abs().max(pair[1].abs())));
        };
        if self.grid.is_some() {
            return Ok(faces.iter().flatten().map(|v| v.abs()).fold(0.0, f64::max));
        }
        let mut maximum: f64 = 0.0;
        for choice in 0..8 {
            let u = tuple(faces, choice);
            for axis in axes {
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
        match self.law {
            InertialLaw::Isotropic(c) => {
                if c != 0.0 {
                    for choice in 0..8 {
                        let side = (choice >> a) & 1;
                        row[side][0] += 0.125 * (c * speed(tuple(faces, choice))?);
                    }
                }
            }
            InertialLaw::Principal { coefficients, axes } => {
                for choice in 0..8 {
                    let u = tuple(faces, choice);
                    let side = (choice >> a) & 1;
                    for (c, r) in coefficients.into_iter().zip(axes) {
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
            }
        }
        for value in row.iter().flatten() {
            finite("simple.forchheimer_row", *value)?;
        }
        Ok(row)
    }

    fn power_density(&self, faces: [[f64; 2]; 3]) -> Result<f64, ChtError> {
        let mut power = 0.0;
        match self.law {
            InertialLaw::Isotropic(c) => {
                if c != 0.0 {
                    for choice in 0..8 {
                        let q = speed(tuple(faces, choice))?;
                        power += 0.125 * (c * q) * q * q;
                    }
                }
            }
            InertialLaw::Principal { coefficients, axes } => {
                for choice in 0..8 {
                    let u = tuple(faces, choice);
                    for (c, r) in coefficients.into_iter().zip(axes) {
                        if c == 0.0 {
                            continue;
                        }
                        let q = dot(r, u);
                        power += 0.125 * c * q.abs() * q * q;
                    }
                }
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
            match self.law {
                InertialLaw::Isotropic(c) => {
                    if c != 0.0 {
                        force += 0.125 * (c * speed(u)?) * u[a];
                    }
                }
                InertialLaw::Principal { coefficients, axes } => {
                    for (c, r) in coefficients.into_iter().zip(axes) {
                        if c == 0.0 {
                            continue;
                        }
                        let q = dot(r, u);
                        force += 0.125 * c * q.abs() * q * r[a];
                    }
                }
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
                "coupled Forchheimer cells need closed normal boundaries or a clear-fluid layer before pressure/inlet openings",
            ));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 2e-12 * expected.abs().max(1.0),
            "{actual:e} != {expected:e}"
        );
    }

    #[test]
    fn g0_vector_law_uses_euclidean_speed_and_rotates_covariantly() {
        let law = ForchheimerResistance::try_isotropic(2.0).unwrap();
        assert_eq!(law.principal_kg_m4(), None);
        assert_eq!(law.material_axes_in_world(), None);
        assert_eq!(law.isotropic_kg_m4(), Some(2.0));
        assert_eq!(law.grid_coefficients(), None);
        let u = [3.0, 4.0, 0.0];
        let rotate = |v: [f64; 3]| [0.6 * v[0] - 0.8 * v[1], 0.8 * v[0] + 0.6 * v[1], v[2]];
        let original = u.map(|value| [value; 2]);
        let rotated = rotate(u).map(|value| [value; 2]);
        let expected = [30.0, 40.0, 0.0];
        near(law.maximum_calibrated_speed(original).unwrap(), 5.0);
        near(law.power_density(original).unwrap(), 250.0);
        near(law.power_density(rotated).unwrap(), 250.0);
        for a in 0..3 {
            for side in 0..2 {
                near(
                    2.0 * law.face_force_density(original, a, side).unwrap(),
                    expected[a],
                );
                near(
                    2.0 * law.face_force_density(rotated, a, side).unwrap(),
                    rotate(expected)[a],
                );
            }
        }
        let separable = ForchheimerResistance::try_new(
            [2.0; 3],
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        )
        .unwrap();
        near(separable.power_density(original).unwrap(), 182.0);
        assert_ne!(separable.power_density(rotated).unwrap(), 182.0);
    }

    #[test]
    fn g0_native_force_is_the_derivative_of_the_positive_cubic_potential() {
        let law = ForchheimerResistance::try_isotropic(7.0).unwrap();
        let faces = [[-3.0, 2.0], [1.0, -4.0], [0.25, -0.5]];
        // Independent quadrature and Euclidean norm, without the native rows.
        let potential = |f: [[f64; 2]; 3]| {
            let mut total = 0.0;
            for x in f[0] {
                for y in f[1] {
                    for z in f[2] {
                        total += 7.0 * (x * x + y * y + z * z).powf(1.5) / 24.0;
                    }
                }
            }
            total
        };
        near(law.power_density(faces).unwrap(), 3.0 * potential(faces));
        let mut work = 0.0;
        for a in 0..3 {
            let row = law.row(faces, a).unwrap();
            for side in 0..2 {
                let force = law.face_force_density(faces, a, side).unwrap();
                near(row[side][0] * faces[a][side] + row[side][1], force);
                let mut plus = faces;
                let mut minus = faces;
                plus[a][side] += 1e-5;
                minus[a][side] -= 1e-5;
                let derivative = (potential(plus) - potential(minus)) / 2e-5;
                assert!((derivative - force).abs() < 1e-7);
                work += force * faces[a][side];
            }
        }
        near(work, law.power_density(faces).unwrap());
        let reversed = faces.map(|pair| pair.map(|value| -value));
        near(law.power_density(reversed).unwrap(), work);
        // Opposing faces cancel in a cell average but keep their real drag.
        let cancelling = [[-0.6, 0.6], [-0.8, 0.8], [0.0; 2]];
        near(law.maximum_calibrated_speed(cancelling).unwrap(), 1.0);
        near(law.power_density(cancelling).unwrap(), 7.0);
    }

    #[test]
    fn g0_zero_vector_coefficient_and_admission_are_exact() {
        let zero = ForchheimerResistance::try_isotropic(0.0).unwrap();
        let faces = [[-3.0, 2.0], [1.0, -4.0], [0.25, -0.5]];
        assert!(zero.is_zero());
        assert_eq!(zero.power_density(faces).unwrap(), 0.0);
        for a in 0..3 {
            assert_eq!(zero.row(faces, a).unwrap(), [[0.0; 2]; 2]);
        }
        let law = ForchheimerResistance::try_isotropic(7.0).unwrap();
        assert_eq!(law.power_density([[0.0; 2]; 3]).unwrap(), 0.0);
        for invalid in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(ForchheimerResistance::try_isotropic(invalid).is_err());
        }
        assert!(law.maximum_calibrated_speed([[f64::NAN; 2]; 3]).is_err());
        near(speed([3e-200, 4e-200, 0.0]).unwrap() / 1e-200, 5.0);
        near(speed([3e200, 4e200, 0.0]).unwrap() / 1e200, 5.0);
        let tiny = f64::from_bits(1);
        let law = ForchheimerResistance::try_isotropic(tiny).unwrap();
        let power = law.power_density([[1e100; 2], [0.0; 2], [0.0; 2]]).unwrap();
        near(power / (tiny * 1e100 * 1e100 * 1e100), 1.0);
    }
}
