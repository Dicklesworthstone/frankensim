//! Work-conjugate coordinates for a sustained directional Coulomb ellipse.
//!
//! For physical tangential force `f = A q`, the conjugate displacement is
//! `s' = A^T s`. The ellipse becomes the disk `|q| <= mu_max * fn`, and
//! `f dot s = q dot s'`. This changes no inertia, compliance or dissipation.
//! The receiving owner must use the complete transformed physical trace and
//! retain its original trajectory and work admission.

use super::ContactLawError;
use crate::frictional_multicontact::elliptic::EllipticFriction;

/// An invertible tangent map; normal coordinates are unchanged.
///
/// Equal coefficients use exact identity, including the frictionless case.
/// A single zero coefficient is deliberately refused: a free-slip axis needs
/// a lower-rank law, and must not acquire fictitious traction through an inverse.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EllipticForceCoordinates {
    force_map: [[f64; 2]; 2],
    inverse_force_map: [[f64; 2]; 2],
    coefficient: f64,
    minimum_scale: f64,
}

impl EllipticForceCoordinates {
    /// Prepare the physical ellipse `R diag(mu_1, mu_2) unit_disk * fn`.
    /// The first principal direction is normalized without overflowing its norm.
    pub fn prepare(friction: EllipticFriction) -> Result<Self, ContactLawError> {
        let mu = friction.coefficients_mu;
        let axis = friction.principal_axis;
        if mu.iter().any(|v| !v.is_finite() || *v < 0.0) || axis.iter().any(|v| !v.is_finite()) {
            return Err(ContactLawError::Invalid(
                "finite nonnegative ellipse coefficients and finite axis required",
            ));
        }
        let scale = axis[0].abs().max(axis[1].abs());
        if scale == 0.0 {
            return Err(ContactLawError::Invalid(
                "nonzero ellipse principal axis required",
            ));
        }
        let identity = [[1.0, 0.0], [0.0, 1.0]];
        if mu[0] == mu[1] {
            return Ok(Self {
                force_map: identity,
                inverse_force_map: identity,
                coefficient: mu[0],
                minimum_scale: 1.0,
            });
        }
        let coefficient = mu[0].max(mu[1]);
        let a = mu[0] / coefficient;
        let b = mu[1] / coefficient;
        if a == 0.0 || b == 0.0 || !(1.0 / a).is_finite() || !(1.0 / b).is_finite() {
            return Err(ContactLawError::Invalid(
                "sustained ellipse requires two representable positive principal coefficients",
            ));
        }
        let scaled = [axis[0] / scale, axis[1] / scale];
        let norm = scaled[0].hypot(scaled[1]);
        let u = [scaled[0] / norm, scaled[1] / norm];
        let v = [-u[1], u[0]];
        // Keep the principal directions as separate columns. Forming a
        // symmetric R D R^T matrix would erase a small positive axis through
        // cancellation when the coefficients are strongly directional.
        let force_map = [[a * u[0], b * v[0]], [a * u[1], b * v[1]]];
        let inverse_force_map = [[u[0] / a, u[1] / a], [v[0] / b, v[1] / b]];
        if force_map
            .iter()
            .chain(&inverse_force_map)
            .flatten()
            .any(|v| !v.is_finite())
        {
            return Err(ContactLawError::NonFinite("ellipse force coordinate map"));
        }
        Ok(Self {
            force_map,
            inverse_force_map,
            coefficient,
            minimum_scale: a.min(b),
        })
    }

    /// Map from solver tangent forces to physical tangent forces.
    pub fn force_map(self) -> [[f64; 2]; 2] {
        self.force_map
    }
    /// Coulomb disk coefficient in solver coordinates.
    pub fn coefficient(self) -> f64 {
        self.coefficient
    }
    /// Smallest singular value, for conservative physical path tolerances.
    pub fn minimum_scale(self) -> f64 {
        self.minimum_scale
    }

    /// Map a normal/two-tangent physical force into solver coordinates [N].
    pub fn force_to_solver(self, force: [f64; 3]) -> Result<[f64; 3], ContactLawError> {
        apply(self.inverse_force_map, force)
    }
    /// Map a normal/two-tangent solver force into physical coordinates [N].
    pub fn force_to_physical(self, force: [f64; 3]) -> Result<[f64; 3], ContactLawError> {
        apply(self.force_map, force)
    }
    /// Map physical gap/slip or velocity into conjugate solver coordinates.
    pub fn path_to_solver(self, path: [f64; 3]) -> Result<[f64; 3], ContactLawError> {
        apply(transpose(self.force_map), path)
    }
    /// Recover physical gap/slip from conjugate solver coordinates.
    pub fn path_to_physical(self, path: [f64; 3]) -> Result<[f64; 3], ContactLawError> {
        apply(transpose(self.inverse_force_map), path)
    }
}

fn transpose(map: [[f64; 2]; 2]) -> [[f64; 2]; 2] {
    [[map[0][0], map[1][0]], [map[0][1], map[1][1]]]
}

fn apply(map: [[f64; 2]; 2], value: [f64; 3]) -> Result<[f64; 3], ContactLawError> {
    if value.iter().any(|v| !v.is_finite()) {
        return Err(ContactLawError::Invalid(
            "finite force/path coordinates required",
        ));
    }
    if map == [[1.0, 0.0], [0.0, 1.0]] {
        return Ok(value);
    }
    let result = [
        value[0],
        map[0][0].mul_add(value[1], map[0][1] * value[2]),
        map[1][0].mul_add(value[1], map[1][1] * value[2]),
    ];
    if result.iter().any(|v| !v.is_finite()) {
        return Err(ContactLawError::NonFinite(
            "ellipse force/path transformation",
        ));
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn work_conjugacy_and_ellipse_support_in_a_rotated_frame() {
        let friction = EllipticFriction {
            coefficients_mu: [0.8, 0.2],
            principal_axis: [3.0, 4.0],
        };
        let c = EllipticForceCoordinates::prepare(friction).unwrap();
        let slip = [0.0, 0.13, -0.07];
        let transformed = c.path_to_solver(slip).unwrap();
        let norm = transformed[1].hypot(transformed[2]);
        let q = [
            5.0,
            -4.0 * transformed[1] / norm,
            -4.0 * transformed[2] / norm,
        ];
        let force = c.force_to_physical(q).unwrap();
        let work = force[1] * slip[1] + force[2] * slip[2];
        assert!((work + 4.0 * norm).abs() < 1e-15);
        let principal_force = [
            0.6 * force[1] + 0.8 * force[2],
            -0.8 * force[1] + 0.6 * force[2],
        ];
        assert!(((principal_force[0] / 4.0).hypot(principal_force[1] / 1.0) - 1.0).abs() < 1e-15);
        let recovered = c.force_to_solver(force).unwrap();
        for i in 0..3 {
            assert!((q[i] - recovered[i]).abs() < 3e-15);
        }
        assert!(
            force[1] * slip[2] - force[2] * slip[1] != 0.0,
            "directional friction generally is not antiparallel to slip"
        );
    }

    #[test]
    fn isotropic_is_exact_identity_and_singular_axes_are_not_regularized() {
        for mu in [0.0, 0.4] {
            let c = EllipticForceCoordinates::prepare(EllipticFriction {
                coefficients_mu: [mu, mu],
                principal_axis: [f64::MAX, f64::MAX],
            })
            .unwrap();
            assert_eq!(c.force_map(), [[1.0, 0.0], [0.0, 1.0]]);
            assert_eq!(
                c.force_to_solver([1.0, 2.0, -3.0]).unwrap(),
                [1.0, 2.0, -3.0]
            );
        }
        for mu in [[0.4, 0.0], [0.0, 0.4], [-0.1, 0.4], [f64::NAN, 0.4]] {
            assert!(
                EllipticForceCoordinates::prepare(EllipticFriction {
                    coefficients_mu: mu,
                    principal_axis: [1.0, 0.0]
                })
                .is_err()
            );
        }
    }
}
