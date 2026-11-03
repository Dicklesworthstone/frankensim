//! Sustained, unregularized normal/Coulomb laws in force and displacement units.
//!
//! A receiving mechanical owner supplies the endpoint gap and the actual
//! tangential displacement accumulated under a held contact force. This module
//! evaluates normal configuration complementarity and the isotropic Coulomb
//! maximum-dissipation inclusion. It does not construct a proxy body, convert
//! impulses to forces, integrate motion, or claim global uniqueness.
//!
//! The positive natural-map scales condition the numerical equations in N/m.
//! They introduce no spring, stored energy, compliance, or material parameter.
//! Physical admission still requires the receiving owner to check its original
//! trajectory, contact mode, work, momentum, and energy balances independently.

use fs_material::cohesive::friction::project_tangent_disk;
use std::fmt;

/// Coupled held-force solve through the original mechanical response.
pub mod solve;

#[cfg(test)]
mod tests;

/// One contact's force and the path produced by that force.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactPoint {
    /// Applied force in N, ordered normal, tangent 1, tangent 2.
    ///
    /// Positive normal force opens the gap. Tangential force opposes slip.
    /// Signed normal trial values are admitted during the nonlinear solve.
    pub force_n: [f64; 3],
    /// Endpoint normal gap in m; nonnegative is geometrically admissible.
    pub gap_m: f64,
    /// Actual accumulated tangential displacement in m in the force's frame.
    pub slip_increment_m: [f64; 2],
    /// Finite nonnegative isotropic Coulomb coefficient.
    pub friction_coefficient: f64,
}

/// Numerical force/displacement conversion scales, not physical stiffnesses.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NaturalMapScales {
    /// Strictly positive finite normal scale in N/m.
    pub normal_n_per_m: f64,
    /// Strictly positive finite common tangent scale in N/m.
    ///
    /// A common scale preserves the Euclidean disk's rotational symmetry.
    pub tangent_n_per_m: f64,
}

impl NaturalMapScales {
    pub(crate) fn validate(self) -> Result<(), ContactLawError> {
        if !self.normal_n_per_m.is_finite()
            || self.normal_n_per_m <= 0.0
            || !self.tangent_n_per_m.is_finite()
            || self.tangent_n_per_m <= 0.0
        {
            return Err(ContactLawError::Invalid(
                "natural-map scales must be finite and strictly positive in N/m",
            ));
        }
        Ok(())
    }
}

/// A force residual and one fixed generalized derivative at its base point.
///
/// Both derivative matrices use normal, tangent 1, tangent 2 ordering. Each
/// Krylov action uses this same selection, including at a contact-mode kink.
/// The path derivative acts on endpoint gap and accumulated tangent slip; it
/// is composed with the original mechanical response by the receiving owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactLinearization {
    /// Natural-map residual in N.
    pub residual_n: [f64; 3],
    /// Derivative with respect to force, holding the mechanical path fixed.
    pub force_derivative: [[f64; 3]; 3],
    /// Derivative with respect to [gap, slip 1, slip 2] in m, holding force fixed.
    pub path_derivative: [[f64; 3]; 3],
}

/// A malformed contact input or unrepresentable numerical operation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ContactLawError {
    /// Inputs violate the declared domain.
    Invalid(&'static str),
    /// A derived residual, projection, or derivative is not finite.
    NonFinite(&'static str),
}

impl fmt::Display for ContactLawError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "invalid sustained-contact input: {reason}"),
            Self::NonFinite(operation) => {
                write!(f, "non-finite sustained-contact operation: {operation}")
            }
        }
    }
}

impl std::error::Error for ContactLawError {}

/// Evaluate the unregularized force/displacement natural map.
///
/// The normal equation is `fn - max(fn - kn * gap, 0) = 0`. Its zeros are
/// `fn >= 0`, `gap >= 0`, and `fn * gap = 0`. The tangent equation is
/// `ft - project_disk(ft - kt * slip, mu * max(fn, 0)) = 0`, using the
/// material owner's shared Euclidean disk projection. Its zeros enforce
/// sticking within the disk or a boundary force opposite nonzero slip.
///
/// At a normal tie this selects the inactive derivative; at a disk boundary
/// it selects the exterior derivative. At the collapsed disk's zero apex it
/// selects zero projection derivative. These choices remain fixed throughout
/// the linearized solve; projecting a Krylov direction would not be linear.
///
/// # Errors
/// Refuses non-finite inputs, negative friction, nonpositive scales, and
/// nonrepresentable projections, residuals, or derivative coefficients.
pub fn linearize(
    point: ContactPoint,
    scales: NaturalMapScales,
) -> Result<ContactLinearization, ContactLawError> {
    scales.validate()?;
    if point.force_n.iter().any(|value| !value.is_finite())
        || !point.gap_m.is_finite()
        || point
            .slip_increment_m
            .iter()
            .any(|value| !value.is_finite())
        || !point.friction_coefficient.is_finite()
        || point.friction_coefficient < 0.0
    {
        return Err(ContactLawError::Invalid(
            "finite force, gap, slip and nonnegative finite friction required",
        ));
    }

    let normal_trial = (-scales.normal_n_per_m).mul_add(point.gap_m, point.force_n[0]);
    let tangent_trial = [
        (-scales.tangent_n_per_m).mul_add(point.slip_increment_m[0], point.force_n[1]),
        (-scales.tangent_n_per_m).mul_add(point.slip_increment_m[1], point.force_n[2]),
    ];
    let capacity = point.friction_coefficient * point.force_n[0].max(0.0);
    if !normal_trial.is_finite() {
        return Err(ContactLawError::NonFinite("normal trial"));
    }
    let projected = project_tangent_disk(tangent_trial, capacity)
        .map_err(|_| ContactLawError::NonFinite("tangent disk projection"))?;
    let length = tangent_trial[0].hypot(tangent_trial[1]);
    let mut result = ContactLinearization {
        residual_n: [0.0; 3],
        force_derivative: [[0.0; 3]; 3],
        path_derivative: [[0.0; 3]; 3],
    };

    if normal_trial > 0.0 {
        // Algebraically fn - (fn - kn * gap), without cancellation at a
        // loaded equilibrium where the force can dwarf the residual.
        result.residual_n[0] = scales.normal_n_per_m * point.gap_m;
        result.path_derivative[0][0] = scales.normal_n_per_m;
    } else {
        result.residual_n[0] = point.force_n[0];
        result.force_derivative[0][0] = 1.0;
    }
    if length < capacity {
        // Preserve tiny physical slip residuals inside a heavily loaded disk.
        for i in 0..2 {
            result.residual_n[i + 1] = scales.tangent_n_per_m * point.slip_increment_m[i];
            result.path_derivative[i + 1][i + 1] = scales.tangent_n_per_m;
        }
    } else if length > 0.0 {
        let unit = [tangent_trial[0] / length, tangent_trial[1] / length];
        let fraction = capacity / length;
        let capacity_derivative = if point.force_n[0] > 0.0 {
            point.friction_coefficient
        } else {
            0.0
        };
        for i in 0..2 {
            result.residual_n[i + 1] = point.force_n[i + 1] - projected[i];
            result.force_derivative[i + 1][0] = -capacity_derivative * unit[i];
            for (j, unit_j) in unit.iter().enumerate() {
                let identity = if i == j { 1.0 } else { 0.0 };
                let projection_derivative = fraction * (identity - unit[i] * unit_j);
                result.force_derivative[i + 1][j + 1] = identity - projection_derivative;
                result.path_derivative[i + 1][j + 1] =
                    scales.tangent_n_per_m * projection_derivative;
            }
        }
    } else {
        // The only remaining case is the collapsed disk at its apex.
        for i in 1..3 {
            result.residual_n[i] = point.force_n[i];
            result.force_derivative[i][i] = 1.0;
        }
    }
    if result.residual_n.iter().any(|value| !value.is_finite())
        || result
            .force_derivative
            .iter()
            .chain(&result.path_derivative)
            .flatten()
            .any(|value| !value.is_finite())
    {
        return Err(ContactLawError::NonFinite("residual or fixed derivative"));
    }
    Ok(result)
}

impl ContactLinearization {
    /// Apply this fixed derivative to force and path directions.
    ///
    /// The path direction is the original mechanical response derivative in
    /// metres per the direction's parameter, not a sampled projection branch.
    ///
    /// # Errors
    /// Refuses non-finite directions, externally modified derivative entries,
    /// or an unrepresentable result.
    pub fn action(
        &self,
        force_direction_n: [f64; 3],
        path_direction_m: [f64; 3],
    ) -> Result<[f64; 3], ContactLawError> {
        if force_direction_n
            .iter()
            .chain(&path_direction_m)
            .chain(self.force_derivative.iter().flatten())
            .chain(self.path_derivative.iter().flatten())
            .any(|value| !value.is_finite())
        {
            return Err(ContactLawError::Invalid(
                "finite directions and derivative entries required",
            ));
        }
        let mut result = [0.0; 3];
        for (row, value) in result.iter_mut().enumerate() {
            for column in 0..3 {
                *value =
                    self.force_derivative[row][column].mul_add(force_direction_n[column], *value);
                *value =
                    self.path_derivative[row][column].mul_add(path_direction_m[column], *value);
            }
            if !value.is_finite() {
                return Err(ContactLawError::NonFinite("fixed derivative action"));
            }
        }
        Ok(result)
    }
}
