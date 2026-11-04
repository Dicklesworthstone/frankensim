//! Linear elasticity: isotropic and orthotropic small-strain laws. The
//! tangent IS the stiffness (state-free), so consistency is exact by
//! construction — these are the baseline the FD gate calibrates itself
//! against.

use crate::{MaterialAdmissibility, MaterialError, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::{Ambition, ModelCard, ValidityDomain};

/// Source-derived reversible thermal response of the original isotropic law.
pub mod thermal;

/// Wave and hydrostatic response from the original admitted elastic tensor.
pub mod acoustic;

/// Mandel stiffness from the independent bulk/shear pair. Keep the original
/// moduli: an E/nu round trip loses bulk response near incompressibility and
/// shear response near the auxetic limit. Admission belongs to the caller.
pub(crate) fn bulk_shear_stiffness(bulk: f64, shear: f64) -> Tangent6 {
    let lambda = bulk - (2.0 / 3.0) * shear;
    let mut stiffness = [[0.0; 6]; 6];
    for (i, row) in stiffness.iter_mut().enumerate().take(3) {
        row[..3].fill(lambda);
        row[i] = bulk + (4.0 / 3.0) * shear;
    }
    for (i, row) in stiffness.iter_mut().enumerate().skip(3) {
        row[i] = 2.0 * shear;
    }
    stiffness
}

/// Isotropic linear elasticity (E, ν).
#[derive(Debug, Clone, PartialEq)]
pub struct IsotropicElastic {
    /// Young's modulus (Pa).
    pub youngs: f64,
    /// Poisson ratio.
    pub poisson: f64,
    /// Calibrated strain magnitude bound (validity domain).
    pub strain_limit: f64,
}

impl IsotropicElastic {
    /// Construct with physical-bounds checking.
    ///
    /// # Errors
    /// [`MaterialError::Parameters`] for E ≤ 0 or ν outside (−1, 0.5).
    pub fn new(youngs: f64, poisson: f64, strain_limit: f64) -> Result<Self, MaterialError> {
        if !(youngs > 0.0 && youngs.is_finite()) {
            return Err(MaterialError::Parameters {
                what: format!("E = {youngs} must be positive"),
            });
        }
        if !(poisson > -1.0 && poisson < 0.5) {
            return Err(MaterialError::Parameters {
                what: format!("nu = {poisson} outside (-1, 0.5)"),
            });
        }
        Ok(IsotropicElastic {
            youngs,
            poisson,
            strain_limit,
        })
    }

    /// Lamé parameters (λ, μ).
    #[must_use]
    pub fn lame(&self) -> (f64, f64) {
        let e = self.youngs;
        let nu = self.poisson;
        let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let mu = e / (2.0 * (1.0 + nu));
        (lambda, mu)
    }

    /// The constant stiffness (Voigt, tensor shear convention).
    #[must_use]
    pub fn stiffness(&self) -> Tangent6 {
        let (lambda, mu) = self.lame();
        let mut c = [[0.0f64; 6]; 6];
        for (i, row) in c.iter_mut().enumerate().take(3) {
            for slot in row.iter_mut().take(3) {
                *slot = lambda;
            }
            row[i] += 2.0 * mu;
        }
        for (i, row) in c.iter_mut().enumerate().skip(3) {
            row[i] = 2.0 * mu; // tensor shear: σ_xy = 2μ ε_xy
        }
        c
    }
}

fn apply(c: &Tangent6, strain: &Voigt) -> Voigt {
    let mut s = [0.0f64; 6];
    for (row, out) in c.iter().zip(s.iter_mut()) {
        // Plain σ = C·ε matvec. This crate uses the TENSOR-shear Voigt
        // convention, so the factor-2 lives in the 2μ shear entries of C (built
        // in `isotropic`/`orthotropic`), NOT here — the shear strain columns
        // must NOT be doubled. Do NOT "fix" this into `2·row[i]·strain[i]` for
        // i ≥ 3; that reintroduces a factor-2 error on every shear stress.
        *out = row[0] * strain[0]
            + row[1] * strain[1]
            + row[2] * strain[2]
            + row[3] * strain[3]
            + row[4] * strain[4]
            + row[5] * strain[5];
    }
    s
}

impl SmallStrainLaw for IsotropicElastic {
    type State = ();

    fn initial_state(&self) {}

    fn stress(&self, strain: &Voigt, (): &Self::State) -> Voigt {
        apply(&self.stiffness(), strain)
    }

    fn tangent(&self, _strain: &Voigt, (): &Self::State) -> Tangent6 {
        self.stiffness()
    }

    fn update_state(&self, _strain: &Voigt, (): &Self::State) {}

    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            has_stored_energy: true,
            dissipation_nonnegative: true, // zero dissipation
            polyconvex: Some(true),        // quadratic energy
            tangent_symmetric: true,
            failure_envelope: "none declared (linear law; validity ends at strain_limit)",
        }
    }

    fn card(&self) -> ModelCard {
        ModelCard::new(
            "material.isotropic-elastic",
            "0.1.0",
            Ambition::Solid,
            vec![
                "linear kinematics (small strain)".to_string(),
                "isotropy".to_string(),
                "no rate effects".to_string(),
            ],
            ValidityDomain::unconstrained().with("strain-magnitude", 0.0, self.strain_limit),
            vec!["silently linear beyond yield of the real material".to_string()],
            0.01,
        )
    }
}

/// Orthotropic linear elasticity (engineering constants, principal axes
/// aligned with the coordinate frame).
#[derive(Debug, Clone, PartialEq)]
pub struct OrthotropicElastic {
    /// Young's moduli (E1, E2, E3).
    pub e: [f64; 3],
    /// Poisson ratios (ν12, ν13, ν23).
    pub nu: [f64; 3],
    /// Shear moduli (G12, G23, G31).
    pub g: [f64; 3],
    /// Calibrated strain magnitude bound.
    pub strain_limit: f64,
}

impl OrthotropicElastic {
    /// Construct with positivity and positive-definite compliance checks.
    ///
    /// # Errors
    /// [`MaterialError::Parameters`] on non-positive moduli or an
    /// indefinite compliance.
    pub fn new(
        e: [f64; 3],
        nu: [f64; 3],
        g: [f64; 3],
        strain_limit: f64,
    ) -> Result<Self, MaterialError> {
        if e.iter()
            .chain(g.iter())
            .any(|&v| !(v > 0.0 && v.is_finite()))
        {
            return Err(MaterialError::Parameters {
                what: "all moduli must be positive".to_string(),
            });
        }
        let law = OrthotropicElastic {
            e,
            nu,
            g,
            strain_limit,
        };
        // Congruence by diag(sqrt(E_i)) removes stress units before testing
        // definiteness. Raw 1/E minors cube the stress scale and can underflow
        // to zero or overflow to infinity for the same physical compliance.
        let (_, lower) = law.normal_cholesky();
        if lower.iter().flatten().any(|value| !value.is_finite())
            || (0..3).any(|i| lower[i][i] <= 0.0)
        {
            return Err(MaterialError::Parameters {
                what: "orthotropic compliance is not positive definite \
                       (Poisson ratios thermodynamically inadmissible)"
                    .to_string(),
            });
        }
        if law
            .stiffness()
            .iter()
            .flatten()
            .any(|value| !value.is_finite())
        {
            return Err(MaterialError::Parameters {
                what: "orthotropic stiffness is not representable".to_string(),
            });
        }
        Ok(law)
    }

    // S = D^-1 R D^-1, D = diag(sqrt(E_i)), R = L L^T with unit diagonal.
    // Per-axis normalization retains orthotropy without forming E_i*E_j or
    // dimensioned compliance determinants. No eigenvalue floor is introduced.
    fn normal_cholesky(&self) -> ([f64; 3], [[f64; 3]; 3]) {
        let roots = self.e.map(fs_math::det::sqrt);
        let coupling = |nu: f64, i: usize, j: usize| {
            let ratio = roots[j] / roots[i];
            if ratio.is_finite() {
                -nu * ratio
            } else {
                // A representable coupling can combine a tiny Poisson ratio
                // with an overflowing modulus-root ratio.
                -(nu * roots[j]) / roots[i]
            }
        };
        let a = coupling(self.nu[0], 0, 1);
        let b = coupling(self.nu[1], 0, 2);
        let c = coupling(self.nu[2], 1, 2);
        let minor2 = (1.0 - a) * (1.0 + a);
        let remainder = (-a).mul_add(b, c);
        let determinant = minor2.mul_add((1.0 - b) * (1.0 + b), -remainder * remainder);
        let l11 = fs_math::det::sqrt(minor2);
        (
            roots,
            [
                [1.0, 0.0, 0.0],
                [a, l11, 0.0],
                [b, remainder / l11, fs_math::det::sqrt(determinant / minor2)],
            ],
        )
    }

    /// The constant stiffness: inverse of the normal-block compliance +
    /// diagonal shear.
    #[must_use]
    pub fn stiffness(&self) -> Tangent6 {
        let (roots, lower) = self.normal_cholesky();
        let mut inverse_lower = [[0.0; 3]; 3];
        for column in 0..3 {
            for row in column..3 {
                let mut value = if row == column { 1.0 } else { 0.0 };
                for k in column..row {
                    value = (-lower[row][k]).mul_add(inverse_lower[k][column], value);
                }
                inverse_lower[row][column] = value / lower[row][row];
            }
        }
        let mut c6 = [[0.0f64; 6]; 6];
        for i in 0..3 {
            for j in i..3 {
                let inverse = (j..3).fold(0.0, |sum, k| {
                    inverse_lower[k][i].mul_add(inverse_lower[k][j], sum)
                });
                let modulus = if i == j {
                    self.e[i]
                } else {
                    roots[i] * roots[j]
                };
                let stiffness = modulus * inverse;
                c6[i][j] = stiffness;
                c6[j][i] = stiffness;
            }
        }
        let [g12, g23, g31] = self.g;
        c6[3][3] = 2.0 * g12;
        c6[4][4] = 2.0 * g23;
        c6[5][5] = 2.0 * g31;
        c6
    }
}

impl SmallStrainLaw for OrthotropicElastic {
    type State = ();

    fn initial_state(&self) {}

    fn stress(&self, strain: &Voigt, (): &Self::State) -> Voigt {
        apply(&self.stiffness(), strain)
    }

    fn tangent(&self, _strain: &Voigt, (): &Self::State) -> Tangent6 {
        self.stiffness()
    }

    fn update_state(&self, _strain: &Voigt, (): &Self::State) {}

    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            has_stored_energy: true,
            dissipation_nonnegative: true,
            polyconvex: Some(true),
            tangent_symmetric: true,
            failure_envelope: "none declared (linear law; validity ends at strain_limit)",
        }
    }

    fn card(&self) -> ModelCard {
        ModelCard::new(
            "material.orthotropic-elastic",
            "0.1.0",
            Ambition::Solid,
            vec![
                "linear kinematics".to_string(),
                "material axes aligned with coordinates".to_string(),
            ],
            ValidityDomain::unconstrained().with("strain-magnitude", 0.0, self.strain_limit),
            vec!["misaligned-axis use requires a rotation the caller owns".to_string()],
            0.02,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64, tolerance: f64) {
        assert!(actual.is_finite() && expected.is_finite());
        assert!(
            (actual - expected).abs() <= tolerance * actual.abs().max(expected.abs()),
            "{actual:e} != {expected:e}"
        );
    }

    #[test]
    fn g0_orthotropic_stress_inverts_engineering_compliance() {
        let law = OrthotropicElastic::new(
            [140e9, 10e9, 7e9],
            [0.31, 0.24, 0.35],
            [5e9, 3e9, 4e9],
            0.02,
        )
        .unwrap();
        // Independent engineering-compliance response to all six stresses.
        let stress = [2e7, -1e7, 3e6, 4e6, -5e6, 6e6];
        let [e1, e2, e3] = law.e;
        let [nu12, nu13, nu23] = law.nu;
        let strain = [
            (stress[0] - nu12 * stress[1] - nu13 * stress[2]) / e1,
            -nu12 * stress[0] / e1 + (stress[1] - nu23 * stress[2]) / e2,
            -nu13 * stress[0] / e1 - nu23 * stress[1] / e2 + stress[2] / e3,
            stress[3] / (2.0 * law.g[0]),
            stress[4] / (2.0 * law.g[1]),
            stress[5] / (2.0 * law.g[2]),
        ];
        for (actual, expected) in law.stress(&strain, &()).into_iter().zip(stress) {
            close(actual, expected, 8e-15);
        }
        let stiffness = law.stiffness();
        for i in 0..6 {
            for j in 0..6 {
                assert_eq!(stiffness[i][j].to_bits(), stiffness[j][i].to_bits());
            }
        }
        for poisson in [-0.8, 0.0, 0.3, 0.4999] {
            let isotropic = IsotropicElastic::new(1.3e9, poisson, 0.02).unwrap();
            let shear = isotropic.lame().1;
            let orthotropic =
                OrthotropicElastic::new([1.3e9; 3], [poisson; 3], [shear; 3], 0.02).unwrap();
            for (actual, expected) in orthotropic
                .stiffness()
                .iter()
                .flatten()
                .zip(isotropic.stiffness().iter().flatten())
            {
                close(*actual, *expected, 3e-12);
            }
        }
    }

    #[test]
    fn g3_orthotropic_stiffness_preserves_the_stress_scale() {
        let e = [140.0, 10.0, 7.0];
        let g = [5.0, 3.0, 4.0];
        let nu = [0.31, 0.24, 0.35];
        let reference = OrthotropicElastic::new(e, nu, g, 0.02).unwrap();
        let strain = [0.001, -0.002, 0.0003, 0.004, -0.001, 0.002];
        for scale in [1e-240, 1e-120, 1e-20, 1.0, 1e20, 1e120, 1e240] {
            let law = OrthotropicElastic::new(e.map(|v| scale * v), nu, g.map(|v| scale * v), 0.02)
                .unwrap();
            for (actual, expected) in law
                .stiffness()
                .iter()
                .flatten()
                .zip(reference.stiffness().iter().flatten())
            {
                close(actual / scale, *expected, 2e-14);
            }
            for (actual, expected) in law
                .stress(&strain, &())
                .into_iter()
                .zip(reference.stress(&strain, &()))
            {
                close(actual / scale, expected, 2e-14);
            }
        }
    }

    #[test]
    fn g3_orthotropic_axis_normalization_retains_large_modulus_contrast() {
        // This dimensionless SPD compliance is independent of the three
        // modulus scales. Its inverse is an analytic reference, not a solve
        // through the implementation being tested.
        let e: [f64; 3] = [1e-240, 1e240, 1.0];
        let coupling = [[0.0, 0.2, 0.1], [0.2, 0.0, 0.3], [0.1, 0.3, 0.0]];
        let cofactors = [[0.91, 0.23, 0.16], [0.23, 0.99, 0.32], [0.16, 0.32, 0.96]];
        for axes in [
            [0, 1, 2],
            [0, 2, 1],
            [1, 0, 2],
            [1, 2, 0],
            [2, 0, 1],
            [2, 1, 0],
        ] {
            let moduli = axes.map(|i| e[i]);
            let poisson = |i: usize, j: usize| {
                coupling[axes[i]][axes[j]] * moduli[i].sqrt() / moduli[j].sqrt()
            };
            let law = OrthotropicElastic::new(
                moduli,
                [poisson(0, 1), poisson(0, 2), poisson(1, 2)],
                [1.0; 3],
                0.02,
            )
            .unwrap();
            let stiffness = law.stiffness();
            for i in 0..3 {
                for j in 0..3 {
                    let scale = moduli[i].sqrt() * moduli[j].sqrt();
                    close(
                        stiffness[i][j] / scale,
                        cofactors[axes[i]][axes[j]] / 0.848,
                        2e-14,
                    );
                }
            }
        }
    }

    #[test]
    fn g0_orthotropic_admission_does_not_depend_on_stress_units() {
        for scale in [1e-240, 1.0, 1e240] {
            for nu in [
                [0.5; 3],
                [0.9; 3],
                [-1.0; 3],
                [1.0, 0.0, 0.0],
                [f64::NAN; 3],
            ] {
                assert!(OrthotropicElastic::new([scale; 3], nu, [scale; 3], 0.02).is_err());
            }
        }
        assert!(OrthotropicElastic::new([1.0; 3], [0.0; 3], [f64::MAX; 3], 0.02).is_err());
    }
}
