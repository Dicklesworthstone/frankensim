//! Linear elasticity: isotropic and orthotropic small-strain laws. The
//! tangent IS the stiffness (state-free), so consistency is exact by
//! construction — these are the baseline the FD gate calibrates itself
//! against.

use crate::{MaterialAdmissibility, MaterialError, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::{Ambition, ModelCard, ValidityDomain};

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
        // Tensor-Voigt contraction: shear strain columns count twice.
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
    /// Construct with positivity + compliance-definiteness spot checks.
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
        // Thermodynamic requirement: the 3×3 normal-block compliance must
        // be positive definite (leading principal minors > 0).
        let s = law.normal_compliance();
        let m1 = s[0][0];
        let m2 = s[0][0] * s[1][1] - s[0][1] * s[1][0];
        let m3 = s[0][0] * (s[1][1] * s[2][2] - s[1][2] * s[2][1])
            - s[0][1] * (s[1][0] * s[2][2] - s[1][2] * s[2][0])
            + s[0][2] * (s[1][0] * s[2][1] - s[1][1] * s[2][0]);
        if !(m1 > 0.0 && m2 > 0.0 && m3 > 0.0) {
            return Err(MaterialError::Parameters {
                what: "orthotropic compliance is not positive definite \
                       (Poisson ratios thermodynamically inadmissible)"
                    .to_string(),
            });
        }
        Ok(law)
    }

    fn normal_compliance(&self) -> [[f64; 3]; 3] {
        let [e1, e2, e3] = self.e;
        let [nu12, nu13, nu23] = self.nu;
        // Symmetry: ν21/E2 = ν12/E1 etc.
        [
            [1.0 / e1, -nu12 / e1, -nu13 / e1],
            [-nu12 / e1, 1.0 / e2, -nu23 / e2],
            [-nu13 / e1, -nu23 / e2, 1.0 / e3],
        ]
    }

    /// The constant stiffness: inverse of the normal-block compliance +
    /// diagonal shear.
    #[must_use]
    pub fn stiffness(&self) -> Tangent6 {
        let s = self.normal_compliance();
        // Invert the 3×3 by cofactors.
        let det = s[0][0] * (s[1][1] * s[2][2] - s[1][2] * s[2][1])
            - s[0][1] * (s[1][0] * s[2][2] - s[1][2] * s[2][0])
            + s[0][2] * (s[1][0] * s[2][1] - s[1][1] * s[2][0]);
        let inv = |r: usize, c: usize| -> f64 {
            let (r1, r2) = ((r + 1) % 3, (r + 2) % 3);
            let (c1, c2) = ((c + 1) % 3, (c + 2) % 3);
            (s[c1][r1] * s[c2][r2] - s[c1][r2] * s[c2][r1]) / det
        };
        let mut c6 = [[0.0f64; 6]; 6];
        for (i, row) in c6.iter_mut().enumerate().take(3) {
            for (j, slot) in row.iter_mut().enumerate().take(3) {
                *slot = inv(i, j);
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
