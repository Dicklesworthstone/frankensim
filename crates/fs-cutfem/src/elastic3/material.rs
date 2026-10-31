//! Prepared reference tensors shared by CutFEM bulk and boundary operators.
//!
//! The native material owner admits the complete Mandel tensor, density and
//! identity. This module only caches its engineering-shear representation and
//! a conservative stiffness scale. No material is resolved or fitted here.

use std::sync::Arc;

use fs_material::IsotropicElastic;
use fs_material::elastic::acoustic::ElasticAcousticState;

use super::ElasticityError3;
use crate::HexCell;

/// Complete reference materials for active background cells.
///
/// All tensors use the simulation's global frame. Density and material identity
/// remain attached to the original immutable state. Design stiffness scales do
/// not modify this reference state or its density. Each returned state applies
/// throughout its cell; material interfaces must therefore follow the supplied
/// cell assignment. The material's explicit strain ceiling is retained, but is
/// not applied to algebraic shape-function strains during element assembly.
#[derive(Clone, Copy)]
pub enum ElasticMaterialField3<'a> {
    /// One state shared by all active cells; derived coefficients are prepared
    /// once and then shared without copying matrices at each cell or quadrature.
    Uniform(&'a ElasticAcousticState),
    /// Select a complete state from the background cell's physical bounds.
    /// The callback must be pure; an error aborts construction. It receives
    /// active cells, rather than a precomputed array in a guessed active order.
    ByCell(&'a dyn Fn(HexCell) -> Result<ElasticAcousticState, ElasticityError3>),
}

/// A cache of already admitted material data, shared cheaply for uniform fields.
#[derive(Clone, PartialEq)]
pub(super) struct PreparedMaterial3 {
    state: ElasticAcousticState,
    engineering_stiffness: [[f64; 6]; 6],
    stiffness_bound: f64,
}

/// Keep the legacy operation tree while general material paths use cached C.
#[derive(Clone, PartialEq)]
pub(super) enum ReferenceMaterial3 {
    Legacy([f64; 2]),
    Source(Arc<PreparedMaterial3>),
}

impl ReferenceMaterial3 {
    pub(super) fn legacy(material: &IsotropicElastic) -> Result<Self, ElasticityError3> {
        let card = IsotropicElastic::new(material.youngs, material.poisson, material.strain_limit)
            .map_err(|_| ElasticityError3::Invalid("invalid isotropic material"))?;
        let (lambda, mu) = card.lame();
        if !lambda.is_finite()
            || !mu.is_finite()
            || mu <= 0.0
            || !card.strain_limit.is_finite()
            || card.strain_limit <= 0.0
            || !((lambda + 2.0 * mu) / mu <= 4.0)
        {
            return Err(ElasticityError3::Invalid(
                "material outside finite compressible regime",
            ));
        }
        Ok(Self::Legacy([lambda, mu]))
    }

    pub(super) fn from_state(state: ElasticAcousticState) -> Result<Self, ElasticityError3> {
        // For engineering strain and physical stress the work-conjugate map
        // is C_eng = D C_M D. Scaling only columns would double shear stresses.
        let factors = [
            1.0,
            1.0,
            1.0,
            core::f64::consts::FRAC_1_SQRT_2,
            core::f64::consts::FRAC_1_SQRT_2,
            core::f64::consts::FRAC_1_SQRT_2,
        ];
        let mandel = state.stiffness_mandel_pa();
        let engineering_stiffness: [[f64; 6]; 6] = std::array::from_fn(|i| {
            std::array::from_fn(|j| mandel[i][j] * factors[i] * factors[j])
        });
        let mut stiffness_bound = 0.0_f64;
        for row in mandel {
            let row_sum: f64 = row.iter().map(|entry| entry.abs()).sum();
            if !row_sum.is_finite() {
                return Err(ElasticityError3::Invalid(
                    "material stiffness bound overflow",
                ));
            }
            stiffness_bound = stiffness_bound.max(row_sum);
        }
        if stiffness_bound <= 0.0
            || engineering_stiffness
                .iter()
                .flatten()
                .any(|v| !v.is_finite())
            || (0..6).any(|i| engineering_stiffness[i][i] <= 0.0)
        {
            return Err(ElasticityError3::Invalid(
                "material engineering stiffness is unrepresentable",
            ));
        }
        Ok(Self::Source(Arc::new(PreparedMaterial3 {
            state,
            engineering_stiffness,
            stiffness_bound,
        })))
    }

    /// Physical stress [xx, yy, zz, xy, yz, zx] of the basis N e_component.
    /// This is a linear operator action, not admission of a physical strain.
    pub(super) fn basis_stress(&self, g: [f64; 3], component: usize) -> [f64; 6] {
        match self {
            Self::Legacy([lambda, mu]) => {
                let mut stress = [lambda * g[component]; 6];
                stress[component] += 2.0 * mu * g[component];
                stress[3] = mu
                    * match component {
                        0 => g[1],
                        1 => g[0],
                        _ => 0.0,
                    };
                stress[4] = mu
                    * match component {
                        1 => g[2],
                        2 => g[1],
                        _ => 0.0,
                    };
                stress[5] = mu
                    * match component {
                        0 => g[2],
                        2 => g[0],
                        _ => 0.0,
                    };
                stress
            }
            Self::Source(material) => {
                // Each displacement basis has only three nonzero engineering
                // strains. Avoid a full 6x6 matvec in stress/traction hot loops.
                let (columns, strain) = match component {
                    0 => ([0, 3, 5], [g[0], g[1], g[2]]),
                    1 => ([1, 3, 4], [g[1], g[0], g[2]]),
                    2 => ([2, 4, 5], [g[2], g[1], g[0]]),
                    _ => panic!("displacement component must lie in 0..3"),
                };
                std::array::from_fn(|i| {
                    let row = &material.engineering_stiffness[i];
                    row[columns[0]] * strain[0]
                        + row[columns[1]] * strain[1]
                        + row[columns[2]] * strain[2]
                })
            }
        }
    }

    pub(super) fn traction(&self, g: [f64; 3], component: usize, n: [f64; 3]) -> [f64; 3] {
        match self {
            Self::Legacy([lambda, mu]) => {
                // Preserve the original Lamé traction operation order.
                let dn: f64 = g.iter().zip(n).map(|(g, n)| g * n).sum();
                std::array::from_fn(|d| {
                    lambda * g[component] * n[d]
                        + mu * (n[component] * g[d] + if component == d { dn } else { 0.0 })
                })
            }
            Self::Source(_) => {
                let s = self.basis_stress(g, component);
                [
                    s[0] * n[0] + s[3] * n[1] + s[5] * n[2],
                    s[3] * n[0] + s[1] * n[1] + s[4] * n[2],
                    s[5] * n[0] + s[4] * n[1] + s[2] * n[2],
                ]
            }
        }
    }

    pub(super) fn nitsche_modulus(&self) -> f64 {
        match self {
            Self::Legacy([lambda, mu]) => lambda + 2.0 * mu,
            Self::Source(material) => material.stiffness_bound,
        }
    }

    pub(super) fn ghost_modulus(&self) -> f64 {
        match self {
            Self::Legacy([_, mu]) => *mu,
            Self::Source(material) => material.stiffness_bound,
        }
    }

    pub(super) fn source_state(&self) -> Option<&ElasticAcousticState> {
        match self {
            Self::Legacy(_) => None,
            Self::Source(material) => Some(&material.state),
        }
    }

    pub(super) fn legacy_lame(&self) -> Option<[f64; 2]> {
        match self {
            Self::Legacy(lame) => Some(*lame),
            Self::Source(_) => None,
        }
    }
}

/// Unpublished construction input; uniform preparation happens only once.
pub(super) enum MaterialProvider3<'a> {
    Legacy(ReferenceMaterial3),
    Uniform(ReferenceMaterial3),
    ByCell(&'a dyn Fn(HexCell) -> Result<ElasticAcousticState, ElasticityError3>),
}

impl<'a> MaterialProvider3<'a> {
    pub(super) fn legacy(material: &IsotropicElastic) -> Result<Self, ElasticityError3> {
        ReferenceMaterial3::legacy(material).map(Self::Legacy)
    }

    pub(super) fn from_field(field: ElasticMaterialField3<'a>) -> Result<Self, ElasticityError3> {
        match field {
            ElasticMaterialField3::Uniform(state) => {
                ReferenceMaterial3::from_state(state.clone()).map(Self::Uniform)
            }
            ElasticMaterialField3::ByCell(callback) => Ok(Self::ByCell(callback)),
        }
    }

    pub(super) fn at(&self, bounds: HexCell) -> Result<ReferenceMaterial3, ElasticityError3> {
        match self {
            Self::Legacy(material) | Self::Uniform(material) => Ok(material.clone()),
            Self::ByCell(callback) => ReferenceMaterial3::from_state(callback(bounds)?),
        }
    }
}

/// Engineering strain [xx, yy, zz, 2xy, 2yz, 2zx] of N e_component.
pub(super) fn engineering_strain(g: [f64; 3], component: usize) -> [f64; 6] {
    match component {
        0 => [g[0], 0.0, 0.0, g[1], 0.0, g[2]],
        1 => [0.0, g[1], 0.0, g[0], g[2], 0.0],
        2 => [0.0, 0.0, g[2], 0.0, g[1], g[0]],
        _ => panic!("displacement component must lie in 0..3"),
    }
}
