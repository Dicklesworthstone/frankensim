//! Receipted reference-ratio ingress, with no inferred missing constants.
//! L/R/T are material axes 1/2/3; the source's major Poisson set governs.
//! A plane result is deliberately not a fabricated three-dimensional law.

use super::*;
use crate::plane_stress_viscosity::MaterialPlane;

/// Longitudinal reference modulus [Pa], at the source's stated condition.
pub const REFERENCE_YOUNG_MODULUS_PROPERTY: &str = "young_modulus_longitudinal";
/// E_R/E_L and E_T/E_L, in material-axis order 2,3.
pub const REFERENCE_YOUNG_RATIO_PROPERTIES: [&str; 2] = ["er_over_el", "et_over_el"];
/// Shear ratios in canonical G12,G23,G31 order.
pub const REFERENCE_SHEAR_RATIO_PROPERTIES: [&str; 3] =
    ["glr_over_el", "grt_over_el", "glt_over_el"];
/// Source major Poisson ratios in nu12,nu13,nu23 order.
pub const REFERENCE_MAJOR_POISSON_PROPERTIES: [&str; 3] = ["nu_lr", "nu_lt", "nu_rt"];
/// Independently measured reciprocal ratios in nu21,nu31,nu32 order.
pub const REFERENCE_MINOR_POISSON_PROPERTIES: [&str; 3] = ["nu_rl", "nu_tl", "nu_tr"];
/// Identity domain binding the source bundle, requested law and derivation rule.
pub const REFERENCE_RATIO_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-material.orthotropic-reference-ratio.v1";

/// Source wood directions; these labels do not infer a spatial orientation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WoodAxis {
    /// Longitudinal, along grain.
    L,
    /// Radial, across growth rings.
    R,
    /// Tangential, tangent to growth rings.
    T,
}
/// L/R/T map to material axes 1/2/3, respectively.
pub const REFERENCE_RATIO_AXES: [WoodAxis; 3] = [WoodAxis::L, WoodAxis::R, WoodAxis::T];

/// Consumer-selected dimensional support; no unrequested modulus is invented.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OrthotropicRequirementSet {
    /// Density and all nine independent three-dimensional constants.
    Full3d,
    /// Density and the four constants governing this ordered material plane.
    Plane(MaterialPlane),
}
/// Explicit rule for the source's independently measured Poisson tables.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PoissonReciprocityRule {
    /// Use nu_LR,nu_LT,nu_RT; expose minor-table discrepancies without averaging.
    MajorRatios,
}
/// Derivation metadata; derived constants are not fabricated material claims.
#[derive(Clone, Debug, PartialEq)]
pub struct ReferenceRatioDerivation {
    /// Source directions assigned to material axes 1/2/3.
    pub axes: [WoodAxis; 3],
    /// Caller-selected reciprocity interpretation.
    pub rule: PoissonReciprocityRule,
    /// Signed relative minor/reciprocal-major minus one, in RL,TL,TR order.
    /// None means the minor is absent or its pair is not requested by the set.
    /// A zero/zero pair has zero residual; unrepresentable residuals refuse.
    pub reciprocity_residuals: [Option<f64>; 3],
    /// Binds the exact source receipts, set, strain limit, axes, rule and residuals.
    pub identity: ContentHash,
}

/// Elastic law confined to one material plane, in engineering shear coordinates.
#[derive(Clone, Debug, PartialEq)]
pub struct OrthotropicPlaneElastic {
    /// Ordered plane: e[0] and e[1] refer to its first and second material axes.
    pub plane: MaterialPlane,
    /// In-plane Young's moduli [Pa].
    pub e: [f64; 2],
    /// Minus the second/first-axis strain ratio under uniaxial first-axis stress.
    pub nu: f64,
    /// In-plane engineering shear modulus [Pa].
    pub g: f64,
    /// Caller-declared calibrated strain magnitude bound.
    pub strain_limit: f64,
}
impl OrthotropicPlaneElastic {
    /// Reduced Q in [sigma_x,sigma_y,tau_xy] = Q[epsilon_x,epsilon_y,gamma_xy].
    #[must_use]
    pub fn stiffness(&self) -> [f64; 9] {
        let denominator = 1.0 - self.nu * (self.nu * self.e[1] / self.e[0]);
        let cross = self.nu * self.e[1] / denominator;
        [
            self.e[0] / denominator,
            cross,
            0.0,
            cross,
            self.e[1] / denominator,
            0.0,
            0.0,
            0.0,
            self.g,
        ]
    }
}
/// A plane state retains its source bundle without asserting out-of-plane data.
#[derive(Clone, Debug, PartialEq)]
pub struct OrthotropicPlaneStatePoint {
    resolved: ResolvedMaterialStatePoint,
    density_kg_m3: f64,
    law: OrthotropicPlaneElastic,
    derivation: ReferenceRatioDerivation,
}
impl OrthotropicPlaneStatePoint {
    /// Complete evidence-bearing bundle for this plane's requested constants.
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedMaterialStatePoint {
        &self.resolved
    }
    /// Density at the queried state [kg/m3].
    #[must_use]
    pub const fn density_kg_m3(&self) -> f64 {
        self.density_kg_m3
    }
    /// Admitted in-plane law.
    #[must_use]
    pub const fn law(&self) -> &OrthotropicPlaneElastic {
        &self.law
    }
    /// Explicit ratio and Poisson-rule provenance for the admitted plane.
    #[must_use]
    pub const fn derivation(&self) -> &ReferenceRatioDerivation {
        &self.derivation
    }
}
/// Full and plane requirements produce distinct laws, never partial 3-D tensors.
#[derive(Clone, Debug, PartialEq)]
pub enum ReferenceRatioOrthotropicStatePoint {
    /// Complete principal-axis three-dimensional law.
    Full3d(OrthotropicElasticStatePoint),
    /// Only the explicitly requested plane is characterized.
    Plane(OrthotropicPlaneStatePoint),
}
impl ReferenceRatioOrthotropicStatePoint {
    /// Original claims and all selected usage receipts, including queried minors.
    #[must_use]
    pub fn resolved(&self) -> &ResolvedMaterialStatePoint {
        match self {
            Self::Full3d(s) => &s.resolved,
            Self::Plane(s) => &s.resolved,
        }
    }
    /// Density at the common query point [kg/m3].
    #[must_use]
    pub fn density_kg_m3(&self) -> f64 {
        match self {
            Self::Full3d(s) => s.density_kg_m3,
            Self::Plane(s) => s.density_kg_m3,
        }
    }
    /// Explicit derivation; canonical states manually wrapped in Full3d retain None.
    #[must_use]
    pub fn derivation(&self) -> Option<&ReferenceRatioDerivation> {
        match self {
            Self::Full3d(s) => s.derivation(),
            Self::Plane(s) => Some(&s.derivation),
        }
    }
    /// Three-dimensional law, only when the complete set was resolved.
    #[must_use]
    pub fn law(&self) -> Option<&OrthotropicElastic> {
        match self {
            Self::Full3d(s) => Some(s.law()),
            Self::Plane(_) => None,
        }
    }
}

/// Resolve source ratios atomically at one supported state, retaining receipts.
/// Missing properties refuse by name. Independently measured minor Poisson
/// ratios are resolved when present for a requested pair, recorded and never
/// averaged into the declared major set. No frequency or moisture correction
/// is introduced. A Full3d request cannot fall back to a plane result.
pub fn resolve_orthotropic_elastic_state_point_from_reference_ratios(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    strain_limit: f64,
    set: OrthotropicRequirementSet,
    rule: PoissonReciprocityRule,
) -> Result<ReferenceRatioOrthotropicStatePoint, MaterialStatePointError> {
    if !strain_limit.is_finite() || strain_limit <= 0.0 {
        return Err(MaterialStatePointError::InvalidDerived {
            quantity: "orthotropic_linear_strain_limit",
        });
    }
    let (young_needed, pair_needed, plane_axes, code) = match set {
        OrthotropicRequirementSet::Full3d => ([true; 3], [true; 3], None, 0_u8),
        OrthotropicRequirementSet::Plane(p) => match p {
            MaterialPlane::Axes12 => (
                [true, true, false],
                [true, false, false],
                Some((p, 0, 1, 0)),
                1,
            ),
            MaterialPlane::Axes21 => (
                [true, true, false],
                [true, false, false],
                Some((p, 1, 0, 0)),
                2,
            ),
            MaterialPlane::Axes13 => (
                [true, false, true],
                [false, true, false],
                Some((p, 0, 2, 1)),
                3,
            ),
            MaterialPlane::Axes23 => (
                [false, true, true],
                [false, false, true],
                Some((p, 1, 2, 2)),
                4,
            ),
        },
    };
    let mut requirements = vec![ScalarPropertyRequirement::try_new(
        REFERENCE_YOUNG_MODULUS_PROPERTY,
        Pressure::DIMS,
        ScalarAdmissibility::StrictlyPositive,
    )?];
    for (index, name) in REFERENCE_YOUNG_RATIO_PROPERTIES.iter().enumerate() {
        if young_needed[index + 1] {
            requirements.push(ScalarPropertyRequirement::try_new(
                *name,
                Dims::NONE,
                ScalarAdmissibility::StrictlyPositive,
            )?);
        }
    }
    // Pair order nu12,nu13,nu23 differs from shear order G12,G23,G31.
    for (pair, shear) in [0, 2, 1].into_iter().enumerate() {
        if pair_needed[pair] {
            requirements.push(ScalarPropertyRequirement::try_new(
                REFERENCE_SHEAR_RATIO_PROPERTIES[shear],
                Dims::NONE,
                ScalarAdmissibility::StrictlyPositive,
            )?);
            requirements.push(ScalarPropertyRequirement::try_new(
                REFERENCE_MAJOR_POISSON_PROPERTIES[pair],
                Dims::NONE,
                ScalarAdmissibility::Finite,
            )?);
        }
    }
    requirements.push(ScalarPropertyRequirement::try_new(
        DENSITY_PROPERTY,
        Density::DIMS,
        ScalarAdmissibility::StrictlyPositive,
    )?);
    // Report absent elastic ratios before density, independent of the generic
    // resolver's alphabetical query order (e.g. both rosewood packs lack both).
    for requirement in &requirements {
        if card.claims().claims_for(requirement.name()).is_empty() {
            return Err(MaterialStatePointError::Query {
                property: requirement.name().into(),
                source: MatDbError::UnknownProperty {
                    property: requirement.name().into(),
                },
            });
        }
    }
    for (pair, name) in REFERENCE_MINOR_POISSON_PROPERTIES.iter().enumerate() {
        if pair_needed[pair] && !card.claims().claims_for(name).is_empty() {
            requirements.push(ScalarPropertyRequirement::try_new(
                *name,
                Dims::NONE,
                ScalarAdmissibility::Finite,
            )?);
        }
    }
    let resolved = resolve_material_state_point(card, point, &requirements, selection)?;
    let value = |name| {
        resolved
            .property(name)
            .map(ResolvedScalarProperty::value_si)
    };
    let el = value(REFERENCE_YOUNG_MODULUS_PROPERTY).expect("required reference modulus");
    let mut e = [None; 3];
    for (index, needed) in young_needed.into_iter().enumerate() {
        if needed {
            let modulus = if index == 0 {
                el
            } else {
                el * value(REFERENCE_YOUNG_RATIO_PROPERTIES[index - 1])
                    .expect("required Young ratio")
            };
            if !modulus.is_finite() || modulus <= 0.0 {
                return Err(MaterialStatePointError::InvalidDerived {
                    quantity: "orthotropic_modulus",
                });
            }
            e[index] = Some(modulus);
        }
    }
    let mut g = [None; 3];
    let nu = REFERENCE_MAJOR_POISSON_PROPERTIES.map(value);
    for (pair, shear) in [0, 2, 1].into_iter().enumerate() {
        if pair_needed[pair] {
            let modulus =
                el * value(REFERENCE_SHEAR_RATIO_PROPERTIES[shear]).expect("required shear ratio");
            if !modulus.is_finite() || modulus <= 0.0 {
                return Err(MaterialStatePointError::InvalidDerived {
                    quantity: "orthotropic_modulus",
                });
            }
            g[shear] = Some(modulus);
        }
    }
    let mut residuals = [None; 3];
    for (pair, (i, j)) in [(0, 1), (0, 2), (1, 2)].into_iter().enumerate() {
        if let Some(minor) = value(REFERENCE_MINOR_POISSON_PROPERTIES[pair]) {
            let expected = nu[pair].expect("requested major") * e[j].expect("requested modulus")
                / e[i].expect("requested modulus");
            let residual = if expected == 0.0 && minor == 0.0 {
                0.0
            } else {
                minor / expected - 1.0
            };
            if !residual.is_finite() {
                return Err(MaterialStatePointError::InvalidDerived {
                    quantity: "orthotropic_reciprocity_residual",
                });
            }
            residuals[pair] = Some(residual);
        }
    }
    let mut identity = DomainHasher::new(REFERENCE_RATIO_IDENTITY_DOMAIN);
    identity.update(resolved.identity().as_bytes());
    identity.update(&[code, 0, 1, 2, 0]); // set, L/R/T axes, MajorRatios rule
    identity.update(&strain_limit.to_bits().to_le_bytes());
    for r in residuals {
        identity.update(&[u8::from(r.is_some())]);
        if let Some(v) = r {
            identity.update(&v.to_bits().to_le_bytes());
        }
    }
    let derivation = ReferenceRatioDerivation {
        axes: REFERENCE_RATIO_AXES,
        rule,
        reciprocity_residuals: residuals,
        identity: identity.finalize(),
    };
    let density_kg_m3 = value(DENSITY_PROPERTY).expect("required density");
    let invalid = |_| MaterialStatePointError::InvalidDerived {
        quantity: "orthotropic_compliance",
    };
    if let Some((plane, i, j, pair)) = plane_axes {
        let e = [
            e[i].expect("requested modulus"),
            e[j].expect("requested modulus"),
        ];
        let mut major = nu[pair].expect("requested major");
        if plane == MaterialPlane::Axes21 {
            major = major * e[0] / e[1];
        }
        let law = OrthotropicPlaneElastic {
            plane,
            e,
            nu: major,
            g: g[[0, 2, 1][pair]].expect("requested shear"),
            strain_limit,
        };
        let denominator = 1.0 - law.nu * (law.nu * law.e[1] / law.e[0]);
        if !denominator.is_finite()
            || denominator <= 0.0
            || law.stiffness().iter().any(|v| !v.is_finite())
        {
            return Err(invalid(()));
        }
        Ok(ReferenceRatioOrthotropicStatePoint::Plane(
            OrthotropicPlaneStatePoint {
                resolved,
                density_kg_m3,
                law,
                derivation,
            },
        ))
    } else {
        let law = OrthotropicElastic::new(
            e.map(|x| x.expect("full modulus")),
            nu.map(|x| x.expect("full major")),
            g.map(|x| x.expect("full shear")),
            strain_limit,
        )
        .map_err(|_| invalid(()))?;
        if law.stiffness().iter().flatten().any(|v| !v.is_finite()) {
            return Err(invalid(()));
        }
        Ok(ReferenceRatioOrthotropicStatePoint::Full3d(
            OrthotropicElasticStatePoint {
                resolved,
                density_kg_m3,
                law,
                derivation: Some(derivation),
            },
        ))
    }
}
