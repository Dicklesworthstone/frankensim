//! Uniform circular specimens lowered to the existing prestressed-string solver.
//!
//! This is the homogeneous isotropic Euler–Bernoulli/Kirchhoff–Carrier rung:
//! `A = pi r²`, `I = pi r⁴/4`, `mu = rho A`, `EA = E A`, `EI = E I`.
//! Length and radius describe the specimen at the resolved material state.
//! Tension, support motion, modal budget, and declared losses remain mechanical
//! inputs. Rebinding holds geometry and tension fixed; it is not a thermal
//! evolution, fixed-mass comparison, imposed-extension solve, or pitch tuner.

use fs_blake3::{ContentHash, DomainHasher};
use fs_material::state_point::{
    DENSITY_PROPERTY, ResolvedMaterialStatePoint, YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Density, Dims, Pressure};
use fs_scenario::PrestressedString;

use crate::acoustic_realize::AcousticRealizeError;

/// One material-bound specimen plus its independently authored mechanical state.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedStringSpecimen {
    string: PrestressedString,
    material: ResolvedMaterialStatePoint,
    area_m2: f64,
    second_moment_m4: f64,
    mass_kg: f64,
    specimen_identity: ContentHash,
}

impl ResolvedStringSpecimen {
    /// Ready for `AcousticAssembly::string` and the existing string realizer.
    #[must_use]
    pub const fn string(&self) -> PrestressedString {
        self.string
    }

    /// Original resolved bundle, including exact card/state identities and receipts.
    #[must_use]
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }

    /// Uniform cross-sectional area [m²].
    #[must_use]
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }

    /// Area second moment about either transverse centroidal axis [m⁴].
    #[must_use]
    pub const fn second_moment_m4(&self) -> f64 {
        self.second_moment_m4
    }

    /// Mass over the complete declared length [kg].
    #[must_use]
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }

    /// Identity of geometry and material resolution, excluding loading and solver
    /// options. It is not an identity of the complete acoustic scenario.
    #[must_use]
    pub const fn specimen_identity(&self) -> ContentHash {
        self.specimen_identity
    }
}

/// Derive a uniform circular string's mass, stiffness, and observer width together.
///
/// `state` must resolve positive `density` [kg/m³] and `young_modulus` [Pa]
/// using `fs_material::state_point::resolve_material_state_point`. The same
/// bundle API drives disc elasticity; this string rung needs no Poisson ratio,
/// yield stress, or unrelated optical/thermal data. Extra resolved properties
/// and all original evidence remain available through the returned specimen.
///
/// The template supplies length, tension, support/polarization choices, modal
/// budget, and declared loss parameters. Its `lin_density_kg_m`,
/// `axial_stiffness_n`, `bending_stiffness_n_m2`, and `width_m` are replaced by
/// the derived values, never merged into material-card authority. Losses are
/// not inferred from material names; the existing realizer's loss and radiation
/// approximations still apply. The diameter supplies its compact strip-radiator
/// width; this binding does not establish circular-wire radiation accuracy.
/// No plasticity, anisotropy, or melting is claimed.
///
/// # Errors
/// Missing or dimensionally wrong required properties, nonpositive/nonfinite
/// geometry/material values, invalid mechanical inputs, or unrepresentable
/// derived quantities refuse before publishing a specimen.
pub fn with_uniform_circular_material_state(
    mut string: PrestressedString,
    radius_m: f64,
    state: &ResolvedMaterialStatePoint,
) -> Result<ResolvedStringSpecimen, AcousticRealizeError> {
    let refuse = |what| AcousticRealizeError::InvalidDescription { what };
    let positive = |v: f64| v.is_finite() && v > 0.0;
    let nonnegative = |v: f64| v.is_finite() && v >= 0.0;
    if !positive(string.length_m) || !positive(radius_m) {
        return Err(refuse(
            "circular string length and radius must be finite and positive",
        ));
    }
    if !positive(string.tension_n)
        || string.n_modes == 0
        || !nonnegative(string.damping_ratio)
        || !nonnegative(string.polarization_detune)
        || string
            .rayleigh
            .is_some_and(|r| !nonnegative(r.alpha_per_s) || !nonnegative(r.beta_s))
    {
        return Err(refuse(
            "string tension, modal budget, or declared loss parameters are invalid",
        ));
    }
    let density = required_positive(state, DENSITY_PROPERTY, Density::DIMS)?;
    let young = required_positive(state, YOUNG_MODULUS_PROPERTY, Pressure::DIMS)?;
    let radius_squared = radius_m * radius_m;
    let area_m2 = core::f64::consts::PI * radius_squared;
    let second_moment_m4 = area_m2 * (0.25 * radius_squared);
    string.lin_density_kg_m = density * area_m2;
    string.axial_stiffness_n = young * area_m2;
    string.bending_stiffness_n_m2 = young * second_moment_m4;
    string.width_m = 2.0 * radius_m;
    let mass_kg = string.lin_density_kg_m * string.length_m;
    if ![
        area_m2,
        second_moment_m4,
        string.lin_density_kg_m,
        string.axial_stiffness_n,
        string.bending_stiffness_n_m2,
        string.width_m,
        mass_kg,
    ]
    .into_iter()
    .all(positive)
    {
        return Err(refuse(
            "circular string geometry/material products overflow or underflow",
        ));
    }
    let mut identity = DomainHasher::new("org.frankensim.fs-couple.circular-string-specimen.v1");
    identity.update(state.identity().as_bytes());
    identity.update(&string.length_m.to_bits().to_le_bytes());
    identity.update(&radius_m.to_bits().to_le_bytes());
    Ok(ResolvedStringSpecimen {
        string,
        material: state.clone(),
        area_m2,
        second_moment_m4,
        mass_kg,
        specimen_identity: identity.finalize(),
    })
}

fn required_positive(
    state: &ResolvedMaterialStatePoint,
    key: &str,
    dims: Dims,
) -> Result<f64, AcousticRealizeError> {
    let Some(property) = state.property(key) else {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "circular string material needs density and Young's modulus",
        });
    };
    let value = property.value_si();
    if property.requirement().dims() != dims || !value.is_finite() || value <= 0.0 {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "circular string density and Young's modulus must be positive SI quantities",
        });
    }
    Ok(value)
}
