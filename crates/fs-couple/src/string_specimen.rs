//! Uniform circular specimens lowered to the existing prestressed-string solver.
//!
//! This is the homogeneous isotropic Euler–Bernoulli/Kirchhoff–Carrier rung:
//! `A = pi r²`, `I = pi r⁴/4`, `mu = rho A`, `EA = E A`, `EI = E I`.
//! Length and radius describe the specimen at the resolved material state.
//! Tension, support motion, modal budget, and declared losses remain mechanical
//! inputs. Rebinding can hold tension or axial extension fixed, or solve the
//! tension for a declared undamped fundamental. Geometry is the current loaded
//! geometry; no thermal evolution or fixed-mass comparison is implied.

use fs_blake3::{ContentHash, DomainHasher};
use fs_material::state_point::{
    DENSITY_PROPERTY, ResolvedMaterialStatePoint, YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Density, Dims, Pressure};
use fs_scenario::PrestressedString;

use crate::acoustic_realize::AcousticRealizeError;

/// Mechanical prestress prescribed independently of material-card authority.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum StringPrestress {
    /// Applied tensile force [N], independent of the material's stiffness.
    FixedTension(f64),
    /// Linear elastic tension `EA * (L - L0) / L0` at the declared loaded
    /// length `L`. Both transverse endpoints must be stationary (pinned, not
    /// rotationally clamped). Cross-section and density describe this loaded
    /// state; Poisson contraction is not solved.
    FixedExtension {
        /// Stress-free length `L0` [m] at the same material state.
        stress_free_length_m: f64,
        /// Caller-declared upper bound on positive engineering strain.
        /// This is an applicability assumption, not a measured yield limit.
        linear_strain_limit: f64,
    },
    /// Solve tension for the primary, undamped pinned-pinned bending mode [Hz].
    /// Includes bending stiffness. It does not tune a damped/nonlinear pressure
    /// waveform or the independently detuned second polarization.
    TargetFundamentalHz(f64),
}

/// One material-bound specimen plus its independently authored mechanical state.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedStringSpecimen {
    string: PrestressedString,
    prestress: StringPrestress,
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

    /// Original mechanical prescription, separate from resolved material claims.
    #[must_use]
    pub const fn prestress(&self) -> StringPrestress {
        self.prestress
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
    string: PrestressedString,
    radius_m: f64,
    state: &ResolvedMaterialStatePoint,
) -> Result<ResolvedStringSpecimen, AcousticRealizeError> {
    with_uniform_circular_material_and_prestress(
        string,
        radius_m,
        state,
        StringPrestress::FixedTension(string.tension_n),
    )
}

/// Bind current circular geometry and material, then resolve one mechanical
/// prestress prescription. The template's independent mass, stiffness, width,
/// and tension are replaced; its supports, losses and mode budget are retained.
///
/// For the pinned-pinned Euler–Bernoulli model, `omega² = (T k² + EI k⁴)/mu`
/// with `k = pi/L`. Inverting gives `T = mu (2 L f)² - EI (pi/L)²`.
/// A target at or below the zero-tension bending frequency is outside the
/// existing tensile-string solver. See [FHWA-HRT-14-049, chapter 2][fhwa],
/// equations of motion for a taut string with flexural stiffness.
///
/// [fhwa]: https://www.fhwa.dot.gov/publications/research/infrastructure/structures/bridge/14049/003.cfm
///
/// Fixed extension uses a homogeneous linear elastic law and the supplied
/// stress-free length at this state. It does not infer thermal strain, yield,
/// creep, transverse contraction, or conservation across different specimens.
/// The resulting specimen identity still excludes loading; `prestress()` keeps
/// the authored prescription without merging it into the material receipts.
///
/// # Errors
/// Refuses invalid geometry/material/template values, slack/compressive or
/// unrepresentable tension, strain beyond its declared limit, and extension or
/// pitch prescriptions with the realizer's moving-end basis.
pub fn with_uniform_circular_material_and_prestress(
    mut string: PrestressedString,
    radius_m: f64,
    state: &ResolvedMaterialStatePoint,
    prestress: StringPrestress,
) -> Result<ResolvedStringSpecimen, AcousticRealizeError> {
    let refuse = |what| AcousticRealizeError::InvalidDescription { what };
    let positive = |v: f64| v.is_finite() && v > 0.0;
    let nonnegative = |v: f64| v.is_finite() && v >= 0.0;
    if !positive(string.length_m) || !positive(radius_m) {
        return Err(refuse(
            "circular string length and radius must be finite and positive",
        ));
    }
    if string.n_modes == 0
        || !nonnegative(string.damping_ratio)
        || !nonnegative(string.polarization_detune)
        || string
            .rayleigh
            .is_some_and(|r| !nonnegative(r.alpha_per_s) || !nonnegative(r.beta_s))
    {
        return Err(refuse(
            "string modal budget or declared loss parameters are invalid",
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
    string.tension_n = resolve_prestress(string, prestress)?;
    let mut identity = DomainHasher::new("org.frankensim.fs-couple.circular-string-specimen.v1");
    identity.update(state.identity().as_bytes());
    identity.update(&string.length_m.to_bits().to_le_bytes());
    identity.update(&radius_m.to_bits().to_le_bytes());
    Ok(ResolvedStringSpecimen {
        string,
        prestress,
        material: state.clone(),
        area_m2,
        second_moment_m4,
        mass_kg,
        specimen_identity: identity.finalize(),
    })
}

fn resolve_prestress(
    string: PrestressedString,
    prestress: StringPrestress,
) -> Result<f64, AcousticRealizeError> {
    let refuse = |what| AcousticRealizeError::InvalidDescription { what };
    let positive = |v: f64| v.is_finite() && v > 0.0;
    if string.moving_end && !matches!(prestress, StringPrestress::FixedTension(_)) {
        return Err(refuse(
            "extension and target pitch require the stationary-end sine basis",
        ));
    }
    let tension = match prestress {
        StringPrestress::FixedTension(tension) => tension,
        StringPrestress::FixedExtension {
            stress_free_length_m,
            linear_strain_limit,
        } => {
            if !positive(stress_free_length_m) || !positive(linear_strain_limit) {
                return Err(refuse(
                    "stress-free length and strain limit must be positive",
                ));
            }
            let strain = (string.length_m - stress_free_length_m) / stress_free_length_m;
            if !positive(strain) || strain > linear_strain_limit {
                return Err(refuse(
                    "string extension is outside its tensile strain domain",
                ));
            }
            string.axial_stiffness_n * strain
        }
        StringPrestress::TargetFundamentalHz(frequency_hz) => {
            if !positive(frequency_hz) {
                return Err(refuse("target fundamental frequency must be positive Hz"));
            }
            let speed = 2.0 * (string.length_m * frequency_hz);
            let wave_number = core::f64::consts::PI / string.length_m;
            string.lin_density_kg_m * speed.powi(2)
                - string.bending_stiffness_n_m2 * wave_number.powi(2)
        }
    };
    if !positive(tension) {
        return Err(refuse(
            "prestress must resolve to finite positive tension; slack or bending-only states need another solver",
        ));
    }
    Ok(tension)
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
