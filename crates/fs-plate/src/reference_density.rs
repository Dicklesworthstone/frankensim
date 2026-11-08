//! Reference-area density on the original flat-plate kinetic fields.
//!
//! A positive P1 density is integrated jointly with the chosen displacement
//! and rotation fields. Taking its mean before the Gram integral would lose
//! mass placement even when the total mass happened to be correct.

use crate::{PlateChart, PlateError, RotaryMass, TransverseMass};

/// Centered, fixed-reference thickness mass moments at one in-plane point.
///
/// `areal_mass_kg_m2 = integral rho(z) dz` weights transverse velocity.
/// `second_moment_kg = integral rho(z) z² dz` weights the two rotations.
/// The declared inventory is symmetric about the plate reference surface;
/// its first thickness moment is zero. Independent moments retain the
/// different rotary inertia of dense face layers and a dense core even
/// when their total mass is equal.
///
/// These are mass moments, not an effective density or an elastic section.
/// They may include symmetric nonstructural mass outside the elastic thickness.
/// No thickness profile, stiffness or density evolution is inferred.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlateSectionMassMoments {
    areal_mass_kg_m2: f64,
    second_moment_kg: f64,
}

impl PlateSectionMassMoments {
    /// Retain positive finite zeroth and second thickness mass moments.
    ///
    /// # Errors
    /// Nonfinite or nonpositive moments refuse. A zero-thickness mass sheet
    /// with no rotary inertia is outside this finite-thickness plate model.
    pub fn try_new(areal_mass_kg_m2: f64, second_moment_kg: f64) -> Result<Self, PlateError> {
        if [areal_mass_kg_m2, second_moment_kg]
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(PlateError::BadSection {
                what: "plate thickness mass moments must be finite and positive",
            });
        }
        Ok(Self {
            areal_mass_kg_m2,
            second_moment_kg,
        })
    }

    /// Zeroth thickness moment, `integral rho dz` [kg/m²].
    #[must_use]
    pub fn areal_mass_kg_m2(self) -> f64 {
        self.areal_mass_kg_m2
    }

    /// Second thickness moment, `integral rho z² dz` [kg].
    #[must_use]
    pub fn second_moment_kg(self) -> f64 {
        self.second_moment_kg
    }
}

/// One authoritative reference inventory; selecting a new field replaces it.
#[derive(Debug, Clone)]
pub(super) enum ReferenceMass {
    Density(Vec<[f64; 3]>),
    SectionMoments(Vec<[PlateSectionMassMoments; 3]>),
}

impl ReferenceMass {
    pub(super) fn validate(&self, elements: usize) -> Result<(), PlateError> {
        match self {
            Self::Density(density) => validate(elements, density),
            Self::SectionMoments(moments) => validate_moments(elements, moments),
        }
    }

    pub(super) fn density(&self) -> Option<&[[f64; 3]]> {
        match self {
            Self::Density(density) => Some(density),
            Self::SectionMoments(_) => None,
        }
    }

    pub(super) fn section_moments(&self) -> Option<&[[PlateSectionMassMoments; 3]]> {
        match self {
            Self::Density(_) => None,
            Self::SectionMoments(moments) => Some(moments),
        }
    }
}

impl PlateChart {
    /// Retain reference density [kg/m³] at each triangle's three corners, in
    /// the exact order of `mesh.tris`. Shared-node values may differ between
    /// triangles, so material jumps do not require averaging across interfaces.
    /// This field replaces the section's scalar density in EVERY chart mass
    /// assembly method; section stiffness/thickness and beam density are separate.
    ///
    /// The plate is symmetric through its supplied, elementwise constant
    /// thickness. This models in-plane grading, not an eccentric laminate.
    /// Coordinates, connectivity and thickness describe the same reference
    /// inventory; changing them requires updating that inventory explicitly.
    /// Constitutive or temperature updates must not resample current density
    /// into this field while claiming a closed body's mass is conserved.
    ///
    /// # Errors
    /// Incomplete coverage or any nonfinite/nonpositive corner density refuses.
    pub fn with_reference_density(mut self, density: Vec<[f64; 3]>) -> Result<Self, PlateError> {
        validate(self.mesh.tris.len(), &density)?;
        self.reference_mass = Some(ReferenceMass::Density(density));
        Ok(self)
    }

    /// Original retained corner densities, or `None` when the current section
    /// field supplies scalar density or explicit section moments are selected.
    /// Assembly rechecks coverage because the
    /// chart's mesh remains publicly editable.
    #[must_use]
    pub fn reference_density(&self) -> Option<&[[f64; 3]]> {
        self.reference_mass
            .as_ref()
            .and_then(ReferenceMass::density)
    }

    /// Retain independent centered thickness mass moments at each triangle's
    /// three corners, in `mesh.tris` order. Both moments interpolate linearly
    /// in-plane; jumps between triangles remain explicit.
    ///
    /// This replaces the scalar/P1 density inventory in EVERY chart mass
    /// assembly method. Elastic thickness/stiffness and stiffener mass remain
    /// independently owned. The original displacement and rotation fields
    /// consume their respective moments without averaging them together.
    /// Coordinates and moments describe a fixed reference inventory; editing
    /// the geometry requires explicitly updating that inventory.
    ///
    /// The flat plate has no membrane coordinates, so eccentric first moments
    /// and membrane/bending kinetic coupling require the shell model.
    ///
    /// # Errors
    /// Incomplete coverage or unrepresentable moments refuse. Assembly also
    /// checks coverage after any public mesh edit and rejects mass overflow.
    pub fn with_reference_section_mass(
        mut self,
        moments: Vec<[PlateSectionMassMoments; 3]>,
    ) -> Result<Self, PlateError> {
        validate_moments(self.mesh.tris.len(), &moments)?;
        self.reference_mass = Some(ReferenceMass::SectionMoments(moments));
        Ok(self)
    }

    /// Retained section moments, or `None` when scalar/P1 density is selected.
    #[must_use]
    pub fn reference_section_mass(&self) -> Option<&[[PlateSectionMassMoments; 3]]> {
        self.reference_mass
            .as_ref()
            .and_then(ReferenceMass::section_moments)
    }
}

pub(super) fn validate(elements: usize, density: &[[f64; 3]]) -> Result<(), PlateError> {
    if density.len() != elements {
        return Err(PlateError::SectionCount {
            expected: elements,
            actual: density.len(),
        });
    }
    if density
        .iter()
        .flatten()
        .any(|rho| !rho.is_finite() || *rho <= 0.0)
    {
        return Err(PlateError::BadSection {
            what: "plate reference density must be finite and positive at every corner",
        });
    }
    Ok(())
}

fn validate_moments(
    elements: usize,
    moments: &[[PlateSectionMassMoments; 3]],
) -> Result<(), PlateError> {
    if moments.len() != elements {
        return Err(PlateError::SectionCount {
            expected: elements,
            actual: moments.len(),
        });
    }
    // Values are immutable and admitted by PlateSectionMassMoments::try_new.
    Ok(())
}

/// Separate the density scale before products and polynomial contractions.
pub(super) fn normalized(density: [f64; 3]) -> (f64, [f64; 3]) {
    let scale = density.into_iter().fold(0.0_f64, f64::max);
    (scale, density.map(|rho| rho / scale))
}

/// Exact mean of `rho(lambda) * lambda^powers` on the unit-area triangle.
/// Density is normalized to at most one. The P1 factor adds one degree:
/// `int(rho*lambda^a)/A = 2 prod(a!)/(sum(a)+2)! *
///                       sum_i rho_i*(a_i+1)/(sum(a)+3)`.
/// All consumers have total degree at most six before this P1 factor.
pub(super) fn moment(density: [f64; 3], powers: [usize; 3]) -> f64 {
    const FACT: [f64; 9] = [1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0];
    let degree: usize = powers.iter().sum();
    let mean = 2.0 * powers.into_iter().map(|p| FACT[p]).product::<f64>() / FACT[degree + 2];
    let weighted = (0..3)
        .map(|i| density[i] * ((powers[i] + 1) as f64 / (degree + 3) as f64))
        .sum::<f64>();
    mean * weighted
}

/// Positive P1 row sums. They preserve inventory and its first in-plane
/// moments; lumping does not reproduce the consistent quadratic kinetic form.
pub(super) fn lumped(scale: f64, density: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| {
        let mut powers = [0; 3];
        powers[i] = 1;
        scale * moment(density, powers)
    })
}

#[allow(clippy::too_many_arguments)]
pub(super) fn panel_mass(
    x: &[f64; 3],
    y: &[f64; 3],
    density: [f64; 3],
    thickness: f64,
    area: f64,
    transverse: TransverseMass,
    rotary: RotaryMass,
) -> Result<[f64; 81], PlateError> {
    let (scale, density) = normalized(density);
    let translation_scale = scale * thickness * area;
    let rotation_scale = scale * thickness * thickness * thickness / 12.0 * area;
    kinetic_mass(
        x,
        y,
        translation_scale,
        density,
        rotation_scale,
        density,
        transverse,
        rotary,
    )
}

pub(super) fn section_panel_mass(
    x: &[f64; 3],
    y: &[f64; 3],
    moments: [PlateSectionMassMoments; 3],
    area: f64,
    transverse: TransverseMass,
    rotary: RotaryMass,
) -> Result<[f64; 81], PlateError> {
    let (translation, translation_weights) = normalized(moments.map(|m| m.areal_mass_kg_m2));
    let (rotation, rotation_weights) = normalized(moments.map(|m| m.second_moment_kg));
    kinetic_mass(
        x,
        y,
        translation * area,
        translation_weights,
        rotation * area,
        rotation_weights,
        transverse,
        rotary,
    )
}

#[allow(clippy::too_many_arguments)]
fn kinetic_mass(
    x: &[f64; 3],
    y: &[f64; 3],
    translation_scale: f64,
    translation_weights: [f64; 3],
    rotation_scale: f64,
    rotation_weights: [f64; 3],
    transverse: TransverseMass,
    rotary: RotaryMass,
) -> Result<[f64; 81], PlateError> {
    let bad = || PlateError::BadSection {
        what: "graded plate kinetic inertia is not representable",
    };
    if [translation_scale, rotation_scale]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return Err(bad());
    }
    let mut mass = [0.0; 81];
    match transverse {
        TransverseMass::Lumped => {
            let diagonal = lumped(translation_scale, translation_weights);
            for i in 0..3 {
                mass[9 * (3 * i) + 3 * i] = diagonal[i];
            }
        }
        TransverseMass::Linear => {
            for i in 0..3 {
                for j in i..3 {
                    let mut powers = [0; 3];
                    powers[i] += 1;
                    powers[j] += 1;
                    let value = translation_scale * moment(translation_weights, powers);
                    mass[9 * (3 * i) + 3 * j] = value;
                    mass[9 * (3 * j) + 3 * i] = value;
                }
            }
        }
        TransverseMass::EdgeCubic => {
            mass = crate::edge_cubic::edge_cubic_mass_graded(
                x,
                y,
                translation_scale,
                translation_weights,
            );
        }
    }
    match rotary {
        RotaryMass::Lumped => {
            let diagonal = lumped(rotation_scale, rotation_weights);
            for i in 0..3 {
                for component in [1, 2] {
                    let dof = 3 * i + component;
                    mass[9 * dof + dof] += diagonal[i];
                }
            }
        }
        RotaryMass::Dkt => {
            let rotation =
                crate::dkt_kinetic::rotary_mass_graded(x, y, rotation_scale, rotation_weights)?;
            for (value, contribution) in mass.iter_mut().zip(rotation) {
                *value += contribution;
            }
        }
    }
    if mass.iter().any(|x| !x.is_finite()) || (0..9).any(|i| mass[9 * i + i] <= 0.0) {
        return Err(bad());
    }
    Ok(mass)
}
