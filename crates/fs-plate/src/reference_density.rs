//! Reference-area density on the original flat-plate kinetic fields.
//!
//! A positive P1 density is integrated jointly with the chosen displacement
//! and rotation fields. Taking its mean before the Gram integral would lose
//! mass placement even when the total mass happened to be correct.

use crate::{PlateChart, PlateError, RotaryMass, TransverseMass};

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
        self.reference_density = Some(density);
        Ok(self)
    }

    /// Original retained corner densities, or `None` when the current section
    /// field supplies scalar density. Assembly rechecks coverage because the
    /// chart's mesh remains publicly editable.
    #[must_use]
    pub fn reference_density(&self) -> Option<&[[f64; 3]]> {
        self.reference_density.as_deref()
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
    let bad = || PlateError::BadSection {
        what: "graded plate kinetic inertia is not representable",
    };
    let (scale, density) = normalized(density);
    let translation_scale = scale * thickness * area;
    let rotation_scale = scale * thickness * thickness * thickness / 12.0 * area;
    if [translation_scale, rotation_scale]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return Err(bad());
    }
    let mut mass = [0.0; 81];
    match transverse {
        TransverseMass::Lumped => {
            let diagonal = lumped(translation_scale, density);
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
                    let value = translation_scale * moment(density, powers);
                    mass[9 * (3 * i) + 3 * j] = value;
                    mass[9 * (3 * j) + 3 * i] = value;
                }
            }
        }
        TransverseMass::EdgeCubic => {
            mass = crate::edge_cubic::edge_cubic_mass_graded(x, y, translation_scale, density);
        }
    }
    match rotary {
        RotaryMass::Lumped => {
            let diagonal = lumped(rotation_scale, density);
            for i in 0..3 {
                for component in [1, 2] {
                    let dof = 3 * i + component;
                    mass[9 * dof + dof] += diagonal[i];
                }
            }
        }
        RotaryMass::Dkt => {
            let rotation = crate::dkt_kinetic::rotary_mass_graded(x, y, rotation_scale, density)?;
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
