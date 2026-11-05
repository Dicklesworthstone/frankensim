//! Consistent inertia of the SAME quadratic rotations used by DKT bending.
//!
//! For a symmetric section, the kinetic form is rho*h*int(w_dot^2) plus
//! rho*h^3/12*int(beta_x_dot^2 + beta_y_dot^2). The internal DKT rotations
//! are minus the physical slopes; their common sign cancels in this Gram
//! matrix. No shear stiffness or independent rotation coordinate is added.

use crate::{
    AssemblyOptions, PlateChart, PlateError, PlateModel, Stiffener, StiffenerMass, TransverseMass,
    assemble_with_sections_inner, edge_coeffs, h_vectors,
};

/// Panel through-thickness inertia, independent of transverse displacement
/// interpolation and stiffener inertia.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum RotaryMass {
    /// Original nodal slope lumping, rho*h^3*area/36 on each slope DOF.
    #[default]
    Lumped,
    /// Exact area integral of the original quadratic DKT rotation field.
    /// This retains inter-node and displacement/slope kinetic coupling.
    /// It differs from gradients of the optional edge-cubic displacement.
    Dkt,
}

impl PlateChart {
    /// Assemble the original stiffness with independent transverse, rotary
    /// and stiffener kinetic fields. `RotaryMass::Dkt` uses the SAME rotations
    /// as bending stiffness, integrated exactly through polynomial degree four.
    /// It replaces the lumped panel slope term directly during assembly; no
    /// already assembled inertia is subtracted or added twice.
    ///
    /// The resulting `PlateModel.m` is consumed unchanged by modal, harmonic
    /// and transient operators. Existing assembly methods keep their defaults.
    /// This is a Rayleigh/Kirchhoff plate with rotary inertia, not a shear-
    /// deformable plate or a three-dimensional short-wavelength material law.
    ///
    /// # Errors
    /// Returns [`PlateError`] for invalid mesh, sections, supports or beams,
    /// or unrepresentable consistent rotary inertia.
    pub fn assemble_with_kinetic_mass(
        &self,
        stiffeners: &[Stiffener],
        opts: &AssemblyOptions,
        transverse_mass: TransverseMass,
        rotary_mass: RotaryMass,
        stiffener_mass: StiffenerMass,
    ) -> Result<PlateModel, PlateError> {
        assemble_with_sections_inner(
            &self.mesh,
            self.section_field(),
            &self.boundary_nodes,
            stiffeners,
            opts,
            transverse_mass,
            rotary_mass,
            stiffener_mass,
        )
    }
}

/// Six positive degree-four triangle points, with weights summing to one.
/// Each quadratic rotation product is degree four, so no spatial sampling
/// error is introduced by the kinetic integral. A positive Gram sum also
/// avoids the cancellation of signed polynomial-moment coefficients.
const ROTATION_RULE: [([f64; 3], f64); 6] = {
    const A: f64 = 0.445_948_490_915_964_9;
    const B: f64 = 0.091_576_213_509_770_74;
    const WA: f64 = 0.223_381_589_678_011_47;
    const WB: f64 = 0.109_951_743_655_321_87;
    [
        ([1.0 - 2.0 * A, A, A], WA),
        ([A, 1.0 - 2.0 * A, A], WA),
        ([A, A, 1.0 - 2.0 * A], WA),
        ([1.0 - 2.0 * B, B, B], WB),
        ([B, 1.0 - 2.0 * B, B], WB),
        ([B, B, 1.0 - 2.0 * B], WB),
    ]
};

/// Exact kinetic Gram matrix [kg*m^2 in slope/slope entries]. Geometry has
/// already passed the original DKT admission. `scale` is rho*h^3*area/12.
pub(super) fn rotary_mass(x: &[f64; 3], y: &[f64; 3], scale: f64) -> Result<[f64; 81], PlateError> {
    let bad = || PlateError::BadSection {
        what: "consistent DKT rotary inertia is not representable",
    };
    if !scale.is_finite() || scale <= 0.0 {
        return Err(bad());
    }
    let cf = edge_coeffs(x, y);
    let mut mass = [0.0; 81];
    for (l, weight) in ROTATION_RULE {
        let shape = [
            l[0] * (2.0 * l[0] - 1.0),
            l[1] * (2.0 * l[1] - 1.0),
            l[2] * (2.0 * l[2] - 1.0),
            4.0 * l[1] * l[2],
            4.0 * l[2] * l[0],
            4.0 * l[0] * l[1],
        ];
        let (hx, hy) = h_vectors(&cf, &shape);
        for r in 0..9 {
            for c in r..9 {
                mass[9 * r + c] += weight * (hx[r] * hx[c] + hy[r] * hy[c]);
            }
        }
    }
    for r in 0..9 {
        for c in r..9 {
            let value = scale * mass[9 * r + c];
            if !value.is_finite() || (r == c && value <= 0.0) {
                return Err(bad());
            }
            mass[9 * r + c] = value;
            mass[9 * c + r] = value;
        }
    }
    Ok(mass)
}

#[cfg(test)]
#[path = "dkt_kinetic_tests.rs"]
mod tests;
