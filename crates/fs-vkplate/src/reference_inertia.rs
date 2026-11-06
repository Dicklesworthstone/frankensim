//! Spatial reference inventory in the original thin-plate grid owner.
//!
//! The field is sampled once at the existing reference-plane grid nodes. Its
//! nodal quadrature is the same `h²` used by the native transverse displacement,
//! force and kinetic-energy owners. It is a spatial discretization of the
//! staircase domain, not exact cut-cell integration of the analytic circle.
//! The positive samples remain in `Body::mass`; no callback is retained or
//! sampled at a deformed position, temperature or later time.

use crate::{Body, PlaneStressViscosity, PlateSpec};

/// Original density on the plate's reference-plane coordinates `[x,y]` [m].
///
/// This explicitly replaces only `PlateSpec::density`. The caller supplies a
/// coherent density, thickness and constitutive description; selecting this
/// field does not infer elastic, loss or thermal properties from its values.
/// Every sampled density and resulting nodal mass must be finite and positive.
/// A field may jump between grid nodes, with no averaging across that interface.
#[derive(Clone, Copy)]
pub enum PlateReferenceDensity<'a> {
    /// Reference mass per volume [kg/m³]. The existing local thickness supplies
    /// the volume measure: `mass_i = rho0(x_i,y_i) * thickness_i * h²`.
    Volumetric(&'a dyn Fn([f64; 2]) -> f64),
    /// Already integrated reference mass per area [kg/m²], for example a film
    /// or layered inventory: `mass_i = sigma0(x_i,y_i) * h²`. Thickness is still
    /// required by the original material/geometry model and is not multiplied
    /// into this inventory again. This thin-plate lane does not resolve rotary
    /// or first-moment inertia of layers; native MITC shells own those fields.
    Areal(&'a dyn Fn([f64; 2]) -> f64),
}

impl Body {
    /// Prepare the original plate with an explicit spatial reference inventory.
    /// Sampling happens exactly once per admitted grid node and retains only
    /// the existing mass array. Clones, warm steps, contacts and radiation reuse
    /// that array without calling or retaining `density`.
    ///
    /// `s.t60` is ignored: its scalar uniform-density dispersion does not define
    /// a graded surface's modal loss. Surface material loss is absent on this
    /// entry; authored support/muffle dashpots remain. Use
    /// [`Self::plate_with_reference_density_and_viscosity`] for sourced bending
    /// loss. Stiffness, prestress and nonlinear coupling retain their original
    /// operators. No geometry or mass renormalization is performed.
    ///
    /// The receiving instrument's existing `stability()` admission must be
    /// checked at its actual clock: it uses these masses with the original
    /// stiffness and loss operators. Density or total mass alone is not a
    /// stability certificate. Existing power-iteration and continuum/model
    /// limits of that admission remain in force.
    ///
    /// # Errors
    /// Original geometry/operator refusals, invalid sampled density/thickness,
    /// or nonpositive, underflowed or nonfinite nodal mass. No partial body is
    /// returned and no numerical mass floor replaces a light inventory.
    pub fn plate_with_reference_density(
        name: &str,
        s: &PlateSpec<'_>,
        density: PlateReferenceDensity<'_>,
    ) -> Result<Self, String> {
        Self::plate_impl(name, s, None, false, Some(density))
    }

    /// The same original reference inventory with a sourced plane-stress
    /// Kelvin–Voigt bending law. The actual mass-normalized viscosity bound and
    /// source frequency domain enter the existing instrument clock admission.
    /// The field changes inertia only; the supplied loss law is independent.
    ///
    /// # Errors
    /// The refusals of [`Self::plate_with_reference_density`] or the existing
    /// sourced viscosity constructor, including absent bending capability.
    pub fn plate_with_reference_density_and_viscosity(
        name: &str,
        s: &PlateSpec<'_>,
        density: PlateReferenceDensity<'_>,
        source: &PlaneStressViscosity,
        angle_rad: f64,
    ) -> Result<Self, String> {
        if !s.bending {
            return Err("bending viscosity requires the plate bending operator".into());
        }
        Self::plate_impl(name, s, Some((source, angle_rad)), false, Some(density))
    }
}

pub(super) fn validate_grid(s: &PlateSpec<'_>) -> Result<(), String> {
    // This is representability admission for the original integer grid, not
    // a new mesh-resolution policy or a bound on allocated process memory.
    let extent = s.radius / s.grid;
    if !s.radius.is_finite()
        || s.radius <= 0.0
        || !s.hole.is_finite()
        || s.hole < 0.0
        || s.hole >= s.radius
        || !s.grid.is_finite()
        || s.grid <= 0.0
        || !extent.is_finite()
        || extent.ceil() > 2_147_483_643.0
    {
        return Err("plate requires finite positive radius/grid, 0 <= hole < radius, and representable grid coordinates".into());
    }
    Ok(())
}

pub(super) fn nodal_masses(
    name: &str,
    s: &PlateSpec<'_>,
    xy: &[[f64; 2]],
    thickness: &[f64],
    field: Option<PlateReferenceDensity<'_>>,
) -> Result<Vec<f64>, String> {
    if field.is_none() && (!s.density.is_finite() || s.density <= 0.0) {
        return Err(format!(
            "{name}: reference volumetric density must be finite and positive"
        ));
    }
    let h = s.grid;
    let mut masses = Vec::with_capacity(xy.len());
    for (node, (&position, &t)) in xy.iter().zip(thickness).enumerate() {
        if !t.is_finite() || t <= 0.0 {
            return Err(format!(
                "{name}: reference thickness at node {node} must be finite and positive"
            ));
        }
        let (density, areal) = match field {
            None => (s.density, false),
            Some(PlateReferenceDensity::Volumetric(sample)) => (sample(position), false),
            Some(PlateReferenceDensity::Areal(sample)) => (sample(position), true),
        };
        if !density.is_finite() || density <= 0.0 {
            return Err(format!(
                "{name}: reference density at node {node} must be finite and positive"
            ));
        }
        // Keep the scalar volumetric constructor's multiplication order.
        let mass = if areal {
            density * h * h
        } else {
            density * t * h * h
        };
        if !mass.is_finite() || mass <= 0.0 {
            return Err(format!(
                "{name}: reference mass at node {node} is nonfinite or underflowed"
            ));
        }
        masses.push(mass);
    }
    Ok(masses)
}

/// Retain the original power-iteration arithmetic on ordinary vectors. A
/// finite operator image can have an overflowing squared norm when physical
/// inertia is very light; normalizing it by infinity would erase that mode.
/// Scaling is needed only when squared storage overflows or underflows.
pub(super) fn spectral_norm(values: &[f64]) -> f64 {
    let squared = values.iter().map(|a| a * a).sum::<f64>();
    if squared.is_finite() && squared >= f64::MIN_POSITIVE {
        return squared.sqrt();
    }
    let mut scale = 0.0_f64;
    for &value in values {
        if !value.is_finite() {
            return f64::INFINITY;
        }
        scale = scale.max(value.abs());
    }
    if scale == 0.0 {
        return 0.0;
    }
    let squared = values
        .iter()
        .map(|value| (value / scale).powi(2))
        .sum::<f64>();
    scale * squared.sqrt()
}

#[cfg(test)]
mod tests;
