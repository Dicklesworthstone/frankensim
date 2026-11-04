//! Source-resolved liquid properties for the native Newtonian film model.
//!
//! The same material owner already supplies momentum and capillarity to the
//! free-surface LBM. This adapter borrows that complete state, retaining its
//! original viscosity basis, query, property receipts and interface context.
//! It requests no caloric properties and performs no new material selection.

use fs_material::liquid::surface::ResolvedLiquidSurfaceState;

use super::{FilmError, FilmMaterial, NusseltFilm};

/// A frozen Newtonian liquid binding for repeated local film evaluations.
///
/// The material owner admits its source conditions. This adapter additionally
/// rejects states requiring evolving shear/pressure rheology, because the
/// Nusselt/Benney model has one constant Newtonian viscosity. Temperature or
/// composition changes require a newly admitted source state; this binding
/// supplies no heat evolution, phase-stability or source-accuracy certificate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialFilmLiquid<'a> {
    material: &'a ResolvedLiquidSurfaceState,
}

impl<'a> MaterialFilmLiquid<'a> {
    /// Bind the exact state used by the existing material free-surface owner.
    /// No properties, source receipts or source coordinates are cloned.
    ///
    /// # Errors
    /// Refuses a state that requires an evolving constitutive owner instead of
    /// silently treating its sampled apparent viscosity as Newtonian.
    pub fn try_new(material: &'a ResolvedLiquidSurfaceState) -> Result<Self, FilmError> {
        if material.rheology_axes().is_some() || material.rheology_pressure_axis().is_some() {
            return Err(FilmError::InvalidInput(
                "Newtonian film cannot freeze an evolving rheology source",
            ));
        }
        Ok(Self { material })
    }

    /// Original source owner, including every retained receipt and query value.
    #[must_use]
    pub const fn material(self) -> &'a ResolvedLiquidSurfaceState {
        self.material
    }

    /// SI coefficients for the existing film and vessel evaluators.
    ///
    /// `None` in the admitted source means capillarity was explicitly neglected;
    /// only that declaration becomes zero. Missing active surface data already
    /// refuses in the material owner. Retain this binding beside numerical
    /// results when their source identity is required.
    #[must_use]
    pub const fn coefficients(self) -> FilmMaterial {
        FilmMaterial {
            density_kg_m3: self.material.density_kg_m3(),
            kinematic_viscosity_m2_s: self.material.kinematic_viscosity_m2_s(),
            surface_tension_n_m: match self.material.surface_tension_n_m() {
                Some(value) => value,
                None => 0.0,
            },
        }
    }

    /// Evaluate the native Nusselt base flow with these SAME source values.
    /// All physical inputs and numerical refusals belong to [`NusseltFilm`].
    /// This fixed-cost call allocates nothing and performs no source queries.
    pub fn film(
        self,
        flow_per_width_m2_s: f64,
        tangential_acceleration_m_s2: f64,
        normal_acceleration_m_s2: f64,
    ) -> Result<NusseltFilm, FilmError> {
        NusseltFilm::new(
            self.coefficients(),
            flow_per_width_m2_s,
            tangential_acceleration_m_s2,
            normal_acceleration_m_s2,
        )
    }
}
