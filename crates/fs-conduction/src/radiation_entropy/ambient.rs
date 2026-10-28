//! Material entropy exchanged by the ORIGINAL ambient radiation law.
//!
//! This is the sum of material reservoir storage rates in stationary radiative
//! exchange, not the entropy flux of the photon field. In a black isothermal
//! bath the field has no net storage; its emission/absorption irreversibility
//! is included in q*(1/T_bath-1/T_surface). Reflected and emitted occupation
//! numbers cannot be treated as independent additive entropy spectra.

use crate::{AmbientRadiationPatch, ConductionError, require_finite};

/// Endpoint rates per original emitting area, W/(m² K), except the heat field.
/// This is an observation of the existing source-bound patch; it neither changes
/// a physical state nor adds a second heat or entropy debit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialRadiationEntropyExchange {
    /// Original net radiative heat flux, W/m², positive from surface to bath.
    pub heat_out_w_m2: f64,
    /// Material storage contribution -q/T_surface, W/(m² K).
    pub surface_storage_rate_w_m2_k: f64,
    /// Black reservoir storage contribution +q/T_bath, W/(m² K).
    pub ambient_storage_rate_w_m2_k: f64,
    /// Total surface-plus-bath production, W/(m² K). Finite-temperature
    /// transfer is nonnegative; equal temperatures give an exact zero.
    pub generation_w_m2_k: f64,
    /// Rounded sum of both storage rates minus stable production, W/(m² K).
    /// This arithmetic discrepancy is NOT additional entropy production.
    pub arithmetic_residual_w_m2_k: f64,
}

impl AmbientRadiationPatch {
    /// Measure material entropy using the SAME source-selected nonlinear heat.
    /// Source-temperature/domain refusals propagate unchanged. The bath is the
    /// original black isothermal view-factor-one reservoir, not a surface query.
    ///
    /// This preserves the owner's area-mean patch approximation. Multiplying by
    /// area and duration supplies endpoint quadrature, not an exact transient
    /// entropy change. For a nonuniform surface use local spatial temperatures;
    /// an area-mean temperature is not a substitute for integrating q/T locally.
    /// Photon emission entropy is a different observable; do not add it to this
    /// closed surface-plus-bath production or interpret it as q/T_surface.
    pub fn material_entropy_exchange(
        &self, surface_temperature_k: f64,
    ) -> Result<MaterialRadiationEntropyExchange, ConductionError> {
        let heat = self.heat_flux_w_m2(surface_temperature_k)?;
        let ambient = self.ambient_temperature_k();
        let surface_rate = require_finite("radiating surface entropy rate", -heat/surface_temperature_k)?;
        let ambient_rate = require_finite("radiative reservoir entropy rate", heat/ambient)?;
        // Difference of temperatures is resolved before division; subtracting
        // rounded reciprocal temperatures can erase near-equilibrium production.
        let generation = require_finite("radiative entropy production",
            (heat/surface_temperature_k)*((surface_temperature_k-ambient)/ambient))?;
        if generation < 0.0
            || (surface_temperature_k != ambient && (heat == 0.0 || generation == 0.0))
        {
            return Err(ConductionError::Radiation {
                surface: self.region().into(),
                what: "nonzero radiative entropy production is unrepresentable or nonpassive".into(),
                fix: "use a representable source-temperature and emissivity scale".into(),
            });
        }
        let residual = require_finite("radiative entropy arithmetic residual",
            (surface_rate+ambient_rate)-generation)?;
        Ok(MaterialRadiationEntropyExchange {
            heat_out_w_m2: heat,
            surface_storage_rate_w_m2_k: surface_rate,
            ambient_storage_rate_w_m2_k: ambient_rate,
            generation_w_m2_k: generation,
            arithmetic_residual_w_m2_k: residual,
        })
    }
}
