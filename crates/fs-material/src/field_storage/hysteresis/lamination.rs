//! Classical eddy-current reaction in insulated thin laminations.
//!
//! For uniform in-plane B, Faraday's law gives E(y) = -y dB/dt.
//! Integrating sigma E² across thickness d gives
//! P = volume * sigma * d²/12 * (dB/dt)². With winding voltage
//! V = N A dB/dt, this is a passive shunt across the ORIGINAL core port,
//! R = 12 N² A rho/(path_length * d²). Its current reacts on the driving
//! circuit; its native resistor work is deposited once in the conducting core.
//!
//! This is the uniform-flux, zero skin-effect term of lamination
//! homogenisation; see Gyselinck et al., "Time-domain finite-element
//! homogenisation of laminated iron cores with net circulating currents"
//! (2015), https://orbi.uliege.be/handle/2268/191457 . Insulation prevents net
//! current between sheets. Through-thickness diffusion, edge/fringing currents,
//! excess/domain-wall loss and insulation damage require their own owners.
//! It must not accompany another resolved eddy-current loss in the SAME volume.
//! The original hysteresis law must describe quasistatic pinning rather than a
//! frequency-dependent total-loss fit. No core-loss coefficient is fitted here.

use super::UniformHystereticCore;
use crate::conductor::{ConductorError, ResolvedConductor, ResolvedConductorState};

/// Physical sheet dimensions and an explicit geometric admission envelope.
/// The magnetic area/path and signed turns come from the original core.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThinLaminationGeometry {
    pub thickness_m: f64,
    /// In-plane transverse width of each conducting sheet, not stack height.
    pub width_m: f64,
    /// Caller-admitted upper bound on d/min(width,path), strictly below one.
    /// This geometric bound is not a skin-effect or physical-error certificate.
    pub maximum_thickness_ratio: f64,
}

/// Source-bound shunt reduction, with physical geometry retained separately
/// from its volume-preserving native transport coordinates.
#[derive(Clone, Debug, PartialEq)]
pub struct LaminatedCoreEddyCurrent {
    core: UniformHystereticCore,
    geometry: ThinLaminationGeometry,
    equivalent_conductor: ResolvedConductor,
}

fn invalid(quantity: &'static str) -> ConductorError {
    ConductorError::InvalidInput { quantity }
}

impl LaminatedCoreEddyCurrent {
    /// Bind the SAME source material and physical core to a lamination shunt.
    /// The native conductor retains the original independent rho/sigma basis,
    /// source keys, validity and uncertainty. Anisotropic, field/frequency-
    /// dependent or phase-changing transport is outside this scalar solid law.
    ///
    /// The caller must separately admit a low-frequency uniform-B regime: the
    /// sheet thickness must be small against magnetic skin depth. This method
    /// does not infer a bandwidth from timestep or from a DC material query.
    pub fn new(
        core: &UniformHystereticCore,
        conductor: ResolvedConductorState,
        geometry: ThinLaminationGeometry,
    ) -> Result<Self, ConductorError> {
        let area = core.linkage_per_t / core.turns;
        let path = core.current_per_field * core.turns;
        let ThinLaminationGeometry {
            thickness_m: d,
            width_m: width,
            maximum_thickness_ratio: limit,
        } = geometry;
        if !d.is_finite()
            || d <= 0.0
            || !width.is_finite()
            || width <= 0.0
            || !limit.is_finite()
            || limit <= 0.0
            || limit >= 1.0
            || d / width > limit
            || d / path > limit
            || d * width > area
        {
            return Err(invalid(
                "positive thin-lamination geometry inside its declared aspect envelope and core area",
            ));
        }
        if conductor.liquid_mass_fraction_axis().is_some()
            || !conductor
                .source_property()
                .is_static_at_fixed_external_conditions()
        {
            return Err(invalid(
                "lamination requires a static scalar solid Ohmic source",
            ));
        }
        if let Some(source) = core.source() {
            if source.material().card_identity() != conductor.material().card_identity()
                || source.material().query() != conductor.material().query()
            {
                return Err(invalid(
                    "magnetic fit and eddy transport require the same core card and complete query",
                ));
            }
        }
        // These are REDUCED transport coordinates, not a physical wire shape.
        // L_e/A_e = 12 N² A/(l d²), while L_e*A_e = original solid volume.
        // Thus the native material/enthalpy owner keeps exactly one inventory.
        let length = 12.0_f64.sqrt() * core.turns.abs() * area / d;
        let cross_section = core.volume_m3 / length;
        let equivalent_conductor = conductor.into_uniform_conductor(length, cross_section)?;
        Ok(Self {
            core: core.clone(),
            geometry,
            equivalent_conductor,
        })
    }

    pub const fn geometry(&self) -> ThinLaminationGeometry {
        self.geometry
    }
    pub const fn turns(&self) -> f64 {
        self.core.turns
    }
    pub fn magnetic_area_m2(&self) -> f64 {
        self.core.linkage_per_t / self.core.turns
    }
    pub fn path_length_m(&self) -> f64 {
        self.core.current_per_field * self.core.turns
    }
    pub const fn volume_m3(&self) -> f64 {
        self.core.volume_m3
    }

    /// Native source-backed shunt; its L/A are equivalent transport coordinates.
    /// They must not be reused as actual wire or thermal-contact geometry.
    pub const fn equivalent_conductor(&self) -> &ResolvedConductor {
        &self.equivalent_conductor
    }

    /// Check the original magnetic fit and physical geometry without requiring
    /// its retained source query to remain at the initial temperature.
    pub fn matches_core(&self, core: &UniformHystereticCore) -> bool {
        self.core.law == core.law
            && self.core.turns.to_bits() == core.turns.to_bits()
            && self.core.linkage_per_t.to_bits() == core.linkage_per_t.to_bits()
            && self.core.current_per_field.to_bits() == core.current_per_field.to_bits()
            && self.core.volume_m3.to_bits() == core.volume_m3.to_bits()
            && self.core.source().map(|s| s.material().card_identity())
                == core.source().map(|s| s.material().card_identity())
    }

    /// Adopt the original native material owner's freshly resolved conductor.
    /// Geometry and source pins stay fixed; only actual source state changes.
    /// This performs no material lookup, loss integration or history advance.
    pub fn with_equivalent_conductor(
        &self,
        conductor: &ResolvedConductor,
    ) -> Result<Self, ConductorError> {
        let original = &self.equivalent_conductor;
        if conductor.length_m().to_bits() != original.length_m().to_bits()
            || conductor.area_m2().to_bits() != original.area_m2().to_bits()
            || conductor.source_basis() != original.source_basis()
            || conductor.source_key() != original.source_key()
            || conductor.material().card_identity() != original.material().card_identity()
            || conductor.source_property().answer().receipt.selected
                != original.source_property().answer().receipt.selected
            || !conductor
                .source_property()
                .is_static_at_fixed_external_conditions()
        {
            return Err(invalid(
                "lamination shunt must retain original geometry, material and electrical source",
            ));
        }
        Ok(Self {
            equivalent_conductor: conductor.clone(),
            ..self.clone()
        })
    }
}
