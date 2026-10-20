//! Source-bound material substitution for the existing plate/contact solver.
//!
//! Geometry, installed tension, T60 loss, supports, rim friction and snare-wire
//! laws remain independent declarations. A material swap never retunes pitch,
//! rescales geometry to preserve mass, or claims a new measured decay time.

use fs_matdb::{NormalizedMaterialCardPack, QueryPoint};
use fs_material::state_point::{
    IsotropicElasticStatePoint, MaterialPropertySelection, resolve_isotropic_elastic_state_point,
};
use fs_qty::{
    QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::sync::Arc;

/// Aggregate import budget; the canonical material-pack decoder owns the wire format.
pub const MAX_MATERIAL_PACK_BYTES: usize = 16 * 1024 * 1024;

/// The primary and secondary vibrating surfaces and one spherical impactor.
/// The complete resolved bundles retain material/card identities, query points,
/// selected claims and usage receipts. Imported packs additionally stay alive.
#[derive(Clone, Debug, PartialEq)]
pub struct KitMaterials {
    primary: IsotropicElasticStatePoint,
    secondary: IsotropicElasticStatePoint,
    striker: IsotropicElasticStatePoint,
    striker_volume_m3: f64,
    striker_tip_radius_m: f64,
    impact_loss_s_m: f64,
    packs: Vec<Arc<NormalizedMaterialCardPack>>,
}

impl KitMaterials {
    /// Bind already-resolved canonical material data. The striker's effective
    /// volume is a geometric/inertial declaration: its moving mass is rho V.
    /// Hertz contact uses its tip radius and the two resolved elastic moduli;
    /// Hunt-Crossley loss is separately declared, never inferred from density.
    pub fn try_new(
        primary: IsotropicElasticStatePoint,
        secondary: IsotropicElasticStatePoint,
        striker: IsotropicElasticStatePoint,
        striker_volume_m3: f64,
        striker_tip_radius_m: f64,
        impact_loss_s_m: f64,
    ) -> Result<Self, String> {
        for (name, value) in [
            ("striker effective volume", striker_volume_m3),
            ("striker tip radius", striker_tip_radius_m),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(format!("{name} must be finite and positive (SI)"));
            }
        }
        if !impact_loss_s_m.is_finite() || impact_loss_s_m < 0.0 {
            return Err("impact loss must be finite and nonnegative [s/m]".into());
        }
        let mass = striker.density_kg_m3() * striker_volume_m3;
        if !mass.is_finite() || mass <= 0.0 {
            return Err("material and effective striker volume give invalid mass".into());
        }
        let result = Self {
            primary,
            secondary,
            striker,
            striker_volume_m3,
            striker_tip_radius_m,
            impact_loss_s_m,
            packs: Vec::new(),
        };
        // Refuse a nonrepresentable contact before building any operators.
        let crate::Law::Power { k, .. } = result.impact_law();
        if !k.is_finite() || k <= 0.0 {
            return Err("material pair and tip radius give invalid Hertz stiffness".into());
        }
        Ok(result)
    }

    /// Decode ordinary canonical material-card packs and resolve at one exact
    /// point. An empty secondary pack explicitly reuses the primary material.
    /// Conflicting claims refuse under SingleClaimOnly; native callers can
    /// resolve pinned claims themselves and use `try_new`.
    pub fn from_card_packs(
        primary: &[u8],
        secondary: &[u8],
        striker: &[u8],
        point: &QueryPoint,
        striker_volume_m3: f64,
        striker_tip_radius_m: f64,
        impact_loss_s_m: f64,
    ) -> Result<Self, String> {
        if primary
            .len()
            .saturating_add(secondary.len())
            .saturating_add(striker.len())
            > MAX_MATERIAL_PACK_BYTES
        {
            return Err("kit material packs exceed the aggregate 16 MiB import budget".into());
        }
        let load = |bytes: &[u8]| {
            NormalizedMaterialCardPack::from_bytes(bytes)
                .map(Arc::new)
                .map_err(|e| e.to_string())
        };
        let primary = load(primary)?;
        let secondary = if secondary.is_empty() {
            Arc::clone(&primary)
        } else {
            load(secondary)?
        };
        let striker = load(striker)?;
        let resolve = |pack: &NormalizedMaterialCardPack| {
            resolve_isotropic_elastic_state_point(
                pack.card(),
                point,
                MaterialPropertySelection::SingleClaimOnly,
            )
            .map_err(|e| e.to_string())
        };
        let mut result = Self::try_new(
            resolve(&primary)?,
            resolve(&secondary)?,
            resolve(&striker)?,
            striker_volume_m3,
            striker_tip_radius_m,
            impact_loss_s_m,
        )?;
        result.packs = vec![primary, secondary, striker];
        Ok(result)
    }

    /// Browser import at a declared absolute material temperature. The axis
    /// name must match the source card; no alias or temperature extrapolation
    /// is invented. Native `from_card_packs` accepts a general typed query point.
    pub fn from_packs_at_temperature(
        primary: &[u8],
        secondary: &[u8],
        striker: &[u8],
        axis: &str,
        temperature_k: f64,
        striker_volume_m3: f64,
        striker_tip_radius_m: f64,
        impact_loss_s_m: f64,
    ) -> Result<Self, String> {
        if axis.trim().is_empty() || axis.len() > 128 {
            return Err(
                "material temperature axis must name a source coordinate (at most 128 bytes)"
                    .into(),
            );
        }
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err("absolute material temperature must be finite and positive [K]".into());
        }
        let point = QueryPoint::new()
            .with_quantity(
                axis,
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::AbsoluteTemperature,
                    ValueForm::Static,
                )),
                temperature_k,
            )
            .map_err(|e| e.to_string())?;
        Self::from_card_packs(
            primary,
            secondary,
            striker,
            &point,
            striker_volume_m3,
            striker_tip_radius_m,
            impact_loss_s_m,
        )
    }

    pub const fn primary(&self) -> &IsotropicElasticStatePoint {
        &self.primary
    }
    pub const fn secondary(&self) -> &IsotropicElasticStatePoint {
        &self.secondary
    }
    pub const fn striker(&self) -> &IsotropicElasticStatePoint {
        &self.striker
    }
    pub fn source_packs(&self) -> &[Arc<NormalizedMaterialCardPack>] {
        &self.packs
    }
    pub fn striker_mass_kg(&self) -> f64 {
        self.striker.density_kg_m3() * self.striker_volume_m3
    }
    pub const fn striker_tip_radius_m(&self) -> f64 {
        self.striker_tip_radius_m
    }

    pub(crate) fn impact_law(&self) -> crate::Law {
        crate::kit::hertz_pair(
            self.striker_tip_radius_m,
            (
                self.striker.young_modulus_pa(),
                self.striker.poisson_ratio(),
            ),
            (
                self.primary.young_modulus_pa(),
                self.primary.poisson_ratio(),
            ),
            self.impact_loss_s_m,
        )
    }

    /// Compact physical metadata; the complete usage receipts remain accessible
    /// through the three resolved bundles. Hash text contains only hex digits.
    pub fn info_json(&self) -> String {
        let surface = |s: &IsotropicElasticStatePoint| {
            format!(
                "{{\"stateIdentity\":\"{}\",\"cardIdentity\":\"{}\",\"densityKgM3\":{},\"youngPa\":{},\"poisson\":{}}}",
                s.resolved().identity(),
                s.resolved().card_identity(),
                s.density_kg_m3(),
                s.young_modulus_pa(),
                s.poisson_ratio()
            )
        };
        format!(
            "{{\"primary\":{},\"secondary\":{},\"striker\":{},\"strikerVolumeM3\":{},\"strikerMassKg\":{},\"strikerTipRadiusM\":{},\"impactLossSM\":{},\"importedPacks\":{}}}",
            surface(&self.primary),
            surface(&self.secondary),
            surface(&self.striker),
            self.striker_volume_m3,
            self.striker_mass_kg(),
            self.striker_tip_radius_m,
            self.impact_loss_s_m,
            self.packs.len()
        )
    }
}
