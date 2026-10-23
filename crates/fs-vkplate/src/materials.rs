//! Source-bound material substitution for the existing plate/contact solver.
//!
//! Geometry, installed tension, supports, rim friction and snare-wire laws
//! remain independent declarations. Source-bound bending viscosity explicitly
//! replaces T60 loss when selected. Drum heads built from imported cards need
//! that, or a declared stated-unknown film loss; they never fall back to the
//! preset T60. A material swap never retunes pitch,
//! rescales geometry to preserve mass, or claims a new measured decay time.

use fs_matdb::{NormalizedMaterialCardPack, QueryPoint};
use fs_material::plane_stress_viscosity::{
    MaterialPlane, PLANE_STRESS_VISCOSITY_PROPERTIES, PlaneStressViscosity,
};
use fs_material::state_point::{
    IsotropicElasticStatePoint, MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_isotropic_elastic_state_point, resolve_material_state_point,
};
use fs_qty::{
    DynViscosity, Pressure, QuantitySpec,
    semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm},
};
use std::sync::Arc;

/// Aggregate import budget; the canonical material-pack decoder owns the wire format.
pub const MAX_MATERIAL_PACK_BYTES: usize = 16 * 1024 * 1024;

/// A drum built from imported head-film cards that make no loss claim and
/// declare none refuses with this; the recording-fitted T60 is not a fallback.
pub const FILM_LOSS_REFUSAL: &str = "head film loss: the imported head-film cards make no loss claim; select a sourced bending loss, or declare the film loss a STATED UNKNOWN entered as 0 (KitMaterials::with_stated_unknown_film_loss); the recording-fitted head T60 is not a fallback";

/// Carried by every drum built under that declaration, in its notes and
/// `info_json`. No edge or hoop loss law exists either; it is entered as 0 too.
pub const FILM_LOSS_STATED_UNKNOWN_RECEIPT: &str = "head film and edge loss: STATED UNKNOWN, entered as 0; omits the measured high-mode excess loss; T60s above the first two modes are upper bounds";

/// Directional in-plane moduli a film card may carry beside its isotropic
/// one; the isotropic head's error is their spread, stated in the receipt.
pub const IN_PLANE_MODULI: [&str; 2] = [
    "young_modulus_machine_direction",
    "young_modulus_transverse_direction",
];

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
    bending_loss: Option<[(PlaneStressViscosity, f64); 2]>,
    film_loss_stated_unknown: bool,
    /// Machine- and transverse-direction moduli [Pa] of each surface's card
    /// at the same state point, when the card carries both.
    in_plane_moduli: [Option<(f64, f64)>; 2],
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
            bending_loss: None,
            film_loss_stated_unknown: false,
            in_plane_moduli: [None, None],
        };
        // Refuse a nonrepresentable contact before building any operators.
        let crate::Law::Power { k, .. } = result.impact_law();
        if !k.is_finite() || k <= 0.0 {
            return Err("material pair and tip radius give invalid Hertz stiffness".into());
        }
        Ok(result)
    }

    /// Select actual reduced plane-stress curvature-rate dissipation on both
    /// surfaces. Each law must describe the same card, coordinates and axis
    /// meanings as that surface's elastic source. Orientations are relative to
    /// the plate grid; no viscous tensor is inferred from elastic Poisson ratio.
    pub fn with_bending_viscosity(
        mut self,
        primary: PlaneStressViscosity,
        secondary: PlaneStressViscosity,
        angles_rad: [f64; 2],
    ) -> Result<Self, String> {
        if self.film_loss_stated_unknown {
            return Err(
                "a sourced bending loss and a stated-unknown film loss are exclusive".into(),
            );
        }
        for (elastic, viscous) in [(&self.primary, &primary), (&self.secondary, &secondary)] {
            let (a, b) = (elastic.resolved(), viscous.source());
            let same_point = a.query_point().len() == b.query_point().len()
                && a.query_point()
                    .iter()
                    .zip(b.query_point())
                    .all(|((an, av), (bn, bv))| an == bn && av.to_bits() == bv.to_bits());
            if a.card_identity() != b.card_identity()
                || !same_point
                || a.axis_quantities() != b.axis_quantities()
            {
                return Err("bending viscosity and elasticity must retain the same material card and exact state point".into());
            }
            for property in a.properties() {
                if !property.is_constant_at_fixed_temperature() {
                    return Err("bending-loss elasticity and density must be constants or temperature-only curves at fixed positive absolute temperature".into());
                }
                property
                    .angular_frequency_band_rad_s()
                    .map_err(|e| e.to_string())?;
            }
        }
        primary
            .rotated_matrix_pa_s(angles_rad[0])
            .map_err(|e| e.to_string())?;
        secondary
            .rotated_matrix_pa_s(angles_rad[1])
            .map_err(|e| e.to_string())?;
        self.bending_loss = Some([(primary, angles_rad[0]), (secondary, angles_rad[1])]);
        Ok(self)
    }

    /// Browser data route for Kelvin--Voigt bending. Frequency is a typed Hz
    /// source coordinate, not a frequency at which to freeze a dispersive curve.
    /// The actual internal clock is checked by the kit builder after admission.
    pub fn from_packs_with_bending_loss(
        primary: &[u8],
        secondary: &[u8],
        striker: &[u8],
        axis: &str,
        temperature_k: f64,
        striker_volume_m3: f64,
        striker_tip_radius_m: f64,
        impact_loss_s_m: f64,
        frequency_axis: &str,
        frequency_hz: f64,
        angles_rad: [f64; 2],
    ) -> Result<Self, String> {
        if frequency_axis.trim().is_empty()
            || frequency_axis.len() > 128
            || frequency_axis == axis
            || !frequency_hz.is_finite()
            || frequency_hz <= 0.0
        {
            return Err("bending viscosity requires a distinct source frequency axis and positive finite frequency [Hz]".into());
        }
        let point = temperature_point(axis, temperature_k)?
            .with_quantity(
                frequency_axis,
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::Frequency(FrequencyConvention::Cyclic),
                    ValueForm::Static,
                )),
                frequency_hz,
            )
            .map_err(|e| e.to_string())?;
        let materials = Self::from_card_packs(
            primary,
            secondary,
            striker,
            &point,
            striker_volume_m3,
            striker_tip_radius_m,
            impact_loss_s_m,
        )?;
        let requirements = PLANE_STRESS_VISCOSITY_PROPERTIES
            .iter()
            .map(|name| {
                ScalarPropertyRequirement::try_new(
                    *name,
                    DynViscosity::DIMS,
                    ScalarAdmissibility::Finite,
                )
                .map_err(|e| e.to_string())
            })
            .collect::<Result<Vec<_>, _>>()?;
        let resolve = |pack: &NormalizedMaterialCardPack| -> Result<PlaneStressViscosity, String> {
            let source = resolve_material_state_point(
                pack.card(),
                &point,
                &requirements,
                MaterialPropertySelection::SingleClaimOnly,
            )
            .map_err(|e| e.to_string())?;
            PlaneStressViscosity::from_material(
                source,
                MaterialPlane::Axes12,
                PLANE_STRESS_VISCOSITY_PROPERTIES,
            )
            .map_err(|e| e.to_string())
        };
        let primary = resolve(&materials.packs[0])?;
        let secondary = resolve(&materials.packs[1])?;
        materials.with_bending_viscosity(primary, secondary, angles_rad)
    }

    pub fn primary_bending_loss(&self) -> Option<(&PlaneStressViscosity, f64)> {
        self.bending_loss.as_ref().map(|v| (&v[0].0, v[0].1))
    }
    pub fn secondary_bending_loss(&self) -> Option<(&PlaneStressViscosity, f64)> {
        self.bending_loss.as_ref().map(|v| (&v[1].0, v[1].1))
    }

    /// Declare that the head films' loss is unknown and enter it as zero; the
    /// drum then carries [`FILM_LOSS_STATED_UNKNOWN_RECEIPT`]. Without this or
    /// a sourced bending loss, a drum built from these cards refuses with
    /// [`FILM_LOSS_REFUSAL`]. Cymbal and hi-hat losses are not affected.
    pub fn with_stated_unknown_film_loss(mut self) -> Result<Self, String> {
        if self.bending_loss.is_some() {
            return Err(
                "a sourced bending loss and a stated-unknown film loss are exclusive".into(),
            );
        }
        self.film_loss_stated_unknown = true;
        Ok(self)
    }
    pub const fn film_loss_is_stated_unknown(&self) -> bool {
        self.film_loss_stated_unknown
    }

    /// Broadband time stepping cannot use a source characterized only at its
    /// query frequency. Require both selected surfaces' coefficient domains to
    /// cover DC through the selected INTERNAL Nyquist frequency, including
    /// elastic coefficients; narrower or incompatible sources refuse.
    pub(crate) fn validate_bending_loss_clock(&self, rate: f64) -> Result<(), String> {
        let Some(loss) = &self.bending_loss else {
            return Ok(());
        };
        let omega = std::f64::consts::PI * rate;
        for (elastic, viscous) in [(&self.primary, &loss[0].0), (&self.secondary, &loss[1].0)] {
            let mut band = viscous.omega_band_rad_s();
            for property in elastic.resolved().properties() {
                if let Some((lo, hi)) = property
                    .angular_frequency_band_rad_s()
                    .map_err(|e| e.to_string())?
                {
                    band = (band.0.max(lo), band.1.min(hi));
                }
            }
            if !omega.is_finite() || rate <= 0.0 || band.0 > 0.0 || band.1 < omega {
                return Err(format!(
                    "bending-loss sources must cover 0 through internal Nyquist {} Hz; admitted intersection is [{}, {}] rad/s",
                    rate / 2.0,
                    band.0,
                    band.1
                ));
            }
        }
        Ok(())
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
        // Optional: a card without both directional claims states none.
        let directional = |pack: &NormalizedMaterialCardPack| -> Option<(f64, f64)> {
            let requirements = IN_PLANE_MODULI
                .iter()
                .map(|name| {
                    ScalarPropertyRequirement::try_new(
                        *name,
                        Pressure::DIMS,
                        ScalarAdmissibility::StrictlyPositive,
                    )
                })
                .collect::<Result<Vec<_>, _>>()
                .ok()?;
            let resolved = resolve_material_state_point(
                pack.card(),
                point,
                &requirements,
                MaterialPropertySelection::SingleClaimOnly,
            )
            .ok()?;
            Some((
                resolved.property(IN_PLANE_MODULI[0])?.value_si(),
                resolved.property(IN_PLANE_MODULI[1])?.value_si(),
            ))
        };
        result.in_plane_moduli = [directional(&primary), directional(&secondary)];
        result.packs = vec![primary, secondary, striker];
        Ok(result)
    }

    /// One line per surface whose card carries directional moduli: the
    /// isotropic head's reduction and its spread, the size of its error.
    pub fn isotropic_reduction_notes(&self) -> Vec<String> {
        [("batter", &self.primary, self.in_plane_moduli[0]), ("resonant", &self.secondary, self.in_plane_moduli[1])]
            .into_iter()
            .filter_map(|(label, surface, moduli)| {
                let (md, td) = moduli?;
                let (e, mean) = (surface.young_modulus_pa(), 0.5 * (md + td));
                let rule = if (e - mean).abs() <= 1e-9 * mean { "mean(MD, TD)" } else { "the card's isotropic claim" };
                Some(format!(
                    "isotropic {label} head: E = {rule} = {:.1} MPa; MD/TD = {:.1} / {:.1} MPa (+/- {:.1} %)",
                    e * 1e-6,
                    md * 1e-6,
                    td * 1e-6,
                    50.0 * (md - td).abs() / mean
                ))
            })
            .collect()
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
        let point = temperature_point(axis, temperature_k)?;
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
        let loss = self.bending_loss.as_ref().map_or_else(|| "null".into(), |loss| {
            let info = |(material, angle): &(PlaneStressViscosity, f64)| {
                let matrix = material.matrix_pa_s().iter().map(f64::to_string).collect::<Vec<_>>().join(",");
                let (lo, hi) = material.omega_band_rad_s();
                format!("{{\"stateIdentity\":\"{}\",\"matrixPaS\":[{}],\"omegaBandRadS\":[{},{}],\"angleRad\":{}}}", material.source().identity(), matrix, lo, hi, angle)
            };
            format!("{{\"model\":\"kelvin-voigt-plane-stress\",\"primary\":{},\"secondary\":{}}}", info(&loss[0]), info(&loss[1]))
        });
        // The receipt is a fixed ASCII constant without quotes or backslashes.
        let film = if self.film_loss_stated_unknown {
            format!("\"{FILM_LOSS_STATED_UNKNOWN_RECEIPT}\"")
        } else {
            "null".into()
        };
        format!(
            "{{\"primary\":{},\"secondary\":{},\"striker\":{},\"strikerVolumeM3\":{},\"strikerMassKg\":{},\"strikerTipRadiusM\":{},\"impactLossSM\":{},\"importedPacks\":{},\"bendingLoss\":{},\"filmLossStatedUnknown\":{}}}",
            surface(&self.primary),
            surface(&self.secondary),
            surface(&self.striker),
            self.striker_volume_m3,
            self.striker_mass_kg(),
            self.striker_tip_radius_m,
            self.impact_loss_s_m,
            self.packs.len(),
            loss,
            film
        )
    }
}

fn temperature_point(axis: &str, temperature_k: f64) -> Result<QueryPoint, String> {
    if axis.trim().is_empty() || axis.len() > 128 {
        return Err(
            "material temperature axis must name a source coordinate (at most 128 bytes)".into(),
        );
    }
    if !temperature_k.is_finite() || temperature_k <= 0.0 {
        return Err("absolute material temperature must be finite and positive [K]".into());
    }
    QueryPoint::new()
        .with_quantity(
            axis,
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            temperature_k,
        )
        .map_err(|e| e.to_string())
}
