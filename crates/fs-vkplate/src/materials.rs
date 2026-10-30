//! Source-bound material substitution for the existing plate/contact solver.
//!
//! Geometry, installed tension, supports, rim friction and snare-wire laws
//! remain independent declarations; a sourced snare-wire steel
//! (`with_snare_wire`) makes the wires coils of that steel. Source-bound
//! bending viscosity explicitly replaces T60 loss when selected. Drum heads
//! built from imported cards need that, the cards' film loss curve, or a
//! declared stated-unknown film loss; they never fall back to the
//! preset T60. A material swap never retunes pitch,
//! rescales geometry to preserve mass, or claims a new measured decay time.

use fs_matdb::{MatDbError, NormalizedMaterialCardPack, QueryPoint};
use fs_material::plane_stress_viscosity::{
    MaterialPlane, PLANE_STRESS_VISCOSITY_PROPERTIES, PlaneStressViscosity,
};
use fs_material::state_point::{
    IsotropicElasticStatePoint, MaterialPropertySelection, MaterialStatePointError,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_isotropic_elastic_state_point,
    resolve_material_state_point,
};
use fs_qty::{
    Dims, DynViscosity, Pressure, QuantitySpec,
    semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm},
};
use std::sync::Arc;

/// Aggregate import budget; the canonical material-pack decoder owns the wire format.
pub const MAX_MATERIAL_PACK_BYTES: usize = 16 * 1024 * 1024;

/// A drum built from imported head-film cards that make no loss claim and
/// declare none refuses with this; the recording-fitted T60 is not a fallback.
pub const FILM_LOSS_REFUSAL: &str = "head film loss: none is selected; select a sourced bending loss or the head cards' film loss curve (KitMaterials::with_film_loss_curve), or declare the film loss a STATED UNKNOWN entered as 0 (KitMaterials::with_stated_unknown_film_loss); the recording-fitted head T60 is not a fallback";

/// Carried by every drum built under that declaration, in its notes and
/// `info_json`. No edge or hoop loss law exists either; it is entered as 0 too.
pub const FILM_LOSS_STATED_UNKNOWN_RECEIPT: &str = "head film and edge loss: STATED UNKNOWN, entered as 0; omits the measured high-mode excess loss; T60s above the first two modes are upper bounds";

/// Directional in-plane moduli a film card may carry beside its isotropic
/// one; the isotropic head's error is their spread, stated in the receipt.
pub const IN_PLANE_MODULI: [&str; 2] = [
    "young_modulus_machine_direction",
    "young_modulus_transverse_direction",
];

/// Claims that make a striker card's tip transversely isotropic about its
/// grain (axis L, Young's modulus `young_modulus` = E_L): E_T/E_L, G_LT/E_L,
/// nu_LT and the transverse-plane Poisson ratio. A card carries all or none.
pub const TIP_ANISOTROPY: [&str; 4] = [
    "et_over_el",
    "glt_over_el",
    "nu_lt",
    "poisson_ratio_transverse_plane",
];

/// The stated limitation of the kit's Hertz contact, printed with every tip.
pub const HERTZ_HALF_SPACE_CAVEAT: &str = "the head side of the Hertz pair is the film as an elastic half-space; a tensioned head film (0.08 to 0.36 mm in this kit) is not one, and that error may exceed the tip bracket";

/// Delafargue and Ulm (2004) indentation modulus of a transversely isotropic
/// solid indented along its symmetry axis (here the grain):
/// `M3 = 2 sqrt((C11 C33 - C13^2) / C11 / (1/C44 + 2 / (sqrt(C11 C33) + C13)))`.
/// It reduces to E / (1 - nu^2) for an isotropic solid.
pub fn axial_indentation_modulus(
    e_l: f64,
    e_p: f64,
    g_lp: f64,
    nu_lp: f64,
    nu_p: f64,
) -> Result<f64, String> {
    // Compliances with axis 3 = L; the plane (1, 2) is isotropic.
    let s = [
        [1.0 / e_p, -nu_p / e_p, -nu_lp / e_l],
        [-nu_p / e_p, 1.0 / e_p, -nu_lp / e_l],
        [-nu_lp / e_l, -nu_lp / e_l, 1.0 / e_l],
    ];
    let minor =
        |r0: usize, r1: usize, c0: usize, c1: usize| s[r0][c0] * s[r1][c1] - s[r0][c1] * s[r1][c0];
    let det =
        s[0][0] * minor(1, 2, 1, 2) - s[0][1] * minor(1, 2, 0, 2) + s[0][2] * minor(1, 2, 0, 1);
    let (c11, c13, c33) = (
        minor(1, 2, 1, 2) / det,
        minor(0, 1, 1, 2) / det,
        minor(0, 1, 0, 1) / det,
    );
    let m = 2.0
        * ((c11 * c33 - c13 * c13) / c11 / (1.0 / g_lp + 2.0 / ((c11 * c33).sqrt() + c13))).sqrt();
    if !(det > 0.0 && c11 > 0.0 && c33 > 0.0 && g_lp > 0.0 && m.is_finite() && m > 0.0) {
        return Err(
            "striker tip: the transversely isotropic constants are not positive definite".into(),
        );
    }
    Ok(m)
}

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
    /// The striker card's TIP_ANISOTROPY values at the same state point.
    tip_anisotropy: Option<[f64; 4]>,
    /// The head cards' film loss curve (`with_film_loss_curve`).
    film_curve: Option<FilmCurve>,
    /// Snare-wire steel and its shear modulus (`with_snare_wire`).
    snare_wire: Option<crate::coil::CoilSteel>,
}

/// The film claim a film loss reads: the loss tangent as a curve over a typed
/// cyclic-frequency axis, at the head's temperature.
pub const FILM_LOSS_TANGENT: &str = "loss_tangent_tensile";

/// Carried in place of FILM_LOSS_STATED_UNKNOWN_RECEIPT by a drum whose film
/// loss comes from its cards: no edge or hoop loss law exists, so that loss
/// is still entered as 0.
pub const EDGE_LOSS_STATED_UNKNOWN_RECEIPT: &str = "head edge, hoop, shell and mount loss: STATED UNKNOWN, entered as 0; omits the measured high-mode excess loss, so T60s above the first two modes stay upper bounds";

/// Where the head cards' FILM_LOSS_TANGENT curve is read.
#[derive(Clone, Debug, PartialEq)]
struct FilmCurve {
    temperature_axis: String,
    frequency_axis: String,
    temperature_k: f64,
    /// Each head card's frequency band [Hz], from its claim's validity.
    band_hz: [(f64, f64); 2],
}

/// One head card's FILM_LOSS_TANGENT at `frequency_hz`, through the canonical
/// resolver (its validity and its linear interpolation between knots).
fn film_curve_tan_delta(
    pack: &NormalizedMaterialCardPack,
    curve: &FilmCurve,
    frequency_hz: f64,
) -> Result<f64, String> {
    let point = temperature_frequency_point(
        &curve.temperature_axis,
        curve.temperature_k,
        &curve.frequency_axis,
        frequency_hz,
    )?;
    let requirement = ScalarPropertyRequirement::try_new(
        FILM_LOSS_TANGENT,
        Dims::NONE,
        ScalarAdmissibility::StrictlyPositive,
    )
    .map_err(|e| e.to_string())?;
    resolve_material_state_point(
        pack.card(),
        &point,
        &[requirement],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .map_err(|e| {
        format!("film loss curve: the head card does not resolve {FILM_LOSS_TANGENT} at {frequency_hz} Hz: {e}")
    })?
    .property(FILM_LOSS_TANGENT)
    .map(|p| p.value_si())
    .ok_or_else(|| format!("film loss curve: {FILM_LOSS_TANGENT} did not resolve"))
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
            tip_anisotropy: None,
            film_curve: None,
            snare_wire: None,
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
        if self.film_loss_stated_unknown || self.film_curve.is_some() {
            return Err(
                "a sourced bending loss is exclusive with a stated-unknown film loss and a film loss curve"
                    .into(),
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
        let point = temperature_frequency_point(axis, temperature_k, frequency_axis, frequency_hz)?;
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
        if self.film_curve.is_some() {
            return Err("a film loss curve and a stated-unknown film loss are exclusive".into());
        }
        self.film_loss_stated_unknown = true;
        Ok(self)
    }
    pub const fn film_loss_is_stated_unknown(&self) -> bool {
        self.film_loss_stated_unknown
    }

    /// Declare the head films' loss hysteretic, with the loss tangent each
    /// head's card gives as its FILM_LOSS_TANGENT curve over `frequency_axis`
    /// (a typed cyclic frequency) at `temperature_k` on `temperature_axis`.
    /// The kit fits each head's generalized-Maxwell law to that curve alone
    /// (`film::fit`) and puts it on the head's bending operator
    /// (`Body::attach_film_maxwell`), so every grid mode carries the curve's
    /// tan delta at its own frequency times its bending share. The edge loss stays a
    /// stated unknown entered as zero (EDGE_LOSS_STATED_UNKNOWN_RECEIPT).
    /// A head card without exactly one such claim, or one that does not
    /// resolve inside its own band at this temperature, refuses.
    pub fn with_film_loss_curve(
        mut self,
        temperature_axis: &str,
        frequency_axis: &str,
        temperature_k: f64,
    ) -> Result<Self, String> {
        if self.bending_loss.is_some() {
            return Err("a sourced bending loss and a film loss curve are exclusive".into());
        }
        if self.film_loss_stated_unknown {
            return Err("a stated-unknown film loss and a film loss curve are exclusive".into());
        }
        if self.packs.len() < 2 {
            return Err("a film loss curve reads the head cards; these materials have none".into());
        }
        let mut band_hz = [(0.0, 0.0); 2];
        for (slot, pack) in band_hz.iter_mut().zip(&self.packs[..2]) {
            let claims = pack.card().claims_for(FILM_LOSS_TANGENT);
            let [(_, claim)] = claims.as_slice() else {
                return Err(format!(
                    "film loss curve: a head card must carry exactly one {FILM_LOSS_TANGENT} claim, not {}",
                    claims.len()
                ));
            };
            *slot = claim
                .validity
                .bound(frequency_axis)
                .filter(|&(lo, hi)| lo.is_finite() && hi.is_finite() && lo > 0.0 && hi > lo)
                .ok_or_else(|| {
                    format!(
                        "film loss curve: the head card's {FILM_LOSS_TANGENT} has no positive band on its {frequency_axis} axis"
                    )
                })?;
        }
        let curve = FilmCurve {
            temperature_axis: temperature_axis.into(),
            frequency_axis: frequency_axis.into(),
            temperature_k,
            band_hz,
        };
        // Each card answers at both ends of its band here, or the binding refuses.
        for (pack, (lo, hi)) in self.packs[..2].iter().zip(band_hz) {
            film_curve_tan_delta(pack, &curve, lo)?;
            film_curve_tan_delta(pack, &curve, hi)?;
        }
        self.film_curve = Some(curve);
        Ok(self)
    }

    pub const fn has_film_loss_curve(&self) -> bool {
        self.film_curve.is_some()
    }

    /// The batter (0) or resonant (1) head card's film loss-tangent band [Hz],
    /// an error without a curve.
    pub fn film_band_hz(&self, head: usize) -> Result<(f64, f64), String> {
        self.film_curve
            .as_ref()
            .map(|curve| curve.band_hz[head.min(1)])
            .ok_or_else(|| "film loss curve: none is declared".to_string())
    }

    /// The batter (0) or resonant (1) head card's film loss tangent at
    /// `frequency_hz`: `Ok(None)` outside its band, an error without a curve.
    pub fn film_tan_delta(&self, head: usize, frequency_hz: f64) -> Result<Option<f64>, String> {
        let curve = self
            .film_curve
            .as_ref()
            .ok_or("film loss curve: none is declared")?;
        if !frequency_hz.is_finite() {
            return Err(format!("film loss curve: frequency {frequency_hz} Hz"));
        }
        let (lo, hi) = curve.band_hz[head.min(1)];
        if frequency_hz < lo || frequency_hz > hi {
            return Ok(None);
        }
        film_curve_tan_delta(&self.packs[head.min(1)], curve, frequency_hz).map(Some)
    }

    /// Bind the snare wires' steel: the wire card (density, E, expansion, heat
    /// capacity, conductivity and its excess bending loss tangent) at
    /// `steel_point` and the shear card's modulus at `shear_point`, each at
    /// its own card's coordinates (`crate::coil::CoilSteel::resolve`). A drum
    /// with wires then builds each strand as a coil of its declared geometry.
    pub fn with_snare_wire(
        mut self,
        steel: &[u8],
        steel_point: &QueryPoint,
        shear: &[u8],
        shear_point: &QueryPoint,
    ) -> Result<Self, String> {
        if steel.len().saturating_add(shear.len()) > MAX_MATERIAL_PACK_BYTES {
            return Err("snare wire packs exceed the aggregate 16 MiB import budget".into());
        }
        let load = |bytes: &[u8]| {
            NormalizedMaterialCardPack::from_bytes(bytes)
                .map(Arc::new)
                .map_err(|e| e.to_string())
        };
        let (steel, shear) = (load(steel)?, load(shear)?);
        self.snare_wire = Some(crate::coil::CoilSteel::resolve(
            &steel,
            steel_point,
            &shear,
            shear_point,
        )?);
        self.packs.extend([steel, shear]);
        Ok(self)
    }

    /// The bound snare-wire steel, if any.
    pub const fn snare_wire(&self) -> Option<&crate::coil::CoilSteel> {
        self.snare_wire.as_ref()
    }

    /// A drum's film-loss receipt with a film loss curve, beside
    /// EDGE_LOSS_STATED_UNKNOWN_RECEIPT.
    pub fn film_loss_note(&self) -> Option<String> {
        self.film_curve.as_ref().map(|curve| {
            let [(b_lo, b_hi), (r_lo, r_hi)] = curve.band_hz;
            format!(
                "head film loss: each head's generalized-Maxwell law on its bending operator, fitted to its card's {FILM_LOSS_TANGENT} curve at {} K alone (batter card {b_lo}-{b_hi} Hz, resonant card {r_lo}-{r_hi} Hz), so every grid mode carries the tan delta at its own frequency times its bending share; tension and the Berger in-plane stretching stay elastic",
                curve.temperature_k
            )
        })
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
        // The tip constants change the contact law. A card carries all of them
        // or none; a partial set, or a claim present but unresolvable here
        // (out of its validity, say), refuses instead of leaving the tip
        // silently isotropic.
        let tip = TIP_ANISOTROPY
            .iter()
            .map(|name| {
                let requirement = ScalarPropertyRequirement::try_new(
                    *name,
                    Dims::NONE,
                    ScalarAdmissibility::Finite,
                )
                .map_err(|e| e.to_string())?;
                match resolve_material_state_point(
                    striker.card(),
                    point,
                    &[requirement],
                    MaterialPropertySelection::SingleClaimOnly,
                ) {
                    Ok(resolved) => resolved
                        .property(name)
                        .map(|p| Some(p.value_si()))
                        .ok_or_else(|| format!("striker tip {name}: not resolved")),
                    Err(MaterialStatePointError::Query {
                        source: MatDbError::UnknownProperty { .. },
                        ..
                    }) => Ok(None),
                    Err(e) => Err(format!("striker tip {name}: {e}")),
                }
            })
            .collect::<Result<Vec<_>, String>>()?;
        result.tip_anisotropy = match tip[..] {
            [Some(et), Some(glt), Some(nu_lt), Some(nu_p)] => {
                axial_indentation_modulus(1.0, et, glt, nu_lt, nu_p)?;
                Some([et, glt, nu_lt, nu_p])
            }
            [None, None, None, None] => None,
            _ => {
                return Err(format!(
                    "striker card carries only part of {TIP_ANISOTROPY:?}"
                ));
            }
        };
        result.packs = vec![primary, secondary, striker];
        Ok(result)
    }

    /// The striker tip's modulus in the Hertz pair: the axial indentation
    /// modulus when its card carries TIP_ANISOTROPY, else E / (1 - nu^2).
    pub fn striker_tip_modulus_pa(&self) -> f64 {
        let (e, nu) = (
            self.striker.young_modulus_pa(),
            self.striker.poisson_ratio(),
        );
        self.tip_anisotropy
            .and_then(|[et, glt, nu_lt, nu_p]| {
                axial_indentation_modulus(e, et * e, glt * e, nu_lt, nu_p).ok()
            })
            .unwrap_or(e / (1.0 - nu * nu))
    }

    /// One line per surface whose card carries directional moduli: the
    /// isotropic head's reduction and its spread, the size of its error.
    pub fn head_reduction_notes(&self) -> Vec<String> {
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

    /// The striker tip's receipt. A tip with TIP_ANISOTROPY prints the
    /// bracket between its axial indentation modulus and the isotropic-E_L
    /// end: E* on the batter film and the Hertz contact-time ratio
    /// (E*_a / E*_i)^(-2/5). Every tip prints HERTZ_HALF_SPACE_CAVEAT.
    pub fn striker_tip_note(&self) -> String {
        let film = self.primary.young_modulus_pa() / (1.0 - self.primary.poisson_ratio().powi(2));
        let e_star = |tip: f64| 1.0 / (1.0 / tip + 1.0 / film);
        let (e, nu) = (
            self.striker.young_modulus_pa(),
            self.striker.poisson_ratio(),
        );
        let isotropic = e / (1.0 - nu * nu);
        if self.tip_anisotropy.is_some() {
            let m = self.striker_tip_modulus_pa();
            format!(
                "striker tip: axial indentation modulus (Delafargue-Ulm 2004, transversely isotropic about the grain from E_L and the card's E_T/E_L, G_LT/E_L, nu_LT, nu_P) {:.3} GPa, E* on the batter film {:.3} GPa; isotropic-E_L end {:.3} GPa, E* {:.3} GPa; contact time {:.3} x the isotropic end; {HERTZ_HALF_SPACE_CAVEAT}",
                m * 1e-9,
                e_star(m) * 1e-9,
                isotropic * 1e-9,
                e_star(isotropic) * 1e-9,
                fs_math::det::pow(e_star(m) / e_star(isotropic), -0.4)
            )
        } else {
            format!(
                "striker tip: isotropic, E / (1 - nu^2) = {:.3} GPa, E* on the batter film {:.3} GPa; {HERTZ_HALF_SPACE_CAVEAT}",
                isotropic * 1e-9,
                e_star(isotropic) * 1e-9
            )
        }
    }

    /// The binding's stated reductions, for a drum's notes: the heads', then
    /// the striker tip's.
    pub fn isotropic_reduction_notes(&self) -> Vec<String> {
        let mut notes = self.head_reduction_notes();
        notes.push(self.striker_tip_note());
        notes
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

    /// Hertz law of the striker tip on the primary film; see
    /// HERTZ_HALF_SPACE_CAVEAT. An anisotropic tip enters through its axial
    /// indentation modulus M, passed as (M, 0) so its compliance is exactly
    /// 1 / M; an isotropic tip keeps its (E, nu) arithmetic unchanged.
    pub(crate) fn impact_law(&self) -> crate::Law {
        let tip = if self.tip_anisotropy.is_some() {
            (self.striker_tip_modulus_pa(), 0.0)
        } else {
            (
                self.striker.young_modulus_pa(),
                self.striker.poisson_ratio(),
            )
        };
        crate::kit::hertz_pair(
            self.striker_tip_radius_m,
            tip,
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
        // Fixed ASCII and numbers only: the axis names are not echoed.
        let curve = self.film_curve.as_ref().map_or_else(
            || "null".into(),
            |c| {
                let [(b_lo, b_hi), (r_lo, r_hi)] = c.band_hz;
                format!(
                    "{{\"property\":\"{FILM_LOSS_TANGENT}\",\"temperatureK\":{},\"bandHz\":[[{b_lo},{b_hi}],[{r_lo},{r_hi}]],\"edgeLoss\":\"{EDGE_LOSS_STATED_UNKNOWN_RECEIPT}\"}}",
                    c.temperature_k
                )
            },
        );
        format!(
            "{{\"primary\":{},\"secondary\":{},\"striker\":{},\"strikerVolumeM3\":{},\"strikerMassKg\":{},\"strikerTipRadiusM\":{},\"impactLossSM\":{},\"importedPacks\":{},\"bendingLoss\":{},\"filmLossStatedUnknown\":{},\"filmLossCurve\":{}}}",
            surface(&self.primary),
            surface(&self.secondary),
            surface(&self.striker),
            self.striker_volume_m3,
            self.striker_mass_kg(),
            self.striker_tip_radius_m,
            self.impact_loss_s_m,
            self.packs.len(),
            loss,
            film,
            curve
        )
    }
}

/// `temperature_point` plus one typed cyclic-frequency [Hz] axis: a source
/// coordinate (a claim valid over a band answers anywhere inside it), not a
/// frequency at which to freeze a dispersive curve.
pub fn temperature_frequency_point(
    axis: &str,
    temperature_k: f64,
    frequency_axis: &str,
    frequency_hz: f64,
) -> Result<QueryPoint, String> {
    temperature_point(axis, temperature_k)?
        .with_quantity(
            frequency_axis,
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::Frequency(FrequencyConvention::Cyclic),
                ValueForm::Static,
            )),
            frequency_hz,
        )
        .map_err(|e| e.to_string())
}

/// A query point on one typed absolute-temperature axis.
pub fn temperature_point(axis: &str, temperature_k: f64) -> Result<QueryPoint, String> {
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
