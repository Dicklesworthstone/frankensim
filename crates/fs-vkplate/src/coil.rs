//! Coiled steel snare wires from sourced data, on the card route (ruling 7373).
//!
//! A snare wire is an open helical coil of round steel wire, pinned at both ends to the end
//! plates. Its bending stiffness, mass and loss here come from a declared geometry, two material
//! cards and the gas record; nothing is fitted to a recording.
//!
//! Geometry ([`CANOPUS_VINTAGE`]): the maker's Japanese page gives wire diameter d 0.5 mm, coil
//! OUTER diameter 1.35 mm and pitch p 3.4 mm. Spring mechanics uses the mean diameter
//! D = 1.35 - 0.5 = 0.85 mm, so the helix angle is a = atan(p / (pi D)) (51.9 degrees).
//!
//! Bending of the coil: a moment about an axis across the coil loads the wire, over a turn, in
//! bending and in torsion. With the wire's EI and GJ and the helix angle a, the turn-averaged
//! compliances give B = 2 sin(a) / (cos^2(a)/GJ + (1 + sin^2(a))/EI), which for a close coil
//! (a -> 0) is 2 p EI GJ / (pi D (EI + GJ)). The coil's mass per unit length is the wire's over
//! sin(a). The wire-bending share of the coil's bending energy is the bending compliance's share,
//! (1 + sin^2 a)/EI over the sum.
//!
//! Loss of a mode at omega whose potential energy is the fraction `share` coil bending (the rest
//! is strainer tension):
//! - the wire's thermoelastic (Zener) bending loss, computed by fs-material from the steel card's
//!   E, density, expansion, heat capacity and conductivity at the wire's radius, plus the card's
//!   measured excess bending loss tangent inside its source band only (a lower bound outside it),
//!   both on the wire-bending share of `share`. Wire torsion is pure shear: it has no
//!   thermoelastic loss, and its anelastic loss is unmeasured and omitted (a lower bound);
//! - the drag of the gas on the moving wire: each wire element moves partly across its own axis
//!   ([`air_loss::oscillating_cylinder`], Stuart's exact law) and partly along it
//!   ([`air_loss::sliding_cylinder`]). For the coil moving across its axis the turn averages of
//!   the squared velocity components are 1 - cos^2(a)/2 across the wire and cos^2(a)/2 along it.
//!   Linear Stokes flow: valid while the wire's amplitude is small against its radius. A rattling
//!   wire exceeds that, so this term is a lower bound. The gas's added mass, about 1.5e-4 of the
//!   steel's, is not applied.
//!
//! Contact with the head ([`turn_contact`]): each turn presses the film as an elliptical Hertz
//! contact of the wire's surface, curved along the helix, its stiffness computed from the two
//! cards. It holds while the contact is small against the film: up to about 0.1 N per turn.
//!
//! Not sourced, and stated in every coil's notes: the strand count, spacing, length, tension,
//! preload and contact loss (see `kit::Wires`). The ends are pinned with no strand-end or cord
//! damping, so the lowest coil modes keep their long material decay.

use crate::GasState;
use crate::air_loss::{oscillating_cylinder, sliding_cylinder};
use fs_matdb::{NormalizedMaterialCardPack, QueryPoint};
use fs_material::section::{
    ConcentricElasticSection, ConcentricMaterialLayer, loss_angle::ConcentricLossAngleSection,
    ted::ConcentricThermoelasticBending,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY, ScalarAdmissibility,
    ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
    YOUNG_MODULUS_PROPERTY, resolve_material_state_point,
};
use fs_qty::{Density, Dims, Pressure};
use std::f64::consts::PI;

/// The steel card's measured excess bending loss tangent (beyond thermoelasticity).
pub const EXCESS_LOSS_PROPERTY: &str = "bending_loss_tangent";
/// The shear card's modulus.
pub const SHEAR_MODULUS_PROPERTY: &str = "shear_modulus";

/// A declared coil: wire diameter, coil outer diameter and pitch [m], and their source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoilGeometry {
    pub wire_diameter_m: f64,
    pub outer_diameter_m: f64,
    pub pitch_m: f64,
    pub source: &'static str,
}

/// Canopus Vintage Snare Wire, the maker's replica of the 1960s Slingerland wire. Its Japanese
/// page states 線径0.5mm、線外径1.35mm、ピッチ3.4mm: wire diameter, coil outer diameter, pitch.
/// The English page (canopusdrums.com/en/product/vintage-snarewire) calls the 1.35 mm a "spiral
/// diameter" without saying which; 線外径 is the outer one.
pub const CANOPUS_VINTAGE: CoilGeometry = CoilGeometry {
    wire_diameter_m: 0.5e-3,
    outer_diameter_m: 1.35e-3,
    pitch_m: 3.4e-3,
    source: "Canopus Vintage Snare Wire, canopusdrums.com/jp/product/vintage-snare-wire/",
};

impl CoilGeometry {
    /// Mean coil diameter D = outer - wire [m].
    #[must_use]
    pub fn mean_diameter_m(&self) -> f64 {
        self.outer_diameter_m - self.wire_diameter_m
    }
    /// Helix angle atan(p / (pi D)) [rad].
    #[must_use]
    pub fn helix_angle_rad(&self) -> f64 {
        (self.pitch_m / (PI * self.mean_diameter_m())).atan()
    }
    /// The principal radii of the wire's surface where a turn rests on a flat, outermost from
    /// the coil's axis [m]: across the wire its own radius d/2; along it the helix's radius of
    /// curvature plus d/2. A helix of radius R = D/2 and lead p/(2 pi) curves at
    /// (R^2 + (p/(2 pi))^2) / R, which for a close coil is R.
    #[must_use]
    pub fn turn_surface_radii_m(&self) -> (f64, f64) {
        let (r, lead, a) = (
            0.5 * self.mean_diameter_m(),
            self.pitch_m / (2.0 * PI),
            0.5 * self.wire_diameter_m,
        );
        (a, (r * r + lead * lead) / r + a)
    }
}

/// The steel's Poisson ratio in the wire's Hertz contact: ASSUMED, neither steel card states one.
/// The film's compliance is about 40 times the steel's, so it moves K by under 0.1 %.
pub const STEEL_POISSON_RATIO_ASSUMED: f64 = 0.3;

/// The wire's Hertz contact on a head film (ruling 7520 L12 (a)): the elliptical contact
/// (`crate::hertz`) of a turn's surface (`CoilGeometry::turn_surface_radii_m`) on the flat
/// film, steel E from the wire card with `STEEL_POISSON_RATIO_ASSUMED` and the film card's E and
/// nu, both half-spaces. For `CANOPUS_VINTAGE` (d 0.5 mm, D 0.85 mm, p 3.4 mm: a 51.9 degree
/// helix) the surface's radii are d/2 = 0.25 mm across the wire and, along it,
/// (R^2 + (p/(2 pi))^2)/R + d/2 = 1.114 + 0.25 = 1.364 mm, so on the 4.903 GPa PET film K is
/// 1.91e8 N/m^1.5 per turn; the wire's radius alone, a point contact, bounds it below at 1.18e8.
/// A wire node of length `node_m` carries `node_m / pitch` turns at one approach, so its
/// stiffness is that many turns'. The film is a half-space only while the contact ellipse is
/// small against its thickness: at the 0.85 mN preload per turn the semi-axes are 2.2 x 6.8 um
/// on 76.2 um of film, in range; near 0.1 N per turn they are 11 x 33 um, and above that load
/// the law is outside its range. Returns the node's law with Hunt-Crossley loss `loss_s_m` and
/// a receipt that states the ellipse at the preload `preload_n_m` [N/m] and at 0.1 N against
/// the film's thickness `film_m`.
pub fn turn_contact(
    geometry: &CoilGeometry,
    steel: &CoilSteel,
    film: (f64, f64),
    film_m: f64,
    node_m: f64,
    preload_n_m: f64,
    loss_s_m: f64,
) -> Result<(crate::Law, String), String> {
    let steel_e = steel.young_modulus_pa();
    let wire = (steel_e, STEEL_POISSON_RATIO_ASSUMED);
    let (across, along) = geometry.turn_surface_radii_m();
    let turn = crate::hertz::EllipticHertz::isotropic(wire, film, across, along)?;
    let point = crate::hertz::EllipticHertz::isotropic(wire, film, across, across)?;
    let turns = node_m / geometry.pitch_m;
    let k = turn.k() * turns;
    if !(k.is_finite() && k > 0.0 && film_m > 0.0 && preload_n_m >= 0.0 && loss_s_m >= 0.0) {
        return Err(format!(
            "snare wire contact: unrepresentable K {k} N/m^1.5, film {film_m} m, preload {preload_n_m} N/m or loss {loss_s_m} s/m"
        ));
    }
    let e_star = 1.0 / ((1.0 - wire.1 * wire.1) / wire.0 + (1.0 - film.1 * film.1) / film.0);
    let preload = preload_n_m * geometry.pitch_m;
    let (pa, pb) = turn.semi_axes_m(preload);
    let (ha, hb) = turn.semi_axes_m(0.1);
    let note = format!(
        "snare wire contact: computed from the cards (Hertz, elliptical, helix): each turn's surface (radius {:.3} mm across the wire, {:.3} mm along its {:.1} deg helix) on the resonant film; steel E {:.4e} Pa (wire card; nu {STEEL_POISSON_RATIO_ASSUMED} assumed), film E {:.4e} Pa, nu {} (film card), E* {:.4e} Pa; K {:.4e} N/m^1.5 per turn (b/a {:.2}; the wire's radius alone, a point contact, gives {:.4e}), {:.3} turns per {:.1} mm node, {:.4e} per node; a half-space while the ellipse is small against the {:.0} um film: semi-axes {:.1} x {:.1} um at the {:.2} mN preload per turn, in range, and {:.0} x {:.0} um at 0.1 N per turn, above which the law is outside its range; Hunt-Crossley loss {loss_s_m} s/m declared",
        across * 1e3,
        along * 1e3,
        geometry.helix_angle_rad().to_degrees(),
        steel_e,
        film.0,
        film.1,
        e_star,
        turn.k(),
        turn.aspect(),
        point.k(),
        turns,
        node_m * 1e3,
        k,
        film_m * 1e6,
        pa * 1e6,
        pb * 1e6,
        preload * 1e3,
        ha * 1e6,
        hb * 1e6,
    );
    Ok((
        crate::Law::Power {
            k,
            alpha: 1.5,
            mu: loss_s_m,
        },
        note,
    ))
}

/// The coil's steel: one resolved state of the wire card (density, E, the three thermal
/// coefficients and its excess bending loss tangent) and the shear modulus from a second card
/// at that card's own state point.
#[derive(Clone, Debug, PartialEq)]
pub struct CoilSteel {
    steel: ResolvedMaterialStatePoint,
    shear: ResolvedMaterialStatePoint,
}

fn requirement(
    name: &str,
    dims: Dims,
    admissibility: ScalarAdmissibility,
) -> Result<ScalarPropertyRequirement, String> {
    ScalarPropertyRequirement::try_new(name, dims, admissibility).map_err(|e| e.to_string())
}

impl CoilSteel {
    /// Resolve the wire card at `steel_point` and the shear card at `shear_point`. A missing,
    /// ambiguous or out-of-domain property refuses, the excess loss tangent included: absent
    /// loss data is not a lossless wire.
    pub fn resolve(
        steel: &NormalizedMaterialCardPack,
        steel_point: &QueryPoint,
        shear: &NormalizedMaterialCardPack,
        shear_point: &QueryPoint,
    ) -> Result<Self, String> {
        let needs = [
            requirement(
                DENSITY_PROPERTY,
                Density::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )?,
            requirement(
                YOUNG_MODULUS_PROPERTY,
                Pressure::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )?,
            requirement(
                LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
                INVERSE_TEMPERATURE_DIMS,
                ScalarAdmissibility::Finite,
            )?,
            requirement(
                SPECIFIC_HEAT_CAPACITY_PROPERTY,
                SPECIFIC_HEAT_CAPACITY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )?,
            requirement(
                THERMAL_CONDUCTIVITY_PROPERTY,
                THERMAL_CONDUCTIVITY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )?,
            requirement(
                EXCESS_LOSS_PROPERTY,
                Dims::NONE,
                ScalarAdmissibility::Finite,
            )?,
        ];
        let steel = resolve_material_state_point(
            steel.card(),
            steel_point,
            &needs,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .map_err(|e| format!("snare wire steel: {e}"))?;
        let shear = resolve_material_state_point(
            shear.card(),
            shear_point,
            &[requirement(
                SHEAR_MODULUS_PROPERTY,
                Pressure::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )?],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .map_err(|e| format!("snare wire shear modulus: {e}"))?;
        Ok(Self { steel, shear })
    }

    /// The wire card's resolved state.
    #[must_use]
    pub const fn steel(&self) -> &ResolvedMaterialStatePoint {
        &self.steel
    }
    /// The shear card's resolved state.
    #[must_use]
    pub const fn shear(&self) -> &ResolvedMaterialStatePoint {
        &self.shear
    }
    /// The wire card's Young's modulus [Pa].
    #[must_use]
    pub fn young_modulus_pa(&self) -> f64 {
        self.steel
            .property(YOUNG_MODULUS_PROPERTY)
            .map_or(f64::NAN, |p| p.value_si())
    }
    fn shear_modulus_pa(&self) -> f64 {
        self.shear
            .property(SHEAR_MODULUS_PROPERTY)
            .map_or(f64::NAN, |p| p.value_si())
    }
}

/// One coiled wire of a declared geometry and steel, in a declared gas.
#[derive(Clone, Debug)]
pub struct CoilWire {
    geometry: CoilGeometry,
    sin: f64,
    cos2: f64,
    bending_n_m2: f64,
    mass_per_m: f64,
    wire_bending_share: f64,
    torsion_share: f64,
    thermoelastic: ConcentricThermoelasticBending,
    excess: ConcentricLossAngleSection,
    gas: GasState,
    steel_note: String,
}

impl CoilWire {
    pub fn new(geometry: &CoilGeometry, steel: &CoilSteel, gas: &GasState) -> Result<Self, String> {
        let (d, dm, p) = (
            geometry.wire_diameter_m,
            geometry.mean_diameter_m(),
            geometry.pitch_m,
        );
        if !(d.is_finite() && d > 0.0 && dm.is_finite() && dm > 0.0 && p.is_finite() && p > 0.0) {
            return Err(format!(
                "snare wire coil needs positive wire diameter, mean diameter (outer - wire) and pitch, got {d} m, {dm} m, {p} m"
            ));
        }
        let fail = |e: MaterialStatePointError| format!("snare wire section: {e}");
        let section = ConcentricElasticSection::try_new(
            0.0,
            vec![ConcentricMaterialLayer {
                outer_radius_m: d / 2.0,
                material: steel.steel.clone(),
            }],
        )
        .map_err(fail)?;
        let excess = section
            .clone()
            .with_loss_angle_bending(EXCESS_LOSS_PROPERTY)
            .map_err(fail)?;
        let thermoelastic = section.clone().with_thermoelastic_bending().map_err(fail)?;
        let g = steel.shear_modulus_pa();
        let ei = section.bending_stiffness_n_m2();
        // Saint-Venant torsion of a solid circle: J = 2 I
        let gj = g * 2.0 * section.second_moment_m4();
        let alpha = geometry.helix_angle_rad();
        let (sin, cos2) = (alpha.sin(), alpha.cos().powi(2));
        let (bend, twist) = ((1.0 + sin * sin) / ei, cos2 / gj);
        let bending_n_m2 = 2.0 * sin / (bend + twist);
        let mass_per_m = section.line_density_kg_m() / sin;
        if !(bending_n_m2.is_finite()
            && bending_n_m2 > 0.0
            && mass_per_m.is_finite()
            && mass_per_m > 0.0)
        {
            return Err(format!(
                "snare wire coil: unrepresentable bending {bending_n_m2} N m^2 or mass {mass_per_m} kg/m"
            ));
        }
        let law = thermoelastic.law();
        let (lo, hi) = excess.omega_band_rad_s();
        let excess_note = format!(
            "excess bending loss tangent {:.3e} inside {:.0}-{:.0} Hz only (thermoelastic alone outside it, a lower bound)",
            excess.hysteretic_stiffness_n_m2() / ei,
            lo / (2.0 * PI),
            hi / (2.0 * PI)
        );
        let steel_note = format!(
            "coil steel: state {} of card {} at {:?}: rho {} kg/m^3, E {:.4e} Pa, thermoelastic strength {:.4e} at {} K; {excess_note}; shear modulus {:.4e} Pa from state {} of card {} at {:?}",
            steel.steel.identity(),
            steel.steel.card_identity(),
            steel.steel.query_point(),
            law.rho,
            law.e,
            law.relaxation_strength(),
            law.t0,
            g,
            steel.shear.identity(),
            steel.shear.card_identity(),
            steel.shear.query_point(),
        );
        Ok(Self {
            geometry: *geometry,
            sin,
            cos2,
            bending_n_m2,
            mass_per_m,
            wire_bending_share: bend / (bend + twist),
            torsion_share: twist / (bend + twist),
            thermoelastic,
            excess,
            gas: *gas,
            steel_note,
        })
    }

    /// Coil bending stiffness B [N m^2].
    #[must_use]
    pub const fn bending_n_m2(&self) -> f64 {
        self.bending_n_m2
    }
    /// Coil mass per unit coil length [kg/m].
    #[must_use]
    pub const fn mass_per_m(&self) -> f64 {
        self.mass_per_m
    }
    /// Fraction of the coil's bending energy that is wire bending (the rest is wire torsion).
    #[must_use]
    pub const fn wire_bending_share(&self) -> f64 {
        self.wire_bending_share
    }

    /// Material part of the decay rate [1/s] of a mode at `omega` whose potential energy is the
    /// fraction `share` coil bending.
    pub fn material_rate(&self, omega: f64, share: f64) -> f64 {
        let ei = self.thermoelastic.elastic().bending_stiffness_n_m2();
        let (lo, hi) = self.excess.omega_band_rad_s();
        let excess = if (lo..=hi).contains(&omega) {
            self.excess.hysteretic_stiffness_n_m2() / ei
        } else {
            0.0
        };
        0.5 * omega
            * share
            * self.wire_bending_share
            * (self.thermoelastic.loss_tangent(omega) + excess)
    }

    /// Gas part of the decay rate [1/s] at `omega`: the turn-averaged drag over twice the
    /// wire's own mass per unit length (both per unit wire length).
    pub fn gas_rate(&self, omega: f64) -> Result<f64, String> {
        let a = self.geometry.wire_diameter_m / 2.0;
        let across = oscillating_cylinder(a, omega, &self.gas)?.resistance_n_s_per_m2;
        let along = sliding_cylinder(a, omega, &self.gas)?;
        let drag = across * (1.0 - self.cos2 / 2.0) + along * self.cos2 / 2.0;
        Ok(drag / (2.0 * self.mass_per_m * self.sin))
    }

    /// Decay rate [1/s] of a mode at `omega` whose potential energy is the fraction `share`
    /// coil bending.
    pub fn decay_rate(&self, omega: f64, share: f64) -> Result<f64, String> {
        Ok(self.material_rate(omega, share) + self.gas_rate(omega)?)
    }

    /// T60 [s] of a free coil wave at `frequency_hz` under `tension` [N], from the continuum
    /// dispersion mu omega^2 = T beta^2 + B beta^4 (for the notes).
    pub fn t60_s(&self, frequency_hz: f64, tension: f64) -> Result<f64, String> {
        let omega = 2.0 * PI * frequency_hz;
        let b = self.bending_n_m2;
        let beta2 = (-tension
            + (tension * tension + 4.0 * b * self.mass_per_m * omega * omega).sqrt())
            / (2.0 * b);
        let share = b * beta2 / (tension + b * beta2);
        Ok(6.907_755 / self.decay_rate(omega, share)?)
    }

    /// Receipts for a piece's notes.
    pub fn notes(&self, tension: f64) -> Result<Vec<String>, String> {
        let g = &self.geometry;
        let t60 = [100.0, 1000.0, 4000.0, 16_000.0]
            .iter()
            .map(|&f| {
                self.t60_s(f, tension)
                    .map(|t| format!("{f:.0} Hz {t:.2} s"))
            })
            .collect::<Result<Vec<_>, String>>()?;
        Ok(vec![
            format!(
                "snare wires: sourced coils ({}): wire {} mm, outer {} mm, so mean {:.3} mm, pitch {} mm, helix {:.1} deg; B {:.4e} N m^2 (wire torsion {:.0} % of its compliance), {:.4e} kg/m",
                g.source,
                g.wire_diameter_m * 1e3,
                g.outer_diameter_m * 1e3,
                g.mean_diameter_m() * 1e3,
                g.pitch_m * 1e3,
                g.helix_angle_rad().to_degrees(),
                self.bending_n_m2,
                100.0 * self.torsion_share,
                self.mass_per_m
            ),
            self.steel_note.clone(),
            format!(
                "coil loss per mode: thermoelastic + excess on the wire-bending share; wire torsion loss unmeasured, omitted (a lower bound); gas drag across and along the wire (exact linear Stokes flow, a lower bound once the wire's amplitude nears its {:.2} mm radius); T60 under the declared tension: {}",
                g.wire_diameter_m * 0.5e3,
                t60.join(", ")
            ),
            COIL_ENDS_NOTE.into(),
        ])
    }
}

/// Every coil states what is still declared and what is absent.
pub const COIL_ENDS_NOTE: &str = "snare wires, declared not sourced: strand count, spacing, length, tension, preload and contact loss (kit::Wires; the card route keeps 88cccbb78's, ruling 7486); ends pinned with no strand-end or cord damping, so the lowest coil modes ring long and 63-125 Hz is too loud";
