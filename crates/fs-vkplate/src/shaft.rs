//! A cymbal struck with the stick's shaft across its edge (ruling 7373 (b)).
//!
//! Stick: the Z5A shape estimate of crates/fs-couple/examples/percussion (the published facts
//! are its 406.4 mm length, 0.560 inch diameter and hickory; the radius stations and the 0.12 m
//! grip are that example's estimates), turned about a pin at the grip. The edge meets the shaft
//! at a declared station. The shaft's contact-effective mass there is the pivot inertia over the
//! lever squared, the rule fs-couple's `StrikerProperties::from_profile` applies at the tip
//! (`kit::STICK_TIP_KG`). The density is the stick card's.
//!
//! Contact: the shaft crosses the edge at right angles, so it is two crossed cylinders: the
//! shaft's radius around its axis and the edge's rounding along it. The wood is indented across
//! its grain. The hickory card's transversely isotropic constants (`materials::TIP_ANISOTROPY`)
//! give the half-space compliance in the two planes through the indentation normal: the plane
//! holding the other transverse axis and the plane holding the grain (Barber, arXiv:1312.2779,
//! eqs 9 and 10; both reduce to (1 - nu^2)/(pi E) when isotropic). Between them the compliance
//! varies as cos 2 phi, the approximation Delafargue and Ulm (2004) found good for orthotropic
//! solids. The bronze adds its isotropic compliance. The pressure on the contact ellipse is
//! Hertzian (Willis 1966), solved by [`crate::hertz::EllipticHertz`], so the gap alone fixes
//! F = K delta^(3/2).
//!
//! Crush (ruling 7476 item 6): at audible stroke speeds that elastic contact's mean pressure
//! exceeds hickory's side hardness, so the wood crushes across its grain. The contact is a
//! [`Crush`]: Hertz until the mean pressure reaches the card's side hardness, read as a mean
//! pressure over the hardness ball's projected area, then crushing at that mean pressure on the
//! Hertz area, unloading from a permanent set (Thornton 1997). The Hertz area is half the edge's
//! truncation of the undeformed shaft, which is the area a hardness test divides by at depth, so at
//! a given crush depth this is the lower of the two readings of the same datum. The card's
//! compression perpendicular to the grain, reported as the stress at the proportional limit with
//! no clearly defined ultimate stress (FPL-GTR-282 chapter 5), is stated with each strike and not
//! used: as the crush pressure it would be a seventh of the hardness pressure.
//!
//! Not modelled: friction along the edge and the stick's bending (normal force only); each hit
//! crushes fresh wood. Unsourced and declared: the stick's angle (the shaft perpendicular to the
//! rim) and where the edge lands on the shaft (`SHOULDER_STATION_M`).

use crate::hertz::EllipticHertz;
use crate::materials::TIP_ANISOTROPY;
use crate::{Crush, Law};
use fs_matdb::{NormalizedMaterialCardPack, QueryPoint};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    YOUNG_MODULUS_PROPERTY, resolve_material_state_point,
};
use fs_qty::{Density, Dims, Force, Pressure};
use std::f64::consts::PI;

/// The Z5A shape estimate: (axial position from the butt, radius) [m], linear between stations
/// (crates/fs-couple/examples/percussion/main.rs). Not a manufacturer drawing.
pub const Z5A_PROFILE: [(f64, f64); 8] = [
    (0.0, 0.0068),
    (0.015, 0.007112),
    (0.280, 0.007112),
    (0.365, 0.0035),
    (0.380, 0.0025),
    (0.393, 0.0045),
    (0.403, 0.0035),
    (0.4064, 0.0),
];
/// The estimate's grip pivot [m from the butt].
pub const Z5A_GRIP_M: f64 = 0.12;
/// The estimate's tip contact [m from the butt], where `kit::STICK_TIP_KG` is integrated.
pub const Z5A_TIP_CONTACT_M: f64 = 0.400;
/// Where a crash edge meets the shaft [m from the butt]: declared, unsourced, the middle of the
/// estimate's taper (0.280 to 0.365 m). The taper's ends span about 18 to 9 g.
pub const SHOULDER_STATION_M: f64 = 0.33;
/// The shaft contact's Hunt-Crossley loss [s/m]: the card route's declared stick loss
/// (`kit::DECLARED_STICK_LOSS_S_M`, 88cccbb78's estimate, ruling 7486), not a card property.
pub const SHAFT_CONTACT_LOSS_S_M: f64 = crate::kit::DECLARED_STICK_LOSS_S_M;
/// The side-hardness test's ball [m]: "the load required to embed a 11.28-mm (0.444-in.) ball to
/// one-half its diameter" (FPL-GTR-282 chapter 5, Hardness). Embedded to half its diameter, its
/// projected area is that of its own diameter.
pub const SIDE_HARDNESS_BALL_DIAMETER_M: f64 = 11.28e-3;
/// The stick card's side hardness [N] and compression perpendicular to the grain [Pa].
pub const SIDE_HARDNESS_PROPERTY: &str = "side_hardness";
pub const COMPRESSION_PERPENDICULAR_PROPERTY: &str = "compression_perpendicular_proportional_limit";
/// Stated with every shaft strike.
pub const SHAFT_UNMODELLED_NOTE: &str = "stick contact, NOT modelled: friction along the edge, the stick's bending, and the crush of a hit before (each hit crushes fresh wood); UNSOURCED and declared: the stick's angle (the shaft perpendicular to the rim) and where the edge lands on it";

/// Radius of a profile at `x` [m], or None outside it.
#[must_use]
pub fn profile_radius_m(profile: &[(f64, f64)], x: f64) -> Option<f64> {
    profile
        .windows(2)
        .find(|w| w[0].0 <= x && x <= w[1].0)
        .map(|w| w[0].1 + (x - w[0].0) / (w[1].0 - w[0].0) * (w[1].1 - w[0].1))
}

/// Pivot inertia per unit density over the lever squared [m^3]: the contact-effective mass is
/// the density times this. Three-point Gauss per segment (exact for the piecewise-linear
/// radius), each disk's transverse inertia r^2/4 included, as fs-couple integrates the tip.
pub fn effective_volume_m3(
    profile: &[(f64, f64)],
    pivot_m: f64,
    contact_m: f64,
) -> Result<f64, String> {
    let inside = |x: f64| matches!((profile.first(), profile.last()), (Some(a), Some(b)) if a.0 <= x && x <= b.0);
    if profile.len() < 2
        || profile.windows(2).any(|w| !(w[1].0 > w[0].0))
        || profile.iter().any(|p| !(p.1 >= 0.0 && p.1.is_finite()))
        || !inside(pivot_m)
        || !inside(contact_m)
        || pivot_m == contact_m
    {
        return Err(format!(
            "stick profile needs ordered stations with finite radii and distinct pivot and contact inside it (pivot {pivot_m} m, contact {contact_m} m)"
        ));
    }
    let root = (3.0_f64 / 5.0).sqrt();
    let mut inertia = 0.0;
    for w in profile.windows(2) {
        let length = w[1].0 - w[0].0;
        for (point, weight) in [(-root, 5.0 / 9.0), (0.0, 8.0 / 9.0), (root, 5.0 / 9.0)] {
            let t = 0.5 * (point + 1.0);
            let x = w[0].0 + t * length;
            let r = w[0].1 + t * (w[1].1 - w[0].1);
            let dv = 0.5 * length * weight * PI * r * r;
            inertia += dv * ((x - pivot_m).powi(2) + 0.25 * r * r);
        }
    }
    Ok(inertia / (contact_m - pivot_m).powi(2))
}

/// Half-space compliances h (point-load deflection F h / r) of a solid transversely isotropic
/// about axis L, indented across it: (the plane through the normal and the other transverse
/// axis, the plane through the normal and L). Inputs as `materials::axial_indentation_modulus`.
pub fn across_grain_compliances(
    e_l: f64,
    e_p: f64,
    g_lp: f64,
    nu_lp: f64,
    nu_p: f64,
) -> Result<(f64, f64), String> {
    // Compliances with axis 3 = L; the plane (1, 2) is isotropic, the normal is axis 1.
    let s = [
        [1.0 / e_p, -nu_p / e_p, -nu_lp / e_l],
        [-nu_p / e_p, 1.0 / e_p, -nu_lp / e_l],
        [-nu_lp / e_l, -nu_lp / e_l, 1.0 / e_l],
    ];
    let minor =
        |r0: usize, r1: usize, c0: usize, c1: usize| s[r0][c0] * s[r1][c1] - s[r0][c1] * s[r1][c0];
    let det =
        s[0][0] * minor(1, 2, 1, 2) - s[0][1] * minor(1, 2, 0, 2) + s[0][2] * minor(1, 2, 0, 1);
    let (c11, c12, c13, c33) = (
        minor(1, 2, 1, 2) / det,
        -minor(0, 2, 1, 2) / det,
        minor(0, 1, 1, 2) / det,
        minor(0, 1, 0, 1) / det,
    );
    // Barber eq 9 with C22 = C11, 2 C66 = C11 - C12; eq 10 with C55 = G_LP.
    let plane = c11 / (PI * (c11 * c11 - c12 * c12));
    let grain = (c33 / (c11 * c33 - c13 * c13) * (1.0 / g_lp + 2.0 / (c13 + (c11 * c33).sqrt())))
        .sqrt()
        / (2.0 * PI);
    if !(det > 0.0
        && c11 > c12.abs()
        && c11 * c33 > c13 * c13
        && g_lp > 0.0
        && plane.is_finite()
        && plane > 0.0
        && grain.is_finite()
        && grain > 0.0)
    {
        return Err("shaft: the transversely isotropic constants are not positive definite".into());
    }
    Ok((plane, grain))
}

/// The shaft of the card's stick striking an isotropic edge.
#[derive(Clone, Debug)]
pub struct ShaftStrike {
    /// Contact-effective mass [kg].
    pub mass_kg: f64,
    /// The elastic contact law, `Law::Power` with the elastic K and `SHAFT_CONTACT_LOSS_S_M`.
    pub law: Law,
    /// The wood's crush across the grain, on that law's Hertz area.
    pub crush: Crush,
    pub contact: EllipticHertz,
    pub notes: Vec<String>,
}

fn requirement(
    name: &str,
    dims: Dims,
    admissibility: ScalarAdmissibility,
) -> Result<ScalarPropertyRequirement, String> {
    ScalarPropertyRequirement::try_new(name, dims, admissibility).map_err(|e| e.to_string())
}

impl ShaftStrike {
    /// The stick card at `point` struck at `station_m` of the Z5A estimate, against a target of
    /// isotropic (E Pa, nu) whose edge is rounded to `edge_radius_m`. The card must carry its
    /// density, E_L, all of `TIP_ANISOTROPY`, its side hardness and its compression
    /// perpendicular to the grain: an across-the-grain contact has no isotropic or elastic
    /// fallback.
    pub fn from_card(
        stick: &[u8],
        point: &QueryPoint,
        station_m: f64,
        edge_radius_m: f64,
        (e_target, nu_target): (f64, f64),
    ) -> Result<Self, String> {
        let pack = NormalizedMaterialCardPack::from_bytes(stick)
            .map_err(|e| format!("shaft stick card: {e}"))?;
        let mut needs = vec![
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
        ];
        for name in TIP_ANISOTROPY {
            needs.push(requirement(name, Dims::NONE, ScalarAdmissibility::Finite)?);
        }
        needs.push(requirement(
            SIDE_HARDNESS_PROPERTY,
            Force::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?);
        needs.push(requirement(
            COMPRESSION_PERPENDICULAR_PROPERTY,
            Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?);
        let wood = resolve_material_state_point(
            pack.card(),
            point,
            &needs,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .map_err(|e| format!("shaft stick card: {e}"))?;
        let value = |name: &str| {
            wood.property(name)
                .map(|p| p.value_si())
                .ok_or_else(|| format!("shaft stick card: {name} not resolved"))
        };
        let (rho, e_l) = (value(DENSITY_PROPERTY)?, value(YOUNG_MODULUS_PROPERTY)?);
        let [et, glt, nu_lt, nu_p] = [
            value(TIP_ANISOTROPY[0])?,
            value(TIP_ANISOTROPY[1])?,
            value(TIP_ANISOTROPY[2])?,
            value(TIP_ANISOTROPY[3])?,
        ];
        let (plane, grain) = across_grain_compliances(e_l, et * e_l, glt * e_l, nu_lt, nu_p)?;
        if !(e_target.is_finite() && e_target > 0.0 && (0.0..0.5).contains(&nu_target)) {
            return Err(format!(
                "shaft target needs E > 0 and 0 <= nu < 0.5, got {e_target} Pa, {nu_target}"
            ));
        }
        let target = (1.0 - nu_target * nu_target) / (PI * e_target);
        let radius = profile_radius_m(&Z5A_PROFILE, station_m)
            .filter(|r| *r > 0.0)
            .ok_or_else(|| format!("shaft station {station_m} m is not on the stick"))?;
        // x2 runs round the shaft, x3 along its axis (the grain), across the edge
        let contact = EllipticHertz::solve(
            plane + target,
            grain + target,
            0.5 / radius,
            0.5 / edge_radius_m,
        )?;
        let mass_kg = rho * effective_volume_m3(&Z5A_PROFILE, Z5A_GRIP_M, station_m)?;
        let hardness = value(SIDE_HARDNESS_PROPERTY)?;
        let perpendicular = value(COMPRESSION_PERPENDICULAR_PROPERTY)?;
        let ball = 0.5 * SIDE_HARDNESS_BALL_DIAMETER_M;
        let crush = Crush::new(contact.area_per_approach(), hardness / (PI * ball * ball))?;
        let modulus = |h: f64| 1e-9 / (PI * h);
        let notes = vec![
            format!(
                "stick: the shaft across the edge at {station_m:.3} m of the Z5A shape estimate (fs-couple examples/percussion), shaft radius {:.2} mm; contact-effective mass {:.2} g about the {Z5A_GRIP_M} m grip (pivot inertia / lever^2) at the density of state {} of card {}, {rho} kg/m^3",
                radius * 1e3,
                mass_kg * 1e3,
                wood.identity(),
                wood.card_identity()
            ),
            format!(
                "stick contact: Hertz of crossed cylinders, the shaft (radius {:.2} mm) across the edge (rounding {:.3} mm); the wood across the grain from E_L and the card's E_T/E_L, G_LT/E_L, nu_LT, nu_P (Barber eqs 9-10 with the Delafargue-Ulm cos 2 phi): {:.3} GPa round the shaft and {:.3} GPa along it, on {:.1} GPa; K {:.4e} N/m^1.5, contact b/a {:.3}; Hunt-Crossley {SHAFT_CONTACT_LOSS_S_M} s/m, the card route's declared stick loss (88cccbb78, ruling 7486), not a card property",
                radius * 1e3,
                edge_radius_m * 1e3,
                modulus(plane),
                modulus(grain),
                modulus(target),
                contact.k(),
                contact.aspect()
            ),
            format!(
                "stick crush across the grain (Thornton 1997 in mean-pressure form): Hertz until the mean pressure reaches {:.1} MPa, the card's side hardness {hardness:.0} N over the projected area of the {:.2} mm hardness ball, from an approach of {:.2} um ({:.2} N); deeper, that mean pressure on the Hertz area ({:.4e} m^2 per m of approach), half the edge's truncation of the shaft; unloading along Hertz from a permanent set. The card's compression perpendicular to the grain, {:.1} MPa at the proportional limit of a plate-loaded block, is not used",
                crush.pressure * 1e-6,
                SIDE_HARDNESS_BALL_DIAMETER_M * 1e3,
                crush.onset(contact.k()) * 1e6,
                crush.force(contact.k(), crush.onset(contact.k())),
                crush.area,
                perpendicular * 1e-6
            ),
            SHAFT_UNMODELLED_NOTE.into(),
        ];
        Ok(Self {
            mass_kg,
            law: Law::Power {
                k: contact.k(),
                alpha: 1.5,
                mu: SHAFT_CONTACT_LOSS_S_M,
            },
            crush,
            contact,
            notes,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64, rel: f64) {
        assert!(
            (actual - expected).abs() <= rel * expected.abs(),
            "{actual} != {expected} (rel {rel})"
        );
    }

    #[test]
    fn across_the_grain_is_isotropic_for_an_isotropic_solid() {
        let e = 12e9;
        for nu in [0.0, 0.25, 0.371, 0.45] {
            let (plane, grain) =
                across_grain_compliances(e, e, e / (2.0 * (1.0 + nu)), nu, nu).unwrap();
            close(plane, (1.0 - nu * nu) / (PI * e), 1e-12);
            close(grain, (1.0 - nu * nu) / (PI * e), 1e-12);
        }
    }

    #[test]
    fn the_tip_station_reproduces_the_kit_tip_mass() {
        // fs-couple printed 0.00724314 kg at 800 kg/m^3 (8 decimals).
        let v = effective_volume_m3(&Z5A_PROFILE, Z5A_GRIP_M, Z5A_TIP_CONTACT_M).unwrap();
        assert!(
            (800.0 * v - crate::kit::STICK_TIP_KG).abs() < 5e-9,
            "{}",
            800.0 * v
        );
        assert!(effective_volume_m3(&Z5A_PROFILE, Z5A_GRIP_M, Z5A_GRIP_M).is_err());
        assert!(effective_volume_m3(&Z5A_PROFILE, Z5A_GRIP_M, 0.5).is_err());
        close(
            profile_radius_m(&Z5A_PROFILE, 0.3225).unwrap(),
            0.005306,
            1e-12,
        );
        assert_eq!(profile_radius_m(&Z5A_PROFILE, 0.41), None);
    }
}
