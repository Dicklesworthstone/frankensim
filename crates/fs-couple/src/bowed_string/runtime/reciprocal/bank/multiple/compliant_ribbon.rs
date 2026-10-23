//! Spatially resolved friction across a COMPLIANT bow-hair ribbon.
//!
//! `UniformBowRibbon` (see `ribbon.rs`) is the rigid limit: a uniform `1/n`
//! load and one common prescribed bow speed, each cell solving only its own
//! slip. Real bow hair is not rigid. Between the rigid bow stick (which moves
//! at the prescribed bow speed) and the string, the hair acts as a longitudinal
//! spring-dashpot in the bowing plane, so the stick speed reaches the string
//! only through a compliant element and the contact point can lag it.
//!
//! Pitteroff & Woodhouse, Acta Acustica 84 (1998) 543-562, "Finite width
//! bowing, Part I", model the ribbon as a homogeneous continuum and ADOPT the
//! spring-dashpot constitutive law (their eq. 1, used throughout as eq. 8):
//!
//! ```text
//!   f(z, t) = -( s * eta(z, t) + d * d/dt eta(z, t) )            (eq. 1)
//! ```
//!
//! with `eta` the hair displacement relative to the stick and `s`, `d` the
//! spring and damping constants DISTRIBUTED per unit length across the ribbon
//! width. The distributed tangential spring follows from the hair material and
//! geometry alone (their eq. 7, `gamma` = the excitation point as a fraction of
//! the hair length, `0.5` at mid-bow):
//!
//! ```text
//!   s = (1 / (gamma (1 - gamma))) * (n_H pi a_H^2 / b) * (E_H / l_H)   (eq. 7)
//! ```
//!
//! Equivalently the LUMPED total spring `S = s b = n_H pi a_H^2 E_H /
//! (l_H gamma (1 - gamma))` is the parallel stiffness of the `n_H` hairs and is
//! independent of the ribbon width `b`; a finite-width calculation spreads it
//! over the contact (P&W p. 547: distributed quantities multiplied by `b` for a
//! point bow, divided by `b` for a finite width). Every number is DATA (hair
//! Young's modulus, radius, count, free length, and the measured distributed
//! damping): swap the material and the compliance changes with NO code change.
//!
//! This module owns NO time loop and NO friction law. It derives the hair
//! compliance from data and lays out the finite-width contact cells, each with
//! its own tangential spring/dashpot. The normal-load distribution (M4.1) and
//! the differential slip (M4.2) EMERGE from the joint `MultiFrictionModalSystem`
//! solve; they are not prescribed here. Compare P&W Fig. 9: a rigid finite bow
//! carries singular edge forces, which the hair compliance spreads into a
//! smooth interior distribution.

use super::super::{ReciprocalStringError, invalid};
use crate::bowed_string::runtime::reciprocal::bank::StringBowLocation;

/// Measured hair-ribbon material and geometry. Every field is DATA, e.g. the
/// matdb pack `bow-hair-horse-pitteroff-woodhouse-acustica1998` (Pitteroff &
/// Woodhouse 1998). Nothing here is derived or defaulted.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HairRibbonProperties {
    /// Real part of the hair Young's modulus `E_H` [Pa] (bowing-plane
    /// longitudinal). P&W: 7 GPa (static / high-frequency) to 9 GPa (dynamic).
    pub youngs_modulus_pa: f64,
    /// Radius `a_H` of one hair [m] (half the measured diameter 0.16-0.25 mm).
    pub hair_radius_m: f64,
    /// Number of hairs `n_H` in immediate contact with the string (P&W: <= 50).
    pub hairs_in_contact: f64,
    /// Free hair length `l_H` between frog and tip [m] (P&W violin: 0.65 m).
    pub free_length_m: f64,
    /// Excitation point `gamma` as a fraction of the free length (mid-bow 0.5).
    pub excitation_fraction: f64,
    /// Measured damping `d` distributed per unit ribbon width [(kg/s)/m], from
    /// the single-hair decay tests (`d = n_H D_hair / b`; P&W violin: 1000).
    /// Zero models an undamped hair.
    pub distributed_damping_kg_s_m: f64,
    /// Hair mass density `rho_H` [kg/m^3] (P&W: about 1100). Retained as the
    /// provenance of the longitudinal wave speed cross-check `E_H = c_H^2 rho`.
    pub density_kg_m3: f64,
}

impl HairRibbonProperties {
    /// All data finite, positive, and `gamma` strictly interior.
    fn check(&self) -> Result<(), ReciprocalStringError> {
        let positive = [
            self.youngs_modulus_pa,
            self.hair_radius_m,
            self.hairs_in_contact,
            self.free_length_m,
            self.density_kg_m3,
        ];
        if positive.iter().any(|v| !v.is_finite() || *v <= 0.0)
            || !self.distributed_damping_kg_s_m.is_finite()
            || self.distributed_damping_kg_s_m < 0.0
            || !self.excitation_fraction.is_finite()
            || !(0.0 < self.excitation_fraction && self.excitation_fraction < 1.0)
        {
            return Err(invalid(
                "hair-ribbon properties need positive finite E, radius, count, length and density, \
                 a nonnegative finite distributed damping and an interior excitation fraction",
            ));
        }
        Ok(())
    }

    /// Lumped tangential spring `S = n_H pi a_H^2 E_H / (l_H gamma (1 - gamma))`
    /// [N/m], the parallel stiffness of the hairs in contact (P&W eq. 7 times
    /// the ribbon width). Independent of the ribbon width.
    #[must_use]
    fn lumped_tangential_stiffness_n_m(&self) -> f64 {
        let cross_section =
            self.hairs_in_contact * core::f64::consts::PI * self.hair_radius_m * self.hair_radius_m;
        cross_section * self.youngs_modulus_pa
            / (self.free_length_m * self.excitation_fraction * (1.0 - self.excitation_fraction))
    }
}

/// A compliant, finite-width bow-hair ribbon: the midpoint contact quadrature
/// of `UniformBowRibbon` plus a per-cell tangential spring/dashpot derived from
/// hair data. Contact zero is the original bound centre; other cells follow
/// from left to right. This is the compliant counterpart of the rigid
/// `UniformBowRibbon`; the normal-load distribution and the differential slip
/// are resolved by the joint contact solver, not stored here.
#[derive(Clone, Debug, PartialEq)]
pub struct CompliantBowRibbon {
    centre: StringBowLocation,
    width_m: f64,
    stations: Vec<f64>,
    properties: HairRibbonProperties,
    lumped_stiffness_n_m: f64,
    lumped_damping_kg_s: f64,
}

impl CompliantBowRibbon {
    /// Lay out `cells` midpoint-quadrature contacts over `width_m` along a
    /// string of length `length_m`, carrying the hair compliance from
    /// `properties`. The layout rules match `UniformBowRibbon` exactly so the
    /// point-bow limit and the rigid ribbon share a quadrature.
    ///
    /// # Errors
    /// Non-physical width/length, an even or out-of-range cell count, a ribbon
    /// reaching the string ends, indistinct stations, or invalid hair data.
    pub fn new(
        centre: StringBowLocation,
        length_m: f64,
        width_m: f64,
        cells: usize,
        properties: HairRibbonProperties,
    ) -> Result<Self, ReciprocalStringError> {
        properties.check()?;
        // Odd quadrature retains the already bound physical centre as port zero.
        if !width_m.is_finite()
            || width_m <= 0.0
            || !length_m.is_finite()
            || length_m <= 0.0
            || !(1..=63).contains(&cells)
            || cells % 2 == 0
        {
            return Err(invalid(
                "compliant ribbon needs positive finite width and an odd cell count in 1..=63",
            ));
        }
        let half = 0.5 * width_m / length_m;
        let left = centre.station_fraction - half;
        let right = centre.station_fraction + half;
        if !left.is_finite()
            || !right.is_finite()
            || left <= 0.0
            || right >= 1.0
            || !(left < centre.station_fraction && centre.station_fraction < right)
        {
            return Err(invalid(
                "the complete finite ribbon must lie strictly inside its actual string",
            ));
        }
        let step = width_m / length_m / cells as f64;
        let mid = cells / 2;
        let mut stations = Vec::with_capacity(cells);
        stations.push(centre.station_fraction);
        let mut previous = left;
        for cell in 0..cells {
            let station = centre.station_fraction + (cell as f64 - mid as f64) * step;
            if !(previous < station && station < right) {
                return Err(invalid(
                    "compliant ribbon cells must have distinct representable interior stations",
                ));
            }
            previous = station;
            if cell != mid {
                stations.push(station);
            }
        }
        let lumped_stiffness_n_m = properties.lumped_tangential_stiffness_n_m();
        let lumped_damping_kg_s = properties.distributed_damping_kg_s_m * width_m;
        Ok(Self {
            centre,
            width_m,
            stations,
            properties,
            lumped_stiffness_n_m,
            lumped_damping_kg_s,
        })
    }

    /// Actual bound ribbon centre on its string.
    #[must_use]
    pub const fn centre(&self) -> StringBowLocation {
        self.centre
    }
    /// Physical contact width along the string [m], not a tonal parameter.
    #[must_use]
    pub const fn width_m(&self) -> f64 {
        self.width_m
    }
    /// Actual contact stations, in the joint solver's port order.
    #[must_use]
    pub fn stations(&self) -> &[f64] {
        &self.stations
    }
    /// Measured hair data this ribbon derived its compliance from.
    #[must_use]
    pub const fn properties(&self) -> &HairRibbonProperties {
        &self.properties
    }
    /// Equal midpoint-quadrature share of the nominal normal force. The solved
    /// per-cell normal load departs from this when a co-located normal contact
    /// resolves it; this is only the quadrature weight.
    #[must_use]
    pub fn load_fraction(&self) -> f64 {
        1.0 / self.stations.len() as f64
    }
    /// Lumped tangential hair spring `S` [N/m] (P&W eq. 7 times the width): the
    /// total parallel stiffness of the hairs, independent of ribbon width.
    #[must_use]
    pub const fn lumped_tangential_stiffness_n_m(&self) -> f64 {
        self.lumped_stiffness_n_m
    }
    /// Lumped tangential hair dashpot `D = d b` [kg/s] (P&W p. 548).
    #[must_use]
    pub const fn lumped_tangential_damping_kg_s(&self) -> f64 {
        self.lumped_damping_kg_s
    }
    /// Distributed tangential spring `s = S / b` [(N/m)/m] (P&W eq. 7). The
    /// width-independent per-unit-length compliance the finite-width model uses.
    #[must_use]
    pub fn distributed_tangential_stiffness_n_m_per_m(&self) -> f64 {
        self.lumped_stiffness_n_m / self.width_m
    }
    /// Tangential spring of ONE contact cell [N/m]: the lumped spring shared
    /// equally over the midpoint cells, so the cells in parallel recover `S`.
    #[must_use]
    pub fn cell_tangential_stiffness_n_m(&self) -> f64 {
        self.lumped_stiffness_n_m / self.stations.len() as f64
    }
    /// Tangential dashpot of ONE contact cell [kg/s], the lumped `D` shared
    /// equally over the midpoint cells.
    #[must_use]
    pub fn cell_tangential_damping_kg_s(&self) -> f64 {
        self.lumped_damping_kg_s / self.stations.len() as f64
    }
    /// Longitudinal (bowing-plane) hair force on ONE cell for a relative
    /// displacement `eta` [m] and relative velocity `eta_dot` [m/s] between the
    /// string contact point and the rigid stick (P&W eq. 1/8):
    /// `f = -(k eta + c eta_dot)`.
    #[must_use]
    pub fn cell_tangential_force_n(&self, eta_m: f64, eta_dot_m_s: f64) -> f64 {
        -(self.cell_tangential_stiffness_n_m() * eta_m
            + self.cell_tangential_damping_kg_s() * eta_dot_m_s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn violin_hair() -> HairRibbonProperties {
        // Pitteroff & Woodhouse 1998 model violin bow (p. 546 / Table I).
        HairRibbonProperties {
            youngs_modulus_pa: 7.0e9,
            hair_radius_m: 0.19e-3 / 2.0,
            hairs_in_contact: 50.0,
            free_length_m: 0.65,
            excitation_fraction: 0.5,
            distributed_damping_kg_s_m: 1000.0,
            density_kg_m3: 1100.0,
        }
    }

    fn centre() -> StringBowLocation {
        StringBowLocation {
            string: 0,
            station_fraction: 0.12,
        }
    }

    #[test]
    fn lumped_spring_reproduces_pitteroff_woodhouse_violin() {
        // P&W quote S = 6 x 10^4 N/m (1 sig. fig.); eq. 7 gives 6.1e4.
        let ribbon = CompliantBowRibbon::new(centre(), 0.33, 0.010, 5, violin_hair()).unwrap();
        let s = ribbon.lumped_tangential_stiffness_n_m();
        assert!((60_000.0..62_000.0).contains(&s), "S = {s} N/m off eq.7");
        assert!(
            (s / 6.0e4 - 1.0).abs() <= 0.05,
            "S = {s} N/m off the paper's 6e4"
        );
        // Distributed s = S / b = 6.1e6 (N/m)/m (paper: 6e6).
        assert!((ribbon.distributed_tangential_stiffness_n_m_per_m() / 6.0e6 - 1.0).abs() <= 0.05);
        // Lumped D = d b = 1000 * 0.010 = 10 kg/s (paper: 10).
        assert!((ribbon.lumped_tangential_damping_kg_s() - 10.0).abs() <= 1e-9);
    }

    #[test]
    fn lumped_spring_is_width_independent_and_cells_recover_it() {
        let narrow = CompliantBowRibbon::new(centre(), 0.33, 0.004, 5, violin_hair()).unwrap();
        let wide = CompliantBowRibbon::new(centre(), 0.33, 0.010, 5, violin_hair()).unwrap();
        // S is a parallel stiffness of the hairs, not a function of the width.
        assert!(
            (narrow.lumped_tangential_stiffness_n_m() - wide.lumped_tangential_stiffness_n_m())
                .abs()
                <= 1e-6
        );
        // The per-cell springs in parallel recover the lumped spring exactly.
        for ribbon in [&narrow, &wide] {
            let cells = ribbon.stations().len() as f64;
            assert!(
                (ribbon.cell_tangential_stiffness_n_m() * cells
                    - ribbon.lumped_tangential_stiffness_n_m())
                .abs()
                    <= 1e-6
            );
        }
    }

    #[test]
    fn point_bow_cell_carries_the_whole_lumped_spring() {
        // A single-cell ribbon is the point bow: all stiffness at the centre.
        let point = CompliantBowRibbon::new(centre(), 0.33, 0.002, 1, violin_hair()).unwrap();
        assert_eq!(point.stations().len(), 1);
        assert!(
            (point.cell_tangential_stiffness_n_m() - point.lumped_tangential_stiffness_n_m()).abs()
                <= 1e-9
        );
        assert!((point.load_fraction() - 1.0).abs() <= 1e-12);
    }

    #[test]
    fn stiffer_or_more_numerous_hair_raises_the_spring() {
        let base = violin_hair();
        let base_s = CompliantBowRibbon::new(centre(), 0.33, 0.010, 5, base)
            .unwrap()
            .lumped_tangential_stiffness_n_m();
        let mut stiffer = base;
        stiffer.youngs_modulus_pa = 9.0e9; // dynamic modulus
        let stiffer_s = CompliantBowRibbon::new(centre(), 0.33, 0.010, 5, stiffer)
            .unwrap()
            .lumped_tangential_stiffness_n_m();
        assert!(stiffer_s > base_s);
        // Scaling is exactly linear in E and in hair count (eq. 7).
        assert!((stiffer_s / base_s - 9.0 / 7.0).abs() <= 1e-9);
    }

    #[test]
    fn constitutive_force_opposes_displacement_and_velocity() {
        let ribbon = CompliantBowRibbon::new(centre(), 0.33, 0.010, 5, violin_hair()).unwrap();
        let k = ribbon.cell_tangential_stiffness_n_m();
        let c = ribbon.cell_tangential_damping_kg_s();
        assert!((ribbon.cell_tangential_force_n(1e-6, 0.0) + k * 1e-6).abs() <= 1e-12);
        assert!((ribbon.cell_tangential_force_n(0.0, 2.0) + c * 2.0).abs() <= 1e-12);
        assert!(ribbon.cell_tangential_force_n(0.0, 0.0).abs() <= 0.0);
    }

    #[test]
    fn invalid_hair_data_and_geometry_refuse() {
        let good = violin_hair();
        // cap+1 cells and even cells refuse.
        assert!(CompliantBowRibbon::new(centre(), 0.33, 0.010, 65, good).is_err());
        assert!(CompliantBowRibbon::new(centre(), 0.33, 0.010, 4, good).is_err());
        // a ribbon wider than its interior reach refuses.
        assert!(CompliantBowRibbon::new(centre(), 0.33, 0.30, 5, good).is_err());
        // non-physical hair data refuses.
        for bad in [
            HairRibbonProperties {
                youngs_modulus_pa: 0.0,
                ..good
            },
            HairRibbonProperties {
                hair_radius_m: -1e-4,
                ..good
            },
            HairRibbonProperties {
                hairs_in_contact: 0.0,
                ..good
            },
            HairRibbonProperties {
                excitation_fraction: 0.0,
                ..good
            },
            HairRibbonProperties {
                excitation_fraction: 1.0,
                ..good
            },
            HairRibbonProperties {
                distributed_damping_kg_s_m: -1.0,
                ..good
            },
        ] {
            assert!(CompliantBowRibbon::new(centre(), 0.33, 0.010, 5, bad).is_err());
        }
    }
}
