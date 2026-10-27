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

use super::super::{ReciprocalStringBank, StringBowLocation};
use super::{MaterialBowContact, MultiBowStringBank};
use crate::bowed_string::runtime::reciprocal::{ReciprocalStringError, invalid, poll};
use crate::modal_acoustic_time::ModalAcousticState;
use crate::render::schedule::force::coupled::contact::friction::{
    FrictionDrive, multiple::MultiFrictionConfig,
};
use fs_exec::CancelGate;

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
    /// Ribbon width `b_ref` [m] at which `distributed_damping_kg_s_m` was stated
    /// (P&W violin: the 10 mm ribbon). The LUMPED total damping `D = d * b_ref`
    /// is then width-free, exactly as the lumped spring `S` is; a finite-width
    /// calculation spreads it over the contact without rescaling the total. Data.
    pub reference_width_m: f64,
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
            self.reference_width_m,
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

    /// Lumped tangential damping `D = d * b_ref` [kg/s], the width-free total
    /// dashpot of the hairs in contact (the distributed `d` times the ribbon
    /// width at which it was stated). Independent of the finite-width contact's
    /// own width, exactly as the lumped spring `S` is (frankensim-o7qbp: the
    /// earlier `d * contact_width` made the total damping scale with the width).
    #[must_use]
    fn lumped_tangential_damping_kg_s(&self) -> f64 {
        self.distributed_damping_kg_s_m * self.reference_width_m
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
        let lumped_damping_kg_s = properties.lumped_tangential_damping_kg_s();
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

    /// Split one total ribbon drive into the per-cell prescribed drives: the
    /// common bow-stick speed and the equal `1/n` share of the total normal
    /// load (uniform normal, P&W Part I). The hair compliance is NOT here; it is
    /// the per-cell series drive bound by `with_compliant_bow_ribbon`. Output
    /// storage is caller-owned and must cover every spatial cell.
    ///
    /// # Errors
    /// A non-finite or negatively loaded drive, or an output length mismatch.
    pub fn project_drive(
        &self,
        drive: FrictionDrive,
        output: &mut [FrictionDrive],
    ) -> Result<(), ReciprocalStringError> {
        if output.len() != self.stations.len()
            || !drive.speed_m_s.is_finite()
            || !drive.normal_force_n.is_finite()
            || drive.normal_force_n < 0.0
        {
            return Err(invalid(
                "ribbon drive requires finite speed, nonnegative total load and one output per cell",
            ));
        }
        output.fill(FrictionDrive {
            speed_m_s: drive.speed_m_s,
            normal_force_n: drive.normal_force_n * self.load_fraction(),
        });
        Ok(())
    }
}

impl ReciprocalStringBank {
    /// Resolve a COMPLIANT finite-width bow-hair ribbon on the existing bow's
    /// string. Like `with_uniform_bow_ribbon` it adds one jointly solved
    /// friction port per contact cell on the original string, but then binds a
    /// per-cell series spring-dashpot drive (M4.4 bow-hair compliance) whose
    /// stiffness and damping come from the hair DATA. The bow stick still moves
    /// at the prescribed common speed; each cell reaches the string only through
    /// its hair, so the contact-point velocities differ across the width and
    /// differential slip EMERGES (Pitteroff & Woodhouse 1998, Part I). The
    /// prescribed normal load is the uniform `1/n` share (Part I); a per-cell
    /// resolved normal distribution (Part III) is the existing contact-load
    /// composition and is not inferred here.
    ///
    /// # Errors
    /// An already bound interface, a cell count over the joint budget, a
    /// non-physical ribbon geometry or hair data, or a friction / series-drive
    /// refusal from the underlying owner.
    pub fn with_compliant_bow_ribbon(
        self,
        width_m: f64,
        cells: usize,
        hair: HairRibbonProperties,
        joint: MultiFrictionConfig,
        gate: &CancelGate,
    ) -> Result<MultiBowStringBank, ReciprocalStringError> {
        poll(gate)?;
        if self.bow_interface().is_some() {
            return Err(invalid(
                "bind ribbon interfaces per cell after resolving the contact area",
            ));
        }
        if cells > joint.max_ports {
            return Err(invalid(
                "ribbon cells exceed the declared joint contact budget",
            ));
        }
        let centre = self.source.options.bow;
        // Rank refusal (owner ruling 2026-10-23): the bowed string must resolve
        // the contact quadrature. With fewer than four modal degrees of freedom
        // per cell the joint compliance matrix over the closely spaced cells is
        // rank-deficient and the differential slip is not trustworthy, so refuse
        // rather than return an unreliable finite-width solve.
        if self.cards[centre.string].mode_count < 4 * cells {
            return Err(invalid(
                "rank-deficient ribbon: the bowed string needs at least four modes per contact cell",
            ));
        }
        let ribbon = CompliantBowRibbon::new(
            centre,
            self.cards[centre.string].length_m,
            width_m,
            cells,
            hair,
        )?;
        let law = self.system.port().law.clone();
        let friction = self.source.options.friction;
        let additional = ribbon.stations[1..]
            .iter()
            .map(|&station_fraction| MaterialBowContact {
                location: StringBowLocation {
                    string: centre.string,
                    station_fraction,
                },
                law: law.clone(),
                friction,
            })
            .collect();
        let bank = self.with_additional_bows(additional, joint, gate)?;
        // One series spring-dashpot per contact cell; the hairs in parallel are
        // the ribbon's lumped stiffness S, so each of the n cells carries S/n.
        let (k, c) = (
            ribbon.cell_tangential_stiffness_n_m(),
            ribbon.cell_tangential_damping_kg_s(),
        );
        let drives = (0..ribbon.stations.len()).map(|i| (i, k, c)).collect();
        let super::MultiBowStringBank {
            source,
            cards,
            system,
            stop_stations,
            bow_stations,
        } = bank;
        let system = system
            .with_series_compliant_drives(drives, gate)
            .map_err(ReciprocalStringError::Friction)?;
        Ok(super::MultiBowStringBank {
            source,
            cards,
            system,
            stop_stations,
            bow_stations,
        })
    }
}

impl MultiBowStringBank {
    /// Replace a string component's modal state before playback, for example to
    /// begin the string from a known Helmholtz motion rather than from rest.
    /// The common clock and every source binding are unchanged; any later
    /// contact history starts fresh.
    ///
    /// # Errors
    /// A sample already rendered, an unknown string index, or a state the modal
    /// basis or budgets reject.
    pub fn restore_string_states(
        &mut self,
        string: usize,
        states: &[ModalAcousticState],
    ) -> Result<(), ReciprocalStringError> {
        self.system
            .restore_component_states(string, states)
            .map_err(ReciprocalStringError::Friction)
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
            reference_width_m: 1.0e-2,
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

    #[test]
    fn project_drive_shares_load_and_keeps_common_speed() {
        let ribbon = CompliantBowRibbon::new(centre(), 0.33, 0.010, 5, violin_hair()).unwrap();
        let mut out = vec![
            FrictionDrive {
                speed_m_s: 0.0,
                normal_force_n: 0.0
            };
            5
        ];
        ribbon
            .project_drive(
                FrictionDrive {
                    speed_m_s: 0.2,
                    normal_force_n: 1.0,
                },
                &mut out,
            )
            .unwrap();
        for d in &out {
            assert!((d.speed_m_s - 0.2).abs() <= 0.0);
            assert!((d.normal_force_n - 0.2).abs() <= 1e-12); // 1.0 / 5 cells
        }
        // A wrong output length or a negative total load refuses.
        let mut short = vec![
            FrictionDrive {
                speed_m_s: 0.0,
                normal_force_n: 0.0
            };
            4
        ];
        assert!(
            ribbon
                .project_drive(
                    FrictionDrive {
                        speed_m_s: 0.2,
                        normal_force_n: 1.0
                    },
                    &mut short
                )
                .is_err()
        );
        assert!(
            ribbon
                .project_drive(
                    FrictionDrive {
                        speed_m_s: 0.2,
                        normal_force_n: -1.0
                    },
                    &mut out
                )
                .is_err()
        );
    }

    #[test]
    fn lumped_spring_and_damping_are_width_invariant() {
        // frankensim-o7qbp: the lumped hair spring S and the lumped hair damping D
        // are properties of the hairs in contact, not of how wide the finite-width
        // contact is laid out. Build the same hair over very different ribbon
        // widths and require S and D bit-identical; also check the violin lumped D
        // equals the P&W 10 kg/s (d = 1000 (kg/s)/m times b_ref = 10 mm).
        let hair = violin_hair();
        let wide = CompliantBowRibbon::new(centre(), 0.33, 1.0e-2, 5, hair).unwrap();
        let narrow = CompliantBowRibbon::new(centre(), 0.33, 1.0e-4, 5, hair).unwrap();
        assert_eq!(
            wide.lumped_stiffness_n_m, narrow.lumped_stiffness_n_m,
            "lumped spring moved with the ribbon width"
        );
        assert_eq!(
            wide.lumped_damping_kg_s, narrow.lumped_damping_kg_s,
            "lumped damping moved with the ribbon width"
        );
        assert!(
            (wide.lumped_damping_kg_s - 10.0).abs() <= 1e-9,
            "violin lumped D {} != 10 kg/s",
            wide.lumped_damping_kg_s
        );
    }
}
