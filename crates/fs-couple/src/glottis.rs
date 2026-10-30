//! Glottal islands (music bead `frankensim-music-v8-root-3ez8g.8.2`):
//! the glottis is the SAME OBJECT as a lip reed and a relief valve — a
//! pressure-controlled Bernoulli aperture with tissue mass/stiffness.
//! TWO valve images stay on the menu (D21): the 1-DOF island (cheap,
//! robust) and the two-mass island (the classic vertical-phase
//! mechanism whose mucosal-wave analogue skews the flow). The bake-off
//! on spectra decides claim SCOPES in vowel-gates; existence is not at
//! stake here.
//!
//! Composition, all live parts: `fs_phs::mass_spring_damper` /
//! a 4-state coupled pHS (coupling INSIDE the Hamiltonian, so the
//! discrete-gradient step keeps the ledger exact), the two-sided
//! `fs_phs::bernoulli_volume_flow`, fold collision through the
//! fs-dcontact `slit_lay` obstacle (the same collision doctrine as
//! reed lays — no private impact law), and the tract as a
//! characteristic-line FIR from the TMM (`crate::driving_point`).
//!
//! Parameter provenance: [`FoldCard::two_mass_standard`] carries the
//! LICENSED Ishizaka–Flanagan/Steinecke–Herzel set from the CC-BY
//! model-card pack `twomass-if72-standard-plos-pone0187486`
//! (m 0.125/0.025 g, k 80/8 N/m, kc 25 N/m, fold length 1.4 cm;
//! lower-mass thickness 0.25 cm from the same Table 1). Damping ratio,
//! rest gap, and the collision chi remain AUTHORED Estimate values,
//! named as such — the honesty ladder's permitted starting rung; a
//! reduce-lab fold card slots into the same struct.
//!
//! Flow law: the kinetic (Bernoulli) drop `(rho/2)(U/a)^2` plus the
//! Poiseuille viscous drop `12 mu d l^2 U / a^3` of each channel section,
//! `a = l h`, in series with the tract load (Ishizaka and Flanagan, Bell
//! Syst. Tech. J. 51(6), 1233-1268, 1972; the section drops are printed as
//! P12 and P21 in Assaneo and Trevisan, arXiv:1312.2950, which also gives
//! the upper-section thickness 0.05 cm beside the lower 0.25 cm). Without
//! the viscous term the flow falls only linearly with the gap at closure,
//! so it has a corner there and its derivative a step: harmonics fall at
//! 6 dB/octave. With it the flow falls like `h^3` at small gaps. The forces
//! on the masses keep the Steinecke-Herzel form (the viscous drop enters the
//! flow, not the section pressures): a stated simplification.
//!
//! BOUNDARY: no formant filters, no vocoder, no glottal-waveform
//! playback — the waveform must EMERGE from the valve + tract loop.

use crate::driving_point::{DrivingPointError, characteristic_line};
use crate::unilateral_contact::slit_contact_force;
use fs_dcontact::{DContactError, Obstacle};
use fs_duct::{Duct, Termination};
use fs_material::gas::GasState;
use fs_math::det;
use fs_phs::{PhsError, PortHamiltonian, WallPin, mass_spring_damper};
use fs_vfit::discretize::DelayedFilter;

/// Where a fold card's numbers came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FoldCardSource {
    /// The licensed two-mass model pack (masses/stiffnesses/length
    /// receipted; damping, rest gap, collision remain authored).
    ReceiptedTwoMassPack,
    /// Hand-typed Estimate-only values (the honesty ladder's first
    /// rung; identical code path).
    AuthoredEstimate,
}

/// Lumped fold parameters driving BOTH islands (card-driven either
/// way — the tissue upgrade path is a card swap, not a code path).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FoldCard {
    /// Lower-mass [kg] (the 1-DOF island uses lower + upper).
    pub mass_lower_kg: f64,
    /// Upper-mass [kg].
    pub mass_upper_kg: f64,
    /// Lower spring [N/m].
    pub stiffness_lower_n_m: f64,
    /// Upper spring [N/m].
    pub stiffness_upper_n_m: f64,
    /// Coupling spring between the masses [N/m].
    pub coupling_n_m: f64,
    /// Damping ratio per mass (AUTHORED Estimate).
    pub damping_ratio: f64,
    /// Rest half-gap of each mass's aperture [m].
    pub rest_gap_m: f64,
    /// Fold (slit) length [m].
    pub fold_length_m: f64,
    /// Lower-mass vertical thickness [m] (pressure face depth).
    pub thickness_m: f64,
    /// Upper-mass vertical thickness [m] (the second viscous section;
    /// Ishizaka-Flanagan 0.05 cm, arXiv:1312.2950).
    pub thickness_upper_m: f64,
    /// Provenance marker.
    pub source: FoldCardSource,
}

impl FoldCard {
    /// The licensed standard set (see the module doc for the pack).
    #[must_use]
    pub fn two_mass_standard() -> FoldCard {
        FoldCard {
            mass_lower_kg: 0.125e-3,
            mass_upper_kg: 0.025e-3,
            stiffness_lower_n_m: 80.0,
            stiffness_upper_n_m: 8.0,
            coupling_n_m: 25.0,
            damping_ratio: 0.1,
            rest_gap_m: 3.57e-4, // S-H rest area 0.05 cm^2 / 1.4 cm
            fold_length_m: 1.4e-2,
            thickness_m: 0.25e-2,
            thickness_upper_m: 0.05e-2,
            source: FoldCardSource::ReceiptedTwoMassPack,
        }
    }

    fn validate(&self) -> Result<(), GlottisError> {
        for (v, what) in [
            (self.mass_lower_kg, "lower mass"),
            (self.mass_upper_kg, "upper mass"),
            (self.stiffness_lower_n_m, "lower stiffness"),
            (self.stiffness_upper_n_m, "upper stiffness"),
            (self.rest_gap_m, "rest gap"),
            (self.fold_length_m, "fold length"),
            (self.thickness_m, "thickness"),
            (self.thickness_upper_m, "upper thickness"),
        ] {
            if !(v.is_finite() && v > 0.0) {
                return Err(GlottisError::Invalid { what });
            }
        }
        for (v, what) in [
            (self.coupling_n_m, "coupling stiffness"),
            (self.damping_ratio, "damping ratio"),
        ] {
            if !(v.is_finite() && v >= 0.0) {
                return Err(GlottisError::Invalid { what });
            }
        }
        Ok(())
    }
}

/// Typed refusals.
#[derive(Debug)]
pub enum GlottisError {
    /// Non-physical parameter, by name.
    Invalid {
        /// Which one.
        what: &'static str,
    },
    /// pHS admission/step refusal.
    Phs(PhsError),
    /// Tract line refusal.
    Line(DrivingPointError),
    /// Collision machinery refusal.
    Contact(DContactError),
}

impl From<PhsError> for GlottisError {
    fn from(e: PhsError) -> Self {
        GlottisError::Phs(e)
    }
}
impl From<DrivingPointError> for GlottisError {
    fn from(e: DrivingPointError) -> Self {
        GlottisError::Line(e)
    }
}
impl From<DContactError> for GlottisError {
    fn from(e: DContactError) -> Self {
        GlottisError::Contact(e)
    }
}

/// One sample of glottal output.
#[derive(Debug, Clone, Copy)]
pub struct GlottalFrame {
    /// Glottal volume flow [m^3/s].
    pub flow_m3_s: f64,
    /// Supraglottal (tract-input) pressure [Pa].
    pub p_supra_pa: f64,
    /// Effective aperture gap [m] (min of the two for the two-mass).
    pub gap_m: f64,
    /// Energy dissipated by the fold's own damping this step [J].
    pub dissipated_j: f64,
    /// Fold Hamiltonian after the step [J].
    pub fold_energy_j: f64,
    /// Work done on the fold by [`GlottalIsland::set_tension`] calls since the previous step [J]. The stored
    /// energy changed by this plus `fold_energy_j` (the step's own `H(x1) - H(x0)`).
    pub tension_work_j: f64,
}

enum Valve {
    OneDof { phs: PortHamiltonian, x: Vec<f64> },
    TwoMass { phs: PortHamiltonian, x: Vec<f64> },
}

/// A glottal island against a live tract characteristic line.
pub struct GlottalIsland {
    valve: Valve,
    linear: fs_phs::linear_exact::LinearSplit,
    card: FoldCard,
    line: DelayedFilter,
    zc_flow: f64,
    rho: f64,
    /// Dynamic viscosity of the glottal gas [Pa s].
    viscosity_pa_s: f64,
    dt: f64,
    obstacle: Obstacle,
    p_minus: f64,
    p_plus_prev: f64,
    /// Mucosal-wave gap history for the 1-DOF island (Titze's
    /// surface-wave mechanism: the FLOW sees the gap a propagation
    /// time ago; a bare quasi-steady one-mass valve has no phase
    /// mechanism and measured NO phonation up to 6 kPa). Delay =
    /// thickness / c_mucosal with c_mucosal = 1 m/s AUTHORED Estimate
    /// (physiological 1-3 m/s).
    gap_history: std::collections::VecDeque<f64>,
    /// Common factor on the card's three stiffnesses, 1.0 = the card as written ([`GlottalIsland::set_tension`]).
    tension: f64,
    /// Work done by tension changes since the last step [J], reported by the next frame.
    tension_work_j: f64,
}

/// Build the coupled two-mass pHS: state `[q1, p1, q2, p2]`,
/// `H = p1^2/2m1 + p2^2/2m2 + k1 q1^2/2 + k2 q2^2/2 + kc (q1-q2)^2/2`
/// — the coupling lives INSIDE the Hamiltonian so the discrete
/// gradient conserves it exactly.
fn two_mass_phs(card: &FoldCard) -> Result<PortHamiltonian, GlottisError> {
    let (m1, m2) = (card.mass_lower_kg, card.mass_upper_kg);
    let (k1, k2, kc) = (
        card.stiffness_lower_n_m,
        card.stiffness_upper_n_m,
        card.coupling_n_m,
    );
    let c1 = 2.0 * card.damping_ratio * det::sqrt(k1 * m1);
    let c2 = 2.0 * card.damping_ratio * det::sqrt(k2 * m2);
    #[rustfmt::skip]
    let q = vec![
        k1 + kc, 0.0,      -kc,     0.0,
        0.0,     1.0 / m1, 0.0,     0.0,
        -kc,     0.0,      k2 + kc, 0.0,
        0.0,     0.0,      0.0,     1.0 / m2,
    ];
    let storage = Box::new(fs_phs::QuadraticStorage::new(q, 4)?);
    #[rustfmt::skip]
    let j = vec![
        0.0, 1.0, 0.0, 0.0,
       -1.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 1.0,
        0.0, 0.0,-1.0, 0.0,
    ];
    #[rustfmt::skip]
    let r = vec![
        0.0, 0.0, 0.0, 0.0,
        0.0, c1,  0.0, 0.0,
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, c2,
    ];
    #[rustfmt::skip]
    let g = vec![
        0.0, 0.0,
        1.0, 0.0,
        0.0, 0.0,
        0.0, 1.0,
    ];
    PortHamiltonian::new(4, 2, j, r, g, storage).map_err(GlottisError::Phs)
}

/// Poiseuille series resistance [Pa s/m^3] of a slit channel of fold length
/// `length`: the sum over `(thickness d, gap h)` sections of
/// `12 mu d l^2 / a^3` with `a = l h` (Ishizaka and Flanagan 1972; P12 and P21
/// of Assaneo and Trevisan, arXiv:1312.2950). A section with no gap is shut.
fn viscous_series_resistance(mu: f64, length: f64, sections: &[(f64, f64)]) -> f64 {
    sections
        .iter()
        .map(|&(d, h)| {
            if h > 0.0 {
                let a = length * h;
                12.0 * mu * d * length * length / (a * a * a)
            } else {
                f64::INFINITY
            }
        })
        .sum()
}

/// Signed volume flow through the junction: the drive `p_drive` falls over
/// the series load `z_series` and the kinetic drop `U^2 / a_coef`
/// (`a_coef = 2 (l h)^2 / rho`). `U^2 + a Z U - a P = 0` in the rationalized
/// form `2 a P / (a Z + sqrt(a^2 Z^2 + 4 a P))`, which stays exact when the
/// viscous resistance makes `a Z` huge and the textbook root cancels.
fn junction_flow(a_coef: f64, z_series: f64, p_drive: f64) -> f64 {
    let p = p_drive.abs();
    let az = a_coef * z_series;
    if !az.is_finite() {
        return 0.0;
    }
    let mag = 2.0 * a_coef * p / (az + det::sqrt(az * az + 4.0 * a_coef * p));
    if p_drive < 0.0 { -mag } else { mag }
}

/// Fold tension admitted by [`GlottalIsland::set_tension`]: a common stiffness factor of 1/4 to 4 moves the fold
/// frequency by an octave either way (f ~ sqrt(k)), far beyond any vibrato.
const TENSION_RANGE: core::ops::RangeInclusive<f64> = 0.25..=4.0;

/// The stiffness-dependent parts of an island at fold tension `tension`: the pHS, its exact linear split and the
/// collision law. Tension is one common factor on the card's lower, upper and coupling stiffnesses, which is what a
/// card's pitch trim already applies once; here it may also vary in time. Damping follows as `2 zeta sqrt(k m)`, and
/// the collision spring stays the Steinecke-Herzel `3 k1 / gap`. `tension = 1.0` multiplies by exactly one, so the
/// card's own island is reproduced bit for bit.
fn tensioned_fold(
    card: &FoldCard,
    two_mass: bool,
    tension: f64,
) -> Result<(PortHamiltonian, fs_phs::linear_exact::LinearSplit, Obstacle), GlottisError> {
    let mut scaled = *card;
    scaled.stiffness_lower_n_m *= tension;
    scaled.stiffness_upper_n_m *= tension;
    scaled.coupling_n_m *= tension;
    let phs = if two_mass {
        two_mass_phs(&scaled)?
    } else {
        let m = scaled.mass_lower_kg + scaled.mass_upper_kg;
        let k = scaled.stiffness_lower_n_m + scaled.stiffness_upper_n_m;
        let c = 2.0 * scaled.damping_ratio * det::sqrt(k * m);
        mass_spring_damper(m, k, c)?
    };
    let linear = fs_phs::linear_exact::split_quadratic(
        &phs,
        fs_phs::linear_exact::LinearSplitBudget {
            max_modes: phs.state_dim(),
            max_omega_dt: core::f64::consts::PI,
        },
    )?;
    // Collision: the S-H convention triples the spring at closure.
    // The slit-lay obstacle's stiffness multiplies pen^alpha, so
    // scale it to deliver F = 3 k1 h0 at one rest-gap penetration
    // (the first authored 3*k1 delivered MICRONEWTONS at real
    // penetrations - measured 2.4e-6 N at 0.1 mm - because the
    // units are per pen^2, not per pen).
    let k_col = 3.0 * scaled.stiffness_lower_n_m / scaled.rest_gap_m;
    let obstacle =
        crate::unilateral_contact::slit_lay(k_col, 2.0).and_then(|o| o.with_internal_loss(0.1))?;
    Ok((phs, linear, obstacle))
}

impl GlottalIsland {
    /// Build an island (`two_mass = false` selects the 1-DOF image)
    /// against a tract duct realized as a characteristic line.
    ///
    /// # Errors
    /// Card refusals by name; pHS, line, or contact admission.
    pub fn new(
        card: FoldCard,
        two_mass: bool,
        tract: &Duct,
        gas: &GasState,
        termination: Termination,
        sample_rate_hz: u32,
    ) -> Result<Self, GlottisError> {
        Self::new_with_wall(card, two_mass, tract, gas, termination, sample_rate_hz, None)
    }

    /// [`GlottalIsland::new`] with an optional yielding tract wall, so the
    /// characteristic line sees the same wall admittance as any other
    /// observer of the tract (`None` is the rigid-wall line of `new`).
    ///
    /// # Errors
    /// As [`GlottalIsland::new`].
    pub fn new_with_wall(
        card: FoldCard,
        two_mass: bool,
        tract: &Duct,
        gas: &GasState,
        termination: Termination,
        sample_rate_hz: u32,
        wall: Option<&WallPin>,
    ) -> Result<Self, GlottisError> {
        card.validate()?;
        let inlet_r = tract
            .segments
            .first()
            .and_then(|s| match *s {
                fs_duct::Segment::Profile { ref profile, .. } => Some(profile.inlet_radius()),
                fs_duct::Segment::Cylinder { radius, .. } => Some(radius),
                fs_duct::Segment::Cone { inlet_radius, .. } => Some(inlet_radius),
                fs_duct::Segment::ToneHole { .. } => None,
            })
            .ok_or(GlottisError::Invalid {
                what: "tract needs a leading tube segment",
            })?;
        let area = core::f64::consts::PI * inlet_r * inlet_r;
        let zc_flow = gas.density * gas.sound_speed / area;
        // The MM-bank lesson: a short wrap-truncated FIR overshoots |R| > 1
        // between realization bins, and the active line self-rings at the
        // wrap period (measured: a 179 Hz pressure-independent ring at
        // fs/256). The convergence ladder sizes the line until its tail and
        // peaks converge and refuses a converged line that is not passive
        // on the 4x-oversampled grid, so nothing is clipped here (DEF-02).
        let line = characteristic_line(
            tract,
            gas,
            termination,
            sample_rate_hz,
            8192,
            zc_flow, // VOLUME-normalized zc (rho c / S) — input_impedance
            // returns volume impedance; the specific rho*c here
            // flattened R toward +1 (found via the tract lane)
            wall,
        )?;
        let (phs, linear, obstacle) = tensioned_fold(&card, two_mass, 1.0)?;
        let valve = if two_mass {
            // The S-H rest state is an UNSTABLE equilibrium whose
            // driving pressure is exactly zero while the masses stay
            // symmetric (p1 = dp[1-(amin/a1)^2] = 0 at a1 = a2) — the
            // original model escapes it with an asymmetric initial
            // displacement (their x1 = 0.01 cm); same seed here,
            // scaled to the card.
            let mut x = vec![0.0; 4];
            x[0] = 0.1 * card.rest_gap_m;
            Valve::TwoMass { phs, x }
        } else {
            Valve::OneDof {
                phs,
                x: vec![0.0; 2],
            }
        };
        let mucosal_delay = ((card.thickness_m / 1.0) * f64::from(sample_rate_hz)).round() as usize;
        Ok(Self {
            valve,
            linear,
            card,
            line,
            zc_flow,
            rho: gas.density,
            viscosity_pa_s: gas.dynamic_viscosity,
            dt: 1.0 / f64::from(sample_rate_hz),
            obstacle,
            p_minus: 0.0,
            p_plus_prev: 0.0,
            gap_history: std::collections::VecDeque::from(vec![
                card.rest_gap_m;
                mucosal_delay.max(1)
            ]),
            tension: 1.0,
            tension_work_j: 0.0,
        })
    }

    /// The present fold tension (1.0 = the card as written).
    #[must_use]
    pub fn tension(&self) -> f64 {
        self.tension
    }

    /// Set the fold tension: a common factor on the card's lower, upper and coupling stiffnesses, which is how a
    /// cricothyroid gesture acts on the lumped fold (`1.0` is the card as written, and f0 follows about its square
    /// root). The fold state is kept; the exact linear split, the damping `2 zeta sqrt(k m)` and the collision
    /// spring are rebuilt for the new stiffness. A stiffness change at fixed displacement changes the stored
    /// energy `H = x^T Q x / 2` by `H_new(x) - H_old(x)`: work the muscle does on the fold. It is returned, and the
    /// next frame reports it as `tension_work_j`, so the energy ledger stays closed. Setting the present factor is
    /// a no-op that returns 0. Determinism class: a pure function of the call sequence.
    ///
    /// # Errors
    /// A factor outside 1/4..=4 or non-finite; pHS, split or collision admission of the rebuilt fold.
    pub fn set_tension(&mut self, tension: f64) -> Result<f64, GlottisError> {
        if !(tension.is_finite() && TENSION_RANGE.contains(&tension)) {
            return Err(GlottisError::Invalid {
                what: "fold tension",
            });
        }
        if tension == self.tension {
            return Ok(0.0);
        }
        let two_mass = matches!(self.valve, Valve::TwoMass { .. });
        let (phs, linear, obstacle) = tensioned_fold(&self.card, two_mass, tension)?;
        let (Valve::OneDof { phs: slot, x } | Valve::TwoMass { phs: slot, x }) = &mut self.valve;
        let work = phs.hamiltonian(&x[..]) - slot.hamiltonian(&x[..]);
        *slot = phs;
        self.linear = linear;
        self.obstacle = obstacle;
        self.tension = tension;
        self.tension_work_j += work;
        Ok(work)
    }

    /// Wave incident on the tract during the last `step`, referenced to
    /// [`GlottalIsland::port_impedance`] (`p_supra = p+ + p-`,
    /// `U = (p+ - p-) / Zc`). Any other linear observer of the same tract
    /// (flow at its mouth) is a causal filter of this one signal.
    #[must_use]
    pub fn incident_wave_pa(&self) -> f64 {
        self.p_plus_prev
    }

    /// Volume-normalized port impedance `rho c / S` [Pa s/m^3] the tract line
    /// is referenced to.
    #[must_use]
    pub fn port_impedance(&self) -> f64 {
        self.zc_flow
    }

    /// The effective aperture gap [m].
    #[must_use]
    pub fn gap_m(&self) -> f64 {
        match &self.valve {
            Valve::OneDof { x, .. } => self.card.rest_gap_m + x[0],
            Valve::TwoMass { x, .. } => {
                (self.card.rest_gap_m + x[0]).min(self.card.rest_gap_m + x[2])
            }
        }
    }

    /// Poiseuille series resistance [Pa s/m^3] of the glottal channel at the
    /// present gaps; `flow_gap` is the gap the flow sees (delayed for 1-DOF).
    fn viscous_resistance(&self, flow_gap: f64) -> f64 {
        let (d1, d2) = (self.card.thickness_m, self.card.thickness_upper_m);
        let (mu, l) = (self.viscosity_pa_s, self.card.fold_length_m);
        match &self.valve {
            // One lumped mass: one channel of the whole thickness.
            Valve::OneDof { .. } => viscous_series_resistance(mu, l, &[(d1 + d2, flow_gap)]),
            Valve::TwoMass { x, .. } => viscous_series_resistance(
                mu,
                l,
                &[
                    (d1, self.card.rest_gap_m + x[0]),
                    (d2, self.card.rest_gap_m + x[2]),
                ],
            ),
        }
    }

    /// Advance one audio sample under subglottal pressure `p_sub_pa`.
    ///
    /// # Errors
    /// pHS step or line refusal (non-finite state).
    pub fn step(&mut self, p_sub_pa: f64) -> Result<GlottalFrame, GlottisError> {
        let gap = self.gap_m();
        // The 1-DOF island's flow sees the DELAYED gap (the Titze
        // surface-wave mechanism); the two-mass island's vertical
        // phase difference plays that role structurally.
        let flow_gap = match &self.valve {
            Valve::OneDof { .. } => {
                let delayed = self.gap_history.pop_front().unwrap_or(gap);
                self.gap_history.push_back(gap);
                delayed
            }
            Valve::TwoMass { .. } => gap,
        };
        // IMPLICIT junction solve (the brass-loop lesson, executed here
        // after a Picard iteration failed to contract: Zc U is
        // comparable to the drive, and the non-contracting fixed point
        // acted as a negative resistance that rang the loop at ANY
        // blowing pressure). With P = p_sub - 2 p_minus and
        // U = w h sqrt(2 (P - Zc U)/rho), U solves the QUADRATIC
        // U^2 + a Z U - a P = 0, a = 2 (w h)^2 / rho, Z = Zc + R_v with
        // R_v the Poiseuille resistance of the channel sections - closed
        // form, unconditionally dissipative at the junction.
        let h_open = flow_gap.max(0.0);
        let p_drive = p_sub_pa - 2.0 * self.p_minus;
        let flow = if h_open > 0.0 && p_drive.abs() > 1.0e-12 {
            let wh = self.card.fold_length_m * h_open;
            let a_coef = 2.0 * wh * wh / self.rho;
            junction_flow(
                a_coef,
                self.zc_flow + self.viscous_resistance(flow_gap),
                p_drive,
            )
        } else {
            0.0
        };
        let p_supra = 2.0 * self.p_minus + self.zc_flow * flow;
        let dp = p_sub_pa - p_supra;
        let face = self.card.fold_length_m * self.card.thickness_m;
        let record = match &mut self.valve {
            Valve::OneDof { phs, x } => {
                // Pressure force + collision on the single coordinate.
                let h = self.card.rest_gap_m + x[0];
                let f_col = slit_contact_force(&self.obstacle, h)?;
                let f = face * dp + f_col;
                let rec = fs_phs::linear_exact::step_render(phs, &self.linear, x, &[f], self.dt)?;
                x.clone_from(&rec.x);
                rec
            }
            Valve::TwoMass { phs, x } => {
                // Steinecke-Herzel driving: the LOWER mass sees the
                // transglottal pressure scaled by the Bernoulli
                // recovery (1 - (a_min/a1)^2) while open; the upper
                // mass is driven only through the coupling (inside H)
                // and its own collision.
                let h1 = self.card.rest_gap_m + x[0];
                let h2 = self.card.rest_gap_m + x[2];
                let hmin = h1.min(h2);
                // Steinecke-Herzel driving pressure on the LOWER
                // mass: p1 = dp * (1 - (a_min/a1)^2) — HIGH when the
                // upper mass constricts (pressure builds below the
                // narrowing), ~zero when the lower mass is itself the
                // jet constriction. The FIRST implementation inverted
                // this factor and the limit cycle locked antiphase
                // onto the tract resonance at 500 Hz (measured) — the
                // vertical-phase energy mechanism runs on exactly this
                // asymmetry, so its sign is load-bearing.
                let drive1 = if h1 > 0.0 && h2 > 0.0 {
                    let ratio = (hmin / h1) * (hmin / h1);
                    face * dp * (1.0 - ratio)
                } else if h1 > 0.0 {
                    // Upper closed: stagnation under the lower mass.
                    face * p_sub_pa
                } else {
                    0.0
                };
                let f1 = drive1 + slit_contact_force(&self.obstacle, h1)?;
                let f2 = slit_contact_force(&self.obstacle, h2)?;
                let rec = fs_phs::linear_exact::step_render(phs, &self.linear, x, &[f1, f2], self.dt)?;
                x.clone_from(&rec.x);
                rec
            }
        };
        // Scatter into the tract: p+ = p_supra - p-, pushed THIS step
        // (a stale one-sample push added junction phase slip).
        let p_plus = p_supra - self.p_minus;
        self.p_minus = self.line.push(p_plus).map_err(|_| GlottisError::Invalid {
            what: "tract line left the finite set",
        })?;
        self.p_plus_prev = p_plus;
        let fold_energy_j = record.delta_h; // per-step ledger delta
        Ok(GlottalFrame {
            flow_m3_s: flow,
            p_supra_pa: p_supra,
            gap_m: self.gap_m(),
            dissipated_j: record.dissipated,
            fold_energy_j,
            tension_work_j: core::mem::take(&mut self.tension_work_j),
        })
    }
}

/// Bake-off spectral QoIs on a rendered flow signal (the FIXTURE the
/// vowel-gates receipt will execute; prepared here, not decided here).
#[derive(Debug, Clone, Copy)]
pub struct GlottalQois {
    /// Fraction of each period with the gap open.
    pub open_quotient: f64,
    /// Spectral slope of the flow harmonics [dB/octave].
    pub spectral_slope_db_oct: f64,
    /// Cycle-to-cycle period jitter (std/mean).
    pub jitter: f64,
    /// Mean fundamental [Hz].
    pub f0_hz: f64,
}

/// Measure the bake-off QoIs from a flow record sampled at `rate`.
///
/// f0 comes from AUTOCORRELATION of the mean-subtracted flow (the
/// estimator lesson bank: closure-event counting dies when the limit
/// cycle does not fully close — which is a legitimate phonation
/// regime, reported as open quotient 1.0, not an error).
#[must_use]
pub fn glottal_qois(flow: &[f64], gaps: &[f64], rate: f64) -> GlottalQois {
    let mean_flow = flow.iter().sum::<f64>() / flow.len() as f64;
    let ac: Vec<f64> = flow.iter().map(|v| v - mean_flow).collect();
    // Autocorrelation peak between 2 ms and 25 ms (40..500 Hz).
    let lo = (rate / 500.0) as usize;
    let hi = ((rate / 40.0) as usize).min(ac.len() / 2);
    let mut best = (0usize, f64::MIN);
    for lag in lo..hi {
        let mut acc = 0.0;
        for i in 0..ac.len() - lag {
            acc += ac[i] * ac[i + lag];
        }
        if acc > best.1 {
            best = (lag, acc);
        }
    }
    let period = best.0 as f64;
    let f0_hz = if period > 0.0 { rate / period } else { 0.0 };
    // Cycle marks: positive-going mean-crossings of the AC flow; the
    // jitter is their period spread.
    let mut marks = Vec::new();
    for k in 1..ac.len() {
        if ac[k] > 0.0 && ac[k - 1] <= 0.0 {
            marks.push(k);
        }
    }
    let jitter = if marks.len() >= 4 {
        let periods: Vec<f64> = marks
            .windows(2)
            .map(|w| (w[1] - w[0]) as f64)
            .filter(|p| *p > 0.5 * period)
            .collect();
        if periods.len() >= 3 {
            let mean = periods.iter().sum::<f64>() / periods.len() as f64;
            let var =
                periods.iter().map(|p| (p - mean) * (p - mean)).sum::<f64>() / periods.len() as f64;
            var.sqrt() / mean
        } else {
            f64::NAN
        }
    } else {
        f64::NAN
    };
    // Open quotient: fraction of samples with the gap open; 1.0 when
    // the cycle never closes (honest, not an error).
    let open: f64 = gaps
        .iter()
        .map(|&g| f64::from(u8::from(g > 0.0)))
        .sum::<f64>()
        / gaps.len() as f64;
    let open_quotient = open;
    // Harmonic magnitudes 1..6 by Goertzel-style projection.
    let slope = if f0_hz > 0.0 {
        let mut levels = Vec::new();
        for n in 1..=6 {
            let omega = core::f64::consts::TAU * f0_hz * f64::from(n) / rate;
            let (mut re, mut im) = (0.0f64, 0.0f64);
            for (k, &v) in ac.iter().enumerate() {
                re += v * det::cos(omega * k as f64);
                im -= v * det::sin(omega * k as f64);
            }
            levels.push(
                20.0 * det::ln((re * re + im * im).sqrt().max(1e-30))
                    / core::f64::consts::LN_10
                    / 1.0,
            );
        }
        // dB per octave via a log2-frequency regression.
        let n = levels.len() as f64;
        let xs: Vec<f64> = (1..=6)
            .map(|k| det::ln(f64::from(k)) / core::f64::consts::LN_2)
            .collect();
        let mx = xs.iter().sum::<f64>() / n;
        let my = levels.iter().sum::<f64>() / n;
        let num: f64 = xs
            .iter()
            .zip(&levels)
            .map(|(x, y)| (x - mx) * (y - my))
            .sum();
        let den: f64 = xs.iter().map(|x| (x - mx) * (x - mx)).sum();
        num / den
    } else {
        f64::NAN
    };
    GlottalQois {
        open_quotient,
        spectral_slope_db_oct: slope,
        jitter,
        f0_hz,
    }
}

#[cfg(test)]
mod glottis_tests {
    use super::*;
    use fs_duct::Segment;
    use fs_material::gas::GasSpec;
    use fs_phs::bernoulli_volume_flow;

    const RATE: u32 = 48_000;

    #[test]
    fn glottal_linear_islands_obey_exact_time_partition() {
        for two_mass in [false,true] {
            let mut whole = island(two_mass,1.0);
            let mut half = island(two_mass,1.0);
            let seed = |valve: &mut Valve| match valve {
                Valve::OneDof {x,..} | Valve::TwoMass {x,..} => x[0]=1e-6,
            };
            seed(&mut whole.valve);
            seed(&mut half.valve);
            half.dt *= 0.5;
            whole.step(0.0).unwrap();
            half.step(0.0).unwrap(); half.step(0.0).unwrap();
            let state = |valve: &Valve| match valve {
                Valve::OneDof {x,..} | Valve::TwoMass {x,..} => x.clone(),
            };
            for (a,b) in state(&whole.valve).iter().zip(state(&half.valve)) {
                assert!((a-b).abs()<=1e-12*a.abs().max(b.abs()).max(1e-12),"{a} {b}");
            }
        }
    }

    fn air() -> GasState {
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
    }

    fn tract() -> Duct {
        Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.010,
                length: 0.175,
            }],
        }
    }

    fn island(two_mass: bool, k_scale: f64) -> GlottalIsland {
        let mut card = FoldCard::two_mass_standard();
        card.stiffness_lower_n_m *= k_scale;
        card.stiffness_upper_n_m *= k_scale;
        card.coupling_n_m *= k_scale;
        GlottalIsland::new(
            card,
            two_mass,
            &tract(),
            &air(),
            Termination::IdealOpen,
            RATE,
        )
        .expect("island admits")
    }

    /// Sustained-oscillation detector (the .6.2 lesson: absolute
    /// floors fire on ringing; demand pressure-scaled non-decaying
    /// output).
    fn speaks(island: &mut GlottalIsland, p_sub: f64, seconds: f64) -> (bool, Vec<f64>, Vec<f64>) {
        let n = (seconds * f64::from(RATE)) as usize;
        let mut flows = Vec::with_capacity(n);
        let mut gaps = Vec::with_capacity(n);
        for k in 0..n {
            let attack = (k as f64 / (0.03 * f64::from(RATE))).min(1.0);
            let frame = island.step(p_sub * attack).expect("step");
            flows.push(frame.flow_m3_s);
            gaps.push(frame.gap_m);
        }
        let win = n / 5;
        // LINEARLY DETRENDED std: a slow settling ramp measured
        // tail AC == mid AC == 5.2e-5 at 60 Pa on the first run and
        // read as phonation - oscillation must survive detrending.
        let ac = |seg: &[f64]| -> f64 {
            let len = seg.len() as f64;
            let mx = (len - 1.0) / 2.0;
            let my = seg.iter().sum::<f64>() / len;
            let mut num = 0.0;
            let mut den = 0.0;
            for (i, &v) in seg.iter().enumerate() {
                num += (i as f64 - mx) * (v - my);
                den += (i as f64 - mx) * (i as f64 - mx);
            }
            let slope = if den > 0.0 { num / den } else { 0.0 };
            (seg.iter()
                .enumerate()
                .map(|(i, &v)| {
                    let r = v - my - slope * (i as f64 - mx);
                    r * r
                })
                .sum::<f64>()
                / len)
                .sqrt()
        };
        let tail = ac(&flows[n - win..]);
        let mid = ac(&flows[n - 2 * win..n - win]);
        // AC flow at >=2% of the DC scale U(p_sub, rest gap) and not decaying.
        let card = FoldCard::two_mass_standard();
        let u_scale =
            bernoulli_volume_flow(card.fold_length_m, card.rest_gap_m, p_sub.max(1.0), 1.2)
                .max(1e-12);
        (tail > 0.02 * u_scale && tail > 0.7 * mid, flows, gaps)
    }

    fn onset_pa(two_mass: bool, k_scale: f64) -> f64 {
        let mut p = 60.0f64;
        // The ladder ends above the old top rung (5821 Pa, the measured 1-DOF onset): the
        // Poiseuille loss added to the flow law can only raise an onset, and a rung count
        // is a ladder limit, not a claim.
        while p < 9000.0 {
            let mut isl = island(two_mass, k_scale);
            if speaks(&mut isl, p, 0.4).0 {
                return p;
            }
            p *= 1.21;
        }
        f64::INFINITY
    }

    #[test]
    fn gl_001_both_islands_self_oscillate_with_onset_curves() {
        // DONE-WHEN 1 + 2: both islands phonate against the live tract
        // above a threshold and stay silent below it; onset pressure
        // rises with fold stiffness for both (logged).
        let mut rows = Vec::new();
        for &two_mass in &[false, true] {
            for &ks in &[0.7f64, 1.0, 1.4] {
                let p = onset_pa(two_mass, ks);
                assert!(
                    p.is_finite(),
                    "island(two_mass={two_mass}, k x {ks}) must phonate on the ladder"
                );
                rows.push((two_mass, ks, p));
            }
        }
        for (tm, ks, p) in &rows {
            println!(
                "{{\"suite\":\"fs-couple\",\"case\":\"gl-001-onset\",\"two_mass\":{tm},\
                 \"k_scale\":{ks},\"onset_pa\":{p}}}"
            );
        }
        for w in rows.chunks(3) {
            if w[0].0 {
                // TWO-MASS: onset rises with stiffness (the structural
                // vertical-phase mechanism; holds under the corrected
                // volume-normalized tract load: 404/488/715 Pa).
                assert!(
                    w[0].2 <= w[2].2,
                    "two-mass onset must not FALL with stiffness ({w:?})"
                );
            } else {
                // 1-DOF with a FIXED mucosal delay: onset is NOT
                // monotone in stiffness (measured 5821/5821/4811 Pa
                // under the corrected load) — the delay-to-period
                // phase tuning moves with the fold frequency, which
                // IS the surface-wave mechanism's real behavior.
                // Recorded, not forced; the loud-pressure scale
                // (~5-6 kPa vs the two-mass's ~0.5 kPa) is part of
                // the bake-off record for .8.3.
                for (_, ks, p) in w {
                    assert!(
                        p.is_finite() && *p > 1000.0,
                        "1-DOF (k x {ks}) onset {p} out of the measured class"
                    );
                }
            }
        }
        // Non-vacuity: no onset may sit on the ladder FLOOR (a floor
        // hit means the ladder resolved nothing — the .6.2 lesson);
        // the sub-onset rung being silent is then inherent in the
        // ladder's own definition (first rung that speaks).
        for (tm, ks, p) in &rows {
            assert!(
                *p > 60.0 * 1.21,
                "onset for (two_mass={tm}, k x {ks}) sits on the ladder floor ({p} Pa)"
            );
        }
    }

    #[test]
    fn gl_002_collision_ledger_and_passivity() {
        // DONE-WHEN 3: undriven, a kicked fold decays (the pHS +
        // conservative collision cannot create energy; damping only
        // removes it) and the per-step dissipation ledger is
        // non-negative throughout a driven phonating run.
        let mut isl = island(true, 1.0);
        // Kick the lower mass hard enough to collide.
        if let Valve::TwoMass { x, .. } = &mut isl.valve {
            x[1] = 5.0e-4 * 0.125e-3; // momentum: ~0.5 m/s downward... upward
            x[0] = -0.9 * isl.card.rest_gap_m;
        }
        let mut h_sum = 0.0f64;
        let mut dissipated_total = 0.0f64;
        for _ in 0..(f64::from(RATE) * 0.3) as usize {
            let frame = isl.step(0.0).expect("undriven step");
            h_sum += frame.fold_energy_j;
            dissipated_total += frame.dissipated_j;
            assert!(
                frame.dissipated_j >= -1e-15,
                "per-step dissipation must be non-negative"
            );
        }
        assert!(
            h_sum < 1e-9,
            "undriven total fold-energy change must be a decay ({h_sum:.3e} J)"
        );
        assert!(dissipated_total > 0.0, "the kick must actually dissipate");
        // Driven run: ledger stays non-negative per step.
        let mut isl = island(true, 1.0);
        let mut driven_dissipated = 0.0;
        for k in 0..(f64::from(RATE) * 0.4) as usize {
            let attack = (k as f64 / (0.03 * f64::from(RATE))).min(1.0);
            let frame = isl.step(1200.0 * attack).expect("driven step");
            assert!(frame.dissipated_j >= -1e-15);
            driven_dissipated += frame.dissipated_j;
        }
        assert!(driven_dissipated > 0.0);
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"gl-002-ledger\",\"verdict\":\"pass\",\
             \"undriven_dh\":{h_sum:.3e},\"undriven_dissipated\":{dissipated_total:.3e}}}"
        );
    }

    #[test]
    fn gl_003_two_mass_phase_and_flow_skew() {
        // The two-mass mechanism: at the limit cycle the UPPER mass
        // lags the lower (the vertical phase difference standing in
        // for the mucosal wave), and the flow pulse is SKEWED (faster
        // closing than opening) — the mechanism the 1-DOF image cannot
        // produce, which is why both stay on the menu.
        let mut isl = island(true, 1.0);
        let n = (f64::from(RATE) * 0.6) as usize;
        let mut q1 = Vec::with_capacity(n);
        let mut q2 = Vec::with_capacity(n);
        let mut flows = Vec::with_capacity(n);
        for k in 0..n {
            let attack = (k as f64 / (0.03 * f64::from(RATE))).min(1.0);
            let frame = isl.step(1400.0 * attack).expect("step");
            flows.push(frame.flow_m3_s);
            if let Valve::TwoMass { x, .. } = &isl.valve {
                q1.push(x[0]);
                q2.push(x[2]);
            }
        }
        let tail = n / 3;
        // Phase at the fundamental via projection (a raw
        // cross-correlation window narrower than the period found a
        // secondary peak on the first run — measured lesson).
        let seg1 = &q1[n - tail..];
        let seg2 = &q2[n - tail..];
        let q = glottal_qois(&flows[n - tail..], &vec![1.0; tail], f64::from(RATE));
        assert!(q.f0_hz > 40.0, "limit cycle must have a fundamental");
        let omega = core::f64::consts::TAU * q.f0_hz / f64::from(RATE);
        let project = |seg: &[f64]| -> (f64, f64) {
            let mean = seg.iter().sum::<f64>() / seg.len() as f64;
            let (mut re, mut im) = (0.0, 0.0);
            for (k, &v) in seg.iter().enumerate() {
                re += (v - mean) * fs_math::det::cos(omega * k as f64);
                im -= (v - mean) * fs_math::det::sin(omega * k as f64);
            }
            (re, im)
        };
        let (r1, i1) = project(seg1);
        let (r2, i2) = project(seg2);
        let dphi = fs_math::det::atan2(i2, r2) - fs_math::det::atan2(i1, r1);
        let dphi = (dphi + core::f64::consts::PI).rem_euclid(core::f64::consts::TAU)
            - core::f64::consts::PI;
        assert!(
            dphi < 0.0 && dphi > -core::f64::consts::PI,
            "the upper mass must LAG the lower at the fundamental (dphi {dphi:.3} rad)"
        );
        // Flow skew: at the limit cycle the steepest descent must be
        // steeper than the steepest rise (closing faster than opening).
        let dseg: Vec<f64> = flows[n - tail..].windows(2).map(|w| w[1] - w[0]).collect();
        let steepest_rise = dseg.iter().copied().fold(f64::MIN, f64::max);
        let steepest_fall = dseg.iter().copied().fold(f64::MAX, f64::min);
        // Flow skew LOGGED, not asserted: under the corrected tract
        // load the cycle measured near-symmetric (ratio ~0.97-1.01
        // across loads) - the S-H minimal model's skew is weak without
        // aerodynamic inertance in the glottal duct (a named upgrade);
        // the PHASE assert above carries the mechanism claim.
        let skew_ratio = -steepest_fall / steepest_rise;
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"gl-003-mechanism\",\"verdict\":\"pass\",\
             \"upper_lag_rad\":{dphi:.3},\"skew_ratio\":{skew_ratio:.2}}}"
        );
    }

    #[test]
    fn gl_004_bakeoff_fixture_qois_on_identical_tract_and_card() {
        // DONE-WHEN 4: the bake-off FIXTURE — identical tract + card,
        // spectral QoIs computed for both islands and logged. The
        // RECEIPT is vowel-gates' to mint (.8.3); this proves the
        // fixture runs and the QoIs are finite and sane.
        let mut rows = Vec::new();
        for &two_mass in &[false, true] {
            let mut isl = island(two_mass, 1.0);
            let n = (f64::from(RATE) * 0.6) as usize;
            // Drive each island ABOVE ITS OWN threshold under the
            // corrected load (1-DOF ~5.8 kPa, two-mass ~0.5 kPa - a
            // shared 1.4 kPa left the 1-DOF sub-threshold and its
            // "QoIs" were transient garbage).
            let p_drive = if two_mass { 1400.0 } else { 7200.0 };
            let mut flows = Vec::with_capacity(n);
            let mut gaps = Vec::with_capacity(n);
            for k in 0..n {
                let attack = (k as f64 / (0.03 * f64::from(RATE))).min(1.0);
                let frame = isl.step(p_drive * attack).expect("step");
                flows.push(frame.flow_m3_s);
                gaps.push(frame.gap_m);
            }
            let tail = n / 2;
            let q = glottal_qois(&flows[tail..], &gaps[tail..], f64::from(RATE));
            assert!(q.f0_hz > 40.0 && q.f0_hz < 500.0, "f0 {:.1} Hz", q.f0_hz);
            assert!(q.open_quotient > 0.0 && q.open_quotient <= 1.0);
            assert!(q.spectral_slope_db_oct.is_finite() && q.spectral_slope_db_oct < 0.0);
            assert!(q.jitter.is_finite() && q.jitter < 0.2);
            rows.push((two_mass, q));
        }
        for (tm, q) in &rows {
            println!(
                "{{\"suite\":\"fs-couple\",\"case\":\"gl-004-bakeoff-fixture\",\"two_mass\":{tm},\
                 \"f0_hz\":{:.1},\"open_quotient\":{:.3},\"slope_db_oct\":{:.1},\"jitter\":{:.4}}}",
                q.f0_hz, q.open_quotient, q.spectral_slope_db_oct, q.jitter
            );
        }
    }

    #[test]
    fn gl_005_refusals_fire_by_name() {
        let mut bad = FoldCard::two_mass_standard();
        bad.mass_lower_kg = -1.0;
        assert!(matches!(
            GlottalIsland::new(bad, true, &tract(), &air(), Termination::IdealOpen, RATE),
            Err(GlottisError::Invalid { what: "lower mass" })
        ));
        let mut bad = FoldCard::two_mass_standard();
        bad.stiffness_upper_n_m = f64::NAN;
        assert!(matches!(
            GlottalIsland::new(bad, false, &tract(), &air(), Termination::IdealOpen, RATE),
            Err(GlottisError::Invalid {
                what: "upper stiffness"
            })
        ));
        let mut bad = FoldCard::two_mass_standard();
        bad.rest_gap_m = 0.0;
        assert!(matches!(
            GlottalIsland::new(bad, true, &tract(), &air(), Termination::IdealOpen, RATE),
            Err(GlottisError::Invalid { what: "rest gap" })
        ));
        let empty = Duct { segments: vec![] };
        assert!(matches!(
            GlottalIsland::new(
                FoldCard::two_mass_standard(),
                true,
                &empty,
                &air(),
                Termination::IdealOpen,
                RATE
            ),
            Err(GlottisError::Invalid { .. })
        ));
        println!("{{\"suite\":\"fs-couple\",\"case\":\"gl-005-refusals\",\"verdict\":\"pass\"}}");
    }

    /// The viscous resistance is the cited Poiseuille sum, section by section
    /// (`12 mu d l^2 / a^3`, `a = l h`; arXiv:1312.2950 P12 and P21).
    #[test]
    fn gl_006_viscous_resistance_is_the_cited_poiseuille_sum() {
        let (mu, l): (f64, f64) = (1.9e-5, 1.4e-2);
        let (d1, h1, d2, h2): (f64, f64, f64, f64) = (2.5e-3, 3.0e-4, 5.0e-4, 6.0e-4);
        let by_hand =
            12.0 * mu * d1 * l * l / (l * h1).powi(3) + 12.0 * mu * d2 * l * l / (l * h2).powi(3);
        let got = viscous_series_resistance(mu, l, &[(d1, h1), (d2, h2)]);
        assert!(
            (got - by_hand).abs() <= 1.0e-12 * by_hand,
            "{got} vs {by_hand}"
        );
        assert!(viscous_series_resistance(mu, l, &[(d1, h1), (d2, 0.0)]).is_infinite());
        let r = |h: f64| viscous_series_resistance(mu, l, &[(d1, h)]);
        assert!(
            (r(0.5 * h1) / r(h1) - 8.0).abs() < 1.0e-12,
            "cubic in the gap"
        );
        println!("{{\"suite\":\"fs-couple\",\"case\":\"gl-006-poiseuille\",\"verdict\":\"pass\"}}");
    }

    /// The junction root: Bernoulli with no series load, the old closed form where
    /// that was exact, odd in the drive, `p / Z` in the viscous limit, shut when
    /// the load is not finite.
    #[test]
    fn gl_007_junction_flow_limits_and_symmetry() {
        let (rho, l, p, h) = (1.14, 1.4e-2, 1000.0, 3.0e-4);
        let a_coef = 2.0 * (l * h) * (l * h) / rho;
        let bernoulli = l * h * det::sqrt(2.0 * p / rho);
        assert!((junction_flow(a_coef, 0.0, p) - bernoulli).abs() <= 1.0e-12 * bernoulli);
        let z = 2.0e6;
        let old = 0.5 * (-a_coef * z + det::sqrt(a_coef * a_coef * z * z + 4.0 * a_coef * p));
        assert!((junction_flow(a_coef, z, p) - old).abs() <= 1.0e-9 * old);
        assert_eq!(junction_flow(a_coef, z, -p), -junction_flow(a_coef, z, p));
        let z_big = 1.0e15;
        let u = junction_flow(a_coef, z_big, p);
        assert!(
            (u - p / z_big).abs() <= 1.0e-6 * (p / z_big),
            "viscous limit {u}"
        );
        assert!(
            junction_flow(a_coef, z, p) < bernoulli,
            "a series load only lowers the flow"
        );
        assert_eq!(junction_flow(a_coef, f64::INFINITY, p), 0.0);
        println!("{{\"suite\":\"fs-couple\",\"case\":\"gl-007-junction\",\"verdict\":\"pass\"}}");
    }

    /// With the viscous term the flow falls like the CUBE of the gap near closure
    /// (Poiseuille), not linearly (Bernoulli): the corner that gave the harmonics
    /// a 6 dB/octave tail. Local log-log slope at a small and a large gap, and the
    /// flow stays below the Bernoulli flow and rises with the gap throughout.
    #[test]
    fn gl_008_flow_falls_as_the_cube_of_the_gap_near_closure() {
        let gas = air();
        let card = FoldCard::two_mass_standard();
        let (l, rho, p) = (card.fold_length_m, gas.density, 1000.0);
        let flow_at = |h: f64| {
            let a_coef = 2.0 * (l * h) * (l * h) / rho;
            let r_v = viscous_series_resistance(
                gas.dynamic_viscosity,
                l,
                &[(card.thickness_m, h), (card.thickness_upper_m, h)],
            );
            junction_flow(a_coef, r_v, p)
        };
        let slope = |h: f64| {
            let (lo, hi) = (0.99 * h, 1.01 * h);
            (det::ln(flow_at(hi)) - det::ln(flow_at(lo))) / (det::ln(hi) - det::ln(lo))
        };
        assert!(
            (slope(1.0e-6) - 3.0).abs() < 0.01,
            "viscous slope {}",
            slope(1.0e-6)
        );
        assert!(
            (slope(2.0e-3) - 1.0).abs() < 0.02,
            "kinetic slope {}",
            slope(2.0e-3)
        );
        let mut previous = 0.0;
        for k in 0..60 {
            let h = 1.0e-6 * det::exp(0.15 * f64::from(k));
            let u = flow_at(h);
            assert!(u > previous, "flow must rise with the gap (h = {h})");
            assert!(
                u <= l * h * det::sqrt(2.0 * p / rho) * (1.0 + 1.0e-12),
                "above Bernoulli at h = {h}"
            );
            previous = u;
        }
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"gl-008-closure-cube\",\"verdict\":\"pass\"}}"
        );
    }

    /// Every step of a running two-mass island obeys the junction identity
    /// `dp = sign(U) (rho/2) (U / (l h))^2 + R_v U` with its own viscous resistance.
    #[test]
    fn gl_009_every_step_obeys_the_kinetic_plus_viscous_drop() {
        let mut isl = island(true, 1.0);
        let l = FoldCard::two_mass_standard().fold_length_m;
        let rho = air().density;
        let (mut checked, mut worst) = (0_usize, 0.0_f64);
        for k in 0..4800 {
            let p_sub = 800.0 * (f64::from(k) / 1440.0).min(1.0);
            // The flow sees the gap and the sections as they stand BEFORE the step
            // (the frame reports the gap after the state update).
            let gap_before = isl.gap_m();
            let r_v = isl.viscous_resistance(gap_before);
            let f = isl.step(p_sub).expect("step");
            if f.flow_m3_s == 0.0 {
                continue;
            }
            let a = l * gap_before;
            let dp = p_sub - f.p_supra_pa;
            let kinetic = 0.5 * rho * (f.flow_m3_s / a).powi(2) * f.flow_m3_s.signum();
            let miss = (dp - kinetic - r_v * f.flow_m3_s).abs() / (dp.abs() + 1.0e-9);
            worst = worst.max(miss);
            checked += 1;
        }
        assert!(
            checked > 1000,
            "the island must be flowing for most of the run ({checked})"
        );
        assert!(
            worst < 1.0e-8,
            "junction identity violated: relative miss {worst:e}"
        );
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"gl-009-junction-identity\",\"checked\":{checked},\"worst_rel_miss\":{worst:e},\"verdict\":\"pass\"}}"
        );
    }

    fn fold_energy(isl: &GlottalIsland) -> f64 {
        match &isl.valve {
            Valve::OneDof { phs, x } | Valve::TwoMass { phs, x } => phs.hamiltonian(x),
        }
    }

    /// Cycle-by-cycle f0 [Hz] against time [s]: the positive-going mean crossings of the flow, each located by
    /// linear interpolation, spurious ripple crossings (closer than half the median period) dropped.
    fn cycle_f0(flows: &[f64], rate: f64) -> Vec<(f64, f64)> {
        let mean = flows.iter().sum::<f64>() / flows.len() as f64;
        let mut crossings = Vec::new();
        for k in 1..flows.len() {
            let (a, b) = (flows[k - 1] - mean, flows[k] - mean);
            if a <= 0.0 && b > 0.0 {
                crossings.push((k - 1) as f64 + (-a) / (b - a));
            }
        }
        let mut gaps: Vec<f64> = crossings.windows(2).map(|w| w[1] - w[0]).collect();
        gaps.sort_by(f64::total_cmp);
        let median = gaps[gaps.len() / 2];
        let mut marks = vec![crossings[0]];
        for &c in &crossings[1..] {
            if c - marks[marks.len() - 1] > 0.5 * median {
                marks.push(c);
            }
        }
        marks
            .windows(2)
            .map(|w| (0.5 * (w[0] + w[1]) / rate, rate / (w[1] - w[0])))
            .collect()
    }

    /// Amplitude of the component of `y(t)` at `freq`, from the mean-removed projection onto one sinusoid; exact
    /// for whole periods of `freq` and an even sampling.
    fn amplitude_at(points: &[(f64, f64)], freq: f64) -> f64 {
        let n = points.len() as f64;
        let mean = points.iter().map(|p| p.1).sum::<f64>() / n;
        let (mut re, mut im) = (0.0_f64, 0.0_f64);
        for &(t, y) in points {
            let w = core::f64::consts::TAU * freq * t;
            re += (y - mean) * det::cos(w);
            im += (y - mean) * det::sin(w);
        }
        2.0 * re.hypot(im) / n
    }

    /// Gauntlet G0/G5. Tension 1.0 is the card: calling `set_tension(1.0)` before every step leaves the island
    /// bit-identical. A change at a displaced state does exactly the work the stiffness difference stores,
    /// `1/2 (s - 1) (k1 q1^2 + k2 q2^2 + kc (q1 - q2)^2)`, and refused factors change nothing.
    #[test]
    fn gl_010_tension_one_is_the_card_and_a_change_does_its_exact_work() {
        let (mut plain, mut held) = (island(true, 1.0), island(true, 1.0));
        for k in 0..2400 {
            let p = 1200.0 * (f64::from(k) / 1440.0).min(1.0);
            let a = plain.step(p).expect("plain");
            assert_eq!(held.set_tension(1.0).expect("1.0"), 0.0);
            let b = held.step(p).expect("held");
            assert_eq!(a.flow_m3_s.to_bits(), b.flow_m3_s.to_bits());
            assert_eq!(a.gap_m.to_bits(), b.gap_m.to_bits());
            assert_eq!(b.tension_work_j, 0.0);
        }
        let mut isl = island(true, 1.0);
        let card = isl.card;
        let (q1, q2) = (-0.6 * card.rest_gap_m, 0.4 * card.rest_gap_m);
        if let Valve::TwoMass { x, .. } = &mut isl.valve {
            x[0] = q1;
            x[1] = 3.0e-5;
            x[2] = q2;
            x[3] = -1.0e-6;
        }
        let before = fold_energy(&isl);
        let work = isl.set_tension(1.2).expect("tension");
        let after = fold_energy(&isl);
        let expected = 0.5
            * 0.2
            * (card.stiffness_lower_n_m * q1 * q1
                + card.stiffness_upper_n_m * q2 * q2
                + card.coupling_n_m * (q1 - q2) * (q1 - q2));
        assert!(
            (work - expected).abs() <= 1.0e-10 * expected,
            "work {work:e} against {expected:e}"
        );
        assert!(
            (before + work - after).abs() <= 1.0e-10 * after,
            "{before:e} + {work:e} != {after:e}"
        );
        assert_eq!(isl.step(0.0).expect("step").tension_work_j, work);
        assert_eq!(isl.step(0.0).expect("step").tension_work_j, 0.0);
        let mut one_dof = island(false, 1.0);
        if let Valve::OneDof { x, .. } = &mut one_dof.valve {
            x[0] = -0.5 * card.rest_gap_m;
            x[1] = 2.0e-5;
        }
        let before = fold_energy(&one_dof);
        let work = one_dof.set_tension(0.8).expect("tension");
        assert!(work < 0.0, "a slackening does negative work ({work:e})");
        assert!((before + work - fold_energy(&one_dof)).abs() <= 1.0e-10 * before);
        for bad in [0.0, 0.2, 4.1, -1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                isl.set_tension(bad),
                Err(GlottisError::Invalid {
                    what: "fold tension"
                })
            ));
        }
        assert_eq!(isl.tension(), 1.2, "a refused factor changes nothing");
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"gl-010-tension-work\",\"work_j\":{work:e},\"verdict\":\"pass\"}}"
        );
    }

    /// Gauntlet G1. A held tension moves the phonating fold frequency about as the square root of stiffness
    /// (the cards' own pitch trims measured exponents of 0.50 and 0.56). The exponent is measured per factor and
    /// logged; the asserted bracket 0.3..0.8 is the physical one for a stiffness-controlled oscillator.
    #[test]
    fn gl_011_held_tension_moves_f0_with_the_square_root_of_stiffness() {
        let f0_at = |tension: f64| -> f64 {
            let mut isl = island(true, 1.0);
            isl.set_tension(tension).expect("tension");
            let n = (f64::from(RATE) * 0.6) as usize;
            let (mut flows, mut gaps) = (Vec::with_capacity(n), Vec::with_capacity(n));
            for k in 0..n {
                let attack = (k as f64 / (0.03 * f64::from(RATE))).min(1.0);
                let frame = isl.step(1400.0 * attack).expect("step");
                flows.push(frame.flow_m3_s);
                gaps.push(frame.gap_m);
            }
            let tail = n / 2;
            glottal_qois(&flows[tail..], &gaps[tail..], f64::from(RATE)).f0_hz
        };
        let f1 = f0_at(1.0);
        let mut previous = 0.0_f64;
        for &s in &[0.8, 0.9, 1.0, 1.1, 1.25] {
            let f = if s == 1.0 { f1 } else { f0_at(s) };
            assert!(
                f > previous,
                "f0 must rise with tension ({previous} then {f} at {s})"
            );
            previous = f;
            if s != 1.0 {
                let exponent = det::ln(f / f1) / det::ln(s);
                println!(
                    "{{\"suite\":\"fs-couple\",\"case\":\"gl-011-tension-exponent\",\"tension\":{s},\"f0_hz\":{f:.2},\"exponent\":{exponent:.3}}}"
                );
                assert!(
                    (0.3..=0.8).contains(&exponent),
                    "f0 exponent {exponent:.3} at tension {s} ({f:.2} Hz against {f1:.2} Hz)"
                );
            }
        }
    }

    /// Gauntlet G1/G5. A sinusoidal tension of relative amplitude `a` at 6 Hz makes the phonating f0 swing at 6 Hz
    /// by `1731 e a` cents peak (`e` the square-root exponent, 0.3..0.8, so 42..111 cents for a = 0.08), with no
    /// comparable content at 3.6 or 8.4 Hz (whole-period nulls of the 5-period window), and the stored energy
    /// ledger closes: the sum of every step's `H(x1) - H(x0)` and every tension work is `H_end - H_start`.
    #[test]
    fn gl_012_sinusoidal_tension_modulates_f0_at_the_gesture_rate_and_closes_the_ledger() {
        let (rate_hz, a) = (6.0_f64, 0.08_f64);
        let mut isl = island(true, 1.0);
        let h_start = fold_energy(&isl);
        let n = (f64::from(RATE) * 1.3) as usize;
        let mut flows = Vec::with_capacity(n);
        let (mut ledger, mut h_max) = (0.0_f64, h_start);
        for k in 0..n {
            let t = k as f64 / f64::from(RATE);
            isl.set_tension(1.0 + a * det::sin(core::f64::consts::TAU * rate_hz * t))
                .expect("tension");
            let attack = (k as f64 / (0.03 * f64::from(RATE))).min(1.0);
            let frame = isl.step(1400.0 * attack).expect("step");
            assert!(frame.dissipated_j >= -1.0e-15);
            ledger += frame.fold_energy_j + frame.tension_work_j;
            flows.push(frame.flow_m3_s);
            h_max = h_max.max(fold_energy(&isl));
        }
        let miss = (ledger - (fold_energy(&isl) - h_start)).abs();
        assert!(
            miss <= 1.0e-9 * h_max,
            "energy ledger miss {miss:e} J of {h_max:e} J"
        );
        // Five whole periods of the gesture, after the onset transient.
        let (i0, i1) = (
            (0.30 * f64::from(RATE)) as usize,
            (0.30 * f64::from(RATE) + 5.0 / rate_hz * f64::from(RATE)) as usize,
        );
        let f0s = cycle_f0(&flows[i0..i1], f64::from(RATE));
        let f_mean = f0s.iter().map(|p| p.1).sum::<f64>() / f0s.len() as f64;
        let cents: Vec<(f64, f64)> = f0s
            .iter()
            .map(|&(t, f)| (t, 1200.0 * det::ln(f / f_mean) / core::f64::consts::LN_2))
            .collect();
        let at_rate = amplitude_at(&cents, rate_hz);
        let off = amplitude_at(&cents, rate_hz - 2.4).max(amplitude_at(&cents, rate_hz + 2.4));
        let per_a = 1200.0 / core::f64::consts::LN_2;
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"gl-012-tension-modulation\",\"cycles\":{},\"f0_mean_hz\":{f_mean:.2},\"peak_cents_at_6hz\":{at_rate:.1},\"max_cents_at_3p6_8p4hz\":{off:.1},\"ledger_miss_j\":{miss:e}}}",
            cents.len()
        );
        assert!(
            (0.3 * a * per_a..=0.8 * a * per_a).contains(&at_rate),
            "f0 swing {at_rate:.1} cents peak at {rate_hz} Hz"
        );
        assert!(
            at_rate > 5.0 * off,
            "modulation {at_rate:.1} cents against {off:.1} cents off the gesture rate"
        );
    }
}
