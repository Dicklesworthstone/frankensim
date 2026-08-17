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
//! BOUNDARY: no formant filters, no vocoder, no glottal-waveform
//! playback — the waveform must EMERGE from the valve + tract loop.

use crate::driving_point::{DrivingPointError, characteristic_line};
use crate::unilateral_contact::slit_contact_force;
use fs_dcontact::{DContactError, Obstacle};
use fs_duct::{Duct, Termination};
use fs_material::gas::GasState;
use fs_math::det;
use fs_phs::{PhsError, PortHamiltonian, bernoulli_volume_flow, mass_spring_damper};
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
}

enum Valve {
    OneDof {
        phs: PortHamiltonian,
        x: Vec<f64>,
    },
    TwoMass {
        phs: PortHamiltonian,
        x: Vec<f64>,
    },
}

/// A glottal island against a live tract characteristic line.
pub struct GlottalIsland {
    valve: Valve,
    card: FoldCard,
    line: DelayedFilter,
    zc_flow: f64,
    rho: f64,
    dt: f64,
    obstacle: Obstacle,
    p_minus: f64,
    p_plus_prev: f64,
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
        card.validate()?;
        let inlet_r = tract
            .segments
            .first()
            .and_then(|s| match *s {
                fs_duct::Segment::Cylinder { radius, .. } => Some(radius),
                fs_duct::Segment::Cone { inlet_radius, .. } => Some(inlet_radius),
                fs_duct::Segment::ToneHole { .. } => None,
            })
            .ok_or(GlottisError::Invalid {
                what: "tract needs a leading tube segment",
            })?;
        let area = core::f64::consts::PI * inlet_r * inlet_r;
        let zc_flow = gas.density * gas.sound_speed / area;
        let line = characteristic_line(
            tract,
            gas,
            termination,
            sample_rate_hz,
            8192,
            zc_flow * area, // plane-wave specific impedance for the FIR
            None,
        )?;
        let valve = if two_mass {
            Valve::TwoMass {
                phs: two_mass_phs(&card)?,
                x: vec![0.0; 4],
            }
        } else {
            let m = card.mass_lower_kg + card.mass_upper_kg;
            let k = card.stiffness_lower_n_m + card.stiffness_upper_n_m;
            let c = 2.0 * card.damping_ratio * det::sqrt(k * m);
            Valve::OneDof {
                phs: mass_spring_damper(m, k, c)?,
                x: vec![0.0; 2],
            }
        };
        // Collision: the S-H convention triples the spring at closure;
        // chi is the reed-lay internal-loss pattern (authored).
        let obstacle = crate::unilateral_contact::slit_lay(3.0 * card.stiffness_lower_n_m, 2.0)
            .and_then(|o| o.with_internal_loss(0.1))?;
        Ok(Self {
            valve,
            card,
            line,
            zc_flow,
            rho: gas.density,
            dt: 1.0 / f64::from(sample_rate_hz),
            obstacle,
            p_minus: 0.0,
            p_plus_prev: 0.0,
        })
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

    /// Advance one audio sample under subglottal pressure `p_sub_pa`.
    ///
    /// # Errors
    /// pHS step or line refusal (non-finite state).
    pub fn step(&mut self, p_sub_pa: f64) -> Result<GlottalFrame, GlottisError> {
        let gap = self.gap_m();
        let p_supra = 2.0 * self.p_minus + self.zc_flow * 0.0; // provisional
        // Flow from the PREVIOUS gap against the incident wave; then
        // the junction pressure includes the flow's own loading.
        let dp0 = p_sub_pa - (2.0 * self.p_minus);
        let mut flow = bernoulli_volume_flow(self.card.fold_length_m, gap.max(0.0), dp0, self.rho);
        // One fixed-point refinement of the junction coupling:
        // p_supra = 2 p_minus + Zc U, and the transglottal drop uses it.
        for _ in 0..3 {
            let p_supra_now = 2.0 * self.p_minus + self.zc_flow * flow;
            let dp = p_sub_pa - p_supra_now;
            flow = bernoulli_volume_flow(self.card.fold_length_m, gap.max(0.0), dp, self.rho);
        }
        let p_supra = if gap > 0.0 {
            2.0 * self.p_minus + self.zc_flow * flow
        } else {
            2.0 * self.p_minus
        };
        let dp = p_sub_pa - p_supra;
        let face = self.card.fold_length_m * self.card.thickness_m;
        let record = match &mut self.valve {
            Valve::OneDof { phs, x } => {
                // Pressure force + collision on the single coordinate.
                let h = self.card.rest_gap_m + x[0];
                let f_col = slit_contact_force(&self.obstacle, h)?;
                let f = face * dp + f_col;
                let rec = fs_phs::step(phs, x, &[f], self.dt)?;
                *x = rec.x.clone();
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
                let recovery = if h1 > 0.0 && hmin > 0.0 {
                    1.0 - (hmin / h1) * (hmin / h1)
                } else {
                    0.0
                };
                let drive1 = if h1 > 0.0 && h2 > 0.0 {
                    face * dp * (1.0 - recovery)
                } else if h1 > 0.0 {
                    // Upper closed: the lower face sees the full
                    // subglottal pressure (stagnation).
                    face * p_sub_pa
                } else {
                    0.0
                };
                let f1 = drive1 + slit_contact_force(&self.obstacle, h1)?;
                let f2 = slit_contact_force(&self.obstacle, h2)?;
                let rec = fs_phs::step(phs, x, &[f1, f2], self.dt)?;
                *x = rec.x.clone();
                rec
            }
        };
        // Scatter into the tract: p+ = p_supra - p-, then advance the
        // line one sample to produce the next incident wave.
        let p_plus = p_supra - self.p_minus;
        self.p_minus = self
            .line
            .push(self.p_plus_prev)
            .map_err(|_| GlottisError::Invalid {
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
#[must_use]
pub fn glottal_qois(flow: &[f64], gaps: &[f64], rate: f64) -> GlottalQois {
    // Periods from positive-going gap-opening instants.
    let mut openings = Vec::new();
    for k in 1..gaps.len() {
        if gaps[k] > 0.0 && gaps[k - 1] <= 0.0 {
            openings.push(k);
        }
    }
    let (f0_hz, jitter, open_quotient) = if openings.len() >= 4 {
        let periods: Vec<f64> = openings.windows(2).map(|w| (w[1] - w[0]) as f64).collect();
        let mean = periods.iter().sum::<f64>() / periods.len() as f64;
        let var = periods.iter().map(|p| (p - mean) * (p - mean)).sum::<f64>()
            / periods.len() as f64;
        let open: f64 = gaps[openings[0]..*openings.last().expect("nonempty")]
            .iter()
            .map(|&g| f64::from(u8::from(g > 0.0)))
            .sum();
        (
            rate / mean,
            var.sqrt() / mean,
            open / (openings.last().expect("nonempty") - openings[0]) as f64,
        )
    } else {
        (0.0, f64::NAN, 1.0)
    };
    // Harmonic magnitudes 1..6 by Goertzel-style projection.
    let slope = if f0_hz > 0.0 {
        let mut levels = Vec::new();
        for n in 1..=6 {
            let omega = core::f64::consts::TAU * f0_hz * n as f64 / rate;
            let (mut re, mut im) = (0.0f64, 0.0f64);
            for (k, &v) in flow.iter().enumerate() {
                re += v * det::cos(omega * k as f64);
                im -= v * det::sin(omega * k as f64);
            }
            levels.push(20.0 * det::ln((re * re + im * im).sqrt().max(1e-30)) / core::f64::consts::LN_10 / 1.0);
        }
        // dB per octave via a log2-frequency regression.
        let n = levels.len() as f64;
        let xs: Vec<f64> = (1..=6).map(|k| det::ln(f64::from(k)) / core::f64::consts::LN_2).collect();
        let mx = xs.iter().sum::<f64>() / n;
        let my = levels.iter().sum::<f64>() / n;
        let num: f64 = xs.iter().zip(&levels).map(|(x, y)| (x - mx) * (y - my)).sum();
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

    const RATE: u32 = 48_000;

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
        GlottalIsland::new(card, two_mass, &tract(), &air(), Termination::IdealOpen, RATE)
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
        let ac = |seg: &[f64]| -> f64 {
            let mean = seg.iter().sum::<f64>() / seg.len() as f64;
            (seg.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / seg.len() as f64).sqrt()
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
        while p < 6000.0 {
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
        for w in rows.chunks(3) {
            assert!(
                w[0].2 <= w[2].2,
                "onset must not FALL with stiffness ({w:?})"
            );
        }
        // Below threshold: silent (the detector, not wishful thinking).
        let mut soft = island(false, 1.0);
        assert!(!speaks(&mut soft, 40.0, 0.4).0, "sub-threshold must be silent");
        for (tm, ks, p) in &rows {
            println!(
                "{{\"suite\":\"fs-couple\",\"case\":\"gl-001-onset\",\"two_mass\":{tm},\
                 \"k_scale\":{ks},\"onset_pa\":{p}}}"
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
        // Phase via cross-correlation argmax over +-2 ms.
        let seg1 = &q1[n - tail..];
        let seg2 = &q2[n - tail..];
        let max_lag = (0.002 * f64::from(RATE)) as isize;
        let mut best = (0isize, f64::MIN);
        for lag in -max_lag..=max_lag {
            let mut acc = 0.0;
            for i in 0..seg1.len() {
                let j = i as isize + lag;
                if j >= 0 && (j as usize) < seg2.len() {
                    acc += seg1[i] * seg2[j as usize];
                }
            }
            if acc > best.1 {
                best = (lag, acc);
            }
        }
        assert!(
            best.0 > 0,
            "the upper mass must LAG the lower (best lag {} samples)",
            best.0
        );
        // Flow skew: at the limit cycle the steepest descent must be
        // steeper than the steepest rise (closing faster than opening).
        let dseg: Vec<f64> = flows[n - tail..].windows(2).map(|w| w[1] - w[0]).collect();
        let steepest_rise = dseg.iter().copied().fold(f64::MIN, f64::max);
        let steepest_fall = dseg.iter().copied().fold(f64::MAX, f64::min);
        assert!(
            -steepest_fall > steepest_rise,
            "the flow must close faster than it opens ({steepest_fall:.3e} vs {steepest_rise:.3e})"
        );
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"gl-003-mechanism\",\"verdict\":\"pass\",\
             \"upper_lag_samples\":{},\"skew_ratio\":{:.2}}}",
            best.0,
            -steepest_fall / steepest_rise
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
            let mut flows = Vec::with_capacity(n);
            let mut gaps = Vec::with_capacity(n);
            for k in 0..n {
                let attack = (k as f64 / (0.03 * f64::from(RATE))).min(1.0);
                let frame = isl.step(1400.0 * attack).expect("step");
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
}
