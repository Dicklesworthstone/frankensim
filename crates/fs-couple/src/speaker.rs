//! Thiele–Small driver + sealed cabinet (music bead
//! `frankensim-music-v8-root-3ez8g.9.4`): the last stage of the
//! electric chain — amp voltage in, radiated pressure out. Everything
//! is existing parts composed (D23): the Bl coupling is an
//! ANTISYMMETRIC J entry (a gyrator IS a Dirac coupling — voltage =
//! Bl·velocity and force = Bl·current come from one off-diagonal pair,
//! power-exact by construction); the suspension and moving mass are
//! the msd storage; the sealed box is EXTRA quadratic stiffness
//! `Sd² ρc² / V` on the cone coordinate (which is exactly why the
//! classic TS resonance-shift algebra holds); the cabinet panel is one
//! fs-plate-sourced mode coupled through the SHARED cavity term
//! `(Sd x − Sp x_p)² ρc²/2V` — the "cabinet sound" as physics,
//! X-Consist with mode truncation DISCLOSED (one panel mode in v1).
//!
//! CARD PROVENANCE: [`TsCard::datasheet_class_6p5`] is an AUTHORED
//! datasheet-class parameter set (Estimate; a manufacturer-sheet
//! ingest upgrades it without changing the runtime shape). The
//! radiation load is the low-`ka` baffled-piston approximation
//! (frequency-independent resistance + added mass — disclosed; the
//! full piston impedance is the recorded refinement). EXCURSION
//! HONESTY: cone travel beyond the card's `x_max` REFUSES — the
//! suspension nonlinearity is not modeled in v1 and silent linear
//! extrapolation presented as forte is exactly the lie this refusal
//! exists to prevent. FORBIDDEN and absent: cabinet IR packs — the
//! cabinet is charts and cards, never a convolution.

use fs_phs::{PhsError, PortHamiltonian, QuadraticStorage, StepRecord};

/// Thiele–Small parameter card (datasheet-class, with provenance).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TsCard {
    /// Voice-coil DC resistance [ohm].
    pub re_ohm: f64,
    /// Voice-coil inductance [H].
    pub le_h: f64,
    /// Force factor Bl [T m].
    pub bl: f64,
    /// Moving mass (cone + coil, without air load) [kg].
    pub mms_kg: f64,
    /// Suspension compliance [m/N].
    pub cms_m_per_n: f64,
    /// Mechanical loss [N s/m].
    pub rms_n_s_m: f64,
    /// Effective piston area [m^2].
    pub sd_m2: f64,
    /// Linear-excursion limit [m].
    pub x_max_m: f64,
}

impl TsCard {
    /// An authored 6.5-inch-woofer-class card (Estimate: typical
    /// datasheet magnitudes, not a specific product).
    #[must_use]
    pub fn datasheet_class_6p5() -> TsCard {
        TsCard {
            re_ohm: 6.4,
            le_h: 0.5e-3,
            bl: 7.5,
            mms_kg: 12.0e-3,
            cms_m_per_n: 1.2e-3,
            rms_n_s_m: 1.2,
            sd_m2: 1.37e-2,
            x_max_m: 4.0e-3,
        }
    }

    fn validate(&self) -> Result<(), SpeakerError> {
        for (v, what) in [
            (self.re_ohm, "Re"),
            (self.le_h, "Le"),
            (self.bl, "Bl"),
            (self.mms_kg, "Mms"),
            (self.cms_m_per_n, "Cms"),
            (self.rms_n_s_m, "Rms"),
            (self.sd_m2, "Sd"),
            (self.x_max_m, "x_max"),
        ] {
            if !(v.is_finite() && v > 0.0) {
                return Err(SpeakerError::Invalid { what });
            }
        }
        Ok(())
    }
}

/// One cabinet-panel mode (fs-plate-sourced frequency/mass/area).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PanelMode {
    /// Modal frequency in vacuo [Hz].
    pub frequency_hz: f64,
    /// Modal mass [kg].
    pub mass_kg: f64,
    /// Effective volume-displacement area of the mode [m^2].
    pub area_m2: f64,
    /// Modal damping ratio.
    pub damping_ratio: f64,
}

/// Typed refusals.
#[derive(Debug)]
pub enum SpeakerError {
    /// Bad card/box parameter, by name.
    Invalid {
        /// What.
        what: &'static str,
    },
    /// Cone excursion left the card's linear range — the suspension
    /// nonlinearity is NOT modeled; refusing beats lying.
    ExcursionExceeded {
        /// Excursion at refusal [m].
        excursion_m: f64,
        /// The card's limit [m].
        x_max_m: f64,
    },
    /// Underlying pHS refusal.
    Phs(PhsError),
}

impl From<PhsError> for SpeakerError {
    fn from(e: PhsError) -> Self {
        SpeakerError::Phs(e)
    }
}

/// The composed driver (+ optional sealed box and panel mode).
pub struct TsDriver {
    phs: PortHamiltonian,
    x: Vec<f64>,
    card: TsCard,
    rho: f64,
    /// Index of the cone-displacement coordinate.
    ix_cone: usize,
    has_panel: bool,
}

/// Air constants for the acoustic terms (20 C).
const RHO: f64 = 1.204;
const C_SOUND: f64 = 343.0;

impl TsDriver {
    /// Compose the driver. `box_volume_m3 = None` is free air;
    /// `panel = Some(..)` adds one cabinet-panel mode coupled through
    /// the shared cavity (requires a box).
    ///
    /// # Errors
    /// Card/box refusals by name; pHS admission.
    pub fn new(
        card: TsCard,
        box_volume_m3: Option<f64>,
        panel: Option<PanelMode>,
    ) -> Result<Self, SpeakerError> {
        card.validate()?;
        if let Some(v) = box_volume_m3
            && !(v.is_finite() && v > 0.0)
        {
            return Err(SpeakerError::Invalid {
                what: "box volume must be positive",
            });
        }
        if panel.is_some() && box_volume_m3.is_none() {
            return Err(SpeakerError::Invalid {
                what: "a panel mode needs a box to couple through",
            });
        }
        // Low-ka baffled-piston load: added mass 8 rho a^3 / 3 and a
        // small frequency-independent resistance class (disclosed
        // approximation).
        let a_eff = (card.sd_m2 / core::f64::consts::PI).sqrt();
        let m_air = 8.0 * RHO * a_eff.powi(3) / 3.0;
        let m_tot = card.mms_kg + m_air;
        let r_rad = RHO * C_SOUND * card.sd_m2 * 0.05;
        let k_susp = 1.0 / card.cms_m_per_n;
        // State: [phi_Le, x_cone, p_cone] (+ [x_p, p_p] with a panel).
        let n = if panel.is_some() { 5 } else { 3 };
        let mut q = vec![0.0; n * n];
        q[0] = 1.0 / card.le_h;
        q[n + 1] = k_susp;
        q[2 * n + 2] = 1.0 / m_tot;
        if let Some(vb) = box_volume_m3 {
            let k_box_scale = RHO * C_SOUND * C_SOUND / vb;
            if let Some(p) = panel {
                if !(p.frequency_hz > 0.0 && p.mass_kg > 0.0 && p.area_m2 > 0.0) {
                    return Err(SpeakerError::Invalid {
                        what: "panel mode values must be positive",
                    });
                }
                // Cavity term (Sd x - Sp xp)^2 * rho c^2 / 2V couples
                // cone and panel; panel's own stiffness/mass added.
                q[n + 1] += k_box_scale * card.sd_m2 * card.sd_m2;
                q[n + 3] -= k_box_scale * card.sd_m2 * p.area_m2;
                q[3 * n + 1] -= k_box_scale * card.sd_m2 * p.area_m2;
                let w_p = core::f64::consts::TAU * p.frequency_hz;
                q[3 * n + 3] = p.mass_kg * w_p * w_p + k_box_scale * p.area_m2 * p.area_m2;
                q[4 * n + 4] = 1.0 / p.mass_kg;
            } else {
                q[n + 1] += k_box_scale * card.sd_m2 * card.sd_m2;
            }
        }
        let mut j = vec![0.0; n * n];
        let set = |j: &mut Vec<f64>, r_: usize, c: usize, v: f64| {
            j[r_ * n + c] += v;
            j[c * n + r_] -= v;
        };
        // Gyrator: phi row gets -Bl * v (col p), p row gets +Bl * i.
        set(&mut j, 0, 2, -card.bl);
        // Kinematics: x_dot = v; p row gets -dH/dx by antisymmetry.
        set(&mut j, 1, 2, 1.0);
        let mut r = vec![0.0; n * n];
        r[0] = card.re_ohm;
        r[2 * n + 2] = card.rms_n_s_m + r_rad;
        if let Some(p) = panel {
            set(&mut j, 3, 4, 1.0);
            let w_p = core::f64::consts::TAU * p.frequency_hz;
            r[4 * n + 4] = 2.0 * p.damping_ratio * p.mass_kg * w_p;
        }
        let mut g = vec![0.0; n];
        g[0] = 1.0;
        let storage = Box::new(QuadraticStorage::new(q, n)?);
        let phs = PortHamiltonian::new(n, 1, j, r, g, storage)?;
        Ok(Self {
            phs,
            x: vec![0.0; n],
            card,
            rho: RHO,
            ix_cone: 1,
            has_panel: panel.is_some(),
        })
    }

    /// One step under terminal voltage `u_v`; returns the record and
    /// the radiated far-field pressure factor `rho Sd a / (4 pi r)` at
    /// `r = 1 m` (compact-source observer).
    ///
    /// # Errors
    /// [`SpeakerError::ExcursionExceeded`] past the card's linear
    /// range; pHS step refusals.
    pub fn step(&mut self, u_v: f64, dt: f64) -> Result<(StepRecord, f64), SpeakerError> {
        let v_before = self.x[2]; // momentum, for acceleration estimate
        let rec = fs_phs::step(&self.phs, &self.x, &[u_v], dt)?;
        self.x.clone_from(&rec.x);
        let excursion = self.x[self.ix_cone];
        if excursion.abs() > self.card.x_max_m {
            return Err(SpeakerError::ExcursionExceeded {
                excursion_m: excursion,
                x_max_m: self.card.x_max_m,
            });
        }
        // a = dv/dt from the momentum update (Mms+air mass folded in).
        let a_eff = (self.card.sd_m2 / core::f64::consts::PI).sqrt();
        let m_tot = self.card.mms_kg + 8.0 * self.rho * a_eff.powi(3) / 3.0;
        let accel = (self.x[2] - v_before) / (m_tot * dt.max(1e-300));
        let p_1m = self.rho * self.card.sd_m2 * accel / (4.0 * core::f64::consts::PI);
        Ok((rec, p_1m))
    }

    /// Voice-coil current [A].
    #[must_use]
    pub fn coil_current_a(&self) -> f64 {
        self.x[0] / self.card.le_h
    }

    /// Cone excursion [m].
    #[must_use]
    pub fn excursion_m(&self) -> f64 {
        self.x[self.ix_cone]
    }

    /// Whether a panel mode is composed in.
    #[must_use]
    pub fn has_panel(&self) -> bool {
        self.has_panel
    }
}

#[cfg(test)]
mod speaker_tests {
    use super::*;
    use fs_math::det;

    const RATE: f64 = 96_000.0;

    /// Steady-state |Z| = |V/I| at frequency `f` by driving sinusoids.
    fn impedance_at(driver_factory: &dyn Fn() -> TsDriver, f: f64, amp_v: f64) -> f64 {
        let mut d = driver_factory();
        let dt = 1.0 / RATE;
        let n = (24.0 / f / dt) as usize;
        let (mut i_re, mut i_im, mut count) = (0.0f64, 0.0f64, 0usize);
        for k in 0..n {
            let ph = core::f64::consts::TAU * f * k as f64 * dt;
            let u = amp_v * det::sin(ph);
            d.step(u, dt).expect("step");
            if k > 2 * n / 3 {
                let i = d.coil_current_a();
                i_re += i * det::sin(ph);
                i_im += i * det::cos(ph);
                count += 1;
            }
        }
        let i_amp = 2.0 * (i_re * i_re + i_im * i_im).sqrt() / count as f64;
        amp_v / i_amp
    }

    fn card() -> TsCard {
        TsCard::datasheet_class_6p5()
    }

    fn fs_free_hz() -> f64 {
        let c = card();
        let a_eff = (c.sd_m2 / core::f64::consts::PI).sqrt();
        let m_tot = c.mms_kg + 8.0 * RHO * a_eff.powi(3) / 3.0;
        1.0 / (core::f64::consts::TAU * (m_tot * c.cms_m_per_n).sqrt())
    }

    fn peak_frequency(factory: &dyn Fn() -> TsDriver, lo: f64, hi: f64) -> f64 {
        let mut best = (lo, 0.0f64);
        let mut f = lo;
        while f <= hi {
            let z = impedance_at(factory, f, 0.5);
            if z > best.1 {
                best = (f, z);
            }
            f *= 1.03;
        }
        best.0
    }

    #[test]
    fn ts_001_sealed_box_shift_matches_the_ts_algebra() {
        // The classic: fc = fs sqrt(1 + Vas/Vb) with
        // Vas = rho c^2 Cms Sd^2 — the analytic oracle from the same
        // card numbers, never from the stepped system.
        let c = card();
        let vas = RHO * C_SOUND * C_SOUND * c.cms_m_per_n * c.sd_m2 * c.sd_m2;
        let vb = 0.020;
        let fs = fs_free_hz();
        let fc_expected = fs * (1.0 + vas / vb).sqrt();
        let free = || TsDriver::new(card(), None, None).expect("free");
        let boxed = || TsDriver::new(card(), Some(0.020), None).expect("boxed");
        let fs_meas = peak_frequency(&free, 20.0, 120.0);
        let fc_meas = peak_frequency(&boxed, 30.0, 250.0);
        let rel_fs = (fs_meas - fs).abs() / fs;
        let rel_fc = (fc_meas - fc_expected).abs() / fc_expected;
        assert!(rel_fs < 0.04, "free-air fs {fs_meas:.1} vs {fs:.1} Hz");
        assert!(
            rel_fc < 0.04,
            "sealed fc {fc_meas:.1} vs TS algebra {fc_expected:.1} Hz (Vas {vas:.4} m^3)"
        );
        assert!(
            fc_meas > 1.2 * fs_meas,
            "the box must stiffen the resonance"
        );
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"ts-001-box-shift\",\"fs\":{fs_meas:.1},\
             \"fc\":{fc_meas:.1},\"fc_analytic\":{fc_expected:.1},\"vas_l\":{:.1}}}",
            vas * 1e3
        );
    }

    #[test]
    fn ts_002_impedance_curve_shows_the_motional_peak() {
        // The electrical impedance vs the analytic TS curve
        // Z = Re + jwLe + Bl^2 / (Rms+Rrad + j(w M - 1/(w Cms))):
        // compared at frequencies across the motional peak.
        let c = card();
        let a_eff = (c.sd_m2 / core::f64::consts::PI).sqrt();
        let m_tot = c.mms_kg + 8.0 * RHO * a_eff.powi(3) / 3.0;
        let r_rad = RHO * C_SOUND * c.sd_m2 * 0.05;
        let factory = || TsDriver::new(card(), None, None).expect("driver");
        let fs = fs_free_hz();
        for f_ratio in [0.5f64, 1.0, 2.0, 6.0] {
            let f = f_ratio * fs;
            let w = core::f64::consts::TAU * f;
            let mech_im = w * m_tot - 1.0 / (w * c.cms_m_per_n);
            let mech = fs_math::c64::C64::new(c.rms_n_s_m + r_rad, mech_im);
            let motional = fs_math::c64::C64::new(c.bl * c.bl, 0.0) * mech.recip();
            let z_analytic = (fs_math::c64::C64::new(c.re_ohm, w * c.le_h) + motional).abs();
            let z = impedance_at(&factory, f, 0.5);
            let rel = (z - z_analytic).abs() / z_analytic;
            assert!(
                rel < 0.05,
                "|Z| at {f_ratio}x fs: {z:.2} vs analytic {z_analytic:.2} (rel {rel:.3})"
            );
            println!(
                "{{\"suite\":\"fs-couple\",\"case\":\"ts-002-impedance\",\"f_ratio\":{f_ratio},\
                 \"z\":{z:.2},\"z_analytic\":{z_analytic:.2}}}"
            );
        }
    }

    #[test]
    fn ts_003_panel_mode_colors_the_response_and_is_disclosed() {
        // The cabinet-panel contrast: with one panel mode coupled
        // through the cavity, the electrical impedance near the panel
        // resonance differs measurably from the rigid box; clamping
        // the panel (omitting it) removes the feature. X-Consist with
        // one retained mode, disclosed.
        let panel = PanelMode {
            frequency_hz: 180.0,
            mass_kg: 0.35,
            area_m2: 0.06,
            damping_ratio: 0.02,
        };
        let rigid = || TsDriver::new(card(), Some(0.020), None).expect("rigid");
        let flexed = || TsDriver::new(card(), Some(0.020), Some(panel)).expect("panel");
        let mut worst = 0.0f64;
        let mut far = 0.0f64;
        for f in [150.0f64, 165.0, 178.0, 182.0, 195.0, 400.0, 800.0] {
            let zr = impedance_at(&rigid, f, 0.5);
            let zf = impedance_at(&flexed, f, 0.5);
            let dev = (zf - zr).abs() / zr;
            if (f - panel.frequency_hz).abs() < 20.0 {
                worst = worst.max(dev);
            } else if f > 300.0 {
                far = far.max(dev);
            }
            println!(
                "{{\"suite\":\"fs-couple\",\"case\":\"ts-003-panel\",\"f\":{f},\"z_rigid\":{zr:.2},\
                 \"z_panel\":{zf:.2},\"dev\":{dev:.4}}}"
            );
        }
        assert!(
            worst > 5.0 * far.max(1e-6),
            "the panel must color the response NEAR its resonance \
             (near dev {worst:.4} vs far {far:.4})"
        );
    }

    #[test]
    fn ts_004_excursion_refuses_past_the_card_limit() {
        // Forte honesty: driving hard enough to exceed x_max REFUSES
        // (the suspension nonlinearity is not modeled); a gentle drive
        // never trips it.
        let mut d = TsDriver::new(card(), None, None).expect("driver");
        let dt = 1.0 / RATE;
        let fs = fs_free_hz();
        let mut refused = false;
        for k in 0..(2.0 * RATE / fs) as usize * 40 {
            let u = 60.0 * det::sin(core::f64::consts::TAU * fs * k as f64 * dt);
            match d.step(u, dt) {
                Ok(_) => {}
                Err(SpeakerError::ExcursionExceeded {
                    excursion_m,
                    x_max_m,
                }) => {
                    assert!(excursion_m.abs() > x_max_m);
                    refused = true;
                    break;
                }
                Err(e) => panic!("unexpected {e:?}"),
            }
        }
        assert!(refused, "a 60 V drive at resonance must exceed 4 mm");
        let mut gentle = TsDriver::new(card(), None, None).expect("driver");
        for k in 0..(RATE * 0.2) as usize {
            let u = 0.5 * det::sin(core::f64::consts::TAU * fs * k as f64 * dt);
            gentle.step(u, dt).expect("gentle stays linear");
        }
        assert!(gentle.excursion_m().abs() <= card().x_max_m);
        println!("{{\"suite\":\"fs-couple\",\"case\":\"ts-004-excursion\",\"verdict\":\"pass\"}}");
    }

    #[test]
    fn ts_005_refusals_and_ledger() {
        let mut bad = card();
        bad.bl = -1.0;
        assert!(matches!(
            TsDriver::new(bad, None, None),
            Err(SpeakerError::Invalid { what: "Bl" })
        ));
        assert!(matches!(
            TsDriver::new(card(), Some(-1.0), None),
            Err(SpeakerError::Invalid { .. })
        ));
        assert!(matches!(
            TsDriver::new(
                card(),
                None,
                Some(PanelMode {
                    frequency_hz: 180.0,
                    mass_kg: 0.3,
                    area_m2: 0.05,
                    damping_ratio: 0.02
                })
            ),
            Err(SpeakerError::Invalid { .. })
        ));
        // Ledger: driven steps keep the supply audit clean.
        let mut d = TsDriver::new(card(), Some(0.02), None).expect("driver");
        let dt = 1.0 / RATE;
        let mut worst = 0.0f64;
        for k in 0..4000 {
            let u = 2.0 * det::sin(core::f64::consts::TAU * 100.0 * k as f64 * dt);
            let (rec, _) = d.step(u, dt).expect("step");
            worst = worst.max((rec.delta_h + rec.dissipated - rec.supplied).abs());
        }
        assert!(worst < 1.0e-9, "supply defect {worst:.3e}");
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"ts-005-ledger\",\"worst_defect\":{worst:.3e}}}"
        );
    }
}
