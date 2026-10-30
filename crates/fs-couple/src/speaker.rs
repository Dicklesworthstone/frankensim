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
//!
//! [`MagneticTsDriver`] additionally accepts bounded caller-supplied Bl(x) and
//! L(x) profiles. It retains this same cone/cabinet model and uses one magnetic
//! energy for nonlinear force, reluctance force and reciprocal back EMF. The
//! original [`TsDriver`] keeps its exact linear held-input propagation.

use fs_material::gas::{GasSpec, GasState};
use fs_phs::linear_exact::{
    LinearPropagator, LinearSplit, LinearSplitBudget, split_quadratic, step_split,
};
use fs_phs::{PhsError, PortHamiltonian, QuadraticStorage, StepRecord};

mod magnetic;
pub use magnetic::thermal;
pub use magnetic::{MagneticCoilSample, MagneticTsDriver, VoiceCoilMagnetics};

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
    /// The shared material owner rejected an inconsistent gas snapshot.
    Gas(fs_material::MaterialError),
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
    linear: LinearSplit,
    x: Vec<f64>,
    card: TsCard,
    exterior_gas: GasState,
    cavity_gas: Option<GasState>,
    moving_mass_kg: f64,
    /// Index of the cone-displacement coordinate.
    ix_cone: usize,
    has_panel: bool,
}

/// Declared room state for the driver's convenience constructor.
const SPEAKER_AIR_TEMPERATURE_K: f64 = 293.15;
const SPEAKER_AIR_PRESSURE_PA: f64 = 101_325.0;

fn declared_air() -> GasState {
    GasState::try_new(
        &GasSpec::dry_air_ussa1976(),
        SPEAKER_AIR_TEMPERATURE_K,
        SPEAKER_AIR_PRESSURE_PA,
    )
    .expect("declared dry-air room state is within GasState's validity window")
}

// Shared physical assembly; the two drivers differ only in their magnetic
// constitutive law and time integrator, never in cone/cavity/panel mechanics.
struct DriverAssembly {
    n: usize,
    q: Vec<f64>,
    j: Vec<f64>,
    r: Vec<f64>,
    g: Vec<f64>,
    moving_mass_kg: f64,
}

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
        Self::new_with_gas(card, box_volume_m3, panel, declared_air())
    }

    /// Compose the same driver in an explicitly supplied gas state.
    ///
    /// Density and sound speed from this state set radiation added mass,
    /// radiation resistance, sealed-box stiffness, panel/cavity coupling and
    /// observer pressure. This is the data-only material/composition swap path;
    /// callers may pass a source-resolved or mixed gas without changing driver
    /// code. The convenience constructor [`Self::new`] retains its declared
    /// room dry-air state.
    ///
    /// # Errors
    /// The same card/box/panel and pHS refusals as [`Self::new`], plus the
    /// shared gas-state coherence check. Changing a public pressure,
    /// temperature or composition-derived coefficient requires reevaluating
    /// the material law before constructing the driver.
    pub fn new_with_gas(
        card: TsCard,
        box_volume_m3: Option<f64>,
        panel: Option<PanelMode>,
        gas: GasState,
    ) -> Result<Self, SpeakerError> {
        Self::new_with_media(card, box_volume_m3, panel, gas, box_volume_m3.map(|_| gas))
    }

    /// Compose with independently resolved exterior and sealed-cavity gases.
    ///
    /// Exterior density and sound speed supply the original radiation load
    /// and microphone pressure. Only the cavity's bulk modulus `rho*c*c`
    /// supplies sealed-volume storage and reciprocal cone/panel coupling.
    /// Changing a sealed fill therefore cannot replace the surrounding gas.
    /// The existing compact radiation and spatially uniform, adiabatic cavity
    /// approximations remain; these are fixed gas snapshots, not gas transport,
    /// heat exchange, molecular relaxation or finite-amplitude compression.
    ///
    /// A cavity gas is required exactly when a box volume is supplied. Both
    /// regions must have the same static pressure because this model is
    /// linearized about zero displacement and contains no pressure-bias load.
    /// Temperature and composition may differ. Source authority stays with
    /// the caller's original material models; both complete snapshots are kept.
    ///
    /// # Errors
    /// Refuses incoherent media, missing/extra cavity gas, unequal static
    /// pressures and the original card, box, panel and pHS admission failures.
    pub fn new_with_media(
        card: TsCard,
        box_volume_m3: Option<f64>,
        panel: Option<PanelMode>,
        exterior_gas: GasState,
        cavity_gas: Option<GasState>,
    ) -> Result<Self, SpeakerError> {
        let DriverAssembly {
            n,
            q,
            j,
            r,
            g,
            moving_mass_kg,
        } = assemble_driver(
            card,
            box_volume_m3,
            panel,
            &exterior_gas,
            cavity_gas.as_ref(),
        )?;
        let storage = Box::new(QuadraticStorage::new(q, n)?);
        let phs = PortHamiltonian::new(n, 1, j, r, g, storage)?;
        let linear = split_quadratic(
            &phs,
            LinearSplitBudget {
                max_modes: n,
                max_omega_dt: core::f64::consts::PI,
            },
        )?;
        Ok(Self {
            phs,
            linear,
            x: vec![0.0; n],
            card,
            exterior_gas,
            cavity_gas,
            moving_mass_kg,
            ix_cone: 1,
            has_panel: panel.is_some(),
        })
    }

    /// Exact full coupled quadratic propagation under held terminal voltage.
    /// No midpoint phase warping is introduced by the coil/cone/panel owner.
    /// The record uses independently integrated loss and step-average current.
    /// This dense reference path is not an allocation-free audio callback.
    /// Returns the record and
    /// the radiated far-field pressure factor `rho Sd a / (4 pi r)` at
    /// `r = 1 m` (compact-source observer).
    ///
    /// # Errors
    /// [`SpeakerError::ExcursionExceeded`] past the card's linear
    /// range; pHS step refusals.
    pub fn step(&mut self, u_v: f64, dt: f64) -> Result<(StepRecord, f64), SpeakerError> {
        let v_before = self.x[2]; // momentum, for acceleration estimate
        let rec = step_split(
            &self.phs,
            &self.linear,
            LinearPropagator::ExactZoh,
            &self.x,
            &[u_v],
            dt,
        )?;
        self.x.clone_from(&rec.x);
        let excursion = self.x[self.ix_cone];
        if excursion.abs() > self.card.x_max_m {
            return Err(SpeakerError::ExcursionExceeded {
                excursion_m: excursion,
                x_max_m: self.card.x_max_m,
            });
        }
        // a = dv/dt from the momentum update (Mms+air mass folded in).
        let accel = (self.x[2] - v_before) / (self.moving_mass_kg * dt.max(1e-300));
        let p_1m =
            self.exterior_gas.density * self.card.sd_m2 * accel / (4.0 * core::f64::consts::PI);
        Ok((rec, p_1m))
    }

    /// Original surrounding medium used by radiation loading and observation.
    #[must_use]
    pub const fn exterior_gas(&self) -> &GasState {
        &self.exterior_gas
    }

    /// Original sealed fill used by cavity compression, absent in free air.
    #[must_use]
    pub const fn cavity_gas(&self) -> Option<&GasState> {
        self.cavity_gas.as_ref()
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

fn assemble_driver(
    card: TsCard,
    box_volume_m3: Option<f64>,
    panel: Option<PanelMode>,
    exterior_gas: &GasState,
    cavity_gas: Option<&GasState>,
) -> Result<DriverAssembly, SpeakerError> {
    card.validate()?;
    exterior_gas
        .validate_coherence()
        .map_err(SpeakerError::Gas)?;
    if box_volume_m3.is_some() != cavity_gas.is_some() {
        return Err(SpeakerError::Invalid {
            what: "a sealed box requires its own cavity gas; free air has no cavity gas",
        });
    }
    if let Some(gas) = cavity_gas {
        gas.validate_coherence().map_err(SpeakerError::Gas)?;
        if gas.pressure != exterior_gas.pressure {
            return Err(SpeakerError::Invalid {
                what: "cavity and exterior need equal static pressure; pressure-bias equilibrium is not modeled",
            });
        }
    }
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
    let m_air = 8.0 * exterior_gas.density * a_eff.powi(3) / 3.0;
    let m_tot = card.mms_kg + m_air;
    let r_rad = exterior_gas.density * exterior_gas.sound_speed * card.sd_m2 * 0.05;
    let k_susp = 1.0 / card.cms_m_per_n;
    // State: [phi_Le, x_cone, p_cone] (+ [x_p, p_p] with a panel).
    let n = if panel.is_some() { 5 } else { 3 };
    let mut q = vec![0.0; n * n];
    q[0] = 1.0 / card.le_h;
    q[n + 1] = k_susp;
    q[2 * n + 2] = 1.0 / m_tot;
    if let (Some(vb), Some(gas)) = (box_volume_m3, cavity_gas) {
        let k_box_scale = gas.density * gas.sound_speed * gas.sound_speed / vb;
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
    Ok(DriverAssembly {
        n,
        q,
        j,
        r,
        g,
        moving_mass_kg: m_tot,
    })
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
        let m_tot = c.mms_kg + 8.0 * declared_air().density * a_eff.powi(3) / 3.0;
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
        let air = declared_air();
        let vas =
            air.density * air.sound_speed * air.sound_speed * c.cms_m_per_n * c.sd_m2 * c.sd_m2;
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
        let air = declared_air();
        let a_eff = (c.sd_m2 / core::f64::consts::PI).sqrt();
        let m_tot = c.mms_kg + 8.0 * air.density * a_eff.powi(3) / 3.0;
        let r_rad = air.density * air.sound_speed * c.sd_m2 * 0.05;
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
            let u = 2.0 * det::sin(core::f64::consts::TAU * 100.0 * f64::from(k) * dt);
            let (rec, _) = d.step(u, dt).expect("step");
            worst = worst.max((rec.delta_h + rec.dissipated - rec.supplied).abs());
        }
        assert!(worst < 1.0e-9, "supply defect {worst:.3e}");
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"ts-005-ledger\",\"worst_defect\":{worst:.3e}}}"
        );
    }

    #[test]
    fn ts_006_admitted_gas_sets_every_acoustic_coefficient() {
        let cold = GasState::try_new(&GasSpec::dry_air_ussa1976(), 273.15, 101_325.0)
            .expect("admitted cold dry air");
        let warm = GasState::try_new(&GasSpec::dry_air_ussa1976(), 313.15, 101_325.0)
            .expect("admitted warm dry air");
        assert!(cold.density > warm.density, "cold air is denser");
        assert!(warm.sound_speed > cold.sound_speed, "warm air is faster");

        let c = card();
        let a_eff = (c.sd_m2 / core::f64::consts::PI).sqrt();
        let box_volume_m3 = 0.020;
        for gas in [cold, warm] {
            let driver = TsDriver::new_with_gas(c, Some(box_volume_m3), None, gas).expect("driver");
            let n = driver.phs.state_dim();
            assert_eq!(n, 3, "the no-panel driver has three states");

            let expected_m_tot = c.mms_kg + 8.0 * gas.density * a_eff.powi(3) / 3.0;
            let mut unit_cone_momentum = vec![0.0; n];
            unit_cone_momentum[2] = 1.0;
            let momentum_effort = driver.phs.effort(&unit_cone_momentum);
            assert_eq!(
                (1.0 / expected_m_tot).to_bits(),
                momentum_effort[2].to_bits(),
                "admitted density sets the cone mass"
            );

            let (_, r, _) = driver.phs.structure();
            let expected_radiation = gas.density * gas.sound_speed * c.sd_m2 * 0.05;
            let actual_radiation = r[2 * n + 2] - c.rms_n_s_m;
            assert!(
                (actual_radiation - expected_radiation).abs() <= expected_radiation.abs() * 1.0e-14,
                "admitted density and sound speed set radiation loss"
            );

            let mut unit_cone_displacement = vec![0.0; n];
            unit_cone_displacement[driver.ix_cone] = 1.0;
            let effort = driver.phs.effort(&unit_cone_displacement);
            let expected_box_stiffness =
                gas.density * gas.sound_speed * gas.sound_speed * c.sd_m2 * c.sd_m2 / box_volume_m3;
            let actual_box_stiffness = effort[driver.ix_cone] - 1.0 / c.cms_m_per_n;
            assert!(
                (actual_box_stiffness - expected_box_stiffness).abs()
                    <= expected_box_stiffness.abs() * 1.0e-14,
                "admitted density and sound speed set box stiffness"
            );
            assert_eq!(
                driver.exterior_gas().density.to_bits(),
                gas.density.to_bits(),
                "observer pressure retains the admitted density"
            );
        }

        let default_driver = TsDriver::new(c, Some(box_volume_m3), None).expect("default");
        assert_eq!(
            default_driver.exterior_gas().density.to_bits(),
            declared_air().density.to_bits(),
            "the public convenience wrapper retains its declared-air provenance"
        );
    }

    #[test]
    fn ts_007_exact_held_voltage_is_invariant_under_time_partition() {
        // G1: the complete coil/cone/cavity/panel state, with noncommuting
        // damping/coupling, must obey the constant-input exponential semigroup.
        let panel = PanelMode {
            frequency_hz: 5_000.0,
            mass_kg: 0.3,
            area_m2: 0.05,
            damping_ratio: 0.02,
        };
        for rate in [48_000, 96_000, 192_000] {
            let dt = 1.0 / f64::from(rate);
            let mut full = TsDriver::new(card(), Some(0.02), Some(panel)).unwrap();
            let mut half = TsDriver::new(card(), Some(0.02), Some(panel)).unwrap();
            let mut supplied_full = 0.0;
            let mut supplied_half = 0.0;
            let mut loss_full = 0.0;
            let mut loss_half = 0.0;
            for _ in 0..20 {
                let (a, _) = full.step(0.01, dt).unwrap();
                supplied_full += a.supplied;
                loss_full += a.dissipated;
                for _ in 0..2 {
                    let (b, _) = half.step(0.01, 0.5 * dt).unwrap();
                    supplied_half += b.supplied;
                    loss_half += b.dissipated;
                }
            }
            for (a, b) in full.x.iter().zip(&half.x) {
                let scale = a.abs().max(b.abs()).max(f64::MIN_POSITIVE);
                assert!(
                    (a - b).abs() <= 1e-12 * scale,
                    "held-voltage state partition differs at {rate} Hz"
                );
            }
            assert!((supplied_full - supplied_half).abs() <= 1e-12 * supplied_full.abs());
            assert!((loss_full - loss_half).abs() <= 1e-12 * loss_full.abs());
        }
    }

    #[test]
    fn gas_changes_must_reach_the_same_electrical_mechanical_acoustic_state() {
        // G0: stale pressure or transport must not silently coexist with the
        // old density/sound speed in cabinet storage and radiated pressure.
        let original = declared_air();
        let corruptions: [fn(&mut GasState); 5] = [
            |gas| gas.pressure *= 2.0,
            |gas| gas.temperature *= 1.1,
            |gas| gas.density *= 2.0,
            |gas| gas.sound_speed *= 1.1,
            |gas| gas.dynamic_viscosity *= 2.0,
        ];
        for corrupt in corruptions {
            let mut stale = original;
            corrupt(&mut stale);
            assert!(matches!(
                TsDriver::new_with_gas(card(), Some(0.02), None, stale),
                Err(SpeakerError::Gas(_))
            ));
        }

        let changed = GasState::try_new(
            &GasSpec::dry_air_ussa1976(),
            original.temperature,
            2.0 * original.pressure,
        )
        .unwrap();
        let driver = TsDriver::new_with_gas(card(), Some(0.02), None, changed).unwrap();
        assert_eq!(driver.exterior_gas().density, 2.0 * original.density);
    }

    #[test]
    fn separate_gases_drive_only_their_own_cavity_and_exterior_operators() {
        // G1/G3: independently check adiabatic pV storage, including its
        // reciprocal panel reaction, and exterior inertia/radiation scaling.
        let air = declared_air();
        let helium = GasState::try_new(
            &fs_material::gas_species::GasSpecies::Helium.spec(),
            air.temperature,
            air.pressure,
        )
        .unwrap();
        let c = card();
        let panel = PanelMode {
            frequency_hz: 90.0,
            mass_kg: 0.3,
            area_m2: 0.05,
            damping_ratio: 0.02,
        };
        let volume = 0.02;
        let baseline =
            TsDriver::new_with_media(c, Some(volume), Some(panel), air, Some(air)).unwrap();
        let mut filled =
            TsDriver::new_with_media(c, Some(volume), Some(panel), air, Some(helium)).unwrap();
        let outside =
            TsDriver::new_with_media(c, Some(volume), Some(panel), helium, Some(air)).unwrap();
        let displacement = [0.0, 2e-4, 0.0, -1e-4, 0.0];
        assert_eq!(
            baseline.phs.effort(&displacement),
            outside.phs.effort(&displacement)
        );
        assert_eq!(baseline.phs.structure(), filled.phs.structure());
        assert_eq!(baseline.moving_mass_kg, filled.moving_mass_kg);
        assert!(outside.moving_mass_kg < baseline.moving_mass_kg);
        assert_eq!(filled.exterior_gas(), &air);
        assert_eq!(filled.cavity_gas(), Some(&helium));

        let effort = filled.phs.effort(&displacement);
        let compression = c.sd_m2 * displacement[1] - panel.area_m2 * displacement[3];
        let k = helium.gamma * helium.pressure / volume;
        let expected_cone = displacement[1] / c.cms_m_per_n + k * c.sd_m2 * compression;
        let expected_panel =
            panel.mass_kg * (core::f64::consts::TAU * panel.frequency_hz).powi(2) * displacement[3]
                - k * panel.area_m2 * compression;
        assert!((effort[1] - expected_cone).abs() < 2e-14 * expected_cone.abs());
        assert!((effort[3] - expected_panel).abs() < 2e-14 * expected_panel.abs());
        assert!(effort[1] > baseline.phs.effort(&displacement)[1]);

        let mut baseline = baseline;
        let mut distinct_pressure = false;
        for _ in 0..120 {
            let previous_momentum = filled.x[2];
            let (record, pressure) = filled.step(0.02, 1.0 / 48_000.0).unwrap();
            let (_, reference) = baseline.step(0.02, 1.0 / 48_000.0).unwrap();
            let expected_pressure =
                air.density * c.sd_m2 * (record.x[2] - previous_momentum) * 48_000.0
                    / (filled.moving_mass_kg * 4.0 * core::f64::consts::PI);
            assert!(
                (pressure - expected_pressure).abs() < 1e-12 * expected_pressure.abs().max(1e-30)
            );
            assert!((record.delta_h + record.dissipated - record.supplied).abs() < 1e-12);
            distinct_pressure |= (pressure - reference).abs() > 1e-9;
        }
        assert!(
            distinct_pressure,
            "the fill must change solved motion and external pressure"
        );
    }

    #[test]
    fn equal_media_preserve_the_complete_legacy_driver() {
        let gas = declared_air();
        for volume in [None, Some(0.02)] {
            let mut old = TsDriver::new_with_gas(card(), volume, None, gas).unwrap();
            let mut explicit =
                TsDriver::new_with_media(card(), volume, None, gas, volume.map(|_| gas)).unwrap();
            for voltage in [0.1, 0.2, 0.0, -0.1] {
                let (a, pa) = old.step(voltage, 1.0 / 48_000.0).unwrap();
                let (b, pb) = explicit.step(voltage, 1.0 / 48_000.0).unwrap();
                assert_eq!(a.x, b.x);
                assert_eq!(
                    (a.supplied, a.dissipated, pa),
                    (b.supplied, b.dissipated, pb)
                );
            }
        }
    }

    #[test]
    fn sealed_media_require_coherent_states_and_a_supported_static_equilibrium() {
        let gas = declared_air();
        assert!(TsDriver::new_with_media(card(), Some(0.02), None, gas, None).is_err());
        assert!(TsDriver::new_with_media(card(), None, None, gas, Some(gas)).is_err());
        let biased = GasState::try_new(
            &GasSpec::dry_air_ussa1976(),
            gas.temperature,
            2.0 * gas.pressure,
        )
        .unwrap();
        assert!(TsDriver::new_with_media(card(), Some(0.02), None, gas, Some(biased)).is_err());
        let mut stale = gas;
        stale.sound_speed *= 1.1;
        for (exterior, cavity) in [(stale, gas), (gas, stale)] {
            assert!(matches!(
                TsDriver::new_with_media(card(), Some(0.02), None, exterior, Some(cavity)),
                Err(SpeakerError::Gas(_)),
            ));
        }
    }
}
