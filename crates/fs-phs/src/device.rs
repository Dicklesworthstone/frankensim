//! Device cards + nonlinear circuit islands (music bead
//! `frankensim-music-v8-root-3ez8g.9.3`): distortion as a CONSTITUTIVE
//! property of a device at an operating point — a triode's transfer
//! curvature, a diode pair's clipping knee — never a waveshaper
//! aesthetic. The split-circuit image factorizes the LINEAR network
//! (solved by the [`crate::circuit`] descriptor machinery) from small
//! nonlinear DEVICE ISLANDS that iterate against it (D14: Newton on
//! the islands only, with ANALYTIC Jacobians — these laws are smooth).
//!
//! CARDS: typed parameters + named provenance + a validity region, in
//! the matdb pattern (locally stored v1; the schema does not preclude
//! matdb migration). The triode is the Koren-class parametric model
//! with the widely published 12AX7 parameter set (authored
//! transcription of Koren's openly published values — Estimate); the
//! diode pair is Shockley-class with series resistance.
//!
//! SPLIT vs MONOLITH: the relaxation here (linear DAE step, island
//! Newton, a fixed number of Gauss–Seidel sweeps with disclosed
//! residuals) is the BUDGET-SHAPED image; the
//! full-DAE-with-devices monolith is the AUTHORITY image, and the
//! bake-off between them belongs to electric-gates (.9.5) — recorded,
//! not decided here.
//!
//! Admission: the DC BIAS solve runs first, and a mis-biased stage
//! (no convergence, or an operating point outside the card's validity
//! region) REFUSES with a named diagnostic — never silent garbage.

use crate::det;

/// Koren-class triode card with provenance and validity.
#[derive(Debug, Clone, PartialEq)]
pub struct TriodeCard {
    /// Amplification factor mu.
    pub mu: f64,
    /// Exponent `ex`.
    pub ex: f64,
    /// Plate-current scale `kg1`.
    pub kg1: f64,
    /// Knee sharpness `kp`.
    pub kp: f64,
    /// Grid-bias curvature `kvb`.
    pub kvb: f64,
    /// Validity: maximum plate voltage the source curves cover [V].
    pub plate_v_max: f64,
    /// Validity: grid range the source curves cover [V] (negative).
    pub grid_v_min: f64,
    /// Provenance string.
    pub source: &'static str,
}

impl TriodeCard {
    /// The widely published Koren 12AX7 parameter set (authored
    /// transcription; Estimate).
    #[must_use]
    pub fn koren_12ax7() -> TriodeCard {
        TriodeCard {
            mu: 100.0,
            ex: 1.4,
            kg1: 1060.0,
            kp: 600.0,
            kvb: 300.0,
            plate_v_max: 400.0,
            grid_v_min: -6.0,
            source: "Koren improved vacuum-tube models (published parameter set for the \
                     12AX7); authored transcription, Estimate",
        }
    }

    /// Plate current [A] for plate/grid voltages (cathode-referenced),
    /// with the analytic partial derivatives (dIp/dVp, dIp/dVg).
    #[must_use]
    pub fn plate_current(&self, v_p: f64, v_g: f64) -> (f64, f64, f64) {
        if v_p <= 0.0 {
            return (0.0, 0.0, 0.0);
        }
        // E1 = (Vp/kp) ln(1 + exp(kp (1/mu + Vg/sqrt(kvb + Vp^2)))).
        let s = (self.kvb + v_p * v_p).sqrt();
        let arg = self.kp * (1.0 / self.mu + v_g / s);
        // Numerically safe softplus.
        let (softplus, sigmoid) = if arg > 30.0 {
            (arg, 1.0)
        } else if arg < -30.0 {
            (det::exp(arg), det::exp(arg))
        } else {
            let e = det::exp(arg);
            ((1.0 + e).ln(), e / (1.0 + e))
        };
        let e1 = v_p / self.kp * softplus;
        if e1 <= 0.0 {
            return (0.0, 0.0, 0.0);
        }
        let ip = 2.0 * e1.powf(self.ex) / self.kg1;
        let dip_de1 = 2.0 * self.ex * e1.powf(self.ex - 1.0) / self.kg1;
        // dE1/dVp and dE1/dVg.
        let ds_dvp = v_p / s;
        let darg_dvp = -self.kp * v_g * ds_dvp / (s * s);
        let de1_dvp = softplus / self.kp + v_p / self.kp * sigmoid * darg_dvp;
        let de1_dvg = v_p / self.kp * sigmoid * self.kp / s;
        (ip, dip_de1 * de1_dvp, dip_de1 * de1_dvg)
    }
}

/// Shockley diode-pair card (antiparallel) with series resistance.
#[derive(Debug, Clone, PartialEq)]
pub struct DiodePairCard {
    /// Saturation current [A].
    pub is_a: f64,
    /// Ideality factor.
    pub n_ideal: f64,
    /// Thermal voltage [V].
    pub vt_v: f64,
    /// Series resistance [ohm].
    pub rs_ohm: f64,
    /// Validity: maximum junction current the law covers [A].
    pub i_max_a: f64,
    /// Provenance.
    pub source: &'static str,
}

impl DiodePairCard {
    /// Silicon-class small-signal pair (authored Estimate).
    #[must_use]
    pub fn silicon_class() -> DiodePairCard {
        DiodePairCard {
            is_a: 2.5e-9,
            n_ideal: 1.75,
            vt_v: 0.02585,
            rs_ohm: 1.0,
            i_max_a: 0.05,
            source: "silicon small-signal class values (1N4148-like magnitudes); authored, \
                     Estimate",
        }
    }

    /// Pair current and its derivative at junction voltage `v`.
    #[must_use]
    pub fn current(&self, v: f64) -> (f64, f64) {
        let x = v / (self.n_ideal * self.vt_v);
        let x = x.clamp(-60.0, 60.0);
        let i = 2.0 * self.is_a * x.sinh();
        let di = 2.0 * self.is_a * x.cosh() / (self.n_ideal * self.vt_v);
        (i, di)
    }
}

/// Typed refusals.
#[derive(Debug)]
pub enum DeviceError {
    /// Bad parameter, by name.
    Invalid {
        /// What.
        what: &'static str,
    },
    /// The DC bias admission failed to converge.
    BiasDidNotConverge {
        /// Residual at stall.
        residual: f64,
    },
    /// The operating point (or a runtime excursion) left the card's
    /// validity region.
    OutsideValidity {
        /// Which bound.
        what: &'static str,
        /// The offending value.
        value: f64,
    },
    /// Circuit-side refusal.
    Circuit(crate::circuit::CircuitError),
}

impl From<crate::circuit::CircuitError> for DeviceError {
    fn from(e: crate::circuit::CircuitError) -> Self {
        DeviceError::Circuit(e)
    }
}

/// Per-block solver telemetry (the disclosure the discipline demands).
#[derive(Debug, Clone, Copy, Default)]
pub struct IslandTelemetry {
    /// Island Newton iterations this block.
    pub iterations: usize,
    /// Worst island residual [A].
    pub worst_residual_a: f64,
    /// Operating-point drift from the bias point [V].
    pub operating_drift_v: f64,
}

/// A single-triode gain stage: supply + plate resistor from the
/// linear DAE, the triode as an island, fixed cathode bias.
pub struct TriodeStage {
    card: TriodeCard,
    dae: crate::circuit::CircuitDae,
    x: Vec<f64>,
    supply_v: f64,
    bias_v: f64,
    plate_node_coord: usize,
    rl_ohm: f64,
    ip_last: f64,
    bias_plate_v: f64,
    /// Fixed island sweeps per sample (deterministic budget).
    sweeps: usize,
    telemetry: IslandTelemetry,
}

impl TriodeStage {
    /// Build and BIAS-ADMIT the stage: supply `supply_v` through
    /// `rl_ohm` into the plate; grid at `bias_v` (negative). The DC
    /// solve is Newton on `(supply − Vp)/RL = Ip(Vp, bias)` with the
    /// card's analytic derivative; failure or an out-of-validity
    /// operating point refuses.
    ///
    /// # Errors
    /// [`DeviceError`] as documented.
    pub fn new(
        card: TriodeCard,
        supply_v: f64,
        rl_ohm: f64,
        bias_v: f64,
    ) -> Result<Self, DeviceError> {
        if !(supply_v.is_finite() && supply_v > 0.0) {
            return Err(DeviceError::Invalid {
                what: "supply must be positive",
            });
        }
        if !(rl_ohm.is_finite() && rl_ohm > 0.0) {
            return Err(DeviceError::Invalid {
                what: "plate resistor must be positive",
            });
        }
        if bias_v > 0.0 || bias_v < card.grid_v_min {
            return Err(DeviceError::OutsideValidity {
                what: "grid bias outside the card's validity",
                value: bias_v,
            });
        }
        // DC bias Newton.
        let mut v_p = supply_v * 0.6;
        let mut converged = false;
        let mut residual = f64::INFINITY;
        for _ in 0..80 {
            let (ip, dip_dvp, _) = card.plate_current(v_p, bias_v);
            let f = (supply_v - v_p) / rl_ohm - ip;
            residual = f.abs();
            if residual < 1.0e-12 * (supply_v / rl_ohm) {
                converged = true;
                break;
            }
            let df = -1.0 / rl_ohm - dip_dvp;
            v_p -= f / df;
            v_p = v_p.clamp(1.0, supply_v);
        }
        if !converged {
            return Err(DeviceError::BiasDidNotConverge { residual });
        }
        if v_p > card.plate_v_max {
            return Err(DeviceError::OutsideValidity {
                what: "bias plate voltage above the card's curves",
                value: v_p,
            });
        }
        // Linear side as a real DAE: supply source 1->0 (port 0),
        // RL 1->2, island current drawn at node 2 (port 1).
        let graph = crate::circuit::CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, crate::circuit::Branch::VoltageSource { port: 0 }),
                (1, 2, crate::circuit::Branch::Resistor { ohms: rl_ohm }),
                (2, 0, crate::circuit::Branch::CurrentSource { port: 1 }),
            ],
            transformers: vec![],
        };
        let dae = crate::circuit::assemble_circuit(&graph)?;
        let (ip0, _, _) = card.plate_current(v_p, bias_v);
        let x0 = vec![0.0; dae.system.state_dim()];
        // THE DESCRIPTOR READ LAW (found by probe): algebraic efforts
        // are MIDPOINTS of the stored coordinate across a step, so the
        // stored x oscillates around the true potential (the dt=0
        // solve from zero lands at exactly 2x). Node potentials must
        // be read as midpoints; u = -ip draws the plate current.
        let x = dae.consistent_initial_state(&x0, &[supply_v, -ip0])?;
        let plate_node_coord = dae.node_potential_index[1];
        Ok(Self {
            card,
            dae,
            x,
            supply_v,
            bias_v,
            plate_node_coord,
            rl_ohm,
            ip_last: ip0,
            bias_plate_v: v_p,
            sweeps: 3,
            telemetry: IslandTelemetry::default(),
        })
    }

    /// Bias-point plate voltage [V].
    #[must_use]
    pub fn bias_plate_v(&self) -> f64 {
        self.bias_plate_v
    }

    /// Last telemetry.
    #[must_use]
    pub fn telemetry(&self) -> IslandTelemetry {
        self.telemetry
    }

    /// One audio sample: grid input `v_in` rides the bias; the island
    /// and the linear DAE relax for a FIXED number of sweeps
    /// (deterministic budget; residual disclosed in telemetry).
    ///
    /// # Errors
    /// Validity excursions and circuit refusals.
    pub fn step(&mut self, v_in: f64, dt: f64) -> Result<f64, DeviceError> {
        let v_g = self.bias_v + v_in;
        if v_g < self.card.grid_v_min {
            return Err(DeviceError::OutsideValidity {
                what: "grid drive below the card's validity",
                value: v_g,
            });
        }
        let mut iterations = 0usize;
        let mut worst = 0.0f64;
        let mut x_next = self.x.clone();
        for _ in 0..self.sweeps {
            let rec = crate::step_descriptor(
                &self.dae.system,
                &self.x,
                &[self.supply_v, -self.ip_last],
                dt,
            )
            .map_err(crate::circuit::CircuitError::from)?;
            // Midpoint read (see the constructor's descriptor-read law).
            let v_p = f64::midpoint(self.x[self.plate_node_coord], rec.x[self.plate_node_coord]);
            x_next = rec.x;
            if v_p > self.card.plate_v_max {
                return Err(DeviceError::OutsideValidity {
                    what: "plate voltage above the card's curves",
                    value: v_p,
                });
            }
            // Island Newton with the LOAD-LINE-AWARE denominator:
            // the plain fixed point has gain -RL*g_p (about -1.7 at
            // this operating point) and DIVERGES — measured, the plate
            // relaxed to 790 V. dVp/dIp = -RL from the linear side, so
            // the Newton step on ip is res / (1 + RL * dIp/dVp).
            let (ip, dip_dvp, _) = self.card.plate_current(v_p.max(0.0), v_g);
            let res = ip - self.ip_last;
            worst = worst.max(res.abs());
            let denom = 1.0 + self.rl_ohm * dip_dvp;
            self.ip_last += res / denom.max(1.0e-6);
            iterations += 1;
        }
        let v_p = f64::midpoint(self.x[self.plate_node_coord], x_next[self.plate_node_coord]);
        self.x = x_next;
        self.telemetry = IslandTelemetry {
            iterations,
            worst_residual_a: worst,
            operating_drift_v: v_p - self.bias_plate_v,
        };
        Ok(v_p)
    }
}

/// A diode-pair clipper: series R from the input, the pair as an
/// island to ground (v1 memoryless output node; a load C is the
/// circuit consumer's addition).
pub struct DiodeClipper {
    card: DiodePairCard,
    r_ohm: f64,
}

impl DiodeClipper {
    /// Build the clipper.
    ///
    /// # Errors
    /// Bad parameters by name.
    pub fn new(card: DiodePairCard, r_ohm: f64) -> Result<Self, DeviceError> {
        if !(r_ohm.is_finite() && r_ohm > 0.0) {
            return Err(DeviceError::Invalid {
                what: "series resistance must be positive",
            });
        }
        for (v, what) in [
            (card.is_a, "Is"),
            (card.n_ideal, "ideality"),
            (card.vt_v, "Vt"),
            (card.rs_ohm, "Rs"),
        ] {
            if !(v.is_finite() && v > 0.0) {
                return Err(DeviceError::Invalid { what });
            }
        }
        Ok(Self { card, r_ohm })
    }

    /// Output voltage for input `v_in`: Newton with the analytic
    /// Jacobian on `(v_in − v)/R = i_pair(v − i R s)` (fixed 40-cap,
    /// state-relative tolerance; returns the solution, iterations,
    /// and the residual — deterministic and disclosed).
    ///
    /// # Errors
    /// [`DeviceError::OutsideValidity`] past the card's current cap.
    pub fn solve(&self, v_in: f64) -> Result<(f64, usize, f64), DeviceError> {
        let mut v = v_in.clamp(-1.0, 1.0);
        let mut residual = f64::INFINITY;
        let mut iters = 0usize;
        for k in 0..40 {
            iters = k + 1;
            let (i, di) = self.card.current(v);
            let f = (v_in - v) / self.r_ohm - i;
            residual = f.abs();
            if residual < 1.0e-14 * (v_in.abs() / self.r_ohm).max(1.0e-12) {
                break;
            }
            let df = -1.0 / self.r_ohm - di;
            let step = f / df;
            v -= step.clamp(-0.3, 0.3);
        }
        let (i, _) = self.card.current(v);
        if i.abs() > self.card.i_max_a {
            return Err(DeviceError::OutsideValidity {
                what: "junction current above the card's cap",
                value: i,
            });
        }
        Ok((v, iters, residual))
    }
}

#[cfg(test)]
mod device_tests {
    use super::*;

    /// Pins the descriptor READ LAW this module depends on: algebraic
    /// efforts are midpoints of the stored coordinate, so the dt=0
    /// consistency solve from a zero seed stores 2x the true potential
    /// and every potential read must be a midpoint across the step.
    #[test]
    fn dv_000_midpoint_read_law() {
        let card = TriodeCard::koren_12ax7();
        let (supply, rl, bias) = (300.0, 100.0e3, -1.5);
        // Re-run the bias Newton by hand.
        let mut v_p = supply * 0.6;
        for _ in 0..80 {
            let (ip, dip_dvp, _) = card.plate_current(v_p, bias);
            let f = (supply - v_p) / rl - ip;
            if f.abs() < 1.0e-15 {
                break;
            }
            v_p -= f / (-1.0 / rl - dip_dvp);
            v_p = v_p.clamp(1.0, supply);
        }
        let (ip0, _, _) = card.plate_current(v_p, bias);
        println!("bias vp {v_p:.3}, ip0 {ip0:.6e}");
        let graph = crate::circuit::CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, crate::circuit::Branch::VoltageSource { port: 0 }),
                (1, 2, crate::circuit::Branch::Resistor { ohms: rl }),
                (2, 0, crate::circuit::Branch::CurrentSource { port: 1 }),
            ],
            transformers: vec![],
        };
        let dae = crate::circuit::assemble_circuit(&graph).expect("assemble");
        let x = dae
            .consistent_initial_state(&vec![0.0; dae.system.state_dim()], &[supply, -ip0])
            .expect("ics");
        // Midpoints against the zero seed ARE the physical potentials:
        // v1 = the supply exactly, v2 = the bias plate voltage from the
        // independent Newton. Raw reads are 2x and MUST NOT be used.
        let v1_mid = 0.5 * x[dae.node_potential_index[0]];
        let v2_mid = 0.5 * x[dae.node_potential_index[1]];
        assert!((v1_mid - supply).abs() < 1.0e-9, "v1 mid {v1_mid}");
        assert!(
            (v2_mid - v_p).abs() < 1.0e-6 * v_p.abs(),
            "v2 mid {v2_mid} vs bias {v_p}"
        );
        assert!(
            (x[dae.node_potential_index[0]] - supply).abs() > 100.0,
            "raw read unexpectedly physical; the midpoint law changed"
        );
    }

    const RATE: f64 = 96_000.0;

    #[test]
    fn dv_001_stage_gain_matches_the_small_signal_derivation() {
        // The stage's measured small-signal gain vs the INDEPENDENT
        // derivation gain = mu RL / (RL + rp), with rp = 1/(dIp/dVp)
        // read from the card law by finite differences IN THE TEST
        // (a different route than the stage's own analytic Jacobian).
        let card = TriodeCard::koren_12ax7();
        let (supply, rl, bias) = (300.0, 100.0e3, -1.5);
        let mut stage = TriodeStage::new(card.clone(), supply, rl, bias).expect("bias admits");
        let vp0 = stage.bias_plate_v();
        // FD transconductance and plate resistance at the bias point.
        let h = 1.0e-4;
        let (ip_p, _, _) = card.plate_current(vp0 + h, bias);
        let (ip_m, _, _) = card.plate_current(vp0 - h, bias);
        let g_p = (ip_p - ip_m) / (2.0 * h);
        let (ig_p, _, _) = card.plate_current(vp0, bias + h);
        let (ig_m, _, _) = card.plate_current(vp0, bias - h);
        let gm = (ig_p - ig_m) / (2.0 * h);
        let rp = 1.0 / g_p;
        let gain_expected = gm * rp * rl / (rl + rp);
        // Measured: tiny sinusoid, steady-state output amplitude.
        let dt = 1.0 / RATE;
        let f = 200.0;
        let n = (12.0 / f / dt) as usize;
        let mut peak = 0.0f64;
        for k in 0..n {
            let vin = 1.0e-3 * crate::det::sin(core::f64::consts::TAU * f * k as f64 * dt);
            let vp = stage.step(vin, dt).expect("step");
            if k > 2 * n / 3 {
                peak = peak.max((vp - vp0).abs());
            }
        }
        let gain = peak / 1.0e-3;
        let rel = (gain - gain_expected).abs() / gain_expected;
        assert!(
            rel < 0.05,
            "stage gain {gain:.1} vs small-signal {gain_expected:.1} (rel {rel:.3})"
        );
        println!(
            "{{\"suite\":\"fs-phs\",\"case\":\"dv-001-gain\",\"gain\":{gain:.1},\
             \"analytic\":{gain_expected:.1},\"vp0\":{vp0:.1}}}"
        );
    }

    #[test]
    fn dv_002_even_harmonics_rise_before_odd() {
        // THE EMERGENT DISTORTION GATE: the single-ended triode's
        // transfer curvature puts H2 above H3 at moderate drive, both
        // rising with level — from the card + bias alone, never
        // shaped. Ladder logged.
        let card = TriodeCard::koren_12ax7();
        let mut rows = Vec::new();
        for &amp in &[0.05f64, 0.2, 0.6] {
            let mut stage = TriodeStage::new(card.clone(), 300.0, 100.0e3, -1.5).expect("bias");
            let dt = 1.0 / RATE;
            let f = 200.0;
            let n = (30.0 / f / dt) as usize;
            let mut out = Vec::with_capacity(n);
            for k in 0..n {
                let vin = amp * crate::det::sin(core::f64::consts::TAU * f * k as f64 * dt);
                out.push(stage.step(vin, dt).expect("step"));
            }
            let tail = &out[n / 2..];
            let mean = tail.iter().sum::<f64>() / tail.len() as f64;
            let line = |h: usize| -> f64 {
                let w = core::f64::consts::TAU * f * h as f64 / RATE;
                let (mut re, mut im) = (0.0, 0.0);
                for (k, &v) in tail.iter().enumerate() {
                    re += (v - mean) * crate::det::cos(w * k as f64);
                    im -= (v - mean) * crate::det::sin(w * k as f64);
                }
                (re * re + im * im).sqrt()
            };
            let (h1, h2, h3) = (line(1), line(2), line(3));
            rows.push((amp, h2 / h1, h3 / h1));
            println!(
                "{{\"suite\":\"fs-phs\",\"case\":\"dv-002-ladder\",\"amp\":{amp},\
                 \"h2_rel\":{:.5},\"h3_rel\":{:.5}}}",
                h2 / h1,
                h3 / h1
            );
        }
        for (amp, h2, h3) in &rows {
            assert!(
                h2 > h3,
                "at drive {amp}: H2 ({h2:.5}) must lead H3 ({h3:.5}) — the single-ended \
                 signature"
            );
        }
        assert!(rows[2].1 > rows[0].1, "H2 must rise with drive");
        assert!(rows[2].2 > rows[0].2, "H3 must rise with drive");
    }

    #[test]
    fn dv_003_diode_clipper_matches_a_method_diverse_solve() {
        // The island's Newton vs a BISECTION solve of the same scalar
        // equation in-test (method diversity), plus the knee sanity:
        // far above the knee the output compresses logarithmically.
        let clipper = DiodeClipper::new(DiodePairCard::silicon_class(), 4.7e3).expect("clipper");
        let card = DiodePairCard::silicon_class();
        for &vin in &[0.05f64, 0.3, 1.0, 3.0] {
            let (v_newton, iters, res) = clipper.solve(vin).expect("newton");
            // Bisection oracle.
            let f = |v: f64| (vin - v) / 4.7e3 - card.current(v).0;
            let (mut lo, mut hi) = (0.0f64, vin.min(1.0));
            for _ in 0..200 {
                let mid = f64::midpoint(lo, hi);
                if f(mid) > 0.0 {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            let v_bisect = f64::midpoint(lo, hi);
            assert!(
                (v_newton - v_bisect).abs() < 1.0e-9,
                "vin {vin}: Newton {v_newton:.9} vs bisection {v_bisect:.9}"
            );
            println!(
                "{{\"suite\":\"fs-phs\",\"case\":\"dv-003-clipper\",\"vin\":{vin},\
                 \"vout\":{v_newton:.5},\"iters\":{iters},\"residual\":{res:.2e}}}"
            );
        }
        // Compression: 3 V in produces far less than 3x the 1 V output.
        let v1 = clipper.solve(1.0).expect("1V").0;
        let v3 = clipper.solve(3.0).expect("3V").0;
        assert!(v3 < 1.5 * v1, "the knee must compress ({v3:.3} vs {v1:.3})");
    }

    #[test]
    fn dv_004_refusals_fire_by_name() {
        let card = TriodeCard::koren_12ax7();
        // Grid bias outside validity.
        assert!(matches!(
            TriodeStage::new(card.clone(), 300.0, 100.0e3, -20.0),
            Err(DeviceError::OutsideValidity { .. })
        ));
        assert!(matches!(
            TriodeStage::new(card.clone(), 300.0, 100.0e3, 1.0),
            Err(DeviceError::OutsideValidity { .. })
        ));
        // Supply above the card's plate curves at bias -> validity.
        assert!(matches!(
            TriodeStage::new(card.clone(), 900.0, 1.0e3, -6.0),
            Err(DeviceError::OutsideValidity { .. })
        ));
        // Runtime grid excursion refuses.
        let mut stage = TriodeStage::new(card, 300.0, 100.0e3, -1.5).expect("bias");
        assert!(matches!(
            stage.step(-10.0, 1.0 / RATE),
            Err(DeviceError::OutsideValidity { .. })
        ));
        // Clipper current cap.
        let clipper = DiodeClipper::new(
            DiodePairCard {
                i_max_a: 1.0e-6,
                ..DiodePairCard::silicon_class()
            },
            10.0,
        )
        .expect("clipper");
        assert!(matches!(
            clipper.solve(5.0),
            Err(DeviceError::OutsideValidity { .. })
        ));
        println!("{{\"suite\":\"fs-phs\",\"case\":\"dv-004-refusals\",\"verdict\":\"pass\"}}");
    }

    #[test]
    fn dv_005_islands_are_deterministic_with_fixed_budgets() {
        // Bitwise repeatability under the fixed sweep budget, and the
        // telemetry discloses iterations/residual/drift.
        let run = || -> (Vec<u64>, IslandTelemetry) {
            let mut stage =
                TriodeStage::new(TriodeCard::koren_12ax7(), 300.0, 100.0e3, -1.5).expect("bias");
            let dt = 1.0 / RATE;
            let mut bits = Vec::new();
            for k in 0..2000 {
                let vin = 0.3 * crate::det::sin(core::f64::consts::TAU * 200.0 * f64::from(k) * dt);
                bits.push(stage.step(vin, dt).expect("step").to_bits());
            }
            (bits, stage.telemetry())
        };
        let (a, ta) = run();
        let (b, _) = run();
        assert_eq!(a, b, "island runs must be bitwise identical");
        assert!(ta.iterations > 0 && ta.worst_residual_a.is_finite());
        println!(
            "{{\"suite\":\"fs-phs\",\"case\":\"dv-005-determinism\",\"iters\":{},\
             \"residual\":{:.2e},\"drift_v\":{:.3}}}",
            ta.iterations, ta.worst_residual_a, ta.operating_drift_v
        );
    }
}
