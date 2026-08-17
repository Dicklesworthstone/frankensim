//! Faraday pickup port (music bead `frankensim-music-v8-root-3ez8g.9.2`):
//! where string motion becomes voltage. `EMF = −dΦ/dt` with the flux a
//! WINDOW-WEIGHTED functional of string state — the pickup senses a
//! weighted segment of the string, not a point, and that aperture
//! window is what makes bridge-vs-neck voicing physical.
//!
//! D20 AS STRUCTURE: the pickup NEVER re-discretizes the string. It
//! owns pose + per-mode gains only — `emf()` is a pure `&self` read of
//! the string owner's modal velocities. There is no step, no state, no
//! copy of string coordinates anywhere in this module.
//!
//! Provenance honesty (v1): the height profile of the effective flux
//! gain and the window shape are AUTHORED Estimate models; the
//! MAGNETOSTATIC BAKE (an offline lab minting the position-dependent B
//! map from magnet geometry, the lips/jet authority-lab pattern) is
//! the RECORDED FOLLOW-UP that upgrades provenance without changing
//! the runtime shape. The flux map is LINEARIZED — Estimate for large
//! excursions, disclosed; saturation/asymmetry enter only with the
//! bake or a measured card.

use fs_math::det;

/// Pickup pose along a string of speaking length `L` — geometry and
/// setup parameters, never magic numbers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PickupPose {
    /// Station as a fraction of speaking length, strictly in (0, 1).
    pub station_frac: f64,
    /// Pole-to-string height [m].
    pub height_m: f64,
    /// Aperture (sensing window width) as a fraction of speaking
    /// length, in (0, 0.5].
    pub aperture_frac: f64,
}

/// Typed refusals.
#[derive(Debug, PartialEq, Eq)]
pub enum PickupError {
    /// Bad pose, by name.
    Invalid {
        /// What.
        what: &'static str,
    },
}

/// A lumped Faraday pickup bound to a modal string basis.
#[derive(Debug, Clone, PartialEq)]
pub struct Pickup {
    pose: PickupPose,
    /// Per-mode window-weighted gains `g_k = ∫ w(x) φ_k(x) dx`.
    gains: Vec<f64>,
    /// Effective flux gain κ(h) [V s / m] (authored Estimate profile).
    flux_gain: f64,
}

/// Authored height profile for the effective flux gain: a dipole-like
/// falloff `κ = κ0 / (1 + h/h0)^2` with `κ0 = 40 mV·s/m` at contact
/// scale and `h0 = 3 mm` (Estimate; the magnetostatic bake replaces
/// this).
fn flux_gain_of_height(height_m: f64) -> f64 {
    let kappa0 = 4.0e-2;
    let h0 = 3.0e-3;
    kappa0 / (1.0 + height_m / h0).powi(2)
}

impl Pickup {
    /// Bind a pickup to the first `n_modes` sinusoidal modes of a
    /// string (mode shapes `sin(k π x / L)` on the unit-length
    /// coordinate — the pickup reads whatever image OWNS those modal
    /// coordinates; it never builds its own).
    ///
    /// The window is a raised cosine over
    /// `[station − aperture/2, station + aperture/2]`, normalized to
    /// unit integral, so the point-pickup limit is exact as the
    /// aperture shrinks.
    ///
    /// # Errors
    /// [`PickupError::Invalid`] on a station outside (0,1), an
    /// aperture outside (0, 0.5], a non-positive height, or zero
    /// modes.
    pub fn bind(pose: PickupPose, n_modes: usize) -> Result<Self, PickupError> {
        if !(pose.station_frac.is_finite() && pose.station_frac > 0.0 && pose.station_frac < 1.0)
        {
            return Err(PickupError::Invalid {
                what: "station must lie strictly inside the speaking length",
            });
        }
        if !(pose.aperture_frac.is_finite()
            && pose.aperture_frac > 0.0
            && pose.aperture_frac <= 0.5)
        {
            return Err(PickupError::Invalid {
                what: "aperture must lie in (0, 0.5] of the speaking length",
            });
        }
        if !(pose.height_m.is_finite() && pose.height_m > 0.0) {
            return Err(PickupError::Invalid {
                what: "height must be positive",
            });
        }
        if n_modes == 0 {
            return Err(PickupError::Invalid {
                what: "a pickup needs at least one mode to read",
            });
        }
        // g_k = ∫ w(x) sin(k π x) dx over the window, 512-point
        // midpoint rule (the window is smooth; the gains are logged
        // and cross-checked analytically in the gates).
        let half = 0.5 * pose.aperture_frac;
        let lo = (pose.station_frac - half).max(0.0);
        let hi = (pose.station_frac + half).min(1.0);
        let steps = 512usize;
        let dx = (hi - lo) / steps as f64;
        let mut gains = vec![0.0f64; n_modes];
        let mut norm = 0.0f64;
        for s in 0..steps {
            let x = lo + (s as f64 + 0.5) * dx;
            // Raised cosine centered on the station.
            let w = 0.5
                + 0.5
                    * det::cos(
                        core::f64::consts::PI * (x - pose.station_frac) / half.max(1e-300),
                    );
            norm += w * dx;
            for (k, g) in gains.iter_mut().enumerate() {
                *g += w * det::sin((k + 1) as f64 * core::f64::consts::PI * x) * dx;
            }
        }
        for g in &mut gains {
            *g /= norm;
        }
        Ok(Self {
            pose,
            gains,
            flux_gain: flux_gain_of_height(pose.height_m),
        })
    }

    /// The pose.
    #[must_use]
    pub fn pose(&self) -> PickupPose {
        self.pose
    }

    /// The per-mode window gains (logged provenance).
    #[must_use]
    pub fn mode_gains(&self) -> &[f64] {
        &self.gains
    }

    /// The effective flux gain [V s / m].
    #[must_use]
    pub fn flux_gain_v_s_per_m(&self) -> f64 {
        self.flux_gain
    }

    /// The EMF for the string owner's modal velocities — a pure read;
    /// the linearized Faraday law `−κ Σ g_k v_k`.
    ///
    /// # Errors
    /// [`PickupError::Invalid`] when the velocity slice is shorter
    /// than the bound mode count.
    pub fn emf_v(&self, modal_velocities: &[f64]) -> Result<f64, PickupError> {
        if modal_velocities.len() < self.gains.len() {
            return Err(PickupError::Invalid {
                what: "velocity slice shorter than the bound mode count",
            });
        }
        let mut acc = 0.0;
        for (g, v) in self.gains.iter().zip(modal_velocities) {
            acc += g * v;
        }
        Ok(-self.flux_gain * acc)
    }

    /// JSON log line: pose, flux gain, and every mode gain.
    #[must_use]
    pub fn debug_line(&self) -> String {
        format!(
            "{{\"suite\":\"fs-couple\",\"case\":\"pickup-binding\",\"station\":{:.4},\
             \"height_m\":{:.4e},\"aperture\":{:.4},\"flux_gain\":{:.4e},\"gains\":{:?}}}",
            self.pose.station_frac,
            self.pose.height_m,
            self.pose.aperture_frac,
            self.flux_gain,
            self.gains
        )
    }
}

#[cfg(test)]
mod pickup_tests {
    use super::*;

    const RATE: f64 = 48_000.0;

    /// Exact-ZOH modal string rotors (the STRING OWNER in these
    /// fixtures): returns per-sample modal velocities for a pluck.
    struct ModalString {
        rot: Vec<(f64, f64)>,
        state: Vec<(f64, f64)>,
        omega: Vec<f64>,
    }

    impl ModalString {
        fn pluck(n_modes: usize, f0: f64, pluck_frac: f64) -> ModalString {
            let dt = 1.0 / RATE;
            let mut rot = Vec::new();
            let mut state = Vec::new();
            let mut omega = Vec::new();
            for k in 1..=n_modes {
                let w = core::f64::consts::TAU * f0 * k as f64;
                let zeta = 1.0e-4_f64;
                let wd = w * (1.0 - zeta * zeta).sqrt();
                let decay = det::exp(-zeta * w * dt);
                rot.push((decay * det::cos(wd * dt), decay * det::sin(wd * dt)));
                // Pluck shape: q_k(0) ~ sin(k pi x_p)/k^2.
                let q0 = det::sin(k as f64 * core::f64::consts::PI * pluck_frac)
                    / (k * k) as f64;
                state.push((q0, 0.0));
                omega.push(w);
            }
            ModalString { rot, state, omega }
        }

        /// Advance one sample; return modal velocities `q̇_k`.
        fn step(&mut self) -> Vec<f64> {
            for ((c, s), st) in self.rot.iter().zip(self.state.iter_mut()) {
                let (re, im) = *st;
                *st = (c * re - s * im, s * re + c * im);
            }
            self.state
                .iter()
                .zip(&self.omega)
                .map(|(&(_, im), &w)| -w * im)
                .collect()
        }
    }

    fn spectrum_lines(signal: &[f64], f0: f64, n_lines: usize) -> Vec<f64> {
        // Per-harmonic projection magnitudes (the string is modal, so
        // lines at k*f0 carry everything).
        (1..=n_lines)
            .map(|k| {
                let omega = core::f64::consts::TAU * f0 * k as f64 / RATE;
                let (mut re, mut im) = (0.0f64, 0.0f64);
                for (n, &v) in signal.iter().enumerate() {
                    re += v * det::cos(omega * n as f64);
                    im -= v * det::sin(omega * n as f64);
                }
                (re * re + im * im).sqrt()
            })
            .collect()
    }

    #[test]
    fn pk_001_station_voicing_matches_the_window_weighted_prediction() {
        // The comb physics: per-mode EMF line amplitudes must match
        // the ANALYTIC prediction |kappa * g_k * omega_k * q_k(0)|
        // computed in-test from the closed-form window integral (an
        // independent recomputation, not the module's own gains), and
        // the bridge pickup must be brighter than the neck pickup.
        let n_modes = 10usize;
        let f0 = 110.0;
        let mut centroids = Vec::new();
        for &(station, name) in &[(0.08f64, "bridge"), (0.35f64, "neck")] {
            let pose = PickupPose {
                station_frac: station,
                height_m: 3.0e-3,
                aperture_frac: 0.04,
            };
            let pickup = Pickup::bind(pose, n_modes).expect("bind");
            println!("{}", pickup.debug_line());
            let mut string = ModalString::pluck(n_modes, f0, 0.22);
            let n = (0.5 * RATE) as usize;
            let mut emf = Vec::with_capacity(n);
            for _ in 0..n {
                let v = string.step();
                emf.push(pickup.emf_v(&v).expect("emf"));
            }
            let lines = spectrum_lines(&emf, f0, n_modes);
            // Analytic prediction per mode (independent integral: the
            // same raised-cosine window integrated by Simpson on a
            // DIFFERENT grid density).
            let half = 0.5 * pose.aperture_frac;
            let analytic_gain = |k: usize| -> f64 {
                let steps = 1001usize;
                let lo = station - half;
                let hi = station + half;
                let h = (hi - lo) / (steps - 1) as f64;
                let f = |x: f64| -> (f64, f64) {
                    let w = 0.5
                        + 0.5
                            * det::cos(core::f64::consts::PI * (x - station) / half);
                    (w, w * det::sin(k as f64 * core::f64::consts::PI * x))
                };
                let mut wsum = 0.0;
                let mut gsum = 0.0;
                for i in 0..steps {
                    let x = lo + i as f64 * h;
                    let coef = if i == 0 || i == steps - 1 {
                        1.0
                    } else if i % 2 == 1 {
                        4.0
                    } else {
                        2.0
                    };
                    let (w, g) = f(x);
                    wsum += coef * w;
                    gsum += coef * g;
                }
                gsum / wsum
            };
            for k in 1..=n_modes {
                let q0 = det::sin(k as f64 * core::f64::consts::PI * 0.22) / (k * k) as f64;
                let expected =
                    (pickup.flux_gain_v_s_per_m() * analytic_gain(k) * core::f64::consts::TAU
                        * f0
                        * k as f64
                        * q0)
                        .abs();
                let measured = lines[k - 1] / (0.5 * RATE * 0.5); // projection norm
                if expected > 1.0e-8 {
                    let rel = (measured - expected).abs() / expected;
                    assert!(
                        rel < 0.05,
                        "{name} mode {k}: line {measured:.4e} vs analytic {expected:.4e} \
                         (rel {rel:.3})"
                    );
                }
            }
            let total: f64 = lines.iter().sum();
            let centroid: f64 = lines
                .iter()
                .enumerate()
                .map(|(i, &m)| (i + 1) as f64 * f0 * m)
                .sum::<f64>()
                / total;
            centroids.push(centroid);
            println!(
                "{{\"suite\":\"fs-couple\",\"case\":\"pk-001-{name}\",\"centroid_hz\":{centroid:.1}}}"
            );
        }
        assert!(
            centroids[0] > 1.2 * centroids[1],
            "bridge must be brighter than neck ({:.1} vs {:.1} Hz)",
            centroids[0],
            centroids[1]
        );
    }

    #[test]
    fn pk_002_two_pickups_comb_when_mixed() {
        // Mixing two pickups sums the per-mode gains BEFORE the
        // magnitude: modes where sin(k pi x_a) and sin(k pi x_b) have
        // opposite signs partially cancel. Assert the mixed line
        // spectrum matches |g_a + g_b| per mode, and exhibit one mode
        // where the mix is SMALLER than either alone (the comb notch).
        let n_modes = 10usize;
        let f0 = 110.0;
        let a = Pickup::bind(
            PickupPose {
                station_frac: 0.15,
                height_m: 3.0e-3,
                aperture_frac: 0.03,
            },
            n_modes,
        )
        .expect("a");
        let b = Pickup::bind(
            PickupPose {
                station_frac: 0.45,
                height_m: 3.0e-3,
                aperture_frac: 0.03,
            },
            n_modes,
        )
        .expect("b");
        let mut string = ModalString::pluck(n_modes, f0, 0.22);
        let n = (0.5 * RATE) as usize;
        let (mut ea, mut eb, mut mix) = (Vec::new(), Vec::new(), Vec::new());
        for _ in 0..n {
            let v = string.step();
            let va = a.emf_v(&v).expect("a");
            let vb = b.emf_v(&v).expect("b");
            ea.push(va);
            eb.push(vb);
            mix.push(va + vb);
        }
        let la = spectrum_lines(&ea, f0, n_modes);
        let lb = spectrum_lines(&eb, f0, n_modes);
        let lm = spectrum_lines(&mix, f0, n_modes);
        let mut notch_found = false;
        for k in 0..n_modes {
            let ga = a.mode_gains()[k];
            let gb = b.mode_gains()[k];
            let denom = (ga.abs() + gb.abs()).max(1e-12);
            let expected_ratio = (ga + gb).abs() / denom;
            let measured_ratio = lm[k] / (la[k] + lb[k]).max(1e-30);
            assert!(
                (measured_ratio - expected_ratio).abs() < 0.05,
                "mode {}: mix ratio {measured_ratio:.3} vs gain-sum prediction \
                 {expected_ratio:.3}",
                k + 1
            );
            if lm[k] < 0.8 * la[k].min(lb[k]) && la[k] > 1e-9 {
                notch_found = true;
            }
        }
        assert!(notch_found, "the mix must carve at least one comb notch");
        println!("{{\"suite\":\"fs-couple\",\"case\":\"pk-002-comb\",\"verdict\":\"pass\"}}");
    }

    #[test]
    fn pk_003_point_limit_and_no_state_duplication() {
        // Aperture -> 0 limit: the gains approach point sampling
        // sin(k pi x0) exactly (the window normalization is what makes
        // this true — a real limit check, not a smoke test). And D20:
        // the pickup is a pure reader — same state slice, same EMF,
        // twice, from two pickups, with the string owner untouched.
        let n_modes = 6usize;
        let x0 = 0.27;
        let tight = Pickup::bind(
            PickupPose {
                station_frac: x0,
                height_m: 3.0e-3,
                aperture_frac: 1.0e-3,
            },
            n_modes,
        )
        .expect("tight");
        for k in 1..=n_modes {
            let point = det::sin(k as f64 * core::f64::consts::PI * x0);
            let g = tight.mode_gains()[k - 1];
            assert!(
                (g - point).abs() < 5.0e-4,
                "mode {k}: tight-aperture gain {g:.6} vs point {point:.6}"
            );
        }
        // Pure-reader law: identical reads, no interior mutation.
        let v = vec![0.1, -0.2, 0.05, 0.0, 0.3, -0.1];
        let e1 = tight.emf_v(&v).expect("read 1");
        let e2 = tight.emf_v(&v).expect("read 2");
        assert!((e1 - e2).abs() == 0.0, "a reader cannot drift");
        // And the EMF recomputes from the SAME slice with the public
        // gains — no hidden copy of string state can exist.
        let manual: f64 = -tight.flux_gain_v_s_per_m()
            * tight
                .mode_gains()
                .iter()
                .zip(&v)
                .map(|(g, vv)| g * vv)
                .sum::<f64>();
        assert!((e1 - manual).abs() < 1e-18);
        println!("{{\"suite\":\"fs-couple\",\"case\":\"pk-003-point-limit\",\"verdict\":\"pass\"}}");
    }

    #[test]
    fn pk_004_refusals_fire_by_name() {
        let ok = PickupPose {
            station_frac: 0.2,
            height_m: 3.0e-3,
            aperture_frac: 0.05,
        };
        for (pose, what) in [
            (
                PickupPose {
                    station_frac: 0.0,
                    ..ok
                },
                "station must lie strictly inside the speaking length",
            ),
            (
                PickupPose {
                    station_frac: 1.2,
                    ..ok
                },
                "station must lie strictly inside the speaking length",
            ),
            (
                PickupPose {
                    aperture_frac: 0.0,
                    ..ok
                },
                "aperture must lie in (0, 0.5] of the speaking length",
            ),
            (
                PickupPose {
                    height_m: -1.0,
                    ..ok
                },
                "height must be positive",
            ),
        ] {
            assert_eq!(
                Pickup::bind(pose, 4).unwrap_err(),
                PickupError::Invalid { what }
            );
        }
        assert!(Pickup::bind(ok, 0).is_err());
        let p = Pickup::bind(ok, 4).expect("bind");
        assert!(p.emf_v(&[0.0, 0.0]).is_err(), "short slice refuses");
        println!("{{\"suite\":\"fs-couple\",\"case\":\"pk-004-refusals\",\"verdict\":\"pass\"}}");
    }

    #[test]
    fn pk_005_emf_drives_the_dae_circuit_end_to_end() {
        // The soft integration with .9.1's circuit DAE, live in the
        // same program increment: the pickup EMF drives a series R
        // into a load C (volume pot + cable — the classic treble
        // rolloff). The load voltage must track the RC-filtered EMF
        // (analytic first-order response at two test frequencies) and
        // the circuit's supply audit must hold throughout.
        use fs_phs::circuit::{Branch, CircuitGraph, assemble_circuit};
        let (r_ohm, c_f) = (250.0e3, 4.7e-10);
        let f_cut = 1.0 / (core::f64::consts::TAU * r_ohm * c_f);
        let graph = CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 2, Branch::Resistor { ohms: r_ohm }),
                (2, 0, Branch::Capacitor { farads: c_f }),
            ],
            transformers: vec![],
        };
        let dae = assemble_circuit(&graph).expect("assemble");
        let dt = 1.0 / RATE;
        for &f_ratio in &[0.25f64, 4.0] {
            let f = f_ratio * f_cut;
            let n_modes = 1usize;
            let pickup = Pickup::bind(
                PickupPose {
                    station_frac: 0.2,
                    height_m: 3.0e-3,
                    aperture_frac: 0.05,
                },
                n_modes,
            )
            .expect("bind");
            let mut string = ModalString::pluck(n_modes, f, 0.22);
            // Undamped-ish single mode: EMF is a near-pure sinusoid.
            let mut x = dae
                .consistent_initial_state(&vec![0.0; dae.system.state_dim()], &[0.0])
                .expect("ics");
            let n = (10.0 / f / dt) as usize;
            let mut emf_amp = 0.0f64;
            let mut out_amp = 0.0f64;
            let mut worst_defect = 0.0f64;
            for k in 0..n {
                let v = string.step();
                let emf = pickup.emf_v(&v).expect("emf");
                let (rec, defect) = dae.step_audited(&x, &[emf], dt).expect("step");
                x = rec.x;
                worst_defect = worst_defect.max(defect);
                if k > 2 * n / 3 {
                    emf_amp = emf_amp.max(emf.abs());
                    out_amp = out_amp.max(x[dae.charge_index[0]].abs() / c_f);
                }
            }
            let expected = emf_amp / (1.0 + (f / f_cut).powi(2)).sqrt();
            let rel = (out_amp - expected).abs() / expected;
            assert!(
                rel < 0.05,
                "RC load at {f_ratio}x cutoff: out {out_amp:.4e} vs analytic {expected:.4e} \
                 (rel {rel:.3})"
            );
            assert!(worst_defect < 1.0e-12, "supply audit {worst_defect:.3e}");
            println!(
                "{{\"suite\":\"fs-couple\",\"case\":\"pk-005-dae\",\"f_ratio\":{f_ratio},\
                 \"out\":{out_amp:.4e},\"analytic\":{expected:.4e},\"defect\":{worst_defect:.2e}}}"
            );
        }
    }
}
