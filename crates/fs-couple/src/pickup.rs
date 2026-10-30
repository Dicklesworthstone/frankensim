//! Faraday pickup port (music bead `frankensim-music-v8-root-3ez8g.9.2`):
//! where string motion becomes voltage. `EMF = −dΦ/dt` with the flux a
//! WINDOW-WEIGHTED functional of string state — the pickup senses a
//! weighted segment of the string, not a point, and that aperture
//! window is what makes bridge-vs-neck voicing physical.
//!
//! D20 AS STRUCTURE: the pickup NEVER re-discretizes the string. It
//! owns pose, per-mode gains and an optional immutable linkage map.
//! Reads borrow the string owner's coordinates; no string state or time
//! integrator is duplicated here.
//!
//! Provenance honesty: `bind` retains the authored height-dependent flux
//! gain and aperture window; its flux law is linearized and remains an
//! estimate for large excursions. `bind_flux_map` instead accepts a bounded,
//! caller-supplied gap-to-linkage map from a measurement or field bake.
//! It computes nonlinear Faraday voltage and reciprocal current force from
//! the same map, without granting its inputs experimental authority.
//! [`field`] compiles that same map from original material states, magnet
//! geometry and winding turns through the native planar field solver.
//!
//! This is a one-coordinate aperture reduction: displacement is positive
//! toward the pole and gap = rest height - window-weighted displacement.
//! The map describes current-independent permanent-magnet linkage, with
//! winding turns and polarity already included. Self-inductance, static
//! magnet attraction, eddy losses and hysteresis belong to their own owners.
//! The finite-step API closes the transducer work identity for a shared
//! step current; it does not make an arbitrary partitioned solve passive.

use fs_math::det;

/// Material-resolved magnetic fields compiled into the same Faraday observer.
pub mod field;

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

/// A sample and tangent of a measured or field-baked linkage map.
/// Slopes are with respect to increasing pole-to-string gap, not displacement.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PickupFluxKnot {
    /// Positive pole-to-string gap [m].
    pub gap_m: f64,
    /// Signed flux linkage [Wb turns = V s], including the coil's turns.
    pub linkage_wb: f64,
    /// Derivative of linkage with respect to gap [V s / m].
    pub slope_wb_per_m: f64,
}

/// Bounded C1 cubic Hermite linkage map. No extrapolation or slope fitting.
///
/// The caller owns the geometry, material, temperature, interpolation error
/// and validity of the supplied samples and slopes. This type establishes
/// numerical interpolation and reciprocity, not source or field authority.
#[derive(Debug, Clone, PartialEq)]
pub struct PickupFluxMap {
    knots: Vec<PickupFluxKnot>,
    /// Nonconstant coefficients in normalized segment coordinate [0, 1].
    coefficients: Vec<[f64; 3]>,
}

impl PickupFluxMap {
    /// Admit ordered finite samples with finite, representable interpolants.
    ///
    /// # Errors
    /// Refuses fewer than two knots, nonpositive or unordered gaps, nonfinite
    /// inputs, or coefficients that overflow. No monotonicity is invented.
    pub fn new(knots: Vec<PickupFluxKnot>) -> Result<Self, PickupError> {
        if knots.len() < 2
            || knots.iter().any(|k| {
                !k.gap_m.is_finite()
                    || k.gap_m <= 0.0
                    || !k.linkage_wb.is_finite()
                    || !k.slope_wb_per_m.is_finite()
            })
        {
            return Err(PickupError::Invalid {
                what: "linkage map needs at least two finite knots at positive gaps",
            });
        }
        let mut coefficients = Vec::with_capacity(knots.len() - 1);
        for pair in knots.windows(2) {
            let [a, b] = [pair[0], pair[1]];
            let width = b.gap_m - a.gap_m;
            if width <= 0.0 {
                return Err(PickupError::Invalid {
                    what: "linkage map gaps must strictly increase",
                });
            }
            let delta = b.linkage_wb - a.linkage_wb;
            let m0 = width * a.slope_wb_per_m;
            let m1 = width * b.slope_wb_per_m;
            let c = [m0, 3.0 * delta - 2.0 * m0 - m1, -2.0 * delta + m0 + m1];
            if c.iter().any(|v| !v.is_finite()) {
                return Err(PickupError::Invalid {
                    what: "linkage interpolant is unrepresentable",
                });
            }
            coefficients.push(c);
        }
        Ok(Self {
            knots,
            coefficients,
        })
    }

    fn segment(&self, gap: f64) -> Result<usize, PickupError> {
        if !gap.is_finite()
            || gap < self.knots[0].gap_m
            || gap > self.knots[self.knots.len() - 1].gap_m
        {
            return Err(PickupError::Invalid {
                what: "string gap lies outside the supplied linkage map",
            });
        }
        Ok(self
            .knots
            .partition_point(|k| k.gap_m <= gap)
            .saturating_sub(1)
            .min(self.coefficients.len() - 1))
    }

    fn sample(&self, gap: f64) -> Result<(f64, f64), PickupError> {
        let i = self.segment(gap)?;
        let a = self.knots[i];
        let b = self.knots[i + 1];
        let width = b.gap_m - a.gap_m;
        let t = (gap - a.gap_m) / width;
        let [c1, c2, c3] = self.coefficients[i];
        let (value, slope) = if gap == b.gap_m {
            (b.linkage_wb, b.slope_wb_per_m)
        } else {
            (
                a.linkage_wb + t * (c1 + t * (c2 + t * c3)),
                (c1 + t * (2.0 * c2 + 3.0 * t * c3)) / width,
            )
        };
        if !value.is_finite() || !slope.is_finite() {
            return Err(PickupError::Invalid {
                what: "linkage evaluation is unrepresentable",
            });
        }
        Ok((value, slope))
    }

    /// Analytic divided difference, including crossings of map knots.
    /// Cancelling the constant linkage algebraically avoids subtracting large
    /// DC flux values for a very small string motion.
    fn secant(&self, a: f64, b: f64) -> Result<f64, PickupError> {
        let lo = a.min(b);
        let hi = a.max(b);
        let first = self.segment(a)?;
        let last = self.segment(b)?;
        if a == b {
            return self.sample(a).map(|(_, slope)| slope);
        }
        let mut slope = 0.0;
        for i in first.min(last)..=first.max(last) {
            let knot = self.knots[i].gap_m;
            let width = self.knots[i + 1].gap_m - knot;
            let left = lo.max(knot);
            let right = hi.min(self.knots[i + 1].gap_m);
            if right <= left {
                continue;
            }
            let t0 = (left - knot) / width;
            let t1 = (right - knot) / width;
            let [c1, c2, c3] = self.coefficients[i];
            let local = (c1 + c2 * (t0 + t1) + c3 * (t0 * t0 + t0 * t1 + t1 * t1)) / width;
            slope += ((right - left) / (hi - lo)) * local;
        }
        if !slope.is_finite() {
            return Err(PickupError::Invalid {
                what: "linkage divided difference is unrepresentable",
            });
        }
        Ok(slope)
    }
}

/// A lumped Faraday pickup bound to a modal string basis.
#[derive(Debug, Clone, PartialEq)]
pub struct Pickup {
    pose: PickupPose,
    /// Per-mode window-weighted gains `g_k = ∫ w(x) φ_k(x) dx`.
    gains: Vec<f64>,
    /// Effective flux gain κ(h) [V s / m] at the rest position.
    flux_gain: f64,
    flux_map: Option<PickupFluxMap>,
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
        if !(pose.station_frac.is_finite() && pose.station_frac > 0.0 && pose.station_frac < 1.0) {
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
                    * det::cos(core::f64::consts::PI * (x - pose.station_frac) / half.max(1e-300));
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
            flux_map: None,
        })
    }

    /// Bind an explicit nonlinear linkage map to the same aperture/mode port.
    /// The map must describe this aperture and remain valid for its weighted
    /// displacement reduction. No spatially resolved field solve is implied.
    ///
    /// # Errors
    /// In addition to [`Self::bind`] refusals, the rest gap must lie in the map.
    pub fn bind_flux_map(
        pose: PickupPose,
        n_modes: usize,
        map: PickupFluxMap,
    ) -> Result<Self, PickupError> {
        let mut pickup = Self::bind(pose, n_modes)?;
        pickup.flux_gain = -map.sample(pose.height_m)?.1;
        pickup.flux_map = Some(map);
        Ok(pickup)
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
    /// than the bound mode count, nonfinite state, or a mapped pickup (which
    /// requires [`Self::emf_at_v`] and the actual positions).
    pub fn emf_v(&self, modal_velocities: &[f64]) -> Result<f64, PickupError> {
        if self.flux_map.is_some() {
            return Err(PickupError::Invalid {
                what: "mapped pickup requires modal positions",
            });
        }
        finite(-self.flux_gain * self.project(modal_velocities)?)
    }

    fn project(&self, values: &[f64]) -> Result<f64, PickupError> {
        if values.len() < self.gains.len() {
            return Err(PickupError::Invalid {
                what: "modal slice shorter than the bound mode count",
            });
        }
        let mut acc = 0.0;
        for (g, v) in self.gains.iter().zip(values) {
            acc += g * v;
        }
        finite(acc)
    }

    fn linkage_and_gain(&self, q: &[f64]) -> Result<(f64, f64), PickupError> {
        self.linkage_and_gain_at_displacement(self.project(q)?)
    }

    fn linkage_and_gain_at_displacement(
        &self,
        displacement: f64,
    ) -> Result<(f64, f64), PickupError> {
        finite(displacement)?;
        if let Some(map) = &self.flux_map {
            let (linkage, slope) = map.sample(self.pose.height_m - displacement)?;
            Ok((linkage, -slope))
        } else {
            Ok((finite(self.flux_gain * displacement)?, self.flux_gain))
        }
    }

    /// Signed linkage [V s] at the owner's modal displacements [m].
    /// The authored linear model uses zero linkage at the rest position.
    ///
    /// # Errors
    /// Refuses short/nonfinite state or motion outside the supplied map.
    pub fn flux_linkage_wb(&self, q: &[f64]) -> Result<f64, PickupError> {
        self.linkage_and_gain(q).map(|(linkage, _)| linkage)
    }

    /// Faraday EMF `-dλ/dt` [V] at the actual displacements and velocities.
    /// Positive coil current is defined so `emf * current` leaves the
    /// mechanical transducer and enters the electrical circuit.
    ///
    /// # Errors
    /// Refuses invalid modal state, unsupported gap, or nonfinite voltage.
    pub fn emf_at_v(&self, q: &[f64], velocities: &[f64]) -> Result<f64, PickupError> {
        self.emf_at_displacement(self.project(q)?, self.project(velocities)?)
    }

    fn emf_at_displacement(&self, displacement: f64, velocity: f64) -> Result<f64, PickupError> {
        let gain = self.linkage_and_gain_at_displacement(displacement)?.1;
        finite(-gain * finite(velocity)?)
    }

    /// Reciprocal generalized force `i ∂λ/∂q_k` [N] on the string modes.
    /// Writes the bound modes into caller-owned storage without allocation;
    /// trailing entries are unchanged. `Σ F_k v_k + emf * i = 0`.
    ///
    /// # Errors
    /// Refuses invalid state/current, unsupported gap, short output or overflow
    /// before modifying any force entry.
    pub fn modal_force_n(
        &self,
        q: &[f64],
        current_a: f64,
        out: &mut [f64],
    ) -> Result<(), PickupError> {
        self.write_force(self.linkage_and_gain(q)?.1, current_a, out)
    }

    fn write_force(&self, gain: f64, current: f64, out: &mut [f64]) -> Result<(), PickupError> {
        let scale = finite(gain * current)?;
        if out.len() < self.gains.len() || self.gains.iter().any(|g| !(scale * g).is_finite()) {
            return Err(PickupError::Invalid {
                what: "modal force output is short or unrepresentable",
            });
        }
        for (force, g) in out.iter_mut().zip(&self.gains) {
            *force = scale * g;
        }
        Ok(())
    }

    fn step_gain(&self, q0: &[f64], q1: &[f64]) -> Result<(f64, f64), PickupError> {
        let d0 = self.project(q0)?;
        let d1 = self.project(q1)?;
        let gain = if let Some(map) = &self.flux_map {
            -map.secant(self.pose.height_m - d0, self.pose.height_m - d1)?
        } else {
            self.flux_gain
        };
        Ok((finite(d1 - d0)?, gain))
    }

    /// Average EMF `-[λ(q1)-λ(q0)]/dt` for a finite displacement step.
    /// Evaluates the Hermite divided difference analytically, including knot
    /// crossings, without cancellation against a large constant flux offset.
    ///
    /// # Errors
    /// Refuses unsupported endpoint states, nonpositive/nonfinite dt or overflow.
    pub fn emf_over_step_v(&self, q0: &[f64], q1: &[f64], dt_s: f64) -> Result<f64, PickupError> {
        if !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(PickupError::Invalid {
                what: "pickup step duration must be positive and finite",
            });
        }
        let (displacement, gain) = self.step_gain(q0, q1)?;
        finite(-gain * displacement / dt_s)
    }

    /// Discrete-gradient force paired with [`Self::emf_over_step_v`].
    /// For the SAME held/step current, `Σ F_k Δq_k = -emf * i * dt` up to
    /// roundoff. Use these forces in the string owner's coupled solve; a
    /// prescribed trajectory alone does not acquire electrical back-reaction.
    ///
    /// # Errors
    /// Refuses invalid endpoints/current/output before changing any force entry.
    pub fn modal_step_force_n(
        &self,
        q0: &[f64],
        q1: &[f64],
        current_a: f64,
        out: &mut [f64],
    ) -> Result<(), PickupError> {
        self.write_force(self.step_gain(q0, q1)?.1, current_a, out)
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

fn finite(value: f64) -> Result<f64, PickupError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(PickupError::Invalid {
            what: "pickup state or response is nonfinite or unrepresentable",
        })
    }
}

#[cfg(test)]
mod nonlinear_pickup_tests {
    use super::*;

    const HEIGHT: f64 = 0.003;
    fn linkage(gap: f64, cubic: f64) -> f64 {
        let x = gap - HEIGHT;
        0.001 - 0.02 * x + 4.0 * x * x + cubic * x * x * x
    }
    fn mapped(cubic: f64) -> Pickup {
        let knots = [0.001, 0.0026, 0.0037, 0.005].map(|gap| {
            let x = gap - HEIGHT;
            PickupFluxKnot {
                gap_m: gap,
                linkage_wb: linkage(gap, cubic),
                slope_wb_per_m: -0.02 + 8.0 * x + 3.0 * cubic * x * x,
            }
        });
        Pickup::bind_flux_map(
            PickupPose {
                station_frac: 0.3,
                height_m: HEIGHT,
                aperture_frac: 0.04,
            },
            3,
            PickupFluxMap::new(knots.to_vec()).unwrap(),
        )
        .unwrap()
    }
    fn dot(a: &[f64], b: &[f64]) -> f64 {
        a.iter().zip(b).map(|(x, y)| x * y).sum()
    }

    #[test]
    fn supplied_linkage_obeys_faraday_and_loaded_lenz_reaction() {
        let pickup = mapped(-1000.0);
        let q = [0.0008, -0.0002, 0.0001];
        let v = [0.14, -0.03, 0.02];
        let gap = HEIGHT - dot(pickup.mode_gains(), &q);
        assert!((pickup.flux_linkage_wb(&q).unwrap() - linkage(gap, -1000.0)).abs() < 1e-17);
        let emf = pickup.emf_at_v(&q, &v).unwrap();
        let h = 1e-7;
        let speed = dot(pickup.mode_gains(), &v);
        let faraday =
            -(linkage(gap - h * speed, -1000.0) - linkage(gap + h * speed, -1000.0)) / (2.0 * h);
        assert!((emf - faraday).abs() < 2e-12, "{emf} vs {faraday}");
        for polarity in [-1.0, 1.0] {
            let velocity = v.map(|x| polarity * x);
            let emf = pickup.emf_at_v(&q, &velocity).unwrap();
            let resistance = 1000.0;
            let current = emf / resistance;
            let mut forces = [0.0; 3];
            pickup.modal_force_n(&q, current, &mut forces).unwrap();
            let mechanical_power = dot(&forces, &velocity);
            assert!(
                mechanical_power < 0.0,
                "a resistive load must damp either direction"
            );
            assert!((mechanical_power + emf * current).abs() < 1e-20);
        }
        assert!(
            pickup.emf_v(&v).is_err(),
            "nonlinear maps require positions"
        );
    }

    #[test]
    fn finite_excursion_crosses_knots_with_reciprocal_work_and_zero_step_limit() {
        let pickup = mapped(-1000.0);
        let q0 = [-0.0014, 0.0004, -0.0002];
        let q1 = [0.0013, -0.0005, 0.0003];
        let dt = 0.001;
        let current = 0.004;
        let gap0 = HEIGHT - dot(pickup.mode_gains(), &q0);
        let gap1 = HEIGHT - dot(pickup.mode_gains(), &q1);
        assert!(gap0 > 0.0037 && gap1 < 0.0026, "cross two internal knots");
        let emf = pickup.emf_over_step_v(&q0, &q1, dt).unwrap();
        let expected = -(linkage(gap1, -1000.0) - linkage(gap0, -1000.0)) / dt;
        assert!((emf - expected).abs() < 1e-14);
        let mut forces = [0.0; 3];
        pickup
            .modal_step_force_n(&q0, &q1, current, &mut forces)
            .unwrap();
        let dq = std::array::from_fn::<_, 3, _>(|i| q1[i] - q0[i]);
        assert!((dot(&forces, &dq) + emf * current * dt).abs() < 1e-20);
        assert_eq!(pickup.emf_over_step_v(&q1, &q0, dt).unwrap(), -emf);
        let mut reverse = [0.0; 3];
        pickup
            .modal_step_force_n(&q1, &q0, current, &mut reverse)
            .unwrap();
        assert_eq!(forces, reverse);
        pickup
            .modal_step_force_n(&q0, &q0, current, &mut forces)
            .unwrap();
        pickup.modal_force_n(&q0, current, &mut reverse).unwrap();
        assert_eq!(forces, reverse);
        assert_eq!(pickup.emf_over_step_v(&q0, &q0, dt).unwrap(), 0.0);
    }

    #[test]
    fn curvature_generates_the_predicted_second_harmonic_without_a_waveshaper() {
        let pickup = mapped(0.0);
        let amplitude = 0.0005;
        let omega = std::f64::consts::TAU * 100.0;
        let samples = 2048;
        let mut second = 0.0;
        for i in 0..samples {
            let phase = std::f64::consts::TAU * i as f64 / samples as f64;
            let q = [
                amplitude * det::cos(phase) / pickup.mode_gains()[0],
                0.0,
                0.0,
            ];
            let v = [
                -amplitude * omega * det::sin(phase) / pickup.mode_gains()[0],
                0.0,
                0.0,
            ];
            second += pickup.emf_at_v(&q, &v).unwrap() * det::sin(2.0 * phase);
        }
        second *= 2.0 / samples as f64;
        let analytic = 4.0 * amplitude * amplitude * omega;
        assert!((second - analytic).abs() < 1e-14, "{second} vs {analytic}");
    }

    #[test]
    fn tiny_motion_keeps_its_emf_despite_a_large_constant_linkage() {
        let map = PickupFluxMap::new(vec![
            PickupFluxKnot {
                gap_m: 1.0,
                linkage_wb: 1e12 - 1.0,
                slope_wb_per_m: 1.0,
            },
            PickupFluxKnot {
                gap_m: 3.0,
                linkage_wb: 1e12 + 1.0,
                slope_wb_per_m: 1.0,
            },
        ])
        .unwrap();
        let pickup = Pickup::bind_flux_map(
            PickupPose {
                station_frac: 0.3,
                height_m: 2.0,
                aperture_frac: 0.04,
            },
            1,
            map,
        )
        .unwrap();
        let q0 = [0.0];
        let q1 = [1e-8];
        assert_eq!(
            pickup.flux_linkage_wb(&q0).unwrap(),
            pickup.flux_linkage_wb(&q1).unwrap()
        );
        let expected = pickup.mode_gains()[0] * q1[0];
        assert_eq!(pickup.emf_over_step_v(&q0, &q1, 1.0).unwrap(), expected);
    }

    #[test]
    fn unsupported_motion_and_nonfinite_inputs_refuse_before_writing_force() {
        let pickup = mapped(0.0);
        let mut force = [9.0; 3];
        for q in [[0.01, 0.0, 0.0], [f64::NAN, 0.0, 0.0]] {
            assert!(pickup.emf_at_v(&q, &[0.1; 3]).is_err());
            assert!(pickup.modal_force_n(&q, 0.001, &mut force).is_err());
            assert_eq!(force, [9.0; 3]);
        }
        assert!(
            pickup
                .modal_force_n(&[0.0; 3], f64::NAN, &mut force)
                .is_err()
        );
        assert_eq!(force, [9.0; 3]);
        assert!(pickup.emf_at_v(&[0.0; 3], &[f64::INFINITY; 3]).is_err());
        assert!(pickup.emf_over_step_v(&[0.0; 3], &[0.0; 3], 0.0).is_err());
        assert!(PickupFluxMap::new(vec![]).is_err());
        let knot = PickupFluxKnot {
            gap_m: HEIGHT,
            linkage_wb: 0.001,
            slope_wb_per_m: 0.0,
        };
        assert!(PickupFluxMap::new(vec![knot, knot]).is_err());
        assert!(
            PickupFluxMap::new(vec![
                knot,
                PickupFluxKnot {
                    gap_m: 0.004,
                    slope_wb_per_m: f64::INFINITY,
                    ..knot
                }
            ])
            .is_err()
        );
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
                // Near-undamped: at 1e-4 the mode-3 amplitude decays
                // ~10% across the projection window and the line reads
                // 5% low (measured, exactly at the band edge).
                let zeta = 1.0e-5_f64;
                let wd = w * (1.0 - zeta * zeta).sqrt();
                let decay = det::exp(-zeta * w * dt);
                rot.push((decay * det::cos(wd * dt), decay * det::sin(wd * dt)));
                // Pluck shape: q_k(0) ~ sin(k pi x_p)/k^2.
                let q0 = det::sin(k as f64 * core::f64::consts::PI * pluck_frac) / (k * k) as f64;
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
                    let w = 0.5 + 0.5 * det::cos(core::f64::consts::PI * (x - station) / half);
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
                let expected = (pickup.flux_gain_v_s_per_m()
                    * analytic_gain(k)
                    * core::f64::consts::TAU
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
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"pk-003-point-limit\",\"verdict\":\"pass\"}}"
        );
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
