//! Quasi-steady planar conduction through two source-resolved k(T) layers.
//!
//! This constitutive flux closure does not integrate thermal storage. The
//! existing enthalpy owner receives its one signed heat transfer. Each layer
//! is homogeneous, isotropic and source-free internally, at fixed composition
//! and phase. Gray radiation may act between the two solved surface traces;
//! distributed generation, radiation and storage inside the material layers
//! require their own resolved equations.

use fs_thermal::enthalpy::{DepositionError, series_contact_conductance_w_k};

mod gray_gap;
pub use gray_gap::{GrayLayerConductionResponse, GrayLayerGap};

mod surface_exchange;
pub use surface_exchange::{SurfaceExchangeState, SurfaceLayerConductionResponse};

/// Maximum retained temperature/conductivity knots in one constitutive chart.
pub const MAX_CONDUCTIVITY_KNOTS: usize = 4096;

/// An explicitly declared, strictly positive piecewise-linear conductivity.
/// Knots are [absolute temperature in K, conductivity in W/(m K)]. Source
/// uncertainty is not inferred. No interpolation outside the knot span occurs.
#[derive(Clone, Debug, PartialEq)]
pub struct ConductivityCurve {
    knots: Vec<[f64; 2]>,
}

impl ConductivityCurve {
    /// Admit an ordered, bounded chart; zero conductivity is outside this
    /// conducting-layer model, not silently replaced with a small positive k.
    pub fn new(knots: Vec<[f64; 2]>) -> Result<Self, DepositionError> {
        if !(2..=MAX_CONDUCTIVITY_KNOTS).contains(&knots.len()) {
            return Err(DepositionError::Budget(
                "conductivity curve needs 2..=4096 knots",
            ));
        }
        if knots.iter().flatten().any(|v| !v.is_finite() || *v <= 0.0)
            || knots.windows(2).any(|p| p[0][0] >= p[1][0])
        {
            return Err(invalid(
                "conductivity knots need increasing positive kelvin and positive finite k",
            ));
        }
        Ok(Self { knots })
    }

    /// Retained, unmodified declared knots.
    pub fn knots(&self) -> &[[f64; 2]] {
        &self.knots
    }

    /// Inclusive absolute-temperature support.
    pub fn temperature_range_k(&self) -> [f64; 2] {
        [self.knots[0][0], self.knots[self.knots.len() - 1][0]]
    }

    fn segment(&self, t: f64) -> Result<usize, DepositionError> {
        let [lo, hi] = self.temperature_range_k();
        if !t.is_finite() || t < lo || t > hi {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(self
            .knots
            .partition_point(|k| k[0] <= t)
            .min(self.knots.len() - 1)
            - 1)
    }

    fn on_segment(&self, i: usize, t: f64) -> Result<f64, DepositionError> {
        let [a, ka] = self.knots[i];
        let [b, kb] = self.knots[i + 1];
        if t == a {
            return Ok(ka);
        }
        if t == b {
            return Ok(kb);
        }
        let fraction = (t - a) / (b - a);
        positive((1.0 - fraction) * ka + fraction * kb)
    }

    /// The original linear chart at one supported temperature; no allocation.
    pub fn conductivity_w_m_k(&self, t: f64) -> Result<f64, DepositionError> {
        self.on_segment(self.segment(t)?, t)
    }

    /// Conductivity and its temperature derivative on the same original
    /// segment, in W/(m K) and W/(m K²). Interior knots select the right
    /// derivative; the upper endpoint selects the left derivative. The
    /// residual still uses the continuous original conductivity chart.
    pub fn conductivity_and_temperature_derivative(
        &self,
        t: f64,
    ) -> Result<(f64, f64), DepositionError> {
        let i = self.segment(t)?;
        let [a, ka] = self.knots[i];
        let [b, kb] = self.knots[i + 1];
        let derivative = finite((kb - ka) / (b - a))?;
        Ok((self.on_segment(i, t)?, derivative))
    }

    /// Signed Kirchhoff difference integral k(T) dT, in W/m. Integrate crossed
    /// segments directly rather than subtracting large absolute primitives.
    /// The trapezoid is exact for this linear k chart, not a quadrature claim
    /// about an underlying unsampled material. Reversal is exactly antisymmetric.
    /// At most 4095 segments are visited; successful evaluation allocates nothing.
    pub fn integral_w_m(&self, initial: f64, final_k: f64) -> Result<f64, DepositionError> {
        let first = self.segment(initial)?;
        let last = self.segment(final_k)?;
        if initial == final_k {
            return Ok(0.0);
        }
        let (lo, hi, first, last, sign) = if initial < final_k {
            (initial, final_k, first, last, 1.0)
        } else {
            (final_k, initial, last, first, -1.0)
        };
        let mut sum = 0.0;
        let mut correction = 0.0;
        for i in first..=last {
            let a = lo.max(self.knots[i][0]);
            let b = hi.min(self.knots[i + 1][0]);
            if a < b {
                let mean = 0.5 * self.on_segment(i, a)? + 0.5 * self.on_segment(i, b)?;
                let term = positive(mean * (b - a))?;
                let y = term - correction;
                let next = finite(sum + y)?;
                correction = (next - sum) - y;
                sum = next;
            }
        }
        Ok(sign * positive(sum)?)
    }

    /// Exact mean of the original piecewise-affine conductivity over a
    /// supported temperature interval. Equal endpoints return the point law.
    ///
    /// Accumulate positive segment means by their temperature widths. Convex
    /// updates avoid an overflowing or underflowing dimensional Kirchhoff
    /// integral when its mean is still representable. No source knot is
    /// skipped, no table is copied, and reversal uses identical arithmetic.
    pub fn mean_conductivity_w_m_k(
        &self,
        initial: f64,
        final_k: f64,
    ) -> Result<f64, DepositionError> {
        let first = self.segment(initial)?;
        let last = self.segment(final_k)?;
        if initial == final_k {
            return self.on_segment(first, initial);
        }
        let (lo, hi, first, last) = if initial < final_k {
            (initial, final_k, first, last)
        } else {
            (final_k, initial, last, first)
        };
        let mut mean: Option<f64> = None;
        for i in first..=last {
            let a = lo.max(self.knots[i][0]);
            let b = hi.min(self.knots[i + 1][0]);
            if a < b {
                let ka = self.on_segment(i, a)?;
                let kb = self.on_segment(i, b)?;
                let segment_mean = 0.5_f64.mul_add(kb - ka, ka);
                mean = Some(match mean {
                    None => segment_mean,
                    Some(previous) => {
                        let fraction = (b - a) / (b - lo);
                        (segment_mean - previous).mul_add(fraction, previous)
                    }
                });
            }
        }
        positive(mean.ok_or(DepositionError::Arithmetic)?)
    }

    /// Invert a signed Kirchhoff difference on this same conductivity chart.
    /// Positive integral raises temperature. Each crossed affine segment has
    /// a scaled quadratic inverse; no iterative solve or extrapolation occurs.
    /// Relative integration retains small differences near a large temperature.
    pub fn temperature_after_integral_w_m(
        &self,
        initial: f64,
        integral: f64,
    ) -> Result<f64, DepositionError> {
        self.inverse_integral_w_m(initial, integral)
            .map(|(temperature, _)| temperature)
    }

    // Retain the inverse's ORIGINAL relative temperature displacement when
    // checking its integral. Reintegrating from the rounded absolute endpoint
    // loses the small drop of a thin, highly conducting layer; k/d then
    // amplifies that Kelvin rounding into a false heat-continuity failure.
    // The returned integral is recomputed from the solved affine-segment
    // displacement, not copied from the requested value or relaxed by a
    // temperature-dependent heat tolerance.
    pub(super) fn inverse_integral_w_m(
        &self,
        initial: f64,
        integral: f64,
    ) -> Result<(f64, f64), DepositionError> {
        self.segment(initial)?;
        if !integral.is_finite() {
            return Err(invalid("conductivity integral must be finite"));
        }
        if integral == 0.0 {
            return Ok((initial, 0.0));
        }
        let heating = integral > 0.0;
        let span = self.temperature_range_k();
        let boundary = if heating { span[1] } else { span[0] };
        let target = integral.abs();
        if let Ok(available) = self.integral_w_m(initial, boundary) {
            if target > available.abs() {
                return Err(DepositionError::MaterialDomain);
            }
            if target == available.abs() {
                return Ok((boundary, available));
            }
        }
        let mut index = if heating {
            self.segment(initial)?
        } else {
            self.knots
                .partition_point(|k| k[0] < initial)
                .saturating_sub(1)
        };
        let mut current = initial;
        let mut cumulative = 0.0;
        let mut correction = 0.0;
        for _ in 0..self.knots.len() - 1 {
            let end = self.knots[if heating { index + 1 } else { index }][0];
            let a = self.on_segment(index, current)?;
            let b = self.on_segment(index, end)?;
            let mean = 0.5 * a + 0.5 * b;
            let width = (end - current).abs();
            let amount = mean * width;
            if amount.is_nan() || amount <= 0.0 {
                return Err(DepositionError::Arithmetic);
            }
            let corrected = amount - correction;
            let next = cumulative + corrected;
            if target <= next {
                if target == next {
                    return Ok((end, integral.signum() * next));
                }
                let energy_fraction = if amount.is_finite() {
                    (target - cumulative) / amount
                } else {
                    ((target - cumulative) / mean) / width
                };
                let scale = a.max(b);
                let (scaled_a, scaled_b) = (a / scale, b / scale);
                let root = ((1.0 - energy_fraction) * scaled_a * scaled_a
                    + energy_fraction * scaled_b * scaled_b)
                    .sqrt();
                let fraction = energy_fraction * (scaled_a + scaled_b) / (scaled_a + root);
                let displacement = fraction * (end - current);
                let t = current + displacement;
                if !fraction.is_finite()
                    || !(0.0..=1.0).contains(&fraction)
                    || !t.is_finite()
                    || t == initial
                    || t < current.min(end)
                    || t > current.max(end)
                {
                    return Err(DepositionError::Arithmetic);
                }
                // The mean over the SOLVED fraction of this original segment
                // is k_start + fraction*(k_end-k_start)/2. Its convex form
                // avoids overflowing positive endpoint conductivities.
                let half_fraction = 0.5 * fraction;
                let mean = positive((1.0 - half_fraction) * a + half_fraction * b)?;
                let partial = positive(mean * displacement.abs())?;
                let actual = positive(finite(cumulative + (partial - correction))?)?;
                return Ok((t, integral.signum() * actual));
            }
            correction = (next - cumulative) - corrected;
            cumulative = next;
            current = end;
            if current == boundary {
                break;
            }
            if heating {
                index += 1;
            } else {
                index -= 1;
            }
        }
        Err(DepositionError::MaterialDomain)
    }
}

/// Explicit scalar work and heat-continuity allowances, independent of spatial
/// model error or the consuming thermal integrator's energy tolerance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayerConductionConfig {
    pub max_iterations: usize,
    pub absolute_heat_tolerance_w: f64,
    pub relative_heat_tolerance: f64,
}

/// Complete physical interface solution. Positive power leaves a and enters b.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayerConductionResponse {
    pub power_a_to_b_w: f64,
    /// Secant lowering for an existing G*(Ta-Tb) link. This is not its Newton
    /// derivative when k depends on temperature; use the derivatives below.
    pub conductance_w_k: f64,
    /// Temperatures on the a/b sides of the contact, in original layer order.
    pub interface_temperature_k: [f64; 2],
    /// [d(power a->b)/dTa, d(power a->b)/dTb] with charts/geometry fixed.
    pub endpoint_derivative_w_k: [f64; 2],
    /// d(power a->b)/dR'' with endpoint temperatures/charts/geometry fixed,
    /// in W²/(m² K). This is the right derivative at perfect contact, R''=0.
    /// A state-dependent contact law chains this with its own dR''/dh.
    pub contact_resistance_derivative_w2_m2_k: f64,
    /// Hot-layer outgoing minus cold-layer incoming power, independent of a/b order.
    /// Each inverse retains its relative affine-segment drop for this audit;
    /// reconstructing heat from rounded absolute interface temperatures can
    /// lose accuracy for thin conducting layers.
    pub layer_heat_residual_w: f64,
    /// Diagnostic temperature-jump residual divided by resistance, W. This is
    /// ill-conditioned near perfect contact, NOT an energy defect. Acceptance
    /// uses the Kelvin residual and explicit representation allowance below.
    pub contact_heat_residual_w: f64,
    /// Actual (T_interface_a-T_interface_b)-R*Q/A, in K.
    pub contact_temperature_residual_k: f64,
    /// Floating-point representation allowance for the contact jump, K;
    /// separate from the caller's physical heat-continuity tolerance.
    pub contact_roundoff_allowance_k: f64,
    pub iterations: usize,
}

/// Two original k(T) charts and fixed planar layer/contact geometry.
/// The thermal network still owns mass, enthalpy, clocks and source accounting.
#[derive(Clone, Copy, Debug)]
pub struct TwoLayerConduction<'a> {
    curves: [&'a ConductivityCurve; 2],
    distances_m: [f64; 2],
    area_m2: f64,
    contact_resistance_m2_k_w: f64,
}

// All interface laws use the SAME supported Kirchhoff inverses, heat checks
// and safeguarded scalar-flux iteration. Only their physical jump relation
// differs. Fluxes and temperatures below are ordered hot/cold, not A/B.
struct LayerFluxTrial {
    flux_w_m2: f64,
    surface_temperature_k: [f64; 2],
    layer_flux_w_m2: [f64; 2],
    iterations: usize,
}

enum LayerFluxUpdate<R> {
    Accepted(R),
    Continue {
        decreasing_residual: f64,
        slope_magnitude: f64,
    },
    // General passive sources can have a nonmonotone differential response.
    // Keep the supported sign bracket without taking an invalid Newton step.
    Bisect {
        decreasing_residual: f64,
    },
}

impl<'a> TwoLayerConduction<'a> {
    pub fn new(
        curves: [&'a ConductivityCurve; 2],
        distances_m: [f64; 2],
        area_m2: f64,
        contact_resistance_m2_k_w: f64,
    ) -> Result<Self, DepositionError> {
        if distances_m.iter().any(|d| !d.is_finite() || *d <= 0.0)
            || !area_m2.is_finite()
            || area_m2 <= 0.0
            || !contact_resistance_m2_k_w.is_finite()
            || contact_resistance_m2_k_w < 0.0
        {
            return Err(invalid(
                "conducting layers require positive lengths/area and nonnegative contact resistance",
            ));
        }
        Ok(Self {
            curves,
            distances_m,
            area_m2,
            contact_resistance_m2_k_w,
        })
    }

    /// Solve the steady layer flux using exact chart integrals and a safeguarded
    /// scalar Newton solve. Each chart covers only its own physical layer, from
    /// its bulk endpoint to its solved contact face. The two supports need not
    /// cover the opposite bulk endpoint or overlap across a finite contact jump.
    /// A source-limited flux bracket refuses an unsupported contact solution;
    /// no trial extrapolates either chart or infers an internal phase transition.
    ///
    /// Contact resistance is area-specific; its temperature jump is j*R, not
    /// Q*R. Source-free layers share one heat flux and contact adds no storage.
    /// Derivatives follow the accepted physical equations, not Newton iterations.
    pub fn evaluate(
        &self,
        temperatures_k: [f64; 2],
        config: LayerConductionConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<LayerConductionResponse, DepositionError> {
        checkpoint()?;
        self.admit_evaluation(temperatures_k, config)?;
        let [ta, tb] = temperatures_k;
        if ta == tb {
            let k = [
                self.curves[0].conductivity_w_m_k(ta)?,
                self.curves[1].conductivity_w_m_k(tb)?,
            ];
            let g = series_contact_conductance_w_k(
                self.area_m2,
                self.distances_m,
                k,
                self.contact_resistance_m2_k_w,
            )?;
            checkpoint()?;
            return Ok(LayerConductionResponse {
                power_a_to_b_w: 0.0,
                conductance_w_k: g,
                interface_temperature_k: [ta, tb],
                endpoint_derivative_w_k: [g, -g],
                contact_resistance_derivative_w2_m2_k: 0.0,
                layer_heat_residual_w: 0.0,
                contact_heat_residual_w: 0.0,
                contact_temperature_residual_k: 0.0,
                contact_roundoff_allowance_k: 0.0,
                iterations: 0,
            });
        }
        let (hot, cold, sign) = if ta > tb { (0, 1, 1.0) } else { (1, 0, -1.0) };
        let (th, tc) = (temperatures_k[hot], temperatures_k[cold]);
        let (kh, kc) = (self.curves[hot], self.curves[cold]);
        let (dh, dc) = (self.distances_m[hot], self.distances_m[cold]);
        let r = self.contact_resistance_m2_k_w;
        let maximum_flux = if r > 0.0 {
            Some(positive((th - tc) / r)?)
        } else {
            None
        };
        self.solve_interface(
            temperatures_k,
            config,
            maximum_flux,
            &mut checkpoint,
            |trial| {
                let [u, v] = trial.surface_temperature_k;
                let [jh, jc] = trial.layer_flux_w_m2;
                let jump = finite(r * jh)?;
                let qh = finite(self.area_m2 * jh)?;
                let qc = finite(self.area_m2 * jc)?;
                let residual = finite(qh - qc)?;
                let jump_residual = finite((u - v) - jump)?;
                let jump_roundoff =
                    finite(8.0 * f64::EPSILON * u.abs().max(v.abs()).max(jump.abs()))?;
                let contact_residual = if r == 0.0 {
                    0.0
                } else {
                    finite(self.area_m2 * ((u - v) / r - jh))?
                };
                let tolerance = finite(
                    config.absolute_heat_tolerance_w
                        + config.relative_heat_tolerance * qh.abs().max(qc.abs()),
                )?;
                // Dividing a representational Kelvin error by a vanishing contact
                // resistance would incorrectly reject the perfect-contact limit.
                // Keep that diagnostic separate; never relax the layer heat gate.
                let jump_tolerance = finite(jump_roundoff + (tolerance / self.area_m2) * r)?;
                if residual.abs() <= tolerance && jump_residual.abs() <= jump_tolerance {
                    let mut interface = [0.0; 2];
                    interface[hot] = u;
                    interface[cold] = v;
                    let ki = [
                        self.curves[0].conductivity_w_m_k(interface[0])?,
                        self.curves[1].conductivity_w_m_k(interface[1])?,
                    ];
                    let g_interface =
                        series_contact_conductance_w_k(self.area_m2, self.distances_m, ki, r)?;
                    let derivatives = [
                        positive(g_interface * (self.curves[0].conductivity_w_m_k(ta)? / ki[0]))?,
                        -positive(g_interface * (self.curves[1].conductivity_w_m_k(tb)? / ki[1]))?,
                    ];
                    // Differentiating the common-flux/contact equations gives
                    // dQ/dR'' = -Q/(d_a/k(Tia) + R'' + d_b/k(Tib)). Use the
                    // already accepted flux density to avoid overflowing Q*G
                    // before division by area, or forming a second interface solve.
                    let contact_derivative = -sign * positive(jh * g_interface)?;
                    // Require a representable nonzero flux even if a loose user
                    // tolerance would otherwise accept an endpoint of the bracket.
                    let g = positive(qh / (th - tc))?;
                    return Ok(LayerFluxUpdate::Accepted(LayerConductionResponse {
                        power_a_to_b_w: sign * qh,
                        conductance_w_k: g,
                        interface_temperature_k: interface,
                        endpoint_derivative_w_k: derivatives,
                        contact_resistance_derivative_w2_m2_k: contact_derivative,
                        layer_heat_residual_w: residual,
                        contact_heat_residual_w: sign * contact_residual,
                        contact_temperature_residual_k: sign * jump_residual,
                        contact_roundoff_allowance_k: jump_roundoff,
                        iterations: trial.iterations,
                    }));
                }
                let gap = finite((u - v) - finite(r * trial.flux_w_m2)?)?;
                let differential_resistance = positive(
                    positive(dh / kh.conductivity_w_m_k(u)?)?
                        + positive(dc / kc.conductivity_w_m_k(v)?)?
                        + r,
                )?;
                Ok(LayerFluxUpdate::Continue {
                    decreasing_residual: gap,
                    slope_magnitude: differential_resistance,
                })
            },
        )
    }

    fn admit_evaluation(
        &self,
        temperatures_k: [f64; 2],
        config: LayerConductionConfig,
    ) -> Result<(), DepositionError> {
        if config.max_iterations == 0
            || !config.absolute_heat_tolerance_w.is_finite()
            || config.absolute_heat_tolerance_w < 0.0
            || !config.relative_heat_tolerance.is_finite()
            || config.relative_heat_tolerance < 0.0
        {
            return Err(invalid(
                "invalid layer-conduction heat tolerance or iteration budget",
            ));
        }
        for (curve, t) in self.curves.into_iter().zip(temperatures_k) {
            curve.conductivity_w_m_k(t)?;
        }
        Ok(())
    }

    fn solve_interface<R>(
        &self,
        temperatures_k: [f64; 2],
        config: LayerConductionConfig,
        maximum_flux_w_m2: Option<f64>,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
        evaluate: impl FnMut(LayerFluxTrial) -> Result<LayerFluxUpdate<R>, DepositionError>,
    ) -> Result<R, DepositionError> {
        self.solve_interface_with_support(
            temperatures_k,
            config,
            maximum_flux_w_m2,
            None,
            checkpoint,
            evaluate,
        )
    }

    fn solve_interface_with_support<R>(
        &self,
        temperatures_k: [f64; 2],
        config: LayerConductionConfig,
        maximum_flux_w_m2: Option<f64>,
        surface_temperature_bounds_k: Option<[[f64; 2]; 2]>,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
        mut evaluate: impl FnMut(LayerFluxTrial) -> Result<LayerFluxUpdate<R>, DepositionError>,
    ) -> Result<R, DepositionError> {
        let (hot, cold) = if temperatures_k[0] > temperatures_k[1] {
            (0, 1)
        } else {
            (1, 0)
        };
        let (th, tc) = (temperatures_k[hot], temperatures_k[cold]);
        let (kh, kc) = (self.curves[hot], self.curves[cold]);
        let (dh, dc) = (self.distances_m[hot], self.distances_m[cold]);
        // Each original source covers only its own physical layer, including
        // when a resistive/radiative interface bridges disjoint source spans.
        let hot_limit = kh.temperature_range_k()[0].max(tc);
        let cold_limit = kc.temperature_range_k()[1].min(th);
        let (hot_limit, cold_limit, hot_near, cold_near) =
            if let Some(bounds) = surface_temperature_bounds_k {
                (
                    hot_limit.max(bounds[hot][0]),
                    cold_limit.min(bounds[cold][1]),
                    th.min(bounds[hot][1]),
                    tc.max(bounds[cold][0]),
                )
            } else {
                (hot_limit, cold_limit, th, tc)
            };
        if hot_limit > hot_near || cold_near > cold_limit {
            return Err(DepositionError::MaterialDomain);
        }
        let hot_integral = kh.integral_w_m(hot_limit, th)?;
        let cold_integral = kc.integral_w_m(tc, cold_limit)?;
        if surface_temperature_bounds_k.is_none() && (hot_integral == 0.0 || cold_integral == 0.0) {
            return Err(DepositionError::MaterialDomain);
        }
        let hot_flux_limit = if hot_integral == 0.0 {
            0.0
        } else {
            positive(hot_integral / dh)?
        };
        let cold_flux_limit = if cold_integral == 0.0 {
            0.0
        } else {
            positive(cold_integral / dc)?
        };
        let mut hi = hot_flux_limit.min(cold_flux_limit);
        if let Some(limit) = maximum_flux_w_m2 {
            hi = hi.min(limit);
        }
        // A source interval can exclude either bulk temperature. Map BOTH
        // interval ends through the same Kirchhoff integrals; q=0 must not be
        // queried when the original surface source does not support it.
        let (hot_minimum, cold_minimum) = if surface_temperature_bounds_k.is_some() {
            let hot = kh.integral_w_m(hot_near, th)?;
            let cold = kc.integral_w_m(tc, cold_near)?;
            (
                if hot == 0.0 { 0.0 } else { positive(hot / dh)? },
                if cold == 0.0 {
                    0.0
                } else {
                    positive(cold / dc)?
                },
            )
        } else {
            (0.0, 0.0)
        };
        let mut lo = hot_minimum.max(cold_minimum);
        if lo > hi {
            return Err(DepositionError::MaterialDomain);
        }
        // The supported upper endpoint distinguishes a supported root from
        // source exhaustion, including a root exactly at a source boundary.
        let mut flux = hi;
        for iteration in 0..config.max_iterations {
            checkpoint()?;
            // Keep source endpoints exact across division/multiplication.
            let (u, hot_integral) = if flux == hot_flux_limit {
                (hot_limit, -hot_integral)
            } else if surface_temperature_bounds_k.is_some() && flux == hot_minimum {
                (hot_near, kh.integral_w_m(th, hot_near)?)
            } else {
                kh.inverse_integral_w_m(th, -finite(flux * dh)?)?
            };
            let (v, cold_integral) = if flux == cold_flux_limit {
                (cold_limit, cold_integral)
            } else if surface_temperature_bounds_k.is_some() && flux == cold_minimum {
                (cold_near, kc.integral_w_m(tc, cold_near)?)
            } else {
                kc.inverse_integral_w_m(tc, finite(flux * dc)?)?
            };
            if let Some(bounds) = surface_temperature_bounds_k {
                if !(bounds[hot][0]..=bounds[hot][1]).contains(&u)
                    || !(bounds[cold][0]..=bounds[cold][1]).contains(&v)
                {
                    return Err(DepositionError::MaterialDomain);
                }
            }
            let update = evaluate(LayerFluxTrial {
                flux_w_m2: flux,
                surface_temperature_k: [u, v],
                layer_flux_w_m2: [finite(-hot_integral / dh)?, finite(cold_integral / dc)?],
                iterations: iteration + 1,
            })?;
            if surface_temperature_bounds_k.is_some() {
                checkpoint()?;
            }
            let (gap, slope) = match update {
                LayerFluxUpdate::Accepted(response) => {
                    checkpoint()?;
                    return Ok(response);
                }
                LayerFluxUpdate::Continue {
                    decreasing_residual,
                    slope_magnitude,
                } => (
                    finite(decreasing_residual)?,
                    Some(positive(slope_magnitude)?),
                ),
                LayerFluxUpdate::Bisect {
                    decreasing_residual,
                } => (finite(decreasing_residual)?, None),
            };
            if iteration == 0 && gap > 0.0 {
                return Err(DepositionError::MaterialDomain);
            }
            if surface_temperature_bounds_k.is_some() {
                if iteration == 0 {
                    if lo == hi {
                        return Err(DepositionError::MaterialDomain);
                    }
                    // The lower source endpoint can be strictly positive. It
                    // must establish the other sign before any interior trial.
                    flux = lo;
                    continue;
                }
                if iteration == 1 && gap < 0.0 {
                    return Err(DepositionError::MaterialDomain);
                }
            }
            if gap > 0.0 {
                lo = flux;
            } else {
                hi = flux;
            }
            let newton = match slope {
                Some(slope) => flux + gap / slope,
                None => f64::NAN,
            };
            let width = hi - lo;
            // Periodic bisection guarantees bracket contraction without
            // rejecting a rapidly converging Newton step near one bracket end.
            let next = if newton.is_finite() && newton > lo && newton < hi && iteration % 8 != 7 {
                newton
            } else {
                lo + 0.5 * width
            };
            if next == flux || next <= lo || next >= hi {
                return Err(DepositionError::Arithmetic);
            }
            flux = next;
        }
        Err(DepositionError::Budget("layer-conduction iterations"))
    }
}

fn invalid(message: &'static str) -> DepositionError {
    DepositionError::InvalidInput(message)
}
fn finite(value: f64) -> Result<f64, DepositionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic)
    }
}
fn positive(value: f64) -> Result<f64, DepositionError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic)
    }
}

#[cfg(test)]
mod interval_mean_tests {
    use super::ConductivityCurve;

    #[test]
    fn source_interval_mean_retains_knots_endpoints_and_original_support() {
        let curve =
            ConductivityCurve::new(vec![[300.0, 2.0], [350.0, 6.0], [400.0, 4.0], [500.0, 8.0]])
                .unwrap();
        // Independent segment areas: 125 + 250 + 412.5, over 150 kelvin.
        let mean = curve.mean_conductivity_w_m_k(325.0, 475.0).unwrap();
        assert!((mean - 5.25).abs() < 2e-15);
        assert_eq!(mean, curve.mean_conductivity_w_m_k(475.0, 325.0).unwrap());
        assert!((mean * 150.0 - curve.integral_w_m(325.0, 475.0).unwrap()).abs() < 2e-12);
        assert_eq!(curve.mean_conductivity_w_m_k(350.0, 350.0).unwrap(), 6.0);
        assert!(curve.mean_conductivity_w_m_k(299.0, 400.0).is_err());
        assert!(curve.mean_conductivity_w_m_k(300.0, f64::NAN).is_err());
    }

    #[test]
    fn representable_mean_survives_unrepresentable_dimensional_integral() {
        let high = ConductivityCurve::new(vec![[1.0, 1e308], [2.0, 1e308], [4.0, 1e308]]).unwrap();
        assert!(high.integral_w_m(1.0, 4.0).is_err());
        assert_eq!(high.mean_conductivity_w_m_k(1.0, 4.0).unwrap(), 1e308);
        let low =
            ConductivityCurve::new(vec![[1e-100, 1e-300], [2e-100, 1e-300], [3e-100, 1e-300]])
                .unwrap();
        assert!(low.integral_w_m(1e-100, 3e-100).is_err());
        assert_eq!(low.mean_conductivity_w_m_k(1e-100, 3e-100).unwrap(), 1e-300);
    }
}
