//! Quasi-steady planar conduction through two source-resolved k(T) layers.
//!
//! This constitutive flux closure does not integrate thermal storage. The
//! existing enthalpy owner receives its one signed heat transfer. Each layer
//! is homogeneous, isotropic and source-free internally, at fixed composition
//! and phase; radiation, distributed generation and storage inside the layers
//! require their own resolved equations.

use fs_thermal::enthalpy::{DepositionError, series_contact_conductance_w_k};

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
    /// Hot-layer outgoing minus cold-layer incoming power, independent of a/b order.
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
    /// scalar Newton solve. Both charts must cover the full endpoint temperature
    /// interval. This deliberately conservative support requirement avoids any
    /// extrapolated trial or silently inferred phase transition within a layer.
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
        for curve in self.curves {
            for t in temperatures_k {
                curve.conductivity_w_m_k(t)?;
            }
        }
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
        let (mut lo, mut hi) = (tc, th);
        let mut u = lo + 0.5 * (hi - lo);
        for iteration in 0..config.max_iterations {
            checkpoint()?;
            let jh = finite(kh.integral_w_m(u, th)? / dh)?;
            let jump = finite(r * jh)?;
            let v = finite(u - jump)?;
            // An interface below the cold boundary is only a bracket sign, not
            // a material evaluation. It implies the physical root lies right.
            if v < tc {
                lo = u;
                let next = lo + 0.5 * (hi - lo);
                if next == u {
                    return Err(DepositionError::Arithmetic);
                }
                u = next;
                continue;
            }
            let jc = finite(kc.integral_w_m(tc, v)? / dc)?;
            let qh = finite(self.area_m2 * jh)?;
            let qc = finite(self.area_m2 * jc)?;
            let residual = finite(qh - qc)?;
            let jump_residual = finite((u - v) - jump)?;
            let jump_roundoff = finite(8.0 * f64::EPSILON * u.abs().max(v.abs()).max(jump.abs()))?;
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
                // Require a representable nonzero flux even if a loose user
                // tolerance would otherwise accept an endpoint of the bracket.
                let g = positive(qh / (th - tc))?;
                checkpoint()?;
                return Ok(LayerConductionResponse {
                    power_a_to_b_w: sign * qh,
                    conductance_w_k: g,
                    interface_temperature_k: interface,
                    endpoint_derivative_w_k: derivatives,
                    layer_heat_residual_w: residual,
                    contact_heat_residual_w: sign * contact_residual,
                    contact_temperature_residual_k: sign * jump_residual,
                    contact_roundoff_allowance_k: jump_roundoff,
                    iterations: iteration + 1,
                });
            }
            if residual > 0.0 {
                lo = u;
            } else {
                hi = u;
            }
            let a = positive(kh.conductivity_w_m_k(u)? / dh)?;
            let b = positive(kc.conductivity_w_m_k(v)? / dc)?;
            let derivative = positive(a + b * (1.0 + r * a))?;
            let newton = u + (jh - jc) / derivative;
            let width = hi - lo;
            // Periodic bisection guarantees bracket contraction without
            // rejecting a rapidly converging Newton step near one bracket end.
            let next = if newton.is_finite() && newton > lo && newton < hi && iteration % 8 != 7 {
                newton
            } else {
                lo + 0.5 * width
            };
            if next == u || next <= lo || next >= hi {
                return Err(DepositionError::Arithmetic);
            }
            u = next;
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
