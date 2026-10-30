//! Ordered planar stacks using the original conductivity charts and heat owner.
//!
//! Shooting in one shared heat flux avoids a dense interface-temperature solve.
//! Every layer uses its exact piecewise-affine Kirchhoff integral; each contact
//! contributes its explicit area-specific temperature jump. No internal layer
//! storage, phase change, radiation, generation or lateral heat flow is inferred.

use super::{ConductivityCurve, LayerConductionConfig};
use fs_thermal::enthalpy::DepositionError;

/// Fixed response/workspace bound; successful evaluation allocates nothing.
pub const MAX_CONDUCTION_LAYERS: usize = 64;

/// One accepted flux and the original ordered layer-face temperatures.
#[derive(Clone, Debug, PartialEq)]
pub struct LayerStackResponse {
    pub power_a_to_b_w: f64,
    /// Secant coefficient for the original conservative thermal link.
    pub conductance_w_k: f64,
    /// [dQ/dTa, dQ/dTb], differentiating physical equations at fixed geometry.
    pub endpoint_derivative_w_k: [f64; 2],
    pub max_layer_heat_residual_w: f64,
    /// Contact-jump residual in kelvin, never divided by a tiny resistance.
    pub max_contact_temperature_residual_k: f64,
    /// Representation allowance, separate from physical heat-continuity error.
    pub max_contact_roundoff_allowance_k: f64,
    pub iterations: usize,
    layer_count: usize,
    temperatures: [[f64; 2]; MAX_CONDUCTION_LAYERS],
}

impl LayerStackResponse {
    /// Each pair is [left face, right face] in the input's original layer order.
    pub fn layer_temperatures_k(&self) -> &[[f64; 2]] {
        &self.temperatures[..self.layer_count]
    }
}

/// Borrowed, source-free, one-dimensional layers with explicitly aligned scalar k.
/// There are N thicknesses and N-1 contact resistances, all in physical order.
#[derive(Clone, Copy, Debug)]
pub struct LayerStackConduction<'a> {
    curves: &'a [ConductivityCurve],
    thicknesses: &'a [f64],
    contacts: &'a [f64],
    area: f64,
}

struct Shot {
    temperatures: [[f64; 2]; MAX_CONDUCTION_LAYERS],
    last_flux: f64,
    last_index: usize,
    // Sensitivities of the last layer's hot inlet: dT/dTh and -dT/dj.
    inlet_temperature_slope: f64,
    inlet_flux_slope: f64,
    last_inlet_k: f64,
}

impl<'a> LayerStackConduction<'a> {
    pub fn new(
        curves: &'a [ConductivityCurve],
        thicknesses_m: &'a [f64],
        contact_resistances_m2_k_w: &'a [f64],
        area_m2: f64,
    ) -> Result<Self, DepositionError> {
        let n = curves.len();
        if n == 0 || n > MAX_CONDUCTION_LAYERS {
            return Err(DepositionError::Budget(
                "conducting stack needs 1..=64 layers",
            ));
        }
        if thicknesses_m.len() != n
            || contact_resistances_m2_k_w.len() != n - 1
            || thicknesses_m.iter().any(|x| !x.is_finite() || *x <= 0.0)
            || contact_resistances_m2_k_w
                .iter()
                .any(|x| !x.is_finite() || *x < 0.0)
            || !area_m2.is_finite()
            || area_m2 <= 0.0
        {
            return Err(invalid(
                "stack geometry needs positive lengths/area and N-1 nonnegative contacts",
            ));
        }
        Ok(Self {
            curves,
            thicknesses: thicknesses_m,
            contacts: contact_resistances_m2_k_w,
            area: area_m2,
        })
    }

    /// Solve a single monotone flux equation with safeguarded Newton/bisection.
    /// All charts must cover the entire endpoint interval. A trial that would
    /// cross the cold boundary supplies a bracket sign without extrapolation.
    /// Checkpoints run per layer and iteration; failed/cancelled trials publish
    /// no response. Numerical work is bounded by layers, knots and caller budget.
    pub fn evaluate(
        &self,
        temperatures_k: [f64; 2],
        config: LayerConductionConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<LayerStackResponse, DepositionError> {
        checkpoint()?;
        if config.max_iterations == 0
            || !config.absolute_heat_tolerance_w.is_finite()
            || config.absolute_heat_tolerance_w < 0.0
            || !config.relative_heat_tolerance.is_finite()
            || config.relative_heat_tolerance < 0.0
        {
            return Err(invalid("invalid stack heat tolerance or iteration budget"));
        }
        let [ta, tb] = temperatures_k;
        let n = self.curves.len();
        for curve in self.curves {
            checkpoint()?;
            curve.conductivity_w_m_k(ta)?;
            curve.conductivity_w_m_k(tb)?;
        }
        if ta == tb {
            let mut resistance = 0.0;
            let mut correction = 0.0;
            for (i, curve) in self.curves.iter().enumerate() {
                checkpoint()?;
                add(
                    &mut resistance,
                    &mut correction,
                    positive(self.thicknesses[i] / curve.conductivity_w_m_k(ta)?)?,
                )?;
                if i < n - 1 {
                    add(&mut resistance, &mut correction, self.contacts[i])?;
                }
            }
            let g = positive(self.area / resistance)?;
            checkpoint()?;
            return Ok(LayerStackResponse {
                power_a_to_b_w: 0.0,
                conductance_w_k: g,
                endpoint_derivative_w_k: [g, -g],
                max_layer_heat_residual_w: 0.0,
                max_contact_temperature_residual_k: 0.0,
                max_contact_roundoff_allowance_k: 0.0,
                iterations: 0,
                layer_count: n,
                temperatures: [[ta, tb]; MAX_CONDUCTION_LAYERS],
            });
        }
        if n == 1 {
            let power =
                finite(self.area * (self.curves[0].integral_w_m(tb, ta)? / self.thicknesses[0]))?;
            let g = positive(power / (ta - tb))?;
            let derivatives = [
                positive(
                    self.area * (self.curves[0].conductivity_w_m_k(ta)? / self.thicknesses[0]),
                )?,
                -positive(
                    self.area * (self.curves[0].conductivity_w_m_k(tb)? / self.thicknesses[0]),
                )?,
            ];
            let mut temperatures = [[0.0; 2]; MAX_CONDUCTION_LAYERS];
            temperatures[0] = [ta, tb];
            checkpoint()?;
            return Ok(LayerStackResponse {
                power_a_to_b_w: power,
                conductance_w_k: g,
                endpoint_derivative_w_k: derivatives,
                max_layer_heat_residual_w: 0.0,
                max_contact_temperature_residual_k: 0.0,
                max_contact_roundoff_allowance_k: 0.0,
                iterations: 0,
                layer_count: 1,
                temperatures,
            });
        }
        let forward = ta > tb;
        let (hot, cold, sign) = if forward {
            (ta, tb, 1.0)
        } else {
            (tb, ta, -1.0)
        };
        // Every individual layer supplies an upper flux bound. Total contact
        // resistance supplies another; no guessed effective conductivity enters.
        let mut upper = f64::INFINITY;
        for (i, curve) in self.curves.iter().enumerate() {
            checkpoint()?;
            upper = upper.min(positive(
                curve.integral_w_m(cold, hot)? / self.thicknesses[i],
            )?);
        }
        let mut contact_sum = 0.0;
        let mut correction = 0.0;
        for &r in self.contacts {
            add(&mut contact_sum, &mut correction, r)?;
        }
        if contact_sum > 0.0 {
            upper = upper.min(positive((hot - cold) / contact_sum)?);
        }
        let mut lower = 0.0;
        let mut flux = positive(0.5 * upper)?;
        for iteration in 0..config.max_iterations {
            checkpoint()?;
            let shot = self.shoot(flux, hot, cold, forward, &mut checkpoint)?;
            let Some(shot) = shot else {
                upper = flux;
                let next = lower + 0.5 * (upper - lower);
                if next <= lower || next >= upper {
                    return Err(DepositionError::Arithmetic);
                }
                flux = next;
                continue;
            };
            let power = positive(self.area * flux)?;
            let last_power = positive(self.area * shot.last_flux)?;
            let residual = finite(power - last_power)?;
            let tolerance = finite(
                config.absolute_heat_tolerance_w
                    + config.relative_heat_tolerance * power.max(last_power),
            )?;
            let last_resistance = positive(self.thicknesses[shot.last_index] / shot.last_inlet_k)?;
            let differential_resistance = positive(last_resistance + shot.inlet_flux_slope)?;
            if residual.abs() <= tolerance {
                let (heat_error, contact_error, contact_roundoff, contacts_admit) =
                    self.residuals(&shot.temperatures, sign * power, tolerance, &mut checkpoint)?;
                if heat_error <= tolerance && contacts_admit {
                    let hot_slope = positive(
                        (self.area / differential_resistance) * shot.inlet_temperature_slope,
                    )?;
                    let cold_slope = positive(
                        (self.area / differential_resistance)
                            * (self.curves[shot.last_index].conductivity_w_m_k(cold)?
                                / shot.last_inlet_k),
                    )?;
                    let derivatives = if forward {
                        [hot_slope, -cold_slope]
                    } else {
                        [cold_slope, -hot_slope]
                    };
                    checkpoint()?;
                    return Ok(LayerStackResponse {
                        power_a_to_b_w: sign * power,
                        conductance_w_k: positive(power / (hot - cold))?,
                        endpoint_derivative_w_k: derivatives,
                        max_layer_heat_residual_w: heat_error,
                        max_contact_temperature_residual_k: contact_error,
                        max_contact_roundoff_allowance_k: contact_roundoff,
                        iterations: iteration + 1,
                        layer_count: n,
                        temperatures: shot.temperatures,
                    });
                }
            }
            if residual > 0.0 {
                upper = flux;
            } else {
                lower = flux;
            }
            // Inverse residual derivative L_last/k_in /(L_last/k_in + b).
            // This bounded ratio avoids forming a potentially huge derivative.
            let newton =
                flux - (flux - shot.last_flux) * (last_resistance / differential_resistance);
            let next =
                if newton.is_finite() && newton > lower && newton < upper && iteration % 8 != 7 {
                    newton
                } else {
                    lower + 0.5 * (upper - lower)
                };
            if next == flux || next <= lower || next >= upper {
                return Err(DepositionError::Arithmetic);
            }
            flux = next;
        }
        Err(DepositionError::Budget("layer-stack conduction iterations"))
    }

    /// March from the hot boundary, inverting N-1 original layer integrals.
    /// The final layer ends at the exact cold boundary and supplies the residual.
    fn shoot(
        &self,
        flux: f64,
        hot: f64,
        cold: f64,
        forward: bool,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<Shot>, DepositionError> {
        let n = self.curves.len();
        let mut temperatures = [[0.0; 2]; MAX_CONDUCTION_LAYERS];
        let mut inlet = hot;
        let mut a = 1.0;
        let mut b = 0.0;
        for step in 0..n {
            checkpoint()?;
            let i = if forward { step } else { n - 1 - step };
            let curve = &self.curves[i];
            if inlet <= cold {
                return Ok(None);
            }
            let inlet_k = curve.conductivity_w_m_k(inlet)?;
            if step == n - 1 {
                temperatures[i] = if forward {
                    [inlet, cold]
                } else {
                    [cold, inlet]
                };
                return Ok(Some(Shot {
                    temperatures,
                    last_flux: positive(curve.integral_w_m(cold, inlet)? / self.thicknesses[i])?,
                    last_index: i,
                    inlet_temperature_slope: a,
                    inlet_flux_slope: b,
                    last_inlet_k: inlet_k,
                }));
            }
            let requested = positive(flux * self.thicknesses[i])?;
            if requested > curve.integral_w_m(cold, inlet)? {
                return Ok(None);
            }
            let outlet = curve.temperature_after_integral_w_m(inlet, -requested)?;
            if outlet < cold {
                return Ok(None);
            }
            temperatures[i] = if forward {
                [inlet, outlet]
            } else {
                [outlet, inlet]
            };
            let outlet_k = curve.conductivity_w_m_k(outlet)?;
            let ratio = positive(inlet_k / outlet_k)?;
            a = positive(a * ratio)?;
            b = positive(b * ratio + self.thicknesses[i] / outlet_k)?;
            let contact = self.contacts[if forward { i } else { i - 1 }];
            let jump = finite(flux * contact)?;
            inlet = finite(outlet - jump)?;
            if inlet < cold {
                return Ok(None);
            }
            b = positive(b + contact)?;
        }
        Err(DepositionError::Arithmetic)
    }

    fn residuals(
        &self,
        temperatures: &[[f64; 2]; MAX_CONDUCTION_LAYERS],
        power: f64,
        tolerance: f64,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(f64, f64, f64, bool), DepositionError> {
        let mut heat_error: f64 = 0.0;
        let mut contact_error: f64 = 0.0;
        let mut roundoff: f64 = 0.0;
        let mut contacts_admit = true;
        for (i, curve) in self.curves.iter().enumerate() {
            checkpoint()?;
            let [left, right] = temperatures[i];
            let actual =
                finite(self.area * (curve.integral_w_m(right, left)? / self.thicknesses[i]))?;
            heat_error = heat_error.max(finite(actual - power)?.abs());
            if i < self.contacts.len() {
                let next = temperatures[i + 1][0];
                let jump = finite((power / self.area) * self.contacts[i])?;
                let error = finite((right - next) - jump)?.abs();
                let allowance =
                    finite(8.0 * f64::EPSILON * right.abs().max(next.abs()).max(jump.abs()))?;
                let limit = finite(allowance + (tolerance / self.area) * self.contacts[i])?;
                contacts_admit &= error <= limit;
                contact_error = contact_error.max(error);
                roundoff = roundoff.max(allowance);
            }
        }
        Ok((heat_error, contact_error, roundoff, contacts_admit))
    }
}

fn add(sum: &mut f64, correction: &mut f64, term: f64) -> Result<(), DepositionError> {
    let corrected = term - *correction;
    let next = finite(*sum + corrected)?;
    *correction = (next - *sum) - corrected;
    *sum = next;
    Ok(())
}
fn invalid(reason: &'static str) -> DepositionError {
    DepositionError::InvalidInput(reason)
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
