//! General passive exchange between the two original material surfaces.

use super::*;

/// One original A-to-B surface flux and its analytic partial derivatives.
///
/// The secant represents `j = G * (Tsa - Tsb)`, with finite G >= 0. It
/// includes the complete interface transfer and excludes both material layers,
/// any surface heat inventory and independent heat sources. The two partials
/// differentiate j with respect to the ACTUAL ordered surface temperatures,
/// holding the source's other coordinates fixed. They may be signed even
/// though the heat itself is passive. `source` retains the original result
/// needed by the caller, including any explicit bulk-energy derivatives.
#[derive(Clone, Debug, PartialEq)]
pub struct SurfaceExchangeState<F> {
    pub flux_a_to_b_w_m2: f64,
    pub conductance_w_m2_k: f64,
    pub surface_derivatives_w_m2_k: [f64; 2],
    pub source: F,
}

/// One conservative material-to-material transfer after eliminating both faces.
#[derive(Clone, Debug, PartialEq)]
pub struct SurfaceLayerConductionResponse<F> {
    pub power_a_to_b_w: f64,
    /// Passive native link secant. This is not the full nonlinear derivative.
    pub conductance_w_k: f64,
    /// Actual material surface temperatures, in the original A/B order.
    pub interface_temperature_k: [f64; 2],
    /// [dQ/dTa, dQ/dTb], with the source's other coordinates held fixed.
    pub endpoint_derivative_w_k: [f64; 2],
    /// A/D multiplying a source's explicit dj/de at fixed surface temperatures.
    /// This already includes both surface responses; do not add them again.
    pub source_flux_derivative_scale_m2: f64,
    /// Hot-layer outgoing minus cold-layer incoming power, independent of order.
    pub layer_heat_residual_w: f64,
    /// Original A-to-B material power minus the original source power.
    pub surface_heat_residual_w: f64,
    /// Total callback trials, including the supported interval endpoints.
    pub iterations: usize,
    /// The accepted callback payload, retained without another source query.
    pub source: F,
}

fn admit_surface_exchange<F>(
    temperatures_k: [f64; 2],
    state: &SurfaceExchangeState<F>,
) -> Result<(), DepositionError> {
    let flux = state.flux_a_to_b_w_m2;
    let g = state.conductance_w_m2_k;
    if !flux.is_finite()
        || !g.is_finite()
        || g < 0.0
        || state
            .surface_derivatives_w_m2_k
            .iter()
            .any(|d| !d.is_finite())
    {
        return Err(DepositionError::MaterialDomain);
    }
    let difference = temperatures_k[0] - temperatures_k[1];
    let represented = finite(g * difference)?;
    if (g == 0.0 && flux != 0.0)
        || (difference == 0.0 && flux != 0.0)
        || (g > 0.0
            && difference != 0.0
            && (flux == 0.0 || represented == 0.0 || flux.signum() != difference.signum()))
        || (flux - represented).abs() > 8.0 * f64::EPSILON * flux.abs().max(represented.abs())
    {
        return Err(invalid(
            "surface exchange needs its passive heat/secant or explicit zero insulation",
        ));
    }
    Ok(())
}

impl TwoLayerConduction<'_> {
    /// Close the same exact material integrals against an original surface law.
    ///
    /// The callback receives the solved faces in ORIGINAL A/B order, including
    /// reversed heat flow. It returns their one passive transfer, both signed
    /// analytic surface partials and a payload retained at convergence. Every
    /// source query is bracketed by cancellation checks. Neither the callback
    /// nor this constitutive solve changes physical history or advances time.
    ///
    /// Each inclusive surface-temperature interval is an ORIGINAL source
    /// support, independent of either bulk temperature. Both interval ends are
    /// mapped into the common Kirchhoff flux coordinate, then intersected with
    /// the original material charts and passive temperature hull. An interval
    /// may exclude the bulk endpoints or contain just an exact supported face.
    /// No callback is evaluated outside its supplied support. Source errors
    /// propagate; bounds alone do not authorize extrapolation over data gaps.
    ///
    /// The original safeguarded scalar driver handles the supported upper and
    /// lower endpoints and every interior trial, each counted against
    /// `max_iterations`. A nonpositive differential D requests bisection at an
    /// intermediate trial. An accepted root requires D > 0 to eliminate both
    /// zero-storage faces locally stably. The bracket and local gate do not
    /// claim global uniqueness or enumerate arbitrary nonlinear source roots.
    ///
    /// Authored contact resistance must be zero: the callback owns the ENTIRE
    /// interface. Both homogeneous source-free material half cells retain their
    /// exact k(T) integrals. Their original thermal owner supplies finite-body
    /// heat capacity, storage, time integration and entropy accounting. No
    /// allocation, conductivity resampling or second scalar solver is added.
    pub fn evaluate_surface_exchange<F>(
        &self,
        temperatures_k: [f64; 2],
        surface_temperature_bounds_k: [[f64; 2]; 2],
        config: LayerConductionConfig,
        mut evaluate: impl FnMut([f64; 2]) -> Result<SurfaceExchangeState<F>, DepositionError>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<SurfaceLayerConductionResponse<F>, DepositionError> {
        checkpoint()?;
        self.admit_evaluation(temperatures_k, config)?;
        if self.contact_resistance_m2_k_w != 0.0 {
            return Err(invalid(
                "surface exchange must replace zero authored contact resistance",
            ));
        }
        for [lo, hi] in surface_temperature_bounds_k {
            if !lo.is_finite() || lo <= 0.0 || !hi.is_finite() || hi < lo {
                return Err(invalid(
                    "surface exchange needs positive ordered source temperature bounds",
                ));
            }
        }
        let [ta, tb] = temperatures_k;
        let (hot, cold, sign) = if ta > tb { (0, 1, 1.0) } else { (1, 0, -1.0) };
        let (th, tc) = (temperatures_k[hot], temperatures_k[cold]);
        self.solve_interface_with_support(
            temperatures_k,
            config,
            None,
            Some(surface_temperature_bounds_k),
            &mut checkpoint,
            |trial| {
                let mut interface = [0.0; 2];
                interface[hot] = trial.surface_temperature_k[0];
                interface[cold] = trial.surface_temperature_k[1];
                let state = evaluate(interface)?;
                admit_surface_exchange(interface, &state)?;
                let [jh, jc] = trial.layer_flux_w_m2;
                let qh = finite(self.area_m2 * jh)?;
                let qc = finite(self.area_m2 * jc)?;
                let power = if qh == 0.0 { 0.0 } else { sign * qh };
                let source_power = finite(self.area_m2 * state.flux_a_to_b_w_m2)?;
                let layer_residual = finite(qh - qc)?;
                let surface_residual = finite(power - source_power)?;
                let tolerance = finite(
                    config.absolute_heat_tolerance_w
                        + config.relative_heat_tolerance
                            * qh.abs().max(qc.abs()).max(source_power.abs()),
                )?;
                let ki = [
                    self.curves[0].conductivity_w_m_k(interface[0])?,
                    self.curves[1].conductivity_w_m_k(interface[1])?,
                ];
                let resistance = [
                    positive(self.distances_m[0] / ki[0])?,
                    positive(self.distances_m[1] / ki[1])?,
                ];
                let [a, b] = state.surface_derivatives_w_m2_k;
                let denominator =
                    finite(a.mul_add(resistance[0], (-b).mul_add(resistance[1], 1.0)))?;
                // A loose heat tolerance cannot turn a conducting interface
                // into exact insulation or accept heat in the wrong direction.
                let same_transfer = if trial.flux_w_m2 == 0.0 {
                    state.flux_a_to_b_w_m2 == 0.0
                } else {
                    qh > 0.0 && qc > 0.0 && sign * state.flux_a_to_b_w_m2 > 0.0
                };
                if same_transfer
                    && layer_residual.abs() <= tolerance
                    && surface_residual.abs() <= tolerance
                {
                    if denominator <= 0.0 {
                        return Err(invalid(
                            "surface exchange has a singular or locally unstable surface balance",
                        ));
                    }
                    let derivatives = [
                        finite(
                            (self.area_m2 * (a / denominator))
                                * (self.curves[0].conductivity_w_m_k(ta)? / ki[0]),
                        )?,
                        finite(
                            (self.area_m2 * (b / denominator))
                                * (self.curves[1].conductivity_w_m_k(tb)? / ki[1]),
                        )?,
                    ];
                    let conductance = if ta == tb {
                        // At equilibrium the original passive interface
                        // secant, not a potentially signed partial, defines
                        // the frozen decay coefficient. Scale the series form
                        // without requiring a representable reciprocal HTC.
                        let g = state.conductance_w_m2_k;
                        if g == 0.0 {
                            0.0
                        } else {
                            let bulk = positive(resistance[0] + resistance[1])?;
                            let ratio = g * bulk;
                            let density = if ratio <= 1.0 {
                                g / (1.0 + ratio)
                            } else {
                                let bulk_g = positive(bulk.recip())?;
                                bulk_g / (1.0 + bulk_g / g)
                            };
                            positive(self.area_m2 * density)?
                        }
                    } else if qh == 0.0 {
                        0.0
                    } else {
                        positive(qh / (th - tc))?
                    };
                    return Ok(LayerFluxUpdate::Accepted(SurfaceLayerConductionResponse {
                        power_a_to_b_w: power,
                        conductance_w_k: conductance,
                        interface_temperature_k: interface,
                        endpoint_derivative_w_k: derivatives,
                        source_flux_derivative_scale_m2: positive(self.area_m2 / denominator)?,
                        layer_heat_residual_w: layer_residual,
                        surface_heat_residual_w: surface_residual,
                        iterations: trial.iterations,
                        source: state.source,
                    }));
                }
                let gap = finite(sign * state.flux_a_to_b_w_m2 - trial.flux_w_m2)?;
                if denominator > 0.0 {
                    Ok(LayerFluxUpdate::Continue {
                        decreasing_residual: gap,
                        slope_magnitude: denominator,
                    })
                } else {
                    Ok(LayerFluxUpdate::Bisect {
                        decreasing_residual: gap,
                    })
                }
            },
        )
    }
}

#[cfg(test)]
#[path = "surface_exchange/tests.rs"]
mod tests;
