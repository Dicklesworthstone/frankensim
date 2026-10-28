//! Positive piecewise-linear conductivity in planar steady stacks.
//!
//! A layer obeys q L = integral(k(T), T_right..T_left). Exact trapezoidal
//! segment integrals and a scaled quadratic inverse preserve each original
//! source knot. One bounded monotone shooting solve determines common flux.
//! Every table must support the entire reservoir interval; this conservative
//! admission deliberately refuses some narrower per-layer domains that a
//! separately justified solve might support. No source extrapolation occurs.

use super::{LayeredWall, checkpoint, finite, invalid, positive};
use crate::{ConductionError, ConductivityTable};
use fs_exec::Cx;

/// Explicit work and numerical residual limits for one nonlinear wall solve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NonlinearWallOptions {
    /// Total shooting trials, including initial bracketing, in 1..=256.
    pub max_iterations: usize,
    /// Maximum final cold-reservoir temperature defect, K.
    pub temperature_tolerance_k: f64,
    /// Absolute constitutive flux residual tolerance, W/m2; nonnegative.
    pub flux_absolute_tolerance_w_per_m2: f64,
    /// Relative constitutive flux residual tolerance; nonnegative.
    pub flux_relative_tolerance: f64,
}

impl NonlinearWallOptions {
    fn validate(self) -> Result<(), ConductionError> {
        if self.max_iterations == 0 || self.max_iterations > 256 {
            return Err(invalid("nonlinear wall iterations must lie in 1..=256"));
        }
        positive(
            "nonlinear wall temperature tolerance",
            self.temperature_tolerance_k,
        )?;
        let absolute = self.flux_absolute_tolerance_w_per_m2;
        let relative = self.flux_relative_tolerance;
        if !absolute.is_finite()
            || !relative.is_finite()
            || absolute < 0.0
            || relative < 0.0
            || (absolute == 0.0 && relative == 0.0)
        {
            return Err(invalid(
                "nonlinear wall flux tolerances must be finite, nonnegative and not both zero",
            ));
        }
        Ok(())
    }

    fn flux_tolerance(self, flux: f64) -> Result<f64, ConductionError> {
        finite(
            "nonlinear wall flux tolerance",
            self.flux_relative_tolerance
                .mul_add(flux.abs(), self.flux_absolute_tolerance_w_per_m2),
        )
    }
}

/// A layer's nominal steady response at the solved boundary pair.
#[derive(Debug, Clone, PartialEq)]
pub struct NonlinearLayerTemperatures {
    /// Physical left face temperature, K.
    pub left_k: f64,
    /// Physical right face temperature, K.
    pub right_k: f64,
    /// Secant resistance (T_left - T_right)/q, m2 K/W.
    /// At equal reservoir temperatures this is L/k(T).
    pub effective_resistance_m2_k_per_w: f64,
    /// Integral of the admitted conductivity law divided by thickness, W/m2.
    pub recomputed_flux_w_per_m2: f64,
    /// Recomputed flux minus the common solved flux, W/m2.
    pub flux_residual_w_per_m2: f64,
}

/// Completed steady nonlinear heat flow with all original material sources.
/// Effective R and U belong to this boundary pair; they are not a tangent
/// conductance, a temperature-independent material property or an uncertainty.
#[derive(Debug, Clone, PartialEq)]
pub struct NonlinearLayeredWallSolution {
    /// Unchanged layer tables, geometry, contacts and boundary assumptions.
    pub wall: LayeredWall,
    /// Physical left, right reservoir temperatures, K.
    pub reservoir_temperatures_k: [f64; 2],
    /// Common planar heat-transfer area, m2.
    pub area_m2: f64,
    /// Signed left-to-right flux density, W/m2.
    pub heat_flux_w_per_m2: f64,
    /// Signed left-to-right total heat rate, W.
    pub heat_rate_w: f64,
    /// Secant inverse resistance for this boundary pair, W/(m2 K).
    pub effective_transmittance_w_per_m2_k: f64,
    /// Boundary temperature difference divided by common flux, m2 K/W.
    pub effective_resistance_m2_k_per_w: f64,
    /// Results in the original physical layer order.
    pub layer_temperatures: Vec<NonlinearLayerTemperatures>,
    /// Maximum absolute layer/contact/film constitutive flux defect, W/m2.
    pub maximum_flux_residual_w_per_m2: f64,
    /// Reconstructed cold reservoir minus its prescribed temperature, K.
    pub cold_boundary_residual_k: f64,
    /// Actual shooting trials performed; zero for equal temperatures.
    pub iterations: usize,
    /// Final signed flux enclosure from monotone shooting, W/m2.
    /// A numerical bracket, not a bound including material uncertainty.
    pub flux_bracket_w_per_m2: [f64; 2],
    /// Original caller-declared work and residual limits.
    pub options: NonlinearWallOptions,
}

struct ConductivityIntegral {
    knots: Vec<(f64, f64)>,
}

impl ConductivityIntegral {
    fn new(
        cx: &Cx<'_>,
        table: &ConductivityTable,
        low: f64,
        high: f64,
    ) -> Result<Self, ConductionError> {
        table.span().check(low)?;
        table.span().check(high)?;
        let mut knots = vec![(low, table.eval(low)?)];
        for (index, &(t, k)) in table.knots().iter().enumerate() {
            if index % 128 == 0 {
                checkpoint(cx, index)?;
            }
            if t > low && t < high {
                knots.push((t, k));
            }
        }
        if high != low {
            knots.push((high, table.eval(high)?));
        }
        for (index, pair) in knots.windows(2).enumerate() {
            if index % 128 == 0 {
                checkpoint(cx, index)?;
            }
            positive(
                "conductivity integral segment",
                segment_area(pair[1].0 - pair[0].0, pair[0].1, pair[1].1)?,
            )?;
        }
        Ok(Self { knots })
    }

    fn value(&self, segment: usize, temperature: f64) -> f64 {
        let (t0, k0) = self.knots[segment];
        let (t1, k1) = self.knots[segment + 1];
        ((temperature - t0) / (t1 - t0)).mul_add(k1 - k0, k0)
    }

    fn integral(&self, cx: &Cx<'_>, low: f64, high: f64) -> Result<f64, ConductionError> {
        if low == high {
            return Ok(0.0);
        }
        let mut sum = 0.0;
        for (index, pair) in self.knots.windows(2).enumerate() {
            if index % 128 == 0 {
                checkpoint(cx, index)?;
            }
            let a = low.max(pair[0].0);
            let b = high.min(pair[1].0);
            if a < b {
                sum = finite(
                    "conductivity integral",
                    sum + segment_area(b - a, self.value(index, a), self.value(index, b))?,
                )?;
            }
        }
        positive("nonzero conductivity integral", sum)?;
        Ok(sum)
    }

    // Returns None if this positive heat flow crosses the cold reservoir.
    // No unsupported temperature is evaluated even in rejected shooting trials.
    fn cool(&self, cx: &Cx<'_>, start: f64, amount: f64) -> Result<Option<f64>, ConductionError> {
        if amount == 0.0 {
            return Ok(Some(start));
        }
        let mut remaining = amount;
        for (index, pair) in self.knots.windows(2).enumerate().rev() {
            if index % 128 == 0 {
                checkpoint(cx, index)?;
            }
            let high = start.min(pair[1].0);
            let low = pair[0].0;
            if high <= low {
                continue;
            }
            let k_high = self.value(index, high);
            let k_low = pair[0].1;
            let width = high - low;
            let capacity = segment_area(width, k_high, k_low)?;
            if remaining > capacity {
                remaining -= capacity;
                continue;
            }
            if remaining == capacity {
                return Ok(Some(low));
            }
            // Scale conductivity and energy before taking the discriminant:
            // a*y + (b-a)*y^2/2 = z, y = cooling distance / segment width.
            // The rationalized positive root avoids subtracting nearby roots.
            let scale = k_high.max(k_low);
            let a = k_high / scale;
            let b = k_low / scale;
            let z = (remaining / capacity) * (0.5 * a + 0.5 * b);
            positive("scaled segment energy", z)?;
            let discriminant = (2.0 * (b - a)).mul_add(z, a * a);
            if !discriminant.is_finite() || discriminant < 0.0 {
                return Err(invalid(
                    "conductivity inverse discriminant is not representable",
                ));
            }
            let denominator = a + fs_math::det::sqrt(discriminant);
            positive("conductivity inverse denominator", denominator)?;
            let fraction = (2.0 * z) / denominator;
            if !fraction.is_finite() || fraction <= 0.0 || fraction > 1.0 {
                return Err(invalid("conductivity inverse left its admitted segment"));
            }
            let temperature = (-width).mul_add(fraction, high);
            if temperature < low || temperature > high {
                return Err(invalid(
                    "conductivity inverse temperature left its source interval",
                ));
            }
            return Ok(Some(temperature));
        }
        Ok(None)
    }
}

fn segment_area(width: f64, k0: f64, k1: f64) -> Result<f64, ConductionError> {
    finite(
        "linear conductivity segment integral",
        width * (0.5 * k0 + 0.5 * k1),
    )
}

struct Trial {
    faces: Vec<[f64; 2]>,
    cold_reservoir: f64,
}

// A linear film/contact temperature drop, or a proved cold-boundary crossing.
fn resistance_drop(
    temperature: f64,
    flux: f64,
    resistance: f64,
    cold: f64,
) -> Result<Option<f64>, ConductionError> {
    if resistance == 0.0 {
        return Ok(Some(temperature));
    }
    if flux > (temperature - cold) / resistance {
        return Ok(None);
    }
    let result = (-flux).mul_add(resistance, temperature);
    if !result.is_finite() || result < cold {
        return Ok(None);
    }
    Ok(Some(result))
}

impl LayeredWall {
    /// Solve positive piecewise-linear k(T) without freezing any source curve.
    ///
    /// All original source knots are integrated exactly in real arithmetic.
    /// Each source table must cover the complete positive reservoir interval,
    /// at its original fixed moisture, pressure and other query coordinates.
    /// This stronger domain admission is explicit; the existing constant solve
    /// continues to check only actual layer-face temperatures.
    ///
    /// A bounded monotone flux search starts from the hotter reservoir; it never
    /// evaluates an out-of-domain trial temperature. A candidate is returned
    /// only when both the temperature and constitutive-flux tolerances pass.
    /// Work is bounded by the declared trials and existing stack knot/layer caps.
    /// Cancellation is checked per layer and at most every 128 source segments.
    /// The result is nominal steady heat, not a transient or material-error bound.
    pub fn solve_temperature_dependent(
        &self,
        cx: &Cx<'_>,
        reservoir_temperatures_k: [f64; 2],
        area_m2: f64,
        options: NonlinearWallOptions,
    ) -> Result<NonlinearLayeredWallSolution, ConductionError> {
        checkpoint(cx, 0)?;
        options.validate()?;
        positive("wall area", area_m2)?;
        for &temperature in &reservoir_temperatures_k {
            positive("reservoir absolute temperature", temperature)?;
        }
        let forward = reservoir_temperatures_k[0] >= reservoir_temperatures_k[1];
        let cold = reservoir_temperatures_k[0].min(reservoir_temperatures_k[1]);
        let hot = reservoir_temperatures_k[0].max(reservoir_temperatures_k[1]);
        let mut laws = Vec::with_capacity(self.layers.len());
        let mut minimum_resistance =
            self.surface_resistances_m2_k_per_w[0] + self.surface_resistances_m2_k_per_w[1];
        for &resistance in &self.contact_resistances_m2_k_per_w {
            minimum_resistance =
                finite("minimum wall resistance", minimum_resistance + resistance)?;
        }
        for (index, layer) in self.layers.iter().enumerate() {
            checkpoint(cx, index)?;
            let law = ConductivityIntegral::new(cx, &layer.conductivity, cold, hot)?;
            let mut maximum = 0.0f64;
            for (knot, point) in law.knots.iter().enumerate() {
                if knot % 128 == 0 {
                    checkpoint(cx, knot)?;
                }
                maximum = maximum.max(point.1);
            }
            let resistance = layer.thickness_m / maximum;
            positive("minimum layer resistance", resistance)?;
            minimum_resistance =
                finite("minimum wall resistance", minimum_resistance + resistance)?;
            laws.push(law);
        }
        positive("minimum wall resistance", minimum_resistance)?;
        if hot == cold {
            return self.finish_nonlinear(
                cx,
                &laws,
                reservoir_temperatures_k,
                area_m2,
                0.0,
                Trial {
                    faces: vec![[hot; 2]; self.layers.len()],
                    cold_reservoir: cold,
                },
                0,
                [0.0; 2],
                options,
            );
        }

        let mut low_flux = 0.0;
        let mut high_flux = (hot - cold) / minimum_resistance;
        positive("nonlinear wall flux upper bound", high_flux)?;
        let mut bracketed = false;
        let mut trial_flux = high_flux;
        let mut last_temperature_residual = hot - cold;
        let mut last_flux_residual = None;
        for iteration in 1..=options.max_iterations {
            checkpoint(cx, iteration)?;
            match self.shoot(cx, &laws, hot, cold, forward, trial_flux)? {
                None => {
                    high_flux = trial_flux;
                    bracketed = true;
                }
                Some(trial) => {
                    last_temperature_residual = trial.cold_reservoir - cold;
                    if last_temperature_residual <= options.temperature_tolerance_k
                        && (bracketed || last_temperature_residual == 0.0)
                    {
                        let signed_flux = if forward { trial_flux } else { -trial_flux };
                        let signed_bracket = if forward {
                            [low_flux, high_flux]
                        } else {
                            [-high_flux, -low_flux]
                        };
                        let candidate = self.finish_nonlinear(
                            cx,
                            &laws,
                            reservoir_temperatures_k,
                            area_m2,
                            signed_flux,
                            trial,
                            iteration,
                            signed_bracket,
                            options,
                        )?;
                        last_flux_residual = Some(candidate.maximum_flux_residual_w_per_m2);
                        if candidate.maximum_flux_residual_w_per_m2
                            <= options.flux_tolerance(signed_flux)?
                        {
                            return Ok(candidate);
                        }
                    }
                    low_flux = trial_flux;
                    if !bracketed {
                        high_flux = finite("nonlinear wall bracketing flux", high_flux * 2.0)?;
                    }
                }
            }
            let next_flux = 0.5 * low_flux + 0.5 * high_flux;
            if next_flux <= low_flux || next_flux >= high_flux {
                break;
            }
            trial_flux = next_flux;
        }
        Err(ConductionError::Conductivity {
            what: format!(
                "bounded nonlinear wall solve did not meet tolerances: cold-boundary residual                  {last_temperature_residual:e} K (limit {:e}), constitutive flux residual                  {last_flux_residual:?} W/m2 (None: no temperature-admissible candidate); shooting budget {}",
                options.temperature_tolerance_k, options.max_iterations,
            ),
        })
    }

    fn shoot(
        &self,
        cx: &Cx<'_>,
        laws: &[ConductivityIntegral],
        hot: f64,
        cold: f64,
        forward: bool,
        flux: f64,
    ) -> Result<Option<Trial>, ConductionError> {
        let first_surface = if forward { 0 } else { 1 };
        let Some(mut temperature) = resistance_drop(
            hot,
            flux,
            self.surface_resistances_m2_k_per_w[first_surface],
            cold,
        )?
        else {
            return Ok(None);
        };
        let mut faces = vec![[0.0; 2]; self.layers.len()];
        for position in 0..self.layers.len() {
            checkpoint(cx, position)?;
            let index = if forward {
                position
            } else {
                self.layers.len() - 1 - position
            };
            let amount = finite(
                "layer conductivity integral",
                flux * self.layers[index].thickness_m,
            )?;
            positive("layer conductivity integral", amount)?;
            let Some(cooled) = laws[index].cool(cx, temperature, amount)? else {
                return Ok(None);
            };
            faces[index] = if forward {
                [temperature, cooled]
            } else {
                [cooled, temperature]
            };
            temperature = cooled;
            if position + 1 < self.layers.len() {
                let contact = if forward { index } else { index - 1 };
                let Some(next) = resistance_drop(
                    temperature,
                    flux,
                    self.contact_resistances_m2_k_per_w[contact],
                    cold,
                )?
                else {
                    return Ok(None);
                };
                temperature = next;
            }
        }
        let Some(cold_reservoir) = resistance_drop(
            temperature,
            flux,
            self.surface_resistances_m2_k_per_w[1 - first_surface],
            cold,
        )?
        else {
            return Ok(None);
        };
        Ok(Some(Trial {
            faces,
            cold_reservoir,
        }))
    }

    #[allow(clippy::too_many_arguments)]
    fn finish_nonlinear(
        &self,
        cx: &Cx<'_>,
        laws: &[ConductivityIntegral],
        reservoir_temperatures_k: [f64; 2],
        area_m2: f64,
        flux: f64,
        trial: Trial,
        iterations: usize,
        flux_bracket_w_per_m2: [f64; 2],
        options: NonlinearWallOptions,
    ) -> Result<NonlinearLayeredWallSolution, ConductionError> {
        let mut layer_temperatures = Vec::with_capacity(self.layers.len());
        let mut maximum_residual = 0.0f64;
        let mut equal_temperature_resistance =
            self.surface_resistances_m2_k_per_w[0] + self.surface_resistances_m2_k_per_w[1];
        for &contact in &self.contact_resistances_m2_k_per_w {
            equal_temperature_resistance = finite(
                "isothermal wall resistance",
                equal_temperature_resistance + contact,
            )?;
        }
        for (index, (layer, face)) in self.layers.iter().zip(&trial.faces).enumerate() {
            checkpoint(cx, index)?;
            let [left, right] = *face;
            let magnitude =
                laws[index].integral(cx, left.min(right), left.max(right))? / layer.thickness_m;
            let recomputed = if left >= right { magnitude } else { -magnitude };
            let residual = finite("nonlinear layer flux residual", recomputed - flux)?;
            maximum_residual = maximum_residual.max(residual.abs());
            let effective_resistance = if flux == 0.0 {
                layer.thickness_m / layer.conductivity.eval(left)?
            } else {
                (left - right) / flux
            };
            positive("effective layer resistance", effective_resistance)?;
            equal_temperature_resistance = finite(
                "isothermal wall resistance",
                equal_temperature_resistance + effective_resistance,
            )?;
            layer_temperatures.push(NonlinearLayerTemperatures {
                left_k: left,
                right_k: right,
                effective_resistance_m2_k_per_w: effective_resistance,
                recomputed_flux_w_per_m2: recomputed,
                flux_residual_w_per_m2: residual,
            });
        }
        let mut check_resistance = |jump: f64, resistance: f64| -> Result<(), ConductionError> {
            if resistance > 0.0 {
                let defect = finite("contact or film flux residual", jump / resistance - flux)?;
                maximum_residual = maximum_residual.max(defect.abs());
            }
            Ok(())
        };
        check_resistance(
            reservoir_temperatures_k[0] - trial.faces[0][0],
            self.surface_resistances_m2_k_per_w[0],
        )?;
        for index in 0..self.contact_resistances_m2_k_per_w.len() {
            check_resistance(
                trial.faces[index][1] - trial.faces[index + 1][0],
                self.contact_resistances_m2_k_per_w[index],
            )?;
        }
        check_resistance(
            trial
                .faces
                .last()
                .ok_or_else(|| invalid("missing nonlinear wall face"))?[1]
                - reservoir_temperatures_k[1],
            self.surface_resistances_m2_k_per_w[1],
        )?;
        let delta_t = reservoir_temperatures_k[0] - reservoir_temperatures_k[1];
        let resistance = if flux == 0.0 {
            equal_temperature_resistance
        } else {
            delta_t / flux
        };
        positive("effective wall resistance", resistance)?;
        let transmittance = 1.0 / resistance;
        positive("effective wall transmittance", transmittance)?;
        let heat = finite("nonlinear wall heat rate", flux * area_m2)?;
        if flux != 0.0 && heat == 0.0 {
            return Err(invalid("nonlinear wall heat rate underflowed"));
        }
        let cold = reservoir_temperatures_k[0].min(reservoir_temperatures_k[1]);
        let cold_residual = finite("cold boundary residual", trial.cold_reservoir - cold)?;
        checkpoint(cx, self.layers.len())?;
        Ok(NonlinearLayeredWallSolution {
            wall: self.clone(),
            reservoir_temperatures_k,
            area_m2,
            heat_flux_w_per_m2: flux,
            heat_rate_w: heat,
            effective_transmittance_w_per_m2_k: transmittance,
            effective_resistance_m2_k_per_w: resistance,
            layer_temperatures,
            maximum_flux_residual_w_per_m2: maximum_residual,
            cold_boundary_residual_k: cold_residual,
            iterations,
            flux_bracket_w_per_m2,
            options,
        })
    }
}
