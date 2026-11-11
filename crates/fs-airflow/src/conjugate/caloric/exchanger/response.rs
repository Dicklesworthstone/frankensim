//! Analytic inlet-enthalpy and conductance response of the same finite exchanger.
//!
//! At the admitted nominal heat Q, integrate I_h,c = integral dq/(C_h,c gap^2).
//! The original implicit UA equation gives dQ/dh_h = m_h I_h/F_Q,
//! dQ/dh_c = -m_c I_c/F_Q and dQ/dUA = 1/F_Q. Here F_Q=1/gap_end
//! for co-current and F_Q=1/gap_end+I_c for counterflow. Every derivative
//! uses original Cp at actual fluid states. No finite differences, frozen
//! average Cp, extra nonlinear solves or generic differentiation system enter.

use super::*;

/// Accuracy of the dimensionless stream-enthalpy response. Quadrature work
/// retains ExchangerConfig's panel cap independently of the primal solve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExchangerResponseConfig {
    pub absolute_enthalpy_gain_tolerance: f64,
    pub relative_enthalpy_gain_tolerance: f64,
}

impl Default for ExchangerResponseConfig {
    fn default() -> Self {
        Self {
            absolute_enthalpy_gain_tolerance: 1e-10,
            relative_enthalpy_gain_tolerance: 1e-8,
        }
    }
}

/// Fixed-mass, fixed-source response at one admitted nominal operating point.
///
/// Input coordinates are each fluid's own relative specific enthalpy [J/kg].
/// Pressure, composition and mass flow stay fixed. A source endpoint supplies
/// an admitted-side tangent, never permission to perturb outside that source.
/// Outlet changes follow dh1_out=dh1_in-dQ/m1 and dh2_out=dh2_in+dQ/m2.
///
/// Own-stream gains are retained separately from heat response so a weak
/// thermal anchor is not lost by subtracting two nearly equal numbers.
/// Allowances cover response quadrature and arithmetic at the admitted Q;
/// primal Q/UA accuracy remains separately reported by ExchangerState. They
/// are numerical allowances, not interval or material-uncertainty certificates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExchangerInletResponse {
    /// d Q_first_to_second / d h_first,in [kg/s], nonnegative.
    pub heat_from_first_inlet_enthalpy_kg_per_s: f64,
    /// d Q_first_to_second / d h_second,in [kg/s], nonpositive.
    pub heat_from_second_inlet_enthalpy_kg_per_s: f64,
    /// d Q_first_to_second / d UA [K], signed like the inlet temperature gap.
    pub heat_from_conductance_k: f64,
    /// d h_first,out / d h_first,in, with the other inlet and UA fixed.
    pub first_outlet_enthalpy_gain: f64,
    /// d h_second,out / d h_second,in, with the other inlet and UA fixed.
    pub second_outlet_enthalpy_gain: f64,
    /// Absolute allowance on both the retained own gain and its heat-transfer
    /// fraction, and on the identity they obey, loss + gain = 1: that identity
    /// holds to the loss error plus the gain error, so the allowance is that
    /// sum. Multiplication by that stream's mass flow gives kg/s allowance.
    pub first_enthalpy_gain_error_allowance: f64,
    pub second_enthalpy_gain_error_allowance: f64,
    pub conductance_response_error_allowance_k: f64,
}

impl CaloricHeatExchanger {
    /// Solve the original exchanger once, then integrate its analytic response.
    ///
    /// The ordinary solve path performs no response work. Constant Cp, zero UA
    /// and the equal-inlet-temperature limit use closed forms; source Cp at the
    /// common equilibrium temperature gives its exact first-order response even
    /// for variable laws. Nonlinear response uses the same original source
    /// breakpoints, interior-pinch checks and bounded panel geometry as the
    /// primal UA integration. Cancellation or any response gate returns no pair.
    pub fn solve_with_inlet_response(
        &self,
        cx: &Cx<'_>,
        config: ExchangerConfig,
        accuracy: ExchangerResponseConfig,
    ) -> Result<(ExchangerState, ExchangerInletResponse), ExchangerError> {
        poll(cx)?;
        accuracy.validate()?;
        let state = self.solve(cx, config)?;
        let first_cp = self
            .first
            .caloric
            .specific_heat_cp_j_kg_k(self.first.inlet_temperature_k)?;
        let second_cp = self
            .second
            .caloric
            .specific_heat_cp_j_kg_k(self.second.inlet_temperature_k)?;
        let first_rate = positive(
            self.first.mass_flow_kg_per_s * first_cp,
            "first response capacity",
        )?;
        let second_rate = positive(
            self.second.mass_flow_kg_per_s * second_cp,
            "second response capacity",
        )?;
        let driving = self.first.inlet_temperature_k - self.second.inlet_temperature_k;
        let response = if self.conductance_w_per_k == 0.0 || driving == 0.0 {
            self.constant_response(first_rate, second_rate, driving)?
        } else {
            let first_is_hot = driving > 0.0;
            let (hot, cold) = if first_is_hot {
                (&self.first, &self.second)
            } else {
                (&self.second, &self.first)
            };
            let pair = Pair {
                hot: SourcePath::new(cx, hot, cold.inlet_temperature_k, -1.0, config)?,
                cold: SourcePath::new(cx, cold, hot.inlet_temperature_k, 1.0, config)?,
                arrangement: self.arrangement,
                driving: driving.abs(),
            };
            if pair.constant_capacities().is_some() {
                self.constant_response(first_rate, second_rate, driving)?
            } else {
                let heat = state.heat_rate_first_to_second_w.abs();
                let (hot_start, cold_start) = pair.point(heat, 0.0)?;
                let (hot_end, cold_end) = pair.point(heat, heat)?;
                let start_gap = positive(
                    hot_start.temperature - cold_start.temperature,
                    "response start gap",
                )?;
                let end_gap = positive(
                    hot_end.temperature - cold_end.temperature,
                    "response end gap",
                )?;
                let moments = pair.response_moments(cx, heat, end_gap, config, accuracy)?;
                let scale = hot_start.temperature
                    + cold_start.temperature
                    + hot_end.temperature
                    + cold_end.temperature;
                let (losses, gains, errors, conductance, conductance_error) =
                    moment_response(self.arrangement, start_gap, end_gap, scale, moments)?;
                let order = if first_is_hot { [0, 1] } else { [1, 0] };
                ExchangerInletResponse {
                    heat_from_first_inlet_enthalpy_kg_per_s: self.first.mass_flow_kg_per_s
                        * losses[order[0]],
                    heat_from_second_inlet_enthalpy_kg_per_s: -self.second.mass_flow_kg_per_s
                        * losses[order[1]],
                    heat_from_conductance_k: if first_is_hot {
                        conductance
                    } else {
                        -conductance
                    },
                    first_outlet_enthalpy_gain: gains[order[0]],
                    second_outlet_enthalpy_gain: gains[order[1]],
                    first_enthalpy_gain_error_allowance: errors[order[0]],
                    second_enthalpy_gain_error_allowance: errors[order[1]],
                    conductance_response_error_allowance_k: conductance_error,
                }
            }
        };
        response.check(self, accuracy)?;
        poll(cx)?;
        Ok((state, response))
    }

    fn constant_response(
        &self,
        first_rate: f64,
        second_rate: f64,
        driving: f64,
    ) -> Result<ExchangerInletResponse, ExchangerError> {
        let (losses, gains, conductance) = if self.conductance_w_per_k == 0.0 {
            ([0.0, 0.0], [1.0, 1.0], driving)
        } else {
            let transfer = constant_heat(
                self.conductance_w_per_k,
                1.0,
                first_rate,
                second_rate,
                self.arrangement,
            )?;
            let minimum = first_rate.min(second_rate);
            let maximum = first_rate.max(second_rate);
            let ratio = minimum / maximum;
            let ntu = self.conductance_w_per_k / minimum;
            let losses = [transfer / first_rate, transfer / second_rate];
            let (gains, response) = match self.arrangement {
                FlowArrangement::CoCurrent => {
                    let attenuation = det::exp(-ntu * (1.0 + ratio));
                    let first_weight = (first_rate / maximum) / (1.0 + ratio);
                    let second_weight = (second_rate / maximum) / (1.0 + ratio);
                    (
                        [
                            first_weight + second_weight * attenuation,
                            second_weight + first_weight * attenuation,
                        ],
                        attenuation,
                    )
                }
                FlowArrangement::Counterflow => {
                    let difference = (maximum - minimum) / maximum;
                    if difference == 0.0 {
                        let retained = if ntu.is_infinite() {
                            0.0
                        } else {
                            1.0 / (1.0 + ntu)
                        };
                        ([retained, retained], retained * retained)
                    } else {
                        let exponent = ntu * difference;
                        let attenuation = det::exp(-exponent);
                        let exchanged = -det::expm1(-exponent);
                        let retained = difference / (difference + ratio * exchanged);
                        let small_gain = retained * attenuation;
                        let gains = if first_rate <= second_rate {
                            [small_gain, retained]
                        } else {
                            [retained, small_gain]
                        };
                        (gains, retained * retained * attenuation)
                    }
                }
            };
            (losses, gains, driving * response)
        };
        let errors = [0, 1].map(|i| 64.0 * f64::EPSILON * (losses[i].abs() + gains[i].abs() + 1.0));
        Ok(ExchangerInletResponse {
            heat_from_first_inlet_enthalpy_kg_per_s: self.first.mass_flow_kg_per_s * losses[0],
            heat_from_second_inlet_enthalpy_kg_per_s: -self.second.mass_flow_kg_per_s * losses[1],
            heat_from_conductance_k: conductance,
            first_outlet_enthalpy_gain: gains[0],
            second_outlet_enthalpy_gain: gains[1],
            first_enthalpy_gain_error_allowance: errors[0],
            second_enthalpy_gain_error_allowance: errors[1],
            conductance_response_error_allowance_k: 64.0
                * f64::EPSILON
                * (driving.abs() + conductance.abs()),
        })
    }
}

impl ExchangerResponseConfig {
    fn validate(self) -> Result<(), ExchangerError> {
        positive(
            self.absolute_enthalpy_gain_tolerance,
            "absolute exchanger response tolerance",
        )?;
        if self.absolute_enthalpy_gain_tolerance >= 1.0
            || !(0.0..1.0).contains(&self.relative_enthalpy_gain_tolerance)
        {
            return Err(ExchangerError::Invalid("exchanger response tolerances"));
        }
        Ok(())
    }

    fn tolerance(self, value: f64) -> f64 {
        self.absolute_enthalpy_gain_tolerance + self.relative_enthalpy_gain_tolerance * value.abs()
    }
}

impl ExchangerInletResponse {
    pub(super) fn check(
        self,
        exchanger: &CaloricHeatExchanger,
        accuracy: ExchangerResponseConfig,
    ) -> Result<(), ExchangerError> {
        let losses = [
            self.heat_from_first_inlet_enthalpy_kg_per_s / exchanger.first.mass_flow_kg_per_s,
            -self.heat_from_second_inlet_enthalpy_kg_per_s / exchanger.second.mass_flow_kg_per_s,
        ];
        let gains = [
            self.first_outlet_enthalpy_gain,
            self.second_outlet_enthalpy_gain,
        ];
        let errors = [
            self.first_enthalpy_gain_error_allowance,
            self.second_enthalpy_gain_error_allowance,
        ];
        for i in 0..2 {
            if !losses[i].is_finite()
                || !gains[i].is_finite()
                || !errors[i].is_finite()
                || errors[i] < 0.0
                || losses[i] < -errors[i]
                || gains[i] < -errors[i]
                || losses[i] > 1.0 + errors[i]
                || gains[i] > 1.0 + errors[i]
                || (losses[i] + gains[i] - 1.0).abs() > errors[i] + 32.0 * f64::EPSILON
                || errors[i]
                    > accuracy
                        .tolerance(losses[i])
                        .min(accuracy.tolerance(gains[i]))
            {
                return Err(ExchangerError::Unresolved(
                    "exchanger enthalpy-response gate",
                ));
            }
        }
        let driving = exchanger.first.inlet_temperature_k - exchanger.second.inlet_temperature_k;
        let tolerance = accuracy.absolute_enthalpy_gain_tolerance * driving.abs()
            + accuracy.relative_enthalpy_gain_tolerance * self.heat_from_conductance_k.abs();
        if !self.heat_from_conductance_k.is_finite()
            || !self.conductance_response_error_allowance_k.is_finite()
            || self.conductance_response_error_allowance_k < 0.0
            || self.conductance_response_error_allowance_k > tolerance
        {
            return Err(ExchangerError::Unresolved(
                "exchanger conductance-response gate",
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Copy)]
struct Moments {
    values: [f64; 2],
    errors: [f64; 2],
}

impl Pair<'_> {
    fn response_moments(
        &self,
        cx: &Cx<'_>,
        heat: f64,
        end_gap: f64,
        config: ExchangerConfig,
        accuracy: ExchangerResponseConfig,
    ) -> Result<Moments, ExchangerError> {
        let breaks = self.source_breaks(heat, config)?;
        let absolute = positive(
            accuracy.absolute_enthalpy_gain_tolerance / end_gap / 16.0,
            "response quadrature tolerance",
        )?;
        let relative = accuracy.relative_enthalpy_gain_tolerance / 16.0;
        let mut values = [Dd::ZERO; 2];
        let mut errors = [Dd::ZERO; 2];
        let mut visited = 0usize;
        let mut pending = Vec::new();
        for interval in breaks.windows(2) {
            pending.push((interval[0], interval[1]));
            while let Some((left, right)) = pending.pop() {
                poll(cx)?;
                visited += 1;
                if visited > config.max_quadrature_panels {
                    return Err(ExchangerError::WorkLimit);
                }
                let panel = self.response_panel(heat, left, right)?;
                let tolerance = absolute * ((right - left) / heat);
                if (0..2).all(|i| panel.errors[i] <= tolerance + relative * panel.values[i]) {
                    for i in 0..2 {
                        values[i] = values[i] + Dd::from_f64(panel.values[i]);
                        errors[i] = errors[i] + Dd::from_f64(panel.errors[i]);
                    }
                } else {
                    let middle = left + 0.5 * (right - left);
                    if middle <= left || middle >= right {
                        return Err(ExchangerError::Unresolved(
                            "exchanger response quadrature precision",
                        ));
                    }
                    pending.push((middle, right));
                    pending.push((left, middle));
                }
            }
        }
        let values = [values[0].to_f64(), values[1].to_f64()];
        let errors = [errors[0].to_f64(), errors[1].to_f64()];
        for i in 0..2 {
            positive(values[i], "exchanger response moment")?;
            nonnegative(errors[i], "exchanger response allowance")?;
        }
        Ok(Moments { values, errors })
    }

    fn response_panel(&self, heat: f64, left: f64, right: f64) -> Result<Moments, ExchangerError> {
        let geometry = self.panel_geometry(heat, left, right)?;
        let [v1, v2, v3, v4] = geometry.derivatives;
        // Normalized derivatives of gap^-2, with its common gap_min^-2
        // factored out. Positive terms bound every sign of original Cp slope.
        let gap = [
            1.0,
            2.0 * v1,
            6.0 * v1 * v1 + 2.0 * v2,
            24.0 * v1 * v1 * v1 + 18.0 * v1 * v2 + 2.0 * v3,
            120.0 * v1 * v1 * v1 * v1
                + 144.0 * v1 * v1 * v2
                + 18.0 * v2 * v2
                + 24.0 * v1 * v3
                + 2.0 * v4,
        ];
        let middle = left + 0.5 * (right - left);
        let slopes = [
            self.hot.slope(middle).abs(),
            self.cold.slope(self.cold_progress(heat, middle)).abs(),
        ];
        let mut values = [0.0; 2];
        let mut errors = [0.0; 2];
        for i in 0..2 {
            let a = geometry.alpha[i];
            // Derivatives of 1/C: (1,a,3a^2,15a^3,105a^4)/C_min.
            let capacity = [1.0, a, 3.0 * a * a, 15.0 * a * a * a, 105.0 * a * a * a * a];
            let mut fourth = 0.0;
            for k in 0..5 {
                if capacity[k] != 0.0 && gap[4 - k] != 0.0 {
                    fourth += [1.0, 4.0, 6.0, 4.0, 1.0][k] * capacity[k] * gap[4 - k];
                }
            }
            let scale =
                (geometry.width / geometry.capacity_min[i]) / geometry.minimum / geometry.minimum;
            let truncation = (scale / 2880.0) * fourth;
            let terms = [0, 1, 2].map(|j| {
                (geometry.width / geometry.capacities[i][j]) / geometry.gaps[j] / geometry.gaps[j]
            });
            let value = terms[0] / 6.0 + (2.0 / 3.0) * terms[1] + terms[2] / 6.0;
            let roundoff = 128.0
                * f64::EPSILON
                * value
                * (1.0
                    + 2.0 * geometry.temperature_scale / geometry.minimum
                    + geometry.temperature_scale * (slopes[i] / geometry.capacity_min[i]));
            if !value.is_finite() || value <= 0.0 || truncation.is_nan() || roundoff.is_nan() {
                return Err(ExchangerError::Unresolved(
                    "exchanger response panel arithmetic",
                ));
            }
            values[i] = value;
            errors[i] = truncation + roundoff;
        }
        Ok(Moments { values, errors })
    }
}

type ResponseParts = ([f64; 2], [f64; 2], [f64; 2], f64, f64);

fn moment_response(
    arrangement: FlowArrangement,
    start: f64,
    end: f64,
    temperature_scale: f64,
    moments: Moments,
) -> Result<ResponseParts, ExchangerError> {
    let gap_error = 64.0 * f64::EPSILON * temperature_scale;
    let relative = 64.0 * f64::EPSILON * (1.0 + temperature_scale / start.min(end));
    let x = moments.values[0] * end;
    let y = moments.values[1] * end;
    let ex = moments.errors[0] * end + relative * x.abs();
    let ey = moments.errors[1] * end + relative * y.abs();
    let ratio = end / start;
    let ratio_error = 2.0 * relative * ratio.abs();
    // `errors[i]` bounds loss[i] + gain[i] - 1, which the response gate enforces and which is
    // exactly zero in exact arithmetic. The loss and the gain are built from moment errors that
    // are independent of each other, so the identity holds to the loss error PLUS the gain error:
    // the larger of the two is not a bound on their sum (at the first measured refusal the
    // identity missed by 1.4125e-10 against allowances of 8.12e-11 and 8.14e-11).
    let (losses, gains, errors, conductance, conductance_error) = match arrangement {
        FlowArrangement::CoCurrent => (
            [x, y],
            [ratio + y, ratio + x],
            [ex + (ratio_error + ey), ey + (ratio_error + ex)],
            end,
            gap_error,
        ),
        FlowArrangement::Counterflow => {
            let denominator = 1.0 + y;
            let error = ey + 8.0 * f64::EPSILON * denominator;
            let lower = positive(denominator - error, "response denominator resolution")?;
            let losses = [x / denominator, y / denominator];
            // I_h-I_c=1/gap_end-1/gap_start. These complements
            // retain weak thermal anchors without subtracting a near-unit loss.
            let gains = [ratio / denominator, 1.0 / denominator];
            let loss_errors = [
                (ex + losses[0].abs() * error) / lower,
                (ey + losses[1].abs() * error) / lower,
            ];
            let gain_errors = [
                (ratio_error + gains[0].abs() * error) / lower,
                gains[1].abs() * error / lower,
            ];
            let errors = [0, 1].map(|i| {
                loss_errors[i]
                    + gain_errors[i]
                    + 16.0 * f64::EPSILON * losses[i].abs().max(gains[i].abs())
            });
            let conductance = end / denominator;
            (
                losses,
                gains,
                errors,
                conductance,
                (gap_error + conductance.abs() * error) / lower,
            )
        }
    };
    if !x.is_finite() || !y.is_finite() || !ratio.is_finite() {
        return Err(ExchangerError::Unresolved(
            "exchanger response normalization",
        ));
    }
    Ok((losses, gains, errors, conductance, conductance_error))
}
