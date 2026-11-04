//! Frequency tangents of the original gas, wall and boundary equations.
//! Source loss squares supply independent bilinear heat observations; no loss
//! is inferred by subtracting storage or by symmetrizing the rounded generator.
use super::*;

/// Original port quantities differentiated with respect to angular frequency.
/// All fields except the ID are per rad/s, in the named original physical units.
#[derive(Clone, Debug, PartialEq)]
pub struct DgHarmonicPortDerivative {
    pub id: PortId,
    pub inflow_m3_s: C64,
    pub pressure_pa: C64,
    pub mean_input_power_w: f64,
}

/// Analytic frequency derivative of every retained physical observation.
/// Every field is per rad/s. Signed energy/heat derivatives are slopes, not
/// negative physical storage or dissipation. Coefficients and geometry stay fixed.
#[derive(Clone, Debug, PartialEq)]
pub struct DgHarmonicDerivative {
    pub energy_coordinates: Vec<C64>,
    pub pressure_pa: Vec<C64>,
    pub volume_flow_m3_s: Vec<C64>,
    pub ports: Vec<DgHarmonicPortDerivative>,
    pub mean_stored_energy_j: f64,
    pub mean_wall_energy_j: f64,
    pub wall_dissipation_w: f64,
    pub molecular_dissipation_w: f64,
    pub matched_absorption_w: f64,
    pub input_power_w: f64,
    pub power_residual_w: f64,
    /// Backward residual of (sI-L)x' = BU' + i*x. The primary equations are
    /// independently checked too; neither residual certifies near-pole forward error.
    pub solve_relative_residual: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct DgHarmonicFrequencyResponse {
    pub response: DgHarmonicResponse,
    pub derivative: DgHarmonicDerivative,
}

fn admit_work<M: DgMedium>(
    view: DgContinuous<'_, M>,
    maximum: usize,
    factor: bool,
) -> Result<(), DuctError> {
    let n = view.dimension();
    let work = n.checked_mul(n).and_then(|square| {
        view.ports()
            .len()
            .checked_mul(2)
            .and_then(|ports| ports.checked_add(24 + if factor { n } else { 0 }))
            .and_then(|terms| square.checked_mul(terms))
    });
    if n == 0 || n > MAX_DG_PROPAGATOR_DOFS || work.is_none_or(|work| work > maximum) {
        return Err(invalid(
            "complete DG frequency tangent exceeds its explicit work budget",
        ));
    }
    Ok(())
}

impl<M: DgMedium> DgDuct<M> {
    /// Original continuous response and analytic derivative per rad/s, sharing
    /// ONE LU. Differentiating (sI-L)x=BU gives (sI-L)x'=BU'+i*x under exp(-i omega t).
    /// Optional port-flow derivatives may name any original port; omitted ports
    /// have zero slope. All gas, viscothermal, molecular and compliant-wall states
    /// remain, with independently differentiated physical storage and heat.
    ///
    /// Work is admitted as n²*(n+2*ports+24) before allocation. The original
    /// material frequency bands, numerical residual gates and cancellation apply.
    /// A singular isolated duct still refuses; larger coupled solves can use the
    /// borrowed derivative observer without inverting its isolated operator.
    /// No accepted state, time or material coefficient changes.
    #[allow(clippy::too_many_arguments)]
    pub fn harmonic_response_with_frequency_derivative(
        &self,
        omega_rad_s: f64,
        inflows: &[(PortId, C64)],
        inflow_derivatives: Option<&[(PortId, C64)]>,
        maximum_work: usize,
        relative_tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgHarmonicFrequencyResponse, DuctError> {
        poll(cancelled)?;
        let view = self.continuous_operator();
        view.admit_frequency(omega_rad_s)?;
        admit_tolerance(relative_tolerance)?;
        admit_work(view, maximum_work, true)?;
        let flows = view.resolve_flows(inflows)?;
        let slopes = inflow_derivatives.unwrap_or(&[]);
        let derivative_flows = view.resolve_flows(slopes)?;
        let (state, factor, scale) = view.solve_harmonic_state(omega_rad_s, &flows, cancelled)?;
        let response = view.observe_harmonic(
            omega_rad_s,
            &state,
            inflows,
            maximum_work,
            relative_tolerance,
            cancelled,
        )?;
        let mut derivative = Vec::with_capacity(state.len());
        for (i, &x) in state.iter().enumerate() {
            poll(cancelled)?;
            let rhs = view.forcing(i, &derivative_flows) + C64::new(0.0, 1.0) * x;
            let scaled = C64::new(rhs.re / scale, rhs.im / scale);
            if !finite(rhs)
                || !finite(scaled)
                || (rhs.re != 0.0 && scaled.re == 0.0)
                || (rhs.im != 0.0 && scaled.im == 0.0)
            {
                return Err(invalid("DG frequency tangent forcing is not representable"));
            }
            derivative.push(scaled);
        }
        factor
            .solve_with_cancel(&mut derivative, &mut *cancelled)
            .map_err(lu_error)?;
        let derivative = view.observe_harmonic_frequency_derivative(
            &response,
            &derivative,
            slopes,
            maximum_work,
            relative_tolerance,
            cancelled,
        )?;
        Ok(DgHarmonicFrequencyResponse {
            response,
            derivative,
        })
    }
}

impl<M: DgMedium> DgContinuous<'_, M> {
    /// Observe a complete original response/tangent supplied by a joint solver.
    /// This does not invert the isolated DG operator, so a valid loaded tangent
    /// at its bare resonance remains admissible. Both original primary equations
    /// and their differentiated equations are checked against the actual port
    /// flows; response diagnostic fields are not trusted as a residual or heat
    /// witness. The same original mass maps and source loss squares map all states.
    ///
    /// Derivatives are per rad/s at fixed coefficients. Omitted port-flow slopes
    /// are zero. Work n²*(2*ports+24) is admitted before response arrays; all bands
    /// and cancellation remain in force. Signed heat slopes do not create a new
    /// material law, spatial-accuracy certificate or transient decay estimate.
    pub fn observe_harmonic_frequency_derivative(
        &self,
        response: &DgHarmonicResponse,
        derivative: &[C64],
        inflow_derivatives: &[(PortId, C64)],
        maximum_work: usize,
        relative_tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgHarmonicDerivative, DuctError> {
        poll(cancelled)?;
        self.admit_frequency(response.omega_rad_s)?;
        admit_tolerance(relative_tolerance)?;
        admit_work(*self, maximum_work, false)?;
        let n = self.dimension();
        let state = &response.energy_coordinates;
        if state.len() != n
            || derivative.len() != n
            || response.ports.len() != self.ports().len()
            || state.iter().chain(derivative).any(|&x| !finite(x))
        {
            return Err(invalid(
                "DG frequency observation requires every finite original state and tangent",
            ));
        }
        let inflows: Vec<_> = response
            .ports
            .iter()
            .map(|port| (port.id, port.inflow_m3_s))
            .collect();
        let flows = self.resolve_flows(&inflows)?;
        let slopes = self.resolve_flows(inflow_derivatives)?;
        let norms = [state.as_slice(), derivative].map(|values| {
            values
                .iter()
                .map(|value| value.abs())
                .fold(0.0_f64, f64::max)
        });
        let s = C64::new(0.0, -response.omega_rad_s);
        let (mut norm_a, mut norm_b, mut norm_r) = (0.0_f64, [0.0_f64; 2], [0.0_f64; 2]);
        for i in 0..n {
            poll(cancelled)?;
            let forcing = [
                self.forcing(i, &flows),
                self.forcing(i, &slopes) + C64::new(0.0, 1.0) * state[i],
            ];
            let mut residual = [-forcing[0], -forcing[1]];
            let mut row_norm = 0.0;
            for j in 0..n {
                let a = C64::from_re(-self.duct.l_tilde[i * n + j])
                    + if i == j { s } else { C64::ZERO };
                residual[0] = residual[0] + a * state[j];
                residual[1] = residual[1] + a * derivative[j];
                row_norm += a.abs();
            }
            if !row_norm.is_finite() || residual.iter().chain(&forcing).any(|&value| !finite(value))
            {
                return Err(invalid("DG frequency original-equation residual overflow"));
            }
            norm_a = norm_a.max(row_norm);
            for k in 0..2 {
                norm_b[k] = norm_b[k].max(forcing[k].abs());
                norm_r[k] = norm_r[k].max(residual[k].abs());
            }
        }
        let mut relative = [0.0; 2];
        for k in 0..2 {
            let scale = norm_a * norms[k] + norm_b[k];
            relative[k] = if scale == 0.0 { 0.0 } else { norm_r[k] / scale };
            if !scale.is_finite() || !relative[k].is_finite() || relative[k] > relative_tolerance {
                return Err(invalid(
                    "DG primary or frequency-tangent equations exceed residual tolerance",
                ));
            }
        }
        let (real, imag): (Vec<_>, Vec<_>) = state.iter().map(|value| (value.re, value.im)).unzip();
        let (d_real, d_imag): (Vec<_>, Vec<_>) =
            derivative.iter().map(|value| (value.re, value.im)).unzip();
        let fields = |a: &[f64], b: &[f64]| {
            let (p, q) = self.duct.from_energy_coords(a);
            let (pi, qi) = self.duct.from_energy_coords(b);
            (
                p.into_iter()
                    .zip(pi)
                    .map(|(re, im)| C64::new(re, im))
                    .collect::<Vec<_>>(),
                q.into_iter()
                    .zip(qi)
                    .map(|(re, im)| C64::new(re, im))
                    .collect::<Vec<_>>(),
            )
        };
        let (pressure, flow) = fields(&real, &imag);
        let (pressure_pa, volume_flow_m3_s) = fields(&d_real, &d_imag);
        let wall_dissipation_w = self.duct.wall_loss_pair::<true>(&real, &d_real)
            + self.duct.wall_loss_pair::<true>(&imag, &d_imag);
        let molecular_dissipation_w = self.duct.molecular_loss_pair::<true>(&real, &d_real)
            + self.duct.molecular_loss_pair::<true>(&imag, &d_imag);
        let mut matched_absorption_w = 0.0;
        for (boundary, node, gas, area) in [
            (
                self.duct.inlet,
                0,
                &self.duct.gas,
                self.duct.inlet_area_m2(),
            ),
            (
                self.duct.outlet,
                self.duct.field_dofs - 1,
                &self.duct.outlet_gas,
                self.duct.outlet_area_m2(),
            ),
        ] {
            if boundary == DgBoundary::Matched {
                let z = gas.density_kg_m3() * gas.sound_speed_m_s() / area;
                matched_absorption_w += 0.5
                    * ((pressure[node].conj() * pressure_pa[node]).re / z
                        + z * (flow[node].conj() * volume_flow_m3_s[node]).re);
            }
        }
        let mut ports = Vec::with_capacity(self.ports().len());
        let mut input_power_w = 0.0;
        for (i, (id, column)) in self.ports().iter().enumerate() {
            poll(cancelled)?;
            let project = |values: &[C64]| {
                column
                    .iter()
                    .zip(values)
                    .fold(C64::ZERO, |sum, (&b, &x)| sum + x.scale(b))
            };
            let p = project(state);
            let dp = project(derivative);
            let power = 0.5 * (dp * flows[i].conj() + p * slopes[i].conj()).re;
            input_power_w += power;
            ports.push(DgHarmonicPortDerivative {
                id: *id,
                inflow_m3_s: slopes[i],
                pressure_pa: dp,
                mean_input_power_w: power,
            });
        }
        let storage = |start: usize| {
            0.5 * state[start..]
                .iter()
                .zip(&derivative[start..])
                .map(|(&value, &slope)| (value.conj() * slope).re)
                .sum::<f64>()
        };
        let mean_stored_energy_j = storage(0);
        let mean_wall_energy_j = storage(2 * self.duct.field_dofs);
        let power_residual_w =
            wall_dissipation_w + molecular_dissipation_w + matched_absorption_w - input_power_w;
        if pressure_pa
            .iter()
            .chain(&volume_flow_m3_s)
            .any(|&value| !finite(value))
            || ports
                .iter()
                .any(|port| !finite(port.pressure_pa) || !port.mean_input_power_w.is_finite())
            || [
                mean_stored_energy_j,
                mean_wall_energy_j,
                wall_dissipation_w,
                molecular_dissipation_w,
                matched_absorption_w,
                input_power_w,
                power_residual_w,
            ]
            .iter()
            .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "DG frequency physical field, storage or power derivative overflow",
            ));
        }
        poll(cancelled)?;
        Ok(DgHarmonicDerivative {
            energy_coordinates: derivative.to_vec(),
            pressure_pa,
            volume_flow_m3_s,
            ports,
            mean_stored_energy_j,
            mean_wall_energy_j,
            wall_dissipation_w,
            molecular_dissipation_w,
            matched_absorption_w,
            input_power_w,
            power_residual_w,
            solve_relative_residual: relative[1],
        })
    }
}

#[cfg(test)]
mod tests;
