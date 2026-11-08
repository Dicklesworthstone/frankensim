//! Same-factor second frequency derivatives of the original loaded pencil.

use super::super::{
    HarmonicResponseError, MatrixLayout, PreparedHarmonicResponse, finite, invalid,
    parameter::PowerSum, poll, scale, scale_complex,
};
use super::{
    FrequencyAction, FrequencyOrder, HarmonicFrequencyControls, HarmonicFrequencyResponse, add,
    vector_exponent,
    work::{ComplexSum, harmonic_mean_power_curvature_with_cancel, velocity_curvature},
};
use crate::quadratic_fraction;
use fs_exec::Cx;
use fs_math::c64::C64;

/// Total additional resources for a response and both frequency derivatives.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicFrequencyCurvatureControls {
    /// Numerical payload beyond the prepared owner and caller arrays.
    /// Conservative bound: 8*n*p + 1536*n + 768*p bytes.
    pub maximum_additional_payload_bytes: usize,
    /// Work score 256*(n*n+n*p+p*p)+1024*(n+p), including all three solves,
    /// physical observations and original-equation audits. Unloaded diagonal
    /// owners use 1280*n. Compact loading substitutes n*n*p for n*n and adds
    /// 256*k*k+1024*k for the existing k-dimensional bordered factor.
    pub maximum_work: u64,
}

/// Original load-port second derivatives per (rad/s)^2.
#[derive(Clone, Debug)]
pub struct HarmonicPortFrequencyCurvatureResponse {
    /// B^T v_omegaomega, with the same fixed original port map.
    pub flow_second_derivative: Vec<C64>,
    /// Z flow_omegaomega + 2 Z_omega flow_omega + Z_omegaomega flow.
    pub effort_second_derivative: Vec<C64>,
    /// Second derivative of signed Re(flow^H effort)/2.
    pub average_absorbed_power_second_derivative: f64,
}

/// Complete physical response, slope and curvature at one driving frequency.
#[derive(Clone, Debug)]
pub struct HarmonicFrequencyCurvatureResponse {
    /// Original motion, work and complete first derivatives.
    pub frequency_response: HarmonicFrequencyResponse,
    /// d^2 q/domega^2 in the original retained coordinates.
    pub displacement_second_derivative: Vec<C64>,
    /// -i(2*q_omega + omega*q_omegaomega), including the DC kinematic term.
    pub velocity_second_derivative: Vec<C64>,
    /// Second derivatives in the original power-conjugate port coordinates.
    pub port_load: Option<HarmonicPortFrequencyCurvatureResponse>,
    /// Second derivative of original source work [W per (rad/s)^2].
    pub average_input_power_second_derivative: f64,
    /// Second derivative of original symmetric viscous work [W per (rad/s)^2].
    pub average_viscous_power_second_derivative: f64,
    /// Source minus independently differentiated viscous and load-port work.
    pub power_balance_residual_second_derivative: f64,
    /// Maximum same-factor solve and published-coordinate second-equation
    /// backward residual. This is not a forward-error bound near a pole.
    pub relative_second_derivative_residual: f64,
}

impl PreparedHarmonicResponse {
    /// Evaluate the exact second frequency derivative with the existing factor.
    ///
    /// Under exp(-i*omega*t), A*q_omegaomega = F_omegaomega
    /// - A_omegaomega*q - 2*A_omega*q_omega, with
    /// -A_omegaomega = 2*M + 2*i*B*Z_omega*B^T + i*omega*B*Z_omegaomega*B^T.
    /// M,C,K and B remain fixed. Source derivatives are per rad/s; None declares
    /// that derivative zero. Loaded owners require BOTH complete finite complex
    /// impedance derivatives; explicit zero matrices declare zero derivatives.
    /// Unloaded owners require neither. Derivatives are supplied by the physical
    /// load owner, never inferred by perturbing its frequency sample.
    ///
    /// The primal, slope and curvature each receive independent binary amplitude
    /// scaling. All original mutual, gyroscopic and nonproportional terms remain;
    /// no new dense pencil/factor, frequency perturbation or line-shape assumption
    /// is introduced. Resource, cancellation and arithmetic refusals publish no
    /// partial result and leave this immutable prepared owner reusable.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_with_frequency_curvature(
        &self,
        force: &[C64],
        force_derivative: Option<&[C64]>,
        force_second_derivative: Option<&[C64]>,
        impedance_derivative: Option<&[C64]>,
        impedance_second_derivative: Option<&[C64]>,
        controls: HarmonicFrequencyCurvatureControls,
        cx: &Cx<'_>,
    ) -> Result<HarmonicFrequencyCurvatureResponse, HarmonicResponseError> {
        poll(cx)?;
        if force_second_derivative
            .is_some_and(|v| v.len() != self.n || v.iter().any(|&z| !finite(z)))
        {
            return Err(invalid(
                "one finite second force frequency derivative per coordinate required",
            ));
        }
        let ports = match (&self.port_load, impedance_second_derivative) {
            (None, None) => 0,
            (Some(load), Some(values))
                if values.len() == load.ports * load.ports && values.iter().all(|&z| finite(z)) =>
            {
                load.ports
            }
            _ => {
                return Err(invalid(
                    "loaded curvature requires a complete finite second impedance derivative; unloaded curvature requires none",
                ));
            }
        };
        admit(self, ports, controls)?;
        let first = self.solve_with_frequency_derivative(
            force,
            force_derivative,
            impedance_derivative,
            HarmonicFrequencyControls {
                maximum_additional_payload_bytes: controls.maximum_additional_payload_bytes,
                maximum_work: controls.maximum_work,
            },
            cx,
        )?;
        let (state, state_exponent) = normalize_motion(self, &first.response.displacement, cx)?;
        let (tangent, tangent_exponent) =
            normalize_motion(self, &first.displacement_derivative, cx)?;
        let state_action = self.frequency_action_order(
            &state,
            FrequencyOrder::Second {
                slope: impedance_derivative,
                curvature: impedance_second_derivative,
            },
            cx,
        )?;
        let tangent_action = self.frequency_action(&tangent, impedance_derivative, cx)?;
        // Convert the scaled-frequency actions to D^-1 times their PHYSICAL
        // force. The factor of two in -2*A_omega*q_omega stays in the exponent.
        let state_action_exponent = state_exponent + state_action.exponent - 2 * self.time_exponent;
        let tangent_action_exponent =
            tangent_exponent + tangent_action.exponent - self.time_exponent + 1;
        let mut rhs_exponent = i32::MIN;
        for i in 0..self.n {
            poll(cx)?;
            for (value, exponent) in [
                (state_action.value[i], state_action_exponent),
                (tangent_action.value[i], tangent_action_exponent),
                (
                    force_second_derivative.map_or(C64::ZERO, |v| v[i]),
                    -self.coordinate_exponents[i],
                ),
            ] {
                for component in [value.re, value.im] {
                    if component != 0.0 {
                        rhs_exponent = rhs_exponent.max(quadratic_fraction(component).1 + exponent);
                    }
                }
            }
        }
        if rhs_exponent == i32::MIN {
            rhs_exponent = 0;
        }
        let mut rhs = Vec::with_capacity(self.n);
        for i in 0..self.n {
            poll(cx)?;
            let mut value = scale_complex(
                force_second_derivative.map_or(C64::ZERO, |v| v[i]),
                -self.coordinate_exponents[i] - rhs_exponent,
            )?;
            let mut correction = C64::ZERO;
            add(
                &mut value,
                &mut correction,
                scale_complex(state_action.value[i], state_action_exponent - rhs_exponent)?,
            );
            add(
                &mut value,
                &mut correction,
                scale_complex(
                    tangent_action.value[i],
                    tangent_action_exponent - rhs_exponent,
                )?,
            );
            if !finite(value) {
                return Err(invalid("frequency curvature source is not representable"));
            }
            rhs.push(value);
        }
        let (curvature, solve_residual) = self.solve_normalized(&rhs, cx)?;
        let mut displacement_second_derivative = Vec::with_capacity(self.n);
        let mut velocity_second_derivative = Vec::with_capacity(self.n);
        for (i, &value) in curvature.iter().enumerate() {
            poll(cx)?;
            let q2 = scale_complex(
                value,
                rhs_exponent - self.coordinate_exponents[i] - 2 * self.time_exponent,
            )?;
            let v2 = velocity_curvature(
                self.angular_frequency_rad_s,
                first.displacement_derivative[i],
                q2,
            )?;
            displacement_second_derivative.push(q2);
            velocity_second_derivative.push(v2);
        }
        let source_residual = self.frequency_curvature_residual(
            &displacement_second_derivative,
            force_second_derivative,
            rhs_exponent,
            &state_action,
            state_action_exponent,
            &tangent_action,
            tangent_action_exponent,
            cx,
        )?;
        let mut input = PowerSum::default();
        let mut viscous = PowerSum::default();
        for i in 0..self.n {
            poll(cx)?;
            input.half_inner(
                force_second_derivative.map_or(C64::ZERO, |v| v[i]),
                first.response.velocity[i],
                1.0,
            )?;
            input.half_inner(
                force_derivative.map_or(C64::ZERO, |v| v[i]),
                first.velocity_derivative[i],
                2.0,
            )?;
            input.half_inner(force[i], velocity_second_derivative[i], 1.0)?;
            for j in self.layout.columns(i, self.n) {
                let c = self.viscous_entry(i, j);
                viscous.half_inner(velocity_second_derivative[i], first.response.velocity[j], c)?;
                // Keep the factor two separate from a possibly huge original C.
                viscous.half_inner(
                    first.velocity_derivative[i],
                    first.velocity_derivative[j],
                    c,
                )?;
                viscous.half_inner(
                    first.velocity_derivative[i],
                    first.velocity_derivative[j],
                    c,
                )?;
                viscous.half_inner(first.response.velocity[i], velocity_second_derivative[j], c)?;
            }
        }
        let input = input.value()?;
        let viscous = viscous.value()?;
        let port_load = match (
            &self.port_load,
            &first.response.port_load,
            &first.port_load,
            impedance_derivative,
            impedance_second_derivative,
        ) {
            (Some(load), Some(primal), Some(slope), Some(dz), Some(ddz)) => {
                let mut flow = Vec::with_capacity(ports);
                for a in 0..ports {
                    poll(cx)?;
                    let mut sum = ComplexSum::default();
                    for (i, &v2) in velocity_second_derivative.iter().enumerate() {
                        sum.product(C64::from_re(load.force_map[i * ports + a]), v2, 1.0)?;
                    }
                    flow.push(sum.value()?);
                }
                let mut effort = Vec::with_capacity(ports);
                for a in 0..ports {
                    poll(cx)?;
                    let mut sum = ComplexSum::default();
                    for b in 0..ports {
                        sum.product(load.impedance[a * ports + b], flow[b], 1.0)?;
                        sum.product(dz[a * ports + b], slope.flow_derivative[b], 2.0)?;
                        sum.product(ddz[a * ports + b], primal.flow[b], 1.0)?;
                    }
                    effort.push(sum.value()?);
                }
                let power = harmonic_mean_power_curvature_with_cancel(
                    &primal.flow,
                    &slope.flow_derivative,
                    &flow,
                    &primal.effort,
                    &slope.effort_derivative,
                    &effort,
                    || cx.checkpoint().is_err(),
                )?;
                Some(HarmonicPortFrequencyCurvatureResponse {
                    flow_second_derivative: flow,
                    effort_second_derivative: effort,
                    average_absorbed_power_second_derivative: power,
                })
            }
            _ => None,
        };
        let load_power = port_load
            .as_ref()
            .map_or(0.0, |v| v.average_absorbed_power_second_derivative);
        let mut balance = PowerSum::default();
        for value in [input, -viscous, -load_power] {
            balance.add_scalar(value)?;
        }
        let balance = balance.value()?;
        poll(cx)?;
        Ok(HarmonicFrequencyCurvatureResponse {
            frequency_response: first,
            displacement_second_derivative,
            velocity_second_derivative,
            port_load,
            average_input_power_second_derivative: input,
            average_viscous_power_second_derivative: viscous,
            power_balance_residual_second_derivative: balance,
            relative_second_derivative_residual: solve_residual.max(source_residual),
        })
    }

    // Reconstruct A*q''-F''+A''*q+2*A'*q' from PUBLISHED physical motion.
    // Keep the source and both actions separate, with their uncancelled row
    // envelopes, rather than auditing only the already summed curvature RHS.
    #[allow(clippy::too_many_arguments)]
    fn frequency_curvature_residual(
        &self,
        displacement: &[C64],
        force: Option<&[C64]>,
        rhs_exponent: i32,
        state: &FrequencyAction,
        state_exponent: i32,
        tangent: &FrequencyAction,
        tangent_exponent: i32,
        cx: &Cx<'_>,
    ) -> Result<f64, HarmonicResponseError> {
        let mut physical = Vec::with_capacity(self.n);
        for (i, &q2) in displacement.iter().enumerate() {
            poll(cx)?;
            physical.push(scale_complex(
                q2,
                self.coordinate_exponents[i] + 2 * self.time_exponent - rhs_exponent,
            )?);
        }
        let state_exponent = state_exponent - rhs_exponent;
        let tangent_exponent = tangent_exponent - rhs_exponent;
        let mut exponent = vector_exponent(&physical);
        for i in 0..self.n {
            poll(cx)?;
            for (bound, shift) in [
                (state.envelope[i], state_exponent),
                (tangent.envelope[i], tangent_exponent),
            ] {
                if bound != 0.0 {
                    exponent = exponent.max(quadratic_fraction(bound).1 + shift);
                }
            }
            if let Some(force) = force {
                for value in [force[i].re, force[i].im] {
                    if value != 0.0 {
                        exponent = exponent.max(
                            quadratic_fraction(value).1
                                - self.coordinate_exponents[i]
                                - rhs_exponent,
                        );
                    }
                }
            }
        }
        for value in &mut physical {
            *value = scale_complex(*value, -exponent)?;
        }
        let mut maximum = 0.0_f64;
        for i in 0..self.n {
            poll(cx)?;
            let source = scale_complex(
                force.map_or(C64::ZERO, |v| v[i]),
                -self.coordinate_exponents[i] - rhs_exponent - exponent,
            )?;
            let mut residual = -source;
            let mut correction = C64::ZERO;
            add(
                &mut residual,
                &mut correction,
                -scale_complex(state.value[i], state_exponent - exponent)?,
            );
            add(
                &mut residual,
                &mut correction,
                -scale_complex(tangent.value[i], tangent_exponent - exponent)?,
            );
            let mut envelope = source.abs()
                + scale(state.envelope[i], state_exponent - exponent)?
                + scale(tangent.envelope[i], tangent_exponent - exponent)?;
            for j in self.operator_columns(i) {
                let coefficient = self.operator_entry(i, j)?;
                add(&mut residual, &mut correction, coefficient * physical[j]);
                envelope += coefficient.abs() * physical[j].abs();
            }
            let relative = if envelope == 0.0 {
                0.0
            } else {
                residual.abs() / envelope
            };
            if !finite(residual) || !envelope.is_finite() || !relative.is_finite() {
                return Err(invalid(
                    "original frequency curvature row audit is not representable",
                ));
            }
            maximum = maximum.max(relative);
        }
        if maximum > self.maximum_relative_residual {
            return Err(HarmonicResponseError::Residual {
                measured: maximum,
                limit: self.maximum_relative_residual,
            });
        }
        Ok(maximum)
    }
}

fn normalize_motion(
    owner: &PreparedHarmonicResponse,
    values: &[C64],
    cx: &Cx<'_>,
) -> Result<(Vec<C64>, i32), HarmonicResponseError> {
    let mut exponent = i32::MIN;
    for (i, &value) in values.iter().enumerate() {
        poll(cx)?;
        for component in [value.re, value.im] {
            if component != 0.0 {
                exponent = exponent.max(
                    quadratic_fraction(component).1
                        + owner.coordinate_exponents[i]
                        + 2 * owner.time_exponent,
                );
            }
        }
    }
    if exponent == i32::MIN {
        exponent = 0;
    }
    let mut normalized = Vec::with_capacity(values.len());
    for (i, &value) in values.iter().enumerate() {
        poll(cx)?;
        normalized.push(scale_complex(
            value,
            owner.coordinate_exponents[i] + 2 * owner.time_exponent - exponent,
        )?);
    }
    Ok((normalized, exponent))
}

fn admit(
    owner: &PreparedHarmonicResponse,
    p: usize,
    controls: HarmonicFrequencyCurvatureControls,
) -> Result<(), HarmonicResponseError> {
    let (n, p) = (owner.n as u128, p as u128);
    let bytes = usize::try_from(8 * n * p + 1536 * n + 768 * p)
        .map_err(|_| invalid("frequency curvature payload extent overflow"))?;
    if bytes > controls.maximum_additional_payload_bytes {
        return Err(HarmonicResponseError::Budget {
            quantity: "frequency curvature additional payload bytes",
            requested: bytes,
            cap: controls.maximum_additional_payload_bytes,
        });
    }
    let k = owner
        .compact_port_factor
        .as_ref()
        .map_or(0, |f| f.dimension()) as u128;
    let structural = if k != 0 {
        n * n * p
    } else if owner.layout == MatrixLayout::Diagonal {
        n
    } else {
        n * n
    };
    let work = 256 * (structural + n * p + p * p + k * k) + 1024 * (n + p + k);
    let requested =
        u64::try_from(work).map_err(|_| invalid("frequency curvature work extent overflow"))?;
    if requested > controls.maximum_work {
        return Err(HarmonicResponseError::WorkBudget {
            requested,
            cap: controls.maximum_work,
        });
    }
    Ok(())
}
