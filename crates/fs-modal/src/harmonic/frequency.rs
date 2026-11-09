//! Analytic frequency tangents of the complete loaded mechanical pencil.
//!
//! A q = F and A = K-w²M-iw(C+B Z(w) Bᵀ) give
//! A q_w = F_w + (2wM+iC+iB Z Bᵀ+iwB Z_w Bᵀ)q.
//! The derivative is applied to a vector, never assembled as another dense
//! operator. The original factor, binary scaling and source coordinates remain.

use super::parameter::PowerSum;
use super::{
    HarmonicResponse, HarmonicResponseError, MatrixLayout, PreparedHarmonicResponse, finite,
    invalid, poll, scale, scale_complex,
};
use crate::quadratic_fraction;
use fs_exec::Cx;
use fs_math::c64::C64;

mod action;
use action::{FrequencyAction, FrequencyOrder, FrequencyPencil, FrequencyPort};
mod curvature;
pub use curvature::{
    HarmonicFrequencyCurvatureControls, HarmonicFrequencyCurvatureResponse,
    HarmonicPortFrequencyCurvatureResponse,
};
mod work;
pub use work::{
    harmonic_mean_power_curvature_with_cancel, harmonic_mean_power_parameter_jet_with_cancel,
    harmonic_quadratic_power_curvature_with_cancel,
    harmonic_quadratic_power_derivative_with_cancel,
    harmonic_quadratic_power_parameter_jet_with_cancel,
};

/// Additional resources for one response plus its same-factor frequency tangent.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicFrequencyControls {
    /// Added numerical payload, excluding the already admitted prepared owner,
    /// caller arrays and allocator metadata. Bound: 8*n*p + 512*n + 384*p bytes.
    pub maximum_additional_payload_bytes: usize,
    /// Work score 64*(n²+n*p+p²)+256*(n+p), including both solves, original
    /// physical observations and the independently reconstructed tangent audit.
    /// An exactly diagonal unloaded modal owner uses 320*n instead.
    /// Compact loading uses 64*(n*n*p+n*p+p*p+k*k)+256*(n+p+k), including
    /// streamed original-entry audits and the k-dimensional bordered solves.
    /// Coupled LU/compact owners multiply this score by three, reserving at
    /// most two same-factor residual corrections per solve. Scalar-only
    /// owners retain the original score. The factorization was admitted
    /// when preparing the owner.
    pub maximum_work: u64,
}

/// Physical port derivatives per rad/s, in the original conjugate port units.
#[derive(Clone, Debug)]
pub struct HarmonicPortFrequencyResponse {
    /// Bᵀ v_w; the real physical port map stays fixed.
    pub flow_derivative: Vec<C64>,
    /// Z_w flow + Z flow_w, including the complete dispersive load.
    pub effort_derivative: Vec<C64>,
    /// Independent derivative of Re(flowᴴ effort)/2. It is signed boundary
    /// work per rad/s, not a fitted material loss or stored-energy estimate.
    pub average_absorbed_power_derivative: f64,
}

/// Original response and its continuous derivative with respect to physical w.
#[derive(Clone, Debug)]
pub struct HarmonicFrequencyResponse {
    /// Same original physical response as `PreparedHarmonicResponse::solve`.
    pub response: HarmonicResponse,
    /// Displacement derivative per rad/s, in the original retained coordinates.
    pub displacement_derivative: Vec<C64>,
    /// -i(q+w*q_w), retaining the kinematic derivative even at DC.
    pub velocity_derivative: Vec<C64>,
    /// Physical load-port derivatives, if the prepared owner contains a load.
    pub port_load: Option<HarmonicPortFrequencyResponse>,
    /// d Re(Fᴴ v)/2, including the authored force slope [W per rad/s].
    pub average_input_power_derivative: f64,
    /// Independently evaluated d Re(vᴴ D v)/2, `D=(C+Cᵀ)/2` [W per rad/s].
    /// Gyroscopic coupling remains in motion but contributes no heat slope.
    pub average_viscous_power_derivative: f64,
    /// Input minus original structural and load-port power slopes [W per rad/s].
    pub power_balance_residual_derivative: f64,
    /// Maximum of the same-factor tangent solve and the independently checked
    /// differentiated equation reconstructed from published physical motion.
    /// This is a backward check, not a near-pole forward-error bound.
    pub relative_derivative_residual: f64,
}

impl PreparedHarmonicResponse {
    /// Evaluate the actual frequency tangent using this frequency's cached LU.
    /// All derivatives are per rad/s under exp(-i*w*t). M,C,K and the physical
    /// port map stay fixed; `force_derivative` is the complete authored dF/dw,
    /// and None explicitly holds the source phasors fixed.
    ///
    /// A loaded owner REQUIRES `impedance_derivative` with one finite complex
    /// dZ/dw per original row-major port entry. Supply a zero matrix to declare
    /// a locally constant impedance. An unloaded owner requires None. The load
    /// owner supplies its actual causal constitutive/radiation derivative at
    /// this frequency; this method neither differentiates a sample nor infers
    /// causality, broadband passivity, material heat or load storage from Z.
    /// The slope need not be Hermitian or positive semidefinite.
    ///
    /// Full mutual and nonproportional terms, source slopes and port dispersion
    /// remain. No frequency perturbation, additional dense factor or diagonal
    /// modal approximation is used. Exact antiresonances and zero forcing have
    /// ordinary state tangents; a true singular loaded pole still refuses.
    /// The added payload/work limits are checked before allocation. Cancellation,
    /// arithmetic or source-row residual refusal publishes no partial result and
    /// leaves the reusable owner unchanged.
    pub fn solve_with_frequency_derivative(
        &self,
        force: &[C64],
        force_derivative: Option<&[C64]>,
        impedance_derivative: Option<&[C64]>,
        controls: HarmonicFrequencyControls,
        cx: &Cx<'_>,
    ) -> Result<HarmonicFrequencyResponse, HarmonicResponseError> {
        poll(cx)?;
        if force_derivative.is_some_and(|v| v.len() != self.n || v.iter().any(|&z| !finite(z))) {
            return Err(invalid(
                "one finite force frequency derivative per coordinate required",
            ));
        }
        let ports = match (&self.port_load, impedance_derivative) {
            (None, None) => 0,
            (Some(load), Some(slope))
                if slope.len() == load.ports * load.ports && slope.iter().all(|&z| finite(z)) =>
            {
                load.ports
            }
            _ => {
                return Err(invalid(
                    "loaded response requires its complete finite impedance slope; unloaded response requires none",
                ));
            }
        };
        let core = self
            .compact_port_factor
            .as_ref()
            .map_or(0, |factor| factor.dimension());
        admit(
            self.n,
            ports,
            self.layout,
            core,
            self.factor.is_some() || self.compact_port_factor.is_some(),
            controls,
        )?;
        let (solution, force_exponent, relative) = self.solve_force(force, cx)?;
        let response = self.observe(force, &solution, force_exponent, relative, cx)?;
        let state_exponent = vector_exponent(&solution);
        let mut normalized = Vec::with_capacity(self.n);
        for &z in &solution {
            normalized.push(scale_complex(z, -state_exponent)?);
        }
        let action = self.frequency_action(&normalized, impedance_derivative, cx)?;
        let action_exponent =
            force_exponent + state_exponent + action.exponent - self.time_exponent;

        // Each derivative RHS gets its own amplitude exponent. A large source
        // slope must not erase the primal force, or a high-Q tangent its drive.
        let mut rhs_exponent = i32::MIN;
        for i in 0..self.n {
            poll(cx)?;
            for value in [action.value[i].re, action.value[i].im] {
                if value != 0.0 {
                    rhs_exponent = rhs_exponent.max(quadratic_fraction(value).1 + action_exponent);
                }
            }
            if let Some(slope) = force_derivative {
                for value in [slope[i].re, slope[i].im] {
                    if value != 0.0 {
                        rhs_exponent = rhs_exponent
                            .max(quadratic_fraction(value).1 - self.coordinate_exponents[i]);
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
            let source = force_derivative.map_or(C64::ZERO, |v| v[i]);
            let value = scale_complex(source, -self.coordinate_exponents[i] - rhs_exponent)?
                + scale_complex(action.value[i], action_exponent - rhs_exponent)?;
            if !finite(value) {
                return Err(invalid("frequency tangent source is not representable"));
            }
            rhs.push(value);
        }
        let (tangent, solve_residual) = self.solve_normalized(&rhs, cx)?;
        let mut displacement_derivative = Vec::with_capacity(self.n);
        let mut velocity_derivative = Vec::with_capacity(self.n);
        for (i, &value) in tangent.iter().enumerate() {
            poll(cx)?;
            let dq = scale_complex(
                value,
                rhs_exponent - self.coordinate_exponents[i] - 2 * self.time_exponent,
            )?;
            let q = response.displacement[i];
            // FMA also avoids overflowing w*dq before a representable sum.
            let dv = C64::new(
                self.angular_frequency_rad_s.mul_add(dq.im, q.im),
                -self.angular_frequency_rad_s.mul_add(dq.re, q.re),
            );
            if !finite(dv) {
                return Err(invalid(
                    "physical velocity frequency derivative is not representable",
                ));
            }
            displacement_derivative.push(dq);
            velocity_derivative.push(dv);
        }
        let source_residual = self.frequency_residual(
            &response.displacement,
            &displacement_derivative,
            force_derivative,
            impedance_derivative,
            force_exponent,
            state_exponent,
            rhs_exponent,
            cx,
        )?;
        // Frequency and material-parameter tangents must use the same signed
        // physical-work owner. Reactive/source-phase terms can individually
        // overflow even when their total is exactly zero. Contract the original
        // damping entries before rounding, including their full nullspaces.
        let mut input = PowerSum::default();
        let mut viscous = PowerSum::default();
        for i in 0..self.n {
            poll(cx)?;
            let df = force_derivative.map_or(C64::ZERO, |v| v[i]);
            input.half_inner(df, response.velocity[i], 1.0)?;
            input.half_inner(force[i], velocity_derivative[i], 1.0)?;
            for j in self.layout.columns(i, self.n) {
                let c = self.viscous_entry(i, j);
                viscous.half_inner(velocity_derivative[i], response.velocity[j], c)?;
                viscous.half_inner(response.velocity[i], velocity_derivative[j], c)?;
            }
        }
        let input = input.value()?;
        let viscous = viscous.value()?;
        let port_load = match (&self.port_load, &response.port_load, impedance_derivative) {
            (Some(load), Some(primal), Some(slope)) => {
                let mut flow_derivative = vec![C64::ZERO; ports];
                for a in 0..ports {
                    poll(cx)?;
                    let mut correction = C64::ZERO;
                    for (i, &dv) in velocity_derivative.iter().enumerate() {
                        add(
                            &mut flow_derivative[a],
                            &mut correction,
                            dv.scale(load.force_map[i * ports + a]),
                        );
                    }
                }
                let mut effort_derivative = vec![C64::ZERO; ports];
                let mut power = PowerSum::default();
                for a in 0..ports {
                    poll(cx)?;
                    let mut correction = C64::ZERO;
                    for b in 0..ports {
                        add(
                            &mut effort_derivative[a],
                            &mut correction,
                            load.impedance[a * ports + b] * flow_derivative[b],
                        );
                        add(
                            &mut effort_derivative[a],
                            &mut correction,
                            slope[a * ports + b] * primal.flow[b],
                        );
                    }
                    power.half_inner(flow_derivative[a], primal.effort[a], 1.0)?;
                    power.half_inner(primal.flow[a], effort_derivative[a], 1.0)?;
                }
                let power = power.value()?;
                if !power.is_finite()
                    || flow_derivative
                        .iter()
                        .chain(&effort_derivative)
                        .any(|&z| !finite(z))
                {
                    return Err(invalid(
                        "physical port frequency tangent or work is not representable",
                    ));
                }
                Some(HarmonicPortFrequencyResponse {
                    flow_derivative,
                    effort_derivative,
                    average_absorbed_power_derivative: power,
                })
            }
            _ => None,
        };
        let load_power = port_load
            .as_ref()
            .map_or(0.0, |p| p.average_absorbed_power_derivative);
        let mut balance = PowerSum::default();
        for value in [input, -viscous, -load_power] {
            balance.add_scalar(value)?;
        }
        let balance = balance.value()?;
        if ![input, viscous, balance].into_iter().all(f64::is_finite) {
            return Err(invalid(
                "frequency derivative of physical power is not representable",
            ));
        }
        poll(cx)?;
        Ok(HarmonicFrequencyResponse {
            response,
            displacement_derivative,
            velocity_derivative,
            port_load,
            average_input_power_derivative: input,
            average_viscous_power_derivative: viscous,
            power_balance_residual_derivative: balance,
            relative_derivative_residual: solve_residual.max(source_residual),
        })
    }

    // Apply -A_omega to a PUBLISHED physical displacement without changing
    // its physical source or the immutable binary coordinate/time maps.
    // The caller admits complete state/load jets and total resources before
    // this internal action. Its vector and uncancelled row envelope are both
    // D^-1 times physical force, divided by 2^the returned exponent.
    pub(super) fn physical_frequency_action(
        &self,
        displacement: &[C64],
        impedance_derivative: Option<&[C64]>,
        cx: &Cx<'_>,
    ) -> Result<(Vec<C64>, Vec<f64>, i32), HarmonicResponseError> {
        let (state, state_exponent) = curvature::normalize_motion(self, displacement, cx)?;
        let action = self.frequency_action(&state, impedance_derivative, cx)?;
        let exponent = state_exponent + action.exponent - self.time_exponent;
        Ok((action.value, action.envelope, exponent))
    }

    fn frequency_action(
        &self,
        state: &[C64],
        impedance_derivative: Option<&[C64]>,
        cx: &Cx<'_>,
    ) -> Result<FrequencyAction, HarmonicResponseError> {
        self.frequency_action_order(state, FrequencyOrder::First(impedance_derivative), cx)
    }

    fn frequency_action_order(
        &self,
        state: &[C64],
        order: FrequencyOrder<'_>,
        cx: &Cx<'_>,
    ) -> Result<FrequencyAction, HarmonicResponseError> {
        let port_scaling = self
            .port_load
            .as_ref()
            .map(|load| load.scaling(&self.coordinate_exponents, cx))
            .transpose()?;
        let port = self
            .port_load
            .as_ref()
            .zip(port_scaling.as_ref())
            .map(|(load, scaling)| FrequencyPort {
                ports: load.ports,
                force_map: &load.force_map,
                impedance: &load.impedance,
                exponents: &scaling.exponents,
                minimum_time_exponent: scaling.minimum_time_exponent,
            });
        FrequencyPencil {
            n: self.n,
            layout: self.layout,
            scaled_frequency: self.scaled_frequency,
            time_exponent: self.time_exponent,
            coordinate_exponents: &self.coordinate_exponents,
            scaled_mass: &self.scaled_mass,
            scaled_damping: &self.scaled_damping,
            port,
        }
        .apply(state, order, || poll(cx))
    }

    // Audit A*q_w-F_w-(2wM+iC+iBZBᵀ+iwBZ_wBᵀ)q from the PUBLISHED physical
    // motion, not the potentially cancelling assembled tangent RHS. One common
    // amplitude scaling precedes row contractions and uncancelled row bounds.
    #[allow(clippy::too_many_arguments)]
    fn frequency_residual(
        &self,
        displacement: &[C64],
        displacement_derivative: &[C64],
        force_derivative: Option<&[C64]>,
        impedance_derivative: Option<&[C64]>,
        force_exponent: i32,
        state_exponent: i32,
        rhs_exponent: i32,
        cx: &Cx<'_>,
    ) -> Result<f64, HarmonicResponseError> {
        let mut state = Vec::with_capacity(self.n);
        let mut tangent = Vec::with_capacity(self.n);
        for i in 0..self.n {
            poll(cx)?;
            let coordinate = self.coordinate_exponents[i] + 2 * self.time_exponent;
            state.push(scale_complex(
                displacement[i],
                coordinate - force_exponent - state_exponent,
            )?);
            tangent.push(scale_complex(
                displacement_derivative[i],
                coordinate - rhs_exponent,
            )?);
        }
        let action = self.frequency_action(&state, impedance_derivative, cx)?;
        let action_exponent =
            force_exponent + state_exponent + action.exponent - self.time_exponent - rhs_exponent;
        let mut exponent = vector_exponent(&tangent);
        for i in 0..self.n {
            if action.envelope[i] != 0.0 {
                exponent = exponent.max(quadratic_fraction(action.envelope[i]).1 + action_exponent);
            }
            if let Some(slope) = force_derivative {
                for v in [slope[i].re, slope[i].im] {
                    if v != 0.0 {
                        exponent = exponent.max(
                            quadratic_fraction(v).1 - self.coordinate_exponents[i] - rhs_exponent,
                        );
                    }
                }
            }
        }
        for q in &mut tangent {
            *q = scale_complex(*q, -exponent)?;
        }
        let mut maximum = 0.0_f64;
        for i in 0..self.n {
            poll(cx)?;
            let source = scale_complex(
                force_derivative.map_or(C64::ZERO, |v| v[i]),
                -self.coordinate_exponents[i] - rhs_exponent - exponent,
            )?;
            let action_value = scale_complex(action.value[i], action_exponent - exponent)?;
            let mut residual = -source - action_value;
            let mut correction = C64::ZERO;
            let mut envelope =
                source.abs() + scale(action.envelope[i], action_exponent - exponent)?;
            for j in self.operator_columns(i) {
                let dq = tangent[j];
                let coefficient = self.operator_entry(i, j)?;
                add(&mut residual, &mut correction, coefficient * dq);
                envelope += coefficient.abs() * dq.abs();
            }
            let relative = if envelope == 0.0 {
                0.0
            } else {
                residual.abs() / envelope
            };
            if !finite(residual) || !envelope.is_finite() || !relative.is_finite() {
                return Err(invalid(
                    "original frequency tangent row audit is not representable",
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

fn admit(
    n: usize,
    p: usize,
    layout: MatrixLayout,
    core: usize,
    coupled_factor: bool,
    controls: HarmonicFrequencyControls,
) -> Result<(), HarmonicResponseError> {
    let (nn, pp) = (n as u128, p as u128);
    let bytes = 8 * nn * pp + 512 * nn + 384 * pp;
    let bytes = usize::try_from(bytes).map_err(|_| invalid("frequency payload extent overflow"))?;
    if bytes > controls.maximum_additional_payload_bytes {
        return Err(HarmonicResponseError::Budget {
            quantity: "frequency tangent additional payload bytes",
            requested: bytes,
            cap: controls.maximum_additional_payload_bytes,
        });
    }
    let kk = core as u128;
    let structural_entries = if core != 0 {
        nn * nn * pp
    } else if layout == MatrixLayout::Diagonal {
        nn
    } else {
        nn * nn
    };
    let work = 64 * (structural_entries + nn * pp + pp * pp + kk * kk) + 256 * (nn + pp + kk);
    let work = if coupled_factor { 3 * work } else { work };
    let requested = u64::try_from(work).map_err(|_| invalid("frequency work extent overflow"))?;
    if requested > controls.maximum_work {
        return Err(HarmonicResponseError::WorkBudget {
            requested,
            cap: controls.maximum_work,
        });
    }
    Ok(())
}

fn vector_exponent(values: &[C64]) -> i32 {
    values
        .iter()
        .flat_map(|v| [v.re, v.im])
        .filter(|v| *v != 0.0)
        .map(|v| quadratic_fraction(v).1)
        .max()
        .unwrap_or(0)
}

fn add(sum: &mut C64, correction: &mut C64, term: C64) {
    let adjusted = term - *correction;
    let next = *sum + adjusted;
    *correction = (next - *sum) - adjusted;
    *sum = next;
}
