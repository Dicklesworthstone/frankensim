//! Analytic frequency tangents of the ORIGINAL joint hereditary equations.

use super::{
    ViscoJointHarmonicError, ViscoJointHarmonicFactor, ViscoJointHarmonicResponse, ViscoPoleError,
    cap, divide, finite, invalid, poll, ratio,
};
use fs_exec::Cx;
use fs_math::c64::C64;
use fs_modal::harmonic::resonance::HarmonicTransferSample;
use fs_modal::harmonic::{HarmonicResponseError, harmonic_mean_power_derivative_with_cancel};

/// Additional resources for one frequency tangent of an admitted joint factor.
#[derive(Clone, Copy, Debug)]
pub struct ViscoJointHarmonicFrequencyBudget {
    /// Conservative added numerical payload: 384*n + 256*p bytes, where n is
    /// the complete original coordinate count and p the total load-port count.
    /// Excludes the separately admitted factor/primal response and allocator metadata.
    pub maximum_additional_payload_bytes: usize,
    /// Conservative added work: 1024*(n*n+n*p+p*p)+512*(n+p). Includes the
    /// tangent solve, up to three same-LU corrections, original-row audits,
    /// load observations and power derivatives. No new matrix is factored.
    pub maximum_work: u64,
}

/// Derivatives of original boundary variables per rad/s, in authored port units.
#[derive(Clone, Debug)]
pub struct ViscoJointHarmonicLoadDerivative {
    /// Every original load state and algebraic reaction; empty for a K/C/M load.
    pub state_derivative: Vec<C64>,
    /// d(B^T*u)/domega, with the original fixed attachment map.
    pub displacement_derivative: Vec<C64>,
    /// d(-i*omega*B^T*u)/domega, including the explicit kinematic term.
    pub flow_derivative: Vec<C64>,
    /// Complete original effort derivative, including direct/derivative terms.
    pub effort_derivative: Vec<C64>,
    /// Signed d Re(flow_a^H*effort_a)/2 per rad/s, in the original port order.
    /// A port exchange, including a negative one, is not a material heat channel.
    pub average_absorbed_port_power_derivative: Vec<f64>,
    /// Complete signed work contraction before rounding individual ports [W s].
    /// No storage or dissipation witness is inferred for a generic realization.
    pub average_absorbed_power_derivative: f64,
}

/// Full admitted forced state and its analytic frequency tangent.
#[derive(Clone, Debug)]
pub struct ViscoJointHarmonicFrequencyResponse {
    /// Original full-cycle material admission, state equations and independent
    /// physical energy accounts are identical to the ordinary harmonic solve.
    pub response: ViscoJointHarmonicResponse,
    /// dx/domega in the exact order of response.physical_vector: displacement,
    /// every original Mandel Maxwell strain, then load states and reactions.
    pub physical_vector_derivative: Vec<C64>,
    /// du/domega [m s], in the original free-DOF order.
    pub displacement_derivative: Vec<C64>,
    /// dv/domega = -i*(u+omega*du/domega) [m].
    pub velocity_derivative: Vec<C64>,
    pub port_load: Option<ViscoJointHarmonicLoadDerivative>,
    pub dynamic_load: Option<ViscoJointHarmonicLoadDerivative>,
    pub state_space_load: Option<ViscoJointHarmonicLoadDerivative>,
    /// Original coefficient-scaled residual of the differentiated force rows.
    pub mechanical_derivative_residual: f64,
    /// Original Maxwell-strain tangent-equation residual; zero without memory.
    pub memory_derivative_residual: f64,
    /// Original load-state tangent residual, including EACH algebraic constraint
    /// and each unqualified state-space row independently. Not a forward bound.
    pub dynamic_load_derivative_residual: f64,
}

impl ViscoJointHarmonicFrequencyResponse {
    /// Feed a fixed real displacement receiver directly to the shared actual
    /// peak/notch and half-power refinement owner. Signed interference survives.
    pub fn displacement_transfer(
        &self,
        receiver: &[f64],
    ) -> Result<HarmonicTransferSample, HarmonicResponseError> {
        HarmonicTransferSample::project_real(
            receiver,
            &self.response.displacement_m,
            &self.displacement_derivative,
        )
    }

    /// Observe physical velocity with its complete kinematic frequency slope.
    pub fn velocity_transfer(
        &self,
        receiver: &[f64],
    ) -> Result<HarmonicTransferSample, HarmonicResponseError> {
        HarmonicTransferSample::project_real(
            receiver,
            &self.response.velocity_m_s,
            &self.velocity_derivative,
        )
    }
}

impl ViscoJointHarmonicFactor<'_, '_> {
    /// Differentiate the complete original hereditary solid and loaded boundary.
    ///
    /// With s=-i*omega and P(s)=P0+s*P1+s²*P2, solve
    /// P(s)*x_omega=f_omega+i*P1*x+2*omega*P2*x using this SAME retained LU.
    /// All tensor-Prony memory, Kelvin viscosity, inertia, mutual port terms,
    /// original independent load states and algebraic reactions remain present.
    /// This requires no standalone impedance and stays valid at a bare-load pole
    /// whenever the complete joint system is nonsingular.
    ///
    /// `force_derivative_n_s=None` fixes the authored peak force; Some supplies
    /// its exact slope per rad/s in original free coordinates. Material laws,
    /// reference geometry, supports, load realization and attachment map are
    /// fixed. Their frequency dependence comes from the original evolution
    /// equations, never from a frozen effective stiffness/damping or a fit.
    ///
    /// Both original primal and tangent equations pass same-factor iterative
    /// refinement. They share the preparation's operator-evaluation ceiling;
    /// every original algebraic row is gated separately. The physical strain
    /// admission is applied to the primal cycle, not to a derivative as though
    /// it were another cycle. A new frequency needs a new factor and fresh
    /// source/load-band admission. No neighbouring-frequency validity follows.
    pub fn solve_with_frequency_derivative(
        &self,
        force_n: &[C64],
        force_derivative_n_s: Option<&[C64]>,
        budget: ViscoJointHarmonicFrequencyBudget,
        cx: &Cx<'_>,
    ) -> Result<ViscoJointHarmonicFrequencyResponse, ViscoJointHarmonicError> {
        let mut cancelled = || cx.checkpoint().is_err();
        poll(&mut cancelled)?;
        let model = self.model;
        let (n, m) = (model.dimension, model.free_dofs().len());
        if force_derivative_n_s.is_some_and(|f| f.len() != m || f.iter().any(|&x| !finite(x))) {
            return Err(
                invalid("one finite force slope per original free coordinate required").into(),
            );
        }
        let p = model.port_load.as_ref().map_or(0, |l| l.port_count())
            + model.dynamic_load.as_ref().map_or(0, |l| l.port_count())
            + model
                .state_space_load
                .as_ref()
                .map_or(0, |l| l.port_count());
        let (nn, pp) = (n as u128, p as u128);
        cap(
            "joint harmonic frequency additional bytes",
            384 * nn + 256 * pp,
            budget.maximum_additional_payload_bytes as u128,
        )?;
        cap(
            "joint harmonic frequency work",
            1024 * (nn * nn + nn * pp + pp * pp) + 512 * (nn + pp),
            budget.maximum_work as u128,
        )?;
        let mut evaluations = 0;
        let response = self.solve_with_evaluations(force_n, &mut evaluations, cx)?;
        let state = &response.physical_vector;
        let s = C64::new(0.0, -self.omega_rad_s);
        model.argument(s, &mut evaluations, &mut cancelled)?;
        // The zero trial extracts the complete original tangent RHS in the
        // already admitted row scales, without retaining another dense matrix.
        let (_, mut derivative) = self.original_frequency_residual(
            s,
            state,
            &vec![C64::ZERO; n],
            force_derivative_n_s,
            &mut cancelled,
        )?;
        self.lu.solve_with_cancel(&mut derivative, &mut cancelled)?;
        for (i, value) in derivative.iter_mut().enumerate() {
            *value = value.scale(model.column_scales[i]);
        }
        let mut residual = [f64::INFINITY; 3];
        for attempt in 0..4 {
            model.argument(s, &mut evaluations, &mut cancelled)?;
            let (r, mut correction) = self.original_frequency_residual(
                s,
                state,
                &derivative,
                force_derivative_n_s,
                &mut cancelled,
            )?;
            residual = r;
            if r.iter()
                .all(|&x| x <= model.budget.physical_residual_tolerance)
            {
                break;
            }
            if attempt == 3 {
                return Err(ViscoPoleError::PhysicalResidual {
                    mechanical: r[0],
                    memory: r[1],
                    dynamic_load: r[2],
                    tolerance: model.budget.physical_residual_tolerance,
                }
                .into());
            }
            self.lu.solve_with_cancel(&mut correction, &mut cancelled)?;
            for i in 0..n {
                derivative[i] = derivative[i] + correction[i].scale(model.column_scales[i]);
            }
        }
        let displacement_derivative = derivative[..m].to_vec();
        let ds = C64::new(0.0, -1.0);
        let velocity_derivative: Vec<_> = displacement_derivative
            .iter()
            .zip(&response.displacement_m)
            .map(|(&du, &u)| s * du + ds * u)
            .collect();
        if velocity_derivative.iter().any(|&v| !finite(v)) {
            return Err(invalid("joint harmonic velocity frequency derivative overflow").into());
        }
        let port_load = model
            .port_load
            .as_ref()
            .zip(response.port_load.as_ref())
            .map(|(load, primal)| {
                load.observe_frequency_derivative(s, &derivative[..m], primal, &mut cancelled)
            })
            .transpose()?;
        let dynamic_load = model
            .dynamic_load
            .as_ref()
            .zip(response.dynamic_load.as_ref())
            .map(|(load, primal)| {
                load.observe_frequency_derivative(model, s, &derivative, primal, &mut cancelled)
            })
            .transpose()?;
        let state_space_load = model
            .state_space_load
            .as_ref()
            .zip(response.state_space_load.as_ref())
            .map(|(load, primal)| {
                load.observe_frequency_derivative(model, s, &derivative, primal, &mut cancelled)
            })
            .transpose()?;
        poll(&mut cancelled)?;
        Ok(ViscoJointHarmonicFrequencyResponse {
            response,
            physical_vector_derivative: derivative,
            displacement_derivative,
            velocity_derivative,
            port_load,
            dynamic_load,
            state_space_load,
            mechanical_derivative_residual: residual[0],
            memory_derivative_residual: residual[1],
            dynamic_load_derivative_residual: residual[2],
        })
    }

    fn original_frequency_residual(
        &self,
        s: C64,
        state: &[C64],
        derivative: &[C64],
        force_derivative: Option<&[C64]>,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<([f64; 3], Vec<C64>), ViscoPoleError> {
        let model = self.model;
        let (n, m) = (model.dimension, model.free_dofs().len());
        if derivative.iter().any(|&x| !finite(x)) {
            return Err(invalid(
                "joint harmonic original-state frequency derivative overflow",
            ));
        }
        let mut numerators = [0.0_f64; 3];
        let mut denominators = [0.0_f64; 3];
        let algebraic_start = model.individually_checked_load_rows();
        let mut algebraic_residual = 0.0_f64;
        let mut correction = vec![C64::ZERO; n];
        for i in 0..n {
            poll(cancelled)?;
            let rhs = if i < m {
                force_derivative.map_or(C64::ZERO, |f| f[i])
            } else {
                C64::ZERO
            };
            let (mut sum, mut compensation) = (C64::ZERO, C64::ZERO);
            let mut envelope = rhs.abs();
            for j in 0..n {
                let at = i * n + j;
                let p1 = model.coefficients[1][at];
                let p2 = model.coefficients[2][at];
                let a = C64::from_re(model.coefficients[0][at]) + s.scale(p1) + (s * s).scale(p2);
                let da = C64::new(-(self.omega_rad_s * p2) * 2.0, -p1);
                for term in [a * derivative[j], da * state[j]] {
                    let term = term - compensation;
                    let next = sum + term;
                    compensation = (next - sum) - term;
                    sum = next;
                }
                // The ORIGINAL differentiated terms set the residual scale,
                // including source action before cancellation at an antiresonance.
                envelope += model.envelope(at, self.omega_rad_s) * derivative[j].abs()
                    + (p1.abs() + (self.omega_rad_s * p2).abs() * 2.0) * state[j].abs();
            }
            let defect = rhs - sum;
            correction[i] = divide(defect, model.row_scales[i]);
            if !finite(defect) || !finite(correction[i]) || !envelope.is_finite() {
                return Err(invalid(
                    "joint harmonic differentiated original equation overflow",
                ));
            }
            let group = if i < m {
                0
            } else if i < model.body_dimension() {
                1
            } else {
                2
            };
            numerators[group] = numerators[group].max(defect.abs());
            denominators[group] = denominators[group].max(envelope);
            if i >= algebraic_start {
                algebraic_residual = algebraic_residual.max(ratio(defect.abs(), envelope));
            }
        }
        let mut residual: [f64; 3] =
            core::array::from_fn(|i| ratio(numerators[i], denominators[i]));
        residual[2] = residual[2].max(algebraic_residual);
        Ok((residual, correction))
    }
}

impl ViscoJointHarmonicLoadDerivative {
    pub(in crate::linear3::visco::poles) fn from_parts(
        state_derivative: Vec<C64>,
        displacement_derivative: Vec<C64>,
        flow_derivative: Vec<C64>,
        effort_derivative: Vec<C64>,
        flow: &[C64],
        effort: &[C64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Self, ViscoPoleError> {
        let work_error = |error| match error {
            HarmonicResponseError::Cancelled => ViscoPoleError::Cancelled,
            HarmonicResponseError::Invalid(reason) => invalid(reason),
            _ => invalid("joint harmonic load frequency work is not representable"),
        };
        let mut port_power = Vec::with_capacity(flow.len());
        for i in 0..flow.len() {
            port_power.push(
                harmonic_mean_power_derivative_with_cancel(
                    &flow[i..i + 1],
                    &effort[i..i + 1],
                    &flow_derivative[i..i + 1],
                    &effort_derivative[i..i + 1],
                    &mut *cancelled,
                )
                .map_err(work_error)?,
            );
        }
        // Contract the ORIGINAL products again for the total, so rounding
        // individual port exchanges cannot erase cancelling mutual work.
        let sum = harmonic_mean_power_derivative_with_cancel(
            flow,
            effort,
            &flow_derivative,
            &effort_derivative,
            &mut *cancelled,
        )
        .map_err(work_error)?;
        if state_derivative
            .iter()
            .chain(&displacement_derivative)
            .chain(&flow_derivative)
            .chain(&effort_derivative)
            .any(|&x| !finite(x))
            || port_power.iter().any(|x| !x.is_finite())
            || !sum.is_finite()
        {
            return Err(invalid(
                "joint harmonic load frequency observation overflow",
            ));
        }
        Ok(Self {
            state_derivative,
            displacement_derivative,
            flow_derivative,
            effort_derivative,
            average_absorbed_port_power_derivative: port_power,
            average_absorbed_power_derivative: sum,
        })
    }
}
