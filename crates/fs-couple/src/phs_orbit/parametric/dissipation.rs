//! Original internal loss evaluated against the complete pumped effort.
//!
//! The pump changes BOTH the conservative field and the effort seen by the
//! constitutive law. Its state derivative is H0_xx+Q(t); no effective damping,
//! artificial external port or fitted modal coefficient replaces D(x,e).

use super::{
    NonlinearDissipation, ParametricError, ParametricHamiltonian, Result, add_matrix_action,
};
use fs_phs::harmonic::equilibrium::LinearizationError;

pub(super) struct Scratch {
    pub(super) flow: Vec<f64>,
    effort_direction: Vec<f64>,
    flow_direction: Vec<f64>,
}
impl Scratch {
    pub(super) fn new(n: usize) -> Self {
        Self {
            flow: vec![0.0; n],
            effort_direction: vec![0.0; n],
            flow_direction: vec![0.0; n],
        }
    }
}

impl<'a> ParametricHamiltonian<'a> {
    /// Attach the SAME full-state nonlinear resisting law used by the original
    /// pHS time stepper/equilibrium owner. Every evaluation uses
    /// `e=grad H0(x)+Q(t)x` and `xdot=(J-R)e-D(x,e)+Gu(t)`.
    ///
    /// The law is smooth, read-only and has complete independent state/effort
    /// actions `D_x dx + D_e de`; all constitutive history stays in original
    /// state coordinates. Its flow need not lie in the range of G. Pointwise
    /// nonnegative e dot D is admitted by the existing fs-phs owner, preserving
    /// signed roundoff and negative differential slopes. This is not a global
    /// or incremental passivity, stability or physical-validity certificate.
    ///
    /// This consumes the prepared source and adds exactly 3*n scalar scratch
    /// entries. Full native/authored storage tangents and at least n permitted
    /// D actions are required before allocation. No second nonlinear law is
    /// silently substituted; combine physical terms in the original callback.
    /// No-D sources retain their original arithmetic and scratch. Warm adapter
    /// field/tangent evaluation allocates nothing; original callbacks may have
    /// their own allocation behavior.
    ///
    /// HB, adaptive HB, variational shooting, sidebands and longer-period views
    /// all retain this source. Parameter sensitivities retain D_e*e_p when the
    /// constitutive law is fixed; explicit loss-direction queries also accept
    /// the independent fixed-state/fixed-effort D_p. Nonlinear heat objectives
    /// retain both these parameter terms and the full original state chain.
    pub fn with_dissipation(mut self, law: NonlinearDissipation<'a>) -> Result<Self> {
        let n = self.system.state_dim();
        if self.dissipation.is_some() || !self.tangent.available() || law.max_tangent_actions < n {
            return Err(ParametricError::Invalid(
                "one original nonlinear loss law, complete storage tangent and n tangent actions required",
            ));
        }
        self.scratch.get_mut().dissipation = Some(Scratch::new(n));
        self.dissipation = Some(law);
        Ok(self)
    }

    pub(super) fn evaluate_dissipation(
        &self,
        state: &[f64],
        scratch: &mut super::Scratch,
    ) -> Result<f64> {
        let Some(law) = self.dissipation else {
            return Ok(0.0);
        };
        let loss = scratch
            .dissipation
            .as_mut()
            .ok_or(ParametricError::Invalid(
                "prepared nonlinear resistance scratch required",
            ))?;
        law.evaluate_into(state, &scratch.effort, &mut loss.flow, || false)
            .map_err(ParametricError::Dissipation)
    }

    pub(super) fn field_from_effort(
        &self,
        scratch: &super::Scratch,
        out: &mut [f64],
    ) -> Result<()> {
        match &scratch.dissipation {
            Some(loss) => self.system.vector_field_with_dissipation_from_effort_into(
                &scratch.effort,
                &scratch.ports,
                &loss.flow,
                out,
            ),
            None => self
                .system
                .vector_field_from_effort_into(&scratch.effort, &scratch.ports, out),
        }
        .map_err(|_| ParametricError::SourceRefused)
    }

    /// Complete continuous tangent `(J-R)(H0_xx+Q)-D_x-D_e(H0_xx+Q)`.
    /// Nominal effort is retained while each independent direction is applied.
    pub(super) fn dissipative_state_jacobian_into(
        &self,
        time: f64,
        state: &[f64],
        out: &mut [f64],
    ) -> Result<()> {
        let n = self.system.state_dim();
        let mut scratch = self.scratch.borrow_mut();
        self.effort(time, state, &mut scratch, false)?;
        self.evaluate_dissipation(state, &mut scratch)?;
        let law = self.dissipation.ok_or(ParametricError::Invalid(
            "original nonlinear resistance required",
        ))?;
        let super::Scratch {
            q,
            effort,
            ports,
            direction,
            flow,
            dissipation,
            ..
        } = &mut *scratch;
        let loss = dissipation.as_mut().ok_or(ParametricError::Invalid(
            "prepared nonlinear resistance scratch required",
        ))?;
        ports.fill(0.0);
        direction.fill(0.0);
        for column in 0..n {
            direction[column] = 1.0;
            if !self
                .tangent
                .apply_hessian(state, direction, &mut loss.effort_direction)
            {
                return Err(ParametricError::SourceRefused);
            }
            add_matrix_action(q, direction, &mut loss.effort_direction)?;
            loss.flow_direction.fill(f64::NAN);
            if !(law.tangent)(
                state,
                effort,
                direction,
                &loss.effort_direction,
                &mut loss.flow_direction,
            ) || loss.flow_direction.iter().any(|v| !v.is_finite())
            {
                return Err(ParametricError::Dissipation(
                    LinearizationError::DissipationRefused {
                        column: Some(column),
                    },
                ));
            }
            self.system
                .vector_field_with_dissipation_from_effort_into(
                    &loss.effort_direction,
                    ports,
                    &loss.flow_direction,
                    flow,
                )
                .map_err(|_| ParametricError::SourceRefused)?;
            for row in 0..n {
                out[row * n + column] = flow[row];
            }
            direction[column] = 0.0;
        }
        Ok(())
    }

    pub(super) fn reset_failure(&self) {
        self.source_refused.set(false);
        *self.failure.borrow_mut() = None;
    }

    pub(super) fn finish<T>(&self, result: Result<T>) -> Result<T> {
        if let Some(failure) = self.failure.borrow_mut().take() {
            Err(failure)
        } else if self.source_refused.get() {
            Err(ParametricError::SourceRefused)
        } else {
            result
        }
    }
}
