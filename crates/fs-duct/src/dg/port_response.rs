//! Same-step boundary relation from the retained exact exponential.

use super::*;

/// Affine step-mean response at one held-flow port, with all other port inputs
/// zero for this step: `mean_pressure = history_pressure + impedance * flow`.
/// The history contains every accepted acoustic and wall-memory state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DgPortResponse {
    /// Mean pressure with zero held flow [Pa], not end-of-step pressure.
    pub history_pressure_pa: f64,
    /// Exact-ZOH mean pressure per held volume flow [Pa s/m^3].
    /// A zero value is possible at a lossless return period; an implicit valve
    /// requiring a strictly positive load must refuse it, never clamp it.
    pub impedance_pa_s_m3: f64,
}

impl DgPropagator {
    /// Read the affine mean-pressure law used by [`Self::step`] without changing
    /// the state or clock. This permits a nonlinear valve to solve pressure and
    /// flow against the actual current duct, before either participant advances.
    ///
    /// All other ports have zero held input in this relation. Previous inputs at
    /// any port are already present in the retained state. For more than one
    /// simultaneously unknown flow a full coupled port matrix is needed; this
    /// scalar relation must not be applied independently to those ports.
    ///
    /// The dense kernel uses `b^T dt_phi1 x / dt` and `b^T port_phi2 / dt`.
    /// The lossless kernel reads the same modal anchor and stable phi functions
    /// as the original propagator. Neither path fits a reflection, integrates a
    /// second model, or substitutes instantaneous pressure for port work.
    ///
    /// # Errors
    /// Foreign state, undeclared port, or a nonfinite/negative response.
    pub fn port_response(
        &self,
        state: &DgState,
        port: PortId,
    ) -> Result<DgPortResponse, DuctError> {
        self.validate_state(state)?;
        let missing = || DuctError::BadParameter {
            what: "DG mean-port response references an undeclared PortId",
        };
        let (history, impedance) = match &self.kernel {
            PropagatorKernel::DenseSeries {
                dt_phi1,
                port_cols_phi2,
                ports,
                ..
            } => {
                let slot = ports
                    .iter()
                    .position(|(id, _)| *id == port)
                    .ok_or_else(missing)?;
                let b = &ports[slot].1;
                let n = self.total_dofs;
                let mut history = 0.0;
                let mut impedance = 0.0;
                for i in 0..n {
                    let mut integral = 0.0;
                    for j in 0..n {
                        integral += dt_phi1[i * n + j] * state.tilde[j];
                    }
                    history += b[i] * integral;
                    impedance += b[i] * port_cols_phi2[slot][i];
                }
                (history / self.dt, impedance / self.dt)
            }
            PropagatorKernel::SkewModal {
                sigmas,
                u_cols,
                v_cols,
                modal_ports,
            } => {
                let (_, bp, bq) = modal_ports
                    .iter()
                    .find(|(id, _, _)| *id == port)
                    .ok_or_else(missing)?;
                let m = self.field_dofs;
                let anchored = !state.modal_anchor_p.is_empty()
                    && state.modal_anchor_dt_bits == Some(self.dt.to_bits());
                let elapsed = state.unforced_steps as f64 * self.dt;
                let mut history = 0.0;
                let mut impedance = 0.0;
                for k in 0..m {
                    let (p, q) = if anchored {
                        let (s, c) = det::sin_cos(sigmas[k] * elapsed);
                        let (ap, aq) = (state.modal_anchor_p[k], state.modal_anchor_q[k]);
                        (c * ap - s * aq, s * ap + c * aq)
                    } else {
                        let mut p = 0.0;
                        let mut q = 0.0;
                        for i in 0..m {
                            p += u_cols[i * m + k] * state.tilde[i];
                            q += v_cols[i * m + k] * state.tilde[m + i];
                        }
                        (p, q)
                    };
                    let (sinc, cosc, p2_diag, _) = modal_phi_coeffs(sigmas[k] * self.dt);
                    history += bp[k] * (sinc * p - cosc * q) + bq[k] * (cosc * p + sinc * q);
                    impedance += self.dt * p2_diag * (bp[k] * bp[k] + bq[k] * bq[k]);
                }
                (history, impedance)
            }
        };
        if !history.is_finite() || !impedance.is_finite() || impedance < 0.0 {
            return Err(DuctError::BadParameter {
                what: "DG mean-port response must be finite with nonnegative impedance",
            });
        }
        Ok(DgPortResponse {
            history_pressure_pa: history,
            impedance_pa_s_m3: impedance,
        })
    }
}

#[cfg(test)]
mod tests;
