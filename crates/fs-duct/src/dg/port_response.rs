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

/// Coupled step-mean response of every declared held-flow port:
/// mean_pressure[i] = history_pressure[i] + sum_j impedance[i, j] * flow[j].
///
/// Flows are positive into the duct at either end. Rows and columns share the
/// declaration order in ports; off-diagonal entries retain the signed response
/// between ports. Their contribution is required when several boundary flows
/// are unknown in the same step.
#[derive(Debug, Clone, PartialEq)]
pub struct DgMultiportResponse {
    /// Declared ports, in the same order as the propagator's step outputs.
    pub ports: Vec<PortId>,
    /// Mean pressure with every held flow zero [Pa], one entry per port.
    pub history_pressure_pa: Vec<f64>,
    /// Exact-ZOH mean pressure per held volume flow [Pa s/m^3], row-major.
    /// The matrix has ports.len() * ports.len() entries. Diagonal entries are
    /// nonnegative; cross-port entries may be negative. No symmetrization or
    /// passivity clipping is applied.
    pub impedance_pa_s_m3: Vec<f64>,
}

impl DgPropagator {
    /// Read the full affine mean-pressure law used by Self::step, without
    /// advancing any acoustic state, wall memory, modal anchor, or clock.
    ///
    /// A loaded outlet and a source at the inlet must be solved together through
    /// this matrix. Applying Self::port_response independently at each port
    /// would discard their same-step coupling. For a chosen flow vector, the
    /// actual step's port work is dt * sum_i flow[i] * mean_pressure[i].
    ///
    /// Both kernels read the original propagator: the dense path contracts its
    /// retained dt_phi1 and port_cols_phi2; the lossless path uses its modal
    /// anchor and stable phi functions. Diagonal entries preserve the scalar
    /// response's arithmetic. Cross terms are neither clamped nor symmetrized.
    /// A duct with no declared ports returns empty vectors.
    ///
    /// # Errors
    /// Foreign state, nonfinite response, or negative diagonal impedance.
    pub fn multiport_response(&self, state: &DgState) -> Result<DgMultiportResponse, DuctError> {
        self.validate_state(state)?;
        let response = match &self.kernel {
            PropagatorKernel::DenseSeries {
                dt_phi1,
                port_cols_phi2,
                ports,
                ..
            } => {
                let n = self.total_dofs;
                let count = ports.len();
                let mut integral = vec![0.0; n];
                for i in 0..n {
                    for j in 0..n {
                        integral[i] += dt_phi1[i * n + j] * state.tilde[j];
                    }
                }
                let mut history = vec![0.0; count];
                let mut impedance = vec![0.0; count * count];
                for (row, (_, b)) in ports.iter().enumerate() {
                    for i in 0..n {
                        history[row] += b[i] * integral[i];
                        for col in 0..count {
                            impedance[row * count + col] += b[i] * port_cols_phi2[col][i];
                        }
                    }
                    history[row] /= self.dt;
                    for col in 0..count {
                        impedance[row * count + col] /= self.dt;
                    }
                }
                DgMultiportResponse {
                    ports: ports.iter().map(|(id, _)| *id).collect(),
                    history_pressure_pa: history,
                    impedance_pa_s_m3: impedance,
                }
            }
            PropagatorKernel::SkewModal {
                sigmas,
                u_cols,
                v_cols,
                modal_ports,
            } => {
                let m = self.field_dofs;
                let count = modal_ports.len();
                let anchored = !state.modal_anchor_p.is_empty()
                    && state.modal_anchor_dt_bits == Some(self.dt.to_bits());
                let elapsed = state.unforced_steps as f64 * self.dt;
                let mut history = vec![0.0; count];
                let mut impedance = vec![0.0; count * count];
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
                    let (sinc, cosc, p2_diag, p2_off) = modal_phi_coeffs(sigmas[k] * self.dt);
                    for (row, (_, bp, bq)) in modal_ports.iter().enumerate() {
                        history[row] +=
                            bp[k] * (sinc * p - cosc * q) + bq[k] * (cosc * p + sinc * q);
                        for (col, (_, cp, cq)) in modal_ports.iter().enumerate() {
                            let value = if row == col {
                                self.dt * p2_diag * (bp[k] * bp[k] + bq[k] * bq[k])
                            } else {
                                self.dt
                                    * (bp[k] * (p2_diag * cp[k] - p2_off * cq[k])
                                        + bq[k] * (p2_off * cp[k] + p2_diag * cq[k]))
                            };
                            impedance[row * count + col] += value;
                        }
                    }
                }
                DgMultiportResponse {
                    ports: modal_ports.iter().map(|(id, _, _)| *id).collect(),
                    history_pressure_pa: history,
                    impedance_pa_s_m3: impedance,
                }
            }
        };
        let count = response.ports.len();
        if response.history_pressure_pa.iter().any(|p| !p.is_finite())
            || response.impedance_pa_s_m3.iter().any(|z| !z.is_finite())
            || (0..count).any(|i| response.impedance_pa_s_m3[i * count + i] < 0.0)
        {
            return Err(DuctError::BadParameter {
                what: "DG multiport response must be finite with nonnegative diagonal impedance",
            });
        }
        Ok(response)
    }

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
