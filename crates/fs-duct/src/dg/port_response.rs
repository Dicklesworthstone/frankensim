//! Same-step boundary relation from the retained exact exponential.

use super::*;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Mutex, PoisonError};

/// `dt_phi1 x` of the dense kernel for one exact state: the state part of the step integral of the
/// energy coordinates, summed as `sum_j dt_phi1[i, j] x[j]` for `j` ascending.
#[derive(Debug)]
pub(super) struct DenseIntegral {
    dt_bits: u64,
    tilde: Vec<f64>,
    pub(super) int: Vec<f64>,
}

/// A one-entry cache of the last dense product, keyed on the dt and the bits of the state. It is
/// shared by clones of a state (so a trial step on a clone reuses the response's product), is never
/// part of equality, and is separate from the Krylov kernel's cache so that the two kernels cannot
/// read each other's entry. The matrix is a deterministic function of the operator identity (which
/// `validate_state` has already checked) and the dt, so the key needs nothing more.
#[derive(Debug, Default)]
pub(super) struct IntegralCache(Mutex<Option<Arc<DenseIntegral>>>);

impl Clone for IntegralCache {
    fn clone(&self) -> Self {
        let held = self.0.lock().unwrap_or_else(PoisonError::into_inner).clone();
        Self(Mutex::new(held))
    }
}

impl PartialEq for IntegralCache {
    fn eq(&self, _: &Self) -> bool {
        true
    }
}

/// Use of the dense response cache, shared by clones of one propagator.
#[derive(Debug, Default)]
pub(super) struct DenseStats {
    calls: AtomicU64,
    products: AtomicU64,
}

/// How often the dense kernel asked for `dt_phi1 x` and how often it had to form it: `calls` counts
/// every [`DgPropagator::multiport_response`], [`DgPropagator::port_response`] and
/// [`DgPropagator::step`] on a dense kernel; `products` the `n x n` products actually summed. A
/// state evaluated repeatedly (a valve iterating on one step) costs one product.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DgDenseResponseReport {
    /// Requests for the state's product.
    pub calls: u64,
    /// Products formed (cache misses).
    pub products: u64,
}

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
    /// `dt_phi1 x` for `state`, reusing the cached product when the dt and every bit of the state
    /// match. The loop is the one the responses and the step ran on every call, in the same order,
    /// so a hit and a miss deliver the same bits.
    pub(super) fn dense_integral(
        &self,
        dt_phi1: &[f64],
        stats: &DenseStats,
        state: &DgState,
    ) -> Arc<DenseIntegral> {
        stats.calls.fetch_add(1, Ordering::Relaxed);
        {
            let held = state.integral_cache.0.lock().unwrap_or_else(PoisonError::into_inner);
            if let Some(hit) = held.as_ref() {
                if hit.dt_bits == self.dt.to_bits()
                    && hit.tilde.len() == state.tilde.len()
                    && hit.tilde.iter().zip(&state.tilde).all(|(a, b)| a.to_bits() == b.to_bits())
                {
                    return Arc::clone(hit);
                }
            }
        }
        let n = self.total_dofs;
        let mut int = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                int[i] += dt_phi1[i * n + j] * state.tilde[j];
            }
        }
        stats.products.fetch_add(1, Ordering::Relaxed);
        let product = Arc::new(DenseIntegral {
            dt_bits: self.dt.to_bits(),
            tilde: state.tilde.clone(),
            int,
        });
        *state.integral_cache.0.lock().unwrap_or_else(PoisonError::into_inner) =
            Some(Arc::clone(&product));
        product
    }

    /// Use of the dense response cache since this propagator (or the one it was cloned from) was
    /// built; `None` on the lossless and sparse kernels.
    #[must_use]
    pub fn dense_response_report(&self) -> Option<DgDenseResponseReport> {
        match &self.kernel {
            PropagatorKernel::DenseSeries { stats, .. } => Some(DgDenseResponseReport {
                calls: stats.calls.load(Ordering::Relaxed),
                products: stats.products.load(Ordering::Relaxed),
            }),
            _ => None,
        }
    }

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
                ports,
                impedance,
                stats,
                ..
            } => {
                let n = self.total_dofs;
                let count = ports.len();
                // The product of this exact state is formed once; the impedance was formed at build.
                let integral = self.dense_integral(dt_phi1, stats, state);
                let mut history = vec![0.0; count];
                for (row, (_, b)) in ports.iter().enumerate() {
                    for i in 0..n {
                        history[row] += b[i] * integral.int[i];
                    }
                    history[row] /= self.dt;
                }
                DgMultiportResponse {
                    ports: ports.iter().map(|(id, _)| *id).collect(),
                    history_pressure_pa: history,
                    impedance_pa_s_m3: impedance.clone(),
                }
            }
            PropagatorKernel::Krylov(kernel) => {
                // The same affine law from the sparse kernel: history from the cached step integral,
                // impedance from the state-independent matrix formed at build.
                let n = self.total_dofs;
                let solved = kernel.solve_cached(self.dt, state)?;
                let count = kernel.ports.len();
                let mut history = vec![0.0; count];
                for (row, (_, b)) in kernel.ports.iter().enumerate() {
                    for i in 0..n {
                        history[row] += b[i] * solved.int[i];
                    }
                    history[row] /= self.dt;
                }
                DgMultiportResponse {
                    ports: kernel.ports.iter().map(|(id, _)| *id).collect(),
                    history_pressure_pa: history,
                    impedance_pa_s_m3: kernel.impedance.clone(),
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
                ports,
                impedance,
                stats,
                ..
            } => {
                let slot = ports
                    .iter()
                    .position(|(id, _)| *id == port)
                    .ok_or_else(missing)?;
                let b = &ports[slot].1;
                let n = self.total_dofs;
                let integral = self.dense_integral(dt_phi1, stats, state);
                let mut history = 0.0;
                for i in 0..n {
                    history += b[i] * integral.int[i];
                }
                // The diagonal of the build-time matrix is this port's impedance, summed and divided
                // exactly as it was here.
                (history / self.dt, impedance[slot * ports.len() + slot])
            }
            PropagatorKernel::Krylov(kernel) => {
                let slot = kernel
                    .ports
                    .iter()
                    .position(|(id, _)| *id == port)
                    .ok_or_else(missing)?;
                let b = &kernel.ports[slot].1;
                let solved = kernel.solve_cached(self.dt, state)?;
                let mut history = 0.0;
                for i in 0..self.total_dofs {
                    history += b[i] * solved.int[i];
                }
                // The diagonal of the build-time matrix, summed and divided as it was here.
                (history / self.dt, kernel.impedance[slot * kernel.ports.len() + slot])
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
