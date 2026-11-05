//! Implicit work-conjugate port loading of the original spatial time step.
use super::*;
use fs_la::factor::{Lu, lu};

/// Additional connection payload and bounded dense port factorization.
#[derive(Clone, Copy, Debug)]
pub struct ViscoDynamicPortBudget {
    /// At most 128 ports per bounded, non-cancellable dense LU call.
    pub maximum_ports: usize,
    /// Connection matrices and cold/warm scratch, excluding the borrowed body
    /// and its existing sparse solver. Each cold inverse action uses the body's
    /// original iteration and true-residual limits.
    pub maximum_payload_bytes: usize,
}

/// Fixed-clock mean impedance, bound to this exact immutable solid operator.
/// Full mutual and nonsymmetric entries survive. This numerical attachment
/// neither certifies a source's passivity nor infers its storage or heat.
pub struct ViscoDynamicPortOperator<'p, 'a> {
    body: &'p ViscoDynamicOperator<'a>,
    ports: usize,
    force_map: Vec<f64>,
    mean_impedance: Vec<f64>,
    response: Vec<f64>,
    compliance: Vec<f64>,
    schur: Lu,
}

/// Accounts on the actual accepted displacement ramp, not a predicted flow.
#[derive(Debug)]
pub struct ViscoDynamicPortResponse {
    pub mean_flows: Vec<f64>,
    pub mean_efforts: Vec<f64>,
    /// Applied nodal force dotted with actual displacement increment [J].
    pub applied_work_j: f64,
    /// Signed effort-flow work INTO the external load [J]. Existing source
    /// memory may release energy; this is not material heat or inferred storage.
    pub boundary_work_j: f64,
}

/// Original full material/inertia result and its separate boundary-work account.
pub struct ViscoDynamicPortStep {
    /// external_work_j is NET nodal work: applied work minus boundary work.
    pub body: ViscoDynamicStep,
    pub port: ViscoDynamicPortResponse,
}

// Fixed-order compensation keeps signed coherent port cancellations visible.
fn sum(values: impl Iterator<Item = f64>) -> f64 {
    let (mut total, mut correction) = (0.0_f64, 0.0_f64);
    for value in values {
        let next = total + value;
        correction += if total.abs() >= value.abs() {
            (total - next) + value
        } else {
            (value - next) + total
        };
        total = next;
    }
    total + correction
}

impl<'a> ViscoDynamicOperator<'a> {
    /// Prepare e_mean=history+D*flow with flow=B^T*v_mean and body force -B*e.
    /// B is row-major free_dofs*ports in equilibrium().free_dofs order; fixed
    /// supports remain stationary. D has ports*ports entries in original order.
    ///
    /// With A=M+h²*K_mean/2, cache W=(h/2)*A^-1*B, S=B^T*W and LU(I+S*D).
    /// Warm steps use one original sparse solve and one small port solve. D is
    /// never inverted: zero/singular direct terms and dark ports are admitted.
    /// A singular coupled system refuses without regularization. Rebinding
    /// constitutive coefficients requires a new original operator/connection.
    pub fn prepare_port_impedance(
        &self,
        ports: usize,
        force_map: &[f64],
        mean_impedance: &[f64],
        budget: ViscoDynamicPortBudget,
        cx: &Cx<'_>,
    ) -> Result<ViscoDynamicPortOperator<'_, 'a>, ViscoElementError> {
        checkpoint(cx)?;
        let n = self.system.equilibrium.free_dofs.len();
        if n == 0 || ports == 0 || ports > 128 || ports > budget.maximum_ports {
            return Err(invalid("bounded nonempty spatial port dimensions required"));
        }
        let entries = n
            .checked_mul(ports)
            .ok_or_else(|| invalid("port map extent overflow"))?;
        let p = ports as u128;
        let payload = 8
            * (2 * n as u128 * p
                + 8 * p * p
                + 8 * n as u128
                + 2 * self.system.full_dofs as u128
                + 16 * p);
        if payload > budget.maximum_payload_bytes as u128
            || force_map.len() != entries
            || mean_impedance.len() != ports * ports
        {
            return Err(invalid(
                "complete bounded spatial port map and mean impedance required",
            ));
        }
        for chunk in force_map.chunks(1024).chain(mean_impedance.chunks(1024)) {
            checkpoint(cx)?;
            if chunk.iter().any(|v| !v.is_finite()) {
                return Err(invalid(
                    "finite original spatial port coefficients required",
                ));
            }
        }
        let mut response = vec![0.0; entries];
        for a in 0..ports {
            checkpoint(cx)?;
            let rhs: Vec<_> = (0..n)
                .map(|i| finite(0.5 * self.dt_s * force_map[i * ports + a]))
                .collect::<Result<_, _>>()?;
            let scale = euclidean_norm(rhs.iter().copied()).max(f64::MIN_POSITIVE);
            let (column, _) = solve_linear(
                &self.operator,
                &self.preconditioner,
                &rhs,
                finite(scale)?,
                self.config.solve,
                cx,
            )?;
            for i in 0..n {
                response[i * ports + a] = finite(column[i])?;
            }
        }
        let mut compliance = vec![0.0; ports * ports];
        for a in 0..ports {
            for b in 0..ports {
                checkpoint(cx)?;
                compliance[a * ports + b] = finite(sum(
                    (0..n).map(|i| force_map[i * ports + a] * response[i * ports + b])
                ))?;
            }
        }
        let mut matrix = vec![0.0; ports * ports];
        for a in 0..ports {
            checkpoint(cx)?;
            for b in 0..ports {
                matrix[a * ports + b] = finite(
                    (if a == b { 1.0 } else { 0.0 })
                        + sum((0..ports)
                            .map(|j| compliance[a * ports + j] * mean_impedance[j * ports + b])),
                )?;
            }
        }
        let schur =
            lu(&matrix, ports).map_err(|_| invalid("singular coupled spatial port impedance"))?;
        checkpoint(cx)?;
        Ok(ViscoDynamicPortOperator {
            body: self,
            ports,
            force_map: force_map.to_vec(),
            mean_impedance: mean_impedance.to_vec(),
            response,
            compliance,
            schur,
        })
    }
}

impl<'p, 'a> ViscoDynamicPortOperator<'p, 'a> {
    pub const fn operator(&self) -> &'p ViscoDynamicOperator<'a> {
        self.body
    }
    pub const fn port_count(&self) -> usize {
        self.ports
    }
    pub fn force_map(&self) -> &[f64] {
        &self.force_map
    }
    pub fn mean_impedance(&self) -> &[f64] {
        &self.mean_impedance
    }

    pub(super) fn admit(
        &self,
        body: &ViscoDynamicOperator<'_>,
        history: &[f64],
        cx: &Cx<'_>,
    ) -> Result<(), ViscoElementError> {
        checkpoint(cx)?;
        if !core::ptr::eq(self.body, body)
            || history.len() != self.ports
            || history.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "same original dynamic operator and finite port history required",
            ));
        }
        Ok(())
    }

    pub(super) fn correct(
        &self,
        velocity: &mut [f64],
        history: &[f64],
        cx: &Cx<'_>,
    ) -> Result<(), ViscoElementError> {
        let p = self.ports;
        let mut flow = vec![0.0; p];
        for a in 0..p {
            checkpoint(cx)?;
            flow[a] = finite(
                sum(velocity
                    .iter()
                    .enumerate()
                    .map(|(i, v)| self.force_map[i * p + a] * v))
                    - sum((0..p).map(|b| self.compliance[a * p + b] * history[b])),
            )?;
        }
        self.schur.solve(&mut flow);
        let effort = self.efforts(&flow, history)?;
        for (i, value) in velocity.iter_mut().enumerate() {
            if i % 128 == 0 {
                checkpoint(cx)?;
            }
            *value = finite(*value - sum((0..p).map(|a| self.response[i * p + a] * effort[a])))?;
        }
        Ok(())
    }

    fn efforts(&self, flow: &[f64], history: &[f64]) -> Result<Vec<f64>, ViscoElementError> {
        if flow.iter().any(|v| !v.is_finite()) {
            return Err(ViscoElementError::NonFiniteResult);
        }
        (0..self.ports)
            .map(|a| {
                finite(
                    history[a]
                        + sum((0..self.ports)
                            .map(|b| self.mean_impedance[a * self.ports + b] * flow[b])),
                )
            })
            .collect()
    }

    pub(super) fn response(
        &self,
        applied: &[f64],
        old: &[f64],
        new: &[f64],
        history: &[f64],
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, ViscoDynamicPortResponse), ViscoElementError> {
        let p = self.ports;
        let free = &self.body.system.equilibrium.free_dofs;
        let h = self.body.dt_s;
        let mut flow = vec![0.0; p];
        for a in 0..p {
            checkpoint(cx)?;
            flow[a] = finite(sum(free
                .iter()
                .enumerate()
                .map(|(i, &d)| self.force_map[i * p + a] * ((new[d] - old[d]) / h))))?;
        }
        let effort = self.efforts(&flow, history)?;
        let mut net = applied.to_vec();
        for (i, &d) in free.iter().enumerate() {
            if i % 128 == 0 {
                checkpoint(cx)?;
            }
            net[d] = finite(net[d] - sum((0..p).map(|a| self.force_map[i * p + a] * effort[a])))?;
        }
        let applied_work_j = finite(sum(applied
            .iter()
            .zip(new.iter().zip(old))
            .map(|(f, (b, a))| f * (b - a))))?;
        let boundary_work_j = finite(h * sum(flow.iter().zip(&effort).map(|(u, e)| u * e)))?;
        Ok((
            net,
            ViscoDynamicPortResponse {
                mean_flows: flow,
                mean_efforts: effort,
                applied_work_j,
                boundary_work_j,
            },
        ))
    }
}

impl<'p, 'a> ViscoDynamicState<'p, 'a> {
    /// Couple one original mean-impedance interval while retaining all material
    /// histories and the existing midpoint inertia/affine displacement ramp.
    /// The caller owns source history, source-band admission and its candidate
    /// transaction. Refusal preserves this complete structural state. Ports
    /// act on free DOFs; nonzero support velocity requires the separate existing
    /// support-motion path and is refused here.
    pub fn advance_with_port_impedance(
        &mut self,
        port: &ViscoDynamicPortOperator<'_, '_>,
        mean_force_n: &[f64],
        history_effort: &[f64],
        cx: &Cx<'_>,
    ) -> Result<ViscoDynamicPortStep, ViscoElementError> {
        let (body, port) =
            self.advance_impl(mean_force_n, None, None, Some((port, history_effort)), cx)?;
        Ok(ViscoDynamicPortStep {
            body,
            port: port.expect("accepted attached interval has its port account"),
        })
    }
}
