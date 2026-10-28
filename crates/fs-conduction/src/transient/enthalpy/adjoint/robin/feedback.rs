//! Typed attachment of consumer-owned convection laws to the accepted step.

use super::super::EnthalpyFeedback;
use super::{
    EnthalpyAdjointError, EnthalpyStepLinearization, RobinPort, RobinStateFeedback,
    bind_boundary_ports, finite, poll, vector,
};
use crate::assemble::ASSEMBLY_TILE;
use fs_exec::Cx;

impl EnthalpyStepLinearization<'_> {
    /// Include convection references and coefficients that depend on wall means.
    ///
    /// Rows and columns of `reference_derivatives` use the same unique `regions`
    /// order: entry `[i][j]` is `d(reference_i)/d(mean_temperature_j)` in K/K.
    /// `htc_wall_derivatives[i]` is `dh_i/d(mean_temperature_i)` in
    /// W/(m² K²). Supply zero for a constant law. These are derivatives of the
    /// ORIGINAL convection rows, even when radiation shares their faces.
    /// Radiation keeps its separate reservoir and its already-bound tangent.
    ///
    /// The consistent P1 residual adds `-dt*h_i*M_i*1*D_i*W*dT/dh` for
    /// references and `dt*M_i*(T-reference_i)*h'_i*W_i*dT/dh` for coefficients.
    /// The chart slope scales columns, including exact latent-plateau zeros.
    /// Both actions, Jacobi, all history/source/parameter pullbacks and physical
    /// convective-heat objective seeds use the same completed tangent.
    ///
    /// The caller owns these constitutive laws and must bind their values and
    /// derivatives at the accepted field. This method checks original port
    /// ownership, shapes, finite arithmetic and work; it does not validate an
    /// external derivative or re-evaluate a consumer's physical law. No primal
    /// solve, endpoint replacement or derivative of solver iterations occurs.
    ///
    /// One consuming attachment is allowed. At most 64 regions are admitted.
    /// `max_feedback_entries` caps all retained nodal factors: existing
    /// radiation factors plus four vectors per selected convection region.
    /// It is not total RSS; the bounded law matrix and disjoint face records
    /// are separate metadata. A failed attachment returns no modified owner.
    pub fn with_robin_feedback(
        mut self,
        cx: &Cx<'_>,
        regions: &[&str],
        reference_derivatives: &[Vec<f64>],
        htc_wall_derivatives: &[f64],
        max_feedback_entries: usize,
    ) -> Result<Self, EnthalpyAdjointError> {
        poll(cx, 0)?;
        let n = self.masses.len();
        let m = regions.len();
        let entries = m
            .checked_mul(2)
            .and_then(|added| self.feedback.len().checked_add(added))
            .and_then(|terms| terms.checked_mul(2))
            .and_then(|vectors| vectors.checked_mul(n));
        if self.robin_feedback.is_some()
            || m > 64
            || reference_derivatives.len() != m
            || htc_wall_derivatives.len() != m
            || entries.is_none_or(|count| count > max_feedback_entries)
            || entries
                .and_then(|count| count.checked_mul(std::mem::size_of::<f64>()))
                .is_none_or(|bytes| isize::try_from(bytes).is_err())
        {
            return Err(EnthalpyAdjointError::InvalidInput(
                "one Robin feedback attachment with at most 64 matching laws must fit the total factor-entry cap",
            ));
        }
        vector(cx, htc_wall_derivatives, m)?;
        for derivatives in reference_derivatives {
            vector(cx, derivatives, m)?;
        }
        let ports = bind_boundary_ports(cx, self.mesh, &self.convection_boundary, regions)?;
        let mut mean_h = Vec::with_capacity(m);
        for port in &ports {
            let mut weights = vec![0.0; n];
            for (vertices, area) in &port.faces {
                poll(cx, mean_h.len())?;
                for &vertex in vertices {
                    weights[vertex] = finite(weights[vertex] + area / port.area_m2 / 3.0)?;
                }
            }
            mean_h.push(self.temperature_pullback(cx, &weights)?);
        }
        let mut updates = Vec::with_capacity(2 * m);
        for (i, port) in ports.iter().enumerate() {
            let mut reference_left = vec![0.0; n];
            let mut coefficient_left = vec![0.0; n];
            for (vertices, area) in &port.faces {
                poll(cx, i)?;
                for (a, &row) in vertices.iter().enumerate() {
                    if self.prescribed[row] { continue; }
                    reference_left[row] =
                        finite(reference_left[row] - self.dt * port.htc_w_m2_k * (area / 3.0))?;
                    for (b, &column) in vertices.iter().enumerate() {
                        let mass = (area / 12.0) * if a == b { 2.0 } else { 1.0 };
                        coefficient_left[row] = finite(
                            coefficient_left[row]
                                + self.dt
                                    * mass
                                    * (self.primal.temperature[column] - port.reference_k),
                        )?;
                    }
                }
            }
            let mut reference_right = vec![0.0; n];
            for (j, weights) in mean_h.iter().enumerate() {
                let derivative = reference_derivatives[i][j];
                if derivative == 0.0 {
                    continue;
                }
                for (vertex, (&weight, entry)) in
                    weights.iter().zip(&mut reference_right).enumerate()
                {
                    if vertex % ASSEMBLY_TILE == 0 {
                        poll(cx, vertex)?;
                    }
                    *entry = finite(derivative.mul_add(weight, *entry))?;
                }
            }
            let mut coefficient_right = Vec::with_capacity(n);
            for (vertex, &weight) in mean_h[i].iter().enumerate() {
                if vertex % ASSEMBLY_TILE == 0 {
                    poll(cx, vertex)?;
                }
                coefficient_right.push(finite(htc_wall_derivatives[i] * weight)?);
            }
            updates.push(EnthalpyFeedback {
                left: reference_left,
                right: reference_right,
            });
            updates.push(EnthalpyFeedback {
                left: coefficient_left,
                right: coefficient_right,
            });
        }
        // Admit the COMPLETE diagonal. Two valid law terms can cancel an
        // intermediate zero; no partial-update Jacobi gate may reject that.
        let mut inverse_diagonal = Vec::with_capacity(n);
        for vertex in 0..n {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            let mut diagonal = if self.prescribed[vertex] {
                self.masses[vertex]
            } else {
                finite(self.masses[vertex] + self.dt * self.tangent.get(vertex, vertex))?
            };
            for update in self.feedback.iter().chain(&updates) {
                diagonal = finite(update.left[vertex].mul_add(update.right[vertex], diagonal))?;
            }
            if diagonal == 0.0 {
                return Err(EnthalpyAdjointError::InvalidInput(
                    "nonzero complete enthalpy Jacobian diagonal required for Jacobi",
                ));
            }
            inverse_diagonal.push(finite(1.0 / diagonal)?);
        }
        self.feedback.extend(updates);
        self.inverse_diagonal = inverse_diagonal;
        self.robin_feedback = Some(RobinStateFeedback {
            ports,
            reference_derivatives: reference_derivatives.to_vec(),
            htc_wall_derivatives: htc_wall_derivatives.to_vec(),
        });
        poll(cx, n)?;
        Ok(self)
    }
}

impl RobinStateFeedback {
    pub(super) fn add_heat_seed(
        &self,
        cx: &Cx<'_>,
        selected: &[RobinPort],
        heat_weights: &[f64],
        temperature: &[f64],
        seed: &mut [f64],
    ) -> Result<(), EnthalpyAdjointError> {
        let mut wall_seeds = vec![0.0; self.ports.len()];
        for (selected, &heat_weight) in selected.iter().zip(heat_weights) {
            poll(cx, 0)?;
            if heat_weight == 0.0 {
                continue;
            }
            let Some(i) = self
                .ports
                .iter()
                .position(|port| port.name == selected.name)
            else {
                continue;
            };
            let port = &self.ports[i];
            let mut mean = 0.0;
            for (vertices, area) in &port.faces {
                poll(cx, i)?;
                for &vertex in vertices {
                    mean = finite(mean + (area / port.area_m2 / 3.0) * temperature[vertex])?;
                }
            }
            let scale = finite(heat_weight * port.area_m2)?;
            wall_seeds[i] = finite(
                wall_seeds[i] + scale * self.htc_wall_derivatives[i] * (mean - port.reference_k),
            )?;
            for (j, value) in wall_seeds.iter_mut().enumerate() {
                *value =
                    finite(*value - scale * port.htc_w_m2_k * self.reference_derivatives[i][j])?;
            }
        }
        for (port, &weight) in self.ports.iter().zip(&wall_seeds) {
            if weight == 0.0 {
                continue;
            }
            for (vertices, area) in &port.faces {
                poll(cx, 0)?;
                let contribution = finite(weight * (area / port.area_m2 / 3.0))?;
                for &vertex in vertices {
                    seed[vertex] = finite(seed[vertex] + contribution)?;
                }
            }
        }
        poll(cx, selected.len())?;
        Ok(())
    }
}
