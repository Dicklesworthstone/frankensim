//! Independent material heat channels along the original complete trajectory.
use super::*;
use fs_phs::linear_exact::integrated_resistance_loss;

/// The SAME full-state step plus original quadrature-point heat channels.
/// The residual is disclosed, not reassigned to any material as invented heat.
#[derive(Clone, Debug)]
pub struct MaterialHeatStep {
    pub dynamics: StepRecord,
    /// Kelvin dashpot plus all Maxwell dashpot heat, in material_points order [J].
    pub point_heat_j: Vec<f64>,
    /// Sum of independently integrated point channels minus total R heat [J].
    pub attribution_residual_j: f64,
}

impl TensorMaterialHistory<'_> {
    /// Advance while retaining separate physical heat for every full-3D point.
    /// Each PSD resistance channel is integrated along the SAME coupled exact
    /// trajectory by fs-phs; state, clock and histories publish only after all
    /// channels succeed and their sum matches total loss within the explicit
    /// absolute attribution allowance. No work-minus-storage residual supplies heat.
    ///
    /// Plane-stress point heat is outside this method's explicit scope; its
    /// existing dynamics and total heat path remain available through advance.
    pub fn advance_with_material_heat(
        &mut self, force_n: &[f64], dt_s: f64,
        absolute_attribution_tolerance_j: f64, gate: &CancelGate,
    ) -> Result<MaterialHeatStep> {
        poll(gate)?;
        if !absolute_attribution_tolerance_j.is_finite() || absolute_attribution_tolerance_j <= 0.0 {
            return Err(bad("positive finite material heat attribution tolerance required"));
        }
        if !self.body.plane_points.is_empty() {
            return Err(bad("separate point heat requires the full-3D material history"));
        }
        let mut candidate = self.clone();
        let dynamics = candidate.advance(force_n, dt_s, gate)?;
        let n = self.body.coordinates();
        let states = self.body.state_dimension();
        let mut at = 2*n;
        let mut heat = Vec::with_capacity(self.body.points.len());
        let mut total = 0.0;
        for point in &self.body.points {
            poll(gate)?;
            let mut resistance = vec![0.0; states*states];
            for i in 0..n { for j in 0..n {
                resistance[(n+i)*states+n+j] = point.viscosity_matrix()[i*n+j];
            }}
            for branch in point.branches() {
                for _ in branch.rows() {
                    resistance[at*states+at] = branch.relaxation_time_s().recip();
                    at += 1;
                }
            }
            let value = if resistance.iter().all(|x| *x == 0.0) { 0.0 } else {
                integrated_resistance_loss(
                    &self.body.system, &self.body.split, &self.state, force_n, dt_s, &resistance,
                ).map_err(TensorMaterialError::Propagation)?
            };
            if !value.is_finite() || value < 0.0 {
                return Err(bad("material point heat left its finite nonnegative domain"));
            }
            total += value;
            heat.push(value);
        }
        let residual = total-dynamics.dissipated;
        if !total.is_finite() || !residual.is_finite() {
            return Err(bad("material heat attribution overflow"));
        }
        if residual.abs() > absolute_attribution_tolerance_j {
            return Err(bad("material heat attribution exceeds its absolute energy allowance"));
        }
        poll(gate)?;
        *self = candidate;
        Ok(MaterialHeatStep { dynamics, point_heat_j: heat, attribution_residual_j: residual })
    }
}
