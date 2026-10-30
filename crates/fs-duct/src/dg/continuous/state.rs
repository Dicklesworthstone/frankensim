//! Original DG state adoption after a complete externally composed trajectory.
use super::*;

/// Integrated exchange at an original port. For a varying inlet flow, work is
/// integral p(t) Q(t) dt and cannot be reconstructed from the two mean values.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DgTrajectoryPort {
    pub id: PortId,
    pub pressure_impulse_pa_s: f64,
    pub work_into_duct_j: f64,
}
impl<M: DgMedium> DgContinuous<'_, M> {
    /// All original acoustic and wall coordinates, after source identity checks.
    pub fn energy_coordinates<'s>(&self, state: &'s DgState) -> Result<&'s [f64], DuctError> {
        if state.operator_identity != self.duct.operator_identity
            || state.tilde.len() != self.dimension()
        {
            return Err(DuctError::BadParameter {
                what: "continuous DG state belongs to another operator",
            });
        }
        Ok(&state.tilde)
    }

    /// Stage an endpoint from this SAME spatial operator's composed trajectory.
    /// The composing integrator supplies actual integrated pressure, port work
    /// and independently integrated original loss. This method checks their
    /// energy closure; it does not infer heat from a remainder or run a solver.
    /// Accepted history is never mutated. Initial energy, parameter work and
    /// prior cumulative work/loss survive; modal propagation anchors are cleared.
    #[allow(clippy::too_many_arguments)]
    pub fn candidate_from_trajectory(
        &self,
        previous: &DgState,
        endpoint: Vec<f64>,
        dt: f64,
        ports: &[DgTrajectoryPort],
        dissipated_j: f64,
        absolute_tolerance_j: f64,
        relative_tolerance: f64,
    ) -> Result<(DgState, DgStepOutput), DuctError> {
        self.energy_coordinates(previous)?;
        let invalid = |what| DuctError::BadParameter { what };
        if endpoint.len() != self.dimension()
            || endpoint.iter().any(|x| !x.is_finite())
            || !dt.is_finite()
            || dt <= 0.0
            || !dissipated_j.is_finite()
            || dissipated_j < 0.0
            || !absolute_tolerance_j.is_finite()
            || absolute_tolerance_j < 0.0
            || !relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&relative_tolerance)
            || ports.len() != self.ports().len()
        {
            return Err(invalid(
                "continuous DG trajectory dimensions, clock or energy tolerance",
            ));
        }
        let step = previous
            .step_count
            .checked_add(1)
            .ok_or_else(|| invalid("continuous DG accepted clock overflows"))?;
        let mut outputs = Vec::with_capacity(ports.len());
        let mut work = 0.0;
        for (id, column) in self.ports() {
            let mut matching = ports.iter().filter(|p| p.id == *id);
            let port = matching
                .next()
                .ok_or_else(|| invalid("continuous DG trajectory omits an original port"))?;
            if matching.next().is_some()
                || !port.pressure_impulse_pa_s.is_finite()
                || !port.work_into_duct_j.is_finite()
            {
                return Err(invalid(
                    "continuous DG trajectory has duplicate/nonfinite original port exchange",
                ));
            }
            work += port.work_into_duct_j;
            outputs.push(PortOutput {
                id: *id,
                conjugate_rate: port.pressure_impulse_pa_s / dt,
                delta_x: port.pressure_impulse_pa_s,
                end_conjugate: column.iter().zip(&endpoint).map(|(b, x)| b * x).sum(),
            });
        }
        let energy = 0.5 * endpoint.iter().map(|x| x * x).sum::<f64>();
        let before = previous.stored_energy_j();
        let residual = energy - before + dissipated_j - work;
        let tolerance = absolute_tolerance_j
            + relative_tolerance * (energy.abs() + before.abs() + work.abs() + dissipated_j);
        if ![energy, work, residual, tolerance]
            .iter()
            .all(|x| x.is_finite())
            || residual.abs() > tolerance
            || outputs
                .iter()
                .any(|p| !p.conjugate_rate.is_finite() || !p.end_conjugate.is_finite())
        {
            return Err(invalid(
                "continuous DG trajectory fails finite original-energy/port closure",
            ));
        }
        let (p, q) = self.duct.from_energy_coords(&endpoint);
        if p.iter().chain(&q).any(|x| !x.is_finite()) {
            return Err(invalid("continuous DG nodal observation overflows"));
        }
        let mut state = previous.clone();
        state.tilde = endpoint;
        state.p = p;
        state.q = q;
        state.nodal_dirty = false;
        state.modal_anchor_p.clear();
        state.modal_anchor_q.clear();
        state.modal_anchor_dt_bits = None;
        state.unforced_steps = 0;
        state.step_count = step;
        state.cumulative_port_work_j += work;
        state.cumulative_dissipated_j += dissipated_j;
        let balance = energy
            - state.initial_energy_j
            - state.cumulative_port_work_j
            - state.cumulative_parameter_work_j
            + state.cumulative_dissipated_j;
        let scale = state.initial_energy_j.abs()
            + state.cumulative_port_work_j.abs()
            + state.cumulative_parameter_work_j.abs();
        let relative = balance.abs() / scale.max(1e-300);
        if ![
            state.cumulative_port_work_j,
            state.cumulative_dissipated_j,
            balance,
            relative,
        ]
        .iter()
        .all(|x| x.is_finite())
        {
            return Err(invalid("continuous DG cumulative work/energy overflows"));
        }
        let output = DgStepOutput {
            port_outputs: outputs,
            ledger: EnergyLedgerRow {
                step,
                stored_energy_j: energy,
                port_work_step_j: work,
                cumulative_port_work_j: state.cumulative_port_work_j,
                cumulative_parameter_work_j: state.cumulative_parameter_work_j,
                cumulative_dissipated_j: state.cumulative_dissipated_j,
                energy_balance_residual_j: balance,
                relative_energy_drift: relative,
                balance_checked: true,
            },
        };
        Ok((state, output))
    }
}
