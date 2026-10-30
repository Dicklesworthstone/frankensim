//! Source-positive molecular heat on the original complete DG trajectory.
use super::*;
use fs_phs::PreparedStepError;
use fs_phs::linear_exact::{DenseExponentialBudget, PreparedSourceLoss};

fn integral_error(error: PreparedStepError) -> DuctError {
    match error {
        PreparedStepError::Cancelled => DuctError::Cancelled { completed_steps: 0 },
        PreparedStepError::Solver(reason) => DuctError::MolecularHeatUnresolved { reason },
    }
}

/// Read-only molecular loss integral tied to one original DG operator and dt.
/// The full native generator, source rows and port ordering remain unchanged.
/// No physical state, second time step, total-loss allocation or endpoint-power
/// approximation is held here. Warm evaluation is a positive quadratic source
/// contraction; a changed gas, material, geometry or clock needs a new cache.
#[derive(Clone, Debug)]
pub struct DgMolecularHeat {
    operator_identity: ContentHash,
    ports: Vec<PortId>,
    integral: PreparedSourceLoss,
}

impl DgMolecularHeat {
    /// Original complete operator identity, including all material memories.
    pub const fn operator_identity(&self) -> ContentHash {
        self.operator_identity
    }
    /// The same held-input interval as the original native propagator [s].
    pub fn duration(&self) -> f64 {
        self.integral.duration()
    }

    /// Integrate molecular slip squares from this accepted INITIAL state under
    /// the native step's held ports. The supplied state is never advanced.
    /// Thermal wall conduction, viscous drag, solid loss and outlet absorption
    /// are excluded; their destinations cannot be inferred from molecular heat.
    pub fn integrated_heat_j(
        &self,
        state: &DgState,
        ports: &[PortInput],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, DuctError> {
        poll(cancelled)?;
        if state.operator_identity != self.operator_identity {
            return Err(invalid(
                "molecular heat requires the original complete DG state",
            ));
        }
        let inputs = resolve_efforts(ports, &self.ports)?;
        self.integral
            .evaluate_controlled(&state.tilde, &inputs, cancelled)
            .map_err(integral_error)
    }
}

impl DgDuct {
    /// Prepare the exact original molecular dissipation channel for one native
    /// propagator. The caller's explicit 2n-coordinate/byte/work allowance is
    /// additional to propagation admission; at most 128 complete DG states are
    /// admitted by the shared positive integration owner. No state is truncated.
    ///
    /// Each retained reservoir contributes its actual source square
    /// `sqrt(rate)*(z-b*(y+a*P*x_pressure))`. Original spatial pressure projection,
    /// local radial memory and every other coupled coordinate remain present in
    /// the trajectory. No molecular source gives an exactly zero integral.
    pub fn prepare_molecular_heat(
        &self,
        propagator: &DgPropagator,
        budget: DenseExponentialBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgMolecularHeat, DuctError> {
        poll(cancelled)?;
        let n = self.total_dofs();
        let m = self.ports.len();
        if propagator.operator_identity != self.operator_identity
            || n > 128
            || 2 * n > budget.max_dimension
        {
            return Err(invalid(
                "molecular heat requires its native propagator and complete bounded state",
            ));
        }
        let count = self
            .molecular
            .as_ref()
            .map_or(Some(0usize), |memory| {
                memory.modes.iter().try_fold(0usize, |count, mode| {
                    count.checked_add(
                        mode.rates
                            .len()
                            .checked_mul(mode.projection.width(self.field_dofs))?,
                    )
                })
            })
            .ok_or_else(|| invalid("molecular heat source extent overflow"))?;
        let bytes = count
            .checked_mul(n)
            .and_then(|v| v.checked_add(n.checked_mul(m)?))
            .and_then(|v| v.checked_mul(8));
        if count > n || bytes.is_none_or(|v| v > budget.max_dense_bytes) {
            return Err(invalid(
                "molecular heat source rows exceed the complete byte allowance",
            ));
        }
        let mut rows = Vec::with_capacity(count * n);
        if let Some(memory) = &self.molecular {
            for mode in &memory.modes {
                let width = mode.projection.width(self.field_dofs);
                for i in 0..width {
                    for (j, (&b, &rate)) in
                        mode.reservoir_scales.iter().zip(&mode.rates).enumerate()
                    {
                        poll(cancelled)?;
                        let root = det::sqrt(rate);
                        let mut row = vec![0.0; n];
                        row[mode.start + (j + 1) * width + i] = root;
                        row[mode.start + i] = -root * b;
                        for (column, weight) in
                            mode.projection.entries(i, -root * b * mode.pressure_scale)
                        {
                            row[column] += weight;
                        }
                        rows.extend_from_slice(&row);
                    }
                }
            }
        }
        let mut inputs = vec![0.0; n * m];
        for (port, (_, column)) in self.ports.iter().enumerate() {
            poll(cancelled)?;
            for i in 0..n {
                inputs[i * m + port] = column[i];
            }
        }
        let integral = PreparedSourceLoss::new_controlled(
            n,
            m,
            &self.l_tilde,
            &inputs,
            &rows,
            propagator.dt,
            budget,
            cancelled,
        )
        .map_err(integral_error)?;
        Ok(DgMolecularHeat {
            operator_identity: self.operator_identity,
            ports: self.ports.iter().map(|(id, _)| *id).collect(),
            integral,
        })
    }
}
