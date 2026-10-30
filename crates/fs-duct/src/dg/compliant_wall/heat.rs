//! Physical solid heat on the original complete acoustic/wall trajectory.
use super::*;
use fs_phs::PreparedStepError;
use fs_phs::linear_exact::{DenseExponentialBudget, PreparedSourceLoss};

fn integral_error(error: PreparedStepError) -> DuctError {
    match error {
        PreparedStepError::Cancelled => DuctError::Cancelled { completed_steps: 0 },
        PreparedStepError::Solver(reason) => DuctError::CompliantWallHeatUnresolved { reason },
    }
}

/// Read-only solid-material heat integral for one original DG operator and dt.
/// Kelvin viscosity and retained Maxwell arms contribute their own physical
/// source squares. Gas viscosity, gas thermal/molecular relaxation, numerical
/// loss and terminal absorption are excluded. No state is advanced here.
#[derive(Clone, Debug)]
pub struct DgCompliantWallHeat {
    operator_identity: ContentHash,
    ports: Vec<PortId>,
    integral: PreparedSourceLoss,
}

impl DgCompliantWallHeat {
    /// Complete original operator, including gas, wall and all retained memory.
    pub const fn operator_identity(&self) -> ContentHash {
        self.operator_identity
    }

    /// Original native held-input interval [s].
    pub fn duration(&self) -> f64 {
        self.integral.duration()
    }

    /// Integrate only solid Kelvin/Maxwell heat from the accepted INITIAL state
    /// under the same held ports as the native step. Warm evaluation contracts
    /// the prepared positive source form; a changed material, gas, geometry or
    /// duration requires a new preparation. No aggregate deficit is allocated.
    pub fn integrated_heat_j(
        &self,
        state: &DgState,
        ports: &[PortInput],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, DuctError> {
        check_cancelled(cancelled)?;
        if state.operator_identity != self.operator_identity {
            return Err(invalid(
                "solid wall heat requires the original complete DG state",
            ));
        }
        let inputs = resolve_efforts(ports, &self.ports)?;
        self.integral
            .evaluate_controlled(&state.tilde, &inputs, cancelled)
            .map_err(integral_error)
    }
}

fn append_square(
    rows: &mut Vec<f64>,
    n: usize,
    coordinate: usize,
    rate: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), DuctError> {
    check_cancelled(cancelled)?;
    if !rate.is_finite() || rate <= 0.0 || coordinate >= n {
        return Err(invalid("solid wall heat source is unrepresentable"));
    }
    let mut row = vec![0.0; n];
    row[coordinate] = det::sqrt(rate);
    rows.extend_from_slice(&row);
    Ok(())
}

impl DgDuct {
    /// Prepare Kelvin and Maxwell source squares on the COMPLETE original
    /// trajectory. The native generator, state coordinates and held-port order
    /// remain unchanged, including all gas and terminal physics. Each regional
    /// coefficient and arm retains its original element support.
    ///
    /// Requires an attached compliant wall and its native propagator. At most
    /// 128 complete states and the explicit additional dense-exponential budget
    /// are admitted, as in the existing molecular source owner. A purely elastic
    /// wall gives exactly zero heat. This is not gas-to-wall heat partition or a
    /// finite-temperature material update; the receiving thermal owner supplies
    /// that separate transaction and accounts for the transfer once.
    pub fn prepare_compliant_wall_heat(
        &self,
        propagator: &DgPropagator,
        budget: DenseExponentialBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgCompliantWallHeat, DuctError> {
        check_cancelled(cancelled)?;
        let n = self.total_dofs();
        let m = self.ports.len();
        if propagator.operator_identity != self.operator_identity
            || n > 128
            || 2 * n > budget.max_dimension
        {
            return Err(invalid(
                "solid wall heat requires its native propagator and complete bounded state",
            ));
        }
        let wall = self.compliant_wall.as_ref().ok_or_else(|| {
            invalid("solid wall heat requires an original compliant material wall")
        })?;
        let local = self.mesh.order + 1;
        let mut count = 0usize;
        for e in 0..self.elements.len() {
            check_cancelled(cancelled)?;
            if wall.pin_at(e).resistance > 0.0 {
                count = count
                    .checked_add(local)
                    .ok_or_else(|| invalid("solid wall heat source extent overflow"))?;
            }
        }
        for branch in &wall.relaxation {
            count = count
                .checked_add(branch.dofs)
                .ok_or_else(|| invalid("solid wall heat source extent overflow"))?;
        }
        let bytes = count
            .checked_mul(n)
            .and_then(|v| v.checked_add(n.checked_mul(m)?))
            .and_then(|v| v.checked_mul(8));
        if count > n || bytes.is_none_or(|v| v > budget.max_dense_bytes) {
            return Err(invalid(
                "solid wall heat source rows exceed the complete byte allowance",
            ));
        }
        let mut rows = Vec::with_capacity(count * n);
        for e in 0..self.elements.len() {
            let pin = wall.pin_at(e);
            if pin.resistance == 0.0 {
                continue;
            }
            for i in 0..local {
                append_square(
                    &mut rows,
                    n,
                    wall.start + self.field_dofs + e * local + i,
                    pin.resistance / pin.surface_density,
                    cancelled,
                )?;
            }
        }
        for branch in &wall.relaxation {
            for coordinate in branch.start..branch.start + branch.dofs {
                append_square(&mut rows, n, coordinate, branch.rate_per_s, cancelled)?;
            }
        }
        let mut inputs = vec![0.0; n * m];
        for (port, (_, column)) in self.ports.iter().enumerate() {
            check_cancelled(cancelled)?;
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
        Ok(DgCompliantWallHeat {
            operator_identity: self.operator_identity,
            ports: self.ports.iter().map(|(id, _)| *id).collect(),
            integral,
        })
    }
}

#[cfg(test)]
mod tests;
