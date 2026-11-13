//! Physical fluid-viscosity heat on the complete original acoustic trajectory.
use super::*;
use fs_phs::PreparedStepError;
use fs_phs::linear_exact::{DenseExponentialBudget, PreparedSourceLoss};

fn invalid(what: &'static str) -> DuctError {
    DuctError::BadParameter { what }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DuctError> {
    if cancelled() {
        Err(DuctError::Cancelled { completed_steps: 0 })
    } else {
        Ok(())
    }
}

fn integral_error(error: PreparedStepError) -> DuctError {
    match error {
        PreparedStepError::Cancelled => DuctError::Cancelled { completed_steps: 0 },
        PreparedStepError::Solver(reason) => DuctError::ViscousHeatUnresolved { reason },
    }
}

/// Original Poiseuille drag and viscous diffusion heat over one native DG step.
/// Every acoustic, thermal, molecular and structural coordinate remains in the
/// trajectory. Only fluid viscous source squares are observed: no thermal-wall
/// heat partition, molecular loss, solid loss or boundary absorption is inferred.
#[derive(Clone, Debug)]
pub struct DgViscousHeat {
    operator_identity: ContentHash,
    ports: Vec<PortId>,
    integral: PreparedSourceLoss,
}

impl DgViscousHeat {
    /// Complete original operator, including all material memories and ports.
    pub const fn operator_identity(&self) -> ContentHash {
        self.operator_identity
    }

    /// The exact held-input interval of the native propagator [s].
    pub fn duration(&self) -> f64 {
        self.integral.duration()
    }

    /// Integrate only fluid viscous heat from the accepted INITIAL state under
    /// the native step's held ports. This observation advances no state. A new
    /// operator or duration requires a new preparation; warm evaluation is a
    /// positive quadratic contraction, never an allocation of total energy loss.
    pub fn integrated_heat_j(
        &self,
        state: &DgState,
        ports: &[PortInput],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, DuctError> {
        poll(cancelled)?;
        if state.operator_identity != self.operator_identity {
            return Err(invalid(
                "viscous heat requires the original complete DG state",
            ));
        }
        let inputs = resolve_efforts(ports, &self.ports)?;
        self.integral
            .evaluate_controlled(&state.tilde, &inputs, cancelled)
            .map_err(integral_error)
    }
}

fn next_source_row(rows: &mut Vec<f64>, n: usize) -> &mut [f64] {
    let start = rows.len();
    rows.resize(start + n, 0.0);
    &mut rows[start..]
}

impl<M: DgMedium> DgDuct<M> {
    /// Prepare the actual fluid-viscosity source integral for one original
    /// propagator. Uniform circular diffusion retains its exact mass basis;
    /// spatial diffusion retains every original material quadrature projection.
    /// Each source is `drag*x²` or `(sqrt(gain)*x-sqrt(rate)*memory)²` from
    /// the same series operator that advances the flow. Thermal shunts, molecular
    /// reservoirs, solid Kelvin/Maxwell arms and terminals contribute no rows.
    ///
    /// The existing positive integration owner admits at most 128 COMPLETE
    /// states and at most n source rows under explicit dimension/byte/work
    /// allowances. No state or quadrature is truncated to fit those limits.
    /// Absent viscous diffusion returns an exactly zero source. This factory
    /// observes heat; the receiving thermal owner supplies its deposition and
    /// material update and removes the internal transfer from exported loss once.
    pub fn prepare_viscous_heat(
        &self,
        propagator: &DgPropagator,
        budget: DenseExponentialBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgViscousHeat, DuctError> {
        poll(cancelled)?;
        let n = self.total_dofs();
        let m = self.field_dofs;
        let ports = self.ports.len();
        if propagator.operator_identity != self.operator_identity
            || n > 128
            || 2 * n > budget.max_dimension
        {
            return Err(invalid(
                "viscous heat requires its native propagator and complete bounded state",
            ));
        }
        let count = match &self.viscothermal {
            Some(ViscothermalMemory::Uniform(model)) => m.checked_mul(
                model.series_terms().len() + usize::from(model.series_drag_per_s() > 0.0),
            ),
            Some(ViscothermalMemory::Spatial(wall)) => {
                wall.viscous_channels(m)
                    .try_fold(0usize, |count, (_, _, drag, terms, _)| {
                        count.checked_add(terms.len() + usize::from(drag > 0.0))
                    })
            }
            None => Some(0),
        }
        .ok_or_else(|| invalid("viscous heat source extent overflow"))?;
        let bytes = count
            .checked_mul(n)
            .and_then(|v| v.checked_add(n.checked_mul(ports)?))
            .and_then(|v| v.checked_mul(8));
        if count > n || bytes.is_none_or(|v| v > budget.max_dense_bytes) {
            return Err(invalid(
                "viscous heat source rows exceed the complete byte allowance",
            ));
        }
        let mut rows = Vec::with_capacity(count * n);
        match &self.viscothermal {
            Some(ViscothermalMemory::Uniform(model)) => {
                if model.series_drag_per_s() > 0.0 {
                    let root = det::sqrt(model.series_drag_per_s());
                    for i in 0..m {
                        poll(cancelled)?;
                        next_source_row(&mut rows, n)[m + i] = root;
                    }
                }
                for (arm, term) in model.series_terms().iter().enumerate() {
                    for i in 0..m {
                        poll(cancelled)?;
                        let row = next_source_row(&mut rows, n);
                        row[m + i] = det::sqrt(term.gain_per_s);
                        row[2 * m + arm * m + i] = -det::sqrt(term.rate_per_s);
                    }
                }
            }
            Some(ViscothermalMemory::Spatial(wall)) => {
                for (field_start, projection, drag, terms, memory_start) in wall.viscous_channels(m)
                {
                    if drag > 0.0 {
                        poll(cancelled)?;
                        let row = next_source_row(&mut rows, n);
                        let root = det::sqrt(drag);
                        for (i, weight) in projection.iter().enumerate() {
                            row[field_start + i] = root * weight;
                        }
                    }
                    for (arm, term) in terms.iter().enumerate() {
                        poll(cancelled)?;
                        let row = next_source_row(&mut rows, n);
                        let root = det::sqrt(term.gain_per_s);
                        for (i, weight) in projection.iter().enumerate() {
                            row[field_start + i] = root * weight;
                        }
                        row[memory_start + arm] = -det::sqrt(term.rate_per_s);
                    }
                }
            }
            None => {}
        }
        let mut inputs = vec![0.0; n * ports];
        for (port, (_, column)) in self.ports.iter().enumerate() {
            poll(cancelled)?;
            for i in 0..n {
                inputs[i * ports + port] = column[i];
            }
        }
        let integral = PreparedSourceLoss::new_controlled(
            n,
            ports,
            &self.l_tilde,
            &inputs,
            &rows,
            propagator.dt,
            budget,
            cancelled,
        )
        .map_err(integral_error)?;
        Ok(DgViscousHeat {
            operator_identity: self.operator_identity,
            ports: self.ports.iter().map(|(id, _)| *id).collect(),
            integral,
        })
    }
}
