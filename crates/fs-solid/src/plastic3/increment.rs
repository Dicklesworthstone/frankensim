//! One Newton/load-path owner for the native small-strain 3D elements.
//!
//! Element operators own geometry, quadrature and constitutive history. This
//! module owns frozen-history trials, prescribed increments, force scaling,
//! sparse corrections, backtracking and atomic publication of accepted paths.

use super::{
    IncrementRecord, IncrementSettings, SmallStrainTetError, checkpoint, dynamics, norm,
    solve_increment_linear,
};
use fs_exec::Cx;
use fs_material::Voigt;
use fs_sparse::{Coo, Csr};

pub(super) trait IncrementOperator {
    type State: Clone + PartialEq;

    fn initial_states(&self) -> Vec<Self::State>;

    /// Internal force and a bound on element stress forces on free rows.
    fn internal(
        &self,
        states: &[Self::State],
        u: &[f64],
        prescribed: Option<&[bool]>,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, f64), SmallStrainTetError>;

    /// The complete derivative, including inertia when the operator has it.
    fn tangent(
        &self,
        states: &[Self::State],
        u: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Csr, SmallStrainTetError>;

    fn tangent_symmetric(&self) -> bool;

    /// Commit only this converged increment; return changed ELEMENT count.
    /// Final stresses must be evaluated from the previous history, before the
    /// update, so a rate-dependent law cannot apply a second physical interval.
    fn commit(
        &self,
        states: &mut [Self::State],
        u: &[f64],
        final_stress: Option<&mut Vec<Voigt>>,
        cx: &Cx<'_>,
    ) -> Result<usize, SmallStrainTetError>;

    /// A geometry/law domain refusal on a proposal may be repairable by a
    /// shorter step. Initial and tangent refusals always propagate; cancellation
    /// is never recoverable. Existing tetrahedral behavior remains the default.
    fn recoverable_trial_error(&self, _error: &SmallStrainTetError) -> bool {
        false
    }
}

pub(super) struct IncrementControls<'a> {
    pub dofs: usize,
    pub prescribed_m: &'a [(usize, f64)],
    pub nodal_forces_n: &'a [(usize, f64)],
    pub load_path: &'a [f64],
    pub settings: IncrementSettings,
}

pub(super) struct IncrementSolution<S> {
    pub displacement_m: Vec<f64>,
    pub point_stress: Vec<Voigt>,
    pub point_states: Vec<S>,
    pub increments: Vec<IncrementRecord>,
}

/// The wrapper validates topology, budgets, controls and initial array shapes
/// before entry. Candidate state stays private until the entire path succeeds.
#[allow(clippy::too_many_lines)] // shared increment / Newton / backtracking loop
pub(super) fn solve<O: IncrementOperator>(
    operator: &O,
    controls: IncrementControls<'_>,
    initial: Option<(&[f64], &[O::State])>,
    inertia: Option<&dynamics::Inertia<'_>>,
    cx: &Cx<'_>,
) -> Result<IncrementSolution<O::State>, SmallStrainTetError> {
    checkpoint(cx)?;
    let IncrementControls {
        dofs,
        prescribed_m,
        nodal_forces_n,
        load_path,
        settings,
    } = controls;
    let mut is_prescribed = vec![false; dofs];
    for &(dof, _) in prescribed_m {
        is_prescribed[dof] = true;
    }
    let free: Vec<usize> = (0..dofs).filter(|&d| !is_prescribed[d]).collect();
    let mut row_of = vec![usize::MAX; dofs];
    for (row, &dof) in free.iter().enumerate() {
        row_of[dof] = row;
    }
    let (mut u, mut states) = match initial {
        Some((displacement, states)) => (displacement.to_vec(), states.to_vec()),
        None => (vec![0.0; dofs], operator.initial_states()),
    };
    let mut increments = Vec::with_capacity(load_path.len());
    let mut point_stress = Vec::with_capacity(states.len());
    for (increment, &load) in load_path.iter().enumerate() {
        let mut f_ext = vec![0.0f64; dofs];
        for &(dof, value) in nodal_forces_n {
            f_ext[dof] += load * value;
        }
        let target: Vec<(usize, f64)> = prescribed_m
            .iter()
            .map(|&(dof, value)| (dof, load * value))
            .collect();
        if f_ext
            .iter()
            .chain(target.iter().map(|(_, v)| v))
            .any(|v| !v.is_finite())
        {
            return Err(SmallStrainTetError::InvalidInput {
                what: "load factor produces non-finite forces or prescribed displacements".into(),
            });
        }
        let scale_rows = if inertia.is_none() {
            Some(is_prescribed.as_slice())
        } else {
            None
        };
        let external_free_norm = free.iter().fold(0.0_f64, |n, &d| n.hypot(f_ext[d]));
        let mut history = Vec::new();
        let mut fallbacks = 0usize;
        let mut converged = false;
        let (mut f_int, initial_scale) = operator.internal(&states, &u, scale_rows, cx)?;
        let mut static_force_scale = initial_scale.max(external_free_norm);
        let mut elastic_free_norm = if inertia.is_some() {
            free.iter().fold(0.0_f64, |n, &d| n.hypot(f_int[d]))
        } else {
            0.0
        };
        if let Some(inertia) = inertia {
            inertia.add_force(&u, &mut f_int, cx)?;
        }
        for _ in 0..settings.max_newton_iterations {
            checkpoint(cx)?;
            let residual: Vec<f64> = free.iter().map(|&d| f_int[d] - f_ext[d]).collect();
            let pending: Vec<(usize, f64)> = target
                .iter()
                .map(|&(dof, value)| (dof, value - u[dof]))
                .collect();
            let pending_any = pending.iter().any(|&(_, dp)| dp != 0.0);
            let scale = if let Some(inertia) = inertia {
                inertia.force_scale(&u, &free, elastic_free_norm, external_free_norm, cx)?
            } else {
                static_force_scale.max(f64::MIN_POSITIVE)
            };
            let relative = norm(&residual) / scale;
            history.push(relative);
            if !pending_any && relative <= settings.relative_tolerance {
                converged = true;
                break;
            }
            // A converged iterate needs no new constitutive tangent or global
            // sparse assembly, including a fully prescribed body.
            let full = operator.tangent(&states, &u, cx)?;
            let mut dp_full = vec![0.0f64; dofs];
            for &(dof, dp) in &pending {
                dp_full[dof] = dp;
            }
            let mut coo = Coo::new(free.len(), free.len());
            let mut rhs: Vec<f64> = residual.iter().map(|r| -r).collect();
            for (row, &dof) in free.iter().enumerate() {
                let (cols, vals) = full.row(dof);
                for (&c, &v) in cols.iter().zip(vals) {
                    if is_prescribed[c] {
                        rhs[row] = (-v).mul_add(dp_full[c], rhs[row]);
                    } else {
                        coo.push(row, row_of[c], v);
                    }
                }
            }
            let (direction, fallback) = solve_increment_linear(
                &coo.assemble(),
                &rhs,
                increment,
                settings,
                operator.tangent_symmetric(),
                cx,
            )?;
            fallbacks += usize::from(fallback);
            let before = norm(&residual);
            let mut alpha = 1.0f64;
            let mut accepted = None;
            for _ in 0..=settings.max_backtracks {
                let mut trial = u.clone();
                for (row, &dof) in free.iter().enumerate() {
                    trial[dof] = alpha.mul_add(direction[row], trial[dof]);
                }
                for &(dof, dp) in &pending {
                    trial[dof] = alpha.mul_add(dp, trial[dof]);
                }
                let (mut f_trial, trial_scale) =
                    match operator.internal(&states, &trial, scale_rows, cx) {
                        Ok(value) => value,
                        Err(error)
                            if error != SmallStrainTetError::Cancelled
                                && operator.recoverable_trial_error(&error) =>
                        {
                            alpha *= 0.5;
                            continue;
                        }
                        Err(error) => return Err(error),
                    };
                let trial_elastic_free_norm = if inertia.is_some() {
                    free.iter().fold(0.0_f64, |n, &d| n.hypot(f_trial[d]))
                } else {
                    0.0
                };
                if let Some(inertia) = inertia {
                    inertia.add_force(&trial, &mut f_trial, cx)?;
                }
                let after: Vec<f64> = free.iter().map(|&d| f_trial[d] - f_ext[d]).collect();
                if pending_any || norm(&after) < before || alpha < 1e-3 {
                    accepted = Some((trial, f_trial, trial_elastic_free_norm, trial_scale));
                    break;
                }
                alpha *= 0.5;
            }
            let Some((trial, f_trial, trial_elastic_free_norm, trial_scale)) = accepted else {
                return Err(SmallStrainTetError::NewtonStalled { increment, history });
            };
            u = trial;
            f_int = f_trial;
            elastic_free_norm = trial_elastic_free_norm;
            // Rejected proposals cannot loosen the force gate. Retain the
            // largest accepted physical scale through unloading.
            static_force_scale = static_force_scale.max(trial_scale);
        }
        if !converged {
            return Err(SmallStrainTetError::NewtonStalled { increment, history });
        }
        let final_stress = if increment + 1 == load_path.len() {
            Some(&mut point_stress)
        } else {
            None
        };
        let evolving = operator.commit(&mut states, &u, final_stress, cx)?;
        increments.push(IncrementRecord {
            load_factor: load,
            residual_history: history,
            reactions_n: prescribed_m
                .iter()
                .map(|&(dof, _)| f_int[dof] - f_ext[dof])
                .collect(),
            evolving_elements: evolving,
            indefinite_fallbacks: fallbacks,
        });
    }
    checkpoint(cx)?;
    Ok(IncrementSolution {
        displacement_m: u,
        point_stress,
        point_states: states,
        increments,
    })
}
