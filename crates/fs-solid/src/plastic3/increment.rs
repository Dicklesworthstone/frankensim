//! One Newton/load-path owner for native small- and finite-strain elements.
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

/// Borrowed contributions from the original time-integration owner. Material
/// assembly and the Newton loop stay shared; the inertia owner supplies its
/// own momentum force, exact derivative and physical free-row force scale.
pub(crate) trait IncrementInertia {
    fn add_force(
        &self,
        u: &[f64],
        force: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), crate::inertia::InertiaError>;

    fn add_tangent(
        &self,
        matrix: &mut Coo,
        cx: &Cx<'_>,
    ) -> Result<(), crate::inertia::InertiaError>;

    fn force_scale(
        &self,
        u: &[f64],
        free: &[usize],
        elastic_free_norm: f64,
        external_free_norm: f64,
        cx: &Cx<'_>,
    ) -> Result<f64, crate::inertia::InertiaError>;
}

impl IncrementInertia for crate::inertia::BackwardEulerInertia<'_> {
    fn add_force(
        &self,
        u: &[f64],
        force: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), crate::inertia::InertiaError> {
        Self::add_force(self, u, force, cx)
    }

    fn add_tangent(
        &self,
        matrix: &mut Coo,
        cx: &Cx<'_>,
    ) -> Result<(), crate::inertia::InertiaError> {
        Self::add_tangent(self, matrix, cx)
    }

    fn force_scale(
        &self,
        u: &[f64],
        free: &[usize],
        elastic_free_norm: f64,
        external_free_norm: f64,
        cx: &Cx<'_>,
    ) -> Result<f64, crate::inertia::InertiaError> {
        Self::force_scale(self, u, free, elastic_free_norm, external_free_norm, cx)
    }
}

impl IncrementInertia for crate::inertia::MidpointInertia<'_> {
    fn add_force(
        &self,
        u: &[f64],
        force: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), crate::inertia::InertiaError> {
        Self::add_force(self, u, force, cx)
    }

    fn add_tangent(
        &self,
        matrix: &mut Coo,
        cx: &Cx<'_>,
    ) -> Result<(), crate::inertia::InertiaError> {
        Self::add_tangent(self, matrix, cx)
    }

    fn force_scale(
        &self,
        u: &[f64],
        free: &[usize],
        elastic_free_norm: f64,
        external_free_norm: f64,
        cx: &Cx<'_>,
    ) -> Result<f64, crate::inertia::InertiaError> {
        Self::force_scale(self, u, free, elastic_free_norm, external_free_norm, cx)
    }
}

pub(crate) trait IncrementOperator {
    type State: Clone + PartialEq;

    fn initial_states(&self) -> Vec<Self::State>;

    /// Translational components per node, used to remove a common translation
    /// from the connected force scale. Existing spatial elements retain 3-D
    /// indexing; plane sheets use the same Newton owner with two components.
    fn components_per_node(&self) -> usize {
        3
    }

    /// Coordinate used only to set the connected equilibrium-force scale.
    /// Translational operators remove a common translation. Shells retain
    /// absolute director rotations: subtracting a common rotation alone is
    /// not a rigid motion without its compatible translation field.
    fn deformation_coordinate(&self, u: &[f64], row: usize, column: usize) -> f64 {
        let components = self.components_per_node();
        u[column] - u[components * (row / components) + column % components]
    }

    /// Internal force and an upper cap on element stress forces on free rows.
    /// The cap alone must not normalize equilibrium: it can include stresses
    /// orthogonal to the row's work-conjugate strain direction.
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
    /// is never recoverable. Element operators opt in explicitly; the default
    /// propagates every refusal.
    fn recoverable_trial_error(&self, _error: &SmallStrainTetError) -> bool {
        false
    }
}

pub(crate) struct IncrementControls<'a> {
    pub dofs: usize,
    pub prescribed_m: &'a [(usize, f64)],
    pub nodal_forces_n: &'a [(usize, f64)],
    pub load_path: &'a [f64],
    pub settings: IncrementSettings,
}

pub(crate) struct IncrementSolution<S> {
    pub displacement_m: Vec<f64>,
    pub point_stress: Vec<Voigt>,
    pub point_states: Vec<S>,
    pub increments: Vec<IncrementRecord>,
}

/// The wrapper validates topology, budgets, controls and initial array shapes
/// before entry. Candidate state stays private until the entire path succeeds.
pub(crate) fn solve<O: IncrementOperator>(
    operator: &O,
    controls: IncrementControls<'_>,
    initial: Option<(&[f64], &[O::State])>,
    inertia: Option<&dynamics::Inertia<'_>>,
    cx: &Cx<'_>,
) -> Result<IncrementSolution<O::State>, SmallStrainTetError> {
    solve_with_inertia(
        operator,
        controls,
        initial,
        inertia.map(|inertia| inertia as &dyn IncrementInertia),
        cx,
    )
}

/// Same accepted-state, backtracking and reaction owner with an explicitly
/// selected inertial rule. The operator's tangent must include this same
/// inertia. A path-mean material law remains the caller's responsibility.
#[allow(clippy::too_many_lines)] // shared increment / Newton / backtracking loop
pub(crate) fn solve_with_inertia<O: IncrementOperator>(
    operator: &O,
    controls: IncrementControls<'_>,
    initial: Option<(&[f64], &[O::State])>,
    inertia: Option<&dyn IncrementInertia>,
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
    let components = operator.components_per_node();
    if !(1..=6).contains(&components) || !dofs.is_multiple_of(components) {
        return Err(SmallStrainTetError::InvalidInput {
            what: "increment coordinate layout requires one to six components per node".into(),
        });
    }
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
        // The start is the accepted configuration evaluated under this
        // increment's laws. A law over an interval dt sees there no total
        // chord against the interval's whole free-strain change, and may
        // refuse that rate although the increment, which ends on its targets,
        // lies inside its domain. With every coordinate prescribed the targets
        // are the endpoint and the start decides nothing, so such a refusal
        // moves the start to the targets. Every other refusal, and any refusal
        // at the targets, still propagates.
        let (mut f_int, mut stress_force_cap) = match operator.internal(&states, &u, scale_rows, cx)
        {
            Ok(value) => value,
            Err(error)
                if free.is_empty()
                    && error != SmallStrainTetError::Cancelled
                    && operator.recoverable_trial_error(&error)
                    && target.iter().any(|&(dof, value)| value - u[dof] != 0.0) =>
            {
                for &(dof, value) in &target {
                    u[dof] = value;
                }
                operator.internal(&states, &u, scale_rows, cx)?
            }
            Err(error) => return Err(error),
        };
        let mut elastic_free_norm = free.iter().fold(0.0_f64, |n, &d| n.hypot(f_int[d]));
        let mut static_force_scale = elastic_free_norm.max(external_free_norm);
        let mut increment_reference = 0.0_f64;
        // The inertial gate's scale is recomputed from each accepted iterate.
        // Like the static scale it never falls within an increment: an
        // iterate that approaches rest must not shrink its own reference.
        let mut inertial_force_scale = 0.0_f64;
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
                inertial_force_scale = inertial_force_scale.max(inertia.force_scale(
                    &u,
                    &free,
                    elastic_free_norm,
                    external_free_norm,
                    cx,
                )?);
                inertial_force_scale
            } else {
                static_force_scale.max(f64::MIN_POSITIVE)
            };
            let before = norm(&residual);
            let relative = before / scale;
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
            let mut connected_scale = 0.0_f64;
            for (row, &dof) in free.iter().enumerate() {
                let (cols, vals) = full.row(dof);
                let mut deformation = [0.0_f64; 6];
                let mut lifting = 0.0_f64;
                for (&c, &v) in cols.iter().zip(vals) {
                    if is_prescribed[c] {
                        rhs[row] = (-v).mul_add(dp_full[c], rhs[row]);
                        lifting = v.mul_add(dp_full[c], lifting);
                    } else {
                        coo.push(row, row_of[c], v);
                    }
                    if inertia.is_none() && v != 0.0 {
                        let component = c % components;
                        let difference = operator.deformation_coordinate(&u, dof, c);
                        deformation[component] = v.mul_add(difference, deformation[component]);
                    }
                }
                if inertia.is_none() {
                    // Contract each displacement component before taking its
                    // magnitude. This retains cancellation of affine fields
                    // against B-bar's mean-zero correction; a translation has
                    // no deformation scale. Opposing x/y/z elastic forces and
                    // prescribed lifts still provide a physical reference for
                    // a nearly equilibrated state with zero applied load.
                    let row_scale = deformation.iter().map(|v| v.abs()).sum::<f64>();
                    if !row_scale.is_finite() || !lifting.is_finite() {
                        return Err(SmallStrainTetError::InvalidInput {
                            what: "connected equilibrium force scale is not representable".into(),
                        });
                    }
                    connected_scale = connected_scale.hypot(row_scale.max(lifting.abs()));
                }
            }
            if inertia.is_none() {
                if !connected_scale.is_finite() {
                    return Err(SmallStrainTetError::InvalidInput {
                        what: "connected equilibrium force norm is not representable".into(),
                    });
                }
                // A large tangent by itself (e.g. an unused bulk penalty) must
                // not weaken the actual-force gate. Conversely, a large stress
                // unrelated to this equation cannot provide the reference.
                // The prescribed lift is this increment's physical force
                // reference. Once the targets are reached the remaining lift is
                // zero, and a homogeneous field cancels the connected
                // deformation scale, so the largest reference of the increment
                // is kept; the accepted stress bound still caps it.
                increment_reference = increment_reference.max(connected_scale);
                static_force_scale =
                    static_force_scale.max(increment_reference.min(stress_force_cap));
                let relative = before / static_force_scale.max(f64::MIN_POSITIVE);
                *history.last_mut().expect("current residual was recorded") = relative;
                if !pending_any && relative <= settings.relative_tolerance {
                    converged = true;
                    break;
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
                        // A shorter trial can repair a Newton proposal that left
                        // the law's domain. With every coordinate prescribed there
                        // is no proposal: the trial is the targets, the endpoint of
                        // the increment, and the accepted history is frozen until
                        // commit, so a shorter lift can only be refused there again
                        // or stall. Its refusal is the answer and propagates.
                        Err(error)
                            if error != SmallStrainTetError::Cancelled
                                && !free.is_empty()
                                && operator.recoverable_trial_error(&error) =>
                        {
                            alpha *= 0.5;
                            continue;
                        }
                        Err(error) => return Err(error),
                    };
                let trial_elastic_free_norm =
                    free.iter().fold(0.0_f64, |n, &d| n.hypot(f_trial[d]));
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
            stress_force_cap = trial_scale;
            // Rejected proposals cannot loosen the force gate. Retain the
            // largest accepted physical scale through unloading.
            static_force_scale = static_force_scale.max(trial_elastic_free_norm);
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
