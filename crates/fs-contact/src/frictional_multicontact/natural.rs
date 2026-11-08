//! Analytic semismooth action of the non-associated Coulomb natural map.
//!
//! For a sliding projection P(y,r)=r*u, u=y/|y|, its differential is
//! dP = u*dr + (r/|y|)*(I-u*u^T)*dy. Inside the disk, dP=dy.
//! These two tangent axes always form one Euclidean disk. At its boundary the
//! exterior derivative is selected; at the collapsed disk's apex the zero
//! projection derivative is selected. The normal positive-part map similarly
//! uses one fixed branch per Newton point. Every action is LINEAR in direction,
//! including at a kink; a direction-dependent difference stencil is not used.
//!
//! Local branch blocks are cached at the exact numerical iterate. Each Krylov
//! action then needs only W*direction and O(N) block work. This is a primal
//! semismooth linearization, not a classical contact-mode adjoint.
//! Large SPD groups reuse their admitted mobility factor as the right inverse
//! W^-1 S^-1, where S contains the natural-map row steps. It is exact on a fully
//! loaded sticking branch (J=S*W). Separating normal branches retain identity
//! preconditioning; their direct identity rows do not need a mobility inverse.
//!
//! For PSD contact mobility, a cached numerical range inverse of this exact
//! Jacobian keeps Newton corrections out of its unresolved impulse gauge.
//! Generic SVD stays in fs-la. This is only a RIGHT preconditioner: the original
//! Jacobian, residual and all physical gates remain unchanged. A discarded
//! linear direction cannot be hidden from FGMRES's true residual check. The
//! fixed 63-row cap bounds factorization work between cancellation polls.

use super::FrictionalContactError as E;
use crate::multicontact::add_product;
use fs_math::dd::Dd;
use fs_math::eft::{two_prod, two_sum};
use fs_solver::NonlinearProblem;
use fs_solver::orthant::active_set::{PreparedSpdMatrix, SpdOrthantError};
use std::cell::{Cell, RefCell};

#[derive(Clone, Copy)]
pub(super) enum Preconditioner<'a> {
    Identity,
    Mobility(&'a PreparedSpdMatrix),
    Range,
}

struct Block {
    direct: [[f64; 3]; 3],
    velocity: [[f64; 3]; 3],
}

struct Linearization {
    point: Vec<f64>,
    residual: Vec<f64>,
    blocks: Vec<Block>,
    inverse: Option<RangeInverse>,
}

struct RangeInverse {
    scale: f64,
    factor: fs_la::factor::Svd,
    cutoff: f64,
}

pub(super) struct NaturalMap<'a, C> {
    matrix: &'a [f64],
    free: &'a [f64],
    targets: &'a [f64],
    steps: &'a [f64],
    mu: &'a [f64],
    evaluations: &'a Cell<usize>,
    maximum_evaluations: usize,
    preconditioner: Preconditioner<'a>,
    cancelled: C,
    failed: Cell<bool>,
    failure: RefCell<Option<E>>,
    cached: RefCell<Option<Linearization>>,
}

impl<'a, C: Fn() -> bool> NaturalMap<'a, C> {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn new(
        matrix: &'a [f64],
        free: &'a [f64],
        targets: &'a [f64],
        steps: &'a [f64],
        mu: &'a [f64],
        evaluations: &'a Cell<usize>,
        maximum_evaluations: usize,
        preconditioner: Preconditioner<'a>,
        cancelled: C,
    ) -> Self {
        Self {
            matrix,
            free,
            targets,
            steps,
            mu,
            evaluations,
            maximum_evaluations,
            preconditioner,
            cancelled,
            failed: Cell::new(false),
            failure: RefCell::new(None),
            cached: RefCell::new(None),
        }
    }

    pub(super) fn has_failed(&self) -> bool {
        self.failed.get()
    }

    pub(super) fn take_failure(&self) -> Option<E> {
        self.failure.borrow_mut().take()
    }

    fn fail(&self, error: E) {
        if !self.failed.replace(true) {
            *self.failure.borrow_mut() = Some(error);
        }
    }

    fn poll(&self) -> bool {
        if self.has_failed() {
            false
        } else if (self.cancelled)() {
            self.fail(E::Cancelled);
            false
        } else {
            true
        }
    }

    fn begin(&self, point: &[f64], output_len: usize) -> bool {
        if !self.poll() {
            return false;
        }
        if point.len() != self.free.len()
            || output_len != self.free.len()
            || point.iter().any(|v| !v.is_finite())
        {
            self.fail(E::Invalid(
                "finite complete natural-map coordinates required",
            ));
            return false;
        }
        if self.evaluations.get() == self.maximum_evaluations {
            self.fail(E::EvaluationBudget);
            return false;
        }
        self.evaluations.set(self.evaluations.get() + 1);
        true
    }

    fn linearize(&self, point: &[f64]) -> bool {
        if self
            .cached
            .borrow()
            .as_ref()
            .is_some_and(|c| c.point == point)
        {
            return true;
        }
        let d = point.len();
        let mut residual = vec![0.0; d];
        let mut blocks = Vec::with_capacity(self.mu.len());
        for k in 0..self.mu.len() {
            if !self.poll() {
                return false;
            }
            let a = 3 * k;
            let mut v = [0.0; 3];
            for b in 0..3 {
                v[b] = affine_row(
                    self.free[a + b],
                    &self.matrix[(a + b) * d..(a + b + 1) * d],
                    point,
                );
            }
            let normal_trial = point[a] - self.steps[a] * (v[0] - self.targets[k]);
            let trial = [
                point[a + 1] - self.steps[a + 1] * v[1],
                point[a + 2] - self.steps[a + 2] * v[2],
            ];
            let radius = self.mu[k] * point[a].max(0.0);
            let length = trial[0].hypot(trial[1]);
            if v.iter()
                .chain(&trial)
                .chain(&[normal_trial, radius, length])
                .any(|x| !x.is_finite())
            {
                self.fail(E::NonFinite);
                return false;
            }
            let mut direct = [[0.0; 3]; 3];
            let mut velocity = [[0.0; 3]; 3];
            if normal_trial >= 0.0 {
                // Equivalent to p - max(p-r*v,0), without subtracting nearly
                // equal impulses at a loaded equilibrium.
                residual[a] = self.steps[a] * (v[0] - self.targets[k]);
                velocity[0][0] = self.steps[a];
            } else {
                residual[a] = point[a];
                direct[0][0] = 1.0;
            }
            if length < radius {
                for b in 1..3 {
                    residual[a + b] = self.steps[a + b] * v[b];
                    velocity[b][b] = self.steps[a + b];
                }
            } else if length > 0.0 {
                let unit = [trial[0] / length, trial[1] / length];
                let fraction = radius / length;
                let dr = if point[a] >= 0.0 { self.mu[k] } else { 0.0 };
                for i in 0..2 {
                    residual[a + 1 + i] = point[a + 1 + i] - radius * unit[i];
                    direct[i + 1][0] = -dr * unit[i];
                    for j in 0..2 {
                        let identity = if i == j { 1.0 } else { 0.0 };
                        let projection = fraction * (identity - unit[i] * unit[j]);
                        direct[i + 1][j + 1] = identity - projection;
                        velocity[i + 1][j + 1] = projection * self.steps[a + j + 1];
                    }
                }
            } else {
                // Collapsed disk at its apex: one fixed generalized derivative,
                // not a direction-dependent projection of the Krylov vector.
                for b in 1..3 {
                    residual[a + b] = point[a + b];
                    direct[b][b] = 1.0;
                }
            }
            if direct
                .iter()
                .chain(&velocity)
                .flatten()
                .any(|v| !v.is_finite())
            {
                self.fail(E::NonFinite);
                return false;
            }
            blocks.push(Block { direct, velocity });
        }
        if residual.iter().any(|v| !v.is_finite()) {
            self.fail(E::NonFinite);
            return false;
        }
        *self.cached.borrow_mut() = Some(Linearization {
            point: point.to_vec(),
            residual,
            blocks,
            inverse: None,
        });
        true
    }

    fn prepare_inverse(&self, point: &[f64]) -> bool {
        if self.cached.borrow().as_ref().unwrap().inverse.is_some() {
            return true;
        }
        let d = point.len();
        let mut jacobian = vec![0.0; d * d];
        {
            let cached = self.cached.borrow();
            for (k, block) in cached.as_ref().unwrap().blocks.iter().enumerate() {
                if !self.poll() {
                    return false;
                }
                let a = 3 * k;
                for i in 0..3 {
                    for j in 0..d {
                        let mut value = if (a..a + 3).contains(&j) {
                            block.direct[i][j - a]
                        } else {
                            0.0
                        };
                        for b in 0..3 {
                            value += block.velocity[i][b] * self.matrix[(a + b) * d + j];
                        }
                        jacobian[(a + i) * d + j] = value;
                    }
                }
            }
        }
        let scale = jacobian.iter().map(|v| v.abs()).fold(0.0, f64::max);
        if !scale.is_finite() || scale <= 0.0 || jacobian.iter().any(|v| !v.is_finite()) {
            self.fail(E::NonFinite);
            return false;
        }
        for value in &mut jacobian {
            *value /= scale;
        }
        if !self.poll() {
            return false;
        }
        let factor = fs_la::factor::svd_jacobi(&jacobian, d, d);
        if !self.poll() {
            return false;
        }
        if factor
            .u
            .iter()
            .chain(&factor.v)
            .chain(&factor.sigma)
            .any(|v| !v.is_finite())
        {
            self.fail(E::NonFinite);
            return false;
        }
        // Numerical linear-solve rank only; this is not a contact-rank claim or
        // a compliance. Any omitted RHS component still faces the unchanged
        // outer Krylov residual test and the original nonlinear/physical gates.
        let cutoff = 128.0 * f64::EPSILON * d as f64 * factor.sigma[0];
        self.cached.borrow_mut().as_mut().unwrap().inverse = Some(RangeInverse {
            scale,
            factor,
            cutoff,
        });
        true
    }
}

impl<C: Fn() -> bool> NonlinearProblem for NaturalMap<'_, C> {
    fn dimension(&self) -> usize {
        self.free.len()
    }

    fn residual(&self, point: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if self.begin(point, output.len()) && self.linearize(point) {
            output.copy_from_slice(&self.cached.borrow().as_ref().unwrap().residual);
        }
    }

    fn jacobian_apply(&self, point: &[f64], direction: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if !self.begin(point, output.len()) {
            return;
        }
        if direction.len() != point.len() || direction.iter().any(|v| !v.is_finite()) {
            self.fail(E::Invalid("finite complete natural-map direction required"));
            return;
        }
        if !self.linearize(point) {
            return;
        }
        let cached = self.cached.borrow();
        let linearization = cached.as_ref().unwrap();
        let d = point.len();
        for (k, block) in linearization.blocks.iter().enumerate() {
            if !self.poll() {
                output.fill(f64::NAN);
                return;
            }
            let a = 3 * k;
            let mut dv = [0.0; 3];
            for b in 0..3 {
                dv[b] = affine_row(0.0, &self.matrix[(a + b) * d..(a + b + 1) * d], direction);
            }
            let local_direction = [
                direction[a],
                dv[0],
                direction[a + 1],
                dv[1],
                direction[a + 2],
                dv[2],
            ];
            for i in 0..3 {
                let coefficients = [
                    block.direct[i][0],
                    block.velocity[i][0],
                    block.direct[i][1],
                    block.velocity[i][1],
                    block.direct[i][2],
                    block.velocity[i][2],
                ];
                output[a + i] = affine_row(0.0, &coefficients, &local_direction);
            }
        }
        if output.iter().any(|v| !v.is_finite()) {
            self.fail(E::NonFinite);
            output.fill(f64::NAN);
        }
    }

    fn preconditioner_apply(
        &self,
        point: &[f64],
        _outer_iteration: usize,
        _inner_iteration: usize,
        residual: &[f64],
        output: &mut [f64],
    ) {
        output.fill(f64::NAN);
        if !self.begin(point, output.len()) {
            return;
        }
        if residual.len() != point.len() || residual.iter().any(|v| !v.is_finite()) {
            self.fail(E::Invalid(
                "finite complete right-preconditioner input required",
            ));
            return;
        }
        match self.preconditioner {
            Preconditioner::Identity => {
                output.copy_from_slice(residual);
                return;
            }
            Preconditioner::Mobility(factor) => {
                if !self.linearize(point) {
                    return;
                }
                if self
                    .cached
                    .borrow()
                    .as_ref()
                    .unwrap()
                    .blocks
                    .iter()
                    .any(|block| block.direct[0][0] != 0.0)
                {
                    output.copy_from_slice(residual);
                    return;
                }
                let rhs: Vec<_> = residual
                    .iter()
                    .zip(self.steps)
                    .map(|(r, s)| r / s)
                    .collect();
                if rhs.iter().any(|v| !v.is_finite()) {
                    self.fail(E::NonFinite);
                    return;
                }
                match factor.solve(&rhs, || !self.poll()) {
                    Ok(solved) => output.copy_from_slice(&solved),
                    Err(SpdOrthantError::Cancelled) => self.fail(E::Cancelled),
                    Err(SpdOrthantError::NonFinite) => self.fail(E::NonFinite),
                    Err(error) => self.fail(E::MobilityNotAdmitted(error)),
                }
                return;
            }
            Preconditioner::Range => {}
        }
        if !self.linearize(point) || !self.prepare_inverse(point) {
            return;
        }
        let cached = self.cached.borrow();
        let inverse = cached.as_ref().unwrap().inverse.as_ref().unwrap();
        let d = point.len();
        output.fill(0.0);
        for k in 0..d {
            if !self.poll() {
                output.fill(f64::NAN);
                return;
            }
            let sigma = inverse.factor.sigma[k];
            if sigma <= inverse.cutoff {
                continue;
            }
            let mut coefficient = 0.0;
            for (i, &r) in residual.iter().enumerate() {
                coefficient += inverse.factor.u[i * d + k] * (r / inverse.scale);
            }
            coefficient /= sigma;
            for (i, value) in output.iter_mut().enumerate() {
                *value += inverse.factor.v[i * d + k] * coefficient;
            }
        }
        if output.iter().any(|v| !v.is_finite()) {
            self.fail(E::NonFinite);
            output.fill(f64::NAN);
        }
    }
}

// The affine offset must participate in the same accumulation as W*p. Rounding
// the product sum first can hide actual slip and manufacture a sticking root,
// even when the exact stopping impulse is representable. TwoProd/TwoSum retain
// the first-order product and summation roundoff in a scalar correction; the
// remaining accumulation error is second order. This keeps the hot row free of
// both a products allocation and full double-double arithmetic per term. If an
// intermediate overflows, reuse the physical audit's scaled product recovery.
fn affine_row(offset: f64, row: &[f64], point: &[f64]) -> f64 {
    let mut high = offset;
    let mut low = 0.0;
    for (&a, &b) in row.iter().zip(point) {
        let (product, product_error) = two_prod(a, b);
        let (next, sum_error) = two_sum(high, product);
        high = next;
        low += sum_error + product_error;
    }
    let result = high + low;
    if result.is_finite() {
        result
    } else {
        row.iter()
            .zip(point)
            .fold(Dd::from_f64(offset), |sum, (&a, &b)| add_product(sum, a, b))
            .to_f64()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compensated_affine_row_preserves_group_cancellation_and_overflow_recovery() {
        let epsilon = 2.0_f64.powi(-27);
        for dimension in [3, 12, 63] {
            let row = vec![1.0 + epsilon; dimension];
            let point = vec![1.0 - epsilon; dimension];
            // Every rounded product equals one, but their actual sum is
            // dimension*(1-e^2). The residual is exactly representable.
            let expected = -(dimension as f64) * epsilon * epsilon;
            assert_eq!(affine_row(-(dimension as f64), &row, &point), expected);
        }
        assert_eq!(affine_row(-f64::MAX, &[f64::MAX], &[2.0]), f64::MAX);
    }

    #[test]
    fn coupled_sticking_iteration_resolves_slip_hidden_by_affine_cancellation() {
        use fs_solver::{NewtonKrylovConfig, NewtonKrylovState};

        // This well-conditioned SPD mobility couples the two tangent axes.
        // The loaded contact has ample static-friction capacity. Its exact
        // stopping impulse is representable, even though forming W*p before
        // adding free motion rounds away the residual tangential velocity.
        let scale = 2.0_f64.powi(54);
        let w = [1.0, 0.0, 0.0, 0.0, 1.0, 0.5, 0.0, 0.5, 1.0];
        let free = [-2.0 * scale, scale, 2.0 * scale];
        let count = Cell::new(0);
        let problem = NaturalMap::new(
            &w,
            &free,
            &[0.0],
            &[1.0; 3],
            &[2.0],
            &count,
            256,
            Preconditioner::Identity,
            || false,
        );
        let seed = vec![2.0 * scale, 1.0, -2.0 * scale];
        let mut residual = [0.0; 3];
        problem.residual(&seed, &mut residual);
        assert_eq!(residual, [0.0, 1.0, 0.5]);

        let config = NewtonKrylovConfig {
            relative_tolerance: 0.0,
            forcing_minimum: 1e-12,
            forcing_maximum: 1e-12,
            ..Default::default()
        };
        let mut state = NewtonKrylovState::new(&problem, seed, config).unwrap();
        let report = state.run(&problem, 4);
        assert!(report.converged, "{report:?}");
        assert!(!problem.has_failed());
        assert_eq!(state.x[0], 2.0 * scale);
        assert!(state.x[1].abs() < 1e-12, "remaining slip: {:?}", state.x);
        assert_eq!(state.x[2], -2.0 * scale);
        problem.residual(&state.x, &mut residual);
        assert!(residual[1].abs() < 1e-12);
        assert!(residual[2].abs() < 1e-12);
    }

    #[test]
    fn sticking_jacobian_retains_product_roundoff_in_coupled_tangent_motion() {
        let epsilon = 2.0_f64.powi(-27);
        let w = [1.0, 0.0, 0.0, 0.0, 1.0 + epsilon, 1.0, 0.0, 1.0, 2.0];
        let count = Cell::new(0);
        let problem = NaturalMap::new(
            &w,
            &[-2.0, 0.0, 0.0],
            &[0.0],
            &[1.0; 3],
            &[2.0],
            &count,
            16,
            Preconditioner::Identity,
            || false,
        );
        let mut action = [0.0; 3];
        problem.jacobian_apply(&[2.0, 0.0, 0.0], &[0.0, 1.0 - epsilon, -1.0], &mut action);
        // (1+e)*(1-e)-1 = -e^2 exactly. A rounded products vector
        // incorrectly makes this mobility direction a zero-slip direction.
        assert_eq!(action, [0.0, -epsilon * epsilon, -1.0 - epsilon]);
        assert!(!problem.has_failed());
    }

    #[test]
    fn smooth_branch_actions_match_independent_central_differences() {
        let w = [1.0, 0.3, 0.2, 0.3, 1.5, -0.1, 0.2, -0.1, 2.0];
        let p = [2.0, 0.2, -0.1];
        let direction = [0.3, -0.4, 0.7];
        for endpoint in [[0.0, 0.02, -0.01], [0.0, 2.0, 3.0], [3.0, 2.0, -1.0]] {
            let mut free = endpoint;
            for i in 0..3 {
                for j in 0..3 {
                    free[i] -= w[3 * i + j] * p[j];
                }
            }
            let count = Cell::new(0);
            let map = NaturalMap::new(
                &w,
                &free,
                &[0.0],
                &[1.0, 0.6, 0.6],
                &[0.7],
                &count,
                32,
                Preconditioner::Identity,
                || false,
            );
            let mut action = [0.0; 3];
            map.jacobian_apply(&p, &direction, &mut action);
            let h = 1e-5;
            let plus = std::array::from_fn::<_, 3, _>(|i| p[i] + h * direction[i]);
            let minus = std::array::from_fn::<_, 3, _>(|i| p[i] - h * direction[i]);
            let mut r_plus = [0.0; 3];
            let mut r_minus = [0.0; 3];
            map.residual(&plus, &mut r_plus);
            map.residual(&minus, &mut r_minus);
            assert!(!map.has_failed());
            for i in 0..3 {
                let difference = (r_plus[i] - r_minus[i]) / (2.0 * h);
                assert!(
                    (action[i] - difference).abs() < 2e-10,
                    "action {} versus independent difference {}",
                    action[i],
                    difference
                );
            }
        }
    }

    #[test]
    fn boundary_action_is_linear_and_does_not_perturb_the_physical_iterate() {
        let w = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let count = Cell::new(0);
        let map = NaturalMap::new(
            &w,
            &[-2.0, -1.0, 0.0],
            &[0.0],
            &[1.0; 3],
            &[0.5],
            &count,
            32,
            Preconditioner::Identity,
            || false,
        );
        let p = [2.0, 1.0, 0.0];
        let a = [0.1, -0.3, 0.4];
        let b = [-0.7, 0.5, 0.2];
        let sum = std::array::from_fn::<_, 3, _>(|i| a[i] + b[i]);
        let mut ja = [0.0; 3];
        let mut jb = [0.0; 3];
        let mut js = [0.0; 3];
        map.jacobian_apply(&p, &a, &mut ja);
        map.jacobian_apply(&p, &b, &mut jb);
        map.jacobian_apply(&p, &sum, &mut js);
        for i in 0..3 {
            assert!((js[i] - ja[i] - jb[i]).abs() < 1e-14);
        }
        assert_eq!(count.get(), 3);
        assert_eq!(map.cached.borrow().as_ref().unwrap().point, p);
        let mut residual = [1.0; 3];
        map.residual(&p, &mut residual);
        assert_eq!(residual, [0.0; 3]);
    }

    #[test]
    fn range_preconditioner_leaves_an_unfittable_rhs_visible_to_the_original_action() {
        let w = [1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0];
        let count = Cell::new(0);
        let map = NaturalMap::new(
            &w,
            &[-1.0, -0.2, 0.0],
            &[0.0],
            &[1.0, 0.5, 0.5],
            &[0.5],
            &count,
            32,
            Preconditioner::Range,
            || false,
        );
        let point = [1.0, 0.1, 0.0];
        let rhs = [2.0, 3.0, 4.0];
        let mut correction = [0.0; 3];
        map.preconditioner_apply(&point, 0, 0, &rhs, &mut correction);
        assert_eq!(correction, [2.0, 3.0, 0.0]);
        let mut actual = [0.0; 3];
        map.jacobian_apply(&point, &correction, &mut actual);
        assert_eq!(actual, [2.0, 3.0, 0.0]);
        assert_eq!(rhs[2] - actual[2], 4.0);
        assert!(!map.has_failed());
    }

    #[test]
    fn prepared_mobility_removes_loaded_chain_krylov_conditioning() {
        use fs_solver::orthant::active_set::{SpdOrthantConfig, prepare_spd_matrix};
        use fs_solver::{NewtonKrylovConfig, NewtonKrylovState};

        let n = 64;
        let d = 3 * n;
        let mut w = vec![0.0; d * d];
        let mut expected = vec![0.0; d];
        let mut free = vec![0.0; d];
        let mut steps = vec![1.0 / 7.0; d];
        for i in 0..n {
            // Unit-mass solid spheres: translational normal mobility and both
            // translation/spin shear mobility, including the shared neighbor.
            for (axis, diagonal, adjacent) in [(0, 2.0, -1.0), (1, 7.0, 1.5), (2, 7.0, 1.5)] {
                let row = 3 * i + axis;
                w[row * d + row] = diagonal;
                if i + 1 < n {
                    w[row * d + row + 3] = adjacent;
                    w[(row + 3) * d + row] = adjacent;
                }
            }
            // The discrete second difference of these loads is 8/n^2. A
            // uniform approaching rate therefore stops at this exact profile.
            expected[3 * i] = 4.0 * ((i + 1) * (n - i)) as f64 / (n * n) as f64;
            free[3 * i] = -8.0 / (n * n) as f64;
            steps[3 * i] = 0.5;
        }
        let factor = prepare_spd_matrix(
            &w,
            d,
            SpdOrthantConfig {
                max_variables: d,
                ..Default::default()
            },
            || false,
        )
        .unwrap();
        let targets = vec![0.0; n];
        let mu = vec![0.3; n];
        let config = NewtonKrylovConfig {
            relative_tolerance: 0.0,
            linear_restart: 64,
            max_linear_cycles: 1,
            forcing_minimum: 1e-10,
            forcing_maximum: 1e-10,
            ..Default::default()
        };
        let mut work = Vec::new();
        for preconditioner in [Preconditioner::Identity, Preconditioner::Mobility(&factor)] {
            let count = Cell::new(0);
            let problem = NaturalMap::new(
                &w,
                &free,
                &targets,
                &steps,
                &mu,
                &count,
                4096,
                preconditioner,
                || false,
            );
            // Both solves start strictly inside the same loaded/sticking
            // branch, so this compares preconditioning, not branch selection.
            let seed = expected.iter().map(|value| 0.5 * value).collect();
            let mut state = NewtonKrylovState::new(&problem, seed, config).unwrap();
            let report = state.run(&problem, 4);
            assert!(report.converged, "{report:?}");
            assert!(!problem.has_failed());
            for (actual, expected) in state.x.iter().zip(&expected) {
                assert!((actual - expected).abs() < 2e-10);
            }
            // Independent local momentum and midpoint work, using actual
            // solved loads rather than the factor or Newton residual.
            let mut actual_work = 0.0;
            let mut expected_work = 0.0;
            for i in 0..n {
                let pn = state.x[3 * i];
                let after = free[3 * i] + 2.0 * pn
                    - if i > 0 { state.x[3 * (i - 1)] } else { 0.0 }
                    - if i + 1 < n { state.x[3 * (i + 1)] } else { 0.0 };
                assert!(after.abs() < 2e-12);
                actual_work += 0.5 * pn * (free[3 * i] + after);
                expected_work += 0.5 * expected[3 * i] * free[3 * i];
            }
            assert!((actual_work - expected_work).abs() < 2e-11);
            assert!(actual_work < 0.0);
            work.push((
                report
                    .history
                    .iter()
                    .map(|step| step.linear_iterations)
                    .sum::<usize>(),
                count.get(),
            ));
        }
        assert!(work[0].0 >= 16, "identity work: {work:?}");
        assert!(work[1].0 <= 2, "prepared work: {work:?}");
        assert!(4 * work[1].1 < work[0].1, "callback reduction: {work:?}");
    }

    #[test]
    fn mobility_preconditioner_retains_separating_rows_and_cancellation() {
        use fs_solver::orthant::active_set::prepare_spd_matrix;

        let mut w = vec![0.0; 36];
        for (axis, diagonal, adjacent) in [(0, 2.0, -1.0), (1, 7.0, 1.5), (2, 7.0, 1.5)] {
            w[axis * 6 + axis] = diagonal;
            w[(axis + 3) * 6 + axis + 3] = diagonal;
            w[axis * 6 + axis + 3] = adjacent;
            w[(axis + 3) * 6 + axis] = adjacent;
        }
        let factor = prepare_spd_matrix(&w, 6, Default::default(), || false).unwrap();
        let point = [1.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let steps = [0.5, 1.0 / 7.0, 1.0 / 7.0, 0.5, 1.0 / 7.0, 1.0 / 7.0];
        let rhs = [1.0, -2.0, 3.0, 4.0, -5.0, 6.0];
        let count = Cell::new(0);
        let separating = NaturalMap::new(
            &w,
            &[-1.0, 0.0, 0.0, 3.0, 0.0, 0.0],
            &[0.0; 2],
            &steps,
            &[0.3; 2],
            &count,
            32,
            Preconditioner::Mobility(&factor),
            || false,
        );
        let mut output = [0.0; 6];
        separating.preconditioner_apply(&point, 0, 0, &rhs, &mut output);
        assert_eq!(output, rhs);
        assert!(!separating.has_failed());

        let free = [-1.0, 0.0, 0.0, -1.0, 0.0, 0.0];
        let polls = Cell::new(0);
        let loaded = NaturalMap::new(
            &w,
            &free,
            &[0.0; 2],
            &steps,
            &[0.3; 2],
            &count,
            32,
            Preconditioner::Mobility(&factor),
            || {
                polls.set(polls.get() + 1);
                false
            },
        );
        loaded.preconditioner_apply(&point, 0, 0, &rhs, &mut output);
        let total = polls.get();
        let mut action = [0.0; 6];
        loaded.jacobian_apply(&point, &output, &mut action);
        for i in 0..6 {
            assert!((action[i] - rhs[i]).abs() < 2e-13);
        }
        for stop in [1, total / 2, total] {
            let seen = Cell::new(0);
            let map = NaturalMap::new(
                &w,
                &free,
                &[0.0; 2],
                &steps,
                &[0.3; 2],
                &count,
                32,
                Preconditioner::Mobility(&factor),
                || {
                    seen.set(seen.get() + 1);
                    seen.get() == stop
                },
            );
            map.preconditioner_apply(&point, 0, 0, &rhs, &mut output);
            assert!(matches!(map.take_failure(), Some(E::Cancelled)));
            assert!(output.iter().all(|value| value.is_nan()));
        }
    }
}
