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
//!
//! For PSD contact mobility, a cached numerical range inverse of this exact
//! Jacobian keeps Newton corrections out of its unresolved impulse gauge.
//! Generic SVD stays in fs-la. This is only a RIGHT preconditioner: the original
//! Jacobian, residual and all physical gates remain unchanged. A discarded
//! linear direction cannot be hidden from FGMRES's true residual check. The
//! fixed 63-row cap bounds factorization work between cancellation polls.

use super::FrictionalContactError as E;
use fs_solver::NonlinearProblem;
use std::cell::{Cell, RefCell};

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
    range_preconditioner: bool,
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
        range_preconditioner: bool,
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
            range_preconditioner,
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
                v[b] = self.free[a + b]
                    + fs_solver::dot(&self.matrix[(a + b) * d..(a + b + 1) * d], point);
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
                dv[b] = fs_solver::dot(&self.matrix[(a + b) * d..(a + b + 1) * d], direction);
            }
            for i in 0..3 {
                output[a + i] = 0.0;
                for j in 0..3 {
                    output[a + i] +=
                        block.direct[i][j] * direction[a + j] + block.velocity[i][j] * dv[j];
                }
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
        if !self.range_preconditioner {
            output.copy_from_slice(residual);
            return;
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

#[cfg(test)]
mod tests {
    use super::*;

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
                false,
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
            false,
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
            true,
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
}
