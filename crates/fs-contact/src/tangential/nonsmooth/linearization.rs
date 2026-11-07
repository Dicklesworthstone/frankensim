//! Fixed-branch tangent of the original two-region elastic/Coulomb return map.
//!
//! For one yielding region, e = K (u-z), z = lambda C^-2 e and
//! |C^-1 e| = 1. Implicit differentiation gives
//! D = H - (H b)(H b)^T / (b^T H b), where
//! H = (K^-1 + lambda C^-2)^-1 and b = C^-2 e.
//! The elastic region instead has D=K. The published plastic increment
//! recovers lambda from the original return, without evaluating a second law.
//!
//! The cached blocks are -dp/dv = dt^2 S^-1 D S^-1 for
//! S=diag(1,1,ell). Every Krylov action at the same Newton point is linear in
//! its direction, including the fixed elastic selection at a yield boundary.
//! A static-to-kinetic drop is still discontinuous; this local branch tangent
//! is neither a smoothing nor a classical derivative across that transition.

use super::{NonsmoothTangentialError as E, poll};
use fs_exec::CancelGate;
use fs_solver::NonlinearProblem;
use fs_tribo::partial_slip::{PartialSlipLaw, PartialSlipStep};
use std::cell::{Cell, RefCell};

type Block = [[f64; 3]; 3];

pub(super) fn impulse_tangent(
    law: &PartialSlipLaw,
    step: &PartialSlipStep,
    ell: f64,
    dt: f64,
) -> Result<Block, E> {
    let p = law.parameters();
    let hard = p.partial_slip_hardening_fraction;
    let fractions = [
        p.partial_slip_onset_fraction * (1.0 - hard),
        1.0 - p.partial_slip_onset_fraction * (1.0 - hard),
    ];
    let mut block = [[0.0; 3]; 3];
    for (region, weight) in [1.0 - hard, hard].into_iter().enumerate() {
        let stiffness = [
            p.tangential_stiffness_n_per_m * weight,
            p.tangential_stiffness_n_per_m * weight,
            p.torsional_stiffness_nm_per_rad * weight,
        ];
        let plastic = step.region_plastic_increment[region];
        let mut diagonal = stiffness;
        let mut normal = [0.0; 3];
        if plastic.iter().any(|&v| v != 0.0) {
            let capacity = [
                step.capacity.kinetic_force_n,
                step.capacity.kinetic_force_n,
                step.capacity.kinetic_torque_nm,
            ]
            .map(|v| v * fractions[region]);
            let elastic = [
                step.next_state.region_elastic_displacement_m()[region][0],
                step.next_state.region_elastic_displacement_m()[region][1],
                step.next_state.region_elastic_twist_rad()[region],
            ];
            let mut q = [0.0; 3];
            let mut numerator = 0.0;
            let mut denominator = 0.0;
            for j in 0..3 {
                q[j] = stiffness[j] * elastic[j] / capacity[j];
                numerator += q[j] * (plastic[j] * capacity[j]);
                denominator += q[j] * q[j];
            }
            let lambda = numerator / denominator;
            if !lambda.is_finite() || lambda < 0.0 {
                return Err(E::Invalid("unrepresentable regional return multiplier"));
            }
            for j in 0..3 {
                diagonal[j] =
                    stiffness[j] / (1.0 + (lambda / capacity[j]) * (stiffness[j] / capacity[j]));
                normal[j] = diagonal[j].sqrt() * q[j] / capacity[j];
            }
            let norm = normal[0].hypot(normal[1]).hypot(normal[2]);
            if !norm.is_finite() || norm <= 0.0 {
                return Err(E::Invalid("unrepresentable regional return tangent"));
            }
            normal = normal.map(|v| v / norm);
        }
        // Scale square roots, rather than forming dt^2 or ell^2 separately.
        let scales = [dt, dt, dt / ell];
        let root: [f64; 3] = std::array::from_fn(|j| scales[j] * diagonal[j].sqrt());
        for i in 0..3 {
            for j in 0..3 {
                let identity = if i == j { 1.0 } else { 0.0 };
                block[i][j] += root[i] * root[j] * (identity - normal[i] * normal[j]);
            }
        }
    }
    if block.iter().flatten().any(|v| !v.is_finite()) {
        return Err(E::Invalid("unrepresentable contact impulse tangent"));
    }
    Ok(block)
}

struct Linearization {
    point: Vec<f64>,
    residual: Vec<f64>,
    blocks: Vec<Block>,
}

struct WorkVectors {
    contact: Vec<f64>,
    mechanical: Vec<f64>,
}

/// Law evaluations occur only when the exact Newton point changes. A Krylov
/// action costs two mobility-factor products and one 3x3 block per contact;
/// it never allocates a dense contact matrix or perturbs retained history.
pub(super) struct PartialSlipProblem<'a, F> {
    dimension: usize,
    factor: &'a [f64],
    mechanical_dofs: usize,
    gate: &'a CancelGate,
    evaluate: F,
    cached: RefCell<Option<Linearization>>,
    work: RefCell<WorkVectors>,
    failed: Cell<bool>,
    failure: RefCell<Option<E>>,
}

impl<'a, F> PartialSlipProblem<'a, F>
where
    F: Fn(&[f64]) -> Result<(Vec<f64>, Vec<Block>), E>,
{
    pub(super) fn new(
        dimension: usize,
        factor: &'a [f64],
        mechanical_dofs: usize,
        gate: &'a CancelGate,
        evaluate: F,
    ) -> Self {
        Self {
            dimension,
            factor,
            mechanical_dofs,
            gate,
            evaluate,
            cached: RefCell::new(None),
            work: RefCell::new(WorkVectors {
                contact: vec![0.0; dimension],
                mechanical: vec![0.0; mechanical_dofs],
            }),
            failed: Cell::new(false),
            failure: RefCell::new(None),
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

    fn linearize(&self, point: &[f64], output_len: usize) -> Result<(), E> {
        poll(self.gate)?;
        if point.len() != self.dimension
            || output_len != self.dimension
            || point.iter().any(|v| !v.is_finite())
        {
            return Err(E::Invalid(
                "finite complete partial-slip coordinates required",
            ));
        }
        if self
            .cached
            .borrow()
            .as_ref()
            .is_some_and(|v| v.point == point)
        {
            return Ok(());
        }
        let (residual, blocks) = (self.evaluate)(point)?;
        if residual.len() != self.dimension
            || 3 * blocks.len() != self.dimension
            || residual.iter().any(|v| !v.is_finite())
            || blocks.iter().flatten().flatten().any(|v| !v.is_finite())
        {
            return Err(E::Invalid(
                "incomplete or nonfinite partial-slip linearization",
            ));
        }
        *self.cached.borrow_mut() = Some(Linearization {
            point: point.to_vec(),
            residual,
            blocks,
        });
        Ok(())
    }

    fn action(&self, point: &[f64], direction: &[f64], output: &mut [f64]) -> Result<(), E> {
        self.linearize(point, output.len())?;
        if direction.len() != self.dimension || direction.iter().any(|v| !v.is_finite()) {
            return Err(E::Invalid(
                "finite complete partial-slip direction required",
            ));
        }
        let cached = self.cached.borrow();
        let linearization = cached.as_ref().unwrap();
        let mut work = self.work.borrow_mut();
        work.contact.fill(0.0);
        for (contact, block) in linearization.blocks.iter().enumerate() {
            poll(self.gate)?;
            for i in 0..3 {
                for j in 0..3 {
                    work.contact[3 * contact + i] += block[i][j] * direction[3 * contact + j];
                }
            }
        }
        work.mechanical.fill(0.0);
        for row in 0..self.dimension {
            poll(self.gate)?;
            let value = work.contact[row];
            for col in 0..self.mechanical_dofs {
                work.mechanical[col] += self.factor[row * self.mechanical_dofs + col] * value;
            }
        }
        for row in 0..self.dimension {
            poll(self.gate)?;
            let change: f64 = (0..self.mechanical_dofs)
                .map(|col| self.factor[row * self.mechanical_dofs + col] * work.mechanical[col])
                .sum();
            output[row] = direction[row] + change;
        }
        if output.iter().any(|v| !v.is_finite()) {
            return Err(E::Invalid("unrepresentable partial-slip Jacobian action"));
        }
        Ok(())
    }
}

impl<F> NonlinearProblem for PartialSlipProblem<'_, F>
where
    F: Fn(&[f64]) -> Result<(Vec<f64>, Vec<Block>), E>,
{
    fn dimension(&self) -> usize {
        self.dimension
    }

    fn residual(&self, point: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if self.has_failed() {
            return;
        }
        match self.linearize(point, output.len()) {
            Ok(()) => output.copy_from_slice(&self.cached.borrow().as_ref().unwrap().residual),
            Err(e) => self.fail(e),
        }
    }

    fn jacobian_apply(&self, point: &[f64], direction: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if self.has_failed() {
            return;
        }
        if let Err(e) = self.action(point, direction, output) {
            output.fill(f64::NAN);
            self.fail(e);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_tribo::partial_slip::{
        GeneralizedWorkOwnership, NormalPatchAuthority, NormalPatchView, PARTIAL_SLIP_MODEL_ID,
        PartialSlipInterface, PartialSlipKinematics, PartialSlipParameters, PartialSlipState,
        PartialSlipStateKind, TangentFrame,
    };

    const DT: f64 = 0.001;
    const ELL: f64 = 0.01;

    fn law() -> PartialSlipLaw {
        PartialSlipLaw::new(
            PARTIAL_SLIP_MODEL_ID,
            "synthetic/regional-tangent",
            PartialSlipParameters {
                static_mu: 0.8,
                kinetic_mu: 0.4,
                tangential_stiffness_n_per_m: 10000.0,
                torsional_stiffness_nm_per_rad: 100.0,
                torsional_capacity_factor: 0.5,
                partial_slip_onset_fraction: 0.5,
                partial_slip_hardening_fraction: 0.4,
            },
        )
        .unwrap()
    }

    fn response(
        law: &PartialSlipLaw,
        old: &PartialSlipState,
        velocity: [f64; 3],
    ) -> PartialSlipStep {
        let patch = NormalPatchView::new(
            "patch",
            "positive-normal-fixture",
            "synthetic/normal",
            NormalPatchAuthority::SyntheticFixture,
            100.0,
            0.02,
            0.01,
            1e-4,
        )
        .unwrap();
        let interface = PartialSlipInterface::new(
            "a->b",
            "history",
            "synthetic/interface",
            NormalPatchAuthority::SyntheticFixture,
        )
        .unwrap();
        law.advance(
            &patch,
            &interface,
            TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
            PartialSlipKinematics {
                creepage: [velocity[0], velocity[1]],
                rolling_speed_mps: 1.0,
                torsional_spin_rad_per_s: velocity[2] / ELL,
                contact_spin_rad_per_s: 0.0,
                dt_s: DT,
            },
            &GeneralizedWorkOwnership::new("patch", "interval", "x", "y", "spin").unwrap(),
            old,
        )
        .unwrap()
    }

    fn impulse(step: &PartialSlipStep) -> [f64; 3] {
        [
            DT * step.tangent_force_n[0],
            DT * step.tangent_force_n[1],
            DT * step.torsional_moment_nm / ELL,
        ]
    }

    fn close(a: f64, b: f64, tolerance: f64) {
        assert!(
            (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1.0),
            "{a:e} != {b:e}"
        );
    }

    #[test]
    fn mixed_shear_spin_tangents_match_original_return_in_each_region_and_are_passive() {
        let law = law();
        let old = PartialSlipState::zero();
        for (velocity, kind) in [
            ([0.2, -0.1, 0.0001], PartialSlipStateKind::Sticking),
            ([7.0, 1.0, 0.015], PartialSlipStateKind::PartialSlip),
            ([30.0, 8.0, 0.2], PartialSlipStateKind::GrossSlide),
            ([30.0, 8.0, 0.0], PartialSlipStateKind::GrossSlide),
            ([0.0, 0.0, 0.2], PartialSlipStateKind::GrossSlide),
        ] {
            let step = response(&law, &old, velocity);
            assert_eq!(step.state, kind);
            let block = impulse_tangent(&law, &step, ELL, DT).unwrap();
            // Differentiate the original law only as an independent test
            // oracle. These points stay strictly inside their selected modes.
            for column in 0..3 {
                let h = 1e-6;
                let mut plus = velocity;
                let mut minus = velocity;
                plus[column] += h;
                minus[column] -= h;
                let upper = response(&law, &old, plus);
                let lower = response(&law, &old, minus);
                assert_eq!(upper.state, kind);
                assert_eq!(lower.state, kind);
                let upper = impulse(&upper);
                let lower = impulse(&lower);
                for row in 0..3 {
                    close(
                        block[row][column],
                        -(upper[row] - lower[row]) / (2.0 * h),
                        2e-7,
                    );
                    close(block[row][column], block[column][row], 1e-14);
                }
            }
            for direction in [
                [1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
                [1.0, -2.0, 0.3],
                [-0.2, 1.0, 2.0],
            ] {
                let mut quadratic = 0.0;
                for row in 0..3 {
                    for col in 0..3 {
                        quadratic += direction[row] * block[row][col] * direction[col];
                    }
                }
                assert!(
                    quadratic >= -1e-12,
                    "incremental friction is active: {quadratic:e}"
                );
            }
        }
    }

    #[test]
    fn a_paused_kinetic_branch_and_reversal_retain_the_elastic_tangent() {
        let law = law();
        let loaded = response(&law, &PartialSlipState::zero(), [30.0, 8.0, 0.2]);
        assert_eq!(loaded.next_state.sliding_regions(), [true, true]);
        for velocity in [[0.0; 3], [-0.01, -0.01, -0.001]] {
            let step = response(&law, &loaded.next_state, velocity);
            assert_eq!(step.region_plastic_increment, [[0.0; 3]; 2]);
            if velocity == [0.0; 3] {
                assert_eq!(step.next_state.sliding_regions(), [true, true]);
            }
            let block = impulse_tangent(&law, &step, ELL, DT).unwrap();
            for row in 0..3 {
                for col in 0..3 {
                    let expected = if row != col {
                        0.0
                    } else if row == 2 {
                        1.0
                    } else {
                        0.01
                    };
                    close(block[row][col], expected, 1e-14);
                }
            }
        }
    }

    #[test]
    fn exact_yield_point_has_one_linear_krylov_operator_and_one_law_evaluation() {
        let law = law();
        let gate = CancelGate::new();
        let evaluations = Cell::new(0);
        let identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let problem = PartialSlipProblem::new(3, &identity, 3, &gate, |velocity: &[f64]| {
            evaluations.set(evaluations.get() + 1);
            let step = response(
                &law,
                &PartialSlipState::zero(),
                velocity.try_into().unwrap(),
            );
            let residual = velocity
                .iter()
                .zip(impulse(&step))
                .map(|(v, p)| v - p)
                .collect();
            Ok((residual, vec![impulse_tangent(&law, &step, ELL, DT)?]))
        });
        // Weak region: k=6000 N/m, static capacity=24 N, dt=0.001 s.
        // The law selects elastic at equality. Crossing this point drops to
        // kinetic traction; a directional finite difference is not a linear
        // Krylov operator here.
        let point = [4.0, 0.0, 0.0];
        let mut residual = [0.0; 3];
        problem.residual(&point, &mut residual);
        assert_eq!(evaluations.get(), 1);
        for direction in [
            [1.0, 2.0, -3.0],
            [-1.0, -2.0, 3.0],
            [2.0, -1.0, 4.0],
            [0.0; 3],
        ] {
            let mut output = [0.0; 3];
            problem.jacobian_apply(&point, &direction, &mut output);
            for j in 0..3 {
                close(
                    output[j],
                    direction[j] * if j == 2 { 2.0 } else { 1.01 },
                    1e-14,
                );
            }
            assert_eq!(evaluations.get(), 1);
        }
        problem.residual(&point, &mut residual);
        assert_eq!(evaluations.get(), 1);
        problem.residual(&[3.0, 0.0, 0.0], &mut residual);
        assert_eq!(evaluations.get(), 2);
        assert!(!problem.has_failed());
    }
}
