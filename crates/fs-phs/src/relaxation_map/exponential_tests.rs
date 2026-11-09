//! Native exact-linear recovery with the original nonlinear Maxwell owner.
use super::*;
use crate::linear_exact::{LinearSplitBudget, PreparedNonlinearZoh};
use crate::{PreparedStepError, QuadraticStorage};

#[derive(Debug)]
struct PairStrain {
    coefficients: [f64; 3],
    declared_quadratic: bool,
}
impl RelaxationStrainMap for PairStrain {
    fn state_dimension(&self) -> usize {
        4
    }
    fn strain(&self, x: &[f64]) -> f64 {
        let [a, b, c] = self.coefficients;
        a * x[0] * x[0] + b * x[0] * x[2] + c * x[2] * x[2]
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        let [a, b, c] = self.coefficients;
        out.copy_from_slice(&[
            2.0 * a * x[0] + b * x[2],
            0.0,
            b * x[0] + 2.0 * c * x[2],
            0.0,
        ]);
    }
    fn hessian_vector(&self, _: &[f64], d: &[f64], out: &mut [f64]) {
        self.gradient(d, out);
    }
    fn is_quadratic_strain(&self) -> bool {
        self.declared_quadratic
    }
}
const PAIRS: [[f64; 3]; 2] = [[0.5, 0.25, 0.125], [0.25, -0.5, 0.375]];
const STIFFNESS: [f64; 2] = [3.0, 5.0];
const TAU: [f64; 2] = [0.15, 0.04];
const INPUT: [f64; 2] = [0.2, -0.1];

fn budget(n: usize) -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: n,
        max_omega_dt: 0.9 * core::f64::consts::PI,
    }
}
fn mapped(
    coefficients: [f64; 3],
    stiffness: f64,
    tau: f64,
    declared: bool,
) -> MappedRelaxationBranch {
    MappedRelaxationBranch {
        strain: Rc::new(PairStrain {
            coefficients,
            declared_quadratic: declared,
        }),
        stiffness,
        relaxation_time_s: tau,
    }
}
fn system(declared: bool) -> PortHamiltonian {
    let q = vec![
        4.0, 0.0, 0.3, 0.0, 0.0, 1.0, 0.0, 0.0, 0.3, 0.0, 9.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    ];
    let j = vec![
        0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0,
    ];
    let r = vec![
        0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3,
    ];
    let g = vec![0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0];
    PortHamiltonian::new(
        4,
        2,
        j,
        r,
        g,
        Box::new(QuadraticStorage::new(q, 4).unwrap()),
    )
    .unwrap()
    .with_mapped_relaxation_branches(
        (0..2)
            .map(|i| mapped(PAIRS[i], STIFFNESS[i], TAU[i], declared))
            .collect(),
    )
    .unwrap()
}

#[test]
fn fixed_nonlinear_strain_has_exact_monotone_memory_recovery_and_positive_heat() {
    let initial = [0.4, 0.0, -0.2, 0.0, -0.07];
    let epsilon = 0.5 * initial[0] * initial[0]
        + 0.25 * initial[0] * initial[2]
        + 0.125 * initial[2] * initial[2];
    let target = 2.0 * epsilon;
    let gap = target - initial[4];
    for ratio in [1e-5, 0.2, 20.0] {
        let system = PortHamiltonian::new(
            4,
            0,
            vec![0.0; 16],
            vec![0.0; 16],
            vec![],
            Box::new(QuadraticStorage::new(vec![0.0; 16], 4).unwrap()),
        )
        .unwrap()
        .with_mapped_relaxation_branches(vec![mapped(PAIRS[0], 4.0, 0.01, true)])
        .unwrap();
        let part = system.storage().quadratic_part().unwrap();
        assert!(!part.exact);
        assert_eq!(part.matrix[24], 1.0);
        let energy = system.hamiltonian(&initial);
        let mut flow = PreparedNonlinearZoh::new_owned(system, 0.01 * ratio, budget(5)).unwrap();
        let record = flow.step(&initial, &[]).unwrap();
        assert_eq!(&record.x[..4], &initial[..4]);
        let expected = target - gap * det::exp(-ratio);
        let heat = 0.5 * gap * gap * (-det::expm1(-2.0 * ratio));
        assert!((record.x[4] - expected).abs() < 2e-12);
        assert!(record.x[4] > initial[4] && record.x[4] <= target);
        assert!((record.dissipated - heat).abs() < 2e-10 * heat.max(1e-10));
        assert!((record.delta_h + record.dissipated).abs() < 2e-12 * energy);
        assert_eq!(record.supplied, 0.0);
    }
}

// Independent continuous q,p,z,heat,work equations. These do not call a native
// storage/derivative callback or a native time integration owner.
fn rhs(x: [f64; 8]) -> [f64; 8] {
    let mut d = [
        x[1],
        -4.0 * x[0] - 0.3 * x[2] - 0.2 * x[1] + INPUT[0],
        x[3],
        -0.3 * x[0] - 9.0 * x[2] - 0.3 * x[3] + INPUT[1],
        0.0,
        0.0,
        0.2 * x[1] * x[1] + 0.3 * x[3] * x[3],
        INPUT[0] * x[1] + INPUT[1] * x[3],
    ];
    for i in 0..2 {
        let [a, b, c] = PAIRS[i];
        let root = det::sqrt(STIFFNESS[i]);
        let epsilon = a * x[0] * x[0] + b * x[0] * x[2] + c * x[2] * x[2];
        let effort = root * epsilon - x[4 + i];
        d[1] -= root * effort * (2.0 * a * x[0] + b * x[2]);
        d[3] -= root * effort * (b * x[0] + 2.0 * c * x[2]);
        d[4 + i] = effort / TAU[i];
        d[6] += effort * effort / TAU[i];
    }
    d
}
fn continuous(mut x: [f64; 8], time: f64) -> [f64; 8] {
    let h = time / 16384.0;
    for _ in 0..16384 {
        let a = rhs(x);
        let b = rhs(core::array::from_fn(|i| x[i] + 0.5 * h * a[i]));
        let c = rhs(core::array::from_fn(|i| x[i] + 0.5 * h * b[i]));
        let d = rhs(core::array::from_fn(|i| x[i] + h * c[i]));
        for i in 0..8 {
            x[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    x
}

#[test]
fn full_cross_mode_mapped_storage_heat_and_work_converge_in_the_same_exact_owner() {
    let initial = [0.6, 0.1, -0.4, 0.2, 0.07, -0.02];
    let reference = continuous([0.6, 0.1, -0.4, 0.2, 0.07, -0.02, 0.0, 0.0], 0.4);
    let mut errors = Vec::new();
    for steps in [40, 80, 160] {
        let mut system = system(true);
        let part = system.storage().quadratic_part().unwrap();
        assert!(!part.exact);
        assert_eq!(part.matrix[2], 0.3);
        assert_eq!(part.matrix[12], 0.3);
        assert_eq!(part.matrix[28], 1.0);
        assert_eq!(part.matrix[35], 1.0);
        let energy = system.hamiltonian(&initial);
        let observation = system.share_storage();
        let mut flow =
            PreparedNonlinearZoh::new_owned(system, 0.4 / f64::from(steps), budget(6)).unwrap();
        let mut x = initial.to_vec();
        let (mut heat, mut work) = (0.0, 0.0);
        for _ in 0..steps {
            let record = flow.step(&x, &INPUT).unwrap();
            assert!(record.dissipated >= 0.0);
            assert!((record.delta_h + record.dissipated - record.supplied).abs() < 2e-11 * energy);
            heat += record.dissipated;
            work += record.supplied;
            x = record.x;
        }
        assert!(heat > 0.0);
        assert!((observation.hamiltonian(&x) + heat - work - energy).abs() < 5e-10 * energy);
        let error = x
            .iter()
            .zip(reference)
            .map(|(a, b)| (a - b).abs())
            .chain([(heat - reference[6]).abs(), (work - reference[7]).abs()])
            .fold(0.0_f64, f64::max);
        errors.push(error);
    }
    assert!(
        errors[1] < 0.35 * errors[0] && errors[2] < 0.35 * errors[1],
        "nonlinear mapped strain refinement: {errors:?}"
    );
    assert!(errors[2] < 3e-5, "complete response and heat: {errors:?}");
}

#[test]
fn undeclared_maps_refuse_and_failed_native_solves_leave_all_outputs_intact() {
    let opaque = system(false);
    assert!(opaque.storage().quadratic_part().is_none());
    assert!(PreparedNonlinearZoh::new_owned(opaque, 0.01, budget(6)).is_err());
    let initial = [0.6, 0.1, -0.4, 0.2, 0.07, -0.02];
    let mut flow = PreparedNonlinearZoh::new_owned(system(true), 0.01, budget(6)).unwrap();
    let mut next = [7.0; 6];
    let mut output = [9.0; 2];
    assert!(matches!(
        flow.step_into_controlled(&initial, &INPUT, &mut next, &mut output, || true),
        Err(PreparedStepError::Cancelled)
    ));
    assert_eq!(next, [7.0; 6]);
    assert_eq!(output, [9.0; 2]);
    flow.set_max_iterations(0).unwrap();
    assert!(
        flow.step_into(&initial, &INPUT, &mut next, &mut output)
            .is_err()
    );
    assert_eq!(next, [7.0; 6]);
    assert_eq!(output, [9.0; 2]);
    flow.set_max_iterations(50).unwrap();
    let accepted = flow.step(&initial, &INPUT).unwrap();
    let fresh = PreparedNonlinearZoh::new_owned(system(true), 0.01, budget(6))
        .unwrap()
        .step(&initial, &INPUT)
        .unwrap();
    assert_eq!(accepted.x, fresh.x);
    assert_eq!(accepted.dissipated.to_bits(), fresh.dissipated.to_bits());
}
