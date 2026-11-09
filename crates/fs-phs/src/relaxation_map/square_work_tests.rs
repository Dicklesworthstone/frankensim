use super::*;
use fs_math::dd::Dd;

#[derive(Debug)]
struct Strain;
impl RelaxationStrainMap for Strain {
    fn state_dimension(&self) -> usize {
        2
    }
    fn is_quadratic_strain(&self) -> bool {
        true
    }
    fn strain(&self, x: &[f64]) -> f64 {
        0.125 + 1.25 * x[0] - 0.75 * x[1] + x[0] * x[0] + 0.5 * x[0] * x[1] + 0.75 * x[1] * x[1]
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out[0] = 1.25 + 2.0 * x[0] + 0.5 * x[1];
        out[1] = -0.75 + 0.5 * x[0] + 1.5 * x[1];
    }
    fn hessian_vector(&self, _: &[f64], d: &[f64], out: &mut [f64]) {
        out[0] = 2.0 * d[0] + 0.5 * d[1];
        out[1] = 0.5 * d[0] + 1.5 * d[1];
    }
}

// Independent two-word polynomial energy and all three work derivatives.
// It does not call the production strain/remainder/Maxwell-energy callbacks.
fn reference(x: &[f64], z: f64, scale: f64) -> (Dd, [Dd; 3]) {
    let c = Dd::from_f64;
    let (u, v, z, s) = (c(x[0]), c(x[1]), c(z), c(scale));
    let e = c(0.125) + c(1.25) * u - c(0.75) * v + u * u + c(0.5) * u * v + c(0.75) * v * v;
    let r = s * e - z;
    (
        c(0.5) * r * r,
        [
            r * s * (c(1.25) + c(2.0) * u + c(0.5) * v),
            r * s * (-c(0.75) + c(0.5) * u + c(1.5) * v),
            -r,
        ],
    )
}

#[test]
fn original_square_work_preserves_tiny_cubic_terms_and_rounded_midpoint() {
    let c = Dd::from_f64;
    for exponent in [0, 8, 16, 22] {
        let step = 2.0_f64.powi(-exponent);
        let a = [0.0625, -0.03125];
        let b = [a[0] + 0.5 * step, a[1] - 0.25 * step];
        let za = 0.0625;
        let zb = za + 0.125 * step;
        for rounded in [false, true] {
            let mut mid = [(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5];
            let mut zm = (za + zb) * 0.5;
            if rounded {
                mid[0] = mid[0].next_up();
                zm = zm.next_down();
            }
            for scale in [1.0, 3.5] {
                let (ha, _) = reference(&a, za, scale);
                let (hb, _) = reference(&b, zb, scale);
                let (_, g) = reference(&mid, zm, scale);
                let expected = (hb
                    - ha
                    - g[0] * (c(b[0]) - c(a[0]))
                    - g[1] * (c(b[1]) - c(a[1]))
                    - g[2] * (c(zb) - c(za)))
                .to_f64();
                let actual = quadratic_strain_midpoint_energy_remainder(
                    &Strain,
                    scale,
                    &a,
                    &b,
                    &mid,
                    [za, zb, zm],
                )
                .unwrap();
                assert!(
                    (actual - expected).abs() <= 2e-12 * expected.abs() + 2e-30,
                    "exponent={exponent}, rounded={rounded}, {actual:e} != {expected:e}"
                );
            }
        }
    }
}

#[test]
fn equilibrium_energy_uses_the_original_coefficient_without_sqrt_renormalization() {
    let a = [0.0625, -0.03125];
    let b = [0.09375, -0.0390625];
    let mid = [(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5];
    let coefficient = 7.3;
    let (ha, _) = reference(&a, 0.0, 1.0);
    let (hb, _) = reference(&b, 0.0, 1.0);
    let (_, g) = reference(&mid, 0.0, 1.0);
    let c = Dd::from_f64;
    let expected = (c(coefficient)
        * (hb - ha - g[0] * (c(b[0]) - c(a[0])) - g[1] * (c(b[1]) - c(a[1]))))
    .to_f64();
    let actual = coefficient
        * quadratic_strain_midpoint_energy_remainder(&Strain, 1.0, &a, &b, &mid, [0.0; 3]).unwrap();
    assert!((actual - expected).abs() < 3e-14 * expected.abs());

    let coupling = RelaxationCoupling::Mapped {
        strain: Rc::new(Strain),
        scale: 3.5,
    };
    let offsets = [0.01, 0.03, 0.02];
    assert_eq!(
        coupling
            .mapped_remainder(&a, &b, &mid, offsets[0], offsets[1], offsets[2])
            .unwrap()
            .to_bits(),
        quadratic_strain_midpoint_energy_remainder(&Strain, 3.5, &a, &b, &mid, offsets)
            .unwrap()
            .to_bits()
    );
}

#[test]
fn incomplete_nonfinite_or_undeclared_square_work_refuses() {
    #[derive(Debug)]
    struct Opaque;
    impl RelaxationStrainMap for Opaque {
        fn state_dimension(&self) -> usize {
            2
        }
        fn strain(&self, _: &[f64]) -> f64 {
            panic!("undeclared map evaluated")
        }
        fn gradient(&self, _: &[f64], _: &mut [f64]) {
            panic!("undeclared map evaluated")
        }
        fn hessian_vector(&self, _: &[f64], _: &[f64], _: &mut [f64]) {
            panic!("undeclared map evaluated")
        }
    }
    let x = [0.0, 0.0];
    assert!(
        quadratic_strain_midpoint_energy_remainder(&Opaque, 1.0, &x, &x, &x, [0.0; 3]).is_none()
    );
    for scale in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            quadratic_strain_midpoint_energy_remainder(&Strain, scale, &x, &x, &x, [0.0; 3])
                .is_none()
        );
    }
    assert!(
        quadratic_strain_midpoint_energy_remainder(&Strain, 1.0, &x[..1], &x, &x, [0.0; 3])
            .is_none()
    );
    assert!(
        quadratic_strain_midpoint_energy_remainder(&Strain, 1.0, &x, &x[..1], &x, [0.0; 3])
            .is_none()
    );
    assert!(
        quadratic_strain_midpoint_energy_remainder(&Strain, 1.0, &x, &x, &x[..1], [0.0; 3])
            .is_none()
    );
    assert!(
        quadratic_strain_midpoint_energy_remainder(&Strain, 1.0, &x, &x, &x, [0.0, f64::NAN, 0.0])
            .is_none()
    );
    let bad = [f64::INFINITY, 0.0];
    assert!(
        quadratic_strain_midpoint_energy_remainder(&Strain, 1.0, &bad, &x, &x, [0.0; 3]).is_none()
    );
}

#[derive(Debug)]
struct MixedStrain {
    curvature: f64,
    declared: bool,
}
const ROW: [f64; 4] = [1.25, 0.0, -0.75, 0.0];
impl RelaxationStrainMap for MixedStrain {
    fn state_dimension(&self) -> usize {
        4
    }
    fn is_quadratic_strain(&self) -> bool {
        true
    }
    fn zero_constant_linear_strain(&self) -> Option<Vec<f64>> {
        self.declared.then(|| ROW.to_vec())
    }
    fn strain(&self, x: &[f64]) -> f64 {
        ROW.iter().zip(x).map(|(b, x)| b * x).sum::<f64>()
            + self.curvature * (x[0] * x[0] + 0.5 * x[0] * x[2] + 0.75 * x[2] * x[2])
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out.copy_from_slice(&[
            ROW[0] + self.curvature * (2.0 * x[0] + 0.5 * x[2]),
            0.0,
            ROW[2] + self.curvature * (0.5 * x[0] + 1.5 * x[2]),
            0.0,
        ]);
    }
    fn hessian_vector(&self, _: &[f64], d: &[f64], out: &mut [f64]) {
        out.copy_from_slice(&[
            self.curvature * (2.0 * d[0] + 0.5 * d[2]),
            0.0,
            self.curvature * (0.5 * d[0] + 1.5 * d[2]),
            0.0,
        ]);
    }
}
fn base() -> PortHamiltonian {
    let q = vec![
        4.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 9.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    ];
    let j = vec![
        0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0,
    ];
    let r = vec![
        0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3,
    ];
    PortHamiltonian::new(
        4,
        2,
        j,
        r,
        vec![0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0],
        Box::new(crate::QuadraticStorage::new(q, 4).unwrap()),
    )
    .unwrap()
}
fn mixed(curvature: f64, declared: bool) -> PortHamiltonian {
    base()
        .with_mapped_relaxation_branches(vec![MappedRelaxationBranch {
            strain: Rc::new(MixedStrain {
                curvature,
                declared,
            }),
            stiffness: 7.3,
            relaxation_time_s: 0.025,
        }])
        .unwrap()
}
fn budget() -> crate::linear_exact::LinearSplitBudget {
    crate::linear_exact::LinearSplitBudget {
        max_modes: 5,
        max_omega_dt: 0.9 * core::f64::consts::PI,
    }
}

#[test]
fn declared_linear_strain_retains_the_complete_exact_maxwell_block() {
    use crate::linear_exact::PreparedNonlinearZoh;
    let direct = base()
        .with_relaxation_branches(vec![crate::RelaxationBranch {
            projection: ROW.to_vec(),
            stiffness: 7.3,
            relaxation_time_s: 0.025,
        }])
        .unwrap();
    let mapped = mixed(0.0, true);
    let expected = direct.storage().quadratic_part().unwrap();
    let actual = mapped.storage().quadratic_part().unwrap();
    assert_eq!(actual.matrix, expected.matrix);
    assert!(!actual.exact); // declaration supplies a split, not a claim of linearity.
    assert_ne!(actual.matrix[4], 0.0);
    assert_ne!(actual.matrix[14], 0.0);
    let initial = [0.4, 0.1, -0.2, 0.2, 0.03];
    let input = [0.1, -0.2];
    assert!((mapped.hamiltonian(&initial) - direct.hamiltonian(&initial)).abs() < 1e-14);
    for (a, b) in mapped.effort(&initial).iter().zip(direct.effort(&initial)) {
        assert!((a - b).abs() < 2e-14);
    }
    let a = PreparedNonlinearZoh::new_owned(mapped, 0.03, budget())
        .unwrap()
        .step(&initial, &input)
        .unwrap();
    let b = PreparedNonlinearZoh::new_owned(direct, 0.03, budget())
        .unwrap()
        .step(&initial, &input)
        .unwrap();
    for (a, b) in a.x.iter().zip(&b.x) {
        assert!((a - b).abs() < 2e-12);
    }
    assert!((a.dissipated - b.dissipated).abs() < 2e-12);
    assert!((a.supplied - b.supplied).abs() < 2e-12);
    assert!(a.dissipated > 0.0);
}

#[test]
fn declaring_the_source_row_keeps_full_nonlinear_storage_force_and_work() {
    use crate::linear_exact::PreparedNonlinearZoh;
    let mut old = mixed(0.2, false);
    let declared = mixed(0.2, true);
    let reference = mixed(0.0, true).storage().quadratic_part().unwrap();
    let part = declared.storage().quadratic_part().unwrap();
    assert_eq!(part.matrix, reference.matrix);
    let old_q = old.storage().quadratic_part().unwrap();
    assert_eq!(old_q.matrix[0], 4.0);
    assert_eq!(old_q.matrix[12], 9.0);
    assert_eq!(old_q.matrix[24], 1.0);
    assert_eq!(old_q.matrix[4], 0.0);
    let a = [0.4, 0.1, -0.2, 0.2, 0.03];
    let b = [0.4001, 0.1002, -0.2003, 0.2001, 0.0304];
    let mid: Vec<_> = a.iter().zip(b).map(|(a, b)| 0.5 * (a + b)).collect();
    assert_eq!(
        old.hamiltonian(&a).to_bits(),
        declared.hamiltonian(&a).to_bits()
    );
    assert_eq!(old.effort(&a), declared.effort(&a));
    assert_eq!(
        old.storage()
            .midpoint_energy_remainder(&a, &b, &mid)
            .unwrap()
            .to_bits(),
        declared
            .storage()
            .midpoint_energy_remainder(&a, &b, &mid)
            .unwrap()
            .to_bits()
    );
    let physical = old.share_storage();
    let mut flow = PreparedNonlinearZoh::new_owned(declared, 0.001, budget()).unwrap();
    let initial = physical.hamiltonian(&a);
    let result = flow.step(&a, &[0.1, -0.2]).unwrap();
    assert!(result.dissipated > 0.0);
    assert!(
        (physical.hamiltonian(&result.x) - initial + result.dissipated - result.supplied).abs()
            < 2e-12
    );
    // Independent complete material force, including the geometric derivative.
    let root = 7.3_f64.sqrt();
    let strain =
        1.25 * a[0] - 0.75 * a[2] + 0.2 * (a[0] * a[0] + 0.5 * a[0] * a[2] + 0.75 * a[2] * a[2]);
    let stress = root * (root * strain - a[4]);
    let force = physical.effort(&a);
    assert!(
        (force[0] - (4.0 * a[0] + stress * (1.25 + 0.2 * (2.0 * a[0] + 0.5 * a[2])))).abs() < 2e-14
    );
    assert!(
        (force[2] - (9.0 * a[2] + stress * (-0.75 + 0.2 * (0.5 * a[0] + 1.5 * a[2])))).abs()
            < 2e-14
    );
}
