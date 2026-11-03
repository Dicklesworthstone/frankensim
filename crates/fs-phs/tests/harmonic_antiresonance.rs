//! Original-row refinement at an exact coupled antiresonance.
use fs_math::c64::C64;
use fs_phs::{
    PortHamiltonian, QuadraticStorage,
    harmonic::{HarmonicBudget, HarmonicError, HarmonicModel},
};

const N: usize = 7;
const BETA: f64 = 13.7;
const COUPLINGS: [f64; 3] = [173.4, 320.1, 56.7];
const TIMES: [f64; 3] = [1e-4, 1e-3, 0.02];

// A host with three real relaxation states drives a lossless mass/spring.
// q, v, z_0..z_2, relative spring displacement, load momentum are original
// physical coordinates. At omega=1000 the host response is exactly zero;
// finite load energy and its opposing force remain. No load elimination.
fn system() -> PortHamiltonian {
    let mut q = vec![0.0; N * N];
    let mut j = q.clone();
    let mut r = q.clone();
    let mut g = vec![0.0; N];
    q[0] = 2100.3_f64.powi(2);
    q[N + 1] = 1.0;
    j[1] = 1.0;
    j[N] = -1.0;
    r[N + 1] = 2.7;
    g[1] = BETA;
    for k in 0..3 {
        let at = k + 2;
        q[0] += COUPLINGS[k] * COUPLINGS[k];
        q[at] = -COUPLINGS[k];
        q[at * N] = -COUPLINGS[k];
        q[at * N + at] = 1.0;
        r[at * N + at] = 1.0 / TIMES[k];
    }
    q[5 * N + 5] = 1_000_000.0;
    q[6 * N + 6] = 1.0;
    j[5 * N + 6] = 1.0;
    j[6 * N + 5] = -1.0;
    j[N + 5] = BETA;
    j[5 * N + 1] = -BETA;
    PortHamiltonian::new(
        N,
        1,
        j,
        r,
        g,
        Box::new(QuadraticStorage::new(q, N).unwrap()),
    )
    .unwrap()
}

fn budget() -> HarmonicBudget {
    HarmonicBudget {
        max_dimension: N,
        max_dense_bytes: 10_000,
        max_dense_work: 1225,
        relative_residual_tolerance: 1e-10,
    }
}

fn near(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() < 2e-12 * expected.abs(),
        "{actual:?} vs {expected:?}"
    );
}

#[test]
fn g1_exact_host_zero_keeps_load_energy_and_original_componentwise_accuracy() {
    let model = HarmonicModel::from_ode(&system(), budget(), || false).unwrap();
    let factor = model.factor(1000.0, || false).unwrap();
    for force in [C64::new(0.6, -0.8), C64::ONE, C64::new(-0.33, 0.51)] {
        let result = factor.response(&[force], || false).unwrap();
        let spring = force.scale(-1e-6);
        let momentum = C64::new(0.0, -1000.0) * spring;
        near(result.state[5], spring);
        near(result.state[6], momentum);
        assert!(result.state[0].abs() < 1e-14 * spring.abs() / BETA);
        assert!(result.state[1].abs() < 1e-14 * momentum.abs() / BETA);
        for i in 0..3 {
            assert!(result.state[i + 2].abs() < 1e-14 * spring.abs() * COUPLINGS[i] / BETA);
        }
        let energy = 0.25 * (1_000_000.0 * spring.norm_sq() + momentum.norm_sq());
        near(
            C64::from_re(result.mean_stored_energy_j),
            C64::from_re(energy),
        );
        assert!(result.dissipated_power_w < 1e-24 * 1000.0 * energy);
        assert!(result.supplied_power_w.abs() < 1e-14 * 1000.0 * energy);
        assert!(result.relative_residual < 1e-10);
    }
}

#[test]
fn g0_refinement_capacity_and_late_cancellation_do_not_poison_the_shared_factor() {
    let source = system();
    let mut small = budget();
    small.max_dense_work -= 1;
    assert!(matches!(
        HarmonicModel::from_ode(&source, small, || false),
        Err(HarmonicError::Invalid { .. })
    ));
    let model = HarmonicModel::from_ode(&source, budget(), || false).unwrap();
    let factor = model.factor(1000.0, || false).unwrap();
    let input = [C64::new(0.6, -0.8)];
    let mut checks = 0;
    let before = factor
        .response(&input, || {
            checks += 1;
            false
        })
        .unwrap();
    let mut at = 0;
    assert_eq!(
        factor.response(&input, || {
            at += 1;
            at == checks - 1
        }),
        Err(HarmonicError::Cancelled)
    );
    assert_eq!(factor.response(&input, || false).unwrap(), before);
    let zero = factor.response(&[C64::ZERO], || false).unwrap();
    assert_eq!(zero.state, vec![C64::ZERO; N]);
    assert_eq!(zero.relative_residual, 0.0);
}
