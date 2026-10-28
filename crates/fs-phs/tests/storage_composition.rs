//! Independent work checks for composed storage with arbitrary energy datums.
use fs_phs::{
    DescriptorPortHamiltonian, PortHamiltonian, QuadraticStorage, ScalarLaw, SeparableStorage,
    StepWorkspace, Storage, SumStorage, discrete_gradient, step, step_descriptor,
};

struct ShiftedQuartic(f64);
impl Storage for ShiftedQuartic {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        self.0 + 0.25 * x[0].powi(4)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out[0] = x[0].powi(3);
    }
    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        let (da, db, m) = (a[0] - mid[0], b[0] - mid[0], mid[0]);
        Some(
            1.5 * m * m * (db * db - da * da)
                + m * (db.powi(3) - da.powi(3))
                + 0.25 * (db.powi(4) - da.powi(4)),
        )
    }
}

fn storage(datum: f64) -> Box<dyn Storage> {
    Box::new(SumStorage {
        a: (Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap()), 1),
        b: (Box::new(ShiftedQuartic(datum)), 1),
    })
}
fn system(datum: f64) -> PortHamiltonian {
    PortHamiltonian::new(
        2,
        1,
        vec![0.0, -1.0, 1.0, 0.0],
        vec![0.2, 0.0, 0.0, 0.0],
        vec![1.0, 0.0],
        storage(datum),
    )
    .unwrap()
}
fn physical_change(a: &[f64], b: &[f64]) -> f64 {
    0.5 * (b[0] - a[0]) * (b[0] + a[0])
        + 0.25 * (b[1] - a[1]) * (b[1] + a[1]) * (b[1] * b[1] + a[1] * a[1])
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 2e-12 * (1.0 + a.abs().max(b.abs())),
        "{a} != {b}"
    );
}

#[test]
fn quadratic_forces_do_not_subtract_large_endpoint_energies() {
    let law = QuadraticStorage::new(vec![2.0, 0.5, 0.5, 3.0], 2).unwrap();
    let a = [1e8, -2e8];
    let b = [1e8 + 1.0, -2e8 + 2.0];
    let mut expected = [0.0; 2];
    law.gradient(&[1e8 + 0.5, -2e8 + 1.0], &mut expected);
    assert_eq!(discrete_gradient(&law, &a, &b), expected);
    let scalar = QuadraticStorage::new(vec![1.0], 1).unwrap();
    let a = 9_007_199_254_740_992.0;
    assert_eq!(
        scalar.midpoint_energy_remainder(&[a], &[a + 2.0], &[a]),
        Some(2.0)
    );
    let stiff = QuadraticStorage::new(vec![1e308], 1).unwrap();
    let remainder = stiff
        .midpoint_energy_remainder(&[1e-154], &[2e-154], &[1.5e-154])
        .unwrap();
    assert!(remainder.is_finite());
    assert!(remainder.abs() < 1e-14);
}

#[test]
fn direct_sum_keeps_both_constitutive_remainders() {
    let law = storage(1e30);
    let gradient = discrete_gradient(law.as_ref(), &[2.0, 1.0], &[4.0, 3.0]);
    assert_eq!(gradient, [4.0, 9.0]);
    assert_eq!(gradient.iter().sum::<f64>() * 2.0, 26.0);
}

#[test]
fn native_and_prepared_motion_and_work_are_energy_datum_invariant() {
    let reference = system(0.0);
    let shifted = system(1e30);
    let initial = [0.7, 1.0];
    let input = [0.4];
    let expected = step(&reference, &initial, &input, 0.01).unwrap();
    let actual = step(&shifted, &initial, &input, 0.01).unwrap();
    for (a, b) in actual.x.iter().zip(&expected.x) {
        close(*a, *b);
    }
    close(actual.delta_h, physical_change(&initial, &actual.x));
    close(actual.delta_h, expected.delta_h);
    close(actual.supplied, expected.supplied);
    close(actual.dissipated, expected.dissipated);
    close(actual.balance_residual(), 0.0);
    assert!(actual.delta_h.abs() > 1e-4);
    let mut workspace = StepWorkspace::new(&shifted).unwrap();
    let (mut endpoint, mut output) = ([0.0; 2], [0.0]);
    let record = workspace
        .step_into(&shifted, &initial, &input, 0.01, &mut endpoint, &mut output)
        .unwrap();
    for (a, b) in endpoint.iter().zip(&expected.x) {
        close(*a, *b);
    }
    close(record.delta_h, physical_change(&initial, &endpoint));
    close(record.balance_residual(), 0.0);
}

#[test]
fn descriptor_energy_change_uses_the_same_constitutive_difference() {
    let system = DescriptorPortHamiltonian::new(
        2,
        2,
        1,
        vec![0.0, -1.0, 1.0, 0.0],
        vec![0.2, 0.0, 0.0, 0.0],
        vec![1.0, 0.0],
        storage(1e30),
    )
    .unwrap();
    let initial = [0.7, 1.0];
    let record = step_descriptor(&system, &initial, &[0.4], 0.01).unwrap();
    close(record.delta_h, physical_change(&initial, &record.x));
    close(record.balance_residual(), 0.0);
    assert!(record.delta_h.abs() > 1e-4);
}

#[test]
fn unsupported_child_does_not_manufacture_a_complete_remainder() {
    let law = SumStorage {
        a: (Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap()), 1),
        b: (
            Box::new(SeparableStorage::new(vec![ScalarLaw {
                h: |x| 0.25 * x.powi(4),
                dh: |x| x.powi(3),
            }])),
            1,
        ),
    };
    assert_eq!(
        law.midpoint_energy_remainder(&[1.0, 1.0], &[2.0, 2.0], &[1.5, 1.5]),
        None
    );
}

struct ShiftedQuadratic {
    original: QuadraticStorage,
    datum: f64,
}
impl Storage for ShiftedQuadratic {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        self.datum + self.original.hamiltonian(x)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        self.original.gradient(x, out);
    }
    fn quadratic_part(&self) -> Option<fs_phs::linear_exact::QuadraticPart> {
        self.original.quadratic_part()
    }
    fn exact_quadratic_matrix(&self) -> Option<&[f64]> {
        self.original.exact_quadratic_matrix()
    }
    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        self.original.midpoint_energy_remainder(a, b, mid)
    }
}

#[test]
fn complete_exact_flow_keeps_small_physical_work_above_an_arbitrary_energy_datum() {
    use fs_phs::linear_exact::{LinearSplitBudget, split_quadratic, step_render};
    let initial = 0.8_f64;
    let force = 0.4;
    let rate = 2.0;
    let dt = 0.25;
    let equilibrium = force / rate;
    let departure = initial - equilibrium;
    let decay = (-rate * dt).exp();
    let expected_state = equilibrium + departure * decay;
    let first_integral = -(-rate * dt).exp_m1() / rate;
    let second_integral = -(-2.0 * rate * dt).exp_m1() / (2.0 * rate);
    let work = force * (equilibrium * dt + departure * first_integral);
    let heat = rate * (equilibrium * equilibrium * dt
        + 2.0 * equilibrium * departure * first_integral
        + departure * departure * second_integral);
    for datum in [0.0, 1e30, -1e30] {
        let system = PortHamiltonian::new(1, 1, vec![0.0], vec![rate], vec![1.0],
            Box::new(ShiftedQuadratic {
                original: QuadraticStorage::new(vec![1.0], 1).unwrap(), datum,
            })).unwrap();
        let split = split_quadratic(&system, LinearSplitBudget {
            max_modes: 1, max_omega_dt: 1.0,
        }).unwrap();
        let record = step_render(&system, &split, &[initial], &[force], dt).unwrap();
        close(record.x[0], expected_state);
        close(record.supplied, work);
        close(record.dissipated, heat);
        close(record.delta_h, 0.5 * (expected_state - initial) * (expected_state + initial));
        close(record.balance_residual(), 0.0);
        assert!(record.delta_h.abs() > 0.1);
    }
}

#[test]
fn frequency_matching_retains_the_original_remainder_and_its_quadratic_correction() {
    use fs_phs::linear_exact::{LinearPropagator, LinearSplitBudget, split_quadratic, step_split};
    let initial = [0.8, 0.3];
    let dt = 0.1;
    let mut reference = None;
    for datum in [0.0, 1e30, -1e30] {
        let system = PortHamiltonian::new(2, 1, vec![0.0, 1.0, -1.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.2], vec![0.0, 1.0],
            Box::new(ShiftedQuadratic {
                original: QuadraticStorage::new(vec![4.0, 0.0, 0.0, 1.0], 2).unwrap(), datum,
            })).unwrap();
        let split = split_quadratic(&system, LinearSplitBudget {
            max_modes: 2, max_omega_dt: 1.0,
        }).unwrap();
        let record = step_split(&system, &split, LinearPropagator::FrequencyMatched,
            &initial, &[0.4], dt).unwrap();
        if let Some(expected) = &reference {
            let expected: &fs_phs::StepRecord = expected;
            for (actual, expected) in record.x.iter().zip(&expected.x) {
                close(*actual, *expected);
            }
            close(record.delta_h, expected.delta_h);
            close(record.supplied, expected.supplied);
            close(record.dissipated, expected.dissipated);
        }
        let correction = split.quadratic_correction(dt).unwrap();
        let mut modified_change = 0.0;
        for i in 0..2 {
            for j in 0..2 {
                let q = correction[i * 2 + j]
                    + if i == j { if i == 0 { 4.0 } else { 1.0 } } else { 0.0 };
                modified_change += 0.5 * (record.x[i] - initial[i]) * q
                    * (record.x[j] + initial[j]);
            }
        }
        close(record.delta_h, modified_change);
        close(record.balance_residual(), 0.0);
        reference = Some(record);
    }
}
