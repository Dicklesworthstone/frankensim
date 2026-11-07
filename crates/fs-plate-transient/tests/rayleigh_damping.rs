//! G0/G1/G3: prestress must not turn authored Rayleigh loss into an energy source.
use fs_plate::PlateModel;
use fs_plate_transient::{
    NoPlateLoad, PlateDynamics, PlateDynamicsBudget, PlateDynamicsError, PlateDynamicsParameters,
};
use fs_solver::NewtonKrylovConfig;
use fs_sparse::{Coo, Csr, LdltOptions};
use fs_time::galpha::{
    ImplicitSolveConfig, OperatorGeneralizedAlpha, SecondOrderProblem, SecondOrderState,
};

fn matrix(entries: [f64; 9]) -> Csr {
    let mut coo = Coo::new(3, 3);
    for row in 0..3 {
        for column in 0..3 {
            if entries[3 * row + column] != 0.0 {
                coo.push(row, column, entries[3 * row + column]);
            }
        }
    }
    coo.assemble()
}

fn plate(stiffness: [f64; 9], scales: [f64; 3]) -> PlateModel {
    PlateModel {
        k: matrix(core::array::from_fn(|i| {
            stiffness[i] * scales[i / 3] * scales[i % 3]
        })),
        m: matrix(core::array::from_fn(|i| {
            if i / 3 == i % 3 { scales[i / 3].powi(2) } else { 0.0 }
        })),
        dof_map: vec![Some(0), Some(1), Some(2)],
        free: 3,
    }
}

fn parameters(alpha: f64, beta: f64) -> PlateDynamicsParameters {
    PlateDynamicsParameters {
        mass_scale: 1.0,
        stiffness_scale: 1.0,
        mass_damping_per_s: alpha,
        stiffness_damping_s: beta,
    }
}

fn budget(plate: &PlateModel) -> PlateDynamicsBudget {
    PlateDynamicsBudget {
        max_dofs: 3,
        max_nonzeros: plate.m.nnz() + plate.k.nnz(),
        max_full_dofs: 3,
        max_load_parameters: 0,
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

#[test]
fn prestress_with_positive_rayleigh_coefficients_cannot_supply_kinetic_power() {
    let plate = plate([1.0, 2.0, 0.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0], [1.0; 3]);
    // Every C diagonal is nonnegative, but v=[1,-1,0] is an active direction.
    // For alpha=0, beta=1 the old accepted law gave v . (-C v) = +2 W.
    let velocity = [1.0, -1.0, 0.0];
    let mut resistance = [0.0; 3];
    plate.k.spmv(&velocity, &mut resistance);
    assert_eq!(-dot(&velocity, &resistance), 2.0);
    assert!(PlateDynamics::new(
        &plate, parameters(0.0, 1.0), &NoPlateLoad, budget(&plate), &mut || false,
    ).is_err());

    // Zero diagonal rows cannot be discarded when they retain cross coupling.
    let zero_row = self::plate(
        [1.0, 0.1, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0; 3],
    );
    assert!(PlateDynamics::new(
        &zero_row, parameters(0.0, 1.0), &NoPlateLoad, budget(&zero_row), &mut || false,
    ).is_err());
}

#[test]
fn actual_combination_admits_indefinite_stiffness_and_preserves_undamped_directions() {
    // K has eigenvalues 3,-1,-1. C=M+K has eigenvalues 4,0,0.
    let k = [1.0, -2.0, 0.0, -2.0, 1.0, 0.0, 0.0, 0.0, -1.0];
    for scales in [[1.0; 3], [2.0_f64.powi(-200), 2.0_f64.powi(200), 1.0]] {
        let plate = plate(k, scales);
        let dynamics = PlateDynamics::new(
            &plate, parameters(1.0, 1.0), &NoPlateLoad, budget(&plate), &mut || false,
        ).unwrap();
        let null_velocity = [1.0 / scales[0], 1.0 / scales[1], 7.0];
        let mut resistance = [0.0; 3];
        dynamics.damping_apply(&null_velocity, &mut resistance);
        assert_eq!(resistance, [0.0; 3], "admission must add no damping floor");
        let velocity = [1.0 / scales[0], -1.0 / scales[1], 0.0];
        dynamics.damping_apply(&velocity, &mut resistance);
        assert_eq!(dot(&velocity, &resistance), 8.0);
    }

    // Complete cancellation is a genuine conservative limit, including a
    // negative reference stiffness. Admission is about C, not stability of K.
    let plate = plate([-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0], [1.0; 3]);
    let dynamics = PlateDynamics::new(
        &plate, parameters(1.0, 1.0), &NoPlateLoad, budget(&plate), &mut || false,
    ).unwrap();
    let mut force = [0.0; 3];
    dynamics.damping_apply(&[2.0, -3.0, 4.0], &mut force);
    assert_eq!(force, [0.0; 3]);
}

#[test]
fn passive_combination_closes_original_midpoint_work_without_added_loss() {
    let plate = plate([1.0, -2.0, 0.0, -2.0, 1.0, 0.0, 0.0, 0.0, -1.0], [1.0; 3]);
    let dynamics = PlateDynamics::new(
        &plate, parameters(1.0, 1.0), &NoPlateLoad, budget(&plate), &mut || false,
    ).unwrap();
    // Excite only the stable eigenvector. The separate undamped unstable
    // directions of K are unexcited; this test does not claim elastic stability.
    let mut state = SecondOrderState::new(
        0.0, &[0.1, -0.1, 0.0], &[0.2, -0.2, 0.0], &[-1.1, 1.1, 0.0],
    );
    let dt = 0.025;
    let method = OperatorGeneralizedAlpha::new(3, dt, 1.0, ImplicitSolveConfig {
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-13,
            relative_tolerance: 1e-12,
            linear_restart: 3,
            max_linear_cycles: 4,
            forcing_minimum: 1e-14,
            forcing_maximum: 1e-12,
            ..NewtonKrylovConfig::default()
        },
        max_newton_iterations: 8,
    });
    let energy = |state: &SecondOrderState| {
        let mut force = [0.0; 3];
        plate.k.spmv(&state.q, &mut force);
        0.5 * (dot(&state.v, &state.v) + dot(&state.q, &force))
    };
    let initial = energy(&state);
    let mut heat = 0.0;
    for _ in 0..64 {
        let before = state.clone();
        method.step(&mut state, &dynamics, &[0.0; 3]).unwrap();
        let mean_velocity: [f64; 3] = core::array::from_fn(|i| {
            0.5 * before.v[i] + 0.5 * state.v[i]
        });
        // Independent original C=M+K quadratic form, not the admission factor.
        let mut viscous_force = [0.0; 3];
        plate.k.spmv(&mean_velocity, &mut viscous_force);
        for i in 0..3 { viscous_force[i] += mean_velocity[i]; }
        let loss = dt * dot(&mean_velocity, &viscous_force);
        assert!(loss >= 0.0);
        assert!((energy(&before) - energy(&state) - loss).abs() < 2e-11 * initial);
        heat += loss;
    }
    assert!(heat > 0.5 * initial);
    assert!((initial - energy(&state) - heat).abs() < 2e-10 * initial);
}

#[test]
fn sparse_admission_honors_factor_budget_and_skips_unneeded_factors() {
    let plate = plate([1.0, -2.0, 0.0, -2.0, 1.0, 0.0, 0.0, 0.0, -1.0], [1.0; 3]);
    let options = LdltOptions { max_factor_bytes: Some(0), ..LdltOptions::default() };
    for alpha in [0.0, 1.0] {
        PlateDynamics::new_with_mass_options(
            &plate, parameters(alpha, 0.0), &NoPlateLoad, budget(&plate), &options, &mut || false,
        ).unwrap();
    }
    assert!(matches!(PlateDynamics::new_with_mass_options(
        &plate, parameters(1.0, 1.0), &NoPlateLoad, budget(&plate), &options, &mut || false,
    ), Err(PlateDynamicsError::DampingFactor(_))));
    assert!(matches!(PlateDynamics::new(
        &plate, parameters(1.0, 1.0), &NoPlateLoad, budget(&plate), &mut || true,
    ), Err(PlateDynamicsError::Cancelled)));
}
