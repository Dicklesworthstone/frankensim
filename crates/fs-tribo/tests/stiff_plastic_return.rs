//! Manufactured limits of the retained Bingham layer, including shear memory.
//! These tests verify the constitutive update, not measured material parameters.
use fs_tribo::layer_contact::tensor::{
    ShearStiffness, ShearViscosity, TensorLayerRelaxationArm,
};
use fs_tribo::layer_contact::ThirdBodyLayerContact;
use fs_tribo::{FrictionRegime, TangentialLaw, TangentialSolveConfig, TangentialStepInput};

fn input(slip: [f64; 2], dt: f64, compliance: f64) -> TangentialStepInput {
    TangentialStepInput {
        normal_n: 1.0,
        slip_increment_m: slip,
        dt_s: dt,
        contact_temperature_k: Some(300.0),
        nominal_pressure_pa: None,
        compliance_m_per_n: compliance,
        solve: TangentialSolveConfig {
            max_iterations: 64,
            maximum_force_n: 10000.0,
            force_absolute_tolerance_n: 1e-12,
            force_relative_tolerance: 1e-13,
        },
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}

fn tensor_model(
    stiffness: [[f64; 2]; 2],
    yield_force: f64,
    resistance: f64,
) -> ThirdBodyLayerContact {
    ThirdBodyLayerContact::new_tensor(
        ShearStiffness::positive_definite(stiffness).unwrap(),
        vec![],
        yield_force,
        resistance,
        300.0,
        1.0,
    )
    .unwrap()
}

fn isotropic_models(resistance: f64) -> [(&'static str, ThirdBodyLayerContact); 2] {
    [
        (
            "scalar",
            ThirdBodyLayerContact::new(1000.0, 0.5, resistance, 300.0, 1.0).unwrap(),
        ),
        (
            "tensor",
            tensor_model([[1000.0, 0.0], [0.0, 1000.0]], 0.5, resistance),
        ),
    ]
}

fn solve(matrix: [[f64; 2]; 2], rhs: [f64; 2]) -> [f64; 2] {
    let determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    [
        (matrix[1][1] * rhs[0] - matrix[0][1] * rhs[1]) / determinant,
        (matrix[0][0] * rhs[1] - matrix[1][0] * rhs[0]) / determinant,
    ]
}

fn quadratic(matrix: [[f64; 2]; 2], vector: [f64; 2]) -> f64 {
    matrix[0][0] * vector[0] * vector[0]
        + 2.0 * matrix[0][1] * vector[0] * vector[1]
        + matrix[1][1] * vector[1] * vector[1]
}

#[test]
fn g1_stiff_returns_retain_plastic_slip_and_match_series_spring_work() {
    let direction = [0.6, -0.8];
    for resistance in [2.0, 1e-6, 1e-12, 1e-24, 1e-310] {
        for compliance in [0.0, 0.001] {
            let x = input(direction.map(|v| 0.01 * v), 0.01, compliance);
            // Independent series equations: du = (2/k + C) F + dp,
            // F = Y + (c/dt) dp. This stays finite when k*dt/c overflows.
            let series_compliance = 2.0 / 1000.0 + compliance;
            let distance = (0.01 - series_compliance * 0.5)
                / (1.0 + series_compliance * resistance / x.dt_s);
            let force = 0.5 + resistance * (distance / x.dt_s);
            let expected_heat = (0.5 + resistance * (distance / x.dt_s)) * distance;
            for (kind, law) in isotropic_models(resistance) {
                let step = law.advance(&law.initial_state(), &x).unwrap_or_else(|error| {
                    panic!("{kind}, resistance {resistance:e}: {error:?}")
                });
                assert_eq!(step.regime, FrictionRegime::Sliding);
                for (i, axis) in direction.into_iter().enumerate() {
                    close(step.traction_n[i], force * axis, 3e-12);
                    close(step.state.plastic_slip_m()[i], distance * axis, 2e-14);
                    close(step.state.elastic_slip_m()[i], 2.0 * force / 1000.0 * axis, 2e-14);
                }
                close(step.dissipated_j, expected_heat, 1e-14);
                close(step.stored_change_j, 2.0 * force * force / 1000.0, 1e-14);
                // Constitutive force tolerance times plastic distance bounds
                // this residual; no closure deficit is counted as heat.
                close(step.closure_residual_j, 0.0, 2e-14);
            }
        }
    }
}

#[test]
fn g1_stiff_return_remains_continuous_across_the_elastic_yield_boundary() {
    let series_compliance = 2.0 / 1000.0 + 0.001;
    for resistance in [2.0, 1e-24] {
        for offset in [-1e-8, 0.0, 1e-8] {
            let slip = series_compliance * 0.5 * (1.0 + offset);
            let x = input([slip, 0.0], 0.01, 0.001);
            let distance = ((slip - series_compliance * 0.5)
                / (1.0 + series_compliance * resistance / x.dt_s))
                .max(0.0);
            let force = if offset <= 0.0 {
                slip / series_compliance
            } else {
                0.5 + resistance * (distance / x.dt_s)
            };
            for (_, law) in isotropic_models(resistance) {
                let step = law.advance(&law.initial_state(), &x).unwrap();
                close(step.traction_n[0], force, 3e-12);
                close(step.state.plastic_slip_m()[0], distance, 2e-14);
                assert_eq!(step.state.plastic_slip_m()[1], 0.0);
                if offset < 0.0 {
                    assert_eq!(step.regime, FrictionRegime::Sticking);
                    assert_eq!(step.state.plastic_slip_m(), [0.0; 2]);
                } else if offset > 0.0 {
                    assert_eq!(step.regime, FrictionRegime::Sliding);
                    assert!(step.state.plastic_slip_m()[0] > 0.0);
                }
                close(step.closure_residual_j, 0.0, 1e-16);
            }
        }
    }
}

#[test]
fn g1_zero_yield_tensor_is_the_exact_series_maxwell_limit_without_iteration() {
    let stiffness = [[800.0, 70.0], [70.0, 500.0]];
    let determinant = 800.0 * 500.0 - 70.0 * 70.0;
    let x = input([0.003, -0.001], 0.01, 0.0002);
    let compliance = [
        [2.0 * 500.0 / determinant + x.compliance_m_per_n, -2.0 * 70.0 / determinant],
        [-2.0 * 70.0 / determinant, 2.0 * 800.0 / determinant + x.compliance_m_per_n],
    ];
    for resistance in [2.0, 1e-12, 1e-24, 1e-310] {
        let law = tensor_model(stiffness, 0.0, resistance);
        let step = law.advance(&law.initial_state(), &x).unwrap();
        // Solve for dp independently: [I + (c/dt)(2 K^-1 + C I)] dp = du.
        let coefficient = resistance / x.dt_s;
        let matrix = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                coefficient * compliance[i][j] + if i == j { 1.0 } else { 0.0 }
            })
        });
        let plastic = solve(matrix, x.slip_increment_m);
        assert_eq!(step.iterations, 0);
        for i in 0..2 {
            let force = resistance * (plastic[i] / x.dt_s);
            close(step.state.plastic_slip_m()[i], plastic[i], 2e-15);
            close(step.traction_n[i], force, 1e-12 * force.abs() + 1e-320);
        }
        let heat = coefficient * (plastic[0] * plastic[0] + plastic[1] * plastic[1]);
        close(step.dissipated_j, heat, 1e-12 * heat + 1e-320);
        close(step.closure_residual_j, 0.0, 1e-16);
    }
}

#[test]
fn g1_stiff_anisotropic_cycles_keep_memory_and_all_physical_dissipation_channels() {
    let arms = [
        [[400.0, 80.0], [80.0, 65.0]],
        [[64.0, -88.0], [-88.0, 121.0]],
    ];
    let times = [0.02, 0.07];
    let kelvin = [[0.0333, -0.009], [-0.009, 0.0193]];
    let law = ThirdBodyLayerContact::new_tensor(
        ShearStiffness::positive_definite([[800.0, 70.0], [70.0, 500.0]]).unwrap(),
        vec![
            TensorLayerRelaxationArm {
                stiffness: ShearStiffness::from_gram_rows(vec![[20.0, 4.0], [0.0, 7.0]]).unwrap(),
                relaxation_time_s: times[0],
            },
            TensorLayerRelaxationArm {
                stiffness: ShearStiffness::from_gram_rows(vec![[8.0, -11.0]]).unwrap(),
                relaxation_time_s: times[1],
            },
        ],
        0.1,
        1e-20,
        300.0,
        1.0,
    )
    .unwrap()
    .with_kelvin_viscosity(
        ShearViscosity::from_gram_rows(vec![[0.18, -0.07], [0.03, 0.12]]).unwrap(),
    )
    .unwrap();
    let mut state = law.initial_state();
    let (mut work, mut heat) = (0.0, 0.0);
    let (mut plastic_heat, mut relaxation_heat, mut kelvin_heat) = (0.0, 0.0, 0.0);
    for i in 0..160 {
        let angle = 0.29 * i as f64;
        let x = input([0.0005 * angle.cos(), 0.0004 * angle.sin()], 0.0005, 0.0001);
        let step = law.advance(&state, &x).unwrap();
        let plastic: [f64; 2] = core::array::from_fn(|j| {
            step.state.plastic_slip_m()[j] - state.plastic_slip_m()[j]
        });
        let distance = plastic[0].hypot(plastic[1]);
        let channels = law.step_dissipation(&state, &step.state, x.dt_s).unwrap();
        if distance > 0.0 {
            close(
                step.traction_n[0].hypot(step.traction_n[1]),
                0.1 + 1e-20 * (distance / x.dt_s),
                3e-12,
            );
        }
        close(
            channels.plastic_j,
            step.traction_n[0] * plastic[0] + step.traction_n[1] * plastic[1],
            2e-15,
        );
        let expected_relaxation: f64 = (0..2)
            .map(|j| {
                let mid = core::array::from_fn(|axis| {
                    0.5 * (state.relaxation_slip_m()[j][axis]
                        + step.state.relaxation_slip_m()[j][axis])
                });
                x.dt_s / times[j] * quadratic(arms[j], mid)
            })
            .sum();
        let delta = core::array::from_fn(|j| {
            step.state.elastic_slip_m()[j] - state.elastic_slip_m()[j]
        });
        close(channels.relaxation_j, expected_relaxation, 1e-18);
        close(channels.kelvin_j, quadratic(kelvin, delta) / x.dt_s, 1e-18);
        assert!(channels.plastic_j >= 0.0);
        assert!(channels.relaxation_j >= 0.0);
        assert!(channels.kelvin_j >= 0.0);
        close(
            step.dissipated_j,
            channels.plastic_j + channels.relaxation_j + channels.kelvin_j,
            1e-18,
        );
        close(step.closure_residual_j, 0.0, 4e-15);
        work += step.traction_n[0] * step.slip_increment_m[0]
            + step.traction_n[1] * step.slip_increment_m[1];
        heat += step.dissipated_j;
        plastic_heat += channels.plastic_j;
        relaxation_heat += channels.relaxation_j;
        kelvin_heat += channels.kelvin_j;
        state = step.state;
    }
    assert!(plastic_heat > 0.0 && relaxation_heat > 0.0 && kelvin_heat > 0.0);
    close(law.stored_energy_j(&state) + heat, work, 2e-13);
}
