use super::*;

fn controls() -> InertialWallsControl {
    InertialWallsControl {
        midpoint_velocity_bracket_m_s: [[-2.0, 2.0]; 2],
        momentum_tolerance_n_s: 1e-11,
        energy_tolerance_j: 1e-11,
        maximum_iterations: 16,
    }
}

fn walls() -> [Piston; 2] {
    [
        Piston::new(2.0, 0.0).unwrap(),
        Piston::new(4.0, 0.0).unwrap(),
    ]
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-11, "{a:.17e} != {b:.17e}");
}

#[test]
fn g1_coupled_elastic_storage_has_analytic_recoil_and_retains_the_last_trial() {
    // Independent two-port fixture: pressure drops linearly with separation,
    // p_mid = 16 - (8/2) * dt * (w_right - w_left), area=2, dt=1/8.
    // Walls of mass 2 and 4 give p_mid=1024/67, w=[-64/67,32/67].
    // Its exact elastic energy change is p_mid times the negative volume
    // increment. This exercises the shared mechanics, not a replacement EOS.
    let dt = 0.125;
    let mut last = ([0.0; 2], 0);
    let result = solve_coupled(
        walls(),
        Time::new(dt),
        StepControl::default(),
        [InertialWallLoad::default(); 2],
        controls(),
        2.0,
        |speed| {
            last = (speed, last.1 + 1);
            let pressure = 16.0 - 4.0 * dt * (speed[1] - speed[0]);
            let impulse = [2.0 * dt * pressure, -2.0 * dt * pressure];
            let work = [impulse[0] * speed[0], impulse[1] * speed[1]];
            Ok::<_, PistonError>(WallGasExchange {
                gas: last,
                momentum_change_n_s: 0.0,
                energy_change_j: work[0] + work[1],
                boundary_impulse_n_s: impulse,
                boundary_work_j: work,
                displacement_m: speed.map(|w| dt * w),
            })
        },
        |_| false,
    )
    .unwrap();
    close(result.midpoint_velocity_m_s[0], -64.0 / 67.0);
    close(result.midpoint_velocity_m_s[1], 32.0 / 67.0);
    close(result.velocity_after_m_s[0], -128.0 / 67.0);
    close(result.velocity_after_m_s[1], 64.0 / 67.0);
    close(result.combined_momentum_residual_n_s, 0.0);
    close(result.total_energy_residual_j, 0.0);
    assert_eq!(result.gas, last);
    assert_eq!(result.gas_trials, last.1);
    assert!(result.gas_trials >= 4);
}

#[test]
fn g0_zero_momentum_defect_does_not_hide_an_independent_energy_defect() {
    for corrupt_boundary_work in [false, true] {
        let mut calls = 0;
        let error = solve_coupled(
            walls(),
            Time::new(0.1),
            StepControl::default(),
            [InertialWallLoad::default(); 2],
            controls(),
            1.0,
            |_| {
                calls += 1;
                Ok::<_, PistonError>(WallGasExchange {
                    gas: (),
                    momentum_change_n_s: 0.0,
                    energy_change_j: if corrupt_boundary_work { 0.0 } else { 1.0 },
                    boundary_impulse_n_s: [0.0; 2],
                    boundary_work_j: if corrupt_boundary_work {
                        [1.0, -1.0]
                    } else {
                        [0.0; 2]
                    },
                    displacement_m: [0.0; 2],
                })
            },
            |_| false,
        )
        .unwrap_err();
        assert!(
            matches!(error, PistonError::Balance { energy_residual_j, .. } if energy_residual_j == 1.0)
        );
        assert_eq!(calls, 1);
    }
}

#[derive(Debug, PartialEq)]
enum OwnerError {
    Mechanics(PistonError),
    SourceLimit { cell: usize },
}

impl From<PistonError> for OwnerError {
    fn from(value: PistonError) -> Self {
        Self::Mechanics(value)
    }
}

#[test]
fn g4_owner_errors_and_invalid_admission_do_not_publish_or_disguise_a_trial() {
    let mut calls = 0;
    let error = solve_coupled(
        walls(),
        Time::new(0.1),
        StepControl::default(),
        [InertialWallLoad::default(); 2],
        controls(),
        1.0,
        |_| -> Result<WallGasExchange<()>, OwnerError> {
            calls += 1;
            Err(OwnerError::SourceLimit { cell: 7 })
        },
        |_| true,
    )
    .unwrap_err();
    assert_eq!(error, OwnerError::SourceLimit { cell: 7 });
    assert_eq!(calls, 1);
    for area in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        let error = solve_coupled(
            walls(),
            Time::new(0.1),
            StepControl::default(),
            [InertialWallLoad::default(); 2],
            controls(),
            area,
            |_| -> Result<WallGasExchange<()>, PistonError> {
                panic!("invalid area reached the original gas owner")
            },
            |_| false,
        )
        .unwrap_err();
        assert!(matches!(error, PistonError::Invalid(_)));
    }
}
