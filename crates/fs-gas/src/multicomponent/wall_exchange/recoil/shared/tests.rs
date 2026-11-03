//! G0/G1/G3/G4 checks against independent simultaneous solutions.
use super::super::super::WallControl;
use super::*;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}
fn cell(mass: f64, velocity: f64, linear: f64, quadratic: f64) -> SharedRecoilCell {
    SharedRecoilCell {
        density_kg_m3: mass,
        momentum_kg_m2_s: mass * velocity,
        volume_m3: 1.0,
        linear_rate_per_s: linear,
        quadratic_per_m: quadratic,
    }
}
fn control(tolerance: f64) -> SharedRecoilControl {
    SharedRecoilControl {
        absolute_velocity_tolerance_m_s: tolerance,
        absolute_energy_tolerance_j: tolerance,
        relative_tolerance: tolerance,
        maximum_steps: 65_536,
    }
}
fn balance() -> RecoilControl {
    RecoilControl {
        wall: WallControl {
            absolute_energy_tolerance_j: 1e-9,
            relative_energy_tolerance: 1e-11,
            maximum_iterations: 80,
        },
        momentum_tolerance_ns: 1e-9,
    }
}

// Independent exponential of the two-slip star generator, followed by exact
// integration of each edge's squared slip. No pair flow or Strang step is used.
fn linear_reference(
    cells: [SharedRecoilCell; 2],
    support: SlidingSupport,
    dt: f64,
) -> ([f64; 2], f64, [f64; 2], [f64; 2]) {
    let mass = cells.map(|c| c.density_kg_m3 * c.volume_m3);
    let rate = cells.map(|c| c.linear_rate_per_s);
    let velocity = cells.map(|c| c.momentum_kg_m2_s / c.density_kg_m3);
    let slip = velocity.map(|u| u - support.velocity_m_s);
    let r = mass.map(|m| m / support.mass_kg);
    let a = [
        [rate[0] * (1.0 + r[0]), r[1] * rate[1]],
        [r[0] * rate[0], rate[1] * (1.0 + r[1])],
    ];
    let gap = ((a[0][0] - a[1][1]).powi(2) + 4.0 * a[0][1] * a[1][0]).sqrt();
    let slow = (a[0][0] + a[1][1] - gap) * 0.5;
    let fast = (a[0][0] + a[1][1] + gap) * 0.5;
    let a_slip = [
        a[0][0] * slip[0] + a[0][1] * slip[1],
        a[1][0] * slip[0] + a[1][1] * slip[1],
    ];
    let cm = core::array::from_fn::<_, 2, _>(|i| (fast * slip[i] - a_slip[i]) / gap);
    let cp = core::array::from_fn::<_, 2, _>(|i| (a_slip[i] - slow * slip[i]) / gap);
    let final_slip = core::array::from_fn::<_, 2, _>(|i| {
        cm[i] * (-slow * dt).exp() + cp[i] * (-fast * dt).exp()
    });
    let total_mass = support.mass_kg + mass[0] + mass[1];
    let center =
        (support.mass_kg * support.velocity_m_s + mass[0] * velocity[0] + mass[1] * velocity[1])
            / total_mass;
    let final_support = center - (mass[0] * final_slip[0] + mass[1] * final_slip[1]) / total_mass;
    let final_velocity = final_slip.map(|s| final_support + s);
    let exposure = |rate: f64| -(-rate * dt).exp_m1() / rate;
    let heat = core::array::from_fn::<_, 2, _>(|i| {
        mass[i]
            * rate[i]
            * (cm[i] * cm[i] * exposure(2.0 * slow)
                + 2.0 * cm[i] * cp[i] * exposure(slow + fast)
                + cp[i] * cp[i] * exposure(2.0 * fast))
    });
    let work = core::array::from_fn::<_, 2, _>(|i| {
        0.5 * mass[i] * (final_velocity[i] - velocity[i]) * (final_velocity[i] + velocity[i])
            + heat[i]
    });
    (final_velocity, final_support, heat, work)
}

#[test]
fn unequal_linear_star_matches_independent_velocities_heat_and_work_and_refines() {
    let cells = [cell(2.0, 4.0, 0.7, 0.0), cell(5.0, -1.2, 2.3, 0.0)];
    let support = SlidingSupport {
        mass_kg: 3.0,
        velocity_m_s: -0.4,
    };
    let reference = linear_reference(cells, support, 0.9);
    let mut errors = Vec::new();
    let mut steps = Vec::new();
    for tolerance in [1e-4, 1e-7, 1e-9] {
        let mut workspace = SharedRecoilWorkspace::new(2, || false).unwrap();
        let update = workspace
            .advance(
                &cells,
                Time::new(0.9),
                support,
                control(tolerance),
                balance(),
                || false,
            )
            .unwrap();
        let mut error = (update.support.velocity_m_s - reference.1).abs();
        for (i, output) in update.cells.iter().enumerate() {
            let velocity = output.momentum / cells[i].density_kg_m3;
            close(velocity, reference.0[i], 40.0 * tolerance);
            close(output.slip_dissipation_j, reference.2[i], 80.0 * tolerance);
            close(output.work_to_gas_j, reference.3[i], 80.0 * tolerance);
            error = error
                .max((velocity - reference.0[i]).abs())
                .max((output.slip_dissipation_j - reference.2[i]).abs())
                .max((output.work_to_gas_j - reference.3[i]).abs());
        }
        close(update.momentum_residual_ns, 0.0, 1e-11);
        close(update.energy_residual_j, 0.0, 1e-10);
        errors.push(error);
        steps.push(update.accepted_steps);
    }
    assert!(errors[1] < errors[0] / 100.0, "{errors:?}");
    assert!(errors[2] < errors[1] / 20.0, "{errors:?}");
    assert!(steps[2] > steps[1] && steps[1] > steps[0], "{steps:?}");
}

#[test]
fn nonlinear_symmetric_star_has_stationary_support_and_correct_local_heat() {
    for (linear, quadratic) in [(0.0, 1.1), (0.8, 0.0), (0.8, 1.1)] {
        let cells = [
            cell(1.5, 2.0, linear, quadratic),
            cell(1.5, -2.0, linear, quadratic),
        ];
        let support = SlidingSupport {
            mass_kg: 0.7,
            velocity_m_s: 0.0,
        };
        let dt = 0.6;
        let factor = if linear == 0.0 {
            1.0 / (1.0 + quadratic * 2.0 * dt)
        } else {
            (-linear * dt).exp() / (1.0 + quadratic * 2.0 * (-(-linear * dt).exp_m1()) / linear)
        };
        let expected_heat = 0.5 * 1.5 * 4.0 * (1.0 - factor * factor);
        let mut workspace = SharedRecoilWorkspace::new(2, || false).unwrap();
        let result = workspace
            .advance(
                &cells,
                Time::new(dt),
                support,
                control(1e-8),
                balance(),
                || false,
            )
            .unwrap();
        close(result.support.velocity_m_s, 0.0, 5e-8);
        for (i, output) in result.cells.iter().enumerate() {
            close(output.momentum, cells[i].momentum_kg_m2_s * factor, 1e-7);
            close(output.slip_dissipation_j, expected_heat, 1e-7);
            close(output.work_to_gas_j, 0.0, 1e-7);
            assert!(output.slip_dissipation_j >= 0.0);
        }
        close(result.energy_residual_j, 0.0, 2e-11);
    }
}

#[test]
fn initially_comoving_cell_receives_indirect_impulse_and_physical_slip_heat() {
    let cells = [cell(1.0, 0.0, 1.0, 0.0), cell(1.0, 1.0, 1.0, 0.0)];
    let support = SlidingSupport {
        mass_kg: 1.0,
        velocity_m_s: 0.0,
    };
    let dt = 0.7;
    let reference = linear_reference(cells, support, dt);
    let mut workspace = SharedRecoilWorkspace::new(2, || false).unwrap();
    let result = workspace
        .advance(
            &cells,
            Time::new(dt),
            support,
            control(1e-9),
            balance(),
            || false,
        )
        .unwrap();
    assert!(result.cells[0].gas_impulse_ns > 0.0);
    assert!(result.cells[0].slip_dissipation_j > 0.0);
    assert!(result.cells[0].momentum < result.support.velocity_m_s);
    for i in 0..2 {
        close(result.cells[i].momentum, reference.0[i], 1e-8);
        close(result.cells[i].slip_dissipation_j, reference.2[i], 1e-8);
        close(result.cells[i].work_to_gas_j, reference.3[i], 1e-8);
    }
}

#[test]
fn nonlinear_shared_exchange_is_boost_covariant_and_converges_across_order_and_clock_partition() {
    let cells = [
        cell(1.0, 2.0, 0.8, 0.3),
        cell(2.5, -0.5, 1.2, 0.1),
        cell(0.8, 0.9, 0.2, 0.9),
    ];
    let support = SlidingSupport {
        mass_kg: 2.0,
        velocity_m_s: 0.2,
    };
    let mut workspace = SharedRecoilWorkspace::new(3, || false).unwrap();
    let original = workspace
        .advance(
            &cells,
            Time::new(0.3),
            support,
            control(1e-8),
            balance(),
            || false,
        )
        .unwrap();
    let expected = original.cells.to_vec();
    let expected_support = original.support;
    for boost in [-20.0, 50.0] {
        let shifted = cells.map(|mut c| {
            c.momentum_kg_m2_s += c.density_kg_m3 * boost;
            c
        });
        let result = workspace
            .advance(
                &shifted,
                Time::new(0.3),
                SlidingSupport {
                    velocity_m_s: support.velocity_m_s + boost,
                    ..support
                },
                control(1e-8),
                balance(),
                || false,
            )
            .unwrap();
        close(
            result.support.velocity_m_s,
            expected_support.velocity_m_s + boost,
            2e-9,
        );
        for i in 0..3 {
            close(
                result.cells[i].momentum,
                expected[i].momentum + cells[i].density_kg_m3 * boost,
                3e-9,
            );
            close(
                result.cells[i].gas_impulse_ns,
                expected[i].gas_impulse_ns,
                3e-9,
            );
            close(
                result.cells[i].slip_dissipation_j,
                expected[i].slip_dissipation_j,
                2e-9,
            );
            close(
                result.cells[i].work_to_gas_j,
                expected[i].work_to_gas_j + boost * expected[i].gas_impulse_ns,
                2e-7,
            );
        }
    }
    let reversed = [cells[2], cells[1], cells[0]];
    let result = workspace
        .advance(
            &reversed,
            Time::new(0.3),
            support,
            control(1e-8),
            balance(),
            || false,
        )
        .unwrap();
    for i in 0..3 {
        close(result.cells[2 - i].momentum, expected[i].momentum, 4e-8);
        close(
            result.cells[2 - i].slip_dissipation_j,
            expected[i].slip_dissipation_j,
            4e-8,
        );
    }
    let half = workspace
        .advance(
            &cells,
            Time::new(0.15),
            support,
            control(5e-9),
            balance(),
            || false,
        )
        .unwrap();
    let half_cells = half.cells.to_vec();
    let half_support = half.support;
    let next = core::array::from_fn::<_, 3, _>(|i| SharedRecoilCell {
        momentum_kg_m2_s: half_cells[i].momentum,
        ..cells[i]
    });
    let second = workspace
        .advance(
            &next,
            Time::new(0.15),
            half_support,
            control(5e-9),
            balance(),
            || false,
        )
        .unwrap();
    close(
        second.support.velocity_m_s,
        expected_support.velocity_m_s,
        2e-8,
    );
    for i in 0..3 {
        close(second.cells[i].momentum, expected[i].momentum, 3e-8);
        close(
            half_cells[i].slip_dissipation_j + second.cells[i].slip_dissipation_j,
            expected[i].slip_dissipation_j,
            4e-8,
        );
        close(
            half_cells[i].work_to_gas_j + second.cells[i].work_to_gas_j,
            expected[i].work_to_gas_j,
            4e-8,
        );
    }
}

#[test]
fn single_active_edge_preserves_original_exact_pair_and_disconnected_cells() {
    let support = SlidingSupport {
        mass_kg: 0.8,
        velocity_m_s: -0.3,
    };
    for (linear, quadratic) in [(1e-12, 0.0), (0.0, 2.3), (1.7, 2.3), (1e5, 2.3)] {
        let active = SharedRecoilCell {
            density_kg_m3: 1.5,
            momentum_kg_m2_s: 1.65,
            volume_m3: 0.4,
            linear_rate_per_s: linear,
            quadratic_per_m: quadratic,
        };
        let cells = [
            cell(1.0, 1e150, 0.0, 0.0),
            active,
            cell(2.0, -9.0, 0.0, 0.0),
        ];
        let exact = mechanics::advance(
            active.momentum_kg_m2_s,
            active.density_kg_m3,
            active.momentum_kg_m2_s / active.density_kg_m3,
            active.volume_m3,
            0.2,
            support,
            RecoilLaw {
                linear_rate_per_s: linear,
                quadratic_per_m: quadratic,
                gas_heat_fraction: 1.0,
            },
        )
        .unwrap();
        let mut workspace = SharedRecoilWorkspace::new(3, || false).unwrap();
        let pointers = (
            workspace.accepted.as_ptr(),
            workspace.coarse.as_ptr(),
            workspace.fine.as_ptr(),
            workspace.updates.as_ptr(),
        );
        let output = workspace
            .advance(
                &cells,
                Time::new(0.2),
                support,
                control(1e-9),
                balance(),
                || false,
            )
            .unwrap();
        assert_eq!(output.accepted_steps, 1);
        assert_eq!(output.attempted_steps, 1);
        assert_eq!(output.cells[1].momentum, exact.momentum);
        assert_eq!(output.support.velocity_m_s, exact.support_velocity_m_s);
        assert_eq!(output.cells[1].work_to_gas_j, exact.work_to_gas_j);
        assert_eq!(output.cells[1].slip_dissipation_j, exact.dissipation_j);
        assert_eq!(output.cells[0].momentum, cells[0].momentum_kg_m2_s);
        assert_eq!(output.cells[2].momentum, cells[2].momentum_kg_m2_s);
        assert_eq!(output.cells[0].slip_dissipation_j, 0.0);
        assert_eq!(
            pointers,
            (
                workspace.accepted.as_ptr(),
                workspace.coarse.as_ptr(),
                workspace.fine.as_ptr(),
                workspace.updates.as_ptr()
            )
        );
    }
    let mut workspace = SharedRecoilWorkspace::new(0, || false).unwrap();
    let empty = workspace
        .advance(
            &[],
            Time::new(5.0),
            support,
            control(1e-8),
            balance(),
            || false,
        )
        .unwrap();
    assert_eq!(empty.support, support);
    assert_eq!(empty.accepted_steps, 0);
    assert_eq!(empty.energy_residual_j, 0.0);
}

#[test]
fn stiff_shared_relaxation_resolves_rates_or_refuses_instead_of_false_heat_convergence() {
    let cells = [cell(2.0, 4.0, 100.0, 0.0), cell(5.0, -1.2, 100.0, 0.0)];
    let support = SlidingSupport {
        mass_kg: 3.0,
        velocity_m_s: -0.4,
    };
    let reference = linear_reference(cells, support, 0.1);
    let mut workspace = SharedRecoilWorkspace::new(2, || false).unwrap();
    let result = workspace
        .advance(
            &cells,
            Time::new(0.1),
            support,
            control(1e-8),
            balance(),
            || false,
        )
        .unwrap();
    assert!(result.accepted_steps >= 67);
    for i in 0..2 {
        close(
            result.cells[i].momentum / cells[i].density_kg_m3,
            reference.0[i],
            3e-7,
        );
        close(result.cells[i].slip_dissipation_j, reference.2[i], 1e-6);
        close(result.cells[i].work_to_gas_j, reference.3[i], 1e-6);
        assert!((-1.2..=4.0).contains(&(result.cells[i].momentum / cells[i].density_kg_m3)));
    }
    let refused = workspace.advance(
        &cells,
        Time::new(100.0),
        support,
        SharedRecoilControl {
            maximum_steps: 1,
            ..control(1e-3)
        },
        balance(),
        || false,
    );
    assert!(matches!(
        refused,
        Err(WallExchangeError::Invalid(
            "shared support integration work exhausted"
        ))
    ));
}

#[test]
fn cancellation_at_each_boundary_and_input_refusals_leave_exact_retry() {
    let cells = [cell(1.0, 1.0, 0.4, 0.1), cell(2.0, -0.5, 0.7, 0.0)];
    let support = SlidingSupport {
        mass_kg: 3.0,
        velocity_m_s: 0.2,
    };
    let mut workspace = SharedRecoilWorkspace::new(2, || false).unwrap();
    let mut polls = 0;
    let original = workspace
        .advance(
            &cells,
            Time::new(0.01),
            support,
            control(1e-5),
            balance(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let expected = original.cells.to_vec();
    let expected_support = original.support;
    for cut in 0..polls {
        let mut seen = 0;
        let cancelled = workspace.advance(
            &cells,
            Time::new(0.01),
            support,
            control(1e-5),
            balance(),
            || {
                let result = seen == cut;
                seen += 1;
                result
            },
        );
        assert!(
            matches!(cancelled, Err(WallExchangeError::Cancelled)),
            "cut={cut}"
        );
        let retry = workspace
            .advance(
                &cells,
                Time::new(0.01),
                support,
                control(1e-5),
                balance(),
                || false,
            )
            .unwrap();
        assert_eq!(retry.cells, expected.as_slice());
        assert_eq!(retry.support, expected_support);
    }
    assert!(SharedRecoilWorkspace::new(2, || true).is_err());
    let mut invalid = cells;
    invalid[0].quadratic_per_m = -1.0;
    assert!(
        workspace
            .advance(
                &invalid,
                Time::new(0.01),
                support,
                control(1e-5),
                balance(),
                || false
            )
            .is_err()
    );
    assert!(
        workspace
            .advance(
                &cells[..1],
                Time::new(0.01),
                support,
                control(1e-5),
                balance(),
                || false
            )
            .is_err()
    );
    assert!(
        workspace
            .advance(
                &cells,
                Time::new(-0.01),
                support,
                control(1e-5),
                balance(),
                || false
            )
            .is_err()
    );
    assert!(
        workspace
            .advance(
                &cells,
                Time::new(0.01),
                support,
                SharedRecoilControl {
                    maximum_steps: 0,
                    ..control(1e-5)
                },
                balance(),
                || false
            )
            .is_err()
    );
    let zero = workspace
        .advance(
            &cells,
            Time::new(0.0),
            support,
            control(1e-5),
            balance(),
            || false,
        )
        .unwrap();
    assert_eq!(zero.support, support);
    assert_eq!(zero.cells[0].momentum, cells[0].momentum_kg_m2_s);
    assert_eq!(zero.accepted_steps, 0);
}
