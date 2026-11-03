//! Independent invariant-group nonlinear mechanics and extensive exchanges.
use super::*;

fn one_interval() -> SharedRecoilControl {
    SharedRecoilControl {
        maximum_steps: 1,
        ..control(1e-10)
    }
}

fn group(velocity: f64, linear: f64, quadratic: f64) -> [SharedRecoilCell; 3] {
    core::array::from_fn(|i| {
        let density = [1.0, 2.0, 4.0][i];
        SharedRecoilCell {
            density_kg_m3: density,
            momentum_kg_m2_s: density * velocity,
            volume_m3: [0.25, 0.125, 0.125][i],
            linear_rate_per_s: linear,
            quadratic_per_m: quadratic,
        }
    })
}

// Independent Riccati solution for slip plus conserved center velocity.
// Work follows support kinetic energy and heat follows reduced kinetic energy;
// neither reference calls the production pair map or integrates split flows.
fn reference(
    mass: f64,
    velocity: f64,
    support: SlidingSupport,
    linear: f64,
    quadratic: f64,
    dt: f64,
) -> (f64, f64, f64, f64) {
    let combined = mass + support.mass_kg;
    let reduced = mass * support.mass_kg / combined;
    let slip = velocity - support.velocity_m_s;
    let time = dt * mass / reduced;
    let factor = if linear == 0.0 {
        1.0 / (1.0 + quadratic * slip.abs() * time)
    } else {
        (-linear * time).exp()
            / (1.0 + quadratic * slip.abs() * -(-linear * time).exp_m1() / linear)
    };
    let next_slip = slip * factor;
    let center = (mass * velocity + support.mass_kg * support.velocity_m_s) / combined;
    let gas = center + support.mass_kg / combined * next_slip;
    let wall = center - mass / combined * next_slip;
    let heat = 0.5 * reduced * (slip - next_slip) * (slip + next_slip);
    let work =
        -0.5 * support.mass_kg * (wall - support.velocity_m_s) * (wall + support.velocity_m_s);
    (gas, wall, heat, work)
}

#[test]
fn uniform_nonlinear_matches_riccati_heat_work_mass_order_and_boost() {
    let support = SlidingSupport {
        mass_kg: 0.4,
        velocity_m_s: -0.3,
    };
    let mut workspace = SharedRecoilWorkspace::new(4, || false).unwrap();
    for linear in [0.0, 0.8] {
        for velocity in [-1.0, 2.0] {
            let active = group(velocity, linear, 1.7);
            let expected = reference(1.0, velocity, support, linear, 1.7, 0.4);
            for boost in [0.0, -20.0, 50.0] {
                for reversed in [false, true] {
                    let mut input = [active[0], active[1], active[2], cell(9.0, 1e150, 0.0, 0.0)]
                        .map(|mut c| {
                            c.momentum_kg_m2_s += c.density_kg_m3 * boost;
                            c
                        });
                    if reversed {
                        input.reverse();
                    }
                    let result = workspace
                        .advance(
                            &input,
                            Time::new(0.4),
                            SlidingSupport {
                                velocity_m_s: support.velocity_m_s + boost,
                                ..support
                            },
                            one_interval(),
                            balance(),
                            || false,
                        )
                        .unwrap();
                    assert_eq!(result.accepted_steps, 1);
                    assert_eq!(result.attempted_steps, 1);
                    assert_eq!(result.velocity_error_estimate_m_s, 0.0);
                    assert_eq!(result.energy_error_estimate_j, 0.0);
                    close(result.support.velocity_m_s, expected.1 + boost, 2e-12);
                    for (c, out) in input.iter().zip(result.cells) {
                        if c.quadratic_per_m == 0.0 {
                            assert_eq!(out.momentum, c.momentum_kg_m2_s);
                            assert_eq!(out.gas_impulse_ns, 0.0);
                            assert_eq!(out.slip_dissipation_j, 0.0);
                            assert_eq!(out.work_to_gas_j, 0.0);
                            continue;
                        }
                        let mass = c.density_kg_m3 * c.volume_m3;
                        close(out.momentum / c.density_kg_m3, expected.0 + boost, 2e-12);
                        close(out.slip_dissipation_j, mass * expected.2, 2e-12);
                        close(
                            out.work_to_gas_j,
                            mass * expected.3 + boost * out.gas_impulse_ns,
                            3e-11,
                        );
                        close(out.gas_impulse_ns, mass * (expected.0 - velocity), 2e-12);
                    }
                    close(result.momentum_residual_ns, 0.0, 2e-12);
                    close(result.energy_residual_j, 0.0, 3e-11);
                }
            }
        }
    }
}

#[test]
fn uniform_nonlinear_is_invariant_under_clock_and_mass_partition() {
    let support = SlidingSupport {
        mass_kg: 0.4,
        velocity_m_s: -0.3,
    };
    for linear in [0.0, 0.8] {
        let cells = group(2.0, linear, 1.7);
        let mut workspace = SharedRecoilWorkspace::new(3, || false).unwrap();
        let whole = workspace
            .advance(
                &cells,
                Time::new(0.4),
                support,
                one_interval(),
                balance(),
                || false,
            )
            .unwrap();
        let expected = whole.cells.to_vec();
        let endpoint = whole.support;
        let first = workspace
            .advance(
                &cells,
                Time::new(0.2),
                support,
                one_interval(),
                balance(),
                || false,
            )
            .unwrap();
        assert_eq!(first.accepted_steps, 1);
        let first_cells = first.cells.to_vec();
        let midpoint = first.support;
        let next = core::array::from_fn::<_, 3, _>(|i| SharedRecoilCell {
            momentum_kg_m2_s: first_cells[i].momentum,
            ..cells[i]
        });
        let second = workspace
            .advance(
                &next,
                Time::new(0.2),
                midpoint,
                one_interval(),
                balance(),
                || false,
            )
            .unwrap();
        assert_eq!(second.accepted_steps, 1);
        close(second.support.velocity_m_s, endpoint.velocity_m_s, 8e-15);
        for i in 0..3 {
            close(second.cells[i].momentum, expected[i].momentum, 8e-15);
            close(
                first_cells[i].slip_dissipation_j + second.cells[i].slip_dissipation_j,
                expected[i].slip_dissipation_j,
                8e-15,
            );
            close(
                first_cells[i].work_to_gas_j + second.cells[i].work_to_gas_j,
                expected[i].work_to_gas_j,
                8e-15,
            );
        }
        let expected_heat: f64 = expected.iter().map(|c| c.slip_dissipation_j).sum();
        let expected_work: f64 = expected.iter().map(|c| c.work_to_gas_j).sum();
        for count in [1, 8] {
            let partition = vec![cell(1.0 / count as f64, 2.0, linear, 1.7); count];
            let mut repartitioned = SharedRecoilWorkspace::new(count, || false).unwrap();
            let result = repartitioned
                .advance(
                    &partition,
                    Time::new(0.4),
                    support,
                    one_interval(),
                    balance(),
                    || false,
                )
                .unwrap();
            assert_eq!(result.accepted_steps, 1);
            close(result.support.velocity_m_s, endpoint.velocity_m_s, 8e-15);
            for out in result.cells {
                close(out.momentum * count as f64, expected[0].momentum, 8e-15);
                close(out.slip_dissipation_j * count as f64, expected_heat, 8e-15);
                close(out.work_to_gas_j * count as f64, expected_work, 8e-15);
            }
        }
    }
}

#[test]
fn stiff_uniform_nonlinear_finishes_in_one_interval_and_heterogeneity_is_not_merged() {
    let support = SlidingSupport {
        mass_kg: 1e-6,
        velocity_m_s: -0.1,
    };
    let mut workspace = SharedRecoilWorkspace::new(3, || false).unwrap();
    for linear in [0.0, 1e7] {
        let cells = group(0.2, linear, 1e8);
        let expected = reference(1.0, 0.2, support, linear, 1e8, 10.0);
        let result = workspace
            .advance(
                &cells,
                Time::new(10.0),
                support,
                one_interval(),
                balance(),
                || false,
            )
            .unwrap();
        assert_eq!(result.accepted_steps, 1);
        assert_eq!(result.attempted_steps, 1);
        close(result.support.velocity_m_s, expected.1, 3e-11);
        for (c, out) in cells.iter().zip(result.cells) {
            let mass = c.density_kg_m3 * c.volume_m3;
            close(out.momentum / c.density_kg_m3, expected.0, 3e-16);
            close(out.slip_dissipation_j, mass * expected.2, 3e-17);
            close(out.work_to_gas_j, mass * expected.3, 3e-17);
        }
    }
    // A subnormal mass fraction can carry resolved energy even though the
    // aggregate specific energy overflows. All actual endpoints remain finite.
    let tiny = f64::from_bits(1);
    let extreme = [tiny, 1e-12].map(|volume| SharedRecoilCell {
        density_kg_m3: 1.0,
        momentum_kg_m2_s: 1e160,
        volume_m3: volume,
        linear_rate_per_s: 0.0,
        quadratic_per_m: 1e-160,
    });
    let mut extreme_workspace = SharedRecoilWorkspace::new(2, || false).unwrap();
    let result = extreme_workspace
        .advance(
            &extreme,
            Time::new(1.0),
            SlidingSupport {
                mass_kg: 1e-12,
                velocity_m_s: 0.0,
            },
            one_interval(),
            RecoilControl {
                momentum_tolerance_ns: 1e140,
                ..balance()
            },
            || false,
        )
        .unwrap();
    let tiny_scale = tiny * 1e160 * 1e160;
    close(
        result.cells[0].slip_dissipation_j,
        tiny_scale * (2.0 / 9.0),
        1e-15,
    );
    close(
        result.cells[0].work_to_gas_j,
        tiny_scale * (-1.0 / 18.0),
        1e-15,
    );
    assert_eq!(result.accepted_steps, 1);

    // One ULP in any part of the invariant descriptor must keep the bounded
    // general solver. Its one-step budget is intentionally insufficient here.
    for changed in 0..3 {
        let mut cells = group(2.0, 0.8, 1.7);
        let value = match changed {
            0 => &mut cells[0].momentum_kg_m2_s,
            1 => &mut cells[0].linear_rate_per_s,
            _ => &mut cells[0].quadratic_per_m,
        };
        *value = f64::from_bits(value.to_bits() + 1);
        let result = workspace.advance(
            &cells,
            Time::new(0.4),
            SlidingSupport {
                mass_kg: 0.4,
                velocity_m_s: -0.3,
            },
            one_interval(),
            balance(),
            || false,
        );
        assert!(
            matches!(
                result,
                Err(WallExchangeError::Invalid(
                    "shared support integration work exhausted"
                ))
            ),
            "{result:?}"
        );
    }
}

#[test]
fn uniform_nonlinear_retries_exactly_after_every_cancellation_boundary() {
    let cells = group(2.0, 0.8, 1.7);
    let support = SlidingSupport {
        mass_kg: 0.4,
        velocity_m_s: -0.3,
    };
    let mut workspace = SharedRecoilWorkspace::new(3, || false).unwrap();
    let pointers = (
        workspace.accepted.as_ptr(),
        workspace.coarse.as_ptr(),
        workspace.fine.as_ptr(),
        workspace.updates.as_ptr(),
    );
    let mut visits = 0;
    let result = workspace
        .advance(
            &cells,
            Time::new(0.4),
            support,
            one_interval(),
            balance(),
            || {
                visits += 1;
                false
            },
        )
        .unwrap();
    let expected = result.cells.to_vec();
    let endpoint = result.support;
    for cut in 0..visits {
        let mut current = 0;
        let refused = workspace.advance(
            &cells,
            Time::new(0.4),
            support,
            one_interval(),
            balance(),
            || {
                let stop = current == cut;
                current += 1;
                stop
            },
        );
        assert!(matches!(refused, Err(WallExchangeError::Cancelled)));
        let retry = workspace
            .advance(
                &cells,
                Time::new(0.4),
                support,
                one_interval(),
                balance(),
                || false,
            )
            .unwrap();
        assert_eq!(retry.cells, expected.as_slice());
        assert_eq!(retry.support, endpoint);
    }
    assert_eq!(
        pointers,
        (
            workspace.accepted.as_ptr(),
            workspace.coarse.as_ptr(),
            workspace.fine.as_ptr(),
            workspace.updates.as_ptr(),
        )
    );
}
