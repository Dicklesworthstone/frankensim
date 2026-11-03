//! Exact collective G1/G3/G4 oracles, including initially-comoving weak heat.
use super::*;

fn one_interval() -> SharedRecoilControl {
    SharedRecoilControl {
        maximum_steps: 1,
        ..control(1e-10)
    }
}

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(expected != 0.0);
    assert!(
        (actual / expected - 1.0).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn common_rate_matches_high_precision_three_cell_heat_work_and_frame_changes() {
    let cells = [
        cell(1.2, -1.0, 2.0, 0.0),
        cell(0.7, 4.0, 2.0, 0.0),
        cell(2.1, 0.5, 2.0, 0.0),
    ];
    let support = SlidingSupport {
        mass_kg: 0.4,
        velocity_m_s: 0.2,
    };
    // Independent 100-digit decimal integration of the two exponential modes
    // and their squared-slip cross term, rounded only here.
    let velocity = [
        0.11971924595987800357,
        1.6256903055208884868,
        0.57151056382818114853,
    ];
    let heat = [
        1.3657886485022222361,
        3.7318345606674291947,
        0.01932491652379493928,
    ];
    let work = [
        0.77438826721414329548,
        -0.94316130001996089238,
        0.09978045731936073382,
    ];
    let final_support = 0.62045376736086010761;
    let mut workspace = SharedRecoilWorkspace::new(3, || false).unwrap();
    for boost in [0.0, -20.0, 50.0] {
        for reversed in [false, true] {
            let mut input = cells.map(|mut c| {
                c.momentum_kg_m2_s += c.density_kg_m3 * boost;
                c
            });
            if reversed {
                input.reverse();
            }
            let result = workspace
                .advance(
                    &input,
                    Time::new(0.6),
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
            close(result.support.velocity_m_s, final_support + boost, 3e-13);
            for j in 0..3 {
                let i = if reversed { 2 - j } else { j };
                let out = result.cells[j];
                close(
                    out.momentum / cells[i].density_kg_m3,
                    velocity[i] + boost,
                    3e-13,
                );
                close(out.slip_dissipation_j, heat[i], 5e-13);
                close(
                    out.work_to_gas_j,
                    work[i] + boost * out.gas_impulse_ns,
                    1e-11,
                );
            }
            close(result.momentum_residual_ns, 0.0, 2e-12);
            close(result.energy_residual_j, 0.0, 5e-11);
        }
    }
    // Both sides of the short-moment evaluation boundary use the same
    // independent two-by-two matrix exponential and exact squared-slip integral.
    let two = [cell(1.0, 2.0, 1.0, 0.0), cell(2.0, -0.5, 1.0, 0.0)];
    let support = SlidingSupport {
        mass_kg: 3.0,
        velocity_m_s: 0.1,
    };
    let mut workspace = SharedRecoilWorkspace::new(2, || false).unwrap();
    for dt in [0.249999999999, 0.250000000001] {
        let reference = linear_reference(two, support, dt);
        let result = workspace
            .advance(
                &two,
                Time::new(dt),
                support,
                one_interval(),
                balance(),
                || false,
            )
            .unwrap();
        for i in 0..2 {
            close(
                result.cells[i].momentum / two[i].density_kg_m3,
                reference.0[i],
                2e-14,
            );
            close(result.cells[i].slip_dissipation_j, reference.2[i], 5e-14);
            close(result.cells[i].work_to_gas_j, reference.3[i], 5e-14);
        }
    }
}

#[test]
fn initially_comoving_cells_retain_cubic_heat_and_sub_ulp_inertia_contrast() {
    // Independent 130-digit decimal references. Subtracting bulk kinetic
    // energies, modal heat squares, or q-1 in f64 erases these channels.
    let cases = [
        (
            [cell(1.0, 0.0, 1.0, 0.0), cell(2.0, 1.0, 1.0, 0.0)],
            SlidingSupport {
                mass_kg: 3.0,
                velocity_m_s: 0.0,
            },
            1e-8,
            [
                3.33333330000000019444e-17,
                6.66666660000000044444e-9,
                1.48148144814814855556e-25,
                1.481481453703703999999998e-25,
                -6.66666654814814942593e-17,
            ],
        ),
        (
            [cell(0.5, 0.0, 1.0, 0.0), cell(0.5, 1.0, 1.0, 0.0)],
            SlidingSupport {
                mass_kg: 1e20,
                velocity_m_s: 0.0,
            },
            1.0,
            [
                1.32120558828557678404e-21,
                3.16060279414278839201e-21,
                1.01038619942792668915e-42,
                1.44678225105718594661e-42,
                -4.99470501117160060876e-22,
            ],
        ),
    ];
    let mut workspace = SharedRecoilWorkspace::new(2, || false).unwrap();
    for (cells, support, dt, expected) in cases {
        let result = workspace
            .advance(
                &cells,
                Time::new(dt),
                support,
                one_interval(),
                balance(),
                || false,
            )
            .unwrap();
        relative(
            result.cells[0].momentum / cells[0].density_kg_m3,
            expected[0],
            3e-14,
        );
        relative(result.support.velocity_m_s, expected[1], 3e-14);
        relative(result.cells[0].slip_dissipation_j, expected[2], 8e-14);
        relative(result.cells[0].work_to_gas_j, expected[3], 8e-14);
        relative(result.cells[1].work_to_gas_j, expected[4], 8e-14);
        assert_eq!(result.accepted_steps, 1);
    }
}

#[test]
fn fully_relaxed_collective_mode_handles_stiff_phase_without_extra_work_or_heat() {
    let tiny_gas = [cell(0.01, 0.1, 1e308, 0.0), cell(0.01, -0.1, 1e308, 0.0)];
    let heavy = SlidingSupport {
        mass_kg: 1e308,
        velocity_m_s: 0.0,
    };
    let mut workspace = SharedRecoilWorkspace::new(2, || false).unwrap();
    let saturated = workspace
        .advance(
            &tiny_gas,
            Time::new(10.0),
            heavy,
            one_interval(),
            balance(),
            || false,
        )
        .unwrap();
    assert_eq!(saturated.support, heavy);
    for output in saturated.cells {
        assert_eq!(output.momentum, 0.0);
        assert_eq!(output.work_to_gas_j, 0.0);
        relative(output.slip_dissipation_j, 0.00005, 1e-14);
    }
    for boost in [0.0, 100.0] {
        let cells = [
            cell(1.0, -2.0 + boost, 1e200, 0.0),
            cell(2.0, 1.0 + boost, 1e200, 0.0),
            cell(1.0, 1e150, 0.0, 0.0),
        ];
        let support = SlidingSupport {
            mass_kg: 1e-10,
            velocity_m_s: boost,
        };
        let mut workspace = SharedRecoilWorkspace::new(3, || false).unwrap();
        let result = workspace
            .advance(
                &cells,
                Time::new(1e200),
                support,
                one_interval(),
                balance(),
                || false,
            )
            .unwrap();
        assert_eq!(result.accepted_steps, 1);
        close(result.support.velocity_m_s, boost, 1e-14);
        for i in 0..2 {
            close(
                result.cells[i].momentum / cells[i].density_kg_m3,
                boost,
                1e-14,
            );
            close(result.cells[i].slip_dissipation_j, [2.0, 1.0][i], 2e-14);
            close(
                result.cells[i].work_to_gas_j,
                boost * result.cells[i].gas_impulse_ns,
                1e-13,
            );
        }
        assert_eq!(result.cells[2].momentum, cells[2].momentum_kg_m2_s);
        assert_eq!(result.cells[2].work_to_gas_j, 0.0);
        assert_eq!(result.cells[2].slip_dissipation_j, 0.0);
    }
}

#[test]
fn exact_collective_candidate_retries_after_every_cancellation_boundary() {
    let cells = [
        cell(1.0, 2.0, 1.3, 0.0),
        cell(2.0, -1.0, 1.3, 0.0),
        cell(0.3, 0.6, 1.3, 0.0),
    ];
    let support = SlidingSupport {
        mass_kg: 0.4,
        velocity_m_s: -0.2,
    };
    let mut workspace = SharedRecoilWorkspace::new(3, || false).unwrap();
    let pointers = (
        workspace.accepted.as_ptr(),
        workspace.coarse.as_ptr(),
        workspace.fine.as_ptr(),
        workspace.updates.as_ptr(),
    );
    let mut polls = 0;
    let result = workspace
        .advance(
            &cells,
            Time::new(0.3),
            support,
            one_interval(),
            balance(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let expected = result.cells.to_vec();
    let endpoint = result.support;
    for cut in 0..polls {
        let mut visits = 0;
        let result = workspace.advance(
            &cells,
            Time::new(0.3),
            support,
            one_interval(),
            balance(),
            || {
                let stop = visits == cut;
                visits += 1;
                stop
            },
        );
        assert!(matches!(result, Err(WallExchangeError::Cancelled)));
        let retry = workspace
            .advance(
                &cells,
                Time::new(0.3),
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
            workspace.updates.as_ptr()
        )
    );
}
