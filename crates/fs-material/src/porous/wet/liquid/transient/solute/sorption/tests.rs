use super::super::transport::{DiluteSoluteBoundary, transport, transport_with_storage};
use super::*;

fn controls<const N: usize>() -> DiluteSoluteControls<N> {
    DiluteSoluteControls {
        absolute_equation_tolerance_mol_m2: [1e-11; N],
        relative_equation_tolerance: 1e-12,
        amount_balance_tolerance_mol_m2: [1e-11; N],
        water_balance_tolerance_m: 1e-13,
        maximum_cell_species: 131072,
        maximum_storage_bytes: 64 * 1024 * 1024,
    }
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e}");
}

#[allow(clippy::too_many_arguments)]
fn solve<const N: usize>(
    dt: f64,
    before: &[f64],
    after: &[f64],
    widths: &[f64],
    water: &[f64],
    mobile: &[[f64; N]],
    bound: &[[f64; N]],
    laws: &[[LinearSoluteSorption; N]],
    d: &[[f64; N]],
    tau: &[f64],
    boundaries: [DiluteSoluteBoundary<N>; 2],
    ctl: DiluteSoluteControls<N>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<SorbingDiluteSoluteStep<N>> {
    let mut storage = SorptionStorage::new(mobile.len(), dt, bound, laws, ctl, &mut cancelled)?;
    let step = transport_with_storage(
        [0.0, dt],
        dt,
        before,
        after,
        widths,
        water,
        mobile,
        d,
        tau,
        boundaries,
        ctl,
        &mut storage,
        &mut cancelled,
    )?;
    storage.finish(mobile, step, &mut cancelled)
}

#[test]
fn g1_closed_kinetic_exchange_matches_independent_two_pool_relaxation() {
    let (v, b, k) = (0.25, 1.0, 0.7);
    let (a0, s0) = (0.2, 0.1);
    let equilibrium = b * (a0 + s0) / (v + b);
    for dt in [1e-6, 0.3, 1e6] {
        let result = solve(
            dt,
            &[v],
            &[v],
            &[1.0],
            &[0.0; 2],
            &[[a0]],
            &[[s0]],
            &[[LinearSoluteSorption::Kinetic {
                capacity_m: b,
                rate_s_inverse: k,
            }]],
            &[[0.0]],
            &[1.0],
            [DiluteSoluteBoundary::Outflow; 2],
            controls(),
            || false,
        )
        .unwrap();
        let expected_s = equilibrium + (s0 - equilibrium) / (1.0 + dt * k * (1.0 + b / v));
        near(result.sorbed_amount_mol_m2[0][0], expected_s, 1e-15);
        near(
            result.aqueous_amount_mol_m2[0][0],
            a0 + s0 - expected_s,
            1e-15,
        );
        near(
            result.exchange_into_sorbed_mol_m2[0][0],
            expected_s - s0,
            1e-15,
        );
        near(result.report.total.amount_change_mol_m2[0], 0.0, 1e-15);
        assert_eq!(result.report.total.linear_solves, 1);
    }
}

#[test]
fn g1_equilibrium_diffusion_uses_retarded_storage_and_original_wet_resistance() {
    let v = [0.2, 0.4];
    let widths = [0.8, 1.2];
    let b = [0.8, 0.1];
    let c0 = [2.0, 0.5];
    let d = [[0.04], [0.09]];
    let tau = [2.0, 3.0];
    let a0 = [[v[0] * c0[0]], [v[1] * c0[1]]];
    let s0 = [[b[0] * c0[0]], [b[1] * c0[1]]];
    let g = 1.0
        / (0.5 * widths[0] / (v[0] / widths[0] * d[0][0] / tau[0])
            + 0.5 * widths[1] / (v[1] / widths[1] * d[1][0] / tau[1]));
    let cap = [v[0] + b[0], v[1] + b[1]];
    let mean = (cap[0] * c0[0] + cap[1] * c0[1]) / (cap[0] + cap[1]);
    for dt in [0.01, 1.0, 1e4] {
        let result = solve(
            dt,
            &v,
            &v,
            &widths,
            &[0.0; 3],
            &a0,
            &s0,
            &b.map(|capacity_m| [LinearSoluteSorption::Equilibrium { capacity_m }]),
            &d,
            &tau,
            [DiluteSoluteBoundary::Outflow; 2],
            controls(),
            || false,
        )
        .unwrap();
        let gap = (c0[0] - c0[1]) / (1.0 + dt * g * (1.0 / cap[0] + 1.0 / cap[1]));
        let expected = [
            mean + cap[1] / (cap[0] + cap[1]) * gap,
            mean - cap[0] / (cap[0] + cap[1]) * gap,
        ];
        for i in 0..2 {
            near(result.concentration_mol_m3[i][0], expected[i], 3e-13);
            near(result.sorbed_amount_mol_m2[i][0], b[i] * expected[i], 3e-13);
        }
        near(
            result.report.total.diffusive_face_amount_mol_m2[1][0],
            dt * g * gap,
            3e-13,
        );
        near(result.report.total.amount_before_mol_m2[0], 2.25, 1e-15);
        near(result.report.total.amount_change_mol_m2[0], 0.0, 3e-13);
    }
}

#[test]
fn g3_stiff_kinetic_sorption_converges_to_equilibrium_without_changing_face_laws() {
    let run = |kinetic| {
        let laws = [0.8, 0.1].map(|capacity_m| {
            [if kinetic {
                LinearSoluteSorption::Kinetic {
                    capacity_m,
                    rate_s_inverse: 1e14,
                }
            } else {
                LinearSoluteSorption::Equilibrium { capacity_m }
            }]
        });
        solve(
            1.0,
            &[0.2, 0.4],
            &[0.2, 0.4],
            &[0.8, 1.2],
            &[0.0; 3],
            &[[0.4], [0.2]],
            &[[1.6], [0.05]],
            &laws,
            &[[0.04], [0.09]],
            &[2.0, 3.0],
            [DiluteSoluteBoundary::Outflow; 2],
            controls(),
            || false,
        )
        .unwrap()
    };
    let a = run(false);
    let b = run(true);
    for i in 0..2 {
        near(
            a.concentration_mol_m3[i][0],
            b.concentration_mol_m3[i][0],
            1e-14,
        );
        near(
            a.sorbed_amount_mol_m2[i][0],
            b.sorbed_amount_mol_m2[i][0],
            1e-14,
        );
    }
    near(
        a.report.total.diffusive_face_amount_mol_m2[1][0],
        b.report.total.diffusive_face_amount_mol_m2[1][0],
        1e-16,
    );
}

#[test]
fn g3_zero_capacity_recovers_original_transport_and_zero_rate_preserves_isolated_pool() {
    let mobile = [[0.3, 0.0], [0.0, 0.6]];
    let d = [[0.1, 0.2]; 2];
    let bounds = [
        DiluteSoluteBoundary::Reservoir {
            concentration_mol_m3: [0.5, 0.2],
            film_resistance_s_m: [0.0, 3.0],
        },
        DiluteSoluteBoundary::Outflow,
    ];
    let plain = transport(
        [0.0, 1.0],
        1.0,
        &[0.3; 2],
        &[0.3; 2],
        &[1.0; 2],
        &[0.0; 3],
        &mobile,
        &d,
        &[1.0; 2],
        bounds,
        controls(),
        || false,
    )
    .unwrap();
    for isolated in [false, true] {
        let old = if isolated {
            [[1.0, 0.5]; 2]
        } else {
            [[0.0; 2]; 2]
        };
        let law = LinearSoluteSorption::Kinetic {
            capacity_m: if isolated { 2.0 } else { 0.0 },
            rate_s_inverse: if isolated { 0.0 } else { 3.0 },
        };
        let sorbed = solve(
            1.0,
            &[0.3; 2],
            &[0.3; 2],
            &[1.0; 2],
            &[0.0; 3],
            &mobile,
            &old,
            &[[law; 2]; 2],
            &d,
            &[1.0; 2],
            bounds,
            controls(),
            || false,
        )
        .unwrap();
        assert_eq!(sorbed.aqueous_amount_mol_m2, plain.amount_mol_m2);
        assert_eq!(sorbed.concentration_mol_m3, plain.concentration_mol_m3);
        assert_eq!(sorbed.sorbed_amount_mol_m2, old);
        assert_eq!(sorbed.exchange_into_sorbed_mol_m2, vec![[0.0; 2]; 2]);
        assert_eq!(
            sorbed.report.total.advective_face_amount_mol_m2,
            plain.report.advective_face_amount_mol_m2
        );
        assert_eq!(
            sorbed.report.total.diffusive_face_amount_mol_m2,
            plain.report.diffusive_face_amount_mol_m2
        );
    }
}

#[test]
fn g3_uniform_equilibrium_survives_changing_water_and_boundary_advection() {
    let before = [0.25, 0.6, 0.3];
    let flux = [0.02, 0.01, 0.015, 0.015];
    let after = core::array::from_fn::<_, 3, _>(|i| before[i] + flux[i] - flux[i + 1]);
    let b = [0.1, 0.7, 1.2];
    let c = 3.0;
    let result = solve(
        1.0,
        &before,
        &after,
        &[1.0; 3],
        &flux,
        &before.map(|v| [v * c]),
        &b.map(|v| [v * c]),
        &b.map(|capacity_m| [LinearSoluteSorption::Equilibrium { capacity_m }]),
        &[[0.02]; 3],
        &[1.0, 2.0, 3.0],
        [
            DiluteSoluteBoundary::InflowConcentration([c]),
            DiluteSoluteBoundary::Outflow,
        ],
        controls(),
        || false,
    )
    .unwrap();
    for i in 0..3 {
        near(result.concentration_mol_m3[i][0], c, 2e-15);
        near(result.sorbed_amount_mol_m2[i][0], b[i] * c, 2e-15);
        near(result.exchange_into_sorbed_mol_m2[i][0], 0.0, 2e-15);
    }
    near(
        result.report.total.amount_change_mol_m2[0],
        c * (flux[0] - flux[3]),
        2e-15,
    );
}

#[test]
fn g1_kinetic_relaxation_has_first_order_time_refinement_and_keeps_trace_release() {
    let (v, b, k, dt_total) = (0.25, 1.0, 0.7, 1.0);
    let exact = 0.24 + (0.1 - 0.24) * (-3.5_f64).exp();
    let errors = [16, 32, 64].map(|steps| {
        let (mut mobile, mut bound) = ([[0.2]], [[0.1]]);
        for _ in 0..steps {
            let r = solve(
                dt_total / steps as f64,
                &[v],
                &[v],
                &[1.0],
                &[0.0; 2],
                &mobile,
                &bound,
                &[[LinearSoluteSorption::Kinetic {
                    capacity_m: b,
                    rate_s_inverse: k,
                }]],
                &[[0.0]],
                &[1.0],
                [DiluteSoluteBoundary::Outflow; 2],
                controls(),
                || false,
            )
            .unwrap();
            mobile[0] = r.aqueous_amount_mol_m2[0];
            bound[0] = r.sorbed_amount_mol_m2[0];
        }
        (bound[0][0] - exact).abs()
    });
    assert!(
        errors[0] / errors[1] > 1.8 && errors[1] / errors[2] > 1.8,
        "{errors:?}"
    );
    let mut ctl = controls();
    ctl.absolute_equation_tolerance_mol_m2 = [1e-214];
    ctl.amount_balance_tolerance_mol_m2 = [1e-214];
    let result = solve(
        2.0,
        &[0.25],
        &[0.25],
        &[1.0],
        &[0.0; 2],
        &[[0.0]],
        &[[1e-200]],
        &[[LinearSoluteSorption::Kinetic {
            capacity_m: 0.0,
            rate_s_inverse: 0.3,
        }]],
        &[[0.0]],
        &[1.0],
        [DiluteSoluteBoundary::Outflow; 2],
        ctl,
        || false,
    )
    .unwrap();
    near(result.aqueous_amount_mol_m2[0][0] / 1e-200, 0.375, 2e-16);
    near(result.sorbed_amount_mol_m2[0][0] / 1e-200, 0.625, 2e-16);
}

#[test]
fn g0_equilibrium_mismatch_invalid_laws_and_additional_payload_budget_refuse() {
    let mobile = [[0.2]];
    let bound = [[0.1]];
    let call = |law, ctl| {
        solve(
            1.0,
            &[0.25],
            &[0.25],
            &[1.0],
            &[0.0; 2],
            &mobile,
            &bound,
            &[[law]],
            &[[0.0]],
            &[1.0],
            [DiluteSoluteBoundary::Outflow; 2],
            ctl,
            || false,
        )
    };
    assert!(matches!(
        call(
            LinearSoluteSorption::Equilibrium { capacity_m: 1.0 },
            controls()
        ),
        Err(DiluteSoluteError::Sorption { .. })
    ));
    for law in [
        LinearSoluteSorption::Equilibrium { capacity_m: -1.0 },
        LinearSoluteSorption::Kinetic {
            capacity_m: 1.0,
            rate_s_inverse: f64::NAN,
        },
    ] {
        assert!(matches!(
            call(law, controls()),
            Err(DiluteSoluteError::Invalid(_))
        ));
    }
    let mut ctl = controls::<1>();
    ctl.maximum_storage_bytes = transport::admit_extent(1, ctl).unwrap().1;
    assert!(matches!(
        call(
            LinearSoluteSorption::Kinetic {
                capacity_m: 1.0,
                rate_s_inverse: 1.0
            },
            ctl
        ),
        Err(DiluteSoluteError::Workspace { .. })
    ));
    assert_eq!(mobile, [[0.2]]);
    assert_eq!(bound, [[0.1]]);
}

#[test]
fn g4_cancellation_during_transport_exposes_no_partial_two_pool_state() {
    let n = 256;
    let mobile = vec![[0.2]; n];
    let bound = vec![[0.1]; n];
    let old_mobile = mobile.clone();
    let old_bound = bound.clone();
    let mut polls = 0;
    let result = solve(
        1.0,
        &vec![0.3; n],
        &vec![0.3; n],
        &vec![1.0; n],
        &vec![0.0; n + 1],
        &mobile,
        &bound,
        &vec![
            [LinearSoluteSorption::Kinetic {
                capacity_m: 1.0,
                rate_s_inverse: 0.5
            }];
            n
        ],
        &vec![[0.1]; n],
        &vec![1.0; n],
        [DiluteSoluteBoundary::Outflow; 2],
        controls(),
        || {
            polls += 1;
            polls == 30
        },
    );
    assert_eq!(result, Err(DiluteSoluteError::Cancelled));
    assert_eq!(mobile, old_mobile);
    assert_eq!(bound, old_bound);
}
