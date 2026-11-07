//! Independent finite-bath dynamics and complete-window physical accounts.
use super::*;
use crate::solution as aqueous;
use fs_couple::moisture::{
    adaptive::LayeredMoistureAdvanceReport,
    radial::{
        MoistureRadial, MoistureRadialGeometry, RadialMoistureSymmetry,
        layered::{MoistureLayeredRadial, MoistureRadialLayer},
    },
    solution::{FiniteMoistureSolution, SolutionMoistureStep},
};

const VOLUME: f64 = 1e-4;
const MOLES: f64 = 0.5;
const SOLVENT_DENSITY: f64 = 0.018 / aqueous::WATER_VOLUME;

#[derive(Clone, Copy)]
struct Physics {
    capacity: [f64; 2],
    join: f64,
    outer: f64,
}

fn time_controls(maximum: f64) -> MoistureTimeControls {
    MoistureTimeControls {
        initial_step_s: maximum,
        minimum_step_s: 1e-8,
        maximum_step_s: maximum,
        maximum_evaluations: 5000,
        moisture_absolute_tolerance: 1.0,
        activity_absolute_tolerance: 1.0,
        reservoir_water_absolute_tolerance_kg: 1.0,
        relative_tolerance: 0.0,
        water_balance_tolerance_kg: 1e-10,
    }
}

fn curved(
    symmetry: RadialMoistureSymmetry,
    initial: [f64; 2],
    cx: &Cx<'_>,
) -> (MoistureLayeredRadial, Physics) {
    let rho = [300.0, 900.0];
    let slope = [0.2, 0.4];
    let diffusion = [0.001, 0.0005];
    let volume = |a: f64, b: f64| match symmetry {
        RadialMoistureSymmetry::Sphere => {
            4.0 * core::f64::consts::PI * (b.powi(3) - a.powi(3)) / 3.0
        }
        RadialMoistureSymmetry::Cylinder { length_m } => {
            core::f64::consts::PI * length_m * (b * b - a * a)
        }
    };
    let area = |r: f64| match symmetry {
        RadialMoistureSymmetry::Sphere => 4.0 * core::f64::consts::PI * r * r,
        RadialMoistureSymmetry::Cylinder { length_m } => 2.0 * core::f64::consts::PI * length_m * r,
    };
    let resistance = |a: f64, b: f64| match symmetry {
        RadialMoistureSymmetry::Sphere => (1.0 / a - 1.0 / b) / (4.0 * core::f64::consts::PI),
        RadialMoistureSymmetry::Cylinder { length_m } => {
            (b / a).ln() / (2.0 * core::f64::consts::PI * length_m)
        }
    };
    let half = [
        rho[0] * slope[0] * diffusion[0] / resistance(0.015, 0.03),
        rho[1] * slope[1] * diffusion[1] / resistance(0.03, 0.04),
    ];
    let outer_bulk = rho[1] * slope[1] * diffusion[1] / resistance(0.04, 0.05);
    let physics = Physics {
        capacity: [
            rho[0] * slope[0] * volume(0.0, 0.03),
            rho[1] * slope[1] * volume(0.03, 0.05),
        ],
        join: 1.0 / (1.0 / half[0] + 1.0 / (0.03 * area(0.03)) + 1.0 / half[1]),
        outer: 1.0 / (1.0 / outer_bulk + 1.0 / (0.04 * area(0.05))),
    };
    let body = MoistureLayeredRadial::new(
        (0..2)
            .map(|i| {
                let card = material(rho[i], vec![(0.0, 0.0), (1.0, slope[i])], [diffusion[i]; 2]);
                let sorption = isotherm(&card);
                let field = MoistureRadial::new(
                    card,
                    point(0.1),
                    "moisture_diffusivity",
                    MaterialPropertySelection::SingleClaimOnly,
                    MoistureRadialGeometry {
                        symmetry,
                        faces_m: if i == 0 {
                            vec![0.0, 0.03]
                        } else {
                            vec![0.03, 0.05]
                        },
                        maximum_moisture: 0.6,
                        maximum_cells: 2,
                    },
                    vec![slope[i] * initial[i]],
                    cx,
                )
                .unwrap();
                MoistureRadialLayer::new(field, sorption, cx).unwrap()
            })
            .collect(),
        2,
        cx,
    )
    .unwrap()
    .with_transfer(
        [
            MoistureActivityTransfer::IdealContact,
            MoistureActivityTransfer::ConductanceKgM2S(0.04),
        ],
        vec![MoistureActivityTransfer::ConductanceKgM2S(0.03)],
        cx,
    )
    .unwrap();
    (body, physics)
}

fn bath_volume(p: Physics, initial: [f64; 2], state: [f64; 2]) -> f64 {
    VOLUME
        - (0..2)
            .map(|i| p.capacity[i] * (state[i] - initial[i]))
            .sum::<f64>()
            / SOLVENT_DENSITY
}

// Independent conservation eliminates the finite bath as V(a). The resulting
// continuous two-cell ODE is integrated by RK4, with a doubled reference mesh
// checked below. No native BE solve, adaptive driver, or numerical layer kernel
// enters this reference. Exterior and interface conductances are derived above.
fn reference(p: Physics, initial: [f64; 2], duration: f64, count: u32) -> [f64; 2] {
    let rate = |a: [f64; 2]| {
        let volume = bath_volume(p, initial, a);
        assert!(volume > 0.0);
        let activity = (-aqueous::WATER_VOLUME * MOLES / volume).exp();
        let join = p.join * (a[1] - a[0]);
        [
            join / p.capacity[0],
            (-join + p.outer * (activity - a[1])) / p.capacity[1],
        ]
    };
    let h = duration / f64::from(count);
    let mut state = initial;
    for _ in 0..count {
        let k1 = rate(state);
        let k2 = rate(std::array::from_fn(|i| state[i] + 0.5 * h * k1[i]));
        let k3 = rate(std::array::from_fn(|i| state[i] + 0.5 * h * k2[i]));
        let k4 = rate(std::array::from_fn(|i| state[i] + h * k3[i]));
        for i in 0..2 {
            state[i] += h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
        }
    }
    state
}

fn activity_error(actual: &[f64], expected: [f64; 2]) -> f64 {
    actual
        .iter()
        .zip(expected)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max)
}

fn accounts(
    bath: &FiniteMoistureSolution,
    endpoint: &FiniteMoistureSolution,
    before: f64,
    after: f64,
    report: &LayeredMoistureAdvanceReport<SolutionMoistureStep<MoistureLayeredReport>>,
) {
    let r = &report.summary;
    near(
        endpoint.solvent_mass_kg() + after,
        bath.solvent_mass_kg() + before,
        1e-10,
    );
    near(r.water_change_kg, after - before, 1e-11);
    near(r.boundary_water_kg.iter().sum(), after - before, 1e-11);
    near(r.transport_balance_residual_kg, 0.0, 1e-10);
    near(r.water_balance_residual_kg, 0.0, 1e-10);
    assert_eq!(
        r.reservoir_water_kg,
        Some([bath.solvent_mass_kg(), endpoint.solvent_mass_kg()])
    );
    assert_eq!(report.layers.len(), 2);
    assert_eq!(
        report.layers[0].boundary_water_kg[1],
        -report.layers[1].boundary_water_kg[0]
    );
    for (i, layer) in report.layers.iter().enumerate() {
        assert_eq!(layer.id, i as u64);
        near(
            layer.water_change_kg,
            layer.boundary_water_kg.iter().sum(),
            1e-11,
        );
        assert!(layer.maximum_accepted_error_ratio <= 1.0);
    }
    assert!(r.maximum_accepted_error_ratio <= 1.0);
    assert_eq!(
        endpoint.solution().solute_moles(),
        bath.solution().solute_moles()
    );
}

#[test]
fn layered_solution_duration_refines_to_independent_finite_bath_dynamics() {
    with_cx(|cx| {
        let initial = [0.8, 0.1];
        let body = stack(initial, [0.2, 0.4], cx);
        let saved = body.activity().to_vec();
        let p = Physics {
            capacity: CAPACITY,
            join: 1.0 / (1.0 / HALF[0] + 1.0 / HALF[1]),
            outer: HALF[1],
        };
        let exact = reference(p, initial, 1.0, 8192);
        let checked = reference(p, initial, 1.0, 4096);
        for i in 0..2 {
            near(exact[i], checked[i], 2e-12);
        }
        let bath = aqueous::reservoir(VOLUME, vec![MOLES]);
        let original_bath = bath.solution().clone();
        let mut errors = Vec::new();
        for maximum in [0.125, 0.0625, 0.03125] {
            let (next_bath, next, report) = bath
                .advance_layered_adaptive_to(
                    &body,
                    1.0,
                    [false, true],
                    config(),
                    nonlinear(),
                    aqueous::controls(VOLUME),
                    time_controls(maximum),
                    cx,
                )
                .unwrap();
            errors.push(activity_error(next.activity(), exact));
            assert_eq!(report.summary.rejected_intervals, 0);
            assert_eq!(report.summary.interval_s, [0.0, 1.0]);
            assert_eq!(report.summary.accepted_intervals, (1.0 / maximum) as usize);
            assert_eq!(
                report.summary.evaluations,
                3 * report.summary.accepted_intervals
            );
            assert_eq!(
                report.summary.accepted_substeps,
                2 * report.summary.accepted_intervals
            );
            assert_eq!(next.time_s(), 1.0);
            accounts(
                &bath,
                &next_bath,
                body.water_mass_kg(),
                next.water_mass_kg(),
                &report,
            );
            assert!(
                next_bath.solution().water_activity().unwrap()
                    < bath.solution().water_activity().unwrap()
            );
            let expected_volume = bath_volume(p, initial, [next.activity()[0], next.activity()[1]]);
            near(
                next_bath.solution().solvent_volume_m3(),
                expected_volume,
                1e-13,
            );
        }
        for pair in errors.windows(2) {
            assert!(
                pair[0] / pair[1] > 1.8 && pair[0] / pair[1] < 2.2,
                "{errors:?}"
            );
        }
        assert_eq!(body.activity(), saved);
        assert_eq!(body.time_s(), 0.0);
        assert_eq!(bath.solution(), &original_bath);
    });
}

#[test]
fn coated_solution_duration_resolves_curved_contact_and_bath_feedback() {
    with_cx(|cx| {
        for symmetry in [
            RadialMoistureSymmetry::Sphere,
            RadialMoistureSymmetry::Cylinder { length_m: 0.1 },
        ] {
            let initial = [0.8, 0.1];
            let (body, p) = curved(symmetry, initial, cx);
            let saved = body.activity().to_vec();
            let exact = reference(p, initial, 4.0, 8192);
            let checked = reference(p, initial, 4.0, 4096);
            for i in 0..2 {
                near(exact[i], checked[i], 2e-12);
            }
            let bath = aqueous::reservoir(VOLUME, vec![MOLES]);
            let mut errors = Vec::new();
            for maximum in [1.0, 0.5, 0.25] {
                let (next_bath, next, report) = bath
                    .advance_layered_radial_adaptive_to(
                        &body,
                        4.0,
                        [false, true],
                        config(),
                        nonlinear(),
                        aqueous::controls(VOLUME),
                        time_controls(maximum),
                        cx,
                    )
                    .unwrap();
                errors.push(activity_error(next.activity(), exact));
                assert_eq!(report.summary.rejected_intervals, 0);
                assert_eq!(report.summary.interval_s, [0.0, 4.0]);
                assert_eq!(next.time_s(), 4.0);
                accounts(
                    &bath,
                    &next_bath,
                    body.water_mass_kg(),
                    next.water_mass_kg(),
                    &report,
                );
                assert!(report.summary.boundary_water_kg[1] > 0.0);
                assert_eq!(report.summary.boundary_water_kg[0], 0.0);
                assert!(
                    next_bath.solution().water_activity().unwrap()
                        < bath.solution().water_activity().unwrap()
                );
                near(
                    next_bath.solution().solvent_volume_m3(),
                    bath_volume(p, initial, [next.activity()[0], next.activity()[1]]),
                    1e-13,
                );
                for i in 0..2 {
                    assert_eq!(
                        next.layers()[i].radial().dry_masses_kg(),
                        body.layers()[i].radial().dry_masses_kg()
                    );
                    assert_eq!(
                        next.layers()[i].sorption().identity(),
                        body.layers()[i].sorption().identity()
                    );
                }
            }
            for pair in errors.windows(2) {
                assert!(
                    pair[0] / pair[1] > 1.8 && pair[0] / pair[1] < 2.2,
                    "{symmetry:?}: {errors:?}"
                );
            }
            assert_eq!(body.activity(), saved);
            assert_eq!(body.time_s(), 0.0);
            assert_eq!(bath.solution().solvent_volume_m3(), VOLUME);
        }
    });
}

fn closed_reference(p: Physics, initial: [f64; 2], duration: f64) -> [f64; 2] {
    let mean =
        (p.capacity[0] * initial[0] + p.capacity[1] * initial[1]) / (p.capacity[0] + p.capacity[1]);
    let decay = (-p.join * (1.0 / p.capacity[0] + 1.0 / p.capacity[1]) * duration).exp();
    initial.map(|a| mean + (a - mean) * decay)
}

#[test]
fn local_water_error_refines_closed_layers_despite_constant_pure_solvent_activity() {
    with_cx(|cx| {
        let initial = [0.8, 0.1];
        let bath = aqueous::reservoir(VOLUME, vec![0.0]);
        let body = stack(initial, [0.2, 0.4], cx);
        let p = Physics {
            capacity: CAPACITY,
            join: 1.0 / (1.0 / HALF[0] + 1.0 / HALF[1]),
            outer: 0.0,
        };
        let exact = closed_reference(p, initial, 1.0);
        let loose = time_controls(1.0);
        let tight = MoistureTimeControls {
            reservoir_water_absolute_tolerance_kg: 1e-6,
            ..loose
        };
        let (_, a, report_a) = bath
            .advance_layered_adaptive_to(
                &body,
                1.0,
                [false; 2],
                config(),
                nonlinear(),
                aqueous::controls(VOLUME),
                loose,
                cx,
            )
            .unwrap();
        let (next_bath, b, report_b) = bath
            .advance_layered_adaptive_to(
                &body,
                1.0,
                [false; 2],
                config(),
                nonlinear(),
                aqueous::controls(VOLUME),
                tight,
                cx,
            )
            .unwrap();
        assert!(report_b.summary.evaluations > report_a.summary.evaluations);
        assert!(report_b.summary.rejected_intervals > 0);
        assert!(activity_error(b.activity(), exact) < 0.3 * activity_error(a.activity(), exact));
        assert_eq!(report_b.summary.activity, Some([1.0; 2]));
        assert_eq!(report_b.summary.boundary_water_kg, [0.0; 2]);
        assert_eq!(next_bath.solution().solvent_volume_m3(), VOLUME);
        accounts(
            &bath,
            &next_bath,
            body.water_mass_kg(),
            b.water_mass_kg(),
            &report_b,
        );

        for symmetry in [
            RadialMoistureSymmetry::Sphere,
            RadialMoistureSymmetry::Cylinder { length_m: 0.1 },
        ] {
            let (body, p) = curved(symmetry, initial, cx);
            let exact = closed_reference(p, initial, 4.0);
            let loose = time_controls(4.0);
            let tight = MoistureTimeControls {
                reservoir_water_absolute_tolerance_kg: 1e-6,
                ..loose
            };
            let (_, a, report_a) = bath
                .advance_layered_radial_adaptive_to(
                    &body,
                    4.0,
                    [false; 2],
                    config(),
                    nonlinear(),
                    aqueous::controls(VOLUME),
                    loose,
                    cx,
                )
                .unwrap();
            let (next_bath, b, report_b) = bath
                .advance_layered_radial_adaptive_to(
                    &body,
                    4.0,
                    [false; 2],
                    config(),
                    nonlinear(),
                    aqueous::controls(VOLUME),
                    tight,
                    cx,
                )
                .unwrap();
            assert!(report_b.summary.evaluations > report_a.summary.evaluations);
            assert!(report_b.summary.rejected_intervals > 0);
            assert!(
                activity_error(b.activity(), exact) < 0.3 * activity_error(a.activity(), exact)
            );
            assert_eq!(report_b.summary.activity, Some([1.0; 2]));
            assert_eq!(report_b.summary.boundary_water_kg, [0.0; 2]);
            assert_eq!(next_bath.solution().solvent_volume_m3(), VOLUME);
            accounts(
                &bath,
                &next_bath,
                body.water_mass_kg(),
                b.water_mass_kg(),
                &report_b,
            );
            for i in 0..2 {
                assert_eq!(
                    b.layers()[i].radial().dry_masses_kg(),
                    body.layers()[i].radial().dry_masses_kg()
                );
            }
        }
    });
}

#[test]
fn layered_solution_duration_discards_a_successful_prefix_when_work_is_exhausted() {
    with_cx(|cx| {
        let bath = aqueous::reservoir(VOLUME, vec![MOLES]);
        let original_bath = bath.solution().clone();
        let body = stack([0.8, 0.1], [0.2, 0.4], cx);
        let saved = body.activity().to_vec();
        let controls = MoistureTimeControls {
            maximum_evaluations: 3,
            ..time_controls(0.25)
        };
        assert!(matches!(
            bath.advance_layered_adaptive_to(
                &body,
                0.5,
                [false, true],
                config(),
                nonlinear(),
                aqueous::controls(VOLUME),
                controls,
                cx,
            ),
            Err(MoistureAdvanceError::WorkLimit {
                time_s: 0.25,
                accepted_intervals: 1,
                evaluations: 3
            })
        ));
        assert_eq!(body.activity(), saved);
        assert_eq!(body.time_s(), 0.0);
        assert_eq!(bath.solution(), &original_bath);
        for symmetry in [
            RadialMoistureSymmetry::Sphere,
            RadialMoistureSymmetry::Cylinder { length_m: 0.1 },
        ] {
            let (body, _) = curved(symmetry, [0.8, 0.1], cx);
            let saved = body.activity().to_vec();
            assert!(matches!(
                bath.advance_layered_radial_adaptive_to(
                    &body,
                    0.5,
                    [false, true],
                    config(),
                    nonlinear(),
                    aqueous::controls(VOLUME),
                    controls,
                    cx,
                ),
                Err(MoistureAdvanceError::WorkLimit {
                    time_s: 0.25,
                    accepted_intervals: 1,
                    evaluations: 3
                })
            ));
            assert_eq!(body.activity(), saved);
            assert_eq!(body.time_s(), 0.0);
            assert_eq!(bath.solution(), &original_bath);
        }
    });
}
