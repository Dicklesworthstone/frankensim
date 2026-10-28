//! Independent normal-motion roots and complete SI energy accounts.
use super::*;
use fs_couple::tribothermal::normal::dynamics::*;

fn initial_bodies<'a>(
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    masses: [f64; 2],
) -> DepositionSimulation<'a> {
    let events = CELLS
        .into_iter()
        .enumerate()
        .map(|(index, id)| {
            let curve = curves[index].curve();
            let cp = (curve.knots().last().unwrap().specific_enthalpy_j_kg
                - curve.knots()[0].specific_enthalpy_j_kg)
                / 1800.0;
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id,
                    segment_id: index as u64,
                    body_id: index as u64,
                    material_id: MaterialId(index as u64),
                    mass_kg: masses[index],
                    specific_h_j_kg: 100.0 * cp,
                    geometry: CellGeometry {
                        start_m: [index as f64, 0.0, 0.0],
                        end_m: [index as f64 + 1.0, 0.0, 0.0],
                        volume_m3: masses[index],
                        surface_area_m2: 6.0,
                    },
                    boundary: common::boundary(),
                }),
            }
        })
        .collect();
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), curves[0].curve() as &dyn ThermalLaw),
            (MaterialId(1), curves[1].curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 4,
            max_links: 4,
            max_events: 8,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-7,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn motion_controls() -> NormalMotionSolve {
    NormalMotionSolve {
        contact: controls(),
        maximum_motion_trials: 200,
        maximum_newton_iterations: 24,
        difference_step_m: 1e-8,
        displacement_tolerance_m: 2e-11,
        momentum_tolerance_n_s: 1e-10,
        total_energy_tolerance_j: 2e-7,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-13,
            relative_tolerance: 1e-12,
            forcing_minimum: 1e-4,
            linear_restart: 1,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn root(mut low: f64, mut high: f64, f: impl Fn(f64) -> f64) -> f64 {
    assert!(f(low) <= 0.0 && f(high) >= 0.0);
    for _ in 0..90 {
        let x = 0.5 * (low + high);
        if f(x) > 0.0 {
            high = x;
        } else {
            low = x;
        }
    }
    0.5 * (low + high)
}
fn mean(e: f64, a: f64, b: f64) -> f64 {
    if a == b {
        (4.0 / 3.0) * e * a.powf(1.5)
    } else {
        e * (basis(b) - basis(a)) / (b - a)
    }
}
fn near_tight(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() < tol, "{a:e} != {b:e}");
}

#[test]
fn actual_two_body_momentum_matches_independent_hertz_root_and_drive_work() {
    let curves = [
        common::material(1.0, 10.0, 0.0).1,
        common::material(1.0, 20.0, 0.0).1,
    ];
    let masses = [2.0, 3.0];
    let thermal = initial_bodies(&curves, masses);
    let pair = InertialNormalPair::new(
        &thermal,
        bind(&thermal, &curves, model(0.002, 0.0, [0.0; 2])),
        [0.2, -0.1],
    )
    .unwrap();
    let dt = 0.001;
    let drive = [1.0, -2.0];
    let free = 0.002 + dt * 0.3 + 0.5 * dt * dt * (drive[0] / masses[0] - drive[1] / masses[1]);
    let expected = root(0.002, free, |d| {
        d - free + 0.5 * dt * dt * (1.0 / 2.0 + 1.0 / 3.0) * mean(1e6, 0.002, d)
    });
    let force = mean(1e6, 0.002, expected);
    let velocities = [
        0.2 + dt * (drive[0] - force) / 2.0,
        -0.1 + dt * (drive[1] + force) / 3.0,
    ];
    let displacements = [
        0.5 * dt * (0.2 + velocities[0]),
        0.5 * dt * (-0.1 + velocities[1]),
    ];
    let step = pair
        .step_candidate(
            &thermal,
            dt,
            drive,
            motion_controls(),
            free,
            &guess(&thermal),
            other,
            || false,
        )
        .unwrap();
    near_tight(step.pair().contact().endpoint().approach_m, expected, 2e-12);
    for i in 0..2 {
        near_tight(step.pair().velocities_mps()[i], velocities[i], 2e-10);
        near_tight(
            step.report().body_displacements_m[i],
            displacements[i],
            2e-13,
        );
    }
    near_tight(
        step.report().external_work_j,
        drive[0] * displacements[0] + drive[1] * displacements[1],
        1e-12,
    );
    near(
        step.report().contact.work.elastic_work_j,
        1e6 * (basis(expected) - basis(0.002)),
    );
    near(step.report().mechanical_energy_residual_j, 0.0);
    near(step.report().total_energy_residual_j, 0.0);
    assert_eq!(step.report().contact.work.physical_dissipation_j, 0.0);
    assert_eq!(step.pair().contact().temperature_k(), [300.0; 2]);
    assert_eq!(pair.velocities_mps(), [0.2, -0.1]);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn elastic_compression_and_rebound_preserve_energy_momentum_and_time_reversal() {
    let curves = [
        common::material(1.0, 10.0, 0.0).1,
        common::material(1.0, 20.0, 0.0).1,
    ];
    let mut thermal = initial_bodies(&curves, [2.0, 3.0]);
    let mut pair = InertialNormalPair::new(
        &thermal,
        bind(&thermal, &curves, model(0.002, 0.0, [0.0; 2])),
        [0.3, -0.1],
    )
    .unwrap();
    let initial_energy =
        pair.kinetic_energy_j().unwrap() + pair.contact().thermodynamics().internal_energy_j;
    let dt = 0.0002;
    for i in 1..=30 {
        let step = pair
            .step_candidate(
                &thermal,
                i as f64 * dt,
                [0.0; 2],
                motion_controls(),
                pair.contact().endpoint().approach_m,
                &guess(&thermal),
                other,
                || false,
            )
            .unwrap();
        near(step.report().total_energy_residual_j, 0.0);
        (thermal, pair, _) = step.into_parts();
        near_tight(
            2.0 * pair.velocities_mps()[0] + 3.0 * pair.velocities_mps()[1],
            0.3,
            2e-13,
        );
        near_tight(
            pair.kinetic_energy_j().unwrap() + pair.contact().thermodynamics().internal_energy_j,
            initial_energy,
            2e-8,
        );
    }
    assert!(pair.velocities_mps()[0] < pair.velocities_mps()[1]);
    assert!(pair.contact().endpoint().approach_m > 0.002);
    // The undamped isothermal discrete-gradient/midpoint rung is reversible.
    pair = InertialNormalPair::new(
        &thermal,
        pair.contact().clone(),
        pair.velocities_mps().map(|v| -v),
    )
    .unwrap();
    for i in 31..=60 {
        let step = pair
            .step_candidate(
                &thermal,
                i as f64 * dt,
                [0.0; 2],
                motion_controls(),
                pair.contact().endpoint().approach_m,
                &guess(&thermal),
                other,
                || false,
            )
            .unwrap();
        (thermal, pair, _) = step.into_parts();
    }
    near_tight(pair.contact().endpoint().approach_m, 0.002, 2e-10);
    near_tight(pair.velocities_mps()[0], -0.3, 2e-8);
    near_tight(pair.velocities_mps()[1], 0.1, 2e-8);
}

#[test]
fn inertial_hunt_crossley_loss_heats_the_two_bodies_exactly_once() {
    let curves = [
        common::material(1.0, 10.0, 0.0).1,
        common::material(1.0, 20.0, 0.0).1,
    ];
    let thermal = initial_bodies(&curves, [2.0, 3.0]);
    let dt = 0.001;
    let damping = 2.0;
    let pair = InertialNormalPair::new(
        &thermal,
        bind(&thermal, &curves, model(0.002, damping, [0.0; 2])),
        [0.2, -0.1],
    )
    .unwrap();
    let free = 0.002 + 0.3 * dt;
    let expected = root(0.002, free, |d| {
        d - free
            + 0.5
                * dt
                * dt
                * (1.0 / 2.0 + 1.0 / 3.0)
                * mean(1e6, 0.002, d)
                * (1.0 + damping * (d - 0.002) / dt)
    });
    let heat = damping * (expected - 0.002) / dt * 1e6 * (basis(expected) - basis(0.002));
    let step = pair
        .step_candidate(
            &thermal,
            dt,
            [0.0; 2],
            motion_controls(),
            free,
            &guess(&thermal),
            other,
            || false,
        )
        .unwrap();
    near_tight(step.pair().contact().endpoint().approach_m, expected, 2e-12);
    near_tight(
        step.report().contact.work.physical_dissipation_j,
        heat,
        2e-10,
    );
    near_tight(
        step.pair().contact().temperature_k()[0],
        300.0 + 0.4 * heat / 20.0,
        2e-9,
    );
    near_tight(
        step.pair().contact().temperature_k()[1],
        300.0 + 0.6 * heat / 60.0,
        2e-9,
    );
    near_tight(
        step.report().kinetic_energy_change_j + step.report().contact.work.internal_energy_change_j,
        -heat,
        2e-9,
    );
    near(step.report().total_energy_residual_j, 0.0);
}

#[test]
fn reversible_normal_cooling_changes_the_inertial_solution_at_the_same_endpoint() {
    let curves = [
        common::material(1.0, 0.03, 0.0).1,
        common::material(1.0, 0.04, 0.0).1,
    ];
    let thermal = initial_bodies(&curves, [2.0, 3.0]);
    let dt = 0.001;
    let pair = InertialNormalPair::new(
        &thermal,
        bind(&thermal, &curves, model(0.002, 0.0, [-1000.0, -500.0])),
        [0.2, -0.1],
    )
    .unwrap();
    let temperatures = |d| {
        let dg = basis(d) - basis(0.002);
        [
            300.0 * 0.06 / (0.06 + 1000.0 * dg),
            300.0 * 0.12 / (0.12 + 500.0 * dg),
        ]
    };
    let stiffness = |d| {
        let t = temperatures(d);
        1e6 - 1000.0 * (t[0] - 300.0) - 500.0 * (t[1] - 300.0)
    };
    let free = 0.002 + 0.3 * dt;
    let expected = root(0.002, free, |d| {
        d - free + 0.5 * dt * dt * (1.0 / 2.0 + 1.0 / 3.0) * mean(stiffness(d), 0.002, d)
    });
    let step = pair
        .step_candidate(
            &thermal,
            dt,
            [0.0; 2],
            motion_controls(),
            free,
            &guess(&thermal),
            other,
            || false,
        )
        .unwrap();
    near_tight(step.pair().contact().endpoint().approach_m, expected, 2e-12);
    for i in 0..2 {
        near_tight(
            step.pair().contact().temperature_k()[i],
            temperatures(expected)[i],
            2e-8,
        );
    }
    near_tight(
        step.report().contact.modulus.modulus_pa,
        stiffness(expected),
        5e-5,
    );
    assert_eq!(step.report().contact.work.physical_dissipation_j, 0.0);
    assert!(
        step.report()
            .contact
            .reversible_heat_j
            .iter()
            .all(|q| *q < 0.0)
    );
    near(step.report().total_energy_residual_j, 0.0);
}

#[test]
fn late_sources_separation_cancellation_and_budget_refusal_preserve_history() {
    let curves = [
        common::material(1.0, 0.03, 0.0).1,
        common::material(1.0, 0.04, 0.0).1,
    ];
    let thermal = initial_bodies(&curves, [2.0, 3.0]);
    let pair = InertialNormalPair::new(
        &thermal,
        bind(&thermal, &curves, model(0.002, 0.0, [-1000.0, -500.0])),
        [0.2, -0.1],
    )
    .unwrap();
    let run = || {
        pair.step_candidate(
            &thermal,
            0.001,
            [0.0; 2],
            motion_controls(),
            0.002,
            &guess(&thermal),
            other,
            || false,
        )
        .unwrap()
    };
    let first = run();
    assert!(matches!(
        pair.step_candidate(
            &thermal,
            0.001,
            [0.0; 2],
            motion_controls(),
            0.002,
            &guess(&thermal),
            other,
            || true
        ),
        Err(NormalThermalError::Cancelled)
    ));
    let mut limited = motion_controls();
    limited.maximum_motion_trials = 1;
    assert!(matches!(
        pair.step_candidate(
            &thermal,
            0.001,
            [0.0; 2],
            limited,
            0.002,
            &guess(&thermal),
            other,
            || false
        ),
        Err(NormalThermalError::TrialBudget)
    ));
    let failed = pair.step_candidate(
        &thermal,
        0.001,
        [0.0; 2],
        motion_controls(),
        0.002,
        &guess(&thermal),
        |state| {
            if state.temperature_k(CELLS[0]).unwrap() < 299.999 {
                return Err("original thermal source boundary");
            }
            Ok(StepCoefficients {
                sources_w: BTreeMap::new(),
                conductances_w_k: BTreeMap::new(),
                payload: (),
            })
        },
        || false,
    );
    assert!(matches!(
        failed,
        Err(NormalThermalError::Source(
            "original thermal source boundary"
        ))
    ));
    let separating =
        InertialNormalPair::new(&thermal, pair.contact().clone(), [-0.3, 0.3]).unwrap();
    assert!(matches!(
        separating.step_candidate(
            &thermal,
            0.1,
            [0.0; 2],
            motion_controls(),
            0.0,
            &guess(&thermal),
            other,
            || false
        ),
        Err(NormalThermalError::Invalid(_))
    ));
    let replay = run();
    assert_eq!(
        first.pair().velocities_mps(),
        replay.pair().velocities_mps()
    );
    assert_eq!(
        first.pair().contact().temperature_k(),
        replay.pair().contact().temperature_k()
    );
    assert_eq!(
        first.pair().contact().endpoint(),
        replay.pair().contact().endpoint()
    );
    assert_eq!(pair.velocities_mps(), [0.2, -0.1]);
    assert_eq!(pair.contact().endpoint().approach_m, 0.002);
    assert_eq!(thermal.time_s(), 0.0);
}
