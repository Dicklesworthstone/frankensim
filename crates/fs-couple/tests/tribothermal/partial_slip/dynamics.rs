//! G1/G3 references use authored caloric/stiffness/friction data.
//! They check equations independently, not measured contact calibration.
use super::*;
use fs_couple::tribothermal::partial_slip::dynamics::*;
use fs_ga::Vec3;
use fs_mbd_material::{RegionBinding, ResolvedBody, VolumeMoments, resolve_body};
use std::cell::Cell;

fn motion_controls() -> InertialContactSolve {
    InertialContactSolve {
        contact: controls(),
        maximum_motion_trials: 80,
        maximum_newton_iterations: 24,
        difference_step_n_s: 1e-5,
        impulse_tolerance_n_s: 1e-8,
        total_energy_tolerance_j: 3e-7,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-10,
            relative_tolerance: 1e-12,
            forcing_minimum: 1e-4,
            linear_restart: 2,
            ..NewtonKrylovConfig::default()
        },
    }
}

fn body(card: &MaterialCard, side: f64) -> ResolvedBody {
    resolve_body(
        &[RegionBinding {
            key: "whole body",
            geometry: VolumeMoments::cuboid(Vec3::new(side, 1.0, 1.0)).unwrap(),
            material: card,
            point: &common::point(300.0),
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        1,
        || false,
    )
    .unwrap()
}

#[test]
fn planar_momentum_and_actual_midpoint_work_match_independent_elastic_solution() {
    let curves = curves();
    let thermal = initial(&curves, [300.0; 2], 0.0);
    let strain = [0.002, -0.003, 0.0];
    let velocities = [[1.0, 0.2], [-0.4, 0.3]];
    let forces = [[2.0, -1.0], [-1.0, 3.0]];
    let pair = InertialContactPair::new(
        &thermal,
        contact(&thermal, &curves, model(1000.0, zero_slopes()), strain),
        frame(),
        velocities,
    )
    .unwrap();
    let dt = 0.01;
    let step = pair
        .step_candidate(
            &thermal,
            dt,
            forces,
            &work("elastic inertial"),
            motion_controls(),
            [0.0; 2],
            &[0.0; 2],
            other,
            || false,
        )
        .unwrap();
    let r = step.report();
    let mut numerical = 0.0;
    let mut stored = 0.0;
    let mut kinetic = 0.0;
    let mut exterior = 0.0;
    for axis in 0..2 {
        // k=1000; mA=mB=1. e'=e+dt*(w+dt*(FA-FB+2F)/2), F=-k e'.
        let e = (strain[axis]
            + dt * (velocities[0][axis] - velocities[1][axis])
            + 0.5 * dt * dt * (forces[0][axis] - forces[1][axis]))
            / (1.0 + 1000.0 * dt * dt);
        let force = -1000.0 * e;
        near(r.contact.contact.tangent_force_n[axis], force);
        numerical += 500.0 * (e - strain[axis]).powi(2);
        stored += 500.0 * (e * e - strain[axis] * strain[axis]);
        for b in 0..2 {
            let v =
                velocities[b][axis] + dt * (forces[b][axis] + if b == 0 { force } else { -force });
            let dx = 0.5 * dt * (velocities[b][axis] + v);
            near(step.pair().velocities_mps()[b][axis], v);
            near(step.pair().translations_m()[b][axis], dx);
            near(r.body_displacements_m[b][axis], dx);
            kinetic += 0.5 * (v * v - velocities[b][axis].powi(2));
            exterior += forces[b][axis] * dx;
            near(r.body_momentum_residual_n_s[b][axis], 0.0);
        }
    }
    near(r.contact.contact.dissipation.heat_j, 0.0);
    near(
        r.contact.contact.generalized_work.numerical_dissipation_j,
        numerical,
    );
    near(r.contact.contact_internal_energy_change_j, stored);
    near(r.kinetic_energy_change_j, kinetic);
    near(r.external_work_j, exterior);
    near(kinetic + stored + numerical, exterior);
    near(r.total_energy_residual_j, 0.0);
    assert_eq!(step.pair().contact().temperature_k(), [300.0; 2]);
}

#[test]
fn source_mass_mu_temperature_and_microslip_heat_drive_the_same_motion() {
    let (a, ca) = common::material(1.0, 0.05, 0.0);
    let (b, cb) = common::material(1.0, 0.1, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let bodies = [body(&cards[0], 1.0), body(&cards[1], 1.0)];
    let point = context();
    let card = source(&cards, [(200.0, 0.1), (600.0, 0.5)], 600.0);
    let thermal = initial(&curves, [300.0; 2], 0.0);
    let original = loaded(
        &thermal,
        &curves,
        data(&cards, &curves, &card, &point),
        0.2,
        [[0.0; 2]; 2],
        1.0,
    );
    assert!(
        InertialContactPair::from_material_bodies(
            &thermal,
            original.clone(),
            [body(&cards[0], 2.0), bodies[1].clone()],
            frame(),
            [[1.0, 0.0], [0.0; 2]],
        )
        .is_err()
    );
    let pair = InertialContactPair::from_material_bodies(
        &thermal,
        original,
        bodies.clone(),
        frame(),
        [[1.0, 0.0], [0.0; 2]],
    )
    .unwrap();
    let step = pair
        .step_candidate(
            &thermal,
            0.1,
            [[0.0; 2]; 2],
            &work("source inertia"),
            motion_controls(),
            [-0.2, 0.0],
            &[0.0; 2],
            other,
            || false,
        )
        .unwrap();
    // Actual du=.1*(1-mu). Independent two-Jenkins plastic work:
    // D=10*mu*du - (1.5^2/750 + 8.5^2/250)*mu^2
    //   =mu-1.292*mu^2. Ca=.05, fractionA=.25 => mu=.2+.005D.
    let heat = root(0.0, 0.3, |d| {
        let mu = 0.2 + 0.005 * d;
        d - mu + 1.292 * mu * mu
    });
    let mu = 0.2 + 0.005 * heat;
    let r = step.report();
    near(step.pair().velocities_mps()[0][0], 1.0 - mu);
    near(step.pair().velocities_mps()[1][0], mu);
    near(r.relative_midpoint_velocity_mps[0], 1.0 - mu);
    near(r.contact.contact.dissipation.heat_j, heat);
    near(
        r.contact.contact.generalized_work.numerical_dissipation_j,
        0.146 * mu * mu,
    );
    near(r.contact.contact_internal_energy_change_j, 0.146 * mu * mu);
    near(r.contact.thermal.stored_energy_change_j, heat);
    near(step.pair().contact().temperature_k()[0], 300.0 + 5.0 * heat);
    near(step.pair().contact().temperature_k()[1], 300.0 + 7.5 * heat);
    let sample = r.contact.material.as_ref().unwrap();
    near(sample.resolved.state_coordinate("speed").unwrap(), 1.0 - mu);
    near(
        sample
            .resolved
            .property(KINETIC_FRICTION_PROPERTY)
            .unwrap()
            .value_si(),
        mu,
    );
    near(r.kinetic_energy_change_j, -mu + mu * mu);
    near(r.total_energy_residual_j, 0.0);
    assert_eq!(step.pair().material_bodies().unwrap(), &bodies);
    assert!(r.contact_trials > r.contact.contact_trials);
    assert_eq!(
        step.pair().contact().checkpoint().state().accepted_steps(),
        8
    );
}

#[test]
fn reciprocal_elastic_heat_and_inertia_match_independent_finite_capacity_root() {
    let curves = [
        common::material(1.0, 0.05, 0.0).1,
        common::material(1.0, 0.1, 0.0).1,
    ];
    let thermal = initial(&curves, [300.0; 2], 0.0);
    let pair = InertialContactPair::new(
        &thermal,
        contact(&thermal, &curves, model(1000.0, slopes()), [0.01, 0.0, 0.0]),
        frame(),
        [[1.0, 0.0], [0.0; 2]],
    )
    .unwrap();
    let dt = 0.01;
    let step = pair
        .step_candidate(
            &thermal,
            dt,
            [[0.0; 2]; 2],
            &work("heated elastic inertia"),
            motion_controls(),
            [0.0; 2],
            &[0.0; 2],
            |_| {
                Ok::<_, Infallible>(StepCoefficients {
                    sources_w: BTreeMap::from([(DESTINATIONS[0], 100.0)]),
                    conductances_w_k: BTreeMap::new(),
                    payload: (),
                })
            },
            || false,
        )
        .unwrap();
    let temperatures = |e: f64| {
        let difference = e * e - 0.01 * 0.01;
        [
            (0.05 * 300.0 + 1.0) / (0.05 + 0.5 * difference),
            0.1 * 300.0 / (0.1 + 0.25 * difference),
        ]
    };
    let strain = root(0.01, 0.02, |e| {
        let [ta, tb] = temperatures(e);
        let k = 1000.0 - (ta - 300.0) - 0.5 * (tb - 300.0);
        e * (1.0 + k * dt * dt) - 0.02
    });
    let [ta, tb] = temperatures(strain);
    let k = 1000.0 - (ta - 300.0) - 0.5 * (tb - 300.0);
    let force = -k * strain;
    let r = step.report();
    near(step.pair().contact().temperature_k()[0], ta);
    near(step.pair().contact().temperature_k()[1], tb);
    near(step.pair().velocities_mps()[0][0], 1.0 + dt * force);
    near(step.pair().velocities_mps()[1][0], -dt * force);
    near(
        r.contact.contact.next_state.region_elastic_displacement_m()[0][0],
        strain,
    );
    near(r.contact.reversible_heat_j[0], 0.05 * (ta - 300.0) - 1.0);
    near(r.contact.reversible_heat_j[1], 0.1 * (tb - 300.0));
    near(r.contact.contact.dissipation.heat_j, 0.0);
    near(r.contact.external_heat_j, 1.0);
    near(r.total_energy_residual_j, 0.0);
}

#[test]
fn retained_elastic_reversal_converts_storage_to_motion_without_friction_heat_and_refines() {
    let curves = curves();
    let horizon = 0.12;
    let omega = (2000.0_f64).sqrt();
    let exact = [
        0.01 * (omega * horizon).cos(),
        -0.01 * omega * (omega * horizon).sin(),
    ];
    let mut errors = Vec::new();
    for count in [24, 48, 96] {
        let mut thermal = initial(&curves, [300.0; 2], 0.0);
        let mut pair = InertialContactPair::new(
            &thermal,
            contact(
                &thermal,
                &curves,
                model(1000.0, zero_slopes()),
                [0.01, 0.0, 0.0],
            ),
            frame(),
            [[0.0; 2]; 2],
        )
        .unwrap();
        let dt = horizon / f64::from(count);
        let mut numerical = 0.0;
        let mut negative_velocity_seen = false;
        let mut reversal_seen = false;
        let mut e = 0.01;
        let mut w = 0.0;
        for index in 1..=count {
            let time = f64::from(index) * dt;
            let actual_dt = time - pair.contact().time_s();
            let next_e = (e + actual_dt * w) / (1.0 + 1000.0 * actual_dt * actual_dt);
            let next_w = w - 2000.0 * actual_dt * next_e;
            let step = pair
                .step_candidate(
                    &thermal,
                    time,
                    [[0.0; 2]; 2],
                    &work("elastic oscillator"),
                    motion_controls(),
                    [0.0; 2],
                    &[0.0; 2],
                    other,
                    || false,
                )
                .unwrap();
            near(
                step.pair()
                    .contact()
                    .checkpoint()
                    .state()
                    .region_elastic_displacement_m()[0][0],
                next_e,
            );
            let velocity = step.pair().velocities_mps()[0][0] - step.pair().velocities_mps()[1][0];
            near(velocity, next_w);
            negative_velocity_seen |= velocity < 0.0;
            reversal_seen |= negative_velocity_seen && velocity > 0.0;
            near(step.report().contact.contact.dissipation.heat_j, 0.0);
            numerical += step
                .report()
                .contact
                .contact
                .generalized_work
                .numerical_dissipation_j;
            let (next_thermal, next_pair, _) = step.into_parts();
            thermal = next_thermal;
            pair = next_pair;
            e = next_e;
            w = next_w;
        }
        assert!(reversal_seen);
        near(
            pair.kinetic_energy_j().unwrap()
                + pair.contact().thermodynamics().internal_energy_j
                + numerical,
            0.05,
        );
        assert_eq!(pair.contact().temperature_k(), [300.0; 2]);
        errors.push(((e - exact[0]).powi(2) + ((w - exact[1]) / omega).powi(2)).sqrt());
    }
    assert!(errors[0] / errors[1] > 1.7 && errors[0] / errors[1] < 2.3);
    assert!(errors[1] / errors[2] > 1.7 && errors[1] / errors[2] < 2.3);
    assert!(errors[2] < 9e-4);
}

#[test]
fn late_failure_cancellation_budget_and_incompatible_history_leave_pair_unchanged() {
    let curves = curves();
    let thermal = initial(&curves, [300.0; 2], 0.0);
    let old = contact(
        &thermal,
        &curves,
        model(1000.0, zero_slopes()),
        [0.01, 0.0, 0.0],
    );
    assert!(
        InertialContactPair::new(
            &thermal,
            old.clone(),
            frame().rotated(0.2).unwrap(),
            [[0.0; 2]; 2],
        )
        .is_err()
    );
    let pair = InertialContactPair::new(&thermal, old, frame(), [[0.0; 2]; 2]).unwrap();
    let cells = thermal.cells().clone();
    let checkpoint = pair.contact().checkpoint().clone();
    let calls = Cell::new(0);
    let failure = pair.step_candidate(
        &thermal,
        0.01,
        [[0.0; 2]; 2],
        &work("atomic"),
        motion_controls(),
        [0.0; 2],
        &[0.0; 2],
        |_| {
            calls.set(calls.get() + 1);
            if calls.get() >= 3 {
                Err("late source failure")
            } else {
                Ok(other(&thermal).unwrap())
            }
        },
        || false,
    );
    assert!(matches!(
        failure,
        Err(CaloricContactError::Source("late source failure"))
    ));
    let mut scarce = motion_controls();
    scarce.maximum_motion_trials = 1;
    assert!(matches!(
        pair.step_candidate(
            &thermal,
            0.01,
            [[0.0; 2]; 2],
            &work("atomic"),
            scarce,
            [0.0; 2],
            &[0.0; 2],
            other,
            || false,
        ),
        Err(CaloricContactError::TrialBudget)
    ));
    let polls = Cell::new(0);
    assert!(matches!(
        pair.step_candidate(
            &thermal,
            0.01,
            [[0.0; 2]; 2],
            &work("atomic"),
            motion_controls(),
            [0.0; 2],
            &[0.0; 2],
            other,
            || {
                polls.set(polls.get() + 1);
                polls.get() > 20
            },
        ),
        Err(CaloricContactError::Cancelled)
    ));
    assert_eq!(thermal.cells(), &cells);
    assert_eq!(pair.contact().checkpoint(), &checkpoint);
    assert_eq!(pair.velocities_mps(), [[0.0; 2]; 2]);
    assert_eq!(pair.translations_m(), [[0.0; 2]; 2]);
    let run = || {
        pair.step_candidate(
            &thermal,
            0.01,
            [[0.0; 2]; 2],
            &work("atomic"),
            motion_controls(),
            [0.0; 2],
            &[0.0; 2],
            other,
            || false,
        )
        .unwrap()
    };
    let a = run();
    let b = run();
    assert_eq!(a.thermal().cells(), b.thermal().cells());
    assert_eq!(a.pair().velocities_mps(), b.pair().velocities_mps());
    assert_eq!(a.pair().translations_m(), b.pair().translations_m());
    assert_eq!(
        a.pair().contact().checkpoint(),
        b.pair().contact().checkpoint()
    );
    assert!(
        pair.step_candidate(
            a.thermal(),
            0.02,
            [[0.0; 2]; 2],
            &work("wrong history"),
            motion_controls(),
            [0.0; 2],
            &[0.0; 2],
            other,
            || false,
        )
        .is_err()
    );
}
