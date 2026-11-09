//! G1/G3/G4: curved native Hertz paths and the original source/thermal owners.
use super::*;
use fs_couple::friction_heat::hertz::transient::propose_rigid_hertz_transient_heat;
use fs_mbd::tangential_friction::hertz::transient::{
    HertzTransientControls, RigidHertzTransientError, RigidHertzTransientReport,
};

fn path_controls() -> HertzTransientControls {
    HertzTransientControls {
        maximum_step_s: 0.02,
        minimum_step_s: 1e-15,
        maximum_attempts: 100_000,
        maximum_evaluations: 400_001,
        scaled_rate_tolerance_m_s: 1e-11,
        scaled_impulse_tolerance_ns: 1e-10,
        scaled_travel_tolerance_m: 1e-10,
        heat_tolerance_j: 1e-12,
        stopping_time_tolerance_s: 1e-7,
    }
}
fn run_transient<'m>(
    pair: &Pair,
    owner: &SourceThermalCoulomb,
    thermal: &DepositionSimulation<'m>,
    route: &FrictionHeatRoute,
    duration: f64,
    path: HertzTransientControls,
    resolve: impl Fn(&DepositionSimulation<'m>, f64) -> Result<RigidHertzHeatSource<f64>, &'static str>,
    gate: &CancelGate,
) -> Result<
    HeatedRigidHertzCandidate<'m, f64, RigidHertzTransientReport>,
    RigidHertzHeatError<&'static str, RigidHertzTransientError>,
> {
    propose_rigid_hertz_transient_heat(
        &pair.bodies,
        &pair.properties,
        pair.contact(owner),
        thermal,
        route,
        thermal.time_s(),
        duration,
        controls(),
        path,
        &[0.0; 3],
        300.0,
        resolve,
        gate,
    )
}

#[test]
fn curved_slip_spin_heat_and_source_temperature_close_before_and_after_arrest() {
    let pair = Pair::new(0.1, 3.0);
    let owner = law(0);
    let thermal = thermal([300.0; 3], So3::identity());
    let routes = routes(&thermal);
    let initial_relative = 0.5 * (0.1 * 0.1 / mobility(false) + 3.0 * 3.0 / mobility(true));
    for duration in [0.2, 12.0] {
        let candidate = run_transient(
            &pair,
            &owner,
            &thermal,
            &routes[0],
            duration,
            path_controls(),
            source,
            &CancelGate::new(),
        )
        .unwrap();
        let mechanical = &candidate.mechanical;
        let path = &mechanical.transient;
        eprintln!(
            "final Hertz heat trial: duration={duration}, accepted={}, rejected={}, evaluations={}, thermal_trials={}",
            path.accepted_steps, path.rejected_steps, path.evaluations, candidate.trials,
        );
        assert!(!path.exact_ray && path.accepted_steps > 0);
        let remaining = 0.5
            * (mechanical.slip_velocity_m_s[1][0].powi(2) / mobility(false)
                + mechanical.spin_rad_s[1].powi(2) / mobility(true));
        let heat = initial_relative - remaining;
        close(mechanical.dissipated_energy_j, heat);
        close(candidate.heat.dissipated_j, heat);
        close(energy(&candidate.thermal) - energy(&thermal), heat);
        close(
            kinetic(&candidate.bodies_after) - kinetic(&pair.bodies),
            -heat,
        );
        close(candidate.combined_energy_residual_j, 0.0);
        close(mechanical.audit.total_energy_residual_j, 0.0);
        if duration > 1.0 {
            let stop = path.stopping_time_bracket_s.unwrap();
            assert!(stop[0] > 0.0 && stop[1] < duration);
            assert!(stop[1] - stop[0] <= path_controls().stopping_time_tolerance_s);
            close(remaining, 0.0);
            close(heat, initial_relative);
        } else {
            assert!(path.stopping_time_bracket_s.is_none());
            assert!(remaining > 0.0 && heat > 0.0);
        }
        // Independent four-coordinate finite-cell/interface balance. The third
        // cell remains on its original ordinary link, outside the heated face.
        let d0 = duration * K[0] * AREA / A;
        let d1 = duration * K[1] * AREA / A;
        let d02 = duration * AREA / (A / K[0] + A / K[2]);
        let expected = solve_dense(
            vec![
                vec![CAP[0] + d0 + d02, 0.0, -d02, -d0],
                vec![0.0, CAP[1] + d1, 0.0, -d1],
                vec![-d02, 0.0, CAP[2] + d02, 0.0],
                vec![-d0, -d1, 0.0, d0 + d1],
            ],
            vec![0.0, 0.0, 0.0, heat],
        );
        for i in 0..3 {
            close(
                candidate.thermal.temperature_k(CellId(i as u64)).unwrap(),
                300.0 + expected[i],
            );
        }
        close(candidate.source, 300.0 + expected[3]);
        assert!(candidate.source > 300.0);
        close(
            candidate.friction_source.law().parameters().kinetic_mu,
            MU[0] - SLOPE[0] * expected[3],
        );
        assert_eq!(candidate.friction_source.owner_identity(), owner.identity());
        assert_eq!(
            candidate
                .friction_source
                .friction_source()
                .state_coordinate("T"),
            Some(candidate.source)
        );
        let entropy = candidate.heat.entropy.unwrap();
        assert!(entropy.friction_generation_j_k > 0.0);
        assert!(entropy.network_conduction_generation_j_k > 0.0);
    }
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(energy(&thermal), 0.0);
}

#[test]
fn curved_source_cancellation_and_quadrature_refusal_preserve_exact_retry() {
    let pair = Pair::new(0.1, 3.0);
    let saved = pair.bodies;
    let owner = law(0);
    let original = owner.clone();
    let thermal = thermal([300.0; 3], So3::identity());
    let saved_cells = thermal.cells().clone();
    let routes = routes(&thermal);
    let duration = 0.2;
    let reference = run_transient(
        &pair,
        &owner,
        &thermal,
        &routes[0],
        duration,
        path_controls(),
        source,
        &CancelGate::new(),
    )
    .unwrap();
    let calls = Cell::new(0);
    let refused = run_transient(
        &pair,
        &owner,
        &thermal,
        &routes[0],
        duration,
        path_controls(),
        |t, temperature| {
            calls.set(calls.get() + 1);
            if calls.get() == reference.trials {
                Err("late curved thermal source")
            } else {
                source(t, temperature)
            }
        },
        &CancelGate::new(),
    );
    assert!(matches!(
        refused,
        Err(RigidHertzHeatError::Source("late curved thermal source"))
    ));
    assert_eq!(calls.get(), reference.trials);
    let gate = CancelGate::new();
    let calls = Cell::new(0);
    assert!(
        run_transient(
            &pair,
            &owner,
            &thermal,
            &routes[0],
            duration,
            path_controls(),
            |t, temperature| {
                calls.set(calls.get() + 1);
                if calls.get() == reference.trials {
                    gate.request();
                }
                source(t, temperature)
            },
            &gate,
        )
        .is_err()
    );
    assert_eq!(calls.get(), reference.trials);
    let exhausted = run_transient(
        &pair,
        &owner,
        &thermal,
        &routes[0],
        duration,
        HertzTransientControls {
            maximum_evaluations: 1,
            ..path_controls()
        },
        source,
        &CancelGate::new(),
    );
    assert!(matches!(exhausted, Err(RigidHertzHeatError::Mechanical(_))));
    assert_eq!(pair.bodies, saved);
    assert_eq!(owner, original);
    assert_eq!(thermal.cells(), &saved_cells);
    assert_eq!(thermal.time_s(), 0.0);
    let retry = run_transient(
        &pair,
        &owner,
        &thermal,
        &routes[0],
        duration,
        path_controls(),
        source,
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(retry.bodies_after, reference.bodies_after);
    assert_eq!(retry.thermal.cells(), reference.thermal.cells());
    assert_eq!(retry.mechanical, reference.mechanical);
}
