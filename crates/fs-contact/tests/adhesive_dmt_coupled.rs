//! G1/G3/G4: original adhesion and source junction shear on complete mechanics.
#[path = "adhesive_dmt_coupled/common.rs"]
mod common;
use common::*;
use fs_contact::normal_patch::dmt::coupled::*;
use fs_exec::CancelGate;

#[test]
fn g1_analytic_real_area_normal_work_and_contact_entry_exit_match_independent_quadrature() {
    let (n, s) = sources(0.01, 1e7, 6e6);
    let owner = DmtCoupledLoadLine::new(&n, &s, compliance(), config()).unwrap();
    for (a, b) in [
        (-1.4 * Z, 0.6 * Z),
        (0.6 * Z, -1.4 * Z),
        (-0.8 * Z, -0.2 * Z),
        (0.2 * Z, 0.8 * Z),
        (0.3 * Z, 0.3 * Z),
    ] {
        let path = owner.path(a, b).unwrap();
        let mut mean = [0.0; 3];
        let steps = 40_000;
        for i in 0..steps {
            let t = (f64::from(i) + 0.5) / f64::from(steps);
            let r = n.evaluate((1.0 - t) * a + t * b).unwrap();
            for (j, v) in [
                r.normal_force_n,
                r.contact_area_m2,
                r.repulsive_normal_force_n,
            ]
            .into_iter()
            .enumerate()
            {
                mean[j] += v / f64::from(steps);
            }
        }
        close(path.mean_normal_force_n, mean[0], 1e-17, 3e-8);
        close(path.mean_real_area_m2, mean[1], 1e-30, 3e-8);
        close(path.mean_repulsive_force_n, mean[2], 1e-17, 3e-8);
        close(
            path.normal_work_j,
            path.mean_normal_force_n * (b - a),
            1e-30,
            1e-14,
        );
        close(
            path.normal_work_j,
            n.evaluate(b).unwrap().reversible_energy_j - n.evaluate(a).unwrap().reversible_energy_j,
            1e-28,
            2e-13,
        );
        close(
            path.static_capacity_n,
            1e7 * path.mean_real_area_m2,
            1e-30,
            2e-15,
        );
        close(
            path.kinetic_capacity_n,
            6e6 * path.mean_real_area_m2,
            1e-30,
            2e-15,
        );
    }
}

#[test]
fn g1_full_off_diagonal_compliance_preserves_real_static_normal_shear_equilibrium() {
    let (n, s) = sources(0.01, 1e7, 6e6);
    let c = compliance();
    let mut owner = DmtCoupledLoadLine::new(&n, &s, c, config()).unwrap();
    let old = owner.initialize([0.25 * Z, 0.1 * Z, -0.2 * Z]).unwrap();
    let q = [0.35 * Z, old.position_m()[1], old.position_m()[2]];
    let path = owner.path(old.position_m()[0], q[0]).unwrap();
    let force = [
        path.mean_normal_force_n,
        0.4 * path.static_capacity_n,
        -0.2 * path.static_capacity_n,
    ];
    let step = owner
        .trial(old, free(c, q, force), 1e-6, &CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(step.regime, DmtCoupledRegime::Sticking);
    assert_eq!(step.friction_dissipation_j, 0.0);
    for i in 0..3 {
        close(step.state().position_m()[i], q[i], 5e-21, 2e-10);
        close(step.mean_force_n[i], force[i], 1e-16, 3e-8);
    }
    assert!(step.mean_force_n[1] != 0.0 && step.mean_force_n[2] != 0.0);
}

#[test]
fn g1_source_shear_is_integrated_with_normal_attachment_and_detachment() {
    let (n, s) = sources(0.01, 1e7, 6e6);
    let c = compliance();
    let mut owner = DmtCoupledLoadLine::new(&n, &s, c, config()).unwrap();
    let gate = CancelGate::new_clock_free();
    for (a, b) in [
        (-0.4 * Z, 0.5 * Z),
        (0.6 * Z, -0.4 * Z),
        (-0.6 * Z, -0.4 * Z),
    ] {
        let old = owner.initialize([a, 0.0, 0.0]).unwrap();
        let q = [b, 0.3 * Z, 0.4 * Z];
        let path = owner.path(a, b).unwrap();
        let force = [
            path.mean_normal_force_n,
            0.6 * path.kinetic_capacity_n,
            0.8 * path.kinetic_capacity_n,
        ];
        let step = owner.trial(old, free(c, q, force), 1e-6, &gate).unwrap();
        for i in 0..3 {
            close(step.state().position_m()[i], q[i], 5e-21, 2e-9);
            close(step.mean_force_n[i], force[i], 1e-16, 3e-8);
        }
        close(
            step.friction_dissipation_j,
            path.kinetic_capacity_n * 0.5 * Z,
            1e-28,
            3e-9,
        );
        if b < 0.0 {
            assert_eq!(step.state().regime(), DmtCoupledRegime::Released);
        }
        if a < 0.0 && b < 0.0 {
            assert_eq!(step.mean_force_n[1..], [0.0, 0.0]);
            assert!(step.mean_force_n[0] < 0.0);
        }
        let actual = owner
            .advance_actual(
                old,
                step.state().position_m(),
                step.mean_force_n,
                step.regime,
                1e-6,
                &gate,
            )
            .unwrap();
        assert_eq!(actual.state(), step.state());
        close(
            actual.friction_dissipation_j,
            step.friction_dissipation_j,
            1e-30,
            1e-14,
        );
    }
}

#[test]
fn g3_occupied_sliding_uses_kinetic_capture_before_restoring_static_strength() {
    let (n, s) = sources(0.01, 1e7, 4e6);
    let c = [[2e-5, 0.0, 0.0], [0.0, 3e-5, 0.0], [0.0, 0.0, 3e-5]];
    let mut owner = DmtCoupledLoadLine::new(&n, &s, c, config()).unwrap();
    let gate = CancelGate::new_clock_free();
    let initial = owner.initialize([0.4 * Z, 0.0, 0.0]).unwrap();
    let p = owner.path(0.4 * Z, 0.4 * Z).unwrap();
    let moved = [0.4 * Z, 0.1 * Z, 0.0];
    let occupied = owner
        .trial(
            initial,
            free(c, moved, [p.mean_normal_force_n, p.kinetic_capacity_n, 0.0]),
            1e-6,
            &gate,
        )
        .unwrap()
        .state();
    assert_eq!(occupied.regime(), DmtCoupledRegime::Sliding);
    let free_next = free(
        c,
        occupied.position_m(),
        [p.mean_normal_force_n, 0.7 * p.static_capacity_n, 0.0],
    );
    let fresh = owner.initialize(occupied.position_m()).unwrap();
    let initially_moving = owner
        .initialize_with_velocity(occupied.position_m(), [0.01, 0.0])
        .unwrap();
    assert_eq!(initially_moving.regime(), DmtCoupledRegime::Sliding);
    assert_eq!(
        owner
            .trial(initially_moving, free_next, 1e-6, &gate)
            .unwrap()
            .regime,
        DmtCoupledRegime::Sliding
    );
    assert_eq!(
        owner.trial(fresh, free_next, 1e-6, &gate).unwrap().regime,
        DmtCoupledRegime::Sticking
    );
    let continuing = owner.trial(occupied, free_next, 1e-6, &gate).unwrap();
    assert_eq!(continuing.regime, DmtCoupledRegime::Sliding);
    assert!(continuing.state().position_m()[1] > occupied.position_m()[1]);
    let capture = free(
        c,
        continuing.state().position_m(),
        [p.mean_normal_force_n, 0.2 * p.static_capacity_n, 0.0],
    );
    let stopped = owner
        .trial(continuing.state(), capture, 1e-6, &gate)
        .unwrap();
    assert_eq!(stopped.state().regime(), DmtCoupledRegime::Sticking);
    let static_again = free(
        c,
        stopped.state().position_m(),
        [p.mean_normal_force_n, 0.7 * p.static_capacity_n, 0.0],
    );
    assert_eq!(
        owner
            .trial(stopped.state(), static_again, 1e-6, &gate)
            .unwrap()
            .regime,
        DmtCoupledRegime::Sticking
    );
    // The position-only path cannot prove that a dynamic endpoint is at rest.
    // An actual moving endpoint must keep kinetic capacity despite zero slip.
    let still_moving = owner
        .with_endpoint_velocity(stopped.state(), [-0.01, 0.0])
        .unwrap();
    assert_eq!(still_moving.position_m(), stopped.state().position_m());
    assert_eq!(still_moving.regime(), DmtCoupledRegime::Sliding);
    assert_eq!(
        owner
            .trial(still_moving, static_again, 1e-6, &gate)
            .unwrap()
            .regime,
        DmtCoupledRegime::Sliding
    );
    assert_eq!(
        owner
            .with_endpoint_velocity(still_moving, [0.0; 2])
            .unwrap()
            .regime(),
        DmtCoupledRegime::Sticking
    );
}

#[test]
fn g3_tangent_basis_rotation_preserves_force_motion_and_work_with_cross_coupling() {
    let (n, s) = sources(0.01, 1e7, 6e6);
    let c = compliance();
    let gate = CancelGate::new_clock_free();
    let r = [[1.0, 0.0, 0.0], [0.0, 0.6, -0.8], [0.0, 0.8, 0.6]];
    let cr = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            (0..3)
                .flat_map(|k| (0..3).map(move |l| r[i][k] * c[k][l] * r[j][l]))
                .sum()
        })
    });
    let mut a = DmtCoupledLoadLine::new(&n, &s, c, config()).unwrap();
    let mut b = DmtCoupledLoadLine::new(&n, &s, cr, config()).unwrap();
    let old = [0.25 * Z, 0.1 * Z, -0.2 * Z];
    let target = [0.5 * Z, 0.4 * Z, 0.2 * Z];
    let p = a.path(old[0], target[0]).unwrap();
    let f = [
        p.mean_normal_force_n,
        0.6 * p.kinetic_capacity_n,
        0.8 * p.kinetic_capacity_n,
    ];
    let qfree = free(c, target, f);
    let first = a
        .trial(a.initialize(old).unwrap(), qfree, 1e-6, &gate)
        .unwrap();
    let second = b
        .trial(
            b.initialize(mul(r, old)).unwrap(),
            mul(r, qfree),
            1e-6,
            &gate,
        )
        .unwrap();
    for i in 0..3 {
        close(
            second.state().position_m()[i],
            mul(r, first.state().position_m())[i],
            1e-20,
            1e-9,
        );
        close(
            second.mean_force_n[i],
            mul(r, first.mean_force_n)[i],
            1e-16,
            1e-8,
        );
    }
    close(
        first.friction_dissipation_j,
        second.friction_dissipation_j,
        1e-28,
        1e-9,
    );
}

#[test]
fn g3_sliding_can_be_admissible_when_the_hypothetical_all_stick_root_is_not() {
    let (n, s) = sources(0.01, 1e7, 6e6);
    let c = compliance();
    let mut owner = DmtCoupledLoadLine::new(&n, &s, c, config()).unwrap();
    let old = owner.initialize([0.8 * Z, 0.0, 0.0]).unwrap();
    let target = [0.9 * Z, -2.9 * Z, 0.0];
    let path = owner.path(old.position_m()[0], target[0]).unwrap();
    let force = [path.mean_normal_force_n, -path.kinetic_capacity_n, 0.0];
    let qfree = free(c, target, force);
    let result = owner
        .trial(old, qfree, 1e-6, &CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(result.regime, DmtCoupledRegime::Sliding);
    for i in 0..3 {
        close(result.state().position_m()[i], target[i], 1e-21, 1e-9);
        close(result.mean_force_n[i], force[i], 1e-16, 1e-8);
    }
    // The static Schur equation is strictly increasing and still negative at
    // the upper source boundary: no all-stick root is admitted in this case.
    let determinant = c[1][1] * c[2][2] - c[1][2] * c[2][1];
    let inverse = [
        [c[2][2] / determinant, -c[1][2] / determinant],
        [-c[2][1] / determinant, c[1][1] / determinant],
    ];
    let elimination = |v: [f64; 2]| {
        [
            inverse[0][0] * v[0] + inverse[0][1] * v[1],
            inverse[1][0] * v[0] + inverse[1][1] * v[1],
        ]
    };
    let cross = elimination([c[1][0], c[2][0]]);
    let free_tangent = elimination([qfree[1], qfree[2]]);
    let schur = c[0][0] - c[0][1] * cross[0] - c[0][2] * cross[1];
    let effective = qfree[0] - c[0][1] * free_tangent[0] - c[0][2] * free_tangent[1];
    assert!(
        Z - effective
            + schur
                * owner
                    .path(old.position_m()[0], Z)
                    .unwrap()
                    .mean_normal_force_n
            < 0.0
    );
}

#[test]
fn g0_g4_sources_stability_actual_force_and_failed_trials_preserve_history() {
    let (n, s) = sources(0.01, 1e7, 6e6);
    let (other, other_s) = sources(0.02, 1e7, 6e6);
    let c = compliance();
    assert!(matches!(
        DmtCoupledLoadLine::new(&n, &other_s, c, config()),
        Err(DmtCoupledError::StateMismatch)
    ));
    assert!(matches!(
        DmtCoupledLoadLine::new(&n, &s, c.map(|row| row.map(|v| v * 1e8)), config()),
        Err(DmtCoupledError::MechanicalStability)
    ));
    let mut owner = DmtCoupledLoadLine::new(&n, &s, c, config()).unwrap();
    let alien = DmtCoupledLoadLine::new(&other, &other_s, c, config())
        .unwrap()
        .initialize([0.1 * Z, 0.0, 0.0])
        .unwrap();
    let state = owner.initialize([0.2 * Z, 0.0, 0.0]).unwrap();
    let gate = CancelGate::new_clock_free();
    let target = [0.3 * Z, 0.2 * Z, 0.0];
    let p = owner.path(state.position_m()[0], target[0]).unwrap();
    let free = free(
        c,
        target,
        [p.mean_normal_force_n, p.kinetic_capacity_n, 0.0],
    );
    let expected = owner.trial(state, free, 1e-6, &gate).unwrap();
    assert!(matches!(
        owner.trial(alien, free, 1e-6, &gate),
        Err(DmtCoupledError::StateMismatch)
    ));
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert_eq!(
        owner.trial(state, free, 1e-6, &cancelled),
        Err(DmtCoupledError::Cancelled)
    );
    assert!(owner.trial(state, [1.0; 3], 1e-6, &gate).is_err());
    let mut limited = config();
    limited.maximum_evaluations = 1;
    assert!(matches!(
        DmtCoupledLoadLine::new(&n, &s, c, limited)
            .unwrap()
            .trial(state, free, 1e-6, &gate),
        Err(DmtCoupledError::Budget)
    ));
    let mut wrong = expected.mean_force_n;
    wrong[1] *= 2.0;
    assert!(
        owner
            .advance_actual(state, target, wrong, DmtCoupledRegime::Sliding, 1e-6, &gate)
            .is_err()
    );
    assert_eq!(owner.trial(state, free, 1e-6, &gate).unwrap(), expected);
    assert!(std::ptr::eq(owner.normal(), &n) && std::ptr::eq(owner.shear(), &s));
}
