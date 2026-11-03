use super::*;

fn controls() -> AuditControls {
    AuditControls {
        force_tolerance_n: 1.0e-8,
        gap_tolerance_m: 1.0e-9,
        slip_tolerance_m: 1.0e-9,
        velocity_tolerance_m_s: 1.0e-8,
        energy_tolerance_j: 1.0e-8,
    }
}

fn support() -> ContactAuditInput {
    ContactAuditInput {
        id: 9,
        gap_before_m: 0.0,
        gap_after_m: 0.0,
        slip_increment_m: [0.0; 3],
        relative_velocity_before_m_s: [0.0; 3],
        relative_velocity_after_m_s: [0.0; 3],
        force_n: [20.0, 0.0, 0.0],
        friction_coefficient: 0.25,
        prescribed_velocity_m_s: [0.0; 3],
        activation_tolerance_m: 1.0e-9,
    }
}

fn body(
    kinetic_change_j: f64,
    authored_work_j: f64,
    native_total_external_work_j: f64,
) -> BodyWork {
    BodyWork {
        kinetic_change_j,
        bulk_change_j: 0.0,
        cohesive_change_j: 0.0,
        material_loss_j: 0.0,
        authored_work_j,
        prescribed_jump_work_j: 0.0,
        native_total_external_work_j,
        native_energy_defect_j: 0.0,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 1.0e-8,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn loaded_rest_stick_slide_and_breakaway_balance_independent_body_work() {
    // Independent 2 kg, gravity 10 m/s² fixtures at h=.4 s. Each authored
    // load, actual velocity, slip and body work is supplied separately.
    for (load, before, after, tangent_force, slip, loss) in [
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (4.0, 0.0, 0.0, -4.0, 0.0, 0.0),
        (0.0, 3.0, 2.0, -5.0, 1.0, 5.0),
        (6.0, 0.0, 0.2, -5.0, 0.04, 0.2),
    ] {
        let contact = ContactAuditInput {
            force_n: [20.0, tangent_force, 0.0],
            slip_increment_m: [0.0, slip, 0.0],
            relative_velocity_before_m_s: [0.0, before, 0.0],
            relative_velocity_after_m_s: [0.0, after, 0.0],
            ..support()
        };
        let kinetic = after * after - before * before;
        let work = body(kinetic, load * slip, (load + tangent_force) * slip);
        let report = audit(&[contact], &[work], 0.4, controls(), &mut || false).unwrap();
        close(report.friction_loss_j, loss);
        close(report.contact_body_work_j, -loss);
        close(report.energy_residual_j, 0.0);
        close(report.contact_work_residual_j, 0.0);
        close(report.normal_constraint_work_j, 0.0);
    }
}

#[test]
fn cold_original_newton_candidate_passes_independent_physical_support_audit() {
    use fs_contact::sustained::{
        NaturalMapScales,
        solve::{SolveControls, solve},
    };
    let mut solve_controls = SolveControls::default();
    solve_controls.newton.absolute_tolerance = 1.0e-11;
    solve_controls.newton.relative_tolerance = 0.0;
    solve_controls.difference_step_n = 1.0e-3;
    solve_controls.force_residual_tolerance_n = 1.0e-9;
    let scales = NaturalMapScales {
        normal_n_per_m: 25.0,
        tangent_n_per_m: 25.0,
    };
    for (load, before) in [(0.0, 0.0), (4.0, 0.0), (0.0, 3.0), (6.0, 0.0)] {
        // Exact constant-force midpoint mechanics. The generic solver uses
        // its original Newton/FGMRES owner, then the separate physical audit
        // checks the original response records. No native body is imitated.
        let candidate = solve(
            &[0.25],
            &[scales],
            &[[0.0; 3]],
            solve_controls,
            &mut |forces, _| {
                let f = forces[0];
                Ok::<_, &'static str>(vec![[
                    0.04 * (f[0] - 20.0),
                    0.4 * before + 0.04 * (load + f[1]),
                    0.04 * f[2],
                ]])
            },
            &mut || false,
        )
        .unwrap();
        let f = candidate.forces_n[0];
        let path = candidate.paths_m[0];
        let endpoint = [
            0.2 * (f[0] - 20.0),
            before + 0.2 * (load + f[1]),
            0.2 * f[2],
        ];
        let contact = ContactAuditInput {
            gap_after_m: path[0],
            force_n: f,
            slip_increment_m: path,
            relative_velocity_before_m_s: [0.0, before, 0.0],
            relative_velocity_after_m_s: endpoint,
            ..support()
        };
        let kinetic = endpoint.iter().map(|v| v * v).sum::<f64>() - before * before;
        let authored = -20.0 * path[0] + load * path[1];
        let total = (f[0] - 20.0) * path[0] + (load + f[1]) * path[1] + f[2] * path[2];
        let checked = audit(
            &[contact],
            &[body(kinetic, authored, total)],
            0.4,
            controls(),
            &mut || false,
        );
        assert!(
            checked.is_ok(),
            "cold load={load}, speed={before}, force={f:?}, path={path:?}: {checked:?}"
        );
    }
}

#[test]
fn moving_support_work_is_distinct_from_friction_and_normal_constraint_work() {
    // Wall tangent speed +2 m/s, body 0→1 m/s, body displacement .2 m.
    // Relative slip is -.6 m. The wall supplies 4 J: 1 J kinetic + 3 J heat.
    let contact = ContactAuditInput {
        force_n: [20.0, 5.0, 0.0],
        slip_increment_m: [0.0, -0.6, 0.0],
        relative_velocity_before_m_s: [0.0, -2.0, 0.0],
        relative_velocity_after_m_s: [0.0, -1.0, 0.0],
        prescribed_velocity_m_s: [0.0, -2.0, 0.0],
        ..support()
    };
    let report = audit(
        &[contact],
        &[body(1.0, 0.0, 1.0)],
        0.4,
        controls(),
        &mut || false,
    )
    .unwrap();
    close(report.surface_work_j, 4.0);
    close(report.friction_loss_j, 3.0);
    close(report.contact_body_work_j, 1.0);
    close(report.relative_contact_work_j, -3.0);
    close(report.energy_residual_j, 0.0);

    // Upward normal wall and body speed 1 m/s under 20 N gravity. Relative
    // normal displacement is zero; actual wall work balances -8 J gravity.
    let contact = ContactAuditInput {
        prescribed_velocity_m_s: [-1.0, 0.0, 0.0],
        ..support()
    };
    let report = audit(
        &[contact],
        &[body(0.0, -8.0, 0.0)],
        0.4,
        controls(),
        &mut || false,
    )
    .unwrap();
    close(report.surface_work_j, 8.0);
    close(report.contact_body_work_j, 8.0);
    close(report.normal_constraint_work_j, 0.0);
    close(report.friction_loss_j, 0.0);
}

#[test]
fn released_body_keeps_its_original_outward_momentum() {
    // 2 kg released upward under 1 N, h=.4: q=.04 m, v=.2 m/s, ΔK=.04 J.
    let contact = ContactAuditInput {
        gap_after_m: 0.04,
        slip_increment_m: [0.04, 0.0, 0.0],
        relative_velocity_after_m_s: [0.2, 0.0, 0.0],
        force_n: [0.0; 3],
        ..support()
    };
    let report = audit(
        &[contact],
        &[body(0.04, 0.04, 0.04)],
        0.4,
        controls(),
        &mut || false,
    )
    .unwrap();
    close(report.friction_loss_j, 0.0);
    close(report.contact_body_work_j, 0.0);
    close(report.energy_residual_j, 0.0);
}

#[test]
fn resolved_energy_alone_cannot_admit_reversal_or_a_false_static_stop() {
    // Constant -5 N friction acting for 1.5 s would reverse v=3 to -.75.
    // Although the midpoint work balances exactly, this is not an admissible
    // uninterrupted Coulomb sliding interval.
    let contact = ContactAuditInput {
        force_n: [20.0, -5.0, 0.0],
        slip_increment_m: [0.0, 1.6875, 0.0],
        relative_velocity_before_m_s: [0.0, 3.0, 0.0],
        relative_velocity_after_m_s: [0.0, -0.75, 0.0],
        ..support()
    };
    assert!(matches!(
        audit(
            &[contact],
            &[body(-8.4375, 0.0, -8.4375)],
            1.5,
            controls(),
            &mut || false
        ),
        Err(AuditError::UnresolvedReversal { contact_id: 9 })
    ));
    // Zero total slip and zero work can also conceal equal/opposite endpoint
    // velocities; a force inside the disk does not make this a static path.
    let contact = ContactAuditInput {
        force_n: [20.0, -4.0, 0.0],
        relative_velocity_before_m_s: [0.0, 3.0, 0.0],
        relative_velocity_after_m_s: [0.0, -3.0, 0.0],
        ..support()
    };
    assert!(matches!(
        audit(
            &[contact],
            &[body(0.0, 0.0, 0.0)],
            3.0,
            controls(),
            &mut || false
        ),
        Err(AuditError::UnresolvedReversal { .. })
    ));
    // With zero dry-friction capacity, a native elastic turnaround has no
    // Coulomb direction to violate and may keep its original midpoint motion.
    let frictionless = ContactAuditInput {
        friction_coefficient: 0.0,
        force_n: [20.0, 0.0, 0.0],
        ..contact
    };
    assert!(
        audit(
            &[frictionless],
            &[body(0.0, 0.0, 0.0)],
            3.0,
            controls(),
            &mut || false
        )
        .is_ok()
    );
}

#[test]
fn tiny_tension_or_penetration_is_not_hidden_by_acceptance_tolerances() {
    for contact in [
        ContactAuditInput {
            gap_after_m: -1.0e-15,
            ..support()
        },
        ContactAuditInput {
            gap_before_m: -1.0e-15,
            ..support()
        },
        ContactAuditInput {
            force_n: [-1.0e-15, 0.0, 0.0],
            ..support()
        },
    ] {
        assert!(matches!(
            audit(
                &[contact],
                &[body(0.0, 0.0, 0.0)],
                0.4,
                controls(),
                &mut || false
            ),
            Err(AuditError::Contact { tolerance: 0.0, .. })
        ));
    }
}

#[test]
fn actual_kinematics_conjugate_work_and_absolute_body_defects_gate_adoption() {
    let wrong_trace = ContactAuditInput {
        slip_increment_m: [0.0, 0.01, 0.0],
        ..support()
    };
    assert!(matches!(
        audit(
            &[wrong_trace],
            &[body(0.0, 0.0, 0.0)],
            0.4,
            controls(),
            &mut || false
        ),
        Err(AuditError::Contact {
            quantity: "actual midpoint kinematics [m]",
            ..
        })
    ));

    let mut false_receipt = body(0.0, 0.0, 1.0);
    false_receipt.bulk_change_j = 1.0;
    assert!(matches!(
        audit(&[support()], &[false_receipt], 0.4, controls(), &mut || {
            false
        }),
        Err(AuditError::Energy {
            quantity: "actual contact/body work",
            ..
        })
    ));

    // Equal/opposite native defects cannot disappear through aggregation.
    let mut positive = body(0.0, 0.0, 0.0);
    positive.native_energy_defect_j = 1.0e-4;
    let mut negative = positive;
    negative.native_energy_defect_j = -1.0e-4;
    assert!(matches!(
        audit(
            &[support()],
            &[positive, negative],
            0.4,
            controls(),
            &mut || false
        ),
        Err(AuditError::Energy {
            quantity: "absolute native body defects",
            ..
        })
    ));

    // The same rule applies when an original receipt falsely reports zero but
    // independently reconstructed body balances expose cancelling defects.
    let positive = body(0.0, 0.0, 1.0e-4);
    let negative = body(0.0, 0.0, -1.0e-4);
    assert!(matches!(
        audit(
            &[support()],
            &[positive, negative],
            0.4,
            controls(),
            &mut || false
        ),
        Err(AuditError::Energy {
            quantity: "absolute native body defects",
            ..
        })
    ));
}

#[test]
fn normal_constraint_drift_remains_a_defect_and_cancellation_returns_no_receipt() {
    let drift = 1.0e-7;
    let contact = ContactAuditInput {
        gap_before_m: drift,
        slip_increment_m: [-drift, 0.0, 0.0],
        relative_velocity_after_m_s: [-5.0e-7, 0.0, 0.0],
        activation_tolerance_m: 1.0e-6,
        ..support()
    };
    let mut work = body(0.0, 0.0, -20.0 * drift);
    work.bulk_change_j = -20.0 * drift;
    let permissive_velocity = AuditControls {
        velocity_tolerance_m_s: 1.0e-6,
        ..controls()
    };
    assert!(matches!(
        audit(&[contact], &[work], 0.4, permissive_velocity, &mut || false),
        Err(AuditError::Energy {
            quantity: "normal constraint work",
            ..
        })
    ));
    assert!(matches!(
        audit(
            &[support()],
            &[body(0.0, 0.0, 0.0)],
            0.4,
            controls(),
            &mut || true
        ),
        Err(AuditError::Cancelled)
    ));
}
