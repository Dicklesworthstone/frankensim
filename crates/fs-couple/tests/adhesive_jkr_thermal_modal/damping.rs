//! Physical support loss, finite-time heating and joint acceptance.
use super::*;

fn oscillator(k: f64, c: f64, extension: f64, time: f64) -> (f64, f64) {
    let decay = c / (2.0 * MU);
    let omega = (k / MU - decay * decay).sqrt();
    let envelope = extension * (-decay * time).exp();
    (
        envelope * ((omega * time).cos() + decay / omega * (omega * time).sin()),
        -envelope * k / (MU * omega) * (omega * time).sin(),
    )
}

fn oscillator_temperature(k: f64, c: f64, extension: f64, time: f64, capacity: f64) -> f64 {
    // Independently integrate the continuous closed-form oscillator power and
    // exact cooling kernel, without consuming the numerical mechanical trace.
    let n = 8192;
    let h = time / n as f64;
    let mut total = 0.0;
    for i in 0..=n {
        let at = i as f64 * h;
        let velocity = oscillator(k, c, extension, at).1;
        let weight = if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        total += weight * c * velocity * velocity * (-ALPHA * (time - at)).exp();
    }
    total * h / (3.0 * capacity)
}

#[test]
fn original_damped_oscillator_heats_during_motion_and_closes_all_destinations() {
    let law = source();
    let (d, e, _) = scales(&law);
    let k = 0.02 * e / (d * d);
    let c = 0.4 * (k * MU).sqrt();
    let extension = 2.0 * d;
    let initial_contact = law.detached_state(-20.0 * d).unwrap();
    let make = || mechanical_with_damping(&law, initial_contact, 0.0, Some((k, extension, c)));
    let caloric = curves(&law, [1.0; 2]);
    let origin = 7.0 * DT;
    let mut bare = make();
    let mut owner = joint_with_support_heat(
        make(),
        timeline(&law, &caloric, [1.0; 2], origin),
        &caloric,
        origin,
    );
    let initial = owner.stored_energy_j().unwrap();
    let (mut outside, mut external, mut loss) = (0.0, 0.0, 0.0);
    for _ in 0..40 {
        let frame = owner.step(&[0.0; 2]).unwrap().clone();
        assert_eq!(frame.mechanical, *bare.step(&[0.0; 2]).unwrap());
        assert_eq!(states(owner.mechanical()), states(&bare));
        assert_eq!(frame.heat.snap.work.dissipated_work_j(), 0.0);
        let support = frame.support_heat.as_ref().unwrap();
        assert!(support.work.dissipated_work_j() > 0.0);
        close(
            support.work.dissipated_work_j(),
            frame.mechanical.damping_dissipation_j,
            e,
            1e-10,
        );
        outside += support.work.other_work_j();
        loss += support.work.dissipated_work_j();
        external += frame.external_transport_heat_j;
        close(
            owner.stored_energy_j().unwrap() + outside - external,
            initial,
            e,
            3e-7,
        );
        assert!(
            frame.support_transfer_residual_j.abs() <= controls().heat.transfer_energy_tolerance_j
        );
        assert!(frame.combined_energy_residual_j.abs() <= frame.combined_energy_tolerance_j);
    }
    let time = owner.mechanical().state().time_s();
    let (exact_extension, exact_velocity) = oscillator(k, c, extension, time);
    close(
        owner.mechanical().state().contact_state().indentation_m()
            - initial_contact.indentation_m()
            + extension,
        exact_extension,
        d,
        2e-4,
    );
    close(
        owner.mechanical().state().velocity_m_s(),
        exact_velocity,
        d / DT,
        2e-4,
    );
    let heat = oscillator_temperature(k, c, extension, time, release(&law));
    for (side, fraction) in [0.6, 0.1].into_iter().enumerate() {
        close(
            owner.thermal().temperature_k(CellId(side as u64)).unwrap(),
            AMBIENT + fraction * heat,
            1.0,
            5e-5,
        );
    }
    assert!(loss > 0.0 && external < 0.0);
    close(outside, 0.3 * loss, e, 1e-10);
    assert_eq!(owner.thermal().time_s(), origin + time);
}

#[test]
fn viscous_support_and_adhesive_snaps_retain_distinct_heat_and_outside_accounts() {
    let law = source();
    let (d, e, _) = scales(&law);
    let k = 0.5 * e / (d * d);
    let c = 0.1 * (k * MU).sqrt();
    let caloric = curves(&law, [1.0; 2]);
    for (attached, extension, expected) in [
        (false, -4.0 * d, JkrJumpKind::Attach),
        (true, 30.0 * d, JkrJumpKind::Detach),
    ] {
        let state = if attached {
            law.evaluate_attached(0.0).unwrap().state()
        } else {
            law.detached_state(-d).unwrap()
        };
        let mut owner = joint_with_support_heat(
            mechanical_with_damping(&law, state, 0.0, Some((k, extension, c))),
            timeline(&law, &caloric, [1.0; 2], 0.0),
            &caloric,
            0.0,
        );
        let initial = owner.stored_energy_j().unwrap();
        let (mut outside, mut external) = (0.0, 0.0);
        let mut seen = false;
        for _ in 0..80 {
            let frame = owner.step(&[0.0; 2]).unwrap().clone();
            let support = frame.support_heat.as_ref().unwrap();
            outside += frame.heat.snap.work.other_work_j() + support.work.other_work_j();
            external += frame.external_transport_heat_j;
            close(
                owner.stored_energy_j().unwrap() + outside - external,
                initial,
                e,
                8e-7,
            );
            assert_eq!(frame.snap_transfer_residual_j, 0.0);
            seen |= frame
                .heat
                .events
                .iter()
                .any(|event| event.jump.kind == expected);
            if seen {
                assert!(support.work.dissipated_work_j() > 0.0);
                assert!(frame.heat.snap.work.dissipated_work_j() > 0.0);
                break;
            }
        }
        assert!(
            seen,
            "original spring/dashpot trajectory must reach {expected:?}"
        );
        assert!(owner.thermal().energy_j() > 0.0);
    }
}

#[test]
fn support_heat_requires_a_route_and_thermal_refusal_rolls_back_damped_mechanics() {
    let law = source();
    let (d, e, _) = scales(&law);
    let k = 0.02 * e / (d * d);
    let c = 0.4 * (k * MU).sqrt();
    let state = law.detached_state(-20.0 * d).unwrap();
    let make = || mechanical_with_damping(&law, state, 0.0, Some((k, 2.0 * d, c)));
    let caloric = curves(&law, [1.0; 2]);
    assert!(matches!(
        JkrThermalModalSystem::new(
            make(),
            timeline(&law, &caloric, [1.0; 2], 0.0),
            [CellId(0), CellId(1)],
            [&caloric[0], &caloric[1]],
            HeatPartition::new(0.25, 0.5, 0.25).unwrap(),
            Some("outgoing elastic waves".into()),
            0.0,
            controls(),
            &CancelGate::new_clock_free(),
        ),
        Err(JkrThermalModalError::Invalid(_))
    ));
    let narrow = curves(&law, [1.0, 1e-10]);
    let mut failed = joint_with_support_heat(
        make(),
        timeline(&law, &narrow, [1.0, 1e-10], 0.0),
        &narrow,
        0.0,
    );
    let bodies = states(failed.mechanical());
    let token = failed.mechanical().state();
    let stamp = failed.heat_history().stamp();
    let cells = failed.thermal().cells().clone();
    assert!(matches!(
        failed.step(&[0.0; 2]),
        Err(JkrThermalModalError::Thermal(_))
    ));
    assert_eq!(states(failed.mechanical()), bodies);
    assert_eq!(failed.mechanical().state(), token);
    assert_eq!(failed.heat_history().stamp(), stamp);
    assert_eq!(failed.thermal().cells(), &cells);
    assert_eq!(failed.thermal().time_s(), 0.0);
    assert!(failed.last_frame().is_none());
    assert!(failed.mechanical().network().last_frame().is_none());
    let mut owner = joint_with_support_heat(
        make(),
        timeline(&law, &caloric, [1.0; 2], 0.0),
        &caloric,
        0.0,
    );
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(owner.step_under_gate(&[0.0; 2], &gate).is_err());
    assert_eq!(owner.mechanical().samples_rendered(), 0);
    assert_eq!(owner.thermal().time_s(), 0.0);
    let accepted = owner.step(&[0.0; 2]).unwrap().clone();
    let mut fresh = joint_with_support_heat(
        make(),
        timeline(&law, &caloric, [1.0; 2], 0.0),
        &caloric,
        0.0,
    );
    assert_eq!(accepted, *fresh.step(&[0.0; 2]).unwrap());
    assert_eq!(owner.thermal().cells(), fresh.thermal().cells());
    assert_eq!(states(owner.mechanical()), states(fresh.mechanical()));
}
