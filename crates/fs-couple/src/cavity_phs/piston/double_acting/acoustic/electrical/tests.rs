//! Native circuit, two-gas impulse and complete DG histories share actual work.
//! Synthetic material sources exercise numerical composition, not measurements.

use super::*;
use crate::cavity_phs::piston::double_acting::acoustic::tests::{
    assert_load_same, assert_piston_same, caloric_energy, load, motion, near,
};
use crate::cavity_phs::piston::double_acting::electrical::tests::{
    circuit, electrical, input, mechanical, piston, solve,
};

fn faces(loads: &[DgPortLoad], areas: [f64; 2]) -> [AcousticMotionFace<'_>; 2] {
    [0, 1].map(|i| AcousticMotionFace {
        load: &loads[i],
        signed_area_m2: areas[i],
    })
}

fn stored(
    piston: &DoubleActingGasPiston,
    circuit: &crate::material_circuit::MaterialCircuit,
) -> f64 {
    caloric_energy(piston.primary().gas())
        + caloric_energy(piston.opposing())
        + 0.5 * piston.primary().mass_kg() * piston.primary().velocity_m_s().powi(2)
        + 0.5 * circuit.charges_c()[0].powi(2) / 0.01
}

fn assert_same(actual: &AcousticElectricalCandidate, expected: &AcousticElectricalCandidate) {
    assert_piston_same(actual.native().piston(), expected.native().piston());
    assert_eq!(
        actual.native().piston_step(),
        expected.native().piston_step()
    );
    assert_eq!(
        actual.native().circuit().circuit(),
        expected.native().circuit().circuit()
    );
    let a = actual.native().circuit().step();
    let b = expected.native().circuit().step();
    assert_eq!(a.native.x, b.native.x);
    assert_eq!(a.native.y, b.native.y);
    assert_eq!(a.resistors, b.resistors);
    assert_eq!(
        [
            a.stored_before_j,
            a.stored_after_j,
            a.heat_j,
            a.energy_residual_j
        ],
        [
            b.stored_before_j,
            b.stored_after_j,
            b.heat_j,
            b.energy_residual_j
        ]
    );
    assert_eq!(
        actual.native().motion().back_emf_v,
        expected.native().motion().back_emf_v
    );
    assert_eq!(
        actual.native().motion().current_a,
        expected.native().motion().current_a
    );
    assert_eq!(
        actual.native().motion().transfer_residual_j,
        expected.native().motion().transfer_residual_j
    );
    assert_eq!(actual.motion(), expected.motion());
    assert_eq!(
        actual.total_energy_residual_j(),
        expected.total_energy_residual_j()
    );
    for (a, b) in actual.loads().iter().zip(expected.loads()) {
        assert_load_same(a, b);
    }
}

#[test]
fn g1_acoustic_reaction_changes_native_current_and_closes_the_complete_energy_window() {
    let mut p = piston(0.03);
    let mut c = circuit(2.0, 0.2, true);
    let dt = input().duration_s;
    let mut loads = vec![load(0.15, dt), load(0.85, dt)];
    let areas = [0.001, -0.0007];
    let initial_energy = stored(&p, &c);
    let mut external = 0.0;
    let mut resistor_heat = 0.0;
    let mut exterior_work = 0.0;
    for sample in 0..12 {
        let drive = DoublePistonInput {
            heat_j: if sample < 3 { [0.02, 0.01] } else { [0.0; 2] },
            external_force_n: if sample < 4 { 0.3 } else { 0.0 },
            ..input()
        };
        let inputs = [0.0, if sample < 4 { 5.0 } else { 0.0 }];
        let before_p = p.clone();
        let before_c = c.clone();
        let before_loads = loads.clone();
        let candidate = p
            .electroacoustic_candidate(
                &faces(&loads, areas),
                electrical(&c, &inputs, 2.0),
                drive,
                mechanical(),
                solve(),
                motion(),
                p.primary().velocity_m_s(),
                || false,
            )
            .unwrap();
        let native = candidate.native();
        let em = native.motion();
        let sound = candidate.motion();
        let piston_step = native.piston_step();
        let dx = piston_step.displacement_m;
        let velocity = dx / dt;
        near(em.back_emf_v, 2.0 * velocity, 2e-9);
        near(em.force_n, 2.0 * em.current_a, 1e-12);
        near(em.electrical_work_j + em.force_n * dx, 0.0, 1e-10);
        near(
            sound.mechanical_work_j + sound.faces.iter().map(|f| f.inlet_work_j).sum::<f64>(),
            0.0,
            1e-10,
        );
        let gas_force = 0.01
            * (piston_step.midpoint[0].state().pressure - piston_step.midpoint[1].state().pressure);
        near(
            p.primary().mass_kg()
                * (native.piston().primary().velocity_m_s() - p.primary().velocity_m_s()),
            dt * (gas_force + drive.external_force_n + em.force_n + sound.force_n),
            2e-9,
        );
        if sample == 0 {
            let unloaded = p
                .electrical_candidate(
                    electrical(&c, &inputs, 2.0),
                    drive,
                    mechanical(),
                    solve(),
                    p.primary().velocity_m_s(),
                    || false,
                )
                .unwrap();
            assert!(sound.force_n < 0.0);
            assert!(
                native.piston().primary().velocity_m_s()
                    < unloaded.piston().primary().velocity_m_s()
            );
            assert!(em.back_emf_v < unloaded.motion().back_emf_v);
            assert!(em.current_a > unloaded.motion().current_a);
            // Independent midpoint RC law: loading changes emf, hence current.
            near(
                em.current_a - unloaded.motion().current_a,
                (unloaded.motion().back_emf_v - em.back_emf_v) / (2.0 + dt / (2.0 * 0.01)),
                1e-10,
            );
        }
        for i in 0..2 {
            let q = areas[i] * velocity;
            let mut replay = before_loads[i].clone();
            let output = replay.preview_step(q).unwrap().accept();
            assert_eq!(sound.faces[i].duct, output);
            assert_load_same(&candidate.loads()[i], &replay);
            assert_eq!(candidate.loads()[i].state().step_count(), sample + 1);
            near(sound.faces[i].inlet_flow_m3_s, q, 1e-17);
        }
        assert_piston_same(&p, &before_p);
        assert_eq!(c, before_c);
        for (a, b) in loads.iter().zip(&before_loads) {
            assert_load_same(a, b);
        }
        // The reserved back-emf source is internal; only port 1 is supplied work.
        external += dt * inputs[1] * native.circuit().step().native.y[1]
            + drive.heat_j.iter().sum::<f64>()
            + drive.external_force_n * dx;
        resistor_heat += native.circuit().step().heat_j;
        exterior_work += sound.terminal_work_j;
        near(
            candidate.total_energy_residual_j(),
            0.0,
            motion().energy_tolerance_j,
        );
        let (native, next_loads, _, _) = candidate.into_parts();
        let (next_c, next_p, _, _, _, _) = native.into_parts();
        p = next_p;
        c = next_c;
        loads = next_loads;
        let acoustic_energy: f64 = loads
            .iter()
            .map(|l| l.stored_energy_j() + l.cumulative_dissipated_j())
            .sum();
        near(
            stored(&p, &c) - initial_energy + resistor_heat + acoustic_energy + exterior_work,
            external,
            (sample + 1) as f64 * motion().energy_tolerance_j,
        );
    }
    assert!(loads.iter().all(|l| l.stored_energy_j() > 0.0));
}

#[test]
fn g3_zero_area_matches_the_original_electrical_solve_and_continues_acoustic_history() {
    let p = piston(0.03);
    let c = circuit(2.0, 0.2, false);
    let loads = vec![
        load(0.15, input().duration_s),
        load(0.85, input().duration_s),
    ];
    let warm = p
        .electroacoustic_candidate(
            &faces(&loads, [0.001, -0.0007]),
            electrical(&c, &[0.0], 2.0),
            input(),
            mechanical(),
            solve(),
            motion(),
            0.03,
            || false,
        )
        .unwrap();
    let (native, loads, _, _) = warm.into_parts();
    let (c, p, _, _, _, _) = native.into_parts();
    let expected = p
        .electrical_candidate(
            electrical(&c, &[0.0], 2.0),
            input(),
            mechanical(),
            solve(),
            p.primary().velocity_m_s(),
            || false,
        )
        .unwrap();
    let candidate = p
        .electroacoustic_candidate(
            &faces(&loads, [0.0; 2]),
            electrical(&c, &[0.0], 2.0),
            input(),
            mechanical(),
            solve(),
            motion(),
            p.primary().velocity_m_s(),
            || false,
        )
        .unwrap();
    assert_piston_same(candidate.native().piston(), expected.piston());
    assert_eq!(candidate.native().piston_step(), expected.piston_step());
    assert_eq!(
        candidate.native().circuit().circuit(),
        expected.circuit().circuit()
    );
    let actual_circuit = candidate.native().circuit().step();
    let expected_circuit = expected.circuit().step();
    assert_eq!(actual_circuit.native.x, expected_circuit.native.x);
    assert_eq!(actual_circuit.native.y, expected_circuit.native.y);
    assert_eq!(actual_circuit.resistors, expected_circuit.resistors);
    assert_eq!(actual_circuit.heat_j, expected_circuit.heat_j);
    assert_eq!(
        candidate.native().total_energy_residual_j(),
        expected.total_energy_residual_j()
    );
    assert_eq!(candidate.motion().force_n, 0.0);
    assert_eq!(candidate.motion().mechanical_work_j, 0.0);
    for (i, original) in loads.iter().enumerate() {
        let mut free = original.clone();
        let output = free.preview_step(0.0).unwrap().accept();
        assert_eq!(candidate.motion().faces[i].duct, output);
        assert_load_same(&candidate.loads()[i], &free);
        assert_ne!(candidate.loads()[i].state(), original.state());
        assert_ne!(
            candidate.loads()[i].terminal().unwrap().pole_states(),
            original.terminal().unwrap().pole_states()
        );
    }
}

#[test]
fn g4_budget_source_and_final_cancellation_preserve_all_owners_and_retry_suffix() {
    let p = piston(0.03);
    let c = circuit(2.0, 0.2, false);
    let loads = vec![
        load(0.15, input().duration_s),
        load(0.85, input().duration_s),
    ];
    let before_p = p.clone();
    let before_c = c.clone();
    let before_loads = loads.clone();
    let areas = [0.001, -0.0007];
    assert!(
        p.electroacoustic_candidate(
            &faces(&loads, areas),
            ElectricalMotionInput {
                resistors: &[],
                ..electrical(&c, &[0.0], 2.0)
            },
            input(),
            mechanical(),
            solve(),
            motion(),
            0.03,
            || false
        )
        .is_err()
    );
    assert!(
        p.electroacoustic_candidate(
            &faces(&loads, areas),
            electrical(&c, &[0.0], 2.0),
            input(),
            mechanical(),
            ElectricalMotionSolve {
                maximum_trials: 1,
                ..solve()
            },
            motion(),
            0.03,
            || false
        )
        .is_err()
    );
    let mut polls = 0;
    let expected = p
        .electroacoustic_candidate(
            &faces(&loads, areas),
            electrical(&c, &[0.0], 2.0),
            input(),
            mechanical(),
            solve(),
            motion(),
            0.03,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut count = 0;
    let refused = p.electroacoustic_candidate(
        &faces(&loads, areas),
        electrical(&c, &[0.0], 2.0),
        input(),
        mechanical(),
        solve(),
        motion(),
        0.03,
        || {
            count += 1;
            count == polls
        },
    );
    assert!(matches!(refused, Err(AcousticMotionError::Cancelled)));
    assert_eq!(count, polls);
    assert_piston_same(&p, &before_p);
    assert_eq!(c, before_c);
    for (a, b) in loads.iter().zip(&before_loads) {
        assert_load_same(a, b);
    }
    let retried = p
        .electroacoustic_candidate(
            &faces(&loads, areas),
            electrical(&c, &[0.0], 2.0),
            input(),
            mechanical(),
            solve(),
            motion(),
            0.03,
            || false,
        )
        .unwrap();
    assert_same(&retried, &expected);
    let (a, loads_a, _, _) = retried.into_parts();
    let (b, loads_b, _, _) = expected.into_parts();
    let (c_a, p_a, _, _, _, _) = a.into_parts();
    let (c_b, p_b, _, _, _, _) = b.into_parts();
    let next_a = p_a
        .electroacoustic_candidate(
            &faces(&loads_a, areas),
            electrical(&c_a, &[0.0], 2.0),
            input(),
            mechanical(),
            solve(),
            motion(),
            p_a.primary().velocity_m_s(),
            || false,
        )
        .unwrap();
    let next_b = p_b
        .electroacoustic_candidate(
            &faces(&loads_b, areas),
            electrical(&c_b, &[0.0], 2.0),
            input(),
            mechanical(),
            solve(),
            motion(),
            p_b.primary().velocity_m_s(),
            || false,
        )
        .unwrap();
    assert_same(&next_a, &next_b);
}
