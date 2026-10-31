//! The two escape reactions share one original exterior and conservative FV receipts.
use super::*;

const FIRST_FEEDBACK: &str = "reservoir_feedback aperture 0 -20 20 16 16\n";
const SECOND_FEEDBACK: &str = "reservoir_feedback aperture 1 -20 20 16 16\n";
const SOLVE: &str =
    "reservoir_feedback_solve 12 1e-5 1e-8 0 4 8 1e-5 0.5 0.9 1.5 1e-4 0.5 1e-6\n";

fn joint_source() -> String {
    multiple_escape_source(false)
        .replace("clock 96000 1", "clock 8000 12")
        .replace("output 256 20", "output 64 20")
        .replace("hold 64 20", "hold 16 2")
        .replace("hold 128 -20", "hold 32 -2")
        .replace("hold 256 0", "hold 64 0")
        .replace("coupling 64 9216", "coupling 256 442368")
        .replace(
            SECOND_ESCAPE,
            &format!("{SECOND_ESCAPE}{FIRST_FEEDBACK}{SECOND_FEEDBACK}{SOLVE}"),
        )
}

#[test]
fn both_escape_reactions_use_the_full_shared_matrix_and_original_fv_work() {
    let text = joint_source();
    let mut coupled = read(&text, 1).into_renderer();
    let mut baseline = read(
        &text.replace(FIRST_FEEDBACK, "").replace(SECOND_FEEDBACK, "").replace(SOLVE, ""),
        1,
    ).into_renderer();
    assert_eq!(coupled.system().reservoir_radiations().len(), 2);
    assert!(coupled.system().reservoir_radiation().is_none());
    assert!(coupled.system().radiation_solve().is_some());
    assert!(coupled.system().load().independent_mouth_response().is_err());
    let mut changed_flow = false;
    let mut reacted = [false; 2];
    let mut mutual_work = [false; 2];
    for _ in 0..64 {
        let response = coupled.system().load().independent_mouths_response().unwrap();
        assert_eq!(response.history_pressure_pa.len(), 3);
        assert_eq!(response.impedance_pa_s_m3.len(), 9);
        assert_ne!(response.impedance_pa_s_m3[5], 0.0);
        assert_ne!(response.impedance_pa_s_m3[7], 0.0);
        let before_step = coupled.system().accepted_steps();
        let mass = coupled.system().zone().total_mass_kg();
        let energy = coupled.system().zone().total_energy_j();
        baseline.block(&mut [0.0]).unwrap();
        coupled.block(&mut [0.0]).unwrap();
        assert_eq!(coupled.system().accepted_steps(), before_step + 1);
        let frame = coupled.last_frame().unwrap();
        assert!(frame.reservoir_radiation.is_none());
        assert_eq!(frame.reservoir_radiations.len(), 2);
        let dt = coupled.system().load().propagator().dt();
        let flows = [
            frame.inlet_flow_m3_s,
            frame.reservoir_radiations[0].outward_flow_m3_s,
            frame.reservoir_radiations[1].outward_flow_m3_s,
        ];
        let pressures = [
            frame.inlet_pressure_pa,
            frame.reservoir_radiations[0].radiation_pressure_pa,
            frame.reservoir_radiations[1].radiation_pressure_pa,
        ];
        for row in 0..3 {
            let mut expected = response.history_pressure_pa[row];
            let mut scale = expected.abs();
            for (column, flow) in flows.iter().enumerate() {
                let contribution = response.impedance_pa_s_m3[row * 3 + column] * flow;
                expected += contribution;
                scale += contribution.abs();
            }
            assert!((pressures[row] - expected).abs()
                <= 1e-10 * scale.max(1.0), "original shared response row {row}");
        }
        assert!(frame.pressure_residual_pa.abs() < 1e-5);
        for (index, mouth) in frame.reservoir_radiations.iter().enumerate() {
            let resolved = &mouth.resolved_exchange;
            assert_eq!(resolved.aperture, index);
            assert_eq!(Some(resolved), frame.reservoir_exchanges.iter()
                .find(|exchange| exchange.aperture == index));
            assert_eq!(mouth.outward_flow_m3_s, resolved.normal_volume_flow_m3_s);
            assert_eq!(mouth.mean_acoustic_pressure_pa, resolved.mean_acoustic_pressure_pa);
            assert!(mouth.pressure_residual_pa.abs() < 1e-5);
            assert_eq!(mouth.acoustic_work_residual_j,
                resolved.acoustic_pressure_work_out_j - mouth.radiation_work_j);
            assert_eq!(mouth.projection_covariance_j,
                resolved.acoustic_pressure_work_out_j
                    - mouth.mean_acoustic_pressure_pa * (mouth.outward_flow_m3_s * dt));
            assert!((mouth.radiation_work_j
                - mouth.radiation_pressure_pa * (mouth.outward_flow_m3_s * dt)).abs()
                <= 64.0 * f64::EPSILON * mouth.radiation_work_j.abs().max(1e-30));
            reacted[index] |= mouth.radiation_pressure_pa.abs() > 1e-9;
            let other = 2 - index;
            mutual_work[index] |=
                (response.impedance_pa_s_m3[(index + 1) * 3 + other] * flows[other]).abs()
                    > 1e-12;
        }
        let mass_residual = frame.flow_mass_change_kg + frame.reservoir_mass_out_kg
            + frame.faces.iter().map(|face| face.mass_out_kg).sum::<f64>();
        assert!(mass_residual.abs() < 5e-12 * mass);
        let supplied = frame.flow.open_energy_j + frame.flow.open_boundary_energy_j
            + frame.flow.wall_heat_j + frame.flow.source_energy_j;
        assert!((frame.flow_energy_change_j - supplied).abs() < 5e-11 * energy.abs());
        changed_flow |= frame.inlet_flow_m3_s != baseline.last_frame().unwrap().inlet_flow_m3_s;
    }
    assert!(changed_flow, "joint reaction must reach the original resolved source");
    assert!(reacted.into_iter().all(|present| present));
    assert!(mutual_work.into_iter().all(|present| present),
        "both escape equations must retain the other escape's pressure contribution");
}

#[test]
fn plural_feedback_requires_distinct_observed_owners_and_explicit_valid_solver() {
    let text = joint_source();
    for invalid in [
        text.replace(SOLVE, ""),
        text.replace(SECOND_FEEDBACK, FIRST_FEEDBACK),
        text.replace(SECOND_FEEDBACK, "reservoir_feedback aperture 2 -20 20 16 16\n"),
        text.replace(FIRST_ESCAPE, ""),
        text.replace(FIRST_FEEDBACK, ""),
        text.replace(FIRST_FEEDBACK, "").replace(SECOND_FEEDBACK, ""),
        text.replace("reservoir_feedback_solve 12 1e-5", "reservoir_feedback_solve 0 1e-5"),
        text.replace("reservoir_feedback_solve 12 1e-5", "reservoir_feedback_solve 12 0"),
        text.replace("0.5 1e-6\n", "0.5\n"),
        text.replace("reservoir_feedback aperture 0 -20 20", "reservoir_feedback aperture 0 20 -20"),
    ] {
        assert!(ResolvedDgFilePerformance::from_bytes(
            invalid.as_bytes(), 1, Arc::new(CancelGate::new_clock_free())
        ).is_err());
    }
    // Keeping one row preserves the established singleton arithmetic and receipts.
    let mut singleton = read(&text.replace(FIRST_FEEDBACK, "").replace(SOLVE, ""), 1)
        .into_renderer();
    assert!(singleton.system().reservoir_radiations().is_empty());
    assert!(singleton.system().radiation_solve().is_none());
    assert_eq!(singleton.system().reservoir_radiation_aperture(), Some(1));
    singleton.block(&mut [0.0]).unwrap();
    assert!(singleton.last_frame().unwrap().reservoir_radiation.is_some());
    assert!(singleton.last_frame().unwrap().reservoir_radiations.is_empty());
}
