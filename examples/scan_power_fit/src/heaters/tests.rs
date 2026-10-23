use super::*;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

const INPUT: &str = "FS_SCAN_HEATERS_V1
heater 0 hold 2 2 2
heater 1 linear 2 2 2
knot 0 0 0.5 0.5 0.5 2 fit 2
knot 1 0 0.5 0.5 0.5 0 fixed
knot 1 0.5 0.5 0.5 0.5 6 fixed
knot 0 1 0.5 0.5 0.5 5 fit 3
knot 0 2 0.5 0.5 0.5 0 fixed
knot 1 2 0.5 0.5 0.5 0 fixed
";
const LIMITS: HeatingProgramBudget = HeatingProgramBudget {
    max_bytes: 4096, max_sources: 2, max_knots: 6,
    max_fitted_powers: 2, max_authored_controls: 30,
};
fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&gate, &Cx::new(&gate, arena, StreamKey { seed: 76, kernel_id: 51, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic));
    });
}

#[test]
fn independent_heaters_keep_exact_events_control_identity_and_absorbed_power() {
    with_cx(|_, cx| {
        let specimen = ScanSpecimen::synthetic_reference(cx).unwrap();
        let program = ScanHeatingProgram::from_text(cx, &specimen.mesh, 3.0, INPUT, LIMITS).unwrap();
        assert_eq!(program.event_times_s(), &[0.5, 1.0, 2.0, 3.0]);
        assert_eq!(program.end_time_s(), 3.0);
        assert_eq!(program.controls(), &[
            ScanPowerControl { source: 0, value: 0, reference_w: 2.0 },
            ScanPowerControl { source: 0, value: 1, reference_w: 3.0 }]);
        let count: usize = program.sources().iter().map(|s| s.schedule_control_count().unwrap()).sum();
        assert_eq!(program.authored_control_count(), count);
        assert_eq!(count, 30);
        for (i, expected) in [(0, 5.0), (1, 4.0)] {
            let sample = program.sources()[i].at(cx, &specimen.mesh, 1.0).unwrap();
            assert!((sample.delivered_power_w - expected).abs() < 1e-12);
        }
        assert_eq!(program.sources()[0].before(cx, &specimen.mesh, 1.0).unwrap().absorbed_power_w, 2.0);
        for source in program.sources() {
            assert_eq!(source.at(cx, &specimen.mesh, 3.0).unwrap().absorbed_power_w, 0.0);
        }
    });
}

#[test]
fn export_accepts_selected_values_but_never_changes_fixed_heating_or_shutdown() {
    with_cx(|_, cx| {
        let specimen = ScanSpecimen::synthetic_reference(cx).unwrap();
        let program = ScanHeatingProgram::from_text(cx, &specimen.mesh, 2.0, INPUT, LIMITS).unwrap();
        let mut sources = program.sources().to_vec();
        sources[0] = sources[0].with_absorbed_power_values(&[8.0, 16.0, 0.0], 0.0, 3).unwrap();
        let exported = program.export(cx, &sources).unwrap();
        assert_eq!(exported[0].powers_w, vec![8.0, 16.0, 0.0]);
        assert_eq!(exported[1].powers_w, vec![0.0, 6.0, 0.0]);
        assert_eq!(exported[0].fitted, vec![true, true, false]);
        sources[1] = sources[1].with_absorbed_power_values(&[0.0, 7.0, 0.0], 0.0, 3).unwrap();
        assert!(program.export(cx, &sources).is_err());
        sources[1] = program.sources()[1].clone();
        sources[0] = sources[0].with_absorbed_power_values(&[8.0, 16.0, 1e-30], 0.0, 3).unwrap();
        assert!(program.export(cx, &sources).is_err());
        assert!(program.export(cx, &sources[..1]).is_err());
    });
}

#[test]
fn declarations_enforce_source_knot_control_caps_and_cancel_before_physics() {
    with_cx(|gate, cx| {
        let specimen = ScanSpecimen::synthetic_reference(cx).unwrap();
        for budget in [HeatingProgramBudget { max_sources: 1, ..LIMITS },
            HeatingProgramBudget { max_knots: 5, ..LIMITS },
            HeatingProgramBudget { max_fitted_powers: 1, ..LIMITS },
            HeatingProgramBudget { max_authored_controls: 29, ..LIMITS },
            HeatingProgramBudget { max_bytes: 16, ..LIMITS }] {
            assert!(ScanHeatingProgram::from_text(cx, &specimen.mesh, 2.0, INPUT, budget).is_err());
        }
        for text in [INPUT.replace("heater 0", "heater 2"), INPUT.replace("hold", "cubic"),
            INPUT.replace("2 fit 2", "-2 fit 2"), INPUT.replace("fit 2", "fit 0"),
            INPUT.replace("knot 0 0 ", "knot 0 0.25 "), INPUT.replace("knot 0 2 ", "knot 0 1 "),
            INPUT.replace("0 fixed", "0 fit 1"), INPUT.replace("2 fit 2", "NaN fit 2"),
            INPUT.replace("5 fit 3", "5 fixed").replace("2 fit 2", "2 fixed")] {
            assert!(ScanHeatingProgram::from_text(cx, &specimen.mesh, 2.0, &text, LIMITS).is_err(), "{text}");
        }
        assert!(ScanHeatingProgram::from_text(cx, &specimen.mesh, 1.5, INPUT, LIMITS).is_err());
        gate.request();
        assert!(ScanHeatingProgram::from_text(cx, &specimen.mesh, 2.0, INPUT, LIMITS).is_err());
    });
}

#[test]
fn legacy_program_keeps_the_exact_original_source_and_selected_scales() {
    with_cx(|_, cx| {
        let specimen = ScanSpecimen::synthetic_reference(cx).unwrap();
        let program = ScanHeatingProgram::from_specimen(cx, &specimen, [0.0, 4.0]).unwrap();
        assert_eq!(program.sources()[0], specimen.source([0.0, 4.0]).unwrap());
        assert_eq!(program.event_times_s(), &[1.0, 2.0]);
        assert_eq!(program.controls().len(), 2);
        let exported = program.export(cx, program.sources()).unwrap();
        assert_eq!(exported[0].times_s, vec![0.0, 1.0, 2.0]);
        assert_eq!(exported[0].powers_w, vec![0.0, 4.0, 0.0]);
    });
}
