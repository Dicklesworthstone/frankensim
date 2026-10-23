use super::*;
use crate::heaters::HeatingProgramBudget;
use crate::specimen::{ScanSpecimen, SpecimenBudget};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

const MODEL: &str = include_str!("../../examples/specimens/rectangular.scan");
const HEATERS: &str = include_str!("../../examples/specimens/multiple.heaters");
const LIMITS: RequirementBudget = RequirementBudget { max_bytes: 4096, max_declarations: 16,
    max_kkt_dimension: 64, max_peak_samples: 128, max_gradient_entries: 1024, max_transpose_solves: 128 };
const INPUT: &str = "FS_SCAN_REQUIREMENTS_V1
power 0 1 * 14 1
energy 0 4 48 1
state 0.5 0 specific_enthalpy_j_kg * 154 1
peak 1 4 liquid_mass_fraction 0.7 1 0.001 all
slew 1 50 20 1
speed 1 0.1 1
";
fn with_program(f: impl FnOnce(&CancelGate, &Cx<'_>, &ScanSpecimen, &ScanHeatingProgram)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 76, kernel_id: 51, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        let model = ScanSpecimen::from_text(&cx, MODEL, SpecimenBudget {
            max_bytes: 4096, max_vertices: 8, max_elements: 6, max_phase_knots: 4 }).unwrap();
        let program = ScanHeatingProgram::from_text(&cx, &model.mesh, model.end_time_s, HEATERS,
            HeatingProgramBudget { max_bytes: 4096, max_sources: 3, max_knots: 9,
                max_fitted_powers: 3, max_authored_controls: 45 }).unwrap();
        f(&gate, &cx, &model, &program);
    });
}

#[test]
fn native_lowering_preserves_physical_units_and_counts_all_refined_peak_samples() {
    with_program(|_, cx, model, program| {
        let requirements = ScanRequirements::from_text(cx, program, &model.mesh, INPUT, LIMITS).unwrap();
        assert_eq!(requirements.event_times_s(), &[0.5, 1.0, 4.0]);
        let bound = requirements.bind(cx, &[0.5, 1., 2., 2.5, 3., 4.]).unwrap();
        assert_eq!(bound.powers, vec![ScanControlBound { coordinate: 1,
            lower: None, upper: Some(14.), scale: 1. }]);
        assert_eq!(bound.energy[0], ScanEnergyCap { start_time_s: 0., end_time_s: 4., maximum_j: 48., scale_j: 1. });
        assert_eq!(bound.states[0].response.endpoint, 0);
        assert_eq!(bound.states[0].response.observable, ScanObservable::SpecificEnthalpy);
        assert_eq!(bound.peaks[0].response.first_endpoint, 1);
        assert_eq!(bound.peaks[0].response.last_endpoint, 5);
        assert_eq!(bound.peaks[0].response.vertices, (0..8).collect::<Vec<_>>());
        assert_eq!(bound.max_peak_samples, 40);
        assert_eq!(bound.response_budget.max_responses, 2);
        assert_eq!(bound.response_budget.max_gradient_entries, 106);
        assert_eq!(bound.response_budget.max_transpose_solves, 13);
        assert_eq!(bound.actuators.max_power_segments, 3);
        assert_eq!(bound.actuators.path_budget.max_segments, 3);
        assert_eq!(bound.kkt_dimension, 19);
        assert!(requirements.bind(cx, &[0.5, 1.1, 2., 3., 4.]).is_err());
        let none = ScanRequirements::none(program, &model.mesh, LIMITS).bind(cx, &[1.,2.,4.]).unwrap();
        assert_eq!(none.kkt_dimension, 6);
        assert_eq!(none.response_budget.max_responses, 0);
    });
}

#[test]
fn invalid_requirements_never_become_ignored_limits_or_held_signal_chords() {
    with_program(|_, cx, model, program| {
        for row in ["power 0 2 * 1 1", "power 2 0 * 1 1", "power 7 0 * 1 1",
            "power 0 0 * * 1", "power 0 0 2 2 1", "power 0 0 -1 2 1",
            "power 0 0 * 2 1\npower 0 0 * 3 1", "energy 2 1 4 1",
            "energy 0 5 4 1", "energy 0 4 -1 1", "state 0 0 temperature_k * 311 1",
            "state 1 8 temperature_k * 311 1", "state 1 0 liquid_mass_fraction * 60 1",
            "state 1 0 temperature * 311 1", "peak 2 1 temperature_k 311 1 1 all",
            "peak 1 4 temperature_k 311 1 0 all", "peak 1 4 temperature_k 311 1 1 7;0",
            "peak 1 4 temperature_k 311 1 1 0;0", "slew 1 NaN 2 1", "speed 1 1 0",
            "slew 1 2 2 1\nslew 1 3 3 1", "power 0 0 * 2 1 extra", "penalty 100"] {
            let text = format!("FS_SCAN_REQUIREMENTS_V1\n{row}\n");
            let error = ScanRequirements::from_text(cx, program, &model.mesh, &text, LIMITS).err().unwrap();
            assert!(error.to_string().contains("requirement line"), "{row}: {error}");
        }
        for row in ["slew 0 50 50 1", "speed 0 1 1"] {
            let text = format!("FS_SCAN_REQUIREMENTS_V1\n{row}\n");
            let requirements = ScanRequirements::from_text(cx, program, &model.mesh, &text, LIMITS).unwrap();
            assert!(requirements.bind(cx, &[1.,2.,3.,4.]).is_err());
        }
        assert!(ScanRequirements::from_text(cx, program, &model.mesh, "FS_SCAN_REQUIREMENTS_V1", LIMITS).is_err());
    });
}

#[test]
fn native_response_work_and_cancellation_are_admitted_before_any_physical_solve() {
    with_program(|gate, cx, model, program| {
        for budget in [RequirementBudget { max_bytes: 8, ..LIMITS }, RequirementBudget { max_declarations: 5, ..LIMITS }] {
            assert!(ScanRequirements::from_text(cx, program, &model.mesh, INPUT, budget).is_err());
        }
        for budget in [RequirementBudget { max_kkt_dimension: 18, ..LIMITS },
            RequirementBudget { max_peak_samples: 39, ..LIMITS },
            RequirementBudget { max_gradient_entries: 105, ..LIMITS },
            RequirementBudget { max_transpose_solves: 12, ..LIMITS }] {
            let requirements = ScanRequirements::from_text(cx, program, &model.mesh, INPUT, budget).unwrap();
            assert!(requirements.bind(cx, &[0.5,1.,2.,2.5,3.,4.]).is_err());
        }
        let requirements = ScanRequirements::from_text(cx, program, &model.mesh, INPUT, LIMITS).unwrap();
        for times in [vec![], vec![0.,4.], vec![1.,1.,4.], vec![1.,f64::NAN,4.], vec![1.,3.]] {
            assert!(requirements.bind(cx, &times).is_err());
        }
        gate.request();
        assert!(requirements.bind(cx, &[0.5,1.,2.,4.]).is_err());
        assert!(ScanRequirements::from_text(cx, program, &model.mesh, INPUT, LIMITS).is_err());
    });
}
