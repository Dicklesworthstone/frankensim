use super::*;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

const DATA: &str = include_str!("../../examples/specimens/rectangular.scan");
const CAP: SpecimenBudget = SpecimenBudget {
    max_bytes: 8192, max_vertices: 8, max_elements: 6, max_phase_knots: 4,
};
fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&gate, &Cx::new(&gate, arena,
            StreamKey { seed: 76, kernel_id: 51, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic));
    });
}

#[test]
fn supplied_mesh_mass_clock_and_exact_shutdown_reach_original_source() {
    with_cx(|_, cx| {
        let model = ScanSpecimen::from_text(cx, DATA, CAP).unwrap();
        assert_eq!(model.mesh.vertex_count(), 8);
        assert_eq!(model.mesh.element_count(), 6);
        let mass: f64 = (0..6).map(|i| model.mesh.element_volume(i) * model.reference_density_kg_m3).sum();
        assert!((mass - 4.0).abs() < 1e-14);
        assert_eq!((model.switch_time_s, model.end_time_s), (1.5, 4.0));
        let source = model.source([8.0, 20.0]).unwrap();
        for (time, expected) in [(0.0, 8.0), (1.5, 20.0), (4.0, 0.0)] {
            let sampled = source.at(cx, &model.mesh, time).unwrap();
            assert_eq!(sampled.absorbed_power_w, expected);
            assert!((sampled.delivered_power_w - expected).abs() < 1e-12);
        }
        assert_eq!(source.before(cx, &model.mesh, 1.5).unwrap().absorbed_power_w, 8.0);
        assert_eq!(source.before(cx, &model.mesh, 4.0).unwrap().absorbed_power_w, 20.0);
        assert_eq!(source.absorbed_power_values(3).unwrap()[2].to_bits(), 0.0_f64.to_bits());
    });
}

#[test]
fn explicit_robin_and_input_identity_do_not_reuse_adiabatic_defaults() {
    with_cx(|_, cx| {
        let original = ScanSpecimen::from_text(cx, DATA, CAP).unwrap();
        let text = DATA.replace("exterior adiabatic", "exterior robin 7.5 293.15");
        let altered = ScanSpecimen::from_text(cx, &text, CAP).unwrap();
        assert_ne!(original.identity, altered.identity);
        assert_eq!(altered.exterior, ThermalBc::robin(7.5, 293.15).unwrap());
        assert_eq!(original.exterior, ThermalBc::adiabatic());
        assert_eq!(altered.reference_density_kg_m3, original.reference_density_kg_m3);
    });
}

#[test]
fn missing_ambiguous_invalid_and_oversized_inputs_refuse_before_fitting() {
    with_cx(|gate, cx| {
        for text in [
            DATA.replace("FS_SCAN_SPECIMEN_V1", "FS_SCAN_SPECIMEN_V2"),
            DATA.replace("conductivity_w_m_k 2\n", ""),
            format!("{DATA}\nconductivity_w_m_k 3\n"),
            DATA.replace("exterior adiabatic", "exterior robin 0 293.15"),
            DATA.replace("exterior adiabatic", "exterior robin 2 -1"),
            DATA.replace("reference_density_kg_m3 2", "reference_density_kg_m3 NaN"),
            DATA.replace("switch_time_s 1.5", "switch_time_s 4"),
            DATA.replace("half_width_m 3 2 2", "half_width_m 0 2 2"),
            DATA.replace("initial_specific_enthalpy_j_kg 150", "initial_specific_enthalpy_j_kg 901"),
            DATA.replace("phase_knot 200 310 1", "phase_knot 200 310 100"),
            DATA.replace("tet 0 1 3 7", "tet 0 0 3 7"),
            DATA.replace("tet 0 1 3 7", "tet 0 1 3 99"),
            DATA.replace("tet 0 5 1 7", "tet 0 1 3 7"),
            DATA.replace("vertex_m 2 1 1", "vertex_m 2 1 0"),
            format!("{DATA}\nundeclared_field 1\n"),
        ] {
            assert!(ScanSpecimen::from_text(cx, &text, CAP).is_err(), "{text}");
        }
        for budget in [SpecimenBudget { max_bytes: 1, ..CAP }, SpecimenBudget { max_vertices: 7, ..CAP },
            SpecimenBudget { max_elements: 5, ..CAP }, SpecimenBudget { max_phase_knots: 3, ..CAP }]
        { assert!(ScanSpecimen::from_text(cx, DATA, budget).is_err()); }
        let orphan = DATA.replace("tet 0 1 3 7", "vertex_m 4 4 4\ntet 0 1 3 7");
        assert!(ScanSpecimen::from_text(cx, &orphan, SpecimenBudget { max_vertices: 9, ..CAP }).is_err());
        gate.request();
        assert!(ScanSpecimen::from_text(cx, DATA, CAP).is_err());
    });
}
