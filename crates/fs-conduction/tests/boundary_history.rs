//! Scenario interpolation and physical boundary-lift regressions.
mod support;

use fs_conduction::bc::{ThermalBc, ThermalBoundaryBuilder};
use fs_conduction::bc::history::{ThermalBoundaryHistory, ThermalRegionHistory};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::{ConductionError, ConductionMesh, HEAT_FLUX_DIMS, HTC_DIMS, TEMPERATURE_DIMS};
use fs_qty::{Dims, QtyAny};
use fs_scenario::bc::{BcKind, BcValue, BoundaryCondition, Physics};
use fs_scenario::scenario::Environment;
use fs_scenario::signal::{Interp, TimeSignal};
use support::{with_cancelled_cx, with_cx};

fn mesh() -> ConductionMesh {
    let (complex, positions) = box_grid([2, 2, 2], [1.0, 1.0, 1.0]);
    ConductionMesh::new(complex, positions).unwrap()
}
fn ramp(from: f64, to: f64, dims: Dims) -> TimeSignal {
    TimeSignal::Ramp { t_start: 0.0, t_end: 2.0,
        from: QtyAny::new(from, dims), to: QtyAny::new(to, dims) }
}
fn row(name: &str, kind: BcKind, signal: TimeSignal) -> BoundaryCondition {
    BoundaryCondition { region: name.into(), physics: Physics::Thermal, kind,
        value: Some(BcValue::Signal(signal)), compatibility: None, frame: 0 }
}
fn temperature_row(name: &str) -> ThermalRegionHistory {
    ThermalRegionHistory { row: row(name, BcKind::Dirichlet, ramp(300.0, 400.0, TEMPERATURE_DIMS)),
        reference_temperature: None }
}
fn left(mesh: &ConductionMesh) -> ThermalBoundaryBuilder<'_> {
    ThermalBoundaryBuilder::new(mesh)
        .region("left", |f| on_box_face(f.centroid[0], 0.0), ThermalBc::dirichlet(300.0).unwrap())
        .unwrap().adiabatic_remainder()
}

#[test]
fn g0_explicit_time_evaluation_preserves_the_steady_refusal() {
    let env = Environment::earth_lab();
    let row = temperature_row("left").row;
    assert!(ThermalBc::from_scenario_row(&row, &env, None).is_err());
    for (time, value) in [(-1.0, 300.0), (1.0, 350.0), (3.0, 400.0)] {
        assert_eq!(ThermalBc::from_scenario_row_at(&row, &env, None, time).unwrap(),
            ThermalBc::dirichlet(value).unwrap());
    }
    assert!(ThermalBc::from_scenario_row_at(&row, &env, None, f64::NAN).is_err());
}

#[test]
fn g1_all_thermal_histories_rebuild_prescribed_values_without_retagging() {
    let mesh = mesh();
    let env = Environment::earth_lab();
    let history = with_cx(|cx| {
        let builder = left(&mesh)
            .region("right", |f| on_box_face(f.centroid[0], 1.0), ThermalBc::robin(10.0, 290.0).unwrap()).unwrap()
            .region("bottom", |f| on_box_face(f.centroid[1], 0.0), ThermalBc::neumann(0.0).unwrap()).unwrap();
        ThermalBoundaryHistory::new(cx, builder, &env, vec![temperature_row("left"),
            ThermalRegionHistory { row: row("right", BcKind::Robin, ramp(10.0, 30.0, HTC_DIMS)),
                reference_temperature: Some(ramp(290.0, 310.0, TEMPERATURE_DIMS)) },
            ThermalRegionHistory { row: row("bottom", BcKind::Neumann, TimeSignal::Table {
                times: vec![0.0, 1.0, 2.0], values: vec![0.0, -100.0, -200.0],
                dims: HEAT_FLUX_DIMS, interp: Interp::Hold }), reference_temperature: None },
        ], 0.0).unwrap()
    });
    with_cx(|cx| {
        let start = history.at(cx, 0.0).unwrap();
        let sampled = history.at(cx, 1.0).unwrap();
        assert_eq!(sampled.conditions()[0], ThermalBc::dirichlet(350.0).unwrap());
        assert_eq!(sampled.conditions()[1], ThermalBc::robin(20.0, 300.0).unwrap());
        assert_eq!(sampled.conditions()[2], ThermalBc::neumann(-100.0).unwrap());
        assert!(sampled.dirichlet().iter().all(|&(_, t)| t == 350.0));
        assert!(!sampled.dirichlet().is_empty());
        assert_eq!(sampled.region_names(), start.region_names());
        assert_eq!(sampled.adiabatic_remainder_faces(), start.adiabatic_remainder_faces());
        for slot in 0..mesh.boundary().len() { assert_eq!(sampled.region_for(slot), start.region_for(slot)); }
        assert_eq!(history.at(cx, 0.0).unwrap(), start);
        assert!(std::ptr::eq(history.mesh(), &mesh));
    });
}

#[test]
fn g3_shared_dirichlet_vertices_keep_declared_region_precedence() {
    let mesh = mesh();
    let env = Environment::earth_lab();
    let builder = |bottom_t| left(&mesh)
        .region("bottom", |f| on_box_face(f.centroid[1], 0.0), ThermalBc::dirichlet(bottom_t).unwrap())
        .unwrap();
    with_cx(|cx| {
        let history = ThermalBoundaryHistory::new(cx, builder(300.0), &env,
            vec![temperature_row("left"), ThermalRegionHistory {
                row: row("bottom", BcKind::Dirichlet, ramp(300.0, 500.0, TEMPERATURE_DIMS)),
                reference_temperature: None }], 0.0).unwrap();
        let expected = ThermalBoundaryBuilder::new(&mesh)
            .region("left", |f| on_box_face(f.centroid[0], 0.0), ThermalBc::dirichlet(350.0).unwrap()).unwrap()
            .region("bottom", |f| on_box_face(f.centroid[1], 0.0), ThermalBc::dirichlet(400.0).unwrap()).unwrap()
            .adiabatic_remainder().finish().unwrap();
        assert_eq!(history.at(cx, 1.0).unwrap(), expected);
    });
}

#[test]
fn g0_bad_names_duplicate_bindings_dimensions_and_future_coefficients_refuse() {
    let mesh = mesh();
    let env = Environment::earth_lab();
    with_cx(|cx| {
        for rows in [vec![temperature_row("absent")], vec![temperature_row("left"), temperature_row("left")]] {
            assert!(ThermalBoundaryHistory::new(cx, left(&mesh), &env, rows, 0.0).is_err());
        }
        let wrong = ThermalRegionHistory { row: row("left", BcKind::Dirichlet, ramp(1.0, 2.0, HEAT_FLUX_DIMS)),
            reference_temperature: None };
        assert!(ThermalBoundaryHistory::new(cx, left(&mesh), &env, vec![wrong], 0.0).is_err());
        let mut wrong = temperature_row("left");
        wrong.reference_temperature = Some(ramp(290.0, 300.0, TEMPERATURE_DIMS));
        assert!(ThermalBoundaryHistory::new(cx, left(&mesh), &env, vec![wrong], 0.0).is_err());
        let robin = || ThermalBoundaryBuilder::new(&mesh)
            .remainder("air", ThermalBc::robin(10.0, 293.15).unwrap()).unwrap();
        let wrong_reference = ThermalRegionHistory {
            row: row("air", BcKind::Robin, ramp(10.0, 20.0, HTC_DIMS)),
            reference_temperature: Some(ramp(1.0, 2.0, HEAT_FLUX_DIMS)) };
        assert!(ThermalBoundaryHistory::new(cx, robin(), &env, vec![wrong_reference], 0.0).is_err());
        let history = ThermalBoundaryHistory::new(cx, robin(), &env, vec![ThermalRegionHistory {
            row: row("air", BcKind::Robin, ramp(10.0, -10.0, HTC_DIMS)),
            reference_temperature: None }], 0.0).unwrap();
        let original = history.at(cx, 0.0).unwrap();
        assert!(history.at(cx, 2.0).is_err());
        assert_eq!(history.at(cx, 0.0).unwrap(), original);
    });
}

#[test]
fn g0_history_evaluation_honors_cancellation_without_changing_the_template() {
    let mesh = mesh();
    let env = Environment::earth_lab();
    let history = with_cx(|cx| ThermalBoundaryHistory::new(cx, left(&mesh), &env,
        vec![temperature_row("left")], 0.0).unwrap());
    with_cancelled_cx(|cx| assert!(matches!(history.at(cx, 1.0), Err(ConductionError::Cancelled { .. }))));
    with_cx(|cx| assert_eq!(history.at(cx, 1.0).unwrap().conditions()[0], ThermalBc::dirichlet(350.0).unwrap()));
}
