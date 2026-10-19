//! Moving-source power and actual P1-load regressions.
mod support;
use fs_conduction::bc::history::march::moving::{MovingHeatBudget, MovingHeatSource};
use fs_conduction::bc::history::march::moving::motion::ScanMotionLimit;
use fs_conduction::{ConductionError, ConductionMesh, ScalarField, COMPONENT_POWER_DIMS};
use fs_conduction::fixtures::box_grid;
use fs_qty::{Dims, QtyAny};
use fs_scenario::signal::{Interp, TimeSignal};
use support::{with_cancelled_cx, with_cx};

const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
fn constant(value: f64, dims: Dims) -> TimeSignal { TimeSignal::Constant(QtyAny::new(value, dims)) }
fn budget() -> MovingHeatBudget {
    MovingHeatBudget { max_vertices: 10000, max_elements: 10000, min_active_vertices: 1,
        power_relative_tolerance: 1e-12 }
}
fn grid(n: usize) -> ConductionMesh {
    let (complex, positions) = box_grid([n; 3], [1.0; 3]);
    ConductionMesh::new(complex, positions).unwrap()
}
fn source(power: TimeSignal) -> MovingHeatSource {
    MovingHeatSource::new([TimeSignal::Ramp { t_start: 0.0, t_end: 1.0,
        from: QtyAny::new(0.25, LENGTH), to: QtyAny::new(0.75, LENGTH) },
        constant(0.5, LENGTH), constant(0.5, LENGTH)], [0.48; 3], power, 0.0, budget()).unwrap()
}
// Independent consistent P1 source load, not the implementation's row-sum loop.
fn integrated_load(mesh: &ConductionMesh, field: &ScalarField) -> f64 {
    let mut assembled = vec![0.0; mesh.vertex_count()];
    for (e, tet) in mesh.complex().tets.iter().enumerate() {
        for &a in tet {
            for &b in tet {
                assembled[a as usize] += mesh.element_volume(e) / 20.0
                    * if a == b { 2.0 * field.at(b as usize) } else { field.at(b as usize) };
            }
        }
    }
    assembled.iter().sum()
}

#[test]
fn g1_motion_limits_follow_anisotropic_footprints_and_clamped_path_segments() {
    let ramp = |end| TimeSignal::Ramp { t_start: 0.0, t_end: 1.0,
        from: QtyAny::new(0.0, LENGTH), to: QtyAny::new(end, LENGTH) };
    let moving = MovingHeatSource::new([ramp(2.0), ramp(3.0), constant(0.0, LENGTH)],
        [0.4, 0.6, 1.0], constant(1.0, COMPONENT_POWER_DIMS), 0.0, budget()).unwrap();
    let policy = ScanMotionLimit { max_normalized_travel: 0.2, max_chebyshev_coefficients: 0 };
    with_cx(|cx| {
        let end = moving.motion_limited_endpoint(cx, 0.0, 1.0, policy).unwrap();
        assert!((end-0.02).abs() < 1e-14);
        assert!(2.0*end/0.4 + 3.0*end/0.6 <= 0.2);
        assert_eq!(moving.motion_limited_endpoint(cx, -1.0, 2.0, policy).unwrap(), 0.0);
        assert_eq!(moving.motion_limited_endpoint(cx, 1.0, 2.0, policy).unwrap(), 2.0);
    });
    with_cancelled_cx(|cx| assert!(moving.motion_limited_endpoint(cx, 0.0, 1.0, policy).is_err()));
}

#[test]
fn g0_motion_limits_distinguish_held_jumps_and_refuse_unrepresentable_steps() {
    let path = |interp| TimeSignal::Table { times: vec![0.0, 0.125, 0.25],
        values: vec![0.0, 10.0, 0.0], dims: LENGTH, interp };
    let moving = |interp| MovingHeatSource::new([path(interp), constant(0.0, LENGTH), constant(0.0, LENGTH)],
        [0.5; 3], constant(0.0, COMPONENT_POWER_DIMS), 0.0, budget()).unwrap();
    let policy = ScanMotionLimit { max_normalized_travel: 0.25, max_chebyshev_coefficients: 0 };
    with_cx(|cx| {
        let linear = moving(Interp::Linear);
        let end = linear.motion_limited_endpoint(cx, 0.0, 1.0, policy).unwrap();
        assert!((end-0.0015625).abs() < 1e-15);
        let held = moving(Interp::Hold);
        assert_eq!(held.motion_limited_endpoint(cx, 0.0, 1.0, policy).unwrap(), 0.125);
        assert_eq!(held.motion_limited_endpoint(cx, 0.125, 1.0, policy).unwrap(), 0.25);
        assert_eq!(held.motion_limited_endpoint(cx, 0.25, 1.0, policy).unwrap(), 1.0);
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY, 1.01] {
            assert!(linear.motion_limited_endpoint(cx, 0.0, 1.0,
                ScanMotionLimit { max_normalized_travel: bad, ..policy }).is_err());
        }
        assert!(linear.motion_limited_endpoint(cx, 0.01, 0.1,
            ScanMotionLimit { max_normalized_travel: f64::MIN_POSITIVE, ..policy }).is_err());
    });
}

#[test]
fn g1_moving_footprint_preserves_watts_through_the_actual_p1_source_load() {
    with_cx(|cx| {
        let heater = source(constant(12.0, COMPONENT_POWER_DIMS));
        for n in [2, 4, 8] {
            let mesh = grid(n);
            for time in [-1.0, 0.0, 0.25, 0.5, 0.75, 1.0, 2.0] {
                let sample = heater.at(cx, &mesh, time).unwrap();
                assert!((integrated_load(&mesh, &sample.source) - 12.0).abs() < 1e-10);
                assert!((sample.delivered_power_w - 12.0).abs() < 1e-10);
                assert!(sample.active_vertices >= 1);
                assert!((0..mesh.vertex_count()).all(|i| sample.source.at(i) >= 0.0));
            }
        }
        let mesh = grid(4);
        let peak_x = |time| {
            let sampled = heater.at(cx, &mesh, time).unwrap();
            let i = (0..mesh.vertex_count()).max_by(|&i, &j|
                sampled.source.at(i).total_cmp(&sampled.source.at(j))).unwrap();
            mesh.positions()[i][0]
        };
        assert_eq!(peak_x(0.0), 0.25);
        assert_eq!(peak_x(1.0), 0.75);
    });
}

#[test]
fn g3_amplitude_scaling_and_exact_left_limits_preserve_the_declared_scan() {
    let mesh = grid(4);
    with_cx(|cx| {
        let a = source(constant(10.0, COMPONENT_POWER_DIMS)).at(cx, &mesh, 0.3).unwrap();
        let b = source(constant(30.0, COMPONENT_POWER_DIMS)).at(cx, &mesh, 0.3).unwrap();
        for i in 0..mesh.vertex_count() {
            assert!((b.source.at(i) - 3.0 * a.source.at(i)).abs() < 1e-10);
        }
        let heater = source(TimeSignal::Table { times: vec![0.0, 0.125, 0.25],
            values: vec![0.0, 120.0, 0.0], dims: COMPONENT_POWER_DIMS, interp: Interp::Hold });
        assert_eq!(heater.at(cx, &mesh, 0.25).unwrap().delivered_power_w, 0.0);
        assert!((heater.before(cx, &mesh, 0.25).unwrap().delivered_power_w - 120.0).abs() < 1e-10);
        assert_eq!(heater.at(cx, &mesh, 0.3).unwrap().delivered_power_w, 0.0);
    });
}

#[test]
fn g0_bad_units_limits_and_unresolved_nonzero_footprints_are_not_silently_repaired() {
    let center = [constant(0.0, LENGTH), constant(0.0, LENGTH), constant(0.0, LENGTH)];
    for radius in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(MovingHeatSource::new(center.clone(), [radius; 3], constant(1.0, COMPONENT_POWER_DIMS), 0.0, budget()).is_err());
    }
    assert!(MovingHeatSource::new(center.clone(), [0.5; 3], constant(1.0, LENGTH), 0.0, budget()).is_err());
    assert!(MovingHeatSource::new(center, [0.5; 3], constant(-1.0, COMPONENT_POWER_DIMS), 0.0, budget()).is_err());
    let mesh = grid(4);
    let outside = [constant(50.0, LENGTH), constant(50.0, LENGTH), constant(50.0, LENGTH)];
    with_cx(|cx| {
        for watts in [0.0, 1.0] {
            let heater = MovingHeatSource::new(outside.clone(), [0.5; 3],
                constant(watts, COMPONENT_POWER_DIMS), 0.0, budget()).unwrap();
            assert_eq!(heater.at(cx, &mesh, 0.0).is_ok(), watts == 0.0);
        }
        let heater = MovingHeatSource::new(std::array::from_fn(|_| constant(0.5, LENGTH)), [0.5; 3],
            constant(1.0, COMPONENT_POWER_DIMS), 0.0,
            MovingHeatBudget { max_vertices: 1, ..budget() }).unwrap();
        assert!(heater.at(cx, &mesh, 0.0).is_err());
        assert!(source(constant(1.0, COMPONENT_POWER_DIMS)).at(cx, &mesh, f64::NAN).is_err());
    });
}

#[test]
fn g4_cancellation_preserves_reusable_immutable_source() {
    let mesh = grid(4);
    let heater = source(constant(10.0, COMPONENT_POWER_DIMS));
    let before = with_cx(|cx| heater.at(cx, &mesh, 0.5).unwrap());
    with_cancelled_cx(|cx| assert!(matches!(heater.at(cx, &mesh, 0.5), Err(ConductionError::Cancelled { .. }))));
    with_cx(|cx| assert_eq!(heater.at(cx, &mesh, 0.5).unwrap(), before));
}
