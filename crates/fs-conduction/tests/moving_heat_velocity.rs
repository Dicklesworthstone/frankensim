use fs_conduction::bc::history::march::moving::{MovingHeatBudget, MovingHeatSource};
use fs_conduction::bc::history::march::moving::motion::velocity::MovingHeatPathBudget;
use fs_conduction::COMPONENT_POWER_DIMS;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_qty::{Dims, QtyAny};
use fs_scenario::signal::{Interp, TimeSignal};

const LENGTH: Dims = Dims([1,0,0,0,0,0]);
const BUDGET: MovingHeatPathBudget = MovingHeatPathBudget { max_knots: 16, max_segments: 16 };
fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| f(&gate,
        &Cx::new(&gate, arena, StreamKey { seed: 31, kernel_id: 51, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic)));
}
fn constant(value: f64) -> TimeSignal { TimeSignal::Constant(QtyAny::new(value, LENGTH)) }
fn table(times: &[f64], values: &[f64]) -> TimeSignal {
    TimeSignal::Table { times: times.to_vec(), values: values.to_vec(), dims: LENGTH, interp: Interp::Linear }
}
fn source(path: [TimeSignal; 3]) -> MovingHeatSource {
    MovingHeatSource::new(path, [1.0; 3],
        TimeSignal::Constant(QtyAny::new(0.0, COMPONENT_POWER_DIMS)), 0.0,
        MovingHeatBudget { max_vertices: 8, max_elements: 6, min_active_vertices: 1,
            power_relative_tolerance: 1e-11 }).unwrap()
}

#[test]
fn g1_union_clock_keeps_diagonal_and_returning_motion_while_power_is_off() {
    let path = source([table(&[0.,1.,2.], &[-1.,1.,-1.]),
        table(&[0.5,1.5], &[0.,3.]), constant(-7.)]);
    with_cx(|_, cx| {
        let segments = path.path_velocity_segments(cx, BUDGET).unwrap();
        assert_eq!(segments.len(), 4);
        for (i, expected) in [[2.,0.,0.],[2.,3.,0.],[-2.,3.,0.],[-2.,0.,0.]].into_iter().enumerate() {
            assert_eq!(segments[i].velocity_m_s, expected);
            assert_eq!(segments[i].end_time_s - segments[i].start_time_s, 0.5);
            assert_eq!(segments[i].axes[0].unwrap().duration_s, 1.0);
        }
        assert_eq!(segments[1].axes[0].unwrap().value_indices, [0,1]);
        assert_eq!(segments[2].axes[0].unwrap().value_indices, [1,2]);
        assert_eq!(segments[1].axes[1].unwrap().value_weights_s_inv, [-1.,1.]);
        assert!(segments.iter().all(|segment| segment.axes[2].is_none()));
    });
}

#[test]
fn g3_authored_velocity_derivatives_survive_foreign_axis_breakpoints() {
    let times = [&[0.,2.,3.][..], &[0.5,1.5][..], &[1.,4.][..]];
    let values = [vec![-1.,2.,1.], vec![0.,-2.], vec![1.,3.]];
    let make = |v: &[Vec<f64>; 3]| source(std::array::from_fn(|a| table(times[a], &v[a])));
    with_cx(|_, cx| {
        let original = make(&values).path_velocity_segments(cx, BUDGET).unwrap();
        for axis in 0..3 {
            for knot in 0..values[axis].len() {
                let mut plus = values.clone(); let mut minus = values.clone();
                plus[axis][knot] += 1e-5; minus[axis][knot] -= 1e-5;
                let p = make(&plus).path_velocity_segments(cx, BUDGET).unwrap();
                let m = make(&minus).path_velocity_segments(cx, BUDGET).unwrap();
                for (i, segment) in original.iter().enumerate() {
                    let expected = segment.axes[axis].map_or(0., |rate|
                        rate.value_indices.iter().zip(rate.value_weights_s_inv)
                            .filter(|(index, _)| **index == knot).map(|(_, weight)| weight).sum::<f64>());
                    let fd = (p[i].velocity_m_s[axis] - m[i].velocity_m_s[axis]) / 2e-5;
                    assert!((fd - expected).abs() < 1e-9);
                }
            }
        }
    });
}

#[test]
fn g0_narrow_clock_uses_the_right_segment_without_a_rounded_midpoint() {
    let end = 1.0_f64.next_up();
    let path = source([table(&[1.,end], &[0.,1.]), table(&[0.,2.], &[0.,0.]), constant(0.)]);
    with_cx(|_, cx| {
        let segments = path.path_velocity_segments(cx, BUDGET).unwrap();
        assert_eq!(segments.len(), 3);
        assert_eq!(segments[1].start_time_s, 1.0);
        assert_eq!(segments[1].end_time_s, end);
        assert_eq!(segments[1].velocity_m_s[0], 1.0 / (end - 1.0));
        assert_eq!(segments[0].velocity_m_s[0], 0.0);
        assert_eq!(segments[2].velocity_m_s[0], 0.0);
    });
}

#[test]
fn g0_g4_path_admission_refuses_jumps_budgets_and_cancellation() {
    with_cx(|gate, cx| {
        let stationary = source([constant(0.), constant(1.), constant(2.)]);
        assert!(stationary.path_velocity_segments(cx, MovingHeatPathBudget { max_knots:0, max_segments:0 }).unwrap().is_empty());
        let moving = source([table(&[0.,1.,2.], &[0.,1.,0.]), constant(0.), constant(0.)]);
        assert!(moving.path_velocity_segments(cx, MovingHeatPathBudget { max_knots:2, ..BUDGET }).is_err());
        assert!(moving.path_velocity_segments(cx, MovingHeatPathBudget { max_segments:1, ..BUDGET }).is_err());
        let held = source([TimeSignal::Table { times:vec![0.,1.], values:vec![0.,1.],
            dims:LENGTH, interp:Interp::Hold }, constant(0.), constant(0.)]);
        assert!(held.path_velocity_segments(cx, BUDGET).is_err());
        let ramp = source([TimeSignal::Ramp { t_start:0., t_end:2., from:QtyAny::new(-2., LENGTH),
            to:QtyAny::new(4., LENGTH) }, constant(0.), constant(0.)]);
        assert_eq!(ramp.path_velocity_segments(cx, BUDGET).unwrap()[0].velocity_m_s, [3.,0.,0.]);
        gate.request();
        assert!(moving.path_velocity_segments(cx, BUDGET).is_err());
    });
}
