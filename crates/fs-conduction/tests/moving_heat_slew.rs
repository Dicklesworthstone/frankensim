use fs_conduction::bc::history::march::moving::{MovingHeatBudget, MovingHeatSource};
use fs_conduction::COMPONENT_POWER_DIMS;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_qty::{Dims, QtyAny};
use fs_scenario::signal::{Interp, TimeSignal};

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&gate, &Cx::new(&gate, arena,
            StreamKey { seed: 41, kernel_id: 61, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic));
    });
}
fn source(power: TimeSignal) -> MovingHeatSource {
    let length = Dims([1, 0, 0, 0, 0, 0]);
    MovingHeatSource::new(std::array::from_fn(|_| TimeSignal::Constant(QtyAny::new(0.5, length))),
        [1.0; 3], power, 0.0, MovingHeatBudget { max_vertices: 8, max_elements: 6,
            min_active_vertices: 1, power_relative_tolerance: 1e-11 }).unwrap()
}
fn table(values: Vec<f64>, interp: Interp) -> MovingHeatSource {
    source(TimeSignal::Table { times: vec![0.0, 2.0, 2.25, 4.0], values,
        dims: COMPONENT_POWER_DIMS, interp })
}

#[test]
fn g1_rates_cover_short_segments_and_fixed_shutdown_without_time_sampling() {
    with_cx(|_, cx| {
        let s = table(vec![0.0, 10.0, 13.0, 0.0], Interp::Linear);
        let rows = s.absorbed_power_slew_segments(cx, 3).unwrap();
        assert_eq!(rows.len(), 3);
        assert_eq!(rows[0].rate_w_s, 5.0);
        assert_eq!(rows[1].rate_w_s, 12.0);
        assert_eq!(rows[2].rate_w_s, -13.0 / 1.75);
        assert_eq!(rows[2].power_value_indices, [2, 3]);
        assert_eq!(rows[2].endpoint_powers_w[1], 0.0);
        assert!(s.absorbed_power_slew_segments(cx, 2).is_err());
        let constant = source(TimeSignal::Constant(QtyAny::new(7.0, COMPONENT_POWER_DIMS)));
        assert!(constant.absorbed_power_slew_segments(cx, 0).unwrap().is_empty());
        let ramp = source(TimeSignal::Ramp { t_start: 1.0, t_end: 3.0,
            from: QtyAny::new(3.0, COMPONENT_POWER_DIMS), to: QtyAny::new(9.0, COMPONENT_POWER_DIMS) });
        let row = ramp.absorbed_power_slew_segments(cx, 1).unwrap()[0];
        assert_eq!(row.rate_w_s, 3.0);
        assert_eq!(row.power_value_weights_s_inv, [-0.5, 0.5]);
    });
}

#[test]
fn g3_power_derivatives_match_perturbations_on_each_original_segment() {
    with_cx(|_, cx| {
        let values = vec![1.0, 10.0, 13.0, 2.0];
        let rows = table(values.clone(), Interp::Linear).absorbed_power_slew_segments(cx, 3).unwrap();
        for (segment, row) in rows.iter().enumerate() {
            for side in 0..2 {
                let index = row.power_value_indices[side];
                let (mut plus, mut minus) = (values.clone(), values.clone());
                plus[index] += 1e-5; minus[index] -= 1e-5;
                let p = table(plus, Interp::Linear).absorbed_power_slew_segments(cx, 3).unwrap();
                let m = table(minus, Interp::Linear).absorbed_power_slew_segments(cx, 3).unwrap();
                let fd = (p[segment].rate_w_s - m[segment].rate_w_s) / 2e-5;
                assert!((fd - row.power_value_weights_s_inv[side]).abs() < 1e-8);
            }
        }
    });
}

#[test]
fn g0_held_jumps_are_not_misrepresented_as_finite_rates_and_cancel_refuses() {
    with_cx(|gate, cx| {
        let held = table(vec![0.0, 10.0, 13.0, 0.0], Interp::Hold);
        assert!(held.absorbed_power_slew_segments(cx, 100).is_err());
        let linear = table(vec![0.0, 10.0, 13.0, 0.0], Interp::Linear);
        gate.request();
        assert!(linear.absorbed_power_slew_segments(cx, 3).is_err());
    });
}
