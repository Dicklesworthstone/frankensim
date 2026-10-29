use super::*;

#[test]
fn circular_terminal_and_its_regional_gas_reach_flue_feedback_and_microphone_pressure() {
    let input = include_str!("../../../../../examples/flue-dg-slit-chimney.performance")
        .replace("dg_terminal_slit 7 0.004 0.00008 0.001 16 0.1",
            "dg_terminal_circular 7 0.0001 0.002 16 0.1");
    let render = |text: &str| {
        let source = FlueDgFilePerformance::from_bytes(text.as_bytes(), 37,
            &CancelGate::new_clock_free()).unwrap();
        let load = source.renderer().voice().dg().unwrap().load();
        // No order of this three-opening fit is passive unchanged: order 4 fits
        // only to 1.9 %, and orders 8 and 12 fit to ~1e-6 with a dip at DC. The
        // selected fit carries a residue repair inside a tenth of the allowance.
        let record = load.radiation_record().unwrap();
        let report = record.report();
        let selected = report.rungs.iter().find(|r| r.order == report.selected_order).unwrap();
        let change = selected.repair_relative_change.unwrap_or_else(|| panic!("{selected:?}"));
        assert!(!selected.continuous_passive
            && change <= 0.1 * record.options().maximum_relative_error, "{selected:?}");
        assert!(report.continuous.passive && report.sampled.passive);
        let channel = load.terminal_channels()[0].channel().unwrap().circular().unwrap();
        assert_eq!(channel.spec().radius_m, 0.0001);
        assert_eq!(channel.spec().radial_cells, 16);
        assert_eq!(channel.spec().sample_period_s, 1.0 / 96000.0);
        assert_eq!(channel.gas(), load.branches()[0].duct().outlet_gas());
        assert_ne!(channel.gas(), load.duct().inlet_gas());
        let gas = *channel.gas();
        let mut renderer = source.into_renderer();
        let mut pressure = vec![0.0; 2048];
        for block in pressure.chunks_mut(37) { renderer.block(block).unwrap(); }
        let load = renderer.voice().dg().unwrap().load();
        let channel = &load.terminal_channels()[0];
        assert!(channel.cumulative_viscous_dissipation_j() > 0.0);
        assert!(channel.cumulative_numerical_dissipation_j() > 0.0);
        assert!(pressure.iter().all(|value| value.is_finite()));
        assert!(pressure.iter().any(|value| value.abs() > 1e-12));
        (pressure, load.state().clone(), gas)
    };
    let (reference, state, gas) = render(&input);
    let changed = input.replace("N2 0.6 CO2 0.4", "N2 0.9 CO2 0.1");
    let (pressure, changed_state, changed_gas) = render(&changed);
    assert_ne!(gas, changed_gas);
    assert_ne!(state, changed_state, "regional gas must affect the coupled source, not only output metadata");
    assert!(pressure.iter().zip(&reference).any(|(a, b)| (a - b).abs() > 1e-10));
}
