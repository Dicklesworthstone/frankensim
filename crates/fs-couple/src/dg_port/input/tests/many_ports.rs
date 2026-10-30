use super::*;

#[test]
fn native_dg_retains_sixteen_chimneys_and_closed_holes_in_one_source() {
    let mut text = String::from(
        "dg_tube 0.34 0.001 baffled 17 1 8 50 2000 0.01 1e-11 2048\ndg_side_openings 16\n",
    );
    for i in 0..16 {
        text.push_str(&format!(
            "dg_side_chimney {} {} {} 0 0.001 0.004 1 1\ndg_chimney_gas 293.15 101325 0\n",
            i + 2,
            i + 1,
            0.02 * (i + 1) as f64
        ));
    }
    text.push_str("dg_terminal_channels 2 0\ndg_terminal_closed 2\ndg_terminal_closed 17\n");
    text.push_str(&format!(
        "dg_radiation computed-baffled {}\ndg_exterior 293.15 101325 0",
        COMPUTED_CONTROLS.replace("1e-8 64", "1e-8 204")
    ));
    let (source, records) = read_receipted(&text, None).unwrap();
    assert_eq!(source.side_openings.len(), 16);
    assert_eq!(source.chimneys.len(), 16);
    assert_eq!(records.len(), 17);
    assert_eq!(source.terminal_channels.len(), 2);
    assert_eq!(
        source.side_openings[15].0,
        DgSidePort {
            id: PortId(17),
            downstream_element: 16
        }
    );
    let branches = source
        .prepare_chimneys(&CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(branches.len(), 16);
    for (i, branch) in branches.iter().enumerate() {
        assert_eq!(branch.port, PortId((i + 2) as u32));
        assert_eq!(branch.duct.total_length(), 0.004);
        assert!(branch.duct.has_viscothermal_memory());
        assert_eq!(branch.duct.inlet_gas(), records[i].state());
        assert_eq!(
            branch.duct.outlet(),
            if i == 0 || i == 15 {
                DgBoundary::Closed
            } else {
                DgBoundary::Port(PortId(1))
            }
        );
    }
    assert_eq!(
        branches
            .iter()
            .filter(|branch| branch.duct.outlet() == DgBoundary::Closed)
            .count(),
        2
    );
    // The larger count cannot bypass original topology, identity or DOF checks.
    for invalid in [
        text.replace("dg_side_openings 16", "dg_side_openings 33"),
        text.replace("dg_side_chimney 17 16", "dg_side_chimney 2 16"),
        text.replace("dg_side_chimney 17 16", "dg_side_chimney 17 15"),
        text.replace("1e-11 2048", "1e-11 1024"),
    ] {
        assert!(read_receipted(&invalid, None).is_err());
    }
}
