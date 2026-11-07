use super::*;

fn input(block: &str) -> String {
    baffled_tube_input(RADIATION_EXTERIOR)
        .replace("1e-11 128", "1e-11 256")
        .replace("\ndg_radiation", &format!(
            "\ndg_side_openings 1\ndg_side_chimney 7 1 0.01 0 0.0006 0.012 1 2\ndg_chimney_gas mixture 303.15 101325 mole 2 N2 0.6 CO2 0.4\n{block}\ndg_radiation",
        ))
}

const CHANNELS: &str = "dg_terminal_channels 2 16\ndg_terminal_closed 1\ndg_terminal_slit 7 0.004 0.00008 0.001 16 0.1";

#[test]
fn declared_closure_and_slit_keep_source_port_order_and_original_gas_records() {
    let (source, records) = read_receipted(&input(CHANNELS), None).unwrap();
    assert_eq!(source.terminal_channels, vec![
        DgTerminalChannelInput { port: PortId(1), kind: DgTerminalChannelKind::Closed },
        DgTerminalChannelInput { port: PortId(7), kind: DgTerminalChannelKind::Slit {
            width_m: 0.004, height_m: 0.00008, length_m: 0.001,
            transverse_cells: 16, maximum_mach: 0.1,
        } },
    ]);
    assert_eq!(source.maximum_terminal_transverse_cells, 16);
    assert_eq!(source.side_openings[0].0.id, PortId(7));
    assert_eq!(source.chimneys[0].port, PortId(7));
    assert_eq!(source.chimneys[0].gas, records[0]);
    assert_eq!(source.radiation.as_ref().unwrap().1, records[1]);
    assert_ne!(records[0].state(), records[1].state());
}

#[test]
fn closures_need_no_transverse_states_and_omitted_terminals_remain_direct() {
    let source = read(&input("dg_terminal_channels 1 0\ndg_terminal_closed 7")).unwrap();
    assert_eq!(source.maximum_terminal_transverse_cells, 0);
    assert_eq!(source.terminal_channels.len(), 1);
    assert_eq!(source.terminal_channels[0].port, PortId(7));
    assert_eq!(source.side_openings.len() + 1, 2);
    assert!(read(TUBE).unwrap().terminal_channels.is_empty());
    assert!(read(&baffled_tube_input(RADIATION_EXTERIOR)).unwrap().terminal_channels.is_empty());
}

#[test]
fn circular_terminal_is_geometry_data_with_the_original_regional_gas_and_cell_budget() {
    let block = "dg_terminal_channels 2 24\ndg_terminal_slit 1 0.004 0.00008 0.001 8 0.1\ndg_terminal_circular 7 0.0001 0.002 16 0.1";
    let (source, records) = read_receipted(&input(block), None).unwrap();
    assert_eq!(source.maximum_terminal_transverse_cells, 24);
    assert_eq!(source.terminal_channels[1], DgTerminalChannelInput {
        port: PortId(7), kind: DgTerminalChannelKind::Circular {
            radius_m: 0.0001, length_m: 0.002, radial_cells: 16, maximum_mach: 0.1,
        },
    });
    assert_eq!(source.chimneys[0].gas, records[0]);
    for (from, to) in [
        ("channels 2 24", "channels 2 23"),
        ("0.0001 0.002 16 0.1", "0 0.002 16 0.1"),
        ("0.0001 0.002 16 0.1", "0.0003 0.002 16 0.1"),
        ("0.0001 0.002 16 0.1", "0.0001 0.002 1 0.1"),
        ("0.0001 0.002 16 0.1", "0.0001 0.002 16 NaN"),
        ("0.0001 0.002 16 0.1", "0.0001 0.002 16 0.31"),
        ("0.0001 0.002 16 0.1", "0.0001 0.002 16 0.1 ignored"),
    ] {
        assert!(read(&input(&block.replace(from, to))).is_err(), "admitted {to}");
    }
}

#[test]
fn annular_terminal_reads_both_radii_and_preserves_regional_gas_and_cell_budget() {
    let block = "dg_terminal_channels 2 24\ndg_terminal_slit 1 0.004 0.00008 0.001 8 0.1\ndg_terminal_annular 7 0.00004 0.0001 0.002 16 0.1";
    let (source, records) = read_receipted(&input(block), None).unwrap();
    assert_eq!(source.maximum_terminal_transverse_cells, 24);
    assert_eq!(source.terminal_channels[1], DgTerminalChannelInput {
        port: PortId(7), kind: DgTerminalChannelKind::Annular {
            inner_radius_m: 0.00004, outer_radius_m: 0.0001, length_m: 0.002,
            radial_cells: 16, maximum_mach: 0.1,
        },
    });
    assert_eq!(source.chimneys[0].gas, records[0]);
    for (from, to) in [
        ("channels 2 24", "channels 2 23"),
        ("0.00004 0.0001", "0 0.0001"),
        ("0.00004 0.0001", "0.0001 0.0001"),
        ("0.00004 0.0001", "0.0002 0.0001"),
        ("0.00004 0.0001", "0.00004 0.0004"),
        ("0.00004 0.0001", "NaN 0.0001"),
        ("0.002 16 0.1", "0.002 1 0.1"),
        ("0.002 16 0.1", "0.002 4097 0.1"),
        ("0.002 16 0.1", "0.002 16 0.31"),
        ("0.002 16 0.1", "0.002 16 0.1 ignored"),
    ] {
        assert!(read(&input(&block.replace(from, to))).is_err(), "admitted {to}");
    }
}

#[test]
fn terminal_source_refuses_undeclared_ports_duplicate_losses_and_unbounded_states() {
    let valid = input(CHANNELS);
    for (from, to) in [
        ("dg_terminal_channels 2 16", "dg_terminal_channels 0 16"),
        ("dg_terminal_channels 2 16", "dg_terminal_channels 3 16"),
        ("dg_terminal_channels 2 16", "dg_terminal_channels 2 15"),
        ("dg_terminal_channels 2 16", "dg_terminal_channels 2 4097"),
        ("dg_terminal_closed 1", "dg_terminal_closed 0"),
        ("dg_terminal_closed 1", "dg_terminal_closed 7"),
        ("dg_terminal_slit 7", "dg_terminal_slit 99"),
        ("dg_terminal_slit 7", "dg_terminal_unknown 7"),
        ("0.004 0.00008", "0 0.00008"),
        ("0.004 0.00008", "0.004 -0.00008"),
        ("0.001 16 0.1", "0.001 0 0.1"),
        ("0.001 16 0.1", "0.001 1 0.1"),
        ("0.004 0.00008", "0.0001 0.00008"),
        ("0.001 16 0.1", "0.001 16 0"),
        ("0.001 16 0.1", "0.001 16 0.31"),
        ("0.001 16 0.1", "0.001 16 NaN"),
        ("0.001 16 0.1", "0.001 16 0.1 ignored"),
    ] {
        assert!(read_receipted(&valid.replace(from, to), None).is_err(), "admitted {to}");
    }
}

#[test]
fn an_orifice_terminal_is_a_chimney_row_with_the_cards_own_vena_contracta_and_no_transverse_state() {
    let block = "dg_terminal_channels 2 0\ndg_terminal_closed 1\ndg_terminal_orifice 7 0.61";
    let source = read(&input(block)).unwrap();
    assert_eq!(source.maximum_terminal_transverse_cells, 0);
    assert_eq!(source.terminal_channels[1], DgTerminalChannelInput {
        port: PortId(7), kind: DgTerminalChannelKind::Orifice { vena_contracta: 0.61 },
    });
    // Another card's coefficient is just another number on the row.
    let other = read(&input("dg_terminal_channels 1 0\ndg_terminal_orifice 7 1.0956")).unwrap();
    assert_eq!(other.terminal_channels[0].kind, DgTerminalChannelKind::Orifice { vena_contracta: 1.0956 });
    for (from, to) in [
        ("dg_terminal_orifice 7 0.61", "dg_terminal_orifice 1 0.61"),
        ("dg_terminal_orifice 7 0.61", "dg_terminal_orifice 99 0.61"),
        ("dg_terminal_orifice 7 0.61", "dg_terminal_orifice 7 0"),
        ("dg_terminal_orifice 7 0.61", "dg_terminal_orifice 7 -0.61"),
        ("dg_terminal_orifice 7 0.61", "dg_terminal_orifice 7 NaN"),
        ("dg_terminal_orifice 7 0.61", "dg_terminal_orifice 7 2.5"),
        ("dg_terminal_orifice 7 0.61", "dg_terminal_orifice 7"),
        ("dg_terminal_orifice 7 0.61", "dg_terminal_orifice 7 0.61 ignored"),
        ("dg_terminal_closed 1", "dg_terminal_closed 7"),
    ] {
        assert!(read(&input(&block.replace(from, to))).is_err(), "admitted {to}");
    }
}
