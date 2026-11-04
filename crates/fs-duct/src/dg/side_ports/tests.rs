use super::*;
use fs_material::gas::GasSpec;

const SIDE: DgSidePort = DgSidePort { id: PortId(7), downstream_element: 1 };

fn bore(inlet: DgBoundary, outlet: DgBoundary) -> DgDuct {
    let geometry = Duct {
        segments: vec![Segment::Cylinder { radius: 0.001, length: 0.05 }],
    };
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    DgDuct::new(&geometry, &gas, DgMesh { elements: 2, order: 2 }, inlet, outlet).unwrap()
}

fn walls(duct: DgDuct, enabled: bool) -> DgDuct {
    if enabled { duct.with_viscothermal(8, [50.0, 2000.0], 0.01).unwrap() } else { duct }
}

#[test]
fn side_port_preserves_original_operator_and_binds_its_trace_to_state_identity() {
    let original = bore(DgBoundary::Port(PortId(1)), DgBoundary::Port(PortId(2)));
    let empty = original.clone().with_side_ports(&[]).unwrap();
    assert_eq!(original.operator_identity(), empty.operator_identity());
    let duct = original.clone().with_side_ports(&[SIDE]).unwrap();
    assert_eq!(duct.l_tilde, original.l_tilde);
    assert_eq!(duct.g_block, original.g_block);
    assert_eq!(duct.node_positions(), original.node_positions());
    assert_eq!(duct.total_dofs(), original.total_dofs());
    assert_eq!(&duct.ports[..2], original.ports.as_slice());
    assert_eq!(duct.side_ports(), &[SIDE]);
    assert_ne!(duct.operator_identity(), original.operator_identity());
    let source = duct.side_port_source(SIDE.id).unwrap();
    assert_eq!(source.port, SIDE);
    assert_eq!(source.station_m, 0.025);
    assert_eq!(source.bore_area_m2, core::f64::consts::PI * 0.001 * 0.001);
    assert_eq!(source.gas, *original.gas());

    let stepper = duct.propagator(1e-5, 1e-11).unwrap();
    let response = stepper.multiport_response(&duct.zero_state()).unwrap();
    let scale = response.impedance_pa_s_m3.iter().map(|value| value.abs()).fold(0.0, f64::max);
    for row in 0..3 {
        for col in 0..3 {
            assert!((response.impedance_pa_s_m3[3 * row + col]
                - response.impedance_pa_s_m3[3 * col + row]).abs() < 1e-12 * scale);
        }
    }
    let mut foreign = original.zero_state();
    let before = foreign.clone();
    assert!(stepper.step(&mut foreign, &[]).is_err());
    assert_eq!(foreign, before, "changed port basis cannot reinterpret an accepted history");
    for invalid in [
        DgSidePort { id: PortId(1), downstream_element: 1 },
        DgSidePort { id: PortId(8), downstream_element: 0 },
        DgSidePort { id: PortId(8), downstream_element: 2 },
    ] {
        assert!(original.clone().with_side_ports(&[invalid]).is_err());
    }
    assert!(original.clone().with_side_ports(&[SIDE, SIDE]).is_err());
    assert!(duct.clone().with_side_ports(&[
        DgSidePort { id: PortId(8), downstream_element: 1 },
    ]).is_err());
    assert!(original.with_side_ports(&[SIDE; MAX_DG_SIDE_PORTS + 1]).is_err());
    assert!(duct.side_port_source(PortId(1)).is_err());
}

#[test]
fn point_inflow_is_conjugate_to_the_actual_trace_and_retains_wall_memory() {
    for wall_memory in [false, true] {
        let original = bore(DgBoundary::Closed, DgBoundary::Closed);
        let before_wall = walls(original.clone().with_side_ports(&[SIDE]).unwrap(), wall_memory);
        let after_wall = walls(original.with_side_ports(&[]).unwrap(), wall_memory)
            .with_side_ports(&[SIDE]).unwrap();
        assert_eq!(before_wall.operator_identity(), after_wall.operator_identity());
        assert_eq!(before_wall.l_tilde, after_wall.l_tilde);
        assert_eq!(before_wall.ports, after_wall.ports);
        let duct = before_wall;
        let dt = 1e-5;
        let stepper = duct.propagator(dt, 1e-11).unwrap();
        let mut state = duct.project_initial_state(|x| 0.3 * det::cos(30.0 * x), |_| 0.0);
        let integral = duct.to_energy_coords(
            &vec![1.0; duct.field_dofs()], &vec![0.0; duct.field_dofs()],
        );
        let volume = |state: &DgState| integral.iter().zip(&state.tilde)
            .map(|(weight, value)| weight * value).sum::<f64>();
        for flow in [1e-8, -2e-8, 0.0, 3e-8] {
            let old_volume = volume(&state);
            let history = state.clone();
            let relation = stepper.multiport_response(&state).unwrap();
            assert_eq!(state, history);
            assert_eq!(relation.ports, vec![SIDE.id]);
            let expected = relation.history_pressure_pa[0]
                + relation.impedance_pa_s_m3[0] * flow;
            let scale = relation.history_pressure_pa[0].abs()
                + (relation.impedance_pa_s_m3[0] * flow).abs();
            let output = stepper.step(&mut state, &[PortInput { id: SIDE.id, effort: flow }]).unwrap();
            stepper.sync_nodal(&mut state).unwrap();
            let port = output.port_outputs[0];
            assert!((port.conjugate_rate - expected).abs() < 1e-11 * scale.max(1e-30));
            let face_pressure = state.pressure()[SIDE.downstream_element * (duct.mesh().order + 1)];
            assert!((port.end_conjugate - face_pressure).abs()
                < 1e-12 * (port.end_conjugate.abs() + face_pressure.abs()).max(1e-30));
            assert_eq!(output.ledger.port_work_step_j, flow * port.delta_x);
            if !wall_memory {
                // Integrate the original discrete continuity equation. A point
                // side source changes total compressive volume by exactly Q dt.
                let target = old_volume + flow * dt;
                assert!((volume(&state) - target).abs()
                    < 1e-9 * (old_volume.abs() + (flow * dt).abs()).max(1e-30));
                assert!(output.ledger.relative_energy_drift < 1e-9);
            }
        }
        if wall_memory { assert!(state.wall_stored_energy_j() > 0.0); }
    }
}

#[test]
fn original_dg_chimney_and_side_port_exchange_equal_pressure_and_opposite_work() {
    // This direct composition is an oracle for a runtime load owner: neither
    // original propagator is replaced, delayed, or collapsed to a reflection.
    for wall_memory in [false, true] {
        let original = bore(DgBoundary::Port(PortId(1)), DgBoundary::Closed);
        let main = walls(original.clone().with_side_ports(&[SIDE]).unwrap(), wall_memory);
        let unbranched = walls(original, wall_memory);
        let geometry = Duct {
            segments: vec![Segment::Cylinder { radius: 0.001, length: 0.03 }],
        };
        let branch = walls(DgDuct::new(&geometry, main.gas(), DgMesh { elements: 1, order: 2 },
            DgBoundary::Port(PortId(9)), DgBoundary::Closed).unwrap(), wall_memory);
        let dt = 1e-5;
        let a = main.propagator(dt, 1e-11).unwrap();
        let b = branch.propagator(dt, 1e-11).unwrap();
        let bare = unbranched.propagator(dt, 1e-11).unwrap();
        let (mut x, mut y, mut reference) = (main.zero_state(), branch.zero_state(), unbranched.zero_state());
        let mut changed_pressure = 0.0_f64;
        let mut source_work = 0.0;
        for sample in 0..32 {
            let inflow = if sample < 8 { 2e-8 } else { 0.0 };
            let main_relation = a.multiport_response(&x).unwrap();
            let branch_relation = b.port_response(&y, PortId(9)).unwrap();
            assert_eq!(main_relation.ports, vec![PortId(1), SIDE.id]);
            let h = &main_relation.history_pressure_pa;
            let z = &main_relation.impedance_pa_s_m3;
            let outward = (h[1] - branch_relation.history_pressure_pa + z[2] * inflow)
                / (z[3] + branch_relation.impedance_pa_s_m3);
            let out_a = a.step(&mut x, &[
                PortInput { id: PortId(1), effort: inflow },
                PortInput { id: SIDE.id, effort: -outward },
            ]).unwrap();
            let out_b = b.step(&mut y, &[PortInput { id: PortId(9), effort: outward }]).unwrap();
            let out_bare = bare.step(&mut reference, &[PortInput { id: PortId(1), effort: inflow }]).unwrap();
            let pa = out_a.port_outputs[1];
            let pb = out_b.port_outputs[0];
            let scale = h[1].abs() + (z[2] * inflow).abs() + (z[3] * outward).abs()
                + branch_relation.history_pressure_pa.abs()
                + (branch_relation.impedance_pa_s_m3 * outward).abs();
            assert!((pa.conjugate_rate - pb.conjugate_rate).abs() < 1e-10 * scale.max(1e-30));
            let exchanged_work = outward * (pb.delta_x - pa.delta_x);
            assert!(exchanged_work.abs() < 1e-10 * (scale * outward.abs() * dt).max(1e-30));
            source_work += inflow * out_a.port_outputs[0].delta_x;
            let balance = x.stored_energy_j() + y.stored_energy_j()
                + x.cumulative_dissipated_j() + y.cumulative_dissipated_j() - source_work;
            assert!(balance.abs() < 1e-9 * source_work.abs().max(1e-30));
            changed_pressure = changed_pressure.max(
                (out_a.port_outputs[0].conjugate_rate - out_bare.port_outputs[0].conjugate_rate).abs());
        }
        assert!(changed_pressure > 1e-6, "the resolved side volume must affect the driving pressure");
        assert!(y.stored_energy_j() > 0.0);
        if wall_memory { assert!(x.wall_stored_energy_j() > 0.0 && y.wall_stored_energy_j() > 0.0); }
    }
}

#[test]
fn side_port_material_comes_from_its_face_in_the_original_primitive_gas_field() {
    let field = crate::profile::gas_field::tests::field([300.0, 450.0], [0.8, 0.2], 101325.0);
    let geometry = Duct {
        segments: vec![Segment::Cone { inlet_radius: 0.001, outlet_radius: 0.002, length: 0.57 }],
    };
    let duct = DgDuct::new_with_gas_field(&geometry, &field, DgMesh { elements: 2, order: 2 },
        DgBoundary::Closed, DgBoundary::Closed).unwrap().with_side_ports(&[SIDE]).unwrap();
    let source = duct.side_port_source(SIDE.id).unwrap();
    assert_eq!(source.station_m, 0.285);
    assert_eq!(source.gas, *field.state_at(source.station_m).unwrap().state());
    assert_ne!(source.gas, *duct.inlet_gas());
    let radius = 0.001 + 0.001 * source.station_m / 0.57;
    let area = core::f64::consts::PI * radius * radius;
    assert!((source.bore_area_m2 - area).abs() < 1e-14 * area);
}
