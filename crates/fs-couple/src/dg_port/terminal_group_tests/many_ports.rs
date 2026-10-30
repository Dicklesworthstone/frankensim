use super::*;

#[test]
fn sixteen_side_ports_keep_full_mutual_reaction_closure_and_acoustic_work() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let side: Vec<_> = (0..16)
        .map(|i| DgSidePort {
            id: PortId(100 + i as u32),
            downstream_element: i + 1,
        })
        .collect();
    let duct = DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.001,
                length: 0.05,
            }],
        },
        &gas,
        DgMesh {
            elements: 17,
            order: 1,
        },
        DgBoundary::Port(INLET),
        DgBoundary::Port(OUTLET),
    )
    .unwrap()
    .with_side_ports(&side)
    .unwrap();
    let mut ports: Vec<_> = side.iter().map(|port| port.id).collect();
    ports.push(OUTLET);
    let n = ports.len();
    // A complete reciprocal positive matrix, including the sixteenth hole.
    // This manufactured load isolates port algebra from a numerical BEM fit.
    let direct: Vec<_> = (0..n * n)
        .map(|i| {
            let (a, b) = (i / n, i % n);
            (if a == b { 2e7 } else { 0.0 }) + 1e5 * (a + 1) as f64 * (b + 1) as f64
        })
        .collect();
    let source = MimoModel {
        ports: n,
        poles: vec![],
        residues: vec![],
        d: direct.clone(),
        e: vec![0.0; n * n],
        reciprocal: true,
    };
    let terminal = PreparedMimoLoad::new(
        &source,
        DT,
        MimoLoadLimits {
            band_rad_s: [1.0, core::f64::consts::PI / DT],
            sampled_passivity_points: 17,
            max_states: 1,
            max_memory_bytes: 1 << 20,
        },
    )
    .unwrap();
    let mut load = DgPortLoad::new_with_terminal_ports(
        duct,
        DT,
        1e-11,
        MAX_DG_PROPAGATOR_DOFS,
        ports.clone(),
        terminal,
    )
    .unwrap()
    .with_terminal_channels(
        vec![DgTerminalChannelInput {
            port: ports[0],
            kind: DgTerminalChannelKind::Closed,
        }],
        0,
    )
    .unwrap();
    let mut inlet_work = 0.0;
    let mut terminal_work = 0.0;
    let mut last_hole_live = false;
    for sample in 0..48 {
        let flow = if sample < 8 { 1e-9 } else { 0.0 };
        let trial = load.preview_step(flow).unwrap();
        let frames = trial.terminal_frames();
        assert_eq!(frames.len(), 17);
        assert_eq!(
            frames.iter().map(|frame| frame.port).collect::<Vec<_>>(),
            ports
        );
        assert_eq!(frames[0].outward_flow_m3_s, 0.0);
        last_hole_live |= frames[15].outward_flow_m3_s.abs() > 1e-16;
        for (i, frame) in frames.iter().enumerate().skip(1) {
            let expected: f64 = direct[i * n..(i + 1) * n]
                .iter()
                .zip(frames)
                .map(|(z, other)| z * other.outward_flow_m3_s)
                .sum();
            near(
                frame.mean_pressure_pa,
                expected,
                frame.mean_pressure_pa.abs() + expected.abs(),
            );
            near(
                frame.port_work_residual_j,
                0.0,
                frame.port_work_j.abs() + frame.load_work_j.abs(),
            );
        }
        inlet_work += trial.inlet_work_j();
        terminal_work += sum_terminal_work(frames);
        trial.accept();
        let stored = load.state().stored_energy_j();
        let dissipated = load.state().cumulative_dissipated_j();
        near(
            stored + dissipated + terminal_work,
            inlet_work,
            stored.abs() + dissipated.abs() + terminal_work.abs() + inlet_work.abs(),
        );
    }
    assert!(
        last_hole_live,
        "the final side port must remain in the actual flow solve"
    );
}
