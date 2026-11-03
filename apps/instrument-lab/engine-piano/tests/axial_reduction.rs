//! G1/G3: axial forces participate in the actual board reduction, before
//! geometry is projected into its retained mechanical/acoustic basis.
use fs_live_piano::{Description, refmods::{
    board_geometry::{BoardGeometry, motion::{SourceBridgeFrame, SourceBridgePort}, ritz::RitzOptions},
    linear::BoardMode, steinway_d::Preset, steinway_scale,
}};

fn panel() -> Preset {
    // One free node: vertical bending and two independent rocking directions.
    // The horizontal bridge supplies its actual arm and signed string axis.
    // This is a manufactured integration fixture, not a piano material claim.
    let mut text = String::from("frankensim-board-geometry-si-v1\nsource,estimated,three-DOF axial-target panel\nsupport,clamped\npretension,0\ndamping,0.01\n");
    for (i, (x, y)) in [(0.0, 0.5), (0.5, 0.0), (1.0, 0.5), (0.5, 1.0), (0.5, 0.5)].iter().enumerate() {
        text.push_str(&format!("node,{i},{x},{y}\n"));
    }
    for (i, t) in [[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]].iter().enumerate() {
        text.push_str(&format!("triangle,{i},{},{},{},0.08,450,1000000,80000,0.3,60000,0\n", t[0], t[1], t[2]));
    }
    text.push_str("fixed,0\nfixed,1\nfixed,2\nfixed,3\nbridge,69,0,0,0,1\n");
    text.push_str("stiffener,1000000,80000,0.0001,0.0000000008333333333333333,0.000000001,0.01,450,0,4,2\n");
    Preset { geometry: text, obj: String::new() }
}

fn compliance(board: &[BoardMode], row: &[f64]) -> f64 {
    assert_eq!(board.len(), row.len());
    board.iter().zip(row).map(|(mode, a)|
        a * a / (std::f64::consts::TAU * mode.frequency_hz).powi(2)).sum()
}

#[test]
fn axial_source_motion_survives_live_pivoted_and_pod_reduction() {
    let preset = panel();
    let geometry = BoardGeometry::read(&preset.geometry).unwrap();
    let attachment = geometry.bridge_attachments(&[69]).unwrap()[0];
    assert_eq!(attachment.axis_out, [0.0, 1.0, 0.0]);
    let source = geometry.prepare_with_motion(&[69], 1000.0).unwrap();
    assert_eq!((source.free_dofs, source.modes.len()), (3, 3));
    let motion = source.motion.as_ref().unwrap();
    let axial = attachment.axial_row(motion).unwrap();
    let (lateral, _) = motion.project_at(attachment.triangle, attachment.weights,
        attachment.arm_m, [1.0, 0.0, 0.0]).unwrap();
    assert!(axial.iter().enumerate().any(|(i, a)| a.abs() > 1e-8
        && source.modes[i].bridge[48].abs() + lateral[i].abs() < 1e-10 * a.abs()),
        "fixture must contain axial motion invisible to both transverse forces");
    let expected = compliance(&source.modes, &axial);
    assert!(expected > 0.0);
    let course = steinway_scale::courses().unwrap().into_iter().find(|c| c.midi == 69).unwrap();
    let primary = SourceBridgePort { triangle: attachment.triangle, weights: attachment.weights,
        arm_m: attachment.arm_m, direction: [0.0, 0.0, 1.0] };
    let transverse_only = [SourceBridgeFrame { midi: 69, primary,
        secondary: SourceBridgePort { direction: [1.0, 0.0, 0.0], ..primary }, axial: None }];
    for spec in ["3,0,10,30", "pod,3,0,4,30"] {
        let options = RitzOptions::parse(spec).unwrap();
        let missing = geometry.prepare_reduced_with_ports(&[69], 1000.0, true,
            false, false, false, &options, Some(&transverse_only)).unwrap();
        let missing_row = attachment.axial_row(missing.motion.as_ref().unwrap()).unwrap();
        assert!(compliance(&missing.modes, &missing_row) < 1e-10 * expected,
            "{spec}: transverse snapshots must omit this axial response");

        let desc = Description::with_axial_reduction_on(&preset, vec![course], 1000.0, spec).unwrap();
        assert_eq!(desc.board.len(), source.modes.len(), "{spec}");
        assert!(desc.board.len() > missing.modes.len(), "{spec}: axial force must select its missing direction");
        let actual = compliance(&desc.board, &desc.bridge_axial_rows[&69]);
        assert!((actual - expected).abs() < 1e-9 * expected,
            "{spec}: retained axial compliance {actual:e} differs from physical source {expected:e}");
        // The retained surface and axial row belong to one modal basis, and
        // an ordinary linear preparation still selects only its original port.
        assert!(desc.surface.iter().all(|p| p.mode_shape.len() == desc.board.len()));
        assert_eq!(desc.damping.as_ref().unwrap().len(), desc.board.len().pow(2));
        let linear = Description::with_reduction_on(&preset, vec![course], 1000.0, spec).unwrap();
        assert!(linear.bridge_axial_rows.is_empty());
        assert!(linear.board.len() < desc.board.len());
    }
}
