//! Original solid surface geometry in the ordinary source/duct/PCM path.
use super::*;

const SURFACE: &str =
    include_str!("../../../../../../../examples/material_modal_valve/solid-surface-dg.performance");
fn surface_source() -> String {
    SURFACE.replace("audio 96000 96000 2048 1", "audio 96000 96000 512 1")
}

#[test]
fn g1_solid_surface_keeps_actual_nodal_motion_area_and_contact_virtual_work() {
    let built = build(&surface_source(), 37, 1000.0, 20000.0, 400.0);
    let valve = valve(&built);
    let closure = valve.plate_closure().unwrap();
    let source = valve.hereditary_solid().unwrap();
    let surface = source.source().surface_geometry().unwrap();
    let face = surface.faces()[valve.solid_surface_closure().unwrap().contact.lay_triangles[0]];
    assert_eq!(face.nodes, [1, 2, 3]);
    close(face.area_m2, 0.01_f64.powi(2) * 3.0_f64.sqrt() / 2.0, 1e-14);
    close(closure.width_m(), 0.01 * 2.0_f64.sqrt(), 1e-14);
    assert_eq!(closure.contact_nodes(), &[1, 2, 3]);
    assert_eq!(closure.displacement_shapes()[1], [0.0, 0.0]);
    let shape = closure.displacement_shapes();
    // Prescribe an independently known surface displacement at the two free
    // nodes by inverting only their observed 2x2 trace, not the contact owner.
    let gaps = [-2e-6, 3e-6];
    let a = gaps[0] - closure.spec().nodal_rest_gap_m[2];
    let b = gaps[1] - closure.spec().nodal_rest_gap_m[3];
    let det = shape[2][0] * shape[3][1] - shape[2][1] * shape[3][0];
    assert!(det.abs() > 0.1);
    let q = [
        (a * shape[3][1] - b * shape[2][1]) / det,
        (shape[2][0] * b - shape[3][0] * a) / det,
    ];
    let nodal = source
        .expand_nodal(&q, &CancelGate::new_clock_free())
        .unwrap();
    for node in [1, 2, 3] {
        let physical = face
            .outward_normal
            .iter()
            .zip(nodal[node])
            .map(|(n, u)| n * u)
            .sum::<f64>();
        close(
            closure.nodal_gap_m(node, &q).unwrap(),
            closure.spec().nodal_rest_gap_m[node] + physical,
            1e-13,
        );
    }
    let velocity = [-0.02, 0.01];
    let nodal_velocity = source
        .expand_nodal(&velocity, &CancelGate::new_clock_free())
        .unwrap();
    let normal_velocity = |node: usize| {
        face.outward_normal
            .iter()
            .zip(nodal_velocity[node])
            .map(|(normal, value)| normal * value)
            .sum::<f64>()
    };
    close(
        closure.mean_gap_velocity_m_s(&velocity).unwrap(),
        0.5 * (normal_velocity(2) + normal_velocity(3)),
        1e-13,
    );
    close(
        closure.open_area_m2(&q).unwrap(),
        closure.width_m() * gaps[1] * gaps[1] / (2.0 * (gaps[1] - gaps[0])),
        1e-13,
    );
    let probe = closure.probe(&q).unwrap();
    assert_eq!(probe.active_lay_points, 1);
    close(probe.active_lay_area_m2, face.area_m2 / 3.0, 1e-14);
    let k = closure.spec().stiffness_pa_per_m_alpha;
    let energy = face.area_m2 / 3.0 * k * (-gaps[0]).powi(3) / 3.0;
    close(closure.stored_energy_j(&q).unwrap(), energy, 1e-13);
    for mode in 0..2 {
        let mut plus = q;
        let mut minus = q;
        let h = 1e-11;
        plus[mode] += h;
        minus[mode] -= h;
        let derivative = (closure.stored_energy_j(&plus).unwrap()
            - closure.stored_energy_j(&minus).unwrap())
            / (2.0 * h);
        let expected = -face.area_m2 / 3.0 * k * gaps[0].powi(2) * shape[2][mode];
        close(derivative, expected, 1e-8);
    }
}

#[test]
fn g3_native_solid_surface_clearance_and_material_reach_pressure_with_partition_identity() {
    let text = surface_source();
    let built = build(&text, 37, 1000.0, 20000.0, 400.0);
    let old_pressure: Vec<_> = valve(&built)
        .spec()
        .modes
        .iter()
        .map(|m| m.pressure_area_m2)
        .collect();
    let simple = build(&hereditary_source(), 37, 1000.0, 20000.0, 400.0);
    assert_eq!(
        old_pressure,
        valve(&simple)
            .spec()
            .modes
            .iter()
            .map(|m| m.pressure_area_m2)
            .collect::<Vec<_>>()
    );
    let (reference, run) = render(built);
    assert!(reference.iter().any(|p| p.abs() > 1e-15));
    let LoadedAperture::Modal(final_valve) = run.system().aperture() else {
        panic!()
    };
    assert!(final_valve.solid_surface_closure().is_some());
    assert!(
        final_valve
            .material_relaxation()
            .unwrap()
            .memory_sqrt_j()
            .iter()
            .flatten()
            .any(|z| *z != 0.0)
    );
    for (change, c33) in [
        (
            text.replace("solid_gap 0.0000002", "solid_gap 0.000004"),
            400.0,
        ),
        (text.clone(), 2000.0),
    ] {
        let (other, other_run) = render(build(&change, 37, 1000.0, 20000.0, c33));
        assert!(
            reference
                .iter()
                .zip(other)
                .any(|(a, b)| (a - b).abs() > 1e-15)
        );
        assert_ne!(
            run.system().load().state(),
            other_run.system().load().state()
        );
    }
    for block in [1, 512] {
        let (other, other_run) = render(build(&text, block, 1000.0, 20000.0, 400.0));
        assert_eq!(reference, other);
        assert_eq!(
            run.system().load().state(),
            other_run.system().load().state()
        );
        assert_eq!(
            encode_pcm16_wav(&reference, 96000, 1.0).unwrap(),
            encode_pcm16_wav(&other, 96000, 1.0).unwrap()
        );
    }
}

#[test]
fn g0_native_solid_surface_refuses_invented_faces_edges_incomplete_gaps_and_limits() {
    let text = surface_source();
    for altered in [
        text.replace("solid_closure spatial 1 ", "solid_closure spatial 0 "),
        text.replace(
            "0.0001 1048576 manufactured-solid-lay",
            "0.0001 1 manufactured-solid-lay",
        ),
        text.replace("solid_clearances 4", "solid_clearances 3"),
        text.replace("solid_lay_face 1 2 3", "solid_lay_face 0 0 1"),
        text.replace(
            "solid_lay_faces 1\nsolid_lay_face 1 2 3",
            "solid_lay_faces 2\nsolid_lay_face 1 2 3\nsolid_lay_face 1 2 3",
        ),
        text.replace(
            "solid_lay_faces 1\nsolid_lay_face 1 2 3",
            "solid_lay_faces 2\nsolid_lay_face 1 2 3\nsolid_lay_face 0 1 2",
        ),
        text.replace("solid_slit_edge 2 3", "solid_slit_edge 0 2"),
        text.replace(
            "solid_slit_edges 1\nsolid_slit_edge 2 3",
            "solid_slit_edges 2\nsolid_slit_edge 2 3\nsolid_slit_edge 3 2",
        ),
        text.replace("solid_gap 0.0000002", "solid_gap -0.001"),
    ] {
        assert!(
            MaterialModalDgValvePerformance::from_bytes_with_sources(
                altered.as_bytes(),
                37,
                &CancelGate::new_clock_free(),
                None,
                |_| Ok(solid_pack(1000.0, 20000.0, 400.0, "", 100000.0))
            )
            .is_err()
        );
    }
}
