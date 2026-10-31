//! Actual sourced physical rotations, not independently authored modal maps.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::shell::ShellRotationGeometry;

#[test]
fn nodal_force_uses_actual_source_directors_and_refuses_ambiguous_offset_spin() {
    use fs_couple::bowed_string::runtime::reciprocal::bank::shell::ShellPointGeometry;
    with_cx(|cx,gate| {
        let bank = build_shell(cx,gate,700.0,0.02,0.17);
        let source = bank.source();
        let node = 5;
        let director = source.material_shell().unwrap().directors.vectors()[node];
        for arm_m in [[0.0;3],director.map(|v|0.02*v)] {
            let geometry = ShellPointGeometry { shell_component:0,node,direction:[0.0,1.0,0.0],arm_m };
            let point = source.shell_point(&geometry,gate).unwrap();
            let row = fs_plate::shell::point::nodal_point_row(geometry.direction,arm_m).unwrap();
            for (mode,&shape) in source.material_shell_image().unwrap().physical_modes().iter().zip(&point.body.shapes) {
                let expected: f64 = (0..6).map(|i|mode[6*node+i]*row[i]).sum();
                assert!((shape-expected).abs()<1e-12*expected.abs().max(1.0));
            }
        }
        let geometry = ShellPointGeometry { shell_component:0,node,direction:[0.0,1.0,0.0],arm_m:[0.02,0.0,0.0] };
        geometry.validate().unwrap(); // Pure syntax cannot know this source's director.
        assert!(source.shell_point(&geometry,gate).unwrap_err().to_string().contains("director-axis drilling"));
    });
}

#[test]
fn physical_rotation_uses_the_original_material_modes_and_refuses_detached_geometry() {
    with_cx(|cx,gate| {
        let geometry = ShellRotationGeometry { shell_component: 0,triangle: 0,
            barycentric: [0.2,0.3,0.5],axis: [0.6,0.8,0.0] };
        let original = build_shell(cx,gate,700.0,0.02,0.17);
        let dense = build_shell(cx,gate,1400.0,0.02,0.17);
        for bank in [&original,&dense] {
            let source = bank.source();
            let rotation = source.shell_rotation(&geometry,gate).unwrap();
            let image = source.material_shell_image().unwrap();
            let original = source.material_shell().unwrap();
            let triangle = original.specimen.mesh().tris[geometry.triangle];
            let facet = original.specimen.mesh().facet(geometry.triangle).unwrap();
            let rows = fs_plate::shell::point::triangle_rotation_rows(geometry.barycentric,geometry.axis,
                facet.frame,facet.gradient,triangle.map(|node|original.directors.vectors()[node])).unwrap();
            assert_eq!(rotation.component,0);
            assert_eq!(rotation.angular_shapes.len(),image.physical_modes().len());
            for (mode,&angle) in image.physical_modes().iter().zip(&rotation.angular_shapes) {
                let expected: f64 = (0..3).map(|a|(0..6)
                    .map(|i|rows[a][i]*mode[6*triangle[a]+i]).sum::<f64>()).sum();
                assert!((angle-expected).abs() < 1e-12*expected.abs().max(1.0));
            }
        }
        assert_ne!(original.source().shell_rotation(&geometry,gate).unwrap().angular_shapes,
            dense.source().shell_rotation(&geometry,gate).unwrap().angular_shapes);
        for bad in [
            ShellRotationGeometry { shell_component: 1,..geometry.clone() },
            ShellRotationGeometry { triangle: usize::MAX,..geometry.clone() },
            ShellRotationGeometry { axis: [0.0;3],..geometry.clone() },
            ShellRotationGeometry { barycentric: [-0.2,0.3,0.9],..geometry.clone() },
            // Original node zero is clamped in all physical rotations.
            ShellRotationGeometry { barycentric: [1.0,0.0,0.0],..geometry.clone() },
        ] {
            assert!(original.source().shell_rotation(&bad,gate).is_err());
        }
        let cancelled = CancelGate::new_clock_free(); cancelled.request();
        assert!(original.source().shell_rotation(&geometry,&cancelled).is_err());
    });
}
