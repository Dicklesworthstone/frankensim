use super::*;

const SPATIAL: &str =
    include_str!("../../../../../../examples/material_modal_valve/paired-spatial.performance");
const TENSOR: &str = "loss tensor eta11 eta22 eta66 eta12 eta16 eta26 1 1000000";
fn source() -> String {
    SPATIAL.replace(
        "audio 768000 48000 2400 0.005",
        "audio 768000 48000 256 0.005",
    )
}

#[test]
fn all_original_modes_and_spatial_geometry_reach_the_scheduled_pressure_path() {
    let text = source();
    let p = build(&text, 37);
    let a = valve(&p);
    let closure = a.plate_closure().unwrap();
    assert_eq!(closure.spec().nodal_rest_gap_m.len(), 21);
    assert_eq!(closure.spec().lay_triangles.len(), 24);
    assert_eq!(closure.spec().provenance, "manufactured-planar-blade-pair");
    assert!(a.viscosity_matrix_n_s_m().is_some());
    for plate in a.material_plates() {
        for (all, local) in closure
            .displacement_shapes()
            .iter()
            .zip(plate.displacement_shapes())
        {
            assert_eq!(&all[plate.mode_range()], local.as_slice());
        }
    }
    let (pressure, source) = render(p);
    let (partitioned, _) = render(build(&text, 1));
    assert_eq!(pressure, partitioned);
    assert!(pressure.iter().all(|p| p.is_finite()));
    assert!(pressure.iter().any(|p| *p != 0.0));
    let (changed_gap, _) = render(build(
        &text.replace("gap 0.000020000", "gap 0.000030000"),
        37,
    ));
    assert!(
        pressure
            .iter()
            .zip(changed_gap)
            .any(|(a, b)| (a - b).abs() > 1e-15)
    );
    let LoadedAperture::Modal(a) = source.source().system().valve() else {
        panic!("original material modes");
    };
    assert_eq!(a.accepted_steps(), 4096);
    let q: Vec<_> = a.modal_states().iter().map(|s| s.displacement_m).collect();
    assert!(
        a.plate_closure()
            .unwrap()
            .probe(&q)
            .unwrap()
            .open_area_m2
            .is_finite()
    );
}

#[test]
fn regional_memory_and_source_material_changes_survive_spatial_contact() {
    let text = source()
        .replace(
            TENSOR,
            "loss regional-maxwell relaxed 1 boundary_modulus boundary_time",
        )
        .replacen("section blade.matpack", "section altered.matpack", 1);
    let load = |modulus, time| {
        MaterialModalValvePerformance::from_bytes_with_sources(
            text.as_bytes(),
            37,
            &CancelGate::new_clock_free(),
            None,
            |path| {
                Ok(if path == "altered.matpack" {
                    pack_with_relaxation(900.0, false, modulus, time)
                } else {
                    pack(900.0, false)
                })
            },
        )
        .unwrap()
    };
    let (a, system) = render(load(1e8, 0.001));
    let (b, _) = render(load(4e8, 0.0002));
    assert!(a.iter().zip(b).any(|(a, b)| (a - b).abs() > 1e-15));
    let LoadedAperture::Modal(valve) = system.source().system().valve() else {
        panic!("original modes");
    };
    assert!(valve.plate_closure().is_some());
    let memory = valve.material_relaxation().unwrap();
    assert!(memory.is_regional());
    assert!(memory.memory_sqrt_j().iter().flatten().any(|z| *z != 0.0));
}

#[test]
fn spatial_file_refuses_incomplete_or_repeated_physical_interface_data() {
    for text in [
        source().replace("clearances 21", "clearances 20"),
        source().replace("lay_triangle 23", "lay_triangle 22"),
        source().replace("lay_triangle 23", "lay_triangle 24"),
        source().replace("gap 0.000020000", "gap -1"),
    ] {
        assert!(
            MaterialModalValvePerformance::from_bytes_with_sources(
                text.as_bytes(),
                37,
                &CancelGate::new_clock_free(),
                None,
                |_| Ok(pack(900.0, false))
            )
            .is_err()
        );
    }
}
