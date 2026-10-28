//! Authored original FE mouths reach the same shell/cavity/exterior owner.
use super::*;

const OPEN_BODY: &str = include_str!("../../open-shell-cavity.body");
const OPEN_PERFORMANCE: &str = include_str!("../../open-shell-cavity.material-strings");

fn open_source() -> input::Spec {
    input::Spec::parse(OPEN_PERFORMANCE.replace("audio 8000 4096", "audio 8000 1024")
        .as_bytes()).unwrap()
}

#[test]
fn declared_opening_faces_reach_joint_bem_and_actual_flow_observation() {
    // The same original reciprocal P1/P0 source now carries independent body
    // and mouth columns. Native execution still owns qualification: construction,
    // strict matrix admission, causal fit and receiver assertions remain hard gates.
    scoped(|cx,gate| {
        let mut loaded = ready(cx,gate,open_source(),OPEN_BODY,700.0);
        let cavity = loaded.performance.bank().source().bank().spatial_cavity().unwrap();
        assert_eq!(*cavity.gas(),loaded.spec.gas);
        assert_eq!(cavity.opening_count(),1);
        let opening = cavity.opening(0).unwrap();
        assert_eq!(opening.boundary_faces,vec![[0,12,14]]);
        assert_eq!(opening.internal_length_m,0.04);
        assert_eq!(opening.internal_resistance_pa_s_m3,0.0);
        assert_eq!(cavity.positions().len(),43);
        assert_eq!(cavity.tetrahedra().len(),80);
        assert_eq!(cavity.coupling().cavity_modes(),1);
        assert_eq!(loaded.performance.bank().system().mode_count(),6);
        let mut scratch = vec![0.0;cavity.observation_state_len()];
        let Some(BodyRadiationReport::ShellBem(load)) = &loaded.body_radiation else {
            panic!("authored openings must retain their original joint shell BEM load");
        };
        assert_eq!(load.opening_ports,1);
        assert_eq!(load.surface_panels,80);
        assert_eq!(load.reciprocal.unknowns,122);
        assert_eq!(load.reciprocal.angular_refinement,0);
        assert_eq!(load.sweep_terms,10_535_395_000);
        assert_eq!(load.gas,loaded.spec.gas);
        let observer = loaded.performance.material_shell_observer().unwrap();
        assert_eq!(observer.body_modes,1);
        assert_eq!(observer.opening_ports,1);
        assert_eq!(observer.gas,loaded.spec.gas);
        let mut reference = observer.transfer().observer().unwrap();
        let body = loaded.spec.strings.len();
        let mut flow_peak = 0.0_f64;
        let mut pressure_peak = 0.0_f64;
        for _ in 0..loaded.spec.samples {
            let mut pressure = [0.0];
            loaded.performance.block(&mut pressure).unwrap();
            let velocity = loaded.performance.bank().system().components()[body].states()[0]
                .velocity_m_sqrt_kg_per_s;
            let neck = loaded.performance.bank().spatial_cavity_neck_observation(0,&mut scratch).unwrap();
            // The shared observer's second column is physical volume flow,
            // not an unscaled modal neck velocity or a separately driven tone.
            assert_eq!(pressure[0],reference.step(&[velocity,neck.volume_flow_m3_s]).unwrap()[0]);
            assert!(pressure[0].is_finite());
            pressure_peak = pressure_peak.max(pressure[0].abs());
            flow_peak = flow_peak.max(neck.volume_flow_m3_s.abs());
        }
        assert!(flow_peak > 0.0);
        assert!(pressure_peak > 1e-14);
        let frame = loaded.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    });
}

#[test]
fn opening_inertia_and_internal_loss_reserve_the_complete_network_before_pack_io() {
    let profile = body::Profile::parse(OPEN_BODY.as_bytes()).unwrap();
    let sealed = body::Profile::parse(BODY.replace(RESOLUTION,
        "interior-fe-v1 2 180 1 256 80 10000000 1e-8").as_bytes()).unwrap();
    assert_eq!((sealed.additional_body_modes(),sealed.additional_body_connections()),(0,1));
    // Reserve the possible internal damper even when this source declares R=0.
    assert_eq!((profile.additional_body_modes(),profile.additional_body_connections()),(1,2));
    let spec = open_source();
    assert_eq!(spec.network_dimensions_with_body_modes(true,
        profile.additional_body_connections(),profile.additional_body_modes()).unwrap(),(6,4,150));
    profile.preflight(&spec).unwrap();
    scoped(|cx,gate| {
        let mut absent = spec.clone();
        absent.body_radiation = None;
        absent.shell_observer = None;
        let mut narrow = spec.clone();
        let Some(input::BodyRadiation::ShellBem(mut radiation)) = narrow.body_radiation else { unreachable!() };
        radiation.maximum_hz = 120.0;
        narrow.body_radiation = Some(input::BodyRadiation::ShellBem(radiation));
        narrow.shell_observer = None;
        for source in [absent,narrow] {
            let mut reads = 0;
            assert!(source.compile_with_body::<NormalizedPack>(cx,gate,Some(profile.clone()),|_| {
                reads += 1;
                Err("material source must not be read before opening/exterior admission".into())
            }).is_err());
            assert_eq!(reads,0);
        }
    });
}

#[test]
fn opening_records_keep_strict_original_face_counts_and_no_ignored_tail() {
    for (from,to) in [
        ("air-openings 1","air-openings 0"),
        ("air-openings 1","air-openings 9"),
        ("air-opening 0.04 0 1","air-opening 0 0 1"),
        ("air-opening 0.04 0 1","air-opening 0.04 -1 1"),
        ("air-opening 0.04 0 1","air-opening NaN 0 1"),
        ("air-opening 0.04 0 1","air-opening 0.04 0 0"),
        ("air-opening 0.04 0 1","air-opening 0.04 0 1 extra"),
        ("air-opening-face 0 12 14","air-opening-face 0 12 43"),
        ("air-opening-face 0 12 14","air-opening-face 0 12 12"),
        ("air-opening-face 0 12 14","air-opening-face 0 12"),
        ("air-opening-face 0 12 14","air-opening-face 0 12 14 extra"),
    ] {
        assert!(body::Profile::parse(OPEN_BODY.replace(from,to).as_bytes()).is_err(),"{to}");
    }
    let repeated = format!("{}air-opening 0.04 0 1\nair-opening-face 14 0 12\n",
        OPEN_BODY.replace("air-openings 1","air-openings 2"));
    assert!(body::Profile::parse(repeated.as_bytes()).is_err());
    assert!(body::Profile::parse(format!("{OPEN_BODY}ignored\n").as_bytes()).is_err());
}
