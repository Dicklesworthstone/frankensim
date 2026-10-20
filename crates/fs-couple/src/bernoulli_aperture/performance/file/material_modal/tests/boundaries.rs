//! Original material key -> both coupled media -> actual multirate pressure.
use super::*;
use crate::bernoulli_aperture::material_boundary::boundary_closing_force_path;
use crate::bernoulli_aperture::two_sided::observer::schedule::boundaries::source_pressure_path;
use fs_scenario::gesture::{
    GestureEvent, GestureSchedule, GestureTarget, GestureTrack, GestureValue,
};

fn keyed(force: bool) -> String {
    let original = input();
    let header = original.split_once("\nschedule\n").unwrap().0
        .replace("network 3 2 1048576\nduct_node inlet\nduct_node junction\nduct_node reflection 0\n",
            "network 4 3 1048576\nduct_node inlet\nduct_node junction\nduct_node reflection 0\nduct_node driven\n")
        .replace("\npropagation ", "\nduct_section 1 3 0.017 0.003 0.002\npropagation ")
        .replace("\ncompile_limits ", concat!(
            "\ndownstream_boundaries\nboundaries 1\nboundary 3 instrument key-1 0 0\n",
            "boundary_plate key.matpack 293.15 1000 0.025 0.01 0.0003 6 2 0.0002 0.0002 0\n",
            "boundary_mode 0 600 0 0.25 0.1\nboundary_contact 10000000000 2 0.5\n",
            "boundary_slit 0.0003\nboundary_loss elastic\ncontrol_source instrument reed\ncompile_limits "));
    let schedule = GestureSchedule::try_new_v2(
        96000,
        vec![
            GestureTrack {
                id: "breath".into(),
                target: GestureTarget::Parameter(source_pressure_path("instrument", "reed")),
                initial: GestureValue::PressurePa(5.0),
                events: vec![GestureEvent {
                    time_s: 0.004,
                    transition_s: 0.001,
                    value: GestureValue::PressurePa(0.0),
                }],
            },
            GestureTrack {
                id: "key".into(),
                target: GestureTarget::Parameter(boundary_closing_force_path(
                    "instrument",
                    "key-1",
                )),
                initial: GestureValue::ForceN(0.0),
                events: if force {
                    vec![
                        GestureEvent {
                            time_s: 1.0 / 96000.0,
                            transition_s: 0.0,
                            value: GestureValue::ForceN(0.0005),
                        },
                        GestureEvent {
                            time_s: 0.003,
                            transition_s: 0.001,
                            value: GestureValue::ForceN(0.0),
                        },
                    ]
                } else {
                    vec![]
                },
            },
        ],
    )
    .unwrap();
    format!(
        "{header}\nschedule\n{}",
        String::from_utf8(schedule.to_canonical_bytes()).unwrap()
    )
}
fn load(text: &str, block: usize, key_density: f64) -> MaterialModalValvePerformance {
    MaterialModalValvePerformance::from_bytes_with_sources(
        text.as_bytes(),
        block,
        &CancelGate::new_clock_free(),
        None,
        |path| {
            Ok(pack(
                if path == "key.matpack" {
                    key_density
                } else {
                    900.0
                },
                false,
            ))
        },
    )
    .unwrap()
}
fn key(source: &ScheduledTwoSidedPressure) -> &crate::bernoulli_aperture::dynamic::DynamicAperture {
    source
        .source()
        .system()
        .boundary_apertures(ValveSide::Downstream)
        .unwrap()
        .aperture(3)
        .unwrap()
}

fn key_loss_pack(eta: f64, upper: f64) -> NormalizedMaterialCardPack {
    let source = pack(900.0, false);
    let mut claims = source.card().claims().clone();
    let observations = claims.claims_for("density")[0].1.observations.clone();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(
            crate::string_specimen::KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY,
            DynViscosity::DIMS,
        ),
        value: PropertyValue::Scalar { value: eta, dims: DynViscosity::DIMS },
        validity: ValidityDomain::unconstrained()
            .with("T", 290.0, 300.0).with("omega", 0.0, upper),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        provenance: Provenance { source: "manufactured material-key viscosity".into(),
            license: "CC0-1.0".into(), artifact: None },
        observations,
    }).unwrap();
    let normalized = NormalizedPack::new("synthetic-key", "synthetic-v1",
        hash_domain("org.frankensim.test.material-key", &eta.to_le_bytes()),
        "CC0-1.0", claims, vec![], vec![]).unwrap();
    NormalizedMaterialCardPack::new(source.card().id().clone(), normalized).unwrap()
}

fn load_key_loss(text: &str, eta: f64, upper: f64)
    -> Result<MaterialModalValvePerformance, PlateValveInputError>
{
    MaterialModalValvePerformance::from_bytes_with_sources(text.as_bytes(), 37,
        &CancelGate::new_clock_free(), None, |path| Ok(if path == "key.matpack" {
            key_loss_pack(eta, upper)
        } else { pack(900.0, false) }))
}

#[test]
fn material_key_viscosity_reaches_the_original_plate_and_rendered_pressure() {
    let text = keyed(true).replace("boundary_loss elastic", "boundary_loss kelvin-voigt");
    let low = load_key_loss(&text, 2e5, 100000.0).unwrap();
    let high = load_key_loss(&text, 8e5, 100000.0).unwrap();
    // Only the key's original material property changes. Reed modes, gas,
    // geometry and every authored force/pressure sample remain identical.
    assert_eq!(valve(&low).spec(), valve(&high).spec());
    assert_eq!(valve(&low).viscosity_matrix_n_s_m(), valve(&high).viscosity_matrix_n_s_m());
    for (source, eta) in [(&low, 2e5), (&high, 8e5)] {
        let aperture = key(source.renderer());
        let plate = aperture.plate_reduction().unwrap();
        let chart = plate.material_chart().unwrap();
        let material = &chart.regions()[0].input.material;
        let law = &chart.bending_viscosity().unwrap()[0];
        assert_eq!(material.card_identity(), source.boundary_records()[0].1.pack().card().content_hash());
        assert_eq!(material.query_point(), &[("T".into(), 293.15), ("omega".into(), 1000.0)]);
        assert_eq!(law.material_state_identity, Some(material.identity()));
        assert_eq!(law.viscosity_pa_s, eta);
        assert_eq!(plate.options().damping_ratio, 0.0);
        assert_eq!(aperture.spec().damping_ratio, plate.damping_ratio());
        assert!(plate.damping_ratio() > 0.0);
        assert!(aperture.relaxation().is_none());
    }
    let a = key(low.renderer()).plate_reduction().unwrap();
    let b = key(high.renderer()).plate_reduction().unwrap();
    assert_eq!((a.mass_kg(), a.stiffness_n_m(), a.pressure_area_m2()),
        (b.mass_kg(), b.stiffness_n_m(), b.pressure_area_m2()));
    assert!((b.damping_ratio() / a.damping_ratio() - 4.0).abs() < 1e-12);
    let (pa, state_a) = render(low);
    let (pb, state_b) = render(high);
    assert!(pa.iter().chain(&pb).all(|p| p.is_finite()));
    assert!(pa.iter().zip(&pb).any(|(a, b)| (a - b).abs() > 1e-15));
    assert_ne!(key(&state_a).state(), key(&state_b).state());
    assert_eq!(key(&state_a).accepted_steps(), 512 * 16);
    assert_eq!(key(&state_b).accepted_steps(), 512 * 16);
}

#[test]
fn explicit_key_viscosity_preserves_zero_loss_and_refuses_missing_or_invalid_source() {
    let elastic = keyed(true);
    let selected = elastic.replace("boundary_loss elastic", "boundary_loss kelvin-voigt");
    let zero = load_key_loss(&selected, 0.0, 100000.0).unwrap();
    assert!(key(zero.renderer()).plate_reduction().unwrap().material_chart().unwrap()
        .bending_viscosity().is_some());
    assert_eq!(render(zero).0, render(load_key_loss(&elastic, 0.0, 100000.0).unwrap()).0);
    let missing = MaterialModalValvePerformance::from_bytes_with_sources(selected.as_bytes(), 37,
        &CancelGate::new_clock_free(), None, |_| Ok(pack(900.0, false)));
    assert!(matches!(missing, Err(PlateValveInputError::Material { .. })));
    assert!(load_key_loss(&selected, -1.0, 100000.0).is_err());
    // The source query (1000 rad/s) is admitted, but its declared validity
    // does not cover the original 0..600 Hz mode search band.
    assert!(load_key_loss(&selected, 0.0, 3000.0).is_err());
    assert!(load_key_loss(&selected.replace("boundary_loss kelvin-voigt",
        "boundary_loss kelvin-voigt maxwell"), 2e5, 100000.0).is_err());
}

#[test]
fn physical_material_keys_and_subframe_force_events_reach_the_original_modal_microphone() {
    let text = keyed(true);
    let loaded = load(&text, 37, 900.0);
    assert_eq!(loaded.info().modes, 6);
    assert_eq!(loaded.boundary_records().len(), 1);
    assert_eq!(loaded.boundary_records()[0].0, ValveSide::Downstream);
    assert_eq!(loaded.boundary_records()[0].1.node(), 3);
    assert_eq!(
        loaded.boundary_records()[0].1.closing_force(),
        &boundary_closing_force_path("instrument", "key-1")
    );
    let owner = key(loaded.renderer());
    assert_eq!(owner.spec().max_steps, 512 * 16);
    assert_eq!(
        owner.viscous_slit().unwrap().gas(),
        loaded
            .renderer()
            .source()
            .system()
            .boundary_apertures(ValveSide::Downstream)
            .unwrap()
            .gas(3)
    );
    assert!(
        loaded
            .renderer()
            .pending_boundary_controls(ValveSide::Downstream)
            .iter()
            .any(|event| event.sample % 16 == 8),
        "force changes must retain sub-audio-frame timing"
    );
    let (pressure, rendered) = render(loaded);
    assert!(pressure.iter().all(|p| p.is_finite()));
    assert!(pressure.iter().any(|p| p.abs() > 1e-15));
    assert_eq!(key(&rendered).accepted_steps(), 512 * 16);
    assert!(
        rendered
            .pending_boundary_controls(ValveSide::Downstream)
            .is_empty()
    );
    assert_eq!(
        rendered.held_boundary_drives(ValveSide::Downstream)[0].opening_force_n,
        0.0
    );
    assert!(
        rendered
            .last_boundary_frame(ValveSide::Downstream)
            .is_some()
    );
    assert_eq!(pressure, render(load(&text, 1, 900.0)).0);
    assert_eq!(pressure, render(load(&text, 128, 900.0)).0);
    let (unchanged_force, _) = render(load(&keyed(false), 37, 900.0));
    let (material, changed) = render(load(&text, 37, 1800.0));
    assert_ne!(pressure, unchanged_force);
    assert_ne!(pressure, material);
    assert_ne!(key(&rendered).state(), key(&changed).state());
    assert_ne!(
        encode_pcm16_wav(&pressure, 48000, 0.005).unwrap().0,
        encode_pcm16_wav(&material, 48000, 0.005).unwrap().0
    );
}

#[test]
fn keys_share_original_material_cache_and_single_blade_owner() {
    let text = keyed(true).replace("boundary_plate key.matpack", "boundary_plate blade.matpack");
    let mut loads = 0;
    let source = MaterialModalValvePerformance::from_bytes_with_sources(
        text.as_bytes(),
        37,
        &CancelGate::new_clock_free(),
        None,
        |path| {
            assert_eq!(path, "blade.matpack");
            loads += 1;
            Ok(pack(900.0, false))
        },
    )
    .unwrap();
    assert_eq!(loads, 1);
    assert_eq!(
        source.boundary_records()[0].1.pack(),
        source.material_records()[0].pack()
    );
    let first = text.find("\nblade ").unwrap();
    let second = first + 1 + text[first + 1..].find("\nblade ").unwrap();
    let tail = second + text[second..].find("\ndownstream_boundaries\n").unwrap();
    let single = format!("{}{}", &text[..second], &text[tail..]).replace("blades 2", "blades 1");
    let mut source = load(&single, 37, 900.0).into_renderer();
    let LoadedAperture::Modal(valve) = source.source().system().valve() else {
        panic!("original modal valve")
    };
    assert_eq!(valve.material_plates().len(), 1);
    assert_eq!(valve.modal_states().len(), 3);
    source.block(&mut [0.0; 37]).unwrap();
    assert_eq!(key(&source).accepted_steps(), 37 * 16);
}

#[test]
fn incomplete_keys_and_unbound_typed_consumers_refuse_before_playback() {
    let text = keyed(true);
    let gate = CancelGate::new_clock_free();
    for altered in [
        text.replace("boundary 3 instrument", "boundary 2 instrument"),
        text.replace(
            "control_source instrument reed",
            "control_source instrument unbound",
        ),
        text.replace(
            "boundary 3 instrument key-1",
            "boundary 3 instrument unknown-key",
        ),
        text.replace(
            "boundary_plate key.matpack",
            "boundary_plate ../key.matpack",
        ),
        text.replace("compile_limits 100000 8192", "compile_limits 1 8192"),
    ] {
        assert!(
            MaterialModalValvePerformance::from_bytes_with_sources(
                altered.as_bytes(),
                37,
                &gate,
                None,
                |_| Ok(pack(900.0, false))
            )
            .is_err()
        );
    }
    let declared = "downstream_boundaries\nboundaries 1";
    let start = text.find(declared).unwrap();
    let end = start + text[start..].find("control_source").unwrap();
    let missing = format!("{}{}", &text[..start], &text[end..]);
    assert!(
        MaterialModalValvePerformance::from_bytes_with_sources(
            missing.as_bytes(),
            37,
            &gate,
            None,
            |_| Ok(pack(900.0, false))
        )
        .is_err()
    );
}
