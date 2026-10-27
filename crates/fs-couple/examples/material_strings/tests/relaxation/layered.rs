//! G0/G3: source-only bonded spectra reach the complete native performer.
use super::*;
use fs_material::section::relaxation::ConcentricPronySection;

fn selected_layered() -> input::Spec {
    let mut s = selected();
    s.strings[0].section = input::WireSection::Bonded {
        core_radius_m: 0.0002,
        sleeves: vec![input::WireLayer {
            pack: "layered-sleeve.fsmatpk".into(), outer_radius_m: 0.0003,
        }],
    };
    s
}

fn original_section(source: &MaterialStringSource) -> &ConcentricPronySection {
    match source {
        MaterialStringSource::WithPolarizations(p) => original_section(p.source()),
        MaterialStringSource::WithTorsion(t) => original_section(t.bending()),
        MaterialStringSource::ConcentricProny { section, .. } => section,
        _ => panic!("bonded source lost its original per-layer constitutive law"),
    }
}

fn build(cx: &Cx<'_>, gate: &CancelGate, s: input::Spec, sleeve_time: f64) -> Ready {
    s.compile(cx, gate, |name| Ok(source_pack(name,
        if name == "layered-sleeve.fsmatpk" { 40e9 } else { 80e9 },
        if name == "layered-sleeve.fsmatpk" { sleeve_time } else { 0.001 },
        200_000.0, None))).unwrap()
}

#[test]
fn sleeve_spectrum_changes_pressure_and_pcm_without_viscosity_or_backbone_retuning() {
    scoped(|cx, gate| {
        let mut a = build(cx, gate, selected_layered(), 0.002);
        let mut b = build(cx, gate, selected_layered(), 0.008);
        let section = original_section(&a.performance.bank().source().bank().strings()[0].source);
        assert!(section.elastic().layers().iter().all(|layer| layer.material.property(ETA).is_none()));
        assert_eq!(section.spectra()[0].terms, vec![(80e9, 0.001)]);
        assert_eq!(section.spectra()[1].terms, vec![(40e9, 0.002)]);
        assert_eq!(a.performance.bank().system().material_memory().len(), 8);
        assert!(a.performance.bank().system().material_memory().iter().any(|z| *z != 0.0));
        assert_eq!(a.spec.bytes, b.spec.bytes);
        assert_eq!(a.catalog["wire.fsmatpk"].pack.content_hash(),
            b.catalog["wire.fsmatpk"].pack.content_hash());
        assert_eq!(a.catalog["body.fsmatpk"].pack.content_hash(),
            b.catalog["body.fsmatpk"].pack.content_hash());
        for (a, b) in a.performance.bank().system().components().iter()
            .zip(b.performance.bank().system().components())
        {
            assert_eq!(a.modes(), b.modes());
        }
        let p = pressure(&mut a, 1);
        let q = pressure(&mut b, 37);
        assert!(p.iter().any(|p| p.abs() > 1e-12));
        assert!(p.iter().zip(&q).any(|(p, q)| (p - q).abs() > 1e-12));
        assert_eq!(p, pressure(&mut build(cx, gate, selected_layered(), 0.002), 37));
        let encode = |p: &[f64]| fs_couple::pcm_wav::encode_pcm16_wav(
            p, a.spec.rate, a.spec.full_scale).unwrap().0;
        assert_ne!(encode(&p)[44..], encode(&q)[44..]);
        for image in [&a, &b] {
            let frame = image.performance.bank().system().last_frame().unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
    });
}

#[test]
fn bonded_source_memory_survives_gas_torsion_two_planes_and_axial_composition() {
    scoped(|cx, gate| {
        let mut s = selected_layered();
        s.samples = 96;
        s.strings[0].bending_modes = 2;
        s.strings[0].twist_modes = 1;
        s.strings[0].polarizations = Some(
            fs_couple::bowed_string::runtime::reciprocal::bank::polarization::TransverseDirections::try_new(
                [1.0, 0.0], [0.8, 0.6]).unwrap());
        s.cylinder_gas_loading = Some(CylinderFlowResolution {
            minimum_wavenumber_radius: 1e-3, maximum_wavenumber_radius: 30.0, intervals: 2,
        });
        s.longitudinal_waves = Some((LongitudinalStringResolution {
            axial_modes: vec![1], maximum_slope: 0.2, maximum_strain: 0.01,
        }, NonlinearModalConfig { max_iterations: 64, maximum_generalized_force: 1e5,
            force_absolute_tolerance: 1e-11, force_relative_tolerance: 1e-10 }));
        s.forces = None;
        s.schedule = GestureSchedule::try_new(1000, vec![GestureTrack {
            id: "bow".into(), target: GestureTarget::BowStroke { string: 0 },
            initial: GestureValue::Bow { velocity_m_per_s: 0.01, normal_force_n: 0.0001, station: 0.2 },
            events: vec![],
        }]).unwrap();
        let mut a = build(cx, gate, s.clone(), 0.002);
        let mut b = build(cx, gate, s, 0.002);
        let source = a.performance.bank().source().bank();
        assert!(source.cylinder_gas_loading().is_some());
        assert!(source.longitudinal_waves().is_some());
        assert!(source.material_cavity().is_some());
        let section = original_section(&source.strings()[0].source);
        assert_eq!(section.spectra()[1].terms, vec![(40e9, 0.002)]);
        assert_eq!(a.performance.bank().system().material_memory().len(), 8);
        assert!(a.performance.bank().system().material_relaxation_branches().iter()
            .all(|arm| arm.projection.len() == 2 * a.performance.bank().system().mode_count()));
        let p = pressure(&mut a, 1);
        assert_eq!(p, pressure(&mut b, 37));
        assert!(p.iter().any(|p| *p != 0.0));
        assert_eq!(a.performance.bank().system().material_memory(),
            b.performance.bank().system().material_memory());
        assert!(a.performance.bank().system().material_memory().iter().any(|z| *z != 0.0));
        let frame = a.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    });
}

#[test]
fn every_layer_supplies_its_own_supported_data_and_counts_toward_preflight() {
    assert!(input::Spec::parse(include_bytes!("../../relaxing-layered-strings.material-strings")).is_ok());
    scoped(|cx, gate| {
        for (omit, band) in [(Some("bending_relaxation_time"), 200_000.0), (None, 1000.0)] {
            assert!(selected_layered().compile(cx, gate, |name| {
                let sleeve = name == "layered-sleeve.fsmatpk";
                Ok(source_pack(name, 80e9, 0.001,
                    if sleeve { band } else { 200_000.0 }, if sleeve { omit } else { None }))
            }).is_err());
        }
        let mut s = selected_layered();
        s.strings[0].bending_modes = 128;
        s.strings[0].polarizations = Some(
            fs_couple::bowed_string::runtime::reciprocal::bank::polarization::TransverseDirections::try_new(
                [1.0, 0.0], [1.0, 0.0]).unwrap());
        assert!(s.string_memory_states().is_err(), "both physical layers must reserve history");
        let mut calls = 0;
        assert!(s.compile(cx, gate, |name| {
            calls += 1;
            Ok(source_pack(name, 80e9, 0.001, 200_000.0, None))
        }).is_err());
        assert_eq!(calls, 0);
    });
}
