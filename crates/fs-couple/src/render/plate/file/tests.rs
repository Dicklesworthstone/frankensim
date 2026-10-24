use super::*;
use crate::thin_plate::CompactBody;

const SOURCE: &[u8] = include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/examples/plate-mesh.performance"));

fn independent_chart() -> (PlateChart, PlateChartRadiation) {
    let mut nodes = Vec::new();
    for y in [0.0, 0.075, 0.15, 0.225, 0.3] {
        for x in [0.0, 0.08, 0.16, 0.24, 0.32, 0.4] { nodes.push((x,y)); }
    }
    let sections = [
        PlateSection::orthotropic_plane_stress_at_angle(7e9, 7e9, 0.3, 2692307692.307692, 0.003, 700.0, 0.0).unwrap(),
        PlateSection::orthotropic_plane_stress_at_angle(11e9, 1.4e9, 0.3, 0.8e9, 0.004, 600.0, 0.35).unwrap(),
    ];
    let mut triangles = Vec::new();
    let mut assigned = Vec::new();
    for j in 0..4 {
        for i in 0..5 {
            let a = j*6+i;
            triangles.extend([[a,a+1,a+7], [a,a+7,a+6]]);
            assigned.extend([sections[usize::from(i>=3)]; 2]);
        }
    }
    let supports = (0..30).filter(|n| n%6==0 || n%6==5 || *n<6 || *n>=24).collect();
    let mesh = PlateMesh::from_unstructured(nodes, triangles).unwrap();
    let chart = PlateChart::with_boundary_and_regions(mesh, sections[0], supports, Vec::new()).unwrap()
        .with_element_sections(assigned).unwrap();
    let mut weights = vec![0.0; 30];
    weights[14] = 0.7;
    weights[15] = 0.3;
    (chart, PlateChartRadiation { assembly: AssemblyOptions {
        support: EdgeSupport::SimplySupported, pretension: 0.0,
    }, eigenvalue_window: ((core::f64::consts::TAU*5.0).powi(2), (core::f64::consts::TAU*700.0).powi(2)),
        n_modes: 3, damping_ratio: 0.02, unit_force_weights: weights })
}
fn source_air() -> GasState { GasState::try_new_moist_air(293.15, 101_325.0, 0.0).unwrap() }
/// Independent retarded reference: output i hears t_{i+1} - r/c, which lies
/// (1 - f) dt into step i - m for a 1 m delay of m + f samples.
fn direct(mut bodies: Vec<CompactBody>, gas: &GasState) -> Vec<f64> {
    let delay = 1.0 * 48_000.0 / gas.sound_speed;
    let whole = delay.floor() as usize;
    let elapsed = (1.0 - (delay - delay.floor())) * 48_000.0_f64.recip();
    let mut starts = Vec::with_capacity(4801);
    (0..4801).map(|i| {
        let force = if i < 37 { 1.0 } else if i < 1701 { 0.0 } else if i < 1738 { -0.4 } else { 0.0 };
        starts.push((bodies.clone(), force));
        for body in &mut bodies { body.drive(force*body.drive_participation, 1.0/48_000.0).unwrap(); }
        if i < whole { return 0.0; }
        let (start, force) = &starts[i - whole];
        let mut p = 0.0;
        for body in start {
            let mut probe = body.clone();
            let acc = probe.drive(force*probe.drive_participation, elapsed).unwrap();
            p += probe.radiate(acc, gas.density, 1.0);
        }
        p
    }).collect()
}

#[test]
fn source_triangle_sections_and_force_footprint_are_the_real_operator_inputs() {
    let parsed = Parsed::read(SOURCE, 37, None).unwrap();
    let (chart, options) = independent_chart();
    assert_eq!(parsed.chart.mesh.nodes, chart.mesh.nodes);
    assert_eq!(parsed.chart.mesh.tris, chart.mesh.tris);
    assert_eq!(parsed.chart.boundary_nodes, chart.boundary_nodes);
    assert_eq!(parsed.options.unit_force_weights, options.unit_force_weights);
    let fs_plate::PlateSectionField::PerElement(actual) = parsed.chart.section_field() else { panic!("missing assignment"); };
    let fs_plate::PlateSectionField::PerElement(expected) = chart.section_field() else { panic!("missing assignment"); };
    assert_eq!(actual.len(), 40);
    for (a,b) in actual.iter().zip(expected) {
        assert_eq!(a.d, b.d);
        assert_eq!(a.thickness, b.thickness);
        assert_eq!(a.density, b.density);
    }
    assert_ne!(actual[0].d, actual[6].d, "regional material must not become decoration");
    assert_eq!(parsed.events[0].sample, 37);
    assert_eq!(parsed.config.gas, source_air());
    assert_eq!(parsed.gas_record, GasRecord::Legacy(source_air()));
    assert_eq!(parsed.config.listener_m, 1.0);
}

#[test]
fn public_loader_renders_geometry_derived_modes_like_direct_chart_mechanics() {
    let (chart, options) = independent_chart();
    let expected = direct(certified_chart_radiators(&chart, &options).unwrap(), &source_air());
    let onset = (48_000.0 / source_air().sound_speed).floor() as usize;
    assert!(expected[..onset].iter().all(|p| *p == 0.0), "silence before r/c");
    assert!(expected[onset + 37..onset + 1701].iter().any(|p| p.abs() > 1e-9), "real unloaded ringdown");
    for block in [1,37,571] {
        let performance = PlatePerformance::from_bytes(SOURCE, block).unwrap();
        let info = performance.info();
        assert_eq!((info.nodes,info.triangles,info.sections), (30,40,2));
        assert_eq!(info.samples, 4801);
        assert_eq!(info.force_events, 3);
        assert_eq!(info.gas, source_air());
        assert_eq!(info.propagation_delay_s, 1.0 / source_air().sound_speed);
        assert!(info.retained_modes > 0 && info.retained_modes <= 3);
        let mut renderer = performance.into_renderer();
        let mut actual = vec![0.0; 4801];
        for chunk in actual.chunks_mut(block) { renderer.block(chunk).unwrap(); }
        assert_eq!(actual.iter().map(|p| p.to_bits()).collect::<Vec<_>>(),
            expected.iter().map(|p| p.to_bits()).collect::<Vec<_>>());
        assert!(renderer.pending_controls().is_empty());
    }
}

#[test]
fn every_record_is_required_and_ignored_suffixes_and_extra_fields_refuse() {
    let text = std::str::from_utf8(SOURCE).unwrap();
    for (offset, _) in text.match_indices('\n').take(text.lines().count()-1) {
        assert!(Parsed::read(&SOURCE[..offset], 64, None).is_err(), "truncation {offset}");
    }
    for changed in [text.to_string()+"ignored\n", text.replace("observer 1.0", "observer 1.0 extra"),
        text.replace("ambient 293.15 101325 0", "ambient 293.15 101325 0 extra"),
        text.replace("performance-v2", "performance-v3")] {
        assert!(Parsed::read(changed.as_bytes(), 64, None).is_err());
    }
    // Version 1 authored a bare density in the observer row; it is not air by default.
    let v1 = text.replace("performance-v2", "performance-v1")
        .replace("ambient 293.15 101325 0\nobserver 1.0", "observer 1.2 1.0");
    assert!(Parsed::read(v1.as_bytes(), 64, None).is_err());
    assert!(Parsed::read(&[255], 64, None).is_err());
    assert!(Parsed::read(&vec![b'x'; MAX_PLATE_PERFORMANCE_BYTES+1], 64, None).is_err());
}

#[test]
fn bad_geometry_sections_weights_times_and_budgets_refuse_before_reduction() {
    let text = std::str::from_utf8(SOURCE).unwrap();
    for (old,new) in [
        ("nodes 30", "nodes 18446744073709551615"),
        ("sections 2", "sections 65"),
        ("section 0.003", "section -0.003"),
        ("triangle 0 1 7 0", "triangle 0 7 1 0"),
        ("triangle 0 1 7 0", "triangle 0 1 700 0"),
        ("triangle 0 1 7 0", "triangle 0 1 7 9"),
        ("support 1\n", "support 0\n"),
        ("weight 15 0.3", "weight 14 0.3"),
        ("weight 15 0.3", "weight 15 0.4"),
        ("initial_force_n 1.0", "initial_force_n NaN"),
        ("force 37 0.0", "force 4801 0.0"),
        ("force 37 0.0", "force 37 101.0"),
        ("1000000.0 32", "1000000.0 2"),
        ("700.0 3", "24000.0 3"),
        ("audio 48000 4801", "audio 0 4801"),
    ] {
        assert!(text.contains(old));
        assert!(Parsed::read(text.replace(old,new).as_bytes(), 64, None).is_err(), "{old} -> {new}");
    }
    assert!(Parsed::read(SOURCE,0,None).is_err());
    assert!(Parsed::read(SOURCE,65_537,None).is_err());
}

#[test]
fn an_unresolvable_gas_record_refuses_by_name_and_never_becomes_air() {
    let text = std::str::from_utf8(SOURCE).unwrap();
    let ambient = "ambient 293.15 101325 0";
    for (record, material) in [
        ("ambient 293.15 101325", false),
        ("ambient NaN 101325 0", false),
        ("ambient 293.15 101325 0.5 0.5", false),
        ("ambient mixture 293.15 101325 mole 1 Xe 1", true),
        ("ambient mixture 293.15 101325 mole 2 N2 0.5 O2 0.4", true),
        ("ambient thermochemical 293.15 101325 mole 1 N2 1", false),
        ("ambient 293.15 0 0", true),
        ("observer 1.0", false),
    ] {
        let changed = text.replace(ambient, record);
        let error = Parsed::read(changed.as_bytes(), 64, None).err().unwrap_or_else(|| panic!("{record} admitted"));
        assert_eq!(matches!(error, PlatePerformanceError::Gas { line: 3, .. }), material, "{record}: {error}");
        if material { assert!(error.to_string().contains("line 3: ambient gas refused"), "{error}"); }
        assert!(PlatePerformance::from_bytes(changed.as_bytes(), 64).is_err());
    }
    let co2 = text.replace(ambient, "ambient mixture 293.15 101325 mole 1 CO2 1");
    let parsed = Parsed::read(co2.as_bytes(), 64, None).unwrap();
    assert!(matches!(parsed.gas_record, GasRecord::Mixture(..)));
    assert!(parsed.config.gas.sound_speed < 0.8 * source_air().sound_speed);
}

#[test]
fn signed_force_events_retain_source_order_for_same_sample_assignments() {
    let text = std::str::from_utf8(SOURCE).unwrap().replace("force 1738 0.0", "force 37 -0.2");
    let performance = PlatePerformance::from_bytes(text.as_bytes(), 64).unwrap();
    let mut renderer = performance.into_renderer();
    let mut initial = [0.0; 37];
    renderer.block(&mut initial).unwrap();
    renderer.block(&mut [0.0; 1]).unwrap();
    let applied = renderer.applied_controls();
    assert_eq!(applied.len(), 2);
    assert_eq!(applied[0].delta, ControlDelta::SetPlateForce { voice:0, force_n:0.0 });
    assert_eq!(applied[1].delta, ControlDelta::SetPlateForce { voice:0, force_n:-0.2 });
    assert_eq!(renderer.pending_controls()[0].sample, 1701);
}
