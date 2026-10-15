use super::*;

#[test]
fn resolved_geometry_keeps_source_features_and_material_authority() {
    for published in [false, true] {
        for (divisions, step) in [(8, 0.025), (16, 0.0125), (32, 0.02), (48, 0.025)] {
            let preset = build_resolved(divisions, step, published).unwrap();
            super::super::board_geometry::BoardGeometry::read(&preset.geometry).unwrap();
            assert_eq!(preset.geometry.lines().filter(|s| s.starts_with("bridge,")).count(), 88);
            assert!(preset.geometry.contains("source,mixed,"));
            assert!(preset.geometry.contains("support,clamped"));
            assert!(preset.geometry.contains("pretension,0"));
            assert!(preset.geometry.contains("damping,0.015"));
            for rib in 1..=17 {
                assert!(preset.obj.contains(&format!("g sugar_pine_rib_{rib}\n")));
            }
            assert_eq!(preset.geometry, build_resolved(divisions, step, published).unwrap().geometry);
        }
    }
}

#[test]
fn unrefined_geometry_is_unchanged_by_the_optional_constructor() {
    let old = build(8).unwrap();
    let resolved = build_resolved(8, 1.0, false).unwrap();
    assert_eq!(resolved.geometry.lines().skip(1).collect::<Vec<_>>(), old.geometry.lines().collect::<Vec<_>>());
    assert_eq!(resolved.obj.lines().skip(1).collect::<Vec<_>>(), old.obj.lines().collect::<Vec<_>>());
    assert!(build(33).is_err());
    assert!(build_with_rt0425_contacts(33).is_err());
}

#[test]
fn resolution_requires_explicit_finite_physical_inputs() {
    for count in [0, 3, 65, usize::MAX] {
        assert!(build_resolved(count, 0.02, true).is_err());
    }
    for step in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MIN_POSITIVE] {
        assert!(build_resolved(16, step, true).is_err());
    }
}
