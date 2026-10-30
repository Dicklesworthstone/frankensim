//! Synthetic source/geometry regressions, not a measured instrument qualification.
use super::*;
use crate::tests::bridge_interfaces::interface;

const WAVES: &str = "longitudinal-waves 0.2 0.01 64 100000 1e-11 1e-10 2 1 2";
const BEARING: &str = "axial-bearings 1\naxial-bearing 0 0.8 0 11 1 0 0 0 0 0.02 0 interface axial.fsintpk 0.0001 normal_stiffness normal_viscosity";

fn text() -> String {
    include_str!("../../longitudinal-cylinder-bow.material-strings")
        .replace("48000 4096 0.01 128", "48000 96 0.01 128")
        .replace(WAVES, &format!("{WAVES}\n{BEARING}"))
}
fn profile() -> Profile {
    Profile::parse(include_bytes!("../../arched.body")).unwrap()
}
fn loaded(text: &str, gate: &CancelGate, k: f64, outer: &str, body: &str, band: f64) -> input::Spec {
    let mut spec = input::Spec::parse(text.as_bytes()).unwrap();
    spec.axial_bearings.load(spec.temperature, spec.omega, gate,
        |_| Ok(interface(k, 1000.0, outer, body, band))).unwrap();
    spec
}
fn ready(cx: &Cx<'_>, gate: &CancelGate, k: f64) -> Ready {
    loaded(&text(), gate, k, "sleeve.fsmatpk", "body.fsmatpk", 200_000.0)
        .compile_with_body(cx, gate, Some(profile()), |name| Ok(pack(name, 1.0, None))).unwrap()
}

#[test]
fn sourced_axial_bearing_reacts_in_original_bow_body_gas_and_pcm_transaction() {
    scoped(|cx, gate| {
        let mut a = ready(cx, gate, 30e6);
        let source = a.performance.bank().source().bank();
        let longitudinal = source.longitudinal_waves().unwrap();
        assert_eq!(longitudinal.bearings().len(), 1);
        let binding = &longitudinal.bearings()[0];
        assert_eq!(binding.geometry.node, 11);
        assert_eq!(binding.geometry.string_axis, [1.0, 0.0, 0.0]);
        assert_eq!(binding.geometry.contact_arm_m, [0.0, 0.0, 0.02]);
        assert_eq!(binding.interface.stiffness_n_m(), 3000.0);
        assert_eq!(binding.interface.damping_n_s_m(), 0.1);
        assert_eq!(binding.interface.resolved().surface_materials(),
            &source.bridge_surface_materials(0).unwrap());
        assert!(source.cylinder_gas_loading().is_some());
        assert!(source.material_cavity().is_some());
        assert!(source.material_shell_image().is_some());
        let axial_component = longitudinal.strings()[0].axial_component();
        assert_eq!(a.performance.bank().system().components().iter()
            .map(|component| component.modes().len()).collect::<Vec<_>>(), vec![3, 2, 1, 12, 2]);
        let expected = pressure(&mut a, 7);
        assert!(expected.iter().any(|p| p.abs() > 1e-15));
        let system = a.performance.bank().system();
        assert!(system.components()[axial_component].states().iter()
            .any(|state| state.displacement_m_sqrt_kg != 0.0));
        assert_ne!(*system.coupled_frame().unwrap().connection_forces_n.last().unwrap(), 0.0);
        let frame = system.last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        assert!(frame.network_dissipation_j >= 0.0);

        let mut partitioned = ready(cx, gate, 30e6);
        assert_eq!(expected, pressure(&mut partitioned, 31));
        let mut changed = ready(cx, gate, 60e6);
        assert_eq!(a.spec.bytes, changed.spec.bytes);
        assert_eq!(a.performance.bank().cards()[0].fundamental_hz(),
            changed.performance.bank().cards()[0].fundamental_hz());
        assert_ne!(expected, pressure(&mut changed, 13));
    });
}

#[test]
fn axial_bearing_uses_actual_assembly_part_material_and_common_body_basis() {
    scoped(|cx, gate| {
        let text = text().replace("0 0.8 0 11", "0 0.8 1 4");
        let profile = || Profile::parse(include_str!("../../assembled.body")
            .replace("region right @body", "region right right-body.fsmatpk").as_bytes()).unwrap();
        let spec = loaded(&text, gate, 30e6, "sleeve.fsmatpk", "right-body.fsmatpk", 200_000.0);
        let mut ready = spec.compile_with_body(cx, gate, Some(profile()),
            |name| Ok(pack(name, 1.0, None))).unwrap();
        let source = ready.performance.bank().source().bank();
        let image = source.material_shell_assembly_image().unwrap();
        assert_eq!(image.physical_modes().len(), 2);
        assert_eq!(source.longitudinal_waves().unwrap().bearings()[0]
            .interface.resolved().surface_materials()[1].chemistry, "source-pack:right-body.fsmatpk");
        assert!(pressure(&mut ready, 31).iter().any(|p| p.abs() > 1e-15));
        let frame = ready.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        // A tied representative on part zero must not select that part's material.
        let wrong = loaded(&text, gate, 30e6, "sleeve.fsmatpk", "body.fsmatpk", 200_000.0);
        assert!(wrong.compile_with_body(cx, gate, Some(profile()),
            |name| Ok(pack(name, 1.0, None))).is_err());
    });
}

#[test]
fn missing_stale_mismatched_geometry_and_source_domains_refuse_before_output() {
    scoped(|cx, gate| {
        let mut calls = 0;
        let spec = input::Spec::parse(text().as_bytes()).unwrap();
        assert!(spec.compile_with_body(cx, gate, Some(profile()), |name| {
            calls += 1; Ok(pack(name, 1.0, None))
        }).is_err());
        assert_eq!(calls, 0);
        let mut stale = loaded(&text(), gate, 30e6, "sleeve.fsmatpk", "body.fsmatpk", 200_000.0);
        stale.temperature = 300.0;
        assert!(stale.compile_with_body(cx, gate, Some(profile()), |name| {
            calls += 1; Ok(pack(name, 1.0, None))
        }).is_err());
        assert_eq!(calls, 0);
        for (outer, body, band) in [
            ("wire.fsmatpk", "body.fsmatpk", 200_000.0),
            ("sleeve.fsmatpk", "other.fsmatpk", 200_000.0),
            ("sleeve.fsmatpk", "body.fsmatpk", 10_000.0),
        ] {
            let spec = loaded(&text(), gate, 30e6, outer, body, band);
            assert!(spec.compile_with_body(cx, gate, Some(profile()),
                |name| Ok(pack(name, 1.0, None))).is_err());
        }
        for input in [text().replace("0 0.8 0 11", "0 0.8 0 3071"),
            text().replace("0 0.8 0 11", "0 0.8 1 11")] {
            let spec = loaded(&input, gate, 30e6, "sleeve.fsmatpk", "body.fsmatpk", 200_000.0);
            assert!(spec.compile_with_body(cx, gate, Some(profile()),
                |name| Ok(pack(name, 1.0, None))).is_err());
        }
    });
}

#[test]
fn axial_bearing_grammar_preserves_physical_geometry_and_counts_complete_setup() {
    let original = text();
    let spec = input::Spec::parse(original.as_bytes()).unwrap();
    let without = input::Spec::parse(original.replace(&format!("{BEARING}\n"), "").as_bytes()).unwrap();
    let (modes, links, work) = spec.network_dimensions(true).unwrap();
    let (old_modes, old_links, _) = without.network_dimensions(true).unwrap();
    assert_eq!(modes, old_modes);
    assert_eq!(links, old_links+1);
    assert_eq!(work, modes*(links+1)*(links+1));
    for (old, new) in [
        ("axial-bearings 1", "axial-bearings 17"),
        ("0 0.8 0 11", "0 1 0 11"),
        ("0 11 1 0 0", "0 11 2 0 0"),
        ("0 0 0.02 0 interface", "0 0 NaN 0 interface"),
        ("0 interface axial.fsintpk", "0 3000 axial.fsintpk"),
        ("axial.fsintpk", "../axial.fsintpk"),
        ("normal_stiffness normal_viscosity", "normal_stiffness normal_stiffness"),
    ] {
        assert!(input::Spec::parse(original.replace(old, new).as_bytes()).is_err(), "accepted {new}");
    }
    assert!(input::Spec::parse(original.replace(&format!("{WAVES}\n"), "").as_bytes()).is_err());
}
