//! Material-only winding swaps must change the actual reciprocal pressure path.
use super::*;

fn wound_spec() -> input::Spec {
    input::Spec::parse(include_str!("../wound-pluck.material-strings")
        .replace("48000 4096 0.01 128", "48000 256 0.01 128").as_bytes()).unwrap()
}
fn wound_ready(cx: &Cx<'_>, gate: &CancelGate, spec: input::Spec, density: f64, band: f64) -> Result<Ready> {
    spec.compile(cx, gate, |name| Ok(if name == "wrap.fsmatpk" {
        // No wrap bending viscosity is needed or silently borrowed from core.
        pack(name, density, band, Some(ETA))
    } else { pack(name, 7800.0, 200_000.0, None) }))
}
fn specimen(ready: &Ready) -> &fs_couple::material_assembly::wound::ResolvedWoundStringSpecimen {
    let source = &ready.performance.bank().source().bank().strings()[0].source;
    let MaterialStringSource::Wound(specimen) = source else { panic!("original winding source retained") };
    specimen
}

#[test]
fn winding_material_swap_reaches_mass_modes_gas_loading_and_real_pressure() {
    scoped(|cx, gate| {
        let mut light = wound_ready(cx, gate, wound_spec(), 8900.0, 200_000.0).unwrap();
        let mut dense = wound_ready(cx, gate, wound_spec(), 17800.0, 200_000.0).unwrap();
        assert_eq!(light.spec.bytes, dense.spec.bytes);
        assert_eq!(specimen(&light).section().core(), specimen(&dense).section().core());
        assert_ne!(specimen(&light).section().winding_material().identity(), specimen(&dense).section().winding_material().identity());
        assert_eq!(specimen(&light).string().tension_n, 20.0);
        assert_eq!(specimen(&dense).string().tension_n, 20.0);
        assert_eq!(specimen(&light).string().bending_stiffness_n_m2,
            specimen(&dense).string().bending_stiffness_n_m2);
        assert!(dense.performance.bank().cards()[0].fundamental_hz() < light.performance.bank().cards()[0].fundamental_hz());
        let original = pressure(&mut light, 37);
        let changed = pressure(&mut dense, 1);
        assert!(original.iter().zip(&changed).any(|(a, b)| (a - b).abs() > 1e-12));
        assert!(original.iter().any(|p| p.abs() > 1e-10));
        let frame = dense.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        let mut replay = wound_ready(cx, gate, wound_spec(), 8900.0, 200_000.0).unwrap();
        assert_eq!(pressure(&mut replay, 128), original);
    });
}

#[test]
fn wound_extension_uses_spring_plus_core_and_gas_radius_is_the_finished_radius() {
    scoped(|cx, gate| {
        let text = String::from_utf8(wound_spec().bytes).unwrap().replace("tension 20", "extension 0.4998 0.001");
        let ready = wound_ready(cx, gate, input::Spec::parse(text.as_bytes()).unwrap(), 8900.0, 200_000.0).unwrap();
        let source = specimen(&ready);
        let strain = (0.5 - 0.4998) / 0.4998;
        assert_eq!(source.string().tension_n, source.section().axial_stiffness_n() * strain);
        assert!(source.section().axial_stiffness_n() > source.section().core().elastic().axial_stiffness_n());
        assert_eq!(source.string().width_m, 2.0 * source.section().outer_radius_m());
        // Source and mechanical state survive the existing cylinder gas owner.
        let gas_text = text.replace("cavity neck", "cylinder-gas-loading 0.001 30 4\ncavity neck");
        let loaded = wound_ready(cx, gate, input::Spec::parse(gas_text.as_bytes()).unwrap(), 8900.0, 200_000.0).unwrap();
        assert_eq!(specimen(&loaded).section(), source.section());
    });
}

#[test]
fn winding_geometry_source_band_and_unmodeled_torsion_refuse() {
    scoped(|cx, gate| {
        assert!(wound_ready(cx, gate, wound_spec(), 8900.0, 1000.0).is_err());
        let text = String::from_utf8(wound_spec().bytes).unwrap();
        for changed in [text.replace("0.00005 0.000101", "0.00005 0.0001"),
            text.replace("20 4 0 0.8", "20 4 2 0.8")] {
            let spec = input::Spec::parse(changed.as_bytes()).unwrap();
            assert!(wound_ready(cx, gate, spec, 8900.0, 200_000.0).is_err());
        }
        assert!(input::Spec::parse(text.replace("0.00005 0.000101", "0.00004 0.000101").as_bytes()).is_err());
    });
}
