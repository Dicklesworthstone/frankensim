//! Fluid memory composes with regional viscosity, torque, stops and preload.
use super::*;

const WIRE: &str = "string wire.fsmatpk 0.5 0.0003 tension 20 1 1 0.8 3000 0.1\nsection bonded 0.0002 1\nlayer sleeve.fsmatpk 0.0003\npolarizations two 0.6 0.8 1 0";
const FLOW: &str = "cylinder-gas-loading 0.001 30 6";
const STOP: &str = "frankensim-string-stops-v1\nstop support 0 0.43 positive bridge:0 0 5000 1 0.2 10 0.001\n";

fn instrument(flow: bool) -> input::Spec {
    let input = include_str!("../cylinder-flow.material-strings")
        .replace("48000 4096 0.01 128", "48000 128 0.01 128")
        .replace("strings 1", "strings 2")
        .replace(WIRE, &format!("{WIRE}\n{}", WIRE.replace("tension 20", "tension 22")))
        .replace("cylinder-gas-loading 0.001 30 12",
            &format!("rotational-bridges 1\nrotational-bridge 0 0.8 0 1 0.001 1e-7 0\n{FLOW}"));
    let input = if flow { input } else { input.replace(&format!("{FLOW}\n"), "") };
    input::Spec::parse(input.as_bytes()).unwrap()
}

fn ready(cx: &Cx<'_>, gate: &CancelGate, flow: bool, right_loss: f64) -> Ready {
    let profile = include_str!("../torsional-bridge.body")
        .replace("region right @body", "region right right.fsmatpk");
    instrument(flow).compile_with_profiles(cx, gate,
        Some(Profile::parse(profile.as_bytes()).unwrap()),
        Some(crate::stops::Profile::parse(STOP.as_bytes()).unwrap()),
        |name| Ok(if name == "right.fsmatpk" { pack("body.fsmatpk", right_loss, None) }
            else { pack(name, 1.0, None) })).unwrap()
}

#[test]
fn regional_body_torque_and_stop_keep_the_same_sources_when_fluid_coordinates_are_appended() {
    scoped(|cx, gate| {
        let mut loaded = ready(cx, gate, true, 1.0);
        let mut plain = ready(cx, gate, false, 1.0);
        assert_eq!(loaded.spec.network_dimensions(true).unwrap(), (45, 43, 45*44*44));
        assert_eq!(loaded.rotational_bridge_baselines_m, vec![Some(0.2)]);
        assert_eq!(loaded.stop_profile.as_ref().unwrap().bytes, STOP.as_bytes());
        assert_eq!(loaded.performance.bank().system().obstacle_count(), 1);
        let a = loaded.performance.bank().source().bank();
        let b = plain.performance.bank().source().bank();
        assert_eq!(a.strings().len(), 2);
        assert_eq!(a.body_connections(), b.body_connections());
        assert_eq!(a.rotational_body_connections(), b.rotational_body_connections());
        assert_eq!(a.material_chart_viscosity().unwrap().damping().matrix_per_s(),
            b.material_chart_viscosity().unwrap().damping().matrix_per_s());
        assert_eq!(a.material_chart().unwrap().specimen.regions().len(), 2);
        assert_eq!(a.material_cavity().unwrap().cavity.gas(), a.gas());
        assert_eq!(a.cylinder_gas_loading().unwrap().strings().len(), 2);
        assert_eq!(loaded.performance.bank().system().components()[2].modes(),
            plain.performance.bank().system().components()[2].modes());
        let contact = loaded.performance.bank().system().obstacle(0).unwrap();
        assert_eq!(contact.right.component, 2);
        assert_eq!(contact.right.shapes, a.strings()[0].connections[0].body.shapes);
        assert_eq!(contact.left.shapes.len(), 3);
        assert_eq!(contact.left.shapes[2], 0.0, "a radial stop invented torsional loading");
        let components = loaded.performance.bank().system().components();
        assert_eq!(components.iter().map(|c| c.modes().len()).collect::<Vec<_>>(),
            vec![3, 3, 2, 1, 18, 18]);
        assert!(components[4..].iter().flat_map(|c| c.states()).all(|s|
            s.displacement_m_sqrt_kg == 0.0 && s.velocity_m_sqrt_kg_per_s == 0.0));
        assert!(loaded.performance.preload_energy_j().unwrap() > 0.0);
        let p = pressure(&mut loaded, 37);
        let q = pressure(&mut plain, 128);
        assert!(p.iter().any(|p| p.abs() > 1e-14));
        assert!(p.iter().zip(q).any(|(p, q)| (p-q).abs() > 1e-17));
        assert!(loaded.performance.bank().system().components()[4].states().iter()
            .any(|s| s.velocity_m_sqrt_kg_per_s != 0.0));
        let frame = loaded.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        let mut replay = ready(cx, gate, true, 1.0);
        assert_eq!(p, pressure(&mut replay, 1));
    });
}

#[test]
fn one_regions_material_loss_still_changes_actual_pressure_with_gas_memory_and_contact() {
    scoped(|cx, gate| {
        let mut original = ready(cx, gate, true, 1.0);
        let mut changed = ready(cx, gate, true, 8.0);
        assert_eq!(original.spec.bytes, changed.spec.bytes);
        assert_eq!(original.catalog["wire.fsmatpk"].pack.content_hash(),
            changed.catalog["wire.fsmatpk"].pack.content_hash());
        assert_eq!(original.catalog["body.fsmatpk"].pack.content_hash(),
            changed.catalog["body.fsmatpk"].pack.content_hash());
        assert_ne!(original.catalog["right.fsmatpk"].pack.content_hash(),
            changed.catalog["right.fsmatpk"].pack.content_hash());
        let a = original.performance.bank().source().bank();
        let b = changed.performance.bank().source().bank();
        for (a, b) in a.cylinder_gas_loading().unwrap().strings().iter()
            .zip(b.cylinder_gas_loading().unwrap().strings()) {
            assert_eq!(a.flow(), b.flow());
        }
        assert_ne!(a.material_chart_viscosity().unwrap().damping().matrix_per_s(),
            b.material_chart_viscosity().unwrap().damping().matrix_per_s());
        assert_eq!(original.performance.bank().system().components()[2].modes(),
            changed.performance.bank().system().components()[2].modes());
        let p = pressure(&mut original, 128);
        let q = pressure(&mut changed, 37);
        assert!(p.iter().zip(q).any(|(p, q)| (p-q).abs() > 1e-15));
        let frame = changed.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    });
}
