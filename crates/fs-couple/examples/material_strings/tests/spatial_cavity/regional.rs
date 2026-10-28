//! Separate frozen acoustic regions through the original material body owner.
use super::*;
use fs_material::gas::GasState;

const REGIONAL: &str = include_str!("../../regional-cavity.material-strings");
const AIR: &str = "mixture 293.15 101325 mole 2 N2 0.8 O2 0.2";

fn with_interior(text: &str, gas: &str) -> String {
    text.replacen("\nbody ", &format!("\ncavity_gas {gas}\nbody "), 1)
}

fn interior(ready: &Ready)
    -> &fs_couple::bowed_string::runtime::reciprocal::bank::cavity::MaterialShellCavity
{
    ready.performance.bank().source().bank().spatial_cavity().unwrap()
}

#[test]
fn explicit_equal_cavity_gas_preserves_legacy_state_and_block_identity() {
    scoped(|cx, gate| {
        let mut inherited = ready(cx, gate, source(PERFORMANCE), BODY, 700.0);
        let explicit_spec = source(&with_interior(PERFORMANCE, AIR));
        assert_eq!(explicit_spec.cavity_gas(), &explicit_spec.gas);
        assert!(explicit_spec.cavity_gas_record.is_some());
        let mut explicit = ready(cx, gate, explicit_spec, BODY, 700.0);
        assert_eq!(uniform_stiffness(&inherited), uniform_stiffness(&explicit));
        let a = pressure(&mut inherited, 1);
        assert!(a.iter().any(|p| p.abs() > 1e-14));
        assert_eq!(a, pressure(&mut explicit, 37));
        for (a, b) in inherited.performance.bank().system().components().iter()
            .zip(explicit.performance.bank().system().components()) {
            assert_eq!(a.states(), b.states());
        }
    });
}

#[test]
fn independent_interior_composition_reaches_compression_and_accepted_string_pressure() {
    scoped(|cx, gate| {
        let mut air = ready(cx, gate, source(PERFORMANCE), BODY, 700.0);
        let mut helium = ready(cx, gate, source(REGIONAL), BODY, 700.0);
        let outer = helium.performance.bank().source().bank().gas();
        assert_eq!(outer, &air.spec.gas);
        let inside = helium.performance.bank().source().bank().spatial_cavity().unwrap().gas();
        assert_eq!(inside, helium.spec.cavity_gas());
        assert_ne!(inside, outer);
        // Uniform compression is K = rho*c^2/V times the SAME signed wall
        // volume projection. Composition changes that real spring before audio.
        let bulk = |g: &GasState| g.density*g.sound_speed*g.sound_speed;
        let expected = bulk(inside)/bulk(air.spec.cavity_gas());
        assert!((uniform_stiffness(&helium)/uniform_stiffness(&air)-expected).abs()
            < 2e-10*expected);
        assert!((expected-1.0).abs() > 0.05);
        let a = pressure(&mut air, 1);
        let b = pressure(&mut helium, 37);
        assert!(a.iter().zip(&b).any(|(a,b)| (a-b).abs() > 1e-14));
        assert_ne!(air.performance.bank().system().components()[0].states(),
            helium.performance.bank().system().components()[0].states());
        let frame = helium.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    });
}

#[test]
fn interior_temperature_moves_pressure_modes_while_exterior_changes_keep_their_original_medium() {
    scoped(|cx, gate| {
        let air = ready(cx, gate, source(&with_interior(PERFORMANCE, AIR)), BODY, 700.0);
        let warm = ready(cx, gate, source(&with_interior(PERFORMANCE,
            "mixture 313.15 101325 mole 2 N2 0.8 O2 0.2")), BODY, 700.0);
        assert_eq!(air.spec.gas, warm.spec.gas);
        let speed_ratio = warm.spec.cavity_gas().sound_speed/air.spec.cavity_gas().sound_speed;
        assert!(speed_ratio > 1.0);
        let frequency_ratio = interior(&warm).modal_report().frequency_hz[1]
            / interior(&air).modal_report().frequency_hz[1];
        assert!((frequency_ratio-speed_ratio).abs() < 2e-8*speed_ratio);
        // Change only the outside. The interior FE source and its eigenproblem
        // stay exactly the same; the declared atmosphere still reaches drag and
        // exterior pressure through the bank's original gas.
        let text = with_interior(PERFORMANCE, AIR).replacen(
            &format!("\ngas {AIR}"), "\ngas mixture 293.15 101325 mole 1 He 1", 1);
        let outside = ready(cx, gate, source(&text), BODY, 700.0);
        assert_ne!(outside.spec.gas, air.spec.gas);
        assert_eq!(interior(&outside).gas(), interior(&air).gas());
        assert_eq!(interior(&outside).modal_report().frequency_hz,
            interior(&air).modal_report().frequency_hz);
        assert_eq!(uniform_stiffness(&outside), uniform_stiffness(&air));
    });
}

#[test]
fn unsupported_reference_pressure_unused_region_and_compact_neck_refuse_before_material_io() {
    let unequal = REGIONAL.replace("cavity_gas mixture 293.15 101325",
        "cavity_gas mixture 293.15 90000");
    assert!(input::Spec::parse(unequal.as_bytes()).is_err());
    let legacy = String::from_utf8(spec().bytes).unwrap();
    assert!(input::Spec::parse(with_interior(&legacy,
        "mixture 293.15 101325 mole 1 He 1").as_bytes()).is_err());
    scoped(|cx, gate| {
        let mut reads = 0;
        assert!(source(REGIONAL).compile::<NormalizedPack>(cx, gate, |_| {
            reads += 1;
            Err("an unconsumed cavity atmosphere must refuse before material I/O".into())
        }).is_err());
        assert_eq!(reads, 0);
    });
}
