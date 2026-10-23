//! Physical neck geometry and the shared gas drive the actual reciprocal load.
use super::*;
use fs_material::gas_species::input::{GasRecord, ThermochemicalPhase, ThermochemicalSources};

#[path = "../../../src/render/schedule/reed/field_input/tests/support.rs"]
mod support;

const LEGACY: &str = "cavity neck 0.003 0.01 0.02 100 1";
const PHYSICAL: &str = "cavity physical-neck 0.003 0.01 0.02 flanged unflanged 1";

fn source() -> String {
    String::from_utf8(spec().bytes)
        .unwrap()
        .replace(LEGACY, PHYSICAL)
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-13 * expected.abs().max(1e-20),
        "actual {actual:e}, expected {expected:e}"
    );
}

// Independently reconstruct physical translation coefficients from acoustic
// geometry and the existing fluid laws, including the area-squared conversion.
fn check_coefficients(spec: &input::Spec) {
    let cavity = &spec.cavity.as_ref().unwrap().cavity;
    assert_eq!(cavity.loss(), CavityNeckLoss::RadiationAtResonance);
    assert_eq!(cavity.gas(), spec.gas_record.state());
    let g = cavity.geometry();
    let ends = fs_couple::cavity_modal::neck_ends(g);
    let area = core::f64::consts::PI * g.neck_radius_m * g.neck_radius_m;
    let length = g.neck_length_m
        + fs_phs::neck_end_correction(g.neck_radius_m, ends.inner).unwrap()
        + fs_phs::neck_end_correction(g.neck_radius_m, ends.outer).unwrap();
    let mass = spec.gas.density * area * length;
    let stiffness = spec.gas.density * spec.gas.sound_speed.powi(2) * area.powi(2) / g.volume_m3;
    let omega = (stiffness / mass).sqrt();
    let radiation = fs_phs::compact_radiation_impedance(
        spec.gas.density,
        spec.gas.sound_speed,
        g.neck_radius_m,
        omega,
        ends.outer,
    )
    .unwrap()
    .0;
    let wave = fs_duct::segment_wave(
        &spec.gas,
        g.neck_radius_m,
        omega,
        fs_duct::LossModel::Bessel,
    )
    .unwrap();
    let viscothermal =
        (fs_math::c64::C64::new(0.0, -1.0) * wave.wavenumber * wave.characteristic_impedance).re
            * g.neck_length_m;
    assert!(radiation > 0.0 && viscothermal > 0.0);
    close(cavity.neck_mass_kg(), mass);
    close(cavity.compression_stiffness_n_m(), stiffness);
    close(cavity.isolated_omega_rad_s(), omega);
    close(
        cavity.neck_resistance_n_s_m(),
        (radiation + viscothermal) * area.powi(2),
    );
}

#[test]
fn physical_neck_gas_and_each_end_condition_reach_shared_load_and_rendered_pressure() {
    scoped(|cx, gate| {
        let baseline = input::Spec::parse(source().as_bytes()).unwrap();
        check_coefficients(&baseline);
        let mut first = compile(cx, gate, baseline.clone(), 7800.0);
        let actual = first.performance.bank().source().bank();
        let neck = &actual.material_cavity().unwrap().cavity;
        assert_eq!(neck.loss(), CavityNeckLoss::RadiationAtResonance);
        assert_eq!(neck.gas(), actual.gas());
        // This is the real reciprocal neck damper, not an output-only decay.
        assert!(
            actual
                .body_connections()
                .iter()
                .any(|c| c.left.component == 1
                    && c.right.component == 1
                    && c.stiffness_n_m == 0.0
                    && c.damping_n_s_m == neck.neck_resistance_n_s_m())
        );
        let original_resistance = neck.neck_resistance_n_s_m();
        let original_mass = neck.neck_mass_kg();
        let original_omega = neck.isolated_omega_rad_s();
        let expected = pressure(&mut first, 1);
        assert!(expected.iter().any(|p| p.abs() > 1e-12));
        let mut replay = compile(cx, gate, baseline.clone(), 7800.0);
        assert_eq!(expected, pressure(&mut replay, 37));
        for (from, to) in [
            ("gas mixture 293.15 101325", "gas mixture 293.15 60000"),
            ("gas mixture 293.15 101325", "gas mixture 313.15 101325"),
            ("mole 2 N2 0.8 O2 0.2", "mole 3 N2 0.7 O2 0.2 CO2 0.1"),
            ("flanged unflanged", "unflanged unflanged"),
            ("flanged unflanged", "flanged flanged"),
        ] {
            let changed = input::Spec::parse(source().replace(from, to).as_bytes()).unwrap();
            check_coefficients(&changed);
            let new_neck = &changed.cavity.as_ref().unwrap().cavity;
            assert_ne!(new_neck.neck_mass_kg(), original_mass, "{to}");
            assert_ne!(
                new_neck.neck_resistance_n_s_m(),
                original_resistance,
                "{to}"
            );
            // Direction and amount of the gas: with the neck's end geometry
            // fixed, the isolated Helmholtz frequency is c sqrt(A / (V L_eff)),
            // so it moves by exactly the sound-speed ratio. The ratio is 1 at
            // 60 kPa (ideal gas), up at 313 K, and down with CO2.
            if !from.starts_with("flanged") {
                close(
                    new_neck.isolated_omega_rad_s() / original_omega,
                    changed.gas.sound_speed / baseline.gas.sound_speed,
                );
            }
            let mut rendered = compile(cx, gate, changed, 7800.0);
            for name in ["wire.fsmatpk", "body.fsmatpk"] {
                assert_eq!(
                    first.catalog[name].pack.content_hash(),
                    rendered.catalog[name].pack.content_hash()
                );
            }
            let output = pressure(&mut rendered, 128);
            assert!(
                expected
                    .iter()
                    .zip(output)
                    .any(|(a, b)| (a - b).abs() > 1e-12),
                "{to}"
            );
            let frame = rendered.performance.bank().system().last_frame().unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
    });
}

#[test]
fn source_backed_thermochemistry_changes_neck_loss_without_rewriting_performance() {
    let original = support::packs("cavity-atmosphere", 0.028, 3.5, 200.0);
    let replacement = support::packs("cavity-atmosphere", 0.028, 4.5, 200.0);
    let text = source().replace(
        "gas mixture 293.15 101325 mole 2 N2 0.8 O2 0.2",
        "gas thermochemical 293.15 101325 mole 1 cavity-atmosphere 1",
    );
    let a_bindings = [original.binding()];
    let b_bindings = [replacement.binding()];
    let a = input::Spec::parse_with_gas_sources(
        text.as_bytes(),
        Some(&ThermochemicalSources {
            species: &a_bindings,
            water_species: None,
            phase: ThermochemicalPhase::FrozenGas,
        }),
    )
    .unwrap();
    let b = input::Spec::parse_with_gas_sources(
        text.as_bytes(),
        Some(&ThermochemicalSources {
            species: &b_bindings,
            water_species: None,
            phase: ThermochemicalPhase::FrozenGas,
        }),
    )
    .unwrap();
    assert_eq!(a.bytes, b.bytes);
    assert_eq!(a.gas.density, b.gas.density);
    assert_eq!(
        a.cavity.as_ref().unwrap().cavity.neck_mass_kg(),
        b.cavity.as_ref().unwrap().cavity.neck_mass_kg()
    );
    assert_ne!(
        a.cavity.as_ref().unwrap().cavity.neck_resistance_n_s_m(),
        b.cavity.as_ref().unwrap().cavity.neck_resistance_n_s_m()
    );
    check_coefficients(&a);
    check_coefficients(&b);
    let GasRecord::Thermochemical(receipt) = &b.gas_record else {
        panic!("original gas receipt lost");
    };
    assert_eq!(
        receipt.source_packs()[0].nasa9_pack,
        replacement.nasa9.content_hash()
    );
    scoped(|cx, gate| {
        let mut first = compile(cx, gate, a, 7800.0);
        let mut second = compile(cx, gate, b, 7800.0);
        let p = pressure(&mut first, 37);
        let q = pressure(&mut second, 128);
        assert!(p.iter().zip(q).any(|(p, q)| (p - q).abs() > 1e-12));
    });
}

#[test]
fn legacy_resistance_and_flanged_ends_are_unchanged_and_physical_neck_refusals_propagate() {
    for text in [
        String::from_utf8(spec().bytes).unwrap(),
        String::from_utf8(spec().bytes)
            .unwrap()
            .replace("293.15 101325", "293.15 60000"),
    ] {
        let parsed = input::Spec::parse(text.as_bytes()).unwrap();
        let neck = &parsed.cavity.unwrap().cavity;
        assert_eq!(neck.loss(), CavityNeckLoss::Resistance(100.0));
        assert_eq!(neck.geometry().inner_end, fs_scenario::NeckEnd::Flanged);
        assert_eq!(neck.geometry().outer_end, fs_scenario::NeckEnd::Flanged);
        assert_eq!(
            neck.neck_resistance_n_s_m(),
            100.0 * neck.neck_area_m2().powi(2)
        );
    }
    for invalid in [
        "cavity physical-neck 0.003 0.01 0.02 automatic unflanged 1",
        "cavity physical-neck 0.003 0.01 0.02 flanged automatic 1",
        "cavity physical-neck 0.003 0.01 0.02 flanged unflanged 100 1",
        "cavity physical-neck 0.003 0.01 0.02 flanged unflanged",
        "cavity physical-neck 0 0.01 0.02 flanged unflanged 1",
        "cavity physical-neck 0.003 0 0.02 flanged unflanged 1",
        "cavity physical-neck 0.003 0.01 0 flanged unflanged 1",
        "cavity physical-neck 0.003 0.01 0.02 flanged unflanged NaN",
        // Positive geometry outside the original compact-radiation ka ceiling.
        "cavity physical-neck 1e-8 0.01 0.02 flanged unflanged 1",
    ] {
        assert!(
            input::Spec::parse(source().replace(PHYSICAL, invalid).as_bytes()).is_err(),
            "{invalid}"
        );
    }
    let fixture =
        input::Spec::parse(include_bytes!("../physical-cavity.material-strings")).unwrap();
    check_coefficients(&fixture);
}
