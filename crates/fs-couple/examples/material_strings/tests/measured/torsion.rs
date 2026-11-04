//! A measured whole string's own measured torsion, selected by `torsion-measured` (README.md, "Measured torsion").
//! The bank tests use the parent module's manufactured adapter pack with a fixture row bound to it; that row carries
//! the cited D3 numbers but is a test fixture, not a measurement of the fixture string.
use super::*;

const D3: &str = "mansour-2016-kaplan-cello-d3-open";
/// The measured-string refusal, word for word (measured.rs).
const REFUSAL: &str = "measured whole strings supply transverse mass and bending stiffness, not torsion, thermal, axial, \
    or intrinsic-loss laws from outside their own source rows";

/// The parent fixture card with three twist modes.
fn twisted() -> String {
    text().replace("116.3 4 0 0.8 ", "116.3 4 3 0.8 ")
}
/// The cited D3 numbers on a row bound to the fixture specimen at its 0.7 m x 0.8 nut-to-bridge span.
fn fixture_row() -> input::MeasuredTorsionRow {
    input::MeasuredTorsionRow {
        id: "fixture-d3-numbers".into(),
        pack: "measured-frontend-fixture".into(),
        linear_density: "linear_density_dominant".into(),
        speaking_length_m: 0.56,
        fundamental_hz: 758.0,
        impedance_kg_s: 1.8,
        quality_factor: 34.0,
        source: "test fixture carrying the D3 row's numbers, not a measurement".into(),
    }
}
fn with_row(card: &str, row: input::MeasuredTorsionRow) -> input::Spec {
    let mut spec = input::Spec::parse(card.as_bytes()).unwrap();
    spec.strings[0].measured_torsion = Some(row);
    spec
}
fn selected(card: &str, record: &str) -> String {
    card.replace(PREPARED, &format!("{record}\n{PREPARED}"))
}

#[test]
fn the_table_holds_the_one_cited_row() {
    let rows = input::measured_torsion_rows().unwrap();
    assert_eq!(rows.len(), 1);
    let row = &rows[0];
    assert_eq!(
        (row.id.as_str(), row.pack.as_str(), row.linear_density.as_str()),
        (D3, "cello-kaplan-ks510-mansour-2016-percival-2013-scratch", "linear_density_ks510_d3")
    );
    assert_eq!(
        (row.speaking_length_m, row.fundamental_hz, row.impedance_kg_s, row.quality_factor),
        (0.69, 758.0, 1.8, 34.0)
    );
    assert!(row.source.contains("\"For the open cello D3 string studied here the torsional fundamental frequency is \
        taken to be 758 Hz, the characteristic torsional impedance is 1.8 kg/s and a constant Q of 34 is assigned to all \
        torsional modes.\""));
    assert!(row.source.contains("Woodhouse & Loach") && row.source.contains("85(5):734-740, 1999"));
    assert!(row.source.contains("Part 2 names neither the string's make nor its length"));
}

#[test]
fn without_a_row_twist_on_a_measured_string_still_refuses_word_for_word() {
    scoped(|cx, gate| {
        let error = ready(cx, gate, input::Spec::parse(twisted().as_bytes()).unwrap(), 1.0).err().unwrap().to_string();
        assert!(error.contains(REFUSAL), "{error}");
        // The row admits twist only: a measured string's other refused laws still refuse with the same words.
        let mut lossy = with_row(&twisted(), fixture_row());
        lossy.strings[0].tension_loss = Some(input::tension_loss_rows().unwrap().remove(0));
        let error = ready(cx, gate, lossy, 1.0).err().unwrap().to_string();
        assert!(error.contains(REFUSAL), "{error}");
        // A row with no twist modes refuses rather than build an empty bank.
        let error = ready(cx, gate, with_row(&text(), fixture_row()), 1.0).err().unwrap().to_string();
        assert!(error.contains("needs the string row's twist mode count"), "{error}");
    });
}

#[test]
fn a_row_builds_the_twist_bank_at_the_measured_fundamental_impedance_and_q() {
    scoped(|cx, gate| {
        let mut built = ready(cx, gate, with_row(&twisted(), fixture_row()), 1.0).unwrap();
        let bank = built.performance.bank();
        let MaterialStringSource::WithPolarizations(planes) = &bank.source().bank().strings()[0].source else {
            panic!("the two declared planes must survive");
        };
        let MaterialStringSource::WithTorsion(twist) = planes.source() else { panic!("the measured twist must exist") };
        assert!(twist.section().is_none());
        assert!(matches!(twist.bending(), MaterialStringSource::Measured(_)));
        assert_eq!(twist.mode_count(), 3);
        let measured = twist.measured().unwrap();
        assert_eq!(
            (measured.fundamental_hz(), measured.impedance_kg_s(), measured.quality_factor(), measured.support_fraction()),
            (758.0, 1.8, 34.0, 0.8)
        );
        // Two planes of the four bending modes, then the three twist modes.
        let modes = bank.system().components()[0].modes();
        assert_eq!(modes.len(), 11);
        assert_eq!(modes[..4], modes[4..8]);
        for (index, mode) in modes[8..].iter().enumerate() {
            let n = (index + 1) as f64;
            assert_eq!(mode.angular_frequency_rad_s, core::f64::consts::TAU * 758.0 * n);
            // The measured Q on every mode, not the shear law's Q falling as 1/n.
            assert_eq!(mode.damping_ratio, 0.5 / 34.0);
            assert_eq!(mode.pressure_per_modal_velocity, fs_math::c64::C64::ZERO);
        }
        // The impedance the bow actually sees. Its surface map of twist mode n at station x is
        // sin(n pi x / s) / sqrt(m L / 2) on the nut-to-bridge span L = s x 0.7 m, with c = omega_1 L / pi,
        // so Z = m c = 2 omega_1 / (pi A^2) for the amplitude A, whatever L is.
        let shapes = &bank.system().port(0).unwrap().left.shapes;
        assert_eq!(shapes.len(), 11);
        let omega_1 = modes[8].angular_frequency_rad_s;
        for n in 1..=3 {
            let amplitude = shapes[7 + n] / (n as f64 * core::f64::consts::PI * 0.2 / 0.8).sin();
            let impedance = 2.0 * omega_1 / (core::f64::consts::PI * amplitude * amplitude);
            assert!((impedance / 1.8 - 1.0).abs() < 1e-12, "mode {n}: Z {impedance}");
        }
        let p = pressure(&mut built, 64);
        assert!(p.iter().all(|value| value.is_finite()) && p.iter().any(|value| value.abs() > 1e-15));
        let frame = built.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    });
}

#[test]
fn a_row_never_reaches_another_string_a_wrong_span_or_a_stopped_string() {
    scoped(|cx, gate| {
        // The cited D3 row, selected by its id, on the fixture string: another specimen refuses.
        let spec = input::Spec::parse(selected(&twisted(), &format!("torsion-measured {D3}")).as_bytes()).unwrap();
        assert_eq!(spec.strings[0].measured_torsion.as_ref().map(|row| row.id.as_str()), Some(D3));
        let error = ready(cx, gate, spec, 1.0).err().unwrap().to_string();
        assert!(error.contains("was measured on cello-kaplan-ks510-mansour-2016-percival-2013-scratch"), "{error}");
        // The fixture specimen under another line-density property refuses too.
        let mut other = fixture_row();
        other.linear_density = "linear_density_other".into();
        let error = ready(cx, gate, with_row(&twisted(), other), 1.0).err().unwrap().to_string();
        assert!(error.contains("was measured on"), "{error}");
        // A fundamental measured on another speaking length refuses.
        let mut longer = fixture_row();
        longer.speaking_length_m = 0.7;
        let error = ready(cx, gate, with_row(&twisted(), longer), 1.0).err().unwrap().to_string();
        assert!(error.contains("speaking length"), "{error}");
        // The row is the open string's: a stop profile refuses before any stop is validated.
        let stops = crate::stops::Profile::parse(include_bytes!("../../contact-bow.stops")).unwrap();
        let error = with_row(&twisted(), fixture_row()).compile_with_profiles(cx, gate, None, Some(stops), |name| {
            Ok(if name == "cello-g-whole.fsmatpk" { measured_pack(1.0) } else { pack(name, 7800.0, 200_000.0, None) })
        }).err().unwrap().to_string();
        assert!(error.contains("a stop profile changes the twisting span"), "{error}");
    });
    // The record only exists on measured sections, names a table row, and takes exactly one id once.
    for record in ["torsion-measured no-such-row", "torsion-measured", &format!("torsion-measured {D3} extra"),
        &format!("torsion-measured {D3}\ntorsion-measured {D3}")] {
        assert!(input::Spec::parse(selected(&twisted(), record).as_bytes()).is_err(), "accepted {record:?}");
    }
    let homogeneous = include_str!("../../pluck.material-strings").replace("cavity neck", &format!("torsion-measured {D3}\ncavity neck"));
    let error = input::Spec::parse(homogeneous.as_bytes()).err().unwrap().to_string();
    assert!(error.contains("only a section measured string takes it"), "{error}");
}
