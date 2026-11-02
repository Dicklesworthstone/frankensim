//! Shared gas propagation observes the original physical plate trajectory.
use super::*;
use fs_couple::air_path::{GasAbsorptionModel, ObserverFilter};
use fs_couple::render::plate::{CompactPlateVoice, PlateVoiceConfig};
use fs_couple::thin_plate::CompactBody;
use fs_material::iso9613::TerrestrialAir;
use fs_scenario::RadiatingPlate;

fn voice(gas: fs_material::gas::GasState) -> CompactPlateVoice {
    let body = CompactBody::from_radiator(RadiatingPlate {
        area_m2: 0.01,
        mass_kg: 0.25,
        frequency_hz: 2000.0,
        damping_ratio: 0.0,
    })
    .unwrap();
    CompactPlateVoice::from_radiators(
        vec![body],
        0.5,
        PlateVoiceConfig {
            sample_rate_hz: 48_000,
            max_modes: 1,
            nyquist_guard_fraction: 0.9,
            gas,
            listener_m: 30.0,
            maximum_abs_force_n: 1.0,
            maximum_abs_pressure_pa: 10.0,
        },
    )
    .unwrap()
}

#[test]
fn causal_atmosphere_matches_analytical_retarded_acceleration_without_changing_motion() {
    let (gas, air) = TerrestrialAir::moist_air(293.15, 101325.0, 0.5).unwrap();
    let count = 12_000;
    for law in [
        GasAbsorptionModel::StokesKirchhoff,
        GasAbsorptionModel::Iso9613Air {
            relative_humidity: 0.5,
            air,
        },
    ] {
        let filter = ObserverFilter::for_law(law, &gas, None, 30.0, 48_000).unwrap();
        let receipt = *filter.receipt();
        assert!(receipt.receiver_delay_s > 0.0);
        assert!(
            (receipt.receiver_delay_s + receipt.low_frequency_delay_s - 30.0 / gas.sound_speed)
                .abs()
                < 1e-14
        );
        let expected: Vec<_> = (0..count)
            .map(|i| {
                // A resting undamped body under held F has acceleration
                // F/m*cos(omega*t). Evaluate it at EACH actual retarded time,
                // then independently convolve the original atmosphere impulse.
                filter
                    .taps()
                    .iter()
                    .take(i + 1)
                    .enumerate()
                    .map(|(j, tap)| {
                        let time = (i + 1 - j) as f64 / 48_000.0 - receipt.receiver_delay_s;
                        if time < 0.0 {
                            0.0
                        } else {
                            tap * gas.density * 0.01 / (core::f64::consts::TAU * 30.0)
                                * (0.5 / 0.25)
                                * (core::f64::consts::TAU * 2000.0 * time).cos()
                        }
                    })
                    .sum::<f64>()
            })
            .collect();
        let peak = expected.iter().fold(0.0_f64, |m, p| m.max(p.abs()));
        assert!(peak > 1e-6);
        let mut plain = voice(gas);
        let mut unfiltered = vec![0.0; count];
        plain.step_block(&mut unfiltered).unwrap();
        let mut reference: Option<Vec<f64>> = None;
        for block in [1, 37, 571] {
            let mut filtered = voice(gas).with_observer_law(law, None).unwrap();
            assert_eq!(filtered.observer_receipt(), Some(&receipt));
            let mut actual = vec![0.0; count];
            for chunk in actual.chunks_mut(block) {
                filtered.step_block(chunk).unwrap();
            }
            for (a, e) in actual.iter().zip(&expected) {
                assert!((a - e).abs() < 2e-9 * peak + 1e-13, "{a} versus {e}");
            }
            assert_ne!(actual, unfiltered);
            assert_eq!(
                filtered.bodies()[0].volume_velocity().to_bits(),
                plain.bodies()[0].volume_velocity().to_bits()
            );
            let mut a = filtered.bodies()[0].clone();
            let mut b = plain.bodies()[0].clone();
            assert_eq!(
                a.drive(0.0, 1.0 / 48_000.0).unwrap().to_bits(),
                b.drive(0.0, 1.0 / 48_000.0).unwrap().to_bits()
            );
            if let Some(previous) = &reference {
                assert_eq!(&actual, previous);
            } else {
                reference = Some(actual);
            }
        }
    }
}

#[test]
fn lossless_selection_preserves_samples_and_cannot_reset_a_warm_filter() {
    let gas = fs_material::gas::GasState::try_new_moist_air(293.15, 101325.0, 0.5).unwrap();
    let mut plain = voice(gas);
    let mut filtered = voice(gas)
        .with_observer_law(GasAbsorptionModel::Lossless, None)
        .unwrap();
    let mut a = vec![0.0; 5000];
    let mut b = vec![0.0; 5000];
    plain.step_block(&mut a).unwrap();
    filtered.step_block(&mut b).unwrap();
    assert_eq!(a, b);
    assert!(
        filtered
            .with_observer_law(GasAbsorptionModel::Lossless, None)
            .is_err()
    );
    assert!(
        plain
            .with_observer_law(GasAbsorptionModel::Lossless, None)
            .is_err()
    );
    assert!(
        voice(gas)
            .with_observer_law(GasAbsorptionModel::Lossless, None)
            .unwrap()
            .with_observer_law(GasAbsorptionModel::Lossless, None)
            .is_err()
    );
}

#[test]
fn native_plate_file_carries_shared_atmosphere_and_receiver_delay_to_pcm() {
    let source = std::str::from_utf8(SOURCE)
        .unwrap()
        .replace("ambient 293.15 101325 0", "ambient 293.15 101325 0.5")
        .replace(
            "observer 1.0\n",
            "observer 10.0\nradiation-load baffled-piston\nobserver_law iso9613\n",
        );
    let performance = PlatePerformance::from_bytes(source.as_bytes(), 4801).unwrap();
    let info = performance.info();
    let filter = info.observer_filter.unwrap();
    assert!(info.baffled_piston_loading);
    assert!(
        (filter.receiver_delay_s + filter.low_frequency_delay_s - info.propagation_delay_s).abs()
            < 1e-14
    );
    let mut renderer = performance.into_renderer();
    let mut pressure = vec![0.0; 4801];
    renderer.block(&mut pressure).unwrap();
    let (expected, _) = encode_pcm16_wav(&pressure, 48_000, info.full_scale_pa).unwrap();
    let root = scratch();
    let input = root.join("atmosphere.performance");
    std::fs::write(&input, &source).unwrap();
    for block in [37, 571] {
        let output = root.join(format!("atmosphere-{block}.wav"));
        assert_eq!(rendered(&input, &output, block), expected);
        let receipt = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        for fragment in [
            "\"law\":\"iso9613\"",
            "\"receiver_delay_s\":",
            "\"low_frequency_delay_s\":",
            "\"molecular_relaxation_omitted\":false",
        ] {
            assert!(receipt.contains(fragment), "{fragment}");
        }
    }
}

#[test]
fn source_law_admission_refuses_incompatible_gas_and_ignored_fields() {
    let text = std::str::from_utf8(SOURCE).unwrap();
    for row in [
        "observer_law",
        "observer_law automatic",
        "observer_law iso9613 extra",
        "observer_law none\nobserver_law none",
        "observer_law relaxation invented",
    ] {
        let source = text.replace("observer 1.0\n", &format!("observer 1.0\n{row}\n"));
        assert!(
            PlatePerformance::from_bytes(source.as_bytes(), 64).is_err(),
            "{row}"
        );
    }
    let co2 = text.replace(
        "ambient 293.15 101325 0",
        "ambient mixture 293.15 101325 mole 1 CO2 1",
    );
    let declaration = |law: &str| {
        co2.replace(
            "observer 1.0\n",
            &format!("observer 1.0\nobserver_law {law}\n"),
        )
    };
    assert!(PlatePerformance::from_bytes(declaration("iso9613").as_bytes(), 64).is_err());
    assert!(
        PlatePerformance::from_bytes(
            declaration("relaxation iso9613-terrestrial-air").as_bytes(),
            64
        )
        .is_err()
    );
    let classical =
        PlatePerformance::from_bytes(declaration("stokes-kirchhoff").as_bytes(), 64).unwrap();
    let filter = classical.info().observer_filter.unwrap();
    assert_eq!(filter.absorption.law, GasAbsorptionModel::StokesKirchhoff);
    assert!(filter.absorption.relaxation_omitted.is_some());
}
