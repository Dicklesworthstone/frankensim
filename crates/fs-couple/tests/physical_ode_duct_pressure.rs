//! G1/G3: actual LC outlet flow reaches the existing physical microphone.
use fs_couple::acoustic_realize::{
    BoreRealization, BoreRealizationReceipt, RealizeOptions, duct_observer::BaffledDuctObserver,
    realize_assembly_with_options,
};
use fs_couple::air_path::GasAbsorptionModel;
use fs_couple::pcm_wav::baffled::CircularOutletReceiver;
use fs_material::gas::{GasSpec, GasState};
use fs_scenario::{
    AcousticAssembly, AmbientGas, BeatingReed, CylinderSegment, Listener, ViscothermalDuct,
    WaveguideEnd,
};

const RATE: u32 = 16_000;

fn gas(pressure_pa: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure_pa).unwrap()
}

fn input(mass_kg: f64) -> AcousticAssembly {
    AcousticAssembly {
        ambient: AmbientGas::sea_level(),
        string: None,
        pluck: None,
        bow: None,
        blow: None,
        soundboard: None,
        body_modes: vec![],
        plate: None,
        cavity: None,
        obstacles: vec![],
        contact_texture: None,
        duct: Some(ViscothermalDuct {
            segments: vec![CylinderSegment {
                radius_m: 0.004,
                outlet_radius_m: 0.004,
                length_m: 0.3,
            }],
            tone_holes: vec![],
            termination: WaveguideEnd::FlangedOpen,
            wall: None,
        }),
        reed: Some(BeatingReed {
            rest_opening_m: 0.0004,
            width_m: 0.01,
            closing_pressure_pa: 10_000.0,
            blowing_pressure_pa: 100.0,
            attack_s: 0.001,
            mass_kg,
            stiffness_n_m: 100.0,
            damping_ratio: 0.1,
        }),
        listener: Listener {
            distance_m: f64::NAN,
        },
        sample_rate_hz: RATE,
        duration_s: 128.0 / f64::from(RATE),
    }
}

fn receiver(gas: GasState, path_samples: usize) -> CircularOutletReceiver {
    // With one equal-area annulus every ray has this exact travel length,
    // even though the original finite-aperture receiver remains in use.
    let path = path_samples as f64 * gas.sound_speed / f64::from(RATE);
    CircularOutletReceiver {
        position_m: [0.0, 0.0, (path * path - 0.004_f64.powi(2) / 2.0).sqrt()],
        radial_rings: 1,
        angular_points: 8,
        maximum_frequency_hz: 500.0,
    }
}

fn options(gas: GasState, path_samples: usize) -> RealizeOptions<'static> {
    RealizeOptions {
        // Deliberately coarse numerical fixture; its original passive-pitch
        // admission still runs. This test is not a musical accuracy claim.
        bore: BoreRealization::lumped_ladder_with_tolerance(4, 500.0),
        gas: Some(gas),
        observer_absorption: Some(GasAbsorptionModel::Lossless),
        baffled_duct: Some(BaffledDuctObserver {
            segment_gases: None,
            exterior_gas: None,
            material_wall: None,
            receiver: receiver(gas, path_samples),
        }),
        ..RealizeOptions::default()
    }
}

#[test]
fn massive_and_quasistatic_reeds_emit_retarded_outlet_pressure_with_inverse_distance_amplitude() {
    for mass in [0.0, 1e-4] {
        let input = input(mass);
        let gas = gas(101_325.0);
        let near = realize_assembly_with_options(&input, &options(gas, 16)).unwrap();
        let far = realize_assembly_with_options(&input, &options(gas, 32)).unwrap();
        assert!(matches!(
            near.bore_receipt,
            Some(BoreRealizationReceipt::LumpedLadder { cells: 4, .. })
        ));
        let binding = near.duct_observer.as_ref().unwrap();
        assert_eq!(binding.axial_delay_samples, 0);
        assert!(near.pressure_pa[..15].iter().all(|p| *p == 0.0));
        let amplitude = near
            .pressure_pa
            .iter()
            .map(|p| p.abs())
            .fold(0.0_f64, f64::max);
        assert!(amplitude > 0.0 && amplitude.is_finite());
        for i in 0..near.pressure_pa.len() - 16 {
            assert!(
                (far.pressure_pa[i + 16] - 0.5 * near.pressure_pa[i]).abs() < 1e-10 * amplitude
            );
        }
        // The explicitly selected receiver never contains instantaneous local
        // mouth pressure, even though that diagnostic has the same Pa units.
        let mut diagnostic = options(gas, 16);
        diagnostic.baffled_duct = None;
        let mut local = input.clone();
        local.listener.distance_m = 1.0;
        let source = realize_assembly_with_options(&local, &diagnostic).unwrap();
        assert!(source.pressure_pa[..15].iter().any(|p| *p != 0.0));
    }
}

#[test]
fn pressure_changes_reach_the_same_source_and_microphone_without_pressure_normalization() {
    let input = input(0.0);
    let low = realize_assembly_with_options(&input, &options(gas(70_000.0), 16)).unwrap();
    let high = realize_assembly_with_options(&input, &options(gas(140_000.0), 16)).unwrap();
    assert_eq!(low.gas.pressure, 70_000.0);
    assert_eq!(high.gas.pressure, 140_000.0);
    let difference: f64 = low
        .pressure_pa
        .iter()
        .zip(&high.pressure_pa)
        .map(|(a, b)| (a - b).abs())
        .sum();
    assert!(difference > 1e-8);
}

#[test]
fn unsupported_boundary_and_region_physics_are_not_silently_omitted() {
    let gas = gas(101_325.0);
    let mut input = input(0.0);
    let mut options = options(gas, 16);
    options.baffled_duct.as_mut().unwrap().exterior_gas =
        Some(GasState::try_new(&GasSpec::dry_air_ussa1976(), 310.0, gas.pressure).unwrap());
    assert!(realize_assembly_with_options(&input, &options).is_err());
    options.baffled_duct.as_mut().unwrap().exterior_gas = None;
    input.duct.as_mut().unwrap().termination = WaveguideEnd::Closed;
    assert!(realize_assembly_with_options(&input, &options).is_err());
    input.duct.as_mut().unwrap().termination = WaveguideEnd::FlangedOpen;
    input.duct.as_mut().unwrap().segments[0].outlet_radius_m = 0.006;
    // A lone widening cone uses the legacy spherical psi image, whose loaded
    // boundary omits a physical Euler near-field flow. Do not turn its load
    // current into a physical microphone; leave ordinary cone rendering intact.
    assert!(realize_assembly_with_options(&input, &options).is_err());
    options.baffled_duct = None;
    input.listener.distance_m = 1.0;
    let legacy = realize_assembly_with_options(&input, &options).unwrap();
    assert!(legacy.pressure_pa.iter().all(|p| p.is_finite()));
    assert!(legacy.pressure_pa.iter().any(|p| *p != 0.0));
}
