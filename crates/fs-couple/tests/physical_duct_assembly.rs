//! The ordinary assembly observes the same physical bore and material media.
use fs_couple::acoustic_realize::{
    BoreRealization, RealizeOptions, duct_observer::BaffledDuctObserver,
    realize_assembly_with_options,
};
use fs_couple::air_path::{GasAbsorptionModel, ObserverFilter};
use fs_couple::driving_point::regional::RegionalDuct;
use fs_couple::pcm_wav::baffled::CircularOutletReceiver;
use fs_couple::render::ReedBoreVoice;
use fs_couple::thin_plate::PlateBank;
use fs_duct::{Duct, Segment, Termination};
use fs_material::gas::{GasSpec, GasState};
use fs_material::gas_species::GasSpecies;
use fs_scenario::{
    AcousticAssembly, AmbientGas, BeatingReed, CylinderSegment, Listener, ToneHole,
    ViscothermalDuct, VolumeVelocityPulse, WaveguideEnd,
};

fn gas(temperature: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, 100_000.0).unwrap()
}
fn input() -> AcousticAssembly {
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
            segments: [0.2, 0.3]
                .map(|length_m| CylinderSegment {
                    radius_m: 0.0015,
                    outlet_radius_m: 0.0015,
                    length_m,
                })
                .to_vec(),
            tone_holes: vec![],
            termination: WaveguideEnd::FlangedOpen,
            wall: None,
        }),
        reed: Some(BeatingReed {
            rest_opening_m: 0.0004,
            width_m: 0.013,
            closing_pressure_pa: 6000.0,
            blowing_pressure_pa: 2800.0,
            attack_s: 0.002,
            mass_kg: 0.0,
            stiffness_n_m: 0.0,
            damping_ratio: 0.35,
        }),
        listener: Listener { distance_m: 1.0 },
        sample_rate_hz: 8_000,
        duration_s: 0.064,
    }
}
fn binding(gases: Option<&[GasState]>, exterior: Option<GasState>) -> BaffledDuctObserver<'_> {
    BaffledDuctObserver {
        segment_gases: gases,
        exterior_gas: exterior,
        receiver: CircularOutletReceiver {
            position_m: [0.0, 0.0, 0.5],
            radial_rings: 4,
            angular_points: 16,
            maximum_frequency_hz: 500.0,
        },
    }
}
fn options(binding: BaffledDuctObserver<'_>) -> RealizeOptions<'_> {
    RealizeOptions {
        baffled_duct: Some(binding),
        observer_absorption: Some(GasAbsorptionModel::Lossless),
        ..RealizeOptions::default()
    }
}

/// G1/G3: compare the public assembly against the existing pressure/flow and
/// spatial receiver owners, then independently check the earliest causal arrival.
#[test]
fn physical_assembly_matches_original_owners_and_retains_propagation_delay() {
    let mut input = input();
    // Explicit material and receiver inputs must not consult these unused defaults.
    input.ambient.temperature_k = -1.0;
    input.listener.distance_m = f64::NAN;
    let gases = [gas(293.15), gas(310.15)];
    let binding = binding(Some(&gases), Some(gas(283.15)));
    let actual = realize_assembly_with_options(&input, &options(binding)).unwrap();
    let info = actual.duct_observer.as_ref().unwrap();
    assert_eq!(info.segment_gases, gases);
    assert_eq!(actual.gas, binding.exterior_gas.unwrap());
    let physical = Duct {
        segments: [0.2, 0.3]
            .map(|length| Segment::Cylinder {
                radius: 0.0015,
                length,
            })
            .to_vec(),
    };
    let regional = RegionalDuct::new(&physical, &gases)
        .unwrap()
        .with_exterior_gas(actual.gas)
        .unwrap();
    let mut voice = ReedBoreVoice::new_with_regional(
        &regional,
        input.reed.unwrap(),
        Termination::FlangedOpen,
        PlateBank::default(),
        0.5,
        input.sample_rate_hz,
        actual.pressure_pa.len(),
        None,
    )
    .unwrap();
    let mut microphone = regional
        .baffled_outlet(
            input.sample_rate_hz,
            voice.line_realization().n_fft,
            binding.receiver,
            None,
        )
        .unwrap();
    microphone.step(voice.outgoing_pressure_pa()).unwrap();
    for &pressure in &actual.pressure_pa {
        voice.step_block(&mut [0.0]).unwrap();
        assert_eq!(
            pressure,
            microphone.step(voice.outgoing_pressure_pa()).unwrap()
        );
    }
    let first_possible = info.axial_delay_samples + info.exterior_delay_samples.0 - 1;
    assert!(
        actual.pressure_pa[..first_possible]
            .iter()
            .all(|p| *p == 0.0)
    );
    assert!(actual.pressure_pa.iter().all(|p| p.is_finite()));
    assert!(actual.pressure_pa.iter().any(|p| *p != 0.0));
}

/// G3/G5: exterior material changes cannot replace the authored interior or
/// its travel time. Interior temperature/composition reaches actual propagation.
#[test]
fn physical_assembly_keeps_interior_exterior_and_observer_law_independent() {
    let input = input();
    let gases = [gas(293.15); 2];
    let homogeneous =
        realize_assembly_with_options(&input, &options(binding(Some(&gases), None))).unwrap();
    let same =
        realize_assembly_with_options(&input, &options(binding(Some(&gases), Some(gases[1]))))
            .unwrap();
    assert_eq!(homogeneous, same);
    let exterior = gas(323.15);
    let outside =
        realize_assembly_with_options(&input, &options(binding(Some(&gases), Some(exterior))))
            .unwrap();
    assert_ne!(homogeneous.pressure_pa, outside.pressure_pa);
    assert_eq!(outside.duct_observer.as_ref().unwrap().segment_gases, gases);
    assert_eq!(
        homogeneous
            .duct_observer
            .as_ref()
            .unwrap()
            .axial_delay_samples,
        outside.duct_observer.as_ref().unwrap().axial_delay_samples
    );
    let helium = GasState::try_new(&GasSpecies::Helium.spec(), 293.15, 100_000.0).unwrap();
    let other_inside = [gases[0], helium];
    let inside = realize_assembly_with_options(
        &input,
        &options(binding(Some(&other_inside), Some(exterior))),
    )
    .unwrap();
    assert_eq!(inside.gas, outside.gas);
    assert_ne!(inside.pressure_pa, outside.pressure_pa);
    assert!(
        inside.duct_observer.as_ref().unwrap().axial_delay_samples
            < outside.duct_observer.as_ref().unwrap().axial_delay_samples
    );
    let mut selected = options(binding(Some(&other_inside), Some(exterior)));
    selected.observer_absorption = Some(GasAbsorptionModel::StokesKirchhoff);
    let absorbed = realize_assembly_with_options(&input, &selected).unwrap();
    let original = ObserverFilter::stokes_kirchhoff(&exterior, 0.5, input.sample_rate_hz).unwrap();
    assert_eq!(
        absorbed.duct_observer.as_ref().unwrap().filter,
        *original.receipt()
    );
    let receipt = absorbed.duct_observer.as_ref().unwrap().filter;
    assert!(
        (receipt.receiver_delay_s + receipt.low_frequency_delay_s - 0.5 / exterior.sound_speed)
            .abs()
            < 1e-15
    );
    assert_ne!(absorbed.pressure_pa, inside.pressure_pa);
}

/// G0: unsupported participants, incoherent snapshots and conflicting source
/// assignments refuse rather than becoming a partial or differently loaded render.
#[test]
fn physical_assembly_refuses_incomplete_or_inconsistent_physical_bindings() {
    let input = input();
    let gases = [gas(293.15); 2];
    let selected = options(binding(Some(&gases), None));
    for altered in [
        RealizeOptions {
            observer_absorption: None,
            ..selected
        },
        RealizeOptions {
            gas: Some(gases[0]),
            ..selected
        },
        RealizeOptions {
            bore: BoreRealization::lumped_ladder(16),
            ..selected
        },
        RealizeOptions {
            observer_absorption: Some(GasAbsorptionModel::MolecularRelaxation),
            ..selected
        },
    ] {
        assert!(realize_assembly_with_options(&input, &altered).is_err());
    }
    let mut pressure_mismatch = gases;
    pressure_mismatch[1] =
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 80_000.0).unwrap();
    let mut stale = gases;
    stale[1].sound_speed *= 1.1;
    for medium in [&gases[..1], &pressure_mismatch, &stale] {
        assert!(
            realize_assembly_with_options(&input, &options(binding(Some(medium), None))).is_err()
        );
    }
    let mut coupled = input.clone();
    coupled.blow = Some(VolumeVelocityPulse {
        peak_m3_s: 1e-5,
        duration_s: 0.01,
    });
    assert!(realize_assembly_with_options(&coupled, &selected).is_err());
    let mut unflanged = input.clone();
    unflanged.duct.as_mut().unwrap().termination = WaveguideEnd::UnflangedOpen;
    assert!(realize_assembly_with_options(&unflanged, &selected).is_err());
    let mut branched = input;
    branched.duct.as_mut().unwrap().tone_holes.push(ToneHole {
        after_segment: 0,
        radius_m: 0.0005,
        chimney_m: 0.001,
        open_fraction: 1.0,
    });
    assert!(realize_assembly_with_options(&branched, &selected).is_err());
}
