//! Native molecular/thermal composition with manufactured reservoir data.
use super::*;
use fs_couple::porous_heating::feedback::covariance::{
    PorousCovarianceBin, PorousCrossSpectrumHeating,
};
use fs_couple::porous_heating::feedback::spectrum::{PorousSpectrumHeating, PorousSpectrumOptions};
use fs_duct::porous::PressureCovariance;
use fs_material::relaxation::{RelaxationAbsorption, RelaxingMedium, input::read_relaxation};

fn molecular(active: bool) -> (RelaxationAbsorption, RelaxingMedium) {
    // Temperature-dependent collision time and caloric strength, with explicit
    // source support. These numerical data do not qualify a measured gas.
    let text = if active {
        "declared 1 synthetic 1 vibration 1 900 1 power 101.325 300 1 200 500"
    } else {
        "declared 1 synthetic 0"
    };
    let law = read_relaxation(&mut text.split_ascii_whitespace()).unwrap();
    let medium = RelaxingMedium::new(
        GasState::try_new(&SPEC, 300.0, 101325.0).unwrap(),
        vec![(law.species()[0].species.clone(), 1.0)],
    )
    .unwrap();
    (law, medium)
}

#[test]
fn molecular_temperature_gradient_reaches_original_oblique_fields_and_cell_heat() {
    let materials = materials();
    let (law, reference) = molecular(true);
    let source = IsobaricRelaxingPorousGas::try_new(gas(), &law, &reference).unwrap();
    assert!(std::ptr::eq(source.relaxation_law(), &law));
    assert!(std::ptr::eq(source.reference_medium(), &reference));
    let thermal = thermal(events([0.0, CP * 40.0], [0.0; 2]), [&SENSIBLE; 2]);
    let mut drive = drive();
    drive.reference_impedance_pa_s_m3 /= 0.8;
    let kt = 0.6 * drive.angular_frequency_rad_s / reference.state().sound_speed;
    let heating = TemperaturePorousHeating::try_new_oblique(
        &thermal,
        &cells(&materials),
        source,
        drive,
        kt,
        options(),
        || Ok(()),
    )
    .unwrap();
    let evaluated: Vec<_> = [300.0, 340.0]
        .into_iter()
        .map(|temperature| {
            let medium = RelaxingMedium::new(
                GasState::try_new(&SPEC, temperature, 101325.0).unwrap(),
                reference.components().to_vec(),
            )
            .unwrap();
            law.evaluate(&medium).unwrap()
        })
        .collect();
    let original: Vec<_> = (0..2)
        .map(|i| {
            PorousLayer::try_new_relaxing_oblique(
                &materials[i],
                &evaluated[i],
                drive.angular_frequency_rad_s,
                [0.03, 0.05][i],
                AREA,
                drive.reference_impedance_pa_s_m3,
                kt,
            )
            .unwrap()
        })
        .collect();
    let stack = PorousLayerStack::try_new(&original, options(), || false).unwrap();
    let field = stack
        .prepare_excitation(drive.incident_pressure_pa, || false)
        .unwrap();
    for i in 0..2 {
        let state = heating.fluid_states()[i].source();
        assert_eq!(state, &evaluated[i]);
        assert!(std::ptr::eq(
            state,
            heating.acoustic_layers()[i].relaxation_state().unwrap()
        ));
        assert_eq!(
            heating.acoustic_layers()[i].two_port().scattering,
            original[i].two_port().scattering
        );
        assert_eq!(
            heating.cells()[i].dissipation,
            field.layer_excitations()[i]
                .dissipation_between(0.0, original[i].thickness_m())
                .unwrap()
        );
    }
    assert_eq!(heating.dissipation(), field.dissipation());
    assert!(evaluated[1].gas().density < evaluated[0].gas().density);
    assert!(evaluated[1].gas().sound_speed > evaluated[0].gas().sound_speed);
    assert_ne!(evaluated[0].processes(), evaluated[1].processes());
    let plain = TemperaturePorousHeating::try_new_oblique(
        &thermal,
        &cells(&materials),
        gas(),
        drive,
        kt,
        options(),
        || Ok(()),
    )
    .unwrap();
    assert!(
        (heating.dissipation().absorbed_power_w - plain.dissipation().absorbed_power_w).abs()
            > 1e-5 * plain.dissipation().absorbed_power_w
    );
    assert_eq!(
        heating.boundary_power().incident_power_w,
        plain.boundary_power().incident_power_w
    );
}

#[test]
fn molecular_feedback_credits_start_heat_then_refreshes_only_changed_native_states() {
    let materials = materials();
    let (law, reference) = molecular(true);
    let source = IsobaricRelaxingPorousGas::try_new(gas(), &law, &reference).unwrap();
    let initial = thermal(events([0.0; 2], [0.0; 2]), [&PLATEAU, &SENSIBLE]);
    let heating = TemperaturePorousHeating::try_new(
        &initial,
        &cells(&materials),
        source,
        drive(),
        options(),
        || Ok(()),
    )
    .unwrap();
    let next = heating
        .advance_to(&initial, 0.125, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(
        next.work,
        PorousFeedbackWork {
            fluid_evaluations: 1,
            gas_evaluations: 1,
            stack_solves: 1
        }
    );
    close(
        next.step.acoustic.deposited_j,
        0.125 * heating.dissipation().absorbed_power_w,
        3e-15,
    );
    assert!(std::ptr::eq(
        heating.fluid_states()[0].source(),
        next.heating.fluid_states()[0].source()
    ));
    assert!(!std::ptr::eq(
        heating.fluid_states()[1].source(),
        next.heating.fluid_states()[1].source()
    ));
    assert_eq!(
        next.heating.drive().reference_impedance_pa_s_m3,
        heating.drive().reference_impedance_pa_s_m3
    );
    assert_eq!(
        next.heating.drive().incident_pressure_pa,
        heating.drive().incident_pressure_pa
    );
    for i in 0..2 {
        let state = next.heating.fluid_states()[i].source();
        let temperature = next.step.thermal.temperature_k(id(i)).unwrap();
        assert_eq!(state.gas().temperature, temperature);
        let original_medium = RelaxingMedium::new(
            GasState::try_new(&SPEC, temperature, 101325.0).unwrap(),
            reference.components().to_vec(),
        )
        .unwrap();
        assert_eq!(state, &law.evaluate(&original_medium).unwrap());
        assert!(std::ptr::eq(
            state,
            next.heating.acoustic_layers()[i]
                .relaxation_state()
                .unwrap()
        ));
    }
    assert_ne!(
        heating.acoustic_layers()[1].two_port().scattering,
        next.heating.acoustic_layers()[1].two_port().scattering
    );
    assert_eq!(initial.temperature_k(id(1)).unwrap(), 300.0);
    assert_eq!(heating.fluid_states()[1].source().gas().temperature, 300.0);
}

#[test]
fn molecular_spectrum_and_covariance_share_states_and_retain_identical_heat() {
    let materials = materials();
    let (law, reference) = molecular(true);
    let source = IsobaricRelaxingPorousGas::try_new(gas(), &law, &reference).unwrap();
    let initial = thermal(events([0.0; 2], [0.0; 2]), [&PLATEAU, &SENSIBLE]);
    let mut bins = [
        drive(),
        PorousHarmonicDrive {
            angular_frequency_rad_s: core::f64::consts::TAU * 700.0,
            ..drive()
        },
    ];
    for bin in &mut bins {
        bin.reference_impedance_pa_s_m3 /= (1.0_f64 - 0.4 * 0.4).sqrt();
    }
    let kt = bins.map(|bin| 0.4 * bin.angular_frequency_rad_s / reference.state().sound_speed);
    let limits = PorousSpectrumOptions {
        maximum_bins: 2,
        maximum_cell_bins: 4,
        composition: options(),
    };
    let spectrum = PorousSpectrumHeating::try_new_oblique(
        &initial,
        &cells(&materials),
        source,
        &bins,
        &kt,
        limits,
        || Ok(()),
    )
    .unwrap();
    for bin in spectrum.bins() {
        for i in 0..2 {
            assert!(std::ptr::eq(
                bin.fluid_states()[i].source(),
                spectrum.fluid_states()[i].source()
            ));
            assert!(std::ptr::eq(
                bin.fluid_states()[i].source(),
                bin.acoustic_layers()[i].relaxation_state().unwrap()
            ));
        }
    }
    let covariances = bins.map(|bin| PorousCovarianceBin {
        angular_frequency_rad_s: bin.angular_frequency_rad_s,
        area_m2: bin.area_m2,
        reference_impedance_pa_s_m3: bin.reference_impedance_pa_s_m3,
        covariance: PressureCovariance::from_peak_pressure(bin.incident_pressure_pa).unwrap(),
    });
    let covariance = PorousCrossSpectrumHeating::try_new_oblique(
        &initial,
        &cells(&materials),
        source,
        &covariances,
        &kt,
        limits,
        || Ok(()),
    )
    .unwrap();
    close(
        spectrum.power().dissipation.absorbed_power_w,
        covariance.power().dissipation.absorbed_power_w,
        2e-12,
    );
    let a = spectrum
        .advance_to(&initial, 0.25, config(0.25), || Ok(()))
        .unwrap();
    let b = covariance
        .advance_to(&initial, 0.25, config(0.25), || Ok(()))
        .unwrap();
    for work in [a.work, b.work] {
        assert_eq!(
            work,
            PorousFeedbackWork {
                fluid_evaluations: 1,
                gas_evaluations: 1,
                stack_solves: 2
            }
        );
    }
    close(
        a.step.acoustic.deposited_j,
        b.step.acoustic.deposited_j,
        2e-12,
    );
    for i in 0..2 {
        close(
            a.step.thermal.temperature_k(id(i)).unwrap(),
            b.step.thermal.temperature_k(id(i)).unwrap(),
            2e-10,
        );
        let state = a.heating.fluid_states()[i].source();
        for bin in a.heating.bins() {
            assert!(std::ptr::eq(state, bin.fluid_states()[i].source()));
            assert!(std::ptr::eq(
                state,
                bin.acoustic_layers()[i].relaxation_state().unwrap()
            ));
        }
    }
}

#[test]
fn empty_molecular_feedback_is_exact_and_source_refusals_publish_no_thermal_state() {
    let materials = materials();
    let (empty_law, reference) = molecular(false);
    let source = IsobaricRelaxingPorousGas::try_new(gas(), &empty_law, &reference).unwrap();
    let initial = thermal(events([0.0, CP * 20.0], [0.0; 2]), [&SENSIBLE; 2]);
    let empty = TemperaturePorousHeating::try_new(
        &initial,
        &cells(&materials),
        source,
        drive(),
        options(),
        || Ok(()),
    )
    .unwrap();
    let plain = TemperaturePorousHeating::try_new(
        &initial,
        &cells(&materials),
        gas(),
        drive(),
        options(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(empty.dissipation(), plain.dissipation());
    let a = empty
        .advance_to(&initial, 0.125, config(0.125), || Ok(()))
        .unwrap();
    let b = plain
        .advance_to(&initial, 0.125, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(a.step.acoustic, b.step.acoustic);
    assert_eq!(a.step.thermal.cells(), b.step.thermal.cells());
    let (law, reference) = molecular(true);
    let wrong = IsobaricPorousGas::try_new(&SPEC, 2.0 * 101325.0, [280.0, 450.0]).unwrap();
    assert!(IsobaricRelaxingPorousGas::try_new(wrong, &law, &reference).is_err());
    let too_wide = IsobaricPorousGas::try_new(&SPEC, 101325.0, [280.0, 550.0]).unwrap();
    assert!(IsobaricRelaxingPorousGas::try_new(too_wide, &law, &reference).is_err());
    let bounded = IsobaricPorousGas::try_new(&SPEC, 101325.0, [280.0, 300.000001]).unwrap();
    let source = IsobaricRelaxingPorousGas::try_new(bounded, &law, &reference).unwrap();
    let cold = thermal(events([0.0; 2], [0.0; 2]), [&SENSIBLE; 2]);
    let heating = TemperaturePorousHeating::try_new(
        &cold,
        &cells(&materials),
        source,
        drive(),
        options(),
        || Ok(()),
    )
    .unwrap();
    assert!(
        heating
            .advance_to(&cold, 0.125, config(0.125), || Ok(()))
            .is_err()
    );
    assert!(
        heating
            .advance_to(&cold, 0.125, config(0.125), || Err(
                DepositionError::Cancelled
            ))
            .is_err()
    );
    assert_eq!(cold.temperature_k(id(0)).unwrap(), 300.0);
    assert_eq!(heating.fluid_states()[0].source().gas().temperature, 300.0);
}
