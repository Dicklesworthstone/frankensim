use super::*;
use fs_couple::porous_heating::feedback::covariance::angular::{
    DiffuseAngularNode, DiffuseHemisphere, PorousAngularMode,
};

fn midpoint(n: usize) -> Vec<DiffuseAngularNode> {
    (0..n)
        .map(|i| DiffuseAngularNode {
            cosine: (i as f64 + 0.5) / n as f64,
            weight_mu: 1.0 / n as f64,
        })
        .collect()
}
fn angular_mode(mut bin: PorousCovarianceBin, cosine: f64) -> PorousAngularMode {
    let exterior = GasState::try_new(&SPEC, 300.0, 101325.0).unwrap();
    bin.reference_impedance_pa_s_m3 /= cosine;
    PorousAngularMode {
        tangential_wavenumber_rad_m: bin.angular_frequency_rad_s / exterior.sound_speed
            * ((1.0 - cosine) * (1.0 + cosine)).sqrt(),
        bin,
    }
}

#[test]
fn g0_independent_mode_keys_preserve_normal_parity_and_original_frequency_refusal() {
    let materials = materials();
    let initial = thermal(events(), [&SENSIBLE; 2]);
    let routes = routes(&materials);
    let normal_modes = bins().map(|bin| angular_mode(bin, 1.0));
    let plain = PorousCrossSpectrumHeating::try_new(
        &initial,
        &routes,
        gas(),
        &bins(),
        options(),
        || Ok(()),
    )
    .unwrap();
    let normal = PorousCrossSpectrumHeating::try_new_independent_modes(
        &initial,
        &routes,
        gas(),
        &normal_modes,
        options(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(plain.power(), normal.power());
    assert_eq!(plain.cells(), normal.cells());
    let plain = plain
        .advance_to(&initial, 0.25, config(0.125), || Ok(()))
        .unwrap();
    let normal = normal
        .advance_to(&initial, 0.25, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(plain.step.thermal.cells(), normal.step.thermal.cells());
    assert_eq!(plain.step.acoustic, normal.step.acoustic);
    assert_eq!(plain.bin_energy, normal.bin_energy);
    assert_eq!(plain.work, normal.work);

    let modes = [
        angular_mode(bins()[0], 0.3),
        angular_mode(bins()[1], 0.6),
        angular_mode(bins()[0], 0.9),
    ];
    let combined = PorousCrossSpectrumHeating::try_new_independent_modes(
        &initial,
        &routes,
        gas(),
        &modes,
        options(),
        || Ok(()),
    )
    .unwrap();
    for (canonical, original) in [1, 2, 0].into_iter().enumerate() {
        let mode = modes[original];
        assert_eq!(
            combined.bins()[canonical].angular_frequency_rad_s,
            mode.bin.angular_frequency_rad_s
        );
        assert_eq!(
            combined.bins()[canonical].reference_impedance_pa_s_m3,
            mode.bin.reference_impedance_pa_s_m3
        );
        assert_eq!(
            combined.tangential_wavenumber_rad_m(canonical),
            Some(mode.tangential_wavenumber_rad_m)
        );
        let single = PorousCrossSpectrumHeating::try_new_oblique(
            &initial,
            &routes,
            gas(),
            &[mode.bin],
            &[mode.tangential_wavenumber_rad_m],
            options(),
            || Ok(()),
        )
        .unwrap();
        assert_eq!(combined.bin_power(canonical).unwrap(), single.power());
        assert_eq!(combined.bin_cells(canonical).unwrap(), single.cells());
    }
    let reversed = [modes[2], modes[1], modes[0]];
    let reversed = PorousCrossSpectrumHeating::try_new_independent_modes(
        &initial,
        &routes,
        gas(),
        &reversed,
        options(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(combined.power(), reversed.power());
    assert_eq!(combined.cells(), reversed.cells());
    // The old entry point never turns repeated physical frequency into an
    // independence declaration, even when the two tangential modes differ.
    assert!(
        PorousCrossSpectrumHeating::try_new_oblique(
            &initial,
            &routes,
            gas(),
            &[modes[0].bin, modes[2].bin],
            &[
                modes[0].tangential_wavenumber_rad_m,
                modes[2].tangential_wavenumber_rad_m
            ],
            options(),
            || Ok(())
        )
        .is_err()
    );
    assert!(
        PorousCrossSpectrumHeating::try_new_independent_modes(
            &initial,
            &routes,
            gas(),
            &[modes[0]; 2],
            options(),
            || Ok(())
        )
        .is_err()
    );
    let mut duplicate = [normal_modes[0]; 2];
    duplicate[1].tangential_wavenumber_rad_m = -0.0;
    assert!(
        PorousCrossSpectrumHeating::try_new_independent_modes(
            &initial,
            &routes,
            gas(),
            &duplicate,
            options(),
            || Ok(())
        )
        .is_err()
    );
    for kt in [-1.0, f64::NAN, f64::INFINITY] {
        let invalid = PorousAngularMode {
            bin: PorousCovarianceBin {
                covariance: covariance(0.0, 0.0, C64::ZERO),
                ..modes[0].bin
            },
            tangential_wavenumber_rad_m: kt,
        };
        assert!(
            PorousCrossSpectrumHeating::try_new_independent_modes(
                &initial,
                &routes,
                gas(),
                &[invalid],
                options(),
                || Ok(())
            )
            .is_err()
        );
    }
}

#[test]
fn g1_diffuse_quadrature_preserves_normal_flux_and_finite_grazing_pressure() {
    let exterior = GasState::try_new(&SPEC, 300.0, 101325.0).unwrap();
    let irradiance = [12.0, 4.0];
    let diffuse = DiffuseHemisphere::try_new(&exterior, 6000.0, AREA, irradiance).unwrap();
    assert_eq!(diffuse.irradiance_w_m2(), irradiance);
    assert_eq!(diffuse.exterior().sound_speed_m_s(), exterior.sound_speed);
    let modes = diffuse.modes(&midpoint(8), 8, || Ok(())).unwrap();
    let mut power = [0.0; 2];
    let mut variance = [0.0; 2];
    for mode in modes {
        let covariance = mode.bin.covariance.matrix_pa2();
        assert_eq!(covariance[1], C64::ZERO);
        assert_eq!(covariance[2], C64::ZERO);
        for side in 0..2 {
            power[side] += covariance[side * 3].re / mode.bin.reference_impedance_pa_s_m3;
            variance[side] += covariance[side * 3].re;
        }
    }
    for side in 0..2 {
        close(power[side], AREA * irradiance[side], 2e-15);
        close(
            variance[side],
            2.0 * exterior.characteristic_impedance * irradiance[side],
            2e-15,
        );
    }
    let epsilon = 1e-8;
    let nodes = [epsilon, 1.0 - epsilon].map(|cosine| DiffuseAngularNode {
        cosine,
        weight_mu: 0.5,
    });
    let grazing = diffuse.modes(&nodes, 2, || Ok(())).unwrap();
    assert_eq!(
        grazing[0].bin.covariance.matrix_pa2(),
        grazing[1].bin.covariance.matrix_pa2()
    );
    close(
        grazing[0].bin.covariance.matrix_pa2()[0].re / grazing[0].bin.reference_impedance_pa_s_m3,
        AREA * irradiance[0] * epsilon,
        2e-15,
    );
    assert!(grazing[0].bin.reference_impedance_pa_s_m3.is_finite());
    let materials = materials();
    let initial = thermal(events(), [&SENSIBLE; 2]);
    let native = PorousCrossSpectrumHeating::try_new_independent_modes(
        &initial,
        &routes(&materials),
        gas(),
        &grazing,
        options(),
        || Ok(()),
    )
    .unwrap();
    close(
        native.power().incident_power_w,
        AREA * irradiance.iter().sum::<f64>(),
        2e-15,
    );
    assert!(native.power().dissipation.absorbed_power_w.is_finite());

    let dark = DiffuseHemisphere::try_new(&exterior, 6000.0, AREA, [0.0; 2]).unwrap();
    for cosine in [0.0, -1.0, 1.1, f64::NAN] {
        assert!(
            dark.modes(
                &[DiffuseAngularNode {
                    cosine,
                    weight_mu: 1.0
                }],
                1,
                || Ok(())
            )
            .is_err()
        );
    }
    for weight_mu in [0.0, -1.0, f64::NAN, f64::INFINITY, 0.5] {
        assert!(
            dark.modes(
                &[DiffuseAngularNode {
                    cosine: 0.5,
                    weight_mu
                }],
                1,
                || Ok(())
            )
            .is_err()
        );
    }
    assert!(
        dark.modes(
            &[DiffuseAngularNode {
                cosine: 0.9,
                weight_mu: 1.0
            }],
            1,
            || Ok(())
        )
        .is_err()
    );
    assert!(diffuse.modes(&midpoint(8), 7, || Ok(())).is_err());
    assert!(diffuse.modes(&midpoint(1), 4097, || Ok(())).is_err());
    assert!(
        diffuse
            .modes(&midpoint(1), 1, || Err(DepositionError::Cancelled))
            .is_err()
    );
    let mut stale = exterior;
    stale.sound_speed *= 2.0;
    assert!(DiffuseHemisphere::try_new(&stale, 6000.0, AREA, irradiance).is_err());
    // A tiny positive weight must not silently become a dark component even
    // when an alternate product ordering could represent the final variance.
    let large_z = fs_material::fluid_acoustics::FluidAcousticState::try_new(1e200, 1.0).unwrap();
    let tiny = f64::from_bits(1);
    assert!((1e200 * 1e-100) * (2.0 * tiny) > 0.0);
    let tiny_node = [
        DiffuseAngularNode {
            cosine: 0.3,
            weight_mu: tiny,
        },
        DiffuseAngularNode {
            cosine: 0.5,
            weight_mu: 1.0,
        },
    ];
    assert!(
        DiffuseHemisphere::try_new(&large_z, 6000.0, 1.0, [1e-100, 0.0])
            .unwrap()
            .modes(&tiny_node, 2, || Ok(()))
            .is_err()
    );
}

#[test]
fn g3_same_frequency_modes_share_one_native_trajectory_and_atomic_refresh() {
    let materials = materials();
    let initial = thermal(events(), [&PLATEAU, &SENSIBLE]);
    let routes = routes(&materials);
    let exterior = GasState::try_new(&SPEC, 300.0, 101325.0).unwrap();
    let modes: Vec<_> = [450.0, 1100.0]
        .into_iter()
        .flat_map(|frequency| {
            DiffuseHemisphere::try_new(
                &exterior,
                core::f64::consts::TAU * frequency,
                AREA,
                [20.0, 5.0],
            )
            .unwrap()
            .modes(&midpoint(4), 4, || Ok(()))
            .unwrap()
        })
        .collect();
    let options = PorousSpectrumOptions {
        maximum_bins: 8,
        maximum_cell_bins: 16,
        ..options()
    };
    let source = PorousCrossSpectrumHeating::try_new_independent_modes(
        &initial,
        &routes,
        gas(),
        &modes,
        options,
        || Ok(()),
    )
    .unwrap();
    let initial_power = source.power();
    let first = source
        .advance_to(&initial, 0.25, config(0.125), || Ok(()))
        .unwrap();
    let continued = first
        .heating
        .advance_to(&first.step.thermal, 0.5, config(0.125), || Ok(()))
        .unwrap();
    let whole = source
        .advance_to(&initial, 0.5, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(whole.step.thermal.cells(), continued.step.thermal.cells());
    assert_eq!(whole.heating.power(), continued.heating.power());
    assert_eq!(whole.work.gas_evaluations, 4);
    assert_eq!(whole.work.stack_solves, 32);
    assert_eq!(whole.step.thermal.temperature_k(IDS[0]).unwrap(), 300.0);
    close(whole.step.acoustic.incident_j, 0.5 * AREA * 50.0, 3e-15);
    close(
        whole.bin_energy.iter().map(|e| e.deposited_j).sum(),
        whole.step.acoustic.deposited_j,
        3e-15,
    );
    close(
        whole.step.heat_report.external_heat_j,
        whole.step.acoustic.deposited_j + 0.075,
        3e-13,
    );
    for i in 0..modes.len() {
        assert_eq!(
            source.bins()[i].angular_frequency_rad_s,
            whole.heating.bins()[i].angular_frequency_rad_s
        );
        assert_eq!(
            source.bins()[i].reference_impedance_pa_s_m3,
            whole.heating.bins()[i].reference_impedance_pa_s_m3
        );
        assert_eq!(
            source.tangential_wavenumber_rad_m(i),
            whole.heating.tangential_wavenumber_rad_m(i)
        );
        close(
            first.bin_energy[i].deposited_j + continued.bin_energy[i].deposited_j,
            whole.bin_energy[i].deposited_j,
            3e-15,
        );
    }
    assert!(
        source
            .advance_to(&initial, 0.5, config(0.125), || Err(
                DepositionError::Cancelled
            ))
            .is_err()
    );
    assert_eq!(source.power(), initial_power);
    assert_eq!(initial.time_s(), 0.0);
}
