//! G1/G3: the public absorbing-incidence query retains original material laws.
use super::*;
use fs_material::conductor::{ElectricalPropertyBasis, ResolvedConductorState};

#[test]
fn original_debye_and_ohmic_sources_drive_absorbing_incidence_without_changing_old_ports() {
    let outside = material(4.0 * E, M);
    let dielectric = finite_layers::debye(2.0 * E, 3.0 * E, 0.7, 0.002);
    let sigma: f64 = 0.02;
    let electric = [C64::new(0.7, 0.2), C64::new(-0.3, 0.5)];
    for basis in [
        ElectricalPropertyBasis::Conductivity,
        ElectricalPropertyBasis::Resistivity,
    ] {
        let stored = match basis {
            ElectricalPropertyBasis::Conductivity => sigma,
            ElectricalPropertyBasis::Resistivity => sigma.recip(),
        };
        let conductor = ResolvedConductorState::try_from_material_with_basis(
            source(
                &[
                    (EPS_NAME, EPS, 3.0 * E),
                    (MU_NAME, MU, 1.4 * M),
                    (basis.property_name(), basis.dimensions(), stored),
                ],
                ValidityDomain::unconstrained().with("T", 280.0, 320.0),
                QueryPoint::new().with("T", 300.0).unwrap(),
            ),
            basis,
        )
        .unwrap();
        for omega in [1e7, 1e9, 1e11] {
            let x = omega * 1e-9;
            let debye_epsilon = C64::new(
                2.0 * E + 3.0 * E / (1.0 + x * x),
                -3.0 * E * x / (1.0 + x * x),
            );
            let ohmic_epsilon = C64::new(3.0 * E, -sigma / omega);
            for (incident, exit, eps_i, eps_t, mu_i, mu_t) in [
                (
                    DielectricResponse::Debye(&dielectric),
                    DielectricResponse::Lossless(&outside),
                    debye_epsilon,
                    C64::from_re(4.0 * E),
                    M,
                    M,
                ),
                (
                    DielectricResponse::Ohmic(&conductor),
                    DielectricResponse::Debye(&dielectric),
                    ohmic_epsilon,
                    debye_epsilon,
                    1.4 * M,
                    M,
                ),
                (
                    DielectricResponse::Debye(&dielectric),
                    DielectricResponse::Ohmic(&conductor),
                    debye_epsilon,
                    ohmic_epsilon,
                    M,
                    1.4 * M,
                ),
            ] {
                let interface =
                    PlanarDielectricInterface::try_new_with_incident_response(incident, exit)
                        .unwrap();
                assert!(std::ptr::eq(
                    interface.incident_material(),
                    incident.material()
                ));
                assert!(std::ptr::eq(
                    interface.transmitted_material(),
                    exit.material()
                ));
                // The old angle/R/T API still refuses absorbing incident ports.
                assert!(interface.response(omega, 0.8, Polarization::Te).is_err());
                let kx = 0.6 * omega * (E * M).sqrt();
                let result = interface
                    .response_at_transverse_wavenumber(omega, kx, electric)
                    .unwrap();
                complex_near(result.incident_permittivity_f_m, eps_i, 4e-14 * eps_i.abs());
                complex_near(
                    result.transmitted_permittivity_f_m,
                    eps_t,
                    4e-14 * eps_t.abs(),
                );
                let qi = (eps_i.scale(omega * omega * mu_i) - C64::from_re(kx * kx)).sqrt();
                let qt = (eps_t.scale(omega * omega * mu_t) - C64::from_re(kx * kx)).sqrt();
                let yi = [eps_i.scale(omega) / qi, qi.scale(1.0 / (omega * mu_i))];
                let yt = [eps_t.scale(omega) / qt, qt.scale(1.0 / (omega * mu_t))];
                for c in 0..2 {
                    let expected_r = (yi[c] - yt[c]) / (yi[c] + yt[c]);
                    let expected_t = yi[c].scale(2.0) / (yi[c] + yt[c]);
                    complex_near(result.reflected_tangential_e_gain[c], expected_r, 2e-13);
                    complex_near(result.transmitted_tangential_e_gain[c], expected_t, 2e-13);
                    complex_near(
                        result.reflected.electric_v_m[c],
                        expected_r * electric[c],
                        2e-13,
                    );
                }
                let incident_power = result.incident_normal_power_w_m2;
                assert!(incident_power > 0.0);
                assert!(result.reflected_normal_power_w_m2 <= 0.0);
                assert!(result.transmitted_normal_power_w_m2 >= 0.0);
                near(
                    incident_power
                        + result.reflected_normal_power_w_m2
                        + result.interference_normal_power_w_m2,
                    result.transmitted_normal_power_w_m2,
                    3e-13 * incident_power,
                );
            }
        }
    }
}

#[test]
fn public_query_refuses_nonfinite_wave_data_and_retains_zero_illumination() {
    let incident = finite_layers::debye(2.0 * E, 3.0 * E, 1.0, 1.0);
    let outside = material(E, M);
    let interface = PlanarDielectricInterface::try_new_with_incident_response(
        DielectricResponse::Debye(&incident),
        DielectricResponse::Lossless(&outside),
    )
    .unwrap();
    assert!(
        interface
            .response_at_transverse_wavenumber(0.0, 0.0, [C64::ONE; 2])
            .is_err()
    );
    assert!(
        interface
            .response_at_transverse_wavenumber(1e9, f64::NAN, [C64::ONE; 2])
            .is_err()
    );
    assert!(
        interface
            .response_at_transverse_wavenumber(1e9, 0.0, [C64::new(f64::INFINITY, 0.0); 2])
            .is_err()
    );
    let zero = interface
        .response_at_transverse_wavenumber(1e9, 0.0, [C64::ZERO; 2])
        .unwrap();
    for field in [zero.incident, zero.reflected, zero.transmitted] {
        assert_eq!(field.electric_v_m, [C64::ZERO; 3]);
        assert_eq!(field.magnetic_a_m, [C64::ZERO; 3]);
    }
    assert_eq!(zero.incident_normal_power_w_m2, 0.0);
    assert_eq!(zero.reflected_normal_power_w_m2, 0.0);
    assert_eq!(zero.interference_normal_power_w_m2, 0.0);
    assert_eq!(zero.transmitted_normal_power_w_m2, 0.0);
}

#[test]
fn original_absorbing_sources_prepare_coherent_films_and_positive_depth_heating() {
    use fs_em::plane_wave::layers::{
        DielectricLayer, PlanarDielectricLayers,
        absorbing_incidence::PlanarAbsorbingDielectricLayers,
    };
    let incident = finite_layers::debye(2.0 * E, 3.0 * E, 0.7, 0.002);
    let outside = material(4.0 * E, 1.2 * M);
    let conductor = ResolvedConductorState::try_from_material_with_basis(
        source(
            &[
                (EPS_NAME, EPS, 3.0 * E),
                (MU_NAME, MU, 1.4 * M),
                (
                    ElectricalPropertyBasis::Conductivity.property_name(),
                    ElectricalPropertyBasis::Conductivity.dimensions(),
                    0.02,
                ),
            ],
            ValidityDomain::unconstrained(),
            QueryPoint::new(),
        ),
        ElectricalPropertyBasis::Conductivity,
    )
    .unwrap();
    let layers = [
        DielectricLayer {
            response: DielectricResponse::Ohmic(&conductor),
            thickness_m: 0.02,
        },
        DielectricLayer {
            response: DielectricResponse::Debye(&incident),
            thickness_m: 0.08,
        },
    ];
    let electric = [C64::new(0.7, 0.2), C64::new(-0.3, 0.5)];
    for omega in [1e8, 1e9, 1e10] {
        let kx = 0.4 * omega * (E * M).sqrt();
        let source = PlanarAbsorbingDielectricLayers::try_new(
            DielectricResponse::Debye(&incident),
            DielectricResponse::Lossless(&outside),
            omega,
            kx,
            &layers,
        )
        .unwrap();
        assert!(std::ptr::eq(
            source.incident_response().material(),
            incident.material()
        ));
        assert!(std::ptr::eq(source.exit_response().material(), &outside));
        assert!(std::ptr::eq(
            source.layer(0).unwrap().response.material(),
            conductor.material()
        ));
        assert_eq!(source.layer_count(), 2);
        assert_eq!(source.angular_frequency_rad_s(), omega);
        assert_eq!(source.transverse_wavenumber_rad_m(), kx);
        let x = omega * 1e-9;
        complex_near(
            source.incident_permittivity_f_m(),
            C64::new(
                2.0 * E + 3.0 * E / (1.0 + x * x),
                -3.0 * E * x / (1.0 + x * x),
            ),
            1e-24,
        );
        // The new explicit-kx owner does not weaken old normalized-port gates.
        assert!(
            PlanarDielectricLayers::try_new_with_incident_response(
                DielectricResponse::Debye(&incident),
                DielectricResponse::Lossless(&outside),
                omega,
                &layers,
            )
            .is_err()
        );
        let actual = source.illuminate(electric).unwrap();
        assert!(std::ptr::eq(actual.source(), &source));
        let response = actual.response();
        assert!(
            actual
                .absorbed_power_per_layer_w_m2()
                .iter()
                .all(|&p| p > 0.0)
        );
        near(
            actual.absorbed_power_per_layer_w_m2().iter().sum(),
            response.absorbed_power_w_m2,
            1e-16,
        );
        near(
            response.incident_normal_power_w_m2
                + response.reflected_normal_power_w_m2
                + response.interference_normal_power_w_m2,
            response.transmitted_normal_power_w_m2 + response.absorbed_power_w_m2,
            1e-15,
        );
        for (i, layer) in layers.iter().enumerate() {
            let n = 200;
            let dz = layer.thickness_m / f64::from(n);
            let mut integral = 0.0;
            for j in 0..=n {
                let point = actual
                    .layer_field(i, layer.thickness_m * f64::from(j) / f64::from(n))
                    .unwrap();
                let weight = if j == 0 || j == n {
                    1.0
                } else if j % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                integral += weight * point.heating_w_m3;
            }
            near(
                integral * dz / 3.0,
                actual.absorbed_power_per_layer_w_m2()[i],
                2e-10 * response.incident_normal_power_w_m2,
            );
        }
        assert!(actual.layer_field(2, 0.0).is_err());
        assert!(actual.layer_field(0, -0.1).is_err());
    }
}

#[test]
fn original_absorbing_sources_support_phase_controlled_two_sided_film_heating() {
    use fs_em::plane_wave::layers::{
        DielectricLayer, absorbing_incidence::PlanarAbsorbingDielectricLayers,
    };
    let dielectric = finite_layers::debye(2.0 * E, 3.0 * E, 0.7, 0.002);
    let source = DielectricResponse::Debye(&dielectric);
    let layers = [DielectricLayer {
        response: source,
        thickness_m: 0.08,
    }];
    let stack =
        PlanarAbsorbingDielectricLayers::try_new(source, source, 1e9, 0.0, &layers).unwrap();
    let prepared = stack.prepare_two_sided().unwrap();
    assert!(std::ptr::eq(prepared.source(), &stack));
    assert!(std::ptr::eq(
        prepared.source().incident_response().material(),
        dielectric.material()
    ));
    assert!(std::ptr::eq(
        prepared.source().exit_response().material(),
        dielectric.material()
    ));
    for channel in 0..2 {
        let mut amplitude = [C64::ZERO; 2];
        amplitude[channel] = C64::ONE;
        let even = prepared.illuminate(amplitude, amplitude).unwrap();
        let odd = prepared
            .illuminate(amplitude, amplitude.map(|e| -e))
            .unwrap();
        assert!(std::ptr::eq(even.source(), &stack));
        assert!(even.response().absorbed_power_w_m2 > 10.0 * odd.response().absorbed_power_w_m2);
        assert!(odd.response().absorbed_power_w_m2 > 0.0);
        for actual in [&even, &odd] {
            let result = actual.response();
            near(
                result.power.iter().map(|p| p.net_into_stack_w_m2).sum(),
                result.absorbed_power_w_m2,
                2e-15,
            );
            near(
                actual.absorbed_power_per_layer_w_m2()[0],
                result.absorbed_power_w_m2,
                2e-15,
            );
            assert!(actual.layer_field(0, 0.02).unwrap().heating_w_m3 > 0.0);
        }
        let single = prepared.illuminate(amplitude, [C64::ZERO; 2]).unwrap();
        let original = stack.illuminate(amplitude).unwrap();
        near(
            single.response().absorbed_power_w_m2,
            original.response().absorbed_power_w_m2,
            1e-15,
        );
    }
}

#[test]
fn explicit_kx_two_sided_fields_reduce_to_the_original_lossless_port_owner() {
    use fs_em::plane_wave::layers::{
        DielectricLayer, PlanarDielectricLayers,
        absorbing_incidence::PlanarAbsorbingDielectricLayers,
    };
    let incident = material(2.0 * E, 1.4 * M);
    let exit = material(3.0 * E, 0.9 * M);
    let dielectric = finite_layers::debye(2.0 * E, 3.0 * E, 0.7, 0.002);
    let layers = [DielectricLayer {
        response: DielectricResponse::Debye(&dielectric),
        thickness_m: 0.08,
    }];
    let omega = 1e9;
    let cosine = 0.8_f64;
    let kx = omega * (2.0 * E * 1.4 * M).sqrt() * ((1.0 - cosine) * (1.0 + cosine)).sqrt();
    let old = PlanarDielectricLayers::try_new(
        &incident,
        DielectricResponse::Lossless(&exit),
        omega,
        &layers,
    )
    .unwrap();
    let stack = PlanarAbsorbingDielectricLayers::try_new(
        DielectricResponse::Lossless(&incident),
        DielectricResponse::Lossless(&exit),
        omega,
        kx,
        &layers,
    )
    .unwrap();
    let prepared = stack.prepare_two_sided().unwrap();
    let a = C64::new(0.7, 0.2);
    let b = C64::new(0.2, -0.4);
    for (channel, pol) in [(0, Polarization::Tm), (1, Polarization::Te)] {
        let prior = old.prepare_illumination(cosine, pol).unwrap();
        let mut heat = [0.0];
        let expected = prior.evaluate([a, b], &mut heat).unwrap();
        let mut left = [C64::ZERO; 2];
        left[channel] = a;
        let mut right = [C64::ZERO; 2];
        right[channel] = b;
        let actual = prepared.illuminate(left, right).unwrap();
        let response = actual.response();
        for face in 0..2 {
            complex_near(
                response.outgoing[face].electric_v_m[channel],
                expected.outgoing_tangential_e_v_m[face],
                2e-13,
            );
            complex_near(
                response.boundary[face].electric_v_m[channel],
                expected.boundary_tangential_e_v_m[face],
                2e-13,
            );
            let h = if channel == 0 {
                response.boundary[face].magnetic_a_m[1]
            } else {
                -response.boundary[face].magnetic_a_m[0]
            };
            complex_near(h, expected.boundary_conjugate_h_a_m[face], 2e-15);
            assert_eq!(response.power[face].interference_into_stack_w_m2, 0.0);
        }
        near(actual.absorbed_power_per_layer_w_m2()[0], heat[0], 2e-15);
    }
}
