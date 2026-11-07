//! G1/G3: the public absorbing-incidence query retains original material laws.
use super::*;
use fs_material::conductor::{ElectricalPropertyBasis, ResolvedConductorState};

#[test]
fn original_debye_and_ohmic_sources_drive_absorbing_incidence_without_changing_old_ports() {
    let outside = material(4.0 * E, M);
    let dielectric = finite_layers::debye(2.0 * E, 3.0 * E, 0.7, 0.002);
    let sigma = 0.02;
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
