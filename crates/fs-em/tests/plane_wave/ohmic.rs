use super::*;
use fs_em::plane_wave::{
    angular_spectrum::{
        AngularSpectrumGrid, AngularSpectrumLayerLimits, AngularSpectrumLimits,
        PlanarAngularSpectrum,
    },
    layers::{DielectricLayer, PlanarDielectricLayers},
};
use fs_material::conductor::{ElectricalPropertyBasis, ResolvedConductorState};
use fs_math::det;

fn conductor(
    epsilon: f64,
    mu: f64,
    sigma: f64,
    basis: ElectricalPropertyBasis,
) -> ResolvedConductorState {
    let value = match basis {
        ElectricalPropertyBasis::Conductivity => sigma,
        ElectricalPropertyBasis::Resistivity => sigma.recip(),
    };
    ResolvedConductorState::try_from_material_with_basis(
        source(
            &[
                (EPS_NAME, EPS, epsilon),
                (MU_NAME, MU, mu),
                (basis.property_name(), basis.dimensions(), value),
            ],
            ValidityDomain::unconstrained().with("T", 280.0, 320.0),
            QueryPoint::new().with("T", 300.0).unwrap(),
        ),
        basis,
    )
    .unwrap()
}

fn relative(actual: f64, expected: f64, tolerance: f64) {
    near(
        actual,
        expected,
        tolerance * actual.abs().max(expected.abs()).max(1e-300),
    );
}

#[test]
fn ohmic_half_space_preserves_original_basis_skin_depth_and_joule_power() {
    let incident = material(E, M);
    for sigma in [0.02, 5.8e7] {
        for basis in [
            ElectricalPropertyBasis::Resistivity,
            ElectricalPropertyBasis::Conductivity,
        ] {
            let metal = conductor(3.0 * E, 1.4 * M, sigma, basis);
            let owner =
                PlanarDielectricInterface::try_new(&incident, DielectricResponse::Ohmic(&metal))
                    .unwrap();
            assert!(std::ptr::eq(owner.transmitted_material(), metal.material()));
            for omega in [1e6, 1e9] {
                let epsilon = C64::new(3.0 * E, -sigma / omega);
                for cosine in [1.0, 0.7, 0.03] {
                    let ki = omega * det::sqrt(E * M);
                    let kx = ki * det::sqrt((1.0 - cosine) * (1.0 + cosine));
                    // Independent dimensional dispersion relation; no normalized
                    // material contrast or production boundary algebra is reused.
                    let kz =
                        (epsilon.scale(omega * omega * 1.4 * M) - C64::from_re(kx * kx)).sqrt();
                    for pol in [Polarization::Te, Polarization::Tm] {
                        let actual = owner.response(omega, cosine, pol).unwrap();
                        let (yi, yt) = match pol {
                            Polarization::Te => {
                                (ki * cosine / (omega * M), kz.scale(1.0 / (omega * 1.4 * M)))
                            }
                            Polarization::Tm => {
                                (omega * E / (ki * cosine), epsilon.scale(omega) / kz)
                            }
                        };
                        let r = (C64::from_re(yi) - yt) / (C64::from_re(yi) + yt);
                        complex_near(actual.reflected_tangential_e, r, 3e-14);
                        complex_near(
                            actual.transmitted_permittivity_f_m,
                            epsilon,
                            2e-15 * epsilon.abs(),
                        );
                        complex_near(
                            actual.transmitted_normal_wavenumber_rad_m,
                            kz,
                            3e-15 * kz.abs(),
                        );
                        let et = actual.transmitted_tangential_e;
                        let ez = match pol {
                            Polarization::Te => C64::ZERO,
                            Polarization::Tm => -(et.scale(kx) / kz),
                        };
                        // Integral_0^infinity sigma |E(z)|^2/2 dz. The TM
                        // longitudinal field matters for oblique lossy incidence.
                        let heat = sigma * (et.norm_sq() + ez.norm_sq()) / (-4.0 * kz.im);
                        relative(heat / (0.5 * yi), actual.transmitted_power_fraction, 3e-13);
                        near(
                            actual.reflected_power_fraction + actual.transmitted_power_fraction,
                            1.0,
                            2e-14,
                        );
                    }
                    if cosine == 1.0 && sigma == 5.8e7 {
                        let skin = det::sqrt(2.0 / (omega * 1.4 * M * sigma));
                        relative(-1.0 / kz.im, skin, 3e-10);
                    }
                }
            }
        }
    }
}

#[test]
fn ohmic_films_recover_resistive_sheet_and_opaque_conductor_limits() {
    let incident = material(E, M);
    let sigma = 1e6;
    let metal = conductor(E, M, sigma, ElectricalPropertyBasis::Conductivity);
    let omega = 1e6;
    let wave_impedance = det::sqrt(M / E);
    // A freestanding sheet of conductance 2/Z has R=T=1/4 and A=1/2.
    // Resolve it as a finite bulk layer, without inserting a surface model.
    let thickness = 2.0 / (wave_impedance * sigma);
    for pol in [Polarization::Te, Polarization::Tm] {
        let sheet = PlanarDielectricLayers::try_new(
            &incident,
            DielectricResponse::Lossless(&incident),
            omega,
            &[DielectricLayer {
                response: DielectricResponse::Ohmic(&metal),
                thickness_m: thickness,
            }],
        )
        .unwrap();
        assert!(std::ptr::eq(
            sheet.layer(0).unwrap().response.material(),
            metal.material()
        ));
        let r = sheet.response(1.0, pol).unwrap();
        complex_near(r.reflected_tangential_e, C64::from_re(-0.5), 5e-11);
        complex_near(r.transmitted_tangential_e, C64::from_re(0.5), 5e-11);
        near(r.reflected_power_fraction, 0.25, 3e-14);
        near(r.transmitted_power_fraction, 0.25, 3e-14);
        near(r.absorbed_power_fraction, 0.5, 3e-14);

        let skin = det::sqrt(2.0 / (omega * M * sigma));
        let opaque = PlanarDielectricLayers::try_new(
            &incident,
            DielectricResponse::Lossless(&incident),
            omega,
            &[DielectricLayer {
                response: DielectricResponse::Ohmic(&metal),
                thickness_m: 1000.0 * skin,
            }],
        )
        .unwrap()
        .response(0.7, pol)
        .unwrap();
        let half = PlanarDielectricInterface::try_new(&incident, DielectricResponse::Ohmic(&metal))
            .unwrap()
            .response(omega, 0.7, pol)
            .unwrap();
        complex_near(
            opaque.reflected_tangential_e,
            half.reflected_tangential_e,
            3e-14,
        );
        assert_eq!(opaque.transmitted_power_fraction, 0.0);
        relative(
            opaque.absorbed_power_fraction,
            half.transmitted_power_fraction,
            3e-13,
        );
    }
}

#[test]
fn reflected_vector_beam_accounts_for_the_same_local_ohmic_heat_once() {
    let incident = material(E, M);
    let sigma = 0.02;
    let metal = conductor(3.0 * E, M, sigma, ElectricalPropertyBasis::Resistivity);
    let thickness = 0.7;
    let stack = PlanarDielectricLayers::try_new(
        &incident,
        DielectricResponse::Lossless(&incident),
        1e9,
        &[DielectricLayer {
            response: DielectricResponse::Ohmic(&metal),
            thickness_m: thickness,
        }],
    )
    .unwrap();
    let grid = AngularSpectrumGrid {
        width: 4,
        height: 4,
        spacing_m: [2.0, 3.0],
    };
    let mut input = vec![[C64::ZERO; 2]; 16];
    input[2] = [C64::ONE, C64::new(0.2, -0.3)];
    input[11] = [C64::new(0.1, 0.8), C64::new(-0.4, 0.2)];
    let beam = PlanarAngularSpectrum::prepare(
        &stack,
        grid,
        &input,
        AngularSpectrumLimits {
            maximum_grid_points: 16,
            maximum_workspace_bytes: 1 << 20,
            maximum_fft_work: 1 << 16,
        },
        || false,
    )
    .unwrap();
    let layer = beam
        .prepare_layer(
            0,
            AngularSpectrumLayerLimits {
                maximum_layer_mode_pairs: 32,
                maximum_workspace_bytes: 1 << 20,
                maximum_fft_work: 1 << 16,
            },
            || false,
        )
        .unwrap();
    relative(
        layer.absorbed_power_w(),
        beam.power().film_absorbed_w,
        3e-12,
    );
    assert!(beam.power().balance_defect_w().abs() < 3e-12 * beam.power().incident_w);
    let left = layer.field(0.0, || false).unwrap();
    let right = layer.field(thickness, || false).unwrap();
    relative(
        left.field().normal_power_w() - right.field().normal_power_w(),
        layer.absorbed_power_w(),
        4e-12,
    );
    let intervals = 64;
    let dz = thickness / intervals as f64;
    let mut heat_integral = 0.0;
    for i in 0..=intervals {
        let field = layer.field(i as f64 * dz, || false).unwrap();
        for (e, &heat) in field
            .field()
            .electric_v_per_m()
            .iter()
            .zip(field.absorbed_power_density_w_m3())
        {
            relative(
                heat,
                0.5 * sigma * e.iter().map(|v| v.norm_sq()).sum::<f64>(),
                3e-15,
            );
        }
        let weight = if i == 0 || i == intervals {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        heat_integral += weight * field.absorbed_power_per_depth_w_m();
    }
    relative(heat_integral * dz / 3.0, layer.absorbed_power_w(), 3e-6);
}

#[test]
fn shared_debye_and_ohmic_sources_add_polarization_and_joule_heat_once() {
    let incident = material(E, M);
    let sigma = 0.02;
    let tau = 2e-9;
    let basis = ElectricalPropertyBasis::Conductivity;
    let shared = source(
        &[
            (EPS_NAME, EPS, 2.0 * E),
            (MU_NAME, MU, M),
            ("delta", EPS, E),
            ("tau", Time::DIMS, tau),
            (basis.property_name(), basis.dimensions(), sigma),
        ],
        ValidityDomain::unconstrained(),
        QueryPoint::new(),
    );
    let dielectric = ResolvedDebyeCapacitor::try_from_material(
        shared.clone(),
        EPS_NAME,
        &[DebyePoleProperties {
            permittivity_increment: "delta".into(),
            relaxation_time: "tau".into(),
        }],
        0.03,
        0.002,
    )
    .unwrap();
    let conductor = ResolvedConductorState::try_from_material(shared).unwrap();
    let response = DielectricResponse::DebyeWithOhmic {
        dielectric: &dielectric,
        conductor: &conductor,
    };
    let interface = PlanarDielectricInterface::try_new(&incident, response).unwrap();
    assert_eq!(
        response.material().identity(),
        conductor.material().identity()
    );
    for omega in [1e8, 5e8, 1e9] {
        let wt = omega * tau;
        let expected_epsilon = C64::new(
            2.0 * E + E / (1.0 + wt * wt),
            -E * wt / (1.0 + wt * wt) - sigma / omega,
        );
        let impedance = (C64::from_re(M) / expected_epsilon).sqrt();
        let incident_impedance = C64::from_re(det::sqrt(M / E));
        let expected_reflection =
            (impedance - incident_impedance) / (impedance + incident_impedance);
        for pol in [Polarization::Te, Polarization::Tm] {
            let actual = interface.response(omega, 1.0, pol).unwrap();
            complex_near(
                actual.transmitted_permittivity_f_m,
                expected_epsilon,
                3e-15 * expected_epsilon.abs(),
            );
            complex_near(actual.reflected_tangential_e, expected_reflection, 3e-14);
        }
        let thickness = 0.15 / (omega * det::sqrt(E * M));
        let stack = PlanarDielectricLayers::try_new(
            &incident,
            DielectricResponse::Lossless(&incident),
            omega,
            &[DielectricLayer {
                response,
                thickness_m: thickness,
            }],
        )
        .unwrap();
        let field = stack.prepare_layer_field(0, 1.0, Polarization::Te).unwrap();
        // Independent cycle-mean Debye resistor loss plus the original Ohmic
        // heat. Neither term is inferred from subtracting external powers.
        let heat_coefficient = 0.5 * (sigma + omega * E * wt / (1.0 + wt * wt));
        let intervals = 64;
        let dz = thickness / intervals as f64;
        let mut integral = 0.0;
        for i in 0..=intervals {
            let e = field.field(i as f64 * dz).unwrap().electric_v_per_m;
            let weight = if i == 0 || i == intervals {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            integral += weight * heat_coefficient * e.norm_sq();
        }
        relative(integral * dz / 3.0, field.absorbed_power_w_m2(), 3e-7);
    }

    // The same nominal conductivity from a different retained bundle must
    // not borrow the Debye owner's source identity or material conditions.
    let foreign = self::conductor(2.0 * E, M, sigma, basis);
    let mismatched = DielectricResponse::DebyeWithOhmic {
        dielectric: &dielectric,
        conductor: &foreign,
    };
    assert!(PlanarDielectricInterface::try_new(&incident, mismatched).is_err());
    assert!(
        PlanarDielectricLayers::try_new(
            &incident,
            DielectricResponse::Lossless(&incident),
            1e9,
            &[DielectricLayer {
                response: mismatched,
                thickness_m: 0.0
            }]
        )
        .is_err()
    );
}

#[test]
fn ohmic_wave_source_requires_displacement_and_magnetic_data() {
    let incident = material(E, M);
    let basis = ElectricalPropertyBasis::Conductivity;
    for missing in [EPS_NAME, MU_NAME] {
        let properties: Vec<_> = [
            (EPS_NAME, EPS, E),
            (MU_NAME, MU, M),
            (basis.property_name(), basis.dimensions(), 0.02),
        ]
        .into_iter()
        .filter(|&(name, _, _)| name != missing)
        .collect();
        let incomplete = ResolvedConductorState::try_from_material(source(
            &properties,
            ValidityDomain::unconstrained(),
            QueryPoint::new(),
        ))
        .unwrap();
        assert!(
            PlanarDielectricInterface::try_new(&incident, DielectricResponse::Ohmic(&incomplete))
                .is_err()
        );
        // Zero thickness cannot bypass source admission in the existing stack.
        assert!(
            PlanarDielectricLayers::try_new(
                &incident,
                DielectricResponse::Lossless(&incident),
                1e9,
                &[DielectricLayer {
                    response: DielectricResponse::Ohmic(&incomplete),
                    thickness_m: 0.0
                }]
            )
            .is_err()
        );
    }
    let metal = conductor(E, M, 0.02, basis);
    let interface =
        PlanarDielectricInterface::try_new(&incident, DielectricResponse::Ohmic(&metal)).unwrap();
    for omega in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::from_bits(1)] {
        assert!(interface.response(omega, 1.0, Polarization::Te).is_err());
    }
}
