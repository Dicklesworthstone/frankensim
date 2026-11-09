#![cfg(feature = "plane-wave")]
//! G1/G3 frequency derivatives against the original Maxwell solves.

use fs_em::plane_wave::{
    DielectricResponse, PlanarDielectricInterface, Polarization,
    angular_spectrum::{
        AngularSpectrumGrid, AngularSpectrumLimits, AngularSpectrumSide, PlanarAngularSpectrum,
    },
    dispersion::IncidenceFrequencyConstraint,
    layers::{DielectricLayer, PlanarDielectricLayers},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    conductor::{ElectricalPropertyBasis, ResolvedConductorState},
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS as EPS, ELECTRIC_PERMITTIVITY_PROPERTY as EPS_NAME,
        MAGNETIC_PERMEABILITY_DIMS as MU, MAGNETIC_PERMEABILITY_PROPERTY as MU_NAME,
        dielectric::{
            DebyePoleProperties, ResolvedDebyeCapacitor,
            lorentz::{LorentzPoleProperties, ResolvedLorentzMaterial},
        },
    },
    state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point,
    },
};
use fs_math::c64::C64;
use fs_qty::{
    Dims, QuantitySpec, Time,
    semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm},
};

const E: f64 = 1e-11;
const M: f64 = 1e-6;
const OMEGA0: f64 = 1e9;
const SIGMA: f64 = 0.02;

fn source(damping: f64) -> ResolvedMaterialStatePoint {
    source_with_storage(damping, 2.0 * E, M)
}

fn source_with_storage(damping: f64, epsilon: f64, mu: f64) -> ResolvedMaterialStatePoint {
    let angular = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Frequency(FrequencyConvention::Angular), ValueForm::Static,
    ));
    let conductivity = ElectricalPropertyBasis::Conductivity;
    let entries = [
        (EPS_NAME, QuantitySpec::dimensional(EPS), epsilon),
        (MU_NAME, QuantitySpec::dimensional(MU), mu),
        ("delta", QuantitySpec::dimensional(EPS), 4.0 * E),
        ("tau", QuantitySpec::dimensional(Time::DIMS), 1.0 / OMEGA0),
        ("omega0", angular, OMEGA0),
        ("gamma", QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])), damping),
        (conductivity.property_name(), QuantitySpec::dimensional(conductivity.dimensions()), SIGMA),
    ];
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, quantity, value) in entries {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::with_quantity(name, quantity),
            value: PropertyValue::Scalar { value, dims: quantity.dims() },
            validity: ValidityDomain::unconstrained(),
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic bulk dispersion fixture".into(),
                license: "test-only".into(), artifact: None,
            },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_with_quantity(
            name, quantity, ScalarAdmissibility::NonNegative,
        ).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic dielectric".into(), phase: "solid".into(),
        process: "frequency derivative fixture".into(), revision: 0,
    }, claims, vec![]).unwrap();
    resolve_material_state_point(
        &card, &QueryPoint::new(), &requirements, MaterialPropertySelection::SingleClaimOnly,
    ).unwrap()
}

fn debye(material: ResolvedMaterialStatePoint, area: f64, gap: f64) -> ResolvedDebyeCapacitor {
    ResolvedDebyeCapacitor::try_from_material(material, EPS_NAME, &[DebyePoleProperties {
        permittivity_increment: "delta".into(), relaxation_time: "tau".into(),
    }], area, gap).unwrap()
}
fn lorentz(material: ResolvedMaterialStatePoint) -> ResolvedLorentzMaterial {
    ResolvedLorentzMaterial::try_from_material(material, EPS_NAME, &[LorentzPoleProperties {
        permittivity_increment: "delta".into(), angular_frequency: "omega0".into(),
        damping_rate: "gamma".into(),
    }]).unwrap()
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual:.17e} != {expected:.17e}");
}
fn complex_near(actual: C64, expected: C64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual:?} != {expected:?}");
}

// Independent five-point differentiation of the original interface solution.
// No production frequency-derivative helper participates in this oracle.
fn native_derivatives(response: DielectricResponse<'_>, omega: f64) -> [C64; 3] {
    let ambient = source(0.03 * OMEGA0);
    let interface = PlanarDielectricInterface::try_new(&ambient, response).unwrap();
    let step = 1e-4 * omega;
    let values: [C64; 5] = core::array::from_fn(|index| {
        interface.response(omega + (index as f64 - 2.0) * step, 1.0, Polarization::Te)
            .unwrap().transmitted_normal_wavenumber_rad_m
    });
    let first = (values[0] - values[1].scale(8.0) + values[3].scale(8.0) - values[4])
        .scale(1.0 / (12.0 * step));
    let second = (-values[0] + values[1].scale(16.0) - values[2].scale(30.0)
        + values[3].scale(16.0) - values[4]).scale(1.0 / (12.0 * step * step));
    [values[2], first, second]
}

#[test]
fn g1_nondispersive_limit_keeps_phase_group_and_exact_zero_gvd() {
    let material = source(0.03 * OMEGA0);
    let response = DielectricResponse::Lossless(&material);
    let inverse_speed = (2.0 * E * M).sqrt();
    for omega in [1e3, OMEGA0, 1e16] {
        let jet = response.bulk_dispersion(omega).unwrap();
        assert!(std::ptr::eq(jet.source().material(), &material));
        assert_eq!(jet.scaled_permittivity_derivatives_f_m(), [C64::ZERO; 2]);
        assert_eq!(jet.wavenumber_second_derivative_s2_m(), C64::ZERO);
        near(jet.group_slowness_s_m(), inverse_speed, 2e-16 * inverse_speed);
        near(jet.group_velocity_m_s().unwrap(), inverse_speed.recip(), 4e-8);
        near(jet.phase_velocity_m_s().unwrap(), inverse_speed.recip(), 4e-8);
        let propagated = jet.propagation(0.17).unwrap();
        near(propagated.phase_delay_s, 0.17 * inverse_speed, 2e-25);
        near(propagated.group_delay_s, propagated.phase_delay_s, 0.0);
        assert_eq!(propagated.group_delay_dispersion_s2, 0.0);
        near(propagated.transfer.norm_sq(), 1.0, 5e-16);
        let identity = jet.propagation(0.0).unwrap();
        assert_eq!(identity.transfer, C64::ONE);
        assert_eq!(identity.group_delay_s, 0.0);
    }
}

#[test]
fn g1_debye_and_ohmic_jets_match_original_maxwell_frequency_sweeps() {
    let material = source(0.03 * OMEGA0);
    let polarization = debye(material.clone(), 0.02, 0.001);
    let different_geometry = debye(material.clone(), 0.31, 0.06);
    let conductor = ResolvedConductorState::try_from_material_with_basis(
        material, ElectricalPropertyBasis::Conductivity,
    ).unwrap();
    for response in [
        DielectricResponse::Debye(&polarization),
        DielectricResponse::Ohmic(&conductor),
        DielectricResponse::DebyeWithOhmic { dielectric: &polarization, conductor: &conductor },
    ] {
        for omega in [0.05 * OMEGA0, OMEGA0, 20.0 * OMEGA0] {
            let jet = response.bulk_dispersion(omega).unwrap();
            let [k, first, second] = native_derivatives(response, omega);
            complex_near(jet.wavenumber_rad_m(), k, 8e-15 * k.abs());
            complex_near(jet.wavenumber_derivative_s_m(), first, 1e-8 * first.abs());
            complex_near(jet.wavenumber_second_derivative_s2_m(), second, 2e-5 * second.abs());
        }
    }
    for omega in [0.01 * OMEGA0, OMEGA0, 100.0 * OMEGA0] {
        let original = DielectricResponse::Debye(&polarization).bulk_dispersion(omega).unwrap();
        let other = DielectricResponse::Debye(&different_geometry).bulk_dispersion(omega).unwrap();
        assert_eq!(original.scaled_permittivity_derivatives_f_m(), other.scaled_permittivity_derivatives_f_m());
        complex_near(original.wavenumber_derivative_s_m(), other.wavenumber_derivative_s_m(),
            5e-15 * original.wavenumber_derivative_s_m().abs());
        let combined = DielectricResponse::DebyeWithOhmic {
            dielectric: &polarization, conductor: &conductor,
        }.bulk_dispersion(omega).unwrap();
        let bare = original.scaled_permittivity_derivatives_f_m();
        let leaking = combined.scaled_permittivity_derivatives_f_m();
        complex_near(leaking[0] - bare[0], C64::new(0.0, SIGMA / omega), 3e-15 * SIGMA / omega);
        complex_near(leaking[1] - bare[1], C64::new(0.0, -2.0 * SIGMA / omega), 5e-15 * SIGMA / omega);
    }
}

#[test]
fn g1_lorentz_keeps_resonant_anomalous_and_evanescent_dispersion() {
    for damping in [0.0, 0.03 * OMEGA0, 2.0 * OMEGA0] {
        let material = source(damping);
        let polarization = lorentz(material.clone());
        let conductor = ResolvedConductorState::try_from_material_with_basis(
            material, ElectricalPropertyBasis::Conductivity,
        ).unwrap();
        for response in [
            DielectricResponse::Lorentz(&polarization),
            DielectricResponse::LorentzWithOhmic { dielectric: &polarization, conductor: &conductor },
        ] {
            for ratio in [0.2, 0.85, 1.0, 1.05, 1.3, 3.0] {
                if damping == 0.0 && ratio == 1.0 { continue; }
                let omega = ratio * OMEGA0;
                let jet = response.bulk_dispersion(omega).unwrap();
                let [k, first, second] = native_derivatives(response, omega);
                complex_near(jet.wavenumber_rad_m(), k, 1e-14 * k.abs());
                complex_near(jet.wavenumber_derivative_s_m(), first, 1e-7 * first.abs());
                complex_near(jet.wavenumber_second_derivative_s2_m(), second, 2e-5 * second.abs());
            }
        }
        if damping == 0.03 * OMEGA0 {
            let jet = DielectricResponse::Lorentz(&polarization).bulk_dispersion(OMEGA0).unwrap();
            assert!(jet.group_slowness_s_m() < 0.0);
            assert!(jet.group_velocity_m_s().unwrap() < 0.0);
            let propagation = jet.propagation(0.001).unwrap();
            assert!(propagation.group_delay_s < 0.0);
            assert!(propagation.transfer.norm_sq() < 1.0);
        }
        if damping == 0.0 {
            let jet = DielectricResponse::Lorentz(&polarization).bulk_dispersion(1.3 * OMEGA0).unwrap();
            assert_eq!(jet.wavenumber_rad_m().re, 0.0);
            assert!(jet.wavenumber_rad_m().im < 0.0);
            assert_eq!(jet.phase_velocity_m_s(), None);
            assert_eq!(jet.group_velocity_m_s(), None);
        }
    }
}

#[test]
fn g3_exact_carrier_transfer_reaches_original_maxwell_fields() {
    let ambient = source(0.05 * OMEGA0);
    let polarization = lorentz(ambient.clone());
    let relaxation = debye(ambient.clone(), 1.0, 1.0);
    let limits = AngularSpectrumLimits {
        maximum_grid_points: 4, maximum_workspace_bytes: 100_000, maximum_fft_work: 10_000,
    };
    let grid = AngularSpectrumGrid { width: 2, height: 2, spacing_m: [1.0, 1.0] };
    for response in [
        DielectricResponse::Debye(&relaxation), DielectricResponse::Lorentz(&polarization),
    ] {
        for ratio in [0.5, 1.0, 2.0] {
            let omega = ratio * OMEGA0;
            let jet = response.bulk_dispersion(omega).unwrap();
            let stack = PlanarDielectricLayers::try_new(&ambient, response, omega, &[]).unwrap();
            let electric = [C64::new(0.8, 0.1), C64::new(-0.3, 0.7)];
            let beam = PlanarAngularSpectrum::prepare(&stack, grid, &[electric; 4], limits, || false).unwrap();
            let face = beam.field(AngularSpectrumSide::Transmitted, 0.0, || false).unwrap();
            for distance in [0.01, 0.1, 1.0] {
                let transfer = jet.propagation(distance).unwrap().transfer;
                let actual = beam.field(AngularSpectrumSide::Transmitted, distance, || false).unwrap();
                for (initial, propagated) in face.electric_v_per_m().iter().zip(actual.electric_v_per_m()) {
                    for component in 0..3 {
                        complex_near(propagated[component], initial[component] * transfer, 2e-13);
                    }
                }
                for (initial, propagated) in face.magnetic_a_per_m().iter().zip(actual.magnetic_a_per_m()) {
                    for component in 0..3 {
                        complex_near(propagated[component], initial[component] * transfer, 2e-15);
                    }
                }
            }
        }
    }
}

#[test]
fn g0_refuses_invalid_frequencies_singular_poles_and_mismatched_sources() {
    let material = source(0.0);
    let polarization = lorentz(material.clone());
    let jet = DielectricResponse::Lossless(&material).bulk_dispersion(OMEGA0).unwrap();
    for invalid in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(jet.propagation(invalid).is_err());
        assert!(DielectricResponse::Lossless(&material).bulk_dispersion(invalid).is_err());
    }
    assert!(DielectricResponse::Lossless(&material).bulk_dispersion(0.0).is_err());
    assert!(DielectricResponse::Lorentz(&polarization).bulk_dispersion(OMEGA0).is_err());
    let other = ResolvedConductorState::try_from_material_with_basis(
        source(0.05 * OMEGA0), ElectricalPropertyBasis::Conductivity,
    ).unwrap();
    assert!(DielectricResponse::LorentzWithOhmic { dielectric: &polarization, conductor: &other }
        .bulk_dispersion(0.5 * OMEGA0).is_err());
}

#[test]
fn g1_high_frequency_loss_derivatives_keep_resolved_absorption_tails() {
    let material = source(0.03 * OMEGA0);
    let relaxation = debye(material.clone(), 1.0, 1.0);
    let conductor = ResolvedConductorState::try_from_material_with_basis(
        material, ElectricalPropertyBasis::Conductivity,
    ).unwrap();
    let inverse_speed = (2.0 * E * M).sqrt();
    for ratio in [1e6_f64, 1e8, 1e12, 1e50] {
        let omega = ratio * OMEGA0;
        let jet = DielectricResponse::Ohmic(&conductor).bulk_dispersion(omega).unwrap();
        // Independent closed form from k² = mu*(epsilon_inf*omega²-i*sigma*omega).
        // Writing sqrt(1-i*t) = u+i*v gives an imaginary first derivative
        // proportional to t³, even when sqrt(1+t²) rounds to exactly one.
        let t = SIGMA / (2.0 * E * omega);
        let radius = (1.0 + t * t).sqrt();
        let u = ((radius + 1.0) * 0.5).sqrt();
        let first = -inverse_speed * t.powi(3) / (8.0 * u.powi(3) * radius);
        let second = (inverse_speed / omega) * t.powi(3) * (radius + 2.0)
            / (8.0 * u * radius.powi(3));
        near(jet.wavenumber_derivative_s_m().im, first, 5e-14 * first.abs());
        near(jet.wavenumber_second_derivative_s2_m().im, second, 5e-14 * second.abs());
        let propagated = jet.propagation(0.1).unwrap();
        near(propagated.log_amplitude_derivative_s, 0.1 * first, 5e-15 * first.abs());
        near(propagated.log_amplitude_second_derivative_s2, 0.1 * second, 5e-15 * second.abs());

        // Independently expand sqrt(epsilon/epsilon_inf) through omega^-3.
        // With a=(sigma+Delta/tau)/epsilon_inf, b=Delta/(epsilon_inf*tau²),
        // c=Delta/(epsilon_inf*tau³), Q=c+a*b/2+a³/8, Im(k')~-S*Q/omega³
        // and Im(k'')~3*S*Q/omega⁴. Here Delta/epsilon_inf=2, tau*OMEGA0=1:
        // Q/OMEGA0³ = 5 (Debye) and 67/8 (Debye with this source's leakage).
        for (response, coefficient) in [
            (DielectricResponse::Debye(&relaxation), 5.0),
            (DielectricResponse::DebyeWithOhmic {
                dielectric: &relaxation, conductor: &conductor,
            }, 67.0 / 8.0),
        ] {
            let jet = response.bulk_dispersion(omega).unwrap();
            let first = -inverse_speed * (coefficient / ratio) / ratio / ratio;
            let second = -3.0 * first / omega;
            // The omitted relative tail is O(ratio^-2); compare each small
            // imaginary component, not a norm dominated by phase slowness.
            near(jet.wavenumber_derivative_s_m().im, first, 1e-10 * first.abs());
            near(jet.wavenumber_second_derivative_s2_m().im, second, 1e-10 * second.abs());
        }
    }
}

fn differentiate_samples(values: [C64; 5], step: f64) -> [C64; 2] {
    [
        (values[0] - values[1].scale(8.0) + values[3].scale(8.0) - values[4])
            .scale(1.0 / (12.0 * step)),
        (-values[0] + values[1].scale(16.0) - values[2].scale(30.0)
            + values[3].scale(16.0) - values[4]).scale(1.0 / (12.0 * step * step)),
    ]
}

#[test]
fn g1_fresnel_frequency_jets_follow_original_dispersive_boundary_fields() {
    use IncidenceFrequencyConstraint::{FixedAngle, FixedTransverseWavenumber};
    let ambient = source(0.03 * OMEGA0);
    let relaxation = debye(ambient.clone(), 0.2, 0.01);
    let resonance = lorentz(ambient.clone());
    let clear_resonance = lorentz(source(0.0));
    let conductor = ResolvedConductorState::try_from_material_with_basis(
        ambient.clone(), ElectricalPropertyBasis::Conductivity,
    ).unwrap();
    for incident in [DielectricResponse::Lossless(&ambient), DielectricResponse::Lorentz(&clear_resonance)] {
        for transmitted in [
            DielectricResponse::Debye(&relaxation), DielectricResponse::Lorentz(&resonance),
            DielectricResponse::Ohmic(&conductor),
            DielectricResponse::LorentzWithOhmic { dielectric: &resonance, conductor: &conductor },
        ] {
            let interface = PlanarDielectricInterface::try_new_with_incident_response(incident, transmitted).unwrap();
            for omega in [0.3 * OMEGA0, 0.8 * OMEGA0] {
                for cosine in [0.5_f64, 0.9, 1.0] {
                    let ki = interface.response(omega, 1.0, Polarization::Te).unwrap()
                        .incident_normal_wavenumber_rad_m;
                    let k_parallel = ki * ((1.0 - cosine) * (1.0 + cosine)).sqrt();
                    for constraint in [FixedAngle, FixedTransverseWavenumber] {
                        for polarization in [Polarization::Te, Polarization::Tm] {
                            let jet = interface.frequency_dispersion(omega, cosine, polarization, constraint).unwrap();
                            assert_eq!(jet.carrier(), interface.response(omega, cosine, polarization).unwrap());
                            assert!(std::ptr::eq(jet.incident_dispersion().source().material(), incident.material()));
                            assert!(std::ptr::eq(jet.transmitted_dispersion().source().material(), transmitted.material()));
                            let step = 1e-4 * omega;
                            let values: [[C64; 3]; 5] = core::array::from_fn(|index| {
                                let frequency = omega + (index as f64 - 2.0) * step;
                                let c = if constraint == FixedAngle { cosine } else {
                                    let wave_number = interface.response(frequency, 1.0, polarization)
                                        .unwrap().incident_normal_wavenumber_rad_m;
                                    (1.0 - (k_parallel / wave_number).powi(2)).sqrt()
                                };
                                let carrier = interface.response(frequency, c, polarization).unwrap();
                                [carrier.reflected_tangential_e, carrier.transmitted_tangential_e,
                                    carrier.transmitted_conjugate_h]
                            });
                            for (component, amplitude) in [jet.reflected_tangential_e(),
                                jet.transmitted_tangential_e(), jet.transmitted_conjugate_h()].into_iter().enumerate()
                            {
                                let [first, second] = differentiate_samples(values.map(|value| value[component]), step);
                                complex_near(amplitude.derivative_s(), first, 2e-8 / omega);
                                complex_near(amplitude.second_derivative_s2(), second, 3e-6 / omega / omega);
                            }
                            let reflected = jet.reflected_tangential_e();
                            let phase = reflected.logarithmic_dispersion().unwrap();
                            // Independent phase samples referenced to the carrier avoid a 2*pi wrap.
                            let phases = values.map(|value| {
                                let ratio = value[0] / values[2][0];
                                C64::from_re(ratio.im.atan2(ratio.re))
                            });
                            let [first, second] = differentiate_samples(phases, step);
                            near(phase.group_delay_s, -first.re, 3e-8 / omega);
                            near(phase.group_delay_dispersion_s2, -second.re, 8e-6 / omega / omega);
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn g1_total_reflection_delay_distinguishes_angle_and_transverse_frequency_paths() {
    use IncidenceFrequencyConstraint::{FixedAngle, FixedTransverseWavenumber};
    let incident = source_with_storage(0.0, 1.0, 1.0);
    let transmitted = source_with_storage(0.0, 0.25, 1.0);
    let interface = PlanarDielectricInterface::try_new(&incident, DielectricResponse::Lossless(&transmitted)).unwrap();
    let omega = OMEGA0;
    let cosine = 0.5_f64;
    let sine_squared = 1.0 - cosine * cosine;
    let kappa = (sine_squared - 0.25).sqrt();
    for polarization in [Polarization::Te, Polarization::Tm] {
        let fixed_angle = interface.frequency_dispersion(omega, cosine, polarization, FixedAngle).unwrap();
        assert_eq!(fixed_angle.reflected_tangential_e().scaled_derivatives(), [C64::ZERO; 2]);
        assert_eq!(fixed_angle.reflected_tangential_e().logarithmic_dispersion().unwrap().group_delay_s, 0.0);
    }
    let fixed_k = interface.frequency_dispersion(omega, cosine, Polarization::Te, FixedTransverseWavenumber).unwrap();
    let phase = fixed_k.reflected_tangential_e().logarithmic_dispersion().unwrap();
    // Differentiate the independent TIR phase 2 atan(kappa/cosine) at fixed k_parallel.
    let delay = 2.0 * sine_squared / (omega * cosine * kappa);
    let gdd = delay / omega * (-1.0 - 1.0 / cosine.powi(2) + 0.25 / kappa.powi(2));
    near(phase.group_delay_s, delay, 5e-15 * delay);
    near(phase.group_delay_dispersion_s2, gdd, 2e-14 * gdd.abs());
    near(phase.log_amplitude_derivative_s, 0.0, 2e-15 / omega);
    near(phase.log_amplitude_second_derivative_s2, 0.0, 2e-14 / omega / omega);
}

#[test]
fn g0_fresnel_zero_has_amplitude_derivatives_but_no_invented_phase() {
    use IncidenceFrequencyConstraint::{FixedAngle, FixedTransverseWavenumber};
    let incident = source_with_storage(0.0, 2.0, 1.0);
    let matched = PlanarDielectricInterface::try_new(&incident, DielectricResponse::Lossless(&incident)).unwrap();
    let response = matched.frequency_dispersion(OMEGA0, 0.5, Polarization::Te, FixedAngle).unwrap();
    assert_eq!(response.reflected_tangential_e().amplitude(), C64::ZERO);
    assert!(response.reflected_tangential_e().logarithmic_dispersion().is_none());
    assert_eq!(response.transmitted_tangential_e().amplitude(), C64::ONE);
    assert_eq!(response.transmitted_tangential_e().scaled_derivatives(), [C64::ZERO; 2]);

    // Binary-exact Brewster point: epsilon_t/epsilon_i=2, mu_t/mu_i=7/8,
    // cosine=1/2 gives q_t=1 and the two TM Fresnel numerators both one.
    let transmitted = source_with_storage(0.0, 4.0, 0.875);
    let brewster = PlanarDielectricInterface::try_new(&incident, DielectricResponse::Lossless(&transmitted)).unwrap();
    let jet = brewster.frequency_dispersion(OMEGA0, 0.5, Polarization::Tm, FixedTransverseWavenumber).unwrap();
    assert_eq!(jet.reflected_tangential_e().amplitude(), C64::ZERO);
    assert!(jet.reflected_tangential_e().logarithmic_dispersion().is_none());
    assert!(jet.reflected_tangential_e().scaled_derivatives()[0].abs() > 0.1);

    let critical_exit = source_with_storage(0.0, 1.5, 1.0);
    let critical = PlanarDielectricInterface::try_new(&incident, DielectricResponse::Lossless(&critical_exit)).unwrap();
    assert!(critical.response(OMEGA0, 0.5, Polarization::Te).is_ok());
    assert!(critical.frequency_dispersion(OMEGA0, 0.5, Polarization::Te, FixedTransverseWavenumber).is_err());
    for frequency in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matched.frequency_dispersion(frequency, 0.5, Polarization::Te, FixedAngle).is_err());
    }
}

#[test]
fn g3_fixed_transverse_fresnel_jet_reaches_the_original_reflected_beam() {
    let ambient = source(0.03 * OMEGA0);
    let resonance = lorentz(ambient.clone());
    let transmitted = DielectricResponse::Lorentz(&resonance);
    let interface = PlanarDielectricInterface::try_new(&ambient, transmitted).unwrap();
    let omega = 0.7 * OMEGA0;
    let cosine = 0.8_f64;
    let ki = interface.response(omega, 1.0, Polarization::Te).unwrap().incident_normal_wavenumber_rad_m;
    let carrier = [ki * (1.0 - cosine * cosine).sqrt(), 0.0];
    let limits = AngularSpectrumLimits {
        maximum_grid_points: 4, maximum_workspace_bytes: 100_000, maximum_fft_work: 10_000,
    };
    let grid = AngularSpectrumGrid { width: 2, height: 2, spacing_m: [1.0, 1.0] };
    let step = omega * 1e-4;
    let fields = core::array::from_fn(|index| {
        let frequency = omega + (index as f64 - 2.0) * step;
        let stack = PlanarDielectricLayers::try_new(&ambient, transmitted, frequency, &[]).unwrap();
        let beam = PlanarAngularSpectrum::prepare_with_carrier(
            &stack, grid, carrier, &[[C64::ZERO, C64::ONE]; 4], limits, || false,
        ).unwrap();
        beam.field(AngularSpectrumSide::Reflected, 0.0, || false).unwrap().electric_v_per_m()[0][1]
    });
    let jet = interface.frequency_dispersion(omega, cosine, Polarization::Te,
        IncidenceFrequencyConstraint::FixedTransverseWavenumber).unwrap().reflected_tangential_e();
    let [first, second] = differentiate_samples(fields, step);
    complex_near(jet.amplitude(), fields[2], 2e-15);
    complex_near(jet.derivative_s(), first, 2e-8 / omega);
    complex_near(jet.second_derivative_s2(), second, 3e-6 / omega / omega);
}

fn check_layer_frequency_sweep(
    incident: DielectricResponse<'_>,
    exit: DielectricResponse<'_>,
    layers: &[DielectricLayer<'_>],
    omega: f64,
    cosine: f64,
    polarization: Polarization,
    constraint: IncidenceFrequencyConstraint,
) {
    let stack =
        PlanarDielectricLayers::try_new_with_incident_response(incident, exit, omega, layers)
            .unwrap();
    let jet = stack
        .second_order_frequency_dispersion(cosine, polarization, constraint)
        .unwrap();
    let first_order = stack
        .frequency_dispersion(cosine, polarization, constraint)
        .unwrap();
    assert!(std::ptr::eq(jet.source(), &stack));
    assert_eq!(jet.carrier(), stack.response(cosine, polarization).unwrap());
    let reference =
        PlanarDielectricInterface::try_new_with_incident_response(incident, exit).unwrap();
    let ki = reference
        .response(omega, 1.0, polarization)
        .unwrap()
        .incident_normal_wavenumber_rad_m;
    let k_parallel = ki * ((1.0 - cosine) * (1.0 + cosine)).sqrt();
    let step = omega * 1e-4;
    let fields: [[C64; 3]; 5] = core::array::from_fn(|index| {
        let frequency = omega + (index as f64 - 2.0) * step;
        let c = if constraint == IncidenceFrequencyConstraint::FixedAngle {
            cosine
        } else {
            let k = reference
                .response(frequency, 1.0, polarization)
                .unwrap()
                .incident_normal_wavenumber_rad_m;
            (1.0 - (k_parallel / k).powi(2)).sqrt()
        };
        let sampled = PlanarDielectricLayers::try_new_with_incident_response(
            incident, exit, frequency, layers,
        )
        .unwrap()
        .response(c, polarization)
        .unwrap();
        [
            sampled.reflected_tangential_e,
            sampled.transmitted_tangential_e,
            sampled.transmitted_conjugate_h,
        ]
    });
    for (index, amplitude) in [
        jet.reflected_tangential_e(),
        jet.transmitted_tangential_e().unwrap(),
        jet.transmitted_conjugate_h().unwrap(),
    ]
    .into_iter()
    .enumerate()
    {
        let [first, second] = differentiate_samples(fields.map(|value| value[index]), step);
        let tangent = [
            first_order.reflected_tangential_e(),
            first_order.transmitted_tangential_e().unwrap(),
            first_order.transmitted_conjugate_h().unwrap(),
        ][index];
        assert_eq!(amplitude.amplitude(), tangent.amplitude());
        complex_near(
            amplitude.scaled_derivatives()[0],
            tangent.scaled_derivative(),
            8e-12 * tangent.scaled_derivative().abs().max(1.0),
        );
        complex_near(amplitude.derivative_s(), first, 8e-8 / omega);
        complex_near(
            amplitude.second_derivative_s2(),
            second,
            8e-6 / omega / omega,
        );
    }
}

#[test]
fn g1_coherent_film_jets_recover_native_bulk_propagation_and_half_space() {
    use IncidenceFrequencyConstraint::{FixedAngle, FixedTransverseWavenumber};
    let material = source(0.0);
    let resonance = lorentz(material.clone());
    let relaxation = debye(material.clone(), 1.0, 1.0);
    let omega = 0.6 * OMEGA0;
    for response in [
        DielectricResponse::Lossless(&material),
        DielectricResponse::Lorentz(&resonance),
    ] {
        let layers = [
            DielectricLayer {
                response,
                thickness_m: 0.04,
            },
            DielectricLayer {
                response,
                thickness_m: 0.06,
            },
        ];
        let stack = PlanarDielectricLayers::try_new_with_incident_response(
            response, response, omega, &layers,
        )
        .unwrap();
        let bulk = response
            .bulk_dispersion(omega)
            .unwrap()
            .propagation(0.1)
            .unwrap();
        for constraint in [FixedAngle, FixedTransverseWavenumber] {
            for polarization in [Polarization::Te, Polarization::Tm] {
                let jet = stack
                    .second_order_frequency_dispersion(1.0, polarization, constraint)
                    .unwrap();
                let amplitude = jet.transmitted_tangential_e().unwrap();
                let phase = amplitude.logarithmic_dispersion().unwrap();
                complex_near(amplitude.amplitude(), bulk.transfer, 2e-14);
                near(phase.group_delay_s, bulk.group_delay_s, 2e-13 / omega);
                near(
                    phase.group_delay_dispersion_s2,
                    bulk.group_delay_dispersion_s2,
                    2e-12 / omega / omega,
                );
                for value in jet.reflected_tangential_e().scaled_derivatives() {
                    complex_near(value, C64::ZERO, 2e-13);
                }
            }
        }
        let exit = DielectricResponse::Debye(&relaxation);
        let empty =
            PlanarDielectricLayers::try_new_with_incident_response(response, exit, omega, &[])
                .unwrap();
        let interface =
            PlanarDielectricInterface::try_new_with_incident_response(response, exit).unwrap();
        for polarization in [Polarization::Te, Polarization::Tm] {
            for constraint in [FixedAngle, FixedTransverseWavenumber] {
                let actual = empty
                    .second_order_frequency_dispersion(0.7, polarization, constraint)
                    .unwrap();
                let expected = interface
                    .frequency_dispersion(omega, 0.7, polarization, constraint)
                    .unwrap();
                for (a, b) in [
                    (
                        actual.reflected_tangential_e(),
                        expected.reflected_tangential_e(),
                    ),
                    (
                        actual.transmitted_tangential_e().unwrap(),
                        expected.transmitted_tangential_e(),
                    ),
                    (
                        actual.transmitted_conjugate_h().unwrap(),
                        expected.transmitted_conjugate_h(),
                    ),
                ] {
                    for (a, b) in a
                        .scaled_derivatives()
                        .into_iter()
                        .zip(b.scaled_derivatives())
                    {
                        complex_near(a, b, 4e-13);
                    }
                }
            }
        }
    }
}

#[test]
fn g1_resonant_etalon_delay_keeps_multiple_reflection_dwell_time() {
    let ambient = source(0.0);
    let omega = 0.7 * OMEGA0;
    let ki = omega * (2.0 * E * M).sqrt();
    for n in [2.0_f64, 8.0, 32.0] {
        let film = source_with_storage(0.0, 2.0 * E * n * n, M);
        let thickness = std::f64::consts::PI / (n * ki);
        let stack = PlanarDielectricLayers::try_new(
            &ambient,
            DielectricResponse::Lossless(&ambient),
            omega,
            &[DielectricLayer {
                response: DielectricResponse::Lossless(&film),
                thickness_m: thickness,
            }],
        )
        .unwrap();
        // Independent Airy-series phase slope at the first half-wave resonance:
        // tau = (n + 1/n)/2 times the single-pass phase time.
        let expected = 0.5 * (n + n.recip()) * std::f64::consts::PI / omega;
        for polarization in [Polarization::Te, Polarization::Tm] {
            let jet = stack
                .second_order_frequency_dispersion(
                    1.0,
                    polarization,
                    IncidenceFrequencyConstraint::FixedAngle,
                )
                .unwrap();
            let phase = jet
                .transmitted_tangential_e()
                .unwrap()
                .logarithmic_dispersion()
                .unwrap();
            near(jet.carrier().transmitted_power_fraction, 1.0, 4e-13);
            near(phase.group_delay_s, expected, 3e-12 * expected);
            near(phase.group_delay_dispersion_s2, 0.0, 2e-10 / omega / omega);
        }
    }
}

#[test]
fn g1_internal_critical_film_has_regular_delay_and_gdd_on_both_frequency_paths() {
    let ambient = source_with_storage(0.0, 1.0, 1.0);
    let critical = source_with_storage(0.0, 0.75, 1.0);
    let exit = source_with_storage(0.0, 2.25, 1.0);
    let layers = [DielectricLayer {
        response: DielectricResponse::Lossless(&critical),
        thickness_m: 1.3 / OMEGA0,
    }];
    for constraint in [
        IncidenceFrequencyConstraint::FixedAngle,
        IncidenceFrequencyConstraint::FixedTransverseWavenumber,
    ] {
        for polarization in [Polarization::Te, Polarization::Tm] {
            check_layer_frequency_sweep(
                DielectricResponse::Lossless(&ambient),
                DielectricResponse::Lossless(&exit),
                &layers,
                OMEGA0,
                0.5,
                polarization,
                constraint,
            );
        }
    }
    let critical_exit = PlanarDielectricLayers::try_new(
        &ambient,
        DielectricResponse::Lossless(&critical),
        OMEGA0,
        &[],
    )
    .unwrap();
    assert!(critical_exit.response(0.5, Polarization::Te).is_ok());
    assert!(
        critical_exit
            .second_order_frequency_dispersion(
                0.5,
                Polarization::Te,
                IncidenceFrequencyConstraint::FixedTransverseWavenumber
            )
            .is_err()
    );
}

#[test]
fn g1_dispersive_absorbing_coating_jets_follow_full_native_spectral_sweeps() {
    let material = source(0.03 * OMEGA0);
    let clear = lorentz(source(0.0));
    let resonance = lorentz(material.clone());
    let relaxation = debye(material.clone(), 1.0, 1.0);
    let conductor = ResolvedConductorState::try_from_material_with_basis(
        material.clone(),
        ElectricalPropertyBasis::Conductivity,
    )
    .unwrap();
    let layers = [
        DielectricLayer {
            response: DielectricResponse::Debye(&relaxation),
            thickness_m: 0.08,
        },
        DielectricLayer {
            response: DielectricResponse::Lorentz(&resonance),
            thickness_m: 0.17,
        },
        DielectricLayer {
            response: DielectricResponse::Ohmic(&conductor),
            thickness_m: 0.03,
        },
    ];
    for incident in [
        DielectricResponse::Lossless(&material),
        DielectricResponse::Lorentz(&clear),
    ] {
        for omega in [0.4 * OMEGA0, 0.8 * OMEGA0] {
            for cosine in [0.6, 1.0] {
                for polarization in [Polarization::Te, Polarization::Tm] {
                    for constraint in [
                        IncidenceFrequencyConstraint::FixedAngle,
                        IncidenceFrequencyConstraint::FixedTransverseWavenumber,
                    ] {
                        check_layer_frequency_sweep(
                            incident,
                            DielectricResponse::Debye(&relaxation),
                            &layers,
                            omega,
                            cosine,
                            polarization,
                            constraint,
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn g3_film_frequency_jets_preserve_subdivision_and_opaque_limits() {
    let ambient = source(0.03 * OMEGA0);
    let resonance = lorentz(ambient.clone());
    let conductor = ResolvedConductorState::try_from_material_with_basis(
        ambient.clone(),
        ElectricalPropertyBasis::Conductivity,
    )
    .unwrap();
    for (response, thickness) in [
        (DielectricResponse::Lorentz(&resonance), 0.3),
        // Each half remains representable; their product underflows.
        (DielectricResponse::Ohmic(&conductor), 400.0),
        (DielectricResponse::Ohmic(&conductor), 1e5),
    ] {
        let one = [DielectricLayer {
            response,
            thickness_m: thickness,
        }];
        let split = [DielectricLayer {
            response,
            thickness_m: thickness * 0.5,
        }; 2];
        let whole = PlanarDielectricLayers::try_new(
            &ambient,
            DielectricResponse::Lossless(&ambient),
            OMEGA0,
            &one,
        )
        .unwrap();
        let pieces = PlanarDielectricLayers::try_new(
            &ambient,
            DielectricResponse::Lossless(&ambient),
            OMEGA0,
            &split,
        )
        .unwrap();
        for polarization in [Polarization::Te, Polarization::Tm] {
            let a = whole
                .second_order_frequency_dispersion(
                    0.8,
                    polarization,
                    IncidenceFrequencyConstraint::FixedAngle,
                )
                .unwrap();
            let b = pieces
                .second_order_frequency_dispersion(
                    0.8,
                    polarization,
                    IncidenceFrequencyConstraint::FixedAngle,
                )
                .unwrap();
            for (a, b) in [
                (
                    Some(a.reflected_tangential_e()),
                    Some(b.reflected_tangential_e()),
                ),
                (a.transmitted_tangential_e(), b.transmitted_tangential_e()),
                (a.transmitted_conjugate_h(), b.transmitted_conjugate_h()),
            ] {
                assert_eq!(a.is_some(), b.is_some());
                let (Some(a), Some(b)) = (a, b) else {
                    continue;
                };
                for (a, b) in a
                    .scaled_derivatives()
                    .into_iter()
                    .zip(b.scaled_derivatives())
                {
                    complex_near(a, b, 2e-10 * b.abs().max(1.0));
                }
            }
            if thickness >= 400.0 {
                let half = PlanarDielectricInterface::try_new(&ambient, response)
                    .unwrap()
                    .frequency_dispersion(
                        OMEGA0,
                        0.8,
                        polarization,
                        IncidenceFrequencyConstraint::FixedAngle,
                    )
                    .unwrap();
                assert_eq!(a.carrier().transmitted_tangential_e, C64::ZERO);
                assert!(a.transmitted_tangential_e().is_none());
                assert!(a.transmitted_conjugate_h().is_none());
                let first_order = pieces
                    .frequency_dispersion(
                        0.8,
                        polarization,
                        IncidenceFrequencyConstraint::FixedAngle,
                    )
                    .unwrap();
                assert!(first_order.transmitted_tangential_e().is_none());
                for (a, b) in a
                    .reflected_tangential_e()
                    .scaled_derivatives()
                    .into_iter()
                    .zip(half.reflected_tangential_e().scaled_derivatives())
                {
                    complex_near(a, b, 2e-13);
                }
            }
        }
    }
}

#[test]
fn g1_refracted_normal_delay_keeps_snell_geometry_and_frequency_constraint() {
    use IncidenceFrequencyConstraint::{FixedAngle, FixedTransverseWavenumber};
    let incident = source_with_storage(0.0, 2.0 * E, M);
    let transmitted = source_with_storage(0.0, 5.0 * E, 0.8 * M);
    let interface = PlanarDielectricInterface::try_new(
        &incident, DielectricResponse::Lossless(&transmitted),
    ).unwrap();
    let cosine = 0.6_f64;
    let incident_slowness = (2.0 * E * M).sqrt();
    let transmitted_slowness = (5.0 * E * 0.8 * M).sqrt();
    let sine_t = incident_slowness / transmitted_slowness * (1.0 - cosine * cosine).sqrt();
    let cosine_t = (1.0 - sine_t * sine_t).sqrt();
    let normal_slowness = transmitted_slowness * cosine_t;
    let distance = 0.23;
    for omega in [1e4, OMEGA0, 1e15] {
        for polarization in [Polarization::Te, Polarization::Tm] {
            let angle = interface.frequency_dispersion(omega, cosine, polarization, FixedAngle).unwrap();
            let fixed_k = interface.frequency_dispersion(
                omega, cosine, polarization, FixedTransverseWavenumber,
            ).unwrap();
            // Snell's conserved tangential wave vector supplies the carrier.
            near(angle.carrier().transmitted_normal_wavenumber_rad_m.re,
                omega * normal_slowness, 5e-15 * omega * normal_slowness);
            near(angle.transmitted_normal_wavenumber_derivative_s_m().re,
                normal_slowness, 5e-15 * normal_slowness);
            assert_eq!(angle.transmitted_normal_wavenumber_second_derivative_s2_m(), C64::ZERO);
            let expected_first = transmitted_slowness / cosine_t;
            let expected_second = -expected_first * sine_t * sine_t / (omega * cosine_t * cosine_t);
            near(fixed_k.transmitted_normal_wavenumber_derivative_s_m().re,
                expected_first, 6e-15 * expected_first);
            near(fixed_k.transmitted_normal_wavenumber_second_derivative_s2_m().re,
                expected_second, 1e-14 * expected_second.abs());
            let propagated = fixed_k.transmitted_propagation(distance).unwrap();
            near(propagated.group_delay_s, distance * expected_first, 6e-15 * distance * expected_first);
            near(propagated.group_delay_dispersion_s2, distance * expected_second,
                1e-14 * distance * expected_second.abs());
            assert_eq!(propagated.log_amplitude_derivative_s, 0.0);
            assert_eq!(propagated.log_amplitude_second_derivative_s2, 0.0);
        }
    }
}

#[test]
fn g1_refracted_jets_follow_original_dispersive_snell_branches() {
    use IncidenceFrequencyConstraint::{FixedAngle, FixedTransverseWavenumber};
    let ambient = source(0.03 * OMEGA0);
    let incident_resonance = lorentz(source(0.0));
    let relaxation = debye(ambient.clone(), 0.2, 0.01);
    let resonance = lorentz(ambient.clone());
    for incident in [DielectricResponse::Lossless(&ambient), DielectricResponse::Lorentz(&incident_resonance)] {
        for transmitted in [DielectricResponse::Debye(&relaxation), DielectricResponse::Lorentz(&resonance),
            DielectricResponse::Lorentz(&incident_resonance)] {
            let interface = PlanarDielectricInterface::try_new_with_incident_response(incident, transmitted).unwrap();
            for omega in [0.4 * OMEGA0, 0.8 * OMEGA0] {
                for cosine in [0.4_f64, 0.9, 1.0] {
                    let k_parallel = interface.response(omega, 1.0, Polarization::Te).unwrap()
                        .incident_normal_wavenumber_rad_m * (1.0 - cosine * cosine).sqrt();
                    for constraint in [FixedAngle, FixedTransverseWavenumber] {
                        let jet = interface.frequency_dispersion(omega, cosine, Polarization::Tm, constraint).unwrap();
                        let step = 1e-4 * omega;
                        let waves = core::array::from_fn(|index| {
                            let frequency = omega + (index as f64 - 2.0) * step;
                            let c = if constraint == FixedAngle { cosine } else {
                                let k = interface.response(frequency, 1.0, Polarization::Tm).unwrap()
                                    .incident_normal_wavenumber_rad_m;
                                (1.0 - (k_parallel / k).powi(2)).sqrt()
                            };
                            interface.response(frequency, c, Polarization::Tm).unwrap()
                                .transmitted_normal_wavenumber_rad_m
                        });
                        let [first, second] = differentiate_samples(waves, step);
                        complex_near(jet.transmitted_normal_wavenumber_derivative_s_m(), first,
                            3e-8 * waves[2].abs() / omega);
                        complex_near(jet.transmitted_normal_wavenumber_second_derivative_s2_m(), second,
                            8e-6 * waves[2].abs() / omega / omega);
                        if cosine == 1.0 {
                            // The exact normal limit keeps the original bulk
                            // owner's weak-loss and optical cancellations.
                            let bulk = transmitted.bulk_dispersion(omega).unwrap();
                            assert_eq!(jet.transmitted_normal_wavenumber_derivative_s_m(),
                                bulk.wavenumber_derivative_s_m());
                            assert_eq!(jet.transmitted_normal_wavenumber_second_derivative_s2_m(),
                                bulk.wavenumber_second_derivative_s2_m());
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn g1_total_internal_reflection_retains_decay_dispersion_without_ray_delay() {
    use IncidenceFrequencyConstraint::{FixedAngle, FixedTransverseWavenumber};
    let incident = source_with_storage(0.0, 1.0, 1.0);
    let transmitted = source_with_storage(0.0, 0.25, 1.0);
    let interface = PlanarDielectricInterface::try_new(
        &incident, DielectricResponse::Lossless(&transmitted),
    ).unwrap();
    let omega = OMEGA0;
    let kappa = 0.5_f64.sqrt();
    let distance = 2.0 / omega;
    for constraint in [FixedAngle, FixedTransverseWavenumber] {
        let jet = interface.frequency_dispersion(omega, 0.5, Polarization::Tm, constraint).unwrap();
        let (first, second) = if constraint == FixedAngle {
            (-kappa, 0.0)
        } else {
            (0.25 / kappa, 0.25 * 0.75 / (omega * kappa.powi(3)))
        };
        complex_near(jet.transmitted_normal_wavenumber_derivative_s_m(), C64::new(0.0, first), 3e-15);
        complex_near(jet.transmitted_normal_wavenumber_second_derivative_s2_m(),
            C64::new(0.0, second), 3e-15 / omega);
        let propagation = jet.transmitted_propagation(distance).unwrap();
        near(propagation.transfer.re, (-omega * kappa * distance).exp(), 2e-15);
        assert_eq!(propagation.transfer.im, 0.0);
        assert_eq!(propagation.phase_delay_s, 0.0);
        assert_eq!(propagation.group_delay_s, 0.0);
        assert_eq!(propagation.group_delay_dispersion_s2, 0.0);
        near(propagation.log_amplitude_derivative_s, first * distance, 5e-15 * distance);
        near(propagation.log_amplitude_second_derivative_s2, second * distance, 5e-15 * distance / omega);
        let identity = jet.transmitted_propagation(0.0).unwrap();
        assert_eq!(identity.transfer, C64::ONE);
        assert_eq!(identity.group_delay_s, 0.0);
        for invalid_distance in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(jet.transmitted_propagation(invalid_distance).is_err());
        }
    }
}

#[test]
fn g3_refracted_normal_transfer_matches_original_beam_and_composes() {
    let ambient = source(0.03 * OMEGA0);
    let resonance = lorentz(ambient.clone());
    let transmitted = DielectricResponse::Lorentz(&resonance);
    let interface = PlanarDielectricInterface::try_new(&ambient, transmitted).unwrap();
    let omega = 0.7 * OMEGA0;
    let cosine = 0.8_f64;
    let ki = interface.response(omega, 1.0, Polarization::Te).unwrap().incident_normal_wavenumber_rad_m;
    let carrier = [ki * (1.0 - cosine * cosine).sqrt(), 0.0];
    let distance = 0.17;
    let limits = AngularSpectrumLimits {
        maximum_grid_points: 4, maximum_workspace_bytes: 100_000, maximum_fft_work: 10_000,
    };
    let grid = AngularSpectrumGrid { width: 2, height: 2, spacing_m: [1.0, 1.0] };
    let step = omega * 1e-4;
    let fields: [C64; 5] = core::array::from_fn(|index| {
        let frequency = omega + (index as f64 - 2.0) * step;
        let stack = PlanarDielectricLayers::try_new(&ambient, transmitted, frequency, &[]).unwrap();
        let beam = PlanarAngularSpectrum::prepare_with_carrier(
            &stack, grid, carrier, &[[C64::ZERO, C64::ONE]; 4], limits, || false,
        ).unwrap();
        let surface = beam.field(AngularSpectrumSide::Transmitted, 0.0, || false).unwrap();
        let exit = beam.field(AngularSpectrumSide::Transmitted, distance, || false).unwrap();
        exit.electric_v_per_m()[0][1] / surface.electric_v_per_m()[0][1]
    });
    let jet = interface.frequency_dispersion(omega, cosine, Polarization::Te,
        IncidenceFrequencyConstraint::FixedTransverseWavenumber).unwrap();
    let propagation = jet.transmitted_propagation(distance).unwrap();
    complex_near(propagation.transfer, fields[2], 3e-15);
    let logarithms = fields.map(|value| {
        let ratio = value / fields[2];
        C64::new(ratio.abs().ln(), ratio.im.atan2(ratio.re))
    });
    let [first, second] = differentiate_samples(logarithms, step);
    near(propagation.group_delay_s, -first.im, 2e-8 / omega);
    near(propagation.group_delay_dispersion_s2, -second.im, 5e-6 / omega / omega);
    near(propagation.log_amplitude_derivative_s, first.re, 2e-8 / omega);
    near(propagation.log_amplitude_second_derivative_s2, second.re, 5e-6 / omega / omega);
    let part = jet.transmitted_propagation(distance * 0.4).unwrap();
    let rest = jet.transmitted_propagation(distance * 0.6).unwrap();
    complex_near(part.transfer * rest.transfer, propagation.transfer, 3e-15);
    near(part.group_delay_s + rest.group_delay_s, propagation.group_delay_s, 1e-24);
    near(part.group_delay_dispersion_s2 + rest.group_delay_dispersion_s2,
        propagation.group_delay_dispersion_s2, 1e-32);
}
