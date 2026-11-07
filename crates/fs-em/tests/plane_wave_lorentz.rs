#![cfg(feature = "plane-wave")]
//! G1/G3 resonant dielectric scattering from synthetic source cards.
//! These regressions cover constitutive coupling and coherent wave accounting;
//! they do not establish measured-material accuracy or transient Maxwell behavior.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, UncertaintyModel,
};
use fs_qty::{
    QuantitySpec,
    semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm},
};

type Entry = (&'static str, QuantitySpec, f64);

fn angular_frequency() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Frequency(FrequencyConvention::Angular),
        ValueForm::Static,
    ))
}

fn card(entries: &[Entry], specimen: &str) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for &(name, quantity, value) in entries {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Scalar {
                    value,
                    dims: quantity.dims(),
                },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic Lorentz wave fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: specimen.into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

#[cfg(test)]
mod integration {
    use super::*;
    use fs_em::plane_wave::{
        DielectricResponse, PlanarDielectricInterface, Polarization,
        layers::{DielectricLayer, PlanarDielectricLayers},
    };
    use fs_matdb::QueryPoint;
    use fs_material::{
        conductor::{ElectricalPropertyBasis, ResolvedConductorState},
        field_storage::{
            ELECTRIC_PERMITTIVITY_DIMS as EPS, ELECTRIC_PERMITTIVITY_PROPERTY as EPS_NAME,
            MAGNETIC_PERMEABILITY_DIMS as MU, MAGNETIC_PERMEABILITY_PROPERTY as MU_NAME,
            dielectric::lorentz::{LorentzPoleProperties, ResolvedLorentzMaterial},
        },
        state_point::{
            MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
            ScalarPropertyRequirement, resolve_material_state_point,
        },
    };
    use fs_math::c64::C64;
    use fs_qty::Dims;

    const E: f64 = 1e-11;
    const M: f64 = 1e-6;
    const NATURAL: f64 = 1e9;
    const SIGMA: f64 = 0.02;

    fn source(entries: &[Entry], specimen: &str) -> ResolvedMaterialStatePoint {
        let requirements = entries
            .iter()
            .map(|&(name, quantity, _)| {
                ScalarPropertyRequirement::try_with_quantity(
                    name,
                    quantity,
                    ScalarAdmissibility::NonNegative,
                )
                .unwrap()
            })
            .collect::<Vec<_>>();
        resolve_material_state_point(
            &card(entries, specimen),
            &QueryPoint::new(),
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap()
    }

    fn outside() -> ResolvedMaterialStatePoint {
        source(
            &[
                (EPS_NAME, QuantitySpec::dimensional(EPS), E),
                (MU_NAME, QuantitySpec::dimensional(MU), M),
            ],
            "synthetic surrounding dielectric",
        )
    }

    fn resonant_source(damping: f64, specimen: &str) -> ResolvedMaterialStatePoint {
        let basis = ElectricalPropertyBasis::Conductivity;
        source(
            &[
                ("epsilon_inf", QuantitySpec::dimensional(EPS), 2.0 * E),
                ("delta", QuantitySpec::dimensional(EPS), 4.0 * E),
                ("omega0", angular_frequency(), NATURAL),
                (
                    "gamma",
                    QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])),
                    damping,
                ),
                (MU_NAME, QuantitySpec::dimensional(MU), M),
                (
                    basis.property_name(),
                    QuantitySpec::dimensional(basis.dimensions()),
                    SIGMA,
                ),
            ],
            specimen,
        )
    }

    fn lorentz(source: ResolvedMaterialStatePoint) -> ResolvedLorentzMaterial {
        ResolvedLorentzMaterial::try_from_material(
            source,
            "epsilon_inf",
            &[LorentzPoleProperties {
                permittivity_increment: "delta".into(),
                angular_frequency: "omega0".into(),
                damping_rate: "gamma".into(),
            }],
        )
        .unwrap()
    }

    fn close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "{actual:e} != {expected:e}"
        );
    }

    fn complex_close(actual: C64, expected: C64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "{actual:?} != {expected:?}"
        );
    }

    #[test]
    fn resonance_and_independent_ohmic_loss_preserve_coherent_passive_scattering() {
        let ambient = outside();
        let material = resonant_source(0.03 * NATURAL, "resonant specimen A");
        let polarization = lorentz(material.clone());
        let conductor = ResolvedConductorState::try_from_material_with_basis(
            material,
            ElectricalPropertyBasis::Conductivity,
        )
        .unwrap();
        assert_eq!(
            polarization.material().identity(),
            conductor.material().identity()
        );
        for ratio in [0.5, 0.99, 1.0, 1.15, 2.0] {
            let omega = ratio * NATURAL;
            let bare_epsilon = polarization.permittivity_f_m(omega).unwrap();
            if ratio == 1.15 {
                assert!(bare_epsilon.re < 0.0);
            }
            for (response, expected_epsilon) in [
                (DielectricResponse::Lorentz(&polarization), bare_epsilon),
                (
                    DielectricResponse::LorentzWithOhmic {
                        dielectric: &polarization,
                        conductor: &conductor,
                    },
                    C64::new(bare_epsilon.re, bare_epsilon.im - SIGMA / omega),
                ),
            ] {
                let interface = PlanarDielectricInterface::try_new(&ambient, response).unwrap();
                let layer = DielectricLayer {
                    response,
                    thickness_m: 0.12,
                };
                let film = PlanarDielectricLayers::try_new(
                    &ambient,
                    DielectricResponse::Lossless(&ambient),
                    omega,
                    &[layer],
                )
                .unwrap();
                let split = PlanarDielectricLayers::try_new(
                    &ambient,
                    DielectricResponse::Lossless(&ambient),
                    omega,
                    &[DielectricLayer {
                        response,
                        thickness_m: 0.04,
                    }; 3],
                )
                .unwrap();
                for mode in [Polarization::Te, Polarization::Tm] {
                    let half_space = interface.response(omega, 0.7, mode).unwrap();
                    complex_close(
                        half_space.transmitted_permittivity_f_m,
                        expected_epsilon,
                        3e-14 * expected_epsilon.abs(),
                    );
                    assert!(half_space.transmitted_normal_wavenumber_rad_m.im <= 0.0);
                    assert!(half_space.reflected_power_fraction >= 0.0);
                    assert!(half_space.transmitted_power_fraction >= 0.0);
                    // This is flux entering the lossy half-space, not a claim
                    // that this power reaches an observer at infinite depth.
                    close(
                        half_space.reflected_power_fraction + half_space.transmitted_power_fraction,
                        1.0,
                        3e-13,
                    );
                    let full = film.response(0.7, mode).unwrap();
                    let subdivided = split.response(0.7, mode).unwrap();
                    assert!(full.absorbed_power_fraction > 0.0);
                    assert!(full.reflected_power_fraction >= 0.0);
                    assert!(full.transmitted_power_fraction >= 0.0);
                    close(
                        full.reflected_power_fraction
                            + full.transmitted_power_fraction
                            + full.absorbed_power_fraction,
                        1.0,
                        3e-13,
                    );
                    complex_close(
                        full.reflected_tangential_e,
                        subdivided.reflected_tangential_e,
                        2e-12,
                    );
                    complex_close(
                        full.transmitted_tangential_e,
                        subdivided.transmitted_tangential_e,
                        2e-12,
                    );
                    close(
                        full.absorbed_power_fraction,
                        subdivided.absorbed_power_fraction,
                        2e-12,
                    );
                }
            }
        }
    }

    #[test]
    fn equal_coefficients_from_different_complete_material_sources_cannot_be_combined() {
        let ambient = outside();
        let polarization = lorentz(resonant_source(0.03 * NATURAL, "specimen A"));
        let other = lorentz(resonant_source(0.03 * NATURAL, "specimen B"));
        assert_eq!(
            polarization.permittivity_f_m(NATURAL).unwrap(),
            other.permittivity_f_m(NATURAL).unwrap()
        );
        assert_ne!(
            polarization.material().identity(),
            other.material().identity()
        );
        let conductor = ResolvedConductorState::try_from_material_with_basis(
            other.material().clone(),
            ElectricalPropertyBasis::Conductivity,
        )
        .unwrap();
        let mismatched = DielectricResponse::LorentzWithOhmic {
            dielectric: &polarization,
            conductor: &conductor,
        };
        assert!(PlanarDielectricInterface::try_new(&ambient, mismatched).is_err());
        assert!(
            PlanarDielectricLayers::try_new(
                &ambient,
                DielectricResponse::Lossless(&ambient),
                NATURAL,
                &[DielectricLayer {
                    response: mismatched,
                    thickness_m: 0.12
                }],
            )
            .is_err()
        );
    }

    #[test]
    fn lossless_negative_permittivity_reflects_and_tunnels_without_invented_heat() {
        let ambient = outside();
        let polarization = lorentz(resonant_source(0.0, "lossless resonant specimen"));
        let response = DielectricResponse::Lorentz(&polarization);
        let interface = PlanarDielectricInterface::try_new(&ambient, response).unwrap();
        let omega = 1.15 * NATURAL;
        assert!(polarization.permittivity_f_m(omega).unwrap().re < 0.0);
        let film = PlanarDielectricLayers::try_new(
            &ambient,
            DielectricResponse::Lossless(&ambient),
            omega,
            &[DielectricLayer {
                response,
                thickness_m: 0.12,
            }],
        )
        .unwrap();
        for mode in [Polarization::Te, Polarization::Tm] {
            let half_space = interface.response(omega, 1.0, mode).unwrap();
            close(half_space.reflected_power_fraction, 1.0, 2e-13);
            close(half_space.transmitted_power_fraction, 0.0, 2e-13);
            assert!(half_space.transmitted_normal_wavenumber_rad_m.im < 0.0);
            let finite = film.response(1.0, mode).unwrap();
            assert!(
                finite.transmitted_power_fraction > 0.0 && finite.transmitted_power_fraction < 1.0
            );
            close(finite.absorbed_power_fraction, 0.0, 3e-13);
            close(
                finite.reflected_power_fraction + finite.transmitted_power_fraction,
                1.0,
                3e-13,
            );
            assert!(interface.response(NATURAL, 1.0, mode).is_err());
        }
        assert!(
            PlanarDielectricLayers::try_new(
                &ambient,
                DielectricResponse::Lossless(&ambient),
                NATURAL,
                &[DielectricLayer {
                    response,
                    thickness_m: 0.12
                }],
            )
            .is_err()
        );
    }

    mod finite_window {
        use super::*;
        use fs_em::plane_wave::angular_spectrum::{
            AngularSpectrumExposure, AngularSpectrumExposureHeatChannel as Channel,
            AngularSpectrumExposureHeatSpectrumLimits, AngularSpectrumExposureLayerLimits,
            AngularSpectrumExposureLimits, AngularSpectrumGrid,
            AngularSpectrumHeatAreaQuadrature as Areas,
            AngularSpectrumHeatDepthSample as DepthSample,
            AngularSpectrumHeatVolumeQuadrature as Volume, AngularSpectrumIntegratedHeatLimits,
            AngularSpectrumLayerField, AngularSpectrumLayerLimits, AngularSpectrumLimits,
            PlanarAngularSpectrum,
        };
        use fs_fft::harmonic_window::{HarmonicWindow, HarmonicWindowProfile as Profile};

        fn beam<'s, 'm>(
            stack: &'s PlanarDielectricLayers<'m>,
            electric: [C64; 2],
        ) -> PlanarAngularSpectrum<'s, 'm> {
            PlanarAngularSpectrum::prepare(
                stack,
                AngularSpectrumGrid {
                    width: 4,
                    height: 1,
                    spacing_m: [2.0, 3.0],
                },
                &[electric; 4],
                AngularSpectrumLimits {
                    maximum_grid_points: 128,
                    maximum_workspace_bytes: 1 << 20,
                    maximum_fft_work: 1 << 20,
                },
                || false,
            )
            .unwrap()
        }

        fn exposure_limits() -> AngularSpectrumExposureLimits {
            AngularSpectrumExposureLimits {
                maximum_lines: 8,
                maximum_point_line_pairs: 1024,
                maximum_point_frequency_pairs: 1024,
                maximum_workspace_bytes: 2 << 20,
                maximum_fft_work: 1 << 20,
            }
        }

        fn heat_limits() -> AngularSpectrumExposureLayerLimits {
            AngularSpectrumExposureLayerLimits {
                maximum_frequency_groups: 16,
                maximum_frequency_layer_mode_pairs: 4096,
                maximum_point_line_pairs: 4096,
                maximum_point_reservoir_frequency_pairs: 16384,
                maximum_workspace_bytes: 16 << 20,
                maximum_fft_work: 1 << 20,
            }
        }

        fn spectrum_limits() -> AngularSpectrumExposureHeatSpectrumLimits {
            AngularSpectrumExposureHeatSpectrumLimits {
                maximum_products: 4096,
                maximum_generated_terms: 8192,
                maximum_workspace_bytes: 32 << 20,
            }
        }

        // Independent real-time oracle for the original oscillator
        // P'' + gamma P' + omega0^2 P = delta_epsilon omega0^2 E.
        // The phasors come from the actual film-field adapter. Neither the
        // production reservoir transfer nor its harmonic product owner is used.
        // Returns Ohmic heat, polarization heat, full storage rate, energy,
        // and the kinetic contribution to that storage rate, in SI units.
        fn instantaneous(
            fields: &[AngularSpectrumLayerField],
            frequencies: &[f64],
            gamma: f64,
            sigma: f64,
            time: f64,
        ) -> [f64; 5] {
            let mut values = [[0.0; 3]; 7];
            let delta = 4.0 * E;
            for (field, &omega) in fields.iter().zip(frequencies) {
                let (sine, cosine) = (omega * time).sin_cos();
                let real = |z: C64| z.re * cosine - z.im * sine;
                let iw = C64::new(0.0, omega);
                let susceptibility = C64::from_re(delta * NATURAL * NATURAL)
                    / C64::new(NATURAL * NATURAL - omega * omega, gamma * omega);
                for axis in 0..3 {
                    let electric = field.field().electric_v_per_m()[0][axis];
                    let magnetic = field.field().magnetic_a_per_m()[0][axis];
                    let p = electric * susceptibility;
                    let phasors = [
                        electric,
                        iw * electric,
                        magnetic,
                        iw * magnetic,
                        p,
                        iw * p,
                        p.scale(-omega * omega),
                    ];
                    for (value, phasor) in values.iter_mut().zip(phasors) {
                        value[axis] += real(phasor);
                    }
                }
            }
            let [
                electric,
                electric_rate,
                magnetic,
                magnetic_rate,
                p,
                velocity,
                acceleration,
            ] = values;
            let dot =
                |a: [f64; 3], b: [f64; 3]| a.into_iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
            let kinetic_rate = dot(velocity, acceleration) / (delta * NATURAL * NATURAL);
            let storage_rate = 2.0 * E * dot(electric, electric_rate)
                + M * dot(magnetic, magnetic_rate)
                + dot(p, velocity) / delta
                + kinetic_rate;
            let energy = 0.5
                * (2.0 * E * dot(electric, electric)
                    + M * dot(magnetic, magnetic)
                    + dot(p, p) / delta
                    + dot(velocity, velocity) / (delta * NATURAL * NATURAL));
            [
                sigma * dot(electric, electric),
                gamma * dot(velocity, velocity) / (delta * NATURAL * NATURAL),
                storage_rate,
                energy,
                kinetic_rate,
            ]
        }

        fn average(
            fields: &[AngularSpectrumLayerField],
            frequencies: &[f64],
            gamma: f64,
            sigma: f64,
            window: HarmonicWindow,
        ) -> [f64; 5] {
            // Composite Simpson integration in normalized time. The triangular
            // detector's kink is exactly a panel boundary, and its weight is
            // applied once to physical products, never to each field amplitude.
            let count = 1024;
            let mut sum = [0.0; 5];
            for i in 0..=count {
                let u = i as f64 / count as f64;
                let sensitivity = match window.profile() {
                    Profile::Uniform => 1.0,
                    Profile::Triangular => 2.0 * (1.0 - (2.0 * u - 1.0).abs()),
                };
                let weight = if i == 0 || i == count {
                    1.0
                } else if i % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                let value = instantaneous(
                    fields,
                    frequencies,
                    gamma,
                    sigma,
                    window.start_s() + u * window.duration_s(),
                );
                for (sum, value) in sum.iter_mut().zip(value) {
                    *sum += weight * sensitivity * value;
                }
            }
            sum.map(|value| value / (3.0 * count as f64))
        }

        fn check(actual: f64, expected: f64) {
            close(actual, expected, 3e-10 * expected.abs().max(1e-3));
        }

        #[test]
        fn coherent_short_window_lorentz_heat_storage_and_volume_use_original_film_fields() {
            let ambient = outside();
            let gamma = 0.08 * NATURAL;
            let material = resonant_source(gamma, "finite-window resonant film");
            let dielectric = lorentz(material.clone());
            let conductor = ResolvedConductorState::try_from_material_with_basis(
                material,
                ElectricalPropertyBasis::Conductivity,
            )
            .unwrap();
            let frequencies = [0.92 * NATURAL, 1.13 * NATURAL];
            let thickness = 0.12;
            let stacks = frequencies.map(|omega| {
                PlanarDielectricLayers::try_new(
                    &ambient,
                    DielectricResponse::Lossless(&ambient),
                    omega,
                    &[DielectricLayer {
                        response: DielectricResponse::LorentzWithOhmic {
                            dielectric: &dielectric,
                            conductor: &conductor,
                        },
                        thickness_m: thickness,
                    }],
                )
                .unwrap()
            });
            let lines = [
                beam(&stacks[0], [C64::new(0.8, 0.2), C64::new(-0.1, 0.3)]),
                beam(&stacks[1], [C64::new(-0.3, 0.6), C64::new(0.4, -0.1)]),
            ];
            let bank =
                AngularSpectrumExposure::try_new(&lines[0], &lines, exposure_limits(), || false)
                    .unwrap();
            let layer = bank.prepare_layer(0, heat_limits(), || false).unwrap();
            let native_layers = lines
                .iter()
                .map(|line| {
                    line.prepare_layer(
                        0,
                        AngularSpectrumLayerLimits {
                            maximum_layer_mode_pairs: 4096,
                            maximum_workspace_bytes: 16 << 20,
                            maximum_fft_work: 1 << 20,
                        },
                        || false,
                    )
                    .unwrap()
                })
                .collect::<Vec<_>>();
            let depth = 0.5 * thickness;
            let fields = native_layers
                .iter()
                .map(|line| line.field(depth, || false).unwrap())
                .collect::<Vec<_>>();
            let channels = [
                Channel::Ohmic,
                Channel::LorentzPole { index: 0 },
                Channel::Polarization,
                Channel::Total,
            ];
            let spectra = channels.map(|channel| {
                layer
                    .heat_spectrum(depth, 0, channel, spectrum_limits(), || false)
                    .unwrap()
            });
            for time in [-0.3 / NATURAL, 0.17 / NATURAL, 1.21 / NATURAL] {
                let reference = instantaneous(&fields, &frequencies, gamma, SIGMA, time);
                for (spectrum, expected) in spectra.iter().zip([
                    reference[0],
                    reference[1],
                    reference[1],
                    reference[0] + reference[1],
                ]) {
                    check(spectrum.evaluate(time).unwrap(), expected);
                }
            }
            assert!(
                layer
                    .heat_spectrum(
                        depth,
                        0,
                        Channel::DebyePole { index: 0 },
                        spectrum_limits(),
                        || false
                    )
                    .is_err()
            );
            assert!(
                layer
                    .heat_spectrum(
                        depth,
                        0,
                        Channel::LorentzPole { index: 1 },
                        spectrum_limits(),
                        || false
                    )
                    .is_err()
            );

            let depth_samples = (0..=8)
                .map(|i| DepthSample {
                    depth_m: thickness * i as f64 / 8.0,
                    weight_m: thickness / 24.0
                        * if i == 0 || i == 8 {
                            1.0
                        } else if i % 2 == 0 {
                            2.0
                        } else {
                            4.0
                        },
                })
                .collect::<Vec<_>>();
            let volume = Volume {
                depth_span_m: [0.0, thickness],
                depth_samples: &depth_samples,
                transverse: Areas::ContinuousPeriodicCell,
            };
            let load = layer
                .integrated_heat_spectrum(
                    volume,
                    Channel::Total,
                    AngularSpectrumIntegratedHeatLimits {
                        maximum_depth_samples: 128,
                        maximum_point_line_pairs: 1 << 20,
                        maximum_fft_work: 1 << 24,
                        maximum_spatial_point_frequency_pairs: 1 << 20,
                        maximum_reservoir_frequency_pairs: 1 << 16,
                        maximum_generated_terms: 8192,
                        maximum_workspace_bytes: 64 << 20,
                    },
                    || false,
                )
                .unwrap();
            close(
                load.represented_volume_m3_dd().to_f64(),
                24.0 * thickness,
                2e-14,
            );
            let native_volume = depth_samples
                .iter()
                .map(|sample| {
                    native_layers
                        .iter()
                        .map(|line| line.field(sample.depth_m, || false).unwrap())
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            for profile in [Profile::Uniform, Profile::Triangular] {
                let window =
                    HarmonicWindow::with_profile(0.17 / NATURAL, 0.43 / NATURAL, profile).unwrap();
                let expected = average(&fields, &frequencies, gamma, SIGMA, window);
                let observed = layer.field(depth, window, || false).unwrap();
                for q in observed.dissipation() {
                    check(q.ohmic_heat_w_m3, expected[0]);
                    check(q.polarization_heat_w_m3, expected[1]);
                    check(q.total_heat_w_m3, expected[0] + expected[1]);
                    check(q.storage_rate_w_m3, expected[2]);
                }
                assert!(
                    expected[4].abs() > 1e-6,
                    "the test must resolve inertial storage"
                );
                for (spectrum, expected) in spectra.iter().zip([
                    expected[0],
                    expected[1],
                    expected[1],
                    expected[0] + expected[1],
                ]) {
                    check(spectrum.observe_window(window).unwrap(), expected);
                }
                let separate = fields
                    .iter()
                    .zip(frequencies)
                    .map(|(field, omega)| {
                        let q =
                            average(std::slice::from_ref(field), &[omega], gamma, SIGMA, window);
                        q[0] + q[1]
                    })
                    .sum::<f64>();
                assert!(
                    (expected[0] + expected[1] - separate).abs() > 1e-8,
                    "coherent cross-frequency heat must not become independent-line power"
                );
                if profile == Profile::Uniform {
                    let start =
                        instantaneous(&fields, &frequencies, gamma, SIGMA, window.start_s())[3];
                    let end = instantaneous(
                        &fields,
                        &frequencies,
                        gamma,
                        SIGMA,
                        window.start_s() + window.duration_s(),
                    )[3];
                    check(
                        observed.dissipation()[0].storage_rate_w_m3,
                        (end - start) / window.duration_s(),
                    );
                }
                let integrated = depth_samples
                    .iter()
                    .zip(&native_volume)
                    .map(|(sample, fields)| {
                        let q = average(fields, &frequencies, gamma, SIGMA, window);
                        24.0 * sample.weight_m * (q[0] + q[1])
                    })
                    .sum::<f64>();
                check(load.observe_window(window).unwrap(), integrated);
            }
        }

        #[test]
        fn lossless_lorentz_window_keeps_inertial_storage_without_heat_and_charges_its_work() {
            let ambient = outside();
            let dielectric = lorentz(resonant_source(0.0, "lossless finite-window film"));
            let omega = 1.15 * NATURAL;
            let stack = PlanarDielectricLayers::try_new(
                &ambient,
                DielectricResponse::Lossless(&ambient),
                omega,
                &[DielectricLayer {
                    response: DielectricResponse::Lorentz(&dielectric),
                    thickness_m: 0.12,
                }],
            )
            .unwrap();
            let lines = [beam(&stack, [C64::new(0.8, 0.2), C64::new(-0.1, 0.3)])];
            let bank =
                AngularSpectrumExposure::try_new(&lines[0], &lines, exposure_limits(), || false)
                    .unwrap();
            let layer = bank.prepare_layer(0, heat_limits(), || false).unwrap();
            // Four grid points, one frequency pair and 3+2*1 contractions:
            // the kinetic Lorentz store must not be charged as a Debye pole.
            assert_eq!(layer.point_reservoir_frequency_pairs(), 20);
            assert!(
                bank.prepare_layer(
                    0,
                    AngularSpectrumExposureLayerLimits {
                        maximum_point_reservoir_frequency_pairs: 19,
                        ..heat_limits()
                    },
                    || false
                )
                .is_err()
            );
            let original = lines[0]
                .prepare_layer(
                    0,
                    AngularSpectrumLayerLimits {
                        maximum_layer_mode_pairs: 4096,
                        maximum_workspace_bytes: 16 << 20,
                        maximum_fft_work: 1 << 20,
                    },
                    || false,
                )
                .unwrap()
                .field(0.06, || false)
                .unwrap();
            let window = HarmonicWindow::new(0.17 / NATURAL, 0.43 / NATURAL).unwrap();
            let q = layer.field(0.06, window, || false).unwrap().dissipation()[0];
            assert_eq!(q.ohmic_heat_w_m3, 0.0);
            assert_eq!(q.polarization_heat_w_m3, 0.0);
            assert_eq!(q.total_heat_w_m3, 0.0);
            let expected = average(&[original], &[omega], 0.0, 0.0, window);
            check(q.storage_rate_w_m3, expected[2]);
            assert!(expected[4].abs() > 1e-6);
            for channel in [
                Channel::LorentzPole { index: 0 },
                Channel::Polarization,
                Channel::Total,
            ] {
                let spectrum = layer
                    .heat_spectrum(0.06, 0, channel, spectrum_limits(), || false)
                    .unwrap();
                assert_eq!(spectrum.evaluate(window.start_s()).unwrap(), 0.0);
                assert_eq!(spectrum.observe_window(window).unwrap(), 0.0);
            }
        }
    }
}
