//! Synthetic dielectric sources exercise the original Maxwell and source owners.
use super::*;
use fs_em::plane_wave::{
    OpticalMaterialResponse,
    angular_spectrum::{
        AngularSpectrumExposure, AngularSpectrumExposureLimits, AngularSpectrumGrid,
        AngularSpectrumLimits, AngularSpectrumSide, PlanarAngularSpectrum,
    },
    layers::{DielectricLayer, PlanarDielectricLayers},
};
use fs_fft::harmonic_window::HarmonicWindow;
use fs_material::state_point::{
    MaterialStatePointError, VISIBLE_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES,
    VISIBLE_DIELECTRIC_ABSORPTION_QUANTITY, VISIBLE_DIELECTRIC_CAUCHY_A_PROPERTY,
    VISIBLE_DIELECTRIC_CAUCHY_B_M2_PROPERTY, VISIBLE_DIELECTRIC_CAUCHY_C_M4_PROPERTY,
    VISIBLE_DIELECTRIC_IOR_PROPERTIES, VISIBLE_DIELECTRIC_REFERENCE_DISTANCE_M_PROPERTY,
    VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES, VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES,
    VISIBLE_DIELECTRIC_TRANSMITTANCE_PROPERTIES, VISIBLE_DIELECTRIC_WAVELENGTHS_NM,
    VisibleOpticalStatePoint, resolve_visible_optical_state_point,
};
use fs_math::det;
use fs_qty::QuantitySpec;

const C0: f64 = 299_792_458.0;
const WL_FREQUENCY: f64 = core::f64::consts::TAU * C0 * 1e9;

#[derive(Clone, Copy, Debug)]
enum Family {
    Cauchy,
    Tabulated,
    Sellmeier,
}

fn glass_card(family: Family, spectral: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let mut add = |name, quantity: QuantitySpec, value| {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Scalar {
                    value,
                    dims: quantity.dims(),
                },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored dielectric Maxwell oracle".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    match family {
        Family::Cauchy => {
            for (name, dims, value) in [
                (VISIBLE_DIELECTRIC_CAUCHY_A_PROPERTY, Dims::NONE, 1.4),
                (
                    VISIBLE_DIELECTRIC_CAUCHY_B_M2_PROPERTY,
                    Dims([2, 0, 0, 0, 0, 0]),
                    5e-15,
                ),
                (
                    VISIBLE_DIELECTRIC_CAUCHY_C_M4_PROPERTY,
                    Dims([4, 0, 0, 0, 0, 0]),
                    2e-28,
                ),
            ] {
                add(name, QuantitySpec::dimensional(dims), value);
            }
        }
        Family::Tabulated => {
            for (name, wavelength) in VISIBLE_DIELECTRIC_IOR_PROPERTIES
                .into_iter()
                .zip(VISIBLE_DIELECTRIC_WAVELENGTHS_NM)
            {
                add(
                    name,
                    QuantitySpec::dimensional(Dims::NONE),
                    phase_index(family, wavelength),
                );
            }
        }
        Family::Sellmeier => {
            for (name, value) in VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES
                .into_iter()
                .zip([0.8, 0.2, 0.5])
            {
                add(name, QuantitySpec::dimensional(Dims::NONE), value);
            }
            for (name, value) in VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES
                .into_iter()
                .zip([0.01e-12, 0.04e-12, 81e-12])
            {
                add(
                    name,
                    QuantitySpec::dimensional(Dims([2, 0, 0, 0, 0, 0])),
                    value,
                );
            }
        }
    }
    if spectral {
        for (name, wavelength) in VISIBLE_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES
            .into_iter()
            .zip(VISIBLE_DIELECTRIC_WAVELENGTHS_NM)
        {
            add(
                name,
                VISIBLE_DIELECTRIC_ABSORPTION_QUANTITY,
                alpha(wavelength),
            );
        }
    } else {
        for name in VISIBLE_DIELECTRIC_TRANSMITTANCE_PROPERTIES {
            add(name, QuantitySpec::dimensional(Dims::NONE), 0.9);
        }
        add(
            VISIBLE_DIELECTRIC_REFERENCE_DISTANCE_M_PROPERTY,
            QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
            0.01,
        );
    }
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(MU_NAME, MU),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: fs_qty::Temperature::DIMS,
                knots: vec![(280.0, 1.2 * M), (320.0, 1.6 * M)],
                dims: MU,
            },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "authored same-state magnetic oracle".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: format!("authored {family:?}"),
            phase: "solid".into(),
            process: "numerical spectrum".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn phase_index(family: Family, wavelength_nm: f64) -> f64 {
    let u2 = (wavelength_nm / 1000.0).powi(2);
    match family {
        Family::Cauchy => 1.4 + 0.005 / u2 + 0.0002 / (u2 * u2),
        Family::Tabulated => 1.45 + 0.0002 * (wavelength_nm - 380.0),
        Family::Sellmeier => det::sqrt(
            1.0 + 0.8 * u2 / (u2 - 0.01) + 0.2 * u2 / (u2 - 0.04) + 0.5 * u2 / (u2 - 81.0),
        ),
    }
}
fn alpha(wavelength_nm: f64) -> f64 {
    2e5 + 1e3 * (wavelength_nm - 380.0)
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with("T", t).unwrap()
}
fn extend(source: &VisibleOpticalStatePoint, card: &MaterialCard) -> VisibleOpticalStatePoint {
    source
        .with_additional_requirements(
            card,
            &[ScalarPropertyRequirement::try_new(
                MU_NAME,
                MU,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap()],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap()
}

#[test]
fn original_dielectric_families_retain_dispersion_phase_and_coherent_film_heat() {
    let incident = material(1.0 / (M * C0 * C0), M);
    for family in [Family::Cauchy, Family::Tabulated, Family::Sellmeier] {
        let card = glass_card(family, true);
        let original = resolve_visible_optical_state_point(
            &card,
            &point(300.0),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert!(
            OpticalMaterialResponse::try_new(&original).is_err(),
            "no default mu"
        );
        let source = extend(&original, &card);
        let original_name = original.resolved().properties()[0].requirement().name();
        assert!(matches!(
            original.with_additional_requirements(
                &card,
                &[ScalarPropertyRequirement::try_new(
                    original_name,
                    MU,
                    ScalarAdmissibility::StrictlyPositive
                )
                .unwrap(),],
                MaterialPropertySelection::SingleClaimOnly
            ),
            Err(MaterialStatePointError::DuplicateRequirement { .. })
        ));
        for property in original.resolved().properties() {
            assert_eq!(
                source.resolved().property(property.requirement().name()),
                Some(property)
            );
        }
        let optical = OpticalMaterialResponse::try_new(&source).unwrap();
        assert!(std::ptr::eq(optical.material(), source.resolved()));
        near(optical.permeability_h_m(), 1.4 * M, 3e-22);
        for wavelength in [380.0, 550.0, 780.0] {
            let omega = WL_FREQUENCY / wavelength;
            let index = C64::new(
                phase_index(family, wavelength),
                -alpha(wavelength) * wavelength * 1e-9 / (4.0 * core::f64::consts::PI),
            );
            complex_near(optical.complex_index(omega).unwrap(), index, 2e-15);
            let epsilon = (index * index).scale(1.0 / (1.4 * M * C0 * C0));
            complex_near(
                optical.permittivity_f_m(omega).unwrap(),
                epsilon,
                5e-15 * epsilon.abs(),
            );
            for cosine in [1.0, 0.6] {
                let k0 = omega / C0;
                let kz = ((index * index).scale(k0 * k0)
                    - C64::from_re(k0 * k0 * (1.0 - cosine * cosine)))
                .sqrt();
                for pol in [Polarization::Te, Polarization::Tm] {
                    let yi = if pol == Polarization::Te {
                        cosine / (M * C0)
                    } else {
                        1.0 / (M * C0 * cosine)
                    };
                    let yt = if pol == Polarization::Te {
                        kz.scale(1.0 / (omega * 1.4 * M))
                    } else {
                        epsilon.scale(omega) / kz
                    };
                    let r01 = (C64::from_re(yi) - yt) / (C64::from_re(yi) + yt);
                    let response = DielectricResponse::OpticalMaterial(optical);
                    let interface =
                        PlanarDielectricInterface::try_new(&incident, response).unwrap();
                    complex_near(
                        interface
                            .response(omega, cosine, pol)
                            .unwrap()
                            .reflected_tangential_e,
                        r01,
                        3e-14,
                    );
                    let thickness = 420e-9;
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
                    let actual = stack.response(cosine, pol).unwrap();
                    let (s, c) = det::sin_cos(kz.re * thickness);
                    let p = C64::new(c, -s).scale(det::exp(kz.im * thickness));
                    let denominator = C64::ONE - r01 * r01 * p * p;
                    let r = r01 * (C64::ONE - p * p) / denominator;
                    let t = (C64::ONE - r01 * r01) * p / denominator;
                    complex_near(actual.reflected_tangential_e, r, 3e-13);
                    complex_near(actual.transmitted_tangential_e, t, 3e-13);
                    near(
                        actual.absorbed_power_fraction,
                        1.0 - r.norm_sq() - t.norm_sq(),
                        5e-13,
                    );
                    assert!(actual.absorbed_power_fraction > 0.0);
                }
            }
        }
        let next = source.at_query_point(&card, &point(310.0)).unwrap();
        near(
            OpticalMaterialResponse::try_new(&next)
                .unwrap()
                .permeability_h_m(),
            1.5 * M,
            3e-22,
        );
        assert_eq!(source.resolved().state_coordinate("T"), Some(300.0));
        assert!(source.at_query_point(&card, &point(350.0)).is_err());
        let optical = OpticalMaterialResponse::try_new(&source).unwrap();
        assert!(
            optical
                .complex_index((WL_FREQUENCY / 780.0).next_down())
                .is_err()
        );
        assert!(
            optical
                .complex_index((WL_FREQUENCY / 380.0).next_up())
                .is_err()
        );
    }
}

#[test]
fn rgb_absorption_and_hidden_poles_cannot_be_promoted_to_wave_sources() {
    for family in [Family::Cauchy, Family::Sellmeier] {
        let card = glass_card(family, false);
        let source = resolve_visible_optical_state_point(
            &card,
            &point(300.0),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert!(OpticalMaterialResponse::try_new(&extend(&source, &card)).is_err());
    }
    let base = glass_card(Family::Sellmeier, true);
    let mut claims = ClaimSet::new();
    for (_, claim) in base.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES[0] {
            claim.value = PropertyValue::Scalar {
                value: (543.21e-9_f64).powi(2),
                dims: Dims([2, 0, 0, 0, 0, 0]),
            };
        }
        claims.insert_claim(claim).unwrap();
    }
    let card = MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap();
    let source = resolve_visible_optical_state_point(
        &card,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let complete = extend(&source, &card);
    assert!(OpticalMaterialResponse::try_new(&complete).is_err());
    assert!(
        source
            .with_additional_requirements(
                &base,
                &[ScalarPropertyRequirement::try_new(
                    MU_NAME,
                    MU,
                    ScalarAdmissibility::StrictlyPositive
                )
                .unwrap(),],
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
    );
}

#[test]
fn glass_exposure_reuses_original_film_and_exit_sources_and_refuses_substitution() {
    let incident = material(1.0 / (M * C0 * C0), M);
    let card = glass_card(Family::Cauchy, true);
    let original = resolve_visible_optical_state_point(
        &card,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = extend(&original, &card);
    let cloned = source.clone();
    let changed = source.at_query_point(&card, &point(310.0)).unwrap();
    assert_eq!(cloned.resolved().identity(), source.resolved().identity());
    assert_ne!(changed.resolved().identity(), source.resolved().identity());
    let grid = AngularSpectrumGrid {
        width: 1,
        height: 1,
        spacing_m: [550e-9, 650e-9],
    };
    let limits = AngularSpectrumLimits {
        maximum_grid_points: 1,
        maximum_workspace_bytes: 1 << 16,
        maximum_fft_work: 1024,
    };
    let exposure_limits = AngularSpectrumExposureLimits {
        maximum_lines: 2,
        maximum_point_line_pairs: 2,
        maximum_point_frequency_pairs: 4,
        maximum_workspace_bytes: 1 << 18,
        maximum_fft_work: 4096,
    };
    let input = [[C64::ONE, C64::new(0.2, -0.1)]];
    fn stack<'a>(
        incident: &'a ResolvedMaterialStatePoint,
        source: &'a VisibleOpticalStatePoint,
        omega: f64,
        as_exit: bool,
    ) -> PlanarDielectricLayers<'a> {
        let response =
            DielectricResponse::OpticalMaterial(OpticalMaterialResponse::try_new(source).unwrap());
        let film = [DielectricLayer {
            response,
            thickness_m: 420e-9,
        }];
        PlanarDielectricLayers::try_new(
            incident,
            if as_exit {
                response
            } else {
                DielectricResponse::Lossless(incident)
            },
            omega,
            if as_exit { &[] } else { &film },
        )
        .unwrap()
    }
    for as_exit in [false, true] {
        let stacks = [530.0, 580.0]
            .map(|wavelength| stack(&incident, &source, WL_FREQUENCY / wavelength, as_exit));
        let lines = [
            PlanarAngularSpectrum::prepare(&stacks[0], grid, &input, limits, || false).unwrap(),
            PlanarAngularSpectrum::prepare(&stacks[1], grid, &input, limits, || false).unwrap(),
        ];
        let exposure =
            AngularSpectrumExposure::try_new(&lines[0], &lines, exposure_limits, || false).unwrap();
        let reflected = exposure
            .field(
                AngularSpectrumSide::Reflected,
                0.0,
                HarmonicWindow::new(0.0, 2e-15).unwrap(),
                || false,
            )
            .unwrap();
        assert!(reflected.normal_power_w().is_finite() && reflected.normal_power_w() < 0.0);
        for foreign in [&cloned, &changed] {
            let foreign_stack = stack(&incident, foreign, WL_FREQUENCY / 530.0, as_exit);
            let foreign_beam =
                PlanarAngularSpectrum::prepare(&foreign_stack, grid, &input, limits, || false)
                    .unwrap();
            assert!(
                AngularSpectrumExposure::try_new(&foreign_beam, &lines, exposure_limits, || false)
                    .is_err()
            );
        }
    }
}
