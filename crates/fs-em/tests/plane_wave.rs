#![cfg(feature = "plane-wave")]
//! Synthetic source-bound plane-wave checks, not material qualification.
use fs_em::plane_wave::{DielectricResponse, PlanarDielectricInterface, Polarization};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS as EPS, ELECTRIC_PERMITTIVITY_PROPERTY as EPS_NAME,
        MAGNETIC_PERMEABILITY_DIMS as MU, MAGNETIC_PERMEABILITY_PROPERTY as MU_NAME,
        dielectric::{DebyePoleProperties, ResolvedDebyeCapacitor},
    },
    state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point,
    },
};
use fs_math::c64::C64;
use fs_qty::{Dims, Time};

const E: f64 = 1e-11;
const M: f64 = 1e-6;

#[path = "plane_wave/illumination.rs"]
mod illumination;

fn source(
    properties: &[(&str, Dims, f64)],
    validity: ValidityDomain,
    point: QueryPoint,
) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for &(name, dims, value) in properties {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: validity.clone(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic plane-wave constitutive fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
                .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic dielectric".into(),
            phase: "solid".into(),
            process: "analytic scattering fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn material(epsilon: f64, mu: f64) -> ResolvedMaterialStatePoint {
    source(
        &[(EPS_NAME, EPS, epsilon), (MU_NAME, MU, mu)],
        ValidityDomain::unconstrained().with("T", 280.0, 320.0),
        QueryPoint::new().with("T", 300.0).unwrap(),
    )
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a} != {b}");
}
fn complex_near(a: C64, b: C64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:?} != {b:?}");
}

#[test]
fn normal_incidence_and_magnetic_impedance_use_the_original_sources() {
    let incident = material(E, M);
    let transmitted = material(4.0 * E, M);
    let interface =
        PlanarDielectricInterface::try_new(&incident, DielectricResponse::Lossless(&transmitted))
            .unwrap();
    assert!(std::ptr::eq(interface.incident_material(), &incident));
    assert!(std::ptr::eq(interface.transmitted_material(), &transmitted));
    for p in [Polarization::Te, Polarization::Tm] {
        let r = interface.response(1e9, 1.0, p).unwrap();
        complex_near(r.reflected_tangential_e, C64::from_re(-1.0 / 3.0), 1e-14);
        complex_near(r.transmitted_tangential_e, C64::from_re(2.0 / 3.0), 1e-14);
        near(r.reflected_power_fraction, 1.0 / 9.0, 1e-14);
        near(r.transmitted_power_fraction, 8.0 / 9.0, 1e-14);
    }
    // Small positive index must survive the normal-incidence subtraction.
    let unit_medium = material(1.0, 1.0);
    let low_index = material(1e-20, 1.0);
    let low =
        PlanarDielectricInterface::try_new(&unit_medium, DielectricResponse::Lossless(&low_index))
            .unwrap();
    for p in [Polarization::Te, Polarization::Tm] {
        let r = low.response(1.0, 1.0, p).unwrap();
        near(r.transmitted_normal_wavenumber_rad_m.re, 1e-10, 1e-24);
        near(r.transmitted_power_fraction, 4e-10, 1e-18);
    }

    // Equal impedance, different phase velocity. Assuming mu=mu0 would fail.
    let magnetic = material(4.0 * E, 4.0 * M);
    let matched =
        PlanarDielectricInterface::try_new(&incident, DielectricResponse::Lossless(&magnetic))
            .unwrap();
    for p in [Polarization::Te, Polarization::Tm] {
        let r = matched.response(1e9, 1.0, p).unwrap();
        complex_near(r.reflected_tangential_e, C64::ZERO, 1e-14);
        near(r.transmitted_power_fraction, 1.0, 1e-14);
        near(
            r.transmitted_normal_wavenumber_rad_m.re / r.incident_normal_wavenumber_rad_m,
            4.0,
            1e-14,
        );
    }
}

#[test]
fn brewster_zero_and_lossless_poynting_flux_are_polarization_specific() {
    let incident = material(E, M);
    let transmitted = material(4.0 * E, M);
    let interface =
        PlanarDielectricInterface::try_new(&incident, DielectricResponse::Lossless(&transmitted))
            .unwrap();
    let cosine = 1.0 / 5.0_f64.sqrt(); // tan(theta_B)=n_t/n_i=2
    let tm = interface.response(1e9, cosine, Polarization::Tm).unwrap();
    let te = interface.response(1e9, cosine, Polarization::Te).unwrap();
    assert!(tm.reflected_tangential_e.abs() < 2e-15);
    near(tm.transmitted_power_fraction, 1.0, 1e-14);
    near(te.reflected_power_fraction, 0.36, 1e-14);
    for c in [0.01, 0.3, 0.8, 1.0] {
        for p in [Polarization::Te, Polarization::Tm] {
            let r = interface.response(1e9, c, p).unwrap();
            let ki = r.incident_normal_wavenumber_rad_m;
            let kt = r.transmitted_normal_wavenumber_rad_m.re;
            let ratio = match p {
                Polarization::Te => kt / ki,
                Polarization::Tm => 4.0 * ki / kt,
            };
            near(
                r.transmitted_power_fraction,
                r.transmitted_tangential_e.norm_sq() * ratio,
                2e-14,
            );
            near(
                r.reflected_power_fraction + r.transmitted_power_fraction,
                1.0,
                2e-14,
            );
        }
    }
}

#[test]
fn total_internal_reflection_keeps_phase_and_a_decaying_evanescent_field() {
    let incident = material(4.0 * E, M);
    let transmitted = material(E, M);
    let interface =
        PlanarDielectricInterface::try_new(&incident, DielectricResponse::Lossless(&transmitted))
            .unwrap();
    for p in [Polarization::Te, Polarization::Tm] {
        let r = interface.response(1e9, 0.6, p).unwrap();
        near(r.reflected_power_fraction, 1.0, 2e-14);
        near(r.transmitted_power_fraction, 0.0, 1e-14);
        assert_eq!(r.transmitted_normal_wavenumber_rad_m.re, 0.0);
        assert!(r.transmitted_normal_wavenumber_rad_m.im < 0.0);
        assert!(
            r.reflected_tangential_e.im.abs() > 0.4,
            "unit magnitude must not erase the TIR phase"
        );
        complex_near(
            r.transmitted_tangential_e,
            C64::ONE + r.reflected_tangential_e,
            1e-14,
        );
        match p {
            Polarization::Te => assert!(r.reflected_tangential_e.im > 0.0),
            Polarization::Tm => assert!(r.reflected_tangential_e.im < 0.0),
        }
    }
    // Exact dyadic critical angle: (n_t/n_i)^2=3/4 and cos(theta)=1/2.
    let critical_medium = material(3.0 * E, M);
    let critical = PlanarDielectricInterface::try_new(
        &incident,
        DielectricResponse::Lossless(&critical_medium),
    )
    .unwrap();
    for (p, expected) in [(Polarization::Te, 1.0), (Polarization::Tm, -1.0)] {
        let r = critical.response(1e9, 0.5, p).unwrap();
        complex_near(r.reflected_tangential_e, C64::from_re(expected), 1e-14);
        near(r.transmitted_power_fraction, 0.0, 1e-14);
    }
}

#[test]
fn identical_media_remain_transparent_arbitrarily_close_to_grazing() {
    let medium = material(E, M);
    let interface =
        PlanarDielectricInterface::try_new(&medium, DielectricResponse::Lossless(&medium)).unwrap();
    for cosine in [1.0, 0.5, 1e-10, 1e-200] {
        for p in [Polarization::Te, Polarization::Tm] {
            let r = interface.response(1e9, cosine, p).unwrap();
            assert_eq!(r.reflected_tangential_e, C64::ZERO);
            assert_eq!(r.transmitted_tangential_e, C64::ONE);
            near(r.transmitted_power_fraction, 1.0, 1e-14);
        }
    }
}

#[test]
fn original_debye_response_supplies_phase_flux_and_geometry_independent_bulk_epsilon() {
    let incident = material(E, M);
    let state = source(
        &[
            (EPS_NAME, EPS, 2.0 * E),
            (MU_NAME, MU, M),
            ("delta", EPS, 2.0 * E),
            ("tau", Time::DIMS, 1e-9),
        ],
        ValidityDomain::unconstrained().with("T", 280.0, 320.0),
        QueryPoint::new().with("T", 300.0).unwrap(),
    );
    let poles = [DebyePoleProperties {
        permittivity_increment: "delta".into(),
        relaxation_time: "tau".into(),
    }];
    let a = ResolvedDebyeCapacitor::try_from_material(state.clone(), EPS_NAME, &poles, 0.02, 0.001)
        .unwrap();
    let b = ResolvedDebyeCapacitor::try_from_material(state, EPS_NAME, &poles, 0.3, 0.02).unwrap();
    let ia = PlanarDielectricInterface::try_new(&incident, DielectricResponse::Debye(&a)).unwrap();
    let ib = PlanarDielectricInterface::try_new(&incident, DielectricResponse::Debye(&b)).unwrap();
    assert!(std::ptr::eq(ia.transmitted_material(), a.material()));
    for p in [Polarization::Te, Polarization::Tm] {
        for c in [0.2, 0.7, 1.0] {
            let r = ia.response(1e9, c, p).unwrap();
            let other = ib.response(1e9, c, p).unwrap();
            complex_near(r.transmitted_permittivity_f_m, C64::new(3.0 * E, -E), 1e-24);
            complex_near(
                r.reflected_tangential_e,
                other.reflected_tangential_e,
                2e-14,
            );
            if c == 1.0 {
                // Independent 80-digit Cartesian evaluation of sqrt(3-i).
                complex_near(
                    r.reflected_tangential_e,
                    C64::new(-0.28180664470320127, 0.07424789305525832),
                    2e-14,
                );
            }
            assert!(r.transmitted_normal_wavenumber_rad_m.re > 0.0);
            assert!(r.transmitted_normal_wavenumber_rad_m.im < 0.0);
            let ki = r.incident_normal_wavenumber_rad_m;
            let kt = r.transmitted_normal_wavenumber_rad_m;
            let (yi, yt) = match p {
                Polarization::Te => (ki / (1e9 * M), kt.scale(1.0 / (1e9 * M))),
                Polarization::Tm => (1e9 * E / ki, r.transmitted_permittivity_f_m.scale(1e9) / kt),
            };
            let flux = r.transmitted_tangential_e.norm_sq() * yt.re / yi;
            near(r.transmitted_power_fraction, flux, 2e-14);
            near(r.reflected_power_fraction + flux, 1.0, 2e-14);
            assert!(flux > 0.0 && flux < 1.0);
        }
    }
}

#[test]
fn missing_magnetic_source_and_invalid_incidence_refuse() {
    let medium = material(E, M);
    let missing_mu = source(
        &[(EPS_NAME, EPS, E)],
        ValidityDomain::unconstrained(),
        QueryPoint::new(),
    );
    assert!(
        PlanarDielectricInterface::try_new(&medium, DielectricResponse::Lossless(&missing_mu),)
            .is_err()
    );
    // A frequency-limited static scalar is not promoted into this explicit
    // nondispersive law just because its current query happens to be in-band.
    let bounded = source(
        &[(EPS_NAME, EPS, E), (MU_NAME, MU, M)],
        ValidityDomain::unconstrained().with("omega", 0.0, 1e9),
        QueryPoint::new().with("omega", 5e8).unwrap(),
    );
    assert!(
        PlanarDielectricInterface::try_new(&medium, DielectricResponse::Lossless(&bounded),)
            .is_err()
    );
    // Positive finite source values must not underflow into a zero-index medium.
    let extreme_incident = material(1e300, 1e-300);
    let extreme_transmitted = material(1e-300, 1e-300);
    let extreme = PlanarDielectricInterface::try_new(
        &extreme_incident,
        DielectricResponse::Lossless(&extreme_transmitted),
    )
    .unwrap();
    assert!(extreme.response(1e9, 1.0, Polarization::Te).is_err());
    assert!(extreme.response(1e9, 1.0, Polarization::Tm).is_err());
    let interface =
        PlanarDielectricInterface::try_new(&medium, DielectricResponse::Lossless(&medium)).unwrap();
    for (omega, cosine) in [
        (0.0, 1.0),
        (f64::NAN, 1.0),
        (1e9, 0.0),
        (1e9, -0.1),
        (1e9, 1.1),
        (1e9, f64::NAN),
    ] {
        assert!(interface.response(omega, cosine, Polarization::Te).is_err());
    }
}

mod finite_layers {
    use super::*;
    use fs_em::plane_wave::layers::{DielectricLayer, MAX_LAYERS, PlanarDielectricLayers};
    use fs_math::det;
    use std::f64::consts::PI;

    pub(super) fn debye(epsilon_inf: f64, delta: f64, area: f64, gap: f64) -> ResolvedDebyeCapacitor {
        let state = source(
            &[
                (EPS_NAME, EPS, epsilon_inf),
                (MU_NAME, MU, M),
                ("delta", EPS, delta),
                ("tau", Time::DIMS, 1e-9),
            ],
            ValidityDomain::unconstrained().with("T", 280.0, 320.0),
            QueryPoint::new().with("T", 300.0).unwrap(),
        );
        ResolvedDebyeCapacitor::try_from_material(
            state,
            EPS_NAME,
            &[DebyePoleProperties {
                permittivity_increment: "delta".into(),
                relaxation_time: "tau".into(),
            }],
            area,
            gap,
        )
        .unwrap()
    }

    // Independent conventional E/H transfer matrices. These modest-depth
    // fixtures can form cos/sin directly; production must also handle opaque
    // films for which the growing entries here would overflow.
    fn transfer_reference(
        films: &[(C64, f64, f64)],
        exit: (C64, f64),
        cosine: f64,
        polarization: Polarization,
    ) -> (C64, C64, f64) {
        let normal = |epsilon: C64, mu: f64| {
            let index = epsilon.scale(mu);
            // Preserve the known incident/exit q=cosine and grazing c^2;
            // subtracting two rounded unit terms would spoil this oracle.
            let real = if cosine > 0.5 {
                index.re - (1.0 - cosine) * (1.0 + cosine)
            } else {
                (index.re - 1.0) + cosine * cosine
            };
            let mut q = if index == C64::ONE {
                C64::from_re(cosine)
            } else {
                C64::new(real, index.im).sqrt()
            };
            if q.im > 0.0 {
                q = -q;
            }
            q
        };
        let j = C64::new(0.0, 1.0);
        let mut matrix = [[C64::ONE, C64::ZERO], [C64::ZERO, C64::ONE]];
        for &(epsilon, mu, depth) in films {
            let q = normal(epsilon, mu);
            let z = q.scale(depth);
            let (s, c) = det::sin_cos(z.re);
            let ch = 0.5 * (det::exp(z.im) + det::exp(-z.im));
            let sh = 0.5 * (det::exp(z.im) - det::exp(-z.im));
            let cos = C64::new(c * ch, -s * sh);
            let sin_over_q = if q == C64::ZERO {
                C64::from_re(depth)
            } else {
                C64::new(s * ch, c * sh) / q
            };
            let (upper, lower) = match polarization {
                Polarization::Te => (
                    j * sin_over_q.scale(mu),
                    j * q * q * sin_over_q.scale(1.0 / mu),
                ),
                Polarization::Tm => (j * q * q * sin_over_q / epsilon, j * epsilon * sin_over_q),
            };
            let slab = [[cos, upper], [lower, cos]];
            let old = matrix;
            for row in 0..2 {
                for column in 0..2 {
                    matrix[row][column] =
                        old[row][0] * slab[0][column] + old[row][1] * slab[1][column];
                }
            }
        }
        let q_exit = normal(exit.0, exit.1);
        let (yi, yt) = match polarization {
            Polarization::Te => (cosine, q_exit.scale(1.0 / exit.1)),
            Polarization::Tm => (1.0 / cosine, exit.0 / q_exit),
        };
        let u = matrix[0][0] + matrix[0][1] * yt;
        let v = matrix[1][0] + matrix[1][1] * yt;
        let denominator = u.scale(yi) + v;
        let r = (u.scale(yi) - v) / denominator;
        let t = C64::from_re(2.0 * yi) / denominator;
        (r, t, t.norm_sq() * yt.re / yi)
    }

    #[test]
    fn empty_and_zero_thickness_films_keep_the_halfspace_sources_and_response() {
        let incident = material(E, M);
        let exit = material(4.0 * E, M);
        let film = material(7.0 * E, 2.0 * M);
        let layers = [DielectricLayer {
            response: DielectricResponse::Lossless(&film),
            thickness_m: 0.0,
        }];
        let bare =
            PlanarDielectricInterface::try_new(&incident, DielectricResponse::Lossless(&exit))
                .unwrap();
        for films in [&[][..], &layers[..]] {
            let stack = PlanarDielectricLayers::try_new(
                &incident,
                DielectricResponse::Lossless(&exit),
                1e9,
                films,
            )
            .unwrap();
            assert!(std::ptr::eq(stack.incident_material(), &incident));
            assert!(std::ptr::eq(stack.exit_material(), &exit));
            assert_eq!(stack.angular_frequency_rad_s(), 1e9);
            assert_eq!(stack.layer_count(), films.len());
            if !films.is_empty() {
                assert!(std::ptr::eq(
                    stack.layer(0).unwrap().response.material(),
                    &film
                ));
            }
            for c in [1.0, 0.5, 0.01, 1e-100] {
                for pol in [Polarization::Te, Polarization::Tm] {
                    let actual = stack.response(c, pol).unwrap();
                    let expected = bare.response(1e9, c, pol).unwrap();
                    assert_eq!(
                        actual.reflected_tangential_e,
                        expected.reflected_tangential_e
                    );
                    assert_eq!(
                        actual.transmitted_tangential_e,
                        expected.transmitted_tangential_e
                    );
                    assert_eq!(
                        actual.transmitted_power_fraction,
                        expected.transmitted_power_fraction
                    );
                    assert_eq!(actual.absorbed_power_fraction, 0.0);
                }
            }
        }
    }

    #[test]
    fn quarter_wave_antireflection_and_fabry_perot_keep_complex_film_phase() {
        let incident = material(E, M);
        let film = material(4.0 * E, M);
        let substrate = material(16.0 * E, M);
        let k = 1e9 * det::sqrt(E) * det::sqrt(M);
        for pol in [Polarization::Te, Polarization::Tm] {
            let layers = [DielectricLayer {
                response: DielectricResponse::Lossless(&film),
                thickness_m: PI / (4.0 * k),
            }];
            let ar = PlanarDielectricLayers::try_new(
                &incident,
                DielectricResponse::Lossless(&substrate),
                1e9,
                &layers,
            )
            .unwrap()
            .response(1.0, pol)
            .unwrap();
            complex_near(ar.reflected_tangential_e, C64::ZERO, 5e-14);
            complex_near(ar.transmitted_tangential_e, C64::new(0.0, -0.5), 5e-14);
            near(ar.transmitted_power_fraction, 1.0, 5e-14);
            let quarter = PlanarDielectricLayers::try_new(
                &incident,
                DielectricResponse::Lossless(&incident),
                1e9,
                &layers,
            )
            .unwrap()
            .response(1.0, pol)
            .unwrap();
            complex_near(quarter.reflected_tangential_e, C64::from_re(-0.6), 5e-14);
            complex_near(quarter.transmitted_tangential_e, C64::new(0.0, -0.8), 5e-14);
            let half = [DielectricLayer {
                thickness_m: PI / (2.0 * k),
                ..layers[0]
            }];
            let resonant = PlanarDielectricLayers::try_new(
                &incident,
                DielectricResponse::Lossless(&incident),
                1e9,
                &half,
            )
            .unwrap()
            .response(1.0, pol)
            .unwrap();
            complex_near(resonant.reflected_tangential_e, C64::ZERO, 5e-14);
            complex_near(resonant.transmitted_tangential_e, -C64::ONE, 5e-14);
            assert_eq!(resonant.absorbed_power_fraction, 0.0);
        }
    }

    #[test]
    fn grazing_resonances_keep_subdivision_and_distinct_film_continuity() {
        let incident = material(1.0, 1.0);
        let first = material(2.0, 1.0);
        for c in [0.1_f64, 0.001, 1e-7, 1e-8] {
            let depth = PI / (2.0 * det::sqrt(1.0 + c * c));
            let quarter = DielectricLayer {
                response: DielectricResponse::Lossless(&first),
                thickness_m: depth,
            };
            let split = PlanarDielectricLayers::try_new(
                &incident,
                DielectricResponse::Lossless(&incident),
                1.0,
                &[quarter; 2],
            )
            .unwrap();
            let whole = PlanarDielectricLayers::try_new(
                &incident,
                DielectricResponse::Lossless(&incident),
                1.0,
                &[DielectricLayer {
                    thickness_m: 2.0 * depth,
                    ..quarter
                }],
            )
            .unwrap();

            // Different original epsilon AND mu, with matched TE film
            // admittances. Coalescing identical films cannot fix this case.
            let epsilon_second = 2.5 + 1.5 * c * c;
            let second = material(epsilon_second, 2.0);
            let second_normal = det::sqrt((2.0 * epsilon_second - 1.0) + c * c);
            let second_depth = PI / (2.0 * second_normal);
            let distinct = PlanarDielectricLayers::try_new(
                &incident,
                DielectricResponse::Lossless(&incident),
                1.0,
                &[
                    quarter,
                    DielectricLayer {
                        response: DielectricResponse::Lossless(&second),
                        thickness_m: second_depth,
                    },
                ],
            )
            .unwrap();
            for pol in [Polarization::Te, Polarization::Tm] {
                let divided = split.response(c, pol).unwrap();
                let single = whole.response(c, pol).unwrap();
                complex_near(
                    divided.reflected_tangential_e,
                    single.reflected_tangential_e,
                    2e-13,
                );
                complex_near(
                    divided.transmitted_tangential_e,
                    single.transmitted_tangential_e,
                    2e-13,
                );
                near(
                    divided.transmitted_power_fraction,
                    single.transmitted_power_fraction,
                    2e-13,
                );
                near(
                    divided.reflected_power_fraction + divided.transmitted_power_fraction,
                    1.0,
                    2e-13,
                );
                assert_eq!(divided.absorbed_power_fraction, 0.0);

                let actual = distinct.response(c, pol).unwrap();
                let expected = transfer_reference(
                    &[
                        (C64::from_re(2.0), 1.0, depth),
                        (C64::from_re(epsilon_second), 2.0, second_depth),
                    ],
                    (C64::ONE, 1.0),
                    c,
                    pol,
                );
                complex_near(actual.reflected_tangential_e, expected.0, 2e-13);
                complex_near(actual.transmitted_tangential_e, expected.1, 2e-13);
                near(actual.transmitted_power_fraction, expected.2, 2e-13);
                near(
                    actual.reflected_power_fraction + actual.transmitted_power_fraction,
                    1.0,
                    2e-13,
                );
                assert_eq!(actual.absorbed_power_fraction, 0.0);
            }
        }
    }

    #[test]
    fn oblique_magnetic_multilayers_match_independent_transfer_matrices() {
        let incident = material(1.0, 1.0);
        let first = material(4.0, 1.0);
        let second = material(0.75, 1.0);
        let third = material(2.25, 2.0);
        let exit = material(1.44, 0.8);
        let layers = [
            DielectricLayer {
                response: DielectricResponse::Lossless(&first),
                thickness_m: 0.17,
            },
            DielectricLayer {
                response: DielectricResponse::Lossless(&second),
                thickness_m: 0.11,
            },
            DielectricLayer {
                response: DielectricResponse::Lossless(&third),
                thickness_m: 0.09,
            },
        ];
        let stack = PlanarDielectricLayers::try_new(
            &incident,
            DielectricResponse::Lossless(&exit),
            1.0,
            &layers,
        )
        .unwrap();
        for c in [0.2, 0.5, 0.8, 1.0] {
            for pol in [Polarization::Te, Polarization::Tm] {
                let actual = stack.response(c, pol).unwrap();
                let expected = transfer_reference(
                    &[
                        (C64::from_re(4.0), 1.0, 0.17),
                        (C64::from_re(0.75), 1.0, 0.11),
                        (C64::from_re(2.25), 2.0, 0.09),
                    ],
                    (C64::from_re(1.44), 0.8),
                    c,
                    pol,
                );
                complex_near(actual.reflected_tangential_e, expected.0, 2e-13);
                complex_near(actual.transmitted_tangential_e, expected.1, 2e-13);
                near(actual.transmitted_power_fraction, expected.2, 2e-13);
                near(
                    actual.reflected_power_fraction + actual.transmitted_power_fraction,
                    1.0,
                    2e-13,
                );
                assert_eq!(actual.absorbed_power_fraction, 0.0);
            }
        }
    }

    #[test]
    fn finite_debye_loss_keeps_original_law_and_is_separate_from_exit_flux() {
        let incident = material(E, M);
        let exit = material(2.25 * E, M);
        let first = debye(2.0 * E, 2.0 * E, 0.02, 0.001);
        let other = debye(2.0 * E, 2.0 * E, 0.3, 0.02);
        let k = 1e9 * det::sqrt(E) * det::sqrt(M);
        for pol in [Polarization::Te, Polarization::Tm] {
            let mut responses = Vec::new();
            for model in [&first, &other] {
                let films = [DielectricLayer {
                    response: DielectricResponse::Debye(model),
                    thickness_m: 0.4 / k,
                }];
                let stack = PlanarDielectricLayers::try_new(
                    &incident,
                    DielectricResponse::Lossless(&exit),
                    1e9,
                    &films,
                )
                .unwrap();
                assert!(std::ptr::eq(
                    stack.layer(0).unwrap().response.material(),
                    model.material()
                ));
                complex_near(
                    stack.layer_permittivity_f_m(0).unwrap(),
                    C64::new(3.0 * E, -E),
                    1e-24,
                );
                let actual = stack.response(0.7, pol).unwrap();
                let reference = transfer_reference(
                    &[(C64::new(3.0, -1.0), 1.0, 0.4)],
                    (C64::from_re(2.25), 1.0),
                    0.7,
                    pol,
                );
                complex_near(actual.reflected_tangential_e, reference.0, 2e-13);
                complex_near(actual.transmitted_tangential_e, reference.1, 2e-13);
                near(actual.transmitted_power_fraction, reference.2, 2e-13);
                near(
                    actual.absorbed_power_fraction,
                    1.0 - reference.0.norm_sq() - reference.2,
                    2e-13,
                );
                assert!(actual.absorbed_power_fraction > 0.1);
                responses.push(actual);
            }
            complex_near(
                responses[0].reflected_tangential_e,
                responses[1].reflected_tangential_e,
                2e-14,
            );
            let no_film = PlanarDielectricLayers::try_new(
                &incident,
                DielectricResponse::Debye(&first),
                1e9,
                &[],
            )
            .unwrap()
            .response(0.7, pol)
            .unwrap();
            assert_eq!(no_film.absorbed_power_fraction, 0.0);
            assert!(no_film.transmitted_power_fraction > 0.5);
            let opaque = [DielectricLayer {
                response: DielectricResponse::Debye(&first),
                thickness_m: 1e6 / k,
            }];
            let thick = PlanarDielectricLayers::try_new(
                &incident,
                DielectricResponse::Lossless(&exit),
                1e9,
                &opaque,
            )
            .unwrap()
            .response(0.7, pol)
            .unwrap();
            complex_near(
                thick.reflected_tangential_e,
                no_film.reflected_tangential_e,
                2e-14,
            );
            assert_eq!(thick.transmitted_power_fraction, 0.0);
            near(
                thick.absorbed_power_fraction,
                no_film.transmitted_power_fraction,
                2e-14,
            );
        }
    }

    #[test]
    fn weak_ultrathin_debye_absorption_survives_power_subtraction_roundoff() {
        let incident = material(E, M);
        let weak = debye(3.0 * E, 2e-20 * E, 0.02, 0.001);
        let k = 1e9 * det::sqrt(E) * det::sqrt(M);
        let films = [DielectricLayer {
            response: DielectricResponse::Debye(&weak),
            thickness_m: 1e-10 / k,
        }];
        let stack = PlanarDielectricLayers::try_new(
            &incident,
            DielectricResponse::Lossless(&incident),
            1e9,
            &films,
        )
        .unwrap();
        for pol in [Polarization::Te, Polarization::Tm] {
            let actual = stack.response(1.0, pol).unwrap();
            // Integral of omega*(-Im epsilon)*|E|^2 / Y_inc across this
            // thin film. Corrections from its O(d) field are below the band.
            near(actual.absorbed_power_fraction / 1e-30, 1.0, 1e-10);
        }
    }

    #[test]
    fn internal_and_exit_critical_angles_and_opaque_tunnelling_remain_finite() {
        let incident = material(1.0, 1.0);
        let critical = material(0.75, 1.0);
        let dense = material(4.0, 1.0);
        let gap = material(0.25, 1.0);
        for pol in [Polarization::Te, Polarization::Tm] {
            let films = [DielectricLayer {
                response: DielectricResponse::Lossless(&critical),
                thickness_m: 4.0,
            }];
            let at_critical = PlanarDielectricLayers::try_new(
                &incident,
                DielectricResponse::Lossless(&incident),
                1.0,
                &films,
            )
            .unwrap()
            .response(0.5, pol)
            .unwrap();
            let expected =
                transfer_reference(&[(C64::from_re(0.75), 1.0, 4.0)], (C64::ONE, 1.0), 0.5, pol);
            complex_near(at_critical.reflected_tangential_e, expected.0, 2e-13);
            complex_near(at_critical.transmitted_tangential_e, expected.1, 2e-13);
            let films = [DielectricLayer {
                response: DielectricResponse::Lossless(&dense),
                thickness_m: 0.3,
            }];
            let exit_critical = PlanarDielectricLayers::try_new(
                &incident,
                DielectricResponse::Lossless(&critical),
                1.0,
                &films,
            )
            .unwrap()
            .response(0.5, pol)
            .unwrap();
            near(exit_critical.reflected_power_fraction, 1.0, 2e-13);
            assert_eq!(exit_critical.transmitted_power_fraction, 0.0);
            for depth in [0.2, 2.0, 20.0, 1e6] {
                let films = [DielectricLayer {
                    response: DielectricResponse::Lossless(&gap),
                    thickness_m: depth,
                }];
                let actual = PlanarDielectricLayers::try_new(
                    &incident,
                    DielectricResponse::Lossless(&incident),
                    1.0,
                    &films,
                )
                .unwrap()
                .response(0.6, pol)
                .unwrap();
                assert_eq!(actual.absorbed_power_fraction, 0.0);
                if depth < 1e6 {
                    let expected = transfer_reference(
                        &[(C64::from_re(0.25), 1.0, depth)],
                        (C64::ONE, 1.0),
                        0.6,
                        pol,
                    );
                    complex_near(actual.reflected_tangential_e, expected.0, 2e-13);
                    complex_near(actual.transmitted_tangential_e, expected.1, 2e-13);
                    near(actual.transmitted_power_fraction, expected.2, 2e-13);
                    assert!(actual.transmitted_power_fraction > 0.0);
                } else {
                    assert_eq!(actual.transmitted_tangential_e, C64::ZERO);
                    near(actual.reflected_power_fraction, 1.0, 2e-13);
                }
            }
        }
    }

    #[test]
    fn finite_layer_caps_geometry_sources_and_incidence_refuse_explicitly() {
        let incident = material(E, M);
        let valid = DielectricLayer {
            response: DielectricResponse::Lossless(&incident),
            thickness_m: 0.1,
        };
        let too_many = [valid; MAX_LAYERS + 1];
        assert!(
            PlanarDielectricLayers::try_new(
                &incident,
                DielectricResponse::Lossless(&incident),
                1e9,
                &too_many,
            )
            .is_err()
        );
        for thickness in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(
                PlanarDielectricLayers::try_new(
                    &incident,
                    DielectricResponse::Lossless(&incident),
                    1e9,
                    &[DielectricLayer {
                        thickness_m: thickness,
                        ..valid
                    }],
                )
                .is_err()
            );
        }
        let missing_mu = source(
            &[(EPS_NAME, EPS, E)],
            ValidityDomain::unconstrained(),
            QueryPoint::new(),
        );
        assert!(
            PlanarDielectricLayers::try_new(
                &incident,
                DielectricResponse::Lossless(&incident),
                1e9,
                &[DielectricLayer {
                    response: DielectricResponse::Lossless(&missing_mu),
                    thickness_m: 0.0
                }],
            )
            .is_err()
        );
        assert!(
            PlanarDielectricLayers::try_new(
                &incident,
                DielectricResponse::Lossless(&incident),
                0.0,
                &[valid],
            )
            .is_err()
        );
        let stack = PlanarDielectricLayers::try_new(
            &incident,
            DielectricResponse::Lossless(&incident),
            1e9,
            &[valid],
        )
        .unwrap();
        for c in [0.0, -0.1, 1.01, f64::NAN] {
            assert!(stack.response(c, Polarization::Te).is_err());
        }
        let dense = material(4.0 * E, M);
        let contrast = PlanarDielectricLayers::try_new(
            &incident,
            DielectricResponse::Lossless(&incident),
            1e9,
            &[DielectricLayer {
                response: DielectricResponse::Lossless(&dense),
                thickness_m: 0.1,
            }],
        )
        .unwrap();
        for pol in [Polarization::Te, Polarization::Tm] {
            assert!(contrast.response(1e-200, pol).is_err());
        }
    }
}
