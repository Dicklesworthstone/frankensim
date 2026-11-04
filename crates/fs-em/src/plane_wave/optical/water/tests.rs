//! Native source-to-interface/film/beam controls. Permeabilities and material
//! cards here are authored fixtures, not measurements or water qualification.
use super::*;
use crate::plane_wave::{
    DielectricResponse, PlanarDielectricInterface, Polarization,
    angular_spectrum::{
        AngularSpectrumExposure, AngularSpectrumExposureLayerLimits, AngularSpectrumExposureLimits,
        AngularSpectrumGrid, AngularSpectrumLimits, AngularSpectrumSide, PlanarAngularSpectrum,
    },
    layers::{DielectricLayer, PlanarDielectricLayers},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS,
        MAGNETIC_PERMEABILITY_PROPERTY,
    },
    state_point::{
        MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
        resolve_material_state_point,
    },
};
use fs_qty::{Dims, semantic::SemanticType};

const MU: f64 = 1.256_637_062_12e-6;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}; tolerance {tolerance:.3e}"
    );
}

fn complex_near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:?} != {expected:?}; tolerance {tolerance:.3e}"
    );
}

fn phase(angle: f64) -> C64 {
    let (sine, cosine) = angle.sin_cos();
    C64::new(cosine, sine)
}

fn snapshot(
    properties: &[(&str, Dims, f64)],
    point: QueryPoint,
    validity: ValidityDomain,
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
                    source: "authored water-optics Maxwell control".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "explicitly selected ordinary-water optical control".into(),
            phase: "caller-owned homogeneous state".into(),
            process: "source-admission and analytic wave control".into(),
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

fn water_state(temperature: f64, density: f64, mu: f64) -> ResolvedMaterialStatePoint {
    snapshot(
        &[
            (DENSITY_PROPERTY, Density::DIMS, density),
            (
                MAGNETIC_PERMEABILITY_PROPERTY,
                MAGNETIC_PERMEABILITY_DIMS,
                mu,
            ),
        ],
        QueryPoint::new().with("T", temperature).unwrap(),
        ValidityDomain::unconstrained().with("T", 200.0, 1000.0),
    )
}

fn unit_index_port() -> ResolvedMaterialStatePoint {
    snapshot(
        &[
            (
                ELECTRIC_PERMITTIVITY_PROPERTY,
                ELECTRIC_PERMITTIVITY_DIMS,
                1.0 / (MU * VACUUM_SPEED_M_S.powi(2)),
            ),
            (
                MAGNETIC_PERMEABILITY_PROPERTY,
                MAGNETIC_PERMEABILITY_DIMS,
                MU,
            ),
        ],
        QueryPoint::new(),
        ValidityDomain::unconstrained(),
    )
}

fn beam_limits() -> AngularSpectrumLimits {
    AngularSpectrumLimits {
        maximum_grid_points: 32,
        maximum_workspace_bytes: 2_000_000,
        maximum_fft_work: 100_000,
    }
}

// Independent Maxwell admittance for a real normal wave number k_z=k0*q.
// The fixed tangential-E convention applies to both TE and TM.
fn admittance(polarization: Polarization, n: f64, q: f64, mu: f64) -> f64 {
    match polarization {
        Polarization::Te => q / (VACUUM_SPEED_M_S * mu),
        Polarization::Tm => n * n / (q * VACUUM_SPEED_M_S * mu),
    }
}

#[test]
fn water_source_keeps_published_absolute_index_group_delay_and_si_mu() {
    // IAPWS R9-97 Table 3: n(273.15 K, 0.1 MPa, 589 nm)=1.334344.
    // Independent density 999.84 kg/m3 is NISTIR 5078 Table 3 (printed p22),
    // rounded to 0.01 kg/m3. The tolerance includes that rounding; no density
    // was inferred by fitting the expected refractive index.
    // https://iapws.org/technical-guidance/release/Rindex.download
    // https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nistir5078.pdf
    for mu in [MU, 2.3 * MU] {
        let source = water_state(273.15, 999.84, mu);
        let law = WaterOpticalResponse::try_lossless(&source, "T").unwrap();
        assert!(core::ptr::eq(law.source(), &source));
        assert_eq!(law.material().identity(), source.identity());
        assert_eq!(
            (
                law.temperature_k(),
                law.density_kg_m3(),
                law.permeability_h_m()
            ),
            (273.15, 999.84, mu)
        );
        let omega = FREQUENCY_WAVELENGTH_NM / 589.0;
        let optical = law.refractive_sample(omega).unwrap();
        near(optical.phase_index, 1.334344, 3.5e-6);
        let epsilon = law.permittivity_f_m(omega).unwrap();
        assert_eq!(epsilon.im, 0.0);
        assert_eq!(
            law.complex_index(omega).unwrap(),
            C64::from_re(optical.phase_index)
        );
        near(
            epsilon.re * mu * VACUUM_SPEED_M_S.powi(2),
            optical.phase_index.powi(2),
            1e-15,
        );
        assert!(optical.group_index > optical.phase_index + 0.01);
        assert!(
            (epsilon.re * mu * VACUUM_SPEED_M_S.powi(2) - optical.group_index.powi(2)).abs() > 0.02
        );
        // Independent frequency differentiation of k=omega*n(omega)/c.
        let h = 1e-5 * omega;
        let k = |w| w * law.complex_index(w).unwrap().re / VACUUM_SPEED_M_S;
        near(
            (k(omega + h) - k(omega - h)) / (2.0 * h) * VACUUM_SPEED_M_S,
            optical.group_index,
            4e-10,
        );
        for wavelength in WAVELENGTH_RANGE_NM {
            law.permittivity_f_m(FREQUENCY_WAVELENGTH_NM / wavelength)
                .unwrap();
        }
        let low = FREQUENCY_WAVELENGTH_NM / WAVELENGTH_RANGE_NM[1];
        let high = FREQUENCY_WAVELENGTH_NM / WAVELENGTH_RANGE_NM[0];
        for bad in [
            0.0,
            -omega,
            f64::INFINITY,
            f64::NAN,
            f64::from_bits(low.to_bits() - 1),
            f64::from_bits(high.to_bits() + 1),
        ] {
            assert!(law.permittivity_f_m(bad).is_err());
        }
    }
}

#[test]
fn water_state_domain_and_original_quantity_sources_are_required() {
    for (temperature, density) in [
        (260.0, 999.0),
        (774.0, 999.0),
        (300.0, -1.0),
        (300.0, 1061.0),
    ] {
        let source = water_state(temperature, density, MU);
        assert!(WaterOpticalResponse::try_lossless(&source, "T").is_err());
    }
    let properties = [
        (DENSITY_PROPERTY, Density::DIMS, 999.84),
        (
            MAGNETIC_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            MU,
        ),
    ];
    let missing_temperature = snapshot(
        &properties,
        QueryPoint::new(),
        ValidityDomain::unconstrained(),
    );
    assert!(WaterOpticalResponse::try_lossless(&missing_temperature, "T").is_err());
    for kind in [
        QuantityKind::AbsoluteTemperature,
        QuantityKind::TemperatureDifference,
    ] {
        let quantity = QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
        let source = snapshot(
            &properties,
            QueryPoint::new()
                .with_quantity("water-temperature", quantity, 300.0)
                .unwrap(),
            ValidityDomain::unconstrained().with_quantity(
                "water-temperature",
                quantity,
                290.0,
                310.0,
            ),
        );
        assert_eq!(
            WaterOpticalResponse::try_lossless(&source, "water-temperature").is_ok(),
            kind == QuantityKind::AbsoluteTemperature
        );
    }
    let t = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let ambiguous = snapshot(
        &properties,
        QueryPoint::new()
            .with("T", 300.0)
            .unwrap()
            .with_quantity("other-temperature", t, 300.0)
            .unwrap(),
        ValidityDomain::unconstrained()
            .with("T", 290.0, 310.0)
            .with_quantity("other-temperature", t, 290.0, 310.0),
    );
    assert!(WaterOpticalResponse::try_lossless(&ambiguous, "T").is_err());
    for incomplete in [&properties[..1], &properties[1..]] {
        let source = snapshot(
            incomplete,
            QueryPoint::new().with("T", 300.0).unwrap(),
            ValidityDomain::unconstrained(),
        );
        assert!(WaterOpticalResponse::try_lossless(&source, "T").is_err());
    }
    let wrong_density = snapshot(
        &[(DENSITY_PROPERTY, Dims::NONE, 999.84), properties[1]],
        QueryPoint::new().with("T", 300.0).unwrap(),
        ValidityDomain::unconstrained(),
    );
    assert!(WaterOpticalResponse::try_lossless(&wrong_density, "T").is_err());
    let frequency_snapshot = snapshot(
        &properties,
        QueryPoint::new()
            .with("T", 300.0)
            .unwrap()
            .with("omega", 1e15)
            .unwrap(),
        ValidityDomain::unconstrained()
            .with("T", 290.0, 310.0)
            .with("omega", 1e14, 1e16),
    );
    assert!(WaterOpticalResponse::try_lossless(&frequency_snapshot, "T").is_err());
}

#[test]
fn water_interfaces_obey_oblique_fresnel_brewster_reciprocity_and_tir() {
    let port = unit_index_port();
    for (temperature, density) in [(273.15, 999.84), (373.15, 962.93), (473.15, 0.46031)] {
        let source = water_state(temperature, density, MU);
        let water = WaterOpticalResponse::try_lossless(&source, "T").unwrap();
        let forward =
            PlanarDielectricInterface::try_new(&port, DielectricResponse::WaterOptical(water))
                .unwrap();
        let reverse = PlanarDielectricInterface::try_new_with_incident_response(
            DielectricResponse::WaterOptical(water),
            DielectricResponse::Lossless(&port),
        )
        .unwrap();
        for wavelength in [226.5, 589.0, 1013.98] {
            let omega = FREQUENCY_WAVELENGTH_NM / wavelength;
            let n = water.complex_index(omega).unwrap().re;
            let k0 = omega / VACUUM_SPEED_M_S;
            for ci in [1.0_f64, 0.8, 0.35] {
                let qt = (n * n - (1.0 - ci * ci)).sqrt();
                for polarization in [Polarization::Te, Polarization::Tm] {
                    let yi = admittance(polarization, 1.0, ci, MU);
                    let yt = admittance(polarization, n, qt, MU);
                    let r = (yi - yt) / (yi + yt);
                    let t = 2.0 * yi / (yi + yt);
                    let result = forward.response(omega, ci, polarization).unwrap();
                    complex_near(result.reflected_tangential_e, C64::from_re(r), 2e-15);
                    complex_near(result.transmitted_tangential_e, C64::from_re(t), 2e-15);
                    near(
                        result.transmitted_normal_wavenumber_rad_m.re / k0,
                        qt,
                        2e-15,
                    );
                    near(result.reflected_power_fraction, r * r, 2e-15);
                    near(result.transmitted_power_fraction, yt / yi * t * t, 3e-15);
                    near(
                        result.reflected_power_fraction + result.transmitted_power_fraction,
                        1.0,
                        2e-15,
                    );
                    let returning = reverse.response(omega, qt / n, polarization).unwrap();
                    complex_near(returning.reflected_tangential_e, C64::from_re(-r), 4e-15);
                    near(
                        returning.transmitted_normal_wavenumber_rad_m.re / k0,
                        ci,
                        4e-15,
                    );
                }
            }
            let brewster_cosine = 1.0 / (1.0 + n * n).sqrt();
            near(
                forward
                    .response(omega, brewster_cosine, Polarization::Tm)
                    .unwrap()
                    .reflected_power_fraction,
                0.0,
                5e-30,
            );
        }
    }
    let source = water_state(273.15, 999.84, MU);
    let water = WaterOpticalResponse::try_lossless(&source, "T").unwrap();
    let omega = FREQUENCY_WAVELENGTH_NM / 589.0;
    let n = water.complex_index(omega).unwrap().re;
    let ci = 0.4_f64;
    let gamma = (n * n * (1.0 - ci * ci) - 1.0).sqrt();
    let reverse = PlanarDielectricInterface::try_new_with_incident_response(
        DielectricResponse::WaterOptical(water),
        DielectricResponse::Lossless(&port),
    )
    .unwrap();
    for polarization in [Polarization::Te, Polarization::Tm] {
        let yi = admittance(polarization, n, n * ci, MU);
        let yt = match polarization {
            Polarization::Te => C64::new(0.0, -gamma / (VACUUM_SPEED_M_S * MU)),
            Polarization::Tm => C64::new(0.0, 1.0 / (gamma * VACUUM_SPEED_M_S * MU)),
        };
        let expected = (C64::from_re(yi) - yt) / (C64::from_re(yi) + yt);
        let reflected = reverse.response(omega, ci, polarization).unwrap();
        complex_near(reflected.reflected_tangential_e, expected, 2e-15);
        near(reflected.reflected_power_fraction, 1.0, 2e-15);
        near(reflected.transmitted_power_fraction, 0.0, 0.0);
        near(
            reflected.transmitted_normal_wavenumber_rad_m.im / (omega / VACUUM_SPEED_M_S),
            -gamma,
            1e-15,
        );
    }
}

#[test]
fn water_films_match_independent_fabry_perot_and_interior_maxwell_fields() {
    let port = unit_index_port();
    for (temperature, density) in [(273.15, 999.84), (373.15, 962.93)] {
        let source = water_state(temperature, density, MU);
        let water = WaterOpticalResponse::try_lossless(&source, "T").unwrap();
        for wavelength in [226.5, 589.0, 1013.98] {
            let omega = FREQUENCY_WAVELENGTH_NM / wavelength;
            let n = water.complex_index(omega).unwrap().re;
            let ci = 0.8_f64;
            let q = (n * n - (1.0 - ci * ci)).sqrt();
            let kz = omega / VACUUM_SPEED_M_S * q;
            for delta in [
                0.0,
                core::f64::consts::FRAC_PI_2,
                core::f64::consts::PI,
                7.3,
            ] {
                let thickness = delta / kz;
                let stack = PlanarDielectricLayers::try_new(
                    &port,
                    DielectricResponse::Lossless(&port),
                    omega,
                    &[DielectricLayer {
                        response: DielectricResponse::WaterOptical(water),
                        thickness_m: thickness,
                    }],
                )
                .unwrap();
                for polarization in [Polarization::Te, Polarization::Tm] {
                    let yi = admittance(polarization, 1.0, ci, MU);
                    let yw = admittance(polarization, n, q, MU);
                    let r = (yi - yw) / (yi + yw);
                    let p = phase(-delta);
                    let denominator = C64::ONE - (p * p).scale(r * r);
                    let expected_r = (C64::ONE - p * p).scale(r) / denominator;
                    let expected_t = p.scale(1.0 - r * r) / denominator;
                    let result = stack.response(ci, polarization).unwrap();
                    complex_near(result.reflected_tangential_e, expected_r, 4e-15);
                    complex_near(result.transmitted_tangential_e, expected_t, 8e-15);
                    near(result.reflected_power_fraction, expected_r.norm_sq(), 2e-15);
                    near(
                        result.transmitted_power_fraction,
                        expected_t.norm_sq(),
                        3e-15,
                    );
                    near(result.absorbed_power_fraction, 0.0, 0.0);
                    near(
                        result.reflected_power_fraction + result.transmitted_power_fraction,
                        1.0,
                        4e-15,
                    );
                    let field = stack.prepare_layer_field(0, ci, polarization).unwrap();
                    let a = C64::from_re(1.0 + r) / denominator;
                    for z in [0.0, 0.5 * thickness, thickness] {
                        let forward = a * phase(-kz * z);
                        let backward = (a * p * p * phase(kz * z)).scale(-r);
                        let actual = field.field(z).unwrap();
                        complex_near(actual.electric_v_per_m, forward + backward, 1e-14);
                        complex_near(
                            actual.conjugate_magnetic_a_per_m,
                            (forward - backward).scale(yw),
                            4e-17,
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn coherent_water_beams_refract_with_exact_phase_and_vector_poynting_flux() {
    let port = unit_index_port();
    let source = water_state(273.15, 999.84, MU);
    let water = WaterOpticalResponse::try_lossless(&source, "T").unwrap();
    let modes = [
        (2_usize, C64::new(0.2, -0.1), C64::new(0.1, 0.3)),
        (3, C64::new(0.8, 0.2), C64::new(-0.5, 0.4)),
        (4, C64::new(-0.2, 0.2), C64::new(0.3, 0.1)),
    ];
    for wavelength in [226.5, 589.0, 1013.98] {
        let lambda = wavelength * 1e-9;
        let omega = FREQUENCY_WAVELENGTH_NM / wavelength;
        let k0 = omega / VACUUM_SPEED_M_S;
        let n = water.complex_index(omega).unwrap().re;
        let grid = AngularSpectrumGrid {
            width: 16,
            height: 1,
            spacing_m: [0.625 * lambda, lambda],
        };
        let input: Vec<_> = (0..grid.width)
            .map(|j| {
                let mut e = [C64::ZERO; 2];
                for &(bin, ex, ey) in &modes {
                    let p =
                        phase(-core::f64::consts::TAU * bin as f64 * j as f64 / grid.width as f64);
                    e[0] = e[0] + ex * p;
                    e[1] = e[1] + ey * p;
                }
                e
            })
            .collect();
        let stack = PlanarDielectricLayers::try_new(
            &port,
            DielectricResponse::WaterOptical(water),
            omega,
            &[],
        )
        .unwrap();
        let beam =
            PlanarAngularSpectrum::prepare(&stack, grid, &input, beam_limits(), || false).unwrap();
        let area = grid.width as f64 * grid.spacing_m[0] * grid.spacing_m[1];
        let mut incident_power = 0.0;
        let mut transmitted_power = 0.0;
        for &(bin, ex, ey) in &modes {
            let kx = core::f64::consts::TAU * bin as f64 / (grid.width as f64 * grid.spacing_m[0]);
            let qi = (1.0 - (kx / k0).powi(2)).sqrt();
            let qt = (n * n - (kx / k0).powi(2)).sqrt();
            for (p, amplitude) in [(Polarization::Tm, ex), (Polarization::Te, ey)] {
                let yi = admittance(p, 1.0, qi, MU);
                let yt = admittance(p, n, qt, MU);
                let t = 2.0 * yi / (yi + yt);
                incident_power += 0.5 * area * yi * amplitude.norm_sq();
                transmitted_power += 0.5 * area * yt * t * t * amplitude.norm_sq();
            }
        }
        near(
            beam.power().incident_w,
            incident_power,
            incident_power * 5e-15,
        );
        near(
            beam.power().transmitted_w,
            transmitted_power,
            transmitted_power * 5e-15,
        );
        near(beam.power().balance_defect_w(), 0.0, incident_power * 5e-15);
        for z in [0.0, 11.0 * lambda] {
            let field = beam
                .field(AngularSpectrumSide::Transmitted, z, || false)
                .unwrap();
            for j in 0..grid.width {
                let x = j as f64 * grid.spacing_m[0];
                let mut expected_e = [C64::ZERO; 3];
                let mut expected_h = [C64::ZERO; 3];
                for &(bin, ex, ey) in &modes {
                    let kx = core::f64::consts::TAU * bin as f64
                        / (grid.width as f64 * grid.spacing_m[0]);
                    let qi = (1.0 - (kx / k0).powi(2)).sqrt();
                    let qt = (n * n - (kx / k0).powi(2)).sqrt();
                    let te = 2.0 * qi / (qi + qt);
                    let tm = 2.0 * qt / (qt + n * n * qi);
                    let ex = ex.scale(tm);
                    let ey = ey.scale(te);
                    let kz = k0 * qt;
                    let p = phase(-kx * x - kz * z);
                    let e = [ex, ey, ex.scale(-kx / kz)];
                    let h = [
                        ey.scale(-kz / (omega * MU)),
                        ex.scale(n * n / (qt * VACUUM_SPEED_M_S * MU)),
                        ey.scale(kx / (omega * MU)),
                    ];
                    for c in 0..3 {
                        expected_e[c] = expected_e[c] + e[c] * p;
                        expected_h[c] = expected_h[c] + h[c] * p;
                    }
                }
                for c in 0..3 {
                    complex_near(field.electric_v_per_m()[j][c], expected_e[c], 7e-14);
                    complex_near(field.magnetic_a_per_m()[j][c], expected_h[c], 4e-16);
                }
            }
            near(
                field.normal_power_w(),
                transmitted_power,
                transmitted_power * 8e-15,
            );
        }
    }
}

#[test]
fn water_multifrequency_source_identity_is_retained_without_inventing_a_heat_reservoir() {
    let port = unit_index_port();
    let source = water_state(273.15, 999.84, MU);
    let warmer_source = water_state(373.15, 962.93, MU);
    let water = WaterOpticalResponse::try_lossless(&source, "T").unwrap();
    let warmer = WaterOpticalResponse::try_lossless(&warmer_source, "T").unwrap();
    let make_stack = |law, wavelength| {
        PlanarDielectricLayers::try_new(
            &port,
            DielectricResponse::Lossless(&port),
            FREQUENCY_WAVELENGTH_NM / wavelength,
            &[DielectricLayer {
                response: DielectricResponse::WaterOptical(law),
                thickness_m: 2e-6,
            }],
        )
        .unwrap()
    };
    let stacks = [make_stack(water, 500.0), make_stack(water, 650.0)];
    assert!(stacks[0].same_sources_and_geometry(&stacks[1]));
    let changed = make_stack(warmer, 650.0);
    assert!(!stacks[0].same_sources_and_geometry(&changed));
    let grid = AngularSpectrumGrid {
        width: 1,
        height: 1,
        spacing_m: [1e-6; 2],
    };
    let beams: Vec<_> = stacks
        .iter()
        .map(|stack| {
            PlanarAngularSpectrum::prepare(
                stack,
                grid,
                &[[C64::ONE, C64::ZERO]],
                beam_limits(),
                || false,
            )
            .unwrap()
        })
        .collect();
    let exposure = AngularSpectrumExposure::try_new(
        &beams[0],
        &beams,
        AngularSpectrumExposureLimits {
            maximum_lines: 2,
            maximum_point_line_pairs: 8,
            maximum_point_frequency_pairs: 8,
            maximum_workspace_bytes: 2_000_000,
            maximum_fft_work: 100_000,
        },
        || false,
    )
    .unwrap();
    assert_eq!(exposure.frequency_count(), 2);
    let result = exposure.prepare_layer(
        0,
        AngularSpectrumExposureLayerLimits {
            maximum_frequency_groups: 2,
            maximum_frequency_layer_mode_pairs: 8,
            maximum_point_line_pairs: 8,
            maximum_point_reservoir_frequency_pairs: 100,
            maximum_workspace_bytes: 2_000_000,
            maximum_fft_work: 100_000,
        },
        || false,
    );
    match result {
        Err(error) => assert!(
            error
                .to_string()
                .contains("original causal material reservoir"),
            "{error}"
        ),
        Ok(_) => {
            panic!("real phase dispersion cannot silently become an instantaneous energy store")
        }
    }
}
