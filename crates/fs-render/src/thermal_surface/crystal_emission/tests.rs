use super::*;
use crate::conductor::{
    CoatingIorSample, CoatingLayer, CoherentAbsorbingCoating, ConductorDataStatus, ConductorSource,
};
use crate::dielectric::{DielectricIorSample, TabulatedIor};
use crate::thermal_surface::CoatingThermalEmission;
use fs_blake3::ContentHash;
use fs_em::plane_wave::{
    DielectricResponse, UniaxialOpticalResponse, uniaxial::oriented::OrientedUniaxialInterface,
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
        MaterialPropertySelection, ResolvedMaterialStatePoint, SPECTRAL_CONDUCTOR_K_QUANTITY,
        ScalarAdmissibility, ScalarPropertyRequirement, SpectralComplexIndexKnot,
        SpectralComplexIndexStatePoint, resolve_material_state_point,
    },
};
use fs_qty::Dims;

const DEPTH_M: f64 = 4.0e-6;
const LAMBDA_NM: f64 = 10_000.0;

fn temperature_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn resolved(properties: Vec<(PropertyKey, f64)>, temperature_k: f64) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (key, value) in properties {
        requirements.push(
            ScalarPropertyRequirement::try_with_quantity(
                key.name(),
                key.quantity(),
                ScalarAdmissibility::NonNegative,
            )
            .unwrap(),
        );
        let dims = key.dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    temperature_quantity(),
                    200.0,
                    2000.0,
                ),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored polarized thermal-emission fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic anisotropic thermal emitter".into(),
            phase: "solid".into(),
            process: "explicit two-axis optical claims".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity("T", temperature_quantity(), temperature_k)
        .unwrap();
    resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

struct Fixture {
    ordinary: SpectralComplexIndexStatePoint,
    extraordinary: SpectralComplexIndexStatePoint,
    exterior: ResolvedMaterialStatePoint,
}

impl Fixture {
    fn new(
        indices: [[f64; 2]; 2],
        temperature_k: f64,
        exterior_index: f64,
        exterior_mu: f64,
    ) -> Self {
        let bindings = |prefix: &str| {
            [0, 1].map(|i| SpectralComplexIndexKnot {
                wavelength_nm: [8000.0, 14_000.0][i],
                real_index_property: PropertyKey::new(format!("{prefix}_n_{i}"), Dims::NONE),
                extinction_index_property: PropertyKey::with_quantity(
                    format!("{prefix}_k_{i}"),
                    SPECTRAL_CONDUCTOR_K_QUANTITY,
                ),
            })
        };
        let ordinary = bindings("ordinary");
        let extraordinary = bindings("extraordinary");
        let mut properties = vec![(
            PropertyKey::new(MAGNETIC_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS),
            VACUUM_MU_H_M,
        )];
        for (axis, knots) in [&ordinary, &extraordinary].into_iter().enumerate() {
            for knot in knots {
                properties.push((knot.real_index_property.clone(), indices[axis][0]));
                properties.push((knot.extinction_index_property.clone(), indices[axis][1]));
            }
        }
        let material = resolved(properties, temperature_k);
        let ordinary =
            SpectralComplexIndexStatePoint::try_from_resolved(material.clone(), &ordinary).unwrap();
        let extraordinary =
            SpectralComplexIndexStatePoint::try_from_resolved(material, &extraordinary).unwrap();
        let exterior = resolved(
            vec![
                (
                    PropertyKey::new(ELECTRIC_PERMITTIVITY_PROPERTY, ELECTRIC_PERMITTIVITY_DIMS),
                    exterior_index * exterior_index / (exterior_mu * LIGHT_M_S * LIGHT_M_S),
                ),
                (
                    PropertyKey::new(MAGNETIC_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS),
                    exterior_mu,
                ),
            ],
            temperature_k,
        );
        Self {
            ordinary,
            extraordinary,
            exterior,
        }
    }

    fn interface(&self, axis: [f64; 3]) -> OrientedUniaxialInterface<'_> {
        OrientedUniaxialInterface::try_new(
            DielectricResponse::Lossless(&self.exterior),
            UniaxialOpticalResponse::try_new(&self.ordinary, &self.extraordinary, 6.1e-6).unwrap(),
            axis,
        )
        .unwrap()
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3.0e-10 * expected.abs().max(f64::MIN_POSITIVE),
        "actual={actual:e}, expected={expected:e}"
    );
}

fn planck_si(temperature_k: f64, wavelength_nm: f64) -> f64 {
    let wavelength = wavelength_nm * 1.0e-9;
    let h = 6.626_070_15e-34;
    let k = 1.380_649e-23;
    2.0 * h * LIGHT_M_S * LIGHT_M_S
        / (wavelength.powi(5) * (h * LIGHT_M_S / (wavelength * k * temperature_k)).exp_m1())
        * 1.0e-9
}

#[test]
fn original_temperature_and_complex_analyzers_match_reciprocal_maxwell_work() {
    let fixture = Fixture::new([[1.55, 0.08], [2.1, 0.25]], 800.0, 1.0, VACUUM_MU_H_M);
    let entrance = fixture.interface([1.0, 2.0, 3.0]);
    let slab = OrientedUniaxialSlab::try_new(
        &entrance,
        DielectricResponse::Lossless(&fixture.exterior),
        DEPTH_M,
    )
    .unwrap();
    let emitter = UniaxialSlabThermalEmission::try_from_material_temperature(&slab, "T").unwrap();
    assert!(core::ptr::eq(emitter.source(), &slab));
    assert_eq!(emitter.temperature_k(), 800.0);
    let c = 0.7_f64;
    let u = [0.6, 0.8];
    let response = emitter.prepare(LAMBDA_NM, c, [3.0, 4.0]).unwrap();
    assert!(core::ptr::eq(response.source(), &slab));
    assert_eq!(response.temperature_k(), 800.0);
    assert_eq!(response.wavelength_nm(), LAMBDA_NM);
    let stokes = response.stokes().unwrap();
    assert!(stokes.u.abs() > 1.0e-5 * stokes.i);
    close(response.radiance_w_m2_sr_nm().unwrap(), stokes.i);
    let omega = core::f64::consts::TAU * LIGHT_M_S * 1.0e9 / LAMBDA_NM;
    for analyzer in [
        [C64::ONE, C64::ZERO],
        [C64::ZERO, C64::ONE],
        [C64::ONE, C64::new(0.0, 1.0)],
        [C64::new(0.3, -0.2), C64::new(-0.1, 0.7)],
    ] {
        // Independent single-illumination solve, using a signed conserved
        // tangential index rather than the prepared polarization cache.
        let norm = (analyzer[0].norm_sq() + analyzer[1].norm_sq()).sqrt();
        let a = analyzer.map(|v| v.scale(1.0 / norm));
        let electric = [
            (a[0].scale(-u[1]) + a[1].scale(c * u[0])).conj(),
            (a[0].scale(u[0]) + a[1].scale(c * u[1])).conj(),
        ];
        let field = slab
            .prepare_absorbing_field_at_tangential_index_in_plane(
                omega,
                -(1.0 - c * c).sqrt(),
                u,
                electric,
            )
            .unwrap();
        let expected = 0.5
            * planck_si(800.0, LAMBDA_NM)
            * field
                .slab_absorption(1.0, 0.0, DEPTH_M)
                .unwrap()
                .absorbed_power_w()
            / field.response().incident_normal_power_w_m2;
        let actual = response.analyzer_radiance_w_m2_sr_nm(analyzer).unwrap();
        close(actual, expected);
        let cross = a[0] * a[1].conj();
        close(
            actual,
            0.5 * (stokes.i
                + stokes.q * (a[0].norm_sq() - a[1].norm_sq())
                + 2.0 * (stokes.u * cross.re + stokes.v * cross.im)),
        );
        for scale in [1.0e-250, 1.0e250] {
            let phase = C64::new(0.6, 0.8).scale(scale);
            close(
                response
                    .analyzer_radiance_w_m2_sr_nm(analyzer.map(|a| a * phase))
                    .unwrap(),
                actual,
            );
        }
    }
}

#[test]
fn isotropic_limit_matches_existing_coating_emission_and_lossless_film_is_dark() {
    for extinction in [0.0, 0.12] {
        let fixture = Fixture::new([[1.7, extinction]; 2], 800.0, 1.0, VACUUM_MU_H_M);
        let entrance = fixture.interface([1.0, 2.0, 3.0]);
        let slab = OrientedUniaxialSlab::try_new(
            &entrance,
            DielectricResponse::Lossless(&fixture.exterior),
            DEPTH_M,
        )
        .unwrap();
        let emitter =
            UniaxialSlabThermalEmission::try_from_material_temperature(&slab, "T").unwrap();
        let film = CoherentAbsorbingCoating::try_new_in_band(
            core::array::from_fn(|i| {
                CoatingIorSample::try_new_in_band(8000.0 + 750.0 * i as f64, 1.7, extinction)
                    .unwrap()
            }),
            DEPTH_M,
            750.0,
            ConductorSource::try_new(ContentHash([91; 32]), ConductorDataStatus::Representative)
                .unwrap(),
        )
        .unwrap();
        let vacuum = TabulatedIor::try_new_in_band(
            core::array::from_fn(|i| {
                DielectricIorSample::try_new_in_band(8000.0 + 750.0 * i as f64, 1.0).unwrap()
            }),
            750.0,
        )
        .unwrap()
        .into();
        let scalar = CoatingThermalEmission::try_transmitting(
            &[CoatingLayer::Absorbing(&film)],
            &[800.0],
            vacuum,
            false,
            ContentHash([92; 32]),
        )
        .unwrap();
        for cosine in [1.0, 0.63] {
            let response = emitter.prepare(LAMBDA_NM, cosine, [3.0, 4.0]).unwrap();
            let expected = scalar
                .evaluate(LAMBDA_NM, cosine)
                .unwrap()
                .polarized_radiance_w_m2_sr_nm();
            for (analyzer, expected) in [[C64::ONE, C64::ZERO], [C64::ZERO, C64::ONE]]
                .into_iter()
                .zip(expected)
            {
                close(
                    response.analyzer_radiance_w_m2_sr_nm(analyzer).unwrap(),
                    expected,
                );
            }
            if extinction == 0.0 {
                assert_eq!(response.stokes().unwrap(), Stokes::unpolarized(0.0));
            } else {
                let s = response.stokes().unwrap();
                assert!(s.u.abs() < 1.0e-12 * s.i && s.v.abs() < 1.0e-12 * s.i);
            }
        }
    }
}

#[test]
fn positive_weak_polarized_emission_survives_a_rounded_stokes_channel() {
    let fixture = Fixture::new([[1.45, 0.12], [2.1, 1.0e-24]], 800.0, 1.0, VACUUM_MU_H_M);
    let entrance = fixture.interface([1.0, 0.0, 0.0]);
    let slab = OrientedUniaxialSlab::try_new(
        &entrance,
        DielectricResponse::Lossless(&fixture.exterior),
        DEPTH_M,
    )
    .unwrap();
    let response = UniaxialSlabThermalEmission::try_from_material_temperature(&slab, "T")
        .unwrap()
        .prepare(LAMBDA_NM, 1.0, [1.0, 0.0])
        .unwrap();
    let strong = response
        .analyzer_radiance_w_m2_sr_nm([C64::ONE, C64::ZERO])
        .unwrap();
    let weak = response
        .analyzer_radiance_w_m2_sr_nm([C64::ZERO, C64::ONE])
        .unwrap();
    assert!(weak > 0.0 && weak < 1.0e-20 * strong);
    let stokes = response.stokes().unwrap();
    assert_eq!(stokes.i - stokes.q, 0.0);
    let omega = core::f64::consts::TAU * LIGHT_M_S * 1.0e9 / LAMBDA_NM;
    let field = slab
        .prepare_absorbing_field(omega, 1.0, [C64::ONE, C64::ZERO])
        .unwrap();
    close(
        weak,
        0.5 * planck_si(800.0, LAMBDA_NM)
            * field
                .slab_absorption(1.0, 0.0, DEPTH_M)
                .unwrap()
                .absorbed_power_w()
            / field.response().incident_normal_power_w_m2,
    );
}

#[test]
fn rotating_crystal_ray_and_analyzer_frame_preserves_emitted_radiance() {
    let fixture = Fixture::new([[1.55, 0.08], [2.1, 0.25]], 800.0, 1.0, VACUUM_MU_H_M);
    let angle = 0.83_f64;
    let (s, c) = angle.sin_cos();
    let rotate = |v: [f64; 3]| [c * v[0] - s * v[1], s * v[0] + c * v[1], v[2]];
    let first = fixture.interface([1.0, 2.0, 3.0]);
    let second = fixture.interface(rotate([1.0, 2.0, 3.0]));
    let first_slab = OrientedUniaxialSlab::try_new(
        &first,
        DielectricResponse::Lossless(&fixture.exterior),
        DEPTH_M,
    )
    .unwrap();
    let second_slab = OrientedUniaxialSlab::try_new(
        &second,
        DielectricResponse::Lossless(&fixture.exterior),
        DEPTH_M,
    )
    .unwrap();
    let plane = rotate([0.6, 0.8, 0.0]);
    let first = UniaxialSlabThermalEmission::try_from_material_temperature(&first_slab, "T")
        .unwrap()
        .prepare(LAMBDA_NM, 0.7, [0.6, 0.8])
        .unwrap();
    let second = UniaxialSlabThermalEmission::try_from_material_temperature(&second_slab, "T")
        .unwrap()
        .prepare(LAMBDA_NM, 0.7, [plane[0], plane[1]])
        .unwrap();
    for analyzer in [
        [C64::ONE, C64::ZERO],
        [C64::ZERO, C64::ONE],
        [C64::ONE, C64::new(0.2, 0.7)],
    ] {
        close(
            first.analyzer_radiance_w_m2_sr_nm(analyzer).unwrap(),
            second.analyzer_radiance_w_m2_sr_nm(analyzer).unwrap(),
        );
    }
    for (expected, actual) in rotate(first.outgoing_direction())
        .into_iter()
        .zip(second.outgoing_direction())
    {
        assert!((actual - expected).abs() < 2.0e-15);
    }
}

#[test]
fn unsupported_state_direction_spectrum_analyzer_and_observation_medium_refuse() {
    let fixture = Fixture::new([[1.55, 0.08], [2.1, 0.25]], 800.0, 1.0, VACUUM_MU_H_M);
    let entrance = fixture.interface([1.0, 2.0, 3.0]);
    let slab = OrientedUniaxialSlab::try_new(
        &entrance,
        DielectricResponse::Lossless(&fixture.exterior),
        DEPTH_M,
    )
    .unwrap();
    assert!(matches!(
        UniaxialSlabThermalEmission::try_from_material_temperature(&slab, "missing"),
        Err(CrystalEmissionError::TemperatureAxis)
    ));
    for t in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(UniaxialSlabThermalEmission::try_isothermal(&slab, t).is_err());
    }
    let emitter = UniaxialSlabThermalEmission::try_from_material_temperature(&slab, "T").unwrap();
    for cosine in [0.0, -0.5, 1.01, f64::NAN] {
        assert!(emitter.prepare(LAMBDA_NM, cosine, [1.0, 0.0]).is_err());
    }
    for plane in [[0.0, 0.0], [f64::INFINITY, 0.0]] {
        assert!(emitter.prepare(LAMBDA_NM, 0.7, plane).is_err());
    }
    for lambda in [0.0, f64::NAN, 7000.0, 15_000.0] {
        assert!(emitter.prepare(lambda, 0.7, [1.0, 0.0]).is_err());
    }
    let response = emitter.prepare(LAMBDA_NM, 0.7, [1.0, 0.0]).unwrap();
    for analyzer in [[C64::ZERO; 2], [C64::new(f64::NAN, 0.0), C64::ONE]] {
        assert!(matches!(
            response.analyzer_radiance_w_m2_sr_nm(analyzer),
            Err(CrystalEmissionError::Analyzer)
        ));
    }
    // n=1 alone does not identify the vacuum impedance of a magnetic medium.
    for (index, mu) in [(1.2, VACUUM_MU_H_M), (1.0, 2.0 * VACUUM_MU_H_M)] {
        let other = Fixture::new([[1.55, 0.08], [2.1, 0.25]], 800.0, index, mu);
        let entrance = other.interface([1.0, 2.0, 3.0]);
        let slab = OrientedUniaxialSlab::try_new(
            &entrance,
            DielectricResponse::Lossless(&other.exterior),
            DEPTH_M,
        )
        .unwrap();
        let emitter =
            UniaxialSlabThermalEmission::try_from_material_temperature(&slab, "T").unwrap();
        assert!(matches!(
            emitter.prepare(LAMBDA_NM, 0.7, [1.0, 0.0]),
            Err(CrystalEmissionError::ObservationMedium)
        ));
    }
}
