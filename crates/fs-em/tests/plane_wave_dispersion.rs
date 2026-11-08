#![cfg(feature = "plane-wave")]
//! G1/G3 frequency derivatives against the original Maxwell solves.

use fs_em::plane_wave::{
    DielectricResponse, PlanarDielectricInterface, Polarization,
    angular_spectrum::{
        AngularSpectrumGrid, AngularSpectrumLimits, AngularSpectrumSide, PlanarAngularSpectrum,
    },
    layers::PlanarDielectricLayers,
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
    let angular = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Frequency(FrequencyConvention::Angular), ValueForm::Static,
    ));
    let conductivity = ElectricalPropertyBasis::Conductivity;
    let entries = [
        (EPS_NAME, QuantitySpec::dimensional(EPS), 2.0 * E),
        (MU_NAME, QuantitySpec::dimensional(MU), M),
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
