#![cfg(feature = "plane-wave")]
//! G0/G1/G3 native coherent-film photoemission with explicitly synthetic sources.

use fs_em::plane_wave::{
    DielectricResponse, Polarization,
    layers::{
        DielectricLayer, IlluminationBudget, PlanarDielectricLayers,
        photoelectric::{EscapeFace, LayerPhotoelectricError, PreparedLayerPhotoelectric},
    },
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
    },
    photoelectric::source::{PhotoelectricKnot, PhotoelectricMaterial},
    state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point,
    },
};
use fs_math::c64::C64;
use fs_qty::{Dims, Length};

const EPSILON: f64 = 8.854_187_812_8e-12;
const PERMEABILITY: f64 = 1.256_637_062_12e-6;
const ENERGY: Dims = Dims([2, 1, -2, 0, 0, 0]);
const WORK: f64 = 3.0e-19;
const THICKNESS: f64 = 12.0e-9;
const YIELD: f64 = 0.2;
const MEAN_FRACTION: f64 = 0.4;
const CHARGE: f64 = 1.602_176_634e-19;

fn source(properties: &[(&str, Dims, f64)], revision: u32) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for &(name, dims, value) in properties {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: Vec::new(),
                provenance: Provenance {
                    source: "synthetic photoemitting-film fixture, no measured-metal claim".into(),
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
            chemistry: "synthetic optical/photoemission fixture".into(),
            phase: "solid".into(),
            process: "manufactured constant law".into(),
            revision,
        },
        claims,
        Vec::new(),
    )
    .unwrap();
    resolve_material_state_point(
        &card,
        &QueryPoint::new(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn vacuum() -> ResolvedMaterialStatePoint {
    source(&[(EPS_NAME, EPS, EPSILON), (MU_NAME, MU, PERMEABILITY)], 0)
}

fn film(revision: u32) -> (ResolvedConductorState, PhotoelectricMaterial) {
    let basis = ElectricalPropertyBasis::Conductivity;
    let original = source(
        &[
            (EPS_NAME, EPS, 2.0 * EPSILON),
            (MU_NAME, MU, PERMEABILITY),
            (basis.property_name(), basis.dimensions(), 1.0e5),
            ("removal_work", ENERGY, WORK),
            ("emission_film_thickness", Length::DIMS, THICKNESS),
            ("internal_yield_short", Dims::NONE, YIELD),
            ("mean_excess_short", Dims::NONE, MEAN_FRACTION),
            ("internal_yield_long", Dims::NONE, YIELD),
            ("mean_excess_long", Dims::NONE, MEAN_FRACTION),
        ],
        revision,
    );
    let conductor =
        ResolvedConductorState::try_from_material_with_basis(original.clone(), basis).unwrap();
    let knots = [
        (200.0, "internal_yield_short", "mean_excess_short"),
        (2000.0, "internal_yield_long", "mean_excess_long"),
    ]
    .map(|(wavelength_nm, yield_name, mean_name)| PhotoelectricKnot {
        wavelength_nm,
        absorbed_quantum_yield_property: PropertyKey::new(yield_name, Dims::NONE),
        mean_excess_fraction_property: PropertyKey::new(mean_name, Dims::NONE),
    });
    let photoelectric = PhotoelectricMaterial::try_from_resolved(
        original,
        PropertyKey::new("removal_work", ENERGY),
        PropertyKey::new("emission_film_thickness", Length::DIMS),
        &knots,
        1800.0,
    )
    .unwrap();
    (conductor, photoelectric)
}

fn omega(wavelength_m: f64) -> f64 {
    core::f64::consts::TAU * 299_792_458.0 / wavelength_m
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs().max(1e-300),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_native_absorption_drives_emitted_charge_and_all_energy_channels_once() {
    let exterior = vacuum();
    let (metal, source) = film(0);
    let wavelength = 250e-9;
    let stack = PlanarDielectricLayers::try_new(
        &exterior,
        DielectricResponse::Lossless(&exterior),
        omega(wavelength),
        &[DielectricLayer {
            response: DielectricResponse::Ohmic(&metal),
            thickness_m: THICKNESS,
        }],
    )
    .unwrap();
    for polarization in [Polarization::Te, Polarization::Tm] {
        let illumination = stack.prepare_illumination(0.73, polarization).unwrap();
        let area = 0.001;
        let owner = PreparedLayerPhotoelectric::try_new(
            &illumination,
            0,
            &source,
            EscapeFace::EntranceVacuum,
            area,
        )
        .unwrap();
        let inputs = [C64::new(100.0, -30.0), C64::new(-25.0, 15.0)];
        let response = owner.evaluate(inputs).unwrap();
        let emission = response.photoelectric();
        let mut native_layer_power = [0.0];
        let native = illumination
            .evaluate(inputs, &mut native_layer_power)
            .unwrap();
        let power_from_ports = (native.power.incoming_w_m2.iter().sum::<f64>()
            - native.power.outgoing_w_m2.iter().sum::<f64>())
            * area;
        let photon_energy = 6.626_070_15e-34 * 299_792_458.0 / wavelength;
        let rate = power_from_ports * YIELD / photon_energy;
        let mean_kinetic = MEAN_FRACTION * (photon_energy - WORK);
        near(emission.absorbed_power_w, power_from_ports);
        near(emission.absorbed_power_w, native_layer_power[0] * area);
        near(emission.emitted_electron_rate_s, rate);
        near(emission.surface_charging_current_a, rate * CHARGE);
        assert_eq!(
            emission.surface_charging_current_a,
            -emission.emitted_electron_current_a
        );
        near(emission.maximum_kinetic_energy_j, photon_energy - WORK);
        near(emission.mean_kinetic_energy_j, mean_kinetic);
        near(emission.kinetic_power_w, rate * mean_kinetic);
        near(emission.removal_power_w, rate * WORK);
        near(
            emission.heat_power_w,
            power_from_ports - rate * (WORK + mean_kinetic),
        );
        near(
            emission.heat_power_w + emission.kinetic_power_w + emission.removal_power_w,
            emission.absorbed_power_w,
        );
        assert!(std::ptr::eq(response.material(), &source));
        assert!(std::ptr::eq(response.illumination(), &illumination));
        assert_eq!(response.layer_index(), 0);
        assert_eq!(response.escape_face(), EscapeFace::EntranceVacuum);
        near(response.vacuum_wavelength_m(), wavelength);
    }
}

#[test]
fn g3_coherent_phase_changes_photoemission_but_independent_sources_add() {
    let exterior = vacuum();
    let (metal, source) = film(0);
    let stack = PlanarDielectricLayers::try_new(
        &exterior,
        DielectricResponse::Lossless(&exterior),
        omega(250e-9),
        &[DielectricLayer {
            response: DielectricResponse::Ohmic(&metal),
            thickness_m: THICKNESS,
        }],
    )
    .unwrap();
    let illumination = stack.prepare_illumination(1.0, Polarization::Te).unwrap();
    let a =
        PreparedLayerPhotoelectric::try_new(&illumination, 0, &source, EscapeFace::ExitVacuum, 1.0)
            .unwrap();
    let b =
        PreparedLayerPhotoelectric::try_new(&illumination, 0, &source, EscapeFace::ExitVacuum, 2.0)
            .unwrap();
    let v = C64::from_re(100.0);
    let positive = a.evaluate([v, v]).unwrap();
    let negative = a.evaluate([v, -v]).unwrap();
    assert!(
        positive.photoelectric().emitted_electron_rate_s
            > 2.0 * negative.photoelectric().emitted_electron_rate_s
    );
    let independent = a.evaluate_independent_sources([v, v]).unwrap();
    let factors = [[v, C64::ZERO], [C64::ZERO, v]];
    let budget = IlluminationBudget {
        maximum_modes: 2,
        maximum_layer_mode_pairs: 2,
    };
    let factored = a
        .evaluate_coherence_modes(&factors, budget, || false)
        .unwrap();
    near(
        factored.photoelectric().emitted_electron_rate_s,
        independent.photoelectric().emitted_electron_rate_s,
    );
    let coherent = a
        .evaluate_coherence_modes(&[[v, v]], budget, || false)
        .unwrap();
    near(
        coherent.photoelectric().emitted_electron_rate_s,
        positive.photoelectric().emitted_electron_rate_s,
    );
    assert!(
        a.evaluate_coherence_modes(&factors, budget, || true)
            .is_err()
    );
    assert!(
        a.evaluate_coherence_modes(
            &factors,
            IlluminationBudget {
                maximum_modes: 1,
                maximum_layer_mode_pairs: 1,
            },
            || false
        )
        .is_err()
    );
    let left = a.evaluate([v, C64::ZERO]).unwrap();
    let right = a.evaluate([C64::ZERO, v]).unwrap();
    near(
        independent.photoelectric().emitted_electron_rate_s,
        left.photoelectric().emitted_electron_rate_s
            + right.photoelectric().emitted_electron_rate_s,
    );
    let scaled = b.evaluate([v.scale(2.0), v.scale(2.0)]).unwrap();
    near(
        scaled.photoelectric().emitted_electron_rate_s,
        8.0 * positive.photoelectric().emitted_electron_rate_s,
    );
    near(
        scaled.photoelectric().heat_power_w,
        8.0 * positive.photoelectric().heat_power_w,
    );
    assert_eq!(
        scaled.photoelectric().maximum_kinetic_energy_j,
        positive.photoelectric().maximum_kinetic_energy_j
    );
    let dark = a.evaluate([C64::ZERO; 2]).unwrap();
    assert_eq!(dark.photoelectric().absorbed_power_w, 0.0);
    assert_eq!(dark.photoelectric().emitted_electron_rate_s, 0.0);
}

#[test]
fn g1_subthreshold_absorbed_light_heats_without_emitting_electrons() {
    let exterior = vacuum();
    let (metal, source) = film(0);
    let stack = PlanarDielectricLayers::try_new(
        &exterior,
        DielectricResponse::Lossless(&exterior),
        omega(1000e-9),
        &[DielectricLayer {
            response: DielectricResponse::Ohmic(&metal),
            thickness_m: THICKNESS,
        }],
    )
    .unwrap();
    let illumination = stack.prepare_illumination(1.0, Polarization::Tm).unwrap();
    let owner = PreparedLayerPhotoelectric::try_new(
        &illumination,
        0,
        &source,
        EscapeFace::EntranceVacuum,
        1.0,
    )
    .unwrap();
    let result = owner.evaluate([C64::from_re(1e6), C64::ZERO]).unwrap();
    assert!(result.photoelectric().absorbed_power_w > 0.0);
    assert_eq!(result.photoelectric().emitted_electron_rate_s, 0.0);
    assert_eq!(
        result.photoelectric().heat_power_w,
        result.photoelectric().absorbed_power_w
    );
}

#[test]
fn g0_wrong_sources_thicknesses_buried_faces_and_invalid_inputs_refuse() {
    let exterior = vacuum();
    let (metal, source) = film(0);
    let (_, changed_source) = film(1);
    let layer = DielectricLayer {
        response: DielectricResponse::Ohmic(&metal),
        thickness_m: THICKNESS,
    };
    let stack = PlanarDielectricLayers::try_new(
        &exterior,
        DielectricResponse::Lossless(&exterior),
        omega(250e-9),
        &[layer; 3],
    )
    .unwrap();
    let illumination = stack.prepare_illumination(1.0, Polarization::Te).unwrap();
    for face in [EscapeFace::EntranceVacuum, EscapeFace::ExitVacuum] {
        assert!(matches!(
            PreparedLayerPhotoelectric::try_new(&illumination, 1, &source, face, 1.0,),
            Err(LayerPhotoelectricError::Invalid(_))
        ));
    }
    assert!(
        PreparedLayerPhotoelectric::try_new(
            &illumination,
            0,
            &changed_source,
            EscapeFace::EntranceVacuum,
            1.0,
        )
        .is_err()
    );
    for area in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            PreparedLayerPhotoelectric::try_new(
                &illumination,
                0,
                &source,
                EscapeFace::EntranceVacuum,
                area,
            )
            .is_err()
        );
    }
    let owner = PreparedLayerPhotoelectric::try_new(
        &illumination,
        0,
        &source,
        EscapeFace::EntranceVacuum,
        1.0,
    )
    .unwrap();
    assert!(
        owner
            .evaluate([C64::new(f64::NAN, 0.0), C64::ZERO])
            .is_err()
    );
    let wrong_thickness = PlanarDielectricLayers::try_new(
        &exterior,
        DielectricResponse::Lossless(&exterior),
        omega(250e-9),
        &[DielectricLayer {
            thickness_m: 2.0 * THICKNESS,
            ..layer
        }],
    )
    .unwrap();
    let wrong_illumination = wrong_thickness
        .prepare_illumination(1.0, Polarization::Te)
        .unwrap();
    assert!(
        PreparedLayerPhotoelectric::try_new(
            &wrong_illumination,
            0,
            &source,
            EscapeFace::EntranceVacuum,
            1.0,
        )
        .is_err()
    );
}

#[path = "plane_wave_photoelectric/depth.rs"]
mod depth;
