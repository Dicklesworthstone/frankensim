//! Manufactured Kelvin--Voigt layer stresses, independent of section averages.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::section::strength::KelvinVoigtTensileCapacity;
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::{Density, DynViscosity, Pressure};

const LIMIT: &str = "tensile_limit";
const AXIAL: &str = "axial_eta";
const BENDING: &str = "bending_eta";

fn material(name: &str, values: [f64; 4], temperature: f64) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (property, value, dims) in [
        ("density", 8000.0, Density::DIMS),
        ("young_modulus", values[0], Pressure::DIMS),
        (AXIAL, values[1], DynViscosity::DIMS),
        (BENDING, values[2], DynViscosity::DIMS),
        (LIMIT, values[3], Pressure::DIMS),
    ] {
        let validity = ValidityDomain::unconstrained().with("T", 280.0, 400.0);
        let validity = match property {
            AXIAL => validity.with("omega", 20.0, 2000.0),
            BENDING => validity.with("omega", 10.0, 1000.0),
            _ => validity,
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(property, dims),
                value: PropertyValue::Scalar { value, dims },
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: format!("manufactured viscous tensile layer {name}"),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(property, dims, ScalarAdmissibility::Finite)
                .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "synthetic solid".into(),
            process: "manufactured section coefficients".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new()
        .with("T", temperature)
        .unwrap()
        .with("omega", 100.0)
        .unwrap();
    resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn section(layers: &[(f64, [f64; 4])], name: &str, temperature: f64) -> ConcentricElasticSection {
    let layers = layers
        .iter()
        .enumerate()
        .map(|(index, (radius, values))| ConcentricMaterialLayer {
            outer_radius_m: *radius,
            material: material(&format!("{name}-{index}"), *values, temperature),
        })
        .collect();
    ConcentricElasticSection::try_new(0.0, layers).unwrap()
}

fn prepared(section: &ConcentricElasticSection) -> KelvinVoigtTensileCapacity {
    let axial = section.clone().with_axial_viscosity(AXIAL).unwrap();
    let bending = section.clone().with_kelvin_voigt_bending(BENDING).unwrap();
    section
        .axial_tensile_capacity(LIMIT)
        .unwrap()
        .with_kelvin_voigt(&axial, &bending)
        .unwrap()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 4e-13 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn original_layer_coefficients_govern_total_stress_instead_of_averaged_loss() {
    let values = [
        (0.001, [200e9, 2e9, 200e9, 100e6]),
        (0.002, [80e9, 16e9, 8e9, 120e6]),
    ];
    let section = section(&values, "layered", 300.0);
    let axial = section.clone().with_axial_viscosity(AXIAL).unwrap();
    let bending = section.clone().with_kelvin_voigt_bending(BENDING).unwrap();
    let elastic = section.axial_tensile_capacity(LIMIT).unwrap();
    let capacity = elastic.clone().with_kelvin_voigt(&axial, &bending).unwrap();
    assert_eq!(capacity.elastic_capacity(), &elastic);
    assert_eq!(capacity.axial_viscosity_property(), AXIAL);
    assert_eq!(capacity.bending_viscosity_property(), BENDING);
    assert_eq!(
        capacity.axial_viscous_rigidity_n_s(),
        axial.viscous_rigidity_n_s()
    );
    assert_eq!(
        capacity.bending_viscous_stiffness_n_m2_s(),
        bending.viscous_stiffness_n_m2_s()
    );
    assert_eq!(capacity.omega_band_rad_s(), (20.0, 1000.0));

    let (strain, rate) = (1e-4, 0.01);
    let curvature = [0.2, -0.1];
    let curvature_rate = [-0.2, 0.4];
    let total = capacity
        .assess_strain_rates(strain, rate, curvature, curvature_rate)
        .unwrap();
    for (index, &(radius, [young, eta_axial, eta_bending, limit])) in values.iter().enumerate() {
        let coefficients = &capacity.layers()[index];
        assert_eq!(coefficients.young_modulus_pa(), young);
        assert_eq!(coefficients.axial_viscosity_pa_s(), eta_axial);
        assert_eq!(coefficients.bending_viscosity_pa_s(), eta_bending);
        // Independent dimensional stress formula, with distinct axial/bending eta/E.
        let mean = young * strain + eta_axial * rate;
        let b0 = young * curvature[0] + eta_bending * curvature_rate[0];
        let b1 = young * curvature[1] + eta_bending * curvature_rate[1];
        let amplitude = radius * (b0 * b0 + b1 * b1).sqrt();
        let observed = &total.layers()[index];
        near(observed.minimum_axial_stress_pa, mean - amplitude);
        near(observed.maximum_axial_stress_pa, mean + amplitude);
        near(
            observed.tensile_utilization,
            (mean + amplitude).max(0.0) / limit,
        );
        assert_eq!(
            *observed,
            capacity
                .assess_layer_strain_rates(index, strain, rate, curvature, curvature_rate)
                .unwrap()
        );
        assert_eq!(
            capacity.elastic_capacity().section().layers()[index].material,
            section.layers()[index].material
        );
    }
    assert_eq!(total.governing_layer(), 1);
    assert!(total.exceeds_limit());
    let elastic_only = elastic
        .assess_axial_bending(
            section.axial_stiffness_n() * strain,
            curvature.map(|kappa| kappa * section.bending_stiffness_n_m2()),
        )
        .unwrap();
    assert!(
        !elastic_only.exceeds_limit(),
        "viscous coating tension is the missing failure demand"
    );

    let rotated = capacity
        .assess_strain_rates(
            strain,
            rate,
            [-curvature[1], curvature[0]],
            [-curvature_rate[1], curvature_rate[0]],
        )
        .unwrap();
    assert_eq!(rotated, total);
}

#[test]
fn elastic_and_viscous_stress_cancel_before_the_tensile_envelope() {
    let capacity = prepared(&section(
        &[(0.001, [200e9, 2e9, 200e9, 100e6])],
        "cancel",
        300.0,
    ));
    let cancelled = capacity
        .assess_strain_rates(1e-4, -0.01, [0.3, -0.4], [-0.3, 0.4])
        .unwrap();
    assert_eq!(cancelled.layers()[0].minimum_axial_stress_pa, 0.0);
    assert_eq!(cancelled.layers()[0].maximum_axial_stress_pa, 0.0);
    assert_eq!(cancelled.utilization(), 0.0);
    let unloading = capacity
        .assess_strain_rates(0.0, -0.01, [0.0; 2], [0.0; 2])
        .unwrap();
    assert_eq!(unloading.layers()[0].maximum_axial_stress_pa, -20e6);
    assert_eq!(unloading.utilization(), 0.0);
}

#[test]
fn explicit_zero_viscosities_recover_original_elastic_stresses() {
    let section = section(
        &[
            (0.001, [200e9, 0.0, 0.0, 100e6]),
            (0.002, [80e9, 0.0, 0.0, 120e6]),
        ],
        "lossless",
        300.0,
    );
    let capacity = prepared(&section);
    let strain = -1e-4;
    let curvature = [0.12, -0.07];
    let actual = capacity
        .assess_strain_rates(strain, 23.0, curvature, [-19.0, 31.0])
        .unwrap();
    let original = capacity
        .elastic_capacity()
        .assess_axial_bending(
            strain * section.axial_stiffness_n(),
            curvature.map(|value| value * section.bending_stiffness_n_m2()),
        )
        .unwrap();
    for (a, b) in actual.layers().iter().zip(original.layers()) {
        near(a.minimum_axial_stress_pa, b.minimum_axial_stress_pa);
        near(a.maximum_axial_stress_pa, b.maximum_axial_stress_pa);
        near(a.tensile_utilization, b.tensile_utilization);
    }
    // An inactive dashpot must not require representing a needless r*kappa_dot.
    assert!(
        capacity
            .assess_strain_rates(0.0, 0.0, [0.0; 2], [f64::from_bits(1), 0.0])
            .is_ok()
    );
}

#[test]
fn splitting_one_material_into_annuli_preserves_the_viscous_tensile_limit() {
    let source = material("same-material", [200e9, 2e9, 3e9, 100e6], 300.0);
    let make = |radii: &[f64]| {
        ConcentricElasticSection::try_new(
            0.0,
            radii
                .iter()
                .map(|radius| ConcentricMaterialLayer {
                    outer_radius_m: *radius,
                    material: source.clone(),
                })
                .collect(),
        )
        .unwrap()
    };
    let whole = prepared(&make(&[0.002]));
    let split = prepared(&make(&[0.0005, 0.001, 0.002]));
    let a = whole
        .assess_strain_rates(1e-4, 0.01, [0.2, -0.3], [-0.4, 0.5])
        .unwrap();
    let b = split
        .assess_strain_rates(1e-4, 0.01, [0.2, -0.3], [-0.4, 0.5])
        .unwrap();
    assert_eq!(a.layers()[0], b.layers()[2]);
    assert_eq!(a.utilization(), b.utilization());
    assert_eq!(b.governing_layer(), 2);
}

#[test]
fn numerically_equal_foreign_or_stale_sources_do_not_bind() {
    let layers = [(0.001, [200e9, 2e9, 3e9, 100e6])];
    let original = section(&layers, "original", 300.0);
    let axial = original.clone().with_axial_viscosity(AXIAL).unwrap();
    let bending = original.clone().with_kelvin_voigt_bending(BENDING).unwrap();
    for other in [
        section(&layers, "foreign", 300.0),
        section(&layers, "original", 310.0),
    ] {
        let foreign_axial = other.clone().with_axial_viscosity(AXIAL).unwrap();
        let foreign_bending = other.with_kelvin_voigt_bending(BENDING).unwrap();
        let capacity = original.axial_tensile_capacity(LIMIT).unwrap();
        assert!(
            capacity
                .clone()
                .with_kelvin_voigt(&foreign_axial, &bending)
                .is_err()
        );
        assert!(
            capacity
                .with_kelvin_voigt(&axial, &foreign_bending)
                .is_err()
        );
    }
}

#[test]
fn invalid_stress_amplitudes_kinematics_and_unrepresentable_results_refuse() {
    let capacity = prepared(&section(
        &[(0.001, [200e9, 2e9, 3e9, 100e6])],
        "domain",
        300.0,
    ));
    let elastic = capacity.elastic_capacity();
    let compression = elastic.assess_layer_normal_stress(0, -20e6, 5e6).unwrap();
    assert_eq!(compression.maximum_axial_stress_pa, -15e6);
    assert_eq!(compression.tensile_utilization, 0.0);
    assert!(elastic.assess_layer_normal_stress(1, 0.0, 0.0).is_err());
    assert!(elastic.assess_layer_normal_stress(0, 0.0, -1.0).is_err());
    assert!(
        elastic
            .assess_layer_normal_stress(0, f64::MAX, f64::MAX)
            .is_err()
    );
    assert!(
        capacity
            .assess_layer_strain_rates(1, 0.0, 0.0, [0.0; 2], [0.0; 2])
            .is_err()
    );
    assert!(
        capacity
            .assess_strain_rates(f64::MAX, 0.0, [0.0; 2], [0.0; 2])
            .is_err()
    );
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(elastic.assess_layer_normal_stress(0, bad, 0.0).is_err());
        assert!(elastic.assess_layer_normal_stress(0, 0.0, bad).is_err());
        assert!(
            capacity
                .assess_strain_rates(bad, 0.0, [0.0; 2], [0.0; 2])
                .is_err()
        );
        assert!(
            capacity
                .assess_strain_rates(0.0, bad, [0.0; 2], [0.0; 2])
                .is_err()
        );
        assert!(
            capacity
                .assess_strain_rates(0.0, 0.0, [0.0, bad], [0.0; 2])
                .is_err()
        );
        assert!(
            capacity
                .assess_strain_rates(0.0, 0.0, [0.0; 2], [bad, 0.0])
                .is_err()
        );
    }
}
