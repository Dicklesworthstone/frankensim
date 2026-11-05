//! Synthetic G0/G1/G3 source and axis tests; no measured-material claim.
use super::*;
use crate::shear_shell::{
    ResolvedShearShell, ShellRegionMaterial, TRANSVERSE_SHEAR_VISCOSITY_PROPERTIES,
    compile_material_shear_shell,
};
use crate::thin_plate::{PlateRegionMaterial, PlateThicknessConstraint};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::plane_stress_viscosity::PLANE_STRESS_VISCOSITY_PROPERTIES;
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_plate::{PlateRegion, shell::{ShellMesh, mitc::{MitcBudget, ShearShellSection}}};
use fs_qty::{Dims, DynViscosity, Pressure};

fn query() -> QueryPoint {
    QueryPoint::new().with_quantity("T", QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature, ValueForm::Static,
    )), 300.0).unwrap().with("omega", 1000.0).unwrap()
}

/// A typed source on a 3-by-3 flat node grid, with one complete region.
/// Parent runtime tests may reduce this same fixture with their own constraints.
pub(in crate::shear_shell::thermoelastic) fn fixture(
    model: PlateMaterialModel,
) -> (ResolvedShearShell, ShellThermalProperties) {
    fixture_with(model, query(), |_| {})
}

pub(in crate::shear_shell::thermoelastic) fn fixture_with(
    model: PlateMaterialModel,
    point: QueryPoint,
    modify: impl FnOnce(&mut Vec<PropertyClaim>),
) -> (ResolvedShearShell, ShellThermalProperties) {
    let mut values = vec![
        ("density", Density::DIMS, 2700.0),
        ("young_modulus", Pressure::DIMS, 70e9),
        ("poisson_ratio", Dims::NONE, 0.3),
        ("young_modulus_1", Pressure::DIMS, 70e9),
        ("young_modulus_2", Pressure::DIMS, 40e9),
        ("poisson_ratio_12", Dims::NONE, 0.23),
        ("shear_modulus_12", Pressure::DIMS, 18e9),
        ("shear_modulus_13", Pressure::DIMS, 16e9),
        ("shear_modulus_23", Pressure::DIMS, 12e9),
        ("alpha_1", INVERSE_TEMPERATURE_DIMS, 2e-5),
        ("alpha_2", INVERSE_TEMPERATURE_DIMS, -8e-6),
        ("alpha_12", INVERSE_TEMPERATURE_DIMS, 3e-6),
        ("cp", SPECIFIC_HEAT_CAPACITY_DIMS, 900.0),
        ("k_normal", THERMAL_CONDUCTIVITY_DIMS, 100.0),
    ];
    values.extend(PLANE_STRESS_VISCOSITY_PROPERTIES.map(|key| (key, DynViscosity::DIMS, 0.0)));
    values.extend(TRANSVERSE_SHEAR_VISCOSITY_PROPERTIES.map(|key| (key, DynViscosity::DIMS, 0.0)));
    let mut claims: Vec<_> = values.into_iter().map(|(key, dims, value)| {
        let band = if key == "k_normal" { (200.0, 1e6) }
            else if key == "cp" || key.starts_with("alpha_") { (100.0, 2e6) }
            else { (50.0, 1e7) };
        PropertyClaim {
            key: PropertyKey::new(key, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("omega", band.0, band.1),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance { source: "synthetic shell thermal fixture".into(),
                license: "CC0-1.0".into(), artifact: None },
        }
    }).collect();
    modify(&mut claims);
    let requirements: Vec<_> = claims.iter().map(|p|
        ScalarPropertyRequirement::try_with_key(&p.key, ScalarAdmissibility::Finite).unwrap()).collect();
    let mut set = ClaimSet::new();
    for claim in claims { set.insert_claim(claim).unwrap(); }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic shell thermal source".into(), phase: "solid".into(),
        process: "numerical fixture".into(), revision: 0,
    }, set, vec![]).unwrap();
    let material = resolve_material_state_point(&card, &point, &requirements,
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    let mut nodes = Vec::new();
    let mut tris = Vec::new();
    for j in 0..3 { for i in 0..3 { nodes.push([0.1 * i as f64, 0.075 * j as f64, 0.0]); } }
    for j in 0..2 { for i in 0..2 {
        let a = 3*j+i;
        tris.extend([[a, a+1, a+4], [a, a+4, a+3]]);
    } }
    let source = compile_material_shear_shell(ShellMesh::new(nodes, tris).unwrap(),
        vec![ShellRegionMaterial {
            region: PlateRegionMaterial {
                region: PlateRegion { name: "source".into(), triangle_indices: (0..8).collect() },
                material, model, thickness_constraint: PlateThicknessConstraint::FixedThickness(0.002),
                material_angle_rad: 0.17,
            },
            shear_correction: 5.0/6.0, axis1_global: [1.0, 0.0, 0.0],
        }], MitcBudget { max_nodes: 9, max_triangles: 8 }, &CancelGate::new_clock_free()).unwrap();
    let expansion = if model.is_isotropic() {
        ShellThermalExpansionProperties::Isotropic { linear_per_k: "alpha_1".into() }
    } else {
        ShellThermalExpansionProperties::Engineering {
            normal_1_per_k: "alpha_1".into(), normal_2_per_k: "alpha_2".into(),
            shear_12_per_k: "alpha_12".into(),
        }
    };
    (source, ShellThermalProperties { expansion, specific_heat: "cp".into(),
        normal_conductivity: "k_normal".into() })
}

#[test]
fn g1_original_source_sets_coupled_capacity_and_intersected_band() {
    let (source, properties) = fixture(PlateMaterialModel::Orthotropic12);
    let region = &source.regions()[0];
    let retained = region.input().region.material.clone();
    let bound = bind_region(region, &properties, (300.0, 3e6)).unwrap();
    assert_eq!(bound.source_band_rad_s, (300.0, 1e6));
    let h = region.thickness_m();
    let q = region.section().principal_matrices().0.coefficients().map(|a| a/h);
    let alpha = [2e-5, -8e-6, 3e-6];
    let beta: [f64; 3] = core::array::from_fn(|i|
        (0..3).map(|j| q[3*i+j]*alpha[j]).sum());
    let capacity = 2700.0*900.0 - 300.0*(0..3).map(|i| alpha[i]*beta[i]).sum::<f64>();
    let tau = h*h*capacity/(core::f64::consts::PI.powi(2)*100.0);
    assert!((bound.law.first_relaxation_time_s()/tau-1.0).abs() < 2e-14);
    let scale = (300.0*h*h*h/(12.0*capacity)).sqrt();
    for (actual, expected) in bound.law.relaxation_factor().into_iter().zip(beta.map(|b| scale*b)) {
        assert!((actual-expected).abs() <= 3e-14*expected.abs());
    }
    assert_eq!(region.input().region.material, retained, "source and receipts cannot change");
}

#[test]
fn g3_axis_exchange_preserves_shear_and_signed_thermal_channel() {
    let mut factors = Vec::new();
    for model in [PlateMaterialModel::Orthotropic12, PlateMaterialModel::Orthotropic21] {
        let (source, properties) = fixture(model);
        factors.push(bind_region(&source.regions()[0], &properties, (300.0, 1e6))
            .unwrap().law.relaxation_factor());
    }
    for (actual, expected) in factors[1].into_iter().zip([factors[0][1], factors[0][0], factors[0][2]]) {
        assert!((actual-expected).abs() <= 3e-14*expected.abs());
    }
    assert!(factors[0][1] < 0.0, "negative normal thermal coupling must survive");
    assert!(factors[0][2] > 0.0, "Axes21 must not reverse in-plane engineering shear");
    let (zero, properties) = fixture_with(PlateMaterialModel::Orthotropic12, query(), |claims| {
        for p in claims.iter_mut().filter(|p| p.key.name().starts_with("alpha_")) {
            p.value = PropertyValue::Scalar { value: 0.0, dims: INVERSE_TEMPERATURE_DIMS };
        }
    });
    assert_eq!(bind_region(&zero.regions()[0], &properties, (300.0, 1e6))
        .unwrap().law.relaxation_factor(), [0.0; 3]);
}

#[test]
fn g0_missing_units_unheld_curves_temperature_and_model_refuse() {
    let (source, mut properties) = fixture(PlateMaterialModel::Orthotropic12);
    properties.normal_conductivity = "missing".into();
    assert!(bind_region(&source.regions()[0], &properties, (300.0, 1e6)).is_err());
    properties.normal_conductivity = "k_normal".into();
    properties.expansion = ShellThermalExpansionProperties::Isotropic { linear_per_k: "alpha_1".into() };
    assert!(bind_region(&source.regions()[0], &properties, (300.0, 1e6)).is_err());
    for point in [QueryPoint::new().with("omega", 1000.0).unwrap(),
        QueryPoint::new().with("T", 300.0).unwrap().with("omega", 1000.0).unwrap(),
        QueryPoint::new().with_quantity("T", QuantitySpec::dimensional(fs_qty::Temperature::DIMS),
            300.0).unwrap().with("omega", 1000.0).unwrap()] {
        let (source, properties) = fixture_with(PlateMaterialModel::Isotropic, point, |_| {});
        assert!(bind_region(&source.regions()[0], &properties, (300.0, 1e6)).is_err());
    }
    for case in 0..3 {
        let (source, properties) = fixture_with(PlateMaterialModel::Isotropic, query(), |claims| {
            let p = claims.iter_mut().find(|p| p.key.name() == "cp").unwrap();
            match case {
                0 => {
                    p.key = PropertyKey::new("cp", Density::DIMS);
                    p.value = PropertyValue::Scalar { value: 900.0, dims: Density::DIMS };
                }
                1 => {
                    p.value = PropertyValue::Curve { abscissa: "omega".into(),
                        abscissa_dims: fs_qty::Frequency::DIMS,
                        knots: vec![(100.0, 900.0), (2e6, 1000.0)], dims: SPECIFIC_HEAT_CAPACITY_DIMS };
                    p.interpolation = InterpolationPolicy::LinearInside;
                }
                _ => p.value = PropertyValue::Scalar { value: 0.0, dims: SPECIFIC_HEAT_CAPACITY_DIMS },
            }
        });
        assert!(bind_region(&source.regions()[0], &properties, (300.0, 1e6)).is_err());
    }
}

#[test]
fn g0_empty_support_and_nonhomogeneous_section_refuse() {
    let (anchored, thermal) = fixture_with(PlateMaterialModel::Isotropic, query(), |claims| {
        let mut anchor = claims[0].clone();
        anchor.key = PropertyKey::new(crate::shear_shell::CONSTANT_LOSS_REFERENCE_OMEGA,
            crate::shear_shell::CONSTANT_LOSS_REFERENCE_OMEGA_DIMS);
        anchor.value = PropertyValue::Scalar { value: 1000.0,
            dims: crate::shear_shell::CONSTANT_LOSS_REFERENCE_OMEGA_DIMS };
        claims.push(anchor);
    });
    assert!(bind_region(&anchored.regions()[0], &thermal, (300.0, 1e6)).is_err());
    let (source, properties) = fixture(PlateMaterialModel::Isotropic);
    let mut region = source.regions()[0].clone();
    assert!(bind_region(&region, &properties, (2e6, 3e6)).is_err());
    let (a, d) = region.section().principal_matrices();
    region.section = ShearShellSection::new_coupled(
        *a, [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], *d,
        [1e6, 0.0, 1e6], 1.0, region.section().areal_mass_kg_m2(),
        0.0, region.section().rotary_inertia_kg(), 0.0,
    ).unwrap();
    assert!(bind_region(&region, &properties, (300.0, 1e6)).is_err());
}
