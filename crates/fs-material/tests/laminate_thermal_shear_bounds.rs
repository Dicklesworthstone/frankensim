//! G0/G1 numerical oracles for whole-box nonlocal shear capacity admission.
//! Authored material cards establish numerical behavior, not measured fidelity.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::laminate::thermal::shear::bounds::LaminateShearCapacityBound;
use fs_material::laminate::thermal::shear::{
    LaminateShearThermoelasticError, ResolvedLaminateShearThermoelastic,
};
use fs_material::laminate::thermal::source::{
    LaminaElasticSource, LaminaExpansionSource, LaminaSourceInput, LaminateThermalSource,
};
use fs_material::laminate::thermal::thermoelastic::AffineLaminaThermoelastic;
use fs_material::laminate::{LaminateShearModel, ReferenceSurface};
use fs_material::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS, IsotropicElasticParameters,
    MaterialPropertySelection, ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES as CTE,
    ScalarAdmissibility, ScalarPropertyRequirement,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use std::sync::Arc;

const LIMIT: f64 = 0.04;
const BOX: [f64; 2] = [295.0, 305.0];
fn absolute_temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", absolute_temperature(), t)
        .unwrap()
}
fn input(
    index: usize,
    thickness: f64,
    angle: f64,
    isotropic: bool,
    variable: bool,
) -> LaminaSourceInput {
    input_with_rate(index, thickness, angle, isotropic, variable, 1.0)
}
fn input_with_rate(
    index: usize,
    thickness: f64,
    angle: f64,
    isotropic: bool,
    variable: bool,
    rate: f64,
) -> LaminaSourceInput {
    let scale = [1.0, 0.73, 1.38][index];
    let mut coefficients = if isotropic {
        vec![
            (
                "young_modulus",
                Pressure::DIMS,
                1200.0 * scale,
                -2.3 * scale,
            ),
            ("poisson_ratio", Dims::NONE, 0.23, 0.0007),
        ]
    } else {
        vec![
            (
                "young_modulus_1",
                Pressure::DIMS,
                1600.0 * scale,
                -2.0 * scale,
            ),
            (
                "young_modulus_2",
                Pressure::DIMS,
                900.0 * scale,
                0.8 * scale,
            ),
            ("poisson_ratio_12", Dims::NONE, 0.18, 0.0004),
            (
                "shear_modulus_12",
                Pressure::DIMS,
                550.0 * scale,
                -0.9 * scale,
            ),
            (
                "shear_modulus_31",
                Pressure::DIMS,
                480.0 * scale,
                0.3 * scale,
            ),
            (
                "shear_modulus_23",
                Pressure::DIMS,
                370.0 * scale,
                -0.2 * scale,
            ),
        ]
    };
    for coefficient in &mut coefficients {
        coefficient.3 *= if variable { rate } else { 0.0 };
    }
    coefficients.extend([
        (CTE[0], INVERSE_TEMPERATURE_DIMS, 1e-5, 0.0),
        (CTE[1], INVERSE_TEMPERATURE_DIMS, -2e-5, 0.0),
        (CTE[2], INVERSE_TEMPERATURE_DIMS, 3e-5, 0.0),
        (DENSITY_PROPERTY, Density::DIMS, 900.0, 0.0),
    ]);
    let mut claims = ClaimSet::new();
    for (name, dims, value, slope) in coefficients {
        let density = name == DENSITY_PROPERTY;
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if density {
                    PropertyValue::Scalar { value, dims }
                } else {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![
                            (280.0, value - 20.0 * slope),
                            (300.0, value),
                            (360.0, value + 60.0 * slope),
                        ],
                        dims,
                    }
                },
                interpolation: if density {
                    InterpolationPolicy::ConstantWithinValidity
                } else {
                    InterpolationPolicy::LinearInside
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    absolute_temperature(),
                    if density { 300.0 } else { 280.0 },
                    if density { 300.0 } else { 360.0 },
                ),
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored shear box oracle".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let req = |name, dims| {
        ScalarPropertyRequirement::try_new(
            name,
            dims,
            if dims == Pressure::DIMS {
                ScalarAdmissibility::StrictlyPositive
            } else {
                ScalarAdmissibility::Finite
            },
        )
        .unwrap()
    };
    LaminaSourceInput {
        card: Arc::new(
            MaterialCard::assemble(
                MaterialStateId {
                    chemistry: format!("shear-box-{index}-{isotropic}-{variable}"),
                    phase: "solid".into(),
                    process: "numerical fixture".into(),
                    revision: 0,
                },
                claims,
                vec![],
            )
            .unwrap(),
        ),
        reference: point(300.0),
        thickness_m: thickness,
        orientation_rad: angle,
        elasticity: if isotropic {
            LaminaElasticSource::Isotropic(IsotropicElasticParameters::YoungPoisson {
                young_modulus: req("young_modulus", Pressure::DIMS),
                poisson_ratio: req("poisson_ratio", Dims::NONE),
            })
        } else {
            LaminaElasticSource::OrthotropicPlane12
        },
        expansion: LaminaExpansionSource::Orthotropic,
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}
fn source(inputs: Vec<LaminaSourceInput>) -> LaminateThermalSource {
    LaminateThermalSource::try_new(inputs, ReferenceSurface::Midsurface, || false).unwrap()
}
fn coupled_source(inputs: Vec<LaminaSourceInput>, angle: f64) -> LaminateThermalSource {
    LaminateThermalSource::try_new_with_shear_model(
        inputs,
        ReferenceSurface::Midsurface,
        LaminateShearModel::EquilibratedCylindrical {
            bending_axis_rad: angle,
        },
        || false,
    )
    .unwrap()
}
fn laws(source: &LaminateThermalSource, range: [f64; 2]) -> Vec<AffineLaminaThermoelastic> {
    source
        .layers()
        .iter()
        .map(|s| AffineLaminaThermoelastic::new(s, range, LIMIT).unwrap())
        .collect()
}

#[test]
fn one_layer_eliminates_q_and_bounds_the_closed_rational_shear_curvature() {
    let h = 0.04;
    let orthotropic = source(vec![input(0, h, 0.37, false, true)]);
    let current = orthotropic.at_temperatures(&[321.0], || false).unwrap();
    let bound =
        LaminateShearCapacityBound::try_new(&laws(&orthotropic, [280.0, 360.0]), &current, || {
            false
        })
        .unwrap();
    // H = (5h/6) R G R^T is affine despite independently variable Q(T).
    assert_eq!(bound.entropy_hessian_absolute_row_sums_j_m2_k2(), &[0.0]);

    let isotropic = source(vec![input(0, h, 0.0, true, true)]);
    let current = isotropic.at_temperatures(&[300.0], || false).unwrap();
    let bound =
        LaminateShearCapacityBound::try_new(&laws(&isotropic, [280.0, 360.0]), &current, || false)
            .unwrap();
    let upper = bound.entropy_hessian_absolute_row_sums_j_m2_k2()[0];
    assert!(upper.is_finite() && upper > 0.0);
    // Independent G = E/[2(1+nu)], G'' = E nu'^2/(1+nu)^3 - E'nu'/(1+nu)^2.
    for i in 0..=80 {
        let t = 280.0 + f64::from(i);
        let e = 1200.0 - 2.3 * (t - 300.0);
        let u = 1.23 + 0.0007 * (t - 300.0);
        let g_second = e * 0.0007_f64.powi(2) / u.powi(3) + 2.3 * 0.0007 / u.powi(2);
        let entropy_curvature = (5.0 / 6.0) * h * LIMIT * LIMIT * g_second.abs();
        assert!(entropy_curvature <= upper, "{entropy_curvature} > {upper}");
    }
}

#[test]
fn independent_temperature_box_encloses_all_mixed_shear_entropy_rows() {
    let source = source(vec![
        input(0, 0.013, 0.0, false, true),
        input(1, 0.007, 0.29, false, true),
        input(2, 0.021, -0.43, false, true),
    ]);
    let laws = laws(&source, BOX);
    let current = source.at_temperatures(&[300.0; 3], || false).unwrap();
    let bound = LaminateShearCapacityBound::try_new(&laws, &current, || false).unwrap();
    let bounds = bound.entropy_hessian_absolute_row_sums_j_m2_k2();
    assert!(bounds.iter().all(|v| v.is_finite() && *v > 0.0));
    let area = 1.7;
    for corner in 0..8 {
        let t: Vec<_> = (0..3)
            .map(|i| {
                if corner & (1 << i) == 0 {
                    BOX[0]
                } else {
                    BOX[1]
                }
            })
            .collect();
        let native = source.at_temperatures(&t, || false).unwrap();
        let shear = ResolvedLaminateShearThermoelastic::try_new(&laws, &native, || false).unwrap();
        for gamma in [[0.013, -0.009], [-0.02, -0.012], [0.0, 0.025]] {
            let mut row_sums = [0.0; 3];
            for j in 0..3 {
                let mut direction = [0.0; 3];
                direction[j] = 1.0;
                let response = shear
                    .thermal_response(gamma, area, &direction, || false)
                    .unwrap();
                for i in 0..3 {
                    row_sums[i] += response.entropy_derivative_j_k2[i].abs();
                }
            }
            for i in 0..3 {
                assert!(
                    row_sums[i] <= area * bounds[i],
                    "cell {i}: {} > {}",
                    row_sums[i],
                    area * bounds[i]
                );
            }
        }
    }
    // Rebinding a different current temperature does not shrink a whole-box
    // bound into a local-only admission.
    let another = source
        .at_temperatures(&[295.0, 305.0, 302.0], || false)
        .unwrap();
    let same = LaminateShearCapacityBound::try_new(&laws, &another, || false).unwrap();
    assert_eq!(bound, same);
}

#[test]
fn fixed_stack_has_zero_curvature_and_refusals_preserve_source_identity() {
    let source = source(vec![
        input(0, 0.02, 0.0, false, false),
        input(1, 0.01, 0.3, false, false),
    ]);
    let laws = laws(&source, BOX);
    let current = source.at_temperatures(&[299.0, 302.0], || false).unwrap();
    let bound = LaminateShearCapacityBound::try_new(&laws, &current, || false).unwrap();
    assert_eq!(
        bound.entropy_hessian_absolute_row_sums_j_m2_k2(),
        &[0.0, 0.0]
    );
    assert!(matches!(
        LaminateShearCapacityBound::try_new(&laws, &current, || true),
        Err(LaminateShearThermoelasticError::Cancelled)
    ));
    let mut calls = 0;
    assert!(matches!(
        LaminateShearCapacityBound::try_new(&laws, &current, || {
            calls += 1;
            calls == 25
        }),
        Err(LaminateShearThermoelasticError::Cancelled)
    ));
    let alien = crate::source(vec![
        input(0, 0.02, 0.0, false, false),
        input(1, 0.01, 0.3, false, false),
    ]);
    let alien = alien.at_temperatures(&[299.0, 302.0], || false).unwrap();
    assert!(LaminateShearCapacityBound::try_new(&laws, &alien, || false).is_err());
    assert!(LaminateShearCapacityBound::try_new(&laws[..1], &current, || false).is_err());
    let outside = source.at_temperatures(&[310.0, 302.0], || false).unwrap();
    assert!(LaminateShearCapacityBound::try_new(&laws, &outside, || false).is_err());
}

#[test]
fn coupled_box_encloses_mixed_entropy_with_co_rotated_physical_probes() {
    // Independent temperatures and unequal thicknesses retain every A^-1 B
    // coupling. These state evaluations check the analytic bound; they are
    // deliberately not inputs to its certificate.
    let range = [299.5, 300.5];
    let mut reference_rows = Vec::new();
    for rotation in [0.0, 0.43] {
        let source = coupled_source(
            vec![
                input(0, 0.013, 0.10 + rotation, false, true),
                input(1, 0.007, 0.29 + rotation, false, true),
                input(2, 0.021, -0.43 + rotation, false, true),
            ],
            0.17 + rotation,
        );
        let laws = laws(&source, range);
        let current = source.at_temperatures(&[300.0; 3], || false).unwrap();
        let bound = LaminateShearCapacityBound::try_new(&laws, &current, || false).unwrap();
        let rows = bound.entropy_hessian_absolute_row_sums_j_m2_k2();
        assert!(rows.iter().all(|x| x.is_finite() && *x > 0.0));
        let (s, c) = rotation.sin_cos();
        let gamma = [c * 0.013 + s * 0.009, s * 0.013 - c * 0.009];
        for (state, temperatures) in [
            [299.5, 300.5, 299.5],
            [300.5, 299.5, 300.5],
            [299.73, 300.26, 300.08],
            [300.0; 3],
        ]
        .into_iter()
        .enumerate()
        {
            let native = source.at_temperatures(&temperatures, || false).unwrap();
            let shear =
                ResolvedLaminateShearThermoelastic::try_new(&laws, &native, || false).unwrap();
            let mut actual = [0.0; 3];
            let mut mixed = 0.0_f64;
            for j in 0..3 {
                let mut direction = [0.0; 3];
                direction[j] = 1.0;
                let response = shear
                    .thermal_response(gamma, 1.7, &direction, || false)
                    .unwrap();
                for i in 0..3 {
                    actual[i] += response.entropy_derivative_j_k2[i].abs();
                    if i != j {
                        mixed = mixed.max(response.entropy_derivative_j_k2[i].abs());
                    }
                }
            }
            assert!(
                mixed > 1e-12,
                "original cells must retain cross-temperature shear entropy"
            );
            for i in 0..3 {
                assert!(
                    actual[i] <= 1.7 * rows[i],
                    "cell {i}: {} > {}",
                    actual[i],
                    1.7 * rows[i]
                );
            }
            if rotation == 0.0 {
                reference_rows.push(actual);
            } else {
                for i in 0..3 {
                    assert!(
                        (actual[i] - reference_rows[state][i]).abs()
                            <= 1e-10 * actual[i].abs().max(reference_rows[state][i].abs())
                    );
                }
            }
            // Rebinding cannot turn a whole-box certificate into a point test.
            let rebound = LaminateShearCapacityBound::try_new(&laws, &native, || false).unwrap();
            assert_eq!(bound, rebound);
        }
    }
}

#[test]
fn coupled_domain_proof_preserves_affine_limit_and_refuses_unseparated_boxes() {
    let isotropic = coupled_source(vec![input(0, 0.04, 0.0, true, true)], 0.17);
    let current = isotropic.at_temperatures(&[300.0], || false).unwrap();
    let bound =
        LaminateShearCapacityBound::try_new(&laws(&isotropic, [280.0, 360.0]), &current, || false)
            .unwrap();
    let upper = bound.entropy_hessian_absolute_row_sums_j_m2_k2()[0];
    for t in [280.0, 300.0, 327.1, 360.0] {
        let e = 1200.0 - 2.3 * (t - 300.0);
        let u: f64 = 1.23 + 0.0007 * (t - 300.0);
        let g_second = e * 0.0007_f64.powi(2) / u.powi(3) + 2.3 * 0.0007 / u.powi(2);
        // Independent homogeneous H=(5h/6)G, with both shear components
        // at their original principal strain limits.
        assert!((5.0 / 6.0) * 0.04 * LIMIT * LIMIT * g_second.abs() <= upper);
    }
    let fixed = coupled_source(
        vec![
            input(0, 0.02, 0.1, false, false),
            input(1, 0.01, -0.3, false, false),
        ],
        0.21,
    );
    let current = fixed.at_temperatures(&[299.0, 302.0], || false).unwrap();
    let fixed_laws = laws(&fixed, BOX);
    let bound = LaminateShearCapacityBound::try_new(&fixed_laws, &current, || false).unwrap();
    assert_eq!(
        bound.entropy_hessian_absolute_row_sums_j_m2_k2(),
        &[0.0, 0.0]
    );
    assert!(matches!(
        LaminateShearCapacityBound::try_new(&fixed_laws, &current, || true),
        Err(LaminateShearThermoelasticError::Cancelled)
    ));
    let alien = coupled_source(
        vec![
            input(0, 0.02, 0.1, false, false),
            input(1, 0.01, -0.3, false, false),
        ],
        0.21,
    );
    let alien = alien.at_temperatures(&[299.0, 302.0], || false).unwrap();
    assert!(LaminateShearCapacityBound::try_new(&fixed_laws, &alien, || false).is_err());

    let fast = coupled_source(vec![input_with_rate(0, 0.04, 0.7, false, true, 10.0)], 0.17);
    let current = fast.at_temperatures(&[300.0], || false).unwrap();
    let wide = laws(&fast, [280.0, 360.0]);
    ResolvedLaminateShearThermoelastic::try_new(&wide, &current, || false).unwrap();
    // The selected state and original positive affine chart are valid. The
    // wide rotated tensor interval cannot prove all inverse pivots; no sampled
    // repair, clipping or pointwise heat-capacity substitution is permitted.
    assert!(LaminateShearCapacityBound::try_new(&wide, &current, || false).is_err());
    let narrow = laws(&fast, [299.9, 300.1]);
    let bound = LaminateShearCapacityBound::try_new(&narrow, &current, || false).unwrap();
    assert_eq!(bound.entropy_hessian_absolute_row_sums_j_m2_k2(), &[0.0]);
}
