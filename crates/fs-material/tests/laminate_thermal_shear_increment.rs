//! G1/G3 finite-temperature work oracles for the native nonlocal shear model.
//! Authored source cards test the discretization, not empirical materials.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::laminate::thermal::shear::{
    LaminateShearThermoelasticError, ResolvedLaminateShearThermoelastic,
};
use fs_material::laminate::thermal::source::{
    LaminaElasticSource, LaminaExpansionSource, LaminaSourceInput, LaminateThermalSource,
};
use fs_material::laminate::thermal::thermoelastic::AffineLaminaThermoelastic;
use fs_material::laminate::{LaminateShearModel, ReferenceSurface};
use fs_material::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS, MaterialPropertySelection,
    ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES as CTE,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use std::sync::Arc;

const RANGE: [f64; 2] = [280.0, 360.0];
const LIMIT: f64 = 0.05;
const MODELS: [LaminateShearModel; 2] = [
    LaminateShearModel::DirectionalWhitney,
    LaminateShearModel::EquilibratedCylindrical {
        bending_axis_rad: 0.21,
    },
];

fn temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn input(index: usize, thickness_m: f64, orientation_rad: f64) -> LaminaSourceInput {
    let scale = [1.0, 0.7, 1.4][index];
    let slope = [1.0, -0.6, 1.8][index];
    let mut claims = ClaimSet::new();
    for (name, dims, value, gradient) in [
        (
            "young_modulus_1",
            Pressure::DIMS,
            1700.0 * scale,
            -2.1 * slope,
        ),
        (
            "young_modulus_2",
            Pressure::DIMS,
            850.0 * scale,
            0.9 * slope,
        ),
        ("poisson_ratio_12", Dims::NONE, 0.18, 0.0004 * slope),
        (
            "shear_modulus_12",
            Pressure::DIMS,
            520.0 * scale,
            -0.8 * slope,
        ),
        (
            "shear_modulus_31",
            Pressure::DIMS,
            470.0 * scale,
            0.7 * slope,
        ),
        (
            "shear_modulus_23",
            Pressure::DIMS,
            330.0 * scale,
            -0.4 * slope,
        ),
        (CTE[0], INVERSE_TEMPERATURE_DIMS, 1e-4, 0.0),
        (CTE[1], INVERSE_TEMPERATURE_DIMS, -3e-5, 0.0),
        (CTE[2], INVERSE_TEMPERATURE_DIMS, 5e-5, 0.0),
        (DENSITY_PROPERTY, Density::DIMS, 900.0, 0.0),
    ] {
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
                            (280.0, value - 20.0 * gradient),
                            (300.0, value),
                            (360.0, value + 60.0 * gradient),
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
                    temperature(),
                    if density { 300.0 } else { RANGE[0] },
                    if density { 300.0 } else { RANGE[1] },
                ),
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored symmetric shear-coordinate work oracle".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    LaminaSourceInput {
        card: Arc::new(
            MaterialCard::assemble(
                MaterialStateId {
                    chemistry: format!("shear-increment-{index}"),
                    phase: "solid".into(),
                    process: "numerical oracle".into(),
                    revision: 0,
                },
                claims,
                vec![],
            )
            .unwrap(),
        ),
        reference: QueryPoint::new()
            .with_quantity("T", temperature(), 300.0)
            .unwrap(),
        thickness_m,
        orientation_rad,
        elasticity: LaminaElasticSource::OrthotropicPlane12,
        expansion: LaminaExpansionSource::Orthotropic,
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}

fn source() -> LaminateThermalSource {
    source_with_model(LaminateShearModel::DirectionalWhitney)
}

fn source_with_model(model: LaminateShearModel) -> LaminateThermalSource {
    LaminateThermalSource::try_new_with_shear_model(
        vec![
            input(0, 0.013, 0.1),
            input(1, 0.007, 0.55),
            input(2, 0.021, -0.32),
        ],
        ReferenceSurface::Midsurface,
        model,
        || false,
    )
    .unwrap()
}

fn laws(source: &LaminateThermalSource, range: [f64; 2]) -> Vec<AffineLaminaThermoelastic> {
    source
        .layers()
        .iter()
        .map(|source| AffineLaminaThermoelastic::new(source, range, LIMIT).unwrap())
        .collect()
}

fn resolved(
    source: &LaminateThermalSource,
    laws: &[AffineLaminaThermoelastic],
    t: &[f64],
) -> ResolvedLaminateShearThermoelastic {
    ResolvedLaminateShearThermoelastic::try_new(
        laws,
        &source.at_temperatures(t, || false).unwrap(),
        || false,
    )
    .unwrap()
}

fn native(source: &LaminateThermalSource, t: &[f64]) -> [f64; 4] {
    *source
        .at_temperatures(t, || false)
        .unwrap()
        .prepared()
        .section()
        .shear_pa_m()
}

fn close(actual: f64, expected: f64, atol: f64, rtol: f64) {
    assert!(
        (actual - expected).abs() <= atol + rtol * actual.abs().max(expected.abs()),
        "{actual:.16e} != {expected:.16e}"
    );
}

fn bilinear(left: [f64; 2], h: [f64; 4], right: [f64; 2]) -> f64 {
    left[0] * (h[0] * right[0] + h[1] * right[1]) + left[1] * (h[2] * right[0] + h[3] * right[1])
}

#[test]
fn homogeneous_increment_and_reciprocal_heat_have_closed_form() {
    for model in MODELS {
        check_homogeneous_increment(model);
    }
}

fn check_homogeneous_increment(model: LaminateShearModel) {
    let angle: f64 = 0.37;
    let thickness = 0.04;
    let source = LaminateThermalSource::try_new_with_shear_model(
        vec![input(0, thickness, angle)],
        ReferenceSurface::Midsurface,
        model,
        || false,
    )
    .unwrap();
    let laws = laws(&source, RANGE);
    let a = resolved(&source, &laws, &[291.0]);
    let b = resolved(&source, &laws, &[341.0]);
    let step = a.temperature_increment(&b, || false).unwrap();
    let increment = step.stiffness_coordinate_increments_pa_m()[0];
    let (sin, cos) = angle.sin_cos();
    let x = (5.0 / 6.0) * thickness * 0.7 * 50.0;
    let y = (5.0 / 6.0) * thickness * -0.4 * 50.0;
    let expected = [
        cos * cos * x + sin * sin * y,
        cos * sin * (x - y),
        cos * sin * (x - y),
        sin * sin * x + cos * cos * y,
    ];
    for k in 0..4 {
        close(increment[k], expected[k], 2e-13, 2e-12);
        close(step.stiffness_increment_residual_pa_m()[k], 0.0, 2e-14, 0.0);
    }

    let area = 0.7;
    let gamma0 = [0.013, -0.009];
    let gamma1 = [0.015, -0.012];
    let s0 = a.state(gamma0, area, || false).unwrap();
    let s1 = b.state(gamma1, area, || false).unwrap();
    let force_work = 0.5
        * area
        * (0..2)
            .map(|k| {
                (s0.force_resultant_n_m[k] + s1.force_resultant_n_m[k]) * (gamma1[k] - gamma0[k])
            })
            .sum::<f64>();
    let parameter_work = 0.5 * area * bilinear(gamma0, increment, gamma1);
    let caloric = -parameter_work - (341.0 * s1.entropy_j_k[0] - 291.0 * s0.entropy_j_k[0]);
    close(
        caloric,
        force_work - (s1.internal_energy_j - s0.internal_energy_j),
        3e-16,
        5e-12,
    );
}

#[test]
fn mixed_cells_telescope_reverse_and_keep_unchanged_coordinate_exactly_zero() {
    for model in MODELS {
        check_mixed_cells(model);
    }
}

fn check_mixed_cells(model: LaminateShearModel) {
    let source = source_with_model(model);
    let laws = laws(&source, RANGE);
    let t0 = [289.0, 305.0, 331.0];
    let t1 = [345.0, 305.0, 293.0];
    let a = resolved(&source, &laws, &t0);
    let b = resolved(&source, &laws, &t1);
    let step = a.temperature_increment(&b, || false).unwrap();
    let back = b.temperature_increment(&a, || false).unwrap();
    let increments = step.stiffness_coordinate_increments_pa_m();
    assert_eq!(increments[1], [0.0; 4]);
    let first = native(&source, &[t1[0], t0[1], t0[2]]);
    let last = native(&source, &[t0[0], t0[1], t1[2]]);
    for k in 0..4 {
        let expected_first =
            0.5 * (first[k] - a.stiffness_pa_m()[k] + b.stiffness_pa_m()[k] - last[k]);
        close(increments[0][k], expected_first, 2e-12, 2e-11);
        close(
            increments.iter().map(|m| m[k]).sum(),
            b.stiffness_pa_m()[k] - a.stiffness_pa_m()[k],
            2e-12,
            2e-11,
        );
        for i in 0..3 {
            close(
                increments[i][k],
                -back.stiffness_coordinate_increments_pa_m()[i][k],
                2e-12,
                2e-11,
            );
        }
    }
    let unchanged = a.temperature_increment(&a.clone(), || false).unwrap();
    assert_eq!(
        unchanged.stiffness_coordinate_increments_pa_m(),
        &[[0.0; 4]; 3]
    );
    assert_eq!(unchanged.stiffness_increment_residual_pa_m(), &[0.0; 4]);

    // Per-cell reciprocal transfers recover the change in the ONE stack
    // internal energy, including mixed-temperature entropy at the fixed cell.
    let (area, gamma0, gamma1) = (0.7, [0.013, -0.009], [0.015, -0.012]);
    let s0 = a.state(gamma0, area, || false).unwrap();
    let s1 = b.state(gamma1, area, || false).unwrap();
    let force_work = 0.5
        * area
        * (0..2)
            .map(|k| {
                (s0.force_resultant_n_m[k] + s1.force_resultant_n_m[k]) * (gamma1[k] - gamma0[k])
            })
            .sum::<f64>();
    let mut caloric = 0.0;
    for i in 0..3 {
        let entropy_work = t1[i] * s1.entropy_j_k[i] - t0[i] * s0.entropy_j_k[i];
        let transfer = -0.5 * area * bilinear(gamma0, increments[i], gamma1) - entropy_work;
        let reverse = -0.5
            * area
            * bilinear(
                gamma1,
                back.stiffness_coordinate_increments_pa_m()[i],
                gamma0,
            )
            + entropy_work;
        close(transfer, -reverse, 3e-16, 5e-12);
        caloric += transfer;
    }
    close(
        caloric,
        force_work - (s1.internal_energy_j - s0.internal_energy_j),
        3e-16,
        5e-12,
    );
}

#[test]
fn symmetric_coordinate_path_converges_to_source_temperature_derivatives() {
    for model in MODELS {
        check_coordinate_refinement(model);
    }
}

fn check_coordinate_refinement(model: LaminateShearModel) {
    let source = source_with_model(model);
    let laws = laws(&source, RANGE);
    let center = [307.0, 319.0, 326.0];
    let direction = [7.0, -5.0, 11.0];
    let at_center = resolved(&source, &laws, &center);
    let error = |width: f64| {
        let t0 = core::array::from_fn::<_, 3, _>(|i| center[i] - 0.5 * width * direction[i]);
        let t1 = core::array::from_fn::<_, 3, _>(|i| center[i] + 0.5 * width * direction[i]);
        let a = resolved(&source, &laws, &t0);
        let b = resolved(&source, &laws, &t1);
        let step = a.temperature_increment(&b, || false).unwrap();
        let mut maximum = 0.0_f64;
        for (i, increment) in step
            .stiffness_coordinate_increments_pa_m()
            .iter()
            .enumerate()
        {
            for (k, &value) in increment.iter().enumerate() {
                let expected =
                    at_center.stiffness_temperature_derivatives_pa_m_k()[i][k] * direction[i];
                maximum = maximum.max((value / width - expected).abs());
            }
        }
        maximum
    };
    let coarse = error(1.0);
    let fine = error(0.5);
    assert!(
        coarse > 1e-10,
        "fixture must resolve nonlocal path truncation"
    );
    assert!(fine < 0.3 * coarse, "{coarse:e} -> {fine:e}");
}

#[test]
fn foreign_source_different_chart_and_cancellation_refuse_transactionally() {
    let source = source();
    let original_laws = laws(&source, RANGE);
    let a = resolved(&source, &original_laws, &[290.0, 310.0, 330.0]);
    let b = resolved(&source, &original_laws, &[345.0, 307.0, 295.0]);
    let other = LaminateThermalSource::try_new(
        source
            .layers()
            .iter()
            .map(|layer| layer.input().clone())
            .collect(),
        ReferenceSurface::Midsurface,
        || false,
    )
    .unwrap();
    let foreign = resolved(&other, &laws(&other, RANGE), &[345.0, 307.0, 295.0]);
    assert!(a.temperature_increment(&foreign, || false).is_err());
    let other_chart = resolved(
        &source,
        &laws(&source, [285.0, 355.0]),
        &[345.0, 307.0, 295.0],
    );
    assert!(a.temperature_increment(&other_chart, || false).is_err());
    assert!(matches!(
        a.temperature_increment(&b, || true),
        Err(LaminateShearThermoelasticError::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        a.temperature_increment(&b, || {
            polls += 1;
            polls >= 16
        }),
        Err(LaminateShearThermoelasticError::Cancelled)
    ));
    let complete = a.temperature_increment(&b, || false).unwrap();
    assert_eq!(complete.stiffness_coordinate_increments_pa_m().len(), 3);
}
