//! G1/G3 numerical oracles for nonlocal thermal shear. These authored material
//! cards test the numerical model; no empirical laminate validation is claimed.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::laminate::ReferenceSurface;
use fs_material::laminate::thermal::shear::{
    LaminateShearThermoelasticError, ResolvedLaminateShearThermoelastic,
};
use fs_material::laminate::thermal::source::{
    LaminaElasticSource, LaminaExpansionSource, LaminaSourceInput, LaminateThermalSource,
};
use fs_material::laminate::thermal::thermoelastic::AffineLaminaThermoelastic;
use fs_material::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS, MaterialPropertySelection,
    ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES as CTE,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use std::sync::Arc;

const RANGE: [f64; 2] = [280.0, 360.0];
const LIMIT: f64 = 0.05;
const GAMMA: [f64; 2] = [0.013, -0.009];
const AREA: f64 = 0.7;

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", absolute_temperature(), t)
        .unwrap()
}
fn absolute_temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn input(index: usize, thickness: f64, angle: f64) -> LaminaSourceInput {
    let scale = [1.0, 0.7, 1.4][index];
    let slope_scale = [1.0, -0.6, 1.8][index];
    let constants = [
        (
            "young_modulus_1",
            Pressure::DIMS,
            1700.0 * scale,
            -2.1 * slope_scale,
        ),
        (
            "young_modulus_2",
            Pressure::DIMS,
            850.0 * scale,
            0.9 * slope_scale,
        ),
        ("poisson_ratio_12", Dims::NONE, 0.18, 0.0004 * slope_scale),
        (
            "shear_modulus_12",
            Pressure::DIMS,
            520.0 * scale,
            -0.8 * slope_scale,
        ),
        (
            "shear_modulus_31",
            Pressure::DIMS,
            470.0 * scale,
            0.7 * slope_scale,
        ),
        (
            "shear_modulus_23",
            Pressure::DIMS,
            330.0 * scale,
            -0.4 * slope_scale,
        ),
        (CTE[0], INVERSE_TEMPERATURE_DIMS, 1e-4, 0.0),
        (CTE[1], INVERSE_TEMPERATURE_DIMS, -3e-5, 0.0),
        (CTE[2], INVERSE_TEMPERATURE_DIMS, 5e-5, 0.0),
        (DENSITY_PROPERTY, Density::DIMS, 900.0, 0.0),
    ];
    let mut claims = ClaimSet::new();
    for (name, dims, value, slope) in constants {
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
                    if density { 300.0 } else { RANGE[0] },
                    if density { 300.0 } else { RANGE[1] },
                ),
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored nonlocal thermal shear oracle".into(),
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
                    chemistry: format!("thermal-shear-{index}"),
                    phase: "solid".into(),
                    process: "numerical derivative fixture".into(),
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
        elasticity: LaminaElasticSource::OrthotropicPlane12,
        expansion: LaminaExpansionSource::Orthotropic,
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}
fn source(inputs: Vec<LaminaSourceInput>, reference: ReferenceSurface) -> LaminateThermalSource {
    LaminateThermalSource::try_new(inputs, reference, || false).unwrap()
}
fn laws(source: &LaminateThermalSource) -> Vec<AffineLaminaThermoelastic> {
    source
        .layers()
        .iter()
        .map(|s| AffineLaminaThermoelastic::new(s, RANGE, LIMIT).unwrap())
        .collect()
}
fn resolved(
    source: &LaminateThermalSource,
    laws: &[AffineLaminaThermoelastic],
    t: &[f64],
) -> ResolvedLaminateShearThermoelastic {
    let native = source.at_temperatures(t, || false).unwrap();
    ResolvedLaminateShearThermoelastic::try_new(laws, &native, || false).unwrap()
}
fn native_h(source: &LaminateThermalSource, t: &[f64]) -> [f64; 4] {
    *source
        .at_temperatures(t, || false)
        .unwrap()
        .prepared()
        .section()
        .shear_pa_m()
}
fn shift(t: &[f64], direction: &[f64], h: f64) -> Vec<f64> {
    t.iter().zip(direction).map(|(&t, &v)| t + h * v).collect()
}
fn close(a: f64, b: f64, atol: f64, rtol: f64) {
    assert!(
        (a - b).abs() <= atol + rtol * a.abs().max(b.abs()),
        "{a:.16e} != {b:.16e}"
    );
}
fn quad(h: [f64; 4], gamma: [f64; 2], area: f64) -> f64 {
    0.5 * area
        * (gamma[0] * (h[0] * gamma[0] + h[1] * gamma[1])
            + gamma[1] * (h[2] * gamma[0] + h[3] * gamma[1]))
}
fn multilayer_inputs(angle: f64) -> Vec<LaminaSourceInput> {
    [(0.013, 0.1), (0.007, 0.55), (0.021, -0.32)]
        .into_iter()
        .enumerate()
        .map(|(i, (h, a))| input(i, h, a + angle))
        .collect()
}

#[test]
fn homogeneous_rotated_layer_has_closed_shear_potential_and_no_affine_g_capacity() {
    let angle: f64 = 0.37;
    let h = 0.04;
    let t = 321.0;
    let source = source(vec![input(0, h, angle)], ReferenceSurface::Midsurface);
    let laws = laws(&source);
    let shear = resolved(&source, &laws, &[t]);
    let derivative = shear.directional_derivatives(&[1.0], || false).unwrap();
    let (s, c) = angle.sin_cos();
    let rotate = |x: f64, y: f64| {
        [
            c * c * x + s * s * y,
            c * s * (x - y),
            c * s * (x - y),
            s * s * x + c * c * y,
        ]
    };
    let expected = rotate(470.0 + 0.7 * (t - 300.0), 330.0 - 0.4 * (t - 300.0));
    let expected_d = rotate(0.7, -0.4);
    for i in 0..4 {
        close(
            derivative[0][i],
            (5.0 / 6.0) * h * expected[i],
            2e-13,
            5e-13,
        );
        close(
            derivative[1][i],
            (5.0 / 6.0) * h * expected_d[i],
            2e-14,
            5e-12,
        );
        close(derivative[2][i], 0.0, 2e-14, 0.0);
    }
    let state = shear.state(GAMMA, AREA, || false).unwrap();
    let principal = [c * GAMMA[0] + s * GAMMA[1], -s * GAMMA[0] + c * GAMMA[1]];
    let psi = (5.0 / 12.0)
        * h
        * AREA
        * ((470.0 + 0.7 * (t - 300.0)) * principal[0].powi(2)
            + (330.0 - 0.4 * (t - 300.0)) * principal[1].powi(2));
    let entropy =
        -(5.0 / 12.0) * h * AREA * (0.7 * principal[0].powi(2) - 0.4 * principal[1].powi(2));
    close(state.free_energy_j, psi, 1e-16, 2e-12);
    close(state.entropy_j_k[0], entropy, 1e-16, 2e-11);
    close(state.internal_energy_j, psi + t * entropy, 1e-16, 2e-11);
    let capacity = shear
        .thermal_response(GAMMA, AREA, &[1.0], || false)
        .unwrap();
    close(capacity.internal_energy_derivative_j_k, 0.0, 1e-14, 0.0);
    // A mechanical derivative independently returns the native force resultant.
    for component in 0..2 {
        let mut plus = GAMMA;
        let mut minus = GAMMA;
        plus[component] += 1e-6;
        minus[component] -= 1e-6;
        let force = (shear.state(plus, AREA, || false).unwrap().free_energy_j
            - shear.state(minus, AREA, || false).unwrap().free_energy_j)
            / (2e-6 * AREA);
        close(force, state.force_resultant_n_m[component], 1e-11, 1e-10);
    }
}

#[test]
fn unequal_layer_temperatures_have_native_directional_derivatives_and_nonlocal_entropy() {
    let source = source(multilayer_inputs(0.0), ReferenceSurface::Midsurface);
    let laws = laws(&source);
    let t = [302.0, 329.0, 314.0];
    let direction = [0.7, -0.4, 1.2];
    let shear = resolved(&source, &laws, &t);
    let d = shear.directional_derivatives(&direction, || false).unwrap();
    assert_eq!(d[0], native_h(&source, &t));
    let h = 0.04;
    let m2 = native_h(&source, &shift(&t, &direction, -2.0 * h));
    let m1 = native_h(&source, &shift(&t, &direction, -h));
    let p1 = native_h(&source, &shift(&t, &direction, h));
    let p2 = native_h(&source, &shift(&t, &direction, 2.0 * h));
    for k in 0..4 {
        let first = (m2[k] - 8.0 * m1[k] + 8.0 * p1[k] - p2[k]) / (12.0 * h);
        let second =
            (-m2[k] + 16.0 * m1[k] - 30.0 * d[0][k] + 16.0 * p1[k] - p2[k]) / (12.0 * h * h);
        close(d[1][k], first, 2e-10, 2e-7);
        close(d[2][k], second, 2e-9, 2e-4);
        let contracted: f64 = shear
            .stiffness_temperature_derivatives_pa_m_k()
            .iter()
            .zip(direction)
            .map(|(g, v)| g[k] * v)
            .sum();
        close(d[1][k], contracted, 1e-13, 2e-11);
    }
    let state = shear.state(GAMMA, AREA, || false).unwrap();
    for i in 0..3 {
        let mut unit = [0.0; 3];
        unit[i] = 1.0;
        let plus = quad(native_h(&source, &shift(&t, &unit, h)), GAMMA, AREA);
        let minus = quad(native_h(&source, &shift(&t, &unit, -h)), GAMMA, AREA);
        close(
            state.entropy_j_k[i],
            -(plus - minus) / (2.0 * h),
            3e-13,
            2e-6,
        );
    }
    let response = shear
        .thermal_response(GAMMA, AREA, &direction, || false)
        .unwrap();
    let plus = resolved(&source, &laws, &shift(&t, &direction, h))
        .state(GAMMA, AREA, || false)
        .unwrap();
    let minus = resolved(&source, &laws, &shift(&t, &direction, -h))
        .state(GAMMA, AREA, || false)
        .unwrap();
    for i in 0..3 {
        close(
            response.entropy_derivative_j_k2[i],
            (plus.entropy_j_k[i] - minus.entropy_j_k[i]) / (2.0 * h),
            1e-13,
            2e-6,
        );
    }
    close(
        response.internal_energy_derivative_j_k,
        (plus.internal_energy_j - minus.internal_energy_j) / (2.0 * h),
        1e-12,
        2e-6,
    );
    let columns: Vec<_> = (0..3)
        .map(|j| {
            let mut unit = [0.0; 3];
            unit[j] = 1.0;
            shear
                .thermal_response(GAMMA, AREA, &unit, || false)
                .unwrap()
                .entropy_derivative_j_k2
        })
        .collect();
    assert!(
        columns[0][1].abs() > 1e-12,
        "heterogeneous stack must retain a mixed thermal derivative"
    );
    for i in 0..3 {
        for j in 0..3 {
            close(columns[i][j], columns[j][i], 1e-16, 2e-10);
        }
    }
}

#[test]
fn reference_translation_and_axis_exchange_preserve_the_native_approximation() {
    let t = [302.0, 329.0, 314.0];
    let direction = [0.7, -0.4, 1.2];
    let inputs = multilayer_inputs(0.0);
    let midsurface = source(inputs.clone(), ReferenceSurface::Midsurface);
    let offset = source(inputs, ReferenceSurface::Offset(0.008));
    let base = resolved(&midsurface, &laws(&midsurface), &t);
    let shifted = resolved(&offset, &laws(&offset), &t);
    let a = base.directional_derivatives(&direction, || false).unwrap();
    let b = shifted
        .directional_derivatives(&direction, || false)
        .unwrap();
    for row in 0..3 {
        for k in 0..4 {
            close(a[row][k], b[row][k], 2e-12, 2e-10);
        }
    }
    // The native directional Whitney model is covariant under x/y exchange.
    // General off-axis multilayers do not claim a coupled 3-D warping solution.
    let rotated = source(
        multilayer_inputs(core::f64::consts::FRAC_PI_2),
        ReferenceSurface::Midsurface,
    );
    let rotated = resolved(&rotated, &laws(&rotated), &t);
    let r = rotated
        .directional_derivatives(&direction, || false)
        .unwrap();
    for row in 0..3 {
        for (k, expected) in [a[row][3], -a[row][2], -a[row][1], a[row][0]]
            .into_iter()
            .enumerate()
        {
            close(r[row][k], expected, 2e-12, 2e-10);
        }
    }
    let original_state = base.state(GAMMA, AREA, || false).unwrap();
    let rotated_state = rotated
        .state([-GAMMA[1], GAMMA[0]], AREA, || false)
        .unwrap();
    close(
        original_state.free_energy_j,
        rotated_state.free_energy_j,
        2e-14,
        1e-10,
    );
    close(
        original_state.internal_energy_j,
        rotated_state.internal_energy_j,
        2e-14,
        1e-10,
    );
    for i in 0..3 {
        close(
            original_state.entropy_j_k[i],
            rotated_state.entropy_j_k[i],
            2e-14,
            1e-10,
        );
    }
}

#[test]
fn source_identity_bounds_cancellation_and_native_fallback_refuse_before_results() {
    let stack = source(multilayer_inputs(0.0), ReferenceSurface::Midsurface);
    let potentials = laws(&stack);
    let t = [302.0, 329.0, 314.0];
    let current = stack.at_temperatures(&t, || false).unwrap();
    assert!(
        ResolvedLaminateShearThermoelastic::try_new(&potentials[..2], &current, || false).is_err()
    );
    let wrong = vec![
        potentials[1].clone(),
        potentials[0].clone(),
        potentials[2].clone(),
    ];
    assert!(ResolvedLaminateShearThermoelastic::try_new(&wrong, &current, || false).is_err());
    let independent_card = source(multilayer_inputs(0.0), ReferenceSurface::Midsurface);
    assert!(
        ResolvedLaminateShearThermoelastic::try_new(&laws(&independent_card), &current, || false)
            .is_err()
    );
    assert!(matches!(
        ResolvedLaminateShearThermoelastic::try_new(&potentials, &current, || true),
        Err(LaminateShearThermoelasticError::Cancelled)
    ));
    let shear = resolved(&stack, &potentials, &t);
    assert!(shear.directional_derivatives(&[1.0], || false).is_err());
    assert!(
        shear
            .directional_derivatives(&[0.0, f64::NAN, 0.0], || false)
            .is_err()
    );
    assert!(shear.state([LIMIT * 2.0, 0.0], AREA, || false).is_err());
    assert!(shear.state(GAMMA, 0.0, || false).is_err());
    let zero = shear.directional_derivatives(&[0.0; 3], || false).unwrap();
    assert_eq!(zero[1], [0.0; 4]);
    assert_eq!(zero[2], [0.0; 4]);
    let mut calls = 0;
    assert!(matches!(
        shear.thermal_response(GAMMA, AREA, &[1.0; 3], || {
            calls += 1;
            calls == 8
        }),
        Err(LaminateShearThermoelasticError::Cancelled)
    ));
    assert_eq!(
        shear.stiffness_pa_m(),
        current.prepared().section().shear_pa_m()
    );
    // An extreme but finite reference thickness activates the native scalar
    // arithmetic's 5/6 fallback (D² overflows). Differentiation must refuse it,
    // even though the source coefficients and native operator are finite.
    let huge = source(vec![input(0, 1e60, 0.0)], ReferenceSurface::Midsurface);
    let huge_laws = laws(&huge);
    let huge_native = huge.at_temperatures(&[300.0], || false).unwrap();
    assert_eq!(
        *huge_native.prepared().section().shear_factor(),
        [5.0 / 6.0; 2]
    );
    assert!(matches!(
        ResolvedLaminateShearThermoelastic::try_new(&huge_laws, &huge_native, || false),
        Err(LaminateShearThermoelasticError::Invalid(_))
    ));
}
