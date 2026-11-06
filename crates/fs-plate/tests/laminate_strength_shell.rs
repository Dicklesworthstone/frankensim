//! Analytical source-bound shell-strength regressions, not material calibration.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::laminate::ReferenceSurface;
use fs_material::laminate::strength::tsai_wu::{
    PreparedLaminateTsaiWu, TsaiWuInteractionSelection,
};
use fs_material::laminate::strength::{
    LaminaStrengthProperties, LaminaStrengthSelection, PreparedLaminateStrength,
};
use fs_material::laminate::thermal::source::{
    LaminaElasticSource, LaminaExpansionSource, LaminaSourceInput, LaminateThermalSource,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as CTE, MaterialPropertySelection,
};
use fs_plate::PlateError;
use fs_plate::shell::ShellMesh;
use fs_plate::shell::mitc::thermal::LaminateThermalField;
use fs_plate::shell::mitc::thermal::strength::{
    LaminateStrengthCriterion as Criterion, ShellStrengthError,
};
use fs_plate::shell::mitc::{
    MitcBudget, MitcShellModel, NodalDirectors, ShellMass, assemble_mitc3p,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec};
use std::sync::Arc;

const H: f64 = 0.02;
const E1: f64 = 100.0e9;
const E2: f64 = 20.0e9;
const G12: f64 = 10.0e9;
const NAMES: [&str; 5] = ["Xt", "Xc", "Yt", "Yc", "S12"];
const LIMITS: [f64; 5] = [120.0e6, 300.0e6, 30.0e6, 90.0e6, 15.0e6];

fn temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn card(alpha: f64) -> Arc<MaterialCard> {
    let values = [
        (DENSITY_PROPERTY, Density::DIMS, 1600.0),
        ("young_modulus_1", Pressure::DIMS, E1),
        ("young_modulus_2", Pressure::DIMS, E2),
        ("poisson_ratio_12", Dims::NONE, 0.0),
        ("shear_modulus_12", Pressure::DIMS, G12),
        ("shear_modulus_31", Pressure::DIMS, 12.0e9),
        ("shear_modulus_23", Pressure::DIMS, 8.0e9),
        (CTE, INVERSE_TEMPERATURE_DIMS, alpha),
        ("interaction", Dims::NONE, 0.2),
    ];
    let mut claims = ClaimSet::new();
    for (name, dims, value) in values.into_iter().chain(
        NAMES
            .into_iter()
            .zip(LIMITS)
            .map(|(name, value)| (name, Pressure::DIMS, value)),
    ) {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    temperature(),
                    280.0,
                    400.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored shell strength oracle".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    Arc::new(
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic orthotropic shell".into(),
                phase: "solid".into(),
                process: "numerical strength fixture".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
    )
}

fn source(alphas: &[f64], angle: f64, current_k: f64) -> PreparedLaminateStrength {
    let layers = alphas
        .iter()
        .map(|&alpha| LaminaSourceInput {
            card: card(alpha),
            reference: QueryPoint::new()
                .with_quantity("T", temperature(), 300.0)
                .unwrap(),
            thickness_m: H / alphas.len() as f64,
            orientation_rad: angle,
            elasticity: LaminaElasticSource::OrthotropicPlane12,
            expansion: LaminaExpansionSource::Isotropic,
            selection: MaterialPropertySelection::SingleClaimOnly,
        })
        .collect();
    let thermal = LaminateThermalSource::try_new(layers, ReferenceSurface::Midsurface, || false)
        .unwrap()
        .at_temperatures(&vec![current_k; alphas.len()], || false)
        .unwrap();
    let selection = LaminaStrengthSelection {
        properties: LaminaStrengthProperties {
            tensile: ["Xt".into(), "Yt".into()],
            compressive: ["Xc".into(), "Yc".into()],
            shear: "S12".into(),
        },
        selection: MaterialPropertySelection::SingleClaimOnly,
    };
    PreparedLaminateStrength::resolve(thermal, &vec![selection; alphas.len()], || false).unwrap()
}

fn coupled(source: PreparedLaminateStrength) -> PreparedLaminateTsaiWu {
    let n = source.thermal().layers().len();
    PreparedLaminateTsaiWu::from_maximum_stress(
        source,
        &vec![
            TsaiWuInteractionSelection {
                property: "interaction".into(),
                selection: MaterialPropertySelection::SingleClaimOnly,
            };
            n
        ],
        || false,
    )
    .unwrap()
}

fn mesh() -> ShellMesh {
    ShellMesh::new(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap()
}

fn model(field: &LaminateThermalField) -> MitcShellModel {
    let mesh = mesh();
    let directors = NodalDirectors::new(&mesh, vec![[0.0, 0.0, 1.0]; 3], 0.9).unwrap();
    assemble_mitc3p(
        &mesh,
        &directors,
        &[field.section().clone()],
        &[],
        ShellMass::Consistent,
        MitcBudget {
            max_nodes: 3,
            max_triangles: 1,
        },
        1.0e-3,
        || false,
    )
    .unwrap()
}

fn rotate(e: [f64; 3], angle: f64) -> [f64; 3] {
    let (s, c) = angle.sin_cos();
    [
        c * c * e[0] + s * s * e[1] - c * s * e[2],
        s * s * e[0] + c * c * e[1] + c * s * e[2],
        2.0 * c * s * (e[0] - e[1]) + (c * c - s * s) * e[2],
    ]
}

fn displacement(e: [f64; 3], k: [f64; 3]) -> [f64; 18] {
    let mut u = [0.0; 18];
    for (node, [x, y, _]) in mesh().nodes.into_iter().enumerate() {
        u[6 * node] = e[0] * x + 0.5 * e[2] * y;
        u[6 * node + 1] = 0.5 * e[2] * x + e[1] * y;
        u[6 * node + 2] = -0.5 * (k[0] * x * x + k[1] * y * y + k[2] * x * y);
        u[6 * node + 3] = -k[1] * y - 0.5 * k[2] * x;
        u[6 * node + 4] = k[0] * x + 0.5 * k[2] * y;
    }
    u
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2.0e-8 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn rotated_native_membrane_uses_material_strength_axes_and_interaction() {
    let strain = [0.0005, 0.0008, -0.0002];
    let stress = [E1 * strain[0], E2 * strain[1], G12 * strain[2]];
    for (layer_angle, section_angle) in [(0.0, 0.0), (0.31, -0.53), (-0.61, 0.77)] {
        let strength = source(&[0.0], layer_angle, 300.0);
        let interaction = coupled(strength.clone());
        let field =
            LaminateThermalField::try_new(Arc::clone(strength.thermal().prepared()), section_angle)
                .unwrap();
        let model = model(&field);
        let shell = model.prepare_laminate_thermal(&[field], || false).unwrap();
        let u = displacement(rotate(strain, layer_angle + section_angle), [0.0; 3]);
        let independent = strength.layer_strength(0).unwrap().assess(stress).unwrap();
        let a = shell
            .laminate_strength_bound(0, Criterion::MaximumStress(&strength), &u, 12, || false)
            .unwrap();
        near(a.utilization_upper_bound(), independent.utilization());
        // Independent positive root of A r^2 + B r = 1, expressed as utilization.
        let [xt, xc, yt, yc, shear] = LIMITS;
        let [s1, s2, t] = stress;
        let quadratic = s1 * s1 / (xt * xc)
            + s2 * s2 / (yt * yc)
            + t * t / (shear * shear)
            + 0.4 * s1 * s2 / (xt * xc * yt * yc).sqrt();
        let linear = s1 * (1.0 / xt - 1.0 / xc) + s2 * (1.0 / yt - 1.0 / yc);
        let expected = 0.5 * (linear + (linear * linear + 4.0 * quadratic).sqrt());
        let b = shell
            .laminate_strength_bound(0, Criterion::TsaiWu(&interaction), &u, 12, || false)
            .unwrap();
        near(b.utilization_upper_bound(), expected);
        assert!((a.utilization_upper_bound() - b.utilization_upper_bound()).abs() > 0.01);
    }
}

#[test]
fn bending_and_thermal_self_stress_survive_zero_section_resultants() {
    let delta = 0.002;
    let strength = source(&[0.0, 4.0e-5], 0.0, 350.0);
    let field =
        LaminateThermalField::try_new(Arc::clone(strength.thermal().prepared()), 0.0).unwrap();
    let model = model(&field);
    let shell = model.prepare_laminate_thermal(&[field], || false).unwrap();
    let u = displacement(
        [delta / 2.0, delta / 2.0, 0.0],
        [1.5 * delta / H, 1.5 * delta / H, 0.0],
    );
    let response = shell.recover(&u, || false).unwrap();
    for value in response.facets[0]
        .membrane_force_n_m
        .into_iter()
        .chain(response.facets[0].bending_moment_n)
    {
        assert!(
            value.abs() < 1.0e-5,
            "self-equilibrated resultants: {value:e}"
        );
    }
    let bound = shell
        .laminate_strength_bound(0, Criterion::MaximumStress(&strength), &u, 24, || false)
        .unwrap();
    near(
        bound.utilization_upper_bound(),
        E1 * delta / (2.0 * LIMITS[0]),
    );
    assert!(bound.utilization_upper_bound() > 0.8);
    assert!(response.strain_energy_j > 0.0);
}

#[test]
fn whole_facet_envelope_covers_nonuniform_condensed_field_with_thermal_stress() {
    let strength = source(&[1.0e-5, -2.0e-5], 0.37, 340.0);
    let interaction = coupled(strength.clone());
    let field =
        LaminateThermalField::try_new(Arc::clone(strength.thermal().prepared()), -0.21).unwrap();
    let model = model(&field);
    let shell = model.prepare_laminate_thermal(&[field], || false).unwrap();
    let u = core::array::from_fn::<_, 18, _>(|i| 0.0002 * ((i * 7 % 13) as f64 - 6.0));
    for criterion in [
        Criterion::MaximumStress(&strength),
        Criterion::TsaiWu(&interaction),
    ] {
        let bound = shell
            .laminate_strength_bound(0, criterion, &u, 24, || false)
            .unwrap();
        let mut largest = 0.0_f64;
        for i in 0..=12 {
            for j in 0..=12 - i {
                let point = [i as f64 / 12.0, j as f64 / 12.0, (12 - i - j) as f64 / 12.0];
                let actual = shell
                    .laminate_strength_at(0, criterion, &u, point, 4, || false)
                    .unwrap();
                largest = largest.max(actual.utilization());
                assert!(actual.utilization() <= bound.utilization_upper_bound() * (1.0 + 2.0e-12));
            }
        }
        assert!(largest > 0.0);
        assert!(bound.utilization_upper_bound() >= largest);
    }
}

#[test]
fn source_temperature_budget_nonfinite_and_cancellation_are_admitted_before_output() {
    let strength = source(&[1.0e-5], 0.0, 320.0);
    let field =
        LaminateThermalField::try_new(Arc::clone(strength.thermal().prepared()), 0.0).unwrap();
    let model = model(&field);
    let shell = model.prepare_laminate_thermal(&[field], || false).unwrap();
    let criterion = Criterion::MaximumStress(&strength);
    for other in [
        source(&[1.0e-5], 0.0, 320.0),
        strength.at_temperatures(&[330.0], || false).unwrap(),
    ] {
        assert!(matches!(
            shell.laminate_strength_bound(
                0,
                Criterion::MaximumStress(&other),
                &[0.0; 18],
                12,
                || false
            ),
            Err(ShellStrengthError::SourceMismatch)
        ));
    }
    assert!(matches!(
        shell.laminate_strength_bound(0, criterion, &[0.0; 18], 11, || false),
        Err(ShellStrengthError::FaceBudget {
            required: 12,
            maximum: 11
        })
    ));
    assert!(
        shell
            .laminate_strength_bound(0, criterion, &[0.0; 17], 12, || false)
            .is_err()
    );
    assert!(
        shell
            .laminate_strength_bound(0, criterion, &[f64::NAN; 18], 12, || false)
            .is_err()
    );
    assert!(
        shell
            .laminate_strength_at(0, criterion, &[0.0; 18], [-0.1, 0.5, 0.6], 2, || false)
            .is_err()
    );
    assert!(matches!(
        shell.laminate_strength_bound(0, criterion, &[0.0; 18], 12, || true),
        Err(ShellStrengthError::Native(PlateError::Cancelled))
    ));
    let mut calls = 0;
    let before = shell
        .laminate_strength_bound(0, criterion, &[0.0; 18], 12, || {
            calls += 1;
            false
        })
        .unwrap();
    let mut visited = 0;
    assert!(matches!(
        shell.laminate_strength_bound(0, criterion, &[0.0; 18], 12, || {
            visited += 1;
            visited == calls
        }),
        Err(ShellStrengthError::Native(PlateError::Cancelled))
    ));
    let after = shell
        .laminate_strength_bound(0, criterion, &[0.0; 18], 12, || false)
        .unwrap();
    assert_eq!(before, after);
}
