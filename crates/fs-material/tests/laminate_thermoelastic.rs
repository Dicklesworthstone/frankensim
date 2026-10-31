//! Independent compliance and thermodynamic oracles for sourced plane stress.
//! Numerical fixtures only; no experimental material qualification is claimed.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::laminate::thermal::source::{
    LaminaElasticSource, LaminaExpansionSource, LaminaSourceInput, LaminateThermalSource,
};
use fs_material::laminate::thermal::thermoelastic::AffineLaminaThermoelastic;
use fs_material::laminate::{LayerLaw, ReferenceSurface};
use fs_material::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS, IsotropicElasticParameters,
    MaterialPropertySelection, ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES as CTE,
    ScalarAdmissibility, ScalarPropertyRequirement,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use std::sync::Arc;

#[derive(Clone, Copy, Debug)]
enum Kind {
    YoungPoisson,
    BulkShear,
    YoungShear,
    Plane12,
}
const KINDS: [Kind; 4] = [
    Kind::YoungPoisson,
    Kind::BulkShear,
    Kind::YoungShear,
    Kind::Plane12,
];
const RANGE: [f64; 2] = [280.0, 360.0];
const LIMIT: f64 = 0.1;
const V: f64 = 0.007;
const TOTAL: [f64; 5] = [0.014, -0.009, 0.012, 0.006, -0.004];
const ALPHA: [f64; 2] = [1e-4, -4e-5];
const ALPHAT: [f64; 2] = [0.5e-6, 1e-6];

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            t,
        )
        .unwrap()
}
fn coefficients(kind: Kind) -> Vec<(&'static str, Dims, f64, f64)> {
    match kind {
        Kind::YoungPoisson => vec![
            ("young_modulus", Pressure::DIMS, 1200.0, -2.3),
            ("poisson_ratio", Dims::NONE, 0.23, 0.0007),
        ],
        Kind::BulkShear => vec![
            ("measured_bulk", Pressure::DIMS, 1100.0, -1.7),
            ("measured_shear", Pressure::DIMS, 450.0, 0.6),
        ],
        Kind::YoungShear => vec![
            ("measured_young", Pressure::DIMS, 900.0, 0.8),
            ("measured_shear", Pressure::DIMS, 400.0, -0.3),
        ],
        Kind::Plane12 => vec![
            ("young_modulus_1", Pressure::DIMS, 1600.0, -2.0),
            ("young_modulus_2", Pressure::DIMS, 900.0, 0.8),
            ("poisson_ratio_12", Dims::NONE, 0.18, 0.0004),
            ("shear_modulus_12", Pressure::DIMS, 550.0, -0.9),
            ("shear_modulus_31", Pressure::DIMS, 480.0, 0.3),
            ("shear_modulus_23", Pressure::DIMS, 370.0, -0.2),
        ],
    }
}
fn req(name: &str, dims: Dims) -> ScalarPropertyRequirement {
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
}
fn prepare_source(kind: Kind, variable: bool, kink: Option<&str>) -> LaminateThermalSource {
    let mut claims = ClaimSet::new();
    let elastic = coefficients(kind);
    for (name, dims, value, slope) in elastic
        .iter()
        .map(|&(name, dims, value, slope)| (name, dims, value, if variable { slope } else { 0.0 }))
        .chain([
            (CTE[0], INVERSE_TEMPERATURE_DIMS, ALPHA[0], ALPHAT[0]),
            (CTE[1], INVERSE_TEMPERATURE_DIMS, ALPHA[1], ALPHAT[1]),
            (CTE[2], INVERSE_TEMPERATURE_DIMS, 2e-4, 0.0),
            (DENSITY_PROPERTY, Density::DIMS, 900.0, 0.0),
        ])
    {
        let quantity = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        let density = name == DENSITY_PROPERTY;
        let mut knots = vec![
            (280.0, value - 20.0 * slope),
            (300.0, value),
            (360.0, value + 60.0 * slope),
        ];
        if kink == Some(name) {
            knots[1].1 += 0.01 * value.abs().max(1e-10);
        }
        // The reduced potential needs no alpha3 derivative or fictitious E3.
        // Its complete original expansion path is still retained and admitted.
        if name == CTE[2] {
            knots[1].1 *= 1.5;
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if density {
                    PropertyValue::Scalar { value, dims }
                } else {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots,
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
                    quantity,
                    if density { 300.0 } else { 280.0 },
                    if density { 300.0 } else { 360.0 },
                ),
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored plane thermoelastic oracle".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = Arc::new(
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: format!("{kind:?}"),
                phase: "solid".into(),
                process: "numerical derivative fixture".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
    );
    let elasticity = match kind {
        Kind::YoungPoisson => {
            LaminaElasticSource::Isotropic(IsotropicElasticParameters::YoungPoisson {
                young_modulus: req(elastic[0].0, Pressure::DIMS),
                poisson_ratio: req(elastic[1].0, Dims::NONE),
            })
        }
        Kind::BulkShear => LaminaElasticSource::Isotropic(IsotropicElasticParameters::BulkShear {
            bulk_modulus: req(elastic[0].0, Pressure::DIMS),
            shear_modulus: req(elastic[1].0, Pressure::DIMS),
        }),
        Kind::YoungShear => {
            LaminaElasticSource::Isotropic(IsotropicElasticParameters::YoungShear {
                young_modulus: req(elastic[0].0, Pressure::DIMS),
                shear_modulus: req(elastic[1].0, Pressure::DIMS),
            })
        }
        Kind::Plane12 => LaminaElasticSource::OrthotropicPlane12,
    };
    LaminateThermalSource::try_new(
        vec![LaminaSourceInput {
            card,
            reference: point(300.0),
            thickness_m: 0.02,
            orientation_rad: 0.0,
            elasticity,
            expansion: LaminaExpansionSource::Orthotropic,
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        ReferenceSurface::Midsurface,
        || false,
    )
    .unwrap()
}
fn close(a: f64, b: f64, atol: f64, rtol: f64) {
    assert!(
        (a - b).abs() <= atol + rtol * a.abs().max(b.abs()),
        "{a:.16e} != {b:.16e}"
    );
}

// Independent reciprocal/compliance derivatives, NOT the production rational Q
// expressions or its second-order arithmetic.
fn reciprocal(x: [f64; 3]) -> [f64; 3] {
    [
        1.0 / x[0],
        -x[1] / x[0].powi(2),
        2.0 * x[1] * x[1] / x[0].powi(3) - x[2] / x[0].powi(2),
    ]
}
fn ratio(n: f64, nt: f64, d: f64, dt: f64) -> [f64; 3] {
    [
        n / d,
        nt / d - n * dt / d.powi(2),
        2.0 * n * dt * dt / d.powi(3) - 2.0 * nt * dt / d.powi(2),
    ]
}
fn oracle_q(kind: Kind, t: f64) -> ([[f64; 9]; 3], [[f64; 2]; 3]) {
    let c: Vec<_> = coefficients(kind)
        .into_iter()
        .map(|(_, _, v, d)| (v + d * (t - 300.0), d))
        .collect();
    let mut q = [[0.0; 9]; 3];
    let g;
    match kind {
        Kind::Plane12 => {
            let s11 = reciprocal([c[0].0, c[0].1, 0.0]);
            let s22 = reciprocal([c[1].0, c[1].1, 0.0]);
            let s12 = ratio(-c[2].0, -c[2].1, c[0].0, c[0].1);
            let s = core::array::from_fn::<_, 3, _>(|i| [[s11[i], s12[i]], [s12[i], s22[i]]]);
            let determinant = s11[0] * s22[0] - s12[0] * s12[0];
            let a = [
                [s22[0] / determinant, -s12[0] / determinant],
                [-s12[0] / determinant, s11[0] / determinant],
            ];
            let product = |a: [[f64; 2]; 2], b: [[f64; 2]; 2]| {
                core::array::from_fn::<_, 2, _>(|i| {
                    core::array::from_fn::<_, 2, _>(|j| {
                        (0..2).map(|k| a[i][k] * b[k][j]).sum::<f64>()
                    })
                })
            };
            let first = product(product(a, s[1]), a).map(|row| row.map(|x| -x));
            let quadratic = product(product(product(product(a, s[1]), a), s[1]), a);
            let second_term = product(product(a, s[2]), a);
            let second = core::array::from_fn::<_, 2, _>(|i| {
                core::array::from_fn::<_, 2, _>(|j| 2.0 * quadratic[i][j] - second_term[i][j])
            });
            for (order, matrix) in [a, first, second].into_iter().enumerate() {
                q[order][0] = matrix[0][0];
                q[order][1] = matrix[0][1];
                q[order][3] = matrix[1][0];
                q[order][4] = matrix[1][1];
            }
            q[0][8] = c[3].0;
            q[1][8] = c[3].1;
            g = [[c[4].0, c[5].0], [c[4].1, c[5].1], [0.0; 2]];
        }
        _ => {
            let (plus_compliance, minus_compliance) = match kind {
                Kind::YoungPoisson => (
                    ratio(1.0 - c[1].0, -c[1].1, c[0].0, c[0].1),
                    ratio(1.0 + c[1].0, c[1].1, c[0].0, c[0].1),
                ),
                Kind::BulkShear => {
                    let k = reciprocal([c[0].0, c[0].1, 0.0]);
                    let g = reciprocal([c[1].0, c[1].1, 0.0]);
                    (
                        core::array::from_fn(|i| g[i] / 6.0 + 2.0 * k[i] / 9.0),
                        g.map(|x| x / 2.0),
                    )
                }
                Kind::YoungShear => {
                    let e = reciprocal([c[0].0, c[0].1, 0.0]);
                    let g = reciprocal([c[1].0, c[1].1, 0.0]);
                    (
                        core::array::from_fn(|i| 2.0 * e[i] - g[i] / 2.0),
                        g.map(|x| x / 2.0),
                    )
                }
                Kind::Plane12 => unreachable!(),
            };
            let plus = reciprocal(plus_compliance);
            let minus = reciprocal(minus_compliance);
            for order in 0..3 {
                q[order][0] = (plus[order] + minus[order]) / 2.0;
                q[order][4] = q[order][0];
                q[order][1] = (plus[order] - minus[order]) / 2.0;
                q[order][3] = q[order][1];
                q[order][8] = minus[order] / 2.0;
            }
            g = core::array::from_fn(|i| [minus[i] / 2.0; 2]);
        }
    }
    (q, g)
}
fn apply(q: &[f64; 9], g: [f64; 2], e: [f64; 5]) -> [f64; 5] {
    core::array::from_fn(|i| {
        if i < 3 {
            (0..3).map(|j| q[3 * i + j] * e[j]).sum()
        } else {
            g[i - 3] * e[i]
        }
    })
}
fn dot(a: [f64; 5], b: [f64; 5]) -> f64 {
    a.into_iter().zip(b).map(|(x, y)| x * y).sum()
}

#[test]
fn g1_all_four_source_schemas_match_independent_compliance_and_thermodynamic_oracles() {
    for kind in KINDS {
        let source = prepare_source(kind, true, None);
        let law = AffineLaminaThermoelastic::new(&source.layers()[0], RANGE, LIMIT).unwrap();
        assert!(!law.is_temperature_independent_elasticity());
        assert_eq!(law.reference_temperature_k(), 300.0);
        let t = 312.5;
        let native = source.at_temperatures(&[t], || false).unwrap();
        let resolved = law.at_materials(&native.layers()[0]).unwrap();
        assert_eq!(resolved.elasticity(), native.layers()[0].elasticity());
        assert_eq!(resolved.expansion(), native.layers()[0].expansion());
        let (q, g) = oracle_q(kind, t);
        for order in 0..3 {
            for (actual, expected) in resolved.stiffness_derivatives_pa()[order]
                .into_iter()
                .zip(q[order])
            {
                close(actual, expected, 1e-12, 2e-12);
            }
            for (actual, expected) in resolved.transverse_shear_derivatives_pa()[order]
                .into_iter()
                .zip(g[order])
            {
                close(actual, expected, 1e-12, 2e-12);
            }
        }
        match &native.prepared().section().layers()[0].law {
            LayerLaw::PlaneStress { stiffness_pa, .. } => {
                assert_eq!(resolved.stiffness_derivatives_pa()[0], *stiffness_pa)
            }
            other => panic!("unexpected native law {other:?}"),
        }
        let mut e = TOTAL;
        let dt = t - 300.0;
        let alpha = [
            ALPHA[0] + ALPHAT[0] * dt,
            ALPHA[1] + ALPHAT[1] * dt,
            0.0,
            0.0,
            0.0,
        ];
        for i in 0..2 {
            e[i] -= ALPHA[i] * dt + 0.5 * ALPHAT[i] * dt * dt;
        }
        let stress = apply(&q[0], g[0], e);
        let first = apply(&q[1], g[1], e);
        let thermal = apply(&q[0], g[0], alpha);
        let expected_psi = 0.5 * V * dot(e, stress);
        let expected_s = V * (dot(alpha, stress) - 0.5 * dot(e, first));
        let expected_c = -t
            * V
            * (0.5 * dot(e, apply(&q[2], g[2], e)) - 2.0 * dot(alpha, first) + dot(alpha, thermal)
                - dot([ALPHAT[0], ALPHAT[1], 0.0, 0.0, 0.0], stress));
        let state = resolved.state(TOTAL, V).unwrap();
        close(state.free_energy_j, expected_psi, 1e-15, 2e-12);
        close(state.entropy_j_k, expected_s, 1e-15, 2e-12);
        close(
            state.internal_energy_j,
            expected_psi + t * expected_s,
            1e-15,
            2e-12,
        );
        close(state.heat_capacity_j_k, expected_c, 1e-15, 2e-12);
        for i in 0..5 {
            close(state.stress_pa[i], stress[i], 1e-12, 2e-12);
            close(
                state.stress_temperature_derivative_pa_k[i],
                first[i] - thermal[i],
                1e-12,
                2e-12,
            );
        }
        assert!(Arc::ptr_eq(
            &law.source().input().card,
            &source.layers()[0].input().card
        ));
        assert!(
            law.source()
                .reference_density()
                .answer()
                .evidence
                .model
                .validity
                .bounds()["T"]
                == (300.0, 300.0)
        );
    }
}

#[test]
fn g3_engineering_work_maxwell_relation_and_temperature_derivatives_share_one_potential() {
    let source = prepare_source(Kind::YoungPoisson, true, None);
    let law = AffineLaminaThermoelastic::new(&source.layers()[0], RANGE, LIMIT).unwrap();
    let t = 315.0;
    let at = law.at_temperature(t).unwrap();
    let base = at.state(TOTAL, V).unwrap();
    let h = 1e-6;
    for i in 0..5 {
        let mut plus = TOTAL;
        let mut minus = TOTAL;
        plus[i] += h;
        minus[i] -= h;
        let a = at.state(plus, V).unwrap();
        let b = at.state(minus, V).unwrap();
        close(
            (a.free_energy_j - b.free_energy_j) / (2.0 * h),
            V * base.stress_pa[i],
            1e-12,
            1e-8,
        );
        close(
            (a.entropy_j_k - b.entropy_j_k) / (2.0 * h),
            -V * base.stress_temperature_derivative_pa_k[i],
            1e-12,
            1e-8,
        );
    }
    let h = 0.002;
    let plus = law.at_temperature(t + h).unwrap().state(TOTAL, V).unwrap();
    let minus = law.at_temperature(t - h).unwrap().state(TOTAL, V).unwrap();
    close(
        -(plus.free_energy_j - minus.free_energy_j) / (2.0 * h),
        base.entropy_j_k,
        1e-12,
        2e-8,
    );
    close(
        (plus.internal_energy_j - minus.internal_energy_j) / (2.0 * h),
        base.heat_capacity_j_k,
        1e-12,
        2e-8,
    );
    for i in 0..5 {
        close(
            (plus.stress_pa[i] - minus.stress_pa[i]) / (2.0 * h),
            base.stress_temperature_derivative_pa_k[i],
            1e-9,
            2e-8,
        );
    }
    let plane = at
        .in_plane_state([TOTAL[0], TOTAL[1], TOTAL[2]], V)
        .unwrap();
    assert_eq!(
        plane,
        at.state([TOTAL[0], TOTAL[1], TOTAL[2], 0.0, 0.0], V)
            .unwrap()
    );
    assert!(
        base.free_energy_j > plane.free_energy_j,
        "raw transverse shear is distinct from in-plane storage"
    );
}

#[test]
fn g1_fixed_elastic_certificate_and_whole_interval_capacity_bound_preserve_signed_thermodynamics() {
    for kind in KINDS {
        for variable in [false, true] {
            let source = prepare_source(kind, variable, None);
            let law = AffineLaminaThermoelastic::new(&source.layers()[0], RANGE, LIMIT).unwrap();
            assert_eq!(law.is_temperature_independent_elasticity(), !variable);
            let lower = law.minimum_heat_capacity_j_k(V).unwrap();
            assert!(lower < 0.0 && lower.is_finite());
            for t in [280.0, 303.0, 327.0, 360.0] {
                let state = law.at_temperature(t).unwrap();
                for signs in 0..32 {
                    let strain =
                        core::array::from_fn(|i| if signs & (1 << i) == 0 { -0.04 } else { 0.04 });
                    assert!(state.state(strain, V).unwrap().heat_capacity_j_k >= lower);
                }
            }
        }
    }
    let source = prepare_source(Kind::Plane12, false, None);
    let law = AffineLaminaThermoelastic::new(&source.layers()[0], RANGE, LIMIT).unwrap();
    let at = law.at_temperature(320.0).unwrap();
    let response = at.in_plane_state(at.free_strain(), V).unwrap();
    close(response.free_energy_j, 0.0, 0.0, 0.0);
    close(response.entropy_j_k, 0.0, 0.0, 0.0);
    close(response.internal_energy_j, 0.0, 0.0, 0.0);
    assert!(
        response.heat_capacity_j_k < 0.0,
        "fixed-total-strain capacity is not the derivative along a moving stress-free path"
    );
    assert!(
        source.layers()[0]
            .input()
            .card
            .claims_for("young_modulus_3")
            .is_empty()
    );
    assert!(
        source.layers()[0]
            .input()
            .card
            .claims_for("poisson_ratio_13")
            .is_empty()
    );
}

#[test]
fn g0_source_kinks_reference_changes_and_strain_temperature_domains_refuse() {
    for kink in ["young_modulus", CTE[0]] {
        let source = prepare_source(Kind::YoungPoisson, true, Some(kink));
        assert!(AffineLaminaThermoelastic::new(&source.layers()[0], RANGE, LIMIT).is_err());
    }
    let source = prepare_source(Kind::Plane12, true, None);
    let law = AffineLaminaThermoelastic::new(&source.layers()[0], RANGE, LIMIT).unwrap();
    assert!(AffineLaminaThermoelastic::new(&source.layers()[0], [301.0, 360.0], LIMIT).is_err());
    for t in [0.0, 279.0, 361.0, f64::NAN] {
        assert!(law.at_temperature(t).is_err());
    }
    let at = law.at_temperature(320.0).unwrap();
    assert!(at.state([0.2; 5], V).is_err());
    assert!(at.state(TOTAL, 0.0).is_err());
    assert!(law.minimum_heat_capacity_j_k(f64::INFINITY).is_err());
    let mut input = source.layers()[0].input().clone();
    // Density is intentionally only supported at the original reference.
    // A new zero-expansion reference cannot replace that fixed source inventory.
    input.reference = point(310.0);
    assert!(
        LaminateThermalSource::try_new(vec![input], ReferenceSurface::Midsurface, || false)
            .is_err()
    );
    let foreign = prepare_source(Kind::YoungPoisson, true, None);
    let foreign = foreign.at_temperatures(&[320.0], || false).unwrap();
    assert!(law.at_materials(&foreign.layers()[0]).is_err());
}
