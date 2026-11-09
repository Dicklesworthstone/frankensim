//! G0/G1/G3: fixed-state and equilibrium derivatives of native thermal stress.
//! Authored temperature-dependent tensors verify mechanics, not calibration.

use std::ops::ControlFlow;

use fs_blake3::ContentHash;
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::thermal::{ThermalJacobianLimits3, ThermoelasticResponse3};
use fs_cutfem::elastic3::{
    CutElasticity3, ElasticMaterialField3, ElasticityError3, ElasticityFormulation3,
    ElasticityOptions3, MaterialRefreshLimits3,
};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_evidence::ValidityDomain;
use fs_ivl::Interval;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::elastic::thermal::{
    AffineThermoelastic, IsothermalThermoelasticInput, IsothermalThermoelasticState,
    ThermoelasticSymmetry,
};
use fs_material::state_point::{
    MaterialPropertySelection, ORTHOTROPIC_POISSON_RATIO_PROPERTIES,
    ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES, ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES,
    ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES, POISSON_RATIO_PROPERTY, ScalarAdmissibility,
    ScalarPropertyRequirement, YOUNG_MODULUS_PROPERTY, resolve_material_state_point,
};
use fs_qty::{
    Dims, Pressure, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_solver::op::LinearOp;

const POINTS: usize = 100_000;
const TOTAL_STRAIN: f64 = 0.1;
const EXPANSION_DIMS: Dims = Dims([0, 0, 0, -1, 0, 0]);

fn poll() -> ControlFlow<()> {
    ControlFlow::Continue(())
}

struct Top;
impl CutSdf3 for Top {
    fn value(&self, p: [f64; 3]) -> f64 {
        p[2] - 0.73
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2] - 0.73, hi[2] - 0.73)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        let value = if axis == HeightAxis::Z { 1.0 } else { 0.0 };
        Interval::new(value, value)
    }
}

fn control<T>(action: impl FnOnce(&mut QuadratureControl3<'_>) -> T) -> T {
    let mut callback = |_| poll();
    let mut control = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 0,
            ..Default::default()
        },
        &mut callback,
    )
    .unwrap();
    action(&mut control)
}

fn domain() -> HexCell {
    HexCell::try_new([0.0; 3], [1.0; 3]).unwrap()
}

fn law(isotropic: bool, nearly_incompressible: bool) -> AffineThermoelastic {
    let quantity = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let mut claims = ClaimSet::new();
    let mut insert = |name: &str, dims: Dims, base: f64, slope: f64| {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![
                        (280.0, base - 20.0 * slope),
                        (300.0, base),
                        (360.0, base + 60.0 * slope),
                    ],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", quantity, 280.0, 360.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "authored native thermal stress derivative fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    let mut elastic_names = Vec::new();
    if isotropic {
        insert(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1200.0, -0.7);
        insert(
            POISSON_RATIO_PROPERTY,
            Dims::NONE,
            if nearly_incompressible { 0.48 } else { 0.2 },
            0.0001,
        );
        elastic_names.extend([
            (YOUNG_MODULUS_PROPERTY, Pressure::DIMS),
            (POISSON_RATIO_PROPERTY, Dims::NONE),
        ]);
    } else {
        for i in 0..3 {
            for (name, dims, base, slope) in [
                (
                    ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES[i],
                    Pressure::DIMS,
                    [1200.0, 900.0, 650.0][i],
                    [-0.7, 0.4, -0.2][i],
                ),
                (
                    ORTHOTROPIC_POISSON_RATIO_PROPERTIES[i],
                    Dims::NONE,
                    [0.16, 0.12, 0.2][i],
                    [0.0006, -0.0004, 0.0005][i],
                ),
                (
                    ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES[i],
                    Pressure::DIMS,
                    [380.0, 260.0, 310.0][i],
                    [-0.12, 0.08, 0.05][i],
                ),
            ] {
                insert(name, dims, base, slope);
                elastic_names.push((name, dims));
            }
        }
    }
    for (i, name) in ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES
        .into_iter()
        .enumerate()
    {
        insert(
            name,
            EXPANSION_DIMS,
            [0.0001, -0.00005, 0.0002][i],
            [0.000001, -0.0000002, -0.0000008][i],
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic oriented thermal solid".into(),
            phase: "solid".into(),
            process: "stress adjoint verification".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity("T", quantity, 300.0)
        .unwrap();
    let requirements = |names: Vec<(&str, Dims)>| {
        names
            .into_iter()
            .map(|(name, dims)| {
                ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap()
            })
            .collect::<Vec<_>>()
    };
    let elastic = resolve_material_state_point(
        &card,
        &point,
        &requirements(elastic_names),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let expansion = resolve_material_state_point(
        &card,
        &point,
        &requirements(
            ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES
                .map(|name| (name, EXPANSION_DIMS))
                .to_vec(),
        ),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    AffineThermoelastic::new(
        &card,
        &elastic,
        &expansion,
        if isotropic {
            ThermoelasticSymmetry::Isotropic
        } else {
            ThermoelasticSymmetry::OrientedOrthotropic {
                principal_to_target: [[0.8, -0.6, 0.0], [0.6, 0.8, 0.0], [0.0, 0.0, 1.0]],
            }
        },
        ThermoelasticSymmetry::OrientedOrthotropic {
            principal_to_target: [[0.8, 0.0, 0.6], [0.0, 1.0, 0.0], [-0.6, 0.0, 0.8]],
        },
        [280.0, 360.0],
        0.05,
    )
    .unwrap()
}

fn hot(law: &AffineThermoelastic, temperature: f64) -> IsothermalThermoelasticState {
    law.at_temperature(temperature)
        .unwrap()
        .isothermal_elasticity(IsothermalThermoelasticInput {
            reference_density_kg_m3: 3.0,
            mechanical_material_identity: ContentHash([89; 32]),
        })
        .unwrap()
}

fn temperature(cell: HexCell) -> f64 {
    let [x, y, z] = cell.lo();
    316.0 + 3.0 * x + 2.0 * y + z
}

fn cartesian(law: &AffineThermoelastic, formulation: ElasticityFormulation3) -> CutElasticity3 {
    let materials = |cell| Ok(hot(law, temperature(cell)));
    control(|control| {
        CutElasticity3::build_with_material_field_and_embedded_dirichlet(
            domain(),
            [2; 3],
            &Top,
            ElasticMaterialField3::ThermoelasticByCell(&materials),
            &|p| p[0] == 0.0,
            &|_, normal| normal[2] > 0.5,
            ElasticityOptions3 {
                formulation,
                ghost_gamma: 2.7,
                ..Default::default()
            },
            Default::default(),
            Default::default(),
            control,
        )
        .unwrap()
    })
}

fn adaptive(law: &AffineThermoelastic, formulation: ElasticityFormulation3) -> AdaptiveElasticity3 {
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let tree = tree
        .refined(&[*tree.leaves().iter().next().unwrap()], poll)
        .unwrap();
    let materials = |cell| Ok(hot(law, temperature(cell)));
    control(|control| {
        AdaptiveElasticity3::build_with_material_field_and_embedded_dirichlet(
            domain(),
            &tree,
            &Top,
            ElasticMaterialField3::ThermoelasticByCell(&materials),
            &|p| p[0] == 0.0,
            &|_, normal| normal[2] > 0.5,
            ElasticityOptions3 {
                formulation,
                ghost_gamma: 2.7,
                ..Default::default()
            },
            Default::default(),
            Default::default(),
            control,
        )
        .unwrap()
    })
}

fn fields(nodes: &[[f64; 3]]) -> (Vec<f64>, Vec<f64>) {
    let u = nodes
        .iter()
        .flat_map(|&[x, y, z]| {
            [
                0.001 * x * x * (1.0 + y),
                -0.0008 * x * (z * z + 0.2 * y),
                0.0006 * x * y * y,
            ]
        })
        .collect();
    // The algebraic direction itself lies outside the physical strain domain.
    let direction = nodes
        .iter()
        .flat_map(|&[x, y, z]| [1.3 * x * (1.0 + y), -0.7 * x * z, 0.9 * x * y])
        .collect();
    (u, direction)
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

fn close(a: f64, b: f64, relative: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a:e} != {b:e}",
    );
}

fn cotangents(points: &[ThermoelasticResponse3]) -> Vec<[f64; 6]> {
    points
        .iter()
        .map(|point| {
            let p = &point.mechanical;
            let [x, y, z] = p.position;
            [1.0 + x, -0.4 + z, 0.6, 0.3 + 0.2 * y, -0.8, 0.9].map(|value| p.weight * value)
        })
        .collect()
}

fn objective(points: &[ThermoelasticResponse3], derivatives: &[[f64; 6]]) -> f64 {
    assert_eq!(points.len(), derivatives.len());
    points
        .iter()
        .zip(derivatives)
        .map(|(point, d)| dot(&point.mechanical.stress, d))
        .sum()
}

// Both native backgrounds expose these same operations. Expand the identical
// numerical experiment without adding a test-only operator abstraction.
macro_rules! fixed_partials {
    ($op:ident, $law:ident, $lower:ident) => {{
        let scales: Vec<_> = (0..$op.cells()).map(|i| 0.6 + 0.3 * i as f64 / $op.cells() as f64).collect();
        $op.set_scales(&scales).unwrap();
        let (u, direction) = fields($op.nodes());
        let points = $op
            .thermoelastic_response(&u, TOTAL_STRAIN, POINTS, poll)
            .unwrap();
        assert!(
            points
                .iter()
                .any(|p| p.mechanical.free_strain[5].abs() > 1e-5)
        );
        let derivatives = cotangents(&points);
        let partial = $op
            .thermoelastic_stress_pullback(&u, TOTAL_STRAIN, &derivatives, poll)
            .unwrap();
        assert_eq!(partial.displacement.len(), $op.n());
        assert_eq!(partial.scales.len(), $op.cells());
        assert_eq!(partial.temperatures.len(), $op.cells());
        let h = 1e-6;
        let values = [-1.0, 1.0].map(|sign| {
            let displaced: Vec<_> = u
                .iter()
                .zip(&direction)
                .map(|(u, d)| u + sign * h * d)
                .collect();
            objective(
                &$op.thermoelastic_response(&displaced, TOTAL_STRAIN, POINTS, poll)
                    .unwrap(),
                &derivatives,
            )
        });
        close(
            dot(&partial.displacement, &direction),
            (values[1] - values[0]) / (2.0 * h),
            2e-8,
            2e-8,
        );
        for cell in 0..$op.cells() {
            let h = 2e-5;
            let values = [-1.0, 1.0].map(|sign| {
                let mut changed = scales.clone();
                changed[cell] += sign * h;
                $op.set_scales(&changed).unwrap();
                objective(
                    &$op.thermoelastic_response(&u, TOTAL_STRAIN, POINTS, poll)
                        .unwrap(),
                    &derivatives,
                )
            });
            close(
                partial.scales[cell],
                (values[1] - values[0]) / (2.0 * h),
                2e-8,
                2e-8,
            );
        }
        $op.set_scales(&scales).unwrap();
        for (cell, target) in $lower.iter().enumerate() {
            let h = 0.002;
            let values = [-1.0, 1.0].map(|sign| {
                let materials = |bounds: HexCell| {
                    let delta = if bounds.lo() == *target {
                        sign * h
                    } else {
                        0.0
                    };
                    Ok(hot(&$law, temperature(bounds) + delta))
                };
                $op.refresh_material_field(
                    ElasticMaterialField3::ThermoelasticByCell(&materials),
                    MaterialRefreshLimits3::default(),
                    poll,
                )
                .unwrap();
                objective(
                    &$op.thermoelastic_response(&u, TOTAL_STRAIN, POINTS, poll)
                        .unwrap(),
                    &derivatives,
                )
            });
            close(
                partial.temperatures[cell],
                (values[1] - values[0]) / (2.0 * h),
                3e-6,
                2e-8,
            );
        }
    }};
}

#[test]
fn cartesian_stress_partials_keep_native_tensor_temperature_and_scale_dependence() {
    for (isotropic, formulation) in [
        (false, ElasticityFormulation3::Displacement),
        (true, ElasticityFormulation3::Displacement),
        (true, ElasticityFormulation3::BBar),
    ] {
        let law = law(isotropic, formulation == ElasticityFormulation3::BBar);
        let mut op = cartesian(&law, formulation);
        let lower: Vec<_> = op.cell_nodes().iter().map(|n| op.nodes()[n[0]]).collect();
        fixed_partials!(op, law, lower);
    }
}

#[test]
fn adaptive_stress_partials_use_original_hanging_reconstruction_and_transpose() {
    for (isotropic, formulation) in [
        (false, ElasticityFormulation3::Displacement),
        (true, ElasticityFormulation3::Displacement),
        (true, ElasticityFormulation3::BBar),
    ] {
        let law = law(isotropic, formulation == ElasticityFormulation3::BBar);
        let mut op = adaptive(&law, formulation);
        assert!(op.physical_nodes().len() > op.nodes().len());
        let lower: Vec<_> = op
            .cell_nodes()
            .iter()
            .map(|n| op.physical_nodes()[n[0]])
            .collect();
        fixed_partials!(op, law, lower);
    }
}

fn stress_square(points: &[ThermoelasticResponse3]) -> (f64, Vec<[f64; 6]>) {
    let mut value = 0.0;
    let mut derivatives = Vec::with_capacity(points.len());
    for point in points {
        let p = &point.mechanical;
        // The objective, not the pullback, owns tensor shear multiplicity.
        let d = std::array::from_fn(|i| p.weight * p.stress[i] * if i < 3 { 1.0 } else { 2.0 });
        value += 0.5 * dot(&p.stress, &d);
        derivatives.push(d);
    }
    (value, derivatives)
}

macro_rules! equilibrium_gradient {
    ($op:ident, $type:ty, $law:ident, $lower:ident) => {{
        let scales: Vec<_> = (0..$op.cells()).map(|i| 0.65 + 0.2 * i as f64 / $op.cells() as f64).collect();
        $op.set_scales(&scales).unwrap();
        let external = $op.body_load(&|_| [0.03, -0.01, 0.02], poll).unwrap();
        let solve = |op: &$type| {
            let prepared = op.prepare_thermoelastic(POINTS, poll).unwrap();
            let free = op.eigenstrain_load(&prepared, poll).unwrap();
            let rhs: Vec<_> = external.iter().zip(free).map(|(f, b)| f + b).collect();
            op.solve_controlled(&rhs, 2e-12, 30_000, 32, |_| poll())
                .unwrap()
        };
        let state = solve(&$op);
        let u = state.coefficients();
        let points = $op
            .thermoelastic_response(u, TOTAL_STRAIN, POINTS, poll)
            .unwrap();
        let (_, derivatives) = stress_square(&points);
        let partial = $op
            .thermoelastic_stress_pullback(u, TOTAL_STRAIN, &derivatives, poll)
            .unwrap();
        let adjoint = $op
            .solve_controlled(&partial.displacement, 2e-12, 30_000, 32, |_| poll())
            .unwrap();
        let residual = $op
            .thermal_residual_temperature_pullback(
                u,
                adjoint.coefficients(),
                TOTAL_STRAIN,
                POINTS,
                poll,
            )
            .unwrap();
        let gradient: Vec<_> = partial
            .temperatures
            .iter()
            .zip(&residual)
            .map(|(direct, implicit)| direct - implicit)
            .collect();
        let direction: Vec<_> = (0..$op.cells())
            .map(|i| 0.4 - 0.17 * (i % 5) as f64)
            .collect();
        let expected = dot(&gradient, &direction);
        // This fixture must exercise the equilibrium term, not merely the
        // already-tested explicit source partial at a fixed displacement.
        assert!(dot(&residual, &direction).abs() > 1e-6);
        let h = 0.002;
        let values = [-1.0, 1.0].map(|sign| {
            let materials = |bounds: HexCell| {
                let cell = $lower.iter().position(|lo| *lo == bounds.lo()).unwrap();
                Ok(hot(&$law, temperature(bounds) + sign * h * direction[cell]))
            };
            $op.refresh_material_field(
                ElasticMaterialField3::ThermoelasticByCell(&materials),
                MaterialRefreshLimits3::default(),
                poll,
            )
            .unwrap();
            let perturbed = solve(&$op);
            let response = $op
                .thermoelastic_response(perturbed.coefficients(), TOTAL_STRAIN, POINTS, poll)
                .unwrap();
            stress_square(&response).0
        });
        close(expected, (values[1] - values[0]) / (2.0 * h), 8e-5, 2e-7);
    }};
}

#[test]
fn solved_cartesian_stress_gradient_combines_direct_temperature_with_full_residual_adjoint() {
    for (isotropic, formulation) in [
        (false, ElasticityFormulation3::Displacement),
        (true, ElasticityFormulation3::BBar),
    ] {
        let law = law(isotropic, formulation == ElasticityFormulation3::BBar);
        let mut op = cartesian(&law, formulation);
        let lower: Vec<_> = op.cell_nodes().iter().map(|n| op.nodes()[n[0]]).collect();
        equilibrium_gradient!(op, CutElasticity3, law, lower);
    }
}

#[test]
fn solved_adaptive_stress_gradient_preserves_cell_temperatures_ghosts_and_supports() {
    for (isotropic, formulation) in [
        (false, ElasticityFormulation3::Displacement),
        (true, ElasticityFormulation3::BBar),
    ] {
        let law = law(isotropic, formulation == ElasticityFormulation3::BBar);
        let mut op = adaptive(&law, formulation);
        let lower: Vec<_> = op
            .cell_nodes()
            .iter()
            .map(|n| op.physical_nodes()[n[0]])
            .collect();
        equilibrium_gradient!(op, AdaptiveElasticity3, law, lower);
    }
}

#[test]
fn prepared_temperature_jacobian_keeps_oriented_tensor_and_expansion_derivatives() {
    let law = law(false, false);
    let mut op = cartesian(&law, ElasticityFormulation3::Displacement);
    let scales: Vec<_> = (0..op.cells()).map(|i| 0.7 + 0.03 * i as f64).collect();
    op.set_scales(&scales).unwrap();
    let (u, _) = fields(op.nodes());
    let lower: Vec<_> = op.cell_nodes().iter().map(|n| op.nodes()[n[0]]).collect();
    let direction: Vec<_> = (0..op.cells())
        .map(|i| 0.4 - 0.17 * (i % 5) as f64)
        .collect();
    let prepared = op
        .prepare_thermal_residual_jacobian(
            &u,
            TOTAL_STRAIN,
            ThermalJacobianLimits3 {
                max_points: POINTS,
                max_workspace_components: 1_000_000,
                max_action_terms: 1_000_000,
            },
            poll,
        )
        .unwrap();
    let mut tangent = vec![0.0; op.n()];
    op.thermal_residual_temperature_apply_into(&prepared, &direction, &mut tangent, poll)
        .unwrap();
    let h = 0.002;
    let values = [-1.0, 1.0].map(|sign| {
        let materials = |bounds: HexCell| {
            let cell = lower.iter().position(|lo| *lo == bounds.lo()).unwrap();
            Ok(hot(&law, temperature(bounds) + sign * h * direction[cell]))
        };
        op.refresh_material_field(
            ElasticMaterialField3::ThermoelasticByCell(&materials),
            MaterialRefreshLimits3::default(),
            poll,
        )
        .unwrap();
        let mut residual = vec![0.0; op.n()];
        op.apply(&u, &mut residual);
        let source = op.prepare_thermoelastic(POINTS, poll).unwrap();
        for (value, force) in residual
            .iter_mut()
            .zip(op.eigenstrain_load(&source, poll).unwrap())
        {
            *value -= force;
        }
        residual
    });
    for (i, derivative) in tangent.iter().enumerate() {
        close(
            *derivative,
            (values[1][i] - values[0][i]) / (2.0 * h),
            4e-6,
            3e-8,
        );
    }
}

#[test]
fn stress_pullback_admits_only_the_primal_and_refuses_invalid_or_cancelled_work() {
    let law = law(true, true);
    let op = cartesian(&law, ElasticityFormulation3::BBar);
    let (u, direction) = fields(op.nodes());
    let preparation = op.prepare_thermoelastic(POINTS, poll).unwrap();
    let before = op
        .thermoelastic_response(&u, TOTAL_STRAIN, POINTS, poll)
        .unwrap();
    let derivatives = cotangents(&before);
    let mut calls = 0;
    let expected = op
        .thermoelastic_stress_pullback(&u, TOTAL_STRAIN, &derivatives, || {
            calls += 1;
            poll()
        })
        .unwrap();
    let large: Vec<_> = derivatives.iter().map(|d| d.map(|v| 1e9 * v)).collect();
    let algebraic = op
        .thermoelastic_stress_pullback(&u, TOTAL_STRAIN, &large, poll)
        .unwrap();
    for (actual, expected) in algebraic.displacement.iter().zip(&expected.displacement) {
        close(*actual, 1e9 * expected, 2e-12, 2e-3);
    }
    assert!(
        op.thermoelastic_stress_pullback(&direction, TOTAL_STRAIN, &derivatives, poll)
            .is_err()
    );
    assert!(
        op.thermoelastic_stress_pullback(&u[..u.len() - 1], TOTAL_STRAIN, &derivatives, poll)
            .is_err()
    );
    assert!(
        op.thermoelastic_stress_pullback(
            &u,
            TOTAL_STRAIN,
            &derivatives[..derivatives.len() - 1],
            poll
        )
        .is_err()
    );
    assert!(
        op.thermoelastic_stress_pullback(&u, f64::NAN, &derivatives, poll)
            .is_err()
    );
    let mut bad = derivatives.clone();
    bad[0][3] = f64::NAN;
    assert!(
        op.thermoelastic_stress_pullback(&u, TOTAL_STRAIN, &bad, poll)
            .is_err()
    );
    for stop in [1, calls / 2, calls] {
        let mut at = 0;
        assert!(matches!(
            op.thermoelastic_stress_pullback(&u, TOTAL_STRAIN, &derivatives, || {
                at += 1;
                if at == stop {
                    ControlFlow::Break(())
                } else {
                    poll()
                }
            }),
            Err(ElasticityError3::Cancelled),
        ));
    }
    op.validate_eigenstrain(&preparation).unwrap();
    assert_eq!(
        op.thermoelastic_response(&u, TOTAL_STRAIN, POINTS, poll)
            .unwrap(),
        before
    );
}

#[path = "cutfem_thermal_stress_sensitivity/entropy.rs"]
mod entropy;
