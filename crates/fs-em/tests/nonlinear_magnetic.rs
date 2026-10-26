//! G1/G3/G4: source-only saturation, independent P1 solution/energy, refusals.
//! Manufactured laws verify implementation; they are not measured soft-iron data.

use fs_em::{
    planar::{
        MagneticError, PlanarMagnetostatics,
        nonlinear::{
            IsotropicMagneticMaterial, NonlinearMagneticControls, NonlinearMagneticError,
            PlanarNonlinearMagnetostatics,
        },
    },
    transfer::{CellSample, PlanarFrame},
};
use fs_evidence::ValidityDomain;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        MAGNETIC_PERMEABILITY_DIMS, MAGNETIC_PERMEABILITY_PROPERTY,
        nonlinear::{REVERSIBLE_FIELD_STRENGTH_PROPERTY, ReversibleMagneticMaterial},
        permanent::{FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS},
    },
    state_point::{
        MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
        resolve_material_state_point,
    },
};
use fs_qty::{QuantitySpec, Temperature};
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::CheckedCgConfig;

fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            temperature,
        )
        .unwrap()
        .with_quantity("B", QuantitySpec::dimensional(FLUX_DENSITY_DIMS), 0.0)
        .unwrap()
}
fn card(value: PropertyValue, interpolation: InterpolationPolicy) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(REVERSIBLE_FIELD_STRENGTH_PROPERTY, FIELD_STRENGTH_DIMS),
            value,
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    200.0,
                    500.0,
                )
                .with_quantity("B", QuantitySpec::dimensional(FLUX_DENSITY_DIMS), 0.0, 4.0),
            interpolation,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured reversible H(B)".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic ferromagnetic law".into(),
            phase: "solid".into(),
            process: "reversible fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn curve(knots: &[(f64, f64)]) -> MaterialCard {
    card(
        PropertyValue::Curve {
            abscissa: "B".into(),
            abscissa_dims: FLUX_DENSITY_DIMS,
            knots: knots.to_vec(),
            dims: FIELD_STRENGTH_DIMS,
        },
        InterpolationPolicy::LinearInside,
    )
}
fn law(source: &MaterialCard) -> ReversibleMagneticMaterial {
    ReversibleMagneticMaterial::resolve(
        source,
        &point(300.0),
        "B",
        4.0,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn limits() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits {
        max_cells: 10000,
        max_dofs: 30000,
        max_triplets: 100000,
    }
}
fn controls() -> NonlinearMagneticControls {
    NonlinearMagneticControls {
        absolute_residual_tolerance_a_m: 1e-11,
        relative_residual_tolerance: 1e-10,
        max_newton_iterations: 30,
        max_backtracks: 30,
        linear: CheckedCgConfig {
            tolerance: 1e-12,
            max_iterations: 1000,
            max_corrections: 2,
        },
    }
}
fn square(depth: f64) -> TriComplex2 {
    TriComplex2::from_indexed_triangles(
        "nonlinear-square",
        vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5]],
        vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
        Metric2::planar(depth).unwrap(),
    )
    .unwrap()
}
fn boundary(offset: f64) -> Vec<(usize, f64)> {
    (0..4).map(|v| (v, offset)).collect()
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 2e-9 * (1.0 + b.abs()), "{a} != {b}");
}

#[test]
fn exact_source_energy_and_tangent_include_every_saturation_segment() {
    let source = curve(&[(0.0, 0.0), (1.0, 1.0), (4.0, 13.0)]);
    let law = law(&source);
    let response = law.evaluate(1.25).unwrap();
    near(response.field_strength_a_m, 2.0);
    near(response.secant_reluctivity_m_h, 1.6);
    near(response.differential_reluctivity_m_h, 4.0);
    near(response.energy_density_j_m3, 0.875);
    near(response.coenergy_density_j_m3, 1.625);
    near(law.evaluate(0.0).unwrap().secant_reluctivity_m_h, 1.0);
    near(law.evaluate(1.0).unwrap().differential_reluctivity_m_h, 4.0);
    near(law.evaluate(4.0).unwrap().energy_density_j_m3, 21.5);
    let epsilon = 1e-5;
    let derivative = (law.evaluate(1.25 + epsilon).unwrap().energy_density_j_m3
        - law.evaluate(1.25 - epsilon).unwrap().energy_density_j_m3)
        / (2.0 * epsilon);
    near(derivative, response.field_strength_a_m);
    assert_eq!(law.material().card_identity(), source.content_hash());
    assert!(law.source_path().knots.iter().all(|(_, a)| {
        a.receipt.selected
            == law
                .material()
                .property(REVERSIBLE_FIELD_STRENGTH_PROPERTY)
                .unwrap()
                .answer()
                .receipt
                .selected
    }));
}

#[test]
fn changing_only_the_source_curve_changes_the_solved_field_and_energy() {
    let straight = law(&curve(&[(0.0, 0.0), (4.0, 4.0)]));
    let saturated = law(&curve(&[(0.0, 0.0), (1.0, 1.0), (4.0, 13.0)]));
    let mut outputs = Vec::new();
    for material in [straight, saturated] {
        let problem = PlanarNonlinearMagnetostatics::new(
            square(0.3),
            vec![IsotropicMagneticMaterial::Reversible(material); 4],
            limits(),
            &mut || false,
        )
        .unwrap();
        let field = problem
            .solve(
                &[12.0; 4],
                &boundary(0.0),
                &[0.0; 5],
                controls(),
                &mut || false,
            )
            .unwrap();
        let amplitude = field.potential_wb_m()[4];
        // Independent one-dof weak equation: 2*H(2*a)=J/3.
        for h in field.magnetic_field_a_m() {
            near(h[0].hypot(h[1]), 2.0);
        }
        assert!(field.report().residual_norm_a_m <= field.report().residual_target_a_m);
        near(field.cell_energy_j().iter().sum(), field.energy_j());
        near(field.cell_coenergy_j().iter().sum(), field.coenergy_j());
        let frame = PlanarFrame::new([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap();
        let samples = [CellSample {
            triangle: 0,
            position_m: [0.5, 0.1, 0.1],
        }];
        let b = field
            .sample_cells(&frame, &samples, 0.0, 1, &mut || false)
            .unwrap();
        near(b[0][0], 2.0 * amplitude);
        near(b[0][1], 0.0);
        outputs.push((amplitude, field.energy_j(), field.coenergy_j()));
    }
    near(outputs[0].0, 1.0);
    near(outputs[0].1, 0.6);
    near(outputs[1].0, 0.625);
    near(outputs[1].1, 0.2625);
    near(outputs[1].2, 0.4875);
    assert_ne!(outputs[0], outputs[1]);
}

#[test]
fn high_contrast_saturation_backtracks_without_extrapolating() {
    let material = law(&curve(&[(0.0, 0.0), (0.25, 0.25), (4.0, 375.25)]));
    let problem = PlanarNonlinearMagnetostatics::new(
        square(1.0),
        vec![IsotropicMagneticMaterial::Reversible(material); 4],
        limits(),
        &mut || false,
    )
    .unwrap();
    let field = problem
        .solve(
            &[60.0; 4],
            &boundary(0.0),
            &[0.0; 5],
            controls(),
            &mut || false,
        )
        .unwrap();
    near(
        field.potential_wb_m()[4],
        0.5 * (0.25 + (10.0 - 0.25) / 100.0),
    );
    assert!(field.report().backtracks > 0);
}

#[test]
fn data_surface_changes_temperature_response_without_a_solver_change() {
    // H=B at 200 K, H=2B at 500 K. Two declared source triangles form the
    // entire rectangle; the exact T=300 slice crosses their shared diagonal.
    let source = card(
        PropertyValue::Surface2D {
            axes: ["B".into(), "T".into()],
            axis_dims: [FLUX_DENSITY_DIMS, Temperature::DIMS],
            samples: vec![
                ([0.0, 200.0], 0.0),
                ([4.0, 200.0], 4.0),
                ([0.0, 500.0], 0.0),
                ([4.0, 500.0], 8.0),
            ],
            triangles: vec![[0, 1, 3], [0, 3, 2]],
            dims: FIELD_STRENGTH_DIMS,
        },
        InterpolationPolicy::LinearInside,
    );
    let cold = ReversibleMagneticMaterial::resolve(
        &source,
        &point(200.0),
        "B",
        4.0,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let hot = ReversibleMagneticMaterial::resolve(
        &source,
        &point(500.0),
        "B",
        4.0,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    near(cold.evaluate(1.0).unwrap().field_strength_a_m, 1.0);
    near(hot.evaluate(1.0).unwrap().field_strength_a_m, 2.0);
    assert_eq!(
        cold.material().card_identity(),
        hot.material().card_identity()
    );
}

#[test]
fn invalid_laws_and_field_envelopes_refuse() {
    for knots in [
        vec![(0.0, 1.0), (4.0, 5.0)],
        vec![(0.0, 0.0), (1.0, 1.0), (4.0, 0.5)],
        vec![(0.0, 0.0), (1.0, 1.0), (4.0, 1.0)],
    ] {
        assert!(
            ReversibleMagneticMaterial::resolve(
                &curve(&knots),
                &point(300.0),
                "B",
                4.0,
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
        );
    }
    let discrete = card(
        PropertyValue::Curve {
            abscissa: "B".into(),
            abscissa_dims: FLUX_DENSITY_DIMS,
            knots: vec![(0.0, 0.0), (4.0, 4.0)],
            dims: FIELD_STRENGTH_DIMS,
        },
        InterpolationPolicy::TabulatedOnly,
    );
    assert!(
        ReversibleMagneticMaterial::resolve(
            &discrete,
            &point(300.0),
            "B",
            4.0,
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
    let law = law(&curve(&[(0.0, 0.0), (4.0, 4.0)]));
    for b in [-1.0, 4.01, f64::NAN, f64::INFINITY] {
        assert!(law.evaluate(b).is_err());
    }
    let problem = PlanarNonlinearMagnetostatics::new(
        square(1.0),
        vec![IsotropicMagneticMaterial::Reversible(law); 4],
        limits(),
        &mut || false,
    )
    .unwrap();
    assert!(matches!(
        problem.solve(
            &[0.0; 4],
            &boundary(0.0),
            &[0.0, 0.0, 0.0, 0.0, 3.0],
            controls(),
            &mut || false
        ),
        Err(NonlinearMagneticError::FieldRange { .. })
    ));
    let mut controls = controls();
    controls.max_newton_iterations = 0;
    assert!(matches!(
        problem.solve(&[12.0; 4], &boundary(0.0), &[0.0; 5], controls, &mut || {
            false
        }),
        Err(NonlinearMagneticError::Nonconvergence(_))
    ));
    assert!(matches!(
        problem.solve(&[12.0; 4], &boundary(0.0), &[0.0; 5], controls, &mut || {
            true
        }),
        Err(NonlinearMagneticError::Magnetic(MagneticError::Cancelled))
    ));
    assert!(
        problem
            .solve(&[0.0; 4], &[(0, 0.0)], &[0.0; 5], controls, &mut || false)
            .is_err()
    );
}

fn linear_material(mu: f64) -> fs_material::state_point::ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(MAGNETIC_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS),
            value: PropertyValue::Scalar {
                value: mu,
                dims: MAGNETIC_PERMEABILITY_DIMS,
            },
            validity: ValidityDomain::unconstrained(),
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured linear field".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "linear".into(),
            phase: "solid".into(),
            process: "fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_material_state_point(
        &card,
        &QueryPoint::new(),
        &[ScalarPropertyRequirement::try_new(
            MAGNETIC_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

#[test]
fn shared_linear_limit_gauge_and_current_reversal_agree() {
    let materials = vec![linear_material(0.75); 4];
    let linear =
        PlanarMagnetostatics::new(square(0.7), materials.clone(), limits(), &mut || false).unwrap();
    let nonlinear = PlanarNonlinearMagnetostatics::new(
        square(0.7),
        materials
            .into_iter()
            .map(IsotropicMagneticMaterial::Linear)
            .collect(),
        limits(),
        &mut || false,
    )
    .unwrap();
    for (offset, current) in [(0.0, 12.0), (100.0, -12.0)] {
        let reference = linear
            .solve(
                &[current; 4],
                &boundary(offset),
                controls().linear,
                &mut || false,
            )
            .unwrap();
        let result = nonlinear
            .solve(
                &[current; 4],
                &boundary(offset),
                &[offset; 5],
                controls(),
                &mut || false,
            )
            .unwrap();
        for (a, b) in result
            .potential_wb_m()
            .iter()
            .zip(reference.potential_wb_m())
        {
            near(*a, *b);
        }
        near(result.energy_j(), reference.energy_j());
        near(result.coenergy_j(), result.energy_j());
        assert_eq!(result.report().newton_iterations, 1);
    }
}

#[test]
fn oblique_field_with_mixed_materials_solves_the_shared_interface_problem() {
    let n = 8;
    let saturated = law(&curve(&[(0.0, 0.0), (1.0, 1.0), (4.0, 13.0)]));
    let left_slope = 1.0;
    let transverse_slope = 1.0;
    let left_b = 2.0_f64.sqrt();
    let left_h = 4.0 * left_b - 3.0;
    let right_slope = left_h / left_b;
    // Continuity of tangential H imposes nu_left*A_x_left=A_x_right.
    // Tangential A derivative is continuous, hence normal B is continuous.
    let exact = |x: f64, y: f64| {
        left_slope * x.min(0.5) + right_slope * (x - 0.5).max(0.0) + transverse_slope * y
    };
    let mut vertices = Vec::new();
    let mut initial = Vec::new();
    let mut boundary = Vec::new();
    for y in 0..=n {
        for x in 0..=n {
            let xx = x as f64 / n as f64;
            let yy = y as f64 / n as f64;
            let value = exact(xx, yy);
            let outer = x == 0 || y == 0 || x == n || y == n;
            if outer {
                boundary.push((vertices.len(), value));
            }
            initial.push(
                value
                    + if outer {
                        0.0
                    } else {
                        0.02 * (std::f64::consts::PI * xx).sin() * (std::f64::consts::PI * yy).sin()
                    },
            );
            vertices.push([xx, yy]);
        }
    }
    let mut faces = Vec::new();
    let mut materials = Vec::new();
    for y in 0..n {
        for x in 0..n {
            let a = (y * (n + 1) + x) as u32;
            let d = a + (n + 1) as u32;
            faces.extend([[a, a + 1, d + 1], [a, d + 1, d]]);
            let material = if x < n / 2 {
                IsotropicMagneticMaterial::Reversible(saturated.clone())
            } else {
                IsotropicMagneticMaterial::Linear(linear_material(1.0))
            };
            materials.extend([material.clone(), material]);
        }
    }
    let mesh = TriComplex2::from_indexed_triangles(
        "nonlinear-interface",
        vertices,
        faces,
        Metric2::planar(0.4).unwrap(),
    )
    .unwrap();
    let problem =
        PlanarNonlinearMagnetostatics::new(mesh, materials, limits(), &mut || false).unwrap();
    let field = problem
        .solve(
            &vec![0.0; problem.mesh().faces().len()],
            &boundary,
            &initial,
            controls(),
            &mut || false,
        )
        .unwrap();
    for (point, &potential) in problem.mesh().vertices().iter().zip(field.potential_wb_m()) {
        near(potential, exact(point[0], point[1]));
    }
    for (b, h) in field
        .magnetic_flux_density_t()
        .iter()
        .zip(field.magnetic_field_a_m())
    {
        near(b[0], transverse_slope);
        near(h[1], -right_slope);
    }
    let left_energy = 0.5 + (left_b - 1.0) + 2.0 * (left_b - 1.0).powi(2);
    let right_energy = 0.5 * (right_slope * right_slope + 1.0);
    near(field.energy_j(), 0.2 * (left_energy + right_energy));
}

#[test]
fn late_cancellation_never_publishes_a_partial_field() {
    let problem = PlanarNonlinearMagnetostatics::new(
        square(1.0),
        vec![
            IsotropicMagneticMaterial::Reversible(law(&curve(&[
                (0.0, 0.0),
                (1.0, 1.0),
                (4.0, 13.0)
            ])));
            4
        ],
        limits(),
        &mut || false,
    )
    .unwrap();
    let mut calls = 0;
    let reference = problem
        .solve(
            &[12.0; 4],
            &boundary(0.0),
            &[0.0; 5],
            controls(),
            &mut || {
                calls += 1;
                false
            },
        )
        .unwrap();
    let mut repeated_calls = 0;
    assert!(matches!(
        problem.solve(
            &[12.0; 4],
            &boundary(0.0),
            &[0.0; 5],
            controls(),
            &mut || {
                repeated_calls += 1;
                repeated_calls == calls
            }
        ),
        Err(NonlinearMagneticError::Magnetic(MagneticError::Cancelled))
    ));
    let again = problem
        .solve(
            &[12.0; 4],
            &boundary(0.0),
            &[0.0; 5],
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(reference.potential_wb_m(), again.potential_wb_m());
}
