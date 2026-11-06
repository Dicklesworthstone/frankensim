//! Native combined-load observations reuse the same solved stresses and frames.
use super::*;
use fs_material::strength::OrthotropicTsaiWu;
use fs_solid::strength::orthotropic::{
    OrthotropicTsaiWuMaterialBinding as InteractionBinding,
    OrthotropicTsaiWuMaterialField as InteractionField,
};

const STRENGTH: [&str; 9] = [
    "tensile_axis_1",
    "compression_1",
    "tensile_axis_2",
    "compression_2",
    "tensile_axis_3",
    "compression_3",
    "strength_12",
    "strength_23",
    "strength_31",
];
const INTERACTION: [&str; 3] = ["f12", "f23", "f31"];

fn interaction_binding(
    source: &MaterialCard,
    temperature: f64,
    rotation: [[f64; 3]; 3],
) -> InteractionBinding {
    let point = QueryPoint::new().with("T", temperature).unwrap();
    let selection = MaterialPropertySelection::SingleClaimOnly;
    let state = resolve_orthotropic_elastic_state_point(source, &point, selection, 0.01).unwrap();
    let criterion =
        OrthotropicTsaiWu::resolve(source, &point, STRENGTH, INTERACTION, selection).unwrap();
    InteractionBinding::try_from_resolved_orthotropic_state(
        &state,
        criterion,
        rotation,
        ORIENTATION,
    )
    .unwrap()
}

#[test]
fn g1_g3_native_triaxial_tension_crosses_interaction_surface_in_rotated_material_axes() {
    with_cx(false, |cx| {
        let source = card(0);
        for rotation in [IDENTITY, fs_material::tensor::rotation([1., 2., 3.], 0.7)] {
            let tensile = binding(&source, 300., rotation);
            let interaction = interaction_binding(&source, 300., rotation);
            assert_eq!(tensile.material(), interaction.material());
            let response = tet_report(
                &tensile,
                [-600. / E[0], -60. / E[1], -30. / E[2], 0., 0., 0.],
                cx,
            );
            let old = response
                .assess_orthotropic_tensile_strength(Field::Uniform(&tensile), 1, cx)
                .unwrap();
            assert!(!old.reaches_limit());
            let report = response
                .assess_orthotropic_strength(InteractionField::Uniform(&interaction), 1, cx)
                .unwrap();
            let assessment = report.governing_sample().unwrap().assessment();
            // All three normalized tensile demands are .6; Xc=1.5Xt.
            // A=.72 (the interaction coefficients sum to zero), B=.6,
            // hence (B+sqrt(B^2+4A))/2=1.2, despite every axis being below 1.
            assert!((assessment.utilization() - 1.2).abs() < 1e-12);
            assert!(report.reaches_limit());
            for (actual, expected) in assessment
                .material_stress_pa()
                .into_iter()
                .zip([600., 60., 30., 0., 0., 0.])
            {
                assert!((actual - expected).abs() < 1e-9);
            }
            assert_eq!(report.response_identity(), response.identity());
            assert_ne!(report.identity(), old.identity());
            assert_eq!(
                report,
                response
                    .assess_orthotropic_strength(
                        InteractionField::PerElement(&[interaction.clone()]),
                        1,
                        cx
                    )
                    .unwrap()
            );
        }
    });
}

#[test]
fn g1_q1_local_temperature_and_orientation_remain_at_each_native_point() {
    with_cx(false, |cx| {
        let source = card(0);
        let cold = binding(&source, 300., IDENTITY);
        let rotation = fs_material::tensor::rotation([1., 2., 3.], 0.7);
        let hot = binding(&source, 500., rotation);
        let tensile: Vec<_> = (0..8)
            .map(|i| if i < 4 { cold.clone() } else { hot.clone() })
            .collect();
        let interactions: Vec<_> = (0..8)
            .map(|i| {
                interaction_binding(
                    &source,
                    if i < 4 { 300. } else { 500. },
                    if i < 4 { IDENTITY } else { rotation },
                )
            })
            .collect();
        let materials: Vec<_> = tensile.iter().map(|b| b.material().clone()).collect();
        let states: Vec<_> = tensile
            .iter()
            .map(|b| thermal(b, [-600. / E[0], -60. / E[1], -30. / E[2], 0., 0., 0.]))
            .collect();
        let fixed: Vec<_> = (0..24).collect();
        let problem = Hex8ElasticProblem {
            nodes_m: &HEX,
            hexahedra: &[[0, 1, 2, 3, 4, 5, 6, 7]],
            materials: TetMaterialField::Uniform(cold.material()),
            fixed_dofs: &fixed,
            formulation: Formulation::Standard,
            budget: TetAssemblyBudget::standard(),
            stiffness_points_per_dir: 2,
        };
        let response = problem
            .solve_thermal_response(
                QuadraticThermalStrainField::ThermomechanicalAtQuadrature {
                    states: &states,
                    materials: &materials,
                    points_per_element: 8,
                },
                TetStaticSolveConfig::default(),
                FRAME,
                0.01,
                cx,
            )
            .unwrap();
        let report = response
            .assess_orthotropic_strength(InteractionField::PerPoint(&interactions), 8, cx)
            .unwrap();
        assert_eq!(report.samples().len(), 8);
        for (i, sample) in report.samples().iter().enumerate() {
            assert_eq!(sample.quadrature_point(), Some(i));
            assert_eq!(
                sample.limit_source_identity(),
                interactions[i].limit().resolved().identity()
            );
            assert!(
                (sample.assessment().utilization() - if i < 4 { 1.2 } else { 2.4 }).abs() < 1e-12
            );
        }
        assert_eq!(
            response.assess_orthotropic_strength(
                InteractionField::Uniform(&interactions[0]),
                8,
                cx
            ),
            Err(TensileAssessmentError::SourceMismatch {
                element: 0,
                point: Some(4)
            })
        );
        assert_eq!(
            response.assess_orthotropic_strength(InteractionField::PerPoint(&interactions), 7, cx),
            Err(TensileAssessmentError::Budget)
        );
    });
}

#[test]
fn g3_source_order_and_g4_stale_condition_and_cancellation_are_preserved() {
    let source = card(0);
    let tensile = binding(&source, 300., IDENTITY);
    let interaction = interaction_binding(&source, 300., IDENTITY);
    let hot = interaction_binding(&source, 500., IDENTITY);
    let point = QueryPoint::new().with("T", 300.).unwrap();
    let state = resolve_orthotropic_elastic_state_point(
        &source,
        &point,
        MaterialPropertySelection::SingleClaimOnly,
        0.01,
    )
    .unwrap();
    assert!(
        InteractionBinding::try_from_resolved_orthotropic_state(
            &state,
            hot.limit().clone(),
            IDENTITY,
            ORIENTATION
        )
        .is_err()
    );
    let response = with_cx(false, |cx| {
        tet_report(&tensile, [-0.0001, 0., 0., 0., 0., 0.], cx)
    });
    with_cx(false, |cx| {
        let first = response
            .assess_orthotropic_strength(InteractionField::Uniform(&interaction), 1, cx)
            .unwrap();
        let mut reordered = STRENGTH;
        reordered.swap(6, 7);
        let swapped = OrthotropicTsaiWu::try_from_resolved(
            interaction.limit().resolved().clone(),
            reordered,
            INTERACTION,
        )
        .unwrap();
        let swapped = InteractionBinding::try_from_resolved_orthotropic_state(
            &state,
            swapped,
            IDENTITY,
            ORIENTATION,
        )
        .unwrap();
        let second = response
            .assess_orthotropic_strength(InteractionField::Uniform(&swapped), 1, cx)
            .unwrap();
        // Same resolved source bundle and same zero-shear result; different
        // meaning assigned to the shear directions must still change identity.
        assert_eq!(
            first.samples()[0].limit_source_identity(),
            second.samples()[0].limit_source_identity()
        );
        assert_eq!(
            first.samples()[0].assessment(),
            second.samples()[0].assessment()
        );
        assert_ne!(first.identity(), second.identity());
        assert_eq!(
            response.assess_orthotropic_strength(InteractionField::Uniform(&hot), 1, cx),
            Err(TensileAssessmentError::SourceMismatch {
                element: 0,
                point: None
            })
        );
    });
    with_cx(true, |cx| {
        assert_eq!(
            response.assess_orthotropic_strength(InteractionField::Uniform(&interaction), 1, cx),
            Err(TensileAssessmentError::Cancelled)
        );
    });
}
