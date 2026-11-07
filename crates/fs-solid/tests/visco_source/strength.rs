//! Analytic strength observations of actual native hereditary stresses.

use super::*;
use fs_material::strength::{OrthotropicTensileLimit, OrthotropicTsaiWu};
use fs_solid::strength::TensileAssessmentError;

const STRENGTHS: [&str; 9] = ["xt", "xc", "yt", "yc", "zt", "zc", "s12", "s23", "s31"];
const INTERACTIONS: [&str; 3] = ["f12", "f23", "f31"];

fn strength_card(changing_stiffness: bool) -> MaterialCard {
    let original = card(changing_stiffness);
    let mut claims = original.claims().clone();
    for (name, value) in STRENGTHS
        .into_iter()
        .zip([15.0, 30.0, 40.0, 50.0, 60.0, 70.0, 8.0, 9.0, 10.0])
    {
        let mut entry = claim(name, Pressure::DIMS, [value, 0.8 * value, 0.8 * value]);
        // Static strengths depend on supported external temperature, never on
        // the dynamic frequency coordinate of the viscoelastic coefficients.
        entry.validity = ValidityDomain::unconstrained().with("T", 300.0, 400.0);
        claims.insert_claim(entry).unwrap();
    }
    for name in INTERACTIONS {
        let mut entry = claim(name, Dims::NONE, [0.0; 3]);
        entry.validity = ValidityDomain::unconstrained().with("T", 300.0, 400.0);
        claims.insert_claim(entry).unwrap();
    }
    MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap()
}

fn limit(card: &MaterialCard, temperature: f64) -> OrthotropicTensileLimit {
    OrthotropicTensileLimit::resolve(
        card,
        &point(temperature),
        ["xt", "yt", "zt"],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn interaction(card: &MaterialCard) -> OrthotropicTsaiWu {
    OrthotropicTsaiWu::resolve(
        card,
        &point(300.0),
        STRENGTHS,
        INTERACTIONS,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

#[test]
fn actual_maxwell_and_kelvin_stress_cross_strength_then_relax_at_fixed_displacement() {
    with_cx(|cx, _| {
        let card = strength_card(false);
        let source = source(&card);
        let system = tet(&source, &[3], cx);
        let criterion = limit(&card, 300.0);
        let prepared = system
            .prepare_orthotropic_strength(vec![criterion.clone()], cx)
            .unwrap();
        let mut state = system.state(cx).unwrap();
        let mut u = [0.0; 12];
        let (strain, h) = (0.01, 0.1);
        u[3] = strain;
        let response = state.step_internal_force(&u, h, cx).unwrap();
        let before = state.energy();
        let report = prepared.assess(&state, 1, cx).unwrap();
        let sample = report.governing_sample().unwrap();
        let maxwell = 600.0 * strain * phi(h, 0.07);
        let expected = 1200.0 * strain + maxwell + 8.0 * strain / h;
        close(sample.assessment.material_stress_pa()[0], expected);
        close(sample.assessment.utilization(), expected / 15.0);
        close(response.full_force_n[3] / VOLUME, expected);
        assert!(report.reaches_limit());
        assert!(
            1200.0 * strain < 15.0,
            "equilibrium-only checking would miss this crossing"
        );
        assert_eq!(
            sample.material_state_identity,
            source.material.state().identity()
        );
        assert_eq!(
            sample.limit_source_identity,
            criterion.resolved().identity()
        );
        assert_eq!(prepared.assess(&state, 1, cx).unwrap(), report);
        assert_eq!(state.energy(), before);
        assert_eq!(state.displacement_m(), &u);

        // Same displacement, different physical history and actual stress.
        let held = 0.2_f64;
        state.step_internal_force(&u, held, cx).unwrap();
        let relaxed = prepared.assess(&state, 1, cx).unwrap();
        close(
            relaxed.samples()[0].assessment.material_stress_pa()[0],
            1200.0 * strain + maxwell * (-held / 0.07).exp(),
        );
        assert!(!relaxed.reaches_limit());
    });
}

#[test]
fn rotated_material_history_preserves_signed_normals_and_physical_shear_for_tsai_wu() {
    with_cx(|cx, _| {
        let card = strength_card(false);
        let mut source = source(&card);
        let c = 0.5_f64.sqrt();
        source.material_to_mesh = [[c, -c, 0.0], [c, c, 0.0], [0.0, 0.0, 1.0]];
        let system = tet(&source, &(0..12).collect::<Vec<_>>(), cx);
        let criterion = interaction(&card);
        let prepared = system
            .prepare_orthotropic_strength(vec![criterion.clone()], cx)
            .unwrap();
        let tensile = system
            .prepare_orthotropic_strength(vec![limit(&card, 300.0)], cx)
            .unwrap();
        let strain = [0.004, -0.002, 0.001, 0.012, 0.002, -0.001];
        let e = [
            [strain[0], strain[3], strain[5]],
            [strain[3], strain[1], strain[4]],
            [strain[5], strain[4], strain[2]],
        ];
        // Independently manufacture affine world displacement R E R^T x.
        let r = source.material_to_mesh;
        let world: [[f64; 3]; 3] = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                (0..3)
                    .flat_map(|a| (0..3).map(move |b| r[i][a] * e[a][b] * r[j][b]))
                    .sum()
            })
        });
        let u: Vec<_> = NODES
            .iter()
            .flat_map(|x| world.map(|row| row.iter().zip(x).map(|(a, b)| a * b).sum::<f64>()))
            .collect();
        let h = 0.05;
        let mut state = system.state(cx).unwrap();
        state.step_internal_force(&u, h, cx).unwrap();
        let report = prepared.assess(&state, 1, cx).unwrap();
        let actual = &report.samples()[0].assessment;
        let mut expected = strain.map(|v| 1200.0 * v);
        expected[0] += strain[0] * (600.0 * phi(h, 0.07) + 8.0 / h);
        for (a, b) in actual.material_stress_pa().into_iter().zip(expected) {
            close(a, b);
        }
        close(
            actual.utilization(),
            criterion
                .assess_material_stress(&expected)
                .unwrap()
                .utilization(),
        );
        assert!(
            report.reaches_limit(),
            "actual shear exceeds its source strength"
        );
        assert!(!tensile.assess(&state, 1, cx).unwrap().reaches_limit());
    });
}

#[test]
fn original_thermal_quadrature_and_bbar_stresses_are_observed_without_averaging() {
    with_cx(|cx, _| {
        let card = strength_card(false);
        let source = source(&card);
        for formulation in [Formulation::Standard, Formulation::BBar] {
            let system = prepare(
                &source,
                SolidMeshRef::Hex8 {
                    nodes_m: &HEX,
                    elements: &HEXES,
                    formulation,
                    stiffness_points_per_dir: 2,
                },
                &[],
                cx,
            );
            let prepared = system
                .prepare_orthotropic_strength(vec![limit(&card, 300.0)], cx)
                .unwrap();
            let id = system
                .reference_material(0)
                .unwrap()
                .material_state_identity;
            let free: Vec<_> = (0..8).map(|q| -0.0006 * (q + 1) as f64).collect();
            let make_thermal = |q: usize, temperature: f64| {
                TetThermalStrainState::try_new(
                    temperature,
                    350.0,
                    [free[q], 0.0, 0.0, 0.0, 0.0, 0.0],
                    id,
                    ContentHash([73; 32]),
                    0.03,
                )
                .unwrap()
            };
            let thermal: Vec<_> = (0..8).map(|q| make_thermal(q, 300.0)).collect();
            let field = |states| QuadraticThermalStrainField::AtQuadrature {
                states,
                points_per_element: 8,
            };
            let mut state = system.state(cx).unwrap();
            let h = 0.02;
            state
                .step_internal_force_with_thermal_strain(&[0.0; 24], h, field(&thermal), cx)
                .unwrap();
            let before = state.energy();
            let report = prepared.assess(&state, 8, cx).unwrap();
            assert_eq!(report.samples().len(), 8);
            let average = free.iter().sum::<f64>() / 8.0;
            for (q, sample) in report.samples().iter().enumerate() {
                let strain = if formulation == Formulation::BBar {
                    [
                        -(2.0 * free[q] + average) / 3.0,
                        (free[q] - average) / 3.0,
                        (free[q] - average) / 3.0,
                    ]
                } else {
                    [-free[q], 0.0, 0.0]
                };
                let mut expected = strain.map(|v| 1200.0 * v);
                expected[0] += strain[0] * (600.0 * phi(h, 0.07) + 8.0 / h);
                for (a, b) in sample
                    .assessment
                    .material_normal_stresses_pa()
                    .into_iter()
                    .zip(expected)
                {
                    close(a, b);
                }
                assert_eq!((sample.element, sample.quadrature_point), (0, q));
            }
            assert_eq!(report.governing_sample().unwrap().quadrature_point, 7);
            assert_eq!(state.energy(), before);
            assert!(matches!(
                prepared.assess(&state, 7, cx),
                Err(TensileAssessmentError::Budget)
            ));

            let mut changed = thermal.clone();
            changed[7] = make_thermal(7, 301.0);
            state
                .step_internal_force_with_thermal_strain(&[0.0; 24], h, field(&changed), cx)
                .unwrap();
            assert!(
                matches!(
                    prepared.assess(&state, 8, cx),
                    Err(TensileAssessmentError::Invalid(_))
                ),
                "a frozen coefficient query cannot confer strength at another actual temperature"
            );
        }
    });
}

#[test]
fn refresh_rebinds_strength_while_foreign_owner_source_budget_and_cancellation_refuse() {
    with_cx(|cx, gate| {
        let card = strength_card(false);
        let source = source(&card);
        let system = tet(&source, &[3], cx);
        let prepared = system
            .prepare_orthotropic_strength(vec![limit(&card, 300.0)], cx)
            .unwrap();
        let mut state = system.state(cx).unwrap();
        let mut u = [0.0; 12];
        u[3] = 0.005;
        state.step_internal_force(&u, 0.02, cx).unwrap();
        let before = state.energy();
        assert!(
            system
                .prepare_orthotropic_strength::<OrthotropicTensileLimit>(vec![], cx)
                .is_err()
        );
        assert!(
            system
                .prepare_orthotropic_strength(vec![limit(&card, 350.0)], cx)
                .is_err()
        );
        assert!(
            system
                .prepare_orthotropic_strength(vec![limit(&strength_card(true), 300.0)], cx)
                .is_err()
        );
        let independent = tet(&source, &[3], cx);
        assert!(
            prepared
                .assess(&independent.state(cx).unwrap(), 1, cx)
                .is_err()
        );
        assert!(matches!(
            prepared.assess(&state, 0, cx),
            Err(TensileAssessmentError::Budget)
        ));

        let material = TetViscoMaterial::from_source(&source.material, IDENTITY).unwrap();
        let unbound = ViscoElements::new(
            SolidMeshRef::Tet4 {
                nodes_m: &NODES,
                elements: &TETS,
            },
            TetViscoField::Uniform(&material),
            &[],
            TetAssemblyBudget::standard(),
            ViscoStateBudget {
                max_internal_variables: 10_000,
            },
            cx,
        )
        .unwrap();
        assert!(
            unbound
                .prepare_orthotropic_strength(vec![limit(&card, 300.0)], cx)
                .is_err()
        );

        let next = system
            .at_query_points(&[&card], &[point(350.0)], cx)
            .unwrap();
        let refreshed_limit = limit(&card, 350.0);
        let next_prepared = next
            .prepare_orthotropic_strength(vec![refreshed_limit.clone()], cx)
            .unwrap();
        let mut changed = state.rebind_materials(&next, cx).unwrap();
        assert!(prepared.assess(&changed.history, 1, cx).is_err());
        assert!(next_prepared.assess(&state, 1, cx).is_err());
        assert!(
            next_prepared.assess(&changed.history, 1, cx).is_err(),
            "zero-duration source refresh cannot invent a Kelvin endpoint rate"
        );
        assert!(changed.history.step_internal_force(&u, 0.0, cx).is_err());
        assert!(next_prepared.assess(&changed.history, 1, cx).is_err());
        let mut continued_u = u;
        continued_u[3] = 0.006;
        let h = 0.03_f64;
        let force = changed
            .history
            .step_internal_force(&continued_u, h, cx)
            .unwrap();
        let report = next_prepared.assess(&changed.history, 1, cx).unwrap();
        let sample = &report.samples()[0];
        let spring =
            u[3] * phi(0.02, 0.07) * (-h / 0.14).exp() + (continued_u[3] - u[3]) * phi(h, 0.14);
        close(
            sample.assessment.material_stress_pa()[0],
            1200.0 * continued_u[3] + 600.0 * spring + 4.0 * (continued_u[3] - u[3]) / h,
        );
        close(
            sample.assessment.material_stress_pa()[0],
            force.full_force_n[3] / VOLUME,
        );
        close(
            sample.assessment.utilization(),
            sample.assessment.material_stress_pa()[0] / 12.0,
        );
        assert_eq!(
            sample.limit_source_identity,
            refreshed_limit.resolved().identity()
        );
        assert_eq!(state.energy(), before);
        assert_eq!(
            next.equilibrium().total_mass_kg,
            system.equilibrium().total_mass_kg
        );

        gate.request();
        assert!(matches!(
            prepared.assess(&state, 1, cx),
            Err(TensileAssessmentError::Cancelled)
        ));
        assert!(matches!(
            system.prepare_orthotropic_strength(vec![limit(&card, 300.0)], cx),
            Err(TensileAssessmentError::Cancelled)
        ));
        assert_eq!(state.energy(), before);
    });
}
