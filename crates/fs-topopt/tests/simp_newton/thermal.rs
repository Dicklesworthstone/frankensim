//! Prescribed-temperature design through the existing constrained Newton owner.
use super::*;
use fs_ascent::reverse_newton::{
    AugmentedNewtonBudget, AugmentedNewtonOptions, AugmentedNewtonStop, AugmentedNewtonStudy,
};
use fs_blake3::ContentHash;
use fs_opt::reverse::physics::box_map::LogisticBox;
use fs_rep_mesh::TetComplex;
use fs_solid::linear3::{
    TetElasticMaterial, TetMaterialField, TetThermalStrainField, TetThermalStrainState,
};
use fs_topopt::RobustPipeline;
use fs_topopt::elasticity::second_order::thermal::ThermalLoadCase;
use fs_topopt::simp_newton::{MASS_DIMS, SimpNewtonResources};

fn thermal_state(material: &TetElasticMaterial, axial: f64, bound: f64) -> TetThermalStrainState {
    TetThermalStrainState::try_new(
        320.0,
        300.0,
        [0.0, 0.0, axial, 0.0, 0.0, 0.0],
        material.material_state_identity,
        ContentHash([49; 32]),
        bound,
    )
    .unwrap()
}
fn thermal_fixture(
    bound: f64,
) -> (
    DesignPipeline,
    DensityElasticity,
    TetElasticMaterial,
    TetThermalStrainState,
    Vec<f64>,
) {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let positions = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let material = TetElasticMaterial::try_new(6.0, 120.0, 0.0, [47; 32]).unwrap();
    let thermal = thermal_state(&material, 0.01, bound);
    let elasticity = DensityElasticity::with_material_field(
        &mesh,
        &positions,
        TetMaterialField::Uniform(&material),
        &|p| p[2] == 0.0,
    )
    .unwrap();
    let pipeline = DesignPipeline {
        filter: DensityFilter::new(&mesh, &positions, 0.15),
        params: SimpParams {
            e_min: 0.05,
            penal: 3.0,
            beta: 0.0,
            eta: 0.5,
        },
    };
    let mut force = vec![0.0; 12];
    force[11] = 0.06;
    (pipeline, elasticity, material, thermal, force)
}

#[test]
fn thermal_newton_finds_the_interior_optimum_and_an_active_physical_mass_cap() {
    with_cx(|_, cx| {
        for cap_kg in [0.55, 0.8] {
            let (pipeline, elasticity, _, thermal, force) = thermal_fixture(0.1);
            let loads = [ThermalLoadCase {
                force: &force,
                weight: 1.0,
                strains: TetThermalStrainField::Uniform(&thermal),
            }];
            let model = SimpNewtonModel::new_thermal(
                &pipeline,
                &elasticity,
                &loads,
                &[1.0 / 6.0],
                options(),
                budget(),
            )
            .unwrap();
            let raw_c = model.compliance("thermal compliance").unwrap();
            let raw_mass = model.mass("physical mass").unwrap();
            let chart = LogisticBox::new(&[0.05], &[0.95], Dims::NONE).unwrap();
            let c = chart.bind(&raw_c).unwrap();
            let mass = chart.bind(&raw_mass).unwrap();
            let start = chart
                .encode(&[if cap_kg < 0.6 { 0.45 } else { 0.78 }], Some(cx))
                .unwrap();
            let scale = c.value_gradient(&start, Some(cx)).unwrap().value;
            let mut builder = ProblemBuilder::new();
            let x = builder
                .var("bounded raw occupancy", Manifold::Rn { dim: 1 }, Dims::NONE)
                .unwrap();
            let cn = builder
                .pde_residual("thermal compliance", x, true, COMPLIANCE_DIMS)
                .unwrap();
            let mn = builder
                .pde_residual("physical mass", x, true, MASS_DIMS)
                .unwrap();
            let unit_c = builder.konst(scale, COMPLIANCE_DIMS).unwrap();
            let objective = builder.div(cn, unit_c).unwrap();
            builder.objective(objective, Sense::Minimize, 1.0).unwrap();
            let unit_mass = builder.konst(1.0, MASS_DIMS).unwrap();
            let normalized_mass = builder.div(mn, unit_mass).unwrap();
            let limit = builder.konst(cap_kg, Dims::NONE).unwrap();
            let constraint = builder.sub(normalized_mass, limit).unwrap();
            builder
                .constraint(
                    constraint,
                    ConstraintKind::LeZero,
                    "source-density mass allowance",
                )
                .unwrap();
            let problem = builder.finish();
            let oracle = ReverseProblem::new_with_physics(
                &problem,
                ReverseLimits {
                    max_nodes: 128,
                    max_scalar_slots: 1024,
                },
                &[
                    PhysicsBinding {
                        node: cn,
                        model: &c,
                    },
                    PhysicsBinding {
                        node: mn,
                        model: &mass,
                    },
                ],
                None,
            )
            .unwrap();
            let mut study = AugmentedNewtonStudy::new(
                &oracle,
                &start,
                AugmentedNewtonOptions::default(),
                Some(cx),
            )
            .unwrap();
            let result = study
                .run(
                    1e-6,
                    300,
                    AugmentedNewtonBudget {
                        objective_evals: 600,
                        hessian_products: 3000,
                        multiplier_updates: 100,
                    },
                    Some(cx),
                )
                .unwrap();
            assert_eq!(result.stop, AugmentedNewtonStop::Converged, "{result:?}");
            assert!(result.kkt.within_tolerance(1e-6));
            let rho = chart.decode(&result.x, Some(cx)).unwrap()[0];
            // One axial tetrahedral DOF: k=E*V=20 and q=k*alpha=.2.
            // c(m)=F^2/(20m)+2F alpha+20m alpha^2 has its interior minimum
            // at m=F/(20 alpha)=.3. Reference mass is rho_material*V=1 kg.
            let unconstrained = ((0.3_f64 - 0.05) / 0.95).cbrt();
            let expected_rho = unconstrained.min(cap_kg);
            close(&[rho], &[expected_rho], 2e-4);
            let m = 0.05 + 0.95 * rho.powi(3);
            let compliance =
                0.06_f64.powi(2) / (20.0 * m) + 2.0 * 0.06 * 0.01 + 20.0 * m * 0.01_f64.powi(2);
            close(&[result.objective], &[compliance / scale], 2e-9);
            close(&result.constraints, &[rho - cap_kg], 2e-10);
            if cap_kg < unconstrained {
                assert!(result.multipliers[0] > 0.0);
            } else {
                assert!(result.constraints[0] < -0.1);
            }
            assert!(elasticity.moduli.iter().all(|m| *m == 1.0));
        }
    });
}

#[test]
fn thermal_views_keep_original_sources_and_share_cache_budget_and_failed_trial_rollback() {
    with_cx(|gate, cx| {
        let (pipeline, elasticity, material, mut thermal, mut force) = thermal_fixture(0.03);
        let retained_thermal = thermal;
        let retained_force = force.clone();
        let cooling = thermal_state(&material, -0.005, 0.03);
        let mut second = vec![0.0; 12];
        second[11] = 0.03;
        let resources = SimpNewtonResources::new(budget());
        let model = {
            let loads = [
                ThermalLoadCase {
                    force: &force,
                    weight: 1.0,
                    strains: TetThermalStrainField::Uniform(&thermal),
                },
                ThermalLoadCase {
                    force: &second,
                    weight: 0.0,
                    strains: TetThermalStrainField::Uniform(&cooling),
                },
            ];
            SimpNewtonModel::new_thermal_with_resources(
                &pipeline,
                &elasticity,
                &loads,
                &[1.0 / 6.0],
                options(),
                resources.clone(),
            )
            .unwrap()
        };
        // The model owns these source snapshots, exactly as it owns dead forces.
        force[11] = 1.0;
        thermal = thermal_state(&material, 0.02, 0.03);
        assert_eq!(force[11], 1.0);
        assert_eq!(thermal.free_strain_mandel()[2], 0.02);
        let c = model.compliance("thermal").unwrap();
        let volume = model.volume("volume").unwrap();
        let mass = model.mass("mass").unwrap();
        let selected = model
            .load_case_compliance(1, "zero-weight cooling case")
            .unwrap();
        let x = [0.5];
        let accepted = c.value_gradient(&x, Some(cx)).unwrap();
        let prepared = model.work();
        assert_eq!(prepared.solve_reservations, 6);
        assert_eq!(prepared.started_solves, 6);
        let m = 0.05 + 0.95 * 0.5_f64.powi(3);
        let expected = (0.06 + m * 0.2).powi(2) / (20.0 * m);
        close(&[accepted.value], &[expected], 2e-11);
        close(
            &[mass.value_gradient(&x, Some(cx)).unwrap().value],
            &[0.5],
            2e-12,
        );
        close(
            &[volume.value_gradient(&x, Some(cx)).unwrap().value],
            &[0.5],
            2e-12,
        );
        assert_eq!(model.work(), prepared);
        let hv = c.hessian_vector_product(&x, &[0.2], Some(cx)).unwrap();
        assert_eq!(
            model.work().solve_reservations,
            prepared.solve_reservations + 3
        );
        let case_hv = selected
            .hessian_vector_product(&x, &[0.2], Some(cx))
            .unwrap();
        close(&case_hv, &[0.25 * hv[0]], 2e-10);
        let before = model.work();
        // Original total/elastic strain bound rejects this diluted design.
        assert!(matches!(
            c.value_gradient(&[0.1], Some(cx)),
            Err(PhysicsError::InvalidInput(_))
        ));
        let failed = model.work();
        assert_eq!(failed.solve_reservations, before.solve_reservations + 6);
        assert!(failed.started_solves > before.started_solves);
        assert_eq!(c.value_gradient(&x, Some(cx)).unwrap(), accepted);
        assert_eq!(model.work(), failed);
        model.set_budget(ComplianceNewtonBudget {
            solve_reservations: failed.solve_reservations,
            ..budget()
        });
        assert!(matches!(
            mass.value_gradient(&[0.6], Some(cx)),
            Err(PhysicsError::Budget { .. })
        ));
        assert_eq!(model.work(), failed);
        model.set_budget(budget());

        let (stage_pipeline, _, _, _, _) = thermal_fixture(0.03);
        let robust = RobustPipeline {
            filter: stage_pipeline.filter,
            params: SimpParams {
                beta: 2.0,
                ..pipeline.params
            },
            eta_offset: 0.05,
        };
        let stage = SimpNewtonModel::new_thermal_robust_with_resources(
            &robust,
            &elasticity,
            &[ThermalLoadCase {
                force: &retained_force,
                weight: 1.0,
                strains: TetThermalStrainField::Uniform(&retained_thermal),
            }],
            &[1.0 / 6.0],
            options(),
            resources,
        )
        .unwrap();
        let second_stage = stage.compliance("robust continuation stage").unwrap();
        second_stage.value_gradient(&[0.6], Some(cx)).unwrap();
        assert_eq!(stage.work(), model.work());
        assert!(model.work().started_solves > failed.started_solves);
        let before = model.work();
        assert_eq!(
            c.hessian_vector_product(&x, &[0.0], Some(cx)).unwrap(),
            [0.0]
        );
        assert_eq!(model.work(), before);
        gate.request();
        assert_eq!(c.value_gradient(&x, Some(cx)), Err(PhysicsError::Cancelled));
        assert_eq!(model.work(), before);
    });
}


fn transverse_thermal(material: &TetElasticMaterial) -> TetThermalStrainState {
    TetThermalStrainState::try_new(
        320.0, 300.0, [0.01, 0.0, 0.0, 0.0, 0.0, 0.0],
        material.material_state_identity, ContentHash([51; 32]), 0.1,
    ).unwrap()
}

#[test]
fn constrained_newton_limits_stress_in_an_independent_zero_weight_thermal_case() {
    use fs_topopt::elasticity::second_order::thermal::stress::ThermalStressOptions;
    with_cx(|_, cx| {
        let (pipeline, elasticity, material, axial, force) = thermal_fixture(0.1);
        let transverse = transverse_thermal(&material);
        let zero = vec![0.0; elasticity.n()];
        let loads = [
            ThermalLoadCase {
                force: &force, weight: 1.0,
                strains: TetThermalStrainField::Uniform(&axial),
            },
            ThermalLoadCase {
                force: &zero, weight: 0.0,
                strains: TetThermalStrainField::Uniform(&transverse),
            },
        ];
        let model = SimpNewtonModel::new_thermal(
            &pipeline, &elasticity, &loads, &[1.0 / 6.0], options(), budget(),
        ).unwrap();
        // The first case has physical sigma_zz=6F=.36 Pa, independent of m.
        // The second has restrained sigma_xx=-1.2m Pa and ZERO compliance.
        // Pick the allowance so the unnormalized fourth moment is one at m=.2.
        let maximum_m: f64 = 0.2;
        let allowance = (0.36_f64.powi(4) + (1.2 * maximum_m).powi(4)).sqrt().sqrt();
        let raw_c = model.compliance("thermal compliance").unwrap();
        let raw_s = model.thermal_stress_moment(
            ThermalStressOptions { allowable_pa: allowance, power: 4 },
            "physical thermal stress",
        ).unwrap();
        let chart = LogisticBox::new(&[0.05], &[0.95], Dims::NONE).unwrap();
        let c = chart.bind(&raw_c).unwrap();
        let stress = chart.bind(&raw_s).unwrap();
        let start = chart.encode(&[0.7], Some(cx)).unwrap();
        let scale = c.value_gradient(&start, Some(cx)).unwrap().value;
        let mut builder = ProblemBuilder::new();
        let x = builder.var("bounded raw occupancy", Manifold::Rn { dim: 1 }, Dims::NONE).unwrap();
        let cn = builder.pde_residual("thermal compliance", x, true, COMPLIANCE_DIMS).unwrap();
        let sn = builder.pde_residual("physical thermal stress", x, true, Dims::NONE).unwrap();
        let unit = builder.konst(scale, COMPLIANCE_DIMS).unwrap();
        let objective = builder.div(cn, unit).unwrap();
        builder.objective(objective, Sense::Minimize, 1.0).unwrap();
        let one = builder.konst(1.0, Dims::NONE).unwrap();
        let constraint = builder.sub(sn, one).unwrap();
        builder.constraint(constraint, ConstraintKind::LeZero, "all-case physical stress").unwrap();
        let problem = builder.finish();
        let oracle = ReverseProblem::new_with_physics(
            &problem,
            ReverseLimits { max_nodes: 128, max_scalar_slots: 1024 },
            &[PhysicsBinding { node: cn, model: &c }, PhysicsBinding { node: sn, model: &stress }],
            None,
        ).unwrap();
        let mut study = AugmentedNewtonStudy::new(
            &oracle, &start, AugmentedNewtonOptions::default(), Some(cx),
        ).unwrap();
        let result = study.run(
            1e-6, 300,
            AugmentedNewtonBudget {
                objective_evals: 600, hessian_products: 3000, multiplier_updates: 100,
            },
            Some(cx),
        ).unwrap();
        assert_eq!(result.stop, AugmentedNewtonStop::Converged, "{result:?}");
        assert!(result.kkt.within_tolerance(1e-6));
        let rho = chart.decode(&result.x, Some(cx)).unwrap()[0];
        let expected = ((maximum_m - 0.05) / 0.95).cbrt();
        close(&[rho], &[expected], 3e-4);
        let m = 0.05 + 0.95 * rho.powi(3);
        assert!(m < 0.25, "the unconstrained compliance optimum is m=.3");
        let physical_moment = (0.36_f64.powi(4) + (1.2 * m).powi(4)) / allowance.powi(4);
        close(&result.constraints, &[physical_moment - 1.0], 1e-8);
        assert!(result.constraints[0].abs() < 2e-6);
        assert!(result.multipliers[0] > 0.0);
        let mass = model.mass("material inventory").unwrap().value_gradient(&[rho], Some(cx)).unwrap();
        close(&[mass.value], &[rho], 1e-10);
    });
}

#[test]
fn thermal_stress_views_reuse_adjoints_and_keep_work_and_criterion_cache_transactional() {
    use fs_topopt::elasticity::second_order::thermal::stress::ThermalStressOptions;
    with_cx(|gate, cx| {
        let (pipeline, elasticity, material, axial, force) = thermal_fixture(0.03);
        let transverse = transverse_thermal(&material);
        let zero = vec![0.0; elasticity.n()];
        let loads = [
            ThermalLoadCase {
                force: &force, weight: 1.0,
                strains: TetThermalStrainField::Uniform(&axial),
            },
            ThermalLoadCase {
                force: &zero, weight: 0.0,
                strains: TetThermalStrainField::Uniform(&transverse),
            },
        ];
        let model = SimpNewtonModel::new_thermal(
            &pipeline, &elasticity, &loads, &[1.0 / 6.0], options(), budget(),
        ).unwrap();
        let stress = model.thermal_stress_moment(
            ThermalStressOptions { allowable_pa: 0.4, power: 4 }, "stress",
        ).unwrap();
        let second = model.thermal_stress_moment(
            ThermalStressOptions { allowable_pa: 0.8, power: 4 }, "second criterion",
        ).unwrap();
        let x = [0.6];
        let sample = stress.value_gradient(&x, Some(cx)).unwrap();
        let prepared = model.work();
        // Existing preparation: 2 loads + forward/compliance/volume/mass filters.
        // Additional stress: 2 retained adjoints + one transpose filter.
        assert_eq!(prepared.solve_reservations, 9);
        assert_eq!(prepared.started_solves, 9);
        assert_eq!(stress.signature().residual_dims, Dims::NONE);
        assert_eq!(stress.clone().value_gradient(&x, Some(cx)).unwrap(), sample);
        model.compliance("compliance").unwrap().value_gradient(&x, Some(cx)).unwrap();
        model.mass("mass").unwrap().value_gradient(&x, Some(cx)).unwrap();
        assert_eq!(model.work(), prepared);
        let hv = stress.hessian_vector_product(&x, &[0.2], Some(cx)).unwrap();
        assert_eq!(model.work().solve_reservations, prepared.solve_reservations + 6);
        assert_eq!(model.work().started_solves, prepared.started_solves + 6);
        let before = model.work();
        assert_eq!(stress.hessian_vector_product(&x, &[0.0], Some(cx)).unwrap(), [0.0]);
        assert_eq!(model.work(), before);
        model.set_budget(ComplianceNewtonBudget { solve_reservations: before.solve_reservations, ..budget() });
        assert!(matches!(second.value_gradient(&x, Some(cx)), Err(PhysicsError::Budget { .. })));
        assert_eq!(stress.value_gradient(&x, Some(cx)).unwrap(), sample);
        assert_eq!(model.work(), before);
        model.set_budget(budget());
        let rescaled = second.value_gradient(&x, Some(cx)).unwrap();
        close(&[rescaled.value], &[sample.value / 16.0], 1e-11);
        close(&rescaled.gradient, &[sample.gradient[0] / 16.0], 1e-10);
        assert_eq!(model.work().solve_reservations, before.solve_reservations + 3);
        let hs = second.hessian_vector_product(&x, &[0.2], Some(cx)).unwrap();
        close(&hs, &[hv[0] / 16.0], 1e-9);
        // A failed source-domain preparation must retain the previous point and
        // its stress adjoints, as it already does for compliance/mass.
        let before = model.work();
        let invalid = [0.1];
        assert!(second.value_gradient(&invalid, Some(cx)).is_err());
        let failed = model.work();
        assert!(failed.solve_reservations > before.solve_reservations);
        assert_eq!(second.value_gradient(&x, Some(cx)).unwrap(), rescaled);
        assert_eq!(model.work(), failed);
        gate.request();
        assert_eq!(second.value_gradient(&x, Some(cx)), Err(PhysicsError::Cancelled));
        assert_eq!(model.work(), failed);
    });
}
