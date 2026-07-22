use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::{
    ConductionMesh, ConductionProblem, ConductivityModel, InitialGuess, LinearConfig, Nonlinearity,
    ScalarField, SolveConfig, StopRule, ThermalBc, ThermalBoundaryBuilder, solve,
};
use fs_convection::{
    CorrelationId, CorrelationInputs, ThermalConductivity, correlation_catalog, evaluate,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_qty::{Length, Temperature};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x0000_C0DE_C011_0000,
                kernel_id: 52,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn config() -> SolveConfig {
    SolveConfig {
        nonlinearity: Nonlinearity::FixedPoint {
            relaxation: 1.0,
            max_backtracks: 4,
        },
        stop: StopRule {
            residual_rtol: 1.0e-11,
            residual_atol: 1.0e-24,
            step_atol: 0.0,
            max_iterations: 8,
        },
        linear: LinearConfig {
            tolerance: 1.0e-13,
            max_iterations: 20_000,
            restart: 40,
        },
        initial: InitialGuess::DirichletMean,
    }
}

#[test]
fn three_flow_rate_heatsink_slab_keeps_correlation_evidence_at_the_robin_seam() {
    const SOLID_LENGTH: f64 = 0.010;
    const WIDTH: f64 = 0.020;
    const HEIGHT: f64 = 0.020;
    const SOLID_K: f64 = 200.0;
    const T_BASE: f64 = 350.0;
    const T_AMBIENT: f64 = 300.0;
    const FLUID_K: f64 = 0.026;
    const PLATE_LENGTH: f64 = 0.050;

    let model_names = correlation_catalog()
        .into_iter()
        .map(|card| card.model.name)
        .collect::<std::collections::BTreeSet<_>>();
    let mut prior_h = 0.0f64;
    let mut prior_heat = 0.0f64;

    for reynolds in [10_000.0, 40_000.0, 160_000.0] {
        let evaluated = evaluate(
            CorrelationId::FlatPlateLaminarAverage,
            CorrelationInputs::forced(reynolds, 0.7),
        )
        .expect("flow point is inside the laminar plate card");
        let coupling = evaluated
            .robin_boundary(
                ThermalConductivity::new(FLUID_K),
                Length::new(PLATE_LENGTH),
                Temperature::new(T_AMBIENT),
            )
            .expect("correlation lowers to Robin");
        let h = coupling.coefficient().value.value();
        assert!(model_names.contains(&coupling.coefficient().model.cards[0]));
        assert_eq!(
            coupling.coefficient().model.cards,
            [CorrelationId::FlatPlateLaminarAverage.name()]
        );
        assert!(coupling.coefficient().model.in_domain);
        match coupling.boundary_condition() {
            ThermalBc::Robin { htc, t_ref } => {
                assert_eq!(htc.at(0).to_bits(), h.to_bits());
                assert_eq!(t_ref.at(0).to_bits(), T_AMBIENT.to_bits());
            }
            other => panic!("correlation lowered to {other:?}, not Robin"),
        }

        let (complex, positions) = box_grid([4, 2, 2], [SOLID_LENGTH, WIDTH, HEIGHT]);
        let mesh = ConductionMesh::new(complex, positions).expect("heatsink slab mesh");
        let material = ConductivityModel::isotropic_declared(SOLID_K).expect("solid material");
        let source = ScalarField::Uniform(0.0);
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region(
                "base",
                |face| on_box_face(face.centroid[0], 0.0),
                ThermalBc::dirichlet(T_BASE).expect("base"),
            )
            .expect("base region")
            .region(
                "correlation-robin",
                |face| on_box_face(face.centroid[0], SOLID_LENGTH),
                coupling.boundary_condition().clone(),
            )
            .expect("Robin region")
            .adiabatic_remainder()
            .finish()
            .expect("boundary partition");
        let solution = with_cx(|cx| {
            solve(
                cx,
                ConductionProblem {
                    mesh: &mesh,
                    boundary: &boundary,
                    material: &material,
                    source: &source,
                },
                config(),
            )
            .expect("conduction solve")
        });

        let area = WIDTH * HEIGHT;
        let analytic_heat =
            area * SOLID_K * h * (T_BASE - T_AMBIENT) / (SOLID_K + h * SOLID_LENGTH);
        let heat = solution.report.energy.robin_out_w;
        let relative_error = (heat - analytic_heat).abs() / analytic_heat;
        assert!(
            relative_error < 1.0e-8,
            "Re={reynolds} h={h} Q={heat} analytic={analytic_heat} rel={relative_error}"
        );
        assert!(solution.report.energy.relative_closure() < 1.0e-10);
        assert!(h > prior_h);
        assert!(heat > prior_heat);
        prior_h = h;
        prior_heat = heat;

        println!(
            "{{\"suite\":\"fs-convection/robin-seam\",\"Re\":{reynolds},\"Pr\":0.7,\"card\":\"{}\",\"h_w_m2k\":{h},\"in_domain\":true,\"model_discrepancy_rel\":{},\"robin_heat_w\":{heat},\"analytic_heat_w\":{analytic_heat},\"relative_error\":{relative_error},\"energy_closure_rel\":{}}}",
            CorrelationId::FlatPlateLaminarAverage.name(),
            coupling.coefficient().model.discrepancy_rel,
            solution.report.energy.relative_closure(),
        );
    }
}
