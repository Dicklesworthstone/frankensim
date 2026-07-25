//! Transient conduction: the exact uniform-heating identity, steady-state
//! recovery against the crate's own steady path, observed temporal order, and
//! the refusal surface (bead `frankensim-extreal-program-f85xj.5.11`).

mod support;

use fs_conduction::ConductionError;
use fs_conduction::bc::{ThermalBc, ThermalBoundary, ThermalBoundaryBuilder};
use fs_conduction::field::ScalarField;
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::material::{ConductivityModel, ConductivityTable};
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::solve::{
    ConductionProblem, InitialGuess, LinearConfig, Nonlinearity, SolveConfig, StopRule, solve,
};
use fs_conduction::transient::{
    TransientConfig, TransientProblem, VolumetricHeatCapacity, assemble_capacitance, march,
};
use fs_rep_mesh::TetComplex;
use support::{with_cancelled_cx, with_cx};

const K: f64 = 10.0;
const RHO_CP: f64 = 2.0e6; // J/(m^3 K), a plausible solid
const T_COLD: f64 = 300.0;

fn unit_mesh(n: usize) -> ConductionMesh {
    let (complex, positions) = box_grid([n, n, n], [1.0, 1.0, 1.0]);
    let complex = TetComplex::from_tets(positions.len(), complex.tets);
    ConductionMesh::new(complex, positions).expect("unit box mesh")
}

fn linear_config() -> LinearConfig {
    LinearConfig {
        tolerance: 1e-13,
        max_iterations: 60_000,
        restart: 60,
    }
}

fn capacity() -> VolumetricHeatCapacity {
    VolumetricHeatCapacity::declared(RHO_CP).expect("capacity admits")
}

fn material() -> ConductivityModel {
    ConductivityModel::isotropic_declared(K).expect("material")
}

/// Entirely adiabatic: no Dirichlet anywhere, so every vertex is free.
fn adiabatic(mesh: &ConductionMesh) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .adiabatic_remainder()
        .finish()
        .expect("all-adiabatic partition")
}

/// Cold on `x = 0`, adiabatic elsewhere.
fn cold_wall(mesh: &ConductionMesh) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region(
            "cold",
            |face| on_box_face(face.centroid[0], 0.0),
            ThermalBc::dirichlet(T_COLD).expect("cold"),
        )
        .expect("cold region")
        .adiabatic_remainder()
        .finish()
        .expect("partition")
}

// ---------------------------------------------------------------------------
// The exact identity.
// ---------------------------------------------------------------------------

#[test]
fn uniform_adiabatic_heating_is_exact_for_every_scheme_and_step() {
    // With no boundary loss and a uniform source, the body heats uniformly at
    // dT/dt = f / (rho c_p). A spatially uniform field is in the kernel of the
    // stiffness operator, so the K terms vanish identically and the theta
    // method integrates C dT/dt = b EXACTLY — for any theta and any step.
    //
    // That makes this a machine-precision check on the capacitance assembly
    // and the stepping algebra together: if the mass matrix were lumped, or
    // scaled wrongly, or the theta weighting misapplied, the rate would drift.
    let mesh = unit_mesh(3);
    let boundary = adiabatic(&mesh);
    let flux_density = 4.0e5; // W/m^3
    let source = ScalarField::uniform("volumetric source", flux_density).expect("source");
    let expected_rate = flux_density / RHO_CP;

    for (theta, dt, steps) in [(1.0, 0.5, 8), (0.5, 0.5, 8), (1.0, 4.0, 3), (0.75, 1.5, 5)] {
        let config = TransientConfig::new(theta, dt, linear_config()).expect("config");
        let initial = vec![T_COLD; mesh.vertex_count()];
        let solution = with_cx(|cx| {
            march(
                cx,
                TransientProblem {
                    mesh: &mesh,
                    boundary: &boundary,
                    material: &material(),
                    source: &source,
                    capacity: capacity(),
                },
                &config,
                &initial,
                steps,
            )
            .expect("march")
        });

        let elapsed = dt * (steps as f64);
        let expected = T_COLD + expected_rate * elapsed;
        for (vertex, value) in solution.temperature.iter().enumerate() {
            assert!(
                (value - expected).abs() < 1e-8,
                "theta={theta} dt={dt}: vertex {vertex} reached {value}, expected {expected}"
            );
        }
        assert!((solution.time_s - elapsed).abs() < 1e-12);
        assert_eq!(solution.steps.len(), steps);
    }
}

#[test]
fn the_capacitance_matrix_holds_the_exact_heat_content() {
    // Summing every entry of C gives rho c_p V: the discrete heat content of
    // a unit uniform temperature rise. This is the property that makes the
    // uniform-heating identity above exact, pinned directly.
    let mesh = unit_mesh(3);
    let capacitance = with_cx(|cx| assemble_capacitance(cx, &mesh, capacity()).expect("assemble"));

    let ones = vec![1.0f64; mesh.vertex_count()];
    let mut product = vec![0.0f64; mesh.vertex_count()];
    capacitance.spmv(&ones, &mut product);
    let total: f64 = product.iter().sum();

    let volume: f64 = (0..mesh.element_count())
        .map(|e| mesh.element_volume(e))
        .sum();
    let expected = RHO_CP * volume;
    assert!(
        (total - expected).abs() < 1e-6 * expected,
        "total heat capacity {total} != rho*c_p*V {expected}"
    );
}

// ---------------------------------------------------------------------------
// Steady-state recovery: cross-check against the crate's own steady path.
// ---------------------------------------------------------------------------

#[test]
fn marching_to_large_time_recovers_the_steady_solution() {
    // The transient must settle onto exactly what the independently tested
    // steady solver produces. This is the cross-check that keeps the new path
    // from drifting onto its own physics.
    let mesh = unit_mesh(3);
    let boundary = cold_wall(&mesh);
    let source = ScalarField::uniform("volumetric source", 5.0e4).expect("source");

    let steady = with_cx(|cx| {
        solve(
            cx,
            ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material(),
                source: &source,
            },
            SolveConfig {
                nonlinearity: Nonlinearity::FixedPoint {
                    relaxation: 1.0,
                    max_backtracks: 8,
                },
                stop: StopRule {
                    residual_rtol: 1e-12,
                    residual_atol: 1e-24,
                    step_atol: 0.0,
                    max_iterations: 12,
                },
                linear: linear_config(),
                initial: InitialGuess::DirichletMean,
            },
        )
        .expect("steady solve")
    });

    // The diffusive time scale is rho c_p L^2 / k = 2e6 / 10 = 2e5 s; march
    // well past it with a coarse backward-Euler step, which damps toward the
    // fixed point rather than resolving the path.
    let config = TransientConfig::backward_euler(5.0e4, linear_config()).expect("config");
    let initial = vec![T_COLD; mesh.vertex_count()];
    let transient = with_cx(|cx| {
        march(
            cx,
            TransientProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material(),
                source: &source,
                capacity: capacity(),
            },
            &config,
            &initial,
            200,
        )
        .expect("march")
    });

    let mut worst = 0.0f64;
    for (a, b) in transient.temperature.iter().zip(steady.temperature.iter()) {
        worst = worst.max((a - b).abs());
    }
    assert!(
        worst < 1e-6,
        "transient settled {worst} K away from the steady solution"
    );
}

#[test]
fn an_initial_condition_inconsistent_with_the_boundary_is_corrected_not_carried() {
    // A caller handing in a uniform field that disagrees with the Dirichlet
    // wall must not have that disagreement persist: the prescribed vertices
    // are lifted at every step.
    let mesh = unit_mesh(2);
    let boundary = cold_wall(&mesh);
    let source = ScalarField::uniform("volumetric source", 0.0).expect("source");
    let config = TransientConfig::backward_euler(1.0e3, linear_config()).expect("config");

    let initial = vec![T_COLD + 50.0; mesh.vertex_count()];
    let solution = with_cx(|cx| {
        march(
            cx,
            TransientProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material(),
                source: &source,
                capacity: capacity(),
            },
            &config,
            &initial,
            1,
        )
        .expect("march")
    });

    for (vertex, position) in mesh.positions().iter().enumerate() {
        if on_box_face(position[0], 0.0) {
            assert!(
                (solution.temperature[vertex] - T_COLD).abs() < 1e-9,
                "Dirichlet vertex {vertex} holds {} not {T_COLD}",
                solution.temperature[vertex]
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Observed temporal order (self-convergence).
// ---------------------------------------------------------------------------

/// Richardson triple: three solutions at `dt`, `dt/2`, `dt/4` give the
/// observed order without needing an external reference solution.
fn observed_order(theta: f64) -> f64 {
    let mesh = unit_mesh(2);
    let boundary = cold_wall(&mesh);
    let source = ScalarField::uniform("volumetric source", 2.0e5).expect("source");
    let final_time = 4.0e4;

    let at = |steps: usize| -> Vec<f64> {
        let config =
            TransientConfig::new(theta, final_time / (steps as f64), linear_config()).expect("cfg");
        let initial = vec![T_COLD; mesh.vertex_count()];
        with_cx(|cx| {
            march(
                cx,
                TransientProblem {
                    mesh: &mesh,
                    boundary: &boundary,
                    material: &material(),
                    source: &source,
                    capacity: capacity(),
                },
                &config,
                &initial,
                steps,
            )
            .expect("march")
            .temperature
        })
    };

    let coarse = at(8);
    let medium = at(16);
    let fine = at(32);

    let diff = |a: &[f64], b: &[f64]| -> f64 {
        a.iter()
            .zip(b.iter())
            .fold(0.0f64, |acc, (x, y)| acc.max((x - y).abs()))
    };
    let first = diff(&coarse, &medium);
    let second = diff(&medium, &fine);
    assert!(
        second > 0.0,
        "the refinement must actually change the answer"
    );
    (first / second).log2()
}

#[test]
fn backward_euler_is_first_order_and_crank_nicolson_is_second() {
    // Self-convergence, not comparison against an analytic solution: the
    // semi-discrete system has no elementary closed form on a tet mesh, so
    // the honest measurement is a Richardson triple. Bands are deliberately
    // loose — the claim is "these are different orders and each is near its
    // nominal one", not a tight order estimate.
    let euler = observed_order(1.0);
    assert!(
        (0.8..1.3).contains(&euler),
        "backward Euler observed order {euler}, expected near 1"
    );

    let nicolson = observed_order(0.5);
    assert!(
        (1.7..2.3).contains(&nicolson),
        "Crank-Nicolson observed order {nicolson}, expected near 2"
    );
    assert!(
        nicolson > euler + 0.5,
        "Crank-Nicolson ({nicolson}) must be measurably higher order than backward Euler ({euler})"
    );
}

#[test]
fn the_nominal_order_matches_the_declared_scheme() {
    assert_eq!(
        TransientConfig::backward_euler(1.0, linear_config())
            .expect("be")
            .nominal_order(),
        1
    );
    assert_eq!(
        TransientConfig::crank_nicolson(1.0, linear_config())
            .expect("cn")
            .nominal_order(),
        2
    );
}

// ---------------------------------------------------------------------------
// Refusals.
// ---------------------------------------------------------------------------

#[test]
fn a_temperature_dependent_conductivity_is_refused_not_linearized() {
    // Freezing k(T) across a step is a different scheme with its own error
    // behaviour. Adopting it silently would make the observed time order
    // depend on how strongly k varies.
    let mesh = unit_mesh(2);
    let boundary = cold_wall(&mesh);
    let source = ScalarField::uniform("volumetric source", 0.0).expect("source");
    let curve = ConductivityTable::declared_curve(vec![(280.0, 9.0), (360.0, 12.0)])
        .expect("temperature-dependent table");
    let varying = ConductivityModel::isotropic(curve);
    assert!(varying.is_temperature_dependent());

    let config = TransientConfig::backward_euler(1.0, linear_config()).expect("config");
    let initial = vec![T_COLD; mesh.vertex_count()];
    let error = with_cx(|cx| {
        march(
            cx,
            TransientProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &varying,
                source: &source,
                capacity: capacity(),
            },
            &config,
            &initial,
            1,
        )
        .expect_err("k(T) must refuse")
    });
    assert!(matches!(
        error,
        ConductionError::Config {
            parameter: "conductivity",
            ..
        }
    ));
}

#[test]
fn configuration_refuses_unstable_and_degenerate_declarations() {
    for theta in [0.49, 0.0, -1.0, 1.5, f64::NAN] {
        assert!(
            TransientConfig::new(theta, 1.0, linear_config()).is_err(),
            "theta {theta} must refuse: below 0.5 the scheme is only conditionally stable"
        );
    }
    for dt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            TransientConfig::new(1.0, dt, linear_config()).is_err(),
            "dt {dt} must refuse"
        );
    }
    for capacity in [0.0, -1.0, f64::NAN] {
        assert!(
            VolumetricHeatCapacity::declared(capacity).is_err(),
            "capacity {capacity} must refuse: a zero capacity is a steady problem, not a fast material"
        );
    }
    assert!(VolumetricHeatCapacity::declared(RHO_CP).is_ok());
}

#[test]
fn a_mismatched_initial_vector_refuses() {
    let mesh = unit_mesh(2);
    let boundary = cold_wall(&mesh);
    let source = ScalarField::uniform("volumetric source", 0.0).expect("source");
    let config = TransientConfig::backward_euler(1.0, linear_config()).expect("config");
    let error = with_cx(|cx| {
        march(
            cx,
            TransientProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material(),
                source: &source,
                capacity: capacity(),
            },
            &config,
            &[T_COLD, T_COLD],
            1,
        )
        .expect_err("short initial vector")
    });
    assert!(matches!(
        error,
        ConductionError::FieldLength {
            field: "initial temperature",
            ..
        }
    ));
}

// ---------------------------------------------------------------------------
// Determinism and cancellation.
// ---------------------------------------------------------------------------

#[test]
fn marching_is_deterministic_across_independent_runs() {
    let build = || {
        let mesh = unit_mesh(3);
        let boundary = cold_wall(&mesh);
        let source = ScalarField::uniform("volumetric source", 1.0e5).expect("source");
        let config = TransientConfig::crank_nicolson(2.0e3, linear_config()).expect("config");
        let initial = vec![T_COLD; mesh.vertex_count()];
        with_cx(|cx| {
            march(
                cx,
                TransientProblem {
                    mesh: &mesh,
                    boundary: &boundary,
                    material: &material(),
                    source: &source,
                    capacity: capacity(),
                },
                &config,
                &initial,
                6,
            )
            .expect("march")
        })
    };
    assert_eq!(build(), build());
}

#[test]
fn a_cancelled_march_publishes_no_partial_solution() {
    let mesh = unit_mesh(2);
    let boundary = cold_wall(&mesh);
    let source = ScalarField::uniform("volumetric source", 0.0).expect("source");
    let config = TransientConfig::backward_euler(1.0, linear_config()).expect("config");
    let initial = vec![T_COLD; mesh.vertex_count()];
    let error = with_cancelled_cx(|cx| {
        march(
            cx,
            TransientProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material(),
                source: &source,
                capacity: capacity(),
            },
            &config,
            &initial,
            4,
        )
        .expect_err("a cancelled march publishes nothing")
    });
    assert!(matches!(error, ConductionError::Cancelled { .. }));
}

// ---------------------------------------------------------------------------
// The transient adjoint.
// ---------------------------------------------------------------------------

use fs_adjoint::verify_gradient;
use fs_conduction::transient::{FinalStateFunctional, source_scale_gradient};

/// Forward objective as a function of the source scale, for the independent
/// finite-difference check.
fn objective_at(
    mesh: &ConductionMesh,
    boundary: &ThermalBoundary,
    base_density: f64,
    functional: &FinalStateFunctional,
    config: &TransientConfig,
    steps: usize,
    scale: f64,
) -> f64 {
    let source =
        ScalarField::uniform("volumetric source", base_density * scale).expect("scaled source");
    let initial = vec![T_COLD; mesh.vertex_count()];
    let solution = with_cx(|cx| {
        march(
            cx,
            TransientProblem {
                mesh,
                boundary,
                material: &material(),
                source: &source,
                capacity: capacity(),
            },
            config,
            &initial,
            steps,
        )
        .expect("march")
    });
    functional.evaluate(&solution.temperature)
}

#[test]
fn the_transient_adjoint_passes_the_crate_gradient_gate() {
    // Finite differences are the INDEPENDENT check, through the same
    // `fs_adjoint::verify_gradient` gate the steady adjoint uses — not a
    // bespoke comparison whose tolerance I could tune until it passed.
    let mesh = unit_mesh(2);
    let boundary = cold_wall(&mesh);
    let base_density = 2.0e5;
    let steps = 6;
    let config = TransientConfig::backward_euler(2.0e4, linear_config()).expect("config");

    // Probe the vertex farthest from the cold wall, where the source actually
    // moves the answer.
    let probe = mesh
        .positions()
        .iter()
        .enumerate()
        .max_by(|a, b| a.1[0].total_cmp(&b.1[0]))
        .map(|(index, _)| index)
        .expect("a vertex");
    let functional =
        FinalStateFunctional::probe(mesh.vertex_count(), probe).expect("probe functional");

    let source = ScalarField::uniform("volumetric source", base_density).expect("source");
    let gradient = with_cx(|cx| {
        source_scale_gradient(
            cx,
            TransientProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material(),
                source: &source,
                capacity: capacity(),
            },
            &config,
            &functional,
            steps,
        )
        .expect("adjoint gradient")
    });
    assert!(
        gradient > 0.0,
        "more source must raise the probe temperature, got {gradient}"
    );

    let objective = |p: &[f64]| -> f64 {
        objective_at(
            &mesh,
            &boundary,
            base_density,
            &functional,
            &config,
            steps,
            p[0],
        )
    };
    let verdict = verify_gradient(&objective, &[1.0], &[gradient], &[vec![1.0]], 1e-4, 1e-6);
    assert!(
        verdict.pass,
        "the transient adjoint must agree with finite differences: max_rel_err={:e} pairs={:?}",
        verdict.max_rel_err, verdict.pairs
    );
    assert_eq!(
        verdict.informative_directions, 1,
        "the probe direction must carry signal, else the pass is vacuous"
    );
}

#[test]
fn the_adjoint_agrees_with_finite_differences_under_crank_nicolson_too() {
    // The adjoint is derived from the theta-method operators, so it must track
    // theta rather than being right only for backward Euler.
    let mesh = unit_mesh(2);
    let boundary = cold_wall(&mesh);
    let base_density = 1.5e5;
    let steps = 5;
    let config = TransientConfig::crank_nicolson(1.5e4, linear_config()).expect("config");
    let functional = FinalStateFunctional::new(vec![
        1.0 / (mesh.vertex_count() as f64);
        mesh.vertex_count()
    ])
    .expect("mean functional");

    let source = ScalarField::uniform("volumetric source", base_density).expect("source");
    let gradient = with_cx(|cx| {
        source_scale_gradient(
            cx,
            TransientProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material(),
                source: &source,
                capacity: capacity(),
            },
            &config,
            &functional,
            steps,
        )
        .expect("adjoint gradient")
    });

    let objective = |p: &[f64]| -> f64 {
        objective_at(
            &mesh,
            &boundary,
            base_density,
            &functional,
            &config,
            steps,
            p[0],
        )
    };
    let verdict = verify_gradient(&objective, &[1.0], &[gradient], &[vec![1.0]], 1e-4, 1e-6);
    assert!(
        verdict.pass,
        "Crank-Nicolson adjoint disagreed with finite differences: max_rel_err={:e} pairs={:?}",
        verdict.max_rel_err, verdict.pairs
    );
    assert_eq!(verdict.informative_directions, 1);
}

#[test]
fn the_adjoint_is_exact_for_the_uniform_adiabatic_case() {
    // The one case with a closed form: with no loss, the mean temperature
    // rises by (s * f / rho c_p) * t, so d(mean)/ds = f * t / (rho c_p)
    // exactly. Any drift here is the adjoint algebra, not discretization.
    let mesh = unit_mesh(2);
    let boundary = adiabatic(&mesh);
    let base_density = 3.0e5;
    let dt = 1.0;
    let steps = 7;
    let config = TransientConfig::backward_euler(dt, linear_config()).expect("config");
    let functional = FinalStateFunctional::new(vec![
        1.0 / (mesh.vertex_count() as f64);
        mesh.vertex_count()
    ])
    .expect("mean functional");

    let source = ScalarField::uniform("volumetric source", base_density).expect("source");
    let gradient = with_cx(|cx| {
        source_scale_gradient(
            cx,
            TransientProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material(),
                source: &source,
                capacity: capacity(),
            },
            &config,
            &functional,
            steps,
        )
        .expect("adjoint gradient")
    });

    let expected = base_density * dt * (steps as f64) / RHO_CP;
    assert!(
        (gradient - expected).abs() < 1e-9 * expected,
        "adjoint gradient {gradient} != closed form {expected}"
    );
}

#[test]
fn the_adjoint_refuses_a_temperature_dependent_material() {
    // The state-independent-operator argument that removes checkpointing does
    // not hold for k(T), so the adjoint must refuse rather than return a
    // gradient whose derivation no longer applies.
    let mesh = unit_mesh(2);
    let boundary = cold_wall(&mesh);
    let source = ScalarField::uniform("volumetric source", 1.0e5).expect("source");
    let curve = ConductivityTable::declared_curve(vec![(280.0, 9.0), (360.0, 12.0)])
        .expect("varying table");
    let varying = ConductivityModel::isotropic(curve);
    let config = TransientConfig::backward_euler(1.0e3, linear_config()).expect("config");
    let functional = FinalStateFunctional::probe(mesh.vertex_count(), 0).expect("probe functional");

    let error = with_cx(|cx| {
        source_scale_gradient(
            cx,
            TransientProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &varying,
                source: &source,
                capacity: capacity(),
            },
            &config,
            &functional,
            2,
        )
        .expect_err("k(T) must refuse")
    });
    assert!(matches!(
        error,
        ConductionError::Config {
            parameter: "conductivity",
            ..
        }
    ));
}

#[test]
fn functional_admission_refuses_degenerate_declarations() {
    assert!(FinalStateFunctional::new(vec![1.0, f64::NAN]).is_err());
    assert!(FinalStateFunctional::probe(4, 4).is_err(), "out of range");
    assert!(FinalStateFunctional::probe(4, 3).is_ok());

    let mesh = unit_mesh(2);
    let boundary = cold_wall(&mesh);
    let source = ScalarField::uniform("volumetric source", 1.0).expect("source");
    let config = TransientConfig::backward_euler(1.0, linear_config()).expect("config");
    let short = FinalStateFunctional::new(vec![1.0, 1.0]).expect("short functional");
    let error = with_cx(|cx| {
        source_scale_gradient(
            cx,
            TransientProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material(),
                source: &source,
                capacity: capacity(),
            },
            &config,
            &short,
            1,
        )
        .expect_err("mismatched weights")
    });
    assert!(matches!(
        error,
        ConductionError::FieldLength {
            field: "functional weights",
            ..
        }
    ));
}

// ---------------------------------------------------------------------------
// G1: the fs-vvreg Level-A lumped-transient binding.
// ---------------------------------------------------------------------------

use fs_vvreg::thermal_level_a::{ThermalLevelAKind, thermal_level_a_cases};

/// The catalog row this fixture answers: theta/theta0 = exp(-t/tau) at
/// t/tau = 1, admitted only for Biot <= 0.1.
fn lumped_reference() -> (f64, f64) {
    let case = thermal_level_a_cases()
        .iter()
        .find(|case| case.id == "thermal-a-lumped-transient")
        .expect("Level-A lumped row");
    assert_eq!(case.kind, ThermalLevelAKind::AnalyticReference);
    let biot_ceiling = case
        .context
        .iter()
        .find(|entry| entry.name == "biot-number")
        .map(|entry| entry.hi)
        .expect("the row declares a Biot ceiling");
    (case.reference_value_si, biot_ceiling)
}

/// Cool a unit cube from `T0` through Robin faces at the requested Biot
/// number, and return the volume-mean normalized excess at exactly one time
/// constant.
///
/// For a unit cube, `V = 1`, `A = 6`, so `Lc = V/A = 1/6`,
/// `h = 6 k Bi`, and `tau = rho c_p V / (h A)`.
fn normalized_excess_at_one_time_constant(biot: f64, n: usize) -> f64 {
    let mesh = unit_mesh(n);
    let htc = 6.0 * K * biot;
    let area = 6.0f64;
    let tau = RHO_CP / (htc * area);

    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "ambient",
            |_face| true,
            ThermalBc::robin(htc, T_COLD).expect("robin"),
        )
        .expect("robin region")
        .adiabatic_remainder()
        .finish()
        .expect("partition");

    let excess0 = 100.0;
    let initial = vec![T_COLD + excess0; mesh.vertex_count()];
    let source = ScalarField::uniform("volumetric source", 0.0).expect("no source");
    // Resolve time finely so the residual discrepancy is the LUMPED
    // approximation, not the integrator.
    let steps = 400;
    let config =
        TransientConfig::crank_nicolson(tau / (steps as f64), linear_config()).expect("config");

    let solution = with_cx(|cx| {
        march(
            cx,
            TransientProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material(),
                source: &source,
                capacity: capacity(),
            },
            &config,
            &initial,
            steps,
        )
        .expect("march")
    });

    let mean: f64 = solution.temperature.iter().sum::<f64>() / (mesh.vertex_count() as f64);
    (mean - T_COLD) / excess0
}

#[test]
fn the_lumped_decay_matches_the_level_a_row_inside_its_declared_biot_regime() {
    // The catalog row is an approximation with a DECLARED validity context,
    // so the honest binding is not "agrees within a tolerance I chose" but
    // "agrees within the regime the row declares".
    let (reference, biot_ceiling) = lumped_reference();
    assert!(
        (reference - (-1.0f64).exp()).abs() < 1e-12,
        "the row is exp(-1)"
    );

    let at_ceiling = normalized_excess_at_one_time_constant(biot_ceiling, 3);
    let error = (at_ceiling - reference).abs() / reference;
    assert!(
        error < 0.05,
        "at the declared Biot ceiling {biot_ceiling} the solve gave {at_ceiling} \
         against the row's {reference} ({:.2}% off)",
        error * 100.0
    );
}

#[test]
fn the_lumped_discrepancy_is_controlled_by_the_biot_number() {
    // This is what makes the binding evidence rather than a coincidence: the
    // lumped model ignores the internal gradient, so its error must SHRINK as
    // Biot does. A fixture that merely landed inside a tolerance at one Biot
    // could be passing for the wrong reason.
    let (reference, biot_ceiling) = lumped_reference();

    let coarse_regime = (normalized_excess_at_one_time_constant(biot_ceiling, 3) - reference).abs();
    let deep_regime =
        (normalized_excess_at_one_time_constant(biot_ceiling / 4.0, 3) - reference).abs();

    assert!(
        deep_regime < coarse_regime,
        "the lumped discrepancy must shrink with Biot: {deep_regime:e} at Bi/4 \
         is not below {coarse_regime:e} at Bi"
    );
    // Roughly first order in Biot, so a 4x reduction should buy most of a 4x
    // error reduction. Loose, because the constant is geometry dependent.
    assert!(
        deep_regime < coarse_regime * 0.6,
        "quartering Biot barely moved the discrepancy ({coarse_regime:e} -> {deep_regime:e}); \
         that suggests the error is NOT the lumped approximation"
    );
}
