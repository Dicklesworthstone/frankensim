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
