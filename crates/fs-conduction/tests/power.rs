//! Component power maps: exactness, superposition, audit refusal, and the
//! analytic uniform-source cross-check (bead
//! `frankensim-extreal-program-f85xj.5.9`).
//!
//! The load-bearing claim is that the power the discrete operator injects
//! equals the power the caller declared IDENTICALLY — not to discretization
//! order. Most tests here are therefore written at `1e-12` relative, and they
//! would fail if the projection ever went back to smearing a per-element
//! density onto nodes.

mod support;

use fs_conduction::bc::{ThermalBc, ThermalBoundaryBuilder};
use fs_conduction::field::ScalarField;
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::material::ConductivityModel;
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::solve::{
    ConductionProblem, InitialGuess, LinearConfig, Nonlinearity, SolveConfig, StopRule, solve,
};
use fs_conduction::{
    ComponentPower, ConductionError, DEFAULT_POWER_TOLERANCE, PowerMap, PowerUncertainty,
};
use fs_rep_mesh::TetComplex;
use support::with_cx;

fn unit_mesh(n: usize) -> ConductionMesh {
    let (complex, positions) = box_grid([n, n, n], [1.0, 1.0, 1.0]);
    let complex = TetComplex::from_tets(positions.len(), complex.tets);
    ConductionMesh::new(complex, positions).expect("unit box mesh")
}

fn config() -> SolveConfig {
    SolveConfig {
        nonlinearity: Nonlinearity::FixedPoint {
            relaxation: 1.0,
            max_backtracks: 8,
        },
        stop: StopRule {
            residual_rtol: 1e-11,
            residual_atol: 1e-24,
            step_atol: 0.0,
            max_iterations: 12,
        },
        linear: LinearConfig {
            tolerance: 1e-13,
            max_iterations: 60_000,
            restart: 60,
        },
        initial: InitialGuess::DirichletMean,
    }
}

/// Every vertex of the mesh, i.e. a component covering the whole domain.
fn all_vertices(mesh: &ConductionMesh) -> Vec<usize> {
    (0..mesh.vertex_count()).collect()
}

/// Vertices on the `x = 0` face — a compact footprint that shares elements
/// with the rest of the mesh, which is where a naive projection leaks power.
fn low_x_vertices(mesh: &ConductionMesh) -> Vec<usize> {
    mesh.positions()
        .iter()
        .enumerate()
        .filter(|(_, p)| on_box_face(p[0], 0.0))
        .map(|(index, _)| index)
        .collect()
}

// ---------------------------------------------------------------------------
// Exactness: the whole point of the module.
// ---------------------------------------------------------------------------

#[test]
fn delivered_power_equals_declared_power_identically() {
    let mesh = unit_mesh(3);
    let component = ComponentPower::new(
        "cpu-die",
        42.5,
        PowerUncertainty::Unstated,
        low_x_vertices(&mesh),
    )
    .expect("component admits");
    let map = PowerMap::new(vec![component], 42.5).expect("map admits");

    let (_, audit) = map
        .volumetric_source(&mesh, DEFAULT_POWER_TOLERANCE)
        .expect("projection admits");

    assert!(
        (audit.delivered_total_w() - 42.5).abs() < 1e-12,
        "delivered {} != declared 42.5; a partial-footprint component must not \
         leak power into the elements it shares",
        audit.delivered_total_w()
    );
    assert!((audit.component_total_w() - 42.5).abs() < 1e-12);
}

#[test]
fn delivered_power_is_exact_across_mesh_refinement() {
    // If the projection were an approximation, this would converge with h.
    // It must instead be exact at every resolution — that is the difference
    // between a discretization property and an arithmetic identity.
    for n in [2usize, 3, 4, 5] {
        let mesh = unit_mesh(n);
        let component = ComponentPower::new(
            "die",
            7.25,
            PowerUncertainty::Unstated,
            low_x_vertices(&mesh),
        )
        .expect("component admits");
        let map = PowerMap::new(vec![component], 7.25).expect("map admits");
        let (_, audit) = map
            .volumetric_source(&mesh, DEFAULT_POWER_TOLERANCE)
            .expect("projection admits");
        assert!(
            (audit.delivered_total_w() - 7.25).abs() < 1e-12,
            "n={n}: delivered {} drifted from 7.25",
            audit.delivered_total_w()
        );
    }
}

#[test]
fn overlapping_components_superpose_without_losing_power() {
    // Two dies under one spreader share vertices. Overlap is admitted, and
    // the delivered total must still be the sum: the delivered functional is
    // linear in the nodal field, so nothing cancels or double-counts.
    let mesh = unit_mesh(3);
    let shared = low_x_vertices(&mesh);
    let a = ComponentPower::new("die-a", 10.0, PowerUncertainty::Unstated, shared.clone())
        .expect("a admits");
    let b =
        ComponentPower::new("die-b", 15.0, PowerUncertainty::Unstated, shared).expect("b admits");
    let map = PowerMap::new(vec![a, b], 25.0).expect("map admits");

    let (_, audit) = map
        .volumetric_source(&mesh, DEFAULT_POWER_TOLERANCE)
        .expect("projection admits");
    assert!(
        (audit.delivered_total_w() - 25.0).abs() < 1e-12,
        "overlapping components delivered {}",
        audit.delivered_total_w()
    );
    assert_eq!(audit.rows().len(), 2);
}

#[test]
fn a_duplicated_vertex_does_not_take_a_double_share() {
    // A vertex listed twice would otherwise receive twice its share while the
    // declared/summed totals still balanced — a silent misdistribution.
    let mesh = unit_mesh(2);
    let mut doubled = low_x_vertices(&mesh);
    doubled.extend_from_slice(&doubled.clone());
    let component =
        ComponentPower::new("die", 3.0, PowerUncertainty::Unstated, doubled).expect("admits");
    assert_eq!(
        component.vertices().len(),
        low_x_vertices(&mesh).len(),
        "duplicate bindings must collapse"
    );

    let map = PowerMap::new(vec![component], 3.0).expect("map admits");
    let (_, audit) = map
        .volumetric_source(&mesh, DEFAULT_POWER_TOLERANCE)
        .expect("projection admits");
    assert!((audit.delivered_total_w() - 3.0).abs() < 1e-12);
}

// ---------------------------------------------------------------------------
// The analytic cross-check: a whole-domain map must equal the uniform source.
// ---------------------------------------------------------------------------

#[test]
fn a_whole_domain_map_reproduces_the_equivalent_uniform_source_solve() {
    // A component covering every vertex with total power P over volume V is
    // the same physics as ScalarField::Uniform(P/V). Solving both and
    // comparing pins the power path to the crate's existing, independently
    // tested source path instead of leaving it on private arithmetic.
    let mesh = unit_mesh(3);
    let total_power = 12.0;
    let volume = 1.0; // unit box

    let map = PowerMap::new(
        vec![
            ComponentPower::new(
                "whole-domain",
                total_power,
                PowerUncertainty::Unstated,
                all_vertices(&mesh),
            )
            .expect("component admits"),
        ],
        total_power,
    )
    .expect("map admits");
    let (mapped_source, audit) = map
        .volumetric_source(&mesh, DEFAULT_POWER_TOLERANCE)
        .expect("projection admits");
    assert!((audit.delivered_total_w() - total_power).abs() < 1e-12);

    let uniform_source =
        ScalarField::uniform("volumetric source", total_power / volume).expect("uniform source");

    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "cold",
            |face| on_box_face(face.centroid[0], 0.0),
            ThermalBc::dirichlet(300.0).expect("cold"),
        )
        .expect("cold region")
        .adiabatic_remainder()
        .finish()
        .expect("complete partition");

    let material = ConductivityModel::isotropic_declared(10.0).expect("material");
    let solve_with = |source: &ScalarField| {
        with_cx(|cx| {
            solve(
                cx,
                ConductionProblem {
                    mesh: &mesh,
                    boundary: &boundary,
                    material: &material,
                    source,
                },
                config(),
            )
            .expect("solve")
        })
    };

    let mapped = solve_with(&mapped_source);
    let uniform = solve_with(&uniform_source);

    let mapped_t = &mapped.temperature;
    let uniform_t = &uniform.temperature;
    assert_eq!(mapped_t.len(), uniform_t.len());
    for (index, (a, b)) in mapped_t.iter().zip(uniform_t.iter()).enumerate() {
        assert!(
            (a - b).abs() < 1e-9,
            "vertex {index}: mapped {a} vs uniform {b}; the power path must \
             agree with the crate's existing uniform-source path"
        );
    }
}

// ---------------------------------------------------------------------------
// The audit: refusal carries the per-component table.
// ---------------------------------------------------------------------------

#[test]
fn a_declared_total_that_does_not_match_refuses_with_the_component_table() {
    let mesh = unit_mesh(2);
    let vertices = low_x_vertices(&mesh);
    let map = PowerMap::new(
        vec![
            ComponentPower::new("cpu", 30.0, PowerUncertainty::Unstated, vertices.clone())
                .expect("cpu"),
            ComponentPower::new("gpu", 45.0, PowerUncertainty::Unstated, vertices).expect("gpu"),
        ],
        100.0, // the forgotten-component case: 75 declared across rows
    )
    .expect("map admits");

    let error = map
        .volumetric_source(&mesh, DEFAULT_POWER_TOLERANCE)
        .expect_err("a map that does not add up must refuse");
    match error {
        ConductionError::ScenarioRow { what, fix, .. } => {
            assert!(what.contains("100"), "the declared total must be named");
            assert!(what.contains("75"), "the component sum must be named");
            // The whole point of the refusal: the reader can see WHICH rows.
            assert!(what.contains("cpu"), "the table must list cpu:\n{what}");
            assert!(what.contains("gpu"), "the table must list gpu:\n{what}");
            assert!(!fix.is_empty());
        }
        other => panic!("expected a scenario-row refusal, got {other:?}"),
    }
}

#[test]
fn a_matching_total_inside_tolerance_is_admitted() {
    let mesh = unit_mesh(2);
    let vertices = low_x_vertices(&mesh);
    // 0.1 + 0.2 is not exactly 0.3 in binary64; the tolerance exists for
    // exactly this, and a strict equality check here would be wrong.
    let map = PowerMap::new(
        vec![
            ComponentPower::new("a", 0.1, PowerUncertainty::Unstated, vertices.clone()).expect("a"),
            ComponentPower::new("b", 0.2, PowerUncertainty::Unstated, vertices).expect("b"),
        ],
        0.3,
    )
    .expect("map admits");
    assert!(
        map.volumetric_source(&mesh, DEFAULT_POWER_TOLERANCE)
            .is_ok(),
        "a binary64 rounding difference must not be reported as an entry error"
    );
}

#[test]
fn a_zero_tolerance_still_admits_an_exactly_matching_total() {
    let mesh = unit_mesh(2);
    let map = PowerMap::new(
        vec![
            ComponentPower::new("a", 4.0, PowerUncertainty::Unstated, low_x_vertices(&mesh))
                .expect("a"),
        ],
        4.0,
    )
    .expect("map admits");
    assert!(map.volumetric_source(&mesh, 0.0).is_ok());
}

// ---------------------------------------------------------------------------
// Uncertainty aggregation.
// ---------------------------------------------------------------------------

#[test]
fn one_unstated_component_makes_the_total_band_unavailable() {
    // The same rule the series thermal-resistance budget uses: an unstated
    // term makes the total unavailable rather than smaller.
    let mesh = unit_mesh(2);
    let vertices = low_x_vertices(&mesh);
    let banded = |name: &str, watts: f64| {
        ComponentPower::new(
            name,
            watts,
            PowerUncertainty::HalfWidth {
                half_width: 0.5,
                confidence: 0.95,
            },
            vertices.clone(),
        )
        .expect("banded component")
    };

    let all_banded = PowerMap::new(vec![banded("a", 1.0), banded("b", 2.0)], 3.0)
        .expect("map admits")
        .volumetric_source(&mesh, DEFAULT_POWER_TOLERANCE)
        .expect("projection")
        .1;
    assert_eq!(
        all_banded.total_half_width_w(),
        Some(1.0),
        "stated half-widths sum conservatively"
    );

    let one_unstated = PowerMap::new(
        vec![
            banded("a", 1.0),
            ComponentPower::new("b", 2.0, PowerUncertainty::Unstated, vertices).expect("b"),
        ],
        3.0,
    )
    .expect("map admits")
    .volumetric_source(&mesh, DEFAULT_POWER_TOLERANCE)
    .expect("projection")
    .1;
    assert_eq!(
        one_unstated.total_half_width_w(),
        None,
        "an unstated component must not silently shrink the total band"
    );
}

// ---------------------------------------------------------------------------
// Admission refusals.
// ---------------------------------------------------------------------------

#[test]
fn component_admission_refuses_degenerate_declarations() {
    let vertices = vec![0usize, 1, 2];
    assert!(
        ComponentPower::new("  ", 1.0, PowerUncertainty::Unstated, vertices.clone()).is_err(),
        "blank name"
    );
    assert!(
        ComponentPower::new("a", -1.0, PowerUncertainty::Unstated, vertices.clone()).is_err(),
        "a heat sink is a boundary condition, not a negative source"
    );
    for bad in [f64::NAN, f64::INFINITY] {
        assert!(
            ComponentPower::new("a", bad, PowerUncertainty::Unstated, vertices.clone()).is_err(),
            "{bad} watts"
        );
    }
    assert!(
        ComponentPower::new("a", 1.0, PowerUncertainty::Unstated, Vec::new()).is_err(),
        "a component that heats nothing"
    );
    assert!(
        ComponentPower::new(
            "a",
            1.0,
            PowerUncertainty::HalfWidth {
                half_width: -0.1,
                confidence: 0.95
            },
            vertices.clone()
        )
        .is_err(),
        "negative half-width"
    );
    assert!(
        ComponentPower::new(
            "a",
            1.0,
            PowerUncertainty::HalfWidth {
                half_width: 0.1,
                confidence: 1.0
            },
            vertices
        )
        .is_err(),
        "confidence must be strictly inside (0, 1)"
    );
}

#[test]
fn map_admission_refuses_duplicates_and_bad_totals() {
    let vertices = vec![0usize, 1];
    let component = |name: &str| {
        ComponentPower::new(name, 1.0, PowerUncertainty::Unstated, vertices.clone())
            .expect("component")
    };
    assert!(PowerMap::new(Vec::new(), 0.0).is_err(), "empty map");
    assert!(
        PowerMap::new(vec![component("a"), component("a")], 2.0).is_err(),
        "duplicate component name"
    );
    assert!(
        PowerMap::new(vec![component("a")], -1.0).is_err(),
        "negative declared total"
    );
    assert!(
        PowerMap::new(vec![component("a")], f64::NAN).is_err(),
        "non-finite declared total"
    );
}

#[test]
fn a_vertex_outside_the_mesh_refuses_by_name() {
    let mesh = unit_mesh(2);
    let outside = mesh.vertex_count() + 5;
    let map = PowerMap::new(
        vec![
            ComponentPower::new("stray", 1.0, PowerUncertainty::Unstated, vec![outside])
                .expect("component"),
        ],
        1.0,
    )
    .expect("map admits");
    let error = map
        .volumetric_source(&mesh, DEFAULT_POWER_TOLERANCE)
        .expect_err("an out-of-range binding must refuse");
    match error {
        ConductionError::ScenarioRow { region, what, .. } => {
            assert_eq!(region, "stray", "the refusal must name the component");
            assert!(what.contains(&outside.to_string()));
        }
        other => panic!("expected a scenario-row refusal, got {other:?}"),
    }
}

#[test]
fn a_negative_tolerance_refuses() {
    let mesh = unit_mesh(2);
    let map = PowerMap::new(
        vec![
            ComponentPower::new("a", 1.0, PowerUncertainty::Unstated, low_x_vertices(&mesh))
                .expect("a"),
        ],
        1.0,
    )
    .expect("map admits");
    assert!(map.volumetric_source(&mesh, -1e-9).is_err());
}

// ---------------------------------------------------------------------------
// Determinism and reporting.
// ---------------------------------------------------------------------------

#[test]
fn projection_is_deterministic_and_component_order_independent() {
    let mesh = unit_mesh(3);
    let vertices = low_x_vertices(&mesh);
    let build = |reversed: bool| {
        let mut components = vec![
            ComponentPower::new("alpha", 5.0, PowerUncertainty::Unstated, vertices.clone())
                .expect("alpha"),
            ComponentPower::new("beta", 6.0, PowerUncertainty::Unstated, vertices.clone())
                .expect("beta"),
        ];
        if reversed {
            components.reverse();
        }
        PowerMap::new(components, 11.0)
            .expect("map admits")
            .volumetric_source(&mesh, DEFAULT_POWER_TOLERANCE)
            .expect("projection")
    };
    let (forward_field, forward_audit) = build(false);
    let (reversed_field, reversed_audit) = build(true);

    assert_eq!(
        forward_field, reversed_field,
        "declaration order must not change the field"
    );
    assert_eq!(forward_audit, reversed_audit);
    // And the map sorts by name, so the audit rows are in a stable order.
    let names: Vec<&str> = forward_audit.rows().iter().map(|row| row.name()).collect();
    assert_eq!(names, vec!["alpha", "beta"]);
}

#[test]
fn the_audit_table_reports_every_row_and_all_three_totals() {
    let mesh = unit_mesh(2);
    let vertices = low_x_vertices(&mesh);
    let map = PowerMap::new(
        vec![
            ComponentPower::new("cpu", 2.0, PowerUncertainty::Unstated, vertices.clone())
                .expect("cpu"),
            ComponentPower::new("dram", 3.0, PowerUncertainty::Unstated, vertices).expect("dram"),
        ],
        5.0,
    )
    .expect("map admits");
    let (_, audit) = map
        .volumetric_source(&mesh, DEFAULT_POWER_TOLERANCE)
        .expect("projection");

    let table = audit.render_table();
    assert!(table.contains("cpu"));
    assert!(table.contains("dram"));
    assert!(table.contains("declared"));
    assert!(table.contains("delivered"));
    assert_eq!(audit.rows().len(), 2);
    for row in audit.rows() {
        assert!(row.bound_volume_m3() > 0.0);
        assert!(row.density_w_per_m3() > 0.0);
        assert!((row.delivered_w() - row.declared_w()).abs() < 1e-12);
    }
}
