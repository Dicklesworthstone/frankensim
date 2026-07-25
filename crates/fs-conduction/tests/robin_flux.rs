//! Per-named-region Robin heat-rate accounting.
//!
//! `EnergyBalance::robin_out_w` is a whole-domain scalar. A conjugate
//! driver has to hand a SPECIFIC surface's heat to a SPECIFIC air path,
//! so the report also carries the same integral restricted to each
//! declared trace. These tests check that decomposition against
//! closed-form 1-D solutions rather than against the total it came from:
//! a test that only summed the parts back up would pass even if two
//! regions had their heat rates swapped.
//!
//! Gauntlet tier: G1 (closed-form oracle) plus G0 (partition invariants).

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::{
    ConductionMesh, ConductionProblem, ConductionSolution, ConductivityModel, InitialGuess,
    LinearConfig, Nonlinearity, ScalarField, SolveConfig, StopRule, ThermalBc,
    ThermalBoundaryBuilder, solve,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x0000_C0DE_B012_0000,
                kernel_id: 57,
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
        // A purely-Robin slab has no Dirichlet row pinning it, so the
        // operator is worse conditioned than the pinned fixtures
        // elsewhere in this crate. Measured achievable relative residuals
        // here are ~1e-12 (two ends) and ~1.3e-12 (one end); the solver
        // fails closed rather than silently returning a looser answer, so
        // this is set above both. It is still four orders tighter than
        // the 1e-9 the closed-form comparisons need.
        linear: LinearConfig {
            tolerance: 1.0e-11,
            max_iterations: 20_000,
            restart: 40,
        },
        initial: InitialGuess::Uniform(300.0),
    }
}

// ---------------------------------------------------------------------
// The doubly-convected generating slab.
//
//   -k T'' = f  on 0 < x < L,  Robin at BOTH ends, sides adiabatic.
//
//   T(x)  = -f x²/(2k) + c₁ x + c₀
//   q(x)  = -k T'(x) = f x - k c₁
//
// Outward flux at x=L is +q(L); at x=0 it is -q(0) = k c₁. Imposing both
// Robin rows and eliminating c₀ gives the closed form below. The two end
// heat rates are deliberately ~20x apart, so a swapped attribution is
// caught by magnitude and not merely by sign.
// ---------------------------------------------------------------------

const LENGTH: f64 = 0.020;
const WIDTH: f64 = 0.010;
const HEIGHT: f64 = 0.010;
const K_SOLID: f64 = 20.0;
const SOURCE_W_M3: f64 = 5.0e5;
const H_WEAK: f64 = 15.0;
const T_INF_WEAK: f64 = 300.0;
const H_STRONG: f64 = 250.0;
const T_INF_STRONG: f64 = 290.0;

struct Analytic {
    heat_out_weak: f64,
    heat_out_strong: f64,
    wall_weak: f64,
    wall_strong: f64,
    area: f64,
}

fn analytic() -> Analytic {
    let area = WIDTH * HEIGHT;
    let numerator = SOURCE_W_M3 * LENGTH - H_STRONG * (T_INF_WEAK - T_INF_STRONG)
        + H_STRONG * SOURCE_W_M3 * LENGTH * LENGTH / (2.0 * K_SOLID);
    let denominator = K_SOLID + H_STRONG * (LENGTH + K_SOLID / H_WEAK);
    let c1 = numerator / denominator;
    let c0 = T_INF_WEAK + K_SOLID * c1 / H_WEAK;
    Analytic {
        heat_out_weak: area * K_SOLID * c1,
        heat_out_strong: area * (SOURCE_W_M3 * LENGTH - K_SOLID * c1),
        wall_weak: c0,
        wall_strong: -SOURCE_W_M3 * LENGTH * LENGTH / (2.0 * K_SOLID) + c1 * LENGTH + c0,
        area,
    }
}

fn generating_slab(cells_x: usize) -> (ConductionMesh, ConductionSolution) {
    let (complex, positions) = box_grid([cells_x, 2, 2], [LENGTH, WIDTH, HEIGHT]);
    let mesh = ConductionMesh::new(complex, positions).expect("generating slab mesh");
    let material = ConductivityModel::isotropic_declared(K_SOLID).expect("solid material");
    let source = ScalarField::Uniform(SOURCE_W_M3);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "weakly-cooled",
            |face| on_box_face(face.centroid[0], 0.0),
            ThermalBc::robin(H_WEAK, T_INF_WEAK).expect("weakly-cooled Robin"),
        )
        .expect("weakly-cooled region")
        .region(
            "strongly-cooled",
            |face| on_box_face(face.centroid[0], LENGTH),
            ThermalBc::robin(H_STRONG, T_INF_STRONG).expect("strongly-cooled Robin"),
        )
        .expect("strongly-cooled region")
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
        .expect("generating slab solve")
    });
    (mesh, solution)
}

#[test]
fn each_robin_region_matches_its_own_closed_form_heat_rate() {
    let expected = analytic();
    let (_mesh, solution) = generating_slab(16);
    let rows = &solution.report.robin_fluxes;

    assert_eq!(rows.len(), 2, "two Robin regions were declared");
    assert_eq!(rows[0].region, "weakly-cooled");
    assert_eq!(rows[1].region, "strongly-cooled");

    // The fixture must be able to tell the two apart. If the ends carried
    // comparable heat, a swapped attribution would survive every
    // assertion below and this file would prove nothing about which
    // surface owns which watt.
    let separation = (expected.heat_out_strong / expected.heat_out_weak).abs();
    assert!(
        separation > 10.0,
        "asymmetric fixture required: ends differ by only {separation}x"
    );

    for (row, expected_heat, expected_wall, expected_h, expected_ref) in [
        (
            &rows[0],
            expected.heat_out_weak,
            expected.wall_weak,
            H_WEAK,
            T_INF_WEAK,
        ),
        (
            &rows[1],
            expected.heat_out_strong,
            expected.wall_strong,
            H_STRONG,
            T_INF_STRONG,
        ),
    ] {
        let heat_error = (row.heat_rate_w - expected_heat).abs() / expected_heat.abs();
        assert!(
            heat_error < 1.0e-9,
            "{}: heat {} vs closed form {expected_heat}, relative {heat_error}",
            row.region,
            row.heat_rate_w
        );
        let wall_error = (row.mean_wall_temperature_k - expected_wall).abs() / expected_wall;
        assert!(
            wall_error < 1.0e-9,
            "{}: wall {} vs closed form {expected_wall}, relative {wall_error}",
            row.region,
            row.mean_wall_temperature_k
        );
        // The means are AREA-WEIGHTED reductions, `Σ(A_f·x)/ΣA_f`. Even a
        // uniform field does not round-trip bitwise through that — it
        // lands within an ULP — so exact-bit equality is the wrong
        // assertion to make here and this pins the reduction instead.
        let htc_error = (row.mean_htc_w_per_m2_k - expected_h).abs() / expected_h;
        assert!(
            htc_error < 1.0e-15,
            "{}: mean h {} vs uniform {expected_h}, relative {htc_error}",
            row.region,
            row.mean_htc_w_per_m2_k
        );
        let ref_error = (row.mean_reference_temperature_k - expected_ref).abs() / expected_ref;
        assert!(
            ref_error < 1.0e-15,
            "{}: mean T_ref {} vs uniform {expected_ref}, relative {ref_error}",
            row.region,
            row.mean_reference_temperature_k
        );
        let area_error = (row.area_m2 - expected.area).abs() / expected.area;
        assert!(area_error < 1.0e-12, "{}: area {}", row.region, row.area_m2);
        assert!(row.faces > 0);
    }
}

#[test]
fn the_two_end_heat_rates_account_for_the_whole_generated_power() {
    // An independent closure: the volumetric source is the only supply
    // and the sides are adiabatic, so the two convective traces must
    // carry all of it. This crosses the region decomposition against the
    // volume integral, which is assembled by a different loop.
    let (_mesh, solution) = generating_slab(16);
    let generated = solution.report.energy.source_w;
    let removed: f64 = solution
        .report
        .robin_fluxes
        .iter()
        .map(|row| row.heat_rate_w)
        .sum();
    let error = (removed - generated).abs() / generated;
    assert!(
        error < 1.0e-9,
        "traces removed {removed} W of {generated} W generated, relative {error}"
    );
}

#[test]
fn the_decomposition_sums_to_the_whole_domain_total_up_to_summation_order() {
    // `robin_out_w` keeps its historical running-accumulator value; the
    // per-region rows use separate accumulators. The two therefore agree
    // only up to floating-point summation order, and this pins that gap
    // rather than pretending it is zero.
    let (_mesh, solution) = generating_slab(16);
    let total = solution.report.energy.robin_out_w;
    let summed: f64 = solution
        .report
        .robin_fluxes
        .iter()
        .map(|row| row.heat_rate_w)
        .sum();
    let error = (summed - total).abs() / total.abs();
    assert!(
        error < 1.0e-12,
        "sum of parts {summed} vs whole-domain {total}, relative {error}"
    );
}

#[test]
fn only_robin_regions_emit_rows_and_they_follow_declaration_order() {
    const T_BASE: f64 = 350.0;
    const T_AMBIENT: f64 = 300.0;
    const H: f64 = 40.0;

    let (complex, positions) = box_grid([4, 2, 2], [LENGTH, WIDTH, HEIGHT]);
    let mesh = ConductionMesh::new(complex, positions).expect("mixed-row mesh");
    let material = ConductivityModel::isotropic_declared(K_SOLID).expect("solid material");
    let source = ScalarField::Uniform(0.0);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "base",
            |face| on_box_face(face.centroid[0], 0.0),
            ThermalBc::dirichlet(T_BASE).expect("Dirichlet base"),
        )
        .expect("base region")
        .region(
            "vented",
            |face| on_box_face(face.centroid[1], 0.0),
            ThermalBc::neumann(0.0).expect("Neumann vent"),
        )
        .expect("vent region")
        .region(
            "convected",
            |face| on_box_face(face.centroid[0], LENGTH),
            ThermalBc::robin(H, T_AMBIENT).expect("Robin"),
        )
        .expect("convected region")
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
        .expect("mixed-row solve")
    });

    let rows = &solution.report.robin_fluxes;
    assert_eq!(
        rows.len(),
        1,
        "Dirichlet, Neumann and the adiabatic remainder own no convective heat"
    );
    assert_eq!(rows[0].region, "convected");

    // Series resistance through the slab, as an independent value.
    let area = WIDTH * HEIGHT;
    let analytic_heat = area * K_SOLID * H * (T_BASE - T_AMBIENT) / (K_SOLID + H * LENGTH);
    let error = (rows[0].heat_rate_w - analytic_heat).abs() / analytic_heat;
    assert!(
        error < 1.0e-9,
        "convected {} vs series-resistance {analytic_heat}, relative {error}",
        rows[0].heat_rate_w
    );
}

#[test]
fn a_robin_region_owning_no_face_emits_no_row() {
    // Declaring a region whose predicate matches nothing is legal. It
    // contributes nothing to the balance and has no area to average
    // over, so it emits no row rather than a zero-area row whose means
    // would be invented.
    let (complex, positions) = box_grid([4, 2, 2], [LENGTH, WIDTH, HEIGHT]);
    let mesh = ConductionMesh::new(complex, positions).expect("empty-region mesh");
    let material = ConductivityModel::isotropic_declared(K_SOLID).expect("solid material");
    let source = ScalarField::Uniform(SOURCE_W_M3);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "nowhere",
            |face| on_box_face(face.centroid[0], 17.0),
            ThermalBc::robin(H_STRONG, T_INF_STRONG).expect("Robin"),
        )
        .expect("empty region")
        .region(
            "weakly-cooled",
            |face| on_box_face(face.centroid[0], 0.0),
            ThermalBc::robin(H_WEAK, T_INF_WEAK).expect("Robin"),
        )
        .expect("weakly-cooled region")
        .adiabatic_remainder()
        .finish()
        .expect("boundary partition");
    assert_eq!(boundary.region_names(), ["nowhere", "weakly-cooled"]);

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
        .expect("empty-region solve")
    });

    let rows = &solution.report.robin_fluxes;
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].region, "weakly-cooled");
    assert!(rows[0].heat_rate_w.is_finite());
    assert!(!rows[0].mean_wall_temperature_k.is_nan());
}

#[test]
fn the_rows_cover_exactly_the_faces_whose_condition_is_robin() {
    // The decomposition is built on `region_for`. Counting Robin faces
    // independently through `condition_for` checks that no face was
    // dropped from a trace or double-counted into one — the wiring
    // failure a heat-rate comparison alone would hide, because a lost
    // face changes the total and the part together.
    let (mesh, solution) = generating_slab(8);
    let robin_faces = (0..mesh.boundary().len())
        .filter(|&slot| {
            matches!(
                solution_boundary_condition(&mesh, slot),
                Some(BoundaryKind::Robin)
            )
        })
        .count();
    let attributed: usize = solution
        .report
        .robin_fluxes
        .iter()
        .map(|row| row.faces)
        .sum();
    assert_eq!(attributed, robin_faces);
    // Two opposing end faces, each 2x2 cells split into two triangles.
    assert_eq!(robin_faces, 16);
}

enum BoundaryKind {
    Robin,
}

/// Re-derive each face's condition kind from the partition the fixture
/// declared, independently of the solve report.
fn solution_boundary_condition(mesh: &ConductionMesh, slot: usize) -> Option<BoundaryKind> {
    let centroid_x = mesh.boundary()[slot].centroid[0];
    (on_box_face(centroid_x, 0.0) || on_box_face(centroid_x, LENGTH)).then_some(BoundaryKind::Robin)
}
