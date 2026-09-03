//! `frankensim-m4rb` vector-compliance DWR acceptance fixtures.

use fs_cutfem::{
    BoundaryTraction, Circle, CutElasticity, CutFemError, CutSdf, CutStabilizationScaling,
    DesignBoxEdge, EdgeBand, HalfPlane, Quadtree,
};
use fs_dwr::{
    ElasticityDwrEstimate, dorfler, estimate_elasticity_compliance,
    estimate_elasticity_compliance_with_boundary_traction,
};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;

fn material() -> IsotropicElastic {
    IsotropicElastic::new(1.0, 0.3, 10.0).expect("fixture material")
}

fn problem<'a>(
    grid: &'a Quadtree,
    sdf: &'a dyn CutSdf,
    material: &'a IsotropicElastic,
    clamp: Option<&'a dyn Fn(f64, f64) -> bool>,
    traction: Option<&'a dyn Fn(f64, f64) -> [f64; 2]>,
    traction_free_interface: bool,
    ghost_gamma: f64,
) -> CutElasticity<'a> {
    CutElasticity {
        grid,
        sdf,
        material,
        nitsche_beta: 100.0,
        ghost_gamma,
        stabilization_scaling: CutStabilizationScaling::MuScaled,
        quad_depth: 3,
        clamp,
        boundary_traction: traction,
        traction_free_interface,
        solver_tol: 1e-12,
        solver_max_iters: 60_000,
    }
}

#[derive(Debug, Clone, Copy)]
struct AlwaysInside;

impl CutSdf for AlwaysInside {
    fn value(&self, _point: [f64; 2]) -> f64 {
        -1.0
    }

    fn gradient(&self, _point: [f64; 2]) -> [f64; 2] {
        [1.0, 0.0]
    }

    fn enclose(&self, _lo: [f64; 2], _hi: [f64; 2]) -> Interval {
        Interval::point(-1.0)
    }
}

/// Material outside an interior circular hole is negative-inside.
#[derive(Debug, Clone, Copy)]
struct CircularHole {
    center: [f64; 2],
    radius: f64,
}

impl CutSdf for CircularHole {
    fn value(&self, point: [f64; 2]) -> f64 {
        let dx = point[0] - self.center[0];
        let dy = point[1] - self.center[1];
        self.radius - dx.hypot(dy)
    }

    fn gradient(&self, point: [f64; 2]) -> [f64; 2] {
        let dx = point[0] - self.center[0];
        let dy = point[1] - self.center[1];
        let norm = dx.hypot(dy);
        if norm < 1e-300 {
            [1.0, 0.0]
        } else {
            [-dx / norm, -dy / norm]
        }
    }

    fn enclose(&self, lo: [f64; 2], hi: [f64; 2]) -> Interval {
        let dx = (Interval::new(lo[0], hi[0]) - Interval::point(self.center[0])).abs();
        let dy = (Interval::new(lo[1], hi[1]) - Interval::point(self.center[1])).abs();
        Interval::point(self.radius) - (dx * dx + dy * dy).sqrt()
    }
}

fn assert_reconstruction(estimate: &ElasticityDwrEstimate) {
    let cell_sum: f64 = estimate.indicators.values().sum();
    let cell_abs: f64 = estimate.indicators.values().map(|value| value.abs()).sum();
    let face_sum: f64 = estimate.face_indicators.values().sum();
    let scale = estimate.eta_signed.abs().max(1.0);
    assert!(
        (cell_sum - estimate.eta_signed).abs() <= 8.0 * f64::EPSILON * scale,
        "cell allocation reconstructs eta: cells={cell_sum:.17e}, eta={:.17e}",
        estimate.eta_signed
    );
    assert!(
        (estimate.terms.total() - estimate.eta_signed).abs() <= 64.0 * f64::EPSILON * scale,
        "term decomposition reconstructs eta: terms={:.17e}, eta={:.17e}",
        estimate.terms.total(),
        estimate.eta_signed
    );
    assert_eq!(
        cell_abs.to_bits(),
        estimate.eta_abs.to_bits(),
        "eta_abs is exactly the deterministic cell marking mass"
    );
    assert!(
        (face_sum - estimate.terms.ghost).abs() <= 16.0 * f64::EPSILON * face_sum.abs().max(1.0),
        "face indicators reconstruct the coarse consistent-energy correction"
    );
    assert_eq!(estimate.ghost_method.as_str(), "coarse-consistent-energy");
    let cell_term_sums = estimate.cell_terms.values().fold(
        fs_dwr::ElasticityResidualTerms::default(),
        |mut sums, terms| {
            sums.bulk += terms.bulk;
            sums.nitsche += terms.nitsche;
            sums.outer_traction += terms.outer_traction;
            sums.ghost += terms.ghost;
            sums
        },
    );
    for (name, cell_sum, global) in [
        ("bulk", cell_term_sums.bulk, estimate.terms.bulk),
        ("Nitsche", cell_term_sums.nitsche, estimate.terms.nitsche),
        (
            "outer traction",
            cell_term_sums.outer_traction,
            estimate.terms.outer_traction,
        ),
        ("ghost", cell_term_sums.ghost, estimate.terms.ghost),
    ] {
        let term_scale = cell_sum.abs().max(global.abs()).max(1.0);
        assert!(
            (cell_sum - global).abs() <= 64.0 * f64::EPSILON * term_scale,
            "per-cell {name} decomposition must reconstruct its global term: cells={cell_sum:.17e}, global={global:.17e}"
        );
    }
    for (cell, terms) in &estimate.cell_terms {
        let indicator = estimate.indicators[cell];
        let term_scale = indicator.abs().max(terms.total().abs()).max(1.0);
        assert!(
            (terms.total() - indicator).abs() <= 16.0 * f64::EPSILON * term_scale,
            "per-cell terms must reconstruct aggregate indicator at {cell:?}"
        );
    }
}

fn assert_coarse_ghost_energy_matches_operator(
    problem: &CutElasticity<'_>,
    body: &dyn Fn(f64, f64) -> [f64; 2],
    embedded_data: &dyn Fn(f64, f64) -> [f64; 2],
    estimate: &ElasticityDwrEstimate,
) {
    let solution = problem
        .solve(body, embedded_data)
        .expect("coarse ghost-energy witness solve");
    let stabilized = problem
        .assemble(body, embedded_data)
        .expect("stabilized coarse operator");
    let unstabilized_problem = CutElasticity {
        ghost_gamma: 0.0,
        ..*problem
    };
    let unstabilized = unstabilized_problem
        .assemble(body, embedded_data)
        .expect("unstabilized coarse operator");
    let stabilized_action = stabilized.apply_vec(solution.coefficients());
    let unstabilized_action = unstabilized.apply_vec(solution.coefficients());
    let operator_energy: f64 = solution
        .coefficients()
        .iter()
        .zip(stabilized_action.iter().zip(&unstabilized_action))
        .map(|(coefficient, (with_ghost, without_ghost))| {
            coefficient * (with_ghost - without_ghost)
        })
        .sum();
    let scale = operator_energy
        .abs()
        .max(estimate.terms.ghost.abs())
        .max(1.0);
    assert!(
        (operator_energy - estimate.terms.ghost).abs() <= 256.0 * f64::EPSILON * scale,
        "coarse face decomposition must match x^T(K_gamma-K_0)x: faces={:.17e}, operator={operator_energy:.17e}",
        estimate.terms.ghost
    );
}

fn assert_bitwise_equal(a: &ElasticityDwrEstimate, b: &ElasticityDwrEstimate) {
    for (left, right) in [
        (a.eta_signed, b.eta_signed),
        (a.eta_abs, b.eta_abs),
        (a.terms.bulk, b.terms.bulk),
        (a.terms.nitsche, b.terms.nitsche),
        (a.terms.outer_traction, b.terms.outer_traction),
        (a.terms.ghost, b.terms.ghost),
        (a.j_primal, b.j_primal),
        (a.j_enriched, b.j_enriched),
    ] {
        assert_eq!(left.to_bits(), right.to_bits(), "bitwise replay");
    }
    assert_eq!(a.dofs, b.dofs);
    assert_eq!(a.enriched_dofs, b.enriched_dofs);
    assert_eq!(a.ghost_method, b.ghost_method);
    assert_map_bits_equal(&a.indicators, &b.indicators);
    assert_cell_term_bits_equal(&a.cell_terms, &b.cell_terms);
    assert_map_bits_equal(&a.face_indicators, &b.face_indicators);
}

fn assert_cell_term_bits_equal(
    a: &BTreeMap<fs_cutfem::CellKey, fs_dwr::ElasticityResidualTerms>,
    b: &BTreeMap<fs_cutfem::CellKey, fs_dwr::ElasticityResidualTerms>,
) {
    assert_eq!(a.len(), b.len());
    for ((key_a, terms_a), (key_b, terms_b)) in a.iter().zip(b) {
        assert_eq!(key_a, key_b);
        for (left, right) in [
            (terms_a.bulk, terms_b.bulk),
            (terms_a.nitsche, terms_b.nitsche),
            (terms_a.outer_traction, terms_b.outer_traction),
            (terms_a.ghost, terms_b.ghost),
        ] {
            assert_eq!(left.to_bits(), right.to_bits(), "cell term at {key_a:?}");
        }
    }
}

fn assert_map_bits_equal<K: Ord + core::fmt::Debug>(a: &BTreeMap<K, f64>, b: &BTreeMap<K, f64>) {
    assert_eq!(a.len(), b.len());
    for ((key_a, value_a), (key_b, value_b)) in a.iter().zip(b) {
        assert_eq!(key_a, key_b);
        assert_eq!(
            value_a.to_bits(),
            value_b.to_bits(),
            "map value at {key_a:?}"
        );
    }
}

fn assert_effectivity(case: &str, level: u32, estimate: &ElasticityDwrEstimate, reference: f64) {
    let true_error = reference - estimate.j_primal;
    let scale = reference.abs().max(estimate.j_primal.abs()).max(1.0);
    assert!(
        true_error.is_finite() && true_error.abs() > 1e-10 * scale,
        "{case} level {level} reference error must be finite and non-degenerate: {true_error:.6e}"
    );
    let effectivity = estimate.eta_signed / true_error;
    let hierarchy_increment = estimate.j_enriched - estimate.j_primal;
    assert!(
        hierarchy_increment.is_finite() && hierarchy_increment.abs() > 1e-10 * scale,
        "{case} level {level} hierarchy increment must be finite and non-degenerate: {hierarchy_increment:.6e}"
    );
    let residual_to_increment = estimate.eta_signed / hierarchy_increment;
    println!(
        "{{\"suite\":\"fs-dwr/elasticity\",\"case\":\"{case}\",\"level\":{level},\"j_h\":{:.10e},\"j_enriched\":{:.10e},\"j_ref\":{reference:.10e},\"true_error\":{true_error:.10e},\"hierarchy_increment\":{hierarchy_increment:.10e},\"eta_signed\":{:.10e},\"eta_abs\":{:.10e},\"bulk\":{:.10e},\"nitsche\":{:.10e},\"outer_traction\":{:.10e},\"ghost\":{:.10e},\"effectivity\":{effectivity:.8},\"residual_to_increment\":{residual_to_increment:.8},\"ghost_method\":\"{}\"}}",
        estimate.j_primal,
        estimate.j_enriched,
        estimate.eta_signed,
        estimate.eta_abs,
        estimate.terms.bulk,
        estimate.terms.nitsche,
        estimate.terms.outer_traction,
        estimate.terms.ghost,
        estimate.ghost_method.as_str(),
    );
    assert!(
        (0.5..=2.0).contains(&effectivity),
        "{case} level {level} effectivity {effectivity:.6} outside [0.5, 2.0]"
    );
}

#[test]
fn elasticity_estimator_refuses_lattice_cap_before_solving() {
    let grid = Quadtree::with_room(0, 16);
    let sdf = AlwaysInside;
    let material = material();
    let zero = |_: f64, _: f64| [0.0, 0.0];
    let err = estimate_elasticity_compliance(
        &problem(&grid, &sdf, &material, None, None, true, 0.0),
        &zero,
        &zero,
    )
    .expect_err("level-16 grid cannot be enriched safely");
    assert!(matches!(err, CutFemError::InvalidElasticityInput { .. }));
}

#[test]
fn graded_elasticity_dwr_reuses_mixed_level_ghost_patches() {
    let mut grid = Quadtree::with_room(1, 2);
    grid.split((1, 1, 0));
    let sdf = HalfPlane {
        normal: [1.0, 0.0],
        offset: 0.6,
    };
    let material = material();
    let zero = |_: f64, _: f64| [0.0, 0.0];
    let body = |_: f64, _: f64| [0.1, -0.07];
    let cut = problem(&grid, &sdf, &material, None, None, false, 0.5);

    let coarse = cut.solve(&body, &zero).expect("graded coarse solve");
    assert!(
        coarse.ghost_faces().iter().any(|(a, b)| a.0 != b.0),
        "fixture must retain mixed-level coarse ghost evidence"
    );
    let estimate = estimate_elasticity_compliance(&cut, &body, &zero).expect("graded DWR estimate");
    assert_reconstruction(&estimate);
    assert_eq!(
        estimate.face_indicators.keys().copied().collect::<Vec<_>>(),
        coarse.ghost_faces().to_vec(),
        "DWR must retain exactly the operator's canonical ghost-patch keys"
    );
    let mixed_ghost_energy: f64 = estimate
        .face_indicators
        .iter()
        .filter_map(|(&(a, b), &value)| (a.0 != b.0).then_some(value))
        .sum();
    assert!(
        mixed_ghost_energy.is_finite() && mixed_ghost_energy > 0.0,
        "mixed-level shared patches must carry nondegenerate positive DWR ghost energy: {mixed_ghost_energy:e}"
    );
    assert_coarse_ghost_energy_matches_operator(&cut, &body, &zero, &estimate);

    let replay = estimate_elasticity_compliance(&cut, &body, &zero).expect("graded DWR replay");
    assert_bitwise_equal(&estimate, &replay);
}

#[test]
#[allow(
    clippy::float_cmp,
    clippy::too_many_lines,
    reason = "exact dyadic box-boundary fixtures and one end-to-end acceptance narrative"
)]
fn elasticity_compliance_effectivity_family_and_determinism() {
    let material = material();
    let (lambda, mu) = material.lame();
    let zero = |_: f64, _: f64| [0.0, 0.0];
    let clamp_left = |x: f64, _: f64| x == 0.0;

    // Fixture 1: full box, exact quadratic displacement, body force, and
    // exact outer traction. It isolates the bulk and dead-traction signs.
    let full = AlwaysInside;
    let displacement_scale = 0.01;
    let full_body = |_: f64, _: f64| [-2.0 * displacement_scale * (lambda + 2.0 * mu), 0.0];
    let full_traction = |x: f64, y: f64| {
        if x == 1.0 {
            [2.0 * displacement_scale * (lambda + 2.0 * mu), 0.0]
        } else if y == 0.0 {
            [0.0, -2.0 * displacement_scale * lambda * x]
        } else if y == 1.0 {
            [0.0, 2.0 * displacement_scale * lambda * x]
        } else {
            [0.0, 0.0]
        }
    };
    let full_reference = {
        let grid = Quadtree::uniform(6);
        problem(
            &grid,
            &full,
            &material,
            Some(&clamp_left),
            Some(&full_traction),
            true,
            0.0,
        )
        .solve(&full_body, &zero)
        .expect("full-box reference")
        .compliance()
    };
    for level in [3u32, 4] {
        let grid = Quadtree::uniform(level);
        let estimate = estimate_elasticity_compliance(
            &problem(
                &grid,
                &full,
                &material,
                Some(&clamp_left),
                Some(&full_traction),
                true,
                0.0,
            ),
            &full_body,
            &zero,
        )
        .expect("full-box estimate");
        assert_reconstruction(&estimate);
        assert!(estimate.terms.bulk.abs() > 0.0);
        assert!(estimate.terms.outer_traction.abs() > 0.0);
        assert_eq!(estimate.terms.nitsche.to_bits(), 0.0f64.to_bits());
        assert_eq!(estimate.terms.ghost.to_bits(), 0.0f64.to_bits());
        assert!(estimate.face_indicators.is_empty());
        assert_effectivity("full-box", level, &estimate, full_reference);
    }

    // Fixture 2: an all-embedded disk with a quadratic exact displacement
    // zero on the circle. It exercises Nitsche and coarse ghost terms.
    let disk = Circle {
        center: [0.47, 0.53],
        radius: 0.31,
    };
    let disk_body = |_: f64, _: f64| [-2.0 * displacement_scale * (lambda + 3.0 * mu), 0.0];
    let disk_reference = {
        let grid = Quadtree::uniform(6);
        problem(&grid, &disk, &material, None, None, false, 0.5)
            .solve(&disk_body, &zero)
            .expect("disk reference")
            .compliance()
    };
    let mut disk_level_three = None;
    for level in [3u32, 4] {
        let grid = Quadtree::uniform(level);
        let cut = problem(&grid, &disk, &material, None, None, false, 0.5);
        let estimate =
            estimate_elasticity_compliance(&cut, &disk_body, &zero).expect("disk estimate");
        assert_reconstruction(&estimate);
        assert!(estimate.terms.bulk.abs() > 0.0);
        assert!(estimate.terms.nitsche.abs() > 0.0);
        assert!(estimate.terms.ghost.abs() > 0.0);
        assert_eq!(estimate.terms.outer_traction.to_bits(), 0.0f64.to_bits());
        assert!(!estimate.face_indicators.is_empty());
        if level == 3 {
            assert_coarse_ghost_energy_matches_operator(&cut, &disk_body, &zero, &estimate);
        }
        assert_effectivity("embedded-disk", level, &estimate, disk_reference);
        if level == 3 {
            disk_level_three = Some((grid, estimate));
        }
    }
    let (replay_grid, first) = disk_level_three.expect("level-three disk estimate");
    let replay = estimate_elasticity_compliance(
        &problem(&replay_grid, &disk, &material, None, None, false, 0.5),
        &disk_body,
        &zero,
    )
    .expect("deterministic replay");
    assert_bitwise_equal(&first, &replay);

    // Disabling ghost stabilization is exact in both term and evidence shape.
    let zero_ghost_grid = Quadtree::uniform(3);
    let zero_ghost = estimate_elasticity_compliance(
        &problem(&zero_ghost_grid, &disk, &material, None, None, false, 0.0),
        &disk_body,
        &zero,
    )
    .expect("ghost-free disk estimate");
    assert_eq!(zero_ghost.terms.ghost.to_bits(), 0.0f64.to_bits());
    assert!(zero_ghost.face_indicators.is_empty());
    assert_reconstruction(&zero_ghost);

    // Fixture 3: off-grid plate with a traction-free circular hole, left
    // clamp, and a smooth right-edge dead traction. It carries the literal
    // stress/cut-boundary localization claim.
    let hole = CircularHole {
        center: [0.47, 0.53],
        radius: 0.24,
    };
    let hole_body = |_: f64, _: f64| [0.0, 0.0];
    let right_traction = |x: f64, y: f64| {
        if x == 1.0 {
            [0.05 * (PI * y).sin().powi(2), 0.0]
        } else {
            [0.0, 0.0]
        }
    };
    let hole_reference = {
        let grid = Quadtree::uniform(6);
        problem(
            &grid,
            &hole,
            &material,
            Some(&clamp_left),
            Some(&right_traction),
            true,
            0.5,
        )
        .solve(&hole_body, &zero)
        .expect("plate-with-hole reference")
        .compliance()
    };
    let mut localization_estimate = None;
    for level in [3u32, 4] {
        let grid = Quadtree::uniform(level);
        let estimate = estimate_elasticity_compliance(
            &problem(
                &grid,
                &hole,
                &material,
                Some(&clamp_left),
                Some(&right_traction),
                true,
                0.5,
            ),
            &hole_body,
            &zero,
        )
        .expect("plate-with-hole estimate");
        assert_reconstruction(&estimate);
        assert!(estimate.terms.bulk.abs() > 0.0);
        assert_eq!(estimate.terms.nitsche.to_bits(), 0.0f64.to_bits());
        assert!(estimate.terms.outer_traction.abs() > 0.0);
        assert!(estimate.terms.ghost.abs() > 0.0);
        assert_effectivity("plate-with-hole", level, &estimate, hole_reference);
        if level == 4 {
            localization_estimate = Some((grid, estimate));
        }
    }
    let (localization_grid, localization) = localization_estimate.expect("level-four hole");
    assert_top_decile_near_hole(&localization_grid, &hole, &localization);
}

fn assert_top_decile_near_hole(
    grid: &Quadtree,
    hole: &CircularHole,
    estimate: &ElasticityDwrEstimate,
) {
    let mut ranked: Vec<(fs_cutfem::CellKey, f64)> = estimate
        .cell_terms
        .iter()
        .map(|(&cell, terms)| {
            // The dead-load edge residual legitimately localizes where the load
            // is applied.  The separate stress/cut-boundary claim therefore
            // uses the non-traction part of the exact per-cell decomposition;
            // the aggregate marking indicator still retains every term.
            let hole_sensitive = (terms.bulk + terms.nitsche) + terms.ghost;
            (cell, hole_sensitive.abs())
        })
        .collect();
    ranked.sort_by(|left, right| right.1.total_cmp(&left.1).then(left.0.cmp(&right.0)));
    let count = ranked.len().div_ceil(10).max(1);
    let top = &ranked[..count];
    let k = 3.0;
    let localization_mass: f64 = ranked.iter().map(|(_, value)| value).sum();
    assert!(
        localization_mass.is_finite() && localization_mass > 0.0,
        "hole-sensitive localization mass must be finite and nonzero"
    );
    let mut top_mass = 0.0;
    for &(cell, indicator) in top {
        let (lo, hi) = grid.rect(cell);
        let center = [f64::midpoint(lo[0], hi[0]), f64::midpoint(lo[1], hi[1])];
        let h = grid.cell_h(cell);
        let boundary_distance =
            ((center[0] - hole.center[0]).hypot(center[1] - hole.center[1]) - hole.radius).abs();
        let cell_reach = f64::sqrt(2.0) * 0.5 * h;
        assert!(
            boundary_distance <= k * h + cell_reach,
            "top-decile cell {cell:?} lies {:.3} cell widths beyond the k={k} hole band",
            (boundary_distance - cell_reach).max(0.0) / h
        );
        top_mass += indicator;
    }
    let mass_fraction = top_mass / localization_mass;
    assert!(
        mass_fraction.is_finite() && mass_fraction >= 0.35,
        "top indicator decile must carry visible marking mass: {mass_fraction:.3}"
    );
    println!(
        "{{\"suite\":\"fs-dwr/elasticity\",\"case\":\"plate-with-hole-localization\",\"indicator_terms\":\"bulk+nitsche+ghost\",\"top_decile_cells\":{},\"k_cells\":{k},\"top_mass_fraction\":{mass_fraction:.8}}}",
        top.len()
    );
}

#[test]
fn edge_band_traction_disjoint_cut_and_refusal_semantics() {
    let material = material();
    let clamp_left = |x: f64, _y: f64| x == 0.0;
    let body_zero = |_: f64, _: f64| [0.0, 0.0];
    let zero = |_: f64, _: f64| [0.0, 0.0];
    let traction_fn = |_x: f64, y: f64| [0.05 * (PI * y).sin().powi(2), 0.0];

    let grid = Quadtree::uniform(3);

    // Case 1: Cut intersects bottom edge (y = 0), but traction is on right edge (x = 1).
    // Disjoint cut edge must succeed under EdgeBand.
    let bottom_notch = HalfPlane {
        normal: [0.0, -1.0],
        offset: -0.1,
    };
    let right_band = EdgeBand::new(DesignBoxEdge::Right, 0.25, 0.75).expect("valid edge band");
    let prob = problem(
        &grid,
        &bottom_notch,
        &material,
        Some(&clamp_left),
        None,
        true,
        0.5,
    );
    let estimate = estimate_elasticity_compliance_with_boundary_traction(
        &prob,
        &body_zero,
        &zero,
        BoundaryTraction::EdgeBand {
            support: right_band,
            value: &traction_fn,
        },
    )
    .expect("disjoint cut with edge band must succeed");
    assert_reconstruction(&estimate);
    assert!(estimate.terms.outer_traction.abs() > 0.0);
    assert!(estimate.terms.bulk.abs() > 0.0);

    // Case 2: Cut intersects the loaded band directly (x = 1, y in [0.25, 0.75]).
    // Cut through loaded band must fail-closed with InvalidElasticityInput.
    let sqrt2 = 2.0f64.sqrt();
    let right_notch = HalfPlane {
        normal: [1.0 / sqrt2, 1.0 / sqrt2],
        offset: 1.5 / sqrt2,
    };
    let prob_cut = problem(
        &grid,
        &right_notch,
        &material,
        Some(&clamp_left),
        None,
        true,
        0.5,
    );
    let result_cut = estimate_elasticity_compliance_with_boundary_traction(
        &prob_cut,
        &body_zero,
        &zero,
        BoundaryTraction::EdgeBand {
            support: right_band,
            value: &traction_fn,
        },
    );
    assert!(
        matches!(result_cut, Err(CutFemError::InvalidElasticityInput { .. })),
        "cut through loaded edge band must refuse"
    );

    // Case 3: Competing boundary traction sources (both problem field and explicit arg) refuse.
    let prob_competing = problem(
        &grid,
        &AlwaysInside,
        &material,
        Some(&clamp_left),
        Some(&traction_fn),
        true,
        0.5,
    );
    let result_competing = estimate_elasticity_compliance_with_boundary_traction(
        &prob_competing,
        &body_zero,
        &zero,
        BoundaryTraction::EdgeBand {
            support: right_band,
            value: &traction_fn,
        },
    );
    assert!(
        matches!(
            result_competing,
            Err(CutFemError::InvalidElasticityInput { .. })
        ),
        "competing boundary traction sources must refuse"
    );
}

#[test]
fn two_cycle_graded_band_journey_is_deterministic_and_refines() {
    // frankensim-v6dn.2.3: two authentic estimate -> band-refine -> graded
    // vector re-solve cycles through the real public path, with bit-exact
    // deterministic evidence across a full replay of the journey.
    let material = material();
    let disk = Circle {
        center: [0.5, 0.5],
        radius: 0.35,
    };
    let body = |x: f64, y: f64| [8.0 * (y - 0.5), 12.0 * (0.5 - x)];
    let zero_embedded = |_: f64, _: f64| [0.0, 0.0];
    let theta = 0.5;

    let reference = {
        let grid = Quadtree::uniform(5);
        let cut = problem(&grid, &disk, &material, None, None, false, 0.5);
        estimate_elasticity_compliance(&cut, &body, &zero_embedded)
            .expect("uniform reference estimate")
            .j_primal
    };

    let run_journey = || -> Vec<(usize, f64, f64, usize, usize)> {
        let mut grid = Quadtree::with_room(3, 6);
        grid.refine_toward_interface(&disk, 3);
        let mut rows = Vec::new();
        for cycle in 0..2 {
            let cut = problem(&grid, &disk, &material, None, None, false, 0.5);
            let solution = cut.solve(&body, &zero_embedded).expect("graded solve");
            let est = estimate_elasticity_compliance(&cut, &body, &zero_embedded)
                .expect("cycle estimate");
            rows.push((
                solution.dof_count(),
                est.j_primal,
                est.eta_abs,
                grid.leaf_count(),
                cycle,
            ));
            if cycle + 1 < 2 {
                let marked = dorfler(&est.indicators, theta);
                assert!(
                    !marked.is_empty(),
                    "cycle {cycle} must mark a positive real split set"
                );
                for cell in &marked {
                    if grid.is_leaf(*cell) {
                        grid.split(*cell);
                    }
                }
                grid.balance();
                assert!(
                    grid.leaf_count() > rows[cycle].3,
                    "cycle {cycle} refinement must grow the leaf set"
                );
            }
        }
        let levels: BTreeSet<u32> = grid.leaves().map(|(level, _, _)| level).collect();
        assert!(
            levels.len() >= 2,
            "the adapted grid must retain mixed leaves: {levels:?}"
        );
        rows
    };

    let run_a = run_journey();
    let run_b = run_journey();
    assert_eq!(run_a, run_b, "the journey must replay bit-exact");
    assert_eq!(run_a.len(), 2, "exactly two authentic cycles are required");

    let (dofs_initial, j_initial, _, _, _) = run_a[0];
    let (dofs_final, j_final, eta_final, _, _) = run_a[1];
    assert!(
        dofs_final > dofs_initial,
        "adapted re-solves must refine: {dofs_initial} -> {dofs_final}"
    );
    let initial_gap = (j_initial - reference).abs();
    let final_gap = (j_final - reference).abs();
    assert!(
        final_gap <= initial_gap,
        "two cycles must move compliance toward the uniform reference: \
         initial {initial_gap:e}, final {final_gap:e}"
    );
    assert!(
        eta_final.is_finite() && eta_final >= 0.0,
        "final marking mass must stay finite and nonnegative"
    );
    println!(
        "{{\"suite\":\"fs-dwr/elasticity\",\"case\":\"two-cycle-graded-band-journey\",\
         \"reference\":{reference:.12e},\"cycles\":[{{\"dofs\":{dofs_initial},\
         \"j\":{:016x},\"eta\":{:016x}}},{{\"dofs\":{dofs_final},\"j\":{:016x},\
         \"eta\":{:016x}}}],\"initial_gap\":{initial_gap:.6e},\
         \"final_gap\":{final_gap:.6e},\"no_claim\":\"2-D Q1 graded elasticity only; no 3D, \
         higher-order, aggregation, near-incompressibility, matrix-free execution, or certified \
         DWR bound claims\"}}",
        j_initial.to_bits(),
        run_a[0].2.to_bits(),
        j_final.to_bits(),
        run_a[1].2.to_bits(),
    );
}
