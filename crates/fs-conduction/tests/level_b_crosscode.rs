//! Level-B cross-code comparison lane (bead
//! `frankensim-extreal-program-f85xj.4.3`).
//!
//! Every `fs-vvreg` Level-B case is rebuilt from its typed catalog
//! definition, solved with THIS crate's assembly/solver stack, and
//! compared probe-by-probe against the frozen external reference
//! (scikit-fem + SuperLU, pinned by `tools/vvref/uv.lock`). The lane
//! also recomputes the canonical BLAKE3 mesh-identity hash and requires
//! it to match the one the external run recorded, which proves both
//! codes assembled on the bit-identical Kuhn mesh.
//!
//! What agreement means: the two implementations share NO code, but by
//! construction they solve the SAME discrete system (P1, element-mean
//! `k(T)`, consistent source/Robin mass). A within-envelope delta is
//! therefore an independent-implementation check on assembly, boundary
//! handling, and linear/nonlinear solving — it is not a discretization
//! or physical-validation claim (fs-vvreg registers these rows as
//! Estimated-colour `CrossCode`). A delta beyond a declared envelope
//! opens an investigation bead; envelopes are never silently widened.

mod support;

use fs_blake3::hash_bytes;
use fs_conduction::bc::{ThermalBc, ThermalBoundaryBuilder};
use fs_conduction::field::ScalarField;
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::material::{ConductivityModel, ConductivityTable};
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::solve::{
    ConductionProblem, ConductionReport, InitialGuess, LinearConfig, Nonlinearity, SolveConfig,
    StopRule, solve,
};
use fs_vvreg::thermal_level_b::{
    ThermalLevelBBcKind, ThermalLevelBCase, ThermalLevelBMaterial, ThermalLevelBReference,
    ThermalLevelBSource, thermal_level_b_case, thermal_level_b_cases, thermal_level_b_references,
};
use support::{json_escape, with_cx};

/// Canonical mesh-identity bytes, mirrored exactly by
/// `tools/vvref/solve_skfem.py::mesh_blake3`.
fn mesh_identity_hash(counts: [usize; 3], positions: &[[f64; 3]], tets: &[[u32; 4]]) -> String {
    let mut payload = Vec::with_capacity(32 + positions.len() * 24 + tets.len() * 16);
    payload.extend_from_slice(b"frankensim-vvref-mesh-v1\n");
    for count in counts {
        payload.extend_from_slice(&(count as u64).to_le_bytes());
    }
    payload.extend_from_slice(&(positions.len() as u64).to_le_bytes());
    for position in positions {
        for value in position {
            payload.extend_from_slice(&value.to_le_bytes());
        }
    }
    payload.extend_from_slice(&(tets.len() as u64).to_le_bytes());
    for tet in tets {
        for vertex in tet {
            payload.extend_from_slice(&vertex.to_le_bytes());
        }
    }
    hash_bytes(&payload).to_hex()
}

/// P1 barycentric interpolation at an arbitrary point: locate a tet
/// whose barycentric coordinates are all nonnegative (within `1e-12`)
/// and blend its four nodal values.
fn interpolate_at(
    tets: &[[u32; 4]],
    positions: &[[f64; 3]],
    values: &[f64],
    point: [f64; 3],
) -> Option<f64> {
    for tet in tets {
        let a = positions[tet[0] as usize];
        let b = positions[tet[1] as usize];
        let c = positions[tet[2] as usize];
        let d = positions[tet[3] as usize];
        let columns = [
            [b[0] - a[0], c[0] - a[0], d[0] - a[0]],
            [b[1] - a[1], c[1] - a[1], d[1] - a[1]],
            [b[2] - a[2], c[2] - a[2], d[2] - a[2]],
        ];
        let rhs = [point[0] - a[0], point[1] - a[1], point[2] - a[2]];
        let det = |m: &[[f64; 3]; 3]| -> f64 {
            m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
        };
        let d0 = det(&columns);
        if d0.abs() < 1e-300 {
            continue;
        }
        let mut replaced = columns;
        let mut bary = [0.0f64; 3];
        for axis in 0..3 {
            for row in 0..3 {
                replaced[row] = columns[row];
                replaced[row][axis] = rhs[row];
            }
            bary[axis] = det(&replaced) / d0;
            for row in 0..3 {
                replaced[row] = columns[row];
            }
        }
        let bary0 = 1.0 - bary[0] - bary[1] - bary[2];
        let all = [bary0, bary[0], bary[1], bary[2]];
        if all.iter().all(|w| *w >= -1e-12) {
            let mut out = 0.0;
            for (weight, vertex) in all.iter().zip(tet) {
                out += weight * values[*vertex as usize];
            }
            return Some(out);
        }
    }
    None
}

struct CaseRun {
    mesh_hash: String,
    positions: Vec<[f64; 3]>,
    tets: Vec<[u32; 4]>,
    temperature: Vec<f64>,
    report: ConductionReport,
}

fn run_case(case: &ThermalLevelBCase) -> CaseRun {
    let (complex, positions) = box_grid(case.mesh_counts, case.mesh_extent);
    let tets = complex.tets.clone();
    let mesh_hash = mesh_identity_hash(case.mesh_counts, &positions, &tets);
    let mesh = ConductionMesh::new(complex, positions.clone()).expect("mesh");

    let material = match case.material {
        ThermalLevelBMaterial::Isotropic { k } => {
            ConductivityModel::isotropic_declared(k).expect("isotropic material")
        }
        ThermalLevelBMaterial::Tensor { k } => {
            ConductivityModel::constant_tensor(k).expect("tensor material")
        }
        ThermalLevelBMaterial::LinearKt { knots } => ConductivityModel::isotropic(
            ConductivityTable::declared_curve(knots.to_vec()).expect("k(T) curve"),
        ),
    };

    let source = match case.source {
        ThermalLevelBSource::None => ScalarField::Uniform(0.0),
        ThermalLevelBSource::PolyXy { q0 } => {
            let ex = case.mesh_extent[0];
            let ey = case.mesh_extent[1];
            let values = positions
                .iter()
                // Exact op order shared with the deck runner:
                // q0 * (x / ex) * (1.0 - y / ey)
                .map(|p| q0 * (p[0] / ex) * (1.0 - p[1] / ey))
                .collect();
            ScalarField::nodal("level-b-source", positions.len(), values).expect("source field")
        }
    };

    let mut builder = ThermalBoundaryBuilder::new(&mesh);
    for bc in case.bcs {
        let axis = bc.axis;
        let target = if bc.at_max {
            case.mesh_extent[axis]
        } else {
            0.0
        };
        let condition = match bc.kind {
            ThermalLevelBBcKind::Dirichlet { t_k } => {
                ThermalBc::dirichlet(t_k).expect("dirichlet bc")
            }
            ThermalLevelBBcKind::Robin { h, t_inf_k } => {
                ThermalBc::robin(h, t_inf_k).expect("robin bc")
            }
        };
        builder = builder
            .region(
                bc.name,
                |f| on_box_face(f.centroid[axis], target),
                condition,
            )
            .expect("region");
    }
    let boundary = builder.adiabatic_remainder().finish().expect("boundary");

    let config = SolveConfig {
        nonlinearity: Nonlinearity::FixedPoint {
            relaxation: 1.0,
            max_backtracks: 8,
        },
        stop: StopRule {
            residual_rtol: 1e-12,
            residual_atol: 1e-24,
            step_atol: 0.0,
            max_iterations: 60,
        },
        linear: LinearConfig {
            tolerance: 1e-13,
            max_iterations: 60_000,
            restart: 60,
        },
        initial: InitialGuess::DirichletMean,
    };
    let solution = with_cx(|cx| {
        solve(
            cx,
            ConductionProblem { element_materials: None, mesh: &mesh,
            boundary: &boundary,
            material: &material,
            source: &source, },
            config,
        )
    })
    .expect("solve");

    CaseRun {
        mesh_hash,
        positions,
        tets,
        temperature: solution.temperature,
        report: solution.report,
    }
}

fn reference_for(case_id: &str) -> &'static ThermalLevelBReference {
    thermal_level_b_references()
        .expect("committed Level-B manifest must verify")
        .iter()
        .find(|reference| reference.case_id == case_id)
        .expect("verification guarantees a block per catalog case")
}

fn vertex_index(case: &ThermalLevelBCase, grid: [usize; 3]) -> usize {
    let py = case.mesh_counts[1] + 1;
    let pz = case.mesh_counts[2] + 1;
    grid[0] * py * pz + grid[1] * pz + grid[2]
}

#[test]
fn every_level_b_case_matches_its_external_reference_within_envelope() {
    for case in thermal_level_b_cases() {
        let reference = reference_for(case.id);
        let run = run_case(case);

        // Both codes assembled on the bit-identical mesh, or the
        // comparison is meaningless.
        assert_eq!(
            run.mesh_hash, reference.mesh_blake3,
            "{}: mesh identity diverged from the external run",
            case.id
        );

        for probe in &reference.probes {
            let vertex = vertex_index(case, probe.grid);
            let ours = run.temperature[vertex];
            let delta = (ours - probe.temperature_k).abs();
            let verdict = if delta <= case.acceptance_atol_k {
                "pass"
            } else {
                "fail"
            };
            println!(
                "{{\"suite\":\"fs-conduction/level-b-crosscode\",\"case\":\"{}\",\
                 \"probe\":{},\"external\":\"{}\",\"ours_k\":{ours},\
                 \"reference_k\":{},\"delta_k\":{delta},\"envelope_k\":{},\
                 \"verdict\":\"{verdict}\",\
                 \"authority\":\"estimated-cross-code-agreement-not-truth\"}}",
                json_escape(case.id),
                probe.index,
                json_escape(&reference.external_code),
                probe.temperature_k,
                case.acceptance_atol_k,
            );
            assert!(
                delta <= case.acceptance_atol_k,
                "{} probe {}: |{ours} - {}| = {delta} K exceeds the declared \
                 envelope {} K; open an investigation bead instead of widening",
                case.id,
                probe.index,
                probe.temperature_k,
                case.acceptance_atol_k
            );

            // The interpolation utility must agree with the nodal value
            // at a node, so probe extraction and interpolation cannot
            // drift apart.
            let interpolated = interpolate_at(
                &run.tets,
                &run.positions,
                &run.temperature,
                run.positions[vertex],
            )
            .expect("probe node sits inside the mesh");
            assert!(
                (interpolated - ours).abs() <= 1e-9,
                "{} probe {}: interpolation at the node disagrees with the nodal value",
                case.id,
                probe.index
            );
        }

        // The discrete field the externals froze bounds our probes too:
        // parity this tight means our field range cannot leave theirs by
        // more than the envelope.
        for probe in &reference.probes {
            let ours = run.temperature[vertex_index(case, probe.grid)];
            assert!(
                ours >= reference.t_min_k - case.acceptance_atol_k
                    && ours <= reference.t_max_k + case.acceptance_atol_k,
                "{}: probe value {ours} outside the external field range",
                case.id
            );
        }

        // Energy balance is this solver's own honesty check; a lane that
        // compares temperatures but ships an unbalanced solve is not a
        // cross-code result worth citing.
        let closure = run.report.energy.relative_closure();
        assert!(
            closure <= 1e-6,
            "{}: energy closure {closure} too loose to cite",
            case.id
        );
        println!(
            "{{\"suite\":\"fs-conduction/level-b-crosscode\",\"case\":\"{}\",\
             \"mesh_blake3\":\"{}\",\"iterations\":{},\"energy_closure\":{closure},\
             \"verdict\":\"pass\"}}",
            json_escape(case.id),
            run.mesh_hash,
            run.report.iterations,
        );
    }
}

#[test]
fn interpolation_reproduces_an_affine_field_exactly() {
    // Synthetic-field correctness for the interpolation utility itself:
    // P1 interpolation is exact on affine data, so any drift is a bug in
    // the point location or the barycentric blend.
    let case = thermal_level_b_case("thermal-b-block-source-robin-v1").expect("case");
    let (complex, positions) = box_grid(case.mesh_counts, case.mesh_extent);
    let affine = |p: [f64; 3]| 300.0 + 17.0 * p[0] - 4.5 * p[1] + 2.25 * p[2];
    let values: Vec<f64> = positions.iter().map(|p| affine(*p)).collect();
    let probes = [
        [0.0123, 0.0087, 0.0111],
        [0.0034, 0.0491, 0.0002],
        [0.0777, 0.0009, 0.0299],
        [0.04, 0.025, 0.015],
    ];
    for point in probes {
        let interpolated = interpolate_at(&complex.tets, &positions, &values, point)
            .expect("interior point locates");
        let exact = affine(point);
        assert!(
            (interpolated - exact).abs() <= 1e-9,
            "affine reproduction failed at {point:?}: {interpolated} vs {exact}"
        );
    }
    // A point outside the box must refuse, not extrapolate.
    assert!(interpolate_at(&complex.tets, &positions, &values, [1.0, 1.0, 1.0]).is_none());
}

#[test]
fn nonlinear_case_exercises_a_real_picard_iteration_here_too() {
    let case = thermal_level_b_case("thermal-b-kt-nonlinear-slab-v1").expect("case");
    let run = run_case(case);
    assert!(
        run.report.iterations > 1,
        "the k(T) case must not converge in a single fixed-point step"
    );
}
