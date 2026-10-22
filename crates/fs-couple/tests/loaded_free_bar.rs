//! M2.5-W4a G3 (IRP bead `frankensim-irp-v5-instrument-realism-pexfp.105`):
//! loaded modes of a real free aluminium bar in air.
//!
//! Fixture data, never authored:
//! - bar: clean 2024-T4 beam 40 of NASA TM X-2852 (Howell 1974), specimen
//!   `aluminum-beams-howell1974` (length, width, measured thickness);
//! - material: the 2024-T3 matdb set of NASA TN D-6448 (E, nu, rho). The
//!   report states no elastic constants for its T4 beams, so this is a
//!   declared temper substitution;
//! - air: `Medium::air()`, the declared room state (the report states none).
//!
//! The structure is the existing TET10 solid, free-free. Radiation loading is
//! the routed modal bake (`modal_radiation_loading_routed`) on the boundary
//! faces of the same Kuhn mesh, each panel carrying the vertex average of
//! u.n, with the formulation `GeometryPolicy` selects (plain CBIE for this
//! thin solid). The bake keeps its unchanged 1e-12 reciprocity gate, which a
//! single-mode matrix passes by construction; multi-mode collocation data
//! do not (the reciprocity decision on .105).
//!
//! The decay comparison with the measured air-minus-vacuum decrement is a
//! one-sided bound only: that difference holds radiation and viscous air loss
//! together (`unknowns.tsv`: u-air-loss-split).
use fs_bem::helmholtz::{self, Formulation, Medium, with_dense_workers};
use fs_bem::modal_loading::{LoadingRoute, SurfaceModeTraces, modal_radiation_loading_routed};
use fs_bem::panel3d::SpherePanels;
use fs_bem::radiation_policy::GeometryPolicy;
use fs_exec::Cx;
use fs_la::eigen::householder_eigh;
use fs_la::eigen_complex::{lu_complex, with_lu_workers};
use fs_la::factor::cholesky;
use fs_math::c64::C64;
use fs_modal::{ModePair, SliceOptions, slice_window};
use fs_solid::linear3::tet10::{QuadraticTetProblem, TetQuadrature, promote_p1_to_p2};
use fs_solid::linear3::{
    TetAssemblyBudget, TetElasticAssembly, TetElasticMaterial, TetMaterialField,
};
use std::time::Instant;

const DIMENSIONS: &str =
    include_str!("../../../data/instruments/geometry/aluminum-beams-howell1974/dimensions.tsv");
const DECAY: &str = include_str!(
    "../../../data/instruments/descriptions/bars/aluminum-beam-free-decay-howell1974.tsv"
);
const ALLOY: &str =
    include_str!("../../../data/matdb/seed-v1/aluminum-2024-t3-nasa-tn-d6448/properties.tsv");
const BEAM: &str = "beam-40";
const TAU: f64 = core::f64::consts::TAU;

fn dimension(part: &str, quantity: &str) -> f64 {
    DIMENSIONS
        .lines()
        .skip(2)
        .find_map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            (f[0] == part && f[1] == "prism" && f[2] == quantity)
                .then(|| f[3].parse::<f64>().expect("specimen value"))
        })
        .unwrap_or_else(|| panic!("specimen {part}/{quantity} missing"))
}

fn scalar(table: &str, name: &str, column: usize) -> f64 {
    table
        .lines()
        .find_map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            (f.first() == Some(&"scalar") && f.get(column) == Some(&name))
                .then(|| f[4].parse::<f64>().expect("source scalar"))
        })
        .unwrap_or_else(|| panic!("source scalar {name} missing"))
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = fs_exec::CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 105,
                kernel_id: 105,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn workers() -> usize {
    std::thread::available_parallelism().map_or(1, |n| n.get().min(16))
}

fn material() -> TetElasticMaterial {
    let identity = fs_blake3::hash_domain(
        "org.frankensim.105.g3.material.v1",
        format!("{ALLOY}\ntemper-substitution=2024-T3 set for the 2024-T4 beams").as_bytes(),
    );
    TetElasticMaterial::try_new(
        scalar(ALLOY, "density", 3),
        scalar(ALLOY, "young_modulus", 3) * 1e6,
        scalar(ALLOY, "poisson_ratio", 3),
        identity.0,
    )
    .expect("isotropic aluminium")
}

struct Bar {
    extent: [f64; 3],
    counts: [usize; 3],
    nodes: Vec<[f64; 3]>,
    assembly: TetElasticAssembly,
}

/// TET10 free-free bar on a Kuhn grid; promotion keeps the P1 vertex ids.
fn bar(counts: [usize; 3], cx: &Cx<'_>) -> Bar {
    let extent = [
        dimension(BEAM, "length"),
        dimension(BEAM, "width"),
        dimension(BEAM, "thickness"),
    ];
    let (complex, p1) = fs_conduction::fixtures::structured_tets(counts, &|i, j, k| {
        [
            extent[0] * i as f64 / counts[0] as f64,
            extent[1] * j as f64 / counts[1] as f64,
            extent[2] * k as f64 / counts[2] as f64,
        ]
    });
    let tets: Vec<[usize; 4]> = complex.tets.iter().map(|t| t.map(|i| i as usize)).collect();
    // A prism: straight mid-edge nodes already lie on the faces.
    let (nodes, elements) = promote_p1_to_p2(&p1, &tets, None, cx).expect("TET10 promotion");
    let mat = material();
    let assembly = QuadraticTetProblem {
        nodes_m: &nodes,
        elements: &elements,
        materials: TetMaterialField::Uniform(&mat),
        fixed_dofs: &[],
        budget: TetAssemblyBudget::standard(),
        quadrature: TetQuadrature {
            stiffness_points_per_dir: 3,
            mass_points_per_dir: 4,
        },
    }
    .assemble(cx)
    .expect("free-free TET10 assembly");
    Bar {
        extent,
        counts,
        nodes,
        assembly,
    }
}

/// Outward boundary triangles of a Kuhn grid, by P1 vertex index, using the
/// same low/low to high/high quad diagonal as the tetrahedra.
fn boundary(counts: [usize; 3]) -> Vec<[usize; 3]> {
    let index = |q: [usize; 3]| (q[0] * (counts[1] + 1) + q[1]) * (counts[2] + 1) + q[2];
    let mut faces = Vec::new();
    for axis in 0..3 {
        let other: Vec<_> = (0..3).filter(|&a| a != axis).collect();
        for high in [false, true] {
            for i in 0..counts[other[0]] {
                for j in 0..counts[other[1]] {
                    let mut q = [0; 3];
                    q[axis] = if high { counts[axis] } else { 0 };
                    q[other[0]] = i;
                    q[other[1]] = j;
                    let a = index(q);
                    q[other[0]] += 1;
                    let b = index(q);
                    q[other[1]] += 1;
                    let c = index(q);
                    q[other[0]] -= 1;
                    let d = index(q);
                    for face in [[a, b, c], [a, c, d]] {
                        faces.push((face, axis, high));
                    }
                }
            }
        }
    }
    faces
        .into_iter()
        .map(|(mut face, axis, high)| {
            // Orientation from the face's own normal component along `axis`.
            let p = |v: usize| {
                let k = v % (counts[2] + 1);
                let j = v / (counts[2] + 1) % (counts[1] + 1);
                let i = v / ((counts[1] + 1) * (counts[2] + 1));
                [i as f64, j as f64, k as f64]
            };
            let (p0, p1, p2) = (p(face[0]), p(face[1]), p(face[2]));
            let ab: [f64; 3] = core::array::from_fn(|k| p1[k] - p0[k]);
            let ac: [f64; 3] = core::array::from_fn(|k| p2[k] - p0[k]);
            let normal =
                ab[(axis + 1) % 3] * ac[(axis + 2) % 3] - ab[(axis + 2) % 3] * ac[(axis + 1) % 3];
            if (normal > 0.0) != high {
                face.swap(1, 2);
            }
            face
        })
        .collect()
}

/// A BEM surface on the Kuhn grid `counts`, nested in the structural grid,
/// and the map from its P1 vertices to the structural P1 vertex ids.
fn surface(bar: &Bar, counts: [usize; 3]) -> (SpherePanels, Vec<[usize; 3]>) {
    let fine = bar.counts;
    assert!(
        fine.iter().zip(counts).all(|(a, b)| a % b == 0),
        "nested BEM grid"
    );
    let map = |v: usize| {
        let k = v % (counts[2] + 1);
        let j = v / (counts[2] + 1) % (counts[1] + 1);
        let i = v / ((counts[1] + 1) * (counts[2] + 1));
        ((i * fine[0] / counts[0]) * (fine[1] + 1) + j * fine[1] / counts[1]) * (fine[2] + 1)
            + k * fine[2] / counts[2]
    };
    let faces: Vec<[usize; 3]> = boundary(counts).into_iter().map(|f| f.map(map)).collect();
    let triangles = faces.iter().map(|f| f.map(|v| bar.nodes[v])).collect();
    (
        SpherePanels::from_triangles(triangles).expect("closed bar surface"),
        faces,
    )
}

/// Panel normal velocity of a unit modal velocity: the vertex average of u.n.
fn trace(bar: &Bar, panels: &SpherePanels, faces: &[[usize; 3]], phi: &[f64]) -> Vec<f64> {
    let mut full = vec![0.0; 3 * bar.nodes.len()];
    for (row, &dof) in bar.assembly.free_dofs.iter().enumerate() {
        full[dof] = phi[row];
    }
    faces
        .iter()
        .zip(panels.normals())
        .map(|(face, n)| {
            face.iter()
                .map(|&v| (0..3).map(|c| full[3 * v + c] * n[c]).sum::<f64>())
                .sum::<f64>()
                / 3.0
        })
        .collect()
}

fn hz(m: &ModePair) -> f64 {
    m.lambda.max(0.0).sqrt() / TAU
}

/// Every in-vacuo mode of a dense pencil (K, M): M = L L^T, C = L^-1 K L^-T,
/// C = V diag(lambda) V^T by Householder-QL, phi = L^-T v (so phi^T M phi = 1).
/// Ascending, rigid modes included.
fn complete_basis(k: &[f64], m: &[f64], n: usize) -> (Vec<f64>, Vec<Vec<f64>>) {
    let chol = cholesky(m, n).expect("SPD consistent mass");
    let l = |i: usize, j: usize| chol.l(i, j);
    let forward = |b: &mut [f64]| {
        for i in 0..n {
            let mut acc = b[i];
            for j in 0..i {
                acc = l(i, j).mul_add(-b[j], acc);
            }
            b[i] = acc / l(i, i);
        }
    };
    let mut y = vec![0.0; n * n];
    let mut column = vec![0.0; n];
    for c in 0..n {
        for i in 0..n {
            column[i] = k[i * n + c];
        }
        forward(&mut column);
        for i in 0..n {
            y[c * n + i] = column[i];
        }
    }
    let mut cmat = vec![0.0; n * n];
    for c in 0..n {
        // Column c of C = L^-1 Y^T e_c, i.e. L^-1 applied to row c of Y.
        for i in 0..n {
            column[i] = y[i * n + c];
        }
        forward(&mut column);
        for i in 0..n {
            cmat[i * n + c] = column[i];
        }
    }
    for i in 0..n {
        for j in 0..i {
            let avg = f64::midpoint(cmat[i * n + j], cmat[j * n + i]);
            cmat[i * n + j] = avg;
            cmat[j * n + i] = avg;
        }
    }
    let (values, vectors) = householder_eigh(&cmat, n);
    let shapes = (0..n)
        .map(|c| {
            let mut phi: Vec<f64> = (0..n).map(|i| vectors[i * n + c]).collect();
            for i in (0..n).rev() {
                let mut acc = phi[i];
                for j in i + 1..n {
                    acc = l(j, i).mul_add(-phi[j], acc);
                }
                phi[i] = acc / l(i, i);
            }
            phi
        })
        .collect();
    (values, shapes)
}

/// The loaded root of w^2 + i z w - W^2 = 0 continuing the in-vacuo W.
fn single_mode_root(big_w: f64, z: C64) -> C64 {
    let i = C64::new(0.0, 1.0);
    let disc = (C64::from_re(4.0 * big_w * big_w) - z * z).sqrt();
    let root = (disc - i * z).scale(0.5);
    if root.re > 0.0 {
        root
    } else {
        (-disc - i * z).scale(0.5)
    }
}

/// Window route: one routed bake per mode (reciprocal by construction), its
/// impedance taken at the real part of the current loaded frequency.
fn window_pole(
    panels: &SpherePanels,
    shape: &[f64],
    big_w: f64,
    medium: Medium,
) -> (C64, C64, f64) {
    let shapes = vec![shape.to_vec()];
    let mut omega = C64::from_re(big_w);
    let mut z = C64::ZERO;
    let mut min_eig = f64::NAN;
    for _ in 0..6 {
        let bake = with_dense_workers(workers(), || {
            modal_radiation_loading_routed(
                SurfaceModeTraces::ClosedBody {
                    surface: panels,
                    shapes: &shapes,
                },
                LoadingRoute::ClosedBodyGeometryPolicy,
                &[omega.re],
                medium,
                fs_blake3::hash_domain("analytic-fixture", b"105-g3-howell-beam-40"),
                1,
                None,
                || false,
            )
        })
        .expect("single-mode routed bake (reciprocal by construction)");
        let row = &bake.rows()[0];
        z = row.impedance()[0];
        min_eig = row.minimum_resistance_eigenvalue();
        let next = single_mode_root(big_w, z);
        let change = (next - omega).abs() / big_w;
        omega = next;
        if change < 1e-13 {
            break;
        }
    }
    (omega, z, min_eig)
}

/// Complete-basis oracle: the coupled FE+BEM system in ALL N in-vacuo modes of
/// the same FE model (rigid modes included), which is the direct coupled solve
/// in another basis. Z is the same BEM owner's complete N x N matrix; the loaded
/// eigenvalue nearest W_r^2 follows by inverse iteration, the frequency by a
/// fixed point on w^2 = mu(w).
fn complete_pole(
    panels: &SpherePanels,
    traces: &[Vec<f64>],
    lambdas: &[f64],
    target: usize,
    medium: Medium,
) -> (C64, f64) {
    let n = traces.len();
    let policy = GeometryPolicy::new(panels).expect("closed bar");
    let mut omega = C64::from_re(lambdas[target].sqrt());
    let mut defect = 0.0;
    for _ in 0..6 {
        let k = omega.re / medium.sound_speed;
        let formulation = policy.formulation(k).expect("policy");
        let mut z = vec![C64::ZERO; n * n];
        for start in (0..n).step_by(helmholtz::MAX_RADIATION_FIELDS_PER_BATCH) {
            let end = (start + helmholtz::MAX_RADIATION_FIELDS_PER_BATCH).min(n);
            let fields: Vec<Vec<C64>> = traces[start..end]
                .iter()
                .map(|t| t.iter().map(|&v| C64::from_re(v)).collect())
                .collect();
            let refs: Vec<&[C64]> = fields.iter().map(Vec::as_slice).collect();
            let solved = with_dense_workers(workers(), || {
                helmholtz::solve_radiation_batch(panels, k, medium, &refs, formulation)
            })
            .expect("complete-basis BEM batch");
            for (column, solution) in (start..end).zip(solved) {
                for (row, t) in traces.iter().enumerate() {
                    let mut force = C64::ZERO;
                    for ((&phi, &p), &area) in t.iter().zip(&solution.pressure).zip(panels.areas())
                    {
                        force = force + p.scale(phi * area);
                    }
                    z[row * n + column] = force;
                }
            }
        }
        let norm = z.iter().fold(0.0_f64, |s, v| s.hypot(v.abs()));
        defect = (0..n)
            .flat_map(|i| (0..n).map(move |j| (i, j)))
            .fold(0.0_f64, |s, (i, j)| {
                s.hypot((z[i * n + j] - z[j * n + i]).abs())
            })
            / norm;
        // A(w) = diag(Lambda) - i w Z(w); eigenvalue mu nearest w^2.
        let i_omega = C64::new(0.0, 1.0) * omega;
        let shift = omega * omega;
        let mut a = vec![C64::ZERO; n * n];
        for r in 0..n {
            for s in 0..n {
                a[r * n + s] = -(i_omega * z[r * n + s]);
            }
            a[r * n + r] = a[r * n + r] + C64::from_re(lambdas[r]) - shift;
        }
        let lu = with_lu_workers(workers(), || lu_complex(&a, n)).expect("shifted pencil");
        let mut x = vec![C64::ZERO; n];
        x[target] = C64::ONE;
        let mut mu = shift;
        for _ in 0..8 {
            let mut y = x.clone();
            lu.solve(&mut y);
            // Rayleigh-type ratio of the shifted inverse.
            let num = x.iter().fold(C64::ZERO, |s, v| s + v.conj() * *v);
            let den = x
                .iter()
                .zip(&y)
                .fold(C64::ZERO, |s, (u, v)| s + u.conj() * *v);
            let next = shift + num / den;
            let scale = y.iter().fold(0.0_f64, |s, v| s.max(v.abs()));
            x = y.into_iter().map(|v| v.scale(1.0 / scale)).collect();
            let done = (next - mu).abs() < 1e-14 * shift.abs();
            mu = next;
            if done {
                break;
            }
        }
        let next = mu.sqrt();
        let next = if next.re > 0.0 { next } else { -next };
        let change = (next - omega).abs() / omega.abs();
        omega = next;
        if change < 1e-13 {
            break;
        }
    }
    (omega, defect)
}

fn cents(a: f64, b: f64) -> f64 {
    1200.0 * (a / b).log2()
}

/// G3: loaded frequencies and decay rates of the window modes from the routed
/// single-mode bakes match the complete-basis coupled FE+BEM solve within 0.1
/// cent and 2%, on two FE/BEM discretizations; radiation damping >= 0 per
/// mode (the bake's passivity gate, plus a direct sign check).
#[test]
fn loaded_modes_match_coupled_solve() {
    let medium = Medium::air();
    for counts in [[12usize, 2, 1], [16, 4, 1]] {
        with_cx(|cx| {
            let started = Instant::now();
            let bar = bar(counts, cx);
            let n = bar.assembly.free_dofs.len();
            let (lambdas, shapes) = complete_basis(
                &bar.assembly.stiffness.to_dense(),
                &bar.assembly.mass.to_dense(),
                n,
            );
            let (panels, faces) = surface(&bar, counts);
            let traces: Vec<Vec<f64>> = shapes
                .iter()
                .map(|phi| trace(&bar, &panels, &faces, phi))
                .collect();
            let f = |r: usize| lambdas[r].max(0.0).sqrt() / TAU;
            assert_eq!(
                (0..n).filter(|&r| f(r) < 1.0).count(),
                6,
                "free-free: six rigid modes"
            );
            // Window: the elastic modes below 1600 Hz.
            let window: Vec<usize> = (6..n).filter(|&r| f(r) < 1600.0).collect();
            assert!(!window.is_empty());
            for &r in &window {
                let big_w = lambdas[r].sqrt();
                let (w_window, z, min_eig) = window_pole(&panels, &traces[r], big_w, medium);
                let (w_complete, defect) = complete_pole(&panels, &traces, &lambdas, r, medium);
                let df = cents(w_window.re, w_complete.re);
                let decay = -w_window.im;
                let decay_complete = -w_complete.im;
                let dd = (decay - decay_complete).abs() / decay_complete;
                println!(
                    "{{\"suite\":\"fs-couple\",\"case\":\"loaded-modes-match-coupled-solve\",\"inputs\":\"NASA TM X-2852 beam 40; matdb 2024-T3 TN D-6448 (temper substitution); Medium::air\",\"axis\":\"fe-bem-grid\",\"rung\":{counts:?},\"dofs\":{n},\"panels\":{},\"mode\":{r},\"f_vacuo_hz\":{},\"f_window_hz\":{},\"f_complete_hz\":{},\"cents\":{df:e},\"decay_window_1_s\":{decay:e},\"decay_complete_1_s\":{decay_complete:e},\"decay_rel\":{dd:e},\"z\":[{:e},{:e}],\"min_eig_re\":{min_eig:e},\"complete_reciprocity_defect\":{defect:e},\"tol\":{{\"cents\":0.1,\"decay\":0.02}},\"verdict\":\"{}\",\"elapsed_s\":{}}}",
                    panels.areas().len(),
                    big_w / TAU,
                    w_window.re / TAU,
                    w_complete.re / TAU,
                    z.re,
                    z.im,
                    if df.abs() <= 0.1 && dd <= 0.02 {
                        "pass"
                    } else {
                        "fail"
                    },
                    started.elapsed().as_secs_f64()
                );
                assert!(
                    z.re >= 0.0 && decay >= 0.0 && decay_complete >= 0.0,
                    "radiation damping >= 0"
                );
                assert!(df.abs() <= 0.1, "mode {r} at {counts:?}: {df} cents");
                assert!(dd <= 0.02, "mode {r} at {counts:?}: decay {dd}");
            }
        });
    }
}

/// G1 (resistance gated apart from the complex error): the radiation decay
/// rate of the fundamental under a nested BEM panel ladder on a fixed fine
/// TET10 structure, with its observed order and the measured air-minus-vacuum
/// decrement of beam 40 as a one-sided physical bound.
#[test]
fn loaded_free_bar_radiation_decay_ladder() {
    let medium = Medium::air();
    with_cx(|cx| {
        let started = Instant::now();
        let fine = [96usize, 16, 2];
        let bar = bar(fine, cx);
        let a = &bar.assembly;
        // Euler-Bernoulli first free-free root, 4.73004^2 / (2 pi L^2) * sqrt(EI / rho A),
        // only to place the certified window.
        let e = scalar(ALLOY, "young_modulus", 3) * 1e6;
        let rho = scalar(ALLOY, "density", 3);
        let [l, _, t] = bar.extent;
        let beam = 4.730_040_7_f64.powi(2) / (TAU * l * l) * (e * t * t / (12.0 * rho)).sqrt();
        let report = slice_window(
            &a.stiffness,
            &a.mass,
            ((TAU * 0.5 * beam).powi(2), (TAU * 1.5 * beam).powi(2)),
            &SliceOptions {
                max_lanczos: 2048,
                max_restarts: 48,
                ritz_tol: 1e-13,
                mass_diagonal_equilibration: true,
                ..SliceOptions::default()
            },
        )
        .expect("certified window around the fundamental");
        // The fundamental vertical bending mode: the one whose end section
        // moves along z (thickness), not y.
        let fundamental = report
            .modes
            .iter()
            .min_by(|x, y| x.lambda.total_cmp(&y.lambda))
            .expect("a mode in the window");
        let f_vacuo = hz(fundamental);
        let measured_vacuum = scalar(DECAY, "beam-40-frequency-3.0e-4atm", 1);
        let delta_air = scalar(DECAY, "beam-40-log-decrement-1atm", 1)
            - scalar(DECAY, "beam-40-log-decrement-3.0e-4atm", 1);
        let mut rows: Vec<(usize, f64, f64)> = Vec::new();
        for s in [4usize, 2, 1] {
            let counts = [fine[0] / s, fine[1] / s, 1];
            let (panels, faces) = surface(&bar, counts);
            let shape = trace(&bar, &panels, &faces, &fundamental.phi);
            let big_w = fundamental.lambda.sqrt();
            let (omega, z, min_eig) = window_pole(&panels, &shape, big_w, medium);
            let decay = -omega.im;
            let decrement = TAU * decay / omega.re;
            println!(
                "{{\"suite\":\"fs-couple\",\"case\":\"loaded-free-bar-radiation-decay\",\"inputs\":\"NASA TM X-2852 beam 40; matdb 2024-T3 TN D-6448; Medium::air\",\"axis\":\"bem-panels\",\"rung\":{counts:?},\"panels\":{},\"f_vacuo_hz\":{f_vacuo},\"f_measured_vacuum_hz\":{measured_vacuum},\"f_loaded_hz\":{},\"decay_1_s\":{decay:e},\"log_decrement\":{decrement:e},\"measured_air_minus_vacuum_decrement\":{delta_air:e},\"z\":[{:e},{:e}],\"min_eig_re\":{min_eig:e},\"elapsed_s\":{}}}",
                panels.areas().len(),
                omega.re / TAU,
                z.re,
                z.im,
                started.elapsed().as_secs_f64()
            );
            rows.push((panels.areas().len(), decay, omega.re));
        }
        // Successive changes of the decay rate (resistance) and of the loaded
        // frequency, coarse to fine.
        let d = [rows[0].1, rows[1].1, rows[2].1];
        let w = [rows[0].2, rows[1].2, rows[2].2];
        let (c1, c2) = ((d[1] - d[0]).abs() / d[2], (d[2] - d[1]).abs() / d[2]);
        let observed = (c1 / c2).log2();
        let finest_cents = cents(w[2], w[1]).abs();
        let decrement = TAU * d[2] / w[2];
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"loaded-free-bar-radiation-decay-orders\",\"decay_changes\":[{c1:e},{c2:e}],\"observed_order\":{observed},\"finest_frequency_change_cents\":{finest_cents:e},\"radiation_decrement\":{decrement:e},\"measured_air_minus_vacuum_decrement\":{delta_air:e},\"ratio\":{},\"tol\":{{\"decay\":0.02,\"cents\":0.1}}}}",
            decrement / delta_air
        );
        assert!(d.iter().all(|&v| v >= 0.0), "radiation damping >= 0");
        assert!(c2 < c1, "decay ladder must be live: {c1} then {c2}");
        assert!(c2 <= 0.02, "finest decay change {c2} > 2%");
        assert!(
            finest_cents <= 0.1,
            "finest loaded frequency change {finest_cents} cents"
        );
        // One-sided: radiation alone cannot exceed the whole measured air loss.
        assert!(
            decrement <= delta_air,
            "radiation {decrement} > measured air loss {delta_air}"
        );
    });
}
