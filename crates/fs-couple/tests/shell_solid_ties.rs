//! M1.2-W3 (bead pexfp.52): shell-solid face ties (`fs_plate::shell::junction::tie_solids_to_shells`).
//! A TET10 body is glued to a shell skin: every node of its glued P2 faces follows the shell through the
//! skin, u = sum_j lambda_j (u_j + theta_j x d_j), and the tied solid coordinates are eliminated by congruence.
//!
//! Fixture: a free 120 x 120 x 3 mm MITC3+ plate (isotropic E 10 GPa, nu 0.3, 450 kg/m3, 8 x 8 cells), with a
//! 30 x 30 x 20 mm TET10 maple block (E 12 GPa, nu 0.3, 600 kg/m3, 2 x 2 x 2 structured cells) standing on
//! the plate's upper skin over its middle, its bottom face glued.
//! - `g0_shell_solid_tie_virtual_work_identity`: over 1000 fs-rand fields, the projected load times the
//!   combined displacement equals the full loads times the full displacements (1e-13 of the sum of the
//!   magnitudes), and the combined strain energy equals the shell's plus the solid's (1e-12 relative);
//! - `g0_glued_block_keeps_six_rigid_modes`: plate and block together have exactly the six rigid modes of
//!   one free body (a missing tie would leave the block's six, an over-constraint fewer);
//! - the refusals: a block lifted off the skin, a block reaching past the master patch, a support-fixed
//!   glued node, and the tie budgets at cap and cap+1.
use fs_exec::Cx;
use fs_modal::{SliceOptions, slice_window};
use fs_plate::ShellMesh;
use fs_plate::shell::ShellModel;
use fs_plate::shell::junction::{
    ShellComponent, SolidFaceTie, SolidFaceTrace, SolidOnShells, TieBudget, TiedShellModel,
    tie_shells, tie_solids_to_shells,
};
use fs_plate::shell::mitc::{
    MitcBudget, NodalDirectors, ShearShellSection, ShellMass, assemble_mitc3p,
};
use fs_rand::StreamKey;
use fs_solid::linear3::tet10::{QuadraticTetProblem, TET10_EDGES, TetQuadrature, promote_p1_to_p2};
use fs_solid::linear3::{TetAssemblyBudget, TetElasticMaterial, TetMaterialField};
use fs_sparse::Csr;

const SIDE: f64 = 0.12;
const H: f64 = 0.003;
const CELLS: usize = 8;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = fs_exec::CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 52,
                kernel_id: 52,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn plate_mesh() -> ShellMesh {
    let id = |i: usize, j: usize| j * (CELLS + 1) + i;
    let mut nodes = Vec::new();
    for j in 0..=CELLS {
        for i in 0..=CELLS {
            nodes.push([
                SIDE * i as f64 / CELLS as f64 - 0.5 * SIDE,
                SIDE * j as f64 / CELLS as f64 - 0.5 * SIDE,
                0.0,
            ]);
        }
    }
    let mut tris = Vec::new();
    for j in 0..CELLS {
        for i in 0..CELLS {
            let (a, b, c, d) = (id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1));
            if (i + j) % 2 == 0 {
                tris.push([a, b, c]);
                tris.push([a, c, d]);
            } else {
                tris.push([a, b, d]);
                tris.push([b, c, d]);
            }
        }
    }
    ShellMesh::new(nodes, tris).unwrap()
}

fn mitc_plate(mesh: &ShellMesh) -> ShellModel {
    let g = 10e9 / 2.6;
    let section = ShearShellSection::from_orthotropic_solid(
        10e9,
        10e9,
        0.3,
        g,
        g,
        g,
        H,
        450.0,
        0.0,
        5.0 / 6.0,
    )
    .unwrap();
    let directors =
        NodalDirectors::new(mesh, vec![[0.0, 0.0, 1.0]; mesh.nodes.len()], 0.9).unwrap();
    assemble_mitc3p(
        mesh,
        &directors,
        &vec![section; mesh.tris.len()],
        &[],
        ShellMass::Lumped,
        MitcBudget {
            max_nodes: mesh.nodes.len(),
            max_triangles: mesh.tris.len(),
        },
        1e-3,
        || false,
    )
    .unwrap()
    .pencil()
    .clone()
}

struct Block {
    nodes: Vec<[f64; 3]>,
    k: Csr,
    m: Csr,
    dof_index: Vec<[Option<usize>; 3]>,
    /// P2 faces on the bottom plane z = z0.
    bottom: Vec<[usize; 6]>,
}

/// A TET10 maple box [x0, x1] x [y0, y1] x [z0, z0 + 0.02], free except for `fixed` full DOFs.
fn block(x: [f64; 2], y: [f64; 2], z0: f64, fixed: &[usize]) -> Block {
    with_cx(|cx| {
        let (complex, p1) = fs_conduction::fixtures::structured_tets([2, 2, 2], &|i, j, k| {
            [
                x[0] + (x[1] - x[0]) * i as f64 / 2.0,
                y[0] + (y[1] - y[0]) * j as f64 / 2.0,
                z0 + 0.02 * k as f64 / 2.0,
            ]
        });
        let tets: Vec<[usize; 4]> = complex.tets.iter().map(|t| t.map(|i| i as usize)).collect();
        let (nodes, tets10) = promote_p1_to_p2(&p1, &tets, None, cx).unwrap();
        let material = TetElasticMaterial::try_new(600.0, 12e9, 0.3, [52; 32]).unwrap();
        let assembly = QuadraticTetProblem {
            nodes_m: &nodes,
            elements: &tets10,
            materials: TetMaterialField::Uniform(&material),
            fixed_dofs: fixed,
            budget: TetAssemblyBudget::standard(),
            quadrature: TetQuadrature {
                stiffness_points_per_dir: 3,
                mass_points_per_dir: 4,
            },
        }
        .assemble(cx)
        .unwrap();
        let mut dof_index = vec![[None; 3]; nodes.len()];
        for (r, &d) in assembly.free_dofs.iter().enumerate() {
            dof_index[d / 3][d % 3] = Some(r);
        }
        // Boundary faces on z = z0: tetrahedron faces whose three corners lie on the plane, with the
        // midside nodes of their edges (TET10 node 4 + e for edge e of TET10_EDGES).
        let on = |n: usize| (nodes[n][2] - z0).abs() < 1e-12;
        let midside = |t: &[usize; 10], a: usize, b: usize| {
            let e = TET10_EDGES
                .iter()
                .position(|&[p, q]| (p, q) == (a.min(b), a.max(b)))
                .unwrap();
            t[4 + e]
        };
        let mut bottom = Vec::new();
        for t in &tets10 {
            for [a, b, c] in [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]] {
                if on(t[a]) && on(t[b]) && on(t[c]) {
                    bottom.push([
                        t[a],
                        t[b],
                        t[c],
                        midside(t, a, b),
                        midside(t, b, c),
                        midside(t, c, a),
                    ]);
                }
            }
        }
        Block {
            nodes,
            k: assembly.stiffness,
            m: assembly.mass,
            dof_index,
            bottom,
        }
    })
}

fn budget() -> TieBudget {
    TieBudget {
        max_segments: 1_000_000,
        max_slave_dofs: 1_000_000,
        gap_tolerance_m: 1e-6,
        angle_tolerance_rad: 0.1,
    }
}

/// The block glued over the master triangles whose centroid lies inside [-r, r]^2.
fn glue(
    mesh: &ShellMesh,
    shell: &ShellModel,
    b: &Block,
    patch_half_width: f64,
    budget: TieBudget,
) -> (
    Result<fs_plate::shell::junction::SolidTiedShellModel, fs_plate::PlateError>,
    TiedShellModel,
) {
    let thickness = vec![H; mesh.nodes.len()];
    let components = [ShellComponent {
        id: 0,
        mesh,
        model: shell,
        nodal_thickness_m: &thickness,
    }];
    let tied = tie_shells(&components, &[], budget, || false).unwrap();
    let master_tris: Vec<usize> = (0..mesh.tris.len())
        .filter(|&t| {
            let c = mesh.tris[t].iter().fold([0.0; 3], |a, &n| {
                [
                    a[0] + mesh.nodes[n][0] / 3.0,
                    a[1] + mesh.nodes[n][1] / 3.0,
                    0.0,
                ]
            });
            c[0].abs() < patch_half_width && c[1].abs() < patch_half_width
        })
        .collect();
    let solid = SolidOnShells {
        trace: SolidFaceTrace {
            positions_m: &b.nodes,
            dof_index: &b.dof_index,
        },
        k: &b.k,
        m: &b.m,
        ties: vec![SolidFaceTie {
            faces: b.bottom.clone(),
            master_component: 0,
            master_tris,
            side: 1,
        }],
    };
    let result = tie_solids_to_shells(&tied.model, &tied, &components, &[solid], budget, || false);
    (result, tied)
}

fn quadratic(a: &Csr, x: &[f64]) -> f64 {
    let mut y = vec![0.0; x.len()];
    a.spmv(x, &mut y);
    x.iter().zip(&y).map(|(p, q)| p * q).sum()
}

const BLOCK_X: [f64; 2] = [-0.015, 0.015];

#[test]
fn g0_shell_solid_tie_virtual_work_identity() {
    let mesh = plate_mesh();
    let shell = mitc_plate(&mesh);
    let b = block(BLOCK_X, BLOCK_X, 0.5 * H, &[]);
    let (glued, tied) = glue(&mesh, &shell, &b, 0.03, budget());
    let glued = glued.unwrap();
    let (n, n_shell, n_solid) = (glued.model.free, glued.shell_free, b.k.nrows());
    println!(
        "combined {n} coordinates: {n_shell} shell, {} solid kept, {} solid eliminated, {} projections, gap {:.1e} m",
        n - n_shell,
        glued.receipt.slave_dofs,
        glued.receipt.segments,
        glued.receipt.max_gap_m
    );
    let mut stream = StreamKey {
        seed: 52,
        kernel: 52,
        tile: 0,
    }
    .stream();
    let (mut worst_work, mut worst_energy) = (0.0_f64, 0.0_f64);
    for _ in 0..1000 {
        let q: Vec<f64> = (0..n).map(|_| stream.next_normal()).collect();
        let f_shell: Vec<f64> = (0..n_shell).map(|_| stream.next_normal()).collect();
        let f_solid: Vec<f64> = (0..n_solid).map(|_| stream.next_normal()).collect();
        let u_solid = glued.expand_solid(0, &q).unwrap();
        // Virtual work: (f_shell + T^T f_solid) . q == f_shell . q_shell + f_solid . u_solid.
        let mut load = glued.project_solid_load(0, &f_solid).unwrap();
        for (l, f) in load.iter_mut().zip(&f_shell) {
            *l += f;
        }
        let reduced: f64 = load.iter().zip(&q).map(|(a, b)| a * b).sum();
        let terms: Vec<f64> = f_shell
            .iter()
            .zip(&q)
            .map(|(a, b)| a * b)
            .chain(f_solid.iter().zip(&u_solid).map(|(a, b)| a * b))
            .collect();
        let full: f64 = terms.iter().sum();
        let scale: f64 = terms.iter().map(|t| t.abs()).sum();
        worst_work = worst_work.max((reduced - full).abs() / scale);
        // Energy: q^T K q == q_shell^T K_shell q_shell + u_solid^T K_solid u_solid (and the same for M).
        for (combined, shell_part, solid_part) in [
            (&glued.model.k, &tied.model.k, &b.k),
            (&glued.model.m, &tied.model.m, &b.m),
        ] {
            let e = quadratic(combined, &q);
            let parts = quadratic(shell_part, &q[..n_shell]) + quadratic(solid_part, &u_solid);
            worst_energy = worst_energy.max((e - parts).abs() / parts.abs());
        }
    }
    println!(
        "virtual work worst {worst_work:.2e} (of the sum of magnitudes), energy worst {worst_energy:.2e} (relative)"
    );
    assert!(worst_work < 1e-13, "virtual work identity {worst_work:.2e}");
    assert!(worst_energy < 1e-12, "energy identity {worst_energy:.2e}");
}

#[test]
fn g0_glued_block_keeps_six_rigid_modes() {
    let mesh = plate_mesh();
    let shell = mitc_plate(&mesh);
    let b = block(BLOCK_X, BLOCK_X, 0.5 * H, &[]);
    let glued = glue(&mesh, &shell, &b, 0.03, budget()).0.unwrap();
    // Shift so the rigid modes sit at +sigma; count the pencil's eigenvalues in (0, 2 sigma).
    let sigma = 1.0;
    let n = glued.model.free;
    let mut coo = fs_sparse::Coo::new(n, n);
    for r in 0..n {
        let (kc, kv) = glued.model.k.row(r);
        for (&c, &v) in kc.iter().zip(kv) {
            coo.push(r, c, v);
        }
        let (mc, mv) = glued.model.m.row(r);
        for (&c, &v) in mc.iter().zip(mv) {
            coo.push(r, c, sigma * v);
        }
    }
    let shifted = coo.assemble();
    let report = slice_window(
        &shifted,
        &glued.model.m,
        (0.5 * sigma, 1.5 * sigma),
        &SliceOptions::default(),
    )
    .unwrap();
    println!("eigenvalues of K + M in (0.5, 1.5): {}", report.expected);
    assert_eq!(report.expected, 6, "one free body has six rigid modes");
}

#[test]
fn g0_uncovered_or_lifted_solid_face_refuses() {
    let mesh = plate_mesh();
    let shell = mitc_plate(&mesh);
    let on_skin = block(BLOCK_X, BLOCK_X, 0.5 * H, &[]);
    assert!(glue(&mesh, &shell, &on_skin, 0.03, budget()).0.is_ok());
    // 1 mm above the skin.
    let lifted = block(BLOCK_X, BLOCK_X, 0.5 * H + 0.001, &[]);
    assert!(glue(&mesh, &shell, &lifted, 0.03, budget()).0.is_err());
    // The master patch covers only [-15, 15] mm cell centres... reach past it with a block from -15 to +45 mm.
    let past = block([-0.015, 0.045], BLOCK_X, 0.5 * H, &[]);
    assert!(glue(&mesh, &shell, &past, 0.016, budget()).0.is_err());
}

#[test]
fn g0_support_fixed_glued_node_refuses() {
    let mesh = plate_mesh();
    let shell = mitc_plate(&mesh);
    let free = block(BLOCK_X, BLOCK_X, 0.5 * H, &[]);
    let node = free.bottom[0][0];
    let fixed = block(BLOCK_X, BLOCK_X, 0.5 * H, &[3 * node + 2]);
    assert!(glue(&mesh, &shell, &fixed, 0.03, budget()).0.is_err());
}

#[test]
fn g0_solid_tie_budget_cap_and_cap_plus_one() {
    let mesh = plate_mesh();
    let shell = mitc_plate(&mesh);
    let b = block(BLOCK_X, BLOCK_X, 0.5 * H, &[]);
    let full = glue(&mesh, &shell, &b, 0.03, budget()).0.unwrap();
    let (segments, dofs) = (full.receipt.segments, full.receipt.slave_dofs);
    assert_eq!(dofs, 3 * segments);
    let at = |max_segments: usize, max_slave_dofs: usize| {
        glue(
            &mesh,
            &shell,
            &b,
            0.03,
            TieBudget {
                max_segments,
                max_slave_dofs,
                ..budget()
            },
        )
        .0
    };
    assert!(at(segments, dofs).is_ok());
    assert!(at(segments - 1, dofs).is_err());
    assert!(at(segments, dofs - 1).is_err());
}
