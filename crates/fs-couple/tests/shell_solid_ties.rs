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
//!
//! Soundpost fixture (`g2_shell_solid_tie_soundpost_fixture`, orchestrator ruling mail 6174): two free
//! quarter-sawn Sitka spruce plates, 120 x 80 x 3 mm below and 100 x 70 x 3 mm above (L along x, R along
//! y, T through the thickness), joined at their centres by a spruce post (L along its axis), 51 mm long
//! between the skins. The tensor is the FPL-GTR-282 seed pack through the M0.2-W1 adapter; maple has no
//! 3-D pack. The post is SQUARE with the declared round post's area, a = r sqrt(pi) (r = 3.175 mm), so its
//! axial stiffness is the round post's exactly, while its bending second moment a^4/12 is pi/3 = 1.047
//! times the round post's pi r^4/4 (+4.7 %). The fixture verifies a shape-agnostic tie, not instrument
//! geometry: the instrument keeps the round post, and the round-post mesher is an open gap of .52.
//! Nothing gated depends on post bending: only the doubly symmetric class (mirror-even about both vertical
//! planes through the post axis) is compared, and in it the post cannot tilt, translate sideways or twist.
//! - Reference: the all-3-D TET10 quarter of that class (both plates and the post as one conforming Kuhn
//!   mesh, graded towards the post) on 4, 6 and 8 element layers per plate. A Richardson bound above
//!   TOL/3 refuses (an unconverged reference).
//! - Tied: the same quarter as MITC3+ grids and a TET10 post glued by `tie_solids_to_shells`, with the
//!   class's symmetry conditions imposed on the combined coordinates.
//! - Gates, 1 %: the first ten symmetric elastic modes, MAC-paired (>= 0.9) on mid-surface w; and the
//!   static system stiffness under a rim squeeze (F = 1 N per plate, up on the lower rim and down on the
//!   upper, uniform per unit rim length; K = F^2 / work).
//! - Printed, not gated: the post-path stiffness, F over the approach of the two mid-surfaces on the post
//!   axis. It is a local quantity: the plates bend across the post's foot, so the axis approach is less
//!   than the post's own shortening, and both models exceed the bare post's E_L a^2 / L. The shell also
//!   has no through-thickness compliance (spruce E_T is 4.3 % of E_L). Measured 2026-10-23: the tied
//!   model is 12 % below the 3-D one. A consumer that needs this local stiffness must not take it from
//!   the shell tie.
use fs_evidence::ValidityDomain;
use fs_exec::Cx;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::OrthotropicElastic;
use fs_material::state_point::*;
use fs_modal::{ModePair, SliceOptions, slice_window};
use fs_plate::ShellMesh;
use fs_plate::shell::ShellModel;
use fs_plate::shell::junction::{
    ShellComponent, SolidFaceTie, SolidFaceTrace, SolidOnShells, TieBudget, TieReceipt,
    TiedShellModel, tie_shells, tie_solids_to_shells,
};
use fs_plate::shell::mitc::{
    MitcBudget, NodalDirectors, ShearShellSection, ShellMass, assemble_mitc3p,
};
use fs_plate::shell::survey::material_axis_angle;
use fs_qty::{Density, Dims, Pressure};
use fs_rand::StreamKey;
use fs_solid::linear3::tet10::{QuadraticTetProblem, TET10_EDGES, TetQuadrature, promote_p1_to_p2};
use fs_solid::linear3::{TetAssemblyBudget, TetElasticMaterial, TetMaterialField};
use fs_sparse::{Csr, DirectOrdering, LdltOptions, NdOptions, SymbolicLdlt};
use std::collections::BTreeMap;
use std::time::Instant;

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

// ---------------------------------------------------------------------------
// The soundpost fixture (g2_shell_solid_tie_soundpost_fixture)
// ---------------------------------------------------------------------------

const SPRUCE: &str =
    include_str!("../../../data/matdb/seed-v1/spruce-sitka-fpl-gtr282/properties.tsv");
/// The declared round soundpost (bead pexfp.52 readiness): radius 3.175 mm (6.35 mm diameter) and
/// 51 mm between the plates' inner skins.
const POST_R: f64 = 3.175e-3;
const POST_L: f64 = 0.051;
/// Plate half-extents [x along L, y along R]: the lower plate and a smaller upper plate, so that no
/// mode of one plate is degenerate with a mode of the other.
const LOWER: [f64; 2] = [0.060, 0.040];
const UPPER: [f64; 2] = [0.050, 0.035];
/// Base rung, per quarter axis: the post's half-width is one cell, then `OUTER_CELLS` cells grow
/// geometrically from it to the plate edge. `LAYERS` element layers per plate and `POST_CELLS` cells
/// along the post. A rung r splits every base cell r times in every direction.
const OUTER_CELLS: usize = 4;
const LAYERS: usize = 2;
const POST_CELLS: usize = 4;
/// 3-D rungs: 4, 6 and 8 element layers through each plate. The tied model's post is meshed as the
/// finest 3-D rung's post; its shells take the shell rungs.
const RUNGS_3D: [usize; 3] = [2, 3, 4];
const RUNGS_SHELL: [usize; 3] = [2, 4, 8];
/// The tied model's solved band [Hz]: wide enough for eleven doubly symmetric modes.
const TIED_BAND: (f64, f64) = (1.0, 8000.0);
const KAPPA: f64 = 5.0 / 6.0;
/// The spruce pack's only validity state.
const MC_PCT: f64 = 12.0;
const TOL: f64 = 0.01;
const MAC_PAIR: f64 = 0.9;
/// The gated modes: the first ten elastic modes of the doubly symmetric class.
const GATED: usize = 10;
const TAU: f64 = core::f64::consts::TAU;

/// The seed pack's scalar rows as a material card (as in `shell_solid_crosscheck.rs`, pexfp.49).
fn pack_card(tsv: &str, label: &str) -> MaterialCard {
    let mut validity: BTreeMap<&str, (&str, f64, f64)> = BTreeMap::new();
    for line in tsv.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.first() == Some(&"validity") {
            validity.insert(f[1], (f[2], f[3].parse().unwrap(), f[4].parse().unwrap()));
        }
    }
    let mut claims = ClaimSet::new();
    for line in tsv.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.first() != Some(&"scalar") {
            continue;
        }
        let raw: f64 = f[4].parse().expect("source scalar");
        let (dims, value) = match f[5] {
            "1" => (Dims::NONE, raw),
            "MPa" => (Pressure::DIMS, raw * 1e6),
            "GPa" => (Pressure::DIMS, raw * 1e9),
            "kg/m3" => (Density::DIMS, raw),
            unit => panic!("unexpected source unit {unit} on {}", f[3]),
        };
        let mut domain = ValidityDomain::unconstrained();
        if let Some(&(axis, lo, hi)) = validity.get(f[1]) {
            domain = domain.with(axis, lo, hi);
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(f[3], dims),
                value: PropertyValue::Scalar { value, dims },
                validity: domain,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: format!("{label} seed-v1 row {}", f[1]),
                    license: "Work-of-the-US-Government-Public-Use-Permitted".into(),
                    artifact: None,
                },
            })
            .expect("pack claim");
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: label.into(),
            phase: "clear wood".into(),
            process: "seed-v1 pack rows".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .expect("pack card")
}

struct Spruce {
    state: OrthotropicElasticStatePoint,
    law: OrthotropicElastic,
    rho: f64,
}

/// The full 3-D Sitka tensor through the M0.2-W1 reference-ratio adapter; no modulus is authored here.
fn spruce() -> Spruce {
    let card = pack_card(SPRUCE, "spruce-sitka-fpl-gtr282");
    let resolved = resolve_orthotropic_elastic_state_point_from_reference_ratios(
        &card,
        &QueryPoint::new().with("moisture_content", MC_PCT).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        1e-3, // declared small-strain domain, not a measured failure strain
        OrthotropicRequirementSet::Full3d,
        PoissonReciprocityRule::MajorRatios,
    )
    .expect("full 3-D tensor from the seed pack");
    let rho = resolved.density_kg_m3();
    let law = resolved.law().expect("full 3-D law").clone();
    let ReferenceRatioOrthotropicStatePoint::Full3d(state) = resolved else {
        panic!("a Full3d request returned a plane state");
    };
    Spruce { state, law, rho }
}

/// The square post of the round post's area: a = r sqrt(pi).
fn post_side() -> f64 {
    POST_R * core::f64::consts::PI.sqrt()
}

/// The lower plate's lower skin, its upper skin (the post's foot), the post's head (the upper plate's
/// lower skin) and the upper plate's upper skin.
fn z_levels() -> [f64; 4] {
    let foot = 0.5 * H;
    let head = foot + POST_L;
    [-0.5 * H, foot, head, head + H]
}

/// `n` equal cells from `a` to `b`, with exact endpoints.
fn cells(a: f64, b: f64, n: usize) -> Vec<f64> {
    (0..=n)
        .map(|i| {
            if i == n {
                b
            } else {
                a + (b - a) * i as f64 / n as f64
            }
        })
        .collect()
}

/// Every cell of `base` split into `r` equal cells.
fn refine(base: &[f64], r: usize) -> Vec<f64> {
    let mut out = vec![base[0]];
    for w in base.windows(2) {
        out.extend_from_slice(&cells(w[0], w[1], r)[1..]);
    }
    out
}

/// One quarter axis from the post axis to `half` at rung `r`: the post's half-width a/2 is one base
/// cell, then `OUTER_CELLS` base cells of sizes (a/2) q^k, with q fitted so they end at `half`. The
/// post's coordinates are computed the same way on every axis, so the blocks share them bitwise.
fn quarter_axis(half: f64, r: usize) -> Vec<f64> {
    let s0 = 0.5 * post_side();
    let span = (half - s0) / s0;
    let sum = |q: f64| {
        let (mut total, mut power) = (0.0, 1.0);
        for _ in 0..OUTER_CELLS {
            power *= q;
            total += power;
        }
        total
    };
    let (mut lo, mut hi) = (1.0, 4.0);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if sum(mid) < span {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let q = 0.5 * (lo + hi);
    let mut base = vec![0.0, s0];
    let mut power = 1.0;
    for k in 0..OUTER_CELLS {
        power *= q;
        let next = if k + 1 == OUTER_CELLS {
            half
        } else {
            base[base.len() - 1] + s0 * power
        };
        base.push(next);
    }
    refine(&base, r)
}

/// The base-grid points of a quarter axis pair at rung `r` (the MAC sample grid).
fn base_points(xs: &[f64], ys: &[f64], r: usize) -> Vec<[f64; 2]> {
    let mut out = Vec::new();
    for i in (0..xs.len()).step_by(r) {
        for j in (0..ys.len()).step_by(r) {
            out.push([xs[i], ys[j]]);
        }
    }
    out
}

/// P1 vertices merged by exact coordinates (-0 and +0 are one vertex).
#[derive(Default)]
struct Vertices {
    nodes: Vec<[f64; 3]>,
    index: BTreeMap<[u64; 3], usize>,
}

impl Vertices {
    fn key(p: [f64; 3]) -> [u64; 3] {
        p.map(|c| (c + 0.0).to_bits())
    }

    fn at(&mut self, p: [f64; 3]) -> usize {
        let next = self.nodes.len();
        let id = *self.index.entry(Self::key(p)).or_insert(next);
        if id == next {
            self.nodes.push(p);
        }
        id
    }

    fn find(&self, p: [f64; 3]) -> usize {
        self.index[&Self::key(p)]
    }
}

/// Kuhn split of every cell of the tensor grid `xs x ys x zs`: six tetrahedra about the cell's main
/// diagonal. The split is translation invariant, so two blocks that share a grid face share its
/// triangles and the merged mesh is conforming.
fn kuhn_block(
    xs: &[f64],
    ys: &[f64],
    zs: &[f64],
    vertices: &mut Vertices,
    tets: &mut Vec<[usize; 4]>,
) {
    const PATHS: [[usize; 3]; 6] = [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    for i in 0..xs.len() - 1 {
        for j in 0..ys.len() - 1 {
            for k in 0..zs.len() - 1 {
                let corner = |d: [usize; 3]| [xs[i + d[0]], ys[j + d[1]], zs[k + d[2]]];
                for path in PATHS {
                    let mut d = [0; 3];
                    let mut t = [vertices.at(corner(d)); 4];
                    for (step, &axis) in path.iter().enumerate() {
                        d[axis] = 1;
                        t[step + 1] = vertices.at(corner(d));
                    }
                    tets.push(t);
                }
            }
        }
    }
}

/// P2 faces (corners, then the midsides of edges 01, 12, 20) of the tetrahedron faces whose three
/// corners satisfy `on`.
fn p2_faces(
    nodes: &[[f64; 3]],
    tets10: &[[usize; 10]],
    on: impl Fn([f64; 3]) -> bool,
) -> Vec<[usize; 6]> {
    let midside = |t: &[usize; 10], a: usize, b: usize| {
        let e = TET10_EDGES
            .iter()
            .position(|&[p, q]| (p, q) == (a.min(b), a.max(b)))
            .unwrap();
        t[4 + e]
    };
    let mut out = Vec::new();
    for t in tets10 {
        for [a, b, c] in [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]] {
            if on(nodes[t[a]]) && on(nodes[t[b]]) && on(nodes[t[c]]) {
                out.push([
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
    out
}

fn triangle_area(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let (u, w) = (
        [b[0] - a[0], b[1] - a[1], b[2] - a[2]],
        [c[0] - a[0], c[1] - a[1], c[2] - a[2]],
    );
    let n = [
        u[1] * w[2] - u[2] * w[1],
        u[2] * w[0] - u[0] * w[2],
        u[0] * w[1] - u[1] * w[0],
    ];
    0.5 * (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt()
}

/// Solve K u = f with the coordinates `pinned` held at zero (the rigid modes removed, for a
/// self-equilibrated load); `u` is zero there.
fn solve_pinned(
    k: &Csr,
    f: &[f64],
    pinned: &[usize],
    coordinates: Option<&[[f64; 3]]>,
) -> Vec<f64> {
    let n = k.nrows();
    let mut row_of = vec![usize::MAX; n];
    let mut kept = Vec::with_capacity(n);
    for r in 0..n {
        if !pinned.contains(&r) {
            row_of[r] = kept.len();
            kept.push(r);
        }
    }
    let mut coo = fs_sparse::Coo::new(kept.len(), kept.len());
    for (i, &r) in kept.iter().enumerate() {
        let (cols, vals) = k.row(r);
        for (&c, &v) in cols.iter().zip(vals) {
            if row_of[c] != usize::MAX {
                coo.push(i, row_of[c], v);
            }
        }
    }
    let reduced = coo.assemble();
    let coordinates: Option<Vec<[f64; 3]>> =
        coordinates.map(|c| kept.iter().map(|&r| c[r]).collect());
    let ordering = if coordinates.is_some() {
        DirectOrdering::NestedDissection(NdOptions {
            leaf_size: 64,
            refinement_passes: 4,
        })
    } else {
        DirectOrdering::Amd
    };
    let factor = SymbolicLdlt::analyze_with_coordinates(&reduced, ordering, coordinates.as_deref())
        .expect("ordering")
        .factor(&reduced, &LdltOptions::default())
        .expect("LDL^T of the pinned stiffness");
    let x = factor.solve(&kept.iter().map(|&r| f[r]).collect::<Vec<_>>());
    let mut u = vec![0.0; n];
    for (i, &r) in kept.iter().enumerate() {
        u[r] = x[i];
    }
    u
}

fn options() -> SliceOptions {
    SliceOptions {
        max_lanczos: 2048,
        max_restarts: 48,
        ritz_tol: 1e-13,
        mass_diagonal_equilibration: true,
        ..SliceOptions::default()
    }
}

fn window(band: (f64, f64)) -> (f64, f64) {
    ((TAU * band.0).powi(2), (TAU * band.1).powi(2))
}

fn hz(m: &ModePair) -> f64 {
    m.lambda.max(0.0).sqrt() / TAU
}

/// One doubly symmetric mode as the comparison sees it: its frequency and the mid-surface w of both
/// plates (lower, then upper) at the quarter's base grid.
#[derive(Clone, Debug)]
struct Mode {
    hz: f64,
    samples: Vec<f64>,
}

#[derive(Debug, PartialEq, Eq)]
enum PairRefusal {
    Unpaired,
    Ambiguous(Vec<usize>),
}

fn mac(a: &[f64], b: &[f64]) -> f64 {
    let dot = |u: &[f64], v: &[f64]| u.iter().zip(v).map(|(x, y)| x * y).sum::<f64>();
    let ab = dot(a, b);
    ab * ab / (dot(a, a) * dot(b, b))
}

/// The unique candidate whose MAC with `target` reaches [`MAC_PAIR`].
fn pair(target: &[f64], candidates: &[&[f64]]) -> Result<usize, PairRefusal> {
    let hits: Vec<usize> = candidates
        .iter()
        .enumerate()
        .filter(|(_, c)| mac(target, c) >= MAC_PAIR)
        .map(|(i, _)| i)
        .collect();
    match hits.as_slice() {
        [i] => Ok(*i),
        [] => Err(PairRefusal::Unpaired),
        _ => Err(PairRefusal::Ambiguous(hits)),
    }
}

/// Share of `v` in span{x, y}: `v^T P v / v^T v` for the orthogonal projector P.
fn span_fraction(v: &[f64], x: &[f64], y: &[f64]) -> f64 {
    let dot = |u: &[f64], w: &[f64]| u.iter().zip(w).map(|(a, b)| a * b).sum::<f64>();
    let (xx, xy, yy, xv, yv) = (dot(x, x), dot(x, y), dot(y, y), dot(x, v), dot(y, v));
    let det = xx * yy - xy * xy;
    if !(det > 0.0) {
        return 0.0;
    }
    let (a, b) = ((yy * xv - xy * yv) / det, (xx * yv - xy * xv) / det);
    (a * xv + b * yv) / dot(v, v)
}

/// Frequency of `target`'s partner among `other` modes (another rung, or the other model), as in
/// `shell_solid_crosscheck.rs`: a unique MAC match; or, when a near-degenerate pair is rotated into a
/// different basis, the target and exactly one other mode both lie at least [`MAC_PAIR`] in the span
/// of the same two `other` modes, and within that certified cluster they pair by frequency rank.
fn track(target: usize, fine: &[Mode], other: &[Mode]) -> Result<f64, PairRefusal> {
    let candidates: Vec<&[f64]> = other.iter().map(|c| c.samples.as_slice()).collect();
    match pair(&fine[target].samples, &candidates) {
        Ok(k) => return Ok(other[k].hz),
        Err(PairRefusal::Ambiguous(hits)) => return Err(PairRefusal::Ambiguous(hits)),
        Err(PairRefusal::Unpaired) => {}
    }
    let mut cluster = None;
    for a in 0..other.len() {
        for b in a + 1..other.len() {
            let (x, y) = (&other[a].samples, &other[b].samples);
            if span_fraction(&fine[target].samples, x, y) < MAC_PAIR {
                continue;
            }
            let partners: Vec<usize> = (0..fine.len())
                .filter(|&j| j != target && span_fraction(&fine[j].samples, x, y) >= MAC_PAIR)
                .collect();
            if let [j] = partners[..] {
                if cluster.is_some() {
                    return Err(PairRefusal::Ambiguous(vec![a, b]));
                }
                cluster = Some((a, b, j));
            }
        }
    }
    let (a, b, j) = cluster.ok_or(PairRefusal::Unpaired)?;
    let (lo, hi) = if other[a].hz <= other[b].hz {
        (other[a].hz, other[b].hz)
    } else {
        (other[b].hz, other[a].hz)
    };
    Ok(if fine[target].hz <= fine[j].hz {
        lo
    } else {
        hi
    })
}

/// Richardson on three rungs (as in `shell_solid_crosscheck.rs`): the order solves
/// `(h1^p - h2^p)/(h2^p - h3^p) = (f1 - f2)/(f2 - f3)` by bisection and is NaN for a non-monotone
/// ladder. Richardson applies only when the order is asymptotic (3..=5); otherwise the declared bound
/// is the last increment. Returns (order, reference, bound).
fn ladder(f: [f64; 3], factors: [usize; 3]) -> (f64, f64, f64) {
    let h = factors.map(|k| 1.0 / k as f64);
    let (d1, d2) = (f[0] - f[1], f[1] - f[2]);
    let pow = |x: f64, p: f64| fs_math::det::pow(x, p);
    let g = |p: f64| (pow(h[0], p) - pow(h[1], p)) / (pow(h[1], p) - pow(h[2], p));
    let ratio = d1 / d2;
    let order = if d1 > 0.0 && d2 > 0.0 && (g(0.05)..=g(12.0)).contains(&ratio) {
        let (mut lo, mut hi) = (0.05, 12.0);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if g(mid) < ratio {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    } else {
        f64::NAN
    };
    if (3.0..=5.0).contains(&order) {
        let tail = d2 * pow(h[2], order) / (pow(h[1], order) - pow(h[2], order));
        (order, f[2] - tail, tail.abs())
    } else {
        (order, f[2], d2.abs())
    }
}

/// The all-3-D quarter on one rung.
struct Reference3d {
    modes: Vec<Mode>,
    /// Eigenvalues below the window by inertia: the rigid translation along z.
    below: usize,
    /// Rim squeeze: F^2 / work for F = 1 N on each plate [N/m].
    k_system: f64,
    /// Rim squeeze: F over the approach of the two mid-surfaces on the post axis [N/m].
    k_post_path: f64,
    dofs: usize,
    seconds: f64,
}

/// The rim-squeeze load on each plate: F = 1 N in total, up on the lower plate and down on the upper,
/// uniform per unit rim length.
const SQUEEZE: [([f64; 2], f64); 2] = [(LOWER, 1.0), (UPPER, -1.0)];

/// The doubly symmetric quarter (x, y >= 0) of the all-3-D TET10 model: both plates and the post as
/// one conforming Kuhn mesh, u_x = 0 on x = 0 and u_y = 0 on y = 0.
fn reference_3d(
    r: usize,
    plate: &TetElasticMaterial,
    post: &TetElasticMaterial,
    band: (f64, f64),
    cx: &Cx<'_>,
) -> Reference3d {
    let started = Instant::now();
    let z = z_levels();
    let (xl, yl) = (quarter_axis(LOWER[0], r), quarter_axis(LOWER[1], r));
    let (xu, yu) = (quarter_axis(UPPER[0], r), quarter_axis(UPPER[1], r));
    let along_post = cells(z[1], z[2], POST_CELLS * r);
    let mut vertices = Vertices::default();
    let mut tets = Vec::new();
    kuhn_block(
        &xl,
        &yl,
        &cells(z[0], z[1], LAYERS * r),
        &mut vertices,
        &mut tets,
    );
    let post_start = tets.len();
    kuhn_block(&xl[..=r], &yl[..=r], &along_post, &mut vertices, &mut tets);
    let post_end = tets.len();
    kuhn_block(
        &xu,
        &yu,
        &cells(z[2], z[3], LAYERS * r),
        &mut vertices,
        &mut tets,
    );
    // Mid-surface w is the mean of the two skins'.
    let mut samples = Vec::new();
    for (xs, ys, lo, hi) in [(&xl, &yl, z[0], z[1]), (&xu, &yu, z[2], z[3])] {
        for [x, y] in base_points(xs, ys, r) {
            samples.push([vertices.find([x, y, lo]), vertices.find([x, y, hi])]);
        }
    }
    let on_axis = z.map(|zz| vertices.find([0.0, 0.0, zz]));
    let pin = vertices.find([0.0, 0.0, along_post[POST_CELLS * r / 2]]);

    let (nodes, tets10) = promote_p1_to_p2(&vertices.nodes, &tets, None, cx).expect("P2 mesh");
    let materials: Vec<TetElasticMaterial> = (0..tets10.len())
        .map(|t| {
            if (post_start..post_end).contains(&t) {
                post.clone()
            } else {
                plate.clone()
            }
        })
        .collect();
    let mut fixed = Vec::new();
    for (n, p) in nodes.iter().enumerate() {
        if p[0] == 0.0 {
            fixed.push(3 * n);
        }
        if p[1] == 0.0 {
            fixed.push(3 * n + 1);
        }
    }
    let assembly = QuadraticTetProblem {
        nodes_m: &nodes,
        elements: &tets10,
        materials: TetMaterialField::PerElement(&materials),
        fixed_dofs: &fixed,
        budget: TetAssemblyBudget::standard(),
        quadrature: TetQuadrature {
            stiffness_points_per_dir: 3,
            mass_points_per_dir: 4,
        },
    }
    .assemble(cx)
    .expect("TET10 assembly");
    let full = |q: &[f64]| {
        let mut u = vec![0.0; 3 * nodes.len()];
        for (row, &d) in assembly.free_dofs.iter().enumerate() {
            u[d] = q[row];
        }
        u
    };

    // Modes: nested dissection, as for the pexfp.49 TET10 plates.
    let mut opts = options();
    opts.ordering = DirectOrdering::NestedDissection(NdOptions {
        leaf_size: 64,
        refinement_passes: 4,
    });
    let report = slice_window(&assembly.stiffness, &assembly.mass, window(band), &opts)
        .expect("certified 3-D window");
    let modes = report
        .modes
        .iter()
        .map(|m| {
            let u = full(&m.phi);
            Mode {
                hz: hz(m),
                samples: samples
                    .iter()
                    .map(|&[a, b]| 0.5 * (u[3 * a + 2] + u[3 * b + 2]))
                    .collect(),
            }
        })
        .collect();

    // Rim squeeze: the z traction F / (rim length x H) on the quarter's two rim faces of each plate
    // (a quarter of the load); a P2 face's consistent load of a uniform traction is A/3 at each midside.
    let mut f_full = vec![0.0; 3 * nodes.len()];
    for ((half, sign), lo, hi) in [(SQUEEZE[0], z[0], z[1]), (SQUEEZE[1], z[2], z[3])] {
        let traction = sign / (4.0 * (half[0] + half[1]) * H);
        for axis in 0..2 {
            let on = |p: [f64; 3]| p[axis] == half[axis] && p[2] >= lo && p[2] <= hi;
            for face in p2_faces(&nodes, &tets10, on) {
                let area = triangle_area(nodes[face[0]], nodes[face[1]], nodes[face[2]]);
                for &m in &face[3..] {
                    f_full[3 * m + 2] += traction * area / 3.0;
                }
            }
        }
    }
    let f: Vec<f64> = assembly.free_dofs.iter().map(|&d| f_full[d]).collect();
    let pin_row = assembly
        .free_dofs
        .iter()
        .position(|&d| d == 3 * pin + 2)
        .unwrap();
    let coordinates: Vec<[f64; 3]> = assembly.free_dofs.iter().map(|&d| nodes[d / 3]).collect();
    let q = solve_pinned(&assembly.stiffness, &f, &[pin_row], Some(&coordinates));
    let work_quarter: f64 = f.iter().zip(&q).map(|(a, b)| a * b).sum();
    let u = full(&q);
    let w = |n: usize| u[3 * n + 2];
    let approach = 0.5 * (w(on_axis[0]) + w(on_axis[1])) - 0.5 * (w(on_axis[2]) + w(on_axis[3]));
    let seconds = started.elapsed().as_secs_f64();
    println!(
        "3-D rung {r}: {} free DOF, {} tetrahedra, {} modes in the window, {} below, nnz(L) {}, {seconds:.1} s",
        assembly.free_dofs.len(),
        tets10.len(),
        report.modes.len(),
        report.below_low,
        report.stats.factor_nnz_l
    );
    Reference3d {
        modes,
        below: report.below_low,
        k_system: 1.0 / (4.0 * work_quarter),
        k_post_path: 1.0 / approach,
        dofs: assembly.free_dofs.len(),
        seconds,
    }
}

/// Parity of a field under a mirror plane through the post axis (x = 0 or y = 0).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Parity {
    Even,
    Odd,
}

/// A symmetry class: the parity of w under x -> -x and under y -> -y.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Class {
    x: Parity,
    y: Parity,
}

const EVEN_EVEN: Class = Class {
    x: Parity::Even,
    y: Parity::Even,
};

impl Class {
    fn parity(self, axis: usize) -> Parity {
        if axis == 0 { self.x } else { self.y }
    }

    /// Translation components (0 = x, 1 = y, 2 = z) that vanish on the plane `axis` = 0.
    fn fixed_translations(self, axis: usize) -> &'static [usize] {
        match (axis, self.parity(axis)) {
            (0, Parity::Even) => &[0],
            (0, Parity::Odd) => &[1, 2],
            (_, Parity::Even) => &[1],
            (_, Parity::Odd) => &[0, 2],
        }
    }

    /// Rotation components that vanish on that plane. A rotation is a pseudovector: under x -> -x it
    /// maps (rx, ry, rz) to (rx, -ry, -rz), under y -> -y to (-rx, ry, -rz).
    fn fixed_rotations(self, axis: usize) -> &'static [usize] {
        match (axis, self.parity(axis)) {
            (0, Parity::Even) => &[1, 2],
            (0, Parity::Odd) => &[0],
            (_, Parity::Even) => &[0, 2],
            (_, Parity::Odd) => &[1],
        }
    }
}

/// A quarter grid (x, y >= 0) on the quarter axes `xs`, `ys` at height `z`; node (i, j) is
/// `j * xs.len() + i`.
fn shell_quarter(xs: &[f64], ys: &[f64], z: f64) -> ShellMesh {
    let id = |i: usize, j: usize| j * xs.len() + i;
    let mut nodes = Vec::with_capacity(xs.len() * ys.len());
    for &y in ys {
        for &x in xs {
            nodes.push([x, y, z]);
        }
    }
    let mut tris = Vec::with_capacity(2 * (xs.len() - 1) * (ys.len() - 1));
    for b in 0..ys.len() - 1 {
        for a in 0..xs.len() - 1 {
            let (p00, p10, p01, p11) = (id(a, b), id(a + 1, b), id(a, b + 1), id(a + 1, b + 1));
            tris.push([p00, p10, p11]);
            tris.push([p00, p11, p01]);
        }
    }
    ShellMesh::new(nodes, tris).expect("quarter shell grid")
}

/// MITC3+ spruce plate: L along x, consistent mass.
fn mitc_spruce(mesh: &ShellMesh, section: &ShearShellSection) -> ShellModel {
    let n = mesh.nodes.len();
    let directors = NodalDirectors::new(mesh, vec![[0.0, 0.0, 1.0]; n], 0.9).unwrap();
    let sections: Vec<ShearShellSection> = (0..mesh.tris.len())
        .map(|e| {
            let angle = material_axis_angle(&mesh.facet(e).unwrap(), [1.0, 0.0, 0.0]).unwrap();
            section.with_material_angle_rad(angle).unwrap()
        })
        .collect();
    assemble_mitc3p(
        mesh,
        &directors,
        &sections,
        &[],
        ShellMass::Consistent,
        MitcBudget {
            max_nodes: n,
            max_triangles: mesh.tris.len(),
        },
        1e-3,
        || false,
    )
    .expect("MITC3+ assembly")
    .pencil()
    .clone()
}

/// The shell coordinates of `tied` that `class` makes vanish: on x = 0 and on y = 0, the translations
/// and rotations of [`Class::fixed_translations`] and [`Class::fixed_rotations`].
fn shell_symmetry_coordinates(
    tied: &TiedShellModel,
    meshes: &[&ShellMesh],
    class: Class,
) -> Vec<usize> {
    let mut out = Vec::new();
    for (c, mesh) in meshes.iter().enumerate() {
        for (n, p) in mesh.nodes.iter().enumerate() {
            for axis in 0..2 {
                if p[axis] != 0.0 {
                    continue;
                }
                let dofs = class
                    .fixed_translations(axis)
                    .iter()
                    .copied()
                    .chain(class.fixed_rotations(axis).iter().map(|r| 3 + r));
                for d in dofs {
                    out.push(tied.map[c][n][d].expect("a free shell coordinate"));
                }
            }
        }
    }
    out
}

/// The pencil without the coordinates `removed`, and the original index of each kept coordinate
/// (ascending).
fn restrict(k: &Csr, m: &Csr, removed: &[usize]) -> (Csr, Csr, Vec<usize>) {
    let n = k.nrows();
    let mut row_of = vec![usize::MAX; n];
    let mut kept = Vec::with_capacity(n);
    for r in 0..n {
        if removed.binary_search(&r).is_err() {
            row_of[r] = kept.len();
            kept.push(r);
        }
    }
    let part = |a: &Csr| {
        let mut coo = fs_sparse::Coo::new(kept.len(), kept.len());
        for (i, &r) in kept.iter().enumerate() {
            let (cols, vals) = a.row(r);
            for (&c, &v) in cols.iter().zip(vals) {
                if row_of[c] != usize::MAX {
                    coo.push(i, row_of[c], v);
                }
            }
        }
        coo.assemble()
    };
    (part(k), part(m), kept)
}

/// The shell-plus-post model on one shell rung, as the doubly symmetric quarter.
struct Tied {
    modes: Vec<Mode>,
    k_system: f64,
    k_post_path: f64,
    dofs: usize,
    receipt: TieReceipt,
    seconds: f64,
}

/// The doubly symmetric quarter of the tied model: two MITC3+ quarter grids at shell rung `r_s`, and
/// the quarter post (Kuhn cells, post rung `r_p`) glued to the lower plate's upper skin and the upper
/// plate's lower skin by `tie_solids_to_shells`. The class's conditions on x = 0 and y = 0 then remove
/// the vanishing shell coordinates and the post's untied ones; the post's glued nodes on those planes
/// follow their master nodes, which the same conditions hold.
fn tied_model(
    r_s: usize,
    r_p: usize,
    section: &ShearShellSection,
    post: &TetElasticMaterial,
    cx: &Cx<'_>,
) -> Tied {
    let started = Instant::now();
    let z = z_levels();
    let (xl, yl) = (quarter_axis(LOWER[0], r_s), quarter_axis(LOWER[1], r_s));
    let (xu, yu) = (quarter_axis(UPPER[0], r_s), quarter_axis(UPPER[1], r_s));
    let lower = shell_quarter(&xl, &yl, 0.5 * (z[0] + z[1]));
    let upper = shell_quarter(&xu, &yu, 0.5 * (z[2] + z[3]));
    let (k_lower, k_upper) = (mitc_spruce(&lower, section), mitc_spruce(&upper, section));
    let thickness = [vec![H; lower.nodes.len()], vec![H; upper.nodes.len()]];
    let components = [
        ShellComponent {
            id: 0,
            mesh: &lower,
            model: &k_lower,
            nodal_thickness_m: &thickness[0],
        },
        ShellComponent {
            id: 1,
            mesh: &upper,
            model: &k_upper,
            nodal_thickness_m: &thickness[1],
        },
    ];
    let tied = tie_shells(&components, &[], budget(), || false).unwrap();

    let s0 = 0.5 * post_side();
    let mut vertices = Vertices::default();
    let mut tets = Vec::new();
    let cross = cells(0.0, s0, r_p);
    kuhn_block(
        &cross,
        &cross,
        &cells(z[1], z[2], POST_CELLS * r_p),
        &mut vertices,
        &mut tets,
    );
    let (nodes, tets10) = promote_p1_to_p2(&vertices.nodes, &tets, None, cx).expect("P2 post");
    let assembly = QuadraticTetProblem {
        nodes_m: &nodes,
        elements: &tets10,
        materials: TetMaterialField::Uniform(post),
        fixed_dofs: &[],
        budget: TetAssemblyBudget::standard(),
        quadrature: TetQuadrature {
            stiffness_points_per_dir: 3,
            mass_points_per_dir: 4,
        },
    }
    .assemble(cx)
    .expect("post assembly");
    let mut dof_index = vec![[None; 3]; nodes.len()];
    for (row, &d) in assembly.free_dofs.iter().enumerate() {
        dof_index[d / 3][d % 3] = Some(row);
    }
    let footprint = |mesh: &ShellMesh| -> Vec<usize> {
        (0..mesh.tris.len())
            .filter(|&t| {
                let c = mesh.tris[t].iter().fold([0.0; 2], |a, &n| {
                    [a[0] + mesh.nodes[n][0] / 3.0, a[1] + mesh.nodes[n][1] / 3.0]
                });
                c[0] < s0 && c[1] < s0
            })
            .collect()
    };
    let solid = SolidOnShells {
        trace: SolidFaceTrace {
            positions_m: &nodes,
            dof_index: &dof_index,
        },
        k: &assembly.stiffness,
        m: &assembly.mass,
        ties: vec![
            SolidFaceTie {
                faces: p2_faces(&nodes, &tets10, |p| p[2] == z[1]),
                master_component: 0,
                master_tris: footprint(&lower),
                side: 1,
            },
            SolidFaceTie {
                faces: p2_faces(&nodes, &tets10, |p| p[2] == z[2]),
                master_component: 1,
                master_tris: footprint(&upper),
                side: -1,
            },
        ],
    };
    let glued = tie_solids_to_shells(&tied.model, &tied, &components, &[solid], budget(), || {
        false
    })
    .expect("post glued to both plates");

    // The doubly symmetric class.
    let mut removed = shell_symmetry_coordinates(&tied, &[&lower, &upper], EVEN_EVEN);
    for (n, p) in nodes.iter().enumerate() {
        for axis in 0..2 {
            if p[axis] != 0.0 {
                continue;
            }
            for &c in EVEN_EVEN.fixed_translations(axis) {
                let Some(row) = dof_index[n][c] else {
                    continue;
                };
                // An untied post coordinate is one appended combined coordinate; a glued one maps
                // onto shell coordinates only.
                let mut unit = vec![0.0; assembly.free_dofs.len()];
                unit[row] = 1.0;
                let load = glued.project_solid_load(0, &unit).unwrap();
                removed.extend((glued.shell_free..glued.model.free).filter(|&i| load[i] != 0.0));
            }
        }
    }
    removed.sort_unstable();
    removed.dedup();
    let (k, m, kept) = restrict(&glued.model.k, &glued.model.m, &removed);
    let index_of = |c: usize, node: usize| {
        let combined = tied.map[c][node][2].expect("a free w coordinate");
        kept.binary_search(&combined).expect("w is never removed")
    };

    // Modes. The quarter keeps one rigid mode, the translation along z.
    let report = slice_window(&k, &m, window(TIED_BAND), &options()).expect("tied window");
    assert_eq!(
        report.below_low, 1,
        "the symmetric quarter has exactly one rigid mode below the window"
    );
    let sample_rows: Vec<usize> = [(&xl, &yl, 0), (&xu, &yu, 1)]
        .iter()
        .flat_map(|&(xs, ys, c)| {
            let mut rows = Vec::new();
            for i in (0..xs.len()).step_by(r_s) {
                for j in (0..ys.len()).step_by(r_s) {
                    rows.push(index_of(c, j * xs.len() + i));
                }
            }
            rows
        })
        .collect();
    let modes: Vec<Mode> = report
        .modes
        .iter()
        .map(|mode| Mode {
            hz: hz(mode),
            samples: sample_rows.iter().map(|&r| mode.phi[r]).collect(),
        })
        .collect();

    // Rim squeeze: the line load F / rim length on the quarter's rim edges x = X and y = Y (consistent:
    // half of each edge's share at each end), with the lower plate's centre held along z.
    let rim_load = |mesh: &ShellMesh, xs: &[f64], ys: &[f64], half: [f64; 2], sign: f64| {
        let id = |i: usize, j: usize| j * xs.len() + i;
        let (nx, ny) = (xs.len() - 1, ys.len() - 1);
        let line = sign / (4.0 * (half[0] + half[1]));
        let mut f = vec![0.0; 6 * mesh.nodes.len()];
        let mut edge = |a: usize, b: usize| {
            let (p, q) = (mesh.nodes[a], mesh.nodes[b]);
            let length = ((q[0] - p[0]).powi(2) + (q[1] - p[1]).powi(2)).sqrt();
            f[6 * a + 2] += 0.5 * line * length;
            f[6 * b + 2] += 0.5 * line * length;
        };
        for i in 0..nx {
            edge(id(i, ny), id(i + 1, ny));
        }
        for j in 0..ny {
            edge(id(nx, j), id(nx, j + 1));
        }
        f
    };
    let loads = [
        rim_load(&lower, &xl, &yl, SQUEEZE[0].0, SQUEEZE[0].1),
        rim_load(&upper, &xu, &yu, SQUEEZE[1].0, SQUEEZE[1].1),
    ];
    let combined = tied.project_loads(&[&loads[0], &loads[1]]).unwrap();
    let f: Vec<f64> = kept
        .iter()
        .map(|&i| if i < combined.len() { combined[i] } else { 0.0 })
        .collect();
    let u = solve_pinned(&k, &f, &[index_of(0, 0)], None);
    let work_quarter: f64 = f.iter().zip(&u).map(|(a, b)| a * b).sum();
    let approach = u[index_of(0, 0)] - u[index_of(1, 0)];
    let seconds = started.elapsed().as_secs_f64();
    println!(
        "shell rung {r_s} (post rung {r_p}): {} of {} coordinates kept, {} modes in the window, tie gap {:.1e} m, {} eliminated, {seconds:.1} s",
        kept.len(),
        glued.model.free,
        modes.len(),
        glued.receipt.max_gap_m,
        glued.receipt.slave_dofs
    );
    Tied {
        modes,
        k_system: 1.0 / (4.0 * work_quarter),
        k_post_path: 1.0 / approach,
        dofs: kept.len(),
        receipt: glued.receipt,
        seconds,
    }
}

fn json_list(values: &[f64], digits: usize) -> String {
    let items: Vec<String> = values.iter().map(|v| format!("{v:.digits$}")).collect();
    format!("[{}]", items.join(","))
}

#[test]
fn g2_shell_solid_tie_soundpost_fixture() {
    let started = Instant::now();
    let s = spruce();
    let oriented = |axes: [[f64; 3]; 3], label: &[u8]| {
        TetElasticMaterial::try_from_resolved_orthotropic_state(
            &s.state,
            axes,
            fs_blake3::hash_domain("org.frankensim.pexfp52.axes.v1", label),
        )
        .expect("oriented spruce tensor")
    };
    // Quarter-sawn plates: L along x, R along y, T through the thickness. The post: L along its axis, R
    // along x, T along y. The columns of `principal_to_world` are the principal axes in world coordinates.
    let plate = oriented(
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        b"plates x=L,y=R,z=T",
    );
    let post = oriented(
        [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]],
        b"post x=R,y=T,z=L",
    );
    let section = ShearShellSection::from_orthotropic_solid(
        s.law.e[0],
        s.law.e[1],
        s.law.nu[0],
        s.law.g[0],
        s.law.g[2],
        s.law.g[1],
        H,
        s.rho,
        0.0,
        KAPPA,
    )
    .expect("shear shell section");
    let a = post_side();
    println!(
        "square post a = r sqrt(pi) = {:.4} mm (area {:.3} mm2, the round post's); its bending second moment a^4/12 is pi/3 = +4.7 % of the round post's pi r^4/4; E_L {:.3} GPa, E_T {:.3} GPa, rho {:.1} kg/m3; the bare post's axial stiffness E_L a^2 / L = {:.4e} N/m",
        1e3 * a,
        1e6 * a * a,
        1e-9 * s.law.e[0],
        1e-9 * s.law.e[2],
        s.rho,
        s.law.e[0] * a * a / POST_L
    );
    with_cx(|cx| {
        let tied: Vec<Tied> = RUNGS_SHELL
            .iter()
            .map(|&r| tied_model(r, RUNGS_3D[2], &section, &post, cx))
            .collect();
        let finest = tied.last().unwrap();
        assert!(
            finest.modes.len() > GATED,
            "only {} doubly symmetric modes below {} Hz",
            finest.modes.len(),
            TIED_BAND.1
        );
        // The 3-D window reaches 30 % above the tied model's eleventh symmetric mode, so the stiffer
        // coarse rungs still hold every gated mode's partner.
        let band = (1.0, 1.3 * finest.modes[GATED].hz);
        let references: Vec<Reference3d> = RUNGS_3D
            .iter()
            .map(|&r| reference_3d(r, &plate, &post, band, cx))
            .collect();
        for (r, reference) in RUNGS_3D.iter().zip(&references) {
            assert_eq!(
                reference.below, 1,
                "3-D rung {r}: the symmetric quarter has one rigid mode"
            );
        }
        let fine = &references[2].modes;
        assert!(fine.len() >= GATED, "only {} 3-D modes", fine.len());

        // The reference: each of the first ten finest-rung modes tracked by MAC to the coarser rungs.
        let mut worst_bound = 0.0_f64;
        let mut worst_dev = 0.0_f64;
        let (mut f_3d, mut f_tied, mut devs) = (Vec::new(), Vec::new(), Vec::new());
        for (index, mode) in fine.iter().take(GATED).enumerate() {
            let mut rungs = [0.0; 3];
            rungs[2] = mode.hz;
            for rung in 0..2 {
                rungs[rung] = track(index, fine, &references[rung].modes).unwrap_or_else(|e| {
                    panic!(
                        "3-D mode {index} at {:.2} Hz untracked on rung {rung}: {e:?}",
                        mode.hz
                    )
                });
            }
            let (order, reference, bound) = ladder(rungs, RUNGS_3D);
            worst_bound = worst_bound.max(bound / reference);
            let partners: Vec<f64> = tied
                .iter()
                .map(|t| {
                    track(index, fine, &t.modes).unwrap_or_else(|e| {
                        panic!("3-D mode {index} at {reference:.2} Hz has no tied partner: {e:?}")
                    })
                })
                .collect();
            let f = *partners.last().unwrap();
            let dev = f / reference - 1.0;
            worst_dev = worst_dev.max(dev.abs());
            println!(
                "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W3\",\"test\":\"g2_shell_solid_tie_soundpost_fixture\",\"tier\":\"G2\",\"inputs\":{{\"fixture\":\"soundpost\",\"tie\":\"rigid\",\"pack\":\"spruce-sitka-fpl-gtr282\",\"mode\":{index}}},\"rungs\":{{\"tet10_layers\":{:?},\"shell_cells_per_quarter_axis\":{:?}}},\"observed\":{{\"f_3d_rungs\":{},\"order\":{order:.3},\"f_3d\":{reference:.4},\"bound_hz\":{bound:.3e},\"f_tied\":{},\"dev\":{dev:.3e}}},\"tolerance\":{TOL}}}",
                RUNGS_3D.map(|r| LAYERS * r),
                RUNGS_SHELL.map(|r| (OUTER_CELLS + 1) * r),
                json_list(&rungs, 4),
                json_list(&partners, 4)
            );
            f_3d.push(reference);
            f_tied.push(f);
            devs.push(dev);
        }

        // Rim squeeze: the system stiffness is gated; the post-path stiffness, a local quantity the shell
        // tie does not reproduce (see the module doc), is printed.
        let k3: [f64; 3] = [0, 1, 2].map(|i| references[i].k_system);
        let (k_order, k_reference, k_bound) = ladder(k3, RUNGS_3D);
        let k_tied: Vec<f64> = tied.iter().map(|t| t.k_system).collect();
        let k_dev = k_tied.last().unwrap() / k_reference - 1.0;
        let p3: [f64; 3] = [0, 1, 2].map(|i| references[i].k_post_path);
        let (p_order, p_reference, p_bound) = ladder(p3, RUNGS_3D);
        let p_tied: Vec<f64> = tied.iter().map(|t| t.k_post_path).collect();
        let p_dev = p_tied.last().unwrap() / p_reference - 1.0;
        println!(
            "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W3\",\"test\":\"g2_shell_solid_tie_soundpost_fixture\",\"tier\":\"G2\",\"inputs\":{{\"fixture\":\"soundpost\",\"load\":\"rim squeeze, 1 N per plate\"}},\"observed\":{{\"k_system_3d_rungs\":{},\"k_system_order\":{k_order:.3},\"k_system_3d\":{k_reference:.6e},\"k_system_bound\":{k_bound:.3e},\"k_system_tied\":{},\"k_system_dev\":{k_dev:.3e},\"k_post_path_3d_rungs\":{},\"k_post_path_order\":{p_order:.3},\"k_post_path_3d\":{p_reference:.6e},\"k_post_path_bound\":{p_bound:.3e},\"k_post_path_tied\":{},\"k_post_path_dev\":{p_dev:.3e}}},\"tolerance\":{{\"k_system\":{TOL},\"k_post_path\":\"ungated\"}}}}",
            json_list(&k3, 2),
            json_list(&k_tied, 2),
            json_list(&p3, 2),
            json_list(&p_tied, 2)
        );
        let dofs_3d: Vec<usize> = references.iter().map(|r| r.dofs).collect();
        let dofs_tied: Vec<usize> = tied.iter().map(|t| t.dofs).collect();
        let seconds_3d: Vec<f64> = references.iter().map(|r| r.seconds).collect();
        let seconds_tied: Vec<f64> = tied.iter().map(|t| t.seconds).collect();
        let verdict = worst_bound <= TOL / 3.0
            && k_bound / k_reference <= TOL / 3.0
            && worst_dev <= TOL
            && k_dev.abs() <= TOL;
        println!(
            "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W3\",\"test\":\"g2_shell_solid_tie_soundpost_fixture\",\"tier\":\"G2\",\"inputs\":{{\"fixture\":\"soundpost\",\"tie\":\"rigid\",\"pack\":\"spruce-sitka-fpl-gtr282\",\"post\":\"square, a = r sqrt(pi), r {POST_R}, L {POST_L}\"}},\"rungs\":{{\"tet10_layers\":{:?},\"dofs_3d\":{dofs_3d:?},\"dofs_tied\":{dofs_tied:?}}},\"observed\":{{\"f_3d\":{},\"f_tied\":{},\"dev\":{},\"max_abs_dev\":{worst_dev:.3e},\"max_rel_bound_3d\":{worst_bound:.3e},\"k_system_dev\":{k_dev:.3e},\"k_post_path_dev\":{p_dev:.3e},\"tie_gap_m\":{:.1e},\"seconds_3d\":{},\"seconds_tied\":{}}},\"tolerance\":{TOL},\"verdict\":\"{}\",\"elapsed_s\":{:.1}}}",
            RUNGS_3D.map(|r| LAYERS * r),
            json_list(&f_3d, 4),
            json_list(&f_tied, 4),
            json_list(&devs, 5),
            finest.receipt.max_gap_m,
            json_list(&seconds_3d, 1),
            json_list(&seconds_tied, 1),
            if verdict { "pass" } else { "fail" },
            started.elapsed().as_secs_f64()
        );
        // Refusal: an unconverged 3-D reference.
        assert!(
            worst_bound <= TOL / 3.0,
            "3-D reference unconverged: frequency bound {worst_bound:.3e} > TOL/3"
        );
        assert!(
            k_bound / k_reference <= TOL / 3.0,
            "3-D reference unconverged: system stiffness bound {:.3e} > TOL/3",
            k_bound / k_reference
        );
        assert!(
            worst_dev <= TOL,
            "a symmetric mode deviates {worst_dev:.3e} from the 3-D reference"
        );
        assert!(
            k_dev.abs() <= TOL,
            "the system stiffness deviates {k_dev:.3e} from the 3-D reference"
        );
    });
}
