//! G0/G1 composition: outline -> plate -> explicit arch/graduation -> shell.
//! Annular membrane oracle: J0(ka)Y0(kb)-J0(kb)Y0(ka)=0 with fixed edges.
use fs_exec::{CancelGate, Cx};
use fs_mesh::planar::*;
use fs_modal::{ModePair, SliceOptions, slice_window};
use fs_plate::shell::arching::*;
use fs_plate::shell::mitc::*;
use fs_plate::shell::profile::ProfileBudget;
use fs_plate::shell::survey::{IsotropicMaterial, MeshShell, SurveyMaterial};
use fs_plate::{AssemblyOptions, EdgeSupport, PlateMesh, PlateSection, assemble};
use fs_sparse::{Coo, Csr, DirectOrdering, LdltOptions, SymbolicLdlt};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 55,
                kernel_id: 55,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}
fn budget() -> PlanarBudget {
    PlanarBudget {
        max_vertices: 100000,
        max_triangles: 200000,
        max_refinement_steps: 100000,
    }
}
fn options(h: f64) -> PlanarOptions<'static> {
    PlanarOptions {
        maximum_edge_m: h,
        minimum_angle_deg: 20.0,
        size_field: None,
        budget: budget(),
    }
}
fn plate(mesh: &PlanarMesh) -> PlateMesh {
    PlateMesh::from_unstructured(
        mesh.vertices.iter().map(|p| (p[0], p[1])).collect(),
        mesh.triangles.clone(),
    )
    .unwrap()
}
fn grid(f: impl Fn(f64, f64) -> f64) -> BicubicGrid {
    let axis = vec![-2.0, -1.0, 0.0, 1.0, 2.0];
    let values = axis
        .iter()
        .flat_map(|&y| axis.iter().map(|&x| f(x, y)).collect::<Vec<_>>())
        .collect();
    BicubicGrid::new(axis.clone(), axis, values).unwrap()
}
fn profile_budget(mesh: &PlateMesh) -> ProfileBudget {
    ProfileBudget {
        max_nodes: mesh.nodes.len(),
        max_triangles: mesh.tris.len(),
        max_feature_evaluations: 100 * mesh.nodes.len(),
    }
}
fn ring(n: usize, r: f64) -> Vec<[f64; 2]> {
    // Explicit dyadic polygon coordinates permit exact midpoint h-refinement.
    // This quantization is geometry data, bounded by sqrt(2)*2^-31 m,
    // not an unreported projection during refinement.
    (0..n)
        .map(|i| {
            let t = 2.0 * core::f64::consts::PI * i as f64 / n as f64;
            [r * fs_math::det::cos(t), r * fs_math::det::sin(t)]
                .map(|x| (x * 1073741824.0).round() / 1073741824.0)
        })
        .collect()
}
fn annulus(n: usize, h: f64, cx: &Cx<'_>) -> PlanarMesh {
    let outer = ring(n, 1.0);
    let inner = ring(n, 0.5);
    mesh_outline(
        &PlanarOutline {
            outer: Loop::Polyline(&outer),
            holes: vec![Loop::Polyline(&inner)],
        },
        &options(h),
        cx,
    )
    .unwrap()
}

#[test]
fn g0_arching_surface_from_parameters() {
    with_cx(|cx| {
        let planar = annulus(32, 0.3, cx);
        let mesh = plate(&planar);
        let state = ReferenceLoadState::Loaded {
            description: "analytic fixture; explicitly loaded reference, no unloading".into(),
        };
        let arch = ArchingSurface::new(
            ArchingParameters::HeightMap(grid(|x, y| 0.02 * (1.0 - x * x) * (1.0 - y * y))),
            state.clone(),
        )
        .unwrap();
        let graduation = ThicknessField::new(ArchingParameters::HeightMap(grid(|x, y| {
            0.003 + 0.0001 * (x + y)
        })))
        .unwrap();
        let lifted = arch.lift(&mesh, profile_budget(&mesh), || false).unwrap();
        assert_eq!(lifted.reference_state, state);
        assert_eq!(lifted.mesh.tris, planar.triangles);
        let directors = NodalDirectors::new(&lifted.mesh, lifted.directors, 0.5).unwrap();
        assert_eq!(directors.vectors().len(), mesh.nodes.len());
        let thickness = graduation
            .sample_mesh(&mesh, mesh.nodes.len(), || false)
            .unwrap();
        let materials = vec![
            SurveyMaterial::Isotropic(IsotropicMaterial {
                young_pa: 1e7,
                poisson: 0.3,
                density_kg_m3: 600.0
            });
            mesh.tris.len()
        ];
        let shell = MeshShell::new_with_materials(
            lifted.mesh,
            thickness,
            &materials,
            profile_budget(&mesh),
        )
        .unwrap();
        assert_eq!(shell.mesh.tris, mesh.tris);
        assert_eq!(shell.mesh.nodes.len(), mesh.nodes.len());
        // Lift is a graph embedding, so boundary connectivity (including hole)
        // is identical, rather than reconstructed or silently filled.
        assert_eq!(
            mesh.boundary_nodes().len(),
            planar
                .receipt
                .boundary_loops
                .iter()
                .map(Vec::len)
                .sum::<usize>()
        );
        for (p, &(x, y)) in shell.mesh.nodes.iter().zip(&mesh.nodes) {
            assert!((p[2] - 0.02 * (1.0 - x * x) * (1.0 - y * y)).abs() < 1e-12);
        }
        println!(
            "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W6\",\"test\":\"g0_arching_surface_from_parameters\",\"tier\":\"G0\",\"inputs\":{{\"holes\":1,\"arching\":\"bicubic-heightmap\",\"graduation\":\"independent\",\"reference_load_state\":\"loaded\"}},\"vertices\":{},\"triangles\":{},\"verdict\":true}}",
            shell.mesh.nodes.len(),
            shell.mesh.tris.len()
        );
    });
}

fn shifted(k: &Csr, m: &Csr, s: f64) -> Csr {
    let mut out = Coo::new(k.nrows(), k.nrows());
    for i in 0..k.nrows() {
        let (c, v) = k.row(i);
        for (&j, &a) in c.iter().zip(v) {
            out.push(i, j, a);
        }
        let (c, v) = m.row(i);
        for (&j, &a) in c.iter().zip(v) {
            out.push(i, j, -s * a);
        }
    }
    out.assemble()
}
fn first_ten(k: &Csr, m: &Csr) -> Vec<ModePair> {
    let mut hi = 1.0;
    for _ in 0..32 {
        let pencil = shifted(k, m, hi);
        let symbolic = SymbolicLdlt::analyze(&pencil, DirectOrdering::Amd).unwrap();
        let count = symbolic
            .factor(&pencil, &LdltOptions::default())
            .unwrap()
            .inertia()
            .negative;
        if count >= 10 {
            let report = slice_window(
                k,
                m,
                (0.0, hi),
                &SliceOptions {
                    ritz_tol: 1e-12,
                    max_lanczos: 1024,
                    max_restarts: 64,
                    mass_diagonal_equilibration: true,
                    ..SliceOptions::default()
                },
            )
            .unwrap();
            assert!(report.modes.len() >= 10);
            return report.modes.into_iter().take(10).collect();
        }
        hi *= 2.0;
    }
    panic!("explicit first-ten spectral search budget exhausted");
}

#[test]
fn g1_plate_with_hole_modes_converge() {
    with_cx(|cx| {
        let mut planar = annulus(64, 0.25, cx);
        let mut dkt = Vec::new();
        let mut mitc = Vec::new();
        let mut hs = Vec::new();
        for rung in 0..4 {
            let mesh = plate(&planar);
            let boundary = mesh.boundary_nodes();
            let section = PlateSection::isotropic(1.0, 0.3, 0.01, 1.0).unwrap();
            let model = assemble(
                &mesh,
                &section,
                &boundary,
                &[],
                &AssemblyOptions {
                    pretension: 0.0,
                    support: EdgeSupport::Clamped,
                },
            )
            .unwrap();
            let modes = first_ten(&model.k, &model.m);
            dkt.push(modes.iter().map(|m| m.lambda).collect::<Vec<_>>());
            let flat = ArchingSurface::new(
                ArchingParameters::HeightMap(grid(|_, _| 0.0)),
                ReferenceLoadState::Unloaded,
            )
            .unwrap()
            .lift(&mesh, profile_budget(&mesh), || false)
            .unwrap();
            let constraints: Vec<_> = boundary
                .iter()
                .map(|&node| ShellConstraint {
                    node,
                    mask: [true; 6],
                    frame: ConstraintFrame::Global,
                })
                .collect();
            let section = ShearShellSection::from_orthotropic_solid(
                1.0,
                1.0,
                0.3,
                1.0 / 2.6,
                1.0 / 2.6,
                1.0 / 2.6,
                0.01,
                1.0,
                0.0,
                5.0 / 6.0,
            )
            .unwrap();
            let directors = NodalDirectors::new(&flat.mesh, flat.directors, 0.5).unwrap();
            let shell = assemble_mitc3p(
                &flat.mesh,
                &directors,
                &vec![section; mesh.tris.len()],
                &constraints,
                ShellMass::Consistent,
                MitcBudget {
                    max_nodes: mesh.nodes.len(),
                    max_triangles: mesh.tris.len(),
                },
                1e-3,
                || false,
            )
            .unwrap();
            let modes = first_ten(&shell.pencil().k, &shell.pencil().m);
            mitc.push(modes.iter().map(|m| m.lambda).collect::<Vec<_>>());
            hs.push(1.0 / (1usize << rung) as f64);
            println!(
                "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W6\",\"test\":\"g1_plate_with_hole_modes_converge\",\"tier\":\"G1\",\"rung\":{rung},\"geometry_axis\":\"fixed-64-edge-dyadic-annulus\",\"vertices\":{},\"triangles\":{},\"dkt\":{:?},\"mitc3p\":{:?}}}",
                mesh.nodes.len(),
                mesh.tris.len(),
                dkt.last().unwrap(),
                mitc.last().unwrap()
            );
            if rung < 3 {
                planar = planar.refine_uniform(budget(), cx).unwrap();
            }
        }
        for (name, rows) in [("DKT", dkt), ("MITC3+", mitc)] {
            for mode in 0..10 {
                // Fixed geometry; three successive eigenvalue differences
                // remove an unknown continuum limit without using a fixture
                // answer. Nonzero changes make the h axis load-bearing.
                let errors: Vec<_> = rows
                    .windows(2)
                    .map(|p| (p[1][mode] - p[0][mode]).abs())
                    .collect();
                let ladder = fs_mms::RefinementLadder::new(hs[..3].to_vec(), errors).unwrap();
                let verdict = fs_mms::OrderGate { theoretical: 2.0 }
                    .check(
                        &format!("{name}-annulus-mode-{mode}"),
                        fs_mms::LadderSide::Primal,
                        &ladder,
                    )
                    .unwrap();
                println!("{}", verdict.json_line(true));
            }
        }
    });
}

fn annular_roots(a: f64, b: f64, count: usize) -> Vec<f64> {
    use fs_aeroac::bessel::{j0, y0};
    let f = |k: f64| j0(k * a) * y0(k * b) - j0(k * b) * y0(k * a);
    let mut roots = Vec::new();
    let mut lo = 0.01;
    let mut left = f(lo);
    for i in 1..=1000 {
        let hi = 0.01 + 0.05 * i as f64;
        let right = f(hi);
        if left * right < 0.0 {
            let mut l = lo;
            let mut h = hi;
            let mut fl = left;
            for _ in 0..55 {
                let mid = f64::midpoint(l, h);
                let fm = f(mid);
                if fl * fm <= 0.0 {
                    h = mid;
                } else {
                    l = mid;
                    fl = fm;
                }
            }
            let root = f64::midpoint(l, h);
            assert!(f(root).abs() < 1e-12);
            roots.push(root);
            if roots.len() == count {
                return roots;
            }
        }
        lo = hi;
        left = right;
    }
    panic!("annulus root-search budget exhausted");
}

#[test]
fn g1_annular_tensioned_film_matches_bessel_roots() {
    with_cx(|cx| {
        let roots = annular_roots(0.5, 1.0, 2);
        let mut errors = vec![Vec::new(); 2];
        for h in [0.2, 0.1, 0.05] {
            let chord_bound_m = 1.0 - fs_math::det::cos(core::f64::consts::PI / 256.0);
            println!(
                "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W6\",\"stage\":\"annulus-geometry-axis\",\"edges_per_loop\":256,\"fixed_across_h_rungs\":true,\"chord_bound_m\":{chord_bound_m},\"quantization_bound_m\":6.6e-10}}"
            );
            let planar = annulus(256, h, cx);
            let mesh = plate(&planar);
            let section = PlateSection::isotropic(1e3, 0.3, 0.001, 1000.0).unwrap();
            let model = assemble(
                &mesh,
                &section,
                &mesh.boundary_nodes(),
                &[],
                &AssemblyOptions {
                    pretension: 2000.0,
                    support: EdgeSupport::SimplySupported,
                },
            )
            .unwrap();
            for (mode, &k) in roots.iter().enumerate() {
                let expected = 2000.0 * k * k; // T/(rho*h)=2000 m²/s².
                let report = slice_window(
                    &model.k,
                    &model.m,
                    (0.8 * expected, 1.2 * expected),
                    &SliceOptions {
                        ritz_tol: 1e-12,
                        max_lanczos: 1024,
                        max_restarts: 64,
                        ..SliceOptions::default()
                    },
                )
                .unwrap();
                // Pair by the independent Bessel mode SHAPE, not by picking
                // whichever eigenvalue happens to be closest to the answer.
                let mut reference = vec![0.0; model.free];
                for (i, &(x, y)) in mesh.nodes.iter().enumerate() {
                    if let Some(dof) = model.dof_map[3 * i] {
                        let r = fs_math::det::sqrt(x * x + y * y);
                        reference[dof] = fs_aeroac::bessel::j0(k * r)
                            * fs_aeroac::bessel::y0(k * 0.5)
                            - fs_aeroac::bessel::j0(k * 0.5) * fs_aeroac::bessel::y0(k * r);
                    }
                }
                let mut mr = vec![0.0; model.free];
                model.m.spmv(&reference, &mut mr);
                let norm = reference
                    .iter()
                    .zip(&mr)
                    .map(|(a, b)| a * b)
                    .sum::<f64>()
                    .sqrt();
                let overlap = |m: &ModePair| {
                    m.phi.iter().zip(&mr).map(|(a, b)| a * b).sum::<f64>().abs() / norm
                };
                let best = report
                    .modes
                    .iter()
                    .max_by(|a, b| overlap(a).total_cmp(&overlap(b)))
                    .expect("axisymmetric mode in declared window");
                assert!(
                    overlap(best) > 0.95,
                    "mode {mode} is not the axisymmetric Bessel shape"
                );
                let error = (best.lambda / expected - 1.0).abs();
                errors[mode].push(error);
                println!(
                    "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W6\",\"test\":\"g1_annular_tensioned_film_matches_bessel_roots\",\"tier\":\"G1\",\"h_m\":{h},\"geometry_axis_edges\":256,\"quantization_bound_m\":6.6e-10,\"mode\":{mode},\"root\":{k},\"expected_lambda\":{expected},\"lambda\":{},\"shape_overlap\":{},\"relative_error\":{error},\"tolerance\":0.005}}",
                    best.lambda,
                    overlap(best)
                );
            }
        }
        for row in errors {
            assert!(
                row.windows(2).all(|p| p[1] < p[0]),
                "annular errors {row:?}"
            );
            assert!(*row.last().unwrap() < 0.005);
        }
    });
}
