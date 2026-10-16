use super::support::*;
use fs_solid::linear3::{TetAssemblyBudget, TetLinearElasticProblem, TetMaterialField};

pub fn cooks_membrane() {
    let started = std::time::Instant::now();
    with_cx(|cx| {
        let m = material(0.9375, 0.25);
        let mut observed = Vec::new();
        let mut dofs = Vec::new();
        // E'=E/(1-nu²), nu'=nu/(1-nu): these plane-strain constants
        // give the plane-stress (E=1,nu=1/3) operator of sol-006.
        for n in [8, 16, 32] {
            let mut mesh = block([n, n, 1], [1.0; 3]);
            for p in &mut mesh.nodes {
                let [s, t, z] = *p;
                *p = [48.0 * s, 44.0 * s + (44.0 - 28.0 * s) * t, z];
            }
            let (nodes, elements) = promote(&mesh, cx);
            let a = assembly(&nodes, &elements, &m, &[], cx);
            let rhs = surface_load(
                &nodes,
                &elements,
                |p| (p[0] - 48.0).abs() < 1e-12,
                [0.0, 1.0 / 16.0, 0.0],
            );
            assert!((rhs.iter().sum::<f64>() - 1.0).abs() < 1e-12);
            let prescribed: Vec<_> = nodes
                .iter()
                .flat_map(|p| {
                    [
                        if p[0].abs() < 1e-12 { Some(0.0) } else { None },
                        if p[0].abs() < 1e-12 { Some(0.0) } else { None },
                        Some(0.0),
                    ]
                })
                .collect();
            let u = solve(&a.stiffness, &rhs, &prescribed);
            let tip = nodes
                .iter()
                .enumerate()
                .find(|(_, p)| {
                    (p[0] - 48.0).abs() < 1e-12 && (p[1] - 52.0).abs() < 1e-12 && p[2].abs() < 1e-12
                })
                .unwrap()
                .0;
            observed.push(u[3 * tip + 1]);
            dofs.push(prescribed.iter().filter(|x| x.is_none()).count());
        }
        let finest = *observed.last().unwrap();
        let change = (observed[2] - observed[1]).abs() / finest.abs();
        log(
            "g2_cooks_membrane_plane_strain_equivalent",
            "G2",
            "Cook-panel-extrusion",
            &[1.0 / 8.0, 1.0 / 16.0, 1.0 / 32.0],
            &dofs,
            &observed,
            0.002,
            "Cook 1974 via sol-006; exact plane-strain mapping",
            started.elapsed().as_secs_f64(),
        );
        assert!((23.5..=24.5).contains(&finest), "Cook tip {finest}");
        assert!(
            change <= 0.002,
            "Cook finest-pair relative difference {change:e}"
        );

        // TET10 on [32,32,1] and P1 on [64,64,2] have the same nodal grid.
        let mut p1 = block([64, 64, 2], [1.0; 3]);
        for p in &mut p1.nodes {
            let [s, t, z] = *p;
            *p = [48.0 * s, 44.0 * s + (44.0 - 28.0 * s) * t, z];
        }
        let a = TetLinearElasticProblem {
            nodes_m: &p1.nodes,
            tetrahedra: &p1.tets,
            materials: TetMaterialField::Uniform(&m),
            fixed_dofs: &[],
            budget: TetAssemblyBudget::standard(),
        }
        .assemble(cx)
        .unwrap();
        // The quadratic trace fixture returns an exactly integrated load on a
        // promoted copy; restricting it to P1 vertices would be incorrect.
        // P1 right edge is a regular grid: integrate its nodal hat functions
        // on each right-face triangle directly, using the corner weights A/3.
        let mut force = vec![0.0; 3 * p1.nodes.len()];
        let mut faces = std::collections::BTreeMap::<[usize; 3], usize>::new();
        for t in &p1.tets {
            for mut f in [
                [t[0], t[1], t[2]],
                [t[0], t[1], t[3]],
                [t[0], t[2], t[3]],
                [t[1], t[2], t[3]],
            ] {
                f.sort_unstable();
                *faces.entry(f).or_default() += 1;
            }
        }
        for (f, count) in faces {
            if count == 1 && f.iter().all(|&i| (p1.nodes[i][0] - 48.0).abs() < 1e-12) {
                let [a, b, c] = f.map(|i| p1.nodes[i]);
                let area =
                    ((b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])).abs() / 2.0;
                for i in f {
                    force[3 * i + 1] += area / 48.0;
                }
            }
        }
        let prescribed: Vec<_> = p1
            .nodes
            .iter()
            .flat_map(|p| {
                [
                    if p[0].abs() < 1e-12 { Some(0.0) } else { None },
                    if p[0].abs() < 1e-12 { Some(0.0) } else { None },
                    Some(0.0),
                ]
            })
            .collect();
        assert_eq!(
            *dofs.last().unwrap(),
            prescribed.iter().filter(|x| x.is_none()).count()
        );
        let u = solve(&a.stiffness, &force, &prescribed);
        let tip = p1
            .nodes
            .iter()
            .enumerate()
            .find(|(_, p)| {
                (p[0] - 48.0).abs() < 1e-12 && (p[1] - 52.0).abs() < 1e-12 && p[2].abs() < 1e-12
            })
            .unwrap()
            .0;
        println!(
            "{{\"suite\":\"irp-s1\",\"test\":\"g2_cooks_membrane_plane_strain_equivalent\",\"element\":\"tet4\",\"matched_free_dofs\":{},\"tip\":{},\"difference_from_tet10\":{},\"literature_band\":[23.5,24.5]}}",
            dofs[2],
            u[3 * tip + 1],
            u[3 * tip + 1] - finest
        );
    });
}

pub fn macneal_harder() {
    // MacNeal & Harder, Finite Elements in Analysis and Design 1 (1985),
    // 3–20, Fig. 4 and Table 3 (p.10), doi:10.1016/0168-874X(85)90003-4:
    // L=6, width=0.2, thickness=0.1, E=1e7, nu=0.3, unit tip shear.
    // Table 3: in-plane 0.1081, out-of-plane 0.4321. These include shear;
    // clamped-end 3-D effects remain disclosed offsets, never fitted away.
    let started = std::time::Instant::now();
    with_cx(|cx| {
        for (axis, reference) in [(1, 0.1081), (2, 0.4321)] {
            let mut ratios = Vec::new();
            let mut dofs = Vec::new();
            for refinement in [1, 2, 4] {
                let mesh = block([6 * refinement, refinement, refinement], [6.0, 0.2, 0.1]);
                let (nodes, elements) = promote(&mesh, cx);
                let a = assembly(&nodes, &elements, &material(1e7, 0.3), &[], cx);
                let mut traction = [0.0; 3];
                traction[axis] = 1.0 / 0.02;
                let force =
                    surface_load(&nodes, &elements, |p| (p[0] - 6.0).abs() < 1e-12, traction);
                assert!((force.iter().sum::<f64>() - 1.0).abs() < 1e-12);
                let prescribed: Vec<_> = nodes
                    .iter()
                    .flat_map(|p| [if p[0].abs() < 1e-12 { Some(0.0) } else { None }; 3])
                    .collect();
                let u = solve(&a.stiffness, &force, &prescribed);
                let mean: f64 = force.iter().zip(&u).map(|(a, b)| a * b).sum();
                ratios.push(mean / reference);
                dofs.push(prescribed.iter().filter(|x| x.is_none()).count());
            }
            log(
                "g2_macneal_harder_straight_cantilever",
                "G2",
                if axis == 1 {
                    "MacNeal-Harder-in-plane"
                } else {
                    "MacNeal-Harder-out-of-plane"
                },
                &[1.0, 0.5, 0.25],
                &dofs,
                &ratios,
                0.015,
                "MacNeal-Harder 1985 Fig.4 Table 3",
                started.elapsed().as_secs_f64(),
            );
            println!(
                "{{\"suite\":\"irp-s1\",\"beam_load_axis\":{axis},\"reference\":{reference},\"end_effect_relative_offset\":{}}}",
                ratios[2] - 1.0
            );
            assert!(
                (0.985..=1.015).contains(&ratios[2]),
                "straight beam axis {axis}: {ratios:?}"
            );
        }
    });
}

pub fn free_bar() {
    let started = std::time::Instant::now();
    with_cx(|cx| {
        let mut values = Vec::new();
        let mut dofs = Vec::new();
        // Central verification can request one additional diagnostic rung.
        // Every original rung and its original order gate remain active.
        let rungs: &[usize] = if std::env::var_os("IRP_TET10_BAR_DIAGNOSTIC_R8").is_some() {
            &[1, 2, 4, 8]
        } else {
            &[1, 2, 4]
        };
        for &r in rungs {
            let mesh = block([6 * r, r, r], [1.0, 0.1, 0.07]);
            let (nodes, elements) = promote(&mesh, cx);
            let a = assembly(&nodes, &elements, &material(1.0, 0.25), &[], cx);
            if r == 1 {
                // An affine P2 mass needs degree four, already exact for
                // Q.mass_points_per_dir=4. Compare assembled operators too:
                // a correct scalar rule alone would not catch using K's rule
                // accidentally for M inside assembly.
                let m = material(1.0, 0.25);
                let higher = fs_solid::linear3::tet10::QuadraticTetProblem {
                    nodes_m: &nodes,
                    elements: &elements,
                    materials: TetMaterialField::Uniform(&m),
                    fixed_dofs: &[],
                    budget: TetAssemblyBudget::standard(),
                    quadrature: fs_solid::linear3::tet10::TetQuadrature {
                        stiffness_points_per_dir: 5,
                        mass_points_per_dir: 5,
                    },
                }
                .assemble(cx)
                .unwrap();
                for (name, low, high) in [
                    ("stiffness", &a.stiffness, &higher.stiffness),
                    ("mass", &a.mass, &higher.mass),
                ] {
                    let low = low.to_dense();
                    let high = high.to_dense();
                    let scale = high.iter().fold(0.0_f64, |a, b| a.max(b.abs()));
                    let difference = low
                        .iter()
                        .zip(&high)
                        .map(|(a, b)| (a - b).abs())
                        .fold(0.0_f64, f64::max);
                    println!(
                        "{{\"suite\":\"irp-s1\",\"test\":\"g1_free_free_bar_eigen_ladder\",\"quadrature_axis\":\"{name}\",\"base_n\":[3,4],\"higher_n\":[5,5],\"relative_max_entry_difference\":{},\"tolerance\":1e-12}}",
                        difference / scale
                    );
                    assert!(difference <= 1e-12 * scale, "{name} quadrature mismatch");
                }
            }
            let report = fs_modal::slice_window(
                &a.stiffness,
                &a.mass,
                (1e-5, 8.0),
                &fs_modal::SliceOptions {
                    max_lanczos: 512,
                    max_restarts: 24,
                    ..fs_modal::SliceOptions::default()
                },
            )
            .unwrap();
            // Torsion interleaves the two bending families. Classify by
            // area-mean transverse motion on the end section: torsion and
            // axial modes have zero mean, whereas bending translates it.
            let end_y = surface_load(
                &nodes,
                &elements,
                |p| p[0].abs() < 1e-12,
                [0.0, 1.0 / 0.007, 0.0],
            );
            let end_z = surface_load(
                &nodes,
                &elements,
                |p| p[0].abs() < 1e-12,
                [0.0, 0.0, 1.0 / 0.007],
            );
            let bending: Vec<_> = report
                .modes
                .iter()
                .filter(|m| {
                    let mean = |load: &[f64]| {
                        m.phi
                            .iter()
                            .zip(&a.free_dofs)
                            .map(|(u, &i)| u * load[i])
                            .sum::<f64>()
                    };
                    let y = mean(&end_y);
                    let z = mean(&end_z);
                    let rms =
                        (m.phi.iter().map(|u| u * u).sum::<f64>() / nodes.len() as f64).sqrt();
                    (y * y + z * z).sqrt() > 0.2 * rms
                })
                .take(4)
                .collect();
            assert_eq!(
                bending.len(),
                4,
                "window must contain four identified bending modes"
            );
            let row: Vec<_> = bending.iter().map(|m| m.lambda).collect();
            let residuals: Vec<_> = bending.iter().map(|m| m.residual).collect();
            println!(
                "{{\"suite\":\"irp-s1\",\"test\":\"g1_free_free_bar_eigen_ladder\",\"refinement\":{r},\"free_dofs\":{},\"minimum_scaled_jacobian\":{},\"eigenvalues\":{row:?},\"eigenvalue_residual_bounds\":{residuals:?}}}",
                a.free_dofs.len(),
                a.minimum_scaled_jacobian
            );
            log(
                "g1_free_free_bar_eigen_ladder",
                "G1",
                "free-free-prismatic-bar",
                &[1.0 / r as f64],
                &[a.free_dofs.len()],
                &row,
                0.3,
                "Richardson limit; E=rho=1, nu=0.25",
                started.elapsed().as_secs_f64(),
            );
            values.push(row);
            dofs.push(a.free_dofs.len());
        }
        let mut orders = Vec::new();
        for first in 0..values.len() - 2 {
            for mode in 0..4 {
                let d0 = values[first][mode] - values[first + 1][mode];
                let d1 = values[first + 1][mode] - values[first + 2][mode];
                assert!(
                    d0 > 0.0 && d1 > 0.0,
                    "eigenvalue ladder not live/monotone: {values:?}"
                );
                let order = fs_math::det::ln(d0 / d1) / fs_math::det::ln(2.0);
                let limit = values[first + 2][mode] - d1 / ((d0 / d1) - 1.0);
                println!(
                    "{{\"suite\":\"irp-s1\",\"test\":\"g1_free_free_bar_eigen_ladder\",\"mode\":{mode},\"refinements\":{:?},\"order\":{{\"observed\":{order},\"theoretical\":4.0,\"band\":0.3}},\"richardson_limit\":{limit}}}",
                    &rungs[first..first + 3]
                );
                orders.push(order);
            }
        }
        // Euler–Bernoulli roots are a model-form comparison, not the 3-D
        // eigenvalue gate. Both transverse inertias are represented.
        for beta in [4.73004_f64, 7.85320, 10.99561] {
            let root4 = (beta * beta) * (beta * beta);
            println!(
                "{{\"suite\":\"irp-s1\",\"model_form\":\"Euler-Bernoulli-free-free\",\"beta\":{beta},\"lambda_y\":{},\"lambda_z\":{},\"gated\":false}}",
                root4 * 0.01 / 12.0,
                root4 * 0.0049 / 12.0
            );
        }
        assert!(
            orders.iter().all(|p| (p - 4.0).abs() <= 0.3),
            "eigen orders {orders:?}: rungs {values:?}, dofs {dofs:?}"
        );
    });
}
