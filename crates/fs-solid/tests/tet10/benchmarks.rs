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
                let u = solve_conditioned(&a.stiffness, &force, &prescribed);
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

pub fn free_bar(qualification: bool) {
    let started = std::time::Instant::now();
    with_cx(|cx| {
        let mut values = Vec::new();
        let mut bounds = Vec::new();
        let mut dofs = Vec::new();
        // The named qualification test always computes r8 (~35 min / 84,099
        // DOFs on ts1). The per-commit lane compares its live r1/r2/r4 results
        // with the explicitly identified historical diagnostic below.
        let rungs: &[usize] = if qualification {
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
                    ritz_tol: 1e-13,
                    mass_diagonal_equilibration: true,
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
            let row: [f64; 4] = core::array::from_fn(|i| bending[i].lambda);
            let residuals: [f64; 4] = core::array::from_fn(|i| bending[i].residual);
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
                0.02,
                "residual fraction; theoretical order 4; E=rho=1, nu=0.25",
                started.elapsed().as_secs_f64(),
            );
            values.push(row);
            bounds.push(residuals);
            dofs.push(a.free_dofs.len());
        }
        if !qualification {
            // Recorded diagnostic, NOT a computed or qualified r8 result.
            // ts1 logs/20261016T121744Z-aeb0c36ee98b.log:3189-3191;
            // MagentaBadger's .58 ruling at 2026-10-16 12:54 UTC. Its residual
            // is too loose for qualification; the named long lane replaces it.
            values.push([
                0.19762753568134261,
                0.38992091447195243,
                1.416470926951944,
                2.6511407190005807,
            ]);
            bounds.push([
                3.846899815313211e-6,
                1.6497145952195841e-6,
                8.365556416457706e-8,
                1.9814741730026972e-8,
            ]);
            println!(
                "{{\"suite\":\"irp-s1\",\"stage\":\"historical-r8-comparison\",\"source\":\"ts1:logs/20261016T121744Z-aeb0c36ee98b.log\",\"qualification\":false}}"
            );
        }
        assert!(
            bar_trend(&values, &bounds, qualification).is_ok(),
            "bar trend failed: {:?}; values {values:?}; dofs {dofs:?}",
            bar_trend(&values, &bounds, qualification)
        );
        // Multiplying the actual assembled consistent mass by 1+h² scales
        // each generalized eigenvalue and its bound by 1/(1+h²), exactly:
        // no second expensive eigensolve is needed to plant this defect.
        let planted: Vec<[f64; 4]> = values
            .iter()
            .enumerate()
            .map(|(i, row)| {
                let h = 1.0 / (1usize << i) as f64;
                row.map(|x| x / (1.0 + h * h))
            })
            .collect();
        let planted_bounds: Vec<[f64; 4]> = bounds
            .iter()
            .enumerate()
            .map(|(i, row)| {
                let h = 1.0 / (1usize << i) as f64;
                row.map(|x| x / (1.0 + h * h))
            })
            .collect();
        assert!(
            matches!(
                bar_trend(&planted, &planted_bounds, qualification),
                Err("order trend") | Err("Richardson agreement")
            ),
            "mass*(1+h²) falsifier must break the trend or limit agreement"
        );
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
    });
}

fn bar_trend(
    values: &[[f64; 4]],
    bounds: &[[f64; 4]],
    qualification: bool,
) -> Result<(), &'static str> {
    if values.len() != 4
        || bounds.len() != 4
        || values.iter().flatten().any(|x| !x.is_finite())
        || bounds.iter().flatten().any(|x| !x.is_finite() || *x < 0.0)
    {
        return Err("invalid eigen ladder");
    }
    for mode in 0..4 {
        let mut orders = [0.0; 2];
        let mut limits = [0.0; 2];
        let mut uncertainties = [0.0; 2];
        let mut truncation = 0.0;
        for first in 0..2 {
            let [x, y, z] = core::array::from_fn(|j| values[first + j][mode]);
            let [ex, ey, ez] = core::array::from_fn(|j| bounds[first + j][mode]);
            let d0 = x - y;
            let d1 = y - z;
            if !(d0 > 0.0 && d1 > 0.0 && d0 > d1) {
                return Err("order trend");
            }
            if (qualification || first == 0) && ey + ez > 0.02 * d1 {
                return Err("eigen residual fraction");
            }
            orders[first] = fs_math::det::ln(d0 / d1) / fs_math::det::ln(2.0);
            let correction = d1 * d1 / (d0 - d1);
            limits[first] = z - correction;
            // Bound changes of the extrapolant over the entire residual box:
            // d1 in [d1-ey-ez,d1+ey+ez], denominator in [s-es,s+es].
            let s = d0 - d1;
            let es = ex + 2.0 * ey + ez;
            if s <= es || d1 <= ey + ez {
                return Err("unresolved Richardson denominator");
            }
            let upper = (d1 + ey + ez) * (d1 + ey + ez) / (s - es);
            let lower = (d1 - ey - ez) * (d1 - ey - ez) / (s + es);
            uncertainties[first] = ez + (upper - correction).max(correction - lower);
            if first == 0 {
                // The two pre-asymptotic limits have truncation error in
                // addition to eigensolver uncertainty. Use the discrepancy
                // between observed-order and theoretical-p=4 corrections;
                // this is an estimated envelope, not a certified limit.
                truncation = (correction - d1 / 15.0).abs();
            }
        }
        let allowance = truncation + uncertainties[0] + uncertainties[1];
        println!(
            "{{\"suite\":\"irp-s1\",\"stage\":\"bar-order-trend\",\"mode\":{mode},\"orders\":{orders:?},\"theoretical_order\":4.0,\"limits\":{limits:?},\"residual_uncertainties\":{uncertainties:?},\"truncation_estimate\":{truncation},\"agreement_bound\":{allowance},\"qualification\":{qualification}}}"
        );
        if !(orders[1] > orders[0] && orders[1] >= 3.5 && orders[1] <= 4.3) {
            return Err("order trend");
        }
        if (limits[0] - limits[1]).abs() > allowance {
            return Err("Richardson agreement");
        }
    }
    Ok(())
}
