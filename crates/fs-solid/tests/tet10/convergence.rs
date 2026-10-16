use super::support::*;
use fs_feec::highorder::simplex::duffy_quadrature;
use fs_mms::{LadderSide, OrderGate, RefinementLadder};
use fs_solid::linear3::TetElasticMaterial;
use fs_solid::linear3::tet10::{BoundaryProjection, promote_p1_to_p2, tet10_point};

const AMPLITUDE: [f64; 3] = [1.0, 0.7, -0.4];

/// Manufactured scalar bubble and its exact gradient/Hessian; all three
/// displacement components use this bubble with different amplitudes.
fn bubble(p: [f64; 3]) -> (f64, [f64; 3], [[f64; 3]; 3]) {
    let f = p.map(|x| x * (1.0 - x));
    let d = p.map(|x| 1.0 - 2.0 * x);
    let value = f.iter().product();
    let grad =
        core::array::from_fn(|i| d[i] * (0..3).filter(|&k| k != i).map(|k| f[k]).product::<f64>());
    let h = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            if i == j {
                -2.0 * (0..3).filter(|&k| k != i).map(|k| f[k]).product::<f64>()
            } else {
                d[i] * d[j] * f[3 - i - j]
            }
        })
    });
    (value, grad, h)
}

fn strain(g: [[f64; 3]; 3]) -> [f64; 6] {
    let s = core::f64::consts::FRAC_1_SQRT_2;
    [
        g[0][0],
        g[1][1],
        g[2][2],
        (g[0][1] + g[1][0]) * s,
        (g[1][2] + g[2][1]) * s,
        (g[2][0] + g[0][2]) * s,
    ]
}

fn body(p: [f64; 3], m: &TetElasticMaterial) -> [f64; 3] {
    let (_, _, h) = bubble(p);
    let mut force = [0.0; 3];
    for j in 0..3 {
        let dg = core::array::from_fn(|i| core::array::from_fn(|k| AMPLITUDE[i] * h[k][j]));
        let de = strain(dg);
        let ds = m
            .stiffness_mandel_pa()
            .map(|r| r.iter().zip(de).map(|(a, b)| a * b).sum::<f64>());
        let s = core::f64::consts::FRAC_1_SQRT_2;
        let ds_tensor = [
            [ds[0], ds[3] * s, ds[5] * s],
            [ds[3] * s, ds[1], ds[4] * s],
            [ds[5] * s, ds[4] * s, ds[2]],
        ];
        for i in 0..3 {
            force[i] -= ds_tensor[i][j];
        }
    }
    force
}

pub fn mms_orders() {
    let started = std::time::Instant::now();
    with_cx(|cx| {
        let region = fs_mesh::RegionSpec {
            id: fs_mesh::RegionId(58),
            kind: fs_mesh::RegionKind::Solid,
            seed: [0.31, 0.37, 0.41],
            triangles: fs_mesh::box_triangles(0),
        };
        let base = fs_mesh::volumetricize(
            fs_mesh::UnverifiedPlc::new(
                fs_mesh::box_vertices(0.0, 1.0, 0.0, 1.0, 0.0, 1.0),
                vec![region],
            ),
            fs_mesh::VolumetricPolicy::fixture_default("m"),
            cx,
        )
        .unwrap();
        let law = fs_material::OrthotropicElastic::new(
            [12.0, 5.0, 3.0],
            [0.2, 0.15, 0.1],
            [2.0, 1.5, 1.0],
            1.0,
        )
        .unwrap();
        let anisotropic = TetElasticMaterial::try_new_oriented_orthotropic(
            1.0,
            &law,
            [[0.8, -0.6, 0.0], [0.6, 0.8, 0.0], [0.0, 0.0, 1.0]],
            [1; 32],
            [2; 32],
        )
        .unwrap();
        for (name, m) in [
            ("isotropic", material(10.0, 0.25)),
            ("rotated-orthotropic", anisotropic),
        ] {
            let mut mesh = base.labeled().refine_uniform();
            let mut hs = Vec::new();
            let mut l2s = Vec::new();
            let mut ens = Vec::new();
            let mut dofs = Vec::new();
            for rung in 0..3 {
                let p1 = Mesh {
                    nodes: mesh.positions().to_vec(),
                    tets: mesh.tets().iter().map(|t| t.map(|i| i as usize)).collect(),
                };
                let (nodes, elements) = promote(&p1, cx);
                let a = assembly(&nodes, &elements, &m, &[], cx);
                let rhs = volume_load(&nodes, &elements, |p| body(p, &m));
                let fixed: Vec<_> = nodes
                    .iter()
                    .flat_map(|p| {
                        [if p.iter()
                            .any(|x| x.abs() < 1e-12 || (*x - 1.0).abs() < 1e-12)
                        {
                            Some(0.0)
                        } else {
                            None
                        }; 3]
                    })
                    .collect();
                let u = solve(&a.stiffness, &rhs, &fixed);
                let mut l2 = 0.0;
                let mut energy = 0.0;
                for (e, t) in elements.iter().enumerate() {
                    for (q, (lambda, w)) in duffy_quadrature(5).into_iter().enumerate() {
                        let p = tet10_point(t.map(|i| nodes[i]), lambda, 1e-12, e, q).unwrap();
                        let (exact, gradient, _) = bubble(p.position_m);
                        let mut uh = [0.0; 3];
                        let mut gu = [[0.0; 3]; 3];
                        for a in 0..10 {
                            for i in 0..3 {
                                uh[i] += p.shapes[a] * u[3 * t[a] + i];
                                for k in 0..3 {
                                    gu[i][k] += p.gradients[a][k] * u[3 * t[a] + i];
                                }
                            }
                        }
                        for i in 0..3 {
                            l2 += w * p.determinant_m3 * (uh[i] - AMPLITUDE[i] * exact).powi(2);
                            for k in 0..3 {
                                gu[i][k] -= AMPLITUDE[i] * gradient[k];
                            }
                        }
                        let error = strain(gu);
                        let ce = m
                            .stiffness_mandel_pa()
                            .map(|r| r.iter().zip(error).map(|(a, b)| a * b).sum::<f64>());
                        energy += w
                            * p.determinant_m3
                            * error.iter().zip(ce).map(|(a, b)| a * b).sum::<f64>();
                    }
                }
                let h = 1.0 / (2usize << rung) as f64;
                hs.push(h);
                l2s.push(l2.sqrt());
                ens.push(energy.sqrt());
                dofs.push(fixed.iter().filter(|x| x.is_none()).count());
                log(
                    "g1_mms_orders_isotropic_and_rotated_orthotropic",
                    "G1",
                    name,
                    &[h],
                    &[*dofs.last().unwrap()],
                    &[l2.sqrt(), energy.sqrt()],
                    0.2,
                    "analytic",
                    started.elapsed().as_secs_f64(),
                );
                if rung < 2 {
                    mesh = mesh.refine_uniform();
                }
            }
            for (norm, errors, order) in [("L2", l2s, 3.0), ("energy", ens, 2.0)] {
                let ladder = RefinementLadder::new(hs.clone(), errors).unwrap();
                let verdict = OrderGate { theoretical: order }
                    .check(&format!("tet10-{name}-{norm}"), LadderSide::Primal, &ladder)
                    .unwrap();
                println!("{}", verdict.json_line(true));
            }
        }
    });
}

struct Cylinder;
impl BoundaryProjection for Cylinder {
    fn project(&self, mut p: [f64; 3]) -> Result<[f64; 3], fs_solid::linear3::TetElasticError> {
        let radius = fs_math::det::sqrt(p[0] * p[0] + p[1] * p[1]);
        // This fixture has only a centre and an outer ring: radial cap edges
        // have r <= 0.5; outer-side chord midpoints have r > cos(pi/8).
        if radius > 0.9 {
            p[0] /= radius;
            p[1] /= radius;
        }
        Ok(p)
    }
}

fn cylinder(n: usize) -> Mesh {
    let mut nodes = vec![[0.0, 0.0, 0.0]];
    for j in 0..n {
        let angle = core::f64::consts::TAU * j as f64 / n as f64;
        nodes.push([fs_math::det::cos(angle), fs_math::det::sin(angle), 0.0]);
    }
    let upper = nodes.len();
    for i in 0..upper {
        nodes.push([nodes[i][0], nodes[i][1], 1.0]);
    }
    let mut tets = Vec::new();
    for j in 0..n {
        let mut triangle = [0, j + 1, (j + 1) % n + 1];
        triangle.sort_unstable();
        let [a, b, c] = triangle;
        let [aa, bb, cc] = triangle.map(|i| i + upper);
        tets.extend([[a, b, c, cc], [a, b, bb, cc], [a, aa, bb, cc]]);
    }
    Mesh { nodes, tets }
}

pub fn curved_geometry() {
    let started = std::time::Instant::now();
    with_cx(|cx| {
        let mut hs = Vec::new();
        let mut straight = Vec::new();
        let mut curved = Vec::new();
        for n in [8, 16, 32] {
            let p1 = cylinder(n);
            hs.push(1.0 / n as f64);
            for (projection, errors) in [
                (None, &mut straight),
                (Some(&Cylinder as &dyn BoundaryProjection), &mut curved),
            ] {
                let (nodes, elements) =
                    promote_p1_to_p2(&p1.nodes, &p1.tets, projection, cx).unwrap();
                let mut volume = 0.0;
                for (e, t) in elements.iter().enumerate() {
                    for (q, (lambda, w)) in duffy_quadrature(4).into_iter().enumerate() {
                        volume += w * tet10_point(t.map(|i| nodes[i]), lambda, 1e-12, e, q)
                            .unwrap()
                            .determinant_m3;
                    }
                }
                errors.push((volume - core::f64::consts::PI).abs());
            }
        }
        let a = RefinementLadder::new(hs.clone(), straight.clone()).unwrap();
        let b = RefinementLadder::new(hs.clone(), curved.clone()).unwrap();
        let p = OrderGate { theoretical: 2.0 }
            .check("tet10-straight-cylinder-geometry", LadderSide::Primal, &a)
            .unwrap();
        let q = fs_mms::fit_order(&b);
        assert!(q.observed >= 2.8, "curved geometry order {}", q.observed);
        assert!(curved.iter().zip(&straight).all(|(a, b)| a < b));
        log(
            "g1_curved_boundary_projection_reduces_geometry_error",
            "G1",
            "unit-cylinder-geometry",
            &hs,
            &[0; 3],
            &straight,
            0.2,
            "analytic volume pi; geometry axis",
            started.elapsed().as_secs_f64(),
        );
        log(
            "g1_curved_boundary_projection_reduces_geometry_error",
            "G1",
            "projected-unit-cylinder-geometry",
            &hs,
            &[0; 3],
            &curved,
            2.8,
            "analytic volume pi; geometry axis",
            started.elapsed().as_secs_f64(),
        );
        println!("{}", p.json_line(true));
        println!(
            "{{\"suite\":\"irp-s1\",\"geometry_order_projected\":{},\"minimum\":2.8}}",
            q.observed
        );
    });
}
