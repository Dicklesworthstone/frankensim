//! G0/G1/G3: actual thin-solid equilibrium and independent work derivatives.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::hyper::{FibreInvariants, Hyperelastic, HyperelasticModel, Volumetric};
use fs_solid::{HyperProblem, Mesh2, NewtonSettings, Patch, SolidError};

#[path = "hyper_plane_stress/reference.rs"]
mod reference;

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 71,
                kernel_id: 49,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

fn nh(mu: f64, lambda: f64) -> Hyperelastic {
    Hyperelastic::new(HyperelasticModel::NeoHookean { mu, lambda }, 4.0).unwrap()
}

fn problem<'a>(mesh: &'a Mesh2, material: &'a Hyperelastic) -> HyperProblem<'a> {
    HyperProblem {
        mesh,
        material,
        dirichlet: vec![],
        traction: vec![],
        settings: NewtonSettings {
            tol: 2e-9,
            max_iters: 40,
            load_steps: 5,
            max_backtracks: 30,
        },
    }
}

fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol * (1.0 + b.abs()), "{a:e} != {b:e}");
}

#[test]
fn g1_dead_traction_solve_recovers_finite_lateral_and_thickness_contraction() {
    with_cx(|cx, _| {
        let (mu, lambda, axial): (f64, f64, f64) = (2.0, 40.0, 1.3);
        // Independent uniaxial NH reference: Pyy=Pzz=0 with equal lateral
        // stretches t gives mu(t²-1)+lambda ln(axial*t²)=0.
        let mut lo = 1.0 / axial.sqrt();
        let mut hi = 1.0;
        for _ in 0..70 {
            let t = 0.5 * (lo + hi);
            if mu * (t * t - 1.0) + lambda * (axial * t * t).ln() > 0.0 {
                hi = t;
            } else {
                lo = t;
            }
        }
        let lateral = 0.5 * (lo + hi);
        let nominal = mu * (axial - lateral * lateral / axial);
        let law = nh(mu, lambda);
        let left = |_: f64, y: f64| [0.0, (lateral - 1.0) * y];
        let pull = |_: f64, _: f64| [nominal, 0.0];
        let local = PlaneStressSettings::default();
        for mesh in [
            Mesh2::triangles(2.0, 1.0, 2, 2),
            Mesh2::quads(2.0, 1.0, 2, 2),
        ] {
            let mut p = problem(&mesh, &law);
            p.dirichlet = vec![(Patch::Left, &left)];
            p.traction = vec![(Patch::Right, &pull)];
            let (u, report) = p.solve_plane_stress(local, cx).unwrap();
            assert_eq!(report.histories.len(), p.settings.load_steps);
            assert!(
                report
                    .histories
                    .iter()
                    .all(|h| h.last().unwrap() < &p.settings.tol)
            );
            for (v, [x, y]) in u.iter().zip(&mesh.nodes) {
                close(v[0], (axial - 1.0) * x, 2e-8);
                close(v[1], (lateral - 1.0) * y, 2e-8);
            }
            let flat: Vec<f64> = u.iter().flatten().copied().collect();
            for state in p.plane_stress_states(&flat, local, cx).unwrap() {
                close(state.deformation_gradient[8], lateral, 2e-8);
                assert!(state.deformation_gradient[8] < 0.99);
                close(state.piola[0], nominal, 2e-8);
                close(state.piola[3], 0.0, 2e-8);
                assert!(state.relative_transverse_residual <= local.relative_tolerance);
            }
            let (residual, _) = p
                .plane_stress_residual_and_tangent(&flat, 1.0, local, cx)
                .unwrap();
            let reaction: f64 = mesh
                .patch_nodes(Patch::Left)
                .iter()
                .map(|&i| residual[2 * i])
                .sum();
            close(reaction, -nominal, 2e-8);
        }
    });
}

#[test]
fn g0_skew_elements_energy_force_and_tangent_share_the_relaxed_material() {
    with_cx(|cx, _| {
        let local = PlaneStressSettings::default();
        let fibre = Hyperelastic::new(
            HyperelasticModel::HolzapfelGasserOgden {
                c10: 1.0,
                k1: 1.8,
                k2: 1.1,
                dispersion: 0.05,
                fibres: vec![[0.8, 0.0, 0.6]],
                invariants: FibreInvariants::Full,
                volumetric: Volumetric::Quadratic { kappa: 20.0 },
            },
            4.0,
        )
        .unwrap();
        for law in [nh(2.0, 40.0), fibre] {
            for mesh in [
                Mesh2::triangles(1.0, 0.7, 1, 1),
                Mesh2::mapped_quads(2, 1, &|x, y| [x + 0.15 * x * y + 0.1 * y, y + 0.1 * x]),
            ] {
                let traction = |x: f64, y: f64| [0.1 + x * 0.04, -0.03 * y];
                let mut p = problem(&mesh, &law);
                p.traction = vec![(Patch::Right, &traction)];
                let u: Vec<f64> = mesh
                    .nodes
                    .iter()
                    .flat_map(|[x, y]| [0.14 * x + 0.03 * x * y, -0.04 * y + 0.02 * x * y])
                    .collect();
                let (force, tangent) = p
                    .plane_stress_residual_and_tangent(&u, 0.7, local, cx)
                    .unwrap();
                let h = 2e-6;
                for col in 0..u.len() {
                    let mut plus = u.clone();
                    let mut minus = u.clone();
                    plus[col] += h;
                    minus[col] -= h;
                    let ep = p
                        .plane_stress_potential_energy(&plus, 0.7, local, cx)
                        .unwrap();
                    let em = p
                        .plane_stress_potential_energy(&minus, 0.7, local, cx)
                        .unwrap();
                    close((ep - em) / (2.0 * h), force[col], 3e-7);
                    let rp = p
                        .plane_stress_residual_and_tangent(&plus, 0.7, local, cx)
                        .unwrap()
                        .0;
                    let rm = p
                        .plane_stress_residual_and_tangent(&minus, 0.7, local, cx)
                        .unwrap()
                        .0;
                    for row in 0..u.len() {
                        close((rp[row] - rm[row]) / (2.0 * h), tangent.get(row, col), 5e-6);
                        close(tangent.get(row, col), tangent.get(col, row), 5e-11);
                    }
                }
            }
        }
    });
}

#[test]
fn g3_rigid_motion_preserves_energy_and_rotates_internal_forces() {
    with_cx(|cx, _| {
        let local = PlaneStressSettings::default();
        let mesh = Mesh2::mapped_quads(2, 1, &|x, y| [x + 0.12 * y, y + 0.07 * x * y]);
        let law = nh(2.0, 40.0);
        let p = problem(&mesh, &law);
        let u: Vec<f64> = mesh
            .nodes
            .iter()
            .flat_map(|[x, y]| [0.18 * x + 0.02 * x * y, -0.04 * y])
            .collect();
        let (s, c) = 0.64_f64.sin_cos();
        let rotated: Vec<f64> = mesh
            .nodes
            .iter()
            .enumerate()
            .flat_map(|(i, [x, y])| {
                let px = x + u[2 * i];
                let py = y + u[2 * i + 1];
                [c * px - s * py - x + 2.0, s * px + c * py - y - 1.0]
            })
            .collect();
        close(
            p.plane_stress_potential_energy(&u, 0.0, local, cx).unwrap(),
            p.plane_stress_potential_energy(&rotated, 0.0, local, cx)
                .unwrap(),
            2e-11,
        );
        let (force, _) = p
            .plane_stress_residual_and_tangent(&u, 0.0, local, cx)
            .unwrap();
        let (rf, _) = p
            .plane_stress_residual_and_tangent(&rotated, 0.0, local, cx)
            .unwrap();
        for i in 0..mesh.nodes.len() {
            close(rf[2 * i], c * force[2 * i] - s * force[2 * i + 1], 2e-10);
            close(
                rf[2 * i + 1],
                s * force[2 * i] + c * force[2 * i + 1],
                2e-10,
            );
        }
    });
}

#[test]
fn g0_rejects_folded_domains_bad_controls_and_cancellation_without_partial_output() {
    with_cx(|cx, gate| {
        let local = PlaneStressSettings::default();
        let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
        let law = nh(2.0, 40.0);
        let p = problem(&mesh, &law);
        let mut folded = vec![0.0; 8];
        folded[6] = -1.5;
        assert!(matches!(
            p.plane_stress_states(&folded, local, cx),
            Err(SolidError::MaterialRefused { .. })
        ));
        assert!(matches!(
            p.plane_stress_states(&[0.0; 6], local, cx),
            Err(SolidError::InvalidInput { .. })
        ));
        let mut limited = law.clone();
        limited.stretch_limit = 1.1;
        let stretch: Vec<f64> = mesh
            .nodes
            .iter()
            .flat_map(|[x, _]| [0.3 * x, 0.0])
            .collect();
        assert!(matches!(
            problem(&mesh, &limited).plane_stress_states(&stretch, local, cx),
            Err(SolidError::MaterialRefused { .. })
        ));
        let mut bad = problem(&mesh, &law);
        bad.settings.load_steps = 0;
        assert!(matches!(
            bad.solve_plane_stress(local, cx),
            Err(SolidError::InvalidInput { .. })
        ));
        let mut invalid_mesh = mesh.clone();
        invalid_mesh.elems[0][0] = usize::MAX;
        assert!(matches!(
            problem(&invalid_mesh, &law).solve_plane_stress(local, cx),
            Err(SolidError::InvalidInput { .. })
        ));
        let cancel_during_load = |_: f64, _: f64| {
            gate.request();
            [0.0, 0.0]
        };
        let mut cancelled = problem(&mesh, &law);
        cancelled.traction = vec![(Patch::Right, &cancel_during_load)];
        assert!(matches!(
            cancelled.plane_stress_potential_energy(&[0.0; 8], 1.0, local, cx),
            Err(SolidError::Cancelled)
        ));
        assert!(matches!(
            p.solve_plane_stress(local, cx),
            Err(SolidError::Cancelled)
        ));
    });
}

#[test]
fn g1_exact_incompressible_traction_solve_preserves_volume_and_ignores_bulk_penalties() {
    with_cx(|cx, _| {
        let (mu, axial): (f64, f64) = (2.0, 1.4);
        // Independent incompressible NH uniaxial solution: t = axial^-1/2,
        // Pxx = mu (axial - axial^-2), W = mu/2 (axial² + 2/axial - 3).
        let lateral = 1.0 / axial.sqrt();
        let nominal = mu * (axial - 1.0 / (axial * axial));
        let left = |_: f64, y: f64| [0.0, (lateral - 1.0) * y];
        let pull = |_: f64, _: f64| [nominal, 0.0];
        let local = PlaneStressSettings::default();
        for mesh in [
            Mesh2::triangles(2.0, 1.0, 2, 2),
            Mesh2::quads(2.0, 1.0, 2, 2),
        ] {
            let mut first_displacement = None;
            for lambda in [0.0, 40.0, 1e300] {
                let law = nh(mu, lambda);
                let mut p = problem(&mesh, &law);
                p.dirichlet = vec![(Patch::Left, &left)];
                p.traction = vec![(Patch::Right, &pull)];
                let (u, report) = p.solve_incompressible_plane_stress(local, cx).unwrap();
                assert_eq!(report.histories.len(), p.settings.load_steps);
                for (v, [x, y]) in u.iter().zip(&mesh.nodes) {
                    close(v[0], (axial - 1.0) * x, 2e-8);
                    close(v[1], (lateral - 1.0) * y, 2e-8);
                }
                if let Some(first) = &first_displacement {
                    assert_eq!(&u, first, "bulk coefficient entered the exact constraint");
                } else {
                    first_displacement = Some(u.clone());
                }
                let flat: Vec<f64> = u.iter().flatten().copied().collect();
                for state in p
                    .incompressible_plane_stress_states(&flat, local, cx)
                    .unwrap()
                {
                    let f = state.deformation_gradient;
                    close((f[0] * f[4] - f[1] * f[3]) * f[8], 1.0, 4e-15);
                    close(f[8], lateral, 2e-8);
                    close(state.piola[0], nominal, 2e-8);
                    close(state.piola[3], 0.0, 2e-8);
                    // The original NH remainder includes -mu ln J; this fixes
                    // the multiplier gauge used by pressure recovery.
                    close(state.pressure_pa, mu * (lateral * lateral - 1.0), 2e-8);
                    for i in [2, 5, 8] {
                        close(state.full_piola[i], 0.0, 2e-10);
                    }
                }
                let (residual, _) = p
                    .incompressible_plane_stress_residual_and_tangent(&flat, 1.0, local, cx)
                    .unwrap();
                let reaction: f64 = mesh
                    .patch_nodes(Patch::Left)
                    .iter()
                    .map(|&i| residual[2 * i])
                    .sum();
                close(reaction, -nominal, 2e-8);
                let storage = p
                    .incompressible_plane_stress_potential_energy(&flat, 1.0, local, cx)
                    .unwrap()
                    + nominal * (axial - 1.0) * 2.0;
                close(storage, mu * (axial * axial + 2.0 / axial - 3.0), 2e-8);
            }
        }
    });
}

#[test]
fn g0_incompressible_skew_element_potential_force_and_tangent_include_pressure_relaxation() {
    with_cx(|cx, _| {
        let local = PlaneStressSettings::default();
        let fibre = Hyperelastic::new(
            HyperelasticModel::HolzapfelGasserOgden {
                c10: 1.0,
                k1: 1.8,
                k2: 1.1,
                dispersion: 0.05,
                fibres: vec![[0.8, 0.0, 0.6]],
                invariants: FibreInvariants::Full,
                volumetric: Volumetric::Quadratic { kappa: 1e300 },
            },
            4.0,
        )
        .unwrap();
        for law in [nh(2.0, 40.0), fibre] {
            for mesh in [
                Mesh2::triangles(1.0, 0.7, 1, 1),
                Mesh2::mapped_quads(2, 1, &|x, y| [x + 0.15 * x * y + 0.1 * y, y + 0.1 * x]),
            ] {
                let traction = |x: f64, y: f64| [0.1 + x * 0.04, -0.03 * y];
                let mut p = problem(&mesh, &law);
                p.traction = vec![(Patch::Right, &traction)];
                let u: Vec<f64> = mesh
                    .nodes
                    .iter()
                    .flat_map(|[x, y]| [0.14 * x + 0.03 * x * y, -0.04 * y + 0.02 * x * y])
                    .collect();
                let (force, tangent) = p
                    .incompressible_plane_stress_residual_and_tangent(&u, 0.7, local, cx)
                    .unwrap();
                let h = 2e-6;
                for col in 0..u.len() {
                    let mut plus = u.clone();
                    let mut minus = u.clone();
                    plus[col] += h;
                    minus[col] -= h;
                    let ep = p
                        .incompressible_plane_stress_potential_energy(&plus, 0.7, local, cx)
                        .unwrap();
                    let em = p
                        .incompressible_plane_stress_potential_energy(&minus, 0.7, local, cx)
                        .unwrap();
                    close((ep - em) / (2.0 * h), force[col], 3e-7);
                    let rp = p
                        .incompressible_plane_stress_residual_and_tangent(&plus, 0.7, local, cx)
                        .unwrap()
                        .0;
                    let rm = p
                        .incompressible_plane_stress_residual_and_tangent(&minus, 0.7, local, cx)
                        .unwrap()
                        .0;
                    for row in 0..u.len() {
                        close((rp[row] - rm[row]) / (2.0 * h), tangent.get(row, col), 5e-6);
                        close(tangent.get(row, col), tangent.get(col, row), 5e-11);
                    }
                }
            }
        }
    });
}

#[test]
fn g0_incompressible_thickness_domain_and_cancellation_remain_enforced() {
    with_cx(|cx, gate| {
        let local = PlaneStressSettings::default();
        let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
        let mut law = nh(2.0, 40.0);
        law.stretch_limit = 2.0;
        let p = problem(&mesh, &law);
        // The in-plane stretches are allowed, but volume conservation would
        // require F22 = 1/4, outside the unchanged lower bound 1/2.
        let stretch: Vec<f64> = mesh.nodes.iter().flat_map(|[x, y]| [*x, *y]).collect();
        assert!(matches!(
            p.incompressible_plane_stress_states(&stretch, local, cx),
            Err(SolidError::MaterialRefused { .. })
        ));
        gate.request();
        assert!(matches!(
            p.solve_incompressible_plane_stress(local, cx),
            Err(SolidError::Cancelled)
        ));
    });
}

fn motion_affine_field(mesh: &Mesh2, gradient: [f64; 4], translation: [f64; 2]) -> Vec<f64> {
    mesh.nodes
        .iter()
        .flat_map(|[x, y]| {
            [
                translation[0] + gradient[0] * x + gradient[1] * y,
                translation[1] + gradient[2] * x + gradient[3] * y,
            ]
        })
        .collect()
}

fn motion_product(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| {
        (0..3).map(|k| a[3 * (ij / 3) + k] * b[3 * k + ij % 3]).sum()
    })
}

#[test]
fn g1_held_plane_motion_has_physical_relaxed_thickness_rate() {
    with_cx(|cx, _| {
        let (mu, lambda, axial): (f64, f64, f64) = (2.0, 40.0, 1.3);
        // True uniaxial NH path: both lateral stresses vanish, with
        // mu*(t²-1)+lambda*ln(axial*t²)=0. Differentiate this scalar equation.
        let mut lo = 1.0 / axial.sqrt();
        let mut hi = 1.0;
        for _ in 0..70 {
            let t = 0.5 * (lo + hi);
            if mu * (t * t - 1.0) + lambda * (axial * t * t).ln() > 0.0 {
                hi = t;
            } else {
                lo = t;
            }
        }
        let lateral = 0.5 * (lo + hi);
        let axial_rate = 0.23;
        let lateral_rate =
            -lambda * lateral * axial_rate / (2.0 * axial * (mu * lateral * lateral + lambda));
        let nominal = mu * (axial - lateral * lateral / axial);
        let law = nh(mu, lambda);
        for mesh in [
            Mesh2::triangles(2.0, 1.0, 2, 2),
            Mesh2::quads(2.0, 1.0, 2, 2),
        ] {
            let p = problem(&mesh, &law);
            let u = motion_affine_field(&mesh, [axial - 1.0, 0.0, 0.0, lateral - 1.0], [0.0; 2]);
            let v = motion_affine_field(&mesh, [axial_rate, 0.0, 0.0, lateral_rate], [0.0; 2]);
            let fields = p
                .plane_stress_held_material_motion(&u, &v, PlaneStressSettings::default(), cx)
                .unwrap();
            close(fields.iter().map(|f| f.reference_area_m2).sum(), 2.0, 1e-12);
            for field in &fields {
                assert!(field.deformation_gradient_rate_per_s[8] < -0.01);
                let stretch = [axial, lateral, lateral];
                let rates = [axial_rate, lateral_rate, lateral_rate];
                for i in 0..3 {
                    for j in 0..3 {
                        close(
                            field.deformation_gradient_rate_per_s[3 * i + j],
                            if i == j { rates[i] } else { 0.0 },
                            2e-10,
                        );
                        close(
                            field.green_lagrange_strain_rate_per_s[3 * i + j],
                            if i == j { stretch[i] * rates[i] } else { 0.0 },
                            2e-10,
                        );
                        close(
                            field.spatial_strain_rate_per_s[3 * i + j],
                            if i == j { rates[i] / stretch[i] } else { 0.0 },
                            2e-10,
                        );
                    }
                }
                close(
                    field.jacobian_rate_per_s,
                    lateral * lateral * axial_rate + 2.0 * axial * lateral * lateral_rate,
                    2e-10,
                );
                close(field.reference_stress_power_density_w_m3, nominal * axial_rate, 2e-9);
                close(
                    field.reference_stress_power_density_w_m3,
                    field.jacobian * field.spatial_stress_power_density_w_m3,
                    2e-12,
                );
            }
            // The same in-plane motion under plane strain has no thickness rate.
            for field in p.plane_strain_held_material_motion(&u, &v, cx).unwrap() {
                assert_eq!(field.deformation_gradient[8], 1.0);
                assert_eq!(field.deformation_gradient_rate_per_s[8], 0.0);
                assert_eq!(field.spatial_strain_rate_per_s[8], 0.0);
                close(
                    field.reference_stress_power_density_w_m3,
                    field.jacobian * field.spatial_stress_power_density_w_m3,
                    2e-12,
                );
            }
        }
    });
}

#[test]
fn g3_held_planar_rates_are_objective_under_moving_rigid_frames() {
    with_cx(|cx, _| {
        let law = nh(2.0, 3.0);
        let angle = 0.63_f64;
        let (s, c) = angle.sin_cos();
        let q = [c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0];
        let qt = [c, s, 0.0, -s, c, 0.0, 0.0, 0.0, 1.0];
        let spin = [0.0, -0.7, 0.0, 0.7, 0.0, 0.0, 0.0, 0.0, 0.0];
        let f = [1.15, 0.12, 0.0, -0.03, 0.91, 0.0, 0.0, 0.0, 1.0];
        let rate = [0.08, 0.05, 0.0, -0.02, 0.04, 0.0, 0.0, 0.0, 0.0];
        let rotated_f = motion_product(&q, &f);
        let rotated_material_rate = motion_product(&q, &rate);
        let rigid_rate = motion_product(&spin, &rotated_f);
        let rotated_rate: [f64; 9] =
            core::array::from_fn(|i| rotated_material_rate[i] + rigid_rate[i]);
        let in_plane = |a: &[f64; 9]| [a[0], a[1], a[3], a[4]];
        let displacement = |mesh: &Mesh2, a: &[f64; 9], translation| {
            motion_affine_field(mesh, [a[0] - 1.0, a[1], a[3], a[4] - 1.0], translation)
        };
        for mesh in [Mesh2::triangles(1.0, 1.0, 1, 1), Mesh2::quads(1.0, 1.0, 1, 1)] {
            let p = problem(&mesh, &law);
            for relaxed in [false, true] {
                let observe = |u: &[f64], v: &[f64]| {
                    if relaxed {
                        p.plane_stress_held_material_motion(u, v, PlaneStressSettings::default(), cx)
                    } else {
                        p.plane_strain_held_material_motion(u, v, cx)
                    }
                    .unwrap()
                };
                let base = observe(
                    &displacement(&mesh, &f, [0.0; 2]),
                    &motion_affine_field(&mesh, in_plane(&rate), [0.0; 2]),
                );
                let rotated = observe(
                    &displacement(&mesh, &rotated_f, [3.0, -7.0]),
                    &motion_affine_field(&mesh, in_plane(&rotated_rate), [13.0, -17.0]),
                );
                for (before, after) in base.iter().zip(rotated) {
                    let expected_d = motion_product(
                        &motion_product(&q, &before.spatial_strain_rate_per_s),
                        &qt,
                    );
                    for ij in 0..9 {
                        close(
                            after.green_lagrange_strain_rate_per_s[ij],
                            before.green_lagrange_strain_rate_per_s[ij],
                            3e-10,
                        );
                        close(after.spatial_strain_rate_per_s[ij], expected_d[ij], 3e-10);
                    }
                    close(after.jacobian_rate_per_s, before.jacobian_rate_per_s, 3e-10);
                    close(
                        after.reference_stress_power_density_w_m3,
                        before.reference_stress_power_density_w_m3,
                        3e-10,
                    );
                }
                for rigid in observe(
                    &displacement(&mesh, &rotated_f, [3.0, -7.0]),
                    &motion_affine_field(&mesh, in_plane(&rigid_rate), [13.0, -17.0]),
                ) {
                    for value in rigid.spatial_strain_rate_per_s
                        .iter().chain(&rigid.green_lagrange_strain_rate_per_s) {
                        close(*value, 0.0, 3e-10);
                    }
                    close(rigid.reference_stress_power_density_w_m3, 0.0, 3e-10);
                }
            }
        }
    });
}

#[test]
fn g1_anisotropic_plane_motion_retains_relaxed_transverse_shear_rate() {
    with_cx(|cx, _| {
        let law = Hyperelastic::new(
            HyperelasticModel::HolzapfelGasserOgden {
                c10: 1.0,
                k1: 1.8,
                k2: 1.1,
                dispersion: 0.05,
                fibres: vec![[0.8, 0.0, 0.6]],
                invariants: FibreInvariants::Full,
                volumetric: Volumetric::Quadratic { kappa: 20.0 },
            },
            4.0,
        ).unwrap();
        let local = PlaneStressSettings::default();
        for mesh in [Mesh2::triangles(1.0, 0.7, 1, 1), Mesh2::quads(1.0, 0.7, 1, 1)] {
            let p = problem(&mesh, &law);
            let u = motion_affine_field(&mesh, [0.14, 0.03, -0.01, 0.08], [0.0; 2]);
            let v = motion_affine_field(&mesh, [0.11, -0.02, 0.03, -0.04], [0.0; 2]);
            let fields = p.plane_stress_held_material_motion(&u, &v, local, cx).unwrap();
            let step = 2e-6;
            let plus: Vec<_> = u.iter().zip(&v).map(|(u, v)| u + step * v).collect();
            let minus: Vec<_> = u.iter().zip(&v).map(|(u, v)| u - step * v).collect();
            let plus = p.plane_stress_states(&plus, local, cx).unwrap();
            let minus = p.plane_stress_states(&minus, local, cx).unwrap();
            assert!(fields.iter().any(|field| {
                field.deformation_gradient_rate_per_s[2].abs()
                    + field.deformation_gradient_rate_per_s[5].abs() > 1e-8
            }));
            for ((field, plus), minus) in fields.iter().zip(plus).zip(minus) {
                for ij in 0..9 {
                    close(
                        field.deformation_gradient_rate_per_s[ij],
                        (plus.deformation_gradient[ij] - minus.deformation_gradient[ij]) / (2.0 * step),
                        3e-6,
                    );
                }
                close(
                    field.reference_stress_power_density_w_m3,
                    field.jacobian * field.spatial_stress_power_density_w_m3,
                    2e-11,
                );
            }
        }
    });
}

#[test]
fn g0_held_plane_motion_preserves_tiny_rates_and_refuses_invalid_motion() {
    with_cx(|cx, gate| {
        let law = nh(2.0, 3.0);
        let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
        let p = problem(&mesh, &law);
        let u = vec![0.0; 2 * mesh.nodes.len()];
        let tiny = 1e-30;
        let v = motion_affine_field(&mesh, [tiny, 0.0, 0.0, 0.0], [0.0; 2]);
        let local = PlaneStressSettings::default();
        let fields = p.plane_stress_held_material_motion(&u, &v, local, cx).unwrap();
        for field in fields {
            assert!((field.deformation_gradient_rate_per_s[0] / tiny - 1.0).abs() < 1e-13);
            assert!((field.deformation_gradient_rate_per_s[8] / tiny + 3.0 / 7.0).abs() < 1e-13);
        }
        assert!(p.plane_stress_held_material_motion(&u, &v[..v.len() - 1], local, cx).is_err());
        assert!(p.plane_strain_held_material_motion(&u, &[f64::NAN; 8], cx).is_err());
        let inverted = motion_affine_field(&mesh, [-2.0, 0.0, 0.0, 0.0], [0.0; 2]);
        assert!(p.plane_stress_held_material_motion(&inverted, &v, local, cx).is_err());
        assert!(p.plane_strain_held_material_motion(&inverted, &v, cx).is_err());
        gate.request();
        assert!(matches!(
            p.plane_stress_held_material_motion(&u, &v, local, cx),
            Err(SolidError::Cancelled)
        ));
    });
}
