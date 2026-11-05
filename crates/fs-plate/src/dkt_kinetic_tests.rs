use super::*;
use crate::{EdgeSupport, PlateMesh, PlateSection};

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() < 3e-12 * scale.max(actual.abs()).max(expected.abs()),
        "{actual} != {expected}"
    );
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn product(m: &[f64; 81], a: &[f64; 9], b: &[f64; 9]) -> f64 {
    (0..9).map(|r| a[r] * dot(&m[9 * r..9 * r + 9], b)).sum()
}

fn monomial(p: [usize; 2], x: &[f64; 3], y: &[f64; 3]) -> [f64; 9] {
    let at = |n: usize, a: usize, b: usize| x[n].powi(a as i32) * y[n].powi(b as i32);
    core::array::from_fn(|i| match i % 3 {
        0 => at(i / 3, p[0], p[1]),
        1 if p[0] > 0 => p[0] as f64 * at(i / 3, p[0] - 1, p[1]),
        2 if p[1] > 0 => p[1] as f64 * at(i / 3, p[0], p[1] - 1),
        _ => 0.0,
    })
}

fn integral(x: usize, y: usize) -> f64 {
    let factorial = |n: usize| (1..=n).fold(1.0, |p, i| p * i as f64);
    factorial(x) * factorial(y) / factorial(x + y + 2)
}

#[test]
fn rigid_motion_and_all_quadratic_fields_have_exact_physical_rotary_energy() {
    let x = [0.0, 1.0, 0.0];
    let y = [0.0, 0.0, 1.0];
    let powers = [[0, 0], [1, 0], [0, 1], [2, 0], [1, 1], [0, 2]];
    for h in [0.001, 0.009, 0.07] {
        let inertia = 460.0 * h * h * h / 12.0;
        let mass = rotary_mass(&x, &y, 0.5 * inertia).unwrap();
        for a in powers {
            for b in powers {
                let mut expected = 0.0;
                if a[0] * b[0] != 0 {
                    expected += (a[0] * b[0]) as f64 * integral(a[0] + b[0] - 2, a[1] + b[1]);
                }
                if a[1] * b[1] != 0 {
                    expected += (a[1] * b[1]) as f64 * integral(a[0] + b[0], a[1] + b[1] - 2);
                }
                close(
                    product(&mass, &monomial(a, &x, &y), &monomial(b, &x, &y)),
                    inertia * expected,
                    inertia,
                );
            }
        }
    }
}

// Independent geometric DKT construction: at an edge midpoint its tangential
// rotation is the cubic Hermite edge derivative, while the normal rotation is
// the linear mean of the endpoint rotations. Interpolate those six values.
// This oracle does not call the production Batoz geometry/rotation tables.
fn edge_rotation(x: &[f64; 3], y: &[f64; 3], l: [f64; 3], q: &[f64; 9]) -> [f64; 2] {
    let mut beta = [0.0; 2];
    for i in 0..3 {
        let shape = l[i] * (2.0 * l[i] - 1.0);
        beta[0] -= shape * q[3 * i + 1];
        beta[1] -= shape * q[3 * i + 2];
    }
    for (i, j) in [(1, 2), (2, 0), (0, 1)] {
        let length = (x[j] - x[i]).hypot(y[j] - y[i]);
        let t = [(x[j] - x[i]) / length, (y[j] - y[i]) / length];
        let n = [-t[1], t[0]];
        let sum = [q[3 * i + 1] + q[3 * j + 1], q[3 * i + 2] + q[3 * j + 2]];
        let tangent = 1.5 * (q[3 * j] - q[3 * i]) / length - 0.25 * dot(&sum, &t);
        let normal = 0.5 * dot(&sum, &n);
        for axis in 0..2 {
            beta[axis] -= 4.0 * l[i] * l[j] * (tangent * t[axis] + normal * n[axis]);
        }
    }
    beta
}

#[test]
fn exact_quartic_integral_matches_independent_edge_constraints_and_duffy_quadrature() {
    for (x, y) in [
        ([0.0, 1.0, 0.0], [0.0, 0.0, 1.0]),
        ([0.2, 1.4, -0.3], [-0.1, 0.8, 1.2]),
    ] {
        let scale = 0.004;
        let mass = rotary_mass(&x, &y, scale).unwrap();
        let root = (3.0_f64 / 5.0).sqrt();
        let points = [0.5 * (1.0 - root), 0.5, 0.5 * (1.0 + root)];
        let weights = [5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0];
        let mut oracle = [0.0; 81];
        for i in 0..3 {
            for j in 0..3 {
                let u = points[i];
                let v = points[j];
                let l = [(1.0 - u) * (1.0 - v), u, (1.0 - u) * v];
                let beta: [[f64; 2]; 9] = core::array::from_fn(|dof| {
                    let mut q = [0.0; 9];
                    q[dof] = 1.0;
                    edge_rotation(&x, &y, l, &q)
                });
                let weight = scale * 2.0 * (1.0 - u) * weights[i] * weights[j];
                for r in 0..9 {
                    for c in 0..9 {
                        oracle[9 * r + c] += weight * dot(&beta[r], &beta[c]);
                    }
                }
            }
        }
        let perm = [1, 2, 0];
        let permuted = rotary_mass(&perm.map(|i| x[i]), &perm.map(|i| y[i]), scale).unwrap();
        for r in 0..9 {
            for c in 0..9 {
                assert_eq!(mass[9 * r + c], mass[9 * c + r]);
                close(mass[9 * r + c], oracle[9 * r + c], scale);
                close(
                    permuted[9 * r + c],
                    mass[9 * (3 * perm[r / 3] + r % 3) + 3 * perm[c / 3] + c % 3],
                    scale,
                );
            }
        }
    }
}

fn chart() -> PlateChart {
    let mesh = PlateMesh::from_unstructured(
        vec![(0.0, 0.0), (1.0, 0.2), (1.1, 0.9), (-0.1, 0.7)],
        vec![[0, 1, 2], [0, 2, 3]],
    )
    .unwrap();
    let a = PlateSection::isotropic(12e9, 0.25, 0.009, 460.0).unwrap();
    let b = PlateSection::isotropic(8e9, 0.2, 0.006, 530.0).unwrap();
    PlateChart::with_boundary_and_regions(mesh, a, vec![], vec![])
        .unwrap()
        .with_element_sections(vec![a, b])
        .unwrap()
}

fn options(support: EdgeSupport) -> AssemblyOptions {
    AssemblyOptions {
        support,
        pretension: 200.0,
    }
}

fn assert_positive_mass(model: &PlateModel) {
    let n = model.free;
    let mut factor = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..=i {
            close(
                model.m.get(i, j),
                model.m.get(j, i),
                (model.m.get(i, i) * model.m.get(j, j)).sqrt(),
            );
            let value = model.m.get(i, j)
                - (0..j)
                    .map(|k| factor[n * i + k] * factor[n * j + k])
                    .sum::<f64>();
            factor[n * i + j] = if i == j {
                assert!(value > 0.0, "nonpositive kinetic pivot {i}: {value}");
                value.sqrt()
            } else {
                value / factor[n * j + j]
            };
        }
    }
}

#[test]
fn live_assembly_preserves_legacy_pencils_stiffeners_and_supports_for_each_kinetic_choice() {
    let chart = chart();
    let opts = options(EdgeSupport::Clamped);
    let beam = Stiffener {
        nodes: vec![0, 2],
        e: 12e9,
        g: 0.8e9,
        area: 0.0004,
        inertia: 3e-8,
        torsion: 2e-8,
        eccentricity: 0.022,
        density: 500.0,
    };
    for panel in [
        TransverseMass::Lumped,
        TransverseMass::Linear,
        TransverseMass::EdgeCubic,
    ] {
        let bare = chart
            .assemble_with_kinetic_mass(&[], &opts, panel, RotaryMass::Dkt, StiffenerMass::Lumped)
            .unwrap();
        let old_bare = chart
            .assemble_with_mass(&[], &opts, panel, StiffenerMass::Lumped)
            .unwrap();
        assert_positive_mass(&bare);
        let rigid: Vec<_> = chart
            .mesh
            .nodes
            .iter()
            .flat_map(|_| [1.0, 0.0, 0.0])
            .collect();
        let mut mr = vec![0.0; bare.free];
        bare.m.spmv(&rigid, &mut mr);
        let expected_mass: f64 = chart
            .mesh
            .tris
            .iter()
            .enumerate()
            .map(|(i, tri)| {
                let p = tri.map(|n| chart.mesh.nodes[n]);
                let area = 0.5
                    * ((p[1].0 - p[0].0) * (p[2].1 - p[0].1)
                        - (p[2].0 - p[0].0) * (p[1].1 - p[0].1));
                let s = chart.section_field().at(i);
                area * s.density * s.thickness
            })
            .sum();
        close(dot(&rigid, &mr), expected_mass, expected_mass);
        for beam_mass in [
            StiffenerMass::Lumped,
            StiffenerMass::ConsistentHermite,
            StiffenerMass::ConsistentEccentric,
        ] {
            let old = chart
                .assemble_with_mass(std::slice::from_ref(&beam), &opts, panel, beam_mass)
                .unwrap();
            let legacy = chart
                .assemble_with_kinetic_mass(
                    std::slice::from_ref(&beam),
                    &opts,
                    panel,
                    RotaryMass::Lumped,
                    beam_mass,
                )
                .unwrap();
            let full = chart
                .assemble_with_kinetic_mass(
                    std::slice::from_ref(&beam),
                    &opts,
                    panel,
                    RotaryMass::Dkt,
                    beam_mass,
                )
                .unwrap();
            assert_eq!(old.k, legacy.k);
            assert_eq!(old.m, legacy.m);
            assert_eq!(old.k, full.k);
            assert_eq!(old.dof_map, full.dof_map);
            assert_ne!(old.m, full.m);
            assert_eq!(
                full.m,
                chart
                    .assemble_with_kinetic_mass(
                        std::slice::from_ref(&beam),
                        &opts,
                        panel,
                        RotaryMass::Dkt,
                        beam_mass
                    )
                    .unwrap()
                    .m
            );
            assert_positive_mass(&full);
            for r in 0..full.free {
                for c in 0..full.free {
                    close(
                        full.m.get(r, c) - bare.m.get(r, c),
                        old.m.get(r, c) - old_bare.m.get(r, c),
                        0.2,
                    );
                }
            }
            for support in [EdgeSupport::Clamped, EdgeSupport::SimplySupported] {
                let mut fixed = chart.clone();
                fixed.boundary_nodes = vec![0, 3];
                let reduced = fixed
                    .assemble_with_kinetic_mass(
                        std::slice::from_ref(&beam),
                        &options(support),
                        panel,
                        RotaryMass::Dkt,
                        beam_mass,
                    )
                    .unwrap();
                assert_positive_mass(&reduced);
                for (r, rr) in reduced.dof_map.iter().enumerate() {
                    for (c, cc) in reduced.dof_map.iter().enumerate() {
                        if let (Some(rr), Some(cc)) = (rr, cc) {
                            assert_eq!(reduced.m.get(*rr, *cc), full.m.get(r, c));
                            assert_eq!(reduced.k.get(*rr, *cc), full.k.get(r, c));
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn bending_wavelength_dispersion_converges_to_rayleigh_kirchhoff_continuum() {
    let section = PlateSection::isotropic(12e9, 0.25, 0.03, 460.0).unwrap();
    let opts = AssemblyOptions {
        support: EdgeSupport::SimplySupported,
        pretension: 0.0,
    };
    // P1 is the live acoustic-radiator field; cubic is the higher-order option.
    // Each keeps its own transverse interpolation error and the same rotations.
    for panel in [TransverseMass::Linear, TransverseMass::EdgeCubic] {
        let mut coarse_shifts = Vec::new();
        for mode_x in [1, 3] {
            let kx = core::f64::consts::PI * f64::from(mode_x);
            let ky = core::f64::consts::PI / 0.8;
            let k2 = kx * kx + ky * ky;
            let rotary_ratio = section.thickness.powi(2) * k2 / 12.0;
            let lambda = section.d[0] * k2 * k2
                / (section.density * section.thickness * (1.0 + rotary_ratio));
            let mut errors = Vec::new();
            let mut shifts = Vec::new();
            for n in [6, 12, 24] {
                let chart =
                    PlateChart::from_mesh(PlateMesh::rectangle(1.0, 0.8, n, n), section).unwrap();
                let model = chart
                    .assemble_with_kinetic_mass(
                        &[],
                        &opts,
                        panel,
                        RotaryMass::Dkt,
                        StiffenerMass::Lumped,
                    )
                    .unwrap();
                let mut q = vec![0.0; model.free];
                for (node, &(x, y)) in chart.mesh.nodes.iter().enumerate() {
                    let (sx, cx) = (kx * x).sin_cos();
                    let (sy, cy) = (ky * y).sin_cos();
                    for (axis, value) in [sx * sy, kx * cx * sy, ky * sx * cy]
                        .into_iter()
                        .enumerate()
                    {
                        if let Some(i) = model.dof_map[3 * node + axis] {
                            q[i] = value;
                        }
                    }
                }
                let mut action = vec![0.0; model.free];
                model.k.spmv(&q, &mut action);
                let stiffness = dot(&q, &action);
                model.m.spmv(&q, &mut action);
                let inertia = dot(&q, &action);
                let omega = (stiffness / inertia).sqrt();
                errors.push((omega / lambda.sqrt() - 1.0).abs());
                let legacy = chart
                    .assemble_with_mass(&[], &opts, panel, StiffenerMass::Lumped)
                    .unwrap();
                legacy.m.spmv(&q, &mut action);
                let legacy_inertia = dot(&q, &action);
                shifts.push((legacy_inertia / inertia).sqrt() - 1.0);
                // The exact continuum rotary fraction grows with k^2. The actual
                // assembled inertia must converge to it, retaining physical units.
                let continuum_inertia =
                    section.density * section.thickness * 0.8 / 4.0 * (1.0 + rotary_ratio);
                if n == 24 {
                    let tolerance = if panel == TransverseMass::Linear {
                        0.03
                    } else {
                        2e-4
                    };
                    assert!((inertia / continuum_inertia - 1.0).abs() < tolerance);
                }
                if n == 6 && mode_x == 1 {
                    let report = crate::modes(
                        &model,
                        (0.65 * lambda, 1.35 * lambda),
                        &crate::SliceOptions::default(),
                    )
                    .unwrap();
                    assert_eq!(report.expected, 1);
                    assert!((report.modes[0].lambda.sqrt() / lambda.sqrt() - 1.0).abs() < 0.04);
                    let old_report = crate::modes(
                        &legacy,
                        (0.65 * lambda, 1.35 * lambda),
                        &crate::SliceOptions::default(),
                    )
                    .unwrap();
                    assert_eq!(old_report.expected, 1);
                    assert!(
                        (report.modes[0].lambda.sqrt() / old_report.modes[0].lambda.sqrt() - 1.0)
                            .abs()
                            > 1e-6
                    );
                }
            }
            assert!(
                errors[1] < 0.5 * errors[0],
                "{panel:?} mode {mode_x}: {errors:?}"
            );
            assert!(
                errors[2] < 0.4 * errors[1],
                "{panel:?} mode {mode_x}: {errors:?}"
            );
            let tolerance = if panel == TransverseMass::Linear {
                0.015
            } else {
                0.004
            };
            assert!(errors[2] < tolerance, "{panel:?} mode {mode_x}: {errors:?}");
            // Both discretizations approach the same continuum, but lumping the
            // slopes gives a different finite-wavelength frequency. That shift
            // must be live, decay under refinement and grow at shorter wavelength.
            assert!(
                shifts.iter().all(|s| *s > 1e-6),
                "mode {mode_x}: {shifts:?}"
            );
            assert!(shifts[1] < 0.3 * shifts[0] && shifts[2] < 0.3 * shifts[1]);
            coarse_shifts.push(shifts[0]);
        }
        assert!(coarse_shifts[1] > 8.0 * coarse_shifts[0]);
    }
}

#[test]
fn unrepresentable_rotary_inertia_refuses_without_a_mass_floor() {
    for scale in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(rotary_mass(&[0.0, 1.0, 0.0], &[0.0, 0.0, 1.0], scale).is_err());
    }
    let mut chart = chart();
    chart.element_sections.as_mut().unwrap()[0].thickness = f64::MIN_POSITIVE;
    assert!(
        chart
            .assemble_with_kinetic_mass(
                &[],
                &options(EdgeSupport::Clamped),
                TransverseMass::Linear,
                RotaryMass::Dkt,
                StiffenerMass::Lumped
            )
            .is_err()
    );
}
