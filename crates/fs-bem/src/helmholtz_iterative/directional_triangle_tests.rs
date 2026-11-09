use super::super::solve_radiation_with_rule;
use super::*;

const K: f64 = 34.0;
const RADIUS: f64 = 0.01;
fn medium() -> Medium {
    Medium {
        density: 1.2,
        sound_speed: 343.0,
    }
}
fn surface(level: u32) -> SpherePanels {
    let unit = SpherePanels::icosphere(RADIUS, level).unwrap();
    let mut triangles = Vec::new();
    for center in [-1.0, 1.0] {
        triangles.extend(
            unit.triangles()
                .unwrap()
                .iter()
                .map(|t| t.map(|p| [p[0], p[1], p[2] + center])),
        );
    }
    SpherePanels::from_triangles(triangles).unwrap()
}
fn controls(order: usize) -> IterativeOptions {
    let directional = DirectionalOptions {
        order,
        leaf_capacity: 64,
        ..DirectionalOptions::for_tolerance(1e-2).unwrap()
    };
    IterativeOptions {
        tol: 1e-10,
        restart: 40,
        max_iters: 240,
        directional: Some(directional),
        precond: Preconditioner::NearFieldBlockJacobi {
            max_block_panels: 20,
        },
        panel_rule: PanelRule::Triangle {
            max_pairs: 2_000_000,
        },
        ..IterativeOptions::default()
    }
}
fn relative(a: &[C64], b: &[C64]) -> f64 {
    (a.iter()
        .zip(b)
        .map(|(a, b)| (*a - *b).norm_sq())
        .sum::<f64>()
        / b.iter().map(|v| v.norm_sq()).sum::<f64>())
    .sqrt()
}
// Independent manufactured outgoing sources lie inside the closed surfaces.
// Boundary p and grad(p) use the closed-form spherical wave directly, without
// production Green kernels, source samplers, panel quadrature or FMM.
fn exact(point: [f64; 3], k: f64) -> (C64, [C64; 3]) {
    let (mut pressure, mut velocity) = (C64::ZERO, [C64::ZERO; 3]);
    for (center, amplitude) in [
        (-1.0, C64::from_re(RADIUS)),
        (1.0, C64::new(0.4, 0.7).scale(RADIUS)),
    ] {
        let displacement = [point[0], point[1], point[2] - center];
        let r = displacement.iter().map(|v| v * v).sum::<f64>().sqrt();
        let phase = C64::new((k * r).cos(), (k * r).sin());
        let p = amplitude * phase.scale(1.0 / r);
        let radial =
            p * C64::new(-1.0 / r, k) / C64::new(0.0, k * medium().density * medium().sound_speed);
        pressure = pressure + p;
        for axis in 0..3 {
            velocity[axis] = velocity[axis] + radial.scale(displacement[axis] / r);
        }
    }
    (pressure, velocity)
}
fn traces(surface: &SpherePanels, k: f64) -> (Vec<C64>, Vec<C64>) {
    surface
        .centroids()
        .iter()
        .zip(surface.normals())
        .map(|(&point, &normal)| {
            let (p, v) = exact(point, k);
            (
                p,
                v.iter()
                    .zip(normal)
                    .fold(C64::ZERO, |sum, (&v, n)| sum + v.scale(n)),
            )
        })
        .unzip()
}
fn solve(
    owner: &IterativeRadiation<'_>,
    drive: &[C64],
    gate: &CancelGate,
) -> (RadiationSolution, SolveReport, ResourceReceipt) {
    let mut state = owner.start(drive, gate).unwrap();
    owner
        .advance(&mut state, owner.options.max_iters + 1, gate)
        .unwrap();
    owner.finish(&state).unwrap()
}

#[test]
fn directional_triangle_matches_dense_above_low_frequency_limit() {
    let (surface, gate) = (surface(0), CancelGate::new_clock_free());
    let (_, drive) = traces(&surface, K);
    let dense = solve_radiation_with_rule(
        &surface,
        K,
        medium(),
        &drive,
        Formulation::BurtonMiller,
        PanelRule::Triangle {
            max_pairs: 2_000_000,
        },
    )
    .unwrap();
    let mut errors = Vec::new();
    for order in [4, 6, 8] {
        let owner = IterativeRadiation::prepare(
            &surface,
            K,
            medium(),
            Formulation::BurtonMiller,
            &controls(order),
            &gate,
        )
        .unwrap();
        if let Layers::BurtonMillerCloud(cloud) = &owner.layers {
            assert!(
                cloud
                    .double
                    .directional_stats()
                    .unwrap()
                    .high_frequency_pairs
                    > 0
            );
            assert_eq!(cloud.double.directional_stats().unwrap().near_pairs, 0);
        } else {
            panic!("original triangle-cloud owner required");
        }
        let (actual, report, receipt) = solve(&owner, &drive, &gate);
        assert!(report.converged_euclidean());
        assert_eq!(receipt.fmm_order, order);
        errors.push(relative(&actual.pressure, &dense.pressure));
    }
    eprintln!("directional triangle dense agreement, order4/6/8: {errors:?}");
    assert!(errors[1] < errors[0] && errors[2] < errors[1]);
    assert!(errors[2] < 2e-6);
    assert!(matches!(
        IterativeRadiation::prepare(
            &surface,
            K,
            medium(),
            Formulation::BurtonMiller,
            &IterativeOptions {
                directional: None,
                ..controls(8)
            },
            &gate
        ),
        Err(IterativeError::Fmm(FmmError::HelmholtzRegime { .. }))
    ));
}

#[test]
fn directional_triangle_robin_fields_and_radiated_work_converge_to_outgoing_sources() {
    let gate = CancelGate::new_clock_free();
    let (beta, z) = (C64::new(0.2, 0.7), medium().density * medium().sound_speed);
    let a = C64::from_re(RADIUS);
    let b = C64::new(0.4, 0.7).scale(RADIUS);
    let exact_power = 2.0 * core::f64::consts::PI / z
        * (a.norm_sq() + b.norm_sq() + 2.0 * (a * b.conj()).re * (2.0 * K).sin() / (2.0 * K));
    let (mut previous_field, mut previous_power) = (f64::INFINITY, f64::INFINITY);
    for level in 0..=2 {
        let surface = surface(level);
        let (pressure, velocity) = traces(&surface, K);
        let beta_values = vec![beta; pressure.len()];
        let drive = pressure
            .iter()
            .zip(&velocity)
            .map(|(&p, &v)| v + beta * p.scale(1.0 / z))
            .collect::<Vec<_>>();
        let owner = IterativeRadiation::prepare_bc(
            &surface,
            K,
            medium(),
            Formulation::BurtonMiller,
            &controls(8),
            super::super::bc::Admittance::Local {
                beta: &beta_values,
                declared_active: false,
            },
            &gate,
        )
        .unwrap();
        let (actual, report, _) = solve(&owner, &drive, &gate);
        assert!(report.converged_euclidean());
        let field =
            relative(&actual.pressure, &pressure).max(relative(&actual.velocity, &velocity));
        let work = (actual.power_diagnostics.surface_power.re - exact_power).abs() / exact_power;
        eprintln!(
            "high-kD triangle Robin level={level}: panels={}, p/u error={field:e}, radiated-work error={work:e}",
            pressure.len()
        );
        assert!(field < previous_field && work < previous_power);
        (previous_field, previous_power) = (field, work);
    }
    assert!(previous_field < 0.02 && previous_power < 0.04);
}

#[test]
fn directional_triangle_static_partition_low_frequency_bits_and_checkpoint_binding() {
    let (surface, gate) = (surface(0), CancelGate::new_clock_free());
    let options = controls(6);
    let owner = IterativeRadiation::prepare(
        &surface,
        K,
        medium(),
        Formulation::BurtonMiller,
        &options,
        &gate,
    )
    .unwrap();
    let triangles = surface.triangles().unwrap();
    let (points, normals, _, _) = radon_cloud(&surface, triangles);
    let static_op = FmmC64::new_double_far_field_within(
        HelmholtzDouble { k: 0.0 },
        points.clone(),
        normals,
        6,
        64 * RADON7.len(),
        options.budget(),
    )
    .unwrap();
    let static_partition = CloudPartition::new(&static_op, surface.areas().len(), points.len());
    if let Layers::BurtonMillerCloud(cloud) = &owner.layers {
        let dynamic = CloudPartition::new(&cloud.double, surface.areas().len(), points.len());
        assert_ne!(
            dynamic.panels(0),
            static_partition.panels(0),
            "fixture must exercise distinct static/dynamic corrections"
        );
    }
    let k = 1.0;
    let (_, drive) = traces(&surface, k);
    let low = IterativeRadiation::prepare(
        &surface,
        k,
        medium(),
        Formulation::BurtonMiller,
        &options,
        &gate,
    )
    .unwrap();
    let original = IterativeRadiation::prepare(
        &surface,
        k,
        medium(),
        Formulation::BurtonMiller,
        &IterativeOptions {
            directional: None,
            fmm_order: 6,
            leaf_capacity: 64,
            ..options.clone()
        },
        &gate,
    )
    .unwrap();
    let (a, _, _) = solve(&low, &drive, &gate);
    let (b, _, _) = solve(&original, &drive, &gate);
    assert_eq!(a.pressure, b.pressure);
    assert_eq!(a.velocity, b.velocity);
    let mut checkpoint = low.start(&drive, &gate).unwrap();
    assert!(matches!(
        original.advance(&mut checkpoint, 1, &gate),
        Err(IterativeError::Invalid { .. })
    ));
    let tiny = IterativeOptions {
        max_near_entries: 1,
        ..options
    };
    assert!(matches!(
        IterativeRadiation::prepare(
            &surface,
            K,
            medium(),
            Formulation::BurtonMiller,
            &tiny,
            &gate
        ),
        Err(IterativeError::NearFieldBudget { .. })
    ));
    gate.request();
    assert!(matches!(
        owner.start(&drive, &gate),
        Err(IterativeError::Cancelled { .. })
    ));
}
