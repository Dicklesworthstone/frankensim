//! G0/G1/G3: stationary pressure-jump sheet scattering, total Robin velocity,
//! finite-source reciprocity and coherent physical receiver fields.
use fs_bem::helmholtz::bc::{
    Admittance, ForceDipole, IncidentField, IncidentSource, Monopole, PlaneWave, StressQuadrupole,
    SurfaceBc,
};
use fs_bem::helmholtz::{HelmholtzError, Medium};
use fs_bem::near_field::Options;
use fs_bem::reciprocal_radiation::open_screen::{
    OpenOptions, OpenPanels, PreparedScreen, ScreenComplexValueLimits, ScreenScatteringLimits,
};
use fs_math::c64::C64;

fn surface() -> OpenPanels {
    let mut vertices = Vec::new();
    for j in 0..4 {
        for i in 0..4 {
            vertices.push([(i as f64 - 1.5) * 0.25, (j as f64 - 1.5) * 0.25, 0.0]);
        }
    }
    let mut triangles = Vec::new();
    for j in 0..3 {
        for i in 0..3 {
            let a = 4 * j + i;
            triangles.push([a, a + 1, a + 5]);
            triangles.push([a, a + 5, a + 4]);
        }
    }
    OpenPanels::new(vertices, triangles).unwrap()
}

fn options(order: usize) -> OpenOptions {
    OpenOptions {
        edge_level: 0,
        maximum_edge_level: 0,
        triangle_order: order,
        angular_refinement: 1,
        maximum_unknowns: 16,
        maximum_cache_bytes: 4_000_000,
        maximum_work: u64::MAX,
        workers: 1,
    }
}

fn limits(surface: &OpenPanels, cases: usize) -> ScreenScatteringLimits {
    ScreenScatteringLimits {
        maximum_cases: cases,
        maximum_jump_entries: cases * surface.vertices().len(),
        maximum_source_point_products: 1_000_000,
        maximum_work: u64::MAX,
        minimum_source_clearance_m: 1e-4,
    }
}

fn quadrature() -> Options {
    Options {
        relative_tolerance: 1e-10,
        maximum_depth: 14,
        maximum_kernel_evaluations: 20_000_000,
    }
}

fn close(actual: &[C64], expected: &[C64], tolerance: f64) {
    assert_eq!(actual.len(), expected.len());
    let scale = expected.iter().map(|z| z.abs()).fold(0.0, f64::max);
    let error = actual
        .iter()
        .zip(expected)
        .map(|(&a, &b)| (a - b).abs())
        .fold(0.0, f64::max);
    assert!(
        error <= tolerance * scale.max(1e-300),
        "error {error}, scale {scale}"
    );
}

fn plane() -> PlaneWave {
    PlaneWave {
        direction: [0.0, 0.0, 1.0],
        amplitude_pa: C64::new(1.0, 0.3),
    }
}

#[test]
fn rigid_and_robin_screens_balance_total_normal_velocity_and_power() {
    let surface = surface();
    let medium = Medium::air();
    let k = 1.2;
    let owner = PreparedScreen::new(&surface, 2.0, options(4), || false).unwrap();
    let incident = plane();
    let rho_c = medium.density * medium.sound_speed;
    let rho_omega = rho_c * k;
    let n = owner.unknowns();
    // Original H WITHOUT the sheet admittance: H[p]/(i*rho*omega)
    // is the weak normal velocity of the scattered double-layer field.
    let h = owner
        .complex_matrix(
            C64::from_re(k),
            ScreenComplexValueLimits {
                maximum_matrix_entries: n * n,
                maximum_workspace_bytes: 4_000_000,
                maximum_work: u64::MAX,
                maximum_growth_exponent: 8.0,
            },
            || false,
        )
        .unwrap();
    for beta in [C64::ZERO, C64::new(0.65, 0.3)] {
        let beta_nodes = vec![beta; surface.vertices().len()];
        let admittance = if beta == C64::ZERO {
            Admittance::Rigid
        } else {
            Admittance::Local {
                beta: &beta_nodes,
                declared_active: false,
            }
        };
        let bank = owner
            .solve_scattering(
                k,
                medium,
                admittance,
                IncidentField::PlaneWave(incident),
                limits(&surface, 1),
                || false,
            )
            .unwrap();
        let scattered = bank.scattered(0).unwrap();
        assert!(bank.maximum_backward_residual() < 1e-10);
        let jump = scattered.pressure_jump_pa();
        let mut active = vec![C64::ZERO; n];
        for (vertex, dof) in owner.pressure_dof_map().iter().enumerate() {
            if let Some(dof) = *dof {
                active[dof] = jump[vertex];
            }
        }
        let mut incident_normal = vec![C64::ZERO; n];
        let mut admittance_normal = vec![C64::ZERO; n];
        // Independent exact P1 triangle integration. Every triangle in this
        // declared planar fixture has area 0.25^2/2 and normal +z.
        let area = 0.25 * 0.25 / 2.0;
        for ids in surface.triangles() {
            for a in 0..3 {
                if let Some(row) = owner.pressure_dof_map()[ids[a]] {
                    incident_normal[row] =
                        incident_normal[row] + incident.amplitude_pa.scale(area / (3.0 * rho_c));
                    for b in 0..3 {
                        let mass = area / 12.0 * if a == b { 2.0 } else { 1.0 };
                        admittance_normal[row] =
                            admittance_normal[row] + (beta * jump[ids[b]]).scale(mass / rho_c);
                    }
                }
            }
        }
        let mut residual = 0.0_f64;
        let mut scale = 0.0_f64;
        for row in 0..n {
            let velocity = (0..n).fold(C64::ZERO, |sum, col| sum + h[row * n + col] * active[col])
                * C64::new(0.0, -1.0 / rho_omega);
            residual =
                residual.max((velocity + incident_normal[row] + admittance_normal[row]).abs());
            scale = scale
                .max(velocity.abs() + incident_normal[row].abs() + admittance_normal[row].abs());
        }
        assert!(
            residual <= 3e-12 * scale,
            "total Robin velocity residual {}",
            residual / scale
        );
        let drive = vec![-incident.amplitude_pa.scale(1.0 / rho_c); surface.vertices().len()];
        let reference = owner
            .solve(
                k,
                medium,
                SurfaceBc {
                    drive_velocity: &drive,
                    admittance,
                },
                || false,
            )
            .unwrap();
        close(jump, reference.pressure_jump_pa(), 3e-12);
        let emitted = scattered.radiated_power_w().unwrap();
        let absorbed = scattered.absorbed_power_w;
        assert!(emitted > 0.0 && absorbed >= 0.0);
        assert!(
            (emitted + absorbed - scattered.driving_power_w).abs()
                <= 1e-9 * scattered.driving_power_w.abs()
        );
        if beta != C64::ZERO {
            assert!(absorbed > 0.0);
        }
        close(
            &bank.fluid_force_n(0, || false).unwrap(),
            &[C64::ZERO, C64::ZERO, -scattered.force_n()],
            3e-14,
        );
    }

    // No normal illumination means no sheet jump, including for a dissipative
    // sheet. Adding a closed-wall -Y*p_inc term would fail this limit.
    let grazing = PlaneWave {
        direction: [1.0, 0.0, 0.0],
        ..incident
    };
    let beta = vec![C64::new(0.65, 0.3); surface.vertices().len()];
    let bank = owner
        .solve_scattering(
            k,
            medium,
            Admittance::Local {
                beta: &beta,
                declared_active: false,
            },
            IncidentField::PlaneWave(grazing),
            limits(&surface, 1),
            || false,
        )
        .unwrap();
    assert!(
        bank.scattered(0)
            .unwrap()
            .pressure_jump_pa()
            .iter()
            .all(|&p| p == C64::ZERO)
    );
    let points = [[0.2, 0.1, 0.5], [-0.1, 0.2, -0.4]];
    let total = bank
        .prepare_receivers(&points, 1e-5, quadrature(), 1_000_000, || false)
        .unwrap()
        .evaluate(0, || false)
        .unwrap();
    for (i, &point) in points.iter().enumerate() {
        let expected = IncidentField::PlaneWave(grazing)
            .pressure_velocity_at(k, medium, point)
            .unwrap();
        close(&[total.scalar.pressure[i]], &[expected.pressure_pa], 1e-14);
        close(
            &total.particle_velocity_m_s[i],
            &expected.velocity_m_s,
            1e-14,
        );
    }
}

#[test]
fn finite_source_diffraction_is_reciprocal_and_total_receivers_preserve_coherence() {
    let surface = surface();
    let medium = Medium::air();
    let k = 1.2;
    let owner = PreparedScreen::new(&surface, 2.0, options(8), || false).unwrap();
    let a = Monopole {
        position_m: [0.1, 0.05, 0.8],
        volume_velocity_m3_s: C64::new(1e-3, -2e-4),
    };
    let b = Monopole {
        position_m: [-0.2, 0.15, -0.9],
        ..a
    };
    let source_a = [a];
    let source_b = [b];
    let both = [IncidentSource::Monopole(a), IncidentSource::Monopole(b)];
    let mixed = [
        IncidentSource::PlaneWave(plane()),
        IncidentSource::ForceDipole(ForceDipole {
            position_m: [0.2, -0.1, -0.7],
            force_on_fluid_n: [C64::new(0.03, 0.02), C64::ZERO, C64::new(-0.01, 0.04)],
        }),
        IncidentSource::StressQuadrupole(StressQuadrupole {
            position_m: [-0.1, 0.1, 0.9],
            stress_moment_n_m: [
                C64::from_re(0.01),
                C64::ZERO,
                C64::new(0.02, 0.01),
                C64::ZERO,
                C64::new(0.003, -0.002),
                C64::ZERO,
            ],
        }),
    ];
    let fields = [
        IncidentField::Monopoles(&source_a),
        IncidentField::Monopoles(&source_b),
        IncidentField::Sources(&both),
        IncidentField::Sources(&mixed),
    ];
    let beta = vec![C64::new(0.4, 0.2); surface.vertices().len()];
    let admittance = Admittance::Local {
        beta: &beta,
        declared_active: false,
    };
    let bank = owner
        .solve_scattering_many(
            k,
            medium,
            admittance,
            &fields,
            limits(&surface, fields.len()),
            || false,
        )
        .unwrap();
    assert_eq!(bank.cases(), fields.len());
    let ab = bank
        .scattered(0)
        .unwrap()
        .pressure_at(&[b.position_m], 1e-5, quadrature(), || false)
        .unwrap();
    let ba = bank
        .scattered(1)
        .unwrap()
        .pressure_at(&[a.position_m], 1e-5, quadrature(), || false)
        .unwrap();
    close(&ab, &ba, 2e-8);
    let points = [[-0.55, 0.2, 0.5], [0.38, -0.2, -0.4]];
    let receivers = bank
        .prepare_receivers(&points, 1e-5, quadrature(), 1_000_000, || false)
        .unwrap();
    let work = receivers.kernel_evaluations();
    assert_eq!(receivers.points_m(), points);
    let mut total = Vec::new();
    for (case, field) in fields.iter().enumerate() {
        let actual = receivers.evaluate(case, || false).unwrap();
        let scattered = bank
            .scattered(case)
            .unwrap()
            .prepare_receivers(&points, 1e-5, quadrature(), 1_000_000, || false)
            .unwrap()
            .evaluate(|| false)
            .unwrap();
        for (i, &point) in points.iter().enumerate() {
            let incident = field.pressure_velocity_at(k, medium, point).unwrap();
            close(
                &[actual.scalar.pressure[i]],
                &[scattered.scalar.pressure[i] + incident.pressure_pa],
                1e-14,
            );
            close(
                &actual.particle_velocity_m_s[i],
                &core::array::from_fn::<_, 3, _>(|axis| {
                    scattered.particle_velocity_m_s[i][axis] + incident.velocity_m_s[axis]
                }),
                1e-14,
            );
        }
        total.push(actual);
    }
    for i in 0..points.len() {
        close(
            &[total[2].scalar.pressure[i]],
            &[total[0].scalar.pressure[i] + total[1].scalar.pressure[i]],
            3e-12,
        );
        close(
            &total[2].particle_velocity_m_s[i],
            &core::array::from_fn::<_, 3, _>(|axis| {
                total[0].particle_velocity_m_s[i][axis] + total[1].particle_velocity_m_s[i][axis]
            }),
            3e-12,
        );
    }
    let single = owner
        .solve_scattering(
            k,
            medium,
            admittance,
            fields[3],
            limits(&surface, 1),
            || false,
        )
        .unwrap();
    assert_eq!(
        single.scattered(0).unwrap().pressure_jump_pa(),
        bank.scattered(3).unwrap().pressure_jump_pa()
    );
    let replay = receivers.evaluate(3, || false).unwrap();
    assert_eq!(replay.scalar.pressure, total[3].scalar.pressure);
    assert_eq!(replay.particle_velocity_m_s, total[3].particle_velocity_m_s);
    assert_eq!(receivers.kernel_evaluations(), work);
}

#[test]
fn scattering_refuses_invalid_sources_limits_and_cancelled_partial_banks() {
    let surface = surface();
    let medium = Medium::air();
    let owner = PreparedScreen::new(&surface, 2.0, options(4), || false).unwrap();
    let field = IncidentField::PlaneWave(plane());
    let controls = limits(&surface, 1);
    for point in [[0.0; 3], surface.vertices()[0]] {
        let sources = [Monopole {
            position_m: point,
            volume_velocity_m3_s: C64::from_re(1e-3),
        }];
        assert!(
            owner
                .solve_scattering(
                    1.2,
                    medium,
                    Admittance::Rigid,
                    IncidentField::Monopoles(&sources),
                    controls,
                    || false
                )
                .is_err()
        );
    }
    for short in [
        ScreenScatteringLimits {
            maximum_cases: 0,
            ..controls
        },
        ScreenScatteringLimits {
            maximum_jump_entries: surface.vertices().len() - 1,
            ..controls
        },
        ScreenScatteringLimits {
            maximum_source_point_products: 1,
            ..controls
        },
        ScreenScatteringLimits {
            maximum_work: 1,
            ..controls
        },
    ] {
        assert!(
            owner
                .solve_scattering(1.2, medium, Admittance::Rigid, field, short, || false)
                .is_err()
        );
    }
    assert!(
        owner
            .solve_scattering_many(1.2, medium, Admittance::Rigid, &[], controls, || false)
            .is_err()
    );
    let negative = vec![C64::from_re(-0.1); surface.vertices().len()];
    assert!(
        owner
            .solve_scattering(
                1.2,
                medium,
                Admittance::Local {
                    beta: &negative,
                    declared_active: false
                },
                field,
                controls,
                || false
            )
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        owner.solve_scattering(1.2, medium, Admittance::Rigid, field, controls, || {
            polls += 1;
            polls >= 12
        }),
        Err(HelmholtzError::Cancelled)
    ));
    let finite = [Monopole {
        position_m: [0.0, 0.0, -0.7],
        volume_velocity_m3_s: C64::from_re(1e-3),
    }];
    let bank = owner
        .solve_scattering(
            1.2,
            medium,
            Admittance::Rigid,
            IncidentField::Monopoles(&finite),
            controls,
            || false,
        )
        .unwrap();
    assert!(
        bank.prepare_receivers(
            &[finite[0].position_m],
            1e-5,
            quadrature(),
            1_000_000,
            || false
        )
        .is_err()
    );
    let points = [[0.2, 0.1, 0.5]];
    assert!(
        bank.prepare_receivers(&points, 1e-5, quadrature(), 0, || false)
            .is_err()
    );
    let rows = bank
        .prepare_receivers(&points, 1e-5, quadrature(), 1_000_000, || false)
        .unwrap();
    let original = rows.evaluate(0, || false).unwrap();
    assert!(rows.evaluate(1, || false).is_err());
    assert!(bank.fluid_force_n(1, || false).is_err());
    assert!(matches!(
        rows.evaluate(0, || true),
        Err(HelmholtzError::Cancelled)
    ));
    let replay = rows.evaluate(0, || false).unwrap();
    assert_eq!(replay.scalar.pressure, original.scalar.pressure);
    assert_eq!(replay.particle_velocity_m_s, original.particle_velocity_m_s);
}
