//! G0/G1/G3: passive porous-screen mechanical work, coherent heat and radiation.
//! These compare the original P1 source/observer; they do not certify a porous
//! material law or convergence of this small mesh to a continuum screen.
use fs_bem::helmholtz::bc::{Admittance, SurfaceBc};
use fs_bem::helmholtz::{HelmholtzError, Medium};
use fs_bem::reciprocal_radiation::open_screen::{
    OpenOptions, OpenPanels, PreparedScreen, ScreenPortLimits,
};
use fs_math::c64::C64;

fn screen() -> OpenPanels {
    let mut vertices = Vec::new();
    for j in 0..4 {
        for i in 0..4 {
            let x = (i as f64 - 1.5) * 0.025;
            let y = (j as f64 - 1.5) * 0.025;
            vertices.push([x, y, 2.0 * x * y]);
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

fn options() -> OpenOptions {
    OpenOptions {
        edge_level: 0,
        maximum_edge_level: 0,
        triangle_order: 4,
        angular_refinement: 1,
        maximum_unknowns: 16,
        maximum_cache_bytes: 1_000_000,
        maximum_work: u64::MAX,
        workers: 1,
    }
}

fn limits(surface: &OpenPanels) -> ScreenPortLimits {
    ScreenPortLimits {
        maximum_ports: 2,
        maximum_matrix_entries: 4,
        maximum_jump_entries: 3 * surface.vertices().len(),
        maximum_work: u64::MAX,
    }
}

fn traces(surface: &OpenPanels) -> (Vec<Vec<f64>>, Vec<Vec<[f64; 3]>>) {
    let nodal = vec![
        surface
            .vertices()
            .iter()
            .map(|p| 1.0 + 2.0 * p[0])
            .collect::<Vec<_>>(),
        surface
            .vertices()
            .iter()
            .map(|p| 0.55 + 6.0 * p[0] - 4.0 * p[1])
            .collect(),
    ];
    let rows = nodal
        .iter()
        .map(|v| {
            surface
                .triangles()
                .iter()
                .map(|t| t.map(|i| v[i]))
                .collect()
        })
        .collect();
    (nodal, rows)
}

fn area(surface: &OpenPanels, triangle: [usize; 3]) -> f64 {
    let p = triangle.map(|i| surface.vertices()[i]);
    let u: [f64; 3] = core::array::from_fn(|i| p[1][i] - p[0][i]);
    let v: [f64; 3] = core::array::from_fn(|i| p[2][i] - p[0][i]);
    let cross = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    0.5 * cross[0].hypot(cross[1]).hypot(cross[2])
}

fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-300),
        "{actual:?} != {expected:?}"
    );
}

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-300),
        "{actual:.17e} != {expected:.17e}"
    );
}

// Exact cubic barycentric moments, independent of the production Duffy rule.
// Integral lambda_i lambda_j lambda_k is A/10, A/30 or A/60.
fn wall_heat(surface: &OpenPanels, beta: &[C64], jump: &[C64], medium: Medium) -> f64 {
    let mut total = 0.0;
    for &t in surface.triangles() {
        let a = area(surface, t);
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    let divisor = if i == j && j == k {
                        10.0
                    } else if i == j || j == k || i == k {
                        30.0
                    } else {
                        60.0
                    };
                    total += a / divisor * beta[t[i]].re * (jump[t[j]] * jump[t[k]].conj()).re;
                }
            }
        }
    }
    total / (2.0 * medium.density * medium.sound_speed)
}

#[test]
fn porous_screen_ports_keep_the_original_work_pairing_and_separate_heat_from_sound() {
    let surface = screen();
    let medium = Medium {
        density: 1.3,
        sound_speed: 340.0,
    };
    let owner = PreparedScreen::new(&surface, 8.0, options(), || false).unwrap();
    let (nodal, rows) = traces(&surface);
    let refs: Vec<_> = rows.iter().map(Vec::as_slice).collect();
    let mut beta: Vec<_> = surface
        .vertices()
        .iter()
        .map(|p| C64::new(1.5 + 4.0 * p[0], 0.8 - 3.0 * p[1]))
        .collect();
    let original_beta = beta.clone();
    let batch = owner
        .solve_velocity_ports_with_local_admittance(
            4.0,
            medium,
            &refs,
            &beta,
            limits(&surface),
            || false,
        )
        .unwrap();
    assert_eq!(batch.local_admittance(), Some(original_beta.as_slice()));
    beta.fill(C64::new(99.0, -42.0));
    assert_eq!(batch.local_admittance(), Some(original_beta.as_slice()));
    let law = Admittance::Local {
        beta: &original_beta,
        declared_active: false,
    };
    for j in 0..2 {
        let drive: Vec<_> = nodal[j].iter().map(|&v| C64::from_re(v)).collect();
        let direct = owner
            .solve(
                4.0,
                medium,
                SurfaceBc {
                    drive_velocity: &drive,
                    admittance: law,
                },
                || false,
            )
            .unwrap();
        for (&actual, &expected) in batch.solutions()[j]
            .pressure_jump_pa()
            .iter()
            .zip(direct.pressure_jump_pa())
        {
            close(actual, expected, 2e-12);
        }
        for i in 0..2 {
            let mut force = C64::ZERO;
            for &t in surface.triangles() {
                let a = area(&surface, t);
                for r in 0..3 {
                    for s in 0..3 {
                        force = force
                            + direct.pressure_jump_pa()[t[s]].scale(
                                nodal[i][t[r]] * a * if r == s { 1.0 / 6.0 } else { 1.0 / 12.0 },
                            );
                    }
                }
            }
            close(batch.impedance()[2 * i + j], force, 2e-12);
        }
    }
    close(batch.impedance()[1], batch.impedance()[2], 2e-12);
    let rigid = owner
        .solve_velocity_ports(4.0, medium, &refs, limits(&surface), || false)
        .unwrap();
    assert!(
        (batch.impedance()[0] - rigid.impedance()[0]).abs() > 1e-3 * rigid.impedance()[0].abs()
    );

    let rates = [C64::new(1.0, 0.4), C64::new(-0.7, 0.2)];
    let drive: Vec<_> = (0..surface.vertices().len())
        .map(|i| rates[0].scale(nodal[0][i]) + rates[1].scale(nodal[1][i]))
        .collect();
    let direct = owner
        .solve(
            4.0,
            medium,
            SurfaceBc {
                drive_velocity: &drive,
                admittance: law,
            },
            || false,
        )
        .unwrap();
    let jump = batch.pressure_jump(&rates, || false).unwrap();
    for (&actual, &expected) in jump.iter().zip(direct.pressure_jump_pa()) {
        close(actual, expected, 2e-12);
    }
    let input = batch.input_power_w(&rates).unwrap();
    let heat = batch.absorbed_power_w(&rates, || false).unwrap();
    let sound = batch.radiated_power_w(&rates, || false).unwrap();
    relative(input, direct.driving_power_w, 2e-12);
    relative(heat, direct.absorbed_power_w, 2e-12);
    relative(
        heat,
        wall_heat(&surface, &original_beta, &jump, medium),
        2e-12,
    );
    relative(sound, direct.radiated_power_w().unwrap(), 2e-12);
    relative(input, heat + sound, 1e-8);
    assert!(heat > 0.0 && sound > 0.0);
    let incoherent: f64 = rates
        .iter()
        .zip(batch.solutions())
        .map(|(rate, source)| rate.norm_sq() * source.absorbed_power_w)
        .sum();
    assert!(
        (heat - incoherent).abs() > 0.1 * heat,
        "cross-port heat must retain destructive interference"
    );
}

#[test]
fn reactive_and_high_q_sheet_heat_keep_clock_phase_and_zero_admittance_limit() {
    let surface = screen();
    let medium = Medium {
        density: 1.3,
        sound_speed: 340.0,
    };
    let owner = PreparedScreen::new(&surface, 8.0, options(), || false).unwrap();
    let (_, rows) = traces(&surface);
    let refs: Vec<_> = rows.iter().map(Vec::as_slice).collect();
    let rates = [C64::new(0.7, -0.4), C64::new(0.3, 0.6)];
    for conductance in [0.0, 1e-20, 0.3] {
        let beta: Vec<_> = surface
            .vertices()
            .iter()
            .map(|p| C64::new(conductance * (1.0 + 3.0 * p[0]), 0.75 + 4.0 * p[1]))
            .collect();
        let batch = owner
            .solve_velocity_ports_with_local_admittance(
                4.0,
                medium,
                &refs,
                &beta,
                limits(&surface),
                || false,
            )
            .unwrap();
        let jump = batch.pressure_jump(&rates, || false).unwrap();
        let heat = batch.absorbed_power_w(&rates, || false).unwrap();
        relative(heat, wall_heat(&surface, &beta, &jump, medium), 2e-12);
        if conductance == 0.0 {
            assert_eq!(heat, 0.0);
        } else {
            assert!(heat > 0.0);
        }
        for phase in [0.41_f64, 1.71, 2.83] {
            let clock = C64::new(phase.cos(), phase.sin());
            relative(
                batch
                    .absorbed_power_w(&rates.map(|v| v * clock), || false)
                    .unwrap(),
                heat,
                2e-13,
            );
        }
    }
    let zeros = vec![C64::ZERO; surface.vertices().len()];
    let porous = owner
        .solve_velocity_ports_with_local_admittance(
            4.0,
            medium,
            &refs,
            &zeros,
            limits(&surface),
            || false,
        )
        .unwrap();
    let rigid = owner
        .solve_velocity_ports(4.0, medium, &refs, limits(&surface), || false)
        .unwrap();
    assert_eq!(porous.impedance(), rigid.impedance());
    for (a, b) in porous.solutions().iter().zip(rigid.solutions()) {
        assert_eq!(a.pressure_jump_pa(), b.pressure_jump_pa());
    }
    assert_eq!(porous.absorbed_power_w(&rates, || false).unwrap(), 0.0);
    assert_eq!(rigid.absorbed_power_w(&rates, || false).unwrap(), 0.0);
    assert!(rigid.local_admittance().is_none());
}

#[test]
fn passive_sheet_law_and_complete_batch_resources_are_admitted_before_solving() {
    let surface = screen();
    let medium = Medium {
        density: 1.3,
        sound_speed: 340.0,
    };
    let owner = PreparedScreen::new(&surface, 8.0, options(), || false).unwrap();
    let (_, rows) = traces(&surface);
    let refs: Vec<_> = rows.iter().map(Vec::as_slice).collect();
    let mut beta = vec![C64::new(0.3, -0.7); surface.vertices().len()];
    let solve = |beta: &[C64], controls| {
        owner.solve_velocity_ports_with_local_admittance(4.0, medium, &refs, beta, controls, || {
            false
        })
    };
    assert!(solve(&beta[..beta.len() - 1], limits(&surface)).is_err());
    for bad in [
        C64::new(-1e-30, 0.0),
        C64::new(f64::NAN, 0.0),
        C64::new(0.0, f64::INFINITY),
    ] {
        beta[5] = bad;
        assert!(solve(&beta, limits(&surface)).is_err());
    }
    beta[5] = C64::new(0.3, -0.7);
    let mut controls = limits(&surface);
    controls.maximum_jump_entries -= 1;
    assert!(solve(&beta, controls).is_err());
    controls = limits(&surface);
    controls.maximum_work = owner.velocity_port_work_with_local_admittance(2).unwrap() - 1;
    assert!(solve(&beta, controls).is_err());
    assert!(matches!(
        owner.solve_velocity_ports_with_local_admittance(
            4.0,
            medium,
            &refs,
            &beta,
            limits(&surface),
            || true,
        ),
        Err(HelmholtzError::Cancelled)
    ));
    let batch = solve(&beta, limits(&surface)).unwrap();
    assert!(batch.absorbed_power_w(&[C64::ONE], || false).is_err());
    assert!(matches!(
        batch.absorbed_power_w(&[C64::ONE; 2], || true),
        Err(HelmholtzError::Cancelled)
    ));
}
