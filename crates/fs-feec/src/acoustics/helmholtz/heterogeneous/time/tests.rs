use super::*;
use crate::acoustics::helmholtz::{HelmholtzFe, Medium};
use crate::fixtures::masked_cube_grid;
use crate::highorder::simplex::LocalFn;
use fs_modal::SliceOptions;
use fs_rep_mesh::TetComplex;

const MEDIUM: Medium = Medium {
    density_kg_m3: 1.3,
    sound_speed_m_s: 2.0,
};

fn brick(nx: usize, lengths: [f64; 3]) -> (TetComplex, Vec<[f64; 3]>) {
    let (complex, mut positions) = masked_cube_grid(nx, 1, 1, &|_, _, _| true);
    for point in &mut positions {
        point[0] *= lengths[0] / nx as f64;
        point[1] *= lengths[1];
        point[2] *= lengths[2];
    }
    (complex, positions)
}

fn plane(geometry: &HelmholtzFe<'_>, x: f64) -> Vec<[u32; 3]> {
    geometry
        .boundary_faces()
        .into_iter()
        .filter(|face| {
            face.iter()
                .all(|&v| (geometry.positions()[v as usize][0] - x).abs() < 1e-12)
        })
        .collect()
}

fn uniform<'f, 'c>(geometry: &'f HelmholtzFe<'c>) -> HeterogeneousHelmholtzFe<'f, 'c> {
    HeterogeneousHelmholtzFe::new(
        geometry,
        4,
        &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
        &mut || false,
    )
    .unwrap()
}

// Pointwise P2 interpolation in the ORIGINAL hierarchical basis. Its edge
// shape is lambda_a lambda_b, hence the four-times midpoint correction.
fn p2(geometry: &HelmholtzFe<'_>, field: impl Fn([f64; 3]) -> f64) -> Vec<f64> {
    assert_eq!(geometry.order(), 2);
    let mut coefficients = vec![0.0; geometry.ndof()];
    for (element, tet) in geometry.tetrahedra().iter().enumerate() {
        for (dof, shape) in geometry.space.element_dofs(element) {
            coefficients[dof] = match shape {
                LocalFn::Vertex(a) => field(geometry.positions()[tet[a] as usize]),
                LocalFn::Edge(a, b, 0) => {
                    let a = geometry.positions()[tet[a] as usize];
                    let b = geometry.positions()[tet[b] as usize];
                    let midpoint = core::array::from_fn(|axis| f64::midpoint(a[axis], b[axis]));
                    4.0 * field(midpoint) - 2.0 * (field(a) + field(b))
                }
                _ => panic!("only original P2 vertex/edge functions expected"),
            };
        }
    }
    coefficients
}

fn mass_projection(mass: &Csr, mode: &[f64], field: &[f64]) -> f64 {
    let mut product = vec![0.0; field.len()];
    mass.spmv(field, &mut product);
    mode.iter().zip(product).map(|(a, b)| a * b).sum()
}

#[test]
fn uniform_compression_has_actual_pressure_work_and_kinetic_energy() {
    let lengths = [1.0, 0.4, 0.2];
    let volume = lengths.iter().product::<f64>();
    let (complex, positions) = brick(2, lengths);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let h = 0.025;
    let prepared = model
        .prepare_time(h, &[], HeterogeneousTimeLimits::default(), &mut || false)
        .unwrap();
    // Exact polynomial solution: psi=a*x²/2+a*c²*t²/2, p=a*c²*t,
    // v_x=-a*x/rho. The inlet is stationary, the x=L wall supplies inward
    // velocity a*L/rho. The velocity field and kinetic energy stay constant.
    for a in [0.02, 2e-10] {
        let potential = p2(&geometry, |point| 0.5 * a * point[0] * point[0]);
        let mut state = prepared
            .state(&potential, &vec![0.0; geometry.ndof()])
            .unwrap();
        let initial = prepared.energy(&state).unwrap();
        let kinetic = a * a * volume * lengths[0] * lengths[0] / (6.0 * MEDIUM.density_kg_m3);
        assert!((initial.kinetic_j / kinetic - 1.0).abs() < 2e-12);
        assert_eq!(initial.compression_j, 0.0);
        let source = plane(&geometry, 1.0)
            .into_iter()
            .map(|face| (face, [a / MEDIUM.density_kg_m3; 3]))
            .collect::<Vec<_>>();
        let mut total_work = 0.0;
        for index in 1..=4 {
            let next = prepared.step(&state, &source, &mut || false).unwrap();
            let time = index as f64 * h;
            let pressure = a * MEDIUM.sound_speed_m_s.powi(2) * time;
            let compression = volume * pressure * pressure
                / (2.0 * MEDIUM.density_kg_m3 * MEDIUM.sound_speed_m_s.powi(2));
            let energy = prepared.energy(&next.state).unwrap();
            for &p in &next.state.pressure_pa()[..positions.len()] {
                assert!((p / pressure - 1.0).abs() < 2e-10);
            }
            assert!((energy.kinetic_j / kinetic - 1.0).abs() < 2e-10);
            assert!((energy.compression_j / compression - 1.0).abs() < 2e-10);
            assert_eq!(next.work.boundary_dissipation_j, 0.0);
            total_work += next.work.input_work_j;
            assert!((total_work / compression - 1.0).abs() < 2e-10);
            assert!(next.work.balance_residual_j().abs() < 2e-10 * (kinetic + compression));
            state = next.state;
        }
    }
}

#[test]
fn resolved_time_phase_matches_midpoint_and_converges_to_the_physical_mode() {
    let (complex, positions) = brick(4, [1.0, 0.2, 0.2]);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let modes = model
        .modes(
            (0.7, 1.3),
            &SliceOptions {
                mass_diagonal_equilibration: true,
                ..SliceOptions::default()
            },
            &mut || false,
        )
        .unwrap();
    assert_eq!(modes.shapes.len(), 1);
    let shape = &modes.shapes[0];
    let omega = det::sqrt(modes.omega_squared_rad2_s2[0]);
    let duration = 0.15;
    let mut errors = Vec::new();
    for count in [8, 16, 32] {
        let h = duration / count as f64;
        let prepared = model
            .prepare_time(h, &[], HeterogeneousTimeLimits::default(), &mut || false)
            .unwrap();
        let mut state = prepared.state(&vec![0.0; model.ndof()], shape).unwrap();
        let initial_energy = prepared.energy(&state).unwrap().total_j();
        for _ in 0..count {
            let step = prepared.step(&state, &[], &mut || false).unwrap();
            assert_eq!(step.work.boundary_dissipation_j, 0.0);
            assert_eq!(step.work.input_work_j, 0.0);
            state = step.state;
        }
        let angle = count as f64 * 2.0 * det::atan(0.5 * omega * h);
        let actual = mass_projection(model.mass(), shape, state.pressure_pa());
        let potential = mass_projection(model.mass(), shape, state.pressure_potential_pa_s());
        assert!((actual - det::cos(angle)).abs() < 5e-8);
        assert!((potential - det::sin(angle) / omega).abs() < 5e-8);
        assert!((prepared.energy(&state).unwrap().total_j() / initial_energy - 1.0).abs() < 2e-9);
        let error = (actual - det::cos(omega * duration)).abs();
        println!("heterogeneous pulse phase steps={count} error={error:e}");
        errors.push(error);
    }
    assert!(errors[0] > 1e-4);
    for pair in errors.windows(2) {
        assert!(
            (3.7..4.3).contains(&(pair[0] / pair[1])),
            "time phase errors {errors:?}"
        );
    }
}

fn packet_potential(x: f64) -> f64 {
    let q = (x - 0.2) / 0.16;
    if q.abs() < 1.0 {
        0.25 * 0.16 * (1.0 - q * q).powi(4)
    } else {
        0.0
    }
}

fn packet_pressure(x: f64) -> f64 {
    let q = (x - 0.2) / 0.16;
    if q.abs() < 1.0 {
        8.0 * 0.25 * q * (1.0 - q * q).powi(3)
    } else {
        0.0
    }
}

#[test]
fn a_broadband_packet_reflects_at_a_density_interface_with_the_right_sign_and_delay() {
    // c is IDENTICAL on both sides. Only rho (and therefore impedance) jumps.
    // The analytic pressure amplitudes are R=0.6, T=1.6 at x=1/2.
    let reflection = 0.6;
    let transmission = 1.6;
    let duration = 0.5;
    let mut errors = Vec::new();
    for nx in [8, 16, 24] {
        let (complex, positions) = brick(nx, [1.0, 0.04, 0.04]);
        let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
        let model = HeterogeneousHelmholtzFe::new(
            &geometry,
            4,
            &mut |_, point| {
                Ok::<_, AcousticFeError>(Medium {
                    density_kg_m3: if point[0] < 0.5 { 1.0 } else { 4.0 },
                    sound_speed_m_s: 1.0,
                })
            },
            &mut || false,
        )
        .unwrap();
        let steps = 8 * nx;
        let prepared = model
            .prepare_time(
                duration / steps as f64,
                &[
                    ResistiveBoundary {
                        faces: plane(&geometry, 0.0),
                        admittance_m_pa_s: 1.0,
                    },
                    ResistiveBoundary {
                        faces: plane(&geometry, 1.0),
                        admittance_m_pa_s: 0.25,
                    },
                ],
                HeterogeneousTimeLimits::default(),
                &mut || false,
            )
            .unwrap();
        let potential = p2(&geometry, |point| packet_potential(point[0]));
        let pressure = p2(&geometry, |point| packet_pressure(point[0]));
        let mut state = prepared.state(&potential, &pressure).unwrap();
        let initial_energy = prepared.energy(&state).unwrap().total_j();
        let mut heat = 0.0;
        for _ in 0..steps {
            let step = prepared.step(&state, &[], &mut || false).unwrap();
            assert_eq!(step.work.input_work_j, 0.0);
            assert!(step.work.boundary_dissipation_j >= 0.0);
            heat += step.work.boundary_dissipation_j;
            state = step.state;
        }
        let final_energy = prepared.energy(&state).unwrap().total_j();
        assert!((final_energy + heat - initial_energy).abs() < 2e-8 * initial_energy);
        let mut error: f64 = 0.0;
        for (vertex, point) in positions.iter().enumerate() {
            let exact = if point[0] < 0.5 {
                packet_pressure(point[0] - duration)
                    + reflection * packet_pressure(1.0 - point[0] - duration)
            } else {
                transmission * packet_pressure(point[0] - duration)
            };
            error = error.max((state.pressure_pa()[vertex] - exact).abs());
        }
        println!("heterogeneous density pulse nx={nx} error={error:e}");
        errors.push(error);
        if nx == 24 {
            let at = |x: f64| {
                let vertex = positions
                    .iter()
                    .position(|point| {
                        (point[0] - x).abs() < 1e-12 && point[1] == 0.0 && point[2] == 0.0
                    })
                    .unwrap();
                state.pressure_pa()[vertex]
            };
            let incident = packet_pressure(0.25);
            assert!((at(0.25) / incident - reflection).abs() < 0.05);
            assert!((at(0.75) / incident - transmission).abs() < 0.05);
        }
    }
    assert!(
        errors[2] < 0.04 && errors[2] < 0.4 * errors[0],
        "pulse errors {errors:?}"
    );
}

#[test]
fn rigid_uniform_pressure_and_zero_state_do_not_create_heat() {
    let (complex, positions) = brick(2, [1.0, 0.2, 0.2]);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let prepared = model
        .prepare_time(0.02, &[], HeterogeneousTimeLimits::default(), &mut || false)
        .unwrap();
    let zero = prepared.zero_state();
    let zero_next = prepared.step(&zero, &[], &mut || false).unwrap();
    assert!(zero_next.state.pressure_pa().iter().all(|&p| p == 0.0));
    assert_eq!(zero_next.work.balance_residual_j(), 0.0);
    let mut state = prepared
        .state(&p2(&geometry, |_| 2.0), &p2(&geometry, |_| 0.7))
        .unwrap();
    let expected =
        0.5 * 0.7_f64.powi(2) * 0.04 / (MEDIUM.density_kg_m3 * MEDIUM.sound_speed_m_s.powi(2));
    for index in 1..=10 {
        let step = prepared.step(&state, &[], &mut || false).unwrap();
        assert_eq!(step.work.boundary_dissipation_j, 0.0);
        let energy = prepared.energy(&step.state).unwrap();
        assert!((energy.compression_j / expected - 1.0).abs() < 1e-10);
        assert!(energy.kinetic_j.abs() < 1e-10 * expected);
        for &potential in &step.state.pressure_potential_pa_s()[..positions.len()] {
            assert!((potential - (2.0 + index as f64 * 0.02 * 0.7)).abs() < 1e-10);
        }
        state = step.state;
    }
}

#[test]
fn sparse_budgets_cancellation_and_original_state_binding_refuse_without_mutation() {
    let (complex, positions) = brick(1, [1.0; 3]);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let other = uniform(&geometry);
    assert!(matches!(
        model.prepare_time(
            0.01,
            &[],
            HeterogeneousTimeLimits {
                max_factor_bytes: 1,
                ..HeterogeneousTimeLimits::default()
            },
            &mut || false
        ),
        Err(HeterogeneousError::Acoustic(AcousticFeError::Factor(
            fs_sparse::direct::LdltError::FactorBudget { .. }
        )))
    ));
    let prepared = model
        .prepare_time(0.01, &[], HeterogeneousTimeLimits::default(), &mut || false)
        .unwrap();
    let foreign = other
        .prepare_time(0.01, &[], HeterogeneousTimeLimits::default(), &mut || false)
        .unwrap();
    let state = prepared
        .state(&p2(&geometry, |point| point[0]), &p2(&geometry, |_| 0.0))
        .unwrap();
    let before = state.pressure_potential_pa_s().to_vec();
    assert!(foreign.step(&state, &[], &mut || false).is_err());
    assert!(matches!(
        prepared.step(&state, &[], &mut || true),
        Err(HeterogeneousError::Cancelled)
    ));
    assert_eq!(state.pressure_potential_pa_s(), before);
    assert!(state.pressure_pa().iter().all(|&p| p == 0.0));
    let active = [ResistiveBoundary {
        faces: plane(&geometry, 1.0),
        admittance_m_pa_s: -1.0,
    }];
    assert!(
        model
            .prepare_time(
                0.01,
                &active,
                HeterogeneousTimeLimits::default(),
                &mut || false
            )
            .is_err()
    );
}
