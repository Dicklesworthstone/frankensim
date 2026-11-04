use super::*;
use crate::acoustics::helmholtz::HelmholtzFe;
use crate::acoustics::helmholtz::heterogeneous::HarmonicLimits;
use crate::acoustics::helmholtz::heterogeneous::time::tests::{MEDIUM, brick, p2, plane, uniform};
use fs_math::det;

fn restricted(prepared: &PreparedMechanicalTime<'_, '_, '_>, field: &[f64]) -> Vec<f64> {
    prepared
        .wall_dofs(0)
        .unwrap()
        .iter()
        .map(|&i| field[i])
        .collect()
}

#[test]
fn exact_wall_motion_preserves_physical_area_storage_and_source_work() {
    let pin = WallPin {
        surface_density: 0.4,
        stiffness_per_area: 1.7,
        resistance: 0.3,
    };
    for height in [0.1, 0.3] {
        let lengths = [1.0, height, 0.2];
        let area = lengths[1] * lengths[2];
        let (complex, positions) = brick(2, lengths);
        let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
        let model = uniform(&geometry);
        let h = 0.0125;
        let prepared = model
            .prepare_time_with_walls(
                h,
                &[],
                &[MechanicalBoundary {
                    faces: plane(&geometry, 1.0),
                    wall: pin,
                }],
                MechanicalTimeLimits::default(),
                &mut || false,
            )
            .unwrap();
        for amplitude in [1.0, 1e-8] {
            // Exact polynomial solution, not a discrete manufactured RHS:
            // w=w0+v*t, p=r*v+K*w; psi=A*x²/2+B*x+integral(p dt),
            // A=K*v/c², B=-rho*v-A*L. Only x=0 receives a constant source.
            let v = 0.02 * amplitude;
            let w0 = 0.03 * amplitude;
            let a = pin.stiffness_per_area * v / MEDIUM.sound_speed_m_s.powi(2);
            let b = -MEDIUM.density_kg_m3 * v - a;
            let p0 = pin.resistance * v + pin.stiffness_per_area * w0;
            let potential = p2(&geometry, |point| {
                0.5 * a * point[0] * point[0] + b * point[0]
            });
            let pressure = p2(&geometry, |_| p0);
            let displacement = restricted(&prepared, &p2(&geometry, |_| w0));
            let velocity = restricted(&prepared, &p2(&geometry, |_| v));
            let mut state = prepared
                .state_with_walls(&potential, &pressure, &[displacement], &[velocity])
                .unwrap();
            let initial = prepared.energy(&state).unwrap().total_j();
            let source = plane(&geometry, 0.0)
                .into_iter()
                .map(|face| (face, [-b / MEDIUM.density_kg_m3; 3]))
                .collect::<Vec<_>>();
            let mut work = 0.0;
            let mut heat = 0.0;
            for index in 1..=8 {
                let step = prepared.step(&state, &source, &mut || false).unwrap();
                let time = index as f64 * h;
                let w = w0 + v * time;
                let p = pin.resistance * v + pin.stiffness_per_area * w;
                let energy = prepared.energy(&step.state).unwrap();
                let wall_kinetic = 0.5 * pin.surface_density * area * v * v;
                let wall_elastic = 0.5 * pin.stiffness_per_area * area * w * w;
                let fluid_kinetic =
                    area / (2.0 * MEDIUM.density_kg_m3) * (a * a / 3.0 + a * b + b * b);
                let fluid_compression =
                    area * p * p / (2.0 * MEDIUM.density_kg_m3 * MEDIUM.sound_speed_m_s.powi(2));
                for (actual, expected) in [
                    (energy.wall_kinetic_j, wall_kinetic),
                    (energy.wall_elastic_j, wall_elastic),
                    (energy.acoustic.kinetic_j, fluid_kinetic),
                    (energy.acoustic.compression_j, fluid_compression),
                ] {
                    assert!(
                        (actual / expected - 1.0).abs() < 3e-10,
                        "energy actual={actual:e} expected={expected:e}"
                    );
                }
                for &actual in &step.state.pressure_pa()[..positions.len()] {
                    assert!((actual / p - 1.0).abs() < 3e-10);
                }
                assert_eq!(step.work.boundary_dissipation_j, 0.0);
                assert!(
                    (step.work.wall_dissipation_j / (h * pin.resistance * area * v * v) - 1.0)
                        .abs()
                        < 3e-10
                );
                let midpoint_pressure =
                    pin.resistance * v + pin.stiffness_per_area * (w0 + v * (time - 0.5 * h));
                let wall_work = h * area * midpoint_pressure * v;
                assert!((step.work.wall_work_j / wall_work - 1.0).abs() < 3e-10);
                work += step.work.input_work_j;
                heat += step.work.wall_dissipation_j;
                assert!((energy.total_j() + heat - initial - work).abs() < 3e-10 * initial);
                state = step.state;
            }
        }
    }
}

#[test]
fn nonuniform_full_trace_returns_stored_energy_and_only_resistance_makes_heat() {
    let lengths = [1.0, 0.2, 0.3];
    let area = lengths[1] * lengths[2];
    let (complex, positions) = brick(2, lengths);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = uniform(&geometry);
    for resistance in [0.0, 0.7] {
        let pin = WallPin {
            surface_density: 0.6,
            stiffness_per_area: 3.0,
            resistance,
        };
        let prepared = model
            .prepare_time_with_walls(
                0.005,
                &[],
                &[MechanicalBoundary {
                    faces: plane(&geometry, 1.0),
                    wall: pin,
                }],
                MechanicalTimeLimits::default(),
                &mut || false,
            )
            .unwrap();
        assert!(
            prepared
                .wall_dofs(0)
                .unwrap()
                .iter()
                .any(|&i| i >= positions.len())
        );
        let amplitude = 0.01;
        let displacement = restricted(
            &prepared,
            &p2(&geometry, |point| {
                amplitude * (1.0 + point[1] / lengths[1]).powi(2)
            }),
        );
        let velocity = vec![0.0; displacement.len()];
        let zeros = vec![0.0; geometry.ndof()];
        let mut state = prepared
            .state_with_walls(&zeros, &zeros, &[displacement], &[velocity])
            .unwrap();
        let initial = prepared.energy(&state).unwrap();
        // Integral_0^1 (1+y)^4 dy = 31/5. Uniform patch averaging or dropping
        // edge bubbles cannot reproduce this independent fourth-order moment.
        let expected = 0.5 * pin.stiffness_per_area * area * amplitude * amplitude * 31.0 / 5.0;
        assert!((initial.wall_elastic_j / expected - 1.0).abs() < 3e-12);
        assert_eq!(initial.acoustic.total_j(), 0.0);
        assert_eq!(initial.wall_kinetic_j, 0.0);
        let mut heat = 0.0;
        for _ in 0..4 {
            let step = prepared.step(&state, &[], &mut || false).unwrap();
            assert_eq!(step.work.input_work_j, 0.0);
            assert_eq!(step.work.boundary_dissipation_j, 0.0);
            assert!(
                step.work.wall_work_j < 0.0,
                "spring must return energy into the fluid"
            );
            if resistance == 0.0 {
                assert_eq!(step.work.wall_dissipation_j, 0.0);
            } else {
                assert!(step.work.wall_dissipation_j > 0.0);
            }
            assert!(step.work.wall_balance_residual_j().abs() < 1e-10 * expected);
            assert!(step.work.acoustic_balance_residual_j().abs() < 1e-10 * expected);
            heat += step.work.wall_dissipation_j;
            state = step.state;
        }
        let energy = prepared.energy(&state).unwrap();
        assert!(energy.acoustic.compression_j > 0.0);
        assert!((energy.total_j() + heat - expected).abs() < 1e-10 * expected);
        let trace = prepared
            .wall_dofs(0)
            .unwrap()
            .iter()
            .copied()
            .filter(|&i| i < positions.len())
            .map(|i| state.pressure_pa()[i])
            .collect::<Vec<_>>();
        let spread = trace.iter().copied().fold(f64::NEG_INFINITY, f64::max)
            - trace.iter().copied().fold(f64::INFINITY, f64::min);
        assert!(
            spread > 1e-7,
            "the local wall must retain its nonuniform pressure field"
        );
    }
}

#[test]
fn reactive_wall_reflection_phase_agrees_with_original_impedance_and_time_refines() {
    let (complex, positions) = brick(8, [1.0, 0.1, 0.1]);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let inlet = plane(&geometry, 0.0);
    let wall_faces = plane(&geometry, 1.0);
    let z0 = MEDIUM.density_kg_m3 * MEDIUM.sound_speed_m_s;
    let pin = WallPin {
        surface_density: 0.6,
        stiffness_per_area: 6.0,
        resistance: 0.4,
    };
    // Frequencies straddle the original wall's sqrt(K/sigma) resonance.
    for omega in [2.0, 5.0] {
        // WallPin uses exp(-i omega t); FEEC uses exp(+i omega t).
        let impedance = fs_phs::wall_specific_impedance(&pin, omega).unwrap().conj();
        let admittance = impedance.recip();
        let harmonic = model
            .prepare(
                omega / core::f64::consts::TAU,
                &[
                    BoundaryAdmittance {
                        faces: inlet.clone(),
                        admittance_m_pa_s: C64::from_re(1.0 / z0),
                    },
                    BoundaryAdmittance {
                        faces: wall_faces.clone(),
                        admittance_m_pa_s: admittance,
                    },
                ],
                HarmonicLimits::default(),
                &mut || false,
            )
            .unwrap();
        let drive = inlet
            .iter()
            .map(|&face| (face, [C64::from_re(2.0 / z0); 3]))
            .collect::<Vec<_>>();
        let field = harmonic.solve(&drive, &mut || false).unwrap();
        let pressure = field.pressure_pa.iter().map(|p| p.re).collect::<Vec<_>>();
        let potential = field
            .pressure_pa
            .iter()
            .map(|p| p.im / omega)
            .collect::<Vec<_>>();
        let (wall_load, wall_area) = geometry.boundary_load(&wall_faces).unwrap();
        let mean = |p: &[C64]| {
            p.iter()
                .zip(&wall_load)
                .fold(C64::ZERO, |sum, (&p, &b)| sum + p.scale(b / wall_area))
        };
        let phase = C64::new(
            det::cos(omega / MEDIUM.sound_speed_m_s),
            det::sin(omega / MEDIUM.sound_speed_m_s),
        );
        let reflection = (impedance - C64::from_re(z0)) / (impedance + C64::from_re(z0));
        assert!((mean(&field.pressure_pa) * phase - C64::ONE - reflection).abs() < 2e-3);
        let duration = core::f64::consts::FRAC_PI_2 / omega;
        let mut errors = Vec::new();
        for count in [8, 16, 32] {
            let h = duration / count as f64;
            let prepared = model
                .prepare_time_with_walls(
                    h,
                    &[ResistiveBoundary {
                        faces: inlet.clone(),
                        admittance_m_pa_s: 1.0 / z0,
                    }],
                    &[MechanicalBoundary {
                        faces: wall_faces.clone(),
                        wall: pin,
                    }],
                    MechanicalTimeLimits::default(),
                    &mut || false,
                )
                .unwrap();
            let v = prepared
                .wall_dofs(0)
                .unwrap()
                .iter()
                .map(|&i| field.pressure_pa[i] * admittance)
                .collect::<Vec<_>>();
            let displacement = v.iter().map(|v| v.im / omega).collect::<Vec<_>>();
            let velocity = v.iter().map(|v| v.re).collect::<Vec<_>>();
            let mut state = prepared
                .state_with_walls(&potential, &pressure, &[displacement], &[velocity])
                .unwrap();
            for index in 0..count {
                let inward = 2.0 / z0 * det::cos(omega * (index as f64 + 0.5) * h);
                let drive = inlet
                    .iter()
                    .map(|&face| (face, [inward; 3]))
                    .collect::<Vec<_>>();
                state = prepared.step(&state, &drive, &mut || false).unwrap().state;
            }
            // At pi/2 the field equals -Im(P). Together with the initial real
            // field this measures the complex reflected wave, not just energy.
            let reconstructed = pressure
                .iter()
                .zip(state.pressure_pa())
                .map(|(&re, &minus_im)| C64::new(re, -minus_im))
                .collect::<Vec<_>>();
            let error = state
                .pressure_pa()
                .iter()
                .zip(&field.pressure_pa)
                .map(|(&p, &reference)| (p + reference.im).abs())
                .fold(0.0_f64, f64::max);
            println!("mechanical wall phase omega={omega} steps={count} error={error:e}");
            if count == 32 {
                assert!((mean(&reconstructed) * phase - C64::ONE - reflection).abs() < 3e-3);
            }
            errors.push(error);
        }
        assert!(
            errors[0] > 1e-5 && errors[2] < 2e-3,
            "temporal errors {errors:?}"
        );
        for pair in errors.windows(2) {
            assert!(
                (3.5..4.5).contains(&(pair[0] / pair[1])),
                "temporal errors {errors:?}"
            );
        }
    }
}

#[test]
fn empty_mechanical_boundary_preserves_the_original_sparse_evolution() {
    let (complex, positions) = brick(2, [1.0, 0.2, 0.3]);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let boundary = [ResistiveBoundary {
        faces: plane(&geometry, 1.0),
        admittance_m_pa_s: 0.4,
    }];
    let bare = model
        .prepare_time(
            0.01,
            &boundary,
            HeterogeneousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let wrapped = model
        .prepare_time_with_walls(
            0.01,
            &boundary,
            &[],
            MechanicalTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let potential = p2(&geometry, |point| 0.02 * point[0] * point[0]);
    let pressure = p2(&geometry, |point| 0.03 * (1.0 + point[1]));
    let mut a = bare.state(&potential, &pressure).unwrap();
    let mut b = wrapped.state(&potential, &pressure).unwrap();
    let source = plane(&geometry, 0.0)
        .into_iter()
        .map(|face| {
            let values = face.map(|v| 0.01 * (1.0 + positions[v as usize][1]));
            (face, values)
        })
        .collect::<Vec<_>>();
    for _ in 0..4 {
        let next_a = bare.step(&a, &source, &mut || false).unwrap();
        let next_b = wrapped.step(&b, &source, &mut || false).unwrap();
        assert_eq!(next_a.state.pressure_pa(), next_b.state.pressure_pa());
        assert_eq!(
            next_a.state.pressure_potential_pa_s(),
            next_b.state.pressure_potential_pa_s()
        );
        assert_eq!(
            next_a.work.boundary_dissipation_j,
            next_b.work.boundary_dissipation_j
        );
        assert_eq!(next_b.work.wall_work_j, 0.0);
        assert_eq!(next_b.work.wall_dissipation_j, 0.0);
        assert_eq!(
            wrapped.energy(&next_b.state).unwrap().total_j(),
            bare.energy(&next_a.state).unwrap().total_j()
        );
        a = next_a.state;
        b = next_b.state;
    }
}

#[test]
fn wall_admission_limits_overlap_cancellation_and_state_binding_are_explicit() {
    let (complex, positions) = brick(1, [1.0; 3]);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let boundary = MechanicalBoundary {
        faces: plane(&geometry, 1.0),
        wall: WallPin {
            surface_density: 0.4,
            stiffness_per_area: 3.0,
            resistance: 0.1,
        },
    };
    let prepare = |wall: MechanicalBoundary| {
        model.prepare_time_with_walls(
            0.01,
            &[],
            &[wall],
            MechanicalTimeLimits::default(),
            &mut || false,
        )
    };
    let prepared = prepare(boundary.clone()).unwrap();
    let other = prepare(boundary.clone()).unwrap();
    let state = prepared.zero_state();
    assert!(other.step(&state, &[], &mut || false).is_err());
    assert!(matches!(
        prepared.step(&state, &[], &mut || true),
        Err(HeterogeneousError::Cancelled)
    ));
    assert!(
        state
            .wall_displacement_m(0)
            .unwrap()
            .iter()
            .all(|&v| v == 0.0)
    );
    assert!(
        prepared
            .state_with_walls(
                state.pressure_potential_pa_s(),
                state.pressure_pa(),
                &[],
                &[]
            )
            .is_err()
    );
    assert!(matches!(
        model.prepare_time_with_walls(
            0.01,
            &[],
            &[boundary.clone()],
            MechanicalTimeLimits {
                max_wall_dofs: 1,
                ..MechanicalTimeLimits::default()
            },
            &mut || false
        ),
        Err(HeterogeneousError::Limit {
            what: "distributed wall trace degrees of freedom",
            ..
        })
    ));
    assert!(
        model
            .prepare_time_with_walls(
                0.01,
                &[ResistiveBoundary {
                    faces: boundary.faces.clone(),
                    admittance_m_pa_s: 0.5,
                }],
                &[boundary.clone()],
                MechanicalTimeLimits::default(),
                &mut || false
            )
            .is_err()
    );
    let mut invalid = boundary;
    invalid.wall.surface_density = 0.0;
    assert!(matches!(prepare(invalid), Err(HeterogeneousError::Wall(_))));
}
