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
                    relaxation: Vec::new(),
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
                    relaxation: Vec::new(),
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
    reflection_phase(&[]);
}

#[test]
fn hereditary_wall_reflection_retains_both_phase_and_amplitude_under_time_refinement() {
    reflection_phase(&[
        AcousticWallRelaxation {
            stiffness_per_area: 3.0,
            relaxation_time_s: 0.08,
        },
        AcousticWallRelaxation {
            stiffness_per_area: 2.0,
            relaxation_time_s: 0.6,
        },
    ]);
}

// The order of the time error is read on the rms over ALL dofs, not on the max. The dof where the max sits hops from count
// to count (omega = 2, hereditary, counts 8 to 256: dofs 49, 61, 123, 67, 49, 63, with signs +, -, +, -, -, +), so adjacent
// ratios of the max scatter between 3.1 and 15.5 and are not an order measurement; the rms ratios from 16 steps on lie in
// [3.64, 4.58] for both omegas, with and without memory arms (measured). A second-order scheme halves the step to a quarter of
// the error: the band is 4 +- 20 %, and a first-order variant must fall outside it.
const ORDER_BAND: core::ops::Range<f64> = 3.4..4.8;

fn reflection_phase(relaxation: &[AcousticWallRelaxation]) {
    for (omega, max, rms) in reflection_errors(relaxation, &[2.0, 5.0], 0.5, &[16, 32, 64, 128]) {
        assert!(
            max[0] > 1e-5 && max[2] < 2e-3,
            "omega {omega}: temporal max errors {max:?}"
        );
        for pair in rms.windows(2) {
            assert!(
                ORDER_BAND.contains(&(pair[0] / pair[1])),
                "omega {omega}: rms errors {rms:?}, max errors {max:?}"
            );
        }
    }
}

// Planted negative: the same predicate on a first-order variant, the inlet drive held at the START of each step instead of
// its midpoint, must fall outside the band at every pair.
#[test]
fn a_drive_held_at_the_step_start_is_first_order_and_the_rms_order_band_refuses_it() {
    let results = reflection_errors(&[], &[2.0], 0.0, &[16, 32, 64, 128]);
    let (_, max, rms) = &results[0];
    for pair in rms.windows(2) {
        assert!(
            !ORDER_BAND.contains(&(pair[0] / pair[1])),
            "rms errors {rms:?}, max errors {max:?}"
        );
    }
}

// For each omega: the max and the rms over all dofs of |p + Im(reference)| at each count. `drive_sample` places the held
// inlet drive of each step inside the step: 0.5 is the midpoint (second order), 0 the step start (first order).
fn reflection_errors(
    relaxation: &[AcousticWallRelaxation],
    omegas: &[f64],
    drive_sample: f64,
    counts: &[usize],
) -> Vec<(f64, Vec<f64>, Vec<f64>)> {
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
    let mut results = Vec::new();
    for &omega in omegas {
        // WallPin uses exp(-i omega t); FEEC uses exp(+i omega t).
        let mut impedance = fs_phs::wall_specific_impedance(&pin, omega).unwrap().conj();
        // Independent continuous Maxwell transform: each arm contributes
        // k*tau/(1+i*omega*tau) to pressure/velocity impedance.
        for arm in relaxation {
            impedance = impedance
                + C64::from_re(arm.stiffness_per_area * arm.relaxation_time_s)
                    / C64::new(1.0, omega * arm.relaxation_time_s);
        }
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
        let mut rms_errors = Vec::new();
        for &count in counts {
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
                        relaxation: relaxation.to_vec(),
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
            let memories = relaxation
                .iter()
                .map(|arm| {
                    v.iter()
                        .map(|v| {
                            let w = C64::new(v.im / omega, -v.re / omega);
                            (w.scale(det::sqrt(arm.stiffness_per_area))
                                / C64::new(1.0, omega * arm.relaxation_time_s))
                            .re
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            let mut state = prepared
                .state_with_wall_memory(
                    &potential,
                    &pressure,
                    &[displacement],
                    &[velocity],
                    &[memories],
                )
                .unwrap();
            for index in 0..count {
                let inward = 2.0 / z0 * det::cos(omega * (index as f64 + drive_sample) * h);
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
            let rms = (state
                .pressure_pa()
                .iter()
                .zip(&field.pressure_pa)
                .map(|(&p, &reference)| (p + reference.im).powi(2))
                .sum::<f64>()
                / field.pressure_pa.len() as f64)
                .sqrt();
            println!(
                "mechanical wall phase omega={omega} steps={count} error={error:e} rms={rms:e}"
            );
            if count == 32 && drive_sample == 0.5 {
                assert!((mean(&reconstructed) * phase - C64::ONE - reflection).abs() < 3e-3);
            }
            errors.push(error);
            rms_errors.push(rms);
        }
        results.push((omega, errors, rms_errors));
    }
    results
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
        relaxation: Vec::new(),
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

#[test]
fn maxwell_wall_ramp_has_independent_storage_heat_and_area_scaling() {
    let pin = WallPin {
        surface_density: 0.4,
        stiffness_per_area: 1.7,
        resistance: 0.3,
    };
    let arms = vec![
        AcousticWallRelaxation {
            stiffness_per_area: 2.1,
            relaxation_time_s: 0.03,
        },
        AcousticWallRelaxation {
            stiffness_per_area: 1.4,
            relaxation_time_s: 0.6,
        },
    ];
    let viscosity: f64 = arms
        .iter()
        .map(|a| a.stiffness_per_area * a.relaxation_time_s)
        .sum();
    let h = 0.0125;
    for height in [0.1, 0.3] {
        let area = height * 0.2;
        let (complex, positions) = brick(2, [1.0, height, 0.2]);
        let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
        let model = uniform(&geometry);
        let prepared = model
            .prepare_time_with_walls(
                h,
                &[],
                &[MechanicalBoundary {
                    faces: plane(&geometry, 1.0),
                    wall: pin,
                    relaxation: arms.clone(),
                }],
                MechanicalTimeLimits::default(),
                &mut || false,
            )
            .unwrap();
        for amplitude in [1.0, 1e-8] {
            let v = 0.02 * amplitude;
            let w0 = 0.03 * amplitude;
            // Continuous exact solution: each arm already follows a steady
            // strain ramp, e_j=sqrt(k_j)*tau_j*v. Memory grows linearly;
            // physical branch storage is constant while it dissipates k*tau*v².
            let a = pin.stiffness_per_area * v / MEDIUM.sound_speed_m_s.powi(2);
            let b = -MEDIUM.density_kg_m3 * v - a;
            let p0 = (pin.resistance + viscosity) * v + pin.stiffness_per_area * w0;
            let potential = p2(&geometry, |point| 0.5 * a * point[0].powi(2) + b * point[0]);
            let pressure = p2(&geometry, |_| p0);
            let displacement = restricted(&prepared, &p2(&geometry, |_| w0));
            let velocity = restricted(&prepared, &p2(&geometry, |_| v));
            let memories = arms
                .iter()
                .map(|arm| {
                    restricted(
                        &prepared,
                        &p2(&geometry, |_| {
                            det::sqrt(arm.stiffness_per_area) * (w0 - arm.relaxation_time_s * v)
                        }),
                    )
                })
                .collect::<Vec<_>>();
            let mut state = prepared
                .state_with_wall_memory(
                    &potential,
                    &pressure,
                    &[displacement],
                    &[velocity],
                    &[memories],
                )
                .unwrap();
            let expected_memory: f64 = arms
                .iter()
                .map(|arm| {
                    0.5 * area * arm.stiffness_per_area * (arm.relaxation_time_s * v).powi(2)
                })
                .sum();
            let initial = prepared.energy(&state).unwrap();
            assert!((initial.wall_relaxation_j / expected_memory - 1.0).abs() < 2e-12);
            let source = plane(&geometry, 0.0)
                .into_iter()
                .map(|face| (face, [-b / MEDIUM.density_kg_m3; 3]))
                .collect::<Vec<_>>();
            let mut work = 0.0;
            let mut heat = 0.0;
            for step_index in 1..=6 {
                let next = prepared.step(&state, &source, &mut || false).unwrap();
                let time = step_index as f64 * h;
                let w = w0 + v * time;
                let p = (pin.resistance + viscosity) * v + pin.stiffness_per_area * w;
                for &actual in &next.state.pressure_pa()[..positions.len()] {
                    assert!((actual / p - 1.0).abs() < 5e-10);
                }
                for (j, arm) in arms.iter().enumerate() {
                    let expected = restricted(
                        &prepared,
                        &p2(&geometry, |_| {
                            det::sqrt(arm.stiffness_per_area) * (w - arm.relaxation_time_s * v)
                        }),
                    );
                    for (&actual, expected) in next
                        .state
                        .wall_relaxation_memory(0, j)
                        .unwrap()
                        .iter()
                        .zip(expected)
                    {
                        assert!((actual - expected).abs() < 1e-10 * amplitude);
                    }
                }
                let energy = prepared.energy(&next.state).unwrap();
                assert!((energy.wall_relaxation_j / expected_memory - 1.0).abs() < 2e-9);
                let maxwell_heat = h * area * viscosity * v * v;
                assert!(
                    (next.work.wall_relaxation_dissipation_j / maxwell_heat - 1.0).abs() < 2e-9
                );
                assert!(
                    (next.work.wall_dissipation_j
                        / (h * area * (pin.resistance + viscosity) * v * v)
                        - 1.0)
                        .abs()
                        < 2e-9
                );
                work += next.work.input_work_j;
                heat += next.work.wall_dissipation_j;
                assert!(
                    (energy.total_j() + heat - initial.total_j() - work).abs()
                        < 2e-9 * initial.total_j()
                );
                state = next.state;
            }
        }
    }
}

#[test]
fn nonuniform_wall_memory_can_release_energy_at_zero_displacement() {
    let lengths = [1.0, 0.2, 0.3];
    let area = lengths[1] * lengths[2];
    let (complex, positions) = brick(2, lengths);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let k = 2.5;
    let tau = 0.02;
    let h = 0.005;
    let prepared = model
        .prepare_time_with_walls(
            h,
            &[],
            &[MechanicalBoundary {
                faces: plane(&geometry, 1.0),
                wall: WallPin {
                    surface_density: 0.6,
                    stiffness_per_area: 3.0,
                    resistance: 0.0,
                },
                relaxation: vec![AcousticWallRelaxation {
                    stiffness_per_area: k,
                    relaxation_time_s: tau,
                }],
            }],
            MechanicalTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let amplitude = 0.01;
    let memory = restricted(
        &prepared,
        &p2(&geometry, |point| {
            det::sqrt(k) * amplitude * (1.0 + point[1] / lengths[1]).powi(2)
        }),
    );
    let zeros = vec![0.0; geometry.ndof()];
    let wall_zeros = vec![0.0; prepared.wall_dofs(0).unwrap().len()];
    let mut state = prepared
        .state_with_wall_memory(
            &zeros,
            &zeros,
            &[wall_zeros.clone()],
            &[wall_zeros],
            &[vec![memory]],
        )
        .unwrap();
    let initial = prepared.energy(&state).unwrap();
    let expected = 0.5 * k * area * amplitude * amplitude * 31.0 / 5.0;
    assert!((initial.wall_relaxation_j / expected - 1.0).abs() < 3e-12);
    assert_eq!(initial.wall_elastic_j, 0.0);
    assert_eq!(initial.wall_kinetic_j, 0.0);
    assert_eq!(initial.acoustic.total_j(), 0.0);
    let mut heat = 0.0;
    for _ in 0..8 {
        let step = prepared.step(&state, &[], &mut || false).unwrap();
        assert_eq!(step.work.input_work_j, 0.0);
        assert_eq!(step.work.boundary_dissipation_j, 0.0);
        assert_eq!(
            step.work.wall_dissipation_j,
            step.work.wall_relaxation_dissipation_j
        );
        assert!(step.work.wall_relaxation_dissipation_j > 0.0);
        assert!(step.work.wall_balance_residual_j().abs() < 2e-10 * expected);
        assert!(step.work.acoustic_balance_residual_j().abs() < 2e-10 * expected);
        heat += step.work.wall_dissipation_j;
        state = step.state;
    }
    let energy = prepared.energy(&state).unwrap();
    assert!(energy.acoustic.total_j() > 0.0);
    assert!(energy.wall_relaxation_j < expected);
    assert!((energy.total_j() + heat - expected).abs() < 2e-10 * expected);
}

#[test]
fn hereditary_walls_require_complete_valid_memory_and_an_explicit_budget() {
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
        relaxation: vec![
            AcousticWallRelaxation {
                stiffness_per_area: 2.0,
                relaxation_time_s: 0.1
            };
            2
        ],
    };
    let prepared = model
        .prepare_time_with_walls(
            0.01,
            &[],
            &[boundary.clone()],
            MechanicalTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let count = prepared.wall_dofs(0).unwrap().len();
    let zeros = vec![0.0; geometry.ndof()];
    let wall_zeros = vec![vec![0.0; count]];
    assert!(
        prepared
            .state_with_walls(&zeros, &zeros, &wall_zeros, &wall_zeros)
            .is_err()
    );
    assert!(
        prepared
            .state_with_wall_memory(
                &zeros,
                &zeros,
                &wall_zeros,
                &wall_zeros,
                &[vec![vec![0.0; count]]]
            )
            .is_err()
    );
    let memory = vec![vec![vec![0.01; count]; 2]];
    let state = prepared
        .state_with_wall_memory(&zeros, &zeros, &wall_zeros, &wall_zeros, &memory)
        .unwrap();
    assert!(matches!(
        prepared.step(&state, &[], &mut || true),
        Err(HeterogeneousError::Cancelled)
    ));
    assert_eq!(state.wall_relaxation_memory(0, 0).unwrap(), memory[0][0]);
    for maximum in [0, 2 * count - 1, 2 * count] {
        let result = model.prepare_time_with_walls(
            0.01,
            &[],
            &[boundary.clone()],
            MechanicalTimeLimits {
                max_wall_relaxation_dofs: maximum,
                ..MechanicalTimeLimits::default()
            },
            &mut || false,
        );
        assert_eq!(result.is_ok(), maximum == 2 * count);
    }
    for (k, tau) in [
        (0.0, 0.1),
        (-1.0, 0.1),
        (1.0, 0.0),
        (f64::NAN, 0.1),
        (1.0, f64::INFINITY),
    ] {
        let mut invalid = boundary.clone();
        invalid.relaxation[0] = AcousticWallRelaxation {
            stiffness_per_area: k,
            relaxation_time_s: tau,
        };
        assert!(
            model
                .prepare_time_with_walls(
                    0.01,
                    &[],
                    &[invalid],
                    MechanicalTimeLimits::default(),
                    &mut || false
                )
                .is_err()
        );
    }
}
