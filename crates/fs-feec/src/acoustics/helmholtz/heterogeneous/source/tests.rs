//! G1/G3 distributed-source checks against continuum continuity and work.

use super::*;
use crate::acoustics::helmholtz::heterogeneous::time::HeterogeneousTimeLimits;
use crate::acoustics::helmholtz::heterogeneous::time::exact::{ExactTimeLimits, ExactVelocityPort};
use crate::acoustics::helmholtz::heterogeneous::{BoundaryAdmittance, HarmonicLimits};
use crate::acoustics::helmholtz::{HelmholtzFe, Medium};
use crate::fixtures::masked_cube_grid;
use fs_math::det;

const MEDIUM: Medium = Medium {
    density_kg_m3: 1.3,
    sound_speed_m_s: 2.0,
};
const LENGTHS: [f64; 3] = [1.0, 0.4, 0.2];
const VOLUME: f64 = LENGTHS[0] * LENGTHS[1] * LENGTHS[2];

fn geometry(nx: usize) -> (fs_rep_mesh::TetComplex, Vec<[f64; 3]>) {
    let (mesh, mut positions) = masked_cube_grid(nx, 1, 1, &|_, _, _| true);
    for point in &mut positions {
        point[0] *= LENGTHS[0] / nx as f64;
        point[1] *= LENGTHS[1];
        point[2] *= LENGTHS[2];
    }
    (mesh, positions)
}

fn end_faces(fe: &HelmholtzFe<'_>, x: f64) -> Vec<[u32; 3]> {
    fe.boundary_faces()
        .into_iter()
        .filter(|face| {
            face.iter()
                .all(|&v| (fe.positions()[v as usize][0] - x).abs() < 1e-12)
        })
        .collect()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn distributed_monopole_matches_loaded_continuum_duct_and_physical_power() {
    // Uniform s=Q/V between TWO anechoic ends. The independent continuum
    // solution is p=P0[1-exp(-i theta)cos(k(x-L/2))], P0=-i rho c^2 s/omega.
    // It tests spatial source phase, propagation and work, not a copied FE RHS.
    let omega = 1.3 * MEDIUM.sound_speed_m_s / LENGTHS[0];
    let k = omega / MEDIUM.sound_speed_m_s;
    let theta = k * LENGTHS[0] / 2.0;
    let travel = C64::new(det::cos(theta), -det::sin(theta));
    let q = C64::new(0.002, -0.0007);
    let p0 = q * C64::new(
        0.0,
        -MEDIUM.density_kg_m3 * MEDIUM.sound_speed_m_s.powi(2) / (VOLUME * omega),
    );
    let exact = |x: f64| p0 * (C64::ONE - travel.scale(det::cos(k * (x - LENGTHS[0] / 2.0))));
    let mean = p0 * (C64::ONE - travel.scale(det::sin(theta) / theta));
    let expected_power = 0.5 * (mean.conj() * q).re;
    let radiated = LENGTHS[1] * LENGTHS[2] * exact(0.0).norm_sq()
        / (MEDIUM.density_kg_m3 * MEDIUM.sound_speed_m_s);
    close(expected_power, radiated, expected_power * 1e-13);
    let mut errors = Vec::new();
    for nx in [1, 2, 4] {
        let (mesh, positions) = geometry(nx);
        let fe = HelmholtzFe::new(&mesh, &positions, 2, MEDIUM).unwrap();
        let model = HeterogeneousHelmholtzFe::new(
            &fe,
            4,
            &mut |_, _| Ok::<_, crate::acoustics::helmholtz::AcousticFeError>(MEDIUM),
            &mut || false,
        )
        .unwrap();
        let source = model
            .volume_velocity_source(
                4,
                &mut |_, _| Ok::<_, HeterogeneousError>(1.0 / VOLUME),
                &mut || false,
            )
            .unwrap();
        close(source.integrated_shape(), 1.0, 2e-14);
        let mut faces = end_faces(&fe, 0.0);
        faces.extend(end_faces(&fe, LENGTHS[0]));
        let owner = model
            .prepare(
                omega / core::f64::consts::TAU,
                &[BoundaryAdmittance {
                    faces,
                    admittance_m_pa_s: C64::from_re(
                        1.0 / (MEDIUM.density_kg_m3 * MEDIUM.sound_speed_m_s),
                    ),
                }],
                HarmonicLimits::default(),
                &mut || false,
            )
            .unwrap();
        let field = owner
            .solve_with_volume_sources(&[], &[(&source, q)], &mut || false)
            .unwrap();
        let mut error: f64 = 0.0;
        for (element, tet) in fe.tetrahedra().iter().enumerate() {
            let x = tet.iter().map(|&v| positions[v as usize][0]).sum::<f64>() / 4.0;
            let pressure = model
                .pressure_at(&field.pressure_pa, element, [0.25; 4])
                .unwrap();
            error = error.max((pressure - exact(x)).abs() / p0.abs());
        }
        close(
            field.input_power_w,
            field.volume_source_power_w,
            expected_power * 1e-12,
        );
        close(
            field.power_balance_residual_w(),
            0.0,
            expected_power * 1e-10,
        );
        if nx == 4 {
            close(
                field.volume_source_power_w,
                expected_power,
                expected_power * 1e-4,
            );
            assert!(error < 1e-4, "continuum field error={error:e}");
        }
        errors.push(error);
        // Coherent copies must cancel in the original field before power.
        let dark = owner
            .solve_with_volume_sources(&[], &[(&source, q), (&source, -q)], &mut || false)
            .unwrap();
        assert!(dark.pressure_pa.iter().all(|&p| p == C64::ZERO));
        assert_eq!(dark.volume_source_power_w, 0.0);
    }
    // These coarse meshes include pre-asymptotic error, so this ladder does
    // not establish a uniform refinement rate. Require monotonic improvement
    // while retaining the independently anchored fine-grid accuracy above.
    assert!(
        errors[1] < errors[0] && errors[2] < errors[1],
        "h ladder={errors:?}"
    );
}

#[test]
fn sealed_heterogeneous_volume_source_compresses_without_boundary_motion() {
    let (mesh, positions) = geometry(2);
    // Equal c, different density: a source proportional to local beta gives
    // exactly uniform p(t)=q*Kref*t/V even across the material interface.
    let material = |x: [f64; 3]| Medium {
        density_kg_m3: if x[0] < 0.5 { 1.3 } else { 3.9 },
        sound_speed_m_s: 2.0,
    };
    let fe = HelmholtzFe::new(&mesh, &positions, 2, MEDIUM).unwrap();
    let model = HeterogeneousHelmholtzFe::new(
        &fe,
        4,
        &mut |_, x| Ok::<_, crate::acoustics::helmholtz::AcousticFeError>(material(x)),
        &mut || false,
    )
    .unwrap();
    let bulk = MEDIUM.density_kg_m3 * MEDIUM.sound_speed_m_s.powi(2);
    let source = model
        .volume_velocity_source(
            4,
            &mut |_, x| {
                let medium = material(x);
                Ok::<_, HeterogeneousError>(
                    bulk / (VOLUME * medium.density_kg_m3 * medium.sound_speed_m_s.powi(2)),
                )
            },
            &mut || false,
        )
        .unwrap();
    close(source.integrated_shape(), 2.0 / 3.0, 2e-14);
    let dt = 0.03;
    let owner = model
        .prepare_time(dt, &[], HeterogeneousTimeLimits::default(), &mut || false)
        .unwrap();
    for q in [0.002, -0.002, 2e-10] {
        let mut state = owner.zero_state();
        let mut work = 0.0;
        for step in 1..=5 {
            let result = owner
                .step_with_volume_sources(&state, &[], &[(&source, q)], &mut || false)
                .unwrap();
            let time = f64::from(step) * dt;
            let pressure = q * bulk * time / VOLUME;
            for element in 0..fe.tetrahedra().len() {
                let probe = model
                    .probe(element, [0.25; 4], &mut |_, x| {
                        Ok::<_, crate::acoustics::helmholtz::AcousticFeError>(material(x))
                    })
                    .unwrap();
                let field = probe.transient(&result.state).unwrap();
                close(field.pressure_pa, pressure, pressure.abs() * 2e-10);
                for velocity in field.particle_velocity_m_s {
                    close(velocity, 0.0, pressure.abs() * 2e-10);
                }
            }
            close(
                result.work.volume_source_work_j,
                result.work.input_work_j,
                result.work.input_work_j.abs() * 1e-12,
            );
            work += result.work.volume_source_work_j;
            let energy = 0.5 * pressure.powi(2) * VOLUME / bulk * source.integrated_shape();
            close(work, energy, energy * 2e-10);
            close(result.work.boundary_dissipation_j, 0.0, energy * 1e-12);
            state = result.state;
        }
        // Reversing injection extracts the stored compression energy.
        let reversal = owner
            .step_with_volume_sources(&state, &[], &[(&source, -q)], &mut || false)
            .unwrap();
        assert!(reversal.work.volume_source_work_j < 0.0);
    }
}

#[test]
fn volume_sources_share_original_boundary_drive_and_do_not_requery_callbacks() {
    let (mesh, positions) = geometry(1);
    let fe = HelmholtzFe::new(&mesh, &positions, 2, MEDIUM).unwrap();
    let model = HeterogeneousHelmholtzFe::new(
        &fe,
        4,
        &mut |_, _| Ok::<_, crate::acoustics::helmholtz::AcousticFeError>(MEDIUM),
        &mut || false,
    )
    .unwrap();
    let mut samples = 0;
    let source = model
        .volume_velocity_source(
            4,
            &mut |_, _| {
                samples += 1;
                Ok::<_, HeterogeneousError>(1.0 / VOLUME)
            },
            &mut || false,
        )
        .unwrap();
    let sampled = samples;
    let omega = 3.0;
    let boundary = end_faces(&fe, LENGTHS[0]);
    let y = 1.0 / (MEDIUM.density_kg_m3 * MEDIUM.sound_speed_m_s);
    let prepared = model
        .prepare(
            omega / core::f64::consts::TAU,
            &[BoundaryAdmittance {
                faces: boundary.clone(),
                admittance_m_pa_s: C64::from_re(y),
            }],
            HarmonicLimits::default(),
            &mut || false,
        )
        .unwrap();
    let p = C64::new(0.03, -0.02);
    let q = p * C64::new(
        0.0,
        omega * VOLUME / (MEDIUM.density_kg_m3 * MEDIUM.sound_speed_m_s.powi(2)),
    );
    let velocity: Vec<_> = boundary
        .iter()
        .map(|&face| (face, [p.scale(y); 3]))
        .collect();
    let field = prepared
        .solve_with_volume_sources(&velocity, &[(&source, q)], &mut || false)
        .unwrap();
    for element in 0..fe.tetrahedra().len() {
        assert!(
            (model
                .pressure_at(&field.pressure_pa, element, [0.25; 4])
                .unwrap()
                - p)
                .abs()
                < p.abs() * 1e-11
        );
    }
    let absorbed = 0.5 * LENGTHS[1] * LENGTHS[2] * y * p.norm_sq();
    close(field.input_power_w, absorbed, absorbed * 1e-11);
    close(field.volume_source_power_w, 0.0, absorbed * 1e-11);
    assert_eq!(samples, sampled);
    let legacy = prepared.solve(&velocity, &mut || false).unwrap();
    let empty = prepared
        .solve_with_volume_sources(&velocity, &[], &mut || false)
        .unwrap();
    assert_eq!(legacy.pressure_pa, empty.pressure_pa);
    assert_eq!(
        legacy.input_power_w.to_bits(),
        empty.input_power_w.to_bits()
    );
}

#[test]
fn volume_source_admission_and_reuse_are_transactional() {
    let (mesh, positions) = geometry(1);
    let fe = HelmholtzFe::new(&mesh, &positions, 2, MEDIUM).unwrap();
    let build = || {
        HeterogeneousHelmholtzFe::new(
            &fe,
            4,
            &mut |_, _| Ok::<_, crate::acoustics::helmholtz::AcousticFeError>(MEDIUM),
            &mut || false,
        )
        .unwrap()
    };
    let model = build();
    let foreign = build();
    let source = model
        .volume_velocity_source(
            4,
            &mut |_, _| Ok::<_, HeterogeneousError>(1.0 / VOLUME),
            &mut || false,
        )
        .unwrap();
    let foreign_source = foreign
        .volume_velocity_source(
            4,
            &mut |_, _| Ok::<_, HeterogeneousError>(1.0 / VOLUME),
            &mut || false,
        )
        .unwrap();
    let mut calls = 0;
    assert!(
        model
            .volume_velocity_source(
                3,
                &mut |_, _| {
                    calls += 1;
                    Ok::<_, HeterogeneousError>(1.0)
                },
                &mut || false
            )
            .is_err()
    );
    assert_eq!(calls, 0);
    assert!(
        model
            .volume_velocity_source(
                4,
                &mut |_, _| Ok::<_, HeterogeneousError>(f64::NAN),
                &mut || false
            )
            .is_err()
    );
    let named = HeterogeneousError::Limit {
        what: "authored source refused",
        requested: 2,
        maximum: 1,
    };
    assert!(matches!(model.volume_velocity_source(4,
        &mut |_, _| Err::<f64, _>(named.clone()), &mut || false), Err(error) if error == named));
    assert!(matches!(
        model.volume_velocity_source(4, &mut |_, _| Ok::<_, HeterogeneousError>(1.0), &mut || {
            true
        }),
        Err(HeterogeneousError::Cancelled)
    ));
    let owner = model
        .prepare_time(0.01, &[], HeterogeneousTimeLimits::default(), &mut || false)
        .unwrap();
    let state = owner.zero_state();
    for invalid in [&[(&foreign_source, 0.1)][..], &[(&source, f64::NAN)][..]] {
        assert!(
            owner
                .step_with_volume_sources(&state, &[], invalid, &mut || false)
                .is_err()
        );
    }
    assert!(matches!(
        owner.step_with_volume_sources(&state, &[], &[(&source, 0.1)], &mut || true),
        Err(HeterogeneousError::Cancelled)
    ));
    assert!(state.pressure_pa().iter().all(|&p| p == 0.0));
    let a = owner
        .step_with_volume_sources(&state, &[], &[(&source, 0.1)], &mut || false)
        .unwrap();
    let b = owner
        .step_with_volume_sources(&state, &[], &[(&source, 0.1)], &mut || false)
        .unwrap();
    assert_eq!(a.state.pressure_pa(), b.state.pressure_pa());
    assert_eq!(a.work.input_work_j.to_bits(), b.work.input_work_j.to_bits());
}

#[test]
fn exact_volume_injection_preserves_si_compression_and_integrated_work() {
    let (mesh, positions) = geometry(1);
    for medium in [
        Medium {
            density_kg_m3: 1.2,
            sound_speed_m_s: 343.0,
        },
        Medium {
            density_kg_m3: 1000.0,
            sound_speed_m_s: 1480.0,
        },
    ] {
        let fe = HelmholtzFe::new(&mesh, &positions, 2, medium).unwrap();
        let model = HeterogeneousHelmholtzFe::new(
            &fe,
            4,
            &mut |_, _| Ok::<_, crate::acoustics::helmholtz::AcousticFeError>(medium),
            &mut || false,
        )
        .unwrap();
        let source = model
            .volume_velocity_source(
                4,
                &mut |_, _| Ok::<_, HeterogeneousError>(1.0 / VOLUME),
                &mut || false,
            )
            .unwrap();
        let dt = 0.05 / medium.sound_speed_m_s;
        let owner = model
            .prepare_exact_time_with_volume_sources(
                dt,
                &[],
                &[],
                &[&source],
                ExactTimeLimits::default(),
                &mut || false,
            )
            .unwrap();
        assert_eq!(owner.port_count(), 0);
        assert_eq!(owner.volume_source_count(), 1);
        for q in [1e-6, 1e-12] {
            let mut state = owner.zero_state();
            let rate = q * medium.density_kg_m3 * medium.sound_speed_m_s.powi(2) / VOLUME;
            for step in 1..=3 {
                let result = owner
                    .step_with_volume_sources(&state, &[], &[q], &mut || false)
                    .unwrap();
                let time = f64::from(step) * dt;
                let pressure = rate * time;
                for element in 0..fe.tetrahedra().len() {
                    let probe = model
                        .probe(element, [0.25; 4], &mut |_, _| {
                            Ok::<_, crate::acoustics::helmholtz::AcousticFeError>(medium)
                        })
                        .unwrap();
                    let field = probe.transient(&result.state).unwrap();
                    close(field.pressure_pa, pressure, pressure.abs() * 2e-9);
                }
                let integral = rate * dt * (time - dt / 2.0);
                close(
                    result.integrated_volume_source_pressure_pa_s[0],
                    integral,
                    integral.abs() * 2e-9,
                );
                close(
                    result.work.volume_source_work_j,
                    q * integral,
                    (q * integral).abs() * 2e-9,
                );
                close(
                    result.work.balance_residual_j(),
                    0.0,
                    (q * integral).abs() * 2e-9,
                );
                assert!(result.integrated_port_pressure_pa_m2_s.is_empty());
                state = result.state;
            }
        }
        let state = owner.zero_state();
        assert!(owner.step(&state, &[], &mut || false).is_err());
        assert!(
            owner
                .step_with_volume_sources(&state, &[], &[f64::NAN], &mut || false)
                .is_err()
        );
        assert!(matches!(
            owner.step_with_volume_sources(&state, &[], &[0.1], &mut || true),
            Err(HeterogeneousError::Cancelled)
        ));
    }
}

#[test]
fn exact_joint_boundary_and_volume_forcing_retains_signed_source_work() {
    let (mesh, positions) = geometry(1);
    let fe = HelmholtzFe::new(&mesh, &positions, 2, MEDIUM).unwrap();
    let model = HeterogeneousHelmholtzFe::new(
        &fe,
        4,
        &mut |_, _| Ok::<_, crate::acoustics::helmholtz::AcousticFeError>(MEDIUM),
        &mut || false,
    )
    .unwrap();
    let source = model
        .volume_velocity_source(
            4,
            &mut |_, _| Ok::<_, HeterogeneousError>(1.0 / VOLUME),
            &mut || false,
        )
        .unwrap();
    let port = ExactVelocityPort {
        profile: end_faces(&fe, LENGTHS[0])
            .into_iter()
            .map(|face| (face, [1.0; 3]))
            .collect(),
    };
    let dt = 0.025;
    let owner = model
        .prepare_exact_time_with_volume_sources(
            dt,
            &[],
            &[port],
            &[&source],
            ExactTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let a = 0.02;
    // psi=a*x^2/2 + alpha*t^2/2, p=alpha*t. Volume extraction opposes
    // the inward right-wall drive without changing the uniform velocity slope.
    let q = -0.5 * a * VOLUME / MEDIUM.density_kg_m3;
    let alpha = 0.5 * a * MEDIUM.sound_speed_m_s.powi(2);
    let mut potential = vec![0.0; fe.ndof()];
    for (element, tet) in fe.tetrahedra().iter().enumerate() {
        for (dof, shape) in fe.space.element_dofs(element) {
            potential[dof] = match shape {
                crate::highorder::simplex::LocalFn::Vertex(v) => {
                    0.5 * a * positions[tet[v] as usize][0].powi(2)
                }
                crate::highorder::simplex::LocalFn::Edge(v, w, 0) => {
                    -0.5 * a
                        * (positions[tet[v] as usize][0] - positions[tet[w] as usize][0]).powi(2)
                }
                _ => panic!("P2 source fixture"),
            };
        }
    }
    let mut state = owner.state(&potential, &vec![0.0; fe.ndof()]).unwrap();
    for step in 1..=4 {
        let result = owner
            .step_with_volume_sources(
                &state,
                &[a * LENGTHS[0] / MEDIUM.density_kg_m3],
                &[q],
                &mut || false,
            )
            .unwrap();
        let time = f64::from(step) * dt;
        let integral = alpha * dt * (time - 0.5 * dt);
        close(
            result.integrated_volume_source_pressure_pa_s[0],
            integral,
            integral * 2e-9,
        );
        close(
            result.integrated_port_pressure_pa_m2_s[0],
            integral * LENGTHS[1] * LENGTHS[2],
            integral * LENGTHS[1] * LENGTHS[2] * 2e-9,
        );
        close(
            result.work.volume_source_work_j,
            q * integral,
            (q * integral).abs() * 2e-9,
        );
        assert!(result.work.volume_source_work_j < 0.0);
        assert!(result.work.input_work_j > 0.0);
        close(
            result.work.balance_residual_j(),
            0.0,
            result.work.input_work_j * 2e-9,
        );
        for element in 0..fe.tetrahedra().len() {
            let probe = model
                .probe(element, [0.25; 4], &mut |_, _| {
                    Ok::<_, crate::acoustics::helmholtz::AcousticFeError>(MEDIUM)
                })
                .unwrap();
            close(
                probe.transient(&result.state).unwrap().pressure_pa,
                alpha * time,
                alpha * time * 2e-9,
            );
        }
        state = result.state;
    }
}
