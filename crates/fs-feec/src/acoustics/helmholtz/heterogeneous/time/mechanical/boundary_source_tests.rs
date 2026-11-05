use super::super::boundary_source_tests::{quadratic_profile, quadratic_work};
use super::*;
use crate::acoustics::helmholtz::heterogeneous::time::tests::{MEDIUM, brick, plane, uniform};
use crate::acoustics::helmholtz::{AcousticFeError, HelmholtzFe};

#[test]
fn midside_excitation_drives_complete_wall_memory_and_observable_acoustic_fields() {
    let lengths = [1.0, 0.2, 0.3];
    let (complex, positions) = brick(2, lengths);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let inlet = plane(&geometry, 0.0);
    let release = plane(&geometry, 1.0);
    let profiles = quadratic_profile(&geometry, &inlet, |x| {
        let y = x[1] / lengths[1];
        4.0 * y * (1.0 - y)
    });
    let source = model
        .quadratic_velocity_source(&profiles, &mut || false)
        .unwrap();
    let volume = model
        .volume_velocity_source(4, &mut |_, _| Ok::<_, HeterogeneousError>(1.0), &mut || {
            false
        })
        .unwrap();
    let owner = model
        .prepare_time_with_walls_and_pressure_release(
            0.005,
            &[],
            &[MechanicalBoundary {
                faces: inlet,
                wall: WallPin {
                    surface_density: 0.6,
                    stiffness_per_area: 3.0,
                    resistance: 0.2,
                },
                relaxation: vec![AcousticWallRelaxation {
                    stiffness_per_area: 2.5,
                    relaxation_time_s: 0.04,
                }],
            }],
            &release,
            MechanicalTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let probe = model
        .probe(0, [0.25; 4], &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM))
        .unwrap();
    let mut state = owner.zero_state();
    let mut input = 0.0;
    let mut heat = 0.0;
    for _ in 0..8 {
        let next = owner
            .step_with_boundary_sources(
                &state,
                &[],
                &[(&source, 0.02)],
                &[(&volume, -0.001)],
                &mut || false,
            )
            .unwrap();
        let expected_boundary = quadratic_work(
            &geometry,
            &profiles,
            state.pressure_pa(),
            next.state.pressure_pa(),
            0.02,
            0.005,
        );
        let scale = next.work.input_work_j.abs()
            + next.work.boundary_source_work_j.abs()
            + next.work.volume_source_work_j.abs();
        assert!((next.work.boundary_source_work_j - expected_boundary).abs() < 3e-10 * scale);
        assert!(
            (next.work.input_work_j
                - next.work.boundary_source_work_j
                - next.work.volume_source_work_j)
                .abs()
                < 3e-10 * scale
        );
        input += next.work.input_work_j;
        heat += next.work.wall_dissipation_j;
        assert!(
            (owner.energy(&next.state).unwrap().total_j() + heat - input).abs() < 3e-10 * input
        );
        assert!(next.work.wall_relaxation_dissipation_j >= 0.0);
        let observed = probe.transient(next.state.acoustic_state()).unwrap();
        assert!(observed.pressure_pa.is_finite());
        assert!(observed.particle_velocity_m_s.iter().all(|v| v.is_finite()));
        for (dof, &fixed) in owner
            .condensed
            .pressure_release
            .constrained
            .iter()
            .enumerate()
        {
            if fixed {
                assert_eq!(next.state.pressure_pa()[dof], 0.0);
                assert_eq!(next.state.pressure_potential_pa_s()[dof], 0.0);
            }
        }
        state = next.state;
    }
    let energy = owner.energy(&state).unwrap();
    assert!(input > 0.0 && heat > 0.0);
    assert!(
        energy.wall_kinetic_j > 0.0
            && energy.wall_elastic_j > 0.0
            && energy.wall_relaxation_j > 0.0
    );
    assert!(
        state
            .wall_relaxation_memory(0, 0)
            .unwrap()
            .iter()
            .any(|&z| z != 0.0)
    );
    let observed = probe.transient(state.acoustic_state()).unwrap();
    assert!(
        observed.pressure_pa != 0.0 || observed.particle_velocity_m_s.iter().any(|&v| v != 0.0)
    );

    // Even a zero-amplitude source cannot own a pressure-release face.
    let invalid_source = model
        .quadratic_velocity_source(
            &quadratic_profile(&geometry, &release, |_| 0.0),
            &mut || false,
        )
        .unwrap();
    let accepted_pressure = state.pressure_pa().to_vec();
    let accepted_memory = state.wall_relaxation_memory(0, 0).unwrap().to_vec();
    assert!(
        owner
            .step_with_boundary_sources(&state, &[], &[(&invalid_source, 0.0)], &[], &mut || false)
            .is_err()
    );
    assert!(matches!(
        owner.step_with_boundary_sources(&state, &[], &[(&source, 0.02)], &[], &mut || true),
        Err(HeterogeneousError::Cancelled)
    ));
    assert_eq!(state.pressure_pa(), accepted_pressure);
    assert_eq!(state.wall_relaxation_memory(0, 0).unwrap(), accepted_memory);
    let old = owner
        .step_with_volume_sources(&state, &[], &[(&volume, 0.001)], &mut || false)
        .unwrap();
    let explicit_empty = owner
        .step_with_boundary_sources(&state, &[], &[], &[(&volume, 0.001)], &mut || false)
        .unwrap();
    for (&a, &b) in old
        .state
        .pressure_pa()
        .iter()
        .zip(explicit_empty.state.pressure_pa())
    {
        assert_eq!(a.to_bits(), b.to_bits());
    }
    assert_eq!(
        old.work.input_work_j.to_bits(),
        explicit_empty.work.input_work_j.to_bits()
    );
    assert_eq!(explicit_empty.work.boundary_source_work_j, 0.0);
}
