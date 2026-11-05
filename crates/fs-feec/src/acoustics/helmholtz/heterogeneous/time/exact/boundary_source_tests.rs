use super::super::AcousticFeError;
use super::super::tests::MEDIUM;
use super::*;
use crate::acoustics::helmholtz::HelmholtzFe;
use crate::acoustics::helmholtz::heterogeneous::MAX_BOUNDARY_SOURCES;
use crate::fixtures::single_tet;

fn constant<'f, 'c>(geometry: &'f HelmholtzFe<'c>) -> HeterogeneousHelmholtzFe<'f, 'c> {
    HeterogeneousHelmholtzFe::new(
        geometry,
        geometry.order() + 2,
        &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
        &mut || false,
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn equivalent(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    let maximum = expected.iter().map(|value| value.abs()).fold(0.0, f64::max);
    for (&a, &b) in actual.iter().zip(expected) {
        close(a, b, 3e-9 * maximum);
    }
}

#[test]
fn exact_midside_motion_has_analytic_pressure_signed_volume_work_and_robin_heat() {
    let (complex, positions) = single_tet();
    let geometry = HelmholtzFe::new(&complex, &positions, 1, MEDIUM).unwrap();
    let model = constant(&geometry);
    // Every corner is stationary; only edge01's midpoint moves. On z=0 the
    // prescribed shape is g=4 lambda0 lambda1, so a corner-only port erases it.
    let boundary_source = model
        .quadratic_velocity_source(&[([0, 1, 2], [0.0, 0.0, 0.0, 1.0, 0.0, 0.0])], &mut || {
            false
        })
        .unwrap();
    let volume_source = model
        .volume_velocity_source(3, &mut |_, _| Ok::<_, HeterogeneousError>(6.0), &mut || {
            false
        })
        .unwrap();
    let rho = MEDIUM.density_kg_m3;
    let speed = MEDIUM.sound_speed_m_s;
    let duration = 0.6 / speed;
    let admittance = 0.2 / (rho * speed);
    let owner = model
        .prepare_exact_time_with_boundary_sources(
            duration,
            &[ResistiveBoundary {
                faces: vec![[0, 1, 2]],
                admittance_m_pa_s: admittance,
            }],
            &[],
            &[&boundary_source],
            &[&volume_source],
            &[[1, 2, 3]],
            ExactTimeLimits {
                max_dofs: 1,
                ..ExactTimeLimits::default()
            },
            &mut || false,
        )
        .unwrap();
    assert_eq!(owner.port_count(), 1);
    assert_eq!(owner.volume_source_count(), 1);
    // G1: original P1 tetrahedron, with only phi0=1-x-y-z free. The original
    // area A=1/2 appears once: int(g)=A/3=1/6 and int(phi0*g)=2A/15=1/15.
    // The normalized volume source contributes int(phi0/V)=1/4, V=1/6.
    let boundary_load = 1.0 / 15.0;
    let volume_load = 1.0 / 4.0;
    close(boundary_source.integrated_shape_m2(), 1.0 / 6.0, 3e-14);
    close(boundary_source.load[0], boundary_load, 3e-14);
    let mass = 1.0 / (60.0 * rho * speed * speed);
    let stiffness = 1.0 / (2.0 * rho);
    let resistance = admittance / 12.0;
    let alpha = resistance / (2.0 * mass);
    let omega = det::sqrt(stiffness / mass - alpha * alpha);
    let phase = omega * duration;
    let decay = det::exp(-alpha * duration);
    let squared_velocity_integral = (1.0 - decay * decay) / (4.0 * alpha)
        - (alpha
            + decay * decay * (-alpha * det::cos(2.0 * phase) + omega * det::sin(2.0 * phase)))
            / (4.0 * (alpha * alpha + omega * omega));
    for scale in [1.0, 1e-8] {
        let velocity = 0.03 * scale;
        let volume_flow = -0.004 * scale;
        let force = boundary_load * velocity + volume_load * volume_flow;
        let potential =
            force / stiffness * (1.0 - decay * (det::cos(phase) + alpha / omega * det::sin(phase)));
        let pressure_scale = force / (mass * omega);
        let pressure = pressure_scale * decay * det::sin(phase);
        let heat = resistance * pressure_scale * pressure_scale * squared_velocity_integral;
        let energy = 0.5 * (mass * pressure * pressure + stiffness * potential * potential);
        let result = owner
            .step_with_volume_sources(
                &owner.zero_state(),
                &[velocity],
                &[volume_flow],
                &mut || false,
            )
            .unwrap();
        close(
            result.state.pressure_pa()[0],
            pressure,
            3e-9 * pressure_scale.abs(),
        );
        close(
            result.state.pressure_potential_pa_s()[0],
            potential,
            3e-9 * potential.abs(),
        );
        for dof in 1..4 {
            assert_eq!(result.state.pressure_pa()[dof], 0.0);
            assert_eq!(result.state.pressure_potential_pa_s()[dof], 0.0);
        }
        let boundary_integral = boundary_load * potential;
        let volume_integral = volume_load * potential;
        let extraction = volume_flow * volume_integral;
        let input = force * potential;
        close(
            result.integrated_port_pressure_pa_m2_s[0],
            boundary_integral,
            3e-9 * boundary_integral,
        );
        close(
            result.integrated_volume_source_pressure_pa_s[0],
            volume_integral,
            3e-9 * volume_integral,
        );
        close(
            result.work.volume_source_work_j,
            extraction,
            3e-9 * extraction.abs(),
        );
        close(result.work.input_work_j, input, 3e-9 * input);
        close(
            result.work.input_work_j,
            velocity * boundary_integral + extraction,
            3e-9 * input,
        );
        close(result.work.boundary_dissipation_j, heat, 3e-9 * heat);
        close(result.work.final_energy_j, energy, 3e-9 * energy);
        close(result.work.balance_residual_j(), 0.0, 3e-9 * input);
        assert!(result.work.volume_source_work_j < 0.0);
        assert!(result.work.boundary_dissipation_j > 0.0);
    }
}

#[test]
fn exact_prepared_linear_shapes_preserve_port_order_and_refuse_foreign_release_or_overbudget_columns()
 {
    let (complex, positions) = single_tet();
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = constant(&geometry);
    let ports = [
        ExactVelocityPort {
            profile: vec![([0, 1, 2], [1.0, -0.5, 0.25])],
        },
        ExactVelocityPort {
            profile: vec![([0, 2, 3], [0.0, 1.0, -1.0])],
        },
    ];
    let shape = [([0, 2, 3], [0.0, 1.0, -1.0, 0.5, 0.0, -0.5])];
    let prepared = model
        .quadratic_velocity_source(&shape, &mut || false)
        .unwrap();
    // Linear data at all six nodes must reproduce the old linear trace.
    // Its signed area integral is zero, but its original P2 load is nonzero.
    close(prepared.integrated_shape_m2(), 0.0, 3e-14);
    assert!(prepared.load.iter().any(|entry| entry.abs() > 1e-4));
    let boundaries = [ResistiveBoundary {
        faces: vec![[0, 1, 2]],
        admittance_m_pa_s: 0.1,
    }];
    let release = [[1, 2, 3]];
    let mixed = model
        .prepare_exact_time_with_boundary_sources(
            0.03,
            &boundaries,
            &ports[..1],
            &[&prepared],
            &[],
            &release,
            ExactTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let linear = model
        .prepare_exact_time_with_pressure_release(
            0.03,
            &boundaries,
            &ports,
            &[],
            &release,
            ExactTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(mixed.port_count(), 2);
    assert_eq!(mixed.volume_source_count(), 0);
    let initial = mixed.zero_state();
    assert!(mixed.step(&initial, &[0.02], &mut || false).is_err());
    assert!(
        mixed
            .step(&initial, &[0.02, -0.04, 0.0], &mut || false)
            .is_err()
    );
    let actual = mixed.step(&initial, &[0.02, -0.04], &mut || false).unwrap();
    let expected = linear
        .step(&initial, &[0.02, -0.04], &mut || false)
        .unwrap();
    equivalent(actual.state.pressure_pa(), expected.state.pressure_pa());
    equivalent(
        actual.state.pressure_potential_pa_s(),
        expected.state.pressure_potential_pa_s(),
    );
    equivalent(
        &actual.integrated_port_pressure_pa_m2_s,
        &expected.integrated_port_pressure_pa_m2_s,
    );
    let work = |w: ExactTimeWork| {
        [
            w.initial_energy_j,
            w.final_energy_j,
            w.delta_energy_j,
            w.input_work_j,
            w.volume_source_work_j,
            w.boundary_dissipation_j,
        ]
    };
    equivalent(&work(actual.work), &work(expected.work));

    let too_many = vec![&prepared; MAX_BOUNDARY_SOURCES + 1];
    assert!(matches!(
        model.prepare_exact_time_with_boundary_sources(
            0.03,
            &[],
            &[],
            &too_many,
            &[],
            &release,
            ExactTimeLimits::default(),
            &mut || false,
        ),
        Err(HeterogeneousError::Limit {
            what: "simultaneous prepared boundary acoustic sources",
            ..
        })
    ));
    let augmented_count =
        MAX_LINEAR_EXPONENTIAL_BLOCK_DIM / 2 - 2 * mixed.pressure_release.free_dofs.len();
    assert!(augmented_count <= MAX_BOUNDARY_SOURCES);
    let augmented_overflow = vec![&prepared; augmented_count];
    assert!(matches!(
        model.prepare_exact_time_with_boundary_sources(
            0.03,
            &[],
            &[],
            &augmented_overflow,
            &[],
            &release,
            ExactTimeLimits::default(),
            &mut || false,
        ),
        Err(HeterogeneousError::Limit {
            what: "exact transient augmented states",
            ..
        })
    ));
    let another_model = constant(&geometry);
    let foreign = another_model
        .quadratic_velocity_source(&shape, &mut || false)
        .unwrap();
    assert!(
        model
            .prepare_exact_time_with_boundary_sources(
                0.03,
                &[],
                &[],
                &[&foreign],
                &[],
                &release,
                ExactTimeLimits::default(),
                &mut || false,
            )
            .is_err()
    );
    let zero_on_release = model
        .quadratic_velocity_source(&[([3, 2, 1], [0.0; 6])], &mut || false)
        .unwrap();
    assert_eq!(zero_on_release.integrated_shape_m2(), 0.0);
    assert!(
        model
            .prepare_exact_time_with_boundary_sources(
                0.03,
                &[],
                &[],
                &[&zero_on_release],
                &[],
                &release,
                ExactTimeLimits::default(),
                &mut || false,
            )
            .is_err()
    );
}
