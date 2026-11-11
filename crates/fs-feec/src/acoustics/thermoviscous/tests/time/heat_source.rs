use super::*;
use crate::acoustics::helmholtz::boundary_layer::BoundaryLayerGeometry;
use crate::acoustics::helmholtz::boundary_layer::thermal_wall::{
    ConductiveThermalWall, ThermalBacking, ThermalWallDepth, WallThermalTransport,
};
use crate::acoustics::thermoviscous::heat_source::ThermoviscousHeatSource;
use crate::acoustics::thermoviscous::time::{
    PreparedThermoviscousHeatTime, TangentialWallConductivity, ThermoviscousCaloricReservoir,
    ThermoviscousWallDiscretization,
};
use crate::acoustics::thermoviscous::wall::ConductiveWallRegion;

fn sealed(mesh: &TetComplex, positions: &[[f64; 3]], size: [f64; 3]) -> ThermoviscousBoundary {
    let mut boundary = boundaries(mesh, positions, size, false);
    for port in boundary.ports.drain(..) {
        boundary
            .no_slip_faces
            .extend(port.traces.into_iter().map(|(face, _)| face));
    }
    boundary
}

fn prepare<'a, 'c>(
    model: &'a ThermoviscousFe<'c>,
    duration_s: f64,
    sources: &[&ThermoviscousHeatSource<'_, '_>],
) -> PreparedThermoviscousHeatTime<'a, 'c> {
    model
        .prepare_time_with_heat_sources(
            duration_s,
            model.medium().cp_j_kg_k,
            &[],
            None,
            None,
            sources,
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        relative(actual, expected) < 2e-8,
        "{actual:e} != {expected:e}"
    );
}

fn same_values(a: impl Iterator<Item = f64>, b: impl Iterator<Item = f64>) {
    let (a, b): (Vec<_>, Vec<_>) = (a.collect(), b.collect());
    assert_eq!(a.len(), b.len());
    let scale = a
        .iter()
        .chain(&b)
        .map(|x| x.abs())
        .fold(1e-12_f64, f64::max);
    for (&actual, &expected) in a.iter().zip(&b) {
        assert!(
            (actual - expected).abs() < 2e-8 * scale,
            "{actual:e} != {expected:e}"
        );
    }
}

fn same_state(a: &ThermoviscousTimeState, b: &ThermoviscousTimeState) {
    same_values(a.pressure_pa.iter().copied(), b.pressure_pa.iter().copied());
    same_values(
        a.velocity_m_s.iter().flatten().copied(),
        b.velocity_m_s.iter().flatten().copied(),
    );
    same_values(
        a.temperature_k.iter().copied(),
        b.temperature_k.iter().copied(),
    );
    same_values(
        a.molecular_temperature_k.iter().flatten().copied(),
        b.molecular_temperature_k.iter().flatten().copied(),
    );
    same_values(
        a.conductive_wall_temperature_k.iter().flatten().copied(),
        b.conductive_wall_temperature_k.iter().flatten().copied(),
    );
    same_values(
        a.wall_displacement_m_sqrt_kg.iter().copied(),
        b.wall_displacement_m_sqrt_kg.iter().copied(),
    );
    same_values(
        a.wall_velocity_m_sqrt_kg_per_s.iter().copied(),
        b.wall_velocity_m_sqrt_kg_per_s.iter().copied(),
    );
}

#[test]
fn g1_held_sealed_heat_release_has_absolute_pressure_temperature_and_signed_work() {
    let size = [0.004, 0.002, 0.002];
    let volume = size.iter().product::<f64>();
    let (mesh, positions) = box_mesh([1, 1, 1], size);
    let gas = air();
    let model = ThermoviscousFe::new(
        &mesh,
        &positions,
        gas,
        &sealed(&mesh, &positions, size),
        limits(),
        &mut || false,
    )
    .unwrap();
    let integral = -2.5;
    let source = model
        .volume_heat_source(4, &mut |_, _| Ok(integral / volume), &mut || false)
        .unwrap();
    let dt = 1e-4;
    let drive = 0.001;
    let prepared = prepare(&model, dt, &[&source]);
    let legacy = model
        .prepare_time(dt, ThermoviscousTimeLimits::default(), &mut || false)
        .unwrap();
    assert_eq!(legacy.system().port_dim(), 0);
    assert_eq!(prepared.system().port_dim(), 1);
    assert_eq!(prepared.source_count(), 1);
    let zero = prepared.zero_state();
    let original = zero.clone();
    let step = prepared
        .step(&zero, &[], &[], &[drive], &mut || false)
        .unwrap();
    assert_eq!(zero, original);
    let cv = gas.cp_j_kg_k - gas.gas_constant_j_kg_k();
    let temperature = integral * drive * dt / (gas.density_kg_m3 * volume * cv);
    let pressure = gas.pressure_pa / gas.temperature_k * temperature;
    for &actual in &step.field.state.temperature_k[..mesh.vertex_count] {
        close(actual, temperature);
    }
    for &actual in &step.field.state.temperature_k[mesh.vertex_count..] {
        assert!(actual.abs() < 2e-8 * temperature.abs());
    }
    for &actual in &step.field.state.pressure_pa {
        close(actual, pressure);
    }
    assert!(
        step.field
            .state
            .velocity_m_s
            .iter()
            .flatten()
            .all(|v| v.abs() < 1e-10)
    );
    let expected_work =
        0.5 * gas.density_kg_m3 * volume * cv / gas.temperature_k * temperature * temperature;
    close(step.field.work.heat_source_work_j, expected_work);
    close(step.heat_source_work_by_source_j[0], expected_work);
    close(step.field.work.final_energy_j, expected_work);
    close(
        step.mean_heat_source_temperature_ratio[0],
        integral * temperature / (2.0 * gas.temperature_k),
    );
    assert!(step.field.work.dissipated_j().abs() < 2e-8 * expected_work);
    assert!(step.field.work.balance_residual_j().abs() < 2e-8 * expected_work);
    // Reverse the held heat to recover the equilibrium. The source absorbs
    // acoustic work during this interval; this signed work must not be clipped.
    let returned = prepared
        .step(&step.field.state, &[], &[], &[-drive], &mut || false)
        .unwrap();
    close(returned.field.work.heat_source_work_j, -expected_work);
    assert!(
        returned
            .field
            .state
            .temperature_k
            .iter()
            .all(|t| t.abs() < 2e-8 * temperature.abs())
    );
    assert!(
        returned
            .field
            .state
            .pressure_pa
            .iter()
            .all(|p| p.abs() < 2e-8 * pressure.abs())
    );
    assert!(returned.field.work.final_energy_j.abs() < 2e-8 * expected_work);
}

#[test]
fn g1_held_heat_resolves_original_molecular_lag_and_same_trajectory_work() {
    let size = [0.004, 0.002, 0.002];
    let volume = size.iter().product::<f64>();
    let (mesh, positions) = box_mesh([1, 1, 1], size);
    let gas = air();
    let model = ThermoviscousFe::new(
        &mesh,
        &positions,
        gas,
        &sealed(&mesh, &positions, size),
        limits(),
        &mut || false,
    )
    .unwrap();
    let source = model
        .volume_heat_source(4, &mut |_, _| Ok(1.0 / volume), &mut || false)
        .unwrap();
    let reservoir = ThermoviscousCaloricReservoir {
        heat_capacity_j_kg_k: 100.0,
        relaxation_time_s: 2e-4,
    };
    let dt = 2e-4;
    let drive = 0.002;
    let frozen_cp = gas.cp_j_kg_k - reservoir.heat_capacity_j_kg_k;
    let prepared = model
        .prepare_time_with_heat_sources(
            dt,
            frozen_cp,
            &[reservoir],
            None,
            None,
            &[&source],
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let step = prepared
        .step(&prepared.zero_state(), &[], &[], &[drive], &mut || false)
        .unwrap();
    let mass = gas.density_kg_m3 * volume;
    let cv = frozen_cp - gas.gas_constant_j_kg_k();
    let capacity = reservoir.heat_capacity_j_kg_k;
    let equilibrium_cv = cv + capacity;
    let rate = (1.0 + capacity / cv) / reservoir.relaxation_time_s;
    let decay = det::exp(-rate * dt);
    let slope = drive / (mass * cv);
    let difference = slope / rate * (1.0 - decay);
    let equilibrium_temperature = drive * dt / (mass * equilibrium_cv);
    let temperature = equilibrium_temperature + capacity / equilibrium_cv * difference;
    let molecular_temperature = equilibrium_temperature - cv / equilibrium_cv * difference;
    for &actual in &step.field.state.temperature_k[..mesh.vertex_count] {
        close(actual, temperature);
    }
    for &actual in &step.field.state.molecular_temperature_k[0][..mesh.vertex_count] {
        close(actual, molecular_temperature);
    }
    for &actual in &step.field.state.pressure_pa {
        close(actual, gas.pressure_pa / gas.temperature_k * temperature);
    }
    let temperature_integral = drive * dt * dt / (2.0 * mass * equilibrium_cv)
        + capacity / equilibrium_cv * slope / rate * (dt - (1.0 - decay) / rate);
    let expected_work = drive / gas.temperature_k * temperature_integral;
    let expected_heat = mass * capacity / (reservoir.relaxation_time_s * gas.temperature_k)
        * (slope / rate).powi(2)
        * (dt - 2.0 * (1.0 - decay) / rate + (1.0 - decay * decay) / (2.0 * rate));
    let expected_energy = 0.5 * mass / gas.temperature_k
        * (cv * temperature * temperature
            + capacity * molecular_temperature * molecular_temperature);
    close(step.field.work.heat_source_work_j, expected_work);
    close(step.field.work.molecular_relaxation_heat_j, expected_heat);
    close(step.field.work.final_energy_j, expected_energy);
    close(
        step.mean_heat_source_temperature_ratio[0],
        temperature_integral / (dt * gas.temperature_k),
    );
    assert!(step.field.work.balance_residual_j().abs() < 2e-8 * expected_work);
}

#[test]
fn g3_heat_and_wall_force_keep_every_loss_channel_and_obey_clock_subdivision() {
    let size = [0.004, 0.002, 0.002];
    let volume = size.iter().product::<f64>();
    let (mesh, positions) = box_mesh([1, 1, 1], size);
    let gas = air();
    let boundary = sealed(&mesh, &positions, size);
    let amplitude = 1.0 / det::sqrt(gas.density_kg_m3 * volume);
    let mut shape = vec![[0.0; 3]; SimplexSpace::new(&mesh, 2).ndof];
    for (coefficient, position) in shape.iter_mut().take(mesh.vertex_count).zip(&positions) {
        for axis in 0..3 {
            coefficient[axis] = amplitude * (position[axis] / size[axis] - 0.5);
        }
    }
    let structure = ModalAcousticTimeModel::try_new_with_free_coordinates(
        48_000,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: 600.0 * std::f64::consts::TAU,
            damping_ratio: 0.1,
            pressure_per_modal_velocity: C64::ZERO,
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    let model = ThermoviscousFe::new_with_moving_wall(
        &mesh,
        &positions,
        gas,
        ThermoviscousThermodynamics::ideal_gas(gas).unwrap(),
        &boundary,
        MovingWall {
            faces: boundary.no_slip_faces.clone(),
            velocity_shapes: vec![shape],
            structure,
        },
        limits(),
        &mut || false,
    )
    .unwrap()
    .with_conductive_walls(
        &[ConductiveWallRegion {
            faces: boundary.no_slip_faces[..1].to_vec(),
            wall: ConductiveThermalWall {
                transport: WallThermalTransport {
                    conductivity_w_m_k: 0.003,
                    volumetric_heat_capacity_j_m3_k: 2000.0,
                },
                depth: ThermalWallDepth::Slab {
                    thickness_m: 0.00002,
                    backing: ThermalBacking::Isothermal,
                },
                minimum_tangential_scale_m: 0.004,
            },
            geometry: BoundaryLayerGeometry {
                minimum_curvature_radius_m: f64::INFINITY,
                minimum_half_gap_m: 0.0001,
                maximum_layer_ratio: 0.05,
            },
            temperature_k: gas.temperature_k,
        }],
        &mut || false,
    )
    .unwrap();
    let uniform = model
        .volume_heat_source(4, &mut |_, _| Ok(1.0 / volume), &mut || false)
        .unwrap();
    let differential = model
        .volume_heat_source(
            4,
            &mut |_, x| Ok((2.0 * x[0] / size[0] - 1.0) / volume),
            &mut || false,
        )
        .unwrap();
    assert!(differential.integrated_shape().abs() < 2e-14);
    let dt = 2e-5;
    let reservoirs = [ThermoviscousCaloricReservoir {
        heat_capacity_j_kg_k: 50.0,
        relaxation_time_s: 2e-5,
    }];
    let wall_elements = [2];
    let tangential = [TangentialWallConductivity {
        isotropic_w_m_k: 0.003,
    }];
    let build = |duration| {
        model
            .prepare_time_with_heat_sources(
                duration,
                gas.cp_j_kg_k - 50.0,
                &reservoirs,
                Some(ThermoviscousWallDiscretization {
                    normal_elements: &wall_elements,
                    maximum_frequency_hz: 2000.0,
                }),
                Some(&tangential),
                &[&uniform, &differential],
                ThermoviscousTimeLimits::default(),
                &mut || false,
            )
            .unwrap()
    };
    let full = build(dt);
    let half = build(0.5 * dt);
    let initial = full.zero_state();
    let whole = full
        .step(&initial, &[], &[1.0], &[0.003, -0.002], &mut || false)
        .unwrap();
    let first = half
        .step(&initial, &[], &[1.0], &[0.003, -0.002], &mut || false)
        .unwrap();
    let second = half
        .step(
            &first.field.state,
            &[],
            &[1.0],
            &[0.003, -0.002],
            &mut || false,
        )
        .unwrap();
    same_state(&whole.field.state, &second.field.state);
    assert_eq!(whole.mean_heat_source_temperature_ratio.len(), 2);
    assert_eq!(whole.field.mean_wall_velocity_m_sqrt_kg_per_s.len(), 1);
    let channels = |work: &crate::acoustics::thermoviscous::time::ThermoviscousTimeWork| {
        [
            work.fluid_viscous_heat_j,
            work.fluid_conductive_heat_j,
            work.structural_damping_heat_j,
            work.molecular_relaxation_heat_j,
            work.solid_conductive_heat_j,
        ]
    };
    for ((actual, a), b) in channels(&whole.field.work)
        .into_iter()
        .zip(channels(&first.field.work))
        .zip(channels(&second.field.work))
    {
        assert!(actual > 0.0, "{:?}", whole.field.work);
        close(actual, a + b);
    }
    for (source, &actual) in whole.heat_source_work_by_source_j.iter().enumerate() {
        close(
            actual,
            first.heat_source_work_by_source_j[source]
                + second.heat_source_work_by_source_j[source],
        );
    }
    close(
        whole.field.work.supplied_j(),
        first.field.work.supplied_j() + second.field.work.supplied_j(),
    );
    assert!(
        whole.field.work.balance_residual_j().abs() < 2e-8 * whole.field.work.supplied_j().abs()
    );
}

#[test]
fn g0_heat_preparation_and_steps_refuse_foreign_nonfinite_and_cancelled_data() {
    let size = [0.004, 0.002, 0.002];
    let (mesh, positions) = box_mesh([1, 1, 1], size);
    let boundary = sealed(&mesh, &positions, size);
    let make = || {
        ThermoviscousFe::new(&mesh, &positions, air(), &boundary, limits(), &mut || false).unwrap()
    };
    let model = make();
    let other = make();
    let source = model
        .volume_heat_source(4, &mut |_, _| Ok(1.0), &mut || false)
        .unwrap();
    let foreign = other
        .volume_heat_source(4, &mut |_, _| Ok(1.0), &mut || false)
        .unwrap();
    for sources in [&[][..], &[&foreign][..]] {
        assert!(
            model
                .prepare_time_with_heat_sources(
                    1e-4,
                    air().cp_j_kg_k,
                    &[],
                    None,
                    None,
                    sources,
                    ThermoviscousTimeLimits::default(),
                    &mut || false,
                )
                .is_err()
        );
    }
    let too_many = vec![&source; fs_phs::linear_exact::MAX_LINEAR_EXPONENTIAL_BLOCK_DIM];
    assert!(
        model
            .prepare_time_with_heat_sources(
                1e-4,
                air().cp_j_kg_k,
                &[],
                None,
                None,
                &too_many,
                ThermoviscousTimeLimits::default(),
                &mut || false,
            )
            .is_err()
    );
    assert!(matches!(
        model.prepare_time_with_heat_sources(
            1e-4,
            air().cp_j_kg_k,
            &[],
            None,
            None,
            &[&source],
            ThermoviscousTimeLimits::default(),
            &mut || true,
        ),
        Err(ThermoviscousTimeError::Fluid(ThermoviscousError::Cancelled))
    ));
    let prepared = prepare(&model, 1e-4, &[&source]);
    let legacy = model
        .prepare_time(1e-4, ThermoviscousTimeLimits::default(), &mut || false)
        .unwrap();
    let mut initial = prepared.zero_state();
    initial.temperature_k[mesh.vertex_count] = 0.25;
    let original = initial.clone();
    let expected = legacy.step(&initial, &[], &[], &mut || false).unwrap();
    let actual = prepared
        .step(&initial, &[], &[], &[0.0], &mut || false)
        .unwrap();
    same_state(&actual.field.state, &expected.state);
    assert_eq!(actual.field.work.heat_source_work_j, 0.0);
    for drive in [&[][..], &[f64::NAN][..], &[f64::INFINITY][..]] {
        assert!(
            prepared
                .step(&initial, &[], &[], drive, &mut || false)
                .is_err()
        );
        assert_eq!(initial, original);
    }
    let mut polls = 0;
    assert!(matches!(
        prepared.step(&initial, &[], &[], &[1.0], &mut || {
            polls += 1;
            polls == 2
        }),
        Err(ThermoviscousTimeError::Fluid(ThermoviscousError::Cancelled))
    ));
    assert_eq!(initial, original);
}

// The wall shape is antisymmetric about the centre of the cell and the one free velocity dof (the interior diagonal edge) is
// symmetric about it, so their exact cross entries are 0. `congruence` sums each pair in a different order per triangle; on
// this fixture one triangle held 1.98e-17 and the other exactly 0, and the mass inverse check refused the matrix.
#[test]
fn the_moving_wall_mass_is_exactly_symmetric_on_the_breathing_wall_fixture() {
    let size = [0.004, 0.002, 0.002];
    let volume = size.iter().product::<f64>();
    let (mesh, positions) = box_mesh([1, 1, 1], size);
    let gas = air();
    let boundary = sealed(&mesh, &positions, size);
    let amplitude = 1.0 / det::sqrt(gas.density_kg_m3 * volume);
    let mut shape = vec![[0.0; 3]; SimplexSpace::new(&mesh, 2).ndof];
    for (coefficient, position) in shape.iter_mut().take(mesh.vertex_count).zip(&positions) {
        for axis in 0..3 {
            coefficient[axis] = amplitude * (position[axis] / size[axis] - 0.5);
        }
    }
    let structure = ModalAcousticTimeModel::try_new_with_free_coordinates(
        48_000,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: 600.0 * std::f64::consts::TAU,
            damping_ratio: 0.1,
            pressure_per_modal_velocity: C64::ZERO,
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    let model = ThermoviscousFe::new_with_moving_wall(
        &mesh,
        &positions,
        gas,
        ThermoviscousThermodynamics::ideal_gas(gas).unwrap(),
        &boundary,
        MovingWall {
            faces: boundary.no_slip_faces.clone(),
            velocity_shapes: vec![shape],
            structure,
        },
        limits(),
        &mut || false,
    )
    .unwrap();
    let n = model.n;
    assert_eq!(
        (
            model.velocity_count,
            model.pressure_start - 3 * model.velocity_count
        ),
        (1, 1)
    );
    for i in 0..n {
        for j in 0..i {
            assert!(
                model.mass[i * n + j] == model.mass[j * n + i],
                "mass pair ({i}, {j}): {:e} against {:e}",
                model.mass[i * n + j],
                model.mass[j * n + i]
            );
        }
    }
}
