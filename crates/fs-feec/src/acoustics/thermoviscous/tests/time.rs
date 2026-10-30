use super::*;
use crate::acoustics::thermoviscous::moving_wall::MovingWall;
use crate::acoustics::thermoviscous::time::{
    ThermoviscousTimeError, ThermoviscousTimeLimits, ThermoviscousTimeState,
};
use fs_modal::time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};

fn translating<'c>(
    mesh: &'c TetComplex,
    pos: &[[f64; 3]],
    density: f64,
    amplitude: f64,
    natural: f64,
    damping: f64,
) -> ThermoviscousFe<'c> {
    // A manufactured zero-expansion liquid; no measured-water claim.
    let medium = ThermoviscousMedium {
        density_kg_m3: density,
        pressure_pa: 101325.0,
        temperature_k: 300.0,
        cp_j_kg_k: 4000.0,
        dynamic_viscosity_pa_s: 0.001,
        thermal_conductivity_w_m_k: 0.6,
        bulk_viscosity_pa_s: 0.002,
    };
    let thermo = ThermoviscousThermodynamics {
        isothermal_compressibility_pa_inv: 1.0 / 2e9,
        isobaric_expansion_per_k: 0.0,
    };
    let geometry = HelmholtzFe::new(
        mesh,
        pos,
        2,
        Medium {
            density_kg_m3: density,
            sound_speed_m_s: thermo.sound_speed_m_s(medium).unwrap(),
        },
    )
    .unwrap();
    let faces = geometry.boundary_faces();
    let mut phi = vec![[0.0; 3]; SimplexSpace::new(mesh, 2).ndof];
    for face in &faces {
        for &v in face {
            phi[v as usize][0] = amplitude;
        }
    }
    let structure = ModalAcousticTimeModel::try_new_with_free_coordinates(
        48_000,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: natural,
            damping_ratio: damping,
            pressure_per_modal_velocity: C64::ZERO,
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    ThermoviscousFe::new_with_moving_wall(
        mesh,
        pos,
        medium,
        thermo,
        &ThermoviscousBoundary {
            no_slip_faces: faces.clone(),
            isothermal_faces: vec![],
            ports: vec![],
        },
        MovingWall {
            faces,
            velocity_shapes: vec![phi],
            structure,
        },
        limits(),
        &mut || false,
    )
    .unwrap()
}

fn relative(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs().max(1e-30)
}

#[test]
fn g1_free_structure_and_original_liquid_mass_accelerate_together() {
    let size = [0.002, 0.001, 0.001];
    let volume = size.iter().product::<f64>();
    let (mesh, pos) = box_mesh([1, 1, 1], size);
    let amplitude = 1.0 / det::sqrt(1000.0 * volume);
    let dt = 1e-4;
    for density in [1000.0, 7000.0] {
        let model = translating(&mesh, &pos, density, amplitude, 0.0, 0.0);
        let prepared = model
            .prepare_time(dt, ThermoviscousTimeLimits::default(), &mut || false)
            .unwrap();
        let mut initial = prepared.zero_state();
        let acceleration = 1.0 / (1.0 + density * volume * amplitude * amplitude);
        // Exact continuum solution: uniform acceleration, zero strain rate,
        // and a fixed linear pressure gradient. P1/P2 represent it exactly.
        for (p, xyz) in initial.pressure_pa.iter_mut().zip(&pos) {
            *p = -density * amplitude * acceleration * (xyz[0] - 0.5 * size[0]);
        }
        let original = initial.clone();
        let frame = prepared.step(&initial, &[], &[1.0], &mut || false).unwrap();
        assert_eq!(initial, original);
        assert!(
            relative(
                frame.state.wall_velocity_m_sqrt_kg_per_s[0],
                acceleration * dt
            ) < 2e-8
        );
        assert!(
            relative(
                frame.state.wall_displacement_m_sqrt_kg[0],
                0.5 * acceleration * dt * dt
            ) < 2e-8
        );
        assert!(
            relative(
                frame.mean_wall_velocity_m_sqrt_kg_per_s[0],
                0.5 * acceleration * dt
            ) < 2e-8
        );
        for velocity in &frame.state.velocity_m_s[..mesh.vertex_count] {
            assert!(relative(velocity[0], amplitude * acceleration * dt) < 2e-8);
            assert!(velocity[1].abs() + velocity[2].abs() < 2e-8 * velocity[0].abs());
        }
        let expected_work = 0.5 * acceleration * dt * dt;
        assert!(relative(frame.work.structural_work_j, expected_work) < 2e-8);
        assert!(relative(frame.work.delta_energy_j, expected_work) < 2e-8);
        assert!(frame.work.dissipated_j().abs() < 2e-8 * expected_work);
    }
}

#[test]
fn g1_original_temperature_pressure_and_viscous_fields_decay_with_closed_work() {
    let size = [0.004, 0.0008, 0.001];
    let (mesh, pos) = box_mesh([1, 1, 1], size);
    let medium = air();
    let boundary = boundaries(&mesh, &pos, size, true);
    let model =
        ThermoviscousFe::new(&mesh, &pos, medium, &boundary, limits(), &mut || false).unwrap();
    let prepared = model
        .prepare_time(1e-4, ThermoviscousTimeLimits::default(), &mut || false)
        .unwrap();
    let mut state = prepared.zero_state();
    let thermal = model
        .temperature_map
        .iter()
        .position(Option::is_some)
        .unwrap();
    state.temperature_k[thermal] = 1.0;
    let initial_energy = prepared.stored_energy_j(&state).unwrap();
    let mut heat = 0.0;
    let mut pressure_seen = 0.0_f64;
    let mut conductive = 0.0;
    let mut viscous = 0.0;
    for _ in 0..128 {
        let frame = prepared
            .step(&state, &[0.0, 0.0], &[], &mut || false)
            .unwrap();
        assert!(frame.work.final_energy_j <= frame.work.initial_energy_j * (1.0 + 1e-10));
        heat += frame.work.dissipated_j();
        conductive += frame.work.fluid_conductive_heat_j;
        viscous += frame.work.fluid_viscous_heat_j;
        pressure_seen = pressure_seen.max(
            frame
                .state
                .pressure_pa
                .iter()
                .map(|p| p.abs())
                .fold(0.0, f64::max),
        );
        state = frame.state;
    }
    assert!(pressure_seen > 0.0 && conductive > 0.0 && viscous > 0.0);
    let final_energy = prepared.stored_energy_j(&state).unwrap();
    assert!(final_energy < initial_energy);
    assert!((final_energy + heat - initial_energy).abs() < 2e-8 * initial_energy);
}

fn phasor_state(field: &ThermoviscousField, angular_frequency: f64) -> ThermoviscousTimeState {
    ThermoviscousTimeState {
        pressure_pa: field.pressure_pa.iter().map(|p| p.re).collect(),
        velocity_m_s: field.velocity_m_s.iter().map(|v| v.map(|z| z.re)).collect(),
        temperature_k: field.temperature_k.iter().map(|v| v.re).collect(),
        wall_displacement_m_sqrt_kg: field
            .moving_wall
            .modal_velocity_m_sqrt_kg_per_s
            .iter()
            .map(|u| u.im / angular_frequency)
            .collect(),
        wall_velocity_m_sqrt_kg_per_s: field
            .moving_wall
            .modal_velocity_m_sqrt_kg_per_s
            .iter()
            .map(|u| u.re)
            .collect(),
    }
}

#[test]
fn g3_same_original_harmonic_load_is_recovered_as_held_force_clock_refines() {
    let size = [0.002, 0.001, 0.001];
    let (mesh, pos) = box_mesh([1, 1, 1], size);
    let amplitude = 1.0 / det::sqrt(1000.0 * size.iter().product::<f64>());
    let model = translating(
        &mesh,
        &pos,
        1000.0,
        amplitude,
        250.0 * std::f64::consts::TAU,
        0.1,
    );
    let frequency = 100.0;
    let omega = std::f64::consts::TAU * frequency;
    let harmonic = model
        .prepare(frequency, &mut || false)
        .unwrap()
        .response_with_wall_forces(&[], &[C64::ONE], &mut || false)
        .unwrap();
    let initial = phasor_state(&harmonic, omega);
    let expected = initial.wall_velocity_m_sqrt_kg_per_s[0];
    let scale = harmonic.moving_wall.modal_velocity_m_sqrt_kg_per_s[0].abs();
    let mut errors = Vec::new();
    for samples in [100, 400] {
        let dt = 1.0 / (frequency * samples as f64);
        let prepared = model
            .prepare_time(dt, ThermoviscousTimeLimits::default(), &mut || false)
            .unwrap();
        let mut state = initial.clone();
        for k in 0..samples {
            let force = (omega * (k as f64 + 0.5) * dt).cos();
            state = prepared
                .step(&state, &[], &[force], &mut || false)
                .unwrap()
                .state;
        }
        errors.push((state.wall_velocity_m_sqrt_kg_per_s[0] - expected).abs() / scale);
    }
    assert!(errors[1] < 0.2 * errors[0], "{errors:?}");
    assert!(errors[1] < 2e-4, "{errors:?}");
}

#[test]
fn g0_cancellation_bad_trace_capacity_and_frequency_walls_refuse_transactionally() {
    use crate::acoustics::helmholtz::boundary_layer::BoundaryLayerGeometry;
    use crate::acoustics::helmholtz::boundary_layer::thermal_wall::{
        ConductiveThermalWall, ThermalBacking, ThermalWallDepth, WallThermalTransport,
    };
    use crate::acoustics::thermoviscous::wall::ConductiveWallRegion;
    let size = [0.004, 0.0008, 0.001];
    let (mesh, pos) = box_mesh([1, 1, 1], size);
    let medium = air();
    let boundary = boundaries(&mesh, &pos, size, false);
    let model =
        ThermoviscousFe::new(&mesh, &pos, medium, &boundary, limits(), &mut || false).unwrap();
    let tiny = ThermoviscousTimeLimits {
        max_states: 1,
        ..ThermoviscousTimeLimits::default()
    };
    assert!(model.prepare_time(1e-4, tiny, &mut || false).is_err());
    {
        let prepared = model
            .prepare_time(1e-4, ThermoviscousTimeLimits::default(), &mut || false)
            .unwrap();
        let state = prepared.zero_state();
        let original = state.clone();
        let mut calls = 0;
        assert!(matches!(
            prepared.step(&state, &[1.0, 0.0], &[], &mut || {
                calls += 1;
                calls == 2
            }),
            Err(ThermoviscousTimeError::Fluid(ThermoviscousError::Cancelled))
        ));
        assert_eq!(state, original);
        let mut invalid = state.clone();
        let fixed = model.velocity_map.iter().position(Option::is_none).unwrap();
        invalid.velocity_m_s[fixed][0] = 1e-10;
        assert!(
            prepared
                .step(&invalid, &[0.0, 0.0], &[], &mut || false)
                .is_err()
        );
    }
    let conductive = model
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
                        backing: ThermalBacking::Adiabatic,
                    },
                    minimum_tangential_scale_m: 0.004,
                },
                geometry: BoundaryLayerGeometry {
                    minimum_curvature_radius_m: f64::INFINITY,
                    minimum_half_gap_m: 0.0001,
                    maximum_layer_ratio: 0.05,
                },
                temperature_k: medium.temperature_k,
            }],
            &mut || false,
        )
        .unwrap();
    assert!(
        conductive
            .prepare_time(1e-4, ThermoviscousTimeLimits::default(), &mut || false)
            .is_err()
    );
}
