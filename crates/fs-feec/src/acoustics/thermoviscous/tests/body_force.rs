//! G1/G3: original weak force, source work, loaded ports and wall reaction.
use super::super::body_force::MAX_BODY_FORCE_SOURCES;
use super::super::moving_wall::{MovingWall, modal_impedance};
use super::*;
use fs_modal::time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use std::cell::Cell;

fn relative(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance * expected.abs().max(1e-30),
        "actual={actual:?}, expected={expected:?}"
    );
}

#[test]
fn g1_original_force_projection_retains_vector_moments_and_zero_resultant_shapes() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let positions = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let boundary = ThermoviscousBoundary {
        no_slip_faces: vec![],
        isothermal_faces: vec![],
        ports: vec![BoundaryMode::uniform(&mesh.faces)],
    };
    let model =
        ThermoviscousFe::new(&mesh, &positions, air(), &boundary, limits(), &mut || false).unwrap();
    let source = model
        .volume_force_source(
            4,
            &mut |_, x| Ok::<_, ThermoviscousError>([1.0 + x[0], 2.0 * x[1], -x[2]]),
            &mut || false,
        )
        .unwrap();
    for (actual, expected) in
        source
            .integrated_shape()
            .into_iter()
            .zip([5.0 / 24.0, 1.0 / 12.0, -1.0 / 24.0])
    {
        assert!((actual - expected).abs() < 3e-15);
    }
    // Independent tetrahedron moments: integral lambda_i=1/24;
    // integral lambda_i lambda_j=(1+delta_ij)/120. Vertex tests are lambda_i.
    for vertex in 0..4 {
        let index = model.velocity_map[vertex].unwrap();
        let expected = [
            1.0 / 24.0 + if vertex == 1 { 1.0 / 60.0 } else { 1.0 / 120.0 },
            if vertex == 2 { 1.0 / 30.0 } else { 1.0 / 60.0 },
            if vertex == 3 {
                -1.0 / 60.0
            } else {
                -1.0 / 120.0
            },
        ];
        for axis in 0..3 {
            assert!(
                (source.load[axis * model.velocity_count + index] - expected[axis]).abs() < 3e-15
            );
        }
    }
    let balanced = model
        .volume_force_source(
            4,
            &mut |_, x| Ok::<_, ThermoviscousError>([x[0] - 0.25, 0.0, 0.0]),
            &mut || false,
        )
        .unwrap();
    assert!(balanced.integrated_shape().iter().all(|v| v.abs() < 3e-15));
    assert!(balanced.load.iter().any(|v| v.abs() > 1e-4));
}

#[test]
fn g1_curved_force_projection_uses_original_positions_and_variable_measure() {
    use crate::acoustics::helmholtz::QuadraticGeometryControls;
    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 41,
                kernel_id: 7,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
        let original = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let a = 0.6;
        let map = |[x, y, z]: [f64; 3]| [x, y, z * (1.0 + a * x)];
        let positions: Vec<_> = original.into_iter().map(map).collect();
        let midsides: Vec<_> = mesh
            .edges
            .iter()
            .map(|&[i, j]| {
                map(std::array::from_fn(|axis| {
                    f64::midpoint(original[i as usize][axis], original[j as usize][axis])
                }))
            })
            .collect();
        let gas = air();
        let thermo = ThermoviscousThermodynamics::ideal_gas(gas).unwrap();
        let geometry = HelmholtzFe::new_quadratic(
            &mesh,
            &positions,
            &midsides,
            2,
            6,
            Medium {
                density_kg_m3: gas.density_kg_m3,
                sound_speed_m_s: thermo.sound_speed_m_s(gas).unwrap(),
            },
            QuadraticGeometryControls::default(),
            &cx,
        )
        .unwrap();
        // A port needs a planar face (`boundary_geometry` refuses a nonaffine one, pinned by `curved.rs`): the two faces
        // this map bends are no-slip walls, the two that stay planar carry the port. The volume-force integrals below
        // do not depend on the boundary data.
        let (planar, bent): (Vec<_>, Vec<_>) = mesh
            .faces
            .iter()
            .copied()
            .partition(|&face| geometry.boundary_geometry(face).is_ok());
        assert_eq!((planar.len(), bent.len()), (2, 2));
        let boundary = ThermoviscousBoundary {
            no_slip_faces: bent,
            isothermal_faces: vec![],
            ports: vec![BoundaryMode::uniform(&planar)],
        };
        let model =
            ThermoviscousFe::from_geometry(geometry, gas, &boundary, limits(), &mut || false)
                .unwrap();
        let source = model
            .volume_force_source(
                6,
                &mut |_, x| Ok::<_, ThermoviscousError>([1.0, x[2], 0.0]),
                &mut || false,
            )
            .unwrap();
        // x=X,y=Y,z=Z(1+aX), determinant=1+aX.
        let expected = [
            1.0 / 6.0 + a / 24.0,
            1.0 / 24.0 + a / 60.0 + a * a / 360.0,
            0.0,
        ];
        for (actual, expected) in source.integrated_shape().into_iter().zip(expected) {
            assert!((actual - expected).abs() < 5e-14);
        }
    });
}

#[test]
fn g3_force_heat_and_traction_share_original_fields_work_and_loaded_scattering() {
    let size = [0.004, 0.001, 0.0012];
    let volume = size.iter().product::<f64>();
    let (mesh, positions) = box_mesh([1, 2, 2], size);
    let model = ThermoviscousFe::new(
        &mesh,
        &positions,
        air(),
        &boundaries(&mesh, &positions, size, true),
        limits(),
        &mut || false,
    )
    .unwrap();
    let calls = Cell::new(0);
    let source = model
        .volume_force_source(
            4,
            &mut |_, x| {
                calls.set(calls.get() + 1);
                Ok::<_, ThermoviscousError>([
                    (1.0 + x[0] / size[0]) / volume,
                    -x[1] / (size[1] * volume),
                    0.0,
                ])
            },
            &mut || false,
        )
        .unwrap();
    let count = calls.get();
    let heat = model
        .volume_heat_source(4, &mut |_, _| Ok(1.0 / volume), &mut || false)
        .unwrap();
    let force_sources = [(&source, C64::new(3e-6, -2e-6))];
    let heat_sources = [(&heat, C64::new(3e-4, -1e-4))];
    let cp = C64::new(970.0, -45.0);
    let prepared = model
        .prepare_with_heat_capacity(1700.0, cp, &mut || false)
        .unwrap();
    let pressure = [C64::new(1.0, 0.2), C64::new(-0.3, 0.1)];
    let combined = prepared
        .response_with_volume_sources(&pressure, &[], &force_sources, &heat_sources, &mut || false)
        .unwrap();
    let force = prepared
        .response_with_volume_sources(&[C64::ZERO; 2], &[], &force_sources, &[], &mut || false)
        .unwrap();
    let thermal = prepared
        .response_with_heat_sources(&pressure, &[], &heat_sources, &mut || false)
        .unwrap();
    for (actual, a, b) in [
        (
            &combined.pressure_pa,
            &force.pressure_pa,
            &thermal.pressure_pa,
        ),
        (
            &combined.temperature_k,
            &force.temperature_k,
            &thermal.temperature_k,
        ),
    ] {
        let scale = a.iter().chain(b).map(|z| z.abs()).fold(1e-30, f64::max);
        for ((x, y), z) in actual.iter().zip(a).zip(b) {
            assert!((*x - *y - *z).abs() < 2e-8 * scale);
        }
    }
    assert!(force.power.body_force_work_w > 0.0);
    assert_work(&combined);
    assert_work(&force);
    let reference: [f64; 2] = [1.3e8, 0.7e8];
    let incident = [C64::new(0.0001, 0.0002), C64::new(-0.0002, 0.00005)];
    let loaded = model
        .prepare_scattering(1700.0, cp, &reference, &mut || false)
        .unwrap();
    let actual = loaded
        .response_with_volume_sources(&incident, &[], &force_sources, &heat_sources, &mut || false)
        .unwrap();
    let direct = prepared
        .response_with_volume_sources(
            &actual.pressure_drive_pa,
            &[],
            &force_sources,
            &heat_sources,
            &mut || false,
        )
        .unwrap();
    for (a, b) in [
        (&actual.pressure_pa, &direct.pressure_pa),
        (&actual.temperature_k, &direct.temperature_k),
    ] {
        let scale = b.iter().map(|z| z.abs()).fold(1e-30, f64::max);
        assert!(a.iter().zip(b).all(|(x, y)| (*x - *y).abs() < 2e-8 * scale));
    }
    let offset = loaded
        .source_outgoing_with_volume_sources(&[], &force_sources, &heat_sources, &mut || false)
        .unwrap();
    let scattering = loaded.scattering(&mut || false).unwrap();
    let mut port_work = 0.0;
    for i in 0..2 {
        let outgoing = (actual.pressure_drive_pa[i]
            - actual.inward_flow_m3_s[i].scale(reference[i]))
        .scale(0.5 / reference[i].sqrt());
        let expected = (0..2).fold(offset[i], |sum, j| {
            sum + scattering[2 * i + j] * incident[j]
        });
        relative(outgoing, expected, 2e-8);
        port_work += 0.5 * (incident[i].norm_sq() - outgoing.norm_sq());
    }
    let scale = actual.power.input_w.abs()
        + actual.power.body_force_work_w.abs()
        + actual.power.heat_source_work_w.abs();
    assert!(
        (actual.power.input_w
            - port_work
            - actual.power.body_force_work_w
            - actual.power.heat_source_work_w)
            .abs()
            < 2e-8 * scale
    );
    assert_work(&actual);
    assert_eq!(
        calls.get(),
        count,
        "warm solves do not resample force profiles"
    );
}

#[test]
fn g1_body_force_retains_constrained_wall_load_and_actual_added_fluid_mass() {
    let size = [0.002, 0.001, 0.001];
    let (mesh, positions) = box_mesh([2, 2, 2], size);
    let medium = ThermoviscousMedium {
        density_kg_m3: 1000.0,
        pressure_pa: 101325.0,
        temperature_k: 300.0,
        cp_j_kg_k: 4000.0,
        dynamic_viscosity_pa_s: 0.001,
        thermal_conductivity_w_m_k: 0.6,
        bulk_viscosity_pa_s: 0.002,
    };
    let thermodynamics = ThermoviscousThermodynamics {
        isothermal_compressibility_pa_inv: 5e-10,
        isobaric_expansion_per_k: 0.0,
    };
    let volume = size.iter().product::<f64>();
    let mass = medium.density_kg_m3 * volume;
    let shape_scale = 1.0 / mass.sqrt();
    let structure = ModalAcousticTimeModel::try_new(
        48000,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: 20.0,
            damping_ratio: 0.05,
            pressure_per_modal_velocity: C64::ZERO,
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    let geometry = HelmholtzFe::new(
        &mesh,
        &positions,
        2,
        Medium {
            density_kg_m3: medium.density_kg_m3,
            sound_speed_m_s: thermodynamics.sound_speed_m_s(medium).unwrap(),
        },
    )
    .unwrap();
    let faces = geometry.boundary_faces();
    let space = SimplexSpace::new(&mesh, 2);
    let mut phi = vec![[0.0; 3]; space.ndof];
    for face in &faces {
        for &vertex in face {
            phi[vertex as usize] = [shape_scale, 0.0, 0.0];
        }
    }
    let model = ThermoviscousFe::new_with_moving_wall(
        &mesh,
        &positions,
        medium,
        thermodynamics,
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
    .unwrap();
    let source = model
        .volume_force_source(
            4,
            &mut |_, _| Ok::<_, ThermoviscousError>([1.0 / volume, 0.0, 0.0]),
            &mut || false,
        )
        .unwrap();
    let force = C64::new(0.001, 0.0003);
    let frequency = 30.0;
    let omega = std::f64::consts::TAU * frequency;
    let field = model
        .prepare(frequency, &mut || false)
        .unwrap()
        .response_with_volume_sources(&[], &[C64::ZERO], &[(&source, force)], &[], &mut || false)
        .unwrap();
    let wall = model.moving_wall().unwrap();
    let impedance = modal_impedance(
        &wall.structure.modes()[0],
        wall.structure.damping_rate_per_s(0),
        omega,
    );
    let expected = force.scale(shape_scale)
        / (impedance + C64::new(0.0, omega * mass * shape_scale * shape_scale));
    let velocity = field.moving_wall.modal_velocity_m_sqrt_kg_per_s[0];
    relative(velocity, expected, 2e-5);
    // No direct structural drive: physical fluid traction alone balances Z_s*u.
    relative(
        field.moving_wall.fluid_force_n_per_sqrt_kg[0],
        impedance * velocity,
        2e-8,
    );
    let expected_work = 0.5 * (force.conj() * velocity.scale(shape_scale)).re;
    assert!(expected_work > 0.0);
    assert!((field.power.body_force_work_w - expected_work).abs() < 2e-5 * expected_work);
    assert!(field.moving_wall.external_input_w == 0.0);
    let fluid_losses = field.power.dissipated_w() - field.power.structural_damping_w;
    let fluid_balance =
        field.power.body_force_work_w + field.moving_wall.power_into_fluid_w - fluid_losses;
    assert!(fluid_balance.abs() < 2e-8 * field.power.body_force_work_w.abs());
    assert_work(&field);
}

#[test]
fn g4_force_profiles_preserve_errors_identity_limits_and_cancellation() {
    #[derive(Debug)]
    enum SourceError {
        Material,
        Interior(ThermoviscousError),
    }
    impl From<ThermoviscousError> for SourceError {
        fn from(error: ThermoviscousError) -> Self {
            Self::Interior(error)
        }
    }
    let size = [0.004, 0.001, 0.0012];
    let (mesh, positions) = box_mesh([1, 1, 1], size);
    let boundary = boundaries(&mesh, &positions, size, true);
    let make = || {
        ThermoviscousFe::new(&mesh, &positions, air(), &boundary, limits(), &mut || false).unwrap()
    };
    let model = make();
    let other = make();
    let mut law = |_, _| Ok::<_, ThermoviscousError>([1.0, -2.0, 3.0]);
    let source = model
        .volume_force_source(4, &mut law, &mut || false)
        .unwrap();
    let foreign = other
        .volume_force_source(4, &mut law, &mut || false)
        .unwrap();
    let prepared = model.prepare(1000.0, &mut || false).unwrap();
    assert!(
        prepared
            .response_with_volume_sources(
                &[C64::ZERO; 2],
                &[],
                &[(&foreign, C64::ONE)],
                &[],
                &mut || false
            )
            .is_err()
    );
    assert!(
        prepared
            .response_with_volume_sources(
                &[C64::ZERO; 2],
                &[],
                &[(&source, C64::new(f64::NAN, 0.0))],
                &[],
                &mut || false
            )
            .is_err()
    );
    assert!(
        prepared
            .response_with_volume_sources(
                &[C64::ZERO; 2],
                &[],
                &vec![(&source, C64::ONE); MAX_BODY_FORCE_SOURCES + 1],
                &[],
                &mut || false
            )
            .is_err()
    );
    for order in [0, 3, 17] {
        assert!(
            model
                .volume_force_source(order, &mut law, &mut || false)
                .is_err()
        );
    }
    assert!(
        model
            .volume_force_source(
                4,
                &mut |_, _| Ok::<_, ThermoviscousError>([f64::NAN, 0.0, 0.0]),
                &mut || false
            )
            .is_err()
    );
    assert!(matches!(
        model.volume_force_source(
            4,
            &mut |_, _| Err::<[f64; 3], _>(SourceError::Material),
            &mut || false
        ),
        Err(SourceError::Material)
    ));
    assert!(matches!(
        model.volume_force_source(4, &mut |_, _| Ok::<_, SourceError>([0.0; 3]), &mut || true),
        Err(SourceError::Interior(ThermoviscousError::Cancelled))
    ));
    let calls = Cell::new(0);
    let result = model.volume_force_source(
        4,
        &mut |_, _| {
            calls.set(calls.get() + 1);
            Ok::<_, ThermoviscousError>([1.0, 0.0, 0.0])
        },
        &mut || calls.get() != 0,
    );
    assert!(matches!(result, Err(ThermoviscousError::Cancelled)));
    assert_eq!(calls.get(), 1);
    assert!(matches!(
        prepared.response_with_volume_sources(
            &[C64::ZERO; 2],
            &[],
            &[(&source, C64::ONE)],
            &[],
            &mut || true
        ),
        Err(ThermoviscousError::Cancelled)
    ));
}
