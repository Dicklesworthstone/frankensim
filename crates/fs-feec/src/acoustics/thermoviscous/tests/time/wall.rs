use super::*;
use crate::acoustics::helmholtz::boundary_layer::BoundaryLayerGeometry;
use crate::acoustics::helmholtz::boundary_layer::thermal_wall::{
    ConductiveThermalWall, ThermalBacking, ThermalWallDepth, WallThermalTransport,
};
use crate::acoustics::thermoviscous::time::{
    PreparedThermoviscousTime, ThermoviscousCaloricReservoir, ThermoviscousWallDiscretization,
};
use crate::acoustics::thermoviscous::wall::ConductiveWallRegion;

fn region(faces: Vec<[u32; 3]>, temperature: f64, backing: ThermalBacking) -> ConductiveWallRegion {
    ConductiveWallRegion {
        faces,
        wall: ConductiveThermalWall {
            // Manufactured transport, not a measured material assignment.
            transport: WallThermalTransport {
                conductivity_w_m_k: 0.003,
                volumetric_heat_capacity_j_m3_k: 2000.0,
            },
            depth: ThermalWallDepth::Slab {
                thickness_m: 0.00002,
                backing,
            },
            minimum_tangential_scale_m: 0.004,
        },
        geometry: BoundaryLayerGeometry {
            minimum_curvature_radius_m: f64::INFINITY,
            minimum_half_gap_m: 0.0001,
            maximum_layer_ratio: 0.05,
        },
        temperature_k: temperature,
    }
}

fn tetra() -> (TetComplex, [[f64; 3]; 4]) {
    (
        TetComplex::from_tets(4, vec![[0, 1, 2, 3]]),
        [
            [0.0; 3],
            [0.002, 0.0, 0.0],
            [0.0, 0.001, 0.0],
            [0.0, 0.0, 0.001],
        ],
    )
}

fn resolution(elements: &[usize]) -> ThermoviscousWallDiscretization<'_> {
    ThermoviscousWallDiscretization {
        normal_elements: elements,
        maximum_frequency_hz: 2000.0,
    }
}

#[test]
fn g1_uniform_insulated_fluid_molecule_and_solid_keep_the_complete_heat_capacity() {
    let (mesh, pos) = tetra();
    let boundary = ThermoviscousBoundary {
        no_slip_faces: mesh.faces.clone(),
        isothermal_faces: vec![],
        ports: vec![],
    };
    let model = ThermoviscousFe::new(&mesh, &pos, air(), &boundary, limits(), &mut || false)
        .unwrap()
        .with_conductive_walls(
            &[region(
                vec![mesh.faces[0]],
                air().temperature_k,
                ThermalBacking::Adiabatic,
            )],
            &mut || false,
        )
        .unwrap();
    let reservoir = ThermoviscousCaloricReservoir {
        heat_capacity_j_kg_k: 100.0,
        relaxation_time_s: 0.001,
    };
    let prepared = model
        .prepare_time_with_thermal_fields(
            1e-4,
            air().cp_j_kg_k - 100.0,
            &[reservoir],
            resolution(&[2]),
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let mut state = prepared.zero_state();
    state.temperature_k[..4].fill(1.0);
    state.molecular_temperature_k[0][..4].fill(1.0);
    let layout = &prepared.conductive_wall_layouts()[0];
    for depth in state.conductive_wall_temperature_k[0].chunks_mut(layout.trace_dofs().len()) {
        for (value, &dof) in depth.iter_mut().zip(layout.trace_dofs()) {
            *value = if dof < 4 { 1.0 } else { 0.0 };
        }
    }
    let (_, _, area) = model.fe.boundary_geometry(mesh.faces[0]).unwrap();
    let volume = model.geometry.vol_signed[0].abs();
    let expected = (air().density_kg_m3 * air().cp_j_kg_k * volume + 2000.0 * 0.00002 * area)
        / (2.0 * air().temperature_k);
    assert!(relative(prepared.stored_energy_j(&state).unwrap(), expected) < 2e-14);
    let frame = prepared.step(&state, &[], &[], &mut || false).unwrap();
    assert!(relative(frame.work.final_energy_j, expected) < 2e-11);
    assert!(frame.work.dissipated_j().abs() < 2e-10 * expected);
    for (actual, before) in frame.state.temperature_k.iter().zip(&state.temperature_k) {
        assert!((actual - before).abs() < 2e-10);
    }
}

#[test]
fn g1_nonuniform_surface_and_depth_energy_preserves_consistent_interface_cross_mass() {
    use crate::acoustics::helmholtz::tri_quadrature;
    use crate::highorder::simplex::LocalFn;
    let (mesh, pos) = tetra();
    let face = mesh.faces[0];
    let boundary = ThermoviscousBoundary {
        no_slip_faces: mesh.faces.clone(),
        isothermal_faces: vec![],
        ports: vec![],
    };
    let model = ThermoviscousFe::new(&mesh, &pos, air(), &boundary, limits(), &mut || false)
        .unwrap()
        .with_conductive_walls(
            &[region(
                vec![face],
                air().temperature_k,
                ThermalBacking::Adiabatic,
            )],
            &mut || false,
        )
        .unwrap();
    let prepared = model
        .prepare_time_with_thermal_fields(
            1e-5,
            air().cp_j_kg_k,
            &[],
            resolution(&[2]),
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let layout = &prepared.conductive_wall_layouts()[0];
    assert_eq!(layout.trace_dofs().len(), 6);
    let mut state = prepared.zero_state();
    for (i, value) in state.temperature_k.iter_mut().enumerate() {
        *value = ((i * 7 + 3) % 11) as f64 / 8.0 - 0.6;
    }
    for (i, value) in state.conductive_wall_temperature_k[0]
        .iter_mut()
        .enumerate()
    {
        *value = ((i * 5 + 2) % 13) as f64 / 7.0 - 0.8;
    }
    // Independent pointwise P2 x P1 integration, not the assembled matrix.
    let face = canonical(face);
    let mut dofs = face.map(|v| v as usize).to_vec();
    for (a, b) in [(0, 1), (0, 2), (1, 2)] {
        dofs.push(model.space.edge_off + mesh.edges.binary_search(&[face[a], face[b]]).unwrap());
    }
    let shapes = [
        LocalFn::Vertex(0),
        LocalFn::Vertex(1),
        LocalFn::Vertex(2),
        LocalFn::Edge(0, 1, 0),
        LocalFn::Edge(0, 2, 0),
        LocalFn::Edge(1, 2, 0),
    ];
    let (_, _, area) = model.fe.boundary_geometry(face).unwrap();
    let mut solid = 0.0;
    for (bary, weight) in tri_quadrature(3) {
        let values: Vec<_> = shapes
            .iter()
            .map(|s| s.eval([bary[0], bary[1], bary[2], 0.0], 2))
            .collect();
        let at_depth = |depth: usize| -> f64 {
            dofs.iter()
                .zip(&values)
                .map(|(&dof, &shape)| {
                    let value = if depth == 0 {
                        state.temperature_k[dof]
                    } else {
                        let index = layout.trace_dofs().binary_search(&dof).unwrap();
                        state.conductive_wall_temperature_k[0][(depth - 1) * 6 + index]
                    };
                    value * shape
                })
                .sum()
        };
        for depth in 0..2 {
            let a = at_depth(depth);
            let b = at_depth(depth + 1);
            for z in [0.5 - 0.5 / det::sqrt(3.0), 0.5 + 0.5 / det::sqrt(3.0)] {
                let t = (1.0 - z) * a + z * b;
                solid += 2000.0 / air().temperature_k * area * weight * 0.00001 * 0.5 * t * t;
            }
        }
    }
    let mut fluid = 0.0;
    for (i, &ti) in state.temperature_k.iter().enumerate() {
        let row = model.temperature_start + model.temperature_map[i].unwrap();
        for (j, &tj) in state.temperature_k.iter().enumerate() {
            let col = model.temperature_start + model.temperature_map[j].unwrap();
            fluid += 0.5 * ti * model.mass[row * model.n + col] * tj / air().temperature_k.powi(2);
        }
    }
    assert!(relative(prepared.stored_energy_j(&state).unwrap(), fluid + solid) < 2e-13);
}

// Query the ACTUAL prepared pHS generator, using its full source Q/J/R/G.
// This checks the production assembly without a duplicate fluid/wall solver.
fn harmonic_flow(prepared: &PreparedThermoviscousTime<'_, '_>, frequency: f64) -> C64 {
    let system = prepared.system();
    let n = system.state_dim();
    let (j, r, g) = system.structure();
    assert_eq!(system.port_dim(), 1);
    let mut a = vec![C64::ZERO; n * n];
    for col in 0..n {
        let mut basis = vec![0.0; n];
        basis[col] = 1.0;
        let effort = system.effort(&basis);
        for row in 0..n {
            a[row * n + col].re = (0..n)
                .map(|k| (r[row * n + k] - j[row * n + k]) * effort[k])
                .sum();
        }
        a[col * n + col].im = std::f64::consts::TAU * frequency;
    }
    let factor = lu_complex_with_cancel(&a, n, &mut || false).unwrap();
    let mut x: Vec<_> = g.iter().map(|&v| C64::from_re(v)).collect();
    factor.solve(&mut x);
    C64::new(
        system.output(&x.iter().map(|z| z.re).collect::<Vec<_>>())[0],
        system.output(&x.iter().map(|z| z.im).collect::<Vec<_>>())[0],
    )
}

#[test]
fn g3_actual_full_time_generator_converges_to_both_original_harmonic_slab_backings() {
    let (mesh, pos) = tetra();
    let wall_face = mesh.faces[0];
    let ports: Vec<_> = mesh
        .faces
        .iter()
        .copied()
        .filter(|f| *f != wall_face)
        .collect();
    let boundary = ThermoviscousBoundary {
        no_slip_faces: vec![wall_face],
        isothermal_faces: vec![],
        ports: vec![BoundaryMode::uniform(&ports)],
    };
    let frequency = 1200.0;
    for backing in [ThermalBacking::Isothermal, ThermalBacking::Adiabatic] {
        let model = ThermoviscousFe::new(&mesh, &pos, air(), &boundary, limits(), &mut || false)
            .unwrap()
            .with_conductive_walls(
                &[region(vec![wall_face], air().temperature_k, backing)],
                &mut || false,
            )
            .unwrap();
        let expected = model
            .prepare(frequency, &mut || false)
            .unwrap()
            .response(&[C64::ONE], &mut || false)
            .unwrap()
            .inward_flow_m3_s[0];
        let mut errors = Vec::new();
        for elements in [1, 2, 4, 8] {
            let prepared = model
                .prepare_time_with_thermal_fields(
                    1.0 / 48_000.0,
                    air().cp_j_kg_k,
                    &[],
                    resolution(&[elements]),
                    ThermoviscousTimeLimits::default(),
                    &mut || false,
                )
                .unwrap();
            let actual = harmonic_flow(&prepared, frequency);
            errors.push((actual - expected).abs() / expected.abs());
        }
        for pair in errors.windows(2) {
            assert!(pair[1] < 0.4 * pair[0], "{backing:?}: {errors:?}");
        }
        assert!(errors[3] < 2e-4, "{backing:?}: {errors:?}");
    }
}

#[test]
fn g1_finite_solid_molecules_and_original_moving_wall_share_one_state_and_work_window() {
    let size = [0.002, 0.001, 0.001];
    let (mesh, pos) = box_mesh([1, 1, 1], size);
    let density = 1000.0;
    let volume = size.iter().product::<f64>();
    let amplitude = 1.0 / det::sqrt(density * volume);
    let raw = translating(&mesh, &pos, density, amplitude, 0.0, 0.0);
    let face = raw.fe.boundary_faces()[0];
    let model = raw
        .with_conductive_walls(
            &[region(vec![face], 300.0, ThermalBacking::Adiabatic)],
            &mut || false,
        )
        .unwrap();
    let reservoir = ThermoviscousCaloricReservoir {
        heat_capacity_j_kg_k: 100.0,
        relaxation_time_s: 0.001,
    };
    let dt = 1e-4;
    let prepared = model
        .prepare_time_with_thermal_fields(
            dt,
            model.medium.cp_j_kg_k - 100.0,
            &[reservoir],
            resolution(&[2]),
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let mut initial = prepared.zero_state();
    let acceleration = 1.0 / (1.0 + density * volume * amplitude * amplitude);
    for (p, xyz) in initial.pressure_pa.iter_mut().zip(&pos) {
        *p = -density * amplitude * acceleration * (xyz[0] - 0.5 * size[0]);
    }
    initial.molecular_temperature_k[0][0] = 0.4;
    initial.conductive_wall_temperature_k[0][0] = 1.0;
    initial.wall_displacement_m_sqrt_kg[0] = 1e-6;
    let before = initial.clone();
    let step = prepared.step(&initial, &[], &[1.0], &mut || false).unwrap();
    assert_eq!(initial, before);
    assert!(
        relative(
            step.state.wall_velocity_m_sqrt_kg_per_s[0],
            acceleration * dt
        ) < 2e-8
    );
    assert!(
        relative(
            step.state.wall_displacement_m_sqrt_kg[0] - 1e-6,
            0.5 * acceleration * dt * dt
        ) < 2e-8
    );
    assert!(step.work.solid_conductive_heat_j > 0.0);
    assert!(step.work.molecular_relaxation_heat_j > 0.0);
    assert!(step.state.temperature_k.iter().any(|t| *t != 0.0));
    assert!(step.work.balance_residual_j().abs() < 1e-10 * step.work.initial_energy_j);
}

#[test]
fn g0_clamped_front_trace_keeps_interior_solid_coefficients_and_refuses_incomplete_states() {
    let (mesh, pos) = tetra();
    let boundary = ThermoviscousBoundary {
        no_slip_faces: mesh.faces.clone(),
        isothermal_faces: vec![mesh.faces[1]],
        ports: vec![],
    };
    let model = ThermoviscousFe::new(&mesh, &pos, air(), &boundary, limits(), &mut || false)
        .unwrap()
        .with_conductive_walls(
            &[region(
                vec![mesh.faces[0]],
                air().temperature_k,
                ThermalBacking::Isothermal,
            )],
            &mut || false,
        )
        .unwrap();
    let prepared = model
        .prepare_time_with_thermal_fields(
            1e-5,
            air().cp_j_kg_k,
            &[],
            resolution(&[2]),
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let layout = &prepared.conductive_wall_layouts()[0];
    assert_eq!(layout.trace_dofs().len(), 6);
    assert_eq!(layout.depths_m(), &[1e-5]);
    let fixed = layout
        .trace_dofs()
        .iter()
        .position(|&d| model.temperature_map[d].is_none())
        .unwrap();
    let mut state = prepared.zero_state();
    state.conductive_wall_temperature_k[0][fixed] = 1.0;
    let original = state.clone();
    let step = prepared.step(&state, &[], &[], &mut || false).unwrap();
    assert!(step.work.solid_conductive_heat_j > 0.0);
    assert!(step.work.final_energy_j < step.work.initial_energy_j);
    assert_eq!(state, original);
    state.conductive_wall_temperature_k[0].pop();
    assert!(prepared.step(&state, &[], &[], &mut || false).is_err());
    assert!(
        model
            .prepare_time_with_thermal_fields(
                1e-5,
                air().cp_j_kg_k,
                &[],
                resolution(&[0]),
                ThermoviscousTimeLimits::default(),
                &mut || false
            )
            .is_err()
    );
    assert!(
        model
            .prepare_time_with_thermal_fields(
                1e-5,
                air().cp_j_kg_k,
                &[],
                resolution(&[]),
                ThermoviscousTimeLimits::default(),
                &mut || false
            )
            .is_err()
    );
    assert!(
        model
            .prepare_time_with_thermal_fields(
                1e-5,
                air().cp_j_kg_k,
                &[],
                resolution(&[128]),
                ThermoviscousTimeLimits::default(),
                &mut || false
            )
            .is_err()
    );
    let before_cancel = original.clone();
    let mut cancelled = 0;
    assert!(matches!(
        prepared.step(&original, &[], &[], &mut || {
            cancelled += 1;
            cancelled == 2
        }),
        Err(ThermoviscousTimeError::Fluid(ThermoviscousError::Cancelled))
    ));
    assert_eq!(original, before_cancel);
    let mut halfspace = model.conductive_walls()[0].clone();
    halfspace.wall.depth = ThermalWallDepth::HalfSpace;
    let model = model
        .with_conductive_walls(&[halfspace], &mut || false)
        .unwrap();
    assert!(
        model
            .prepare_time_with_thermal_fields(
                1e-5,
                air().cp_j_kg_k,
                &[],
                resolution(&[2]),
                ThermoviscousTimeLimits::default(),
                &mut || false
            )
            .is_err()
    );
}
