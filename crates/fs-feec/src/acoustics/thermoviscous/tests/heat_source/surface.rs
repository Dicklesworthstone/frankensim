//! Physical surface fluxes through the unchanged full thermal/acoustic owners.
use super::*;
use crate::acoustics::helmholtz::QuadraticGeometryControls;
use crate::acoustics::helmholtz::boundary_layer::BoundaryLayerGeometry;
use crate::acoustics::helmholtz::boundary_layer::thermal_wall::{
    ConductiveThermalWall, ThermalBacking, ThermalWallDepth, WallThermalTransport,
};
use crate::acoustics::thermoviscous::moving_wall::MovingWall;
use crate::acoustics::thermoviscous::time::{
    ThermoviscousCaloricReservoir, ThermoviscousTimeLimits, ThermoviscousWallDiscretization,
};
use crate::acoustics::thermoviscous::wall::ConductiveWallRegion;
use fs_modal::time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance * expected.abs().max(1e-25),
        "actual={actual:e}, expected={expected:e}"
    );
}

fn same_values(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    let magnitude = expected.iter().map(|x| x.abs()).fold(1e-25, f64::max);
    assert!(
        actual
            .iter()
            .zip(expected)
            .all(|(a, b)| (a - b).abs() < 2e-8 * magnitude)
    );
}

/// Every boundary face (a face of exactly one tetrahedron) as a no-slip wall.
fn all_walls(mesh: &TetComplex) -> ThermoviscousBoundary {
    let no_slip_faces = mesh
        .faces
        .iter()
        .copied()
        .filter(|face| {
            mesh.tets
                .iter()
                .filter(|tet| face.iter().all(|v| tet.contains(v)))
                .count()
                == 1
        })
        .collect();
    ThermoviscousBoundary {
        no_slip_faces,
        isothermal_faces: vec![],
        ports: vec![],
    }
}

/// `fixtures::single_tet` with one interior vertex (index 4, its centroid)
/// joined to the four corners: four tetrahedra on the same outer boundary.
/// With every outer face a wall, a lone tetrahedron has no free velocity or
/// temperature node, which the thermoviscous owner refuses; the interior
/// vertex supplies them without changing any boundary face or its geometry.
fn centered_tet() -> (TetComplex, Vec<[f64; 3]>) {
    let (_, mut positions) = crate::fixtures::single_tet();
    let centroid = std::array::from_fn(|axis| positions.iter().map(|p| p[axis]).sum::<f64>() / 4.0);
    positions.push(centroid);
    let mesh = TetComplex::from_tets(
        5,
        vec![[4, 1, 2, 3], [0, 4, 2, 3], [0, 1, 4, 3], [0, 1, 2, 4]],
    );
    (mesh, positions)
}

fn constant_shape(model: &ThermoviscousFe<'_>, faces: &[[u32; 3]], density: f64) -> Vec<f64> {
    let mut shape = vec![0.0; model.space.ndof];
    for face in faces {
        for &v in face {
            shape[v as usize] = density;
        }
    }
    shape
}

#[test]
fn g1_surface_flux_retains_midside_moments_signed_shapes_and_strict_admission() {
    let (mesh, positions) = centered_tet();
    let boundary = all_walls(&mesh);
    let make = || {
        ThermoviscousFe::new(&mesh, &positions, air(), &boundary, limits(), &mut || false).unwrap()
    };
    let model = make();
    let face = [0, 1, 2];
    let edge = model.space.edge_off + mesh.edges.binary_search(&[0, 1]).unwrap();
    let mut shape = vec![0.0; model.space.ndof];
    // Hierarchical edge basis is lambda0*lambda1: coefficient FOUR makes
    // a unit midpoint value with zero at every vertex.
    shape[edge] = 4.0;
    let source = model
        .surface_heat_source(&[face], &shape, &mut || false)
        .unwrap();
    close(source.integrated_shape(), 1.0 / 6.0, 2e-14);
    for (vertex, expected) in [(0, 1.0 / 15.0), (1, 1.0 / 15.0), (2, 1.0 / 30.0)] {
        close(
            source.load[model.temperature_map[vertex].unwrap()],
            expected,
            2e-14,
        );
    }
    close(
        source.load[model.temperature_map[edge].unwrap()],
        1.0 / 45.0,
        2e-14,
    );
    let reversed = model
        .surface_heat_source(&[[2, 0, 1]], &shape, &mut || false)
        .unwrap();
    assert_eq!(source.load, reversed.load);
    let mut signed = vec![0.0; model.space.ndof];
    signed[0] = -1.0;
    signed[1] = 1.0;
    let differential = model
        .surface_heat_source(&[face], &signed, &mut || false)
        .unwrap();
    assert!(differential.integrated_shape().abs() < 2e-14);
    assert!(differential.load.iter().any(|x| x.abs() > 1e-3));
    assert!(
        model
            .surface_heat_source(&[], &shape, &mut || false)
            .is_err()
    );
    assert!(
        model
            .surface_heat_source(&[face, [2, 1, 0]], &shape, &mut || false)
            .is_err()
    );
    assert!(
        model
            .surface_heat_source(&[[0, 1, 99]], &shape, &mut || false)
            .is_err()
    );
    assert!(
        model
            .surface_heat_source(&[face], &shape[1..], &mut || false)
            .is_err()
    );
    let mut bad = shape.clone();
    bad[3] = 1.0;
    assert!(
        model
            .surface_heat_source(&[face], &bad, &mut || false)
            .is_err()
    );
    bad[3] = f64::NAN;
    assert!(
        model
            .surface_heat_source(&[face], &bad, &mut || false)
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        model.surface_heat_source(&[face], &shape, &mut || {
            polls += 1;
            polls == 5
        }),
        Err(ThermoviscousError::Cancelled)
    ));
    assert_eq!(
        source.load,
        model
            .surface_heat_source(&[face], &shape, &mut || false)
            .unwrap()
            .load
    );
    let other = make();
    assert!(
        other
            .prepare(1000.0, &mut || false)
            .unwrap()
            .response_with_heat_sources(&[], &[], &[(&source, C64::ZERO)], &mut || false)
            .is_err()
    );
    let isothermal = ThermoviscousFe::new(
        &mesh,
        &positions,
        air(),
        &ThermoviscousBoundary {
            isothermal_faces: boundary.no_slip_faces.clone(),
            ..boundary.clone()
        },
        limits(),
        &mut || false,
    )
    .unwrap();
    assert!(
        isothermal
            .surface_heat_source(&[face], &shape, &mut || false)
            .is_err()
    );
    // Adjacent exact-isothermal faces can eliminate this entire trace even
    // when the selected face itself was not explicitly listed as isothermal.
    let constrained = ThermoviscousFe::new(
        &mesh,
        &positions,
        air(),
        &ThermoviscousBoundary {
            isothermal_faces: boundary
                .no_slip_faces
                .iter()
                .copied()
                .filter(|f| *f != face)
                .collect(),
            ..boundary
        },
        limits(),
        &mut || false,
    )
    .unwrap();
    assert!(
        constrained
            .surface_heat_source(&[face], &shape, &mut || false)
            .is_err()
    );
}

#[test]
fn g1_curved_surface_flux_uses_the_original_variable_measure() {
    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 73,
                kernel_id: 5,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let (mesh, original) = centered_tet();
        let a = 0.6;
        let map = |[x, y, z]: [f64; 3]| [x, y, z * (1.0 + a * x)];
        let positions: Vec<_> = original.iter().copied().map(map).collect();
        let mids: Vec<_> = mesh
            .edges
            .iter()
            .map(|&[i, j]| {
                map(std::array::from_fn(|axis| {
                    f64::midpoint(original[i as usize][axis], original[j as usize][axis])
                }))
            })
            .collect();
        let gas = air();
        let geometry = HelmholtzFe::new_quadratic(
            &mesh,
            &positions,
            &mids,
            2,
            6,
            Medium {
                density_kg_m3: gas.density_kg_m3,
                sound_speed_m_s: ThermoviscousThermodynamics::ideal_gas(gas)
                    .unwrap()
                    .sound_speed_m_s(gas)
                    .unwrap(),
            },
            QuadraticGeometryControls::default(),
            &cx,
        )
        .unwrap();
        let model =
            ThermoviscousFe::from_geometry(geometry, gas, &all_walls(&mesh), limits(), &mut || {
                false
            })
            .unwrap();
        let mut shape = vec![0.0; model.space.ndof];
        shape[model.space.edge_off + mesh.edges.binary_search(&[0, 1]).unwrap()] = 4.0;
        let source = model
            .surface_heat_source(&[[0, 1, 3]], &shape, &mut || false)
            .unwrap();
        // Physical y=0 surface measure is (1+a*x) dX dZ. Its independent
        // simplex integral against 4*lambda0*lambda1 is 1/6+a/15.
        close(source.integrated_shape(), 1.0 / 6.0 + a / 15.0, 3e-13);
        let vertices = (0..mesh.vertex_count)
            .map(|v| source.load[model.temperature_map[v].unwrap()])
            .sum::<f64>();
        close(vertices, source.integrated_shape(), 2e-14);
    });
}

#[test]
fn g1_signed_surface_heating_recovers_fourier_gradient_and_absolute_work() {
    let size = [0.002, 0.001, 0.001];
    let area = size[1] * size[2];
    let volume = size.iter().product::<f64>();
    let gas = air();
    // A declared zero-expansion fluid isolates the exact Fourier boundary law;
    // it is a numerical constitutive fixture, not an air material claim.
    let thermo = ThermoviscousThermodynamics {
        isothermal_compressibility_pa_inv: 1.0 / gas.pressure_pa,
        isobaric_expansion_per_k: 0.0,
    };
    let omega = 35.0;
    let amplitude = C64::new(0.001, 0.0003);
    for cells in [1, 2, 4] {
        let (mesh, positions) = box_mesh([cells, 1, 1], size);
        let boundary = sealed(&mesh, &positions, size);
        let model = ThermoviscousFe::new_with_thermodynamics(
            &mesh,
            &positions,
            gas,
            thermo,
            &boundary,
            limits(),
            &mut || false,
        )
        .unwrap();
        let end = |x: f64| x.abs() < 1e-12 * size[0] || (x - size[0]).abs() < 1e-12 * size[0];
        let faces: Vec<_> = boundary
            .no_slip_faces
            .iter()
            .copied()
            .filter(|face| {
                face.iter().all(|&v| end(positions[v as usize][0]))
                    && face.iter().all(|&v| {
                        (positions[v as usize][0] - positions[face[0] as usize][0]).abs()
                            < 1e-12 * size[0]
                    })
            })
            .collect();
        let mut shape = vec![0.0; model.space.ndof];
        for face in &faces {
            for &v in face {
                shape[v as usize] = if positions[v as usize][0] < size[0] / 2.0 {
                    -1.0 / area
                } else {
                    1.0 / area
                };
            }
        }
        let surface = model
            .surface_heat_source(&faces, &shape, &mut || false)
            .unwrap();
        assert!(surface.integrated_shape().abs() < 3e-14);
        let storage = model
            .volume_heat_source(
                4,
                &mut |_, x| Ok((1.0 + x[0] / size[0]) / volume),
                &mut || false,
            )
            .unwrap();
        // T'=A*(1+x/L), v=p=0. Interior q'=i*omega*rho*Cp*T';
        // inward flux is kappa*dT'/dn, negative at x=0 and positive at x=L.
        let field = model
            .prepare(omega / std::f64::consts::TAU, &mut || false)
            .unwrap()
            .response_with_heat_sources(
                &[],
                &[],
                &[
                    (
                        &surface,
                        amplitude.scale(gas.thermal_conductivity_w_m_k * area / size[0]),
                    ),
                    (
                        &storage,
                        amplitude
                            * C64::new(0.0, omega * gas.density_kg_m3 * gas.cp_j_kg_k * volume),
                    ),
                ],
                &mut || false,
            )
            .unwrap();
        for (v, &actual) in field.temperature_k.iter().enumerate() {
            let expected = if v < mesh.vertex_count {
                amplitude.scale(1.0 + positions[v][0] / size[0])
            } else {
                C64::ZERO
            };
            assert!((actual - expected).abs() < 2e-9 * amplitude.abs());
        }
        assert!(field.pressure_pa.iter().all(|p| p.abs() < 1e-12));
        assert!(field.velocity_m_s.iter().flatten().all(|v| v.abs() < 1e-12));
        let expected_work = gas.thermal_conductivity_w_m_k * area * amplitude.norm_sq()
            / (2.0 * gas.temperature_k * size[0]);
        close(field.power.heat_source_work_w, expected_work, 2e-8);
        close(field.power.conductive_w, expected_work, 2e-8);
        assert_work(&field);
    }
}

#[test]
fn g3_surface_heat_emits_through_original_loaded_acoustic_ports() {
    let size = [0.004, 0.001, 0.0012];
    let (mesh, positions) = box_mesh([1, 2, 2], size);
    let boundary = boundaries(&mesh, &positions, size, false);
    let model =
        ThermoviscousFe::new(&mesh, &positions, air(), &boundary, limits(), &mut || false).unwrap();
    let faces = boundary.no_slip_faces[..1].to_vec();
    let (_, _, area) = model.fe.boundary_geometry(faces[0]).unwrap();
    let shape = constant_shape(&model, &faces, 1.0 / area);
    let source = model
        .surface_heat_source(&faces, &shape, &mut || false)
        .unwrap();
    close(source.integrated_shape(), 1.0, 3e-14);
    // A traction port is not a wall-temperature source owner.
    assert!(
        model
            .surface_heat_source(&[boundary.ports[0].traces[0].0], &shape, &mut || false)
            .is_err()
    );
    let cp = C64::new(970.0, -45.0);
    let reference: [f64; 2] = [1.3e8, 0.7e8];
    let sources = [(&source, C64::new(0.0003, -0.0002))];
    let unloaded = model
        .prepare_with_heat_capacity(1700.0, cp, &mut || false)
        .unwrap();
    let source_only = unloaded
        .response_with_heat_sources(&[C64::ZERO; 2], &[], &sources, &mut || false)
        .unwrap();
    let admittance = model
        .admittance_with_heat_capacity(1700.0, cp, &mut || false)
        .unwrap();
    // Independent physical termination P+R*(Y*P+Q_heat)=0, with no incident wave.
    let mut matrix = vec![C64::ZERO; 4];
    let mut pressure = vec![C64::ZERO; 2];
    for i in 0..2 {
        pressure[i] = source_only.inward_flow_m3_s[i].scale(-reference[i]);
        for j in 0..2 {
            matrix[2 * i + j] = admittance[2 * i + j].scale(reference[i])
                + C64::from_re(if i == j { 1.0 } else { 0.0 });
        }
    }
    lu_complex_with_cancel(&matrix, 2, || false)
        .unwrap()
        .solve_with_cancel(&mut pressure, || false)
        .unwrap();
    let expected = unloaded
        .response_with_heat_sources(&pressure, &[], &sources, &mut || false)
        .unwrap();
    let loaded = model
        .prepare_scattering(1700.0, cp, &reference, &mut || false)
        .unwrap();
    let actual = loaded
        .response_with_sources(&[C64::ZERO; 2], &[], &sources, &mut || false)
        .unwrap();
    same_resolved_fields(&actual, &expected);
    assert_work(&actual);
    let outgoing = loaded
        .source_outgoing(&[], &sources, &mut || false)
        .unwrap();
    let emitted = 0.5 * outgoing.iter().map(|a| a.norm_sq()).sum::<f64>();
    assert!(emitted > 0.0);
    close(
        actual.power.heat_source_work_w,
        actual.power.dissipated_w() + emitted,
        2e-8,
    );
    let zero = loaded
        .response_with_sources(
            &[C64::ZERO; 2],
            &[],
            &[sources[0], (&source, sources[0].1.scale(-1.0))],
            &mut || false,
        )
        .unwrap();
    assert!(zero.temperature_k.iter().all(|t| *t == C64::ZERO));
    assert_eq!(zero.power.heat_source_work_w, 0.0);
}

#[test]
fn g3_surface_heat_retains_moving_wall_molecular_and_solid_fields_in_time() {
    let length: f64 = 0.002;
    let (mesh, reference) = crate::fixtures::single_tet();
    let positions: Vec<_> = reference.iter().map(|p| p.map(|x| x * length)).collect();
    let volume = length.powi(3) / 6.0;
    let gas = air();
    let boundary = all_walls(&mesh);
    let scale = 1.0 / det::sqrt(gas.density_kg_m3 * volume);
    let mut shape = vec![[0.0; 3]; SimplexSpace::new(&mesh, 2).ndof];
    for (coefficient, position) in shape.iter_mut().take(mesh.vertex_count).zip(&positions) {
        for axis in 0..3 {
            coefficient[axis] = scale * (position[axis] / length - 0.25);
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
    let face = [0, 1, 2];
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
            faces: vec![face],
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
    let surface_shape = constant_shape(&model, &[face], 2.0 / length.powi(2));
    let source = model
        .surface_heat_source(&[face], &surface_shape, &mut || false)
        .unwrap();
    let reservoirs = [ThermoviscousCaloricReservoir {
        heat_capacity_j_kg_k: 50.0,
        relaxation_time_s: 2e-5,
    }];
    let frozen_cp = gas.cp_j_kg_k - 50.0;
    let frequency = 800.0;
    let cp = C64::from_re(frozen_cp)
        + C64::new(1.0, std::f64::consts::TAU * frequency * 2e-5)
            .recip()
            .scale(50.0);
    let harmonic = model
        .prepare_with_heat_capacity(frequency, cp, &mut || false)
        .unwrap()
        .response_with_heat_sources(
            &[],
            &[C64::ZERO],
            &[(&source, C64::new(0.001, 0.0003))],
            &mut || false,
        )
        .unwrap();
    assert_work(&harmonic);
    assert!(harmonic.moving_wall.modal_velocity_m_sqrt_kg_per_s[0].abs() > 0.0);
    assert!(harmonic.power.solid_conductive_w > 0.0 && harmonic.power.caloric_w > 0.0);
    assert!(harmonic.power.structural_damping_w > 0.0);
    let depth = [2];
    let prepare = |duration_s| {
        model
            .prepare_time_with_heat_sources(
                duration_s,
                frozen_cp,
                &reservoirs,
                Some(ThermoviscousWallDiscretization {
                    normal_elements: &depth,
                    maximum_frequency_hz: 2000.0,
                }),
                None,
                &[&source],
                ThermoviscousTimeLimits::default(),
                &mut || false,
            )
            .unwrap()
    };
    let dt = 2e-5;
    let full = prepare(dt);
    let half = prepare(0.5 * dt);
    let initial = full.zero_state();
    let saved = initial.clone();
    let whole = full
        .step(&initial, &[], &[0.0], &[0.001], &mut || false)
        .unwrap();
    let first = half
        .step(&initial, &[], &[0.0], &[0.001], &mut || false)
        .unwrap();
    let second = half
        .step(&first.field.state, &[], &[0.0], &[0.001], &mut || false)
        .unwrap();
    assert_eq!(initial, saved);
    let actual = &whole.field.state;
    let expected = &second.field.state;
    for (a, b) in [
        (&actual.pressure_pa, &expected.pressure_pa),
        (&actual.temperature_k, &expected.temperature_k),
        (
            &actual.wall_displacement_m_sqrt_kg,
            &expected.wall_displacement_m_sqrt_kg,
        ),
        (
            &actual.wall_velocity_m_sqrt_kg_per_s,
            &expected.wall_velocity_m_sqrt_kg_per_s,
        ),
    ] {
        same_values(a, b);
    }
    same_values(
        &actual
            .velocity_m_s
            .iter()
            .flatten()
            .copied()
            .collect::<Vec<_>>(),
        &expected
            .velocity_m_s
            .iter()
            .flatten()
            .copied()
            .collect::<Vec<_>>(),
    );
    for (a, b) in actual
        .molecular_temperature_k
        .iter()
        .zip(&expected.molecular_temperature_k)
        .chain(
            actual
                .conductive_wall_temperature_k
                .iter()
                .zip(&expected.conductive_wall_temperature_k),
        )
    {
        same_values(a, b);
    }
    assert!(
        whole.field.state.molecular_temperature_k[0]
            .iter()
            .any(|x| x.abs() > 0.0)
    );
    assert!(
        whole.field.state.conductive_wall_temperature_k[0]
            .iter()
            .any(|x| x.abs() > 0.0)
    );
    assert!(whole.field.state.wall_velocity_m_sqrt_kg_per_s[0].abs() > 0.0);
    close(
        whole.heat_source_work_by_source_j[0],
        first.heat_source_work_by_source_j[0] + second.heat_source_work_by_source_j[0],
        2e-8,
    );
    assert!(
        whole.field.work.fluid_conductive_heat_j > 0.0
            && whole.field.work.solid_conductive_heat_j > 0.0
    );
    assert!(
        whole.field.work.molecular_relaxation_heat_j > 0.0
            && whole.field.work.structural_damping_heat_j > 0.0
    );
    assert!(
        whole.field.work.balance_residual_j().abs() < 2e-8 * whole.field.work.supplied_j().abs()
    );
}
