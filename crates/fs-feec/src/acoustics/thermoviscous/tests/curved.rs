use super::*;
use crate::acoustics::helmholtz::QuadraticGeometryControls;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_modal::time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};

fn run(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 19,
                kernel_id: 95,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

fn acoustic(medium: ThermoviscousMedium, thermo: ThermoviscousThermodynamics) -> Medium {
    Medium {
        density_kg_m3: medium.density_kg_m3,
        sound_speed_m_s: thermo.sound_speed_m_s(medium).unwrap(),
    }
}

fn mapped(
    mesh: &TetComplex,
    original: &[[f64; 3]],
    map: impl Fn([f64; 3]) -> [f64; 3],
) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let positions = original.iter().copied().map(&map).collect();
    let midpoints = mesh
        .edges
        .iter()
        .map(|&[a, b]| {
            map(std::array::from_fn(|i| {
                f64::midpoint(original[a as usize][i], original[b as usize][i])
            }))
        })
        .collect();
    (positions, midpoints)
}

fn coefficients(
    mesh: &TetComplex,
    positions: &[[f64; 3]],
    midpoints: &[[f64; 3]],
    field: impl Fn([f64; 3]) -> f64,
) -> Vec<f64> {
    let mut result: Vec<_> = positions.iter().copied().map(&field).collect();
    for (&[a, b], &point) in mesh.edges.iter().zip(midpoints) {
        // The original edge basis is lambda_a lambda_b, not 4 lambda_a lambda_b.
        result.push(4.0 * (field(point) - 0.5 * (result[a as usize] + result[b as usize])));
    }
    result
}

fn bilinear(matrix: &[f64], left: &[f64], right: &[f64]) -> f64 {
    let n = left.len();
    (0..n)
        .map(|i| left[i] * (0..n).map(|j| matrix[i * n + j] * right[j]).sum::<f64>())
        .sum()
}

fn relative(got: f64, expected: f64, tolerance: f64) {
    assert!(
        (got - expected).abs() <= tolerance * expected.abs().max(1e-30),
        "got {got:e}, expected {expected:e}"
    );
}

#[test]
fn consumed_affine_geometry_preserves_original_mixed_matrices_and_loads() {
    let size = [0.004, 0.001, 0.0012];
    let (mesh, positions) = box_mesh([1, 2, 2], size);
    let gas = air();
    let thermo = ThermoviscousThermodynamics::ideal_gas(gas).unwrap();
    let boundary = boundaries(&mesh, &positions, size, true);
    let original =
        ThermoviscousFe::new(&mesh, &positions, gas, &boundary, limits(), &mut || false).unwrap();
    let geometry = HelmholtzFe::new(&mesh, &positions, 2, acoustic(gas, thermo)).unwrap();
    let consumed =
        ThermoviscousFe::from_geometry(geometry, gas, &boundary, limits(), &mut || false).unwrap();
    assert!(!consumed.is_curved());
    assert_eq!(consumed.geometry_quadrature_order(), 4);
    assert_eq!(consumed.constant, original.constant);
    assert_eq!(consumed.mass, original.mass);
    assert_eq!(consumed.loads, original.loads);
    assert_eq!(consumed.areas, original.areas);
}

#[test]
fn mapped_physical_linear_fields_recover_exact_mixed_stress_storage_and_divergence() {
    run(|cx| {
        let original = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
        let a = 0.6;
        let (positions, midpoints) =
            mapped(&mesh, &original, |[x, y, z]| [x, y, z * (1.0 + a * x)]);
        let volume = 1.0 / 6.0 + a / 24.0;
        let mut gas = air();
        gas.bulk_viscosity_pa_s = 2.7e-5;
        let thermo = ThermoviscousThermodynamics::ideal_gas(gas).unwrap();
        let geometry = HelmholtzFe::new_quadratic(
            &mesh,
            &positions,
            &midpoints,
            2,
            7,
            acoustic(gas, thermo),
            QuadraticGeometryControls::default(),
            cx,
        )
        .unwrap();
        let walls = geometry.boundary_faces();
        let boundary = ThermoviscousBoundary {
            no_slip_faces: walls.clone(),
            isothermal_faces: vec![],
            ports: vec![],
        };
        // Inspect the full operator before its already-existing moving-wall
        // congruence. No constrained boundary field is faked as a free solution.
        let model = ThermoviscousFe::assemble_geometry(
            geometry,
            gas,
            thermo,
            false,
            &boundary,
            &walls.into_iter().collect(),
            limits(),
            &mut || false,
        )
        .unwrap();
        assert!(model.is_curved());
        assert_eq!(model.geometry_quadrature_order(), 7);
        let velocity = [
            coefficients(&mesh, &positions, &midpoints, |[x, y, _]| 2.0 * x - y),
            coefficients(&mesh, &positions, &midpoints, |[_, y, z]| 3.0 * y + z),
            coefficients(&mesh, &positions, &midpoints, |[x, _, z]| -x + 4.0 * z),
        ];
        let thermal = coefficients(&mesh, &positions, &midpoints, |[x, y, z]| {
            1.0 + 2.0 * x - 3.0 * y + z
        });
        // Observe an independently manufactured physical field through the
        // production curved-point factory, without solving a constrained state.
        let point_state = crate::acoustics::thermoviscous::time::ThermoviscousTimeState {
            pressure_pa: vec![2.0; mesh.vertex_count],
            velocity_m_s: (0..model.space.ndof)
                .map(|i| std::array::from_fn(|axis| velocity[axis][i]))
                .collect(),
            temperature_k: thermal.clone(),
            molecular_temperature_k: vec![],
            conductive_wall_temperature_k: vec![],
            wall_displacement_m_sqrt_kg: vec![],
            wall_velocity_m_sqrt_kg_per_s: vec![],
            wall_relaxation_coordinates_sqrt_j: vec![],
        };
        let bary = [0.1, 0.2, 0.3, 0.4];
        let sample = model
            .prepare_field_at(0, bary, &mut || false)
            .unwrap()
            .instantaneous(&point_state, &mut || false)
            .unwrap();
        let [x, y, z] = [bary[1], bary[2], bary[3] * (1.0 + a * bary[1])];
        for (actual, expected) in
            sample
                .velocity_m_s
                .into_iter()
                .zip([2.0 * x - y, 3.0 * y + z, -x + 4.0 * z])
        {
            relative(actual, expected, 3e-12);
        }
        relative(sample.temperature_k, 1.0 + 2.0 * x - 3.0 * y + z, 3e-12);
        for (j, expected) in [2.0, -3.0, 1.0].into_iter().enumerate() {
            relative(sample.temperature_gradient_k_m[j], expected, 3e-12);
            relative(
                sample.heat_flux_w_m2[j],
                -gas.thermal_conductivity_w_m_k * expected,
                3e-12,
            );
        }
        relative(sample.stress.bulk_pa, 9.0 * gas.bulk_viscosity_pa_s, 3e-12);
        sample.energy_flux().unwrap().total_w_m2().unwrap();
        let mut v = vec![0.0; model.n];
        let mut theta = vec![0.0; model.n];
        let mut constant_theta = vec![0.0; model.n];
        let mut p = vec![0.0; model.n];
        for (i, &index) in model.velocity_map.iter().enumerate() {
            for axis in 0..3 {
                v[axis * model.velocity_count + index.unwrap()] = velocity[axis][i];
            }
        }
        for (i, &index) in model.temperature_map.iter().enumerate() {
            theta[model.temperature_start + index.unwrap()] = thermal[i];
            if i < mesh.vertex_count {
                constant_theta[model.temperature_start + index.unwrap()] = 1.0;
            }
        }
        p[model.pressure_start..model.temperature_start].fill(1.0);
        // grad(v)=[[2,-1,0],[0,3,1],[-1,0,4]], div(v)=9,
        // |dev(sym grad(v))|^2=3.5; grad(theta)=(2,-3,1).
        relative(
            bilinear(&model.constant, &v, &v),
            volume * (7.0 * gas.dynamic_viscosity_pa_s + 81.0 * gas.bulk_viscosity_pa_s),
            3e-12,
        );
        relative(
            bilinear(&model.constant, &theta, &theta),
            14.0 * gas.thermal_conductivity_w_m_k * gas.temperature_k * volume,
            3e-12,
        );
        relative(bilinear(&model.constant, &p, &v), 9.0 * volume, 3e-12);
        relative(bilinear(&model.constant, &v, &p), -9.0 * volume, 3e-12);
        relative(
            bilinear(&model.mass, &p, &p),
            thermo.isothermal_compressibility_pa_inv * volume,
            3e-12,
        );
        relative(
            bilinear(&model.mass, &constant_theta, &constant_theta),
            gas.density_kg_m3 * gas.cp_j_kg_k * gas.temperature_k * volume,
            3e-12,
        );
        relative(
            bilinear(&model.mass, &p, &constant_theta),
            -thermo.isobaric_expansion_per_k * gas.temperature_k * volume,
            3e-12,
        );
    });
}

#[test]
fn curved_duct_keeps_work_reciprocity_and_independent_quadrature_refinement() {
    run(|cx| {
        let size = [0.004, 0.001, 0.0012];
        let (mesh, original) = box_mesh([1, 2, 2], size);
        let gas = air();
        let thermo = ThermoviscousThermodynamics::ideal_gas(gas).unwrap();
        // End faces remain affine; both displacement and area variation of the
        // intervening wall are represented by the same exact quadratic map.
        let (positions, midpoints) = mapped(&mesh, &original, |[x, y, z]| {
            [
                x,
                y,
                z * (1.0 + 0.4 * x / size[0]) + 0.8 * size[2] * x * (size[0] - x) / size[0].powi(2),
            ]
        });
        let boundary = boundaries(&mesh, &original, size, true);
        let drive = [C64::new(1.0, 0.2), C64::new(-0.3, 0.1)];
        let mut responses = Vec::new();
        for order in [4, 7] {
            let geometry = HelmholtzFe::new_quadratic(
                &mesh,
                &positions,
                &midpoints,
                2,
                order,
                acoustic(gas, thermo),
                QuadraticGeometryControls::default(),
                cx,
            )
            .unwrap();
            let model =
                ThermoviscousFe::from_geometry(geometry, gas, &boundary, limits(), &mut || false)
                    .unwrap();
            let cp = C64::new(970.0, -45.0);
            let field = model
                .prepare_with_heat_capacity(1700.0, cp, &mut || false)
                .unwrap()
                .response(&drive, &mut || false)
                .unwrap();
            assert_work(&field);
            assert!(
                field.power.viscous_w > 0.0
                    && field.power.conductive_w > 0.0
                    && field.power.caloric_w > 0.0
            );
            let y = model
                .admittance_with_heat_capacity(1700.0, cp, &mut || false)
                .unwrap();
            assert!((y[1] - y[2]).abs() < 2e-9 * y[1].abs());
            responses.push(field.inward_flow_m3_s);
        }
        for (a, b) in responses[0].iter().zip(&responses[1]) {
            assert!((*a - *b).abs() < 2e-5 * b.abs());
        }
    });
}

#[test]
fn curved_moving_wall_retains_full_vector_no_slip_and_reciprocal_work() {
    run(|cx| {
        let length = 0.002;
        let original = [
            [0.0, 0.0, 0.0],
            [length, 0.0, 0.0],
            [0.0, length, 0.0],
            [0.0, 0.0, length],
        ];
        let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
        let (positions, midpoints) = mapped(&mesh, &original, |[x, y, z]| {
            [x, y, z * (1.0 + 0.5 * x / length)]
        });
        let gas = air();
        let thermo = ThermoviscousThermodynamics::ideal_gas(gas).unwrap();
        let geometry = HelmholtzFe::new_quadratic(
            &mesh,
            &positions,
            &midpoints,
            2,
            7,
            acoustic(gas, thermo),
            QuadraticGeometryControls::default(),
            cx,
        )
        .unwrap();
        let faces = geometry.boundary_faces();
        let amplitude =
            1.0 / det::sqrt(gas.density_kg_m3 * length.powi(3) * (1.0 / 6.0 + 0.5 / 24.0));
        let shapes = [
            coefficients(&mesh, &positions, &midpoints, |[x, _, z]| {
                amplitude * (x + 0.5 * z) / length
            }),
            coefficients(&mesh, &positions, &midpoints, |[_, y, _]| {
                amplitude * y / length
            }),
            coefficients(&mesh, &positions, &midpoints, |[x, _, z]| {
                amplitude * (0.25 * x + z) / length
            }),
        ];
        let phi: Vec<_> = (0..shapes[0].len())
            .map(|i| std::array::from_fn(|c| shapes[c][i]))
            .collect();
        let source = ModalAcousticTimeModel::try_new_with_free_coordinates(
            48_000,
            vec![ModalAcousticMode {
                angular_frequency_rad_s: 300.0,
                damping_ratio: 0.02,
                pressure_per_modal_velocity: C64::ZERO,
            }],
            ModalAcousticTimeBudget::audible_reference(),
        )
        .unwrap();
        let boundary = ThermoviscousBoundary {
            no_slip_faces: faces.clone(),
            isothermal_faces: vec![],
            ports: vec![],
        };
        let model = ThermoviscousFe::from_geometry_with_moving_wall(
            geometry,
            gas,
            thermo,
            &boundary,
            MovingWall {
                faces,
                velocity_shapes: vec![phi],
                structure: source,
            },
            limits(),
            &mut || false,
        )
        .unwrap();
        let frequency = 1700.0;
        let force = C64::new(1.0, 0.3);
        let field = model
            .prepare(frequency, &mut || false)
            .unwrap()
            .response_with_wall_forces(&[], &[force], &mut || false)
            .unwrap();
        assert_work(&field);
        let wall = model.moving_wall().unwrap();
        let u = field.moving_wall.modal_velocity_m_sqrt_kg_per_s[0];
        for (actual, shape) in field.velocity_m_s.iter().zip(&wall.velocity_shapes[0]) {
            for c in 0..3 {
                assert!((actual[c] - u.scale(shape[c])).abs() < 1e-12 * u.abs() * amplitude);
            }
        }
        let solid_force = crate::acoustics::thermoviscous::moving_wall::modal_impedance(
            &wall.structure.modes()[0],
            wall.structure.damping_rate_per_s(0),
            std::f64::consts::TAU * frequency,
        ) * u;
        let fluid_force = field.moving_wall.fluid_force_n_per_sqrt_kg[0];
        assert!(
            (force + fluid_force - solid_force).abs() < 2e-8 * force.abs().max(solid_force.abs())
        );
        assert!(field.power.viscous_w > 0.0);
        relative(
            field.moving_wall.power_into_fluid_w,
            field.power.dissipated_w() - field.power.structural_damping_w,
            2e-8,
        );
    });
}

#[test]
fn geometry_admission_refuses_curved_pressure_traces_wrong_basis_and_wrong_fluid() {
    run(|cx| {
        let original = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
        let (positions, midpoints) =
            mapped(&mesh, &original, |[x, y, z]| [x, y, z * (1.0 + 0.5 * x)]);
        let gas = air();
        let thermo = ThermoviscousThermodynamics::ideal_gas(gas).unwrap();
        let boundary = ThermoviscousBoundary {
            no_slip_faces: vec![],
            isothermal_faces: vec![],
            ports: vec![BoundaryMode::uniform(&mesh.faces)],
        };
        let curved = HelmholtzFe::new_quadratic(
            &mesh,
            &positions,
            &midpoints,
            2,
            5,
            acoustic(gas, thermo),
            QuadraticGeometryControls::default(),
            cx,
        )
        .unwrap();
        assert!(
            matches!(ThermoviscousFe::from_geometry(curved, gas, &boundary, limits(), &mut || false),
            Err(ThermoviscousError::Acoustic(AcousticFeError::InvalidInput { what }))
            if what.contains("curved face"))
        );
        for (order, density) in [(1, gas.density_kg_m3), (2, 2.0 * gas.density_kg_m3)] {
            let geometry = HelmholtzFe::new(
                &mesh,
                &original,
                order,
                Medium {
                    density_kg_m3: density,
                    ..acoustic(gas, thermo)
                },
            )
            .unwrap();
            assert!(matches!(
                ThermoviscousFe::from_geometry(geometry, gas, &boundary, limits(), &mut || false),
                Err(ThermoviscousError::InvalidInput(_))
            ));
        }
        let geometry = HelmholtzFe::new(&mesh, &original, 2, acoustic(gas, thermo)).unwrap();
        assert!(matches!(
            ThermoviscousFe::from_geometry(geometry, gas, &boundary, limits(), &mut || true),
            Err(ThermoviscousError::Cancelled)
        ));
    });
}
