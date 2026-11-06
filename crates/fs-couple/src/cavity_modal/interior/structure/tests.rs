use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_feec::acoustics::helmholtz::{HelmholtzFe, Medium};
use fs_feec::fixtures::{masked_cube_grid, single_tet};
use fs_solid::linear3::{
    TetAssemblyBudget, TetElasticMaterial, TetLinearElasticProblem, TetMaterialField,
};
use fs_sparse::Coo;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 81,
                kernel_id: 19,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn controls() -> HarmonicResponseControls {
    HarmonicResponseControls {
        maximum_coordinates: 128,
        maximum_payload_bytes: 16 * 1024 * 1024,
        maximum_relative_residual: 1e-10,
    }
}

fn dense(values: &[f64], n: usize) -> Csr {
    let mut coo = Coo::new(n, n);
    for i in 0..n {
        for j in 0..n {
            if values[i * n + j] != 0.0 {
                coo.push(i, j, values[i * n + j]);
            }
        }
    }
    coo.assemble()
}

fn uniform<'f, 'c>(
    geometry: &'f HelmholtzFe<'c>,
    medium: Medium,
) -> HeterogeneousHelmholtzFe<'f, 'c> {
    HeterogeneousHelmholtzFe::new(
        geometry,
        geometry.order() + 2,
        &mut |_, _| Ok::<_, AcousticFeError>(medium),
        &mut || false,
    )
    .unwrap()
}

fn product(matrix: &Csr, vector: &[C64], row: usize) -> C64 {
    let (columns, values) = matrix.row(row);
    columns
        .iter()
        .zip(values)
        .fold(C64::ZERO, |sum, (&column, &value)| {
            sum + vector[column].scale(value)
        })
}

#[test]
fn original_nonuniform_pressure_and_full_structural_pencil_balance_work() {
    let (complex, points) = single_tet();
    let medium = Medium {
        density_kg_m3: 1.2,
        sound_speed_m_s: 3.0,
    };
    let geometry = HelmholtzFe::new(&complex, &points, 2, medium).unwrap();
    let interior = HeterogeneousHelmholtzFe::new(
        &geometry,
        5,
        &mut |_, p: [f64; 3]| {
            Ok::<_, AcousticFeError>(Medium {
                density_kg_m3: 1.1 + 0.2 * p[0],
                sound_speed_m_s: 3.0 + 0.4 * p[1],
            })
        },
        &mut || false,
    )
    .unwrap();
    let mass = dense(&[2.0, 0.4, 0.4, 1.1], 2);
    let stiffness = dense(&[70.0, -13.0, -13.0, 40.0], 2);
    // Symmetric loss plus a genuine structural gyroscopic part.
    let damping = dense(&[0.7, 0.3, -0.1, 0.5], 2);
    let traces = [
        vec![([0, 1, 2], [0.7, -0.3, 0.9])],
        vec![([2, 0, 1], [1.1, -0.2, 0.8])],
    ];
    let wall = ResistiveBoundary {
        faces: vec![[0, 1, 3]],
        admittance_m_pa_s: 0.08,
    };
    let surface = geometry.boundary_mass(&wall.faces).unwrap();
    with_cx(|cx| {
        let model = StructuralInterior::new(
            &interior,
            StructuralMatrices {
                mass: &mass,
                stiffness: &stiffness,
                damping: Some(&damping),
            },
            &traces,
            &[wall],
            controls(),
            cx,
        )
        .unwrap();
        let force = [C64::new(0.7, -0.2), C64::new(-0.3, 0.4)];
        let flow = model
            .project_inward_velocity(&[(
                [1, 2, 3],
                [
                    C64::new(0.1, 0.03),
                    C64::new(-0.04, 0.07),
                    C64::new(0.02, -0.01),
                ],
            )])
            .unwrap();
        let frequency = 0.7;
        let omega = core::f64::consts::TAU * frequency;
        let prepared = model.prepare(frequency, cx).unwrap();
        let response = prepared.solve(&force, &flow, cx).unwrap();
        assert_eq!(response.pressure_pa.len(), geometry.ndof());
        assert_eq!(response.structural_displacement.len(), 2);
        for row in 0..2 {
            let residual = product(&stiffness, &response.structural_displacement, row)
                - product(&mass, &response.structural_displacement, row).scale(omega * omega)
                + product(&damping, &response.structural_velocity, row)
                - response.fluid_force[row]
                - force[row];
            assert!(residual.abs() < 1e-10, "structural row {row}: {residual:?}");
        }
        for row in 0..geometry.ndof() {
            let residual = product(interior.stiffness(), &response.pressure_potential_pa_s, row)
                - product(interior.mass(), &response.pressure_potential_pa_s, row)
                    .scale(omega * omega)
                + product(&surface, &response.pressure_pa, row).scale(0.08)
                + response.displaced_flow_m3_s[row]
                - flow[row];
            assert!(residual.abs() < 1e-10, "pressure row {row}: {residual:?}");
            let pressure = C64::new(
                omega * response.pressure_potential_pa_s[row].im,
                -omega * response.pressure_potential_pa_s[row].re,
            );
            assert!((pressure - response.pressure_pa[row]).abs() < 1e-13);
        }
        let boundary_work = work(&response.pressure_pa, &response.displaced_flow_m3_s);
        assert!((boundary_work - response.average_fluid_to_structure_power_w).abs() < 1e-12);
        let independently_integrated_wall_power = 0.5
            * (0..geometry.ndof())
                .map(|row| {
                    let y_p = product(&surface, &response.pressure_pa, row).scale(0.08);
                    let p = response.pressure_pa[row];
                    p.re * y_p.re + p.im * y_p.im
                })
                .sum::<f64>();
        assert!(independently_integrated_wall_power > 0.0);
        assert!(
            (independently_integrated_wall_power - response.average_wall_dissipation_w).abs()
                < 1e-12
        );
        let input =
            response.average_structural_input_power_w + response.average_acoustic_input_power_w;
        let losses =
            response.average_structural_dissipation_w + response.average_wall_dissipation_w;
        assert!((input - losses).abs() < 2e-11);
        assert!(response.power_balance_residual_w.abs() < 2e-11);
        // Reusing one whole-system factor preserves complex source superposition.
        let doubled_force = force.map(|z| z.scale(2.0));
        let doubled_flow: Vec<_> = flow.iter().map(|z| z.scale(2.0)).collect();
        let doubled = prepared.solve(&doubled_force, &doubled_flow, cx).unwrap();
        for (a, b) in doubled.pressure_pa.iter().zip(&response.pressure_pa) {
            assert!((*a - b.scale(2.0)).abs() < 1e-11);
        }
    });
}

#[test]
fn full_pressure_field_recovers_closed_cavity_spring_and_added_air_mass() {
    let (complex, points) = masked_cube_grid(1, 1, 1, &|_, _, _| true);
    let medium = Medium {
        density_kg_m3: 1.3,
        sound_speed_m_s: 3.0,
    };
    let geometry = HelmholtzFe::new(&complex, &points, 2, medium).unwrap();
    let interior = uniform(&geometry, medium);
    let piston: Vec<_> = geometry
        .boundary_faces()
        .into_iter()
        .filter(|face| face.iter().all(|&vertex| points[vertex as usize][0] == 1.0))
        .map(|face| (face, [1.0; 3]))
        .collect();
    assert_eq!(piston.len(), 2);
    let mass = dense(&[0.7], 1);
    let stiffness = dense(&[30.0], 1);
    let damping = dense(&[0.4], 1);
    with_cx(|cx| {
        let matrices = || StructuralMatrices {
            mass: &mass,
            stiffness: &stiffness,
            damping: Some(&damping),
        };
        let model =
            StructuralInterior::new(&interior, matrices(), &[piston], &[], controls(), cx).unwrap();
        let omega = 0.05;
        let force = C64::new(1.0, 0.2);
        let zero = vec![C64::ZERO; geometry.ndof()];
        let response = model
            .prepare(omega / core::f64::consts::TAU, cx)
            .unwrap()
            .solve(&[force], &zero, cx)
            .unwrap();
        let reaction = response.fluid_force[0] / response.structural_displacement[0];
        // A=L=V=1. Exact 1D closed-cavity reaction is
        // -rho*c*omega*cot(omega*L/c)*A, whose expansion retains BOTH
        // rho*c^2*A/L stiffness and rho*A*L/3 added inertial mass.
        let bulk = medium.density_kg_m3 * medium.sound_speed_m_s * medium.sound_speed_m_s;
        let added_mass = (reaction.re + bulk) / (omega * omega);
        let expected_mass = medium.density_kg_m3 / 3.0;
        assert!(
            (added_mass / expected_mass - 1.0).abs() < 3e-4,
            "{added_mass}"
        );
        assert!(reaction.im.abs() < 1e-9);
        let phase = omega / medium.sound_speed_m_s;
        let exact_reaction =
            -medium.density_kg_m3 * medium.sound_speed_m_s * omega * fs_math::det::cos(phase)
                / fs_math::det::sin(phase);
        let expected = force / C64::new(30.0 - omega * omega * 0.7 - exact_reaction, -omega * 0.4);
        assert!((response.structural_displacement[0] - expected).abs() / expected.abs() < 2e-8);
        assert!(geometry.ndof() > points.len());
        assert_eq!(response.pressure_pa.len(), geometry.ndof());

        // Zero physical normal motion recovers the original dry structure while
        // retaining, rather than deleting or regularizing, the full fluid space.
        let uncoupled =
            StructuralInterior::new(&interior, matrices(), &[Vec::new()], &[], controls(), cx)
                .unwrap();
        let dry = uncoupled
            .prepare(omega / core::f64::consts::TAU, cx)
            .unwrap()
            .solve(&[force], &zero, cx)
            .unwrap();
        let expected_dry = force / C64::new(30.0 - omega * omega * 0.7, -omega * 0.4);
        assert!((dry.structural_displacement[0] - expected_dry).abs() < 1e-12);
        assert!(dry.pressure_pa.iter().all(|p| p.abs() < 1e-13));
    });
}

#[test]
fn actual_solid_coordinates_give_consistent_face_work_independent_of_winding() {
    let (complex, fluid_points) = single_tet();
    let medium = Medium {
        density_kg_m3: 1.2,
        sound_speed_m_s: 3.0,
    };
    let geometry = HelmholtzFe::new(&complex, &fluid_points, 1, medium).unwrap();
    let interior = uniform(&geometry, medium);
    let solid_points = [
        fluid_points[0],
        fluid_points[1],
        fluid_points[2],
        [0.0, 0.0, -1.0],
    ];
    let tetrahedra = [[0, 1, 2, 3]];
    let material = TetElasticMaterial::try_new(6.0, 120.0, 0.25, [7; 32]).unwrap();
    let fixed = [0, 1, 3, 4, 6, 7, 9, 10, 11];
    let problem = TetLinearElasticProblem {
        nodes_m: &solid_points,
        tetrahedra: &tetrahedra,
        materials: TetMaterialField::Uniform(&material),
        fixed_dofs: &fixed,
        budget: TetAssemblyBudget::standard(),
    };
    with_cx(|cx| {
        let assembly = problem.assemble(cx).unwrap();
        assert_eq!(assembly.free_dofs, [2, 5, 8]);
        let face = [ConformingTetFace {
            fluid: [0, 1, 2],
            solid: [0, 1, 2],
        }];
        let model = StructuralInterior::from_tet_elastic(
            &interior,
            ConformingTetStructure {
                assembly: &assembly,
                nodes_m: &solid_points,
                tetrahedra: &tetrahedra,
                interface: &face,
            },
            None,
            &[],
            controls(),
            cx,
        )
        .unwrap();
        let n = model.structural_dofs() + model.pressure_dofs();
        assert_eq!(model.structural_dofs(), 3);
        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(model.mass[i * n + j], assembly.mass.get(i, j));
                assert_eq!(model.stiffness[i * n + j], assembly.stiffness.get(i, j));
                // Outward fluid normal is -z. Integral(N_i N_j)dA has the
                // consistent triangle values A/6 on diagonal and A/12 off it.
                let expected = -(if i == j { 2.0 } else { 1.0 }) / 24.0;
                assert!((model.coupling()[i * 3 + j] - expected).abs() < 1e-14);
            }
            assert!(model.coupling()[3 * 3 + i].abs() < 1e-14);
        }
        let reversed = [ConformingTetFace {
            fluid: [2, 1, 0],
            solid: [2, 1, 0],
        }];
        let other = StructuralInterior::from_tet_elastic(
            &interior,
            ConformingTetStructure {
                assembly: &assembly,
                nodes_m: &solid_points,
                tetrahedra: &tetrahedra,
                interface: &reversed,
            },
            None,
            &[],
            controls(),
            cx,
        )
        .unwrap();
        for (a, b) in model.coupling().iter().zip(other.coupling()) {
            assert!((a - b).abs() < 1e-14);
        }
        let q = [
            C64::new(0.03, -0.02),
            C64::new(-0.01, 0.04),
            C64::new(0.02, 0.01),
        ];
        let response = model
            .prepare(0.31, cx)
            .unwrap()
            .solve(&q, &vec![C64::ZERO; geometry.ndof()], cx)
            .unwrap();
        assert!(response.relative_residual < 1e-10);
        assert_eq!(
            response.structural_displacement.len(),
            assembly.free_dofs.len()
        );
        assert!(
            (work(&response.pressure_pa, &response.displaced_flow_m3_s)
                - response.average_fluid_to_structure_power_w)
                .abs()
                < 1e-12
        );

        // Matching vertex labels cannot admit a body lying inside the fluid.
        let mut wrong_side = solid_points;
        wrong_side[3][2] = 1.0;
        let rejected = StructuralInterior::from_tet_elastic(
            &interior,
            ConformingTetStructure {
                assembly: &assembly,
                nodes_m: &wrong_side,
                tetrahedra: &tetrahedra,
                interface: &face,
            },
            None,
            &[],
            controls(),
            cx,
        );
        assert!(matches!(rejected, Err(StructuralInteriorError::Invalid(_))));
    });
}

#[test]
fn coupled_factor_handles_a_bare_acoustic_pole_without_inverting_the_fluid() {
    let (complex, points) = single_tet();
    let medium = Medium {
        density_kg_m3: 1.2,
        sound_speed_m_s: 2.0,
    };
    let geometry = HelmholtzFe::new(&complex, &points, 1, medium).unwrap();
    let interior = uniform(&geometry, medium);
    let omega = fs_math::det::sqrt(80.0) * medium.sound_speed_m_s;
    // The reference P1 tetrahedron has this exact nonconstant generalized
    // eigenpair. Its face-volume coupling is nonzero.
    let eigenvector = [-3.0, 1.0, 1.0, 1.0].map(|value| C64::new(value, 0.0));
    for row in 0..4 {
        let bare_residual = product(interior.stiffness(), &eigenvector, row)
            - product(interior.mass(), &eigenvector, row).scale(omega * omega);
        assert!(bare_residual.abs() < 1e-12);
    }
    let mass = dense(&[1.0], 1);
    let stiffness = dense(&[10.0], 1);
    let damping = dense(&[0.8], 1);
    with_cx(|cx| {
        let model = StructuralInterior::new(
            &interior,
            StructuralMatrices {
                mass: &mass,
                stiffness: &stiffness,
                damping: Some(&damping),
            },
            &[vec![([0, 1, 2], [1.0; 3])]],
            &[],
            controls(),
            cx,
        )
        .unwrap();
        let response = model
            .prepare(omega / core::f64::consts::TAU, cx)
            .unwrap()
            .solve(&[C64::new(1.0, 0.0)], &[C64::ZERO; 4], cx)
            .unwrap();
        assert!(response.relative_residual < 1e-10);
        assert!(response.structural_displacement[0].abs() < 1e-10);
        assert!((response.fluid_force[0] + C64::new(1.0, 0.0)).abs() < 1e-8);
        assert!(response.pressure_pa.iter().any(|p| p.abs() > 1.0));
    });
}
