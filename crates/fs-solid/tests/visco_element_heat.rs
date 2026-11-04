//! Original spatial Kelvin/Maxwell loss reaches its physical element account.
//! All material values are synthetic numerical fixtures.

use fs_exec::Cx;
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_solid::linear3::hex::{HEX27_REFERENCE_NODES, HexConnectivity, HexQuadrature};
use fs_solid::linear3::visco::*;
use fs_solid::linear3::{
    TetAssemblyBudget, TetStaticSolveConfig, TetThermalStrainField, TetThermalStrainState,
};

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = fs_exec::CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 771,
                kernel_id: 771,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ));
    });
}
fn diagonal(values: [f64; 6]) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { values[i] } else { 0.0 }))
}
fn close(actual: f64, expected: f64) {
    let tolerance = 3e-11 * actual.abs().max(expected.abs()).max(1e-12);
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}
fn material(e: f64, eta: f64, branch: f64, identity: u8) -> TetViscoMaterial {
    let branches = if branch == 0.0 {
        Vec::new()
    } else {
        vec![PronyBranch::new(diagonal([branch; 6]), 0.15).unwrap()]
    };
    TetViscoMaterial {
        law: TensorProny::new_mandel(diagonal([e; 6]), branches, 0.05)
            .unwrap()
            .with_kelvin_viscosity_mandel(diagonal([eta; 6]))
            .unwrap(),
        density_kg_m3: 1000.0 + f64::from(identity),
        material_state_identity: [identity; 32],
    }
}
fn prepare(
    mesh: SolidMeshRef<'_>,
    materials: &[TetViscoMaterial],
    fixed: &[usize],
    cx: &Cx<'_>,
) -> ViscoElements {
    ViscoElements::new(
        mesh,
        TetViscoField::PerElement(materials),
        fixed,
        TetAssemblyBudget::standard(),
        ViscoStateBudget {
            max_internal_variables: 10_000,
        },
        cx,
    )
    .unwrap()
}
fn two_tets() -> ([[f64; 3]; 8], [[usize; 4]; 2]) {
    // Unequal physical volumes, unrelated to the very different loss strengths.
    (
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [3.0, 0.0, 0.0],
            [5.0, 0.0, 0.0],
            [3.0, 1.0, 0.0],
            [3.0, 0.0, 1.0],
        ],
        [[0, 1, 2, 3], [4, 5, 6, 7]],
    )
}
fn check_sums(step: &ViscoForceStep) {
    close(
        step.element_steps.iter().map(|e| e.energy.stored_j).sum(),
        step.energy.stored_j,
    );
    close(
        step.element_steps.iter().map(|e| e.energy.work_j).sum(),
        step.energy.work_j,
    );
    close(
        step.element_steps
            .iter()
            .map(|e| e.energy.dissipated_j)
            .sum(),
        step.energy.dissipated_j,
    );
    close(
        step.element_steps.iter().map(|e| e.step_work_j).sum(),
        step.step_work_j,
    );
    close(
        step.element_steps
            .iter()
            .map(|e| e.step_mechanical_work_j)
            .sum(),
        step.step_mechanical_work_j,
    );
    close(
        step.element_steps
            .iter()
            .map(|e| e.step_eigenstrain_work_j)
            .sum(),
        step.step_eigenstrain_work_j,
    );
    close(
        step.element_steps.iter().map(|e| e.step_dissipated_j).sum(),
        step.step_dissipated_j,
    );
    for element in &step.element_steps {
        assert!(element.step_dissipated_j >= 0.0);
        let residual =
            element.step_work_j - element.stored_energy_change_j - element.step_dissipated_j;
        let scale = element.step_work_j.abs()
            + element.stored_energy_change_j.abs()
            + element.step_dissipated_j;
        assert!(
            residual.abs() <= 3e-11 * scale.max(1e-12),
            "element energy residual={residual:e}"
        );
        close(
            element.step_work_j,
            element.step_mechanical_work_j + element.step_eigenstrain_work_j,
        );
    }
}

#[test]
fn heterogeneous_heat_uses_actual_material_ramps_and_retains_memory_on_refusal() {
    with_cx(|cx| {
        let (nodes, cells) = two_tets();
        let materials = [
            material(1200.0, 0.0, 0.0, 1),
            material(800.0, 3.0, 400.0, 2),
        ];
        let solid = prepare(
            SolidMeshRef::Tet4 {
                nodes_m: &nodes,
                elements: &cells,
            },
            &materials,
            &[0, 1, 2],
            cx,
        );
        assert_eq!(solid.element_count(), 2);
        assert_eq!(solid.element_reference_mass_kg(2), None);
        let volumes = [1.0 / 6.0, 1.0 / 3.0];
        for e in 0..2 {
            close(solid.element_reference_volume_m3(e).unwrap(), volumes[e]);
            close(
                solid.element_reference_mass_kg(e).unwrap(),
                volumes[e] * materials[e].density_kg_m3,
            );
        }
        close(
            (0..2)
                .map(|e| solid.element_reference_mass_kg(e).unwrap())
                .sum(),
            solid.equilibrium().total_mass_kg,
        );
        let mut state = solid.state(cx).unwrap();
        assert!(core::ptr::eq(state.system(), &solid));
        let mut original = [materials[0].law.state(), materials[1].law.state()];
        for (strain, dt) in [
            ([0.003, -0.002], 0.04),
            ([0.003, -0.002], 0.3),
            ([-0.001, 0.0005], 0.09),
        ] {
            let displacement: Vec<_> = nodes
                .iter()
                .enumerate()
                .flat_map(|(n, p)| {
                    let e = n / 4;
                    [
                        strain[e] * (p[0] - if e == 0 { 0.0 } else { 3.0 }),
                        0.0,
                        0.0,
                    ]
                })
                .collect();
            let step = state.step_internal_force(&displacement, dt, cx).unwrap();
            assert_eq!(step.element_steps.len(), 2);
            check_sums(&step);
            assert_eq!(step.element_steps[0].step_dissipated_j, 0.0);
            assert!(step.element_steps[1].step_dissipated_j > 0.0);
            for e in 0..2 {
                let expected = materials[e]
                    .law
                    .step_with_mean(&mut original[e], [strain[e], 0.0, 0.0, 0.0, 0.0, 0.0], dt)
                    .unwrap();
                close(
                    step.element_steps[e].step_dissipated_j,
                    volumes[e] * expected.step_dissipated,
                );
                close(
                    step.element_steps[e].energy.stored_j,
                    volumes[e] * original[e].stored,
                );
                close(
                    step.element_steps[e].step_work_j,
                    volumes[e] * expected.step_work,
                );
            }
        }
        // A later element refuses after the earlier one has already been evaluated.
        let mut control = state.clone();
        let mut invalid = state.displacement_m().to_vec();
        invalid[15] = 1e4;
        assert!(state.step_internal_force(&invalid, 0.01, cx).is_err());
        let held = state.displacement_m().to_vec();
        let actual = state.step_internal_force(&held, 0.02, cx).unwrap();
        let expected = control.step_internal_force(&held, 0.02, cx).unwrap();
        assert_eq!(actual.element_steps, expected.element_steps);
        assert_eq!(actual.energy, expected.energy);
    });
}

#[test]
fn restrained_expansion_keeps_element_eigenstrain_work_distinct_from_heat() {
    with_cx(|cx| {
        let (nodes, cells) = two_tets();
        let materials = [
            material(1200.0, 2.0, 0.0, 6),
            material(800.0, 3.0, 400.0, 7),
        ];
        let solid = prepare(
            SolidMeshRef::Tet4 {
                nodes_m: &nodes,
                elements: &cells,
            },
            &materials,
            &[],
            cx,
        );
        let free_strain = [0.001, 0.002];
        let thermal: Vec<_> = (0..2)
            .map(|e| {
                TetThermalStrainState::try_new(
                    320.0,
                    300.0,
                    [free_strain[e], 0.0, 0.0, 0.0, 0.0, 0.0],
                    materials[e].material_state_identity,
                    fs_blake3::ContentHash([71; 32]),
                    0.05,
                )
                .unwrap()
            })
            .collect();
        let mut state = solid.state(cx).unwrap();
        let step = state
            .step_internal_force_with_thermal_strain(
                &[0.0; 24],
                0.04,
                TetThermalStrainField::PerElement(&thermal),
                cx,
            )
            .unwrap();
        check_sums(&step);
        for element in &step.element_steps {
            assert_eq!(element.step_mechanical_work_j, 0.0);
            close(element.step_work_j, element.step_eigenstrain_work_j);
            assert!(element.step_eigenstrain_work_j > element.step_dissipated_j);
        }
        // No new expansion or nodal work; only the Maxwell element still heats.
        let hold = state.step_internal_force(&[0.0; 24], 0.2, cx).unwrap();
        check_sums(&hold);
        assert_eq!(hold.element_steps[0].step_dissipated_j, 0.0);
        assert!(hold.element_steps[1].step_dissipated_j > 0.0);
        for element in &hold.element_steps {
            assert_eq!(element.step_eigenstrain_work_j, 0.0);
            assert_eq!(element.step_mechanical_work_j, 0.0);
        }
    });
}

#[test]
fn curved_hex27_heat_integrates_every_original_quadrature_point() {
    with_cx(|cx| {
        // Exact Q2 geometry, with rational physical strains under u_x=b*xi².
        let (length, warp, amplitude, eta, dt) = (0.1, 0.2_f64, 1e-5, 4.0, 0.03);
        let nodes = HEX27_REFERENCE_NODES.map(|p| {
            let [x, y, z] = p.map(f64::from);
            [length * x * (1.0 + warp * y * y), length * y, length * z]
        });
        let cells = [core::array::from_fn(|i| i)];
        let displacement: Vec<_> = HEX27_REFERENCE_NODES
            .iter()
            .flat_map(|p| [amplitude * f64::from(p[0]).powi(2), 0.0, 0.0])
            .collect();
        let materials = [material(1200.0, eta, 0.0, 3)];
        let i0 = 2.0 * warp.sqrt().atan() / warp.sqrt();
        let i2 = (2.0 - i0) / warp;
        let axial = (16.0 / 3.0) * amplitude * amplitude * length * i0;
        let shear = (32.0 / 5.0) * amplitude * amplitude * length * warp * warp * i2;
        let expected_heat = eta * (axial + shear) / dt;
        let mut previous_error = f64::INFINITY;
        for order in [3, 5, 7] {
            let solid = prepare(
                SolidMeshRef::Hex {
                    nodes_m: &nodes,
                    elements: HexConnectivity::Hex27(&cells),
                    quadrature: HexQuadrature {
                        points_per_dir: order,
                    },
                },
                &materials,
                &[],
                cx,
            );
            let step = solid
                .state(cx)
                .unwrap()
                .step_internal_force(&displacement, dt, cx)
                .unwrap();
            assert_eq!(step.element_steps.len(), 1);
            check_sums(&step);
            let error = (step.element_steps[0].step_dissipated_j / expected_heat - 1.0).abs();
            assert!(
                error < previous_error,
                "quadrature={order}, error={error:e}"
            );
            if order == 7 {
                assert!(error < 1e-8, "{error:e}");
            }
            previous_error = error;
            close(
                solid.element_reference_volume_m3(0).unwrap(),
                8.0 * length.powi(3) * (1.0 + warp / 3.0),
            );
        }
    });
}

#[test]
fn sparse_inertial_ringdown_returns_original_element_heat_without_modal_projection() {
    with_cx(|cx| {
        let (nodes, cells) = two_tets();
        let materials = [
            material(1200.0, 0.0, 0.0, 4),
            material(800.0, 3.0, 400.0, 5),
        ];
        let fixed: Vec<_> = (0..24).filter(|d| ![3, 15].contains(d)).collect();
        let solid = prepare(
            SolidMeshRef::Tet4 {
                nodes_m: &nodes,
                elements: &cells,
            },
            &materials,
            &fixed,
            cx,
        );
        let prepared = solid
            .dynamic_operator(
                0.005,
                ViscoDynamicConfig {
                    maximum_steps: 100,
                    maximum_free_dofs: 24,
                    solve: TetStaticSolveConfig {
                        relative_residual_tolerance: 1e-13,
                        maximum_iterations: 40,
                    },
                    absolute_energy_tolerance_j: 1e-13,
                    relative_energy_tolerance: 1e-11,
                },
                cx,
            )
            .unwrap();
        assert!(core::ptr::eq(prepared.system(), &solid));
        let mut velocity = [0.0; 24];
        velocity[3] = 0.001;
        velocity[15] = -0.002;
        let mut state = prepared
            .start(solid.state(cx).unwrap(), &velocity, cx)
            .unwrap();
        assert!(core::ptr::eq(state.operator(), &prepared));
        assert_eq!(state.operator().duration_s(), 0.005);
        let initial = state.kinetic_energy_j();
        let mut heat = [0.0; 2];
        for _ in 0..80 {
            let step = state.advance(&[0.0; 24], cx).unwrap();
            check_sums(&step.constitutive);
            for (total, element) in heat.iter_mut().zip(&step.constitutive.element_steps) {
                *total += element.step_dissipated_j;
            }
            close(
                step.kinetic_energy_j
                    + step.constitutive.energy.stored_j
                    + heat.iter().sum::<f64>(),
                initial,
            );
        }
        assert_eq!(heat[0], 0.0);
        assert!(heat[1] > 0.0);
    });
}
