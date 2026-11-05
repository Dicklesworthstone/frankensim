//! G1/G3: anisotropic pore drag in native momentum, continuity and physical work.

use fs_flux::{
    BrinkmanCell, FluxParams, FluxSystem, MaterialFluxNumerics, MaterialFluxSystem,
    MaterialStokesBoundary, StokesBoundary, StokesPressureReference, TriMesh,
};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::porous::permeability::IntrinsicPermeability;
use fs_solid::Mesh2;
use std::cell::Cell;

fn tensor(angle: f64) -> IntrinsicPermeability {
    let (s, c) = angle.sin_cos();
    IntrinsicPermeability::orthotropic(
        [0.002, 0.0005, 0.003],
        [[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]],
    )
    .unwrap()
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} vs {expected:e}, tolerance {tolerance:e}"
    );
}
fn no_slip(_: [f64; 2]) -> [f64; 2] {
    [0.0; 2]
}

#[test]
fn g1_pressure_driven_anisotropic_darcy_slab_needs_no_artificial_velocity_pin() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 3));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        nu: 0.02,
        ..FluxParams::default()
    };
    let permeability = tensor(0.53);
    let cells = vec![BrinkmanCell::try_new(permeability, 1.0).unwrap(); mesh.tris.len()];
    // A prescribed world pressure gradient drives transverse flow in a
    // rotated medium. Darcy gives u = -K grad(p)/nu independently of FEM.
    let gradient = [-3.0, 1.0];
    let (s, c) = 0.53_f64.sin_cos();
    let kxx = 0.002 * c * c + 0.0005 * s * s;
    let kxy = (0.002 - 0.0005) * s * c;
    let kyy = 0.002 * s * s + 0.0005 * c * c;
    let exact = [
        -(kxx * gradient[0] + kxy * gradient[1]) / params.nu,
        -(kxy * gradient[0] + kyy * gradient[1]) / params.nu,
    ];
    let traction = |x: [f64; 2], normal: [f64; 2]| {
        let pressure = 8.0 + gradient[0] * x[0] + gradient[1] * x[1];
        [-pressure * normal[0], -pressure * normal[1]]
    };
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| StokesBoundary::TotalTraction {
            edge,
            value: &traction,
        })
        .collect();
    let assembly = system
        .assemble_brinkman_mixed(params, &cells, &|_| [0.0; 2], &boundaries)
        .unwrap();
    let flow = assembly.solve().unwrap();
    let (error, divergence) = system.velocity_error(&flow.field.x, &|_| exact);
    assert!(
        error < 2e-10 && divergence < 2e-10,
        "error={error:e}, div={divergence:e}"
    );
    for (pressure, x) in flow.pressure.iter().zip(&mesh.centroids) {
        near(
            *pressure,
            8.0 + gradient[0] * x[0] + gradient[1] * x[1],
            2e-9,
        );
    }
    let power: f64 = assembly
        .cell_darcy_power_per_density(&flow.field.x[..system.n_u])
        .unwrap()
        .iter()
        .sum();
    near(
        power,
        -gradient[0] * exact[0] - gradient[1] * exact[1],
        2e-9,
    );
    let free = vec![BrinkmanCell::free_fluid(); mesh.tris.len()];
    assert!(
        system
            .assemble_brinkman_mixed(params, &free, &|_| [0.0; 2], &boundaries)
            .is_err()
    );
}

#[test]
fn g1_rotated_anisotropic_affine_flow_has_exact_pressure_flux_and_dissipation() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 3));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        nu: 0.02,
        ..FluxParams::default()
    };
    let ratio = 1.7;
    let velocity = |x: [f64; 2]| [0.3 + 0.2 * x[1], -0.1];
    for angle in [0.0, 0.53, std::f64::consts::FRAC_PI_2] {
        let cells = vec![BrinkmanCell::try_new(tensor(angle), ratio).unwrap(); mesh.tris.len()];
        // Independently rotate diag(1/k1,1/k2), retaining both cross terms.
        let (s, c) = angle.sin_cos();
        let rxx = c * c / 0.002 + s * s / 0.0005;
        let rxy = s * c * (1.0 / 0.002 - 1.0 / 0.0005);
        let ryy = s * s / 0.002 + c * c / 0.0005;
        let force = |x| {
            let u = velocity(x);
            [
                params.nu * (rxx * u[0] + rxy * u[1]),
                params.nu * (rxy * u[0] + ryy * u[1]),
            ]
        };
        let traction = |_: [f64; 2], n: [f64; 2]| {
            let shear = params.nu * ratio * 0.2;
            [-3.0 * n[0] + shear * n[1], shear * n[0] - 3.0 * n[1]]
        };
        let boundaries: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(e, _)| mesh.is_boundary(*e))
            .map(|(edge, geometry)| {
                if geometry.normal[0] < -0.5 {
                    StokesBoundary::Velocity {
                        edge,
                        value: &velocity,
                    }
                } else {
                    StokesBoundary::TotalTraction {
                        edge,
                        value: &traction,
                    }
                }
            })
            .collect();
        let assembly = system
            .assemble_brinkman_mixed(params, &cells, &force, &boundaries)
            .unwrap();
        let flow = assembly.solve().unwrap();
        let (error, divergence) = system.velocity_error(&flow.field.x, &velocity);
        assert!(
            error < 2e-10 && divergence < 2e-10,
            "angle={angle}, error={error:e}, divergence={divergence:e}"
        );
        assert_eq!(
            flow.pressure_reference,
            StokesPressureReference::TractionDatum
        );
        for pressure in flow.pressure {
            near(pressure, 3.0, 2e-9);
        }
        let powers = assembly
            .cell_darcy_power_per_density(&flow.field.x[..system.n_u])
            .unwrap();
        assert!(powers.iter().all(|power| *power >= 0.0));
        let expected = params.nu * (rxx * (0.09 + 0.06 + 0.04 / 3.0) - 0.08 * rxy + 0.01 * ryy);
        near(powers.iter().sum(), expected, 2e-9);
    }
}

#[test]
fn g1_brinkman_channel_refines_to_independent_hyperbolic_profile() {
    let (nu, permeability, ratio) = (0.7_f64, 0.05_f64, 1.5_f64);
    let length = (ratio * permeability).sqrt();
    let amplitude = permeability / nu;
    let velocity = |x: [f64; 2]| {
        [
            amplitude * (1.0 - ((x[1] - 0.5) / length).cosh() / (0.5 / length).cosh()),
            0.0,
        ]
    };
    let traction = |x: [f64; 2], n: [f64; 2]| {
        let shear = -nu * ratio * amplitude / length * ((x[1] - 0.5) / length).sinh()
            / (0.5 / length).cosh();
        [shear * n[1], shear * n[0]]
    };
    let mut errors = Vec::new();
    for n in [3, 6] {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, n, n));
        let system = FluxSystem::new(&mesh);
        let cells = vec![
            BrinkmanCell::try_new(
                IntrinsicPermeability::isotropic(permeability).unwrap(),
                ratio
            )
            .unwrap();
            mesh.tris.len()
        ];
        let boundaries: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(e, _)| mesh.is_boundary(*e))
            .map(|(edge, geometry)| {
                if geometry.normal[0].abs() > 0.5 {
                    StokesBoundary::TotalTraction {
                        edge,
                        value: &traction,
                    }
                } else {
                    StokesBoundary::Velocity {
                        edge,
                        value: &no_slip,
                    }
                }
            })
            .collect();
        let flow = system
            .assemble_brinkman_mixed(
                FluxParams {
                    nu,
                    ..FluxParams::default()
                },
                &cells,
                &|_| [1.0, 0.0],
                &boundaries,
            )
            .unwrap()
            .solve()
            .unwrap();
        let (error, divergence) = system.velocity_error(&flow.field.x, &velocity);
        assert!(divergence < 2e-10);
        errors.push(error);
    }
    assert!(errors[0] / errors[1] > 2.5, "channel refinement {errors:?}");
    assert!(errors[1] < 0.0015, "channel L2 {errors:?}");
}

#[test]
fn g0_reaction_is_symmetric_positive_and_matches_independent_power() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        nu: 0.3,
        ..FluxParams::default()
    };
    let porous = vec![BrinkmanCell::try_new(tensor(0.53), 1.0).unwrap(); mesh.tris.len()];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(e, _)| mesh.is_boundary(*e))
        .map(|(edge, _)| StokesBoundary::Velocity {
            edge,
            value: &no_slip,
        })
        .collect();
    let base = system
        .assemble_stokes_mixed(params, &|_| [0.0; 2], &boundaries)
        .unwrap();
    let reaction = system
        .assemble_brinkman_mixed(params, &porous, &|_| [0.0; 2], &boundaries)
        .unwrap();
    let base_matrix = base.matrix().to_dense();
    let porous_matrix = reaction.matrix().to_dense();
    let mut velocity = vec![0.0; system.n_u];
    for (i, u) in velocity.iter_mut().enumerate() {
        if !mesh.is_boundary(i / 2) {
            *u = (i as f64).sin();
        }
    }
    let mut power = 0.0;
    for i in 0..system.n_u {
        for j in 0..system.n_u {
            let drag = porous_matrix[i * system.n + j] - base_matrix[i * system.n + j];
            if !mesh.is_boundary(i / 2) && !mesh.is_boundary(j / 2) {
                near(
                    drag,
                    porous_matrix[j * system.n + i] - base_matrix[j * system.n + i],
                    5e-12,
                );
            }
            power += velocity[i] * drag * velocity[j];
        }
    }
    assert!(power > 0.0);
    let integrated: f64 = reaction
        .cell_darcy_power_per_density(&velocity)
        .unwrap()
        .iter()
        .sum();
    near(power, integrated, integrated * 2e-12);
    let free = vec![BrinkmanCell::free_fluid(); mesh.tris.len()];
    let free_assembly = system
        .assemble_brinkman_mixed(params, &free, &|_| [0.0; 2], &boundaries)
        .unwrap();
    assert_eq!(free_assembly.matrix().to_dense(), base_matrix);
    assert_eq!(free_assembly.rhs(), base.rhs());
    assert_eq!(
        free_assembly
            .cell_darcy_power_per_density(&velocity)
            .unwrap(),
        vec![0.0; mesh.tris.len()]
    );
}

struct CountedMaterial {
    state: FluidTransportState,
    reads: Cell<usize>,
}
impl FluidTransport for CountedMaterial {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.reads.set(self.reads.get() + 1);
        Ok(self.state)
    }
}

#[test]
fn g3_original_material_viscosity_scales_flow_and_work_without_density_confusion() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 3));
    let system = FluxSystem::new(&mesh);
    let cells: Vec<_> = mesh
        .centroids
        .iter()
        .map(|point| {
            if point[0] < 0.5 {
                BrinkmanCell::try_new(tensor(0.53), 2.0).unwrap()
            } else {
                BrinkmanCell::free_fluid()
            }
        })
        .collect();
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(e, _)| mesh.is_boundary(*e))
        .map(|(edge, _)| MaterialStokesBoundary::Velocity {
            edge,
            value: &no_slip,
        })
        .collect();
    let mut results = Vec::new();
    for (rho, mu) in [(1.3, 0.2), (7.0, 0.4)] {
        let material = CountedMaterial {
            state: FluidTransportState::try_new(rho, 1000.0, 0.5, mu).unwrap(),
            reads: Cell::new(0),
        };
        let binding =
            MaterialFluxSystem::new(&system, &material, MaterialFluxNumerics::default()).unwrap();
        let solution = binding
            .solve_brinkman_mixed(&cells, &|x| [-x[1], x[0]], &boundaries)
            .unwrap();
        assert_eq!(material.reads.get(), 1);
        assert_eq!(solution.transport(), &material.state);
        assert_eq!(solution.cells(), cells.as_slice());
        for (power, cell) in solution.cell_darcy_power_w_per_m().iter().zip(&cells) {
            assert!(*power >= 0.0);
            if cell.permeability().is_none() {
                assert_eq!(*power, 0.0);
            }
        }
        results.push(solution);
    }
    let (a, b) = (&results[0], &results[1]);
    assert!(a.darcy_power_w_per_m() > 0.0);
    for (ua, ub) in a.velocity_dofs().iter().zip(b.velocity_dofs()) {
        near(*ub, 0.5 * ua, 2e-11);
    }
    for (pa, pb) in a.pressure_pa().iter().zip(b.pressure_pa()) {
        near(*pb, *pa, 2e-10);
    }
    near(
        b.darcy_power_w_per_m(),
        0.5 * a.darcy_power_w_per_m(),
        2e-11,
    );
}

#[test]
fn g0_seals_wrong_frames_and_unrepresentable_transport_are_refused() {
    assert!(BrinkmanCell::try_new(IntrinsicPermeability::isotropic(0.0).unwrap(), 1.0).is_err());
    for ratio in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(BrinkmanCell::try_new(tensor(0.0), ratio).is_err());
    }
    let r = 0.5_f64.sqrt();
    let tilted = IntrinsicPermeability::orthotropic(
        [1.0, 2.0, 3.0],
        [[r, 0.0, r], [0.0, 1.0, 0.0], [-r, 0.0, r]],
    )
    .unwrap();
    assert!(BrinkmanCell::try_new(tilted, 1.0).is_err());
    // An arbitrarily stiff unrelated y axis cannot conceal xz coupling.
    for scale in [1e-60, 1.0, 1e60] {
        let high_contrast = IntrinsicPermeability::orthotropic(
            [scale, 1e-30 * scale, 2.0 * scale],
            [[r, 0.0, r], [0.0, 1.0, 0.0], [-r, 0.0, r]],
        )
        .unwrap();
        assert!(BrinkmanCell::try_new(high_contrast, 1.0).is_err());
    }
    // A fully rotated isotropic tensor still leaves the xy plane invariant;
    // its off-plane terms cancel rather than being discarded axis by axis.
    let (s, c) = 0.37_f64.sin_cos();
    let (t, d) = 0.71_f64.sin_cos();
    let isotropic = IntrinsicPermeability::orthotropic(
        [7.0; 3],
        [[c * d, s * d, -t], [-s, c, 0.0], [c * t, s * t, d]],
    )
    .unwrap();
    assert!(BrinkmanCell::try_new(isotropic, 1.0).is_ok());
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(e, _)| mesh.is_boundary(*e))
        .map(|(edge, _)| StokesBoundary::Velocity {
            edge,
            value: &no_slip,
        })
        .collect();
    let cells = vec![
        BrinkmanCell::try_new(
            IntrinsicPermeability::isotropic(f64::MIN_POSITIVE).unwrap(),
            1.0
        )
        .unwrap();
        mesh.tris.len()
    ];
    assert!(
        system
            .assemble_brinkman_mixed(
                FluxParams {
                    nu: 100.0,
                    ..FluxParams::default()
                },
                &cells,
                &|_| [0.0; 2],
                &boundaries
            )
            .is_err()
    );
    assert!(
        system
            .assemble_brinkman_mixed(FluxParams::default(), &[], &|_| [0.0; 2], &boundaries)
            .is_err()
    );
    let cells = vec![BrinkmanCell::free_fluid(); mesh.tris.len()];
    assert!(
        system
            .assemble_brinkman_mixed(
                FluxParams::default(),
                &cells,
                &|_| [f64::NAN, 0.0],
                &boundaries
            )
            .is_err()
    );
}

#[test]
fn g0_drag_in_one_component_does_not_anchor_a_disconnected_free_fluid() {
    let mesh = TriMesh::from_mesh2(&Mesh2 {
        nodes: vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [2.0, 0.0],
            [3.0, 0.0],
            [2.0, 1.0],
        ],
        elems: vec![vec![0, 1, 2], vec![3, 4, 5]],
        patches: Vec::new(),
    });
    let system = FluxSystem::new(&mesh);
    let cells = [
        BrinkmanCell::try_new(tensor(0.0), 1.0).unwrap(),
        BrinkmanCell::free_fluid(),
    ];
    let traction = |point: [f64; 2], normal: [f64; 2]| {
        let pressure = if point[0] < 1.5 { 3.0 } else { 7.0 };
        [-pressure * normal[0], -pressure * normal[1]]
    };
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .map(|edge| StokesBoundary::TotalTraction {
            edge,
            value: &traction,
        })
        .collect();
    assert!(matches!(
        system.assemble_brinkman_mixed(FluxParams::default(), &cells, &|_| [0.0; 2], &boundaries),
        Err(fs_flux::StokesBoundaryError::UnanchoredComponent { component: 1 })
    ));
    // Grounding both domains through their own physical drag admits their
    // independently prescribed traction pressures without velocity pins.
    let cells = [cells[0]; 2];
    let solution = system
        .assemble_brinkman_mixed(FluxParams::default(), &cells, &|_| [0.0; 2], &boundaries)
        .unwrap()
        .solve()
        .unwrap();
    assert!(
        solution.field.x[..system.n_u]
            .iter()
            .all(|u| u.abs() < 1e-12)
    );
    assert!((solution.pressure[0] - 3.0).abs() < 1e-11);
    assert!((solution.pressure[1] - 7.0).abs() < 1e-11);
}

#[test]
fn g1_translating_anisotropic_skeleton_preserves_relative_flow_and_work() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 3));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        nu: 0.02,
        ..FluxParams::default()
    };
    let (s, c) = 0.53_f64.sin_cos();
    let kxx = 0.002 * c * c + 0.0005 * s * s;
    let kxy = (0.002 - 0.0005) * s * c;
    let kyy = 0.002 * s * s + 0.0005 * c * c;
    let gradient = [-3.0, 1.0];
    let relative = [
        -(kxx * gradient[0] + kxy * gradient[1]) / params.nu,
        -(kxy * gradient[0] + kyy * gradient[1]) / params.nu,
    ];
    let traction = |x: [f64; 2], n: [f64; 2]| {
        let pressure = 8.0 + gradient[0] * x[0] + gradient[1] * x[1];
        [-pressure * n[0], -pressure * n[1]]
    };
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| StokesBoundary::TotalTraction {
            edge,
            value: &traction,
        })
        .collect();
    for boost in [[0.0, 0.0], [0.7, -0.4], [-0.3, 0.6]] {
        let skeleton = [-0.9 + boost[0], 0.4 + boost[1]];
        let exact = [skeleton[0] + relative[0], skeleton[1] + relative[1]];
        let cell = BrinkmanCell::try_new(tensor(0.53), 1.0)
            .unwrap()
            .with_skeleton_velocity(skeleton)
            .unwrap();
        assert_eq!(cell.skeleton_velocity(), skeleton);
        let assembly = system
            .assemble_brinkman_mixed(
                params,
                &vec![cell; mesh.tris.len()],
                &|_| [0.0; 2],
                &boundaries,
            )
            .unwrap();
        let flow = assembly.solve().unwrap();
        let (error, divergence) = system.velocity_error(&flow.field.x, &|_| exact);
        assert!(error < 2e-10 && divergence < 2e-10);
        for (pressure, x) in flow.pressure.iter().zip(&mesh.centroids) {
            near(
                *pressure,
                8.0 + gradient[0] * x[0] + gradient[1] * x[1],
                2e-9,
            );
        }
        let velocity = &flow.field.x[..system.n_u];
        let powers = assembly.cell_darcy_power_per_density(velocity).unwrap();
        let exchanges = assembly.cell_darcy_exchange_per_density(velocity).unwrap();
        for ((exchange, area), power) in exchanges.iter().zip(&mesh.areas).zip(powers) {
            for a in 0..2 {
                near(exchange.skeleton_force[a], -area * gradient[a], 2e-10);
            }
            near(
                exchange.fluid_power,
                area * (gradient[0] * exact[0] + gradient[1] * exact[1]),
                2e-10,
            );
            near(
                exchange.skeleton_power,
                -area * (gradient[0] * skeleton[0] + gradient[1] * skeleton[1]),
                2e-10,
            );
            near(
                exchange.dissipation,
                -area * (gradient[0] * relative[0] + gradient[1] * relative[1]),
                2e-10,
            );
            near(
                exchange.fluid_power + exchange.skeleton_power + exchange.dissipation,
                0.0,
                2e-13,
            );
            assert_eq!(exchange.dissipation, power);
            if boost == [0.0, 0.0] {
                assert!(exchange.fluid_power > 0.0);
                assert!(exchange.skeleton_power < 0.0);
            }
        }
    }
}

#[test]
fn g1_cellwise_moving_medium_retains_material_force_and_boundary_reactions() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 3));
    let system = FluxSystem::new(&mesh);
    let (density, mu, ratio) = (7.0, 0.2, 1.0);
    let material = CountedMaterial {
        state: FluidTransportState::try_new(density, 1000.0, 0.5, mu).unwrap(),
        reads: Cell::new(0),
    };
    let (s, c) = 0.53_f64.sin_cos();
    let rxx = c * c / 0.002 + s * s / 0.0005;
    let rxy = s * c * (1.0 / 0.002 - 1.0 / 0.0005);
    let ryy = s * s / 0.002 + c * c / 0.0005;
    let velocity = |x: [f64; 2]| [0.3 + 0.2 * x[1], -0.1];
    let cells: Vec<_> = mesh
        .centroids
        .iter()
        .map(|point| {
            if point[0] > 0.7 {
                BrinkmanCell::free_fluid()
            } else {
                BrinkmanCell::try_new(tensor(0.53), ratio)
                    .unwrap()
                    .with_skeleton_velocity([0.8 - 0.3 * point[1], 0.1 + 0.2 * point[0]])
                    .unwrap()
            }
        })
        .collect();
    let force = |point: [f64; 2]| {
        let cell = cells[mesh.containing_cell(point).unwrap()];
        if cell.permeability().is_none() {
            return [0.0; 2];
        }
        let u = velocity(point);
        let frame = cell.skeleton_velocity();
        let relative = [u[0] - frame[0], u[1] - frame[1]];
        [
            mu * (rxx * relative[0] + rxy * relative[1]),
            mu * (rxy * relative[0] + ryy * relative[1]),
        ]
    };
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| MaterialStokesBoundary::Velocity {
            edge,
            value: &velocity,
        })
        .collect();
    let binding =
        MaterialFluxSystem::new(&system, &material, MaterialFluxNumerics::default()).unwrap();
    let solution = binding
        .solve_brinkman_mixed(&cells, &force, &boundaries)
        .unwrap();
    assert_eq!(material.reads.get(), 1);
    let (error, divergence) = system.velocity_error(&solution.discrete().x, &velocity);
    assert!(
        error < 2e-10 && divergence < 2e-10,
        "error={error:e}, divergence={divergence:e}"
    );
    assert!(
        solution
            .pressure_pa()
            .iter()
            .all(|pressure| pressure.abs() < 2e-8)
    );
    for (i, exchange) in solution.cell_darcy_exchange().iter().enumerate() {
        if cells[i].permeability().is_none() {
            assert_eq!(*exchange, fs_flux::MaterialBrinkmanDragExchange::default());
            continue;
        }
        let area = mesh.areas[i];
        let mean_u = velocity(mesh.centroids[i]);
        let frame = cells[i].skeleton_velocity();
        let mean_r = [mean_u[0] - frame[0], mean_u[1] - frame[1]];
        // Exact uniform-triangle second moment, independent of solver quadrature.
        let variance_y = mesh.tris[i]
            .iter()
            .map(|&vertex| (mesh.verts[vertex][1] - mesh.centroids[i][1]).powi(2))
            .sum::<f64>()
            / 12.0;
        let expected_force = [
            area * mu * (rxx * mean_r[0] + rxy * mean_r[1]),
            area * mu * (rxy * mean_r[0] + ryy * mean_r[1]),
        ];
        let expected_loss = area
            * mu
            * (rxx * (mean_r[0].powi(2) + 0.04 * variance_y)
                + 2.0 * rxy * mean_r[0] * mean_r[1]
                + ryy * mean_r[1].powi(2));
        for a in 0..2 {
            near(exchange.skeleton_force_n_per_m[a], expected_force[a], 2e-9);
        }
        near(exchange.dissipation_w_per_m, expected_loss, 2e-9);
        near(
            exchange.skeleton_power_w_per_m,
            expected_force[0] * frame[0] + expected_force[1] * frame[1],
            2e-9,
        );
        near(
            exchange.fluid_power_w_per_m
                + exchange.skeleton_power_w_per_m
                + exchange.dissipation_w_per_m,
            0.0,
            2e-12,
        );
        assert_eq!(
            exchange.dissipation_w_per_m,
            solution.cell_darcy_power_w_per_m()[i]
        );
    }
    assert_eq!(
        solution.darcy_exchange().dissipation_w_per_m,
        solution.darcy_power_w_per_m()
    );
    // The new volume load must also enter the original constrained reaction rows.
    // Constant shear with zero pressure has zero net boundary traction. Omitting
    // D*u_s from those rows would report a spurious reaction on the outer walls.
    let selected: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| edge)
        .collect();
    let load = solution
        .boundary_loads(&selected, [0.0; 2], &[Some(0.0)])
        .unwrap();
    near(load.stress.force_n_per_m[0], 0.0, 2e-8);
    near(load.stress.force_n_per_m[1], 0.0, 2e-8);
    near(load.stress.moment_n, 0.0, 2e-8);
}

#[test]
fn g0_comotion_and_explicit_stationary_skeleton_limits() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        nu: 0.02,
        ..FluxParams::default()
    };
    let stationary = BrinkmanCell::try_new(tensor(0.53), 1.0).unwrap();
    let velocity = |_: [f64; 2]| [0.7, -0.4];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| StokesBoundary::Velocity {
            edge,
            value: &velocity,
        })
        .collect();
    let base = system
        .assemble_brinkman_mixed(
            params,
            &vec![stationary; mesh.tris.len()],
            &|_| [0.0; 2],
            &boundaries,
        )
        .unwrap();
    let explicit = system
        .assemble_brinkman_mixed(
            params,
            &vec![stationary.with_skeleton_velocity([0.0; 2]).unwrap(); mesh.tris.len()],
            &|_| [0.0; 2],
            &boundaries,
        )
        .unwrap();
    assert_eq!(base.matrix().to_dense(), explicit.matrix().to_dense());
    assert_eq!(base.rhs(), explicit.rhs());
    let moving = system
        .assemble_brinkman_mixed(
            params,
            &vec![
                stationary
                    .with_skeleton_velocity(velocity([0.0; 2]))
                    .unwrap();
                mesh.tris.len()
            ],
            &|_| [0.0; 2],
            &boundaries,
        )
        .unwrap();
    assert_eq!(base.matrix().to_dense(), moving.matrix().to_dense());
    let flow = moving.solve().unwrap();
    let (error, divergence) = system.velocity_error(&flow.field.x, &velocity);
    assert!(error < 2e-11 && divergence < 2e-11);
    for exchange in moving
        .cell_darcy_exchange_per_density(&flow.field.x[..system.n_u])
        .unwrap()
    {
        assert!(
            exchange
                .skeleton_force
                .iter()
                .all(|force| force.abs() < 2e-11)
        );
        assert!(exchange.dissipation < 1e-23);
    }
    assert!(moving.cell_darcy_exchange_per_density(&[]).is_err());
    assert!(
        moving
            .cell_darcy_exchange_per_density(&vec![f64::NAN; system.n_u])
            .is_err()
    );
    assert!(
        BrinkmanCell::free_fluid()
            .with_skeleton_velocity([1.0, 0.0])
            .is_err()
    );
    assert_eq!(
        BrinkmanCell::free_fluid()
            .with_skeleton_velocity([0.0; 2])
            .unwrap(),
        BrinkmanCell::free_fluid()
    );
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(stationary.with_skeleton_velocity([0.0, bad]).is_err());
    }
    let overflowing = stationary
        .with_skeleton_velocity([f64::MAX, f64::MAX])
        .unwrap();
    assert!(
        system
            .assemble_brinkman_mixed(
                params,
                &vec![overflowing; mesh.tris.len()],
                &|_| [0.0; 2],
                &boundaries
            )
            .is_err()
    );
}
