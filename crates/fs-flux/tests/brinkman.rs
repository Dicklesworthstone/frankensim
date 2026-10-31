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
    let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .map(|edge| StokesBoundary::TotalTraction {
            edge,
            value: &traction,
        })
        .collect();
    assert!(matches!(
        system.assemble_brinkman_mixed(FluxParams::default(), &cells, &|_| [0.0; 2], &boundaries),
        Err(fs_flux::StokesBoundaryError::DisconnectedMesh)
    ));
}
