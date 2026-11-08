//! Native EM force, exact P1 transfer and independent mechanical trajectories.
//! Synthetic constitutive laws test implementation, not experimental devices.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_couple::spatial_induction_vibration::{
    SpatialInductionSolidMap, SpatialInductionVibration, SpatialInductionVibrationControls,
    SpatialInductionVibrationError,
};
use fs_em::spatial::{
    SpatialMagneticMaterial,
    vector_potential::{
        SpatialVectorMagnetostatics, VectorMagneticBoundary, VectorMagneticControls,
        VectorPotentialFace,
        induction::{
            SpatialInduction, SpatialInductionLimits,
            harmonic::{SpatialHarmonicInductionControls, SpatialHarmonicInductionResponse},
            mechanical::SpatialLorentzLoadLimits,
        },
    },
};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    conductor::{
        ELECTRICAL_CONDUCTIVITY_DIMS, ELECTRICAL_CONDUCTIVITY_PROPERTY, ResolvedConductorState,
    },
    field_storage::{MAGNETIC_PERMEABILITY_DIMS, MAGNETIC_PERMEABILITY_PROPERTY},
    state_point::{
        MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
        resolve_material_state_point,
    },
};
use fs_math::c64::C64;
use fs_solid::{
    dynamics::ExactLinearSolidControls,
    linear3::{TetAssemblyBudget, TetElasticMaterial, TetLinearElasticProblem, TetMaterialField},
};
use fs_solver::CheckedCgConfig;
use fs_sparse::Coo;

fn with_cx(gate: &CancelGate, run: impl FnOnce(&Cx<'_>)) {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 85,
                kernel_id: 86,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance,
        "actual={actual:e}, expected={expected:e}, error={:e}",
        (actual - expected).abs()
    );
}
fn controls() -> SpatialInductionVibrationControls {
    SpatialInductionVibrationControls {
        loads: SpatialLorentzLoadLimits {
            max_vertices: 8,
            max_cells: 6,
        },
        structure: ExactLinearSolidControls {
            maximum_free_dofs: 32,
            maximum_preparation_bytes: 64 * 1024 * 1024,
            maximum_exponential_multiply_terms: 10_000_000_000,
            maximum_relative_operator_error: 1e-12,
            maximum_work_error_j: 1e-9,
            maximum_steps: 1000,
        },
        maximum_transfer_scalars: 10_000,
    }
}
fn with_field(alpha: f64, run: impl FnOnce(&SpatialHarmonicInductionResponse<'_, '_>)) {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value) in [
        (
            MAGNETIC_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            2.0,
        ),
        (
            ELECTRICAL_CONDUCTIVITY_PROPERTY,
            ELECTRICAL_CONDUCTIVITY_DIMS,
            5.0,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured induction vibration law".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let bundle = resolve_material_state_point(
        &card,
        &QueryPoint::new(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let conductor = ResolvedConductorState::try_from_material(bundle.clone()).unwrap();
    let mut positions = Vec::new();
    for z in [0.0, 1.0] {
        for y in [0.0, 1.0] {
            for x in [0.0, 1.0] {
                positions.push([x, y, z]);
            }
        }
    }
    let problem = SpatialVectorMagnetostatics::new(
        positions,
        vec![
            [0, 1, 3, 7],
            [0, 3, 2, 7],
            [0, 2, 6, 7],
            [0, 6, 4, 7],
            [0, 4, 5, 7],
            [0, 5, 1, 7],
        ],
        vec![SpatialMagneticMaterial::Linear(bundle); 6],
        WeightedAssemblyLimits {
            max_cells: 6,
            max_dofs: 100,
            max_triplets: 5000,
        },
        &mut || false,
    )
    .unwrap();
    let boundary = VectorMagneticBoundary {
        potential: problem
            .exterior_faces()
            .into_iter()
            .map(|face| VectorPotentialFace {
                face,
                edge_integrals_wb: [0.0; 3],
            })
            .collect(),
    };
    let cg = CheckedCgConfig {
        tolerance: 1e-12,
        max_iterations: 3000,
        max_corrections: 2,
    };
    let initial = problem
        .solve(
            &[[0.0; 3]; 6],
            &boundary,
            VectorMagneticControls {
                field: cg,
                gauge: cg,
                maximum_gauge_product_terms: 100_000,
                current_tolerance_a: 2e-8,
                flux_tolerance_wb: 2e-9,
                energy_tolerance_j: 2e-8,
            },
            &mut || false,
        )
        .unwrap();
    let state = SpatialInduction::new(
        &initial,
        vec![Some(conductor); 6],
        SpatialInductionLimits {
            maximum_unknowns: 100,
            maximum_conducting_cells: 6,
            maximum_assembly_terms: 100_000,
        },
        2e-9,
        &mut || false,
    )
    .unwrap();
    let velocity: Vec<_> = problem
        .positions_m()
        .iter()
        .map(|p| p.map(|x| alpha * x))
        .collect();
    let response = state
        .solve_harmonic(
            &[[C64::new(3.2, -2.4); 3]; 6],
            &velocity,
            3.0,
            SpatialHarmonicInductionControls {
                relative_tolerance: 2e-11,
                absolute_residual_a: 1e-11,
                restart: 32,
                maximum_cycles: 100,
                maximum_krylov_scalars: 100_000,
                maximum_field_scalars: 100_000,
                current_tolerance_a: 2e-8,
                flux_tolerance_wb: 2e-9,
                power_tolerance_w: 2e-8,
                maximum_cell_magnetic_peclet: 100.0,
                maximum_cell_skin_ratio: 1000.0,
            },
            &mut || false,
        )
        .unwrap();
    run(&response);
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(state.potential_wb(), initial.potential_wb());
}
fn body<'a>(
    nodes: &'a [[f64; 3]],
    tets: &'a [[usize; 4]],
    material: &'a TetElasticMaterial,
    fixed: &'a [usize],
) -> TetLinearElasticProblem<'a> {
    TetLinearElasticProblem {
        nodes_m: nodes,
        tetrahedra: tets,
        materials: TetMaterialField::Uniform(material),
        fixed_dofs: fixed,
        budget: TetAssemblyBudget {
            maximum_nodes: 8,
            maximum_tetrahedra: 6,
            maximum_free_dofs: 32,
            minimum_scaled_jacobian: 1e-12,
        },
    }
}
fn correspondence(vertices: &[usize]) -> SpatialInductionSolidMap<'_> {
    SpatialInductionSolidMap {
        field_vertices: vertices,
        field_cells: &[0],
    }
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn force_at(response: &SpatialHarmonicInductionResponse<'_, '_>, theta: f64) -> f64 {
    // Independent degree-two tetrahedron quadrature of N_1*(J(t) cross B(t))_x.
    let cell = response.cells().iter().find(|cell| cell.cell == 0).unwrap();
    let (sin, cos) = theta.sin_cos();
    let at = |value: C64| value.re * cos + value.im * sin;
    let b = response.flux_density_t()[0].map(at);
    let high = (5.0 + 3.0 * 5.0_f64.sqrt()) / 20.0;
    let low = (5.0 - 5.0_f64.sqrt()) / 20.0;
    let mut force = 0.0;
    for point in 0..4 {
        let shape: [f64; 4] = std::array::from_fn(|i| if i == point { high } else { low });
        let j = std::array::from_fn(|axis| {
            (0..4)
                .map(|i| shape[i] * at(cell.current_density_a_m2[i][axis]))
                .sum()
        });
        force += cell.volume_m3 / 4.0 * shape[1] * cross(j, b)[0];
    }
    force
}
fn force_coefficients(response: &SpatialHarmonicInductionResponse<'_, '_>) -> [f64; 3] {
    let first = force_at(response, 0.0);
    let opposite = force_at(response, std::f64::consts::FRAC_PI_2);
    let mean = 0.5 * (first + opposite);
    [
        mean,
        0.5 * (first - opposite),
        force_at(response, std::f64::consts::FRAC_PI_4) - mean,
    ]
}
fn snapshot(v: &SpatialInductionVibration<'_, '_, '_>) -> Vec<u64> {
    v.solid()
        .displacement_m()
        .iter()
        .chain(v.solid().velocity_mps())
        .chain(v.cosine_force_n())
        .chain(v.sine_force_n())
        .copied()
        .chain([
            v.solid().time_s(),
            v.solid().stored_energy_j(),
            v.solid().cumulative_work_j(),
            v.solid().cumulative_dissipated_energy_j(),
        ])
        .map(f64::to_bits)
        .chain([v.solid().accepted_steps()])
        .collect()
}

#[test]
fn physical_mean_and_twice_frequency_forces_drive_exact_motion_and_independent_work() {
    for alpha in [0.0, 0.35] {
        with_field(alpha, |response| {
            with_cx(&CancelGate::new_clock_free(), |cx| {
                let vertices = [0, 1, 3, 7];
                let nodes = vertices.map(|v| response.problem().positions_m()[v]);
                let fixed: Vec<_> = (0..12).filter(|&dof| dof != 3).collect();
                let material = TetElasticMaterial::try_new(2.0, 12.0, 0.0, [41; 32]).unwrap();
                let problem = body(&nodes, &[[0, 1, 2, 3]], &material, &fixed);
                let mut viscosity = Coo::new(1, 1);
                let damping = 0.04;
                viscosity.push(0, 0, damping);
                let viscosity = viscosity.assemble();
                let (q0, v0) = (0.01, 0.03);
                let mut u = [[0.0; 3]; 4];
                u[1][0] = q0;
                let mut v = [[0.0; 3]; 4];
                v[1][0] = v0;
                let mut vibration = SpatialInductionVibration::prepare(
                    response,
                    &problem,
                    correspondence(&vertices),
                    Some(&viscosity),
                    0.1,
                    &u,
                    &v,
                    controls(),
                    cx,
                )
                .unwrap();
                assert!(std::ptr::eq(vibration.field(), response));
                assert!(std::ptr::eq(vibration.problem(), &problem));
                assert_eq!(vibration.solid().free_dofs(), &[3]);
                assert_eq!(vibration.solid().angular_frequencies_rad_s(), &[0.0, 6.0]);
                let [f0, fr, fi] = force_coefficients(response);
                near(vibration.cosine_force_n()[0], f0, 2e-10);
                near(vibration.cosine_force_n()[1], fr, 2e-10);
                near(vibration.sine_force_n()[1], fi, 2e-10);
                assert!(fr.hypot(fi) > 1e-5);
                if alpha == 0.0 {
                    near(f0, 0.0, 2e-9);
                } else {
                    assert!(f0.abs() > 1e-6);
                }
                // N_1=x-y, V=1/6, rho=2: M=rho*V/10=1/30.
                // E=12, nu=0: K=V*mu*(|grad N_1|²+grad_x²)=3.
                let (mass, stiffness, omega) = (1.0 / 30.0, 3.0, 6.0);
                let amplitude =
                    C64::new(fr, fi) / C64::new(stiffness - mass * omega * omega, -damping * omega);
                let decay = damping / (2.0 * mass);
                let natural = (stiffness / mass - decay * decay).sqrt();
                let a = q0 - f0 / stiffness - amplitude.re;
                let b = (v0 - omega * amplitude.im + decay * a) / natural;
                let trajectory = |time: f64| {
                    let (sn, cs) = (natural * time).sin_cos();
                    let (s, c) = (omega * time).sin_cos();
                    let envelope = (-decay * time).exp();
                    let q = f0 / stiffness
                        + amplitude.re * c
                        + amplitude.im * s
                        + envelope * (a * cs + b * sn);
                    let v = omega * (amplitude.im * c - amplitude.re * s)
                        + envelope
                            * ((natural * b - decay * a) * cs - (natural * a + decay * b) * sn);
                    (q, v, f0 + fr * c + fi * s)
                };
                for _ in 0..13 {
                    let step = vibration.advance(cx).unwrap();
                    let (q, v, _) = trajectory(step.time_s);
                    near(vibration.solid().displacement_m()[0], q, 3e-10);
                    near(vibration.solid().velocity_mps()[0], v, 3e-9);
                }
                let slices = 4096;
                let (mut work, mut heat) = (0.0, 0.0);
                for i in 0..=slices {
                    let (_, v, f) = trajectory(1.3 * i as f64 / slices as f64);
                    let weight = if i == 0 || i == slices {
                        1.0
                    } else if i % 2 == 0 {
                        2.0
                    } else {
                        4.0
                    };
                    work += weight * f * v;
                    heat += weight * damping * v * v;
                }
                work *= 1.3 / (3 * slices) as f64;
                heat *= 1.3 / (3 * slices) as f64;
                near(vibration.solid().cumulative_work_j(), work, 3e-10);
                near(
                    vibration.solid().cumulative_dissipated_energy_j(),
                    heat,
                    3e-10,
                );
                // Joule heating belongs to the unchanged EM source, not the solid's C.
                assert!(response.balance().mean_ohmic_power_w > 0.0);
                let now = vibration.solid().time_s();
                let (_, _, force) = trajectory(now);
                near(
                    vibration.cosine_force_n()[0] + vibration.cosine_force_n()[1],
                    force,
                    2e-10,
                );
            })
        });
    }
}

#[test]
fn twice_frequency_resonance_keeps_secular_growth_and_phase_across_partitions() {
    with_field(0.0, |response| {
        with_cx(&CancelGate::new_clock_free(), |cx| {
            let vertices = [0, 1, 3, 7];
            let nodes = vertices.map(|v| response.problem().positions_m()[v]);
            let fixed: Vec<_> = (0..12).filter(|&dof| dof != 3).collect();
            let material = TetElasticMaterial::try_new(2.0, 4.8, 0.0, [42; 32]).unwrap();
            let problem = body(&nodes, &[[0, 1, 2, 3]], &material, &fixed);
            let prepare = |dt| {
                SpatialInductionVibration::prepare(
                    response,
                    &problem,
                    correspondence(&vertices),
                    None,
                    dt,
                    &[[0.0; 3]; 4],
                    &[[0.0; 3]; 4],
                    controls(),
                    cx,
                )
                .unwrap()
            };
            let [f0, fr, fi] = force_coefficients(response);
            let (mass, omega) = (1.0 / 30.0, 6.0);
            let mut split = prepare(0.07);
            for _ in 0..25 {
                split.advance(cx).unwrap();
            }
            let time = split.solid().time_s();
            let (s, c) = (omega * time).sin_cos();
            let q = f0 / (mass * omega * omega) * (1.0 - c)
                + fr / (2.0 * mass * omega) * time * s
                + fi / (2.0 * mass * omega) * (s / omega - time * c);
            let v = f0 / (mass * omega) * s
                + fr / (2.0 * mass * omega) * (s + omega * time * c)
                + fi / (2.0 * mass) * time * s;
            near(split.solid().displacement_m()[0], q, 3e-9);
            near(split.solid().velocity_mps()[0], v, 3e-8);
            let mut whole = prepare(time);
            whole.advance(cx).unwrap();
            near(whole.solid().displacement_m()[0], q, 3e-9);
            near(whole.solid().velocity_mps()[0], v, 3e-8);
            near(
                whole.solid().cumulative_work_j(),
                split.solid().cumulative_work_j(),
                3e-9,
            );
            near(
                whole.solid().cumulative_work_j(),
                0.5 * mass * (v * v + omega * omega * q * q),
                3e-9,
            );
        })
    });
}

#[test]
fn original_mesh_permutations_preserve_motion_and_refusals_preserve_phase() {
    with_field(0.0, |response| {
        let gate = CancelGate::new_clock_free();
        with_cx(&gate, |cx| {
            let vertices = [0, 1, 3, 7];
            let nodes = vertices.map(|v| response.problem().positions_m()[v]);
            let fixed: Vec<_> = (0..12).filter(|&dof| dof != 3).collect();
            let material = TetElasticMaterial::try_new(2.0, 12.0, 0.0, [43; 32]).unwrap();
            let problem = body(&nodes, &[[0, 1, 2, 3]], &material, &fixed);
            let prepare = |map, control| {
                SpatialInductionVibration::prepare(
                    response,
                    &problem,
                    map,
                    None,
                    0.11,
                    &[[0.0; 3]; 4],
                    &[[0.0; 3]; 4],
                    control,
                    cx,
                )
            };
            let mut ordinary = prepare(correspondence(&vertices), controls()).unwrap();
            let permuted_vertices = [7, 3, 1, 0];
            let permuted_nodes = permuted_vertices.map(|v| response.problem().positions_m()[v]);
            let permuted_fixed: Vec<_> = (0..12).filter(|&dof| dof != 6).collect();
            let permuted_problem =
                body(&permuted_nodes, &[[3, 2, 1, 0]], &material, &permuted_fixed);
            let mut permuted = SpatialInductionVibration::prepare(
                response,
                &permuted_problem,
                correspondence(&permuted_vertices),
                None,
                0.11,
                &[[0.0; 3]; 4],
                &[[0.0; 3]; 4],
                controls(),
                cx,
            )
            .unwrap();
            for _ in 0..4 {
                ordinary.advance(cx).unwrap();
                permuted.advance(cx).unwrap();
                near(
                    ordinary.solid().displacement_m()[0],
                    permuted.solid().displacement_m()[0],
                    2e-11,
                );
                near(
                    ordinary.solid().velocity_mps()[0],
                    permuted.solid().velocity_mps()[0],
                    2e-10,
                );
                near(
                    ordinary.solid().cumulative_work_j(),
                    permuted.solid().cumulative_work_j(),
                    2e-11,
                );
            }
            for map in [
                correspondence(&[0, 0, 3, 7]),
                correspondence(&[0, 2, 3, 7]),
                SpatialInductionSolidMap {
                    field_vertices: &vertices,
                    field_cells: &[1],
                },
            ] {
                assert!(prepare(map, controls()).is_err());
            }
            assert!(
                prepare(
                    correspondence(&vertices),
                    SpatialInductionVibrationControls {
                        maximum_transfer_scalars: 0,
                        ..controls()
                    }
                )
                .is_err()
            );
            let mut moving_support = [[0.0; 3]; 4];
            moving_support[0][0] = 1.0;
            assert!(
                SpatialInductionVibration::prepare(
                    response,
                    &problem,
                    correspondence(&vertices),
                    None,
                    0.1,
                    &moving_support,
                    &[[0.0; 3]; 4],
                    controls(),
                    cx
                )
                .is_err()
            );
            let mut one_step = controls();
            one_step.structure.maximum_steps = 1;
            let mut limited = prepare(correspondence(&vertices), one_step).unwrap();
            limited.advance(cx).unwrap();
            let accepted = snapshot(&limited);
            assert!(limited.advance(cx).is_err());
            assert_eq!(snapshot(&limited), accepted);
            let accepted = snapshot(&ordinary);
            gate.request();
            assert!(matches!(
                ordinary.advance(cx),
                Err(SpatialInductionVibrationError::Cancelled)
            ));
            assert_eq!(snapshot(&ordinary), accepted);
        });
    });
}
