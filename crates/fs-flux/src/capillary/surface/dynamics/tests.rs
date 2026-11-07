use super::*;

fn controls() -> SurfaceDynamicsControls {
    let mut solve = ImplicitSolveConfig::default();
    solve.max_newton_iterations = 32;
    solve.newton.absolute_tolerance = 1e-13;
    solve.newton.relative_tolerance = 1e-12;
    SurfaceDynamicsControls {
        solve,
        maximum_krylov_restart: 32,
        impulse_tolerance_n_s: 1e-11,
        kinematic_tolerance_m: 1e-12,
        energy_tolerance_j: 1e-11,
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}

fn positions(state: &SecondOrderState) -> Vec<Vector> {
    state
        .q
        .chunks_exact(3)
        .map(|v| [v[0], v[1], v[2]])
        .collect()
}

fn momentum(state: &SecondOrderState, masses: &[f64]) -> (Vector, Vector) {
    let mut linear = [0.0; 3];
    let mut angular = [0.0; 3];
    for ((q, v), mass) in state
        .q
        .chunks_exact(3)
        .zip(state.v.chunks_exact(3))
        .zip(masses)
    {
        let p = [mass * v[0], mass * v[1], mass * v[2]];
        let l = cross([q[0], q[1], q[2]], p);
        for d in 0..3 {
            linear[d] += p[d];
            angular[d] += l[d];
        }
    }
    (linear, angular)
}

fn tetrahedron() -> ([Vector; 4], [[usize; 3]; 4], [f64; 4], [Vector; 4]) {
    let q = [
        [1.1, 0.1, -0.1],
        [-0.4, 0.9, 0.2],
        [-0.5, -0.5, 0.4],
        [0.1, -0.2, 1.3],
    ];
    let faces = [[0, 1, 2], [0, 3, 1], [0, 2, 3], [1, 3, 2]];
    let masses = [0.7, 1.1, 1.3, 0.9];
    let v = core::array::from_fn(|i| {
        let spin = cross([0.7, -0.4, 1.8], q[i]);
        [
            spin[0] + 0.2 * q[i][1],
            spin[1] - 0.1 * q[i][2],
            spin[2] + 0.15 * q[i][0],
        ]
    });
    (q, faces, masses, v)
}

#[test]
fn insoluble_coating_breathes_with_its_thermodynamic_tension_and_no_external_supply() {
    use super::super::surfactant::SurfactantLaw;
    let q = [
        [1.0, 0.0, 0.0],
        [-0.5, 3.0_f64.sqrt() / 2.0, 0.0],
        [-0.5, -3.0_f64.sqrt() / 2.0, 0.0],
    ];
    let masses = [0.8; 3];
    let faces = [[0, 1, 2]];
    let model = CapillarySurfaceDynamics::try_new(&masses, &faces, 3, 1).unwrap();
    let amounts = [0.1];
    let coating = SurfactantField::Uniform {
        law: SurfactantLaw::IdealDilute {
            clean_tension_n_m: 0.6,
            thermal_energy_j_mol: 1.0,
            reference_concentration_mol_m2: 0.1,
        },
        amounts_mol: &amounts,
    };
    let inviscid = SurfaceViscosityField::Uniform(SurfaceViscosity::ZERO);
    let force = [[0.0; 3]; 3];
    let speed = 0.2;
    let velocity = q.map(|p| p.map(|x| speed * x));
    let mut state = model
        .initial_state_surfactant(0.0, &q, &velocity, coating, inviscid, &force)
        .unwrap();
    let dt = 0.08;
    let area = 3.0 * 3.0_f64.sqrt() / 4.0;
    // Independently reduce E(r)=sigma0*c*r² + RT*N*(ln(N/(Gamma_ref*c*r²))-1)
    // to its scalar midpoint momentum equation. The area secant, rather than
    // the mean endpoint tension, is required for this finite-step identity.
    let balance = |radius: f64| {
        let ratio = radius * radius;
        let secant = if ratio == 1.0 {
            0.6 - amounts[0] / area
        } else {
            0.6 - amounts[0] * ratio.ln() / (area * (ratio - 1.0))
        };
        2.0 * masses[0] * (radius - 1.0 - dt * speed) / (dt * dt)
            + secant * 3.0_f64.sqrt() * (1.0 + radius) / 4.0
    };
    let (mut lo, mut hi) = (0.95, 1.05);
    assert!(balance(lo) < 0.0 && balance(hi) > 0.0);
    for _ in 0..60 {
        let middle = 0.5 * lo + 0.5 * hi;
        if balance(middle) > 0.0 {
            hi = middle;
        } else {
            lo = middle;
        }
    }
    let radius = 0.5 * lo + 0.5 * hi;
    let final_speed = 2.0 * (radius - 1.0) / dt - speed;
    let report = model
        .advance_surfactant(
            &mut state,
            coating,
            inviscid,
            &force,
            dt,
            controls(),
            || false,
        )
        .unwrap();
    for i in 0..3 {
        for d in 0..3 {
            close(state.q[3 * i + d], radius * q[i][d], 4e-12);
            close(state.v[3 * i + d], final_speed * q[i][d], 4e-12);
        }
    }
    assert_eq!(report.surface.coefficient_supply_j, 0.0);
    assert_eq!(amounts, [0.1]);
    close(
        report.kinetic_change_j + report.surface.energy_change_j,
        0.0,
        5e-13,
    );
}

#[test]
fn nonuniform_conserved_coating_preserves_momentum_and_closes_viscous_surface_energy() {
    use super::super::surfactant::SurfactantLaw;
    let (q, faces, masses, v) = tetrahedron();
    let model = CapillarySurfaceDynamics::try_new(&masses, &faces, 4, 4).unwrap();
    let amounts = [0.02, 0.07, 0.03, 0.05];
    let coating = SurfactantField::Uniform {
        law: SurfactantLaw::Langmuir {
            clean_tension_n_m: 0.3,
            thermal_energy_j_mol: 0.5,
            saturation_concentration_mol_m2: 0.8,
        },
        amounts_mol: &amounts,
    };
    let viscosity = SurfaceViscosityField::Uniform(SurfaceViscosity {
        shear_n_s_m: 0.1,
        dilatational_n_s_m: 0.15,
    });
    let force = [[0.0; 3]; 4];
    let mut state = model
        .initial_state_surfactant(0.0, &q, &v, coating, viscosity, &force)
        .unwrap();
    let initial = kinetic(&masses, &state.v).unwrap()
        + assemble_surfactant_step(&q, &q, &faces, coating)
            .unwrap()
            .energy
            .energy_before_j;
    let (linear, angular) = momentum(&state, &masses);
    let mut losses = 0.0;
    for _ in 0..12 {
        let report = model
            .advance_surfactant(
                &mut state,
                coating,
                viscosity,
                &force,
                0.005,
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(report.surface.coefficient_supply_j, 0.0);
        losses += report.surface_viscous_dissipation_j;
        close(
            report.kinetic_after_j + report.surface.energy_after_j + losses,
            initial,
            5e-11,
        );
        let (p, l) = momentum(&state, &masses);
        for d in 0..3 {
            close(p[d], linear[d], 3e-11);
            close(l[d], angular[d], 3e-11);
        }
    }
    assert!(losses > 0.0);
    assert_eq!(amounts, [0.02, 0.07, 0.03, 0.05]);
    let before = state.clone();
    let overfilled = [100.0; 4];
    let invalid = SurfactantField::Uniform {
        law: SurfactantLaw::Langmuir {
            clean_tension_n_m: 0.3,
            thermal_energy_j_mol: 0.5,
            saturation_concentration_mol_m2: 0.8,
        },
        amounts_mol: &overfilled,
    };
    assert!(matches!(
        model.advance_surfactant(
            &mut state,
            invalid,
            viscosity,
            &force,
            0.005,
            controls(),
            || false
        ),
        Err(SurfaceDynamicsError::Surfactant(_))
    ));
    assert_eq!(state.q, before.q);
    assert_eq!(state.v, before.v);
    assert_eq!(state.a, before.a);
    assert_eq!(state.t, before.t);
    assert_eq!(state.steps, before.steps);
    assert_eq!(state.history.len(), before.history.len());
}

#[test]
fn moving_coating_diffusion_closes_the_complete_energy_and_amount_account() {
    use super::super::surfactant::{
        SurfactantLaw,
        diffusion::{SurfaceDiffusionControls, SurfaceDiffusionEdge},
    };
    use super::surfactant_transport::SurfaceSurfactantTransport;
    let (q, faces, masses, v) = tetrahedron();
    let model = CapillarySurfaceDynamics::try_new(&masses, &faces, 4, 4).unwrap();
    let law = SurfactantLaw::Langmuir {
        clean_tension_n_m: 0.3,
        thermal_energy_j_mol: 0.5,
        saturation_concentration_mol_m2: 0.8,
    };
    let mut amounts = [0.02, 0.07, 0.03, 0.05];
    let initial_amount = amounts.iter().sum::<f64>();
    let coating = SurfactantField::Uniform {
        law,
        amounts_mol: &amounts,
    };
    let viscosity = SurfaceViscosityField::Uniform(SurfaceViscosity {
        shear_n_s_m: 0.1,
        dilatational_n_s_m: 0.15,
    });
    let external = [[0.0; 3]; 4];
    let mut state = model
        .initial_state_surfactant(0.0, &q, &v, coating, viscosity, &external)
        .unwrap();
    let initial_energy = kinetic(&masses, &state.v).unwrap()
        + assemble_surfactant_step(&q, &q, &faces, coating)
            .unwrap()
            .energy
            .energy_before_j;
    let (linear, angular) = momentum(&state, &masses);
    let edges = [
        SurfaceDiffusionEdge {
            faces: [0, 1],
            conductance_m2_s: 0.4,
        },
        SurfaceDiffusionEdge {
            faces: [1, 2],
            conductance_m2_s: 0.2,
        },
        SurfaceDiffusionEdge {
            faces: [2, 3],
            conductance_m2_s: 0.3,
        },
    ];
    let transport = SurfaceSurfactantTransport {
        law,
        edges: &edges,
        controls: SurfaceDiffusionControls {
            maximum_step_s: 0.002,
            maximum_substeps: 8,
            maximum_faces: 4,
            maximum_edges: 3,
            amount_tolerance_mol: 1e-13,
            energy_increase_tolerance_j: 1e-13,
        },
    };
    let mut losses = 0.0;
    for _ in 0..12 {
        let report = model
            .advance_surfactant_diffusive(
                &mut state,
                &mut amounts,
                transport,
                viscosity,
                &external,
                0.005,
                controls(),
                || false,
            )
            .unwrap();
        assert!(report.diffusion_free_energy_dissipation_j > 0.0);
        assert!(report.mechanics.surface_viscous_dissipation_j > 0.0);
        assert_eq!(report.mechanics.surface.coefficient_supply_j, 0.0);
        losses += report.diffusion_free_energy_dissipation_j
            + report.mechanics.surface_viscous_dissipation_j;
        close(
            report.mechanics.kinetic_after_j + report.surface_free_energy_after_j + losses,
            initial_energy,
            8e-11,
        );
        close(amounts.iter().sum(), initial_amount, 1e-13);
        close(report.energy_residual_j, 0.0, 1e-11);
        let positions = positions(&state);
        let stored = assemble_surfactant_step(
            &positions,
            &positions,
            &faces,
            SurfactantField::Uniform {
                law,
                amounts_mol: &amounts,
            },
        )
        .unwrap();
        close(
            stored.energy.energy_after_j,
            report.surface_free_energy_after_j,
            1e-13,
        );
        let (p, l) = momentum(&state, &masses);
        for d in 0..3 {
            close(p[d], linear[d], 3e-11);
            close(l[d], angular[d], 3e-11);
        }
    }

    // Observe the final callback of a successful deterministic trial, then
    // cancel that same trial immediately before it would publish either state.
    let original = state.clone();
    let original_amounts = amounts;
    let mut calls = 0;
    model
        .advance_surfactant_diffusive(
            &mut state,
            &mut amounts,
            transport,
            viscosity,
            &external,
            0.005,
            controls(),
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    let mut retry = original.clone();
    let mut retry_amounts = original_amounts;
    let mut current = 0;
    let result = model.advance_surfactant_diffusive(
        &mut retry,
        &mut retry_amounts,
        transport,
        viscosity,
        &external,
        0.005,
        controls(),
        || {
            current += 1;
            current == calls
        },
    );
    assert!(matches!(result, Err(SurfaceDynamicsError::Cancelled)));
    assert_eq!(retry_amounts, original_amounts);
    assert_eq!(retry.q, original.q);
    assert_eq!(retry.v, original.v);
    assert_eq!(retry.a, original.a);
    assert_eq!(retry.t, original.t);
    assert_eq!(retry.steps, original.steps);
    assert_eq!(retry.history.len(), original.history.len());
}

#[test]
fn empty_diffusion_graph_recovers_the_original_coated_mechanical_trajectory() {
    use super::super::surfactant::{SurfactantLaw, diffusion::SurfaceDiffusionControls};
    use super::surfactant_transport::SurfaceSurfactantTransport;
    let (q, faces, masses, v) = tetrahedron();
    let model = CapillarySurfaceDynamics::try_new(&masses, &faces, 4, 4).unwrap();
    let law = SurfactantLaw::IdealDilute {
        clean_tension_n_m: 0.3,
        thermal_energy_j_mol: 0.5,
        reference_concentration_mol_m2: 0.1,
    };
    let mut amounts = [0.02, 0.07, 0.03, 0.05];
    let inviscid = SurfaceViscosityField::Uniform(SurfaceViscosity::ZERO);
    let external = [[0.0; 3]; 4];
    let coating = SurfactantField::Uniform {
        law,
        amounts_mol: &amounts,
    };
    let mut reference = model
        .initial_state_surfactant(0.0, &q, &v, coating, inviscid, &external)
        .unwrap();
    let mut composed = reference.clone();
    model
        .advance_surfactant(
            &mut reference,
            coating,
            inviscid,
            &external,
            0.005,
            controls(),
            || false,
        )
        .unwrap();
    let report = model
        .advance_surfactant_diffusive(
            &mut composed,
            &mut amounts,
            SurfaceSurfactantTransport {
                law,
                edges: &[],
                controls: SurfaceDiffusionControls {
                    maximum_step_s: 0.002,
                    maximum_substeps: 8,
                    maximum_faces: 4,
                    maximum_edges: 0,
                    amount_tolerance_mol: 0.0,
                    energy_increase_tolerance_j: 0.0,
                },
            },
            inviscid,
            &external,
            0.005,
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(reference.q, composed.q);
    assert_eq!(reference.v, composed.v);
    assert_eq!(reference.a, composed.a);
    assert_eq!(amounts, [0.02, 0.07, 0.03, 0.05]);
    assert_eq!(report.diffusion_free_energy_dissipation_j, 0.0);
    assert_eq!(report.amount_defect_mol, 0.0);
}

fn partition_laws(langmuir: bool) -> [super::super::surfactant::SurfactantLaw; 4] {
    use super::super::surfactant::SurfactantLaw;
    let clean = [0.3, 0.4, 0.5, 0.6];
    let references = [0.4, 0.8, 1.2, 0.6];
    core::array::from_fn(|i| {
        if langmuir {
            SurfactantLaw::Langmuir {
                clean_tension_n_m: clean[i],
                thermal_energy_j_mol: 0.5,
                saturation_concentration_mol_m2: references[i],
            }
        } else {
            SurfactantLaw::IdealDilute {
                clean_tension_n_m: clean[i],
                thermal_energy_j_mol: 0.5,
                reference_concentration_mol_m2: references[i],
            }
        }
    })
}

fn diffusion_controls(
    edges: usize,
) -> super::super::surfactant::diffusion::SurfaceDiffusionControls {
    super::super::surfactant::diffusion::SurfaceDiffusionControls {
        maximum_step_s: 0.002,
        maximum_substeps: 8,
        maximum_faces: 4,
        maximum_edges: edges,
        amount_tolerance_mol: 1e-13,
        energy_increase_tolerance_j: 1e-13,
    }
}

fn adaptive_diffusion_controls(
    edges: usize,
) -> super::super::surfactant::diffusion::adaptive::AdaptiveSurfaceDiffusionControls {
    super::super::surfactant::diffusion::adaptive::AdaptiveSurfaceDiffusionControls {
        initial_step_s: 0.02,
        minimum_step_s: 1e-10,
        maximum_step_s: 0.02,
        absolute_concentration_tolerance_mol_m2: 1e-9,
        relative_tolerance: 1e-7,
        maximum_trials: 1_000,
        maximum_pair_exchanges: 15_000,
        maximum_faces: 4,
        maximum_edges: edges,
        amount_tolerance_mol: 1e-13,
        energy_increase_tolerance_j: 1e-13,
    }
}

#[test]
fn adaptive_partitioned_motion_uses_solved_areas_and_independent_exact_pair_flows() {
    use super::super::surfactant::{SurfactantFace, diffusion::SurfacePartitionDiffusionEdge};
    use super::surfactant_transport::AdaptiveSurfacePartitionedSurfactantTransport;
    let (q, faces, masses, _) = tetrahedron();
    let velocity = q.map(|p| p.map(|x| 0.8 * x));
    let model = CapillarySurfaceDynamics::try_new(&masses, &faces, 4, 4).unwrap();
    let laws = partition_laws(false);
    let original_amounts = [0.02, 0.07, 0.03, 0.05];
    let coating = |amounts: [f64; 4]| {
        core::array::from_fn::<_, 4, _>(|i| SurfactantFace {
            law: laws[i],
            amount_mol: amounts[i],
        })
    };
    let areas = |p: &[Vector]| {
        core::array::from_fn::<_, 4, _>(|i| geometry(faces[i].map(|j| p[j])).unwrap().area)
    };
    // This exact two-reservoir solution is independent of the production pair
    // map and controller. References set chemical capacities, not raw area.
    let exact_pair = |amounts: [f64; 4], area: [f64; 4], duration: f64| {
        let capacity = [0.4 * area[0], 0.8 * area[1]];
        let total = amounts[0] + amounts[1];
        let equilibrium = total * capacity[0] / (capacity[0] + capacity[1]);
        let decay = (-0.7 * (1.0 / capacity[0] + 1.0 / capacity[1]) * duration).exp();
        let first = equilibrium + (amounts[0] - equilibrium) * decay;
        [first, total - first, amounts[2], amounts[3]]
    };
    let storage = |p: &[Vector], amounts: [f64; 4]| {
        coating(amounts)
            .into_iter()
            .zip(areas(p))
            .map(|(face, area)| face.at_area(area).unwrap().free_energy_j)
            .sum::<f64>()
    };
    let inviscid = SurfaceViscosityField::Uniform(SurfaceViscosity::ZERO);
    let external = [[0.0; 3]; 4];
    let initial = model
        .initial_state_surfactant(
            0.0,
            &q,
            &velocity,
            SurfactantField::PerFace(&coating(original_amounts)),
            inviscid,
            &external,
        )
        .unwrap();
    let initial_energy = kinetic(&masses, &initial.v).unwrap() + storage(&q, original_amounts);
    let dt = 0.03;
    let intermediate_amounts = exact_pair(original_amounts, areas(&q), 0.5 * dt);
    let mut reference = initial.clone();
    model
        .advance_surfactant(
            &mut reference,
            SurfactantField::PerFace(&coating(intermediate_amounts)),
            inviscid,
            &external,
            dt,
            controls(),
            || false,
        )
        .unwrap();
    let expected = exact_pair(
        intermediate_amounts,
        areas(&positions(&reference)),
        0.5 * dt,
    );
    let held_area = exact_pair(original_amounts, areas(&q), dt);
    assert!((expected[0] - held_area[0]).abs() > 1e-7);

    let mut state = initial;
    let mut amounts = original_amounts;
    let mut adaptive = adaptive_diffusion_controls(1);
    adaptive.maximum_trials = 2;
    adaptive.maximum_pair_exchanges = 6;
    let report = model
        .advance_partitioned_surfactant_diffusive_adaptive(
            &mut state,
            &mut amounts,
            AdaptiveSurfacePartitionedSurfactantTransport {
                laws: &laws,
                edges: &[SurfacePartitionDiffusionEdge {
                    faces: [0, 1],
                    conductance_mol_s: 0.7,
                }],
                controls: adaptive,
            },
            inviscid,
            &external,
            dt,
            controls(),
            || false,
        )
        .unwrap();
    for (&actual, expected) in amounts.iter().zip(expected) {
        close(actual, expected, 5e-13);
    }
    for (actual, expected) in state
        .q
        .iter()
        .zip(&reference.q)
        .chain(state.v.iter().zip(&reference.v))
    {
        close(*actual, *expected, 5e-12);
    }
    close(amounts.iter().sum(), original_amounts.iter().sum(), 2e-14);
    assert!(amounts.iter().all(|n| *n > 0.0));
    close(
        kinetic(&masses, &state.v).unwrap()
            + storage(&positions(&state), amounts)
            + report.coupled.diffusion_free_energy_dissipation_j,
        initial_energy,
        3e-11,
    );
    assert_eq!(report.total_trials, 2);
    assert_eq!(report.total_evaluated_pair_exchanges, 6);
    assert_eq!(report.diffusion_before.evaluated_pair_exchanges, 3);
    assert_eq!(report.diffusion_after.evaluated_pair_exchanges, 3);
    assert_eq!(
        report.coupled.diffusion_before,
        report.diffusion_before.diffusion
    );
    assert_eq!(
        report.coupled.diffusion_after,
        report.diffusion_after.diffusion
    );
}

#[test]
fn adaptive_partition_equilibrium_moves_without_spurious_chemical_redistribution() {
    use super::super::surfactant::{SurfactantFace, diffusion::SurfacePartitionDiffusionEdge};
    use super::surfactant_transport::AdaptiveSurfacePartitionedSurfactantTransport;
    let (q, faces, masses, _) = tetrahedron();
    let model = CapillarySurfaceDynamics::try_new(&masses, &faces, 4, 4).unwrap();
    let translation = [0.2, -0.1, 0.3];
    let velocity = [translation; 4];
    let inviscid = SurfaceViscosityField::Uniform(SurfaceViscosity::ZERO);
    let edges = core::array::from_fn::<_, 3, _>(|i| SurfacePartitionDiffusionEdge {
        faces: [i, i + 1],
        conductance_mol_s: 1e6,
    });
    for langmuir in [false, true] {
        let laws = partition_laws(langmuir);
        let references = [0.4, 0.8, 1.2, 0.6];
        let original_amounts = core::array::from_fn::<_, 4, _>(|i| {
            geometry(faces[i].map(|j| q[j])).unwrap().area * references[i] * 0.05
        });
        let coating = core::array::from_fn::<_, 4, _>(|i| SurfactantFace {
            law: laws[i],
            amount_mol: original_amounts[i],
        });
        // Balance the actual capillary forces to isolate rigid translation.
        // These unequal concentrations have equal activity/occupancy.
        let external: Vec<_> =
            assemble_surfactant_step(&q, &q, &faces, SurfactantField::PerFace(&coating))
                .unwrap()
                .forces_n
                .into_iter()
                .map(|f| f.map(|x| -x))
                .collect();
        let mut state = model
            .initial_state_surfactant(
                0.0,
                &q,
                &velocity,
                SurfactantField::PerFace(&coating),
                inviscid,
                &external,
            )
            .unwrap();
        let mut amounts = original_amounts;
        let report = model
            .advance_partitioned_surfactant_diffusive_adaptive(
                &mut state,
                &mut amounts,
                AdaptiveSurfacePartitionedSurfactantTransport {
                    laws: &laws,
                    edges: &edges,
                    controls: adaptive_diffusion_controls(3),
                },
                inviscid,
                &external,
                0.02,
                controls(),
                || false,
            )
            .unwrap();
        for i in 0..4 {
            close(amounts[i], original_amounts[i], 3e-14);
            for d in 0..3 {
                close(state.q[3 * i + d], q[i][d] + 0.02 * translation[d], 5e-12);
            }
        }
        close(
            report.coupled.diffusion_free_energy_dissipation_j,
            0.0,
            3e-14,
        );
        assert_eq!(report.total_trials, 2);
        assert_eq!(report.total_evaluated_pair_exchanges, 30);
    }
}

#[test]
fn adaptive_empty_graph_recovers_native_mechanics_without_consuming_pair_or_trial_budget() {
    use super::surfactant_transport::AdaptiveSurfaceSurfactantTransport;
    let (q, faces, masses, velocity) = tetrahedron();
    let model = CapillarySurfaceDynamics::try_new(&masses, &faces, 4, 4).unwrap();
    let law = partition_laws(false)[0];
    let original_amounts = [0.02, 0.07, 0.03, 0.05];
    let coating = SurfactantField::Uniform {
        law,
        amounts_mol: &original_amounts,
    };
    let inviscid = SurfaceViscosityField::Uniform(SurfaceViscosity::ZERO);
    let external = [[0.0; 3]; 4];
    let mut reference = model
        .initial_state_surfactant(0.0, &q, &velocity, coating, inviscid, &external)
        .unwrap();
    let mut state = reference.clone();
    model
        .advance_surfactant(
            &mut reference,
            coating,
            inviscid,
            &external,
            0.005,
            controls(),
            || false,
        )
        .unwrap();
    let mut adaptive = adaptive_diffusion_controls(0);
    adaptive.maximum_trials = 1;
    adaptive.maximum_pair_exchanges = 0;
    let mut amounts = original_amounts;
    let report = model
        .advance_surfactant_diffusive_adaptive(
            &mut state,
            &mut amounts,
            AdaptiveSurfaceSurfactantTransport {
                law,
                edges: &[],
                controls: adaptive,
            },
            inviscid,
            &external,
            0.005,
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(state, reference);
    assert_eq!(amounts, original_amounts);
    assert_eq!(report.total_trials, 0);
    assert_eq!(report.total_evaluated_pair_exchanges, 0);
    assert_eq!(report.coupled.diffusion_free_energy_dissipation_j, 0.0);
}

#[test]
fn adaptive_work_is_shared_across_motion_and_late_refusals_preserve_both_states() {
    use super::super::surfactant::{
        SurfactantFace,
        diffusion::{SurfacePartitionDiffusionEdge, adaptive::AdaptiveSurfaceDiffusionError},
    };
    use super::surfactant_transport::{
        AdaptiveSurfacePartitionedSurfactantTransport, SurfaceDiffusionPhase,
    };
    let (q, faces, masses, velocity) = tetrahedron();
    let model = CapillarySurfaceDynamics::try_new(&masses, &faces, 4, 4).unwrap();
    let laws = partition_laws(false);
    let original_amounts = [0.02, 0.07, 0.03, 0.05];
    let coating = core::array::from_fn::<_, 4, _>(|i| SurfactantFace {
        law: laws[i],
        amount_mol: original_amounts[i],
    });
    let inviscid = SurfaceViscosityField::Uniform(SurfaceViscosity::ZERO);
    let external = [[0.0; 3]; 4];
    let mut initial = model
        .initial_state_surfactant(
            0.0,
            &q,
            &velocity,
            SurfactantField::PerFace(&coating),
            inviscid,
            &external,
        )
        .unwrap();
    model
        .advance_surfactant(
            &mut initial,
            SurfactantField::PerFace(&coating),
            inviscid,
            &external,
            0.005,
            controls(),
            || false,
        )
        .unwrap();
    assert!(!initial.history.is_empty());
    let edges = core::array::from_fn::<_, 3, _>(|i| SurfacePartitionDiffusionEdge {
        faces: [i, i + 1],
        conductance_mol_s: [30.0, 20.0, 40.0][i],
    });
    let adaptive = adaptive_diffusion_controls(3);
    let transport = AdaptiveSurfacePartitionedSurfactantTransport {
        laws: &laws,
        edges: &edges,
        controls: adaptive,
    };
    let mut state = initial.clone();
    let mut amounts = original_amounts;
    let mut calls = 0;
    let report = model
        .advance_partitioned_surfactant_diffusive_adaptive(
            &mut state,
            &mut amounts,
            transport,
            inviscid,
            &external,
            0.02,
            controls(),
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    assert!(report.diffusion_before.rejected_intervals > 0);
    assert!(report.diffusion_after.accepted_intervals > 1);
    for restrict_trials in [true, false] {
        let mut limited = adaptive;
        if restrict_trials {
            limited.maximum_trials = report.total_trials - 1;
        } else {
            limited.maximum_pair_exchanges = report.total_evaluated_pair_exchanges - 1;
        }
        let mut retry = initial.clone();
        let mut retry_amounts = original_amounts;
        let result = model.advance_partitioned_surfactant_diffusive_adaptive(
            &mut retry,
            &mut retry_amounts,
            AdaptiveSurfacePartitionedSurfactantTransport {
                controls: limited,
                ..transport
            },
            inviscid,
            &external,
            0.02,
            controls(),
            || false,
        );
        match result {
            Err(SurfaceDynamicsError::AdaptiveDiffusion {
                phase: SurfaceDiffusionPhase::AfterMotion,
                error: AdaptiveSurfaceDiffusionError::TrialBudget { completed_trials },
            }) if restrict_trials => assert_eq!(completed_trials, limited.maximum_trials),
            Err(SurfaceDynamicsError::AdaptiveDiffusion {
                phase: SurfaceDiffusionPhase::AfterMotion,
                error:
                    AdaptiveSurfaceDiffusionError::PairExchangeBudget {
                        completed,
                        required,
                    },
            }) if !restrict_trials => {
                assert!(completed >= report.diffusion_before.evaluated_pair_exchanges);
                assert!(completed + required > limited.maximum_pair_exchanges);
            }
            other => panic!("expected shared final-half budget refusal, got {other:?}"),
        }
        assert_eq!(retry, initial);
        assert_eq!(retry_amounts, original_amounts);
    }
    let mut retry = initial.clone();
    let mut retry_amounts = original_amounts;
    let mut current = 0;
    let result = model.advance_partitioned_surfactant_diffusive_adaptive(
        &mut retry,
        &mut retry_amounts,
        transport,
        inviscid,
        &external,
        0.02,
        controls(),
        || {
            current += 1;
            current == calls
        },
    );
    assert!(matches!(result, Err(SurfaceDynamicsError::Cancelled)));
    assert_eq!(retry, initial);
    assert_eq!(retry_amounts, original_amounts);
}

#[test]
fn moving_partitioned_coating_preserves_original_energy_amount_and_momentum() {
    use super::super::surfactant::{SurfactantFace, diffusion::SurfacePartitionDiffusionEdge};
    use super::surfactant_transport::SurfacePartitionedSurfactantTransport;
    let (q, faces, masses, v) = tetrahedron();
    let model = CapillarySurfaceDynamics::try_new(&masses, &faces, 4, 4).unwrap();
    let viscosity = SurfaceViscosityField::Uniform(SurfaceViscosity {
        shear_n_s_m: 0.1,
        dilatational_n_s_m: 0.15,
    });
    let external = [[0.0; 3]; 4];
    let edges = core::array::from_fn::<_, 3, _>(|i| SurfacePartitionDiffusionEdge {
        faces: [i, i + 1],
        conductance_mol_s: [0.4, 0.2, 0.3][i],
    });
    for langmuir in [false, true] {
        let laws = partition_laws(langmuir);
        let initial_amounts = [0.02, 0.07, 0.03, 0.05];
        let mut amounts = initial_amounts;
        let coating = core::array::from_fn::<_, 4, _>(|i| SurfactantFace {
            law: laws[i],
            amount_mol: amounts[i],
        });
        let mut state = model
            .initial_state_surfactant(
                0.0,
                &q,
                &v,
                SurfactantField::PerFace(&coating),
                viscosity,
                &external,
            )
            .unwrap();
        // Independently reconstruct physical storage using each ORIGINAL law
        // and actual endpoint area, rather than the composed account alone.
        let storage = |p: &[Vector], n: &[f64]| {
            (0..4)
                .map(|i| {
                    SurfactantFace {
                        law: laws[i],
                        amount_mol: n[i],
                    }
                    .at_area(geometry(faces[i].map(|j| p[j])).unwrap().area)
                    .unwrap()
                    .free_energy_j
                })
                .sum::<f64>()
        };
        let mut stored = storage(&q, &amounts);
        let initial_energy = kinetic(&masses, &state.v).unwrap() + stored;
        let (linear, angular) = momentum(&state, &masses);
        let mut viscous_loss = 0.0;
        let mut diffusion_loss = 0.0;
        for _ in 0..12 {
            let report = model
                .advance_partitioned_surfactant_diffusive(
                    &mut state,
                    &mut amounts,
                    SurfacePartitionedSurfactantTransport {
                        laws: &laws,
                        edges: &edges,
                        controls: diffusion_controls(3),
                    },
                    viscosity,
                    &external,
                    0.005,
                    controls(),
                    || false,
                )
                .unwrap();
            let after = storage(&positions(&state), &amounts);
            close(report.surface_free_energy_before_j, stored, 2e-14);
            close(report.surface_free_energy_after_j, after, 2e-14);
            close(report.surface_free_energy_change_j, after - stored, 2e-14);
            stored = after;
            viscous_loss += report.mechanics.surface_viscous_dissipation_j;
            diffusion_loss += report.diffusion_free_energy_dissipation_j;
            close(
                kinetic(&masses, &state.v).unwrap() + stored + viscous_loss + diffusion_loss,
                initial_energy,
                5e-11,
            );
            close(amounts.iter().sum(), initial_amounts.iter().sum(), 2e-14);
            assert!(amounts.iter().all(|n| *n > 0.0));
            assert_eq!(report.mechanics.surface.coefficient_supply_j, 0.0);
            let (p, l) = momentum(&state, &masses);
            for d in 0..3 {
                close(p[d], linear[d], 3e-11);
                close(l[d], angular[d], 3e-11);
            }
        }
        assert!(viscous_loss > 0.0 && diffusion_loss > 0.0);
        assert!(
            amounts
                .iter()
                .zip(initial_amounts)
                .any(|(&after, before)| (after - before).abs() > 1e-4)
        );
    }
}

#[test]
fn partitioned_uniform_law_recovers_composition_and_late_cancellation_is_atomic() {
    use super::super::surfactant::{
        SurfactantLaw,
        diffusion::{SurfaceDiffusionEdge, SurfacePartitionDiffusionEdge},
    };
    use super::surfactant_transport::{
        SurfacePartitionedSurfactantTransport, SurfaceSurfactantTransport,
    };
    let (q, faces, masses, v) = tetrahedron();
    let model = CapillarySurfaceDynamics::try_new(&masses, &faces, 4, 4).unwrap();
    let law = SurfactantLaw::IdealDilute {
        clean_tension_n_m: 0.3,
        thermal_energy_j_mol: 0.5,
        reference_concentration_mol_m2: 2.0,
    };
    let initial_amounts = [0.02, 0.07, 0.03, 0.05];
    let coating = SurfactantField::Uniform {
        law,
        amounts_mol: &initial_amounts,
    };
    let viscosity = SurfaceViscosityField::Uniform(SurfaceViscosity::ZERO);
    let external = [[0.0; 3]; 4];
    let mut initial = model
        .initial_state_surfactant(0.0, &q, &v, coating, viscosity, &external)
        .unwrap();
    model
        .advance_surfactant(
            &mut initial,
            coating,
            viscosity,
            &external,
            0.005,
            controls(),
            || false,
        )
        .unwrap();
    assert!(!initial.history.is_empty());
    let edges = core::array::from_fn::<_, 3, _>(|i| SurfaceDiffusionEdge {
        faces: [i, i + 1],
        conductance_m2_s: [0.4, 0.2, 0.3][i],
    });
    let partition_edges = edges.map(|edge| SurfacePartitionDiffusionEdge {
        faces: edge.faces,
        conductance_mol_s: 2.0 * edge.conductance_m2_s,
    });
    let laws = [law; 4];
    let transport = SurfacePartitionedSurfactantTransport {
        laws: &laws,
        edges: &partition_edges,
        controls: diffusion_controls(3),
    };
    let mut reference = initial.clone();
    let mut reference_amounts = initial_amounts;
    model
        .advance_surfactant_diffusive(
            &mut reference,
            &mut reference_amounts,
            SurfaceSurfactantTransport {
                law,
                edges: &edges,
                controls: diffusion_controls(3),
            },
            viscosity,
            &external,
            0.005,
            controls(),
            || false,
        )
        .unwrap();
    let mut composed = initial.clone();
    let mut amounts = initial_amounts;
    let mut calls = 0;
    model
        .advance_partitioned_surfactant_diffusive(
            &mut composed,
            &mut amounts,
            transport,
            viscosity,
            &external,
            0.005,
            controls(),
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    assert_eq!(reference, composed);
    assert_eq!(reference_amounts, amounts);

    let mut retry = initial.clone();
    let mut retry_amounts = initial_amounts;
    let mut current = 0;
    let result = model.advance_partitioned_surfactant_diffusive(
        &mut retry,
        &mut retry_amounts,
        transport,
        viscosity,
        &external,
        0.005,
        controls(),
        || {
            current += 1;
            current == calls
        },
    );
    assert!(matches!(result, Err(SurfaceDynamicsError::Cancelled)));
    assert_eq!(retry_amounts, initial_amounts);
    assert_eq!(retry, initial);
}

#[test]
fn empty_partition_graph_recovers_per_face_mechanics_and_bad_law_count_is_atomic() {
    use super::super::surfactant::{SurfactantFace, SurfactantLaw};
    use super::surfactant_transport::SurfacePartitionedSurfactantTransport;
    let (q, faces, masses, v) = tetrahedron();
    let model = CapillarySurfaceDynamics::try_new(&masses, &faces, 4, 4).unwrap();
    let mut laws = partition_laws(false);
    laws[1] = partition_laws(true)[1];
    laws[2] = SurfactantLaw::IdealDilute {
        clean_tension_n_m: 0.5,
        thermal_energy_j_mol: 0.7,
        reference_concentration_mol_m2: 1.2,
    };
    // Disconnected components may use different law families and temperatures.
    let initial_amounts = [0.02, 0.07, 0.03, 0.05];
    let coating = core::array::from_fn::<_, 4, _>(|i| SurfactantFace {
        law: laws[i],
        amount_mol: initial_amounts[i],
    });
    let viscosity = SurfaceViscosityField::Uniform(SurfaceViscosity::ZERO);
    let external = [[0.0; 3]; 4];
    let mut reference = model
        .initial_state_surfactant(
            0.0,
            &q,
            &v,
            SurfactantField::PerFace(&coating),
            viscosity,
            &external,
        )
        .unwrap();
    let initial = reference.clone();
    let mut composed = reference.clone();
    model
        .advance_surfactant(
            &mut reference,
            SurfactantField::PerFace(&coating),
            viscosity,
            &external,
            0.005,
            controls(),
            || false,
        )
        .unwrap();
    let mut amounts = initial_amounts;
    let report = model
        .advance_partitioned_surfactant_diffusive(
            &mut composed,
            &mut amounts,
            SurfacePartitionedSurfactantTransport {
                laws: &laws,
                edges: &[],
                controls: diffusion_controls(0),
            },
            viscosity,
            &external,
            0.005,
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(reference.q, composed.q);
    assert_eq!(reference.v, composed.v);
    assert_eq!(reference.a, composed.a);
    assert_eq!(amounts, initial_amounts);
    assert_eq!(report.diffusion_free_energy_dissipation_j, 0.0);
    let mut retry = initial.clone();
    let result = model.advance_partitioned_surfactant_diffusive(
        &mut retry,
        &mut amounts,
        SurfacePartitionedSurfactantTransport {
            laws: &laws[..3],
            edges: &[],
            controls: diffusion_controls(0),
        },
        viscosity,
        &external,
        0.005,
        controls(),
        || false,
    );
    assert!(matches!(
        result,
        Err(SurfaceDynamicsError::Diffusion(
            SurfaceDiffusionError::InvalidInput(_)
        ))
    ));
    assert_eq!(amounts, initial_amounts);
    assert_eq!(retry, initial);
}

#[test]
fn interfacial_dilatational_viscosity_matches_the_analytical_damped_breathing_mode() {
    let q = [
        [1.0, 0.0, 0.0],
        [-0.5, 3.0_f64.sqrt() / 2.0, 0.0],
        [-0.5, -3.0_f64.sqrt() / 2.0, 0.0],
    ];
    let masses = [0.8; 3];
    let faces = [[0, 1, 2]];
    let surface = CapillarySurfaceDynamics::try_new(&masses, &faces, 3, 1).unwrap();
    let sigma = 0.6;
    let kappa = 0.9;
    let tension = SurfaceTensionField::Uniform(sigma);
    // A uniform dilation has zero deviatoric strain; the shear coefficient
    // cannot affect this mode. The area strain produces c=sqrt(3)*kappa/m.
    let viscosity = SurfaceViscosityField::Uniform(SurfaceViscosity {
        shear_n_s_m: 3.7,
        dilatational_n_s_m: kappa,
    });
    let velocity = q.map(|p| p.map(|x| 0.2 * x));
    let force = [[0.0; 3]; 3];
    let mut state = surface
        .initial_state_viscous(0.0, &q, &velocity, tension, viscosity, &force)
        .unwrap();
    let h = 0.08;
    let k = 3.0_f64.sqrt() * sigma / (2.0 * masses[0]);
    let c = 3.0_f64.sqrt() * kappa / masses[0];
    let denominator = 1.0 + 0.5 * h * c + 0.25 * h * h * k;
    let radius = (1.0 + 0.5 * h * c - 0.25 * h * h * k + h * 0.2) / denominator;
    let speed = ((1.0 - 0.5 * h * c - 0.25 * h * h * k) * 0.2 - h * k) / denominator;
    let report = surface
        .advance_viscous(
            &mut state,
            tension,
            tension,
            viscosity,
            &force,
            h,
            controls(),
            || false,
        )
        .unwrap();
    for i in 0..3 {
        for d in 0..3 {
            close(state.q[3 * i + d], radius * q[i][d], 3e-12);
            close(state.v[3 * i + d], speed * q[i][d], 3e-12);
        }
    }
    let midpoint_speed = 0.5 * (0.2 + speed);
    close(
        report.surface_viscous_dissipation_j,
        3.0 * 3.0_f64.sqrt() * kappa * midpoint_speed * midpoint_speed * h,
        3e-13,
    );
    close(
        report.kinetic_change_j
            + report.surface.energy_change_j
            + report.surface_viscous_dissipation_j,
        0.0,
        3e-13,
    );
    assert!(report.surface_viscous_dissipation_j > 0.0);
}

#[test]
fn viscous_surface_dissipates_energy_while_preserving_linear_and_angular_momentum() {
    let (q, faces, masses, v) = tetrahedron();
    let surface = CapillarySurfaceDynamics::try_new(&masses, &faces, 4, 4).unwrap();
    let tension = SurfaceTensionField::Uniform(0.3);
    let coefficients = [
        SurfaceViscosity {
            shear_n_s_m: 0.2,
            dilatational_n_s_m: 0.4,
        },
        SurfaceViscosity {
            shear_n_s_m: 0.5,
            dilatational_n_s_m: 0.1,
        },
        SurfaceViscosity::ZERO,
        SurfaceViscosity {
            shear_n_s_m: 0.3,
            dilatational_n_s_m: 0.2,
        },
    ];
    let viscosity = SurfaceViscosityField::PerFace(&coefficients);
    let force = [[0.0; 3]; 4];
    let mut state = surface
        .initial_state_viscous(0.0, &q, &v, tension, viscosity, &force)
        .unwrap();
    let (linear, angular) = momentum(&state, &masses);
    let initial = kinetic(&masses, &state.v).unwrap()
        + assemble_surface_tension_step(&q, &q, &faces, tension, tension)
            .unwrap()
            .energy
            .energy_before_j;
    let mut losses = 0.0;
    for _ in 0..20 {
        let report = surface
            .advance_viscous(
                &mut state,
                tension,
                tension,
                viscosity,
                &force,
                0.005,
                controls(),
                || false,
            )
            .unwrap();
        assert!(report.surface_viscous_dissipation_j > 0.0);
        losses += report.surface_viscous_dissipation_j;
        close(
            report.kinetic_after_j + report.surface.energy_after_j + losses,
            initial,
            5e-11,
        );
        let (p, l) = momentum(&state, &masses);
        for d in 0..3 {
            close(p[d], linear[d], 3e-11);
            close(l[d], angular[d], 3e-11);
        }
    }
    assert!(losses > 1e-4);
}

#[test]
fn g1_equilateral_breathing_matches_midpoint_solution_and_second_order_time_refinement() {
    let q = [
        [1.0, 0.0, 0.0],
        [-0.5, 3.0_f64.sqrt() / 2.0, 0.0],
        [-0.5, -3.0_f64.sqrt() / 2.0, 0.0],
    ];
    let mass = [0.8; 3];
    let faces = [[0, 1, 2]];
    let sigma = 0.6;
    let field = SurfaceTensionField::Uniform(sigma);
    let surface = CapillarySurfaceDynamics::try_new(&mass, &faces, 3, 1).unwrap();
    let velocity = q.map(|p| p.map(|v| 0.2 * v));
    let force = [[0.0; 3]; 3];
    let initial = surface
        .initial_state(0.0, &q, &velocity, field, &force)
        .unwrap();
    // A=(3 sqrt(3)/4) r² and three masses give r''+omega² r=0.
    // The original structural midpoint step is its exact Cayley map.
    let square = 3.0_f64.sqrt() * sigma / (2.0 * mass[0]);
    let h = 0.08;
    let factor = 0.25 * h * h * square;
    let expected_r = ((1.0 - factor) + h * 0.2) / (1.0 + factor);
    let expected_v = ((1.0 - factor) * 0.2 - h * square) / (1.0 + factor);
    let mut state = initial.clone();
    let report = surface
        .advance(&mut state, field, field, &force, h, controls(), || false)
        .unwrap();
    for i in 0..3 {
        for d in 0..3 {
            close(state.q[3 * i + d], expected_r * q[i][d], 2e-13);
            close(state.v[3 * i + d], expected_v * q[i][d], 2e-12);
        }
    }
    close(
        report.kinetic_change_j + report.surface.energy_change_j,
        0.0,
        2e-13,
    );
    assert_eq!(state.steps, 1);
    assert_eq!(state.history.len(), 1);
    let time = 0.4;
    let omega = square.sqrt();
    let exact = (omega * time).cos() + 0.2 * (omega * time).sin() / omega;
    let errors: Vec<_> = [8, 16, 32]
        .map(|steps| {
            let mut state = initial.clone();
            for _ in 0..steps {
                surface
                    .advance(
                        &mut state,
                        field,
                        field,
                        &force,
                        time / steps as f64,
                        controls(),
                        || false,
                    )
                    .unwrap();
            }
            assert_eq!(state.history.len(), steps);
            (state.q[0] - exact).abs()
        })
        .into();
    for pair in errors.windows(2) {
        assert!(
            pair[0] / pair[1] > 3.9 && pair[0] / pair[1] < 4.1,
            "{errors:?}"
        );
    }
}

#[test]
fn g1_rotating_nonuniform_surface_preserves_actual_energy_linear_and_angular_momentum() {
    let (q, faces, masses, v) = tetrahedron();
    let surface = CapillarySurfaceDynamics::try_new(&masses, &faces, 4, 4).unwrap();
    let sigma = [0.2, 0.3, 0.5, 0.4];
    let field = SurfaceTensionField::PerFace(&sigma);
    let force = [[0.0; 3]; 4];
    let mut state = surface.initial_state(0.0, &q, &v, field, &force).unwrap();
    let (linear, angular) = momentum(&state, &masses);
    let initial_surface = assemble_surface_tension_step(&q, &q, &faces, field, field).unwrap();
    let initial_energy =
        kinetic(&masses, &state.v).unwrap() + initial_surface.energy.energy_before_j;
    for _ in 0..40 {
        let report = surface
            .advance(&mut state, field, field, &force, 0.005, controls(), || {
                false
            })
            .unwrap();
        close(
            report.kinetic_after_j + report.surface.energy_after_j,
            initial_energy,
            3e-11,
        );
        assert_eq!(report.surface.coefficient_supply_j, 0.0);
        assert_eq!(report.external_work_j, 0.0);
        let (p, l) = momentum(&state, &masses);
        for d in 0..3 {
            close(p[d], linear[d], 2e-11);
            close(l[d], angular[d], 3e-11);
        }
    }
    let edges = |p: &[Vector]| {
        (0..3)
            .map(|i| dot(sub(p[i + 1], p[0]), sub(p[i + 1], p[0])))
            .collect::<Vec<_>>()
    };
    let before = edges(&q);
    let after = edges(&positions(&state));
    assert!(
        before.iter().zip(after).any(|(a, b)| (a - b).abs() > 0.01),
        "the fixture must actually deform, not merely test rigid transport"
    );
}

#[test]
fn g1_prescribed_sigma_change_supplies_free_energy_and_held_force_equilibrium_is_stationary() {
    let q = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let faces = [[0, 1, 2]];
    let masses = [1.0; 3];
    let surface = CapillarySurfaceDynamics::try_new(&masses, &faces, 3, 1).unwrap();
    let zero = SurfaceTensionField::Uniform(0.0);
    let sigma = SurfaceTensionField::Uniform(0.6);
    let no_force = [[0.0; 3]; 3];
    let mut state = surface
        .initial_state(0.0, &q, &no_force, zero, &no_force)
        .unwrap();
    let report = surface
        .advance(&mut state, zero, sigma, &no_force, 0.1, controls(), || {
            false
        })
        .unwrap();
    assert!(report.kinetic_after_j > 0.0);
    assert!(report.surface.area_after_m2 < report.surface.area_before_m2);
    close(
        report.surface.coefficient_supply_j,
        0.3 * (report.surface.area_before_m2 + report.surface.area_after_m2),
        1e-14,
    );
    close(
        report.kinetic_change_j + report.surface.energy_change_j,
        report.surface.coefficient_supply_j,
        2e-12,
    );

    let initial_force = triangle_surface_tension(q, 0.6).unwrap().forces_n;
    let balanced = initial_force.map(|f| f.map(|v| -v));
    let mut equilibrium = surface
        .initial_state(0.0, &q, &no_force, sigma, &balanced)
        .unwrap();
    let initial = equilibrium.clone();
    let report = surface
        .advance(
            &mut equilibrium,
            sigma,
            sigma,
            &balanced,
            0.3,
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(equilibrium.q, initial.q);
    assert_eq!(equilibrium.v, initial.v);
    assert_eq!(report.external_work_j, 0.0);
    assert_eq!(report.kinetic_change_j, 0.0);

    // The zero-tension limit is actual uniformly accelerated nodal inertia.
    // Nonzero exterior work must enter once with the correct sign and units.
    let load = [[0.3, -0.2, 0.1]; 3];
    let mut ballistic = surface
        .initial_state(0.0, &q, &no_force, zero, &load)
        .unwrap();
    let h = 0.2;
    let report = surface
        .advance(&mut ballistic, zero, zero, &load, h, controls(), || false)
        .unwrap();
    for i in 0..3 {
        for d in 0..3 {
            close(
                ballistic.q[3 * i + d],
                q[i][d] + 0.5 * h * h * load[i][d] / masses[i],
                1e-14,
            );
            close(ballistic.v[3 * i + d], h * load[i][d] / masses[i], 1e-13);
        }
    }
    assert!(report.external_work_j > 0.0);
    close(report.kinetic_change_j, report.external_work_j, 2e-14);
    assert_eq!(report.surface.energy_change_j, 0.0);
}

#[test]
fn g3_surface_motion_is_objective_under_finite_rotation_translation_and_checkpoint_replay() {
    let (q, faces, masses, v) = tetrahedron();
    let surface = CapillarySurfaceDynamics::try_new(&masses, &faces, 4, 4).unwrap();
    let rotate = |p: Vector| [p[2], p[1], -p[0]];
    let transformed = q.map(|p| {
        let p = rotate(p);
        [p[0] + 2.0, p[1] - 1.0, p[2] + 3.0]
    });
    let sigma0 = [0.2, 0.3, 0.5, 0.4];
    let sigma1 = [0.25, 0.28, 0.6, 0.35];
    let a = SurfaceTensionField::PerFace(&sigma0);
    let b = SurfaceTensionField::PerFace(&sigma1);
    let force = [[0.0; 3]; 4];
    let mut original = surface.initial_state(0.0, &q, &v, a, &force).unwrap();
    let mut replay = original.clone();
    let mut rotated = surface
        .initial_state(0.0, &transformed, &v.map(rotate), a, &force)
        .unwrap();
    let report = surface
        .advance(&mut original, a, b, &force, 0.03, controls(), || false)
        .unwrap();
    let other = surface
        .advance(&mut rotated, a, b, &force, 0.03, controls(), || false)
        .unwrap();
    surface
        .advance(&mut replay, a, b, &force, 0.03, controls(), || false)
        .unwrap();
    assert_eq!(original, replay);
    for i in 0..4 {
        let p = rotate([
            original.q[3 * i],
            original.q[3 * i + 1],
            original.q[3 * i + 2],
        ]);
        let v = rotate([
            original.v[3 * i],
            original.v[3 * i + 1],
            original.v[3 * i + 2],
        ]);
        for d in 0..3 {
            close(rotated.q[3 * i + d], p[d] + [2.0, -1.0, 3.0][d], 2e-12);
            close(rotated.v[3 * i + d], v[d], 2e-11);
        }
    }
    close(
        report.surface.energy_after_j,
        other.surface.energy_after_j,
        2e-12,
    );
    close(
        report.surface.coefficient_supply_j,
        other.surface.coefficient_supply_j,
        2e-13,
    );
    close(report.energy_residual_j, 0.0, 2e-12);
}

#[test]
fn g4_cancelled_invalid_and_unconverged_steps_preserve_native_history_and_every_state_field() {
    let (q, faces, masses, v) = tetrahedron();
    let surface = CapillarySurfaceDynamics::try_new(&masses, &faces, 4, 4).unwrap();
    let field = SurfaceTensionField::Uniform(0.4);
    let force = [[0.0; 3]; 4];
    let mut state = surface.initial_state(0.0, &q, &v, field, &force).unwrap();
    surface
        .advance(&mut state, field, field, &force, 0.01, controls(), || false)
        .unwrap();
    let checkpoint = state.clone();
    let mut calls = 0;
    let error = surface
        .advance(&mut state, field, field, &force, 0.01, controls(), || {
            calls += 1;
            calls >= 9
        })
        .unwrap_err();
    assert!(matches!(error, SurfaceDynamicsError::Cancelled));
    assert_eq!(state, checkpoint);
    let bad = SurfaceTensionField::Uniform(-0.1);
    assert!(matches!(
        surface.advance(&mut state, field, bad, &force, 0.01, controls(), || false),
        Err(SurfaceDynamicsError::Surface(
            CapillarySurfaceError::InvalidSurfaceTension
        ))
    ));
    assert_eq!(state, checkpoint);
    let mut budget = controls();
    budget.solve.max_newton_iterations = 1;
    assert!(matches!(
        surface.advance(&mut state, field, field, &force, 0.3, budget, || false),
        Err(SurfaceDynamicsError::Time(_))
    ));
    assert_eq!(state, checkpoint);
    let mut exact_gate = controls();
    exact_gate.impulse_tolerance_n_s = 0.0;
    exact_gate.kinematic_tolerance_m = 0.0;
    exact_gate.energy_tolerance_j = 0.0;
    assert!(matches!(
        surface.advance(&mut state, field, field, &force, 0.01, exact_gate, || false),
        Err(SurfaceDynamicsError::Balance { .. })
    ));
    assert_eq!(state, checkpoint);
    assert!(CapillarySurfaceDynamics::try_new(&[1.0, 0.0, 1.0], &[[0, 1, 2]], 3, 1).is_err());
    assert!(CapillarySurfaceDynamics::try_new(&masses, &faces, 3, 4).is_err());
}
