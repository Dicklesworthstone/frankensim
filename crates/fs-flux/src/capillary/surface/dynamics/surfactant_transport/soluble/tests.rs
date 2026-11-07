use super::*;
use fs_material::osmosis::nonideal::VirialSolutionModel;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use std::sync::Arc;

const RT: f64 = R * 300.0;
const SATURATION: f64 = 1e-5;
const VOLUME: f64 = 1e-6;
const REFERENCE: f64 = 2.0;

fn law(saturation: f64) -> SurfactantLaw {
    SurfactantLaw::Langmuir {
        clean_tension_n_m: 0.2,
        thermal_energy_j_mol: RT,
        saturation_concentration_mol_m2: saturation,
    }
}
fn bulk() -> AqueousSolution {
    AqueousSolution::nonideal(
        Arc::new(
            VirialSolutionModel::new(300.0, 18e-6, REFERENCE, vec![20.0; 2], vec![0.0; 4]).unwrap(),
        ),
        VOLUME,
        vec![4e-6, 1e-6],
    )
    .unwrap()
}
fn chemistry(step: f64, faces: usize, edges: usize) -> SolubleSurfaceControls {
    SolubleSurfaceControls {
        maximum_step_s: step,
        maximum_substeps: 100,
        maximum_root_evaluations: 20_000,
        maximum_surface_pair_exchanges: 1000,
        maximum_faces: faces,
        maximum_edges: edges,
        root_tolerance_mol: 1e-20,
        amount_tolerance_mol: 1e-15,
        energy_tolerance_j: 1e-11,
    }
}
fn mechanical() -> SurfaceDynamicsControls {
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
fn inviscid() -> SurfaceViscosityField<'static> {
    SurfaceViscosityField::Uniform(SurfaceViscosity::ZERO)
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; allowance {tolerance:e}"
    );
}
fn same_state(a: &SecondOrderState, b: &SecondOrderState) {
    assert_eq!(a.q, b.q);
    assert_eq!(a.v, b.v);
    assert_eq!(a.a, b.a);
    assert_eq!(a.t, b.t);
    assert_eq!(a.steps, b.steps);
    assert_eq!(a.history, b.history);
}
fn positions(state: &SecondOrderState) -> Vec<Vector> {
    state
        .q
        .chunks_exact(3)
        .map(|q| [q[0], q[1], q[2]])
        .collect()
}

// Independent lattice free energy and finite-bath implicit quadratic.
fn surface_energy(area: f64, amount: f64) -> f64 {
    let capacity = area * SATURATION;
    let theta = amount / capacity;
    0.2 * area + RT * (amount * theta.ln() + (capacity - amount) * (-theta).ln_1p())
}
fn adsorption(area: f64, bulk: f64, original: f64, dt: f64, kd: f64, offset: f64) -> f64 {
    let capacity = area * SATURATION;
    let total = bulk + original;
    let a = kd * (-offset / RT).exp() / (VOLUME * REFERENCE);
    let b = 1.0 + dt * kd + dt * a * (total + capacity);
    let c = original + dt * a * total * capacity;
    2.0 * c / (b + (b * b - 4.0 * dt * a * c).sqrt())
}

#[test]
fn moving_finite_bath_matches_independent_exchange_and_breathing_equations() {
    let q = [
        [1.0, 0.0, 0.0],
        [-0.5, 3.0_f64.sqrt() / 2.0, 0.0],
        [-0.5, -3.0_f64.sqrt() / 2.0, 0.0],
    ];
    let masses = [0.8; 3];
    let faces = [[0, 1, 2]];
    let model = CapillarySurfaceDynamics::try_new(&masses, &faces, 3, 1).unwrap();
    let area = 3.0 * 3.0_f64.sqrt() / 4.0;
    let dt = 0.16;
    let offset = RT * 2.0_f64.ln();
    let kinetics = SolubleSurfaceKinetics {
        bulk_species: 0,
        standard_potential_offset_j_mol: offset,
        desorption_rates_s_inverse: &[3.0],
    };
    let laws = [law(SATURATION)];
    let transport = SolubleSurfaceTransport {
        laws: &laws,
        kinetics,
        edges: &[],
        controls: chemistry(dt, 1, 0),
    };
    for speed in [-0.3, 0.3] {
        let mut bath = bulk();
        let old_bulk = bath.solute_moles()[0];
        let mut amounts = [1e-6];
        let initial = amounts[0];
        let midpoint = adsorption(area, old_bulk, initial, 0.5 * dt, 3.0, offset);
        let middle_bulk = old_bulk + initial - midpoint;
        let residual = |radius: f64| {
            let final_area = area * radius * radius;
            let tension = (surface_energy(final_area, midpoint) - surface_energy(area, midpoint))
                / (final_area - area);
            2.0 * masses[0] * (radius - 1.0 - dt * speed) / (dt * dt)
                + tension * 3.0_f64.sqrt() * (1.0 + radius) / 4.0
        };
        let (mut low, mut high) = if speed > 0.0 {
            (1.001, 1.2)
        } else {
            (0.8, 0.999)
        };
        assert!(residual(low) < 0.0 && residual(high) > 0.0);
        for _ in 0..60 {
            let middle = 0.5 * (low + high);
            if residual(middle) > 0.0 {
                high = middle;
            } else {
                low = middle;
            }
        }
        let radius = 0.5 * (low + high);
        let after_area = area * radius * radius;
        let expected = adsorption(after_area, middle_bulk, midpoint, 0.5 * dt, 3.0, offset);
        let expected_bulk = middle_bulk + midpoint - expected;
        let mut state = model
            .initial_state_surfactant(
                0.0,
                &q,
                &q.map(|v| v.map(|x| speed * x)),
                SurfactantField::Uniform {
                    law: laws[0],
                    amounts_mol: &amounts,
                },
                inviscid(),
                &[[0.0; 3]; 3],
            )
            .unwrap();
        let report = model
            .advance_soluble_surfactant(
                &mut state,
                &mut bath,
                &mut amounts,
                transport,
                inviscid(),
                &[[0.0; 3]; 3],
                dt,
                mechanical(),
                || false,
            )
            .unwrap();
        for i in 0..3 {
            for d in 0..3 {
                close(state.q[3 * i + d], radius * q[i][d], 5e-12);
                close(
                    state.v[3 * i + d],
                    (2.0 * (radius - 1.0) / dt - speed) * q[i][d],
                    5e-11,
                );
            }
        }
        close(amounts[0], expected, 2e-19);
        close(bath.solute_moles()[0], expected_bulk, 2e-19);
        close(
            report.before_motion.surface_amount_after_mol,
            midpoint,
            2e-19,
        );
        assert!(
            (expected - adsorption(area, middle_bulk, midpoint, 0.5 * dt, 3.0, offset)).abs()
                > 1e-9
        );
        let mixing = |n: f64| RT * n * ((n / (VOLUME * REFERENCE)).ln() - 1.0);
        let expected_energy = mixing(expected_bulk) - mixing(old_bulk)
            + surface_energy(after_area, expected)
            - surface_energy(area, initial)
            + offset * (expected - initial);
        close(report.free_energy_change_j, expected_energy, 2e-13);
        close(
            report.mechanics.surface.energy_change_j,
            surface_energy(after_area, midpoint) - surface_energy(area, midpoint),
            2e-13,
        );
        close(report.amount_defect_mol, 0.0, 1e-18);
        close(report.energy_account_residual_j, 0.0, 2e-13);
        close(report.energy_residual_j, 0.0, 2e-12);
        assert!(report.exchange_free_energy_dissipation_j > 0.0);
        assert_eq!(report.total_substeps, 2);
        assert_eq!(report.total_surface_pair_exchanges, 0);
        assert_eq!(bath.solvent_volume_m3(), VOLUME);
        assert_eq!(bath.solute_moles()[1], 1e-6);
    }
}

fn patch() -> ([Vector; 4], [[usize; 3]; 2], [f64; 4], [Vector; 4]) {
    let q = [
        [0.0, 0.0, 0.0],
        [1.2, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.2, 1.0, 0.2],
    ];
    let v = q.map(|p| [0.2 * p[0], 0.1 * p[1], 0.05 * p[2]]);
    (q, [[0, 1, 2], [1, 3, 2]], [0.7, 0.8, 0.9, 1.0], v)
}

#[test]
fn heterogeneous_moving_sites_reuse_native_geometry_diffusion_and_mechanical_loss() {
    let (q, faces, masses, v) = patch();
    let model = CapillarySurfaceDynamics::try_new(&masses, &faces, 4, 2).unwrap();
    let laws = [law(SATURATION), law(2.0 * SATURATION)];
    let edges = [SurfacePartitionDiffusionEdge {
        faces: [0, 1],
        conductance_mol_s: 3e-5,
    }];
    let kinetics = SolubleSurfaceKinetics {
        bulk_species: 0,
        standard_potential_offset_j_mol: RT * 0.7,
        desorption_rates_s_inverse: &[4.0, 0.0],
    };
    let transport = SolubleSurfaceTransport {
        laws: &laws,
        kinetics,
        edges: &edges,
        controls: chemistry(0.02, 2, 1),
    };
    let viscosity = SurfaceViscosityField::Uniform(SurfaceViscosity {
        shear_n_s_m: 0.01,
        dilatational_n_s_m: 0.02,
    });
    let mut amounts = [2e-7, 2e-6];
    let coating = [
        SurfactantFace {
            law: laws[0],
            amount_mol: amounts[0],
        },
        SurfactantFace {
            law: laws[1],
            amount_mol: amounts[1],
        },
    ];
    let original = model
        .initial_state_surfactant(
            0.0,
            &q,
            &v,
            SurfactantField::PerFace(&coating),
            viscosity,
            &[[0.0; 3]; 4],
        )
        .unwrap();
    let mut state = original.clone();
    let mut reference = original;
    let mut bath = bulk();
    let mut reference_bath = bath.clone();
    let mut reference_amounts = amounts;
    let dt = 0.08;
    let before_areas = faces.map(|f| geometry(f.map(|i| q[i])).unwrap().area);
    advance_soluble_surface(
        &mut reference_bath,
        &laws,
        &before_areas,
        &mut reference_amounts,
        kinetics,
        &edges,
        0.5 * dt,
        transport.controls,
        || false,
    )
    .unwrap();
    let middle = std::array::from_fn::<_, 2, _>(|i| SurfactantFace {
        law: laws[i],
        amount_mol: reference_amounts[i],
    });
    model
        .advance_surfactant(
            &mut reference,
            SurfactantField::PerFace(&middle),
            viscosity,
            &[[0.0; 3]; 4],
            dt,
            mechanical(),
            || false,
        )
        .unwrap();
    let moved = positions(&reference);
    let after_areas = faces.map(|f| geometry(f.map(|i| moved[i])).unwrap().area);
    advance_soluble_surface(
        &mut reference_bath,
        &laws,
        &after_areas,
        &mut reference_amounts,
        kinetics,
        &edges,
        0.5 * dt,
        transport.controls,
        || false,
    )
    .unwrap();
    let report = model
        .advance_soluble_surfactant(
            &mut state,
            &mut bath,
            &mut amounts,
            transport,
            viscosity,
            &[[0.0; 3]; 4],
            dt,
            mechanical(),
            || false,
        )
        .unwrap();
    same_state(&state, &reference);
    assert_eq!(amounts, reference_amounts);
    assert_eq!(bath, reference_bath);
    assert!(std::ptr::eq(
        bath.nonideal_model().unwrap(),
        reference_bath.nonideal_model().unwrap()
    ));
    assert!(report.surface_diffusion_free_energy_dissipation_j > 0.0);
    assert!(report.mechanics.surface_viscous_dissipation_j > 0.0);
    assert!(report.exchange_free_energy_dissipation_j > 0.0);
    close(
        report.total_amount_before_mol,
        report.total_amount_after_mol,
        1e-18,
    );
    close(report.energy_residual_j, 0.0, 2e-12);
    assert_eq!(report.total_surface_pair_exchanges, 8);
    assert_ne!(before_areas, after_areas);
}

#[test]
fn sealed_soluble_surface_retains_original_native_motion_and_bulk_inventory() {
    let (q, faces, masses, v) = patch();
    let model = CapillarySurfaceDynamics::try_new(&masses, &faces, 4, 2).unwrap();
    let laws = [law(SATURATION); 2];
    let mut amounts = [1e-6, 2e-6];
    let initial = model
        .initial_state_surfactant(
            0.0,
            &q,
            &v,
            SurfactantField::Uniform {
                law: laws[0],
                amounts_mol: &amounts,
            },
            inviscid(),
            &[[0.0; 3]; 4],
        )
        .unwrap();
    let mut state = initial.clone();
    let mut expected = initial;
    let mut bath = bulk();
    let saved = bath.clone();
    let transport = SolubleSurfaceTransport {
        laws: &laws,
        kinetics: SolubleSurfaceKinetics {
            bulk_species: 0,
            standard_potential_offset_j_mol: RT * 3.0,
            desorption_rates_s_inverse: &[0.0; 2],
        },
        edges: &[],
        controls: chemistry(0.1, 2, 0),
    };
    model
        .advance_surfactant(
            &mut expected,
            SurfactantField::Uniform {
                law: laws[0],
                amounts_mol: &amounts,
            },
            inviscid(),
            &[[0.0; 3]; 4],
            0.1,
            mechanical(),
            || false,
        )
        .unwrap();
    let report = model
        .advance_soluble_surfactant(
            &mut state,
            &mut bath,
            &mut amounts,
            transport,
            inviscid(),
            &[[0.0; 3]; 4],
            0.1,
            mechanical(),
            || false,
        )
        .unwrap();
    same_state(&state, &expected);
    assert_eq!(bath, saved);
    assert_eq!(amounts, [1e-6, 2e-6]);
    assert_eq!(report.chemical_free_energy_change_j, 0.0);
    assert_eq!(report.exchange_free_energy_dissipation_j, 0.0);
    assert_eq!(report.surface_diffusion_free_energy_dissipation_j, 0.0);
}

#[test]
fn final_stage_work_refusal_and_final_cancellation_preserve_all_three_states() {
    let (q, faces, masses, v) = patch();
    let model = CapillarySurfaceDynamics::try_new(&masses, &faces, 4, 2).unwrap();
    let laws = [law(SATURATION); 2];
    let initial_amounts = [1e-6, 2e-6];
    let initial_bulk = bulk();
    let initial = model
        .initial_state_surfactant(
            0.0,
            &q,
            &v,
            SurfactantField::Uniform {
                law: laws[0],
                amounts_mol: &initial_amounts,
            },
            inviscid(),
            &[[0.0; 3]; 4],
        )
        .unwrap();
    let edges = [SurfacePartitionDiffusionEdge {
        faces: [0, 1],
        conductance_mol_s: 2e-5,
    }];
    let transport = SolubleSurfaceTransport {
        laws: &laws,
        kinetics: SolubleSurfaceKinetics {
            bulk_species: 0,
            standard_potential_offset_j_mol: RT * 0.7,
            desorption_rates_s_inverse: &[2.0, 1.0],
        },
        edges: &edges,
        controls: chemistry(0.1, 2, 1),
    };
    let mut calls = 0;
    let baseline = model
        .advance_soluble_surfactant(
            &mut initial.clone(),
            &mut initial_bulk.clone(),
            &mut initial_amounts.clone(),
            transport,
            inviscid(),
            &[[0.0; 3]; 4],
            0.1,
            mechanical(),
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    for budget in 0..4 {
        let mut limited = transport;
        match budget {
            0 => limited.controls.maximum_substeps = baseline.before_motion.substeps,
            1 => {
                limited.controls.maximum_root_evaluations = baseline.before_motion.root_evaluations
            }
            2 => {
                limited.controls.maximum_surface_pair_exchanges =
                    baseline.before_motion.surface_pair_exchanges
            }
            _ => {}
        }
        let mut state = initial.clone();
        let mut bath = initial_bulk.clone();
        let mut amounts = initial_amounts;
        let mut current = 0;
        let error = model
            .advance_soluble_surfactant(
                &mut state,
                &mut bath,
                &mut amounts,
                limited,
                inviscid(),
                &[[0.0; 3]; 4],
                0.1,
                mechanical(),
                || {
                    current += 1;
                    budget == 3 && current == calls
                },
            )
            .unwrap_err();
        if budget == 3 {
            assert!(matches!(
                error,
                SolubleSurfaceDynamicsError::Mechanics(SurfaceDynamicsError::Cancelled)
            ));
        } else {
            assert!(matches!(
                error,
                SolubleSurfaceDynamicsError::Transport {
                    phase: SurfaceDiffusionPhase::AfterMotion,
                    ..
                }
            ));
        }
        same_state(&state, &initial);
        assert_eq!(bath, initial_bulk);
        assert_eq!(amounts, initial_amounts);
    }
}
