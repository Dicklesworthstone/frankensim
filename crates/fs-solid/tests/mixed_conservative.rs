//! Exact P2/P1 volume constraints with conservative finite-strain motion.
//! Physical storage, weak volume, and signed-mass momenta are evaluated without
//! the assembled residual or material discrete-gradient implementation.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_feec::highorder::simplex::duffy_quadrature;
use fs_material::hyper::discrete_gradient::DiscreteGradientSettings;
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_solid::hyper3::mixed::conservative::{
    ConservativeFormulation, ConservativeTetDefinition, ConservativeTetSolid,
};
use fs_solid::hyper3::mixed::{MixedLawField, MixedUpSettings};
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::linear3::tet10::{TET10_EDGES, promote_p1_to_p2};
use fs_solver::nonlinear::NewtonKrylovConfig;

const VOLUME: f64 = 1.0 / 6.0;
const DENSITY: f64 = 6.0;
const WORK_ALLOWANCE: f64 = 2e-9;

fn with_cx(run: impl FnOnce(&Cx<'_>)) {
    let cancel = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            &cancel,
            arena,
            StreamKey {
                seed: 62,
                kernel_id: 23,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

fn mesh(cx: &Cx<'_>) -> (Vec<[f64; 3]>, Vec<[usize; 10]>) {
    promote_p1_to_p2(
        &[[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        &[[0, 1, 2, 3]],
        None,
        cx,
    )
    .unwrap()
}

fn law(kappa: f64) -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::MooneyRivlin {
            c10: 4.0,
            c01: 2.0,
            kappa,
        },
        4.0,
    )
    .unwrap()
}

fn definition<'a>(
    nodes_m: &'a [[f64; 3]],
    elements: &'a [[usize; 10]],
    material: &'a Hyperelastic,
    densities: &'a [f64],
) -> ConservativeTetDefinition<'a> {
    ConservativeTetDefinition {
        nodes_m,
        elements,
        material: MixedLawField::Uniform(material),
        reference_density_kg_m3: densities,
        formulation: ConservativeFormulation::ExactIncompressible {
            maximum_relative_volume_residual: 1e-9,
            maximum_initial_relative_volume_rate_per_s: 1e-10,
            maximum_point_volume_work_defect: 1e-12,
        },
        budget: TetAssemblyBudget::standard(),
        quadrature_points_per_dir: 4,
        maximum_quadrature_points: 64,
        maximum_mass_assembly_bytes: 128 * 1024,
        material_settings: DiscreteGradientSettings {
            minimum_relative_strain_increment: 1e-6,
            absolute_work_tolerance_pa: 1e-10,
            relative_work_tolerance: 1e-12,
        },
    }
}

fn settings() -> MixedUpSettings {
    MixedUpSettings {
        load_steps: 1,
        maximum_newton_iterations: 30,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 2e-11,
            relative_tolerance: 0.0,
            ..NewtonKrylovConfig::default()
        },
        length_scale_m: 1.0,
        stress_scale_pa: 12.0,
        maximum_assembly_bytes: 2 * 1024 * 1024,
        maximum_factor_bytes: 2 * 1024 * 1024,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn initial_state(nodes: &[[f64; 3]]) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
    // This nonaffine triangular map has J=1 pointwise. With current position x,
    // v=t+(Omega+D)x has div_x(v)=tr(D)=0, including its deforming component.
    let displacement: Vec<_> = nodes
        .iter()
        .map(|x| [0.25 * x[1] + 0.4 * x[1] * x[2], 0.0, 0.0])
        .collect();
    let velocity = nodes
        .iter()
        .zip(&displacement)
        .map(|(x, u)| {
            let position = core::array::from_fn(|c| x[c] + u[c] - 0.25);
            let spin = cross([0.3, -0.2, 0.5], position);
            [
                0.2 + spin[0] + 0.08 * position[0],
                -0.15 + spin[1] - 0.06 * position[1],
                0.1 + spin[2] - 0.02 * position[2],
            ]
        })
        .collect();
    (displacement, velocity)
}

// Exact affine P2 mass / (rho*V), integrated from barycentric monomials.
fn mass_fraction(a: usize, b: usize) -> f64 {
    let numerator = match (a < 4, b < 4) {
        (true, true) => {
            if a == b {
                6.0
            } else {
                1.0
            }
        }
        (true, false) => {
            if TET10_EDGES[b - 4].contains(&a) {
                -4.0
            } else {
                -6.0
            }
        }
        (false, true) => return mass_fraction(b, a),
        (false, false) => {
            if a == b {
                32.0
            } else if TET10_EDGES[a - 4]
                .iter()
                .any(|i| TET10_EDGES[b - 4].contains(i))
            {
                16.0
            } else {
                8.0
            }
        }
    };
    numerator / 420.0
}

fn independent_motion(
    nodes: &[[f64; 3]],
    element: &[usize; 10],
    displacement: &[[f64; 3]],
    velocity: &[[f64; 3]],
) -> (f64, [f64; 3], [f64; 3]) {
    let mut momentum = vec![[0.0; 3]; nodes.len()];
    for (a, &i) in element.iter().enumerate() {
        for (b, &j) in element.iter().enumerate() {
            for c in 0..3 {
                momentum[i][c] += DENSITY * VOLUME * mass_fraction(a, b) * velocity[j][c];
            }
        }
    }
    let (mut kinetic, mut linear, mut angular) = (0.0, [0.0; 3], [0.0; 3]);
    for i in 0..nodes.len() {
        let position = core::array::from_fn(|c| nodes[i][c] + displacement[i][c]);
        let spin = cross(position, momentum[i]);
        for c in 0..3 {
            kinetic += 0.5 * velocity[i][c] * momentum[i][c];
            linear[c] += momentum[i][c];
            angular[c] += spin[c];
        }
    }
    (kinetic, linear, angular)
}

fn independent_fields(element: &[usize; 10], u: &[[f64; 3]]) -> (f64, [f64; 4]) {
    let (mut stored, mut volume) = (0.0, [0.0; 4]);
    for (lambda, weight) in duffy_quadrature(4) {
        // This fixture's reference tetrahedron is the unit simplex, so its
        // physical barycentric derivatives and geometry determinant are exact.
        let dl = [[-1.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let mut gradients = [[0.0; 3]; 10];
        for a in 0..4 {
            gradients[a] = dl[a].map(|v| (4.0 * lambda[a] - 1.0) * v);
        }
        for (e, &[a, b]) in TET10_EDGES.iter().enumerate() {
            gradients[4 + e] =
                core::array::from_fn(|j| 4.0 * (lambda[a] * dl[b][j] + lambda[b] * dl[a][j]));
        }
        let mut f = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        for a in 0..10 {
            for i in 0..3 {
                for j in 0..3 {
                    f[3 * i + j] += u[element[a]][i] * gradients[a][j];
                }
            }
        }
        let cof = [
            f[4] * f[8] - f[5] * f[7],
            f[5] * f[6] - f[3] * f[8],
            f[3] * f[7] - f[4] * f[6],
            f[2] * f[7] - f[1] * f[8],
            f[0] * f[8] - f[2] * f[6],
            f[1] * f[6] - f[0] * f[7],
            f[1] * f[5] - f[2] * f[4],
            f[2] * f[3] - f[0] * f[5],
            f[0] * f[4] - f[1] * f[3],
        ];
        let j = f[0] * cof[0] + f[1] * cof[1] + f[2] * cof[2];
        assert!(j > 0.0);
        let i1: f64 = f.iter().map(|v| v * v).sum();
        let i2: f64 = cof.iter().map(|v| v * v).sum();
        // Physical Mooney-Rivlin remainder, without calling either full energy
        // (which contains kappa) or the new two-state material implementation.
        let iso = j.powf(-2.0 / 3.0);
        stored += weight * (4.0 * (iso * i1 - 3.0) + 2.0 * (iso * iso * i2 - 3.0));
        for a in 0..4 {
            volume[a] += weight * lambda[a] * (j - 1.0);
        }
    }
    (stored, volume)
}

fn state_bits(solid: &ConservativeTetSolid<'_>) -> Vec<u64> {
    solid
        .displacement_m()
        .iter()
        .flatten()
        .chain(solid.velocity_mps().iter().flatten())
        .chain(solid.algorithmic_pressure_pa())
        .copied()
        .chain([
            solid.stored_energy_j(),
            solid.kinetic_energy_j(),
            solid.total_energy_j(),
            solid.time_s(),
        ])
        .map(f64::to_bits)
        .collect()
}

#[test]
fn g1_free_exact_body_preserves_weak_volume_energy_and_both_momenta() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(60.0);
        let densities = [DENSITY];
        let (u0, v0) = initial_state(&nodes);
        let mut solid = definition(&nodes, &elements, &material, &densities)
            .prepare(&u0, &v0, cx)
            .unwrap();
        assert_eq!(solid.pressure_nodes(), &[0, 1, 2, 3]);
        let (kinetic0, linear0, angular0) = independent_motion(&nodes, &elements[0], &u0, &v0);
        let (stored0, mut previous_volume) = independent_fields(&elements[0], &u0);
        let h0 = kinetic0 + stored0;
        for volume in previous_volume {
            near(volume, 0.0, 2e-16);
        }
        let length2 = |u: &[[f64; 3]]| {
            (0..3)
                .map(|c| (nodes[1][c] + u[1][c] - nodes[0][c] - u[0][c]).powi(2))
                .sum::<f64>()
        };
        let initial_length2 = length2(&u0);
        let (mut worst_energy, mut worst_linear, mut worst_angular, mut worst_volume, mut strain) =
            (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
        for _ in 0..32 {
            let step = solid
                .advance(0.01, &[], &[], [0.0; 3], settings(), WORK_ALLOWANCE, cx)
                .unwrap();
            assert!(step.reactions_n.is_empty());
            assert_eq!(step.inertia.numerical_loss_j, 0.0);
            let next = &step.solid;
            let (stored, volume) = independent_fields(&elements[0], next.displacement_m());
            let (kinetic, linear, angular) = independent_motion(
                &nodes,
                &elements[0],
                next.displacement_m(),
                next.velocity_mps(),
            );
            let pressure_work: f64 = next
                .algorithmic_pressure_pa()
                .iter()
                .enumerate()
                .map(|(a, p)| p * (volume[a] - previous_volume[a]))
                .sum();
            near(next.stored_energy_j(), stored, 2e-11);
            near(step.inertia.kinetic_after_j, kinetic, 2e-11);
            near(step.energy.pressure_work_j, pressure_work, 2e-11);
            near(step.energy.pressure_work_j, 0.0, WORK_ALLOWANCE);
            near(
                step.energy.material_work_j,
                stored - solid.stored_energy_j(),
                WORK_ALLOWANCE,
            );
            near(step.energy.balance_error_j, 0.0, WORK_ALLOWANCE);
            near(step.energy.cumulative_balance_error_j, 0.0, WORK_ALLOWANCE);
            for a in 0..4 {
                near(step.pressure_residual_m3[a], volume[a], 2e-15);
                worst_volume = worst_volume.max(volume[a].abs() / VOLUME);
            }
            worst_energy = worst_energy.max((stored + kinetic - h0).abs());
            for c in 0..3 {
                worst_linear = worst_linear.max((linear[c] - linear0[c]).abs());
                worst_angular = worst_angular.max((angular[c] - angular0[c]).abs());
            }
            strain = strain.max((length2(next.displacement_m()) - initial_length2).abs());
            previous_volume = volume;
            solid = step.solid;
        }
        assert!(strain > 1e-3, "the free body must change strain");
        near(worst_volume, 0.0, 1e-9);
        near(worst_energy, 0.0, WORK_ALLOWANCE);
        near(worst_linear, 0.0, 1e-10);
        near(worst_angular, 0.0, 1e-10);
        println!(
            "exact free tet: energy_drift={worst_energy:e} J, linear_drift={worst_linear:e}, angular_drift={worst_angular:e}, weak_volume={worst_volume:e}, length_squared_change={strain:e}"
        );
    });
}

#[test]
fn g0_exact_constraint_trajectory_is_bitwise_independent_of_unused_bulk_modulus() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let (ordinary, enormous) = (law(60.0), law(1e300));
        let densities = [DENSITY];
        let (u0, v0) = initial_state(&nodes);
        let mut first = definition(&nodes, &elements, &ordinary, &densities)
            .prepare(&u0, &v0, cx)
            .unwrap();
        let mut second = definition(&nodes, &elements, &enormous, &densities)
            .prepare(&u0, &v0, cx)
            .unwrap();
        assert_eq!(state_bits(&first), state_bits(&second));
        for _ in 0..4 {
            let a = first
                .advance(0.01, &[], &[], [0.0; 3], settings(), WORK_ALLOWANCE, cx)
                .unwrap();
            let b = second
                .advance(0.01, &[], &[], [0.0; 3], settings(), WORK_ALLOWANCE, cx)
                .unwrap();
            assert_eq!(state_bits(&a.solid), state_bits(&b.solid));
            assert_eq!(a.energy, b.energy);
            assert_eq!(a.inertia, b.inertia);
            assert_eq!(a.pressure_residual_m3, b.pressure_residual_m3);
            first = a.solid;
            second = b.solid;
        }
    });
}

#[test]
fn g0_incompatible_initial_volume_and_redundant_pressure_rows_refuse_without_mutation() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(1e300);
        let densities = [DENSITY];
        let definition = definition(&nodes, &elements, &material, &densities);
        let zero = vec![[0.0; 3]; nodes.len()];
        let dilated: Vec<_> = nodes.iter().map(|x| [0.02 * x[0], 0.0, 0.0]).collect();
        assert!(definition.prepare(&dilated, &zero, cx).is_err());
        // At F=I this velocity has div(v)=0.02/s, so each independent
        // normalized P1 weak volume rate is 0.02/4 = 0.005/s.
        assert!(definition.prepare(&zero, &dilated, cx).is_err());
        for rate in [-1.0, f64::NAN, f64::INFINITY] {
            let invalid = ConservativeTetDefinition {
                formulation: ConservativeFormulation::ExactIncompressible {
                    maximum_relative_volume_residual: 1e-9,
                    maximum_initial_relative_volume_rate_per_s: rate,
                    maximum_point_volume_work_defect: 1e-12,
                },
                ..definition
            };
            assert!(invalid.prepare(&zero, &zero, cx).is_err());
        }
        let accepted = definition.prepare(&zero, &zero, cx).unwrap();
        let before = state_bits(&accepted);
        let (_, volume) = independent_fields(&elements[0], &zero);
        assert_eq!(volume, [0.0; 4]);
        // Even an exactly solved residual cannot select a pressure gauge here:
        // the endpoint volume equations have no free displacement directions.
        let closed: Vec<_> = (0..3 * nodes.len()).map(|d| (d, 0.0)).collect();
        let refused =
            accepted.advance(0.01, &closed, &[], [0.0; 3], settings(), WORK_ALLOWANCE, cx);
        assert!(refused.is_err(), "{refused:?}");
        assert_eq!(state_bits(&accepted), before);
    });
}

fn impulse_controls() -> fs_solid::hyper3::mixed::impulse::MixedImpulseControls {
    use fs_solid::hyper3::mixed::impulse::{MixedImpulseControls, SolidImpulseControls};
    MixedImpulseControls {
        balance: SolidImpulseControls {
            solve: fs_solver::CheckedCgConfig {
                tolerance: 1e-11,
                max_iterations: 64,
                max_corrections: 0,
            },
            maximum_workspace_components: 8192,
            linear_tolerance_kg_m_s: 1e-10,
            angular_tolerance_kg_m2_s: 1e-10,
            work_tolerance_j: WORK_ALLOWANCE,
        },
        velocity_scale_m_s: 1.0,
        pressure_impulse_scale_pa_s: 1.0,
        work_scale_j: 1.0,
        maximum_assembly_bytes: 2 * 1024 * 1024,
        maximum_factor_bytes: 2 * 1024 * 1024,
        maximum_relative_volume_rate_per_s: 1e-10,
    }
}

#[test]
fn g1_g3_pressure_impulse_removes_only_volume_loading_and_keeps_rigid_momentum() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(1e300);
        let densities = [DENSITY];
        let gradient = [
            [-1.0, -1.0, -1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let translation = [0.2, -0.3, 0.1];
        for rotated in [false, true] {
            let rotate = |x: [f64; 3]| if rotated { [-x[1], x[0], x[2]] } else { x };
            let displacement: Vec<_> = nodes
                .iter()
                .map(|&x| {
                    let y = rotate(x);
                    [y[0] + 0.3 - x[0], y[1] - 0.2 - x[1], y[2] + 0.7 - x[2]]
                })
                .collect();
            let accepted = definition(&nodes, &elements, &material, &densities)
                .prepare(&displacement, &vec![[0.0; 3]; nodes.len()], cx)
                .unwrap();
            let before = state_bits(&accepted);
            let mut applied = Vec::new();
            for (a, &node) in elements[0].iter().enumerate() {
                // Exact P2 row mass is negative at corners. A uniform pressure
                // impulse integrates to zero there, and V*(grad La+grad Lb) at
                // an edge. These barycentric integrals do not use native B/M.
                let row_mass = DENSITY * VOLUME * if a < 4 { -0.05 } else { 0.2 };
                let volume_load = if a < 4 {
                    [0.0; 3]
                } else {
                    let [i, j] = TET10_EDGES[a - 4];
                    rotate(core::array::from_fn(|c| {
                        VOLUME * (gradient[i][c] + gradient[j][c])
                    }))
                };
                for c in 0..3 {
                    applied.push((3 * node + c, row_mass * translation[c] + volume_load[c]));
                }
            }
            let event = accepted
                .impulse_candidate(&[], &applied, None, [0.0; 3], impulse_controls(), cx)
                .unwrap();
            assert_eq!(state_bits(&accepted), before);
            assert_eq!(event.solid.displacement_m(), accepted.displacement_m());
            assert_eq!(
                event.solid.algorithmic_pressure_pa(),
                accepted.algorithmic_pressure_pa()
            );
            assert_eq!(
                event.solid.stored_energy_j().to_bits(),
                accepted.stored_energy_j().to_bits()
            );
            assert_eq!(event.solid.time_s(), accepted.time_s());
            assert!(event.impulse.momentum.solve.converged_euclidean());
            for velocity in event.solid.velocity_mps() {
                for c in 0..3 {
                    near(velocity[c], translation[c], 2e-10);
                }
            }
            for pressure in &event.impulse.pressure_impulses_pa_s {
                near(*pressure, 1.0, 2e-9);
            }
            let (kinetic, linear, angular) = independent_motion(
                &nodes,
                &elements[0],
                event.solid.displacement_m(),
                event.solid.velocity_mps(),
            );
            let center = rotate([0.25; 3]);
            let expected_angular = cross(
                [center[0] + 0.3, center[1] - 0.2, center[2] + 0.7],
                translation,
            );
            near(
                kinetic,
                0.5 * translation.iter().map(|v| v * v).sum::<f64>(),
                1e-11,
            );
            for c in 0..3 {
                near(linear[c], translation[c], 1e-11);
                near(angular[c], expected_angular[c], 1e-11);
                near(event.impulse.pressure.linear_kg_m_s[c], 0.0, 1e-11);
                near(event.impulse.pressure.angular_kg_m2_s[c], 0.0, 1e-11);
            }
            near(event.impulse.pressure_work_j, 0.0, 1e-11);
            near(event.impulse.external_work_j, kinetic, 1e-11);
            // Continue the original conservative owner: an impulse's supplied
            // energy must enter its cumulative work, without resetting time.
            let step = event
                .solid
                .advance(0.01, &[], &[], [0.0; 3], settings(), WORK_ALLOWANCE, cx)
                .unwrap();
            near(step.energy.cumulative_balance_error_j, 0.0, WORK_ALLOWANCE);
            near(
                step.solid.total_energy_j(),
                event.solid.total_energy_j(),
                WORK_ALLOWANCE,
            );
            assert!(
                accepted
                    .impulse_candidate(&[], &applied, Some(0), [0.0; 3], impulse_controls(), cx)
                    .is_err()
            );
        }
    });
}

#[test]
fn g1_supported_pressure_gauge_has_analytic_reactions_and_incompatible_events_refuse() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(60.0);
        let densities = [DENSITY];
        let definition = definition(&nodes, &elements, &material, &densities);
        let zero = vec![[0.0; 3]; nodes.len()];
        let accepted = definition.prepare(&zero, &zero, cx).unwrap();
        let before = state_bits(&accepted);
        // Only motion along the three edges leaving corner zero is free.
        // Constant pressure is a true null mode on these tangent directions.
        let free: Vec<_> = (0..3)
            .map(|c| {
                let edge = TET10_EDGES
                    .iter()
                    .position(|e| e.contains(&0) && e.contains(&(c + 1)))
                    .unwrap();
                3 * elements[0][4 + edge] + c
            })
            .collect();
        let held: Vec<_> = (0..3 * nodes.len()).filter(|d| !free.contains(d)).collect();
        let applied = [(free[0], 0.1)];
        let event = accepted
            .impulse_candidate(&held, &applied, Some(0), [0.0; 3], impulse_controls(), cx)
            .unwrap();
        for velocity in event.solid.velocity_mps().iter().flatten() {
            near(*velocity, 0.0, 1e-11);
        }
        // B[:,edge01.x]=V/5*[1,-1,0,0], lambda0=0, so lambda1=-3 Pa s.
        for (pressure, expected) in event
            .impulse
            .pressure_impulses_pa_s
            .iter()
            .zip([0.0, -3.0, 0.0, 0.0])
        {
            near(*pressure, expected, 1e-10);
        }
        near(event.impulse.momentum.support.linear_kg_m_s[0], -0.1, 1e-11);
        near(event.impulse.pressure_work_j, 0.0, 1e-11);
        assert!(
            accepted
                .impulse_candidate(&held, &applied, None, [0.0; 3], impulse_controls(), cx)
                .is_err()
        );
        assert!(
            accepted
                .impulse_candidate(&[0, 0], &applied, None, [0.0; 3], impulse_controls(), cx)
                .is_err()
        );
        let mut limited = impulse_controls();
        limited.maximum_factor_bytes = 1;
        assert!(
            accepted
                .impulse_candidate(&held, &applied, Some(0), [0.0; 3], limited, cx)
                .is_err()
        );
        assert_eq!(state_bits(&accepted), before);

        // Initial admission and event accuracy are explicit separate budgets.
        // The event cannot quietly project away a pre-existing expansion rate.
        let relaxed = ConservativeTetDefinition {
            formulation: ConservativeFormulation::ExactIncompressible {
                maximum_relative_volume_residual: 1e-9,
                maximum_initial_relative_volume_rate_per_s: 1.0,
                maximum_point_volume_work_defect: 1e-12,
            },
            ..definition
        };
        let expansion: Vec<_> = nodes.iter().map(|x| x.map(|v| 0.02 * v)).collect();
        let incompatible = relaxed.prepare(&zero, &expansion, cx).unwrap();
        let incompatible_before = state_bits(&incompatible);
        assert!(
            incompatible
                .impulse_candidate(&[], &[], None, [0.0; 3], impulse_controls(), cx)
                .is_err()
        );
        assert_eq!(state_bits(&incompatible), incompatible_before);
    });
}
