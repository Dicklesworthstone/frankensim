//! Independent finite-strain energy, momentum and time-accuracy checks.
//! The nonlinear shear oracle is a scalar quartic oscillator derived from the
//! original Ogden energy; its roots and continuous reference do not use the FE
//! residual or the new material discrete derivative.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_feec::highorder::simplex::duffy_quadrature;
use fs_material::hyper::discrete_gradient::DiscreteGradientSettings;
use fs_material::hyper::{Hyperelastic, HyperelasticModel, OgdenTerm, Volumetric};
use fs_solid::hyper3::mixed::conservative::{
    ConservativeFormulation, ConservativeTetDefinition, ConservativeTetError, ConservativeTetStep,
};
use fs_solid::hyper3::mixed::{MixedLawField, MixedUpError, MixedUpSettings};
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::linear3::tet10::{TET10_EDGES, promote_p1_to_p2, tet10_point};
use fs_solver::nonlinear::NewtonKrylovConfig;

const VOLUME: f64 = 1.0 / 6.0;
const DENSITY: f64 = 6.0;
const MU: f64 = 12.0;
const WORK_ALLOWANCE: f64 = 2e-9;

fn with_cx(cancelled: bool, run: impl FnOnce(&Cx<'_>)) {
    let cancel = CancelGate::new_clock_free();
    if cancelled {
        cancel.request();
    }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            &cancel,
            arena,
            StreamKey {
                seed: 62,
                kernel_id: 22,
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

fn law() -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::MooneyRivlin {
            c10: 4.0,
            c01: 2.0,
            kappa: 60.0,
        },
        4.0,
    )
    .unwrap()
}

fn material_settings() -> DiscreteGradientSettings {
    DiscreteGradientSettings {
        minimum_relative_strain_increment: 1e-6,
        absolute_work_tolerance_pa: 1e-10,
        relative_work_tolerance: 1e-12,
    }
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
        formulation: ConservativeFormulation::FullMaterial,
        budget: TetAssemblyBudget::standard(),
        quadrature_points_per_dir: 4,
        maximum_quadrature_points: 64,
        maximum_mass_assembly_bytes: 128 * 1024,
        material_settings: material_settings(),
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
        stress_scale_pa: MU,
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

fn check_work(step: &ConservativeTetStep<'_>) {
    let energy = &step.energy;
    near(
        energy.dead_load_work_j + energy.support_work_j,
        energy.stored_after_j - energy.stored_before_j + step.inertia.kinetic_change_j,
        WORK_ALLOWANCE,
    );
    near(
        energy.material_work_j,
        energy.stored_after_j - energy.stored_before_j,
        WORK_ALLOWANCE,
    );
    near(
        step.inertia.inertial_work_j,
        step.inertia.kinetic_change_j,
        2e-11,
    );
    assert_eq!(step.inertia.numerical_loss_j, 0.0);
    near(energy.balance_error_j, 0.0, WORK_ALLOWANCE);
    near(energy.cumulative_balance_error_j, 0.0, WORK_ALLOWANCE);
}

// Independent affine P2 mass entries / (rho*V), obtained by integrating
// barycentric monomials. Signed corner-edge entries must remain present.
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

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
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

fn independent_storage(
    nodes: &[[f64; 3]],
    element: &[usize; 10],
    displacement: &[[f64; 3]],
    material: &Hyperelastic,
) -> f64 {
    duffy_quadrature(4)
        .into_iter()
        .enumerate()
        .map(|(q, (lambda, weight))| {
            let point = tet10_point(element.map(|i| nodes[i]), lambda, 1e-10, 0, q).unwrap();
            let mut f = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
            for a in 0..10 {
                for i in 0..3 {
                    for j in 0..3 {
                        f[3 * i + j] += displacement[element[a]][i] * point.gradients[a][j];
                    }
                }
            }
            weight * point.determinant_m3 * material.energy(&f)
        })
        .sum()
}

#[test]
fn g1_free_deforming_tet_preserves_energy_and_both_momenta() {
    with_cx(false, |cx| {
        let (nodes, elements) = mesh(cx);
        let material = law();
        let initial_u: Vec<_> = nodes
            .iter()
            .map(|x| {
                [
                    0.08 * x[0] + 0.48 * x[1] * x[2],
                    -0.06 * x[1] + 0.03 * x[2],
                    0.04 * x[2],
                ]
            })
            .collect();
        let initial_v: Vec<_> = nodes
            .iter()
            .zip(&initial_u)
            .map(|(x, u)| {
                let position = core::array::from_fn(|i| x[i] + u[i] - 0.25);
                let spin = cross([0.3, -0.2, 0.5], position);
                [
                    0.2 + spin[0] + 0.1 * (x[0] - 0.25),
                    -0.15 + spin[1],
                    0.1 + spin[2],
                ]
            })
            .collect();
        let densities = [DENSITY];
        let mut solid = definition(&nodes, &elements, &material, &densities)
            .prepare(&initial_u, &initial_v, cx)
            .unwrap();
        let (kinetic0, linear0, angular0) =
            independent_motion(&nodes, &elements[0], &initial_u, &initial_v);
        let h0 = independent_storage(&nodes, &elements[0], &initial_u, &material) + kinetic0;
        let length2 = |u: &[[f64; 3]]| {
            (0..3)
                .map(|c| (nodes[1][c] + u[1][c] - nodes[0][c] - u[0][c]).powi(2))
                .sum::<f64>()
        };
        let initial_length2 = length2(&initial_u);
        let (mut worst_energy, mut worst_linear, mut worst_angular, mut deformation) =
            (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
        for _ in 0..48 {
            let step = solid
                .advance(0.01, &[], &[], [0.0; 3], settings(), WORK_ALLOWANCE, cx)
                .unwrap();
            check_work(&step);
            assert!(step.reactions_n.is_empty());
            let next = &step.solid;
            let (kinetic, linear, angular) = independent_motion(
                &nodes,
                &elements[0],
                next.displacement_m(),
                next.velocity_mps(),
            );
            let stored =
                independent_storage(&nodes, &elements[0], next.displacement_m(), &material);
            near(next.stored_energy_j(), stored, 2e-11);
            near(step.inertia.kinetic_after_j, kinetic, 2e-11);
            worst_energy = worst_energy.max((stored + kinetic - h0).abs());
            for c in 0..3 {
                worst_linear = worst_linear.max((linear[c] - linear0[c]).abs());
                worst_angular = worst_angular.max((angular[c] - angular0[c]).abs());
            }
            deformation = deformation.max((length2(next.displacement_m()) - initial_length2).abs());
            solid = step.solid;
        }
        assert!(
            deformation > 1e-3,
            "the free body must actually change strain"
        );
        near(worst_energy, 0.0, WORK_ALLOWANCE);
        near(worst_linear, 0.0, 1e-10);
        near(worst_angular, 0.0, 1e-10);
        println!(
            "free nonlinear tet: energy_drift={worst_energy:e} J, linear_drift={worst_linear:e}, angular_drift={worst_angular:e}, length_squared_change={deformation:e}"
        );
    });
}

// The exact nonlinear shear restriction is U(q)=a*q^2+b*q^4. Its continuous
// reference uses a small independent scalar RK4 integration; a doubled grid
// below bounds reference contamination of the measured second-order ladder.
fn scalar_reference(
    intervals: usize,
    duration: f64,
    amplitude: f64,
    mass: f64,
    a: f64,
    b: f64,
) -> (f64, f64) {
    let h = duration / intervals as f64;
    let force = |q: f64| 2.0 * a * q + 4.0 * b * q.powi(3);
    let (mut q, mut v) = (amplitude, 0.0);
    for _ in 0..intervals {
        let k1 = (v, -force(q) / mass);
        let k2 = (v + 0.5 * h * k1.1, -force(q + 0.5 * h * k1.0) / mass);
        let k3 = (v + 0.5 * h * k2.1, -force(q + 0.5 * h * k2.0) / mass);
        let k4 = (v + h * k3.1, -force(q + h * k3.0) / mass);
        q += h * (k1.0 + 2.0 * k2.0 + 2.0 * k3.0 + k4.0) / 6.0;
        v += h * (k1.1 + 2.0 * k2.1 + 2.0 * k3.1 + k4.1) / 6.0;
    }
    (q, v)
}

fn scalar_step(q: f64, v: f64, h: f64, amplitude: f64, mass: f64, a: f64, b: f64) -> (f64, f64) {
    let residual = |next: f64| {
        2.0 * mass / (h * h) * (next - q - h * v) + (next + q) * (a + b * (next * next + q * q))
    };
    let (mut lo, mut hi) = (-2.0 * amplitude, 2.0 * amplitude);
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..90 {
        let mid = 0.5 * (lo + hi);
        if residual(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let next = 0.5 * (lo + hi);
    (next, 2.0 * (next - q) / h - v)
}

#[test]
fn g1_quartic_ogden_shear_matches_independent_steps_and_second_order_time_ladder() {
    with_cx(false, |cx| {
        let (nodes, elements) = mesh(cx);
        let material = Hyperelastic::new(
            HyperelasticModel::Ogden {
                terms: vec![OgdenTerm {
                    mu_pa: MU / 2.0,
                    alpha: 4.0,
                }],
                volumetric: Volumetric::Quadratic { kappa: 10.0 * MU },
            },
            4.0,
        )
        .unwrap();
        let edge = TET10_EDGES
            .iter()
            .position(|e| e.contains(&2) && e.contains(&3))
            .unwrap();
        let node = elements[0][4 + edge];
        let dof = 3 * node;
        let (amplitude, duration) = (0.3, 0.2);
        let mut initial_u = vec![[0.0; 3]; nodes.len()];
        initial_u[node][0] = amplitude;
        let initial_v = vec![[0.0; 3]; nodes.len()];
        let fixed: Vec<_> = (0..3 * nodes.len())
            .filter(|&d| d != dof)
            .map(|d| (d, 0.0))
            .collect();
        // u_x=4*q*y*z gives J=1; sum(lambda^4)-3 = 4*s^2+s^4,
        // s^2=16*q^2*(y^2+z^2). Their simplex means are1/5 and1/15.
        let mass = 8.0 * DENSITY * VOLUME / 105.0;
        let a = 8.0 * MU * VOLUME / 5.0;
        let b = 32.0 * MU * VOLUME / 15.0;
        let potential = |q: f64| a * q * q + b * q.powi(4);
        let reference = scalar_reference(4096, duration, amplitude, mass, a, b);
        let finer_reference = scalar_reference(8192, duration, amplitude, mass, a, b);
        near(reference.0, finer_reference.0, 1e-11);
        near(reference.1, finer_reference.1, 1e-11);
        let mut errors = Vec::new();
        for intervals in [8, 16, 32] {
            let densities = [DENSITY];
            let mut solid = definition(&nodes, &elements, &material, &densities)
                .prepare(&initial_u, &initial_v, cx)
                .unwrap();
            near(solid.stored_energy_j(), potential(amplitude), 2e-12);
            let dt = duration / f64::from(intervals);
            let (mut q, mut v, mut worst_oracle) = (amplitude, 0.0, 0.0_f64);
            for _ in 0..intervals {
                (q, v) = scalar_step(q, v, dt, amplitude, mass, a, b);
                let step = solid
                    .advance(dt, &fixed, &[], [0.0; 3], settings(), WORK_ALLOWANCE, cx)
                    .unwrap();
                check_work(&step);
                let actual_q = step.solid.displacement_m()[node][0];
                let actual_v = step.solid.velocity_mps()[node][0];
                worst_oracle = worst_oracle.max((actual_q - q).abs());
                near(actual_q, q, 5e-10);
                near(actual_v, v, 2e-8);
                near(step.energy.stored_after_j, potential(actual_q), 2e-11);
                near(
                    step.inertia.kinetic_after_j,
                    0.5 * mass * actual_v * actual_v,
                    2e-11,
                );
                near(
                    step.energy.stored_after_j + step.inertia.kinetic_after_j,
                    potential(amplitude),
                    WORK_ALLOWANCE,
                );
                solid = step.solid;
            }
            let actual_q = solid.displacement_m()[node][0];
            let actual_v = solid.velocity_mps()[node][0];
            assert!(
                actual_q < 0.0,
                "the nonlinear mode must cross its reference"
            );
            let error = ((actual_q - reference.0) / amplitude)
                .hypot((actual_v - reference.1) / (amplitude * (2.0 * a / mass).sqrt()));
            errors.push(error);
            println!(
                "quartic P2 shear intervals={intervals}, trajectory_error={error:e}, scalar_oracle_error={worst_oracle:e}"
            );
        }
        for pair in errors.windows(2) {
            let order = (pair[0] / pair[1]).log2();
            assert!(
                (1.8..2.2).contains(&order),
                "nonlinear time order={order:e}, errors={errors:?}"
            );
        }
    });
}

#[test]
fn g1_constant_body_acceleration_has_exact_midpoint_motion_and_work() {
    with_cx(false, |cx| {
        let (nodes, elements) = mesh(cx);
        let material = law();
        let densities = [DENSITY];
        let velocity = [0.2, -0.1, 0.05];
        let acceleration = [0.3, -0.2, 0.1];
        let dt = 0.025;
        let mut solid = definition(&nodes, &elements, &material, &densities)
            .prepare(
                &vec![[0.0; 3]; nodes.len()],
                &vec![velocity; nodes.len()],
                cx,
            )
            .unwrap();
        for count in 1..=5 {
            let time = f64::from(count) * dt;
            let expected_v: [f64; 3] =
                core::array::from_fn(|i| velocity[i] + time * acceleration[i]);
            let expected_u: [f64; 3] =
                core::array::from_fn(|i| time * velocity[i] + 0.5 * time * time * acceleration[i]);
            let step = solid
                .advance(
                    dt,
                    &[],
                    &[],
                    acceleration.map(|a| DENSITY * a),
                    settings(),
                    WORK_ALLOWANCE,
                    cx,
                )
                .unwrap();
            check_work(&step);
            assert!(step.reactions_n.is_empty());
            for (u, v) in step
                .solid
                .displacement_m()
                .iter()
                .zip(step.solid.velocity_mps())
            {
                for c in 0..3 {
                    near(u[c], expected_u[c], 3e-11);
                    near(v[c], expected_v[c], 3e-9);
                }
            }
            near(step.energy.stored_after_j, 0.0, 2e-12);
            near(step.solid.reference_mass_kg(), DENSITY * VOLUME, 2e-14);
            near(
                step.inertia.kinetic_after_j,
                0.5 * DENSITY * VOLUME * expected_v.iter().map(|v| v * v).sum::<f64>(),
                2e-10,
            );
            solid = step.solid;
        }
    });
}

#[test]
fn g0_moving_support_work_and_rejected_candidates_preserve_accepted_state() {
    with_cx(false, |cx| {
        let (nodes, elements) = mesh(cx);
        let material = law();
        let densities = [DENSITY];
        let zero = vec![[0.0; 3]; nodes.len()];
        let solid = definition(&nodes, &elements, &material, &densities)
            .prepare(&zero, &zero, cx)
            .unwrap();
        let (dt, velocity) = (0.02, [0.3, -0.2, 0.1]);
        let supports: Vec<_> = (0..3 * nodes.len())
            .map(|d| (d, 0.5 * dt * velocity[d % 3]))
            .collect();
        let step = solid
            .advance(dt, &supports, &[], [0.0; 3], settings(), WORK_ALLOWANCE, cx)
            .unwrap();
        assert!(step.newton.is_none());
        for &(d, reaction) in &step.reactions_n {
            let fraction = if elements[0][..4].contains(&(d / 3)) {
                -1.0 / 20.0
            } else {
                1.0 / 5.0
            };
            near(
                reaction,
                DENSITY * VOLUME * fraction * velocity[d % 3] / dt,
                2e-10,
            );
        }
        let kinetic = 0.5 * DENSITY * VOLUME * velocity.iter().map(|v| v * v).sum::<f64>();
        near(step.energy.support_work_j, kinetic, 2e-11);
        near(step.inertia.kinetic_after_j, kinetic, 2e-11);
        check_work(&step);

        // Explicitly allow a coarse local strain-midpoint approximation, then
        // demand a tighter assembled work gate. This forces a real energy
        // rejection, independent of whether roundoff happens to be exactly0.
        let nonlinear = Hyperelastic::new(
            HyperelasticModel::Ogden {
                terms: vec![OgdenTerm {
                    mu_pa: MU / 2.0,
                    alpha: 4.0,
                }],
                volumetric: Volumetric::Quadratic { kappa: 10.0 * MU },
            },
            4.0,
        )
        .unwrap();
        let edge = TET10_EDGES
            .iter()
            .position(|e| e.contains(&2) && e.contains(&3))
            .unwrap();
        let dof = 3 * elements[0][4 + edge];
        let mut initial = zero.clone();
        initial[dof / 3][0] = 0.3;
        let mut definition = definition(&nodes, &elements, &nonlinear, &densities);
        definition.material_settings = DiscreteGradientSettings {
            minimum_relative_strain_increment: 1.0,
            absolute_work_tolerance_pa: 100.0,
            relative_work_tolerance: 0.0,
        };
        let accepted = definition.prepare(&initial, &zero, cx).unwrap();
        let prescribed: Vec<_> = (0..3 * nodes.len())
            .map(|d| (d, if d == dof { 0.18 } else { 0.0 }))
            .collect();
        let before = (
            accepted.displacement_m().to_vec(),
            accepted.velocity_mps().to_vec(),
            accepted.stored_energy_j(),
            accepted.time_s(),
        );
        let refused = accepted.advance(0.02, &prescribed, &[], [0.0; 3], settings(), 1e-12, cx);
        assert!(
            matches!(refused, Err(ConservativeTetError::EnergyGate { .. })),
            "{refused:?}"
        );
        with_cx(true, |cancelled| {
            assert!(matches!(
                accepted.advance(
                    0.02,
                    &prescribed,
                    &[],
                    [0.0; 3],
                    settings(),
                    WORK_ALLOWANCE,
                    cancelled
                ),
                Err(ConservativeTetError::Native(MixedUpError::Cancelled))
            ));
        });
        assert_eq!(
            before,
            (
                accepted.displacement_m().to_vec(),
                accepted.velocity_mps().to_vec(),
                accepted.stored_energy_j(),
                accepted.time_s()
            )
        );
    });
}
