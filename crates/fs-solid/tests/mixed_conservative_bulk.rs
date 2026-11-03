//! Original finite-bulk storage in conservative P2/P1 motion.
//! The physical bulk energy, weak volume moments and rational P2 inertia below
//! are independent of the assembled residual and local secant inverse.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_feec::highorder::simplex::duffy_quadrature;
use fs_material::hyper::discrete_gradient::DiscreteGradientSettings;
use fs_material::hyper::discrete_gradient::bulk::DiscreteBulkSettings;
use fs_material::hyper::{Hyperelastic, HyperelasticModel, OgdenTerm, Volumetric};
use fs_solid::hyper3::mixed::conservative::{
    ConservativeFormulation, ConservativeTetDefinition, ConservativeTetSolid,
};
use fs_solid::hyper3::mixed::{MixedLawField, MixedUpSettings};
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::linear3::tet10::{TET10_EDGES, promote_p1_to_p2};
use fs_solver::nonlinear::NewtonKrylovConfig;

const VOLUME: f64 = 1.0 / 6.0;
const DENSITY: f64 = 6.0;
const KAPPA: f64 = 60.0;
const WORK_ALLOWANCE: f64 = 2e-9;

#[path = "mixed_conservative_bulk/impulse.rs"]
mod impulse;

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
                kernel_id: 24,
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

fn law(logarithmic: bool, kappa: f64) -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::Ogden {
            // alpha=2 gives the independently evaluated remainder
            // W_iso=mu/2*(J^(-2/3)*|F|^2-3).
            terms: vec![OgdenTerm {
                mu_pa: 12.0,
                alpha: 2.0,
            }],
            volumetric: if logarithmic {
                Volumetric::LogSquared { kappa }
            } else {
                Volumetric::Quadratic { kappa }
            },
        },
        4.0,
    )
    .unwrap()
}

fn bulk_settings() -> DiscreteBulkSettings {
    DiscreteBulkSettings {
        minimum_volume_ratio: 0.25,
        maximum_volume_ratio: 2.5,
        maximum_iterations: 96,
        absolute_pressure_tolerance_pa: 1e-11,
        relative_pressure_tolerance: 1e-12,
        absolute_work_tolerance_pa: 1e-10,
        relative_work_tolerance: 1e-12,
    }
}

fn definition<'a>(
    nodes_m: &'a [[f64; 3]],
    elements: &'a [[usize; 10]],
    material: MixedLawField<'a>,
    densities: &'a [f64],
) -> ConservativeTetDefinition<'a> {
    ConservativeTetDefinition {
        nodes_m,
        elements,
        material,
        reference_density_kg_m3: densities,
        formulation: ConservativeFormulation::OriginalMaterial {
            bulk: bulk_settings(),
            maximum_relative_volume_residual: 1e-9,
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
    let displacement: Vec<_> = nodes
        .iter()
        .map(|x| [0.04 * x[0] + 0.25 * x[1] + 0.2 * x[1] * x[2], 0.0, 0.0])
        .collect();
    let velocity = nodes
        .iter()
        .zip(&displacement)
        .map(|(x, u)| {
            let p = core::array::from_fn(|c| x[c] + u[c] - 0.25);
            let spin = cross([0.3, -0.2, 0.5], p);
            // Nonzero divergence is physical for the finite-bulk formulation.
            [
                0.2 + spin[0] + 0.08 * p[0],
                -0.15 + spin[1] - 0.04 * p[1],
                0.1 + spin[2] + 0.02 * p[2],
            ]
        })
        .collect();
    (displacement, velocity)
}

fn bulk_energy(logarithmic: bool, kappa: f64, j: f64) -> f64 {
    let strain = if logarithmic { j.ln() } else { j - 1.0 };
    0.5 * kappa * strain * strain
}

fn initial_pressure(logarithmic: bool) -> [f64; 4] {
    // The nonaffine triangular initial map has J=1.04 pointwise. Its original
    // endpoint bulk stress is a continuous P1 field, so the initial mixed
    // state satisfies both constitutive and weak-volume equations.
    let j = 1.04_f64;
    [if logarithmic {
        KAPPA * j.ln() / j
    } else {
        KAPPA * (j - 1.0)
    }; 4]
}

// Exact affine P2 mass / (rho*V), from barycentric monomial integrals.
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
    u: &[[f64; 3]],
    v: &[[f64; 3]],
) -> (f64, [f64; 3], [f64; 3]) {
    let mut momentum = vec![[0.0; 3]; nodes.len()];
    for (a, &i) in element.iter().enumerate() {
        for (b, &j) in element.iter().enumerate() {
            for c in 0..3 {
                momentum[i][c] += DENSITY * VOLUME * mass_fraction(a, b) * v[j][c];
            }
        }
    }
    let (mut kinetic, mut linear, mut angular) = (0.0, [0.0; 3], [0.0; 3]);
    for i in 0..nodes.len() {
        let spin = cross(core::array::from_fn(|c| nodes[i][c] + u[i][c]), momentum[i]);
        for c in 0..3 {
            kinetic += 0.5 * v[i][c] * momentum[i][c];
            linear[c] += momentum[i][c];
            angular[c] += spin[c];
        }
    }
    (kinetic, linear, angular)
}

// Reference simplex geometry is exact. Retained constitutive volumes are
// deliberately independent inputs: mixed weak constraints do not require
// equality to geometric J at every point.
fn independent_fields(
    element: &[usize; 10],
    u: &[[f64; 3]],
    jbar: &[f64],
    logarithmic: bool,
    kappas: &[f64],
) -> (f64, f64, [f64; 4], f64) {
    assert_eq!(jbar.len(), 64);
    let (mut remainder, mut bulk, mut volume, mut mismatch) = (0.0, 0.0, [0.0; 4], 0.0_f64);
    for (q, (lambda, weight)) in duffy_quadrature(4).into_iter().enumerate() {
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
        let j = f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
            + f[2] * (f[3] * f[7] - f[4] * f[6]);
        assert!(j > 0.0);
        remainder +=
            weight * 6.0 * (j.powf(-2.0 / 3.0) * f.iter().map(|v| v * v).sum::<f64>() - 3.0);
        bulk += weight
            * bulk_energy(
                logarithmic,
                if kappas.len() == 1 {
                    kappas[0]
                } else {
                    kappas[q]
                },
                jbar[q],
            );
        for a in 0..4 {
            volume[a] += weight * lambda[a] * (j - jbar[q]);
        }
        mismatch = mismatch.max((j - jbar[q]).abs());
    }
    (remainder, bulk, volume, mismatch)
}

fn state_bits(solid: &ConservativeTetSolid<'_>) -> Vec<u64> {
    solid
        .displacement_m()
        .iter()
        .flatten()
        .chain(solid.velocity_mps().iter().flatten())
        .chain(solid.algorithmic_pressure_pa())
        .chain(solid.bulk_volume_ratios())
        .copied()
        .chain([
            solid.stored_energy_j(),
            solid.bulk_stored_energy_j(),
            solid.kinetic_energy_j(),
            solid.time_s(),
        ])
        .map(f64::to_bits)
        .collect()
}

#[test]
fn g1_original_bulk_free_motion_conserves_physical_energy_and_momenta() {
    with_cx(false, |cx| {
        let (nodes, elements) = mesh(cx);
        let densities = [DENSITY];
        let (u0, v0) = initial_state(&nodes);
        for logarithmic in [false, true] {
            let material = law(logarithmic, KAPPA);
            let mut solid = definition(
                &nodes,
                &elements,
                MixedLawField::Uniform(&material),
                &densities,
            )
            .prepare_with_pressure(&u0, &v0, &initial_pressure(logarithmic), cx)
            .unwrap();
            let (r0, b0, weak0, mismatch0) = independent_fields(
                &elements[0],
                &u0,
                solid.bulk_volume_ratios(),
                logarithmic,
                &[KAPPA],
            );
            let (k0, p0, l0) = independent_motion(&nodes, &elements[0], &u0, &v0);
            near(mismatch0, 0.0, 2e-15);
            for value in weak0 {
                near(value, 0.0, 2e-16);
            }
            near(solid.bulk_stored_energy_j(), b0, 2e-12);
            near(solid.total_energy_j(), r0 + b0 + k0, 2e-12);
            assert!(b0 > 1e-3);
            let (mut h_error, mut p_error, mut l_error, mut volume_error, mut bulk_change) =
                (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
            for _ in 0..32 {
                let before_bulk = solid.bulk_stored_energy_j();
                let before_remainder = solid.stored_energy_j() - before_bulk;
                let step = solid
                    .advance(0.005, &[], &[], [0.0; 3], settings(), WORK_ALLOWANCE, cx)
                    .unwrap();
                let next = &step.solid;
                let (remainder, bulk, weak, _) = independent_fields(
                    &elements[0],
                    next.displacement_m(),
                    next.bulk_volume_ratios(),
                    logarithmic,
                    &[KAPPA],
                );
                let (kinetic, linear, angular) = independent_motion(
                    &nodes,
                    &elements[0],
                    next.displacement_m(),
                    next.velocity_mps(),
                );
                near(next.bulk_stored_energy_j(), bulk, 2e-11);
                near(next.stored_energy_j(), remainder + bulk, 2e-11);
                near(next.kinetic_energy_j(), kinetic, 2e-11);
                near(
                    step.energy.pressure_work_j,
                    bulk - before_bulk,
                    WORK_ALLOWANCE,
                );
                near(
                    step.energy.material_work_j,
                    remainder - before_remainder,
                    WORK_ALLOWANCE,
                );
                near(step.energy.balance_error_j, 0.0, WORK_ALLOWANCE);
                near(step.energy.cumulative_balance_error_j, 0.0, WORK_ALLOWANCE);
                assert_eq!(step.inertia.numerical_loss_j, 0.0);
                for a in 0..4 {
                    near(step.pressure_residual_m3[a], weak[a], 2e-15);
                    volume_error = volume_error.max(weak[a].abs() / VOLUME);
                }
                h_error = h_error.max((remainder + bulk + kinetic - r0 - b0 - k0).abs());
                bulk_change = bulk_change.max((bulk - b0).abs());
                for c in 0..3 {
                    p_error = p_error.max((linear[c] - p0[c]).abs());
                    l_error = l_error.max((angular[c] - l0[c]).abs());
                }
                solid = step.solid;
            }
            near(h_error, 0.0, WORK_ALLOWANCE);
            near(p_error, 0.0, 1e-10);
            near(l_error, 0.0, 1e-10);
            near(volume_error, 0.0, 1e-9);
            assert!(
                bulk_change > 1e-3,
                "finite compressibility must exchange physical bulk energy"
            );
            println!(
                "finite bulk log={logarithmic}: H={h_error:e} J, P={p_error:e}, L={l_error:e}, weak={volume_error:e}, bulk_change={bulk_change:e} J"
            );
        }
    });
}

#[test]
fn g1_log_bulk_retained_history_has_second_order_time_convergence() {
    with_cx(false, |cx| {
        let (nodes, elements) = mesh(cx);
        let densities = [DENSITY];
        let material = law(true, KAPPA);
        let (initial_u, initial_v) = initial_state(&nodes);
        let definition = definition(
            &nodes,
            &elements,
            MixedLawField::Uniform(&material),
            &densities,
        );
        let duration = 0.08;
        let evolve = |intervals: usize| {
            let mut solid = definition
                .prepare_with_pressure(&initial_u, &initial_v, &initial_pressure(true), cx)
                .unwrap();
            let initial_energy = solid.total_energy_j();
            for _ in 0..intervals {
                let step = solid
                    .advance(
                        duration / intervals as f64,
                        &[],
                        &[],
                        [0.0; 3],
                        settings(),
                        WORK_ALLOWANCE,
                        cx,
                    )
                    .unwrap();
                near(step.solid.total_energy_j(), initial_energy, WORK_ALLOWANCE);
                solid = step.solid;
            }
            (
                solid.displacement_m().to_vec(),
                solid.velocity_mps().to_vec(),
                solid.bulk_volume_ratios().to_vec(),
            )
        };
        let reference = evolve(128);
        let mut errors = Vec::new();
        for intervals in [8, 16, 32] {
            let actual = evolve(intervals);
            // Dimensionless RMS trajectory metric with explicit SI scales:
            // 4 cm displacement, 0.5 m/s velocity, 0.05 volume-ratio change.
            let displacement_error: f64 = actual
                .0
                .iter()
                .flatten()
                .zip(reference.0.iter().flatten())
                .map(|(a, b)| ((a - b) / 0.04).powi(2))
                .sum();
            let velocity_error: f64 = actual
                .1
                .iter()
                .flatten()
                .zip(reference.1.iter().flatten())
                .map(|(a, b)| ((a - b) / 0.5).powi(2))
                .sum();
            let volume_error: f64 = actual
                .2
                .iter()
                .zip(&reference.2)
                .map(|(a, b)| ((a - b) / 0.05).powi(2))
                .sum();
            let error = ((displacement_error + velocity_error) / (3 * nodes.len()) as f64
                + volume_error / actual.2.len() as f64)
                .sqrt();
            assert!(error.is_finite() && error > 0.0);
            errors.push(error);
            println!("original log bulk intervals={intervals}, trajectory_error={error:e}");
        }
        for pair in errors.windows(2) {
            let order = (pair[0] / pair[1]).log2();
            assert!(
                (1.8..2.2).contains(&order),
                "retained original log bulk time order={order:e}, errors={errors:?}"
            );
        }
    });
}

#[test]
fn g1_homogeneous_dilation_uses_original_secant_pressure_and_support_work() {
    with_cx(false, |cx| {
        let (nodes, elements) = mesh(cx);
        let densities = [DENSITY];
        let zero = vec![[0.0; 3]; nodes.len()];
        for logarithmic in [false, true] {
            let material = law(logarithmic, KAPPA);
            let mut solid = definition(
                &nodes,
                &elements,
                MixedLawField::Uniform(&material),
                &densities,
            )
            .prepare(&zero, &zero, cx)
            .unwrap();
            let mut before_j = 1.0;
            for stretch in [1.08_f64, 0.94, 1.04] {
                let after_j = stretch.powi(3);
                let pressure = (bulk_energy(logarithmic, KAPPA, after_j)
                    - bulk_energy(logarithmic, KAPPA, before_j))
                    / (after_j - before_j);
                let prescribed: Vec<_> = (0..3 * nodes.len())
                    .map(|d| (d, (stretch - 1.0) * nodes[d / 3][d % 3]))
                    .collect();
                let before_energy = solid.total_energy_j();
                let step = solid
                    .advance(
                        0.05,
                        &prescribed,
                        &[],
                        [0.0; 3],
                        settings(),
                        WORK_ALLOWANCE,
                        cx,
                    )
                    .unwrap();
                // Finite pressure compliance resolves all four pressure DOFs
                // even with every displacement prescribed; no gauge is needed.
                assert!(step.newton.is_some());
                for &p in step.solid.algorithmic_pressure_pa() {
                    near(p, pressure, 2e-8);
                }
                for &j in step.solid.bulk_volume_ratios() {
                    near(j, after_j, 2e-10);
                }
                for &r in &step.pressure_residual_m3 {
                    near(r, 0.0, VOLUME * 1e-9);
                }
                let bulk = VOLUME * bulk_energy(logarithmic, KAPPA, after_j);
                let kinetic = independent_motion(
                    &nodes,
                    &elements[0],
                    step.solid.displacement_m(),
                    step.solid.velocity_mps(),
                )
                .0;
                near(step.solid.bulk_stored_energy_j(), bulk, 2e-10);
                near(step.solid.stored_energy_j(), bulk, 2e-10);
                near(
                    step.energy.pressure_work_j,
                    VOLUME * pressure * (after_j - before_j),
                    WORK_ALLOWANCE,
                );
                near(
                    step.energy.support_work_j,
                    bulk + kinetic - before_energy,
                    WORK_ALLOWANCE,
                );
                near(step.energy.balance_error_j, 0.0, WORK_ALLOWANCE);
                assert_eq!(step.energy.dead_load_work_j, 0.0);
                solid = step.solid;
                before_j = after_j;
            }
        }
    });
}

#[test]
fn g0_quadrature_materials_retain_distinct_bulk_compliance_and_storage() {
    with_cx(false, |cx| {
        let (nodes, elements) = mesh(cx);
        let densities = [DENSITY];
        let kappas: Vec<_> = (0..64)
            .map(|q| 20.0 + f64::from((q * 17 % 101) as u32))
            .collect();
        let materials: Vec<_> = kappas.iter().map(|&k| law(false, k)).collect();
        let zero = vec![[0.0; 3]; nodes.len()];
        let solid = definition(
            &nodes,
            &elements,
            MixedLawField::PerQuadrature(&materials),
            &densities,
        )
        .prepare(&zero, &zero, cx)
        .unwrap();
        let stretch = 1.04;
        let prescribed: Vec<_> = (0..3 * nodes.len())
            .map(|d| (d, (stretch - 1.0) * nodes[d / 3][d % 3]))
            .collect();
        let step = solid
            .advance(
                0.05,
                &prescribed,
                &[],
                [0.0; 3],
                settings(),
                WORK_ALLOWANCE,
                cx,
            )
            .unwrap();
        let mut pressure_work = 0.0;
        for (q, (lambda, weight)) in duffy_quadrature(4).into_iter().enumerate() {
            let p: f64 = (0..4)
                .map(|a| lambda[a] * step.solid.algorithmic_pressure_pa()[a])
                .sum();
            let j = step.solid.bulk_volume_ratios()[q];
            near(j, 1.0 + 2.0 * p / kappas[q], 2e-12);
            pressure_work += weight * p * (j - 1.0);
        }
        let (remainder, bulk, weak, mismatch) = independent_fields(
            &elements[0],
            step.solid.displacement_m(),
            step.solid.bulk_volume_ratios(),
            false,
            &kappas,
        );
        assert!(
            mismatch > 1e-3,
            "heterogeneous weak compliance must not be averaged into pointwise J"
        );
        near(remainder, 0.0, 2e-12);
        near(step.solid.bulk_stored_energy_j(), bulk, 2e-12);
        near(step.energy.pressure_work_j, pressure_work, WORK_ALLOWANCE);
        near(pressure_work, bulk, WORK_ALLOWANCE);
        for a in 0..4 {
            near(weak[a], 0.0, VOLUME * 1e-9);
            near(step.pressure_residual_m3[a], weak[a], 2e-15);
        }
    });
}

#[test]
fn g0_expanded_initial_volume_requires_an_admissible_explicit_pressure_predictor() {
    with_cx(false, |cx| {
        let (nodes, elements) = mesh(cx);
        let densities = [DENSITY];
        let material = law(false, KAPPA);
        let zero = vec![[0.0; 3]; nodes.len()];
        let initial_j = 2.1_f64;
        let stretch = initial_j.cbrt();
        let initial: Vec<_> = nodes
            .iter()
            .map(|x| x.map(|v| (stretch - 1.0) * v))
            .collect();
        let definition = definition(
            &nodes,
            &elements,
            MixedLawField::Uniform(&material),
            &densities,
        );
        // Zero initial endpoint pressure reconstructs Jbar=1, incompatible
        // with the supplied geometric J=2.1. Nothing is projected.
        assert!(definition.prepare(&initial, &zero, cx).is_err());
        let pressure = KAPPA * (initial_j - 1.0);
        let solid = definition
            .prepare_with_pressure(&initial, &zero, &[pressure; 4], cx)
            .unwrap();
        assert_eq!(solid.displacement_m(), initial);
        for &j in solid.bulk_volume_ratios() {
            near(j, initial_j, 2e-14);
        }
        let stored = VOLUME * bulk_energy(false, KAPPA, initial_j);
        near(solid.bulk_stored_energy_j(), stored, 2e-12);
        let prescribed: Vec<_> = initial.iter().flatten().copied().enumerate().collect();
        let step = solid
            .advance(
                0.05,
                &prescribed,
                &[],
                [0.0; 3],
                settings(),
                WORK_ALLOWANCE,
                cx,
            )
            .unwrap();
        for &p in step.solid.algorithmic_pressure_pa() {
            near(p, pressure, 2e-9);
        }
        near(step.solid.bulk_stored_energy_j(), stored, 2e-12);
        near(step.solid.total_energy_j(), stored, 2e-12);
        near(step.energy.pressure_work_j, 0.0, 2e-12);
    });
}

#[test]
fn g0_finite_bulk_domains_budgets_and_cancelled_steps_preserve_accepted_history() {
    with_cx(false, |cx| {
        let (nodes, elements) = mesh(cx);
        let densities = [DENSITY];
        let material = law(true, KAPPA);
        let zero = vec![[0.0; 3]; nodes.len()];
        let velocity: Vec<_> = nodes.iter().map(|x| [0.02 * x[0], 0.0, 0.0]).collect();
        let mut definition = definition(
            &nodes,
            &elements,
            MixedLawField::Uniform(&material),
            &densities,
        );
        definition.formulation = ConservativeFormulation::OriginalMaterial {
            bulk: DiscreteBulkSettings {
                maximum_volume_ratio: 1.1,
                ..bulk_settings()
            },
            maximum_relative_volume_residual: 1e-9,
            maximum_point_volume_work_defect: 1e-12,
        };
        let accepted = definition.prepare(&zero, &velocity, cx).unwrap();
        let before = state_bits(&accepted);
        let outside: Vec<_> = (0..3 * nodes.len())
            .map(|d| (d, 0.1 * nodes[d / 3][d % 3]))
            .collect();
        assert!(
            accepted
                .advance(
                    0.05,
                    &outside,
                    &[],
                    [0.0; 3],
                    settings(),
                    WORK_ALLOWANCE,
                    cx
                )
                .is_err()
        );
        let beyond: Vec<_> = nodes.iter().map(|x| x.map(|v| 0.1 * v)).collect();
        assert!(definition.prepare(&beyond, &zero, cx).is_err());
        let limited = MixedUpSettings {
            maximum_assembly_bytes: 1,
            ..settings()
        };
        assert!(
            accepted
                .advance(0.005, &[], &[], [0.0; 3], limited, WORK_ALLOWANCE, cx)
                .is_err()
        );
        with_cx(true, |cancelled| {
            assert!(
                accepted
                    .advance(
                        0.005,
                        &[],
                        &[],
                        [0.0; 3],
                        settings(),
                        WORK_ALLOWANCE,
                        cancelled
                    )
                    .is_err()
            );
        });
        assert_eq!(state_bits(&accepted), before);
        let coupled = Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 12.0,
                lambda: 60.0,
            },
            4.0,
        )
        .unwrap();
        definition.material = MixedLawField::Uniform(&coupled);
        assert!(definition.prepare(&zero, &zero, cx).is_err());
    });
}
