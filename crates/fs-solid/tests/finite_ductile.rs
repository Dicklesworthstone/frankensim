//! G0/G3: original finite tetrahedra, multiplicative plasticity and damage.
//! Synthetic constants exercise the numerical law; no alloy is qualified here.
#![deny(unsafe_code)]

#[path = "finite_ductile/dynamics.rs"]
mod dynamics;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::ductile::evolution::DuctileJ2Settings;
use fs_material::ductile::{DuctileDomain, TabulatedDuctileLocus};
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::finite::ductile::FiniteDuctileJ2;
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters};
use fs_solid::brittle::TetCrackBandGeometry;
use fs_solid::ductile::{DuctileStepConditions, DuctileTemperatures};
use fs_solid::hyper3::plastic::FinitePlasticStepControls;
use fs_solid::hyper3::plastic::ductile::{FiniteDuctileTetError, FiniteDuctileTetSolid};
use fs_solid::loads3::TetLoadGeometry3;
use fs_solid::{IncrementSettings, TetAssemblyBudget};

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const ONSET: f64 = 0.002;
const GP: f64 = 0.05;

fn with_cx<R>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 37,
                kernel_id: 103,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&gate, &cx)
    })
}

fn model(width: f64, damage_limit: f64) -> FiniteDuctileJ2<TabulatedDuctileLocus> {
    let native = FiniteJ2::new(
        FiniteJ2Parameters {
            shear_modulus_pa: 500.0,
            bulk_modulus_pa: 1000.0,
            yield_stress_pa: 3.0,
            hardening_modulus_pa: 0.0,
        },
        FiniteJ2Limits {
            minimum_jacobian: 0.25,
            maximum_elastic_log_strain: 0.5,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
    )
    .unwrap();
    let locus = TabulatedDuctileLocus::new(
        vec![(-20.0, ONSET), (20.0, ONSET)],
        DuctileDomain {
            triaxiality: [-20.0, 20.0],
            strain_rate_s_inv: [1e-10, 10.0],
            temperature_k: [300.0, 300.0],
        },
    )
    .unwrap();
    FiniteDuctileJ2::new(
        native,
        locus,
        DuctileJ2Settings {
            characteristic_length_m: width,
            post_initiation_plastic_work_j_m2: GP,
            maximum_damage_increment: damage_limit,
            energy_tolerances: EnergyTolerances {
                energy_j_m3: 1e-9,
                yield_pa: 1e-8,
            },
        },
    )
    .unwrap()
}

fn conditions() -> DuctileStepConditions<'static> {
    DuctileStepConditions {
        duration_s: 0.1,
        temperatures_k: DuctileTemperatures::Uniform(300.0),
    }
}
fn controls() -> FinitePlasticStepControls {
    FinitePlasticStepControls {
        equilibrium_work_tolerance_j: 1e-9,
        maximum_endpoint_work_defect_j: 1.0,
        yield_tolerance_pa: 1e-8,
    }
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-9 * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}
fn stretch(e: f64) -> [f64; 9] {
    [
        e.exp(),
        0.0,
        0.0,
        0.0,
        (-0.5 * e).exp(),
        0.0,
        0.0,
        0.0,
        (-0.5 * e).exp(),
    ]
}
fn multiply(a: [f64; 9], b: [f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| {
        (0..3)
            .map(|k| a[3 * (ij / 3) + k] * b[3 * k + ij % 3])
            .sum()
    })
}
fn affine(nodes: &[[f64; 3]], f: [f64; 9], origin: [f64; 3], offset: usize) -> Vec<(usize, f64)> {
    nodes
        .iter()
        .enumerate()
        .flat_map(|(node, x)| {
            let x = core::array::from_fn::<_, 3, _>(|j| x[j] - origin[j]);
            (0..3).map(move |i| {
                (
                    3 * (node + offset) + i,
                    (0..3)
                        .map(|j| (f[3 * i + j] - if i == j { 1.0 } else { 0.0 }) * x[j])
                        .sum(),
                )
            })
        })
        .collect()
}

#[test]
fn affine_finite_plastic_path_loses_native_reaction_and_preserves_work() {
    with_cx(|_, cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let bands = [TetCrackBandGeometry::from_element(&geometry, 0, [1.0, 0.0, 0.0]).unwrap()];
        let material = model(bands[0].band_width_m(), 0.04);
        let materials = [&material];
        let mut solid = FiniteDuctileTetSolid::try_new(
            &geometry,
            &materials,
            &bands,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let mut previous_damage = 0.0;
        for index in 1..=50 {
            let strain = 0.08 * index as f64 / 50.0;
            let supports = affine(&NODES, stretch(strain), [0.0; 3], 0);
            let step = solid
                .advance(
                    &supports,
                    &[],
                    [0.0; 3],
                    conditions(),
                    IncrementSettings::default(),
                    controls(),
                    cx,
                )
                .unwrap();
            let alpha = (strain - 3.0 / 1500.0).max(0.0);
            let g = (-bands[0].band_width_m() * 3.0 * (alpha - ONSET).max(0.0) / GP).exp();
            let tau = (1000.0 * strain).min(2.0);
            let expected_piola = g * tau / strain.exp();
            near(step.solid.element_states()[0].native_state().alpha, alpha);
            near(step.solid.element_states()[0].integrity(), g);
            near(step.solid.element_first_piola_pa()[0][0], expected_piola);
            near(step.equilibrium.reactions_n[3], expected_piola / 6.0);
            let damage = step.solid.element_states()[0].damage();
            assert!(damage >= previous_damage);
            assert!(damage - previous_damage <= 0.04 + 1e-13);
            near(step.energy.balance_error_j, 0.0);
            near(step.energy.equilibrium_work_residual_j, 0.0);
            assert!(
                step.energy.absolute_endpoint_work_defect_j + 1e-14
                    >= step.energy.endpoint_work_defect_j.abs()
            );
            previous_damage = damage;
            solid = step.solid;
        }
        assert!(previous_damage > 0.75);
        assert!(solid.element_states()[0].dissipated_plastic_j_m3() > 0.0);
        assert!(solid.element_states()[0].dissipated_damage_j_m3() > 0.0);
    });
}

#[test]
fn free_degree_solves_the_damaged_finite_tangent_and_restarts_without_time_replay() {
    with_cx(|_, cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let bands = [TetCrackBandGeometry::from_element(&geometry, 0, [1.0, 0.0, 0.0]).unwrap()];
        let material = model(bands[0].band_width_m(), 1.0);
        let materials = [&material];
        let solid = FiniteDuctileTetSolid::try_new(
            &geometry,
            &materials,
            &bands,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        // F=diag(lambda,1,1): q_trial=2 mu log(lambda), tau_xx=K log(lambda)+2Y/3.
        let lambda = 1.02_f64;
        let alpha = (1000.0 * lambda.ln() - 3.0) / 1500.0;
        let g = (-bands[0].band_width_m() * 3.0 * (alpha - ONSET) / GP).exp();
        let piola = g * (1000.0 * lambda.ln() + 2.0) / lambda;
        let supports: Vec<_> = (0..12)
            .filter(|&dof| dof != 3)
            .map(|dof| (dof, 0.0))
            .collect();
        let force = [(3, piola / 6.0)];
        let step = solid
            .advance(
                &supports,
                &force,
                [0.0; 3],
                conditions(),
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        near(step.solid.displacement_m()[1][0], lambda - 1.0);
        near(step.solid.element_first_piola_pa()[0][0], piola);
        near(step.solid.element_states()[0].native_state().alpha, alpha);
        near(step.solid.element_states()[0].integrity(), g);
        assert!(step.solid.element_states()[0].damage() > 0.1);
        assert!(step.equilibrium.residual_history.len() > 1);
        near(step.energy.balance_error_j, 0.0);
        let repeated = step
            .solid
            .advance(
                &supports,
                &force,
                [0.0; 3],
                conditions(),
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(repeated.solid.element_states(), step.solid.element_states());
        near(repeated.energy.internal_endpoint_work_j, 0.0);
        near(repeated.energy.plastic_dissipation_j, 0.0);
        near(repeated.energy.damage_release_j, 0.0);
        assert_eq!(solid.element_states()[0], material.initial_state());
    });
}

#[test]
fn finite_rotations_keep_signed_defects_separate_and_opposing_cells_cannot_cancel_budget() {
    with_cx(|_, cx| {
        let mut nodes = NODES.to_vec();
        nodes.extend(NODES.map(|x| [x[0] + 2.0, x[1], x[2]]));
        let tets = [[0, 1, 2, 3], [4, 5, 6, 7]];
        let geometry =
            TetLoadGeometry3::new(&nodes, &tets, TetAssemblyBudget::standard(), 8, cx).unwrap();
        let bands = [
            TetCrackBandGeometry::from_element(&geometry, 0, [1.0, 0.0, 0.0]).unwrap(),
            TetCrackBandGeometry::from_element(&geometry, 1, [1.0, 0.0, 0.0]).unwrap(),
        ];
        let material = model(bands[0].band_width_m(), 0.04);
        let materials = [&material, &material];
        let mut solid = FiniteDuctileTetSolid::try_new(
            &geometry,
            &materials,
            &bands,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        for index in 1..=40 {
            let f = stretch(0.05 * index as f64 / 40.0);
            let mut supports = affine(&nodes[..4], f, [0.0; 3], 0);
            supports.extend(affine(&nodes[4..], f, [2.0, 0.0, 0.0], 4));
            solid = solid
                .advance(
                    &supports,
                    &[],
                    [0.0; 3],
                    conditions(),
                    IncrementSettings::default(),
                    controls(),
                    cx,
                )
                .unwrap()
                .solid;
        }
        let old = solid.element_states().to_vec();
        let f = stretch(0.05);
        let (s, c) = (core::f64::consts::PI / 3.0).sin_cos();
        let rx = [1.0, 0.0, 0.0, 0.0, c, -s, 0.0, s, c];
        let ry = [0.0, 0.0, 1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0];
        let mut supports = affine(&nodes[..4], multiply(rx, f), [0.0; 3], 0);
        supports.extend(affine(&nodes[4..], multiply(ry, f), [2.0, 0.0, 0.0], 4));
        let rotated = solid
            .advance(
                &supports,
                &[],
                [0.0; 3],
                conditions(),
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        for (before, after) in old.iter().zip(rotated.solid.element_states()) {
            near(after.native_state().alpha, before.native_state().alpha);
            near(after.integrity(), before.integrity());
            near(after.stored_energy_j_m3(), before.stored_energy_j_m3());
        }
        near(rotated.energy.plastic_dissipation_j, 0.0);
        near(rotated.energy.damage_release_j, 0.0);
        near(rotated.energy.endpoint_work_defect_j, 0.0);
        assert!(rotated.energy.absolute_endpoint_work_defect_j > 0.05);
        assert!(
            rotated.solid.element_states()[0].endpoint_work_defect_j_m3()
                < old[0].endpoint_work_defect_j_m3()
        );
        let tight = FinitePlasticStepControls {
            maximum_endpoint_work_defect_j: 0.5 * rotated.energy.absolute_endpoint_work_defect_j,
            ..controls()
        };
        assert!(matches!(
            solid.advance(
                &supports,
                &[],
                [0.0; 3],
                conditions(),
                IncrementSettings::default(),
                tight,
                cx
            ),
            Err(FiniteDuctileTetError::Energy {
                quantity: "sum of absolute element integration defects",
                ..
            })
        ));
        assert_eq!(solid.element_states(), old);
    });
}

#[test]
fn geometry_damage_and_cancellation_refusals_leave_the_original_finite_specimen_unchanged() {
    with_cx(|gate, cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let bands = [TetCrackBandGeometry::from_element(&geometry, 0, [1.0, 0.0, 0.0]).unwrap()];
        let material = model(bands[0].band_width_m(), 0.04);
        let wrong = model(2.0 * bands[0].band_width_m(), 0.04);
        assert!(
            FiniteDuctileTetSolid::try_new(
                &geometry,
                &[&wrong],
                &bands,
                TetAssemblyBudget::standard(),
                cx
            )
            .is_err()
        );
        let materials = [&material];
        let solid = FiniteDuctileTetSolid::try_new(
            &geometry,
            &materials,
            &bands,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let supports = affine(&NODES, stretch(0.08), [0.0; 3], 0);
        assert!(
            solid
                .advance(
                    &supports,
                    &[],
                    [0.0; 3],
                    conditions(),
                    IncrementSettings::default(),
                    controls(),
                    cx
                )
                .is_err()
        );
        assert_eq!(solid.element_states()[0], material.initial_state());
        assert_eq!(solid.displacement_m(), &[[0.0; 3]; 4]);
        let zero_time = DuctileStepConditions {
            duration_s: 0.0,
            ..conditions()
        };
        assert!(
            solid
                .advance(
                    &supports,
                    &[],
                    [0.0; 3],
                    zero_time,
                    IncrementSettings::default(),
                    controls(),
                    cx
                )
                .is_err()
        );
        gate.request();
        assert!(matches!(
            solid.advance(
                &supports,
                &[],
                [0.0; 3],
                conditions(),
                IncrementSettings::default(),
                controls(),
                cx
            ),
            Err(FiniteDuctileTetError::Solve(
                fs_solid::plastic3::SmallStrainTetError::Cancelled
            ))
        ));
        assert_eq!(solid.element_states()[0], material.initial_state());
    });
}
