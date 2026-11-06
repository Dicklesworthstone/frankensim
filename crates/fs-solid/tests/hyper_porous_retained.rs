//! G0/G1/G3 retained pressure-sensitive solids on the original native solver.
//! Independent hydrostatic formulas and rigid motions verify implementation;
//! these synthetic coefficients do not establish measured foam accuracy.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::plastic::finite::porous::{
    FinitePorous, FinitePorousLimits, FinitePorousParameters,
};
use fs_solid::hyper3::plastic::porous::{FinitePorousTetError, FinitePorousTetSolid};
use fs_solid::hyper3::plastic::{
    FinitePlasticLawField, FinitePlasticStepControls, FinitePlasticTetError,
};
use fs_solid::hyper3::{HyperTetError, HyperTetSettings};
use fs_solid::linear3::TetAssemblyBudget;

const MU: f64 = 80.0;
const K: f64 = 160.0;
const Y: f64 = 2.0;
const H: f64 = 8.0;
const A: f64 = 1.2;
const NODES: [[f64; 3]; 4] = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const SUPPORTS: [(usize, f64); 6] = [(0, 0.0), (1, 0.0), (2, 0.0), (4, 0.0), (5, 0.0), (8, 0.0)];

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 42,
                kernel_id: 99,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate)
    })
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}, tolerance {tolerance:.3e}"
    );
}

fn law(pressure_sensitivity: f64) -> FinitePorous {
    FinitePorous::new(
        FinitePorousParameters {
            shear_modulus_pa: MU,
            bulk_modulus_pa: K,
            yield_stress_pa: Y,
            hardening_modulus_pa: H,
            pressure_sensitivity,
        },
        FinitePorousLimits {
            minimum_jacobian: 0.2,
            maximum_elastic_log_strain: 0.5,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
    )
    .unwrap()
}

fn controls() -> FinitePlasticStepControls {
    FinitePlasticStepControls {
        equilibrium_work_tolerance_j: 2e-7,
        maximum_endpoint_work_defect_j: 2.0,
        yield_tolerance_pa: 2e-7,
    }
}

// Independent hydrostatic reduction. beta maps signed mean Kirchhoff stress
// to the normalized ellipse radius and alpha to plastic log-volume loss.
fn hydrostatic(stretch: f64) -> (f64, f64, f64) {
    let beta = A / (1.0 + A * A / 9.0).sqrt();
    let alpha = (-3.0 * K * beta * stretch.ln() - Y) / (H + K * beta * beta);
    let plastic_j = (-beta * alpha).exp();
    let tau = -(Y + H * alpha) / beta;
    (alpha, plastic_j, tau)
}

fn affine(nodes: &[[f64; 3]], f: [f64; 9], shift: [f64; 3]) -> Vec<(usize, f64)> {
    nodes
        .iter()
        .enumerate()
        .flat_map(|(node, x)| {
            (0..3).map(move |i| {
                let current = (0..3).fold(shift[i], |sum, j| f[3 * i + j].mul_add(x[j], sum));
                (3 * node + i, current - x[i])
            })
        })
        .collect()
}

#[test]
fn force_driven_compaction_unloads_and_reloads_with_original_history_and_volume() {
    with_cx(|cx, _| {
        let material = law(A);
        let mut solid = FinitePorousTetSolid::try_new(
            &NODES,
            &TETS,
            FinitePlasticLawField::Uniform(&material),
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let peak = 0.8;
        let (alpha, plastic_j, peak_tau) = hydrostatic(peak);
        let unloaded = plastic_j.cbrt();
        let reloaded = (peak * unloaded).sqrt();
        let mut heat = 0.0;
        for (index, (stretch, tau)) in [
            (peak, peak_tau),
            (unloaded, 0.0),
            (reloaded, 0.5 * peak_tau),
        ]
        .into_iter()
        .enumerate()
        {
            let nominal = tau / stretch;
            // Six free deformation coordinates: three axial stretches and
            // three shears. V0=1/6; reference shape gradients are Cartesian.
            let forces = [(3, nominal / 6.0), (7, nominal / 6.0), (11, nominal / 6.0)];
            let old_state = solid.element_states().to_vec();
            let old_u = solid.displacement_m().to_vec();
            let old_storage = solid.stored_energy_j();
            let step = solid
                .advance(
                    &SUPPORTS,
                    &forces,
                    [0.0; 3],
                    HyperTetSettings::default(),
                    controls(),
                    cx,
                )
                .unwrap();
            assert!(step.equilibrium.steps[0].linear_iterations > 0);
            assert_eq!(solid.element_states(), old_state);
            assert_eq!(solid.displacement_m(), old_u);
            near(solid.stored_energy_j(), old_storage, 0.0);
            near(step.energy.stored_before_j, old_storage, 0.0);
            for i in 0..9 {
                near(
                    step.solid.element_deformation()[0][i],
                    if i % 4 == 0 { stretch } else { 0.0 },
                    3e-8,
                );
                near(
                    step.solid.element_first_piola_pa()[0][i],
                    if i % 4 == 0 { nominal } else { 0.0 },
                    4e-7,
                );
            }
            let state = &step.solid.element_states()[0];
            near(state.alpha, alpha, 3e-8);
            near(material.plastic_jacobian(state).unwrap(), plastic_j, 3e-8);
            let stored = (0.5 * tau * tau / K + 0.5 * H * alpha * alpha) / 6.0;
            let dissipated = if index == 0 { Y * alpha / 6.0 } else { 0.0 };
            near(step.solid.stored_energy_j(), stored, 4e-7);
            near(step.energy.stored_after_j, stored, 4e-7);
            near(step.energy.plastic_dissipation_j, dissipated, 3e-7);
            near(
                step.energy.incremental_potential_j,
                stored + dissipated,
                4e-7,
            );
            near(step.equilibrium.strain_energy_j, stored + dissipated, 4e-7);
            near(step.energy.equilibrium_work_residual_j, 0.0, 2e-7);
            near(
                step.energy.dead_load_endpoint_work_j + step.energy.support_endpoint_work_j,
                stored - old_storage + dissipated + step.energy.endpoint_work_defect_j,
                8e-7,
            );
            if index > 0 {
                assert_eq!(
                    step.solid.element_states(),
                    old_state,
                    "elastic unload/reload must retain compaction"
                );
                near(step.energy.plastic_dissipation_j, 0.0, 0.0);
            }
            heat += step.energy.plastic_dissipation_j;
            solid = step.solid;
        }
        assert!(plastic_j < 0.55);
        assert!(unloaded < 0.82);
        near(heat, Y * alpha / 6.0, 4e-7);
        assert_eq!(solid.reference_nodes_m(), NODES);
        assert_eq!(solid.reference_tetrahedra(), TETS);
    });
}

#[test]
fn heterogeneous_elements_and_rigid_motion_preserve_distinct_plastic_volumes() {
    with_cx(|cx, _| {
        // A shared face joins reference volumes 1/6 and 1/3. The second
        // material is pressure independent and cannot compact plastically.
        let nodes = [NODES[0], NODES[1], NODES[2], NODES[3], [1.0; 3]];
        let tets = [[0, 1, 2, 3], [4, 2, 1, 3]];
        let materials = [law(A), law(0.0)];
        let solid = FinitePorousTetSolid::try_new(
            &nodes,
            &tets,
            FinitePlasticLawField::PerElement(&materials),
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let stretch = 0.85;
        let f = [stretch, 0.0, 0.0, 0.0, stretch, 0.0, 0.0, 0.0, stretch];
        let mut bounded = controls();
        // A finite rotation of a prestressed body has a signed endpoint
        // quadrature error. Permit it explicitly, while checking heat is zero.
        bounded.maximum_endpoint_work_defect_j = 100.0;
        let compacted = solid
            .advance(
                &affine(&nodes, f, [0.0; 3]),
                &[],
                [0.0; 3],
                HyperTetSettings::default(),
                bounded,
                cx,
            )
            .unwrap();
        let (alpha, jp, tau) = hydrostatic(stretch);
        near(compacted.solid.element_states()[0].alpha, alpha, 3e-8);
        near(
            materials[0]
                .plastic_jacobian(&compacted.solid.element_states()[0])
                .unwrap(),
            jp,
            3e-8,
        );
        assert_eq!(
            compacted.solid.element_states()[1],
            materials[1].initial_state()
        );
        near(
            compacted.energy.plastic_dissipation_j,
            Y * alpha / 6.0,
            4e-7,
        );
        let elastic_tau = 3.0 * K * stretch.ln();
        let storage = (0.5 * tau * tau / K + 0.5 * H * alpha * alpha) / 6.0
            + 0.5 * elastic_tau * elastic_tau / K / 3.0;
        near(compacted.solid.stored_energy_j(), storage, 5e-7);
        let (sin, cos) = 0.7_f64.sin_cos();
        let rotation = [cos, -sin, 0.0, sin, cos, 0.0, 0.0, 0.0, 1.0];
        let rotated_f = rotation.map(|value| stretch * value);
        let rotated = compacted
            .solid
            .advance(
                &affine(&nodes, rotated_f, [0.4, -0.7, 1.2]),
                &[],
                [0.0; 3],
                HyperTetSettings::default(),
                bounded,
                cx,
            )
            .unwrap();
        assert_eq!(
            rotated.solid.element_states(),
            compacted.solid.element_states()
        );
        near(rotated.energy.plastic_dissipation_j, 0.0, 0.0);
        near(rotated.solid.stored_energy_j(), storage, 5e-7);
        assert!(rotated.energy.endpoint_work_defect_j < -1.0);
        for element in 0..2 {
            let old_p = compacted.solid.element_first_piola_pa()[element];
            let new_p = rotated.solid.element_first_piola_pa()[element];
            for i in 0..9 {
                let expected = (0..3)
                    .map(|j| rotation[3 * (i / 3) + j] * old_p[3 * j + i % 3])
                    .sum();
                near(new_p[i], expected, 5e-7);
                near(
                    rotated.solid.element_deformation()[element][i],
                    rotated_f[i],
                    3e-8,
                );
            }
        }
    });
}

#[test]
fn refused_work_and_cancellation_leave_compacted_state_retryable() {
    with_cx(|cx, gate| {
        let material = law(A);
        let initial = FinitePorousTetSolid::try_new(
            &NODES,
            &TETS,
            FinitePlasticLawField::Uniform(&material),
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let f = [0.9, 0.0, 0.0, 0.0, 0.9, 0.0, 0.0, 0.0, 0.9];
        let solid = initial
            .advance(
                &affine(&NODES, f, [0.0; 3]),
                &[],
                [0.0; 3],
                HyperTetSettings::default(),
                controls(),
                cx,
            )
            .unwrap()
            .solid;
        let before_state = solid.element_states().to_vec();
        let before_f = solid.element_deformation().to_vec();
        let before_u = solid.displacement_m().to_vec();
        let before_p = solid.element_first_piola_pa().to_vec();
        let before_storage = solid.stored_energy_j();
        let target = affine(
            &NODES,
            [0.8, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 0.8],
            [0.0; 3],
        );
        let mut strict = controls();
        strict.maximum_endpoint_work_defect_j = 0.0;
        let refused = solid
            .advance(
                &target,
                &[],
                [0.0; 3],
                HyperTetSettings::default(),
                strict,
                cx,
            )
            .unwrap_err();
        assert!(matches!(
            refused,
            FinitePorousTetError::Solid(FinitePlasticTetError::EnergyGate {
                quantity: "porous signed endpoint integration defect",
                ..
            })
        ));
        let retry = solid
            .advance(
                &target,
                &[],
                [0.0; 3],
                HyperTetSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        let repeat = solid
            .advance(
                &target,
                &[],
                [0.0; 3],
                HyperTetSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(retry.equilibrium, repeat.equilibrium);
        assert_eq!(retry.energy, repeat.energy);
        assert_eq!(retry.solid.element_states(), repeat.solid.element_states());
        assert!(retry.solid.element_states()[0].alpha > before_state[0].alpha);
        gate.request();
        let cancelled = solid
            .advance(
                &target,
                &[],
                [0.0; 3],
                HyperTetSettings::default(),
                controls(),
                cx,
            )
            .unwrap_err();
        assert!(matches!(
            cancelled,
            FinitePorousTetError::Solid(FinitePlasticTetError::Native(HyperTetError::Cancelled))
        ));
        assert_eq!(solid.element_states(), before_state);
        assert_eq!(solid.element_deformation(), before_f);
        assert_eq!(solid.displacement_m(), before_u);
        assert_eq!(solid.element_first_piola_pa(), before_p);
        near(solid.stored_energy_j(), before_storage, 0.0);
    });
}

#[test]
fn law_count_and_physical_endpoint_controls_are_admitted_before_state_changes() {
    with_cx(|cx, _| {
        let material = law(A);
        let no_laws: [FinitePorous; 0] = [];
        let refused = FinitePorousTetSolid::try_new(
            &NODES,
            &TETS,
            FinitePlasticLawField::PerElement(&no_laws),
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap_err();
        assert!(matches!(
            refused,
            FinitePorousTetError::Solid(FinitePlasticTetError::Native(
                HyperTetError::LawCountMismatch {
                    expected: 1,
                    got: 0
                }
            ))
        ));
        let solid = FinitePorousTetSolid::try_new(
            &NODES,
            &TETS,
            FinitePlasticLawField::Uniform(&material),
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let multiple = HyperTetSettings {
            load_steps: 2,
            ..HyperTetSettings::default()
        };
        assert!(
            solid
                .advance(&SUPPORTS, &[], [0.0; 3], multiple, controls(), cx)
                .is_err()
        );
        let mut nonfinite = controls();
        nonfinite.yield_tolerance_pa = f64::NAN;
        assert!(
            solid
                .advance(
                    &SUPPORTS,
                    &[],
                    [0.0; 3],
                    HyperTetSettings::default(),
                    nonfinite,
                    cx
                )
                .is_err()
        );
        assert_eq!(solid.element_states(), [material.initial_state()]);
        assert_eq!(solid.displacement_m(), [[0.0; 3]; 4]);
    });
}
