//! Retained finite-strain material histories on the native tetrahedral solver.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_material::visco_hyper::{BranchEnergy, FiniteStrainBranch, ViscoHyperelastic};
use fs_solid::hyper3::visco::{
    FiniteViscoStepControls, FiniteViscoTetError, FiniteViscoTetSolid, ViscoelasticLawField,
};
use fs_solid::hyper3::{HyperTetProblem, HyperTetSettings, HyperelasticLawField};
use fs_solid::linear3::TetAssemblyBudget;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 32,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn law(mu: f64, branch_mu: f64, tau_s: f64) -> ViscoHyperelastic {
    let equilibrium = Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu,
            lambda: 2.0 * mu,
        },
        4.0,
    )
    .unwrap();
    let branches = if branch_mu == 0.0 {
        vec![]
    } else {
        vec![FiniteStrainBranch {
            energy: BranchEnergy::NeoHookean { mu_pa: branch_mu },
            tau_s,
        }]
    };
    ViscoHyperelastic::new(equilibrium, branches).unwrap()
}

fn settings() -> HyperTetSettings {
    HyperTetSettings {
        load_steps: 1,
        relative_tolerance: 1e-11,
        linear_tolerance: 1e-12,
        ..HyperTetSettings::default()
    }
}
fn controls() -> FiniteViscoStepControls {
    FiniteViscoStepControls {
        equilibrium_work_tolerance_j: 2e-7,
        maximum_endpoint_work_defect_j: 1e4,
    }
}

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const TET: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const PATCH: [[f64; 3]; 5] = [NODES[0], NODES[1], NODES[2], NODES[3], [0.25, 0.25, 0.25]];
const PATCH_TETS: [[usize; 4]; 4] = [[0, 1, 2, 4], [0, 1, 4, 3], [0, 4, 2, 3], [4, 1, 2, 3]];

fn affine(f: &[f64; 9], point: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| (0..3).map(|j| f[3 * i + j] * point[j]).sum::<f64>() - point[i])
}
fn boundary(f: &[f64; 9]) -> Vec<(usize, f64)> {
    NODES
        .iter()
        .enumerate()
        .flat_map(|(node, &point)| {
            affine(f, point)
                .into_iter()
                .enumerate()
                .map(move |(c, value)| (3 * node + c, value))
        })
        .collect()
}
fn axial_supports() -> Vec<(usize, f64)> {
    (0..12)
        .filter(|&dof| dof != 3)
        .map(|dof| (dof, 0.0))
        .collect()
}
fn close(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * a.abs().max(b.abs()).max(1.0),
        "{a:e} versus {b:e}"
    );
}

#[test]
fn affine_patch_relaxation_matches_material_point_with_free_interior_node() {
    with_cx(|cx| {
        let law = law(10_000.0, 30_000.0, 0.15);
        let mut solid = FiniteViscoTetSolid::try_new(
            &PATCH,
            &PATCH_TETS,
            ViscoelasticLawField::Uniform(&law),
            TetAssemblyBudget::default(),
            cx,
        )
        .unwrap();
        let f = [1.25, 0.03, 0.0, 0.0, 0.92, 0.0, 0.0, 0.0, 0.9];
        let prescribed = boundary(&f);
        let mut material_state = law.initial_state();
        let mut initial_stress = 0.0;
        for (step, dt) in [0.0, 0.05, 0.1, 0.2].into_iter().enumerate() {
            let native = law.step_response(&material_state, &f, dt).unwrap();
            let accepted = solid
                .advance(dt, &prescribed, &[], [0.0; 3], settings(), controls(), cx)
                .unwrap();
            for (node, &point) in PATCH.iter().enumerate() {
                for c in 0..3 {
                    close(
                        accepted.solid.displacement_m()[node][c],
                        affine(&f, point)[c],
                        2e-10,
                    );
                }
            }
            for (state, stress) in accepted
                .solid
                .element_states()
                .iter()
                .zip(accepted.solid.element_first_piola_pa())
            {
                for i in 0..9 {
                    close(stress[i], native.first_piola_pa[i], 2e-9);
                }
                for (actual, expected) in state
                    .inverse_viscous_cauchy_green
                    .iter()
                    .flatten()
                    .zip(native.state.inverse_viscous_cauchy_green.iter().flatten())
                {
                    close(*actual, *expected, 2e-10);
                }
            }
            close(
                accepted.energy.stored_after_j,
                native.state.stored / 6.0,
                2e-10,
            );
            close(
                accepted.energy.dissipated_j,
                native.dissipated_pa / 6.0,
                2e-9,
            );
            assert!(accepted.energy.dissipated_j >= 0.0);
            if step == 0 {
                initial_stress = native.first_piola_pa[0];
            } else {
                assert!(accepted.energy.native_loading_work_j.abs() < 2e-7);
                assert!(accepted.energy.internal_endpoint_work_j.abs() < 2e-7);
            }
            material_state = native.state;
            solid = accepted.solid;
        }
        assert!(solid.element_first_piola_pa()[0][0] < initial_stress * 0.7);
        close(solid.time_s(), 0.35, 1e-14);
    });
}

#[test]
fn maintained_load_creeps_and_books_physical_loss_separately_from_endpoint_error() {
    with_cx(|cx| {
        let law = law(10_000.0, 30_000.0, 0.15);
        let prescribed = axial_supports();
        let force = 900.0;
        let mut solid = FiniteViscoTetSolid::try_new(
            &NODES,
            &TET,
            ViscoelasticLawField::Uniform(&law),
            TetAssemblyBudget::default(),
            cx,
        )
        .unwrap();
        let jump = solid
            .advance(
                0.0,
                &prescribed,
                &[(3, force)],
                [0.0; 3],
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(jump.energy.dissipated_j, 0.0);
        let initial_extension = jump.solid.displacement_m()[1][0];
        assert!(initial_extension > 0.0);
        solid = jump.solid;
        let mut loss = 0.0;
        for _ in 0..20 {
            let before = solid.displacement_m()[1][0];
            let step = solid
                .advance(
                    0.03,
                    &prescribed,
                    &[(3, force)],
                    [0.0; 3],
                    settings(),
                    controls(),
                    cx,
                )
                .unwrap();
            assert!(step.solid.displacement_m()[1][0] > before);
            close(
                step.solid.element_first_piola_pa()[0][0] / 6.0,
                force,
                2e-10,
            );
            close(
                step.energy.dead_load_endpoint_work_j,
                force * (step.solid.displacement_m()[1][0] - before),
                1e-12,
            );
            close(
                step.energy.stored_after_j - step.energy.stored_before_j
                    + step.energy.dissipated_j
                    + step.energy.endpoint_work_defect_j,
                step.energy.dead_load_endpoint_work_j + step.energy.support_endpoint_work_j,
                2e-8,
            );
            assert!(step.energy.dissipated_j > 0.0);
            loss += step.energy.dissipated_j;
            solid = step.solid;
        }
        assert!(solid.displacement_m()[1][0] > 1.4 * initial_extension);
        close(solid.element_states()[0].dissipated / 6.0, loss, 1e-12);
        close(solid.time_s(), 0.6, 1e-14);
    });
}

#[test]
fn zero_branches_recover_native_hyperelastic_equilibrium_and_storage() {
    with_cx(|cx| {
        let law = law(10_000.0, 0.0, 0.15);
        let prescribed = axial_supports();
        let forces = [(3, 900.0)];
        let reference = HyperTetProblem {
            nodes_m: &NODES,
            tetrahedra: &TET,
            material: HyperelasticLawField::Uniform(&law.equilibrium),
            prescribed_m: &prescribed,
            nodal_forces_n: &forces,
            body_force_n_m3: [0.0; 3],
            budget: TetAssemblyBudget::default(),
            settings: settings(),
        }
        .solve(cx)
        .unwrap();
        let initial = FiniteViscoTetSolid::try_new(
            &NODES,
            &TET,
            ViscoelasticLawField::Uniform(&law),
            TetAssemblyBudget::default(),
            cx,
        )
        .unwrap();
        let step = initial
            .advance(
                0.3,
                &prescribed,
                &forces,
                [0.0; 3],
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        for (actual, expected) in step
            .solid
            .displacement_m()
            .iter()
            .flatten()
            .zip(reference.displacement_m.iter().flatten())
        {
            close(*actual, *expected, 1e-12);
        }
        close(step.energy.stored_after_j, reference.strain_energy_j, 1e-12);
        close(
            step.energy.incremental_potential_j,
            reference.strain_energy_j,
            1e-12,
        );
        assert_eq!(step.energy.dissipated_j, 0.0);
    });
}

#[test]
fn heterogeneous_histories_remain_compatible_and_rejected_candidates_leave_no_memory() {
    with_cx(|cx| {
        let soft = law(10_000.0, 30_000.0, 0.03);
        let stiff = law(20_000.0, 50_000.0, 0.3);
        let laws = [soft.clone(), stiff.clone(), soft.clone(), stiff];
        let initial = FiniteViscoTetSolid::try_new(
            &PATCH,
            &PATCH_TETS,
            ViscoelasticLawField::PerElement(&laws),
            TetAssemblyBudget::default(),
            cx,
        )
        .unwrap();
        let state_before = initial.element_states().to_vec();
        let f = [1.18, 0.02, 0.0, 0.0, 0.96, 0.0, 0.0, 0.0, 0.95];
        let prescribed = boundary(&f);
        let strict = FiniteViscoStepControls {
            maximum_endpoint_work_defect_j: 0.0,
            ..controls()
        };
        assert!(matches!(
            initial.advance(0.04, &prescribed, &[], [0.0; 3], settings(), strict, cx),
            Err(FiniteViscoTetError::EnergyGate {
                quantity: "signed endpoint integration defect",
                ..
            })
        ));
        assert_eq!(initial.element_states(), state_before);
        assert_eq!(initial.displacement_m(), &[[0.0; 3]; 5]);
        assert_eq!(initial.time_s(), 0.0);
        let step = initial
            .advance(0.04, &prescribed, &[], [0.0; 3], settings(), controls(), cx)
            .unwrap();
        let replay = initial
            .advance(0.04, &prescribed, &[], [0.0; 3], settings(), controls(), cx)
            .unwrap();
        assert_eq!(step.solid.element_states(), replay.solid.element_states());
        assert_eq!(step.solid.displacement_m(), replay.solid.displacement_m());
        assert_eq!(step.energy, replay.energy);
        let center = step.solid.displacement_m()[4];
        let homogeneous_center = affine(&f, PATCH[4]);
        assert!(
            center
                .iter()
                .zip(homogeneous_center)
                .any(|(a, b)| (a - b).abs() > 1e-4)
        );
        let held = step
            .solid
            .advance(0.06, &prescribed, &[], [0.0; 3], settings(), controls(), cx)
            .unwrap();
        for element in 0..4 {
            let actual = &held.solid.element_states()[element];
            let expected = laws[element]
                .step_response(
                    &step.solid.element_states()[element],
                    &actual.deformation_gradient,
                    0.06,
                )
                .unwrap();
            assert_eq!(*actual, expected.state);
            assert_eq!(
                held.solid.element_first_piola_pa()[element],
                expected.first_piola_pa
            );
        }
        assert!(
            initial
                .advance(-0.1, &prescribed, &[], [0.0; 3], settings(), controls(), cx)
                .is_err()
        );
        assert!(
            FiniteViscoTetSolid::try_new(
                &PATCH,
                &PATCH_TETS,
                ViscoelasticLawField::PerElement(&laws[..3]),
                TetAssemblyBudget::default(),
                cx
            )
            .is_err()
        );
    });
}
