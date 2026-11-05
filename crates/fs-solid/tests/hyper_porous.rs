//! Force-driven finite porous compaction through the existing tetrahedral solver.
//! The original mesh and caller-owned history survive loading/unloading.
//! Synthetic homogeneous coefficients do not establish measured foam accuracy.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::plastic::finite::porous::{
    FinitePorous, FinitePorousLimits, FinitePorousParameters,
};
use fs_solid::hyper3::{HyperTetProblem, HyperTetSettings, HyperelasticLawField};
use fs_solid::linear3::TetAssemblyBudget;

const MU: f64 = 80.0;
const K: f64 = 160.0;
const Y: f64 = 2.0;
const H: f64 = 8.0;
const A: f64 = 1.2;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}, allowance {tolerance:.3e}"
    );
}

#[test]
fn native_force_equilibrium_retains_compaction_and_elastic_reload_without_a_new_solver() {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 42,
                kernel_id: 98,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let nodes = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let tets = [[0, 1, 2, 3]];
        // Remove six rigid motions; all three axial stretches and all three
        // independent shears remain unknowns in the native equilibrium solve.
        let prescribed = [(0, 0.0), (1, 0.0), (2, 0.0), (4, 0.0), (5, 0.0), (8, 0.0)];
        let material = FinitePorous::new(
            FinitePorousParameters {
                shear_modulus_pa: MU,
                bulk_modulus_pa: K,
                yield_stress_pa: Y,
                hardening_modulus_pa: H,
                pressure_sensitivity: A,
            },
            FinitePorousLimits {
                minimum_jacobian: 0.2,
                maximum_elastic_log_strain: 0.5,
                maximum_plastic_log_strain: 2.0,
                maximum_accumulated_plastic_strain: 3.0,
            },
        )
        .unwrap();
        let beta = A / (1.0 + A * A / 9.0).sqrt();
        let peak_stretch = 0.8_f64;
        let alpha = (-3.0 * K * beta * peak_stretch.ln() - Y) / (H + K * beta * beta);
        let plastic_j = (-beta * alpha).exp();
        let unloaded_stretch = (-beta * alpha / 3.0).exp();
        let reload_stretch = (peak_stretch * unloaded_stretch).sqrt();
        let peak_kirchhoff = -(Y + H * alpha) / beta;
        let mut accepted_state = material.initial_state();
        let mut accepted_u = vec![0.0; 12];
        let mut accumulated_heat = 0.0;

        for (step, stretch) in [peak_stretch, unloaded_stretch, reload_stretch]
            .into_iter()
            .enumerate()
        {
            let expected_kirchhoff = match step {
                0 => peak_kirchhoff,
                1 => 0.0,
                _ => 0.5 * peak_kirchhoff,
            };
            let nominal = expected_kirchhoff / stretch;
            // V_ref=1/6 and gradients of the three axial vertices are unit
            // Cartesian vectors. These are dead referential force resultants.
            let forces = [(3, nominal / 6.0), (7, nominal / 6.0), (11, nominal / 6.0)];
            let saved_state = accepted_state.clone();
            let frozen = [material.frozen(&accepted_state).unwrap()];
            let problem = HyperTetProblem {
                nodes_m: &nodes,
                tetrahedra: &tets,
                material: HyperelasticLawField::PerElement(&frozen),
                prescribed_m: &prescribed,
                nodal_forces_n: &forces,
                body_force_n_m3: [0.0; 3],
                budget: TetAssemblyBudget::standard(),
                settings: HyperTetSettings::default(),
            };
            let equilibrium = problem.solve_from(&accepted_u, &cx).unwrap();
            assert!(equilibrium.steps[0].linear_iterations > 0);
            let candidate_u: Vec<_> = equilibrium
                .displacement_m
                .iter()
                .flatten()
                .copied()
                .collect();
            let fields = problem.element_responses(&candidate_u, &cx).unwrap();
            assert_eq!(fields.len(), 1);
            let field = &fields[0];
            let response = material
                .response(&field.deformation_gradient, &accepted_state)
                .unwrap();
            material
                .validate_committed(&field.deformation_gradient, &response.state, 2e-7)
                .unwrap();
            // Native Newton iterations and field observation are frozen-history
            // trials. The caller publishes only after every candidate is admitted.
            assert_eq!(accepted_state, saved_state);
            let expected_f = [stretch, 0.0, 0.0, 0.0, stretch, 0.0, 0.0, 0.0, stretch];
            for i in 0..9 {
                near(field.deformation_gradient[i], expected_f[i], 3e-8);
                let diagonal = i % 4 == 0;
                near(
                    field.first_piola_pa[i],
                    if diagonal { nominal } else { 0.0 },
                    3e-7,
                );
                near(
                    field.cauchy_pa[i],
                    if diagonal {
                        expected_kirchhoff / stretch.powi(3)
                    } else {
                        0.0
                    },
                    5e-7,
                );
            }
            near(field.reference_volume_m3, 1.0 / 6.0, 1e-15);
            near(field.current_volume_m3, stretch.powi(3) / 6.0, 3e-8);
            near(response.state.alpha, alpha, 3e-8);
            near(response.plastic_jacobian, plastic_j, 3e-8);
            let storage = 0.5 * expected_kirchhoff.powi(2) / K + 0.5 * H * alpha * alpha;
            let heat = if step == 0 { Y * alpha } else { 0.0 };
            near(response.stored_energy_j_m3, storage, 4e-7);
            near(response.plastic_dissipation_j_m3, heat, 3e-7);
            near(equilibrium.strain_energy_j, (storage + heat) / 6.0, 4e-7);
            near(field.potential_density_j_m3, storage + heat, 4e-7);
            if step > 0 {
                assert_eq!(
                    response.state, accepted_state,
                    "unload/reload must remain elastic"
                );
                near(response.plastic_increment, 0.0, 0.0);
            }
            accumulated_heat += response.plastic_dissipation_j_m3 / 6.0;
            accepted_state = response.state;
            accepted_u = candidate_u;
        }
        assert!(accepted_state.alpha > 0.5);
        assert!(material.plastic_jacobian(&accepted_state).unwrap() < 0.55);
        near(accumulated_heat, Y * alpha / 6.0, 4e-7);
    });
}
