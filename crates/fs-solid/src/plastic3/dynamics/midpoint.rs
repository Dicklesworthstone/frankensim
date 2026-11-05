//! Centered momentum with the existing P1 constitutive and Newton owners.
//!
//! The retaining material owner must supply a law whose response is the
//! work-conjugate mean stress over the accepted-to-proposed strain path, and
//! whose tangent differentiates that response with respect to the endpoint.
//! `update_state` still returns the accepted ENDPOINT history. Merely using an
//! arbitrary endpoint stress with midpoint inertia does not conserve energy.
//!
//! Loads are mean forces over this interval; prescribed values are endpoint
//! displacements. Reactions are the complete mean internal-plus-inertial force
//! minus applied force. Their work uses the actual displacement increment.
//! The original `MidpointInertia` owns `v_next = 2*du/dt-v_old`, the consistent
//! mass action and the kinetic-work receipt, with zero numerical kinetic loss.
//! Temporal accuracy still requires physical timestep refinement.
//!
//! At material switches a caller may explicitly enable numerical constitutive
//! columns using the original solver's fallible finite-difference adapter.
//! These are Newton corrections only; all forces, accepted histories and work
//! retain the actual material response. Physical-domain refusals still fail.

use super::*;
use crate::inertia::MidpointInertia;
use fs_solver::NonlinearProblem;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};

/// Disposable proposal; accepted displacement, velocity and history are borrowed.
#[derive(Clone, Debug)]
pub(crate) struct MidpointStep<S> {
    /// Endpoint displacement, strain and committed material history. The
    /// `element_stress` entries are the assembled RESPONSE stresses, usually
    /// path means. They are not endpoint stresses; the material owner must
    /// recover those independently before publishing its endpoint state.
    pub(crate) solution: SmallStrainTetSolution<S>,
    pub(crate) velocity_mps: Vec<[f64; 3]>,
    pub(crate) energy: InertialEnergy,
    /// Number of element assemblies that required numerical material columns.
    /// Each uses at most seven actual stress calls at frozen accepted history.
    pub(crate) numerical_tangent_evaluations: usize,
}

impl<L: SmallStrainLaw> SmallStrainTetProblem<'_, L> {
    /// Reuse an immutable mass bound to this reference geometry/inventory and
    /// the native shared Newton solve. The crate-level material owner admits
    /// path-mean stress, history, forces and full energy before publication.
    /// This deliberately is not a public midpoint claim for an arbitrary law.
    /// `difference_step_strain = None` requires every needed analytic tangent;
    /// `Some(step)` opts into bounded numerical columns when it is unavailable.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn advance_midpoint_with_mass(
        &self,
        displacement_m: &[f64],
        velocity_mps: &[[f64; 3]],
        element_states: &[L::State],
        mass: &P1ReferenceMass,
        dt_s: f64,
        difference_step_strain: Option<f64>,
        cx: &Cx<'_>,
    ) -> Result<MidpointStep<L::State>, SmallStrainTetError> {
        checkpoint(cx)?;
        if self.load_path != [1.0] {
            return Err(invalid(
                "one mean load at factor one is required per midpoint interval",
            ));
        }
        if difference_step_strain.is_some_and(|step| !step.is_finite() || step <= 0.0) {
            return Err(invalid("positive finite numerical strain step required"));
        }
        let fallback = difference_step_strain.map(|difference_step_strain| TangentFallback {
            difference_step_strain,
            evaluations: std::cell::Cell::new(0),
        });
        let inertia =
            MidpointInertia::new(mass, displacement_m, velocity_mps.as_flattened(), dt_s, cx)?;
        let solution = self.solve_candidate_with_tangent(
            Some((displacement_m, element_states)),
            Some(&inertia),
            fallback.as_ref(),
            cx,
        )?;
        let (velocity_mps, energy) = inertia.motion(solution.displacement_m.as_flattened(), cx)?;
        checkpoint(cx)?;
        Ok(MidpointStep {
            solution,
            velocity_mps,
            energy,
            numerical_tangent_evaluations: fallback.map_or(0, |value| value.evaluations.get()),
        })
    }
}

/// Assemble six local Newton columns, not a second material law or global
/// finite-difference matrix. Original source/domain and cancellation refusals
/// are preserved, including any refusal at a numerical stencil point.
pub(in crate::plastic3) fn numerical_tangent<L: SmallStrainLaw>(
    law: &L,
    strain: &Voigt,
    state: &L::State,
    element: usize,
    difference_step_strain: f64,
    cx: &Cx<'_>,
) -> Result<Tangent6, SmallStrainTetError> {
    let map_error = |error| match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(what) => SmallStrainTetError::InvalidInput {
            what: format!("element {element} numerical constitutive tangent: {what}"),
        },
    };
    let stress_at = |strain: &Voigt| {
        let stress = law.stress(strain, state);
        if stress.iter().any(|value| !value.is_finite()) {
            return Err(SmallStrainTetError::InvalidInput {
                what: format!("element {element} refused a numerical stress evaluation"),
            });
        }
        Ok(stress)
    };
    checkpoint(cx)?;
    let base_stress = stress_at(strain)?;
    let numerical = FallibleFiniteDifference::new(6, difference_step_strain, |values: &[f64]| {
        checkpoint(cx)?;
        // Cache only this immutable accepted-history/base-strain evaluation.
        // The native difference owner still selects and admits every stencil.
        let stress = if values == strain.as_slice() {
            base_stress
        } else {
            stress_at(&core::array::from_fn(|i| values[i]))?
        };
        Ok(stress.to_vec())
    })
    .map_err(map_error)?;
    let mut tangent = [[0.0; 6]; 6];
    for column in 0..6 {
        let mut direction = [0.0; 6];
        direction[column] = 1.0;
        let mut stress_change = [0.0; 6];
        numerical.jacobian_apply(strain, &direction, &mut stress_change);
        if let Some(error) = numerical.take_failure() {
            return Err(map_error(error));
        }
        for row in 0..6 {
            tangent[row][column] = stress_change[row];
        }
    }
    checkpoint(cx)?;
    Ok(tangent)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    use fs_material::IsotropicElastic;
    use fs_material::brittle::{
        BrittleCrackBand, BrittleCrackBandParameters, BrittleCrackBandSettings,
    };

    const NODES: [[f64; 3]; 4] = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];

    fn problem<'a, L: SmallStrainLaw>(
        law: &'a L,
        prescribed: &'a [(usize, f64)],
    ) -> SmallStrainTetProblem<'a, L> {
        SmallStrainTetProblem {
            nodes_m: &NODES,
            tetrahedra: &TETS,
            laws: SmallStrainLawField::Uniform(law),
            prescribed_m: prescribed,
            nodal_forces_n: &[],
            load_path: &[1.0],
            budget: TetAssemblyBudget::standard(),
            settings: IncrementSettings::default(),
        }
    }

    fn close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "{actual:e} != {expected:e}"
        );
    }

    #[test]
    fn native_inertial_rules_match_their_independent_elastic_mode_and_work() {
        let gate = CancelGate::new_clock_free();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 41,
                    kernel_id: 67,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            let young = 1.2e6;
            let poisson = 0.2;
            let band = BrittleCrackBand::new(
                IsotropicElastic::new(young, poisson, 0.03).unwrap(),
                BrittleCrackBandParameters {
                    tensile_strength_pa: 2000.0,
                    shear_strength_pa: 2000.0,
                    mode_i_fracture_energy_j_m2: 1.0,
                    mode_ii_fracture_energy_j_m2: 2.0,
                    bk_exponent: 1.5,
                },
                BrittleCrackBandSettings {
                    band_width_m: 1.0 / 3.0,
                    normal: [1.0, 0.0, 0.0],
                    maximum_strain: 0.03,
                    energy_tolerance_j_m3: 1e-7,
                },
            )
            .unwrap();
            let endpoint_law = band.small_strain_law();
            let mean_law = band.mean_small_strain_law();
            let fixed: Vec<_> = (0..12)
                .filter(|dof| *dof != 3)
                .map(|dof| (dof, 0.0))
                .collect();
            let old_q = 1e-4;
            let old_v = 0.01;
            let old_state = band
                .trial(&[old_q, 0.0, 0.0, 0.0, 0.0, 0.0], &band.initial_state())
                .unwrap()
                .state;
            let mut old_u = [0.0; 12];
            old_u[3] = old_q;
            let mut velocity = [[0.0; 3]; 4];
            velocity[1][0] = old_v;
            // Independent one-free-DOF reduction: V=1/6, M_ii=rho*V/10,
            // and K_ii=(lambda+2*mu)*V for this reference tetrahedron.
            let modulus = young * (1.0 - poisson) / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
            let stiffness = modulus / 6.0;
            for rho in [6.0, 12.0] {
                let mass =
                    P1ReferenceMass::new([(TETS[0], 1.0 / 6.0)].into_iter(), &[rho], 12, &cx)
                        .unwrap();
                let modal_mass = rho / 60.0;
                for dt in [0.0002, 0.001] {
                    let expected_be =
                        modal_mass * (old_q + dt * old_v) / (modal_mass + dt * dt * stiffness);
                    let expected_midpoint = ((modal_mass - 0.25 * dt * dt * stiffness) * old_q
                        + modal_mass * dt * old_v)
                        / (modal_mass + 0.25 * dt * dt * stiffness);
                    let endpoint = problem(&endpoint_law, &fixed)
                        .advance_backward_euler_with_mass(
                            &old_u,
                            &velocity,
                            &[old_state],
                            &mass,
                            dt,
                            &cx,
                        )
                        .unwrap();
                    let midpoint = problem(&mean_law, &fixed)
                        .advance_midpoint_with_mass(
                            &old_u,
                            &velocity,
                            &[old_state],
                            &mass,
                            dt,
                            None,
                            &cx,
                        )
                        .unwrap();
                    close(endpoint.solution.displacement_m[1][0], expected_be, 1e-14);
                    close(
                        midpoint.solution.displacement_m[1][0],
                        expected_midpoint,
                        1e-14,
                    );
                    close(
                        endpoint.velocity_mps[1][0],
                        (expected_be - old_q) / dt,
                        1e-11,
                    );
                    close(
                        midpoint.velocity_mps[1][0],
                        2.0 * (expected_midpoint - old_q) / dt - old_v,
                        1e-11,
                    );
                    close(
                        endpoint.solution.element_stress[0][0],
                        modulus * expected_be,
                        1e-8,
                    );
                    close(
                        midpoint.solution.element_stress[0][0],
                        0.5 * modulus * (old_q + expected_midpoint),
                        1e-8,
                    );
                    close(
                        endpoint.energy.kinetic_change_j
                            + endpoint.energy.numerical_loss_j
                            + stiffness * expected_be * (expected_be - old_q),
                        0.0,
                        1e-11,
                    );
                    close(
                        midpoint.energy.kinetic_change_j
                            + 0.5
                                * stiffness
                                * (expected_midpoint * expected_midpoint - old_q * old_q),
                        0.0,
                        1e-11,
                    );
                    assert!(endpoint.energy.numerical_loss_j > 0.0);
                    assert_eq!(midpoint.energy.numerical_loss_j, 0.0);
                    assert_eq!(midpoint.numerical_tangent_evaluations, 0);
                    close(
                        midpoint.solution.element_states[0].strain()[0],
                        expected_midpoint,
                        1e-14,
                    );
                    assert_eq!(midpoint.solution.element_states[0].damage(), 0.0);
                }
            }
            assert_eq!(old_state.strain()[0], old_q);
            assert_eq!(old_state.damage(), 0.0);
            assert_eq!(old_u[3], old_q);
            assert_eq!(velocity[1][0], old_v);
        });
    }
}
