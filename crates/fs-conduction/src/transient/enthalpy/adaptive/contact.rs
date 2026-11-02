//! Adaptive source-backed contact through the original enthalpy controller.
//!
//! Every coarse and fine endpoint solves the actual contact/radiation laws
//! together. The controller compares h, T and liquid fraction, commits only a
//! complete accepted pair, and integrates only that pair's external heat.

use fs_solver::Globalization;

use super::*;
use crate::interface::temperature::TemperatureContacts;
use crate::interface::temperature::enthalpy::{
    ContactEnthalpyNewtonConfig, check_coupled_endpoint,
};
use crate::transient::enthalpy::PointwiseBoundary;
use crate::{MaterialPointwiseAmbientRadiation, MaterialPointwiseEnclosureRadiation};

/// Contact controls added to the adaptive state's unchanged endpoint policy.
///
/// Newton/FGMRES controls and energy tolerance come from
/// [`AdaptiveEnthalpyConfig::endpoint`]. That Newton policy must use Armijo line
/// search, zero relative tolerance and an absolute target no looser than the
/// physical contact residual target below. No call widens those gates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AdaptiveContactEnthalpyConfig {
    /// Residual callbacks per endpoint, including initialization and rejected
    /// line-search/source-domain trials. Physical endpoint audits are separate.
    pub max_residual_evaluations_per_endpoint: usize,
    /// Original matching-pair/overlap-triangle contact tangent allowance.
    pub max_face_pairs: usize,
    /// Absolute norm of the complete physical nodal energy residual, J.
    pub residual_tolerance_j: f64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::transient::enthalpy::{EnthalpyBudget, EnthalpyEndpointEvaluation};
    use crate::{ConductionMesh, ConductivityModel, ScalarField, ThermalBoundaryBuilder};
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_blake3::ContentHash;
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
    use fs_rep_mesh::TetComplex;

    #[test]
    fn physical_endpoint_refusal_retries_without_publishing_comparison_history() {
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let gate = CancelGate::new_clock_free();
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 1,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            let mesh = ConductionMesh::new(
                TetComplex::from_tets(4, vec![[0, 1, 2, 3]]),
                vec![[0.; 3], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            )
            .unwrap();
            let curve = EquilibriumEnthalpyPhaseCurve::try_new(
                ContentHash([33; 32]),
                vec![
                    EnthalpyPhaseKnot {
                        specific_enthalpy_j_kg: 0.,
                        temperature_k: 300.,
                        liquid_mass_fraction: 0.,
                        bulk_density_kg_m3: 1.,
                    },
                    EnthalpyPhaseKnot {
                        specific_enthalpy_j_kg: 100.,
                        temperature_k: 400.,
                        liquid_mass_fraction: 1.,
                        bulk_density_kg_m3: 1.,
                    },
                ],
            )
            .unwrap();
            let storage = EnthalpyBackwardEuler::uniform(
                &cx,
                &mesh,
                &curve,
                1.,
                EnthalpyBudget {
                    max_vertices: 4,
                    max_elements: 1,
                },
            )
            .unwrap();
            let boundary = ThermalBoundaryBuilder::new(&mesh)
                .adiabatic_remainder()
                .finish()
                .unwrap();
            let material = ConductivityModel::isotropic_declared(1.).unwrap();
            let source = ScalarField::Uniform(0.);
            let problem = ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                element_materials: None,
                source: &source,
            };
            let endpoint = EnthalpyStepConfig {
                newton: fs_solver::NewtonKrylovConfig {
                    relative_tolerance: 0.,
                    ..Default::default()
                },
                max_newton_iterations: 8,
                energy_tolerance_j: 1e-8,
            };
            let config = AdaptiveEnthalpyConfig {
                enthalpy_absolute_tolerance_j_kg: 1.,
                temperature_absolute_tolerance_k: 1.,
                liquid_fraction_absolute_tolerance: 0.1,
                relative_tolerance: 0.,
                minimum_step_s: 0.1,
                maximum_step_s: 1.,
                endpoint,
            };
            let mut state = storage
                .adaptive_state(&cx, &[50.; 4], 0., 1., config)
                .unwrap();
            let physical = EnthalpyEndpointEvaluation {
                temperature: vec![350.; 4],
                residual_j: vec![1e-8, -1e-8, 0., 0.],
                residual_norm_j: 2_f64.sqrt() * 1e-8,
                free_dofs: 4,
                stored_energy_change_j: 0.,
                source_w: 0.,
                neumann_out_w: 0.,
                dirichlet_in_w: 0.,
                convective_out_w: 0.,
                radiation_out_w: 0.,
                energy_residual_j: 0.,
                convective_robin_fluxes: vec![],
            };
            let controls = ContactEnthalpyNewtonConfig {
                step: endpoint,
                max_residual_evaluations: 32,
                max_face_pairs: 1,
                residual_tolerance_j: 1e-8,
            };
            assert!(matches!(
                check_coupled_endpoint(&physical, controls),
                Err(EnthalpyError::EndpointResidual { .. })
            ));
            let mut calls = 0;
            let mut observed = 0;
            let report = state
                .advance_with(
                    &cx,
                    1.,
                    2,
                    |old, dt| {
                        calls += 1;
                        // Exercise a typed failed independent audit, not a fabricated
                        // Newton stall. All later responses are real equilibrium steps.
                        if calls == 1 {
                            check_coupled_endpoint(&physical, controls)?;
                        }
                        storage.advance(&cx, problem, None, old, dt, endpoint)
                    },
                    |trial| {
                        observed += 1;
                        assert_eq!(trial.start_time_s, 0.);
                        assert_eq!(trial.end_time_s, 0.5);
                        assert_eq!(trial.first_half.specific_enthalpy_j_kg, [50.; 4]);
                        assert_eq!(trial.second_half.specific_enthalpy_j_kg, [50.; 4]);
                        ControlFlow::Continue(())
                    },
                )
                .unwrap();
            assert_eq!(report.status, AdaptiveEnthalpyStatus::TrialLimit);
            assert_eq!(
                (report.rejected, report.accepted, report.endpoint_solves),
                (1, 1, 4)
            );
            assert_eq!((calls, observed), (4, 1));
            assert_eq!(state.time_s(), 0.5);
            assert_eq!(state.specific_enthalpy_j_kg(), [50.; 4]);
            assert_eq!(report.integrated_net_input_j, 0.);
            assert_eq!(report.energy_residual_j, 0.);
        });
    }
}

impl AdaptiveContactEnthalpyConfig {
    fn endpoint(
        self,
        step: EnthalpyStepConfig,
        max_trials: usize,
    ) -> Result<ContactEnthalpyNewtonConfig, EnthalpyError> {
        if self.max_residual_evaluations_per_endpoint == 0
            || self.max_face_pairs == 0
            || !self.residual_tolerance_j.is_finite()
            || self.residual_tolerance_j <= 0.0
            || step.newton.relative_tolerance != 0.0
            || step.newton.absolute_tolerance > self.residual_tolerance_j
            || !matches!(step.newton.globalization, Globalization::LineSearch(_))
        {
            return Err(EnthalpyError::InvalidInput(
                "adaptive contact needs positive work/joule gates, line search, zero endpoint relative tolerance and an absolute target no looser than the physical target",
            ));
        }
        // Every attempted comparison has at most three endpoint solves. Bound
        // rejected/coarse work as well as accepted work before changing state.
        let endpoints = max_trials
            .checked_mul(3)
            .ok_or(EnthalpyError::InvalidInput(
                "adaptive contact work allowance overflows",
            ))?;
        endpoints
            .checked_mul(self.max_residual_evaluations_per_endpoint)
            .and_then(|_| endpoints.checked_mul(step.max_newton_iterations))
            .and_then(|attempts| attempts.checked_mul(step.newton.linear_restart))
            .and_then(|columns| columns.checked_mul(step.newton.max_linear_cycles))
            .ok_or(EnthalpyError::InvalidInput(
                "adaptive contact work allowance overflows",
            ))?;
        Ok(ContactEnthalpyNewtonConfig {
            step,
            max_residual_evaluations: self.max_residual_evaluations_per_endpoint,
            max_face_pairs: self.max_face_pairs,
            residual_tolerance_j: self.residual_tolerance_j,
        })
    }
}

#[derive(Clone, Copy)]
enum ContactRadiation<'a> {
    None,
    Ambient(&'a PointwiseAmbientRadiation<'a>),
    MaterialAmbient(&'a MaterialPointwiseAmbientRadiation<'a>),
    Enclosure(&'a PointwiseEnclosureRadiation<'a>),
    MaterialEnclosure(&'a MaterialPointwiseEnclosureRadiation<'a>),
}

impl<'a> ContactRadiation<'a> {
    fn boundary(self) -> Option<PointwiseBoundary<'a>> {
        match self {
            Self::None => None,
            Self::Ambient(owner) => Some(PointwiseBoundary::Ambient(owner)),
            Self::MaterialAmbient(owner) => Some(PointwiseBoundary::MaterialAmbient(owner)),
            Self::Enclosure(owner) => Some(PointwiseBoundary::Enclosure(owner)),
            Self::MaterialEnclosure(owner) => Some(PointwiseBoundary::MaterialEnclosure(owner)),
        }
    }
}

impl AdaptiveEnthalpyState<'_, '_, '_> {
    /// Advance source-backed contact with adaptive primary-enthalpy time error.
    ///
    /// Contact geometry, pressure and pinned source laws stay fixed throughout
    /// this event window; each endpoint re-queries their actual temperatures.
    /// Sharp latent plateaus retain independent h and liquid-fraction gates.
    /// `max_trials` counts accepted and rejected comparisons; actual attempted
    /// endpoint solves, including comparison work, appear in the shared report.
    ///
    /// Nonlinear, local residual and energy refusals use the original
    /// controller's bounded step reduction. Source-support, binding and
    /// cancellation errors return their original refusal with the last
    /// accepted physical history unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to_with_temperature_contacts(
        &mut self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        contacts: &TemperatureContacts<'_, '_>,
        config: AdaptiveContactEnthalpyConfig,
        end_time_s: f64,
        max_trials: usize,
    ) -> Result<AdaptiveEnthalpyReport, EnthalpyError> {
        self.advance_to_with_temperature_contacts_observed(
            cx,
            problem,
            contacts,
            config,
            end_time_s,
            max_trials,
            |_| ControlFlow::Continue(()),
        )
    }

    /// Observe both accepted fine endpoints, with the same atomic acceptance,
    /// stop/resume and cancellation semantics as [`Self::advance_to_observed`].
    /// Each endpoint includes its independently audited contact heat reports.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to_with_temperature_contacts_observed(
        &mut self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        contacts: &TemperatureContacts<'_, '_>,
        config: AdaptiveContactEnthalpyConfig,
        end_time_s: f64,
        max_trials: usize,
        observe: impl FnMut(EnthalpyAcceptedTrial<'_>) -> ControlFlow<()>,
    ) -> Result<AdaptiveEnthalpyReport, EnthalpyError> {
        self.advance_contact(
            cx,
            problem,
            contacts,
            config,
            ContactRadiation::None,
            end_time_s,
            max_trials,
            observe,
        )
    }

    /// Adapt the same simultaneous contact/enthalpy solve with local P1
    /// ambient radiation. Its prepared owner must retain this exact mesh and
    /// boundary; convection and signed radiative heat remain separate accounts.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to_with_temperature_contacts_and_pointwise_ambient_radiation(
        &mut self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        contacts: &TemperatureContacts<'_, '_>,
        config: AdaptiveContactEnthalpyConfig,
        radiation: &PointwiseAmbientRadiation<'_>,
        end_time_s: f64,
        max_trials: usize,
    ) -> Result<AdaptiveEnthalpyReport, EnthalpyError> {
        self.advance_to_with_temperature_contacts_and_pointwise_ambient_radiation_observed(
            cx,
            problem,
            contacts,
            config,
            radiation,
            end_time_s,
            max_trials,
            |_| ControlFlow::Continue(()),
        )
    }

    /// Stream accepted contact/ambient-radiation fine pairs without retaining
    /// rejected trajectories or integrating their heat into accepted history.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to_with_temperature_contacts_and_pointwise_ambient_radiation_observed(
        &mut self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        contacts: &TemperatureContacts<'_, '_>,
        config: AdaptiveContactEnthalpyConfig,
        radiation: &PointwiseAmbientRadiation<'_>,
        end_time_s: f64,
        max_trials: usize,
        observe: impl FnMut(EnthalpyAcceptedTrial<'_>) -> ControlFlow<()>,
    ) -> Result<AdaptiveEnthalpyReport, EnthalpyError> {
        self.advance_contact(
            cx,
            problem,
            contacts,
            config,
            ContactRadiation::Ambient(radiation),
            end_time_s,
            max_trials,
            observe,
        )
    }

    /// Adapt source-backed contact and ambient emissivity in the same enthalpy
    /// solve. Each candidate refreshes both original material responses; the
    /// radiation tangent retains liquid-fraction feedback on latent plateaus.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to_with_temperature_contacts_and_material_ambient_radiation(
        &mut self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        contacts: &TemperatureContacts<'_, '_>,
        config: AdaptiveContactEnthalpyConfig,
        radiation: &MaterialPointwiseAmbientRadiation<'_>,
        end_time_s: f64,
        max_trials: usize,
    ) -> Result<AdaptiveEnthalpyReport, EnthalpyError> {
        self.advance_to_with_temperature_contacts_and_material_ambient_radiation_observed(
            cx, problem, contacts, config, radiation, end_time_s, max_trials,
            |_| ControlFlow::Continue(()),
        )
    }

    /// Observe accepted contact/material-radiation fine pairs and their actual
    /// nodal emissivity receipts. Rejected work cannot contribute heat or alter
    /// accepted history; all original endpoint and source gates still apply.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to_with_temperature_contacts_and_material_ambient_radiation_observed(
        &mut self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        contacts: &TemperatureContacts<'_, '_>,
        config: AdaptiveContactEnthalpyConfig,
        radiation: &MaterialPointwiseAmbientRadiation<'_>,
        end_time_s: f64,
        max_trials: usize,
        observe: impl FnMut(EnthalpyAcceptedTrial<'_>) -> ControlFlow<()>,
    ) -> Result<AdaptiveEnthalpyReport, EnthalpyError> {
        self.advance_contact(
            cx, problem, contacts, config, ContactRadiation::MaterialAmbient(radiation),
            end_time_s, max_trials, observe,
        )
    }

    /// Adapt local emission, reflected enclosure exchange and source-backed
    /// contact in one enthalpy system. Closed exchange remains internal;
    /// an open enclosure's explicit bath enters accepted external heat once.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to_with_temperature_contacts_and_pointwise_enclosure_radiation(
        &mut self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        contacts: &TemperatureContacts<'_, '_>,
        config: AdaptiveContactEnthalpyConfig,
        radiation: &PointwiseEnclosureRadiation<'_>,
        end_time_s: f64,
        max_trials: usize,
    ) -> Result<AdaptiveEnthalpyReport, EnthalpyError> {
        self.advance_to_with_temperature_contacts_and_pointwise_enclosure_radiation_observed(
            cx,
            problem,
            contacts,
            config,
            radiation,
            end_time_s,
            max_trials,
            |_| ControlFlow::Continue(()),
        )
    }

    /// Observe accepted fine endpoints with their original local-emission and
    /// enclosure reports. Contact and closed radiative exchange never become
    /// fictitious external input, even when the local fluxes are nonzero.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to_with_temperature_contacts_and_pointwise_enclosure_radiation_observed(
        &mut self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        contacts: &TemperatureContacts<'_, '_>,
        config: AdaptiveContactEnthalpyConfig,
        radiation: &PointwiseEnclosureRadiation<'_>,
        end_time_s: f64,
        max_trials: usize,
        observe: impl FnMut(EnthalpyAcceptedTrial<'_>) -> ControlFlow<()>,
    ) -> Result<AdaptiveEnthalpyReport, EnthalpyError> {
        self.advance_contact(
            cx,
            problem,
            contacts,
            config,
            ContactRadiation::Enclosure(radiation),
            end_time_s,
            max_trials,
            observe,
        )
    }

    /// Adapt source contact and evolving reflected emissivity together. Each
    /// trial includes the full h derivative of emission and reflection, even
    /// when melting changes emissivity at constant temperature.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to_with_temperature_contacts_and_material_enclosure_radiation(
        &mut self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        contacts: &TemperatureContacts<'_, '_>,
        config: AdaptiveContactEnthalpyConfig,
        radiation: &MaterialPointwiseEnclosureRadiation<'_>,
        end_time_s: f64,
        max_trials: usize,
    ) -> Result<AdaptiveEnthalpyReport, EnthalpyError> {
        self.advance_to_with_temperature_contacts_and_material_enclosure_radiation_observed(
            cx, problem, contacts, config, radiation, end_time_s, max_trials,
            |_| ControlFlow::Continue(()),
        )
    }

    /// Observe accepted material-enclosure/contact endpoints with their actual
    /// nodal source receipts and independently checked local energy balances.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to_with_temperature_contacts_and_material_enclosure_radiation_observed(
        &mut self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        contacts: &TemperatureContacts<'_, '_>,
        config: AdaptiveContactEnthalpyConfig,
        radiation: &MaterialPointwiseEnclosureRadiation<'_>,
        end_time_s: f64,
        max_trials: usize,
        observe: impl FnMut(EnthalpyAcceptedTrial<'_>) -> ControlFlow<()>,
    ) -> Result<AdaptiveEnthalpyReport, EnthalpyError> {
        self.advance_contact(
            cx, problem, contacts, config, ContactRadiation::MaterialEnclosure(radiation),
            end_time_s, max_trials, observe,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_contact(
        &mut self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        contacts: &TemperatureContacts<'_, '_>,
        controls: AdaptiveContactEnthalpyConfig,
        radiation: ContactRadiation<'_>,
        end: f64,
        max_trials: usize,
        observe: impl FnMut(EnthalpyAcceptedTrial<'_>) -> ControlFlow<()>,
    ) -> Result<AdaptiveEnthalpyReport, EnthalpyError> {
        poll(cx, 0)?;
        let storage = self.storage;
        let config = controls.endpoint(self.config.endpoint, max_trials)?;
        self.advance_with(
            cx,
            end,
            max_trials,
            |old, dt| {
                let (solved, actual) = storage.advance_with_temperature_contacts(
                    cx,
                    problem,
                    contacts,
                    old,
                    dt,
                    config,
                    radiation.boundary(),
                )?;
                let interfaces = Some(actual.frozen_interfaces());
                let candidate = &solved.specific_enthalpy_j_kg;
                let physical = match radiation {
                    ContactRadiation::None => storage.evaluate_endpoint(
                        cx,
                        problem,
                        interfaces,
                        old,
                        candidate,
                        dt,
                        config.step,
                        None,
                    ),
                    ContactRadiation::Ambient(owner) => storage
                        .evaluate_endpoint_with_pointwise_ambient_radiation(
                            cx,
                            problem,
                            interfaces,
                            old,
                            candidate,
                            dt,
                            config.step,
                            owner,
                        ),
                ContactRadiation::MaterialAmbient(owner) => storage
                    .evaluate_endpoint_with_material_ambient_radiation(
                        cx, problem, interfaces, old, candidate, dt, config.step, owner,
                    ),
                    ContactRadiation::Enclosure(owner) => storage
                        .evaluate_endpoint_with_pointwise_enclosure_radiation(
                            cx,
                            problem,
                            interfaces,
                            old,
                            candidate,
                            dt,
                            config.step,
                            owner,
                        ),
                    ContactRadiation::MaterialEnclosure(owner) => storage
                        .evaluate_endpoint_with_material_enclosure_radiation(
                            cx, problem, interfaces, old, candidate, dt, config.step, owner,
                        ),
                }?;
                check_coupled_endpoint(&physical, config)?;
                // All physical gates precede entry into the shared controller.
                // Its contact fluxes and pointwise reports already refer to this
                // exact accepted law; no discarded source-state clone is needed.
                Ok(solved)
            },
            observe,
        )
    }
}
