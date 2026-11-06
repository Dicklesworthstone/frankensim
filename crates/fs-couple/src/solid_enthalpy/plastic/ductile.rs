//! Ductile initiation on the original spatial thermo-mechanical trajectory.
//!
//! The native thermal/plastic owner resolves and balances the complete physical
//! interval before its actual stress, plastic history and caloric temperatures
//! enter the shared initiation criterion. The same duration supplies inertia,
//! material kinetics and ductile plastic-strain rate. Native J2 and complete
//! Chaboche histories remain unchanged; their physical and reversible heat is
//! neither recomputed nor repartitioned here.
//!
//! Original source references and both endpoint material bundles are admitted.
//! A source locus frozen at one temperature cannot qualify another temperature.
//! Authored numerical Johnson--Cook/table laws retain their explicit domains
//! and do not acquire material-card authority from this coupling.
//!
//! Every input is borrowed. A failure or initiation crossing returns none of
//! the proposed heat, observer, mechanics, motion or ductility state. Callers
//! refine that complete original interval to locate initiation. Endpoint-held
//! accumulation remains first order; this gate adds no post-initiation damage,
//! fracture heat, topology change or separated/contacting fragments.

use super::dynamics::{ThermoplasticDynamics, ThermoplasticDynamicsCandidate};
use super::*;
use fs_solid::ductile::{
    DuctileAccepted, DuctileAdvanceError, DuctileHistory, DuctileLawField, DuctileMaterialBinding,
    DuctileSolidError, DuctileStepConditions, DuctileTemperatures, PlasticStrainHistory,
};

impl<'a, S: ThermoplasticSourceFamily<'a>> ThermoplasticSolid<'a, S>
where
    <S::Material as SmallStrainLaw>::State: PlasticStrainHistory,
{
    /// Read temperature from the same admitted caloric object and retained
    /// cell enthalpy. No material database query, ambient-temperature substitute
    /// or independent temperature state is introduced.
    fn ductile_temperatures<L: DuctileMaterialBinding, E>(
        &self,
        history: &DuctileHistory<'_, L>,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, DuctileAdvanceError<ThermoplasticError<E>>> {
        cx.checkpoint().map_err(|_| DuctileSolidError::Cancelled)?;
        if history.points().len() != self.cells.len() {
            return Err(DuctileSolidError::Input(
                "one ductile history per original thermoplastic tetrahedron required",
            )
            .into());
        }
        let laws = history.laws();
        let mut temperatures = Vec::with_capacity(self.cells.len());
        for (point, bound) in self.cells.iter().enumerate() {
            if point % 256 == 0 {
                cx.checkpoint().map_err(|_| DuctileSolidError::Cancelled)?;
            }
            let law = match laws {
                DuctileLawField::Uniform(law) => law,
                DuctileLawField::PerPoint(laws) => &laws[point],
            };
            let source = bound.source.thermal_source();
            if !law.matches_material(source.reference().material()) {
                return Err(DuctileSolidError::Input(
                    "ductile calibration must match the original thermal/mechanical source reference",
                ).into());
            }
            let state = source
                .caloric_source()
                .curve()
                .state_at_specific_enthalpy(bound.accepted.specific_h_j_kg)
                .map_err(|error| {
                    DuctileAdvanceError::Mechanics(ThermoplasticError::Material(
                        FeedbackError::Caloric(error),
                    ))
                })?;
            temperatures.push(state.temperature_k());
        }
        Ok(temperatures)
    }

    fn admit_ductile_history<L: DuctileMaterialBinding, E>(
        &self,
        history: &DuctileHistory<'_, L>,
        duration_s: f64,
        cx: &Cx<'_>,
    ) -> Result<(), DuctileAdvanceError<ThermoplasticError<E>>> {
        let temperatures = self.ductile_temperatures(history, cx)?;
        history.validate_material_history(
            &self.solution.element_states,
            &self.materials,
            DuctileStepConditions {
                duration_s,
                temperatures_k: DuctileTemperatures::PerPoint(&temperatures),
            },
            cx,
        )?;
        Ok(())
    }

    fn check_ductile_candidate<'d, L: DuctileMaterialBinding, E>(
        &self,
        candidate: &Self,
        history: &DuctileHistory<'d, L>,
        duration_s: f64,
        cx: &Cx<'_>,
    ) -> Result<DuctileHistory<'d, L>, DuctileAdvanceError<ThermoplasticError<E>>> {
        let temperatures = candidate.ductile_temperatures(history, cx)?;
        history
            .check_material_tet_candidate(
                &self.solution.element_states,
                &self.materials,
                &candidate.solution,
                &candidate.materials,
                DuctileStepConditions {
                    duration_s,
                    temperatures_k: DuctileTemperatures::PerPoint(&temperatures),
                },
                cx,
            )
            .map_err(DuctileAdvanceError::Ductility)
    }

    /// Solve one original coupled interval, then admit every ductile material
    /// point at its actual caloric endpoint. Accept the complete returned
    /// `mechanics` candidate and `ductility` together. A refusal preserves the
    /// original thermal state, observer, displacement and plastic history.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_ductility<'t, 'd, O: DepositionObserver, L: DuctileMaterialBinding>(
        &self,
        history: &DuctileHistory<'d, L>,
        thermal: &DepositionSimulation<'t>,
        observer: &O,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        thermal_config: StepConfig,
        config: ThermoplasticConfig,
        cx: &Cx<'_>,
    ) -> Result<
        DuctileAccepted<'d, ThermoplasticCandidate<'a, 't, O, S>, L>,
        DuctileAdvanceError<ThermoplasticError<O::Error>>,
    > {
        let duration_s = target_s - self.time_s;
        self.admit_ductile_history(history, duration_s, cx)?;
        let mechanics = self
            .advance(
                thermal,
                observer,
                target_s,
                prescribed_m,
                nodal_forces_n,
                settings,
                thermal_config,
                config,
                cx,
            )
            .map_err(DuctileAdvanceError::Mechanics)?;
        let ductility = self.check_ductile_candidate(mechanics.solid(), history, duration_s, cx)?;
        Ok(DuctileAccepted {
            mechanics,
            ductility,
        })
    }
}

impl<'a, S: ThermoplasticSourceFamily<'a>> ThermoplasticDynamics<'a, S>
where
    <S::Material as SmallStrainLaw>::State: PlasticStrainHistory,
{
    /// Native momentum, reciprocal material heat and initiation share one
    /// physical interval. Kinetic and constitutive numerical losses retain
    /// their original accounts. No candidate state escapes an initiation or
    /// source/domain refusal, including a candidate observer or velocity.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_ductility<'t, 'd, O: DepositionObserver, L: DuctileMaterialBinding>(
        &self,
        history: &DuctileHistory<'d, L>,
        thermal: &DepositionSimulation<'t>,
        observer: &O,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        thermal_config: StepConfig,
        config: ThermoplasticConfig,
        cx: &Cx<'_>,
    ) -> Result<
        DuctileAccepted<'d, ThermoplasticDynamicsCandidate<'a, 't, O, S>, L>,
        DuctileAdvanceError<ThermoplasticError<O::Error>>,
    > {
        let duration_s = target_s - self.time_s();
        self.solid()
            .admit_ductile_history(history, duration_s, cx)?;
        let mechanics = self
            .advance(
                thermal,
                observer,
                target_s,
                prescribed_m,
                nodal_forces_n,
                settings,
                thermal_config,
                config,
                cx,
            )
            .map_err(DuctileAdvanceError::Mechanics)?;
        let ductility = self.solid().check_ductile_candidate(
            mechanics.dynamics().solid(),
            history,
            duration_s,
            cx,
        )?;
        Ok(DuctileAccepted {
            mechanics,
            ductility,
        })
    }
}
