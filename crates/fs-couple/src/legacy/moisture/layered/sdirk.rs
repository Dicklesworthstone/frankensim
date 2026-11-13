//! Second-order evolution of original heterogeneous slab water inventories.
use super::*;

/// Both source-admitted stages of one second-order layered physical interval.
#[derive(Clone, Debug)]
pub struct MoistureLayeredSdirk2Report {
    /// Full original-to-final storage, boundary and interface water accounts.
    /// Each transfer already contains both RK-weighted stage contributions;
    /// surface activities, material source points and residuals describe the
    /// final stage, while interval/storage begin at the original accepted state.
    pub endpoint: MoistureLayeredReport,
    /// Diagnostic state at t0+gamma*dt with gamma*dt times its face rates.
    /// Do not add these water accounts to the complete endpoint account again.
    pub first_stage: MoistureLayeredReport,
}

impl MoistureLayeredSlab {
    /// Alexander SDIRK2 with the original nonlinear storage and interface laws.
    ///
    /// Both stages retain the actual material-specific M(a), source-integrated
    /// D(M), dry cell masses and explicit surface/interface resistances. The
    /// first face rates supply physical kg loads to the original storage of
    /// stage two; neither capacities nor activities are advanced as substitutes
    /// for water. Temperature, geometry and exterior reservoirs stay fixed.
    ///
    /// L-stability does not imply monotonicity. An inadmissible stage returns
    /// ActivityStageRange; reduce the interval rather than clipping the state
    /// or silently replacing this method with backward Euler. This fixed-
    /// geometry method does not add swelling feedback or a temporal certificate.
    pub fn advance_sdirk2(
        &self,
        dt: f64,
        faces: &[MoistureLayeredFace; 2],
        step: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        cx: &Cx<'_>,
    ) -> Result<(Self, MoistureLayeredSdirk2Report), MoistureError> {
        poll(cx)?;
        let reservoir = self.sdirk_reservoir(faces);
        let (candidate, first_stage) = self.with_sdirk_storage(|state| {
            kernel::sdirk2(
                state,
                dt,
                reservoir,
                step,
                nonlinear,
                |implicit_dt, profiles, activity, reference, system, cx| {
                    self.assemble_paths(
                        implicit_dt,
                        reservoir,
                        profiles,
                        activity,
                        reference,
                        None,
                        system,
                        cx,
                    )
                },
                cx,
            )
        })?;
        let next = self.with_profile(
            &candidate.activity,
            &candidate.moisture,
            candidate.report.interval_s[1],
            cx,
        )?;
        poll(cx)?;
        Ok((
            next,
            MoistureLayeredSdirk2Report {
                endpoint: candidate.report,
                first_stage,
            },
        ))
    }

    /// First-stage face loads for the original finite-reservoir driver.
    /// Capture once after its activity root converges, without a second solve.
    pub(in crate::moisture) fn sdirk_load(
        &self,
        dt: f64,
        first: &Self,
        faces: &[MoistureLayeredFace; 2],
        cx: &Cx<'_>,
    ) -> Result<kernel::StageLoad, MoistureError> {
        poll(cx)?;
        if first.definition_identity() != self.definition_identity() {
            return Err(MoistureError::Invalid(
                "layered SDIRK load must retain the original material and geometry",
            ));
        }
        let reservoir = self.sdirk_reservoir(faces);
        self.with_sdirk_storage(|state| {
            kernel::capture_sdirk_load(
                state,
                dt,
                &first.activity,
                first.time_s(),
                reservoir,
                |implicit_dt, profiles, activity, reference, system, cx| {
                    self.assemble_paths(
                        implicit_dt,
                        reservoir,
                        profiles,
                        activity,
                        reference,
                        None,
                        system,
                        cx,
                    )
                },
                cx,
            )
        })
    }

    /// Finish the native second stage using this original storage and the
    /// independently solved endpoint reservoir activity. Reported transfers
    /// include both RK stage weights exactly once.
    pub(in crate::moisture) fn advance_sdirk_stage(
        &self,
        dt: f64,
        faces: &[MoistureLayeredFace; 2],
        known: &kernel::StageLoad,
        step: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        cx: &Cx<'_>,
    ) -> Result<(Self, MoistureLayeredReport), MoistureError> {
        poll(cx)?;
        let reservoir = self.sdirk_reservoir(faces);
        let candidate = self.with_sdirk_storage(|state| {
            kernel::finish_sdirk(
                state,
                dt,
                reservoir,
                step,
                nonlinear,
                known,
                |implicit_dt, profiles, activity, reference, system, cx| {
                    self.assemble_paths(
                        implicit_dt,
                        reservoir,
                        profiles,
                        activity,
                        reference,
                        None,
                        system,
                        cx,
                    )
                },
                cx,
            )
        })?;
        let next = self.with_profile(
            &candidate.activity,
            &candidate.moisture,
            candidate.report.interval_s[1],
            cx,
        )?;
        poll(cx)?;
        Ok((next, candidate.report))
    }

    fn sdirk_reservoir(&self, faces: &[MoistureLayeredFace; 2]) -> [Option<f64>; 2] {
        std::array::from_fn(|side| match faces[side] {
            _ if self.exterior_transfer[side] == MoistureActivityTransfer::Sealed => None,
            MoistureLayeredFace::Sealed => None,
            MoistureLayeredFace::PrescribedRelativeHumidity(a) => Some(a),
        })
    }

    // All stage paths borrow the exact same native material/storage views.
    // Neither a homogeneous field nor a replacement capacity is introduced.
    fn with_sdirk_storage<T>(
        &self,
        operation: impl FnOnce(LayeredStorage<'_>) -> Result<T, MoistureError>,
    ) -> Result<T, MoistureError> {
        let original: Vec<_> = self
            .layers
            .iter()
            .flat_map(|layer| layer.slab.moisture.iter().copied())
            .collect();
        let dry_masses: Vec<_> = self
            .layer_for_cell
            .iter()
            .map(|&layer| self.layers[layer].slab.dry_mass_per_cell_kg)
            .collect();
        let sources: Vec<_> = self
            .layers
            .iter()
            .map(|layer| LayeredSource {
                card: &layer.slab.card,
                transport: &layer.slab.transport,
                property: &layer.slab.diffusivity_property,
                sorption: &layer.sorption,
            })
            .collect();
        operation(LayeredStorage {
            sources: &sources,
            layer_for_cell: &self.layer_for_cell,
            activity: &self.activity,
            original: &original,
            dry_masses: &dry_masses,
            activity_bounds: self.activity_bounds,
            time_s: self.time_s(),
            dry_mass_kg: self.dry_mass_kg(),
            water_mass_kg: self.water_mass_kg(),
        })
    }
}
