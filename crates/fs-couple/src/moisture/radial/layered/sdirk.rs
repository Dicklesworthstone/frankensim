//! The shared layered SDIRK kernel on native heterogeneous curved shells.
use super::*;
use crate::moisture::layered::MoistureLayeredSdirk2Report;

impl MoistureLayeredRadial {
    /// Second-order source-integrated transport in a coated sphere or cylinder.
    ///
    /// Exact radial paths, unequal original dry masses and material-specific
    /// sorption stay in the shared layered activity/storage kernel at both
    /// stages. The final boundary and interface accounts contain both weighted
    /// face transfers; the first-stage account is diagnostic, not an extra debit.
    /// Reservoir activities, temperature and reference geometry stay fixed.
    /// Inadmissible stages refuse atomically with ActivityStageRange instead of
    /// clipping or selecting another time integrator. Ordinary advance is BE.
    pub fn advance_sdirk2(
        &self,
        dt: f64,
        faces: &[MoistureLayeredFace; 2],
        step: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        cx: &Cx<'_>,
    ) -> Result<(Self, MoistureLayeredSdirk2Report), MoistureError> {
        poll(cx)?;
        let reservoir = self.sdirk_reservoir(faces)?;
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
                        system,
                        cx,
                    )
                },
                cx,
            )
        })?;
        let (next, endpoint) = self.accept_sdirk_candidate(candidate, cx)?;
        Ok((
            next,
            MoistureLayeredSdirk2Report {
                endpoint,
                first_stage,
            },
        ))
    }

    /// Physical first-stage face loads, including each unequal shell and the
    /// exact curved interface path. Only the native converged state is used.
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
                "radial layered SDIRK load must retain the original material and geometry",
            ));
        }
        let reservoir = self.sdirk_reservoir(faces)?;
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
                        system,
                        cx,
                    )
                },
                cx,
            )
        })
    }

    /// Complete the native second stage against a separately solved finite
    /// reservoir activity, retaining the original shell storage and geometry.
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
        let reservoir = self.sdirk_reservoir(faces)?;
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
                        system,
                        cx,
                    )
                },
                cx,
            )
        })?;
        self.accept_sdirk_candidate(candidate, cx)
    }

    fn sdirk_reservoir(
        &self,
        faces: &[MoistureLayeredFace; 2],
    ) -> Result<[Option<f64>; 2], MoistureError> {
        if self.boundary_areas[0] == 0.0 && !matches!(faces[0], MoistureLayeredFace::Sealed) {
            return Err(MoistureError::Invalid(
                "a solid radial center has zero area and must remain sealed",
            ));
        }
        Ok(std::array::from_fn(|side| match faces[side] {
            _ if self.exterior_transfer[side] == MoistureActivityTransfer::Sealed => None,
            MoistureLayeredFace::Sealed => None,
            MoistureLayeredFace::PrescribedRelativeHumidity(a) => Some(a),
        }))
    }

    fn with_sdirk_storage<T>(
        &self,
        operation: impl FnOnce(LayeredStorage<'_>) -> Result<T, MoistureError>,
    ) -> Result<T, MoistureError> {
        let original: Vec<_> = self
            .layers
            .iter()
            .flat_map(|layer| layer.radial.moisture.iter().copied())
            .collect();
        let dry_masses: Vec<_> = self
            .layers
            .iter()
            .flat_map(|layer| layer.radial.dry_masses_kg.iter().copied())
            .collect();
        let sources: Vec<_> = self
            .layers
            .iter()
            .map(MoistureRadialLayer::source)
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

    fn accept_sdirk_candidate(
        &self,
        candidate: kernel::LayeredCandidate,
        cx: &Cx<'_>,
    ) -> Result<(Self, MoistureLayeredReport), MoistureError> {
        let mut next = self.clone();
        next.activity = candidate.activity;
        let mut offset = 0;
        for layer in &mut next.layers {
            poll(cx)?;
            let n = layer.radial.moisture.len();
            layer
                .radial
                .moisture
                .copy_from_slice(&candidate.moisture[offset..offset + n]);
            layer.radial.time_s = candidate.report.interval_s[1];
            offset += n;
        }
        poll(cx)?;
        Ok((next, candidate.report))
    }
}
