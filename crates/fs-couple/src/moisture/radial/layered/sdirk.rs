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
        if self.boundary_areas[0] == 0.0 && !matches!(faces[0], MoistureLayeredFace::Sealed) {
            return Err(MoistureError::Invalid(
                "a solid radial center has zero area and must remain sealed",
            ));
        }
        let reservoir = std::array::from_fn(|side| match faces[side] {
            _ if self.exterior_transfer[side] == MoistureActivityTransfer::Sealed => None,
            MoistureLayeredFace::Sealed => None,
            MoistureLayeredFace::PrescribedRelativeHumidity(a) => Some(a),
        });
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
        let (candidate, first_stage) = kernel::sdirk2(
            LayeredStorage {
                sources: &sources,
                layer_for_cell: &self.layer_for_cell,
                activity: &self.activity,
                original: &original,
                dry_masses: &dry_masses,
                activity_bounds: self.activity_bounds,
                time_s: self.time_s(),
                dry_mass_kg: self.dry_mass_kg(),
                water_mass_kg: self.water_mass_kg(),
            },
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
        )?;
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
        Ok((
            next,
            MoistureLayeredSdirk2Report {
                endpoint: candidate.report,
                first_stage,
            },
        ))
    }
}
