//! Borrowed views of the original solved hydraulic owner.
//!
//! A quasi-steady rate is held over its explicitly supplied thermal interval;
//! it carries no initial/final momentum history. An inertial view retains its
//! actual native old, average and new rates. Both share source-law admission.

use super::*;

#[derive(Clone, Copy)]
pub(super) enum GraphFlow<'a> {
    Inertial(&'a InertialGraphReport),
    QuasiSteady {
        hydraulics: &'a GraphSolution,
        duration_s: f64,
    },
}

impl GraphFlow<'_> {
    pub(super) fn duration_s(self) -> f64 {
        match self {
            Self::Inertial(report) => report.duration_s,
            Self::QuasiSteady { duration_s, .. } => duration_s,
        }
    }

    pub(super) fn average(self, branch: usize) -> f64 {
        match self {
            Self::Inertial(report) => report.average_volume_flows_m3_s[branch],
            Self::QuasiSteady { hydraulics, .. } => hydraulics.branches[branch].flow.value(),
        }
    }

    pub(super) fn initial(self, branch: usize) -> f64 {
        match self {
            Self::Inertial(report) => report.initial_flows_m3_s[branch],
            Self::QuasiSteady { .. } => self.average(branch),
        }
    }

    pub(super) fn final_flow(self, branch: usize) -> f64 {
        match self {
            Self::Inertial(report) => report.final_flows_m3_s[branch],
            Self::QuasiSteady { .. } => self.average(branch),
        }
    }

    pub(super) fn envelope(self, branch: usize) -> f64 {
        match self {
            Self::Inertial(report) => report.speed_envelope_m3_s[branch],
            Self::QuasiSteady { .. } => self.average(branch).abs(),
        }
    }

    pub(super) fn pressure_drop(self, branch: usize) -> f64 {
        match self {
            Self::Inertial(report) => report.dissipative_pressure_drops_pa[branch],
            Self::QuasiSteady { hydraulics, .. } => hydraulics.branches[branch]
                .dissipative_pressure_drop
                .value(),
        }
    }

    pub(super) fn gravity_head(self, branch: usize) -> f64 {
        match self {
            Self::Inertial(report) => report.gravity_pressure_heads_pa[branch],
            Self::QuasiSteady { hydraulics, .. } => hydraulics
                .gravity
                .as_ref()
                .map_or(0.0, |gravity| gravity.pressure_heads_pa()[branch]),
        }
    }

    /// Admission precedes all indexed access, including heat-law evaluation.
    pub(super) fn admit(
        self,
        cx: &Cx<'_>,
        core: &BoundGraphCore<'_>,
    ) -> Result<(), SourceGraphError> {
        poll(cx)?;
        let count = core.laws.len();
        if !positive(self.duration_s()) {
            return Err(SourceGraphError::Invalid(
                "native source admission requires the positive thermal interval duration",
            ));
        }
        match self {
            Self::Inertial(report) => {
                if [
                    report.average_volume_flows_m3_s.len(),
                    report.initial_flows_m3_s.len(),
                    report.final_flows_m3_s.len(),
                    report.speed_envelope_m3_s.len(),
                    report.dissipative_pressure_drops_pa.len(),
                    report.gravity_pressure_heads_pa.len(),
                ]
                .into_iter()
                .any(|n| n != count)
                {
                    return Err(SourceGraphError::Invalid(
                        "complete native inertial source report required",
                    ));
                }
            }
            Self::QuasiSteady { hydraulics, .. } => {
                if hydraulics.branches.len() != count
                    || hydraulics.pressures.len() != core.template.node_count()
                    || hydraulics.pressures.iter().any(|p| !p.value().is_finite())
                {
                    return Err(SourceGraphError::Invalid(
                        "complete native quasi-steady hydraulic solution required",
                    ));
                }
                for (actual, expected) in hydraulics.branches.iter().zip(core.template.branches()) {
                    poll(cx)?;
                    if actual.from != expected.from
                        || actual.to != expected.to
                        || actual.loss.name != expected.loss.name
                    {
                        return Err(SourceGraphError::Invalid(
                            "quasi-steady solution changed the original source graph topology",
                        ));
                    }
                }
                let gravity = hydraulics
                    .gravity
                    .as_ref()
                    .ok_or(SourceGraphError::Invalid(
                        "quasi-steady source solution lost its original gravity frame",
                    ))?;
                if gravity.pressure_heads_pa().len() != count
                    || gravity.positions_m() != core.positions.as_slice()
                    || gravity.acceleration_m_s2() != core.acceleration
                    || gravity.reference_transport().density_kg_m3()
                        != core.reference.density_kg_m3()
                    || gravity.maximum_relative_density_departure()
                        != core.limits.maximum_relative_density_departure
                {
                    return Err(SourceGraphError::Invalid(
                        "quasi-steady solution changed the original source gravity frame",
                    ));
                }
                core.inventory
                    .admit_gravity(cx, gravity)
                    .map_err(SourceGraphError::Inventory)?;
            }
        }
        Ok(())
    }
}
