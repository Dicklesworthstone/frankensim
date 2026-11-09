//! Explicit ambient discretization on the original scan history.
//!
//! Declarations and event sampling are shared. Area-mean exchange and local P1
//! emission keep their distinct original solvers and physical reports.

use super::{EnthalpyHistoryRun, HeterogeneousEnthalpyError, RadiationRegionHistory,
    ScanRadiation, ScanningEnthalpyRun, conduction, invalid, patches_before, poll};
use crate::{AmbientRadiationConfig, AmbientRadiationPatch, ConductionError,
    PointwiseAmbientRadiation};
use fs_exec::Cx;

pub(super) enum ScanRadiationBinding {
    MeanAmbient(ScanRadiation),
    PointwiseAmbient(Vec<RadiationRegionHistory>),
}

impl ScanRadiationBinding {
    pub(super) fn patches(&self) -> &[RadiationRegionHistory] {
        match self {
            Self::MeanAmbient(radiation) => &radiation.patches,
            Self::PointwiseAmbient(patches) => patches,
        }
    }

    pub(super) fn mean_config(&self) -> Option<AmbientRadiationConfig> {
        match self {
            Self::MeanAmbient(radiation) => Some(radiation.config),
            Self::PointwiseAmbient(_) => None,
        }
    }

    pub(super) fn is_pointwise(&self) -> bool {
        matches!(self, Self::PointwiseAmbient(_))
    }

    pub(super) fn sample(&self, cx: &Cx<'_>, time_s: f64, before: bool)
        -> Result<Vec<AmbientRadiationPatch>, ConductionError>
    {
        sample(cx, self.patches(), time_s, before)
    }
}

pub(crate) fn sample(cx: &Cx<'_>, histories: &[RadiationRegionHistory], time_s: f64, before: bool)
    -> Result<Vec<AmbientRadiationPatch>, ConductionError>
{
    if before { return patches_before(cx, histories, time_s); }
    let mut patches = Vec::with_capacity(histories.len());
    for (index, history) in histories.iter().enumerate() {
        poll(cx, index)?;
        patches.push(history.at(time_s)?);
    }
    poll(cx, histories.len())?;
    Ok(patches)
}

// Shared by direct scan attachment and native SDIRK history attachment. This
// uses the pointwise owner's natural-boundary composition, including spatial
// Robin data, nonzero flux and nodes shared with adjacent prescribed faces.
pub(crate) fn admit_pointwise_ambient(cx: &Cx<'_>, run: &EnthalpyHistoryRun<'_, '_, '_>,
    histories: &[RadiationRegionHistory]) -> Result<(), HeterogeneousEnthalpyError>
{
    poll(cx, run.accepted_steps()).map_err(conduction)?;
    if run.accepted_steps() != 0 {
        return Err(invalid("pointwise radiation must be attached before the first physical step"));
    }
    let boundary = run.problem.boundary.at(cx, run.time_s()).map_err(conduction)?;
    let patches = sample(cx, histories, run.time_s(), false).map_err(conduction)?;
    let owner = PointwiseAmbientRadiation::new(cx, run.problem.boundary.mesh(), &boundary, &patches)
        .map_err(conduction)?;
    owner.load_and_report(cx, run.temperature()).map_err(conduction)?;
    poll(cx, histories.len()).map_err(conduction)?;
    Ok(())
}

impl<'h, 'm, 'c> ScanningEnthalpyRun<'h, 'm, 'c> {
    /// Attach local P1 ambient emission with independent reservoir histories.
    /// The original pointwise owner integrates epsilon*sigma*(T_h^4-T_a^4)
    /// over each actual trace. Its fixed sourced emissivity and support remain
    /// authoritative; this does not select evolving material emissivity or an
    /// enclosure. Original convection, prescribed flux and bath reactions keep
    /// their own physical accounts alongside radiation.
    ///
    /// Both temperature-iterated and phase-aware Newton contacts remain active.
    /// All event/adaptive trials sample their own boundary and reservoir side;
    /// rejected work publishes no caloric state, time or accepted-step count.
    /// Select one ambient radiation discretization before stepping; an existing
    /// radiation binding refuses replacement.
    pub fn with_pointwise_ambient_radiation(mut self, cx: &Cx<'_>,
        patches: Vec<RadiationRegionHistory>) -> Result<Self, HeterogeneousEnthalpyError>
    {
        poll(cx, self.run.accepted_steps()).map_err(conduction)?;
        if self.radiation.is_some() {
            return Err(invalid("a scan must select one ambient radiation discretization before stepping"));
        }
        admit_pointwise_ambient(cx, &self.run, &patches)?;
        self.radiation = Some(ScanRadiationBinding::PointwiseAmbient(patches));
        Ok(self)
    }

    /// Original local-emission declarations, absent for area-mean or no radiation.
    #[must_use]
    pub fn pointwise_ambient_patches(&self) -> Option<&[RadiationRegionHistory]> {
        match &self.radiation {
            Some(ScanRadiationBinding::PointwiseAmbient(patches)) => Some(patches),
            _ => None,
        }
    }
}
