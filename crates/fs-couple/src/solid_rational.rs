//! Original admitted rational radiation/impedance attached to hereditary solids.
//!
//! The original continuous A/B/C/D/E and port order enter the solid owner's
//! complete polynomial. Accepted ZOH history and the physical source clock are
//! borrowed unchanged; steady harmonic/free-mode analysis neither advances
//! history nor turns it into additional harmonic forcing. Original passivity
//! reports retain their classification and do not supply an energy witness.

use fs_exec::Cx;
use fs_solid::linear3::visco::{
    PreparedViscoPoles, ViscoElements, ViscoJointHarmonicBudget, ViscoJointHarmonicError,
    ViscoJointHarmonicFactor, ViscoPoleBudget, ViscoPoleDynamicBudget, ViscoPoleError,
    ViscoPoleStateSpaceLoad,
};
use fs_vfit::mimo::{
    MimoPassivityReport,
    load::{MimoLoadError, PreparedMimoLoad},
};

/// Separate allowances for the original real realization and the complete
/// physical solid/boundary image. Each owner's numerical payload is additive.
#[derive(Clone, Copy, Debug)]
pub struct SolidMimoResonanceBudget {
    pub maximum_realization_bytes: usize,
    pub poles: ViscoPoleBudget,
    pub load: ViscoPoleDynamicBudget,
}

#[derive(Clone, Debug, PartialEq)]
pub enum SolidMimoResonanceError {
    Invalid(&'static str),
    Source(MimoLoadError),
    Solid(ViscoPoleError),
    Harmonic(ViscoJointHarmonicError),
    Cancelled,
}

impl core::fmt::Display for SolidMimoResonanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-SOLID-MIMO-RESONANCE: {self:?}")
    }
}
impl std::error::Error for SolidMimoResonanceError {}
impl From<MimoLoadError> for SolidMimoResonanceError {
    fn from(error: MimoLoadError) -> Self {
        Self::Source(error)
    }
}
impl From<ViscoPoleError> for SolidMimoResonanceError {
    fn from(error: ViscoPoleError) -> Self {
        Self::Solid(error)
    }
}
impl From<ViscoJointHarmonicError> for SolidMimoResonanceError {
    fn from(error: ViscoJointHarmonicError) -> Self {
        Self::Harmonic(error)
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), SolidMimoResonanceError> {
    cx.checkpoint().map_err(|_| SolidMimoResonanceError::Cancelled)
}

/// Complete original solid/material/rational-boundary resonance image.
/// The original admitted source remains available with its reports and history.
pub struct PreparedSolidMimoResonance<'body, 'load> {
    source: &'load PreparedMimoLoad,
    model: PreparedViscoPoles<'body>,
}

impl<'body, 'load> PreparedSolidMimoResonance<'body, 'load> {
    /// Couple the source's original flow-to-effort impedance through one
    /// row-major free_displacements*source_ports force map. The load receives
    /// map^T*u' and the body receives the opposing map*effort.
    ///
    /// All real/conjugate source coordinates and full mutual/direct/derivative
    /// terms survive. The PreparedMimoLoad proper-model admission currently
    /// makes E exactly zero; this adapter still passes that original term.
    /// No digital mean-direct matrix, fitted damping, pole reduction, inverse
    /// bare impedance or invented boundary storage enters the model.
    ///
    /// Dimensions and the complete map are checked before materializing the
    /// bounded original realization. Source and solid budgets remain distinct.
    /// The resulting generic boundary reports signed work; its complete total
    /// storage and dissipation accounts remain None without a physical witness.
    pub fn new(
        body: &'body ViscoElements,
        source: &'load PreparedMimoLoad,
        force_map: &[f64],
        budget: SolidMimoResonanceBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, SolidMimoResonanceError> {
        checkpoint(cx)?;
        let view = source.continuous_operator();
        let (states, ports) = (view.dimension(), view.ports());
        let mechanical = body.equilibrium().free_dofs.len();
        let coordinates = mechanical
            .checked_add(body.internal_variables())
            .and_then(|n| n.checked_add(states))
            .ok_or(SolidMimoResonanceError::Invalid("complete solid/source extent overflow"))?;
        let entries = mechanical
            .checked_mul(ports)
            .ok_or(SolidMimoResonanceError::Invalid("solid/source port map extent overflow"))?;
        if mechanical == 0 || force_map.len() != entries
            || coordinates > budget.poles.maximum_coordinates
            || states > budget.load.maximum_states || ports > budget.load.maximum_ports
        {
            return Err(SolidMimoResonanceError::Invalid(
                "complete solid/source states, original map and explicit dimensions required",
            ));
        }
        for chunk in force_map.chunks(64) {
            checkpoint(cx)?;
            if chunk.iter().any(|v| !v.is_finite()) {
                return Err(SolidMimoResonanceError::Invalid("finite original force map required"));
            }
        }
        let matrices = view.state_space(
            budget.maximum_realization_bytes,
            &mut || cx.checkpoint().is_err(),
        )?;
        checkpoint(cx)?;
        let band = source.limits().band_rad_s;
        let model = body.prepare_resonance_poles_with_state_space_load(
            ViscoPoleStateSpaceLoad {
                states: matrices.n,
                ports: matrices.ports,
                a: &matrices.a,
                b: &matrices.b,
                c: &matrices.c,
                d: &matrices.d,
                e: &matrices.e,
                force_map,
                band_rad_s: (band[0], band[1]),
            },
            budget.poles,
            budget.load,
            cx,
        )?;
        checkpoint(cx)?;
        Ok(Self { source, model })
    }

    pub fn source(&self) -> &'load PreparedMimoLoad {
        self.source
    }

    /// Original continuous and sampled reports; no stronger certificate added.
    pub fn reports(&self) -> (&MimoPassivityReport, &MimoPassivityReport) {
        self.source.reports()
    }

    /// Shared owner for full original-state poles, contours and observations.
    pub fn model(&self) -> &PreparedViscoPoles<'body> {
        &self.model
    }

    /// Reusable full-state continuous factor. Source and material bands both
    /// apply; the source's accepted history and sample clock remain unchanged.
    pub fn harmonic_factor(
        &self,
        omega_rad_s: f64,
        budget: ViscoJointHarmonicBudget,
        cx: &Cx<'_>,
    ) -> Result<ViscoJointHarmonicFactor<'_, 'body>, SolidMimoResonanceError> {
        checkpoint(cx)?;
        self.source.continuous_operator().admit_frequency(omega_rad_s)?;
        Ok(self.model.harmonic_factor(omega_rad_s, budget, cx)?)
    }
}
