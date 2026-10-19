//! Source-retaining composition with the existing nonlinear mechanics and audio.
//!
//! String bodies remain the original complete prefix. Caller-supplied receivers,
//! contacts, felt and enclosed gas volumes use that SAME mechanical basis and
//! time owner. No receiver response, bridge law or pressure gain is fabricated.

use super::{NonlinearMaterialString, invalid};
use crate::render::plate::impact::{
    ImpactBody, ImpactConfig, ImpactError, ImpactFrame, ImpactSystem,
    PreparedImpactSystem, VolumeSpring,
};
use crate::render::plate::impact::audio::ImpactSource;
use crate::render::plate::impact::felt::FeltPad;
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;

impl NonlinearMaterialString {
    /// Prepare the original nonlinear time owner while retaining this material
    /// source and its physical point maps through subsequent pressure rendering.
    ///
    /// All string coordinates come first: both transverse planes, then any
    /// original torsion. Receiver bodies follow in supplied order. Contact rows,
    /// felt maps and volume areas must address this COMPLETE basis; the existing
    /// owner validates them. A receiver and its reciprocal physical connection
    /// must be supplied explicitly; this operation invents neither one.
    ///
    /// The returned run implements the original `ImpactSource`, so it can be
    /// passed directly to `ImpactPressureRenderer` with solver-derived receiver
    /// velocity projections and radiation in the declared gas.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare(
        self,
        receivers: Vec<ImpactBody>,
        contacts: Vec<Obstacle>,
        pads: Vec<FeltPad>,
        volumes: Vec<VolumeSpring>,
        config: ImpactConfig,
        gate: &CancelGate,
    ) -> Result<MaterialStringRun, ImpactError> {
        if gate.is_requested() {
            return Err(ImpactError::Cancelled);
        }
        if config.dt_s.to_bits() != (1.0 / f64::from(self.sample_rate_hz())).to_bits() {
            return Err(invalid("nonlinear material runtime must use its actual source clock"));
        }
        let mut bodies = self.bodies().to_vec();
        bodies.extend(receivers);
        let system = ImpactSystem::new(bodies, contacts, pads, volumes, config)?
            .prepare_analytic()?;
        if gate.is_requested() {
            return Err(ImpactError::Cancelled);
        }
        Ok(MaterialStringRun { source: self, system })
    }
}

/// One original material source, its connected receivers and accepted mechanics.
/// Storage, dissipation, external work, trials and rollback belong exclusively
/// to the existing prepared ImpactSystem. This wrapper adds no stepping state.
pub struct MaterialStringRun {
    source: NonlinearMaterialString,
    system: PreparedImpactSystem,
}
impl MaterialStringRun {
    /// Original section, atmosphere, physical directions and point/torque maps.
    pub const fn source(&self) -> &NonlinearMaterialString {
        &self.source
    }
    /// Read the accepted full state, string tension, energy and solver diagnostics.
    pub const fn system(&self) -> &PreparedImpactSystem {
        &self.system
    }
    /// Advance the existing shared mechanical owner once, transactionally.
    pub fn step(
        &mut self, forces: &[f64], gate: &CancelGate,
    ) -> Result<ImpactFrame, ImpactError> {
        self.system.step(forces, gate)
    }
}
impl ImpactSource for MaterialStringRun {
    type Frame = ImpactFrame;
    fn state(&self) -> &[f64] {
        ImpactSource::state(&self.system)
    }
    fn mode_count(&self) -> usize {
        ImpactSource::mode_count(&self.system)
    }
    fn samples_rendered(&self) -> u64 {
        ImpactSource::samples_rendered(&self.system)
    }
    fn sample_period_s(&self) -> f64 {
        ImpactSource::sample_period_s(&self.system)
    }
    fn remaining_steps(&self) -> u64 {
        ImpactSource::remaining_steps(&self.system)
    }
    fn maximum_generalized_force(&self) -> f64 {
        ImpactSource::maximum_generalized_force(&self.system)
    }
    fn advance(
        &mut self, forces: &[f64], gate: &CancelGate,
    ) -> Result<Self::Frame, ImpactError> {
        self.step(forces, gate)
    }
}
