//! Native moving receivers and ideal images for the SAME stress-potential Hessian.

use super::PreparedMovingQuadrupoleSample;
use crate::AeroacError;
use crate::curle3d::moving::trajectory::{
    PreparedMovingReceiverSource, PreparedReflectedMovingReceiverSource,
    PreparedReflectedSourceTrajectory,
};
use fs_math::dd::Dd;

impl PreparedMovingReceiverSource<'_> {
    /// Complete accelerating quadrupole pressure at the original moving receiver.
    ///
    /// The same causal pair supplies the admitted root, clearance and clocks.
    /// Source acceleration and jerk come from its original quintic. The stress
    /// tensor and both derivatives still belong to the medium's fixed axes and
    /// returned source time; microphone motion adds no pressure multiplier.
    pub fn prepare_quadrupole_receiver(
        &self,
        receiver_time_s: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedMovingQuadrupoleSample, AeroacError> {
        let propagation = self.prepare_receiver(receiver_time_s, cancelled)?.pressure;
        let (source, jerk) = self
            .source()
            .kinematics_with_jerk_at(propagation.source_time_s())?;
        let mut sample = self
            .observer_at(receiver_time_s)?
            .prepare_kinematic_quadrupole(propagation, source, jerk)?;
        // Receiver and source have independent admitted speed envelopes.
        sample.maximum_mach = sample.maximum_mach.max(self.receiver().maximum_mach());
        Ok(sample)
    }
}

impl PreparedReflectedSourceTrajectory<'_> {
    /// One original ideal image of a moving physical volume-stress source.
    ///
    /// Query the ORIGINAL tensor at this path's own retarded time. The inverse
    /// receiver image supplies the full tensor transformation implicitly; no
    /// source tensor is rotated or phase-shifted a second time. Rigid/release
    /// parity is folded into the pressure coefficients exactly once.
    pub fn prepare_quadrupole_path(
        &self,
        path: usize,
        receiver_time_s: f64,
        maximum_time_residual_s: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedMovingQuadrupoleSample, AeroacError> {
        let propagation =
            self.prepare_path(path, receiver_time_s, maximum_time_residual_s, cancelled)?;
        let (source, jerk) = self
            .trajectory()
            .kinematics_with_jerk_at(propagation.source_time_s())?;
        let (observer, parity, velocity_map) = self
            .reflection_observer()
            .stationary_field_path(path, source.position_m)?;
        Ok(observer
            .prepare_kinematic_quadrupole(propagation, source, jerk)?
            .with_image(parity, velocity_map))
    }
}

impl PreparedReflectedMovingReceiverSource<'_> {
    /// One complete quadrupole image at the original moving microphone event.
    ///
    /// The image's original causal solve and current virtual receiver supply
    /// all near, induction, acceleration and jerk terms. Summing `pressure_dd`
    /// over these paths retains interference before any output rounding.
    /// Tensor histories stay in their original medium-frame axes and clocks.
    pub fn prepare_quadrupole_path(
        &self,
        path: usize,
        receiver_time_s: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedMovingQuadrupoleSample, AeroacError> {
        let propagation = self
            .prepare_path(path, receiver_time_s, cancelled)?
            .pressure;
        let (source, jerk) = self
            .source()
            .kinematics_with_jerk_at(propagation.source_time_s())?;
        let (observer, parity, velocity_map) = self.force_history_path_at(path, receiver_time_s)?;
        let mut sample = observer
            .prepare_kinematic_quadrupole(propagation, source, jerk)?
            .with_image(parity, velocity_map);
        sample.maximum_mach = sample.maximum_mach.max(self.receiver().maximum_mach());
        Ok(sample)
    }
}

impl PreparedMovingQuadrupoleSample {
    fn with_image(mut self, parity: f64, velocity_map: [[f64; 3]; 3]) -> Self {
        for coefficient in self.pressure_coefficients.iter_mut().flatten() {
            *coefficient *= parity;
        }
        // The optional observed-clock-rate query takes PHYSICAL receiver
        // velocity. Pull the virtual ray direction back by H^-T, without
        // pressure parity, before its receiver-velocity dot product.
        self.direction = core::array::from_fn(|i| {
            ((0..3).fold(Dd::ZERO, |sum, j| {
                sum + Dd::from_f64(velocity_map[i][j]) * Dd::from_f64(self.direction[j])
            }) * Dd::from_f64(parity))
            .to_f64()
        });
        self
    }
}
