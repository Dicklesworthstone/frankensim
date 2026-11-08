//! Uniform acceleration of the original condensed MITC reference inventory.

use super::{MitcShellModel, PlateError};

impl MitcShellModel {
    /// Complete global reference-node forces [N] and axial moments [N m] from
    /// an explicit uniform physical acceleration [m/s²].
    ///
    /// The original physical element mass supplies three cached translation
    /// responses during assembly. This retains the same P1 I0/I1/I2 field,
    /// director interpolation, elastic bubble condensation and chosen mass
    /// image. In particular, asymmetric layers retain signed first-moment
    /// torques. Numerical relative-drilling inertia supplies no body force.
    ///
    /// Every original coordinate is returned, including supported ones.
    /// Constrained mass columns still contribute to free loading. Apply the
    /// original [`Self::project_load`] or native static solver to these complete
    /// loads, and add any separately authored dead loads exactly once.
    ///
    /// This is a fixed-reference dead acceleration. No current-density update,
    /// finite-rotation follower load, rotating-frame force or support motion
    /// is inferred. Independent transient director coordinates use the retained
    /// enriched motion owner's own full-inertia load.
    ///
    /// # Errors
    /// Nonfinite acceleration, force/moment overflow or cancellation. Evaluation
    /// is linear in original coordinate count and returns no partial load.
    pub fn body_acceleration_load(
        &self,
        acceleration_m_s2: [f64; 3],
        cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<f64>, PlateError> {
        self.body_acceleration
            .full_load_with_cancellation(acceleration_m_s2, cancelled)
    }

    /// The same physical body load pulled into the original free support axes.
    ///
    /// Oblique supports use the native work-conjugate rotation after the full
    /// mass response, rather than multiplying a reduced mass by free translation.
    /// This preserves `f_world · u_world = f_free · u_free` for admissible motion.
    ///
    /// # Errors
    /// The complete-load refusals, unrepresentable support projection, or
    /// cancellation during either full loading or support projection.
    pub fn project_body_acceleration(
        &self,
        acceleration_m_s2: [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<f64>, PlateError> {
        let full = self.body_acceleration_load(acceleration_m_s2, &mut cancelled)?;
        self.project_load_cancellable(&full, &mut cancelled)
    }
}
