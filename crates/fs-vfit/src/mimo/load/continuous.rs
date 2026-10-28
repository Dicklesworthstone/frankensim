//! Continuous access to an already admitted proper rational load.
//!
//! This retains the original pole realization, port ordering and admission
//! reports. Those reports do not supply a storage matrix for arbitrary fits.

use super::{MimoLoadError, PreparedMimoLoad};
use crate::mimo::MimoStateSpace;

/// Transactional transfer of the original real-coordinate pole history.
pub mod state;

/// Borrowed source owner for x' = A x + B u, y = C x + D u.
///
/// Matrices come from the original MimoModel::state_space. A real pole uses
/// one real coordinate per input port. A conjugate pair uses all real parts
/// followed by all NEGATIVE imaginary parts of the original MimoZoh states;
/// its real A block is [[Re(p), Im(p)], [-Im(p), Re(p)]].
///
/// The original sampled report may be GridOnly. This view neither constructs
/// nor claims a physical storage or dissipation certificate for the terminal.
/// Its exchanged work remains a boundary account unless a separate physical
/// storage owner is supplied.
#[derive(Clone, Copy, Debug)]
pub struct MimoContinuous<'a> {
    load: &'a PreparedMimoLoad,
}

impl PreparedMimoLoad {
    /// Immutable source-bound access to the original continuous realization.
    /// Accepted pole history, runtime clock, limits and reports are unchanged.
    #[must_use]
    pub const fn continuous_operator(&self) -> MimoContinuous<'_> {
        MimoContinuous { load: self }
    }
}

impl MimoContinuous<'_> {
    /// Complete original real dimension, including both coordinates of pairs.
    #[must_use]
    pub fn dimension(&self) -> usize {
        self.load.source.order() * self.load.source.ports
    }

    /// Number of inputs/outputs in the original multiport order.
    #[must_use]
    pub fn ports(&self) -> usize {
        self.load.source.ports
    }

    /// Positive finite frequency inside this original source's declared band.
    /// This checks the continuous source band, not a discrete transfer function.
    pub fn admit_frequency(&self, omega_rad_s: f64) -> Result<(), MimoLoadError> {
        let [low, high] = self.load.limits.band_rad_s;
        if !omega_rad_s.is_finite() || omega_rad_s <= 0.0 || omega_rad_s < low || omega_rad_s > high
        {
            return Err(MimoLoadError(
                "continuous MIMO frequency is outside the admitted source band",
            ));
        }
        Ok(())
    }

    /// Materialize the SAME original real A/B/C/D/E under an explicit byte cap.
    ///
    /// D is the instantaneous physical direct term, not the held-step mean
    /// direct matrix. E is zero by the prepared owner's proper-model admission.
    /// No balancing, state reduction, fitting, timestep conversion or D inverse
    /// occurs. A purely static or singular-D model remains representable.
    ///
    /// The retained matrix payload is 8*(n*n + 2*n*ports + 2*ports*ports) bytes,
    /// excluding vector/allocator overhead and caller storage. Original
    /// admission already bounds n <= 256 and ports <= 128; these counts fit
    /// usize on every supported target. Cancellation is checked at entry and
    /// after the single bounded original construction. Refusal has no effects.
    pub fn state_space(
        &self,
        maximum_dense_bytes: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MimoStateSpace, MimoLoadError> {
        if cancelled() {
            return Err(MimoLoadError("continuous MIMO construction cancelled"));
        }
        let n = self.dimension();
        let p = self.ports();
        let required = (n * n + 2 * n * p + 2 * p * p) * core::mem::size_of::<f64>();
        if required > maximum_dense_bytes {
            return Err(MimoLoadError(
                "original continuous MIMO matrices exceed the explicit byte budget",
            ));
        }
        let original = self.load.source.state_space();
        if original
            .a
            .iter()
            .chain(&original.b)
            .chain(&original.c)
            .chain(&original.d)
            .chain(&original.e)
            .any(|v| !v.is_finite())
        {
            return Err(MimoLoadError(
                "original continuous MIMO realization is not finite",
            ));
        }
        if cancelled() {
            return Err(MimoLoadError("continuous MIMO construction cancelled"));
        }
        Ok(original)
    }
}

#[cfg(test)]
mod tests;
