//! Retain complete dry-interface sources in the EXISTING material performers.
//! Interface physics lives in the same single/joint force solver. These methods
//! neither reassemble a body nor introduce an independent rendering loop.

use super::MultiBowStringBank;
use crate::bowed_string::runtime::reciprocal::{ReciprocalMaterialString, ReciprocalStringError};
use crate::bowed_string::runtime::reciprocal::bank::ReciprocalStringBank;
use crate::render::schedule::force::coupled::contact::friction::interface::DryFrictionBinding;

impl ReciprocalMaterialString {
    /// Attach the original interface source to the identical retained bow law.
    /// It remains in the physical solver through `into_parts` and scheduling.
    /// Material temperature and interface temperature are separate fixed inputs;
    /// no surface-temperature, rosin chemistry or heat-transfer model is inferred.
    pub fn with_bow_interface(mut self, binding: DryFrictionBinding) -> Result<Self, ReciprocalStringError> {
        self.system = self.system.with_interface(binding).map_err(ReciprocalStringError::Friction)?;
        Ok(self)
    }
    /// Complete original interface data used by subsequent bow steps.
    #[must_use]
    pub fn bow_interface(&self) -> Option<&DryFrictionBinding> { self.system.interface() }
}
impl ReciprocalStringBank {
    /// Bind the existing bow on its actual transverse/torsional surface map.
    /// Unplayed strings, receiver/cavity motion, gas, numerical budgets and the
    /// accepted sample clock are untouched. Conversion to several bows preserves
    /// this source as contact zero; scheduling retains the same physical owner.
    /// Supply `binding.card().friction_law().clone()` at construction rather than
    /// independently copying its coefficients. A mismatch or rebind refuses.
    pub fn with_bow_interface(mut self, binding: DryFrictionBinding) -> Result<Self, ReciprocalStringError> {
        self.system = self.system.with_interface(binding).map_err(ReciprocalStringError::Friction)?;
        Ok(self)
    }
    /// Source checked against the actual current bow contact, not ambient air.
    #[must_use]
    pub fn bow_interface(&self) -> Option<&DryFrictionBinding> { self.system.interface() }
}
impl MultiBowStringBank {
    /// Bind an original interface to one actual contact in drive order.
    /// Each bow keeps its own law, pressure/speed domain, nominal area and fixed
    /// temperature. All contacts validate before ONE joint state publication.
    /// Existing motion and original source bindings remain intact. This does not
    /// automatically establish correspondence to bulk material card identities:
    /// the interface producer owns the ordered surface/texture/history assignment.
    pub fn with_bow_interface(mut self, port: usize, binding: DryFrictionBinding)
        -> Result<Self, ReciprocalStringError>
    {
        self.system = self.system.with_interface(port, binding).map_err(ReciprocalStringError::Friction)?;
        Ok(self)
    }
    /// Original interface for this bow, also available through a scheduled bank.
    #[must_use]
    pub fn bow_interface(&self, port: usize) -> Option<&DryFrictionBinding> { self.system.interface(port) }
}
