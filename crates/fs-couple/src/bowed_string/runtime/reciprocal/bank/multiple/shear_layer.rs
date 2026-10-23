//! Explicit material interfacial shear in the retained instrument bank.
//! No new excitation, body, observer, material defaults or rendering loop.

use super::{MultiBowStringBank, ReciprocalStringError};
use fs_exec::CancelGate;
use fs_material::state_point::ResolvedShearLayer;

impl MultiBowStringBank {
    /// Bind isothermal thin shear layers to selected original bow/contact ports.
    ///
    /// Each source retains its actual G, shear yield stress, plastic viscosity,
    /// nominal contact area, thickness and small elastic-shear limit. The
    /// material's internal displacement and physical loss enter the SAME joint
    /// string/body/cavity mechanics and therefore the existing pressure output.
    /// Geometry/material swaps propagate through the existing coupled response;
    /// no spectrogram target, fitted friction curve or audio shaping is used.
    ///
    /// This is a fixed-layer, retained-compression rung. Bind once before
    /// playback; dry-interface cards cannot also control those ports. Opening,
    /// migrating layer history and nonlinear geometric-storage composition
    /// refuse explicitly. Temperature is the declared frozen source state;
    /// conduction, flash heating, wear and measured rosin behavior are not
    /// supplied by this adapter.
    pub fn with_bow_shear_layers(
        mut self,
        bindings: Vec<(usize, ResolvedShearLayer)>,
        gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        self.system = self
            .system
            .with_material_shear_layers(bindings, gate)
            .map_err(ReciprocalStringError::Friction)?;
        Ok(self)
    }

    /// Complete original layer material and geometry, retained through
    /// scheduling and into_parts alongside its constitutive history.
    #[must_use]
    pub fn bow_shear_layer(&self, port: usize) -> Option<&ResolvedShearLayer> {
        self.system.material_shear_layer(port)
    }
}
