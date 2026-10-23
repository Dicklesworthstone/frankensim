//! Explicit material interfacial shear in the retained instrument bank.
//! No new excitation, body, observer, material defaults or rendering loop.

use super::{MultiBowStringBank, ReciprocalStringError};
use fs_exec::CancelGate;
use fs_material::state_point::ResolvedShearLayer;
use fs_thermal::enthalpy::{DepositionSimulation, backward_euler::StepCoefficients};
use crate::render::schedule::force::coupled::contact::friction::{FrictionDrive,
    multiple::layer::thermal::{LayerThermalCell, LayerThermalSolve, LayerThermalFrame, LayerThermalError}};

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

    /// Continue the original instrument and its retained-film thermal cells in
    /// one transaction. Temperature-dependent yield/viscosity change the actual
    /// source traction, string/body/cavity motion and the SAME pressure observer.
    ///
    /// The original source bundle must include constant reference density and
    /// modulus. Geometry/strain history remain fixed and the film must stay solid
    /// and compressed. The existing network owns heat capacity and conductive
    /// partition; this adds neither an authored heat split nor a thermal store.
    /// See the shared modal adapter for endpoint accuracy and moving-grid limits.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_layer_thermal_network<'t, S, E>(
        &mut self,
        thermal: &mut DepositionSimulation<'t>,
        bindings: &[LayerThermalCell<'_>],
        controls: LayerThermalSolve,
        energy_guess_j: &[f64],
        external: &[f64],
        drives: &[FrictionDrive],
        gate: &CancelGate,
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
    ) -> Result<LayerThermalFrame<S>, LayerThermalError<E>> {
        self.system.step_with_layer_thermal_network(thermal, bindings, controls,
            energy_guess_j, external, drives, gate, resolve)
    }

    /// Complete original layer material and geometry, retained through
    /// scheduling and into_parts alongside its constitutive history.
    #[must_use]
    pub fn bow_shear_layer(&self, port: usize) -> Option<&ResolvedShearLayer> {
        self.system.material_shear_layer(port)
    }
}
