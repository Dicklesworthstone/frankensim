//! Explicit material interfacial shear in the retained instrument bank.
//! No new excitation, body, observer, material defaults or rendering loop.

use super::{MultiBowStringBank, ReciprocalStringError};
use fs_exec::CancelGate;
use fs_material::state_point::ResolvedShearLayer;
use fs_thermal::enthalpy::{DepositionSimulation, backward_euler::StepCoefficients};
use crate::render::schedule::force::coupled::contact::friction::{FrictionDrive,
    multiple::layer::thermal::{LayerThermalCell, LayerThermalSolve, LayerThermalFrame, LayerThermalError,
        SeriesDriveHeatDestination}};

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
    /// and migrating layer history refuse explicitly. Nonlinear string strain
    /// and any compliant-drive springs participate in the same staged solve
    /// with the layer's own reversible storage. Temperature is the declared frozen source state;
    /// conduction, flash heating, wear and measured rosin behavior are not
    /// supplied by this adapter.
    pub fn with_bow_shear_layers(
        mut self,
        bindings: Vec<(usize, ResolvedShearLayer)>,
        gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        if !bindings.is_empty() && self.bow_station_motion_limits().is_some() {
            return Err(super::invalid("moving bow stations cannot bind fixed material shear history"));
        }
        self.system = self
            .system
            .with_material_shear_layers(bindings, gate)
            .map_err(ReciprocalStringError::Friction)?;
        Ok(self)
    }

    /// Continue the original instrument and its retained-film thermal cells in
    /// one transaction. Temperature-dependent elasticity/yield/viscosity change the actual
    /// source traction, string/body/cavity motion and the SAME pressure observer.
    ///
    /// The original source bundle must include constant reference density and
    /// constant or affine modulus. Reversible heat and mechanical entropy keep
    /// their own account beside plastic loss. Geometry and the tangent frame
    /// remain fixed, strain history is retained, and the film must stay solid
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
        if self.bow_station_motion_limits().is_some() {
            return Err(LayerThermalError::Invalid("contact heat requires transported bow geometry and material history"));
        }
        self.system.step_with_layer_thermal_network(thermal, bindings, controls,
            energy_guess_j, external, drives, gate, resolve)
    }

    /// Continue original instrument motion, layer thermodynamics and explicitly
    /// routed series-hair heating on the SAME finite thermal network.
    /// The heat is the actual drive dashpot loss; stored hair spring energy is
    /// retained mechanically. Existing conduction can carry that heat into a
    /// material film, changing its traction and this instrument's pressure.
    /// Hair drive coefficients remain their declared fixed values: no material
    /// temperature law or thermal mass is inferred from spring/dashpot data.
    /// Empty layer bindings permit a series-drive-only thermal continuation.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_contact_thermal_network<'t, S, E>(
        &mut self,
        thermal: &mut DepositionSimulation<'t>,
        bindings: &[LayerThermalCell<'_>],
        series_heat: &[SeriesDriveHeatDestination],
        controls: LayerThermalSolve,
        energy_guess_j: &[f64],
        external: &[f64],
        drives: &[FrictionDrive],
        gate: &CancelGate,
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
    ) -> Result<LayerThermalFrame<S>, LayerThermalError<E>> {
        if self.bow_station_motion_limits().is_some() {
            return Err(LayerThermalError::Invalid("contact heat requires transported bow geometry and material history"));
        }
        self.system.step_with_contact_thermal_network(thermal, bindings, series_heat,
            controls, energy_guess_j, external, drives, gate, resolve)
    }

    /// Complete original layer material and geometry, retained through
    /// scheduling and into_parts alongside its constitutive history.
    #[must_use]
    pub fn bow_shear_layer(&self, port: usize) -> Option<&ResolvedShearLayer> {
        self.system.material_shear_layer(port)
    }
}
