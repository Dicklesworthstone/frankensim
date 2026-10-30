//! Explicit material interfacial shear in the retained instrument bank.
//! No new excitation, body, observer, material defaults or rendering loop.

/// Moving dry-contact finite-cell heat on the original bow geometry.
pub mod dry;

use super::{MultiBowStringBank, ReciprocalStringError};
use crate::acoustic_realize::AcousticRealizeError;
use crate::render::RenderError;
use fs_exec::CancelGate;
use fs_material::state_point::ResolvedShearLayer;
use fs_thermal::enthalpy::{DepositionSimulation, backward_euler::StepCoefficients};
use crate::render::schedule::force::coupled::contact::friction::{FrictionDrive,
    multiple::FrictionPortFrame,
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

    /// Attach each declared film to its original STRING material patch.
    /// The opposing straight fiber ribbon may translate and accelerate through
    /// the existing support-motion owner. String station, tangent direction,
    /// film mass/geometry and elastic/plastic/Maxwell history remain retained.
    /// Normal and tangential reactions share the same moving hair receiver;
    /// film-side migration, wear, opening and mass transport remain unsupported.
    pub fn with_bow_shear_layers_attached_to_strings(
        mut self, bindings: Vec<(usize, ResolvedShearLayer)>, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        self.system = self.system.with_material_shear_layers_attached_to_left(bindings, gate)
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
        self.step_with_contact_thermal_network(thermal, bindings, &[], controls,
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
            return self.with_bow_geometry(external, drives, None, None, None, gate,
                |system, step| system.step_with_contact_thermal_geometry(thermal, bindings,
                    series_heat, controls, energy_guess_j, step.external, step.drives,
                    step.geometry, step.endpoints, step.normal_geometry, step.stop_offsets_m,
                    step.config, gate, |endpoint, _ports| resolve(endpoint)));
        }
        self.system.step_with_contact_thermal_network(thermal, bindings, series_heat,
            controls, energy_guess_j, external, drives, gate, resolve)
    }

    /// Combine finite contact heat with independent finger/fret station motion
    /// and normal support translation on this original string/body/cavity bank.
    ///
    /// Positions are in original stop order; omitted controls hold accepted
    /// geometry. Material maps are projected from the original string basis,
    /// and the generic contact owner accounts actuator work. The source station
    /// cursor advances only after mechanics and heat accept together. Bow/film
    /// footprint migration still requires a transported material-history owner.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_contact_thermal_network_and_stops<'t, S, E: core::fmt::Display>(
        &mut self,
        thermal: &mut DepositionSimulation<'t>,
        bindings: &[LayerThermalCell<'_>],
        series_heat: &[SeriesDriveHeatDestination],
        controls: LayerThermalSolve,
        energy_guess_j: &[f64],
        external: &[f64],
        drives: &[FrictionDrive],
        stations: Option<&[f64]>,
        offsets_m: Option<&[f64]>,
        gate: &CancelGate,
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
    ) -> Result<LayerThermalFrame<S>, RenderError> {
        self.step_with_contact_thermal_feedback_and_stops(thermal, bindings, series_heat,
            controls, energy_guess_j, external, drives, stations, offsets_m, gate,
            |endpoint, _ports| resolve(endpoint))
    }

    /// The original retained-film/stop transaction with trial-dependent heat
    /// transfer. Actual candidate port slip determines which opposing material
    /// cells meet the retained film; source evaluation never advances a cursor.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_contact_thermal_feedback_and_stops<'t, S, E: core::fmt::Display>(
        &mut self,
        thermal: &mut DepositionSimulation<'t>,
        bindings: &[LayerThermalCell<'_>],
        series_heat: &[SeriesDriveHeatDestination],
        controls: LayerThermalSolve,
        energy_guess_j: &[f64],
        external: &[f64],
        drives: &[FrictionDrive],
        stations: Option<&[f64]>,
        offsets_m: Option<&[f64]>,
        gate: &CancelGate,
        resolve: impl Fn(&DepositionSimulation<'t>, &[FrictionPortFrame]) -> Result<StepCoefficients<S>, E>,
    ) -> Result<LayerThermalFrame<S>, RenderError> {
        let endpoints = stations.map(|s| self.prepare_stop_station_geometry(s, gate)).transpose()?;
        let failure = |error: LayerThermalError<E>|
            RenderError::Voice(AcousticRealizeError::Nonlinear(error.to_string()));
        let accepted = if self.bow_station_motion_limits().is_some() {
            self.with_bow_geometry(external, drives, None, endpoints.as_deref(), offsets_m, gate,
                |system, step| system.step_with_contact_thermal_geometry(thermal, bindings,
                    series_heat, controls, energy_guess_j, step.external, step.drives,
                    step.geometry, step.endpoints, step.normal_geometry, step.stop_offsets_m,
                    step.config, gate, resolve)).map_err(failure)?
        } else {
            self.system.step_with_contact_thermal_feedback(
                thermal, bindings, series_heat, controls, energy_guess_j, external,
                drives, endpoints.as_deref(), offsets_m, gate, resolve,
            ).map_err(failure)?
        };
        // The shared owner has published mechanics, normal geometry and heat.
        // No fallible operation may separate that sample from its station cursor.
        if let Some(stations) = stations {
            self.stop_stations.as_mut().expect("admitted stop station state")
                .accepted.copy_from_slice(stations);
        }
        Ok(accepted)
    }

    /// Complete original layer material and geometry, retained through
    /// scheduling and into_parts alongside its constitutive history.
    #[must_use]
    pub fn bow_shear_layer(&self, port: usize) -> Option<&ResolvedShearLayer> {
        self.system.material_shear_layer(port)
    }
}
