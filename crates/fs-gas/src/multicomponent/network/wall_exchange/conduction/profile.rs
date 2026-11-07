//! Source-integrated solid profiles on the original gas/wall transactions.

use super::*;
use fs_material::thermal::conduction::{ConductivityCurve, LayerConductionConfig};
use fs_thermal::enthalpy::DepositionError;

#[derive(Clone, Debug, PartialEq)]
pub(super) struct IntegratedWallProfile {
    curve: ConductivityCurve,
    geometry: ReservoirLayerGeometry,
    maximum_solid_conductance_w_k: f64,
}

impl MaterialWallConduction {
    /// Integrate the original piecewise-affine k(T) through the source-free solid
    /// and solve its actual face temperature against the original gas film.
    /// Both retained wall and gas temperatures are the proposed implicit endpoint.
    /// Preparation is shared across clones; no curve is rebuilt during evolution.
    ///
    /// This opt-in requires `ImplicitEndpoint` evaluation. Finite local and
    /// shared walls retain their original material caloric binding; a prescribed
    /// local reservoir supplies external heat without a stored wall inventory.
    /// Branch, plenum and shared gas ports use the same solid/film endpoint law.
    /// Wall-to-wall links require a two-layer law and still refuse this opt-in.
    /// A phase-dependent k(T,phi) cannot be integrated from a lumped wall's one h:
    /// its missing spatial phase profile is not inferred. No storage is added,
    /// and quasi-steady conduction does not resolve transient internal diffusion.
    pub fn with_integrated_profile(mut self) -> Result<Self, WallExchangeError> {
        let ConductionSource::Temperature(source) = &self.source else {
            return Err(WallExchangeError::Invalid(
                "integrated conductivity requires a temperature-only source; spatial phase state is absent",
            ));
        };
        // The canonical source-line owner guarantees continuous piecewise-affine
        // interpolation. Retain every conductivity breakpoint, without adding
        // redundant density/Cp knots or changing the original source selection.
        let line = &source.source_lines()[2];
        let mut knots = Vec::new();
        knots
            .try_reserve_exact(line.knots.len())
            .map_err(|_| WallExchangeError::Invalid("integrated wall profile allocation"))?;
        let mut maximum_k = 0.0_f64;
        for (temperature, answer) in &line.knots {
            let conductivity = answer.evidence.value.value;
            maximum_k = maximum_k.max(conductivity);
            knots.push([*temperature, conductivity]);
        }
        let curve = ConductivityCurve::new(knots)?;
        let geometry = match self.geometry {
            Some(geometry) => geometry,
            None => ReservoirLayerGeometry::planar(self.area.value(), self.length.value())?,
        };
        let maximum_solid_conductance_w_k = maximum_k * geometry.conductance_multiplier_m();
        if !positive(maximum_solid_conductance_w_k) {
            return Err(WallExchangeError::Invalid(
                "integrated wall conductance range",
            ));
        }
        self.integrated_profile = Some(Arc::new(IntegratedWallProfile {
            curve,
            geometry,
            maximum_solid_conductance_w_k,
        }));
        self.evaluation = WallConductionEvaluation::ImplicitEndpoint;
        Ok(self)
    }

    /// Whether the original k(T) is integrated through the solid profile.
    pub fn has_integrated_profile(&self) -> bool {
        self.integrated_profile.is_some()
    }

    pub(crate) fn profile_gas_sensitivity_bound(&self, film: f64) -> Option<f64> {
        self.integrated_profile
            .as_ref()
            .map(|profile| profile.gas_conductance_bound(film))
    }

    pub(crate) fn profile_response(
        &self,
        temperatures_k: [f64; 2],
        film: f64,
        duration: fs_qty::Time,
        control: WallControl,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<WallHeatResponse, WallExchangeError> {
        self.integrated_profile
            .as_ref()
            .ok_or(WallExchangeError::Invalid(
                "integrated wall profile missing",
            ))?
            .evaluate(temperatures_k, film, duration, control, cancelled)
    }
}

impl IntegratedWallProfile {
    /// Global bound on |dP/dTgas|, including a film-dominated or insulating wall.
    /// The original gas caloric inverse uses this to fund its temperature error.
    pub(super) fn gas_conductance_bound(&self, film: f64) -> f64 {
        let solid = self.maximum_solid_conductance_w_k;
        let smaller = film.min(solid);
        smaller / (1.0 + smaller / film.max(solid))
    }

    pub(super) fn evaluate(
        &self,
        temperatures_k: [f64; 2],
        film: f64,
        duration: fs_qty::Time,
        control: WallControl,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<WallHeatResponse, WallExchangeError> {
        if cancelled() {
            return Err(WallExchangeError::Cancelled);
        }
        // A disabled gas film is genuinely insulating; never invert zero or
        // require the disconnected gas temperature to lie in the solid chart.
        if film == 0.0 {
            self.curve.conductivity_w_m_k(temperatures_k[0])?;
            return Ok(WallHeatResponse {
                conductance_w_per_k: 0.0,
                power_to_gas_w: 0.0,
                temperature_derivative_w_per_k: Some([0.0; 2]),
                flux_error_w: 0.0,
            });
        }
        let dt = duration.value();
        let absolute_heat_tolerance_w = if dt > 0.0 {
            (0.0625 * control.absolute_energy_tolerance_j / dt).min(f64::MAX / 16.0)
        } else {
            0.0625 * control.absolute_energy_tolerance_j
        };
        let response = self.geometry.evaluate_material_to_bath(
            &self.curve,
            temperatures_k,
            1.0 / film,
            LayerConductionConfig {
                max_iterations: control.maximum_iterations,
                absolute_heat_tolerance_w,
                relative_heat_tolerance: 0.0625 * control.relative_energy_tolerance,
            },
            || {
                if cancelled() {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            },
        )?;
        Ok(WallHeatResponse {
            conductance_w_per_k: response.conductance_w_k,
            power_to_gas_w: response.power_to_reservoir_w,
            temperature_derivative_w_per_k: Some(response.endpoint_derivative_w_k),
            flux_error_w: response.heat_residual_w.abs(),
        })
    }
}
