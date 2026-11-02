//! Exact piecewise-constant generation in the original planar wall owner.
//!
//! In a layer, `dq/dx = s` and `dT/dx = -q/k`, hence
//! `q(x) = q_left + s*x` and `T(x) = T_left - x*(q_left + s*x/2)/k`.
//! Contacts and films drop temperature by their resistance times the LOCAL
//! flux. Summing those drops closes both reservoir conditions in one scalar
//! equation; imposing a common flux would violate conservation.
//! See MIT 3.044, spring 2013, problem set 2 (internally heated slab):
//! <https://ocw.mit.edu/courses/3-044-materials-processing-spring-2013/ff9d5adcbafc16451f8e329865aba5e7_MIT3_044S13_pset2solns.pdf>.

use super::{LayeredWall, checkpoint, finite, invalid, positive};
use crate::ConductionError;
use fs_exec::Cx;

/// Exact quadratic temperature and linear flux within one homogeneous layer.
#[derive(Debug, Clone, PartialEq)]
pub struct GeneratedLayerTemperatures {
    /// Absolute temperatures on the left and right faces, K.
    pub face_temperatures_k: [f64; 2],
    /// Signed left-to-right heat flux on the two faces, W/m2.
    pub face_fluxes_w_per_m2: [f64; 2],
    /// Spatial mean of the quadratic temperature profile, K.
    pub mean_temperature_k: f64,
    /// Minimum absolute temperature, K.
    pub minimum_temperature_k: f64,
    /// Minimum's position measured from the left face, m.
    pub minimum_position_m: f64,
    /// Maximum absolute temperature, K.
    pub maximum_temperature_k: f64,
    /// Maximum's position measured from the left face, m.
    pub maximum_position_m: f64,
    /// Fourier flux at the layer midpoint, recovered from its rounded faces.
    pub recomputed_midpoint_flux_w_per_m2: f64,
    /// Recovered midpoint flux minus the solved midpoint flux, W/m2.
    pub midpoint_flux_residual_w_per_m2: f64,
}

/// Steady conduction with declared per-layer internal heat production/removal.
#[derive(Debug, Clone, PartialEq)]
pub struct GeneratedLayeredWallSolution {
    /// Unchanged geometry, contacts, films, source tables and receipts.
    pub wall: LayeredWall,
    /// Original left, right reservoir temperatures, K.
    pub reservoir_temperatures_k: [f64; 2],
    /// Original heat-transfer area, m2.
    pub area_m2: f64,
    /// Original uniform source densities in physical layer order, W/m3.
    /// Negative entries declare heat removal, not an inferred mechanism.
    pub volumetric_generation_w_m3: Vec<f64>,
    /// Signed left-to-right flux at the left and right stack ends, W/m2.
    pub boundary_fluxes_w_per_m2: [f64; 2],
    /// Heat delivered OUTWARD to the left and right reservoirs, W.
    /// Either value may be negative when that reservoir supplies heat.
    pub reservoir_heat_rates_w: [f64; 2],
    /// Volume-integrated signed internal generation, W.
    pub generated_heat_rate_w: f64,
    /// Outward reservoir heat minus internal generation, W.
    pub energy_residual_w: f64,
    /// Original layer order, including actual interior extrema.
    pub layer_temperatures: Vec<GeneratedLayerTemperatures>,
    /// Largest absolute recovered layer/film/contact flux defect, W/m2.
    /// A roundoff diagnostic, not an uncertainty certificate.
    pub maximum_flux_residual_w_per_m2: f64,
    /// Reconstructed right reservoir minus its prescribed temperature, K.
    pub right_boundary_residual_k: f64,
}

impl LayeredWall {
    /// Solve `-(k T')' = s` with one uniform source density per layer.
    ///
    /// Conductivity must be temperature independent; each layer may have its
    /// own value and source span. Source densities are explicit finite W/m3
    /// values in original layer order. Positive means production; negative
    /// means removal. No Joule, reaction or refrigeration law is inferred.
    ///
    /// The solution is analytic for this one-dimensional model, with work linear
    /// in the layers and retained knots and O(layers) result storage. No mesh or
    /// nonlinear iteration is needed. Contact jumps use local
    /// flux. Both faces and interior extrema must lie in the material span.
    /// Diagnostics describe returned rounded values. Cancellation or refusal
    /// publishes no partial solution. Nonlinear k(T), temperature-dependent
    /// generation, lateral spreading and transient storage are outside this
    /// constant-coefficient analytical path.
    pub fn solve_with_generation(
        &self,
        cx: &Cx<'_>,
        reservoir_temperatures_k: [f64; 2],
        area_m2: f64,
        volumetric_generation_w_m3: &[f64],
    ) -> Result<GeneratedLayeredWallSolution, ConductionError> {
        checkpoint(cx, 0)?;
        positive("wall area", area_m2)?;
        for &temperature in &reservoir_temperatures_k {
            positive("reservoir absolute temperature", temperature)?;
        }
        if volumetric_generation_w_m3.len() != self.layers.len() {
            return Err(invalid("supply exactly one volumetric source per layer"));
        }
        let mut resistance = Vec::with_capacity(self.layers.len());
        let mut generation = Vec::with_capacity(self.layers.len());
        let mut total_resistance = self.surface_resistances_m2_k_per_w[0];
        let mut cumulative_generation = 0.0;
        let mut generated_drop = 0.0;
        for (index, (layer, &source)) in self
            .layers
            .iter()
            .zip(volumetric_generation_w_m3)
            .enumerate()
        {
            checkpoint(cx, index)?;
            if layer.conductivity.is_temperature_dependent() {
                return Err(ConductionError::Conductivity {
                    what: format!(
                        "layer {:?} varies with temperature; internally heated planar conduction requires constant conductivity",
                        layer.name
                    ),
                });
            }
            finite("volumetric generation", source)?;
            let r = layer.thickness_m / layer.conductivity.knots()[0].1;
            positive("layer area-specific resistance", r)?;
            let g = finite("layer generation per area", source * layer.thickness_m)?;
            if source != 0.0 && g == 0.0 {
                return Err(invalid("nonzero volumetric generation underflowed"));
            }
            generated_drop = finite(
                "generation temperature drop",
                r.mul_add(cumulative_generation + 0.5 * g, generated_drop),
            )?;
            cumulative_generation = finite("cumulative generation", cumulative_generation + g)?;
            total_resistance = finite("total resistance", total_resistance + r)?;
            if let Some(&contact) = self.contact_resistances_m2_k_per_w.get(index) {
                generated_drop = finite(
                    "generation contact drop",
                    contact.mul_add(cumulative_generation, generated_drop),
                )?;
                total_resistance = finite("total resistance", total_resistance + contact)?;
            }
            resistance.push(r);
            generation.push(g);
        }
        let right_film = self.surface_resistances_m2_k_per_w[1];
        total_resistance = finite("total resistance", total_resistance + right_film)?;
        positive("total resistance", total_resistance)?;
        generated_drop = finite(
            "generation film drop",
            right_film.mul_add(cumulative_generation, generated_drop),
        )?;
        let driving_drop = finite(
            "driving temperature drop",
            reservoir_temperatures_k[0] - reservoir_temperatures_k[1] - generated_drop,
        )?;
        let left_flux = finite("left boundary heat flux", driving_drop / total_resistance)?;
        if driving_drop != 0.0 && left_flux == 0.0 {
            return Err(invalid("nonzero boundary heat flux underflowed"));
        }
        let right_flux = finite(
            "right boundary heat flux",
            left_flux + cumulative_generation,
        )?;
        let reservoir_heat_rates_w = [
            finite("left reservoir heat rate", -left_flux * area_m2)?,
            finite("right reservoir heat rate", right_flux * area_m2)?,
        ];
        for (&flux, &heat) in [left_flux, right_flux].iter().zip(&reservoir_heat_rates_w) {
            if flux != 0.0 && heat == 0.0 {
                return Err(invalid("nonzero reservoir heat rate underflowed"));
            }
        }
        let generated_heat_rate_w =
            finite("total generated heat", cumulative_generation * area_m2)?;
        if cumulative_generation != 0.0 && generated_heat_rate_w == 0.0 {
            return Err(invalid("nonzero total generated heat underflowed"));
        }
        let energy_residual_w = finite(
            "energy residual",
            reservoir_heat_rates_w[0] + reservoir_heat_rates_w[1] - generated_heat_rate_w,
        )?;
        let left_film = self.surface_resistances_m2_k_per_w[0];
        let mut temperature = finite(
            "left face temperature",
            (-left_flux).mul_add(left_film, reservoir_temperatures_k[0]),
        )?;
        let mut maximum_flux_residual = 0.0_f64;
        if left_film > 0.0 {
            maximum_flux_residual = flux_defect(
                (reservoir_temperatures_k[0] - temperature) / left_film,
                left_flux,
            )?;
        }
        let mut produced = 0.0;
        let mut layers = Vec::with_capacity(self.layers.len());
        for (index, layer) in self.layers.iter().enumerate() {
            checkpoint(cx, index)?;
            let r = resistance[index];
            let g = generation[index];
            let flux_in = finite("layer inlet flux", left_flux + produced)?;
            produced = finite("accumulated source", produced + g)?;
            let flux_out = finite("layer outlet flux", left_flux + produced)?;
            let midpoint_flux = finite("layer midpoint flux", flux_in + 0.5 * g)?;
            let left = temperature;
            let right = finite("layer right temperature", (-r).mul_add(midpoint_flux, left))?;
            let mean = finite(
                "layer mean temperature",
                (-r).mul_add(0.5 * flux_in + g / 6.0, left),
            )?;
            let (mut minimum, mut minimum_x, mut maximum, mut maximum_x) = if left <= right {
                (left, 0.0, right, layer.thickness_m)
            } else {
                (right, layer.thickness_m, left, 0.0)
            };
            // At most one stationary point; a dimensionless position avoids
            // squaring a potentially large coordinate.
            if (flux_in < 0.0 && flux_out > 0.0) || (flux_in > 0.0 && flux_out < 0.0) {
                let fraction = finite("stationary layer coordinate", -flux_in / g)?;
                if fraction <= 0.0 || fraction >= 1.0 {
                    return Err(invalid(
                        "interior temperature extremum is not representable",
                    ));
                }
                let value = finite(
                    "interior temperature extremum",
                    (-0.5 * r * fraction).mul_add(flux_in, left),
                )?;
                let x = finite("interior extremum position", fraction * layer.thickness_m)?;
                if value < minimum {
                    minimum = value;
                    minimum_x = x;
                }
                if value > maximum {
                    maximum = value;
                    maximum_x = x;
                }
            }
            for value in [minimum, maximum] {
                positive("layer absolute temperature", value)?;
                layer.conductivity.span().check(value)?;
            }
            let recomputed = finite("recovered layer midpoint flux", (left - right) / r)?;
            let residual = finite("midpoint Fourier residual", recomputed - midpoint_flux)?;
            maximum_flux_residual = maximum_flux_residual.max(residual.abs());
            layers.push(GeneratedLayerTemperatures {
                face_temperatures_k: [left, right],
                face_fluxes_w_per_m2: [flux_in, flux_out],
                mean_temperature_k: mean,
                minimum_temperature_k: minimum,
                minimum_position_m: minimum_x,
                maximum_temperature_k: maximum,
                maximum_position_m: maximum_x,
                recomputed_midpoint_flux_w_per_m2: recomputed,
                midpoint_flux_residual_w_per_m2: residual,
            });
            temperature = right;
            if let Some(&contact) = self.contact_resistances_m2_k_per_w.get(index) {
                temperature = finite(
                    "contact downstream temperature",
                    (-contact).mul_add(flux_out, right),
                )?;
                if contact > 0.0 {
                    maximum_flux_residual = maximum_flux_residual
                        .max(flux_defect((right - temperature) / contact, flux_out)?);
                }
            }
        }
        let reconstructed_right = finite(
            "right reservoir reconstruction",
            (-right_film).mul_add(right_flux, temperature),
        )?;
        let right_boundary_residual_k = finite(
            "right boundary residual",
            reconstructed_right - reservoir_temperatures_k[1],
        )?;
        if right_film > 0.0 {
            maximum_flux_residual = maximum_flux_residual.max(flux_defect(
                (temperature - reservoir_temperatures_k[1]) / right_film,
                right_flux,
            )?);
        }
        checkpoint(cx, self.layers.len())?;
        Ok(GeneratedLayeredWallSolution {
            wall: self.clone(),
            reservoir_temperatures_k,
            area_m2,
            volumetric_generation_w_m3: volumetric_generation_w_m3.to_vec(),
            boundary_fluxes_w_per_m2: [left_flux, right_flux],
            reservoir_heat_rates_w,
            generated_heat_rate_w,
            energy_residual_w,
            layer_temperatures: layers,
            maximum_flux_residual_w_per_m2: maximum_flux_residual,
            right_boundary_residual_k,
        })
    }
}

fn flux_defect(recovered: f64, expected: f64) -> Result<f64, ConductionError> {
    finite("film/contact Fourier residual", recovered - expected).map(f64::abs)
}
