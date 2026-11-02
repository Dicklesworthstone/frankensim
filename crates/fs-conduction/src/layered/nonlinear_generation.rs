//! Internally heated stacks with the original positive piecewise-linear k(T).
//!
//! For a planar wall the Kirchhoff potential K(T)=integral(k(T), dT) obeys
//! K''=-s, so its profile is quadratic even when temperature is not. The
//! shared radial caller supplies exact reciprocal-area and source moments.
//! Exact source-segment inversion recovers temperature; a monotone scalar
//! search closes the right reservoir. See Peavy, DOI:10.1177/109719639602000107
//! for the steady one-dimensional Kirchhoff transformation. In all geometries
//! dQ/dr=s*A(r) and Q=-A(r)*dK/dr; geometry owns only those spatial integrals.

use super::geometry::LayerShape;
use super::nonlinear::ConductivityIntegral;
use super::{LayeredWall, NonlinearWallOptions, checkpoint, finite, invalid, positive};
use crate::ConductionError;
use fs_exec::Cx;

/// Explicit material-temperature admission and bounded nonlinear work.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NonlinearGenerationOptions {
    /// Closed positive interval containing every material temperature, K.
    /// Every layer must support the entire interval. This conservative common
    /// interval is caller-declared; it is never inferred by extrapolating a
    /// table. Reservoirs outside it are allowed behind nonzero surface films.
    pub temperature_interval_k: [f64; 2],
    /// Shooting budget and temperature/constitutive-flux residual tolerances.
    pub solver: NonlinearWallOptions,
}

/// A layer's nonlinear temperature profile and exact linear flux profile.
#[derive(Debug, Clone, PartialEq)]
pub struct NonlinearGeneratedLayerTemperatures {
    /// Left and right absolute face temperatures, K.
    pub face_temperatures_k: [f64; 2],
    /// Signed left-to-right flux on both faces, W/m2.
    pub face_fluxes_w_per_m2: [f64; 2],
    /// Minimum temperature over the entire layer, K.
    pub minimum_temperature_k: f64,
    /// Minimum's position from this layer's left face, m.
    pub minimum_position_m: f64,
    /// Maximum temperature over the entire layer, K.
    pub maximum_temperature_k: f64,
    /// Maximum's position from this layer's left face, m.
    pub maximum_position_m: f64,
    /// Integral(k(T), T_right..T_left)/thickness, W/m2.
    /// This is the spatial mean flux; k at a representative temperature is
    /// not substituted for the original conductivity curve.
    pub recomputed_mean_flux_w_per_m2: f64,
    /// Recovered minus solved mean flux, W/m2.
    pub mean_flux_residual_w_per_m2: f64,
}

/// Completed internally heated nonlinear wall, retaining the original sources.
#[derive(Debug, Clone, PartialEq)]
pub struct NonlinearGeneratedLayeredWallSolution {
    /// Unchanged geometry, films, contacts, source tables and source receipts.
    pub wall: LayeredWall,
    /// Original left and right reservoir temperatures, K.
    pub reservoir_temperatures_k: [f64; 2],
    /// Original common heat-transfer area, m2.
    pub area_m2: f64,
    /// Original signed uniform generation per physical layer, W/m3.
    pub volumetric_generation_w_m3: Vec<f64>,
    /// Signed left-to-right flux at both stack ends, W/m2.
    pub boundary_fluxes_w_per_m2: [f64; 2],
    /// Heat delivered outward to the left and right reservoirs, W.
    pub reservoir_heat_rates_w: [f64; 2],
    /// Signed volume-integrated internal generation, W.
    pub generated_heat_rate_w: f64,
    /// Outward reservoir heat minus internal generation, W.
    pub energy_residual_w: f64,
    /// Temperatures and extrema in physical layer order.
    pub layer_temperatures: Vec<NonlinearGeneratedLayerTemperatures>,
    /// Largest layer/contact/film constitutive flux defect, W/m2.
    pub maximum_flux_residual_w_per_m2: f64,
    /// Reconstructed right reservoir minus the prescribed temperature, K.
    pub right_boundary_residual_k: f64,
    /// Actual shooting trials performed.
    pub iterations: usize,
    /// Final numerical bracket for the left boundary flux, W/m2.
    /// Floating-point diagnostic, not a material/roundoff uncertainty bound.
    pub left_flux_bracket_w_per_m2: [f64; 2],
    /// Unchanged caller-declared domain, work and residual limits.
    pub options: NonlinearGenerationOptions,
    laws: Vec<ConductivityIntegral>,
    pub(super) shapes: Vec<LayerShape>,
}

impl NonlinearGeneratedLayeredWallSolution {
    /// Recover the actual nonlinear profile at a position inside one layer.
    /// The retained source segments are inverted without a mesh, interpolation
    /// between face temperatures, or new material queries. Cancellation is
    /// checked before recovery and at most every 128 crossed source segments.
    pub fn temperature_at(
        &self,
        cx: &Cx<'_>,
        layer_index: usize,
        position_m: f64,
    ) -> Result<f64, ConductionError> {
        checkpoint(cx, layer_index)?;
        let layer = self
            .wall
            .layers
            .get(layer_index)
            .ok_or_else(|| invalid("temperature query layer index is out of range"))?;
        if !position_m.is_finite() || position_m < 0.0 || position_m > layer.thickness_m {
            return Err(invalid(
                "temperature query must lie inside its physical layer",
            ));
        }
        let row = &self.layer_temperatures[layer_index];
        if position_m == 0.0 {
            return Ok(row.face_temperatures_k[0]);
        }
        if position_m == layer.thickness_m {
            return Ok(row.face_temperatures_k[1]);
        }
        let drop = finite(
            "interior Kirchhoff drop",
            self.shapes[layer_index].potential_drop(
                row.face_fluxes_w_per_m2[0],
                self.volumetric_generation_w_m3[layer_index] * self.shapes[layer_index].volume(),
                position_m,
            ),
        )?;
        self.laws[layer_index]
            .shift(cx, row.face_temperatures_k[0], drop)?
            .ok_or_else(|| invalid("temperature query left the solved material interval"))
    }
}

enum Shot {
    TooHot,
    TooCold,
    Complete {
        right_residual: f64,
        maximum_flux_residual: f64,
    },
}

impl LayeredWall {
    /// Solve `-(k(T) T')' = s` with one finite signed source per layer.
    ///
    /// Each source segment is integrated and inverted analytically. The
    /// Kirchhoff profile is quadratic and its possible interior extremum is
    /// checked, so supported faces cannot hide an unsupported hot/cold spot.
    /// Contacts and films use their local flux, including reversed flow.
    ///
    /// A monotone shooting search enforces both declared reservoir conditions.
    /// No trial evaluates conductivity outside the declared common interval.
    /// Only a candidate meeting both temperature and constitutive-flux limits
    /// is returned. Work is bounded by the original layer/knot caps and the
    /// caller's shooting limit, with one reused O(layers) trial buffer.
    ///
    /// Generation is prescribed independently of temperature. No reaction,
    /// Joule feedback, lateral spreading, transient storage or uncertainty
    /// certificate is inferred. The constant-conductivity analytical method
    /// remains available without a shooting search or common-domain admission.
    pub fn solve_temperature_dependent_with_generation(
        &self,
        cx: &Cx<'_>,
        reservoir_temperatures_k: [f64; 2],
        area_m2: f64,
        volumetric_generation_w_m3: &[f64],
        options: NonlinearGenerationOptions,
    ) -> Result<NonlinearGeneratedLayeredWallSolution, ConductionError> {
        self.solve_generated_geometry(
            cx,
            reservoir_temperatures_k,
            area_m2,
            volumetric_generation_w_m3,
            options,
            self.layers.iter().map(|layer| LayerShape::planar(layer.thickness_m)).collect(),
            &self.contact_resistances_m2_k_per_w,
            self.surface_resistances_m2_k_per_w,
        )
    }

    // Shared material integral, source-domain/extremum checks and scalar solve.
    // The private radial caller uses q=Q/Aref and R_reference=R_areal*Aref/Aface.
    // Its public report converts these private normalized quantities to actual
    // heat rates and local fluxes. Planar arithmetic stays area-independent.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn solve_generated_geometry(
        &self,
        cx: &Cx<'_>,
        reservoir_temperatures_k: [f64; 2],
        area_m2: f64,
        volumetric_generation_w_m3: &[f64],
        options: NonlinearGenerationOptions,
        shapes: Vec<LayerShape>,
        contacts: &[f64],
        surfaces: [f64; 2],
    ) -> Result<NonlinearGeneratedLayeredWallSolution, ConductionError> {
        checkpoint(cx, 0)?;
        options.solver.validate()?;
        positive("wall area", area_m2)?;
        for &temperature in &reservoir_temperatures_k {
            positive("reservoir absolute temperature", temperature)?;
        }
        let [low, high] = options.temperature_interval_k;
        positive("admitted minimum material temperature", low)?;
        positive("admitted maximum material temperature", high)?;
        if low >= high {
            return Err(invalid(
                "material temperature interval must have positive width",
            ));
        }
        if volumetric_generation_w_m3.len() != self.layers.len() {
            return Err(invalid("supply exactly one volumetric source per layer"));
        }
        for boundary in 0..2 {
            let temperature = reservoir_temperatures_k[boundary];
            if surfaces[boundary] == 0.0
                && (temperature < low || temperature > high)
            {
                return Err(invalid(
                    "a prescribed face lies outside the material interval",
                ));
            }
        }
        let mut laws = Vec::with_capacity(self.layers.len());
        let mut generation = Vec::with_capacity(self.layers.len());
        let mut total_generation = 0.0;
        for (index, (layer, &source)) in self
            .layers
            .iter()
            .zip(volumetric_generation_w_m3)
            .enumerate()
        {
            checkpoint(cx, index)?;
            finite("volumetric generation", source)?;
            let amount = finite("layer generation per reference area", source * shapes[index].volume())?;
            if source != 0.0 && amount == 0.0 {
                return Err(invalid("nonzero volumetric generation underflowed"));
            }
            total_generation = finite("cumulative generation", total_generation + amount)?;
            generation.push(amount);
            laws.push(ConductivityIntegral::new(
                cx,
                &layer.conductivity,
                low,
                high,
            )?);
        }
        // Both first-layer endpoints must lie in [low,high]. Therefore its
        // integral drop lies in [-capacity,capacity], bounding the inlet flux
        // regardless of generation in later layers or external reservoirs.
        let capacity = laws[0].integral(cx, low, high)?;
        let half_width = capacity / shapes[0].resistance();
        positive("inlet flux bracket half-width", half_width)?;
        let center = -shapes[0].weighted_flux(0.0, generation[0]);
        let mut bracket = [
            finite("inlet flux lower bound", center - half_width)?,
            finite("inlet flux upper bound", center + half_width)?,
        ];
        if bracket[0] >= bracket[1] {
            return Err(invalid("inlet flux bracket is not representable"));
        }
        let mut rows = Vec::with_capacity(self.layers.len());
        let mut last_residual = None;
        for iteration in 1..=options.solver.max_iterations {
            checkpoint(cx, iteration)?;
            let flux = 0.5 * bracket[0] + 0.5 * bracket[1];
            rows.clear();
            match self.shoot_generated(
                cx,
                &laws,
                reservoir_temperatures_k,
                &generation,
                &shapes,
                contacts,
                surfaces,
                flux,
                [low, high],
                &mut rows,
            )? {
                Shot::TooHot => bracket[0] = flux,
                Shot::TooCold => bracket[1] = flux,
                Shot::Complete {
                    right_residual,
                    maximum_flux_residual,
                } => {
                    last_residual = Some(right_residual);
                    let flux_scale = rows
                        .iter()
                        .flat_map(|row| row.face_fluxes_w_per_m2)
                        .map(f64::abs)
                        .fold(0.0_f64, f64::max);
                    if right_residual.abs() <= options.solver.temperature_tolerance_k
                        && maximum_flux_residual <= options.solver.flux_tolerance(flux_scale)?
                    {
                        let right_flux = finite("right boundary flux", flux + total_generation)?;
                        let heat = [
                            finite("left reservoir heat rate", -flux * area_m2)?,
                            finite("right reservoir heat rate", right_flux * area_m2)?,
                        ];
                        let generated_heat_rate_w =
                            finite("total generated heat", total_generation * area_m2)?;
                        if (flux != 0.0 && heat[0] == 0.0)
                            || (right_flux != 0.0 && heat[1] == 0.0)
                            || (total_generation != 0.0 && generated_heat_rate_w == 0.0)
                        {
                            return Err(invalid("nonzero heat rate underflowed"));
                        }
                        let energy_residual_w =
                            finite("energy residual", heat[0] + heat[1] - generated_heat_rate_w)?;
                        checkpoint(cx, self.layers.len())?;
                        return Ok(NonlinearGeneratedLayeredWallSolution {
                            wall: self.clone(),
                            reservoir_temperatures_k,
                            area_m2,
                            volumetric_generation_w_m3: volumetric_generation_w_m3.to_vec(),
                            boundary_fluxes_w_per_m2: [flux, right_flux],
                            reservoir_heat_rates_w: heat,
                            generated_heat_rate_w,
                            energy_residual_w,
                            layer_temperatures: rows,
                            maximum_flux_residual_w_per_m2: maximum_flux_residual,
                            right_boundary_residual_k: right_residual,
                            iterations: iteration,
                            left_flux_bracket_w_per_m2: bracket,
                            options,
                            laws,
                            shapes,
                        });
                    }
                    if right_residual > 0.0 {
                        bracket[0] = flux;
                    } else {
                        bracket[1] = flux;
                    }
                }
            }
            let next = 0.5 * bracket[0] + 0.5 * bracket[1];
            if next <= bracket[0] || next >= bracket[1] {
                break;
            }
        }
        Err(ConductionError::Conductivity {
            what: format!(
                "internally heated nonlinear wall found no solution meeting material-domain and residual limits within {} trials; last supported right-boundary residual {last_residual:?} K",
                options.solver.max_iterations,
            ),
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn shoot_generated(
        &self,
        cx: &Cx<'_>,
        laws: &[ConductivityIntegral],
        reservoirs: [f64; 2],
        generation: &[f64],
        shapes: &[LayerShape],
        contacts: &[f64],
        surfaces: [f64; 2],
        inlet_flux: f64,
        interval: [f64; 2],
        rows: &mut Vec<NonlinearGeneratedLayerTemperatures>,
    ) -> Result<Shot, ConductionError> {
        let mut temperature = finite(
            "left face temperature",
            (-inlet_flux).mul_add(surfaces[0], reservoirs[0]),
        )?;
        let mut maximum_flux_residual = 0.0_f64;
        let mut produced = 0.0;
        for (index, (layer, law)) in self.layers.iter().zip(laws).enumerate() {
            checkpoint(cx, index)?;
            if temperature < interval[0] {
                return Ok(Shot::TooCold);
            }
            if temperature > interval[1] {
                return Ok(Shot::TooHot);
            }
            let left = temperature;
            let flux_in = finite("layer inlet flux", inlet_flux + produced)?;
            produced = finite("accumulated source", produced + generation[index])?;
            let flux_out = finite("layer outlet flux", inlet_flux + produced)?;
            let shape = shapes[index];
            let mean_flux = finite(
                "resistance-weighted reference flux",
                shape.weighted_flux(flux_in, generation[index]),
            )?;
            let drop = finite("layer Kirchhoff drop", mean_flux * shape.resistance())?;
            if mean_flux != 0.0 && drop == 0.0 {
                return Err(invalid("nonzero Kirchhoff drop underflowed"));
            }
            let Some(right) = law.shift(cx, left, drop)? else {
                return Ok(if drop > 0.0 {
                    Shot::TooCold
                } else {
                    Shot::TooHot
                });
            };
            let (mut minimum, mut minimum_x, mut maximum, mut maximum_x) = if left <= right {
                (left, 0.0, right, layer.thickness_m)
            } else {
                (right, layer.thickness_m, left, 0.0)
            };
            if (flux_in < 0.0 && flux_out > 0.0) || (flux_in > 0.0 && flux_out < 0.0) {
                let fraction = finite("stationary layer coordinate", -flux_in / generation[index])?;
                if fraction <= 0.0 || fraction >= 1.0 {
                    return Err(invalid(
                        "interior temperature extremum is not representable",
                    ));
                }
                let x = finite("interior extremum position", shape.stationary_position(fraction))?;
                if x <= 0.0 || x >= layer.thickness_m {
                    return Err(invalid("interior extremum position is not representable"));
                }
                let extremum_drop = finite(
                    "interior extremum Kirchhoff drop",
                    shape.potential_drop(flux_in, generation[index], x),
                )?;
                let Some(value) = law.shift(cx, left, extremum_drop)? else {
                    return Ok(if extremum_drop > 0.0 {
                        Shot::TooCold
                    } else {
                        Shot::TooHot
                    });
                };
                if value < minimum {
                    minimum = value;
                    minimum_x = x;
                }
                if value > maximum {
                    maximum = value;
                    maximum_x = x;
                }
            }
            let magnitude = law.integral(cx, left.min(right), left.max(right))? / shape.resistance();
            let recomputed = if left >= right { magnitude } else { -magnitude };
            let residual = finite("mean Fourier residual", recomputed - mean_flux)?;
            maximum_flux_residual = maximum_flux_residual.max(residual.abs());
            rows.push(NonlinearGeneratedLayerTemperatures {
                face_temperatures_k: [left, right],
                face_fluxes_w_per_m2: [flux_in, flux_out],
                minimum_temperature_k: minimum,
                minimum_position_m: minimum_x,
                maximum_temperature_k: maximum,
                maximum_position_m: maximum_x,
                recomputed_mean_flux_w_per_m2: recomputed,
                mean_flux_residual_w_per_m2: residual,
            });
            temperature = right;
            if let Some(&contact) = contacts.get(index) {
                temperature = finite(
                    "contact downstream temperature",
                    (-flux_out).mul_add(contact, right),
                )?;
                if contact > 0.0 {
                    maximum_flux_residual = maximum_flux_residual.max(
                        finite(
                            "contact flux residual",
                            (right - temperature) / contact - flux_out,
                        )?
                        .abs(),
                    );
                }
            }
        }
        let right_flux = finite("right boundary flux", inlet_flux + produced)?;
        let right_residual = finite(
            "right reservoir residual",
            (-right_flux).mul_add(surfaces[1], temperature)
                - reservoirs[1],
        )?;
        for (jump, resistance, flux) in [
            (
                reservoirs[0] - rows[0].face_temperatures_k[0],
                surfaces[0],
                inlet_flux,
            ),
            (
                temperature - reservoirs[1],
                surfaces[1],
                right_flux,
            ),
        ] {
            if resistance > 0.0 {
                maximum_flux_residual = maximum_flux_residual
                    .max(finite("surface film flux residual", jump / resistance - flux)?.abs());
            }
        }
        Ok(Shot::Complete {
            right_residual,
            maximum_flux_residual,
        })
    }
}
