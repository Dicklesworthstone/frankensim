//! Finite moving-counterface conduction around a material-attached shear film.
//!
//! The film stays at its original mechanical/material patch. Its finite mass,
//! elastic/plastic history and plastic heat remain with the existing layer
//! owner. Only the opposite solid's material-coordinate footprint moves. This
//! does not transport film mass or strain history, move a mechanical contact
//! map, or supply an interfacial resistance. Both surfaces remain adiabatic to
//! ambient: changing exposed-area convection is outside this approximation.

use super::*;
use dry::DryThermalCell;

/// Signed internal Fourier transfer at one candidate thermal/mechanical state.
/// Positive watts enter a cell. Plastic heat is absent: its original layer
/// owner deposits that loss once. The small signed sum is exposed roundoff,
/// never an external supply or a replacement for the coupled energy gate.
#[derive(Clone, Debug, PartialEq)]
pub struct RetainedFilmHeatFrame {
    pub sources_w: BTreeMap<CellId, f64>,
    pub internal_transfer_residual_w: f64,
}

/// One fixed original shear-film cell and a sliding finite receiving strip.
/// The producer owns the material correspondence and must derive each footprint
/// from the SAME trial's actual relative displacement. No motion is inferred
/// from a commanded drive, and no receiver history is stored in this adapter.
#[derive(Clone, Debug)]
pub struct RetainedFilmHeat<'a> {
    port: usize,
    film: DryThermalCell<'a>,
    surface: Vec<DryThermalCell<'a>>,
    receiving_cells: BTreeSet<CellId>,
    area_m2: f64,
}

impl<'a> RetainedFilmHeat<'a> {
    /// Bind original film geometry, caloric cell, material and finite receivers.
    /// Existing links between the film and these cells would duplicate this
    /// conduction and refuse. Fixed links to its attached substrate may remain.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        layer: &ResolvedShearLayer,
        binding: LayerThermalCell<'a>,
        film: DryThermalCell<'a>,
        surface: Vec<DryThermalCell<'a>>,
        thermal: &DepositionSimulation<'_>,
        maximum_cells: usize,
    ) -> Result<Self, LayerThermalError> {
        let bad = LayerThermalError::Invalid;
        if surface.is_empty() || maximum_cells == 0 || surface.len() >= maximum_cells {
            return Err(bad(
                "retained film conduction needs bounded finite receiving cells",
            ));
        }
        if film.cell != binding.cell
            || film.material.material.content_hash() != binding.material.content_hash()
            || film.transport.card_identity() != layer.resolved().card_identity()
            || film.material.temperature_axis != binding.temperature_axis
            || !same_size(2.0 * film.half_depth_m, layer.thickness_m())
            || !same_size(
                (film.interval_m[1] - film.interval_m[0]) * film.width_m,
                layer.area_m2(),
            )
        {
            return Err(bad(
                "retained film heat must keep the original layer material, area, thickness and cell",
            ));
        }
        let cell = thermal
            .cells()
            .get(&film.cell)
            .ok_or(bad("retained film thermal cell is inactive"))?;
        if !thermal.uses_material_law(cell.material_id, binding.caloric)
            || thermal.temperature_k(film.cell)? != layer.temperature_kelvin()
        {
            return Err(bad(
                "retained film must start from its accepted caloric owner and temperature",
            ));
        }
        film.conductance::<Infallible>(thermal, layer.area_m2())?;
        let mut ids = BTreeSet::new();
        let mut previous = None;
        for body in &surface {
            if body.cell == film.cell
                || !ids.insert(body.cell)
                || previous.is_some_and(|end| end != body.interval_m[0])
            {
                return Err(bad(
                    "retained film receivers must be distinct contiguous material cells",
                ));
            }
            body.conductance::<Infallible>(thermal, layer.area_m2())?;
            previous = Some(body.interval_m[1]);
        }
        if thermal.links().values().any(|link| {
            link.a == film.cell && ids.contains(&link.b)
                || link.b == film.cell && ids.contains(&link.a)
        }) {
            return Err(bad(
                "an existing link duplicates retained-film counterface conduction",
            ));
        }
        Ok(Self {
            port: binding.port,
            film,
            surface,
            receiving_cells: ids,
            area_m2: layer.area_m2(),
        })
    }

    #[must_use]
    pub const fn port(&self) -> usize {
        self.port
    }
    #[must_use]
    pub const fn film_cell(&self) -> CellId {
        self.film.cell
    }
    #[must_use]
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }

    /// Original half-film and half-body Fourier resistances act in series on
    /// the actual overlap. Uncovered cells receive no heat, but keep their
    /// original finite enthalpy for a later revisit. All source claims refresh
    /// at the candidate cell temperatures without extending their domains.
    ///
    /// Add these signed rates to the existing `StepCoefficients::sources_w` in
    /// the joint thermal trial. The producer must preserve actual increments
    /// when merging other sources and retain the original global energy gate.
    pub fn heat_rates_w(
        &self,
        endpoint: &DepositionSimulation<'_>,
        mut interval_m: [f64; 2],
        width_m: f64,
    ) -> Result<RetainedFilmHeatFrame, LayerThermalError> {
        let bad = LayerThermalError::Invalid;
        if !width_m.is_finite() || width_m <= 0.0 || interval_m.iter().any(|x| !x.is_finite()) {
            return Err(bad(
                "retained film footprint requires finite coordinates and positive width",
            ));
        }
        // Snap representational translation roundoff, not physical coverage.
        for edge in &mut interval_m {
            for boundary in self.surface.iter().flat_map(|body| body.interval_m) {
                if (*edge - boundary).abs() <= 8.0 * f64::EPSILON * edge.abs().max(boundary.abs()) {
                    *edge = boundary;
                    break;
                }
            }
        }
        let length = interval_m[1] - interval_m[0];
        let area_roundoff = 64.0
            * f64::EPSILON
            * width_m
            * interval_m[0].abs().max(interval_m[1].abs()).max(length);
        if length <= 0.0
            || !length.is_finite()
            || (length * width_m - self.area_m2).abs() > area_roundoff
            || interval_m[0] < self.surface[0].interval_m[0]
            || interval_m[1] > self.surface.last().expect("nonempty strip").interval_m[1]
        {
            return Err(bad(
                "moving counterface must cover the original retained film area inside its finite strip",
            ));
        }
        if endpoint.links().values().any(|link| {
            link.a == self.film.cell && self.receiving_cells.contains(&link.b)
                || link.b == self.film.cell && self.receiving_cells.contains(&link.a)
        }) {
            return Err(bad(
                "an existing link duplicates retained-film counterface conduction",
            ));
        }
        let mut sources_w = BTreeMap::new();
        for body in &self.surface {
            if body.width_m < width_m {
                return Err(bad("retained film footprint exceeds receiving strip width"));
            }
            // Admit uncovered cells too: their source/history remains live.
            body.conductance::<Infallible>(endpoint, self.area_m2)?;
            let overlap =
                body.interval_m[1].min(interval_m[1]) - body.interval_m[0].max(interval_m[0]);
            if overlap > 0.0 {
                let area = number(overlap * width_m)?;
                let (film_g, film_t) = self.film.conductance::<Infallible>(endpoint, area)?;
                let (body_g, body_t) = body.conductance::<Infallible>(endpoint, area)?;
                // Harmonic sum without overflowing either reciprocal.
                let (small, large) = (film_g.min(body_g), film_g.max(body_g));
                let conductance = number(small / (1.0 + small / large))?;
                if conductance <= 0.0 {
                    return Err(bad(
                        "positive retained-film conductance is not representable",
                    ));
                }
                let watts = number(conductance * (film_t - body_t))?;
                sources_w.insert(body.cell, watts);
            }
        }
        let mut debit = 0.0;
        for watts in sources_w.values() {
            debit = number(debit + watts)?;
        }
        sources_w.insert(self.film.cell, -debit);
        let (mut residual, mut absolute) = (0.0, 0.0);
        for watts in sources_w.values() {
            residual = number(residual + watts)?;
            absolute = number(absolute + watts.abs())?;
        }
        let roundoff = 32.0 * f64::EPSILON * sources_w.len() as f64 * absolute;
        if !roundoff.is_finite() || residual.abs() > roundoff {
            return Err(bad(
                "retained-film conduction does not conserve internal heat within summation roundoff",
            ));
        }
        Ok(RetainedFilmHeatFrame {
            sources_w,
            internal_transfer_residual_w: residual,
        })
    }
}
