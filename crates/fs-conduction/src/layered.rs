//! Planar, one-dimensional steady heat flow through an ordered material stack.
//!
//! Each layer retains the existing conductivity table and its source receipts.
//! This rung admits only temperature-independent conductivity; it never freezes
//! a varying source curve at a representative temperature. Contacts and surface
//! films are explicit area-specific resistances, including explicit zero for
//! perfect contact or a prescribed surface temperature. There is no volumetric
//! heating, lateral spreading, radiation, moisture transport or transient claim.

use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_matdb::{ClaimSelection, MaterialCard, PropertyKey, QueryPoint};

use crate::{ConductionError, ConductivityTable};

/// Maximum layers in one bounded planar solve.
pub const MAX_LAYERED_WALL_LAYERS: usize = 256;
/// Maximum total retained conductivity knots in one stack.
pub const MAX_LAYERED_WALL_KNOTS: usize = 65_536;

/// One homogeneous planar layer; coordinates and thickness are in metres.
#[derive(Debug, Clone, PartialEq)]
pub struct PlanarThermalLayer {
    name: String,
    thickness_m: f64,
    conductivity: ConductivityTable,
    material_card: Option<ContentHash>,
}

impl PlanarThermalLayer {
    /// Bind geometry to an admitted conductivity table. Source receipts stay
    /// on the table; a table alone does not establish a material-card identity.
    /// The caller is responsible for matching geometric/process conditions
    /// (for example a source's nominal-thickness coordinate) to this specimen.
    pub fn new(
        name: impl Into<String>,
        thickness_m: f64,
        conductivity: ConductivityTable,
    ) -> Result<Self, ConductionError> {
        let name = name.into();
        if name.trim().is_empty() || name.len() > 128 {
            return Err(invalid("layer name must contain 1..=128 UTF-8 bytes"));
        }
        positive("layer thickness", thickness_m)?;
        Ok(Self {
            name,
            thickness_m,
            conductivity,
            material_card: None,
        })
    }

    /// Resolve one exact material-card conductivity over fully declared query
    /// points, then retain both the card identity and all original receipts.
    /// This adds no temperature or moisture support beyond the selected claim.
    #[allow(clippy::too_many_arguments)]
    pub fn from_card(
        name: impl Into<String>,
        thickness_m: f64,
        card: &MaterialCard,
        property: &PropertyKey,
        temperature_axis: &str,
        points: &[QueryPoint],
        selection: ClaimSelection,
    ) -> Result<Self, ConductionError> {
        let conductivity = match selection {
            ClaimSelection::Policy(policy) => ConductivityTable::from_claims_at_query_points(
                card.claims(),
                property,
                temperature_axis,
                points,
                policy,
            )?,
            ClaimSelection::Pinned(claim) => ConductivityTable::from_claims_pinned_at_query_points(
                card.claims(),
                property,
                temperature_axis,
                points,
                claim,
            )?,
        };
        let mut layer = Self::new(name, thickness_m, conductivity)?;
        layer.material_card = Some(card.content_hash());
        Ok(layer)
    }

    /// Stable caller-supplied layer name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Layer thickness along the heat-flow direction, m.
    #[must_use]
    pub const fn thickness_m(&self) -> f64 {
        self.thickness_m
    }

    /// The unchanged material table, including conditions and usage receipts.
    #[must_use]
    pub const fn conductivity(&self) -> &ConductivityTable {
        &self.conductivity
    }

    /// Exact source card when constructed through from_card.
    #[must_use]
    pub const fn material_card(&self) -> Option<ContentHash> {
        self.material_card
    }
}

/// An ordered, equal-area stack with explicit contact and boundary assumptions.
#[derive(Debug, Clone, PartialEq)]
pub struct LayeredWall {
    layers: Vec<PlanarThermalLayer>,
    contact_resistances_m2_k_per_w: Vec<f64>,
    surface_resistances_m2_k_per_w: [f64; 2],
}

impl LayeredWall {
    /// Admit a nonempty stack. Contacts are in left-to-right order and number
    /// exactly one fewer than the layers. Surface resistances are left, right.
    /// All resistances are in m2 K/W, finite and nonnegative. Explicit zero
    /// means perfect contact (interior) or a prescribed face temperature
    /// (surface); no default resistance or convective correlation is supplied.
    pub fn new(
        layers: Vec<PlanarThermalLayer>,
        contact_resistances_m2_k_per_w: Vec<f64>,
        surface_resistances_m2_k_per_w: [f64; 2],
    ) -> Result<Self, ConductionError> {
        if layers.is_empty() || layers.len() > MAX_LAYERED_WALL_LAYERS {
            return Err(invalid("a planar stack requires 1..=256 layers"));
        }
        if contact_resistances_m2_k_per_w.len() != layers.len() - 1 {
            return Err(invalid(
                "supply exactly one contact resistance per adjacent layer pair",
            ));
        }
        let mut names = std::collections::BTreeSet::new();
        let mut knots = 0usize;
        for layer in &layers {
            if !names.insert(layer.name()) {
                return Err(invalid("layer names must be unique within the stack"));
            }
            knots = knots
                .checked_add(layer.conductivity.knots().len())
                .ok_or_else(|| invalid("conductivity knot count overflowed"))?;
            if knots > MAX_LAYERED_WALL_KNOTS {
                return Err(invalid(
                    "the stack exceeds 65536 retained conductivity knots",
                ));
            }
        }
        for &value in contact_resistances_m2_k_per_w
            .iter()
            .chain(surface_resistances_m2_k_per_w.iter())
        {
            if !value.is_finite() || value < 0.0 {
                return Err(invalid(
                    "contact and surface resistances must be finite and nonnegative",
                ));
            }
        }
        Ok(Self {
            layers,
            contact_resistances_m2_k_per_w,
            surface_resistances_m2_k_per_w,
        })
    }

    /// Ordered layers and original source tables.
    #[must_use]
    pub fn layers(&self) -> &[PlanarThermalLayer] {
        &self.layers
    }

    /// Explicit adjacent-layer contact resistances, m2 K/W.
    #[must_use]
    pub fn contact_resistances_m2_k_per_w(&self) -> &[f64] {
        &self.contact_resistances_m2_k_per_w
    }

    /// Explicit left, right surface-film resistances, m2 K/W.
    #[must_use]
    pub const fn surface_resistances_m2_k_per_w(&self) -> [f64; 2] {
        self.surface_resistances_m2_k_per_w
    }

    /// Solve the constant-conductivity planar problem analytically.
    ///
    /// Positive heat flux travels left to right; reservoirs may be supplied in
    /// either temperature order. Every layer's two actual face temperatures
    /// must be inside that layer's source span. All source conditions remain
    /// those used to construct its table. The returned residual re-evaluates
    /// each Fourier law from the rounded face temperatures; it is a numerical
    /// diagnostic, not an uncertainty or discretization-error certificate.
    ///
    /// # Errors
    /// Cancellation, invalid area/absolute temperatures, varying conductivity,
    /// unsupported source temperatures, or nonrepresentable arithmetic. No
    /// partial candidate or mutation is published on a refusal.
    pub fn solve(
        &self,
        cx: &Cx<'_>,
        reservoir_temperatures_k: [f64; 2],
        area_m2: f64,
    ) -> Result<LayeredWallSolution, ConductionError> {
        checkpoint(cx, 0)?;
        positive("wall area", area_m2)?;
        for &temperature in &reservoir_temperatures_k {
            positive("reservoir absolute temperature", temperature)?;
        }

        let mut resistances = Vec::with_capacity(self.layers.len());
        let mut total = self.surface_resistances_m2_k_per_w[0];
        for (index, layer) in self.layers.iter().enumerate() {
            checkpoint(cx, index)?;
            if layer.conductivity.is_temperature_dependent() {
                return Err(ConductionError::Conductivity {
                    what: format!(
                        "layer {:?} varies with temperature; the constant planar solve cannot freeze it",
                        layer.name
                    ),
                });
            }
            let conductivity = layer.conductivity.knots()[0].1;
            let resistance = layer.thickness_m / conductivity;
            positive("layer area-specific resistance", resistance)?;
            resistances.push(resistance);
            total = finite("total resistance", total + resistance)?;
            if let Some(&contact) = self.contact_resistances_m2_k_per_w.get(index) {
                total = finite("total resistance", total + contact)?;
            }
        }
        total = finite(
            "total resistance",
            total + self.surface_resistances_m2_k_per_w[1],
        )?;
        positive("total resistance", total)?;
        let transmittance = 1.0 / total;
        positive("thermal transmittance", transmittance)?;
        let delta_t = reservoir_temperatures_k[0] - reservoir_temperatures_k[1];
        let flux = finite("heat flux", delta_t / total)?;
        let heat = finite("heat rate", flux * area_m2)?;
        if delta_t != 0.0 && (flux == 0.0 || heat == 0.0) {
            return Err(invalid(
                "nonzero heat transfer underflowed its representation",
            ));
        }

        // Evaluate every face from the same reservoir pair and cumulative
        // resistance, avoiding accumulated layer-by-layer temperature drift.
        let face = |r: f64| -> f64 {
            if r == 0.0 {
                reservoir_temperatures_k[0]
            } else if r == total {
                reservoir_temperatures_k[1]
            } else {
                (-delta_t).mul_add(r / total, reservoir_temperatures_k[0])
            }
        };
        let mut cumulative = self.surface_resistances_m2_k_per_w[0];
        let mut layer_results = Vec::with_capacity(self.layers.len());
        let mut maximum_flux_residual = 0.0f64;
        for (index, (layer, &resistance)) in self.layers.iter().zip(&resistances).enumerate() {
            checkpoint(cx, index)?;
            let left = face(cumulative);
            cumulative = finite("cumulative resistance", cumulative + resistance)?;
            let right = face(cumulative);
            for temperature in [left, right] {
                positive("layer face absolute temperature", temperature)?;
                layer.conductivity.span().check(temperature)?;
            }
            let actual_flux = finite(
                "recomputed layer heat flux",
                layer.conductivity.eval(left)? * ((left - right) / layer.thickness_m),
            )?;
            let residual = finite("layer flux residual", actual_flux - flux)?;
            maximum_flux_residual = maximum_flux_residual.max(residual.abs());
            layer_results.push(PlanarLayerTemperatures {
                left_k: left,
                right_k: right,
                resistance_m2_k_per_w: resistance,
                recomputed_flux_w_per_m2: actual_flux,
                flux_residual_w_per_m2: residual,
            });
            if let Some(&contact) = self.contact_resistances_m2_k_per_w.get(index) {
                cumulative = finite("cumulative resistance", cumulative + contact)?;
            }
        }
        let reconstructed_right = (-flux).mul_add(
            self.surface_resistances_m2_k_per_w[1],
            layer_results
                .last()
                .ok_or_else(|| invalid("missing layer result"))?
                .right_k,
        );
        let boundary_residual = finite(
            "right boundary temperature residual",
            reconstructed_right - reservoir_temperatures_k[1],
        )?;
        checkpoint(cx, self.layers.len())?;
        Ok(LayeredWallSolution {
            wall: self.clone(),
            reservoir_temperatures_k,
            area_m2,
            heat_flux_w_per_m2: flux,
            heat_rate_w: heat,
            thermal_transmittance_w_per_m2_k: transmittance,
            total_resistance_m2_k_per_w: total,
            layer_temperatures: layer_results,
            maximum_flux_residual_w_per_m2: maximum_flux_residual,
            right_boundary_residual_k: boundary_residual,
        })
    }
}

/// Both sides of one layer; adjacent rows retain any explicit contact jump.
#[derive(Debug, Clone, PartialEq)]
pub struct PlanarLayerTemperatures {
    /// Temperature on the layer's left face, K.
    pub left_k: f64,
    /// Temperature on the layer's right face, K.
    pub right_k: f64,
    /// Layer thickness divided by its conductivity, m2 K/W.
    pub resistance_m2_k_per_w: f64,
    /// Fourier flux recomputed from the returned rounded face temperatures.
    pub recomputed_flux_w_per_m2: f64,
    /// Recomputed Fourier flux minus the common solved flux, W/m2.
    pub flux_residual_w_per_m2: f64,
}

/// A completed nominal steady solve, retaining geometry, assumptions and sources.
/// No material uncertainty has been propagated through this calculation.
#[derive(Debug, Clone, PartialEq)]
pub struct LayeredWallSolution {
    /// Unchanged input stack, including exact source tables and receipts.
    pub wall: LayeredWall,
    /// Left, right reservoir absolute temperatures, K.
    pub reservoir_temperatures_k: [f64; 2],
    /// Common heat-transfer area, m2.
    pub area_m2: f64,
    /// Signed common left-to-right heat flux density, W/m2.
    pub heat_flux_w_per_m2: f64,
    /// Signed left-to-right total heat flow, W.
    pub heat_rate_w: f64,
    /// Nominal inverse total area-specific resistance, W/(m2 K).
    pub thermal_transmittance_w_per_m2_k: f64,
    /// All layer, contact and surface resistances in series, m2 K/W.
    pub total_resistance_m2_k_per_w: f64,
    /// Layer face temperatures in the original stack order.
    pub layer_temperatures: Vec<PlanarLayerTemperatures>,
    /// Maximum absolute Fourier-law residual from the returned temperatures.
    pub maximum_flux_residual_w_per_m2: f64,
    /// Reconstructed right reservoir minus the prescribed right reservoir, K.
    pub right_boundary_residual_k: f64,
}

fn invalid(message: &str) -> ConductionError {
    ConductionError::Config {
        parameter: "layered-wall",
        what: message.to_string(),
    }
}

fn finite(name: &str, value: f64) -> Result<f64, ConductionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(&format!("{name} is not finite")))
    }
}

fn positive(name: &str, value: f64) -> Result<(), ConductionError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(invalid(&format!(
            "{name} must be finite and strictly positive"
        )))
    }
}

fn checkpoint(cx: &Cx<'_>, at: usize) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled {
        stage: "layered-wall",
        at,
    })
}
