//! SI boundary exchange on the original planar double-population heat solver.
//!
//! The slab's retained density, Cp, dx/dt and signed temperature coordinate
//! supply every conversion. A film coefficient is an explicit physical input;
//! this adapter does not select a correlation or debit a finite reservoir.

use super::{MaterialLbmError, MaterialSlab, ThermalLbm, ThermalLbm3};
use crate::thermal::boundary::{ThermalWall, ThermalWallError, ThermalWallExchange};

/// Prescribed thermal law on one stationary, impermeable wall or surface link.
/// Bottom/top orientation never changes the sign convention: positive flux
/// enters the fluid. Wall area comes from the existing lattice, not this law.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MaterialThermalWall {
    /// Prescribed positive absolute wall temperature [K].
    Temperature(f64),
    /// Exactly zero heat exchange with the exterior.
    Adiabatic,
    /// Prescribed signed inward heat-flux density [W/m²].
    HeatFlux(f64),
    /// Exterior Newton cooling/heating: `q_in = h (T_ambient - T_wall)`.
    /// The wall trace is solved by the native reflection rule; it is not the
    /// temperature of the adjacent fluid cell.
    Convective {
        /// Finite nonnegative film coefficient [W/(m² K)]; zero insulates.
        h_w_m2_k: f64,
        /// Positive absolute temperature of the prescribed reservoir [K].
        ambient_temperature_k: f64,
    },
}

/// Measured exchange on one wall during the last accepted thermal step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MaterialWallExchange {
    /// Net heat entering the represented fluid through this wall [J].
    /// This is the sum of actual incoming-minus-outgoing thermal populations,
    /// converted once using the represented cell volume and heat capacity.
    pub heat_in_j: f64,
    /// Area-mean temperature of the native halfway wall trace [K].
    pub mean_temperature_k: f64,
    /// Smallest local halfway wall temperature [K]. A negative normalized
    /// temperature scale reverses which native endpoint supplies this bound.
    pub min_temperature_k: f64,
    /// Largest local halfway wall temperature [K], with the same signed
    /// coordinate conversion as the mean and minimum.
    pub max_temperature_k: f64,
}

fn invalid(reason: &'static str) -> MaterialLbmError {
    MaterialLbmError::Invalid(reason)
}

fn wall_error(error: ThermalWallError) -> MaterialLbmError {
    match error {
        ThermalWallError::Invalid(reason) => invalid(reason),
    }
}

impl<L> MaterialSlab<L> {
    fn normalized_wall_temperature(&self, temperature_k: f64) -> Result<f64, MaterialLbmError> {
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(invalid(
                "thermal wall and reservoir temperatures must be positive finite kelvin",
            ));
        }
        let difference = temperature_k - self.top_temperature_k;
        let normalized = difference / self.temperature_difference_k;
        if !normalized.is_finite() || (difference != 0.0 && normalized == 0.0) {
            return Err(invalid(
                "wall temperature is not representable in the retained thermal coordinate",
            ));
        }
        Ok(normalized)
    }

    fn thermal_wall_inputs<const N: usize>(
        &self,
        walls: [MaterialThermalWall; N],
    ) -> Result<[ThermalWall; N], MaterialLbmError> {
        let units = self.force_units;
        // Heat-flux density times dt/(rho Cp dx) gives the scalar entering
        // ONE normal lattice link per step. DeltaT is signed for q, while H
        // multiplies a difference in the same temperature coordinate.
        let transfer = |value: f64| -> Result<f64, MaterialLbmError> {
            let lattice = ((value / units.transport.density_kg_m3())
                / units.transport.specific_heat_j_kg_k()) / units.velocity_scale;
            if !value.is_finite() || !lattice.is_finite() || (value != 0.0 && lattice == 0.0) {
                return Err(invalid(
                    "thermal wall transfer is not representable in lattice units",
                ));
            }
            Ok(lattice)
        };
        let mut result = [ThermalWall::Adiabatic; N];
        for (output, input) in result.iter_mut().zip(walls) {
            *output = match input {
                MaterialThermalWall::Temperature(temperature) => {
                    ThermalWall::Temperature(self.normalized_wall_temperature(temperature)?)
                }
                MaterialThermalWall::Adiabatic => ThermalWall::Adiabatic,
                MaterialThermalWall::HeatFlux(flux) => {
                    let heat = transfer(flux)? / self.temperature_difference_k;
                    if !heat.is_finite() || (flux != 0.0 && heat == 0.0) {
                        return Err(invalid(
                            "wall heat flux is not representable in the retained thermal coordinate",
                        ));
                    }
                    ThermalWall::HeatFlux(heat)
                }
                MaterialThermalWall::Convective {
                    h_w_m2_k,
                    ambient_temperature_k,
                } => {
                    if !h_w_m2_k.is_finite() || h_w_m2_k < 0.0 {
                        return Err(invalid(
                            "thermal wall film coefficient must be finite and nonnegative",
                        ));
                    }
                    ThermalWall::Convective {
                        transfer: transfer(h_w_m2_k)?,
                        ambient: self.normalized_wall_temperature(ambient_temperature_k)?,
                    }
                }
            };
            output.validate().map_err(wall_error)?;
        }
        Ok(result)
    }

    fn spatial_thermal_wall_inputs(
        &self, walls: [Vec<MaterialThermalWall>; 2], links: usize,
    ) -> Result<[Vec<ThermalWall>; 2], MaterialLbmError> {
        if links == 0 || walls.iter().any(|wall| wall.len() != links) {
            return Err(invalid("one physical thermal law per wall link is required"));
        }
        let mut result = [Vec::new(), Vec::new()];
        for (output, input) in result.iter_mut().zip(walls) {
            *output = input.into_iter().map(|wall| {
                self.thermal_wall_inputs([wall]).map(|native| native[0])
            }).collect::<Result<_, _>>()?;
        }
        Ok(result)
    }

    fn physical_wall_exchange(
        &self,
        exchange: Option<[ThermalWallExchange; 2]>,
        depth_m: f64,
    ) -> Result<Option<[MaterialWallExchange; 2]>, MaterialLbmError> {
        if !depth_m.is_finite() || depth_m <= 0.0 {
            return Err(invalid(
                "thermal wall reporting needs a positive finite represented depth",
            ));
        }
        let Some(exchange) = exchange else {
            return Ok(None);
        };
        let units = self.force_units;
        let volume = (units.cell_size_m * units.cell_size_m) * depth_m;
        let energy_unit = ((units.transport.density_kg_m3() * volume)
            * units.transport.specific_heat_j_kg_k())
            * self.temperature_difference_k;
        if !volume.is_finite() || volume <= 0.0 || !energy_unit.is_finite() || energy_unit == 0.0 {
            return Err(invalid("thermal wall energy unit is unrepresentable"));
        }
        let convert =
            |wall: ThermalWallExchange| -> Result<MaterialWallExchange, MaterialLbmError> {
                let heat_in_j = wall.heat_in * energy_unit;
                if !heat_in_j.is_finite() || (wall.heat_in != 0.0 && heat_in_j == 0.0) {
                    return Err(invalid(
                        "measured thermal wall heat is unrepresentable in joules",
                    ));
                }
                let [minimum, maximum] = if self.temperature_difference_k > 0.0 {
                    [wall.min_temperature, wall.max_temperature]
                } else {
                    [wall.max_temperature, wall.min_temperature]
                };
                Ok(MaterialWallExchange {
                    heat_in_j,
                    mean_temperature_k: self.temperature_si(wall.mean_temperature)?,
                    min_temperature_k: self.temperature_si(minimum)?,
                    max_temperature_k: self.temperature_si(maximum)?,
                })
            };
        Ok(Some([convert(exchange[0])?, convert(exchange[1])?]))
    }
}

impl MaterialSlab<ThermalLbm> {
    /// Set bottom/top thermal walls in SI using this slab's retained material.
    /// All validation and conversions precede mutation. No populations or
    /// temperature coordinates are reset, so insulation and heating schedules
    /// can be changed without erasing fluid or thermal history.
    ///
    /// These laws apply to ordinary, external-force and checked source steps.
    /// The original small-density-variation and constant-storage assumptions
    /// still apply; a prescribed reservoir is not a finite heat-storage body.
    pub fn set_thermal_walls(
        &mut self,
        walls: [MaterialThermalWall; 2],
    ) -> Result<(), MaterialLbmError> {
        if self.lattice.t_bottom != 1.0 || self.lattice.t_top != 0.0 {
            return Err(invalid(
                "material thermal walls refuse independently changed normalized wall temperatures",
            ));
        }
        let walls = self.thermal_wall_inputs(walls)?;
        self.lattice.set_thermal_walls(walls).map_err(wall_error)
    }

    /// Set local bottom/top surface laws in SI, one entry per x link. Both
    /// arrays must contain `lattice.grid.nx` entries. Shape admission and all
    /// signed temperature/heat conversions precede mutation of the slab.
    /// Local patches use the same conservative reflection as uniform walls.
    pub fn set_spatial_thermal_walls(
        &mut self, walls: [Vec<MaterialThermalWall>; 2],
    ) -> Result<(), MaterialLbmError> {
        if self.lattice.t_bottom != 1.0 || self.lattice.t_top != 0.0 {
            return Err(invalid(
                "material thermal walls refuse independently changed normalized wall temperatures",
            ));
        }
        let walls = self.spatial_thermal_wall_inputs(walls, self.lattice.grid.nx)?;
        self.lattice.set_spatial_thermal_walls(walls).map_err(wall_error)
    }

    /// Last measured wall heat and trace temperatures, with explicit 2-D
    /// extrusion depth [m]. Changing depth scales joules, never temperature or
    /// flow. `None` means no exchange is available after initialization or a
    /// wall-law change. The observation does not advance or modify the solver.
    pub fn last_wall_exchange(
        &self,
        depth_m: f64,
    ) -> Result<Option<[MaterialWallExchange; 2]>, MaterialLbmError> {
        self.physical_wall_exchange(self.lattice.last_wall_exchange(), depth_m)
    }
}

impl MaterialSlab<ThermalLbm3> {
    /// Three-dimensional counterpart with the same captured SI conversion.
    /// D3Q7 has a different thermal relaxation rate, but the physical net wall
    /// heat does not acquire another lattice-weight or relaxation factor.
    pub fn set_thermal_walls(
        &mut self,
        walls: [MaterialThermalWall; 2],
    ) -> Result<(), MaterialLbmError> {
        if self.lattice.t_bottom != 1.0 || self.lattice.t_top != 0.0 {
            return Err(invalid(
                "material thermal walls refuse independently changed normalized wall temperatures",
            ));
        }
        let walls = self.thermal_wall_inputs(walls)?;
        self.lattice.set_thermal_walls(walls).map_err(wall_error)
    }

    /// Set local bottom/top surface laws in SI, in x-fastest (`x + nx*y`)
    /// order. Each array must contain `lattice.grid.nx * lattice.grid.ny`
    /// entries. Every patch uses the captured material and thermal coordinate;
    /// invalid shape or conversion preserves all prior laws and reports.
    pub fn set_spatial_thermal_walls(
        &mut self, walls: [Vec<MaterialThermalWall>; 2],
    ) -> Result<(), MaterialLbmError> {
        if self.lattice.t_bottom != 1.0 || self.lattice.t_top != 0.0 {
            return Err(invalid(
                "material thermal walls refuse independently changed normalized wall temperatures",
            ));
        }
        let links = self.lattice.grid.nx.checked_mul(self.lattice.grid.ny)
            .ok_or(invalid("wall plane overflows"))?;
        let walls = self.spatial_thermal_wall_inputs(walls, links)?;
        self.lattice.set_spatial_thermal_walls(walls).map_err(wall_error)
    }

    /// Last measured wall heat [J] and mean trace temperatures [K], using the
    /// actual physical cell volume `dx³`. No extrusion depth is inferred.
    pub fn last_wall_exchange(&self) -> Result<Option<[MaterialWallExchange; 2]>, MaterialLbmError> {
        self.physical_wall_exchange(self.lattice.last_wall_exchange(), self.force_units.cell_size_m)
    }
}
