//! Conservative transient passive heat coupled to transported turbulence.
//!
//! rho cp T_t = d_y[(kappa + rho cp nu_t/Pr_t) T_y] + q_volume.
//! The prepared mechanical endpoint supplies nu_t, a prescribed flux enters
//! the wall half-cell, and the centerline exterior heat flux is zero.

use super::{
    RansChannelGrid, RansFieldState, RansSolver, RansTransportError, RansTransportReport,
    RansTransportWorkspace, checkpoint, volume,
};

/// Heat inputs held constant over one common physical interval.
#[derive(Clone, Copy, Debug, Default)]
pub struct RansHeatLoad<'a> {
    /// Prescribed wall heat flux [W/m²], positive into the channel.
    pub wall_heat_flux_w_m2: f64,
    /// Optional nodal volumetric heating [W/m³], including both boundary nodes.
    /// The dual-cell integral is the trapezoidal integral on the authored grid.
    /// Signed values permit cooling; a nonpositive temperature endpoint refuses.
    pub volumetric_heating_w_m3: Option<&'a [f64]>,
}

/// Mechanical and sensible-heat budgets for one jointly accepted interval.
/// Thermal quantities are per unit wall area; density and heat capacity are
/// constant. Sensible inventory uses zero kelvin as its reference temperature.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RansThermalTransportReport {
    /// The same mechanical report returned by the mechanical-only advance.
    pub mechanical: RansTransportReport,
    /// Integral rho cp T dy before the step [J/m²].
    pub thermal_energy_before_j_m2: f64,
    /// Integral rho cp T dy at the accepted endpoint [J/m²].
    pub thermal_energy_after_j_m2: f64,
    /// Sum of local rho cp volume (T_new - T_old) increments [J/m²].
    /// This avoids subtracting two large accumulated inventories.
    pub thermal_energy_change_j_m2: f64,
    /// Imposed heat input through the wall [W/m²].
    pub wall_heat_flux_w_m2: f64,
    /// Integral of the prescribed volumetric heating across the half-channel [W/m²].
    pub volumetric_heat_input_w_m2: f64,
    /// Accepted local-increment rate minus wall and volume heat input [W/m²].
    /// This measures numerical balance, including endpoint rounding.
    pub thermal_balance_defect_w_m2: f64,
}

impl RansSolver {
    /// Advance momentum, turbulence and passive temperature atomically.
    ///
    /// Thermal diffusivity is kappa/(rho cp) + nu_t/Pr_t, with nu_t from
    /// the prepared mechanical endpoint. Shared face conductances act on the
    /// same stretched grid, including a finite wall and centerline half-cell.
    /// The wall admits the authored flux and the centerline is insulated.
    /// Backward Euler solves temperature increments, preserving a uniform
    /// unforced temperature exactly. The original nine workspace arrays are
    /// reused; stepping allocates no numerical storage.
    ///
    /// Heat is passive: it does not change material properties or produce
    /// buoyancy. Mechanical dissipation is not deposited as heat. No axial
    /// heat removal is supplied, so nonzero net heating continues to change
    /// the thermal inventory and does not define a thermal steady state.
    ///
    /// # Errors
    /// Mechanical admission/errors pass through. Thermal material parameters
    /// must be positive and finite; sources must be finite and correctly sized.
    /// Unrepresentable operators, inventories or nonpositive temperatures
    /// refuse. Every error and cancellation preserves all accepted fields.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_launder_sharma_with_heat(
        &self,
        grid: &RansChannelGrid,
        state: &mut RansFieldState,
        pressure_gradient_pa_m: f64,
        time_step_s: f64,
        heat: RansHeatLoad<'_>,
        work: &mut RansTransportWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RansThermalTransportReport, RansTransportError> {
        checkpoint(0, &mut cancelled)?;
        let n = grid
            .num_cells
            .checked_add(1)
            .ok_or(RansTransportError::InvalidInput("node count overflow"))?;
        if [self.rho, self.cp, self.k_thermal, self.pr_turbulent]
            .iter()
            .any(|&value| !value.is_finite() || value <= 0.0)
            || !heat.wall_heat_flux_w_m2.is_finite()
            || heat
                .volumetric_heating_w_m3
                .is_some_and(|source| source.len() != n)
        {
            return Err(RansTransportError::InvalidInput(
                "thermal material or heat load",
            ));
        }
        if let Some(source) = heat.volumetric_heating_w_m3 {
            for (i, value) in source.iter().enumerate() {
                checkpoint(i, &mut cancelled)?;
                if !value.is_finite() {
                    return Err(RansTransportError::InvalidInput(
                        "nonfinite volumetric heating",
                    ));
                }
            }
        }
        let capacity = self.rho * self.cp;
        let molecular = self.k_thermal / capacity;
        let wall_flux = heat.wall_heat_flux_w_m2 / capacity;
        if !capacity.is_finite()
            || capacity <= 0.0
            || !molecular.is_finite()
            || molecular <= 0.0
            || !wall_flux.is_finite()
        {
            return Err(RansTransportError::Arithmetic {
                equation: "thermal material",
                node: 0,
            });
        }

        let mechanical = self.prepare_launder_sharma(
            grid,
            state,
            pressure_gradient_pa_m,
            time_step_s,
            work,
            &mut cancelled,
        )?;
        let dt = time_step_s;
        let mut volumetric_input = Sum::default();
        // Solve (V/dt - L) delta_T = L T_old + heat inputs. Retaining only
        // differences on the RHS avoids injecting roundoff into a constant T.
        for i in 0..n {
            checkpoint(i, &mut cancelled)?;
            let v = thermal_volume(grid, i);
            let source = heat.volumetric_heating_w_m3.map_or(0.0, |values| values[i]);
            let storage = v / dt;
            let input = v * (source / capacity);
            if !storage.is_finite() || storage <= 0.0 || !input.is_finite() {
                return Err(RansTransportError::Arithmetic {
                    equation: "thermal storage/source",
                    node: i,
                });
            }
            volumetric_input.add(v * source);
            work.lower[i] = 0.0;
            work.diagonal[i] = storage;
            work.upper[i] = 0.0;
            work.rhs[i] = input;
        }
        work.rhs[0] += wall_flux;
        for i in 0..n - 1 {
            checkpoint(i, &mut cancelled)?;
            let dy = grid.y_nodes[i + 1] - grid.y_nodes[i];
            let eddy_face = 0.5 * work.viscosity[i] + 0.5 * work.viscosity[i + 1];
            let conductance = (molecular + eddy_face / self.pr_turbulent) / dy;
            if !conductance.is_finite() || conductance <= 0.0 {
                return Err(RansTransportError::Arithmetic {
                    equation: "thermal diffusion",
                    node: i,
                });
            }
            // One face contributes equal and opposite fluxes, including the
            // wall-adjacent face. The wall itself has thermal storage.
            let flux = conductance * (state.temp_k[i + 1] - state.temp_k[i]);
            work.diagonal[i] += conductance;
            work.diagonal[i + 1] += conductance;
            work.upper[i] = -conductance;
            work.lower[i + 1] = -conductance;
            work.rhs[i] += flux;
            work.rhs[i + 1] -= flux;
        }
        for i in 0..n {
            checkpoint(i, &mut cancelled)?;
            if !work.diagonal[i].is_finite() || work.diagonal[i] <= 0.0 || !work.rhs[i].is_finite()
            {
                return Err(RansTransportError::Arithmetic {
                    equation: "temperature",
                    node: i,
                });
            }
        }
        work.solve("temperature", &mut cancelled)?;

        let mut before = Sum::default();
        let mut after = Sum::default();
        let mut change = Sum::default();
        for i in 0..n {
            checkpoint(i, &mut cancelled)?;
            let old = state.temp_k[i];
            let next = old + work.rhs[i];
            let cell_capacity = capacity * thermal_volume(grid, i);
            if !next.is_finite()
                || next <= 0.0
                || !cell_capacity.is_finite()
                || cell_capacity <= 0.0
            {
                return Err(RansTransportError::Arithmetic {
                    equation: "positive temperature",
                    node: i,
                });
            }
            before.add(cell_capacity * old);
            after.add(cell_capacity * next);
            // Budget the actual representable endpoint, not an increment that
            // might disappear when added to a large background temperature.
            change.add(cell_capacity * (next - old));
            work.rhs[i] = next;
        }
        let net_input = heat.wall_heat_flux_w_m2 + volumetric_input.value;
        let defect = change.value / dt - net_input;
        if [
            before.value,
            after.value,
            change.value,
            volumetric_input.value,
            net_input,
            defect,
        ]
        .iter()
        .any(|value| !value.is_finite())
        {
            return Err(RansTransportError::Arithmetic {
                equation: "thermal budget",
                node: n - 1,
            });
        }
        let report = RansThermalTransportReport {
            mechanical,
            thermal_energy_before_j_m2: before.value,
            thermal_energy_after_j_m2: after.value,
            thermal_energy_change_j_m2: change.value,
            wall_heat_flux_w_m2: heat.wall_heat_flux_w_m2,
            volumetric_heat_input_w_m2: volumetric_input.value,
            thermal_balance_defect_w_m2: defect,
        };
        checkpoint(0, &mut cancelled)?;
        work.publish_mechanics(state);
        state.temp_k.copy_from_slice(&work.rhs);
        Ok(report)
    }
}

fn thermal_volume(grid: &RansChannelGrid, i: usize) -> f64 {
    if i == 0 {
        0.5 * (grid.y_nodes[1] - grid.y_nodes[0])
    } else {
        volume(grid, i)
    }
}

#[derive(Default)]
struct Sum {
    value: f64,
    compensation: f64,
}

impl Sum {
    fn add(&mut self, value: f64) {
        let increment = value - self.compensation;
        let next = self.value + increment;
        self.compensation = (next - self.value) - increment;
        self.value = next;
    }
}

#[cfg(test)]
mod tests;
