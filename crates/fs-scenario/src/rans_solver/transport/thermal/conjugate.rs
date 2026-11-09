//! Monolithic backward-Euler conduction between finite solid layers and the
//! original transported channel. Exact block elimination supplies a positive
//! wall impedance to the existing fluid thermal operator; no coupling iteration
//! or empirical heat-transfer coefficient is needed.

use super::{
    RansChannelGrid, RansFieldState, RansHeatCoupling, RansHeatLoad, RansSolver,
    RansThermalTransportReport, RansTransportError, RansTransportWorkspace, Sum, checkpoint,
};

/// One cell of a stationary, constant-property conducting wall. Layers are
/// ordered from the exterior face toward the fluid interface.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RansWallLayer {
    /// Cell thickness normal to the wall [m].
    pub thickness_m: f64,
    /// Volumetric sensible-heat capacity rho*cp [J/(m³ K)].
    pub volumetric_heat_capacity_j_m3_k: f64,
    /// Normal thermal conductivity [W/(m K)].
    pub conductivity_w_m_k: f64,
}

impl RansWallLayer {
    fn capacity(self) -> f64 {
        self.thickness_m * self.volumetric_heat_capacity_j_m3_k
    }

    fn half_resistance(self) -> f64 {
        0.5 * (self.thickness_m / self.conductivity_w_m_k)
    }
}

/// Admitted wall material and cell-center temperatures. Its clone contains all
/// state needed for deterministic replay alongside the original fluid state.
#[derive(Clone, Debug, PartialEq)]
pub struct RansConjugateWall {
    layers: Vec<RansWallLayer>,
    temperature_k: Vec<f64>,
}

impl RansConjugateWall {
    /// Admit finite positive properties, temperature, capacity and face
    /// conductances before allocating. The final solid half-cell touches the
    /// original fluid node at y=0, without a contact resistance or wall law.
    pub fn try_new(
        layers: &[RansWallLayer],
        temperature_k: &[f64],
        maximum_layers: usize,
    ) -> Result<Self, RansTransportError> {
        let n = layers.len();
        if n == 0 || temperature_k.len() != n {
            return Err(RansTransportError::InvalidInput(
                "wall layer/temperature shape",
            ));
        }
        if n > maximum_layers {
            return Err(RansTransportError::NodeBudget {
                needed: n,
                maximum: maximum_layers,
            });
        }
        n.checked_mul(4 * core::mem::size_of::<f64>())
            .ok_or(RansTransportError::Allocation)?;
        for (i, layer) in layers.iter().enumerate() {
            if [
                layer.thickness_m,
                layer.volumetric_heat_capacity_j_m3_k,
                layer.conductivity_w_m_k,
                layer.capacity(),
                layer.half_resistance(),
                temperature_k[i],
            ]
            .iter()
            .any(|&x| !x.is_finite() || x <= 0.0)
            {
                return Err(RansTransportError::InvalidInput(
                    "wall material or temperature",
                ));
            }
            let resistance = layer.half_resistance()
                + layers.get(i + 1).map_or(0.0, |next| next.half_resistance());
            let conductance = 1.0 / resistance;
            if !conductance.is_finite() || conductance <= 0.0 {
                return Err(RansTransportError::InvalidInput("wall face conductance"));
            }
        }
        let mut owned_layers = Vec::new();
        owned_layers
            .try_reserve_exact(n)
            .map_err(|_| RansTransportError::Allocation)?;
        owned_layers.extend_from_slice(layers);
        let mut owned_temperature = Vec::new();
        owned_temperature
            .try_reserve_exact(n)
            .map_err(|_| RansTransportError::Allocation)?;
        owned_temperature.extend_from_slice(temperature_k);
        Ok(Self {
            layers: owned_layers,
            temperature_k: owned_temperature,
        })
    }

    #[must_use]
    pub fn layers(&self) -> &[RansWallLayer] {
        &self.layers
    }

    #[must_use]
    pub fn temperatures_k(&self) -> &[f64] {
        &self.temperature_k
    }

    // Adaptive trials own a private wall image; only the completed physical
    // window may copy its admitted temperatures back to the caller's wall.
    pub(in crate::rans_solver::transport) fn temperatures_mut(&mut self) -> &mut [f64] {
        &mut self.temperature_k
    }

    pub(in crate::rans_solver::transport) fn numerical_bytes(&self) -> usize {
        self.layers.capacity() * core::mem::size_of::<RansWallLayer>()
            + self.temperature_k.capacity() * core::mem::size_of::<f64>()
    }

    // Conductance from layer i to i+1, or to the fluid interface at the end.
    fn conductance(&self, i: usize) -> f64 {
        1.0 / (self.layers[i].half_resistance()
            + self
                .layers
                .get(i + 1)
                .map_or(0.0, |next| next.half_resistance()))
    }
}

/// Loads held over one common interval. Positive exterior heat enters the
/// outer solid face; the solid/fluid exchange is solved, never prescribed.
#[derive(Clone, Copy, Debug, Default)]
pub struct RansConjugateHeatLoad<'a> {
    /// Exterior wall heat input [W/m²].
    pub exterior_heat_flux_w_m2: f64,
    /// Optional volumetric heating in each solid cell [W/m³].
    pub wall_volumetric_heating_w_m3: Option<&'a [f64]>,
    /// Optional nodal fluid volumetric heating [W/m³].
    pub fluid_volumetric_heating_w_m3: Option<&'a [f64]>,
    /// Whether original physical mechanical losses enter fluid sensible heat.
    pub coupling: RansHeatCoupling,
    /// Optional streamwise Boussinesq feedback from the fluid temperature.
    pub buoyancy: Option<super::RansBoussinesqForce>,
}

/// Original nine fluid buffers and four scalars per wall cell. Every step
/// reuses these arrays; no full-fluid copy, dense matrix or iteration is added.
pub struct RansConjugateWorkspace {
    fluid: RansTransportWorkspace,
    diagonal: Vec<f64>,
    particular: Vec<f64>,
    response: Vec<f64>,
    temperature: Vec<f64>,
}

impl RansConjugateWorkspace {
    /// Allocate the original fluid scratch and four wall arrays under separate
    /// node/layer ceilings. The material and temperature stay with the caller.
    pub fn try_new(
        grid: &RansChannelGrid,
        wall: &RansConjugateWall,
        maximum_fluid_nodes: usize,
        maximum_wall_layers: usize,
    ) -> Result<Self, RansTransportError> {
        let n = wall.layers.len();
        if n > maximum_wall_layers {
            return Err(RansTransportError::NodeBudget {
                needed: n,
                maximum: maximum_wall_layers,
            });
        }
        n.checked_mul(4 * core::mem::size_of::<f64>())
            .ok_or(RansTransportError::Allocation)?;
        let field = || {
            let mut v = Vec::new();
            v.try_reserve_exact(n)
                .map_err(|_| RansTransportError::Allocation)?;
            v.resize(n, 0.0);
            Ok::<_, RansTransportError>(v)
        };
        Ok(Self {
            fluid: RansTransportWorkspace::try_new(grid, maximum_fluid_nodes)?,
            diagonal: field()?,
            particular: field()?,
            response: field()?,
            temperature: field()?,
        })
    }

    /// Numerical scratch payload, excluding headers and the caller's
    /// accepted wall material/temperature and fluid state.
    #[must_use]
    pub fn numerical_bytes(&self) -> usize {
        self.fluid.numerical_bytes()
            + [
                &self.diagonal,
                &self.particular,
                &self.response,
                &self.temperature,
            ]
            .iter()
            .map(|v| v.capacity() * core::mem::size_of::<f64>())
            .sum::<usize>()
    }
}

/// Independent physical inventories and residuals per unit wall area. The
/// common interface flux is computed from the represented endpoint temperatures.
/// Residuals are reported without repairing either inventory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RansConjugateReport {
    /// Original mechanics and fluid thermal budget with the resolved wall flux.
    pub fluid: RansThermalTransportReport,
    /// Sum of solid cell capacity times opening temperature [J/m²].
    pub wall_energy_before_j_m2: f64,
    /// Sum of solid cell capacity times accepted temperature [J/m²].
    pub wall_energy_after_j_m2: f64,
    /// Sum of local capacity-weighted represented temperature increments [J/m²].
    pub wall_energy_change_j_m2: f64,
    /// Authored exterior input, positive into the solid [W/m²].
    pub exterior_heat_flux_w_m2: f64,
    /// Thickness-integrated authored solid volumetric heating [W/m²].
    pub wall_volumetric_heat_input_w_m2: f64,
    /// Actual solid-to-fluid conductance times endpoint temperature jump [W/m²].
    pub interface_heat_flux_w_m2: f64,
    /// Actual interface flux minus the eliminated fluid boundary flux [W/m²].
    pub interface_equation_defect_w_m2: f64,
    /// Wall storage rate minus exterior and volumetric input plus interface loss.
    pub wall_balance_defect_w_m2: f64,
    /// Fluid-plus-wall storage rate minus all external and physical heat inputs.
    pub combined_thermal_balance_defect_w_m2: f64,
}

// Series admittance ab/(a+b), retaining the smaller finite term without a
// product overflow or cancellation of two nearly equal conductances.
fn series(a: f64, b: f64) -> f64 {
    let low = a.min(b);
    low / (1.0 + low / a.max(b))
}

fn arithmetic(node: usize) -> RansTransportError {
    RansTransportError::Arithmetic {
        equation: "conjugate wall",
        node,
    }
}

impl RansSolver {
    /// Advance transported mechanics, fluid heat and a finite layered wall in
    /// one atomic physical step. The thermal subsystem is one backward-Euler
    /// solve, reduced exactly by a tridiagonal Schur complement. Complexity is
    /// O(fluid nodes + wall layers), including original mechanical work.
    ///
    /// Wall cells use physical half-cell series resistances. The fluid retains
    /// its actual wall-node half-volume, endpoint turbulent conductivity and
    /// optional physical dissipative heating. Only fluid viscosity is supplied
    /// by the turbulence closure; no eddy conductivity enters the solid.
    ///
    /// This is fixed-property, stationary, one-dimensional conduction, with
    /// no temperature feedback into momentum/material, contact resistance,
    /// radiation or phase change. Any error/cancellation preserves both states.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_launder_sharma_with_conjugate_wall(
        &self,
        grid: &RansChannelGrid,
        state: &mut RansFieldState,
        wall: &mut RansConjugateWall,
        pressure_gradient_pa_m: f64,
        time_step_s: f64,
        heat: RansConjugateHeatLoad<'_>,
        work: &mut RansConjugateWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RansConjugateReport, RansTransportError> {
        checkpoint(0, &mut cancelled)?;
        let n = wall.layers.len();
        if n == 0
            || !time_step_s.is_finite()
            || time_step_s <= 0.0
            || !heat.exterior_heat_flux_w_m2.is_finite()
            || state.temp_k.is_empty()
            || [
                &work.diagonal,
                &work.particular,
                &work.response,
                &work.temperature,
            ]
            .iter()
            .any(|v| v.len() != n)
            || heat
                .wall_volumetric_heating_w_m3
                .is_some_and(|v| v.len() != n)
        {
            return Err(RansTransportError::InvalidInput(
                "conjugate state, load or workspace",
            ));
        }
        let reference = state.temp_k[0];
        if !reference.is_finite() || reference <= 0.0 {
            return Err(RansTransportError::InvalidInput(
                "fluid interface temperature",
            ));
        }
        let mut wall_source = Sum::default();
        let mut admittance = 0.0;
        // Solve wall temperature relative to the original fluid interface.
        // A uniform unforced temperature therefore has an exactly zero RHS.
        // Positive prefix admittances avoid subtractive Schur cancellation.
        for i in 0..n {
            checkpoint(i, &mut cancelled)?;
            let storage = wall.layers[i].capacity() / time_step_s;
            let source = heat.wall_volumetric_heating_w_m3.map_or(0.0, |v| v[i])
                * wall.layers[i].thickness_m;
            if !storage.is_finite() || storage <= 0.0 || !source.is_finite() {
                return Err(arithmetic(i));
            }
            wall_source.add(source);
            let mut rhs = storage * (wall.temperature_k[i] - reference) + source;
            let previous_admittance = admittance;
            admittance = storage;
            if i == 0 {
                rhs += heat.exterior_heat_flux_w_m2;
            } else {
                let link = wall.conductance(i - 1);
                admittance += series(link, previous_admittance);
                rhs += (link / work.diagonal[i - 1]) * work.particular[i - 1];
            }
            work.diagonal[i] = admittance + wall.conductance(i);
            work.particular[i] = rhs;
            if !admittance.is_finite()
                || admittance <= 0.0
                || !work.diagonal[i].is_finite()
                || work.diagonal[i] <= 0.0
                || !rhs.is_finite()
            {
                return Err(arithmetic(i));
            }
        }
        let conductance = wall.conductance(n - 1);
        let response = conductance / work.diagonal[n - 1];
        let input = response * work.particular[n - 1];
        let impedance = series(conductance, admittance);
        if !input.is_finite() || !impedance.is_finite() || impedance <= 0.0 {
            return Err(arithmetic(n - 1));
        }
        work.particular[n - 1] /= work.diagonal[n - 1];
        work.response[n - 1] = response;
        for i in (0..n - 1).rev() {
            checkpoint(i, &mut cancelled)?;
            let factor = wall.conductance(i) / work.diagonal[i];
            work.particular[i] =
                work.particular[i] / work.diagonal[i] + factor * work.particular[i + 1];
            work.response[i] = factor * work.response[i + 1];
        }
        let mut fluid = self.prepare_launder_sharma_with_heat(
            grid,
            state,
            pressure_gradient_pa_m,
            time_step_s,
            RansHeatLoad {
                wall_heat_flux_w_m2: input,
                volumetric_heating_w_m3: heat.fluid_volumetric_heating_w_m3,
                coupling: heat.coupling,
                buoyancy: heat.buoyancy,
            },
            impedance,
            &mut work.fluid,
            &mut cancelled,
        )?;
        let delta = work.fluid.rhs[0] - reference;
        let mut before = Sum::default();
        let mut after = Sum::default();
        let mut change = Sum::default();
        for i in 0..n {
            checkpoint(i, &mut cancelled)?;
            let next = reference + (work.particular[i] + work.response[i] * delta);
            if !next.is_finite() || next <= 0.0 {
                return Err(arithmetic(i));
            }
            let capacity = wall.layers[i].capacity();
            before.add(capacity * wall.temperature_k[i]);
            after.add(capacity * next);
            change.add(capacity * (next - wall.temperature_k[i]));
            work.temperature[i] = next;
        }
        // Independently recompute physical interface transfer from the actual
        // representable endpoints, exposing any Schur/rounding discrepancy.
        let interface = conductance * (work.temperature[n - 1] - work.fluid.rhs[0]);
        let interface_defect = interface - fluid.wall_heat_flux_w_m2;
        fluid.wall_heat_flux_w_m2 = interface;
        fluid.thermal_balance_defect_w_m2 = fluid.thermal_energy_change_j_m2 / time_step_s
            - interface
            - fluid.volumetric_heat_input_w_m2
            - fluid.molecular_heat_input_w_m2
            - fluid.turbulent_heat_input_w_m2;
        let wall_defect =
            change.value / time_step_s - heat.exterior_heat_flux_w_m2 - wall_source.value
                + interface;
        let combined = (fluid.thermal_energy_change_j_m2 + change.value) / time_step_s
            - heat.exterior_heat_flux_w_m2
            - wall_source.value
            - fluid.volumetric_heat_input_w_m2
            - fluid.molecular_heat_input_w_m2
            - fluid.turbulent_heat_input_w_m2;
        if [
            before.value,
            after.value,
            change.value,
            wall_source.value,
            interface,
            interface_defect,
            fluid.thermal_balance_defect_w_m2,
            wall_defect,
            combined,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(arithmetic(n - 1));
        }
        let report = RansConjugateReport {
            fluid,
            wall_energy_before_j_m2: before.value,
            wall_energy_after_j_m2: after.value,
            wall_energy_change_j_m2: change.value,
            exterior_heat_flux_w_m2: heat.exterior_heat_flux_w_m2,
            wall_volumetric_heat_input_w_m2: wall_source.value,
            interface_heat_flux_w_m2: interface,
            interface_equation_defect_w_m2: interface_defect,
            wall_balance_defect_w_m2: wall_defect,
            combined_thermal_balance_defect_w_m2: combined,
        };
        checkpoint(0, &mut cancelled)?;
        work.fluid.publish_mechanics(state);
        state.temp_k.copy_from_slice(&work.fluid.rhs);
        wall.temperature_k.copy_from_slice(&work.temperature);
        Ok(report)
    }
}
