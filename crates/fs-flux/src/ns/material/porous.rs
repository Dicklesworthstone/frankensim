//! Molecular source viscosity and physical units for native Brinkman momentum.

use super::{MaterialFluxError, MaterialFluxSystem, MaterialStokesBoundary};
use crate::ns::FluxSolution;
use crate::ns::boundary::{StokesBoundary, StokesPressureReference, StokesSolution};
use crate::ns::porous::BrinkmanCell;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};

/// Accepted native flow, original transport snapshot and actual pore dissipation.
pub struct MaterialBrinkmanSolution {
    flow: StokesSolution,
    transport: FluidTransportState,
    cells: Vec<BrinkmanCell>,
    pressure_pa: Vec<f64>,
    cell_darcy_power_w_per_m: Vec<f64>,
    darcy_power_w_per_m: f64,
}

impl MaterialBrinkmanSolution {
    /// Native BDM1/P0 result, including independent equation residuals.
    #[must_use]
    pub const fn discrete(&self) -> &FluxSolution {
        &self.flow.field
    }
    /// Physical superficial-velocity moments [m/s].
    #[must_use]
    pub fn velocity_dofs(&self) -> &[f64] {
        &self.flow.field.x[..self.flow.field.n_u]
    }
    /// Pressure [Pa] in the declared boundary datum, never an inferred EOS state.
    #[must_use]
    pub fn pressure_pa(&self) -> &[f64] {
        &self.pressure_pa
    }
    /// The traction datum or all-velocity mean-zero gauge.
    #[must_use]
    pub const fn pressure_reference(&self) -> StokesPressureReference {
        self.flow.pressure_reference
    }
    /// Original immutable material snapshot used by both stress and pore drag.
    #[must_use]
    pub const fn transport(&self) -> &FluidTransportState {
        &self.transport
    }
    /// Exact geometry and explicitly supplied effective-viscosity closures.
    #[must_use]
    pub fn cells(&self) -> &[BrinkmanCell] {
        &self.cells
    }
    /// `integral(mu u.K^-1.u)` per triangle, per out-of-plane depth [W/m].
    /// Physical dissipation is reported once; no thermal solver is updated here.
    #[must_use]
    pub fn cell_darcy_power_w_per_m(&self) -> &[f64] {
        &self.cell_darcy_power_w_per_m
    }
    /// Sum of the original cell Darcy-power quadratures [W/m].
    #[must_use]
    pub const fn darcy_power_w_per_m(&self) -> f64 {
        self.darcy_power_w_per_m
    }
}

impl<F: FluidTransport + ?Sized> MaterialFluxSystem<'_, '_, '_, F> {
    /// Solve the existing mixed-boundary momentum operator with intrinsic pore
    /// drag and an explicitly chosen effective-viscosity ratio in every cell.
    ///
    /// Force is N/m³, boundary velocity m/s and TOTAL traction Pa. The same
    /// source-owned molecular viscosity and reference density set Darcy drag,
    /// effective shear stress, pressure conversion and power. Geometry is m²,
    /// not a fluid-dependent hydraulic conductivity or a membrane permeance.
    /// Source queries remain on the original bound material; no gas state or
    /// property fit is manufactured. Supplied cell laws alone carry no measured
    /// material authority, and do not infer saturation or applicability bounds.
    ///
    /// This is a steady, rigid-frame incompressible closure. It does not deposit
    /// drag work into heat or add poromechanical, capillary or gas-slip physics.
    ///
    /// # Errors
    /// Preserves material and native Stokes refusals, rejecting nonfinite loads,
    /// unrepresentable conversions or unconverged flow before returning a result.
    pub fn solve_brinkman_mixed(
        &self,
        cells: &[BrinkmanCell],
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialStokesBoundary<'_>],
    ) -> Result<MaterialBrinkmanSolution, MaterialFluxError> {
        let density = self.transport.density_kg_m3();
        let scaled_tractions: Vec<_> = boundaries
            .iter()
            .map(|boundary| {
                move |point: [f64; 2], normal: [f64; 2]| match boundary {
                    MaterialStokesBoundary::TotalTraction { value, .. } => {
                        let physical = value(point, normal);
                        [physical[0] / density, physical[1] / density]
                    }
                    MaterialStokesBoundary::Velocity { .. } => {
                        unreachable!("only traction edges request stress")
                    }
                }
            })
            .collect();
        let primitive_boundaries: Vec<_> = boundaries
            .iter()
            .zip(&scaled_tractions)
            .map(|(boundary, scaled)| match boundary {
                MaterialStokesBoundary::Velocity { edge, value } => StokesBoundary::Velocity {
                    edge: *edge,
                    value: *value,
                },
                MaterialStokesBoundary::TotalTraction { edge, .. } => {
                    StokesBoundary::TotalTraction {
                        edge: *edge,
                        value: scaled,
                    }
                }
            })
            .collect();
        let assembly = self
            .system
            .assemble_brinkman_mixed(
                self.params,
                cells,
                &|point| {
                    let physical = force_density(point);
                    [physical[0] / density, physical[1] / density]
                },
                &primitive_boundaries,
            )
            .map_err(MaterialFluxError::Boundary)?;
        let flow = assembly.solve().map_err(MaterialFluxError::Boundary)?;
        let cell_darcy_power_w_per_m: Vec<_> = assembly
            .cell_darcy_power_per_density(&flow.field.x[..self.system.n_u])
            .map_err(MaterialFluxError::Boundary)?
            .into_iter()
            .map(|power| density * power)
            .collect();
        let darcy_power_w_per_m: f64 = cell_darcy_power_w_per_m.iter().sum();
        let pressure_pa: Vec<_> = flow
            .pressure
            .iter()
            .map(|pressure| density * pressure)
            .collect();
        if pressure_pa.iter().any(|value| !value.is_finite())
            || cell_darcy_power_w_per_m
                .iter()
                .any(|value| !value.is_finite())
            || !darcy_power_w_per_m.is_finite()
        {
            return Err(MaterialFluxError::Invalid(
                "Brinkman physical pressure or dissipation is not representable",
            ));
        }
        Ok(MaterialBrinkmanSolution {
            flow,
            transport: self.transport,
            cells: cells.to_vec(),
            pressure_pa,
            cell_darcy_power_w_per_m,
            darcy_power_w_per_m,
        })
    }
}
