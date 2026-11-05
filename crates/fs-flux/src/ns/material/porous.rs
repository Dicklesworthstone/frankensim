//! Molecular source viscosity and physical units for native Brinkman momentum.

mod forchheimer;
pub use forchheimer::MaterialForchheimerSolution;

use super::{MaterialBoundaryLoad, MaterialFluxError, MaterialFluxSystem, MaterialStokesBoundary};
use crate::ns::FluxSolution;
use crate::ns::boundary::{StokesBoundary, StokesPressureReference, StokesSolution};
use crate::ns::porous::{BrinkmanCell, BrinkmanDragExchange};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};

/// Physical porous-resistance force and work per metre of out-of-plane depth.
/// The accessor identifies the Darcy or Forchheimer law being accounted for.
///
/// Positive power acts ON the named subsystem. The two signed powers sum
/// to minus the nonnegative relative dissipation. No heat or solid dynamics
/// are updated implicitly; the coupled owner decides where these exchanges go.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MaterialBrinkmanDragExchange {
    /// Reaction exerted by the fluid on the skeleton [N/m].
    pub skeleton_force_n_per_m: [f64; 2],
    /// Signed work of the identified resistance on the fluid [W/m].
    pub fluid_power_w_per_m: f64,
    /// Signed reaction work on the prescribed skeleton [W/m].
    pub skeleton_power_w_per_m: f64,
    /// Nonnegative relative loss [W/m]: `integral(mu w.K^-1.w)` for Darcy
    /// or `integral(rho beta |w|^3)` for Forchheimer, where `w = u-u_s`.
    pub dissipation_w_per_m: f64,
}

impl MaterialBrinkmanDragExchange {
    fn from_kinematic(
        value: BrinkmanDragExchange,
        density: f64,
    ) -> Result<Self, MaterialFluxError> {
        let physical = Self {
            skeleton_force_n_per_m: value.skeleton_force.map(|force| density * force),
            fluid_power_w_per_m: density * value.fluid_power,
            skeleton_power_w_per_m: density * value.skeleton_power,
            dissipation_w_per_m: density * value.dissipation,
        };
        physical.admit()?;
        Ok(physical)
    }

    fn add(&mut self, value: Self) {
        for a in 0..2 {
            self.skeleton_force_n_per_m[a] += value.skeleton_force_n_per_m[a];
        }
        self.fluid_power_w_per_m += value.fluid_power_w_per_m;
        self.skeleton_power_w_per_m += value.skeleton_power_w_per_m;
        self.dissipation_w_per_m += value.dissipation_w_per_m;
    }

    fn admit(self) -> Result<(), MaterialFluxError> {
        if !self
            .skeleton_force_n_per_m
            .iter()
            .all(|force| force.is_finite())
            || !self.fluid_power_w_per_m.is_finite()
            || !self.skeleton_power_w_per_m.is_finite()
            || !self.dissipation_w_per_m.is_finite()
            || self.dissipation_w_per_m < 0.0
        {
            Err(MaterialFluxError::Invalid(
                "Brinkman physical resistance force or work is not representable",
            ))
        } else {
            Ok(())
        }
    }
}

/// Accepted native flow, original transport snapshot and actual pore dissipation.
pub struct MaterialBrinkmanSolution {
    flow: StokesSolution,
    transport: FluidTransportState,
    cells: Vec<BrinkmanCell>,
    pressure_pa: Vec<f64>,
    cell_darcy_power_w_per_m: Vec<f64>,
    darcy_power_w_per_m: f64,
    cell_darcy_exchange: Vec<MaterialBrinkmanDragExchange>,
    darcy_exchange: MaterialBrinkmanDragExchange,
}

impl MaterialBrinkmanSolution {
    /// Physical discrete boundary reactions, retaining pore drag in the same
    /// momentum equation. Mean-zero components need explicit pressure offsets
    /// [Pa]; natural components retain their traction datum with None.
    ///
    /// # Errors
    /// Refuses invalid selections/references or unrepresentable physical loads.
    pub fn boundary_loads(
        &self,
        edges: &[usize],
        origin_m: [f64; 2],
        pressure_offsets_pa: &[Option<f64>],
    ) -> Result<MaterialBoundaryLoad, MaterialFluxError> {
        super::mixed::physical_boundary_load(
            &self.flow,
            self.transport.density_kg_m3(),
            edges,
            origin_m,
            pressure_offsets_pa,
        )
    }
    /// Native BDM1/P0 result, including independent equation residuals.
    #[must_use]
    pub const fn discrete(&self) -> &FluxSolution {
        &self.flow.field
    }
    /// Physical macroscopic-velocity moments [m/s]. Prescribed skeleton motion
    /// uses the one-domain/phi=1 convention documented by [`BrinkmanCell`].
    #[must_use]
    pub fn velocity_dofs(&self) -> &[f64] {
        &self.flow.field.x[..self.flow.field.n_u]
    }
    /// Pressure [Pa] in the declared boundary datum, never an inferred EOS state.
    #[must_use]
    pub fn pressure_pa(&self) -> &[f64] {
        &self.pressure_pa
    }
    /// The traction datum, all-velocity mean-zero gauge, or componentwise mix.
    #[must_use]
    pub const fn pressure_reference(&self) -> StokesPressureReference {
        self.flow.pressure_reference
    }
    /// Local datum in original cell order; independent closed components have
    /// no known relative pressure offset.
    #[must_use]
    pub fn pressure_reference_at(&self, cell: usize) -> Option<StokesPressureReference> {
        self.flow.pressure_reference_per_cell.get(cell).copied()
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
    /// `integral(mu (u-u_s).K^-1.(u-u_s))` per triangle and unit depth [W/m].
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

    /// Original cell quadratures converted with the SAME material density.
    /// Includes skeleton reaction and both signed powers, not only the loss.
    #[must_use]
    pub fn cell_darcy_exchange(&self) -> &[MaterialBrinkmanDragExchange] {
        &self.cell_darcy_exchange
    }

    /// Sum of the original cell force/work exchanges, in deterministic order.
    #[must_use]
    pub const fn darcy_exchange(&self) -> &MaterialBrinkmanDragExchange {
        &self.darcy_exchange
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
    /// This is a steady incompressible closure on frozen geometry. Authored
    /// skeleton velocities use the one-domain/high-porosity approximation
    /// described by [`BrinkmanCell::with_skeleton_velocity`]. Their reaction
    /// work remains separate from relative dissipation: neither is silently
    /// deposited into heat. No solid dynamics, moving interface, finite-porosity
    /// storage correction, poromechanics, capillary or gas-slip law is added.
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
        self.solve_brinkman_mixed_with_history(cells, force_density, boundaries, None)
    }

    /// Advance native Stokes–Brinkman by one fully implicit backward-Euler step.
    ///
    /// Force density [N/m³], total traction [Pa], velocities [m/s] and cell laws
    /// prescribe the endpoint; the timestep is seconds. The immutable bound
    /// transport snapshot supplies molecular viscosity and the SAME density to
    /// inertia, drag, stress and physical reaction recovery. Changing density
    /// between steps requires an explicitly conservative state transfer; a new
    /// material binding alone does not perform that transfer.
    ///
    /// Previous velocity contains exactly the native BDM1 velocity moments and
    /// is never modified. The original pressure constraints, boundary reactions
    /// and checked linear solve are retained. Returned Darcy powers are endpoint
    /// rates, not accumulated heat; numerical time-discretization loss is separate.
    /// Geometry stays fixed and no convection, finite-porosity storage, moving
    /// interface or solid dynamics is inferred from a moving skeleton velocity.
    ///
    /// # Errors
    /// Preserves [`Self::solve_brinkman_mixed`] material and native refusals,
    /// also rejecting invalid history or timestep before sampling callbacks.
    pub fn bdf1_brinkman_mixed_step(
        &self,
        cells: &[BrinkmanCell],
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialStokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
    ) -> Result<MaterialBrinkmanSolution, MaterialFluxError> {
        self.solve_brinkman_mixed_with_history(
            cells,
            force_density,
            boundaries,
            Some((previous_velocity, timestep_s)),
        )
    }

    fn solve_brinkman_mixed_with_history(
        &self,
        cells: &[BrinkmanCell],
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialStokesBoundary<'_>],
        history: Option<(&[f64], f64)>,
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
                    MaterialStokesBoundary::Velocity { .. }
                    | MaterialStokesBoundary::NavierSlip { .. } => {
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
                MaterialStokesBoundary::NavierSlip {
                    edge,
                    value,
                    tangential_resistance_pa_s_per_m,
                } => StokesBoundary::NavierSlip {
                    edge: *edge,
                    value: *value,
                    tangential_resistance_m_s: super::mixed::kinematic_slip_resistance(
                        *tangential_resistance_pa_s_per_m,
                        density,
                    ),
                },
                MaterialStokesBoundary::TotalTraction { edge, .. } => {
                    StokesBoundary::TotalTraction {
                        edge: *edge,
                        value: scaled,
                    }
                }
            })
            .collect();
        let force = |point| {
            let physical = force_density(point);
            [physical[0] / density, physical[1] / density]
        };
        let assembly = match history {
            Some((previous_velocity, timestep_s)) => self.system.assemble_brinkman_mixed_bdf1(
                self.params,
                cells,
                &force,
                &primitive_boundaries,
                previous_velocity,
                timestep_s,
            ),
            None => self.system.assemble_brinkman_mixed(
                self.params,
                cells,
                &force,
                &primitive_boundaries,
            ),
        }
        .map_err(MaterialFluxError::Boundary)?;
        let flow = assembly.solve().map_err(MaterialFluxError::Boundary)?;
        let exchanges = assembly
            .cell_darcy_exchange_per_density(&flow.field.x[..self.system.n_u])
            .map_err(MaterialFluxError::Boundary)?;
        material_brinkman_solution(flow, self.transport, cells, exchanges)
    }
}

// One physical recovery path for all native Brinkman solves. Child closures may
// retain additional resistance channels without replacing this source binding,
// pressure conversion or the original discrete boundary reaction owner.
fn material_brinkman_solution(
    flow: StokesSolution,
    transport: FluidTransportState,
    cells: &[BrinkmanCell],
    exchanges: Vec<BrinkmanDragExchange>,
) -> Result<MaterialBrinkmanSolution, MaterialFluxError> {
    let density = transport.density_kg_m3();
    let cell_darcy_exchange: Vec<_> = exchanges
        .into_iter()
        .map(|value| MaterialBrinkmanDragExchange::from_kinematic(value, density))
        .collect::<Result<_, _>>()?;
    let mut darcy_exchange = MaterialBrinkmanDragExchange::default();
    let cell_darcy_power_w_per_m: Vec<_> = cell_darcy_exchange
        .iter()
        .map(|exchange| {
            darcy_exchange.add(*exchange);
            exchange.dissipation_w_per_m
        })
        .collect();
    darcy_exchange.admit()?;
    let darcy_power_w_per_m = darcy_exchange.dissipation_w_per_m;
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
        transport,
        cells: cells.to_vec(),
        pressure_pa,
        cell_darcy_power_w_per_m,
        darcy_power_w_per_m,
        cell_darcy_exchange,
        darcy_exchange,
    })
}
