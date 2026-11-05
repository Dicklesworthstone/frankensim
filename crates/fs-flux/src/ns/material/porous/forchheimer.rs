//! Original molecular transport and SI work for nonlinear porous momentum.

use super::{MaterialBrinkmanDragExchange, MaterialBrinkmanSolution, material_brinkman_solution};
use crate::material::{MaterialFluxError, MaterialFluxSystem, MaterialStokesBoundary};
use crate::ns::boundary::StokesBoundary;
use crate::ns::porous::{ForchheimerCell, ForchheimerControls};
use fs_material::fluid_transport::FluidTransport;

/// Accepted nonlinear flow with source-bound physical Darcy and inertial work.
/// The original Brinkman result owns pressure, boundary reactions and the Darcy
/// account; the added Forchheimer account is separate and counted exactly once.
pub struct MaterialForchheimerSolution {
    flow: MaterialBrinkmanSolution,
    cells: Vec<ForchheimerCell>,
    cell_forchheimer_exchange: Vec<MaterialBrinkmanDragExchange>,
    forchheimer_exchange: MaterialBrinkmanDragExchange,
    iterations: usize,
    line_search_reductions: usize,
    relative_update: f64,
}

impl MaterialForchheimerSolution {
    /// Complete accepted NONLINEAR flow, original transport snapshot, physical
    /// pressure/boundary reactions, and the separately identified Darcy account.
    #[must_use]
    pub const fn flow(&self) -> &MaterialBrinkmanSolution {
        &self.flow
    }

    /// Original anisotropic permeability, skeleton motion and authored beta.
    #[must_use]
    pub fn cells(&self) -> &[ForchheimerCell] {
        &self.cells
    }

    /// Physical form-drag reactions and signed work per triangle, per unit
    /// out-of-plane depth. Same source density as momentum and pressure.
    #[must_use]
    pub fn cell_forchheimer_exchange(&self) -> &[MaterialBrinkmanDragExchange] {
        &self.cell_forchheimer_exchange
    }

    /// Deterministic cell sum of form-drag force and signed work, in SI units.
    #[must_use]
    pub const fn forchheimer_exchange(&self) -> &MaterialBrinkmanDragExchange {
        &self.forchheimer_exchange
    }

    /// Positive endpoint `integral(rho beta |u-u_s|^3)` [W/m]. This rate is not
    /// a heat increment, skeleton work or backward-Euler numerical dissipation.
    #[must_use]
    pub const fn forchheimer_power_w_per_m(&self) -> f64 {
        self.forchheimer_exchange.dissipation_w_per_m
    }

    /// Newton solves performed; zero also identifies the original zero-beta path.
    #[must_use]
    pub const fn iterations(&self) -> usize {
        self.iterations
    }

    /// Total rejected line-search trials.
    #[must_use]
    pub const fn line_search_reductions(&self) -> usize {
        self.line_search_reductions
    }

    /// Full, unrelaxed final relative velocity correction.
    #[must_use]
    pub const fn relative_update(&self) -> f64 {
        self.relative_update
    }
}

impl<F: FluidTransport + ?Sized> MaterialFluxSystem<'_, '_, '_, F> {
    /// Solve nonlinear steady Brinkman–Forchheimer momentum with the original
    /// once-bound density and molecular viscosity. Force is N/m³, velocity m/s,
    /// total traction Pa and authored isotropic beta 1/m. Physical form drag is
    /// `rho beta |u-u_s| (u-u_s)`, without a factor of 1/2.
    ///
    /// Molecular viscosity sets the original anisotropic Darcy force; the cell's
    /// effective-viscosity ratio only changes shear stress. The original source
    /// is not re-queried or replaced by a fitted effective permeability. A beta
    /// input alone carries no measured calibration authority or physical support
    /// envelope. A directional gas calibration cannot be inferred to apply here.
    ///
    /// The primitive owner's saturated, fixed-domain phi=1 convention, Newton
    /// controls, pure-callback requirement and physical regime limits apply.
    /// All boundary reactions are evaluated from the accepted nonlinear state.
    /// Darcy and Forchheimer force/work remain separate, with no implicit heat.
    ///
    /// # Errors
    /// Preserves native material, physical representability, boundary and
    /// nonlinear convergence refusals; no candidate state is returned on failure.
    pub fn solve_brinkman_forchheimer_mixed(
        &self,
        cells: &[ForchheimerCell],
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialStokesBoundary<'_>],
        initial_velocity: &[f64],
        controls: ForchheimerControls,
    ) -> Result<MaterialForchheimerSolution, MaterialFluxError> {
        self.solve_material_forchheimer(
            cells,
            force_density,
            boundaries,
            initial_velocity,
            controls,
            None,
        )
    }

    /// Implicit native BDF1 mass/history and endpoint nonlinear pore resistance.
    /// History and the once-bound reference density are immutable through the
    /// step. The caller owns the clock and commits the returned velocity only
    /// after success. Changing density requires a conservative state transfer.
    ///
    /// This is unsteady incompressible Stokes on frozen geometry, with the same
    /// phi=1 convention. It adds neither convective/finite-porosity inertia nor
    /// solid dynamics. Reported physical powers exclude BDF1 numerical loss.
    ///
    /// # Errors
    /// Also preserves native timestep/history admission before load callbacks.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf1_brinkman_forchheimer_mixed_step(
        &self,
        cells: &[ForchheimerCell],
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialStokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
        controls: ForchheimerControls,
    ) -> Result<MaterialForchheimerSolution, MaterialFluxError> {
        self.solve_material_forchheimer(
            cells,
            force_density,
            boundaries,
            previous_velocity,
            controls,
            Some(timestep_s),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn solve_material_forchheimer(
        &self,
        cells: &[ForchheimerCell],
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialStokesBoundary<'_>],
        initial_velocity: &[f64],
        controls: ForchheimerControls,
        timestep_s: Option<f64>,
    ) -> Result<MaterialForchheimerSolution, MaterialFluxError> {
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
                    tangential_resistance_m_s: crate::material::mixed::kinematic_slip_resistance(
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
        let solution = if let Some(dt) = timestep_s {
            self.system.bdf1_brinkman_forchheimer_mixed_step(
                self.params,
                cells,
                &force,
                &primitive_boundaries,
                initial_velocity,
                dt,
                controls,
            )
        } else {
            self.system.solve_brinkman_forchheimer_mixed(
                self.params,
                cells,
                &force,
                &primitive_boundaries,
                initial_velocity,
                controls,
            )
        }
        .map_err(MaterialFluxError::Boundary)?;
        let cell_forchheimer_exchange: Vec<_> = solution
            .cell_forchheimer_exchange
            .into_iter()
            .map(|exchange| MaterialBrinkmanDragExchange::from_kinematic(exchange, density))
            .collect::<Result<_, _>>()?;
        let mut forchheimer_exchange = MaterialBrinkmanDragExchange::default();
        for &exchange in &cell_forchheimer_exchange {
            forchheimer_exchange.add(exchange);
        }
        forchheimer_exchange.admit()?;
        let brinkman: Vec<_> = cells.iter().map(|cell| *cell.brinkman()).collect();
        let flow = material_brinkman_solution(
            solution.flow,
            self.transport,
            &brinkman,
            solution.cell_darcy_exchange,
        )?;
        Ok(MaterialForchheimerSolution {
            flow,
            cells: cells.to_vec(),
            cell_forchheimer_exchange,
            forchheimer_exchange,
            iterations: solution.iterations,
            line_search_reductions: solution.line_search_reductions,
            relative_update: solution.relative_update,
        })
    }
}
