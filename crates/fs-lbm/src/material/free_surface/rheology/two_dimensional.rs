//! Source-pinned shear rheology for the original dense D2Q9 free surface.
//!
//! The incoming Guo-corrected stress closes the same bounded constitutive root
//! as sparse material flow. Only the resulting molecular tau is installed in
//! the native conservative candidate; collision, streaming and VOF conversion
//! remain owned by FreeSurface. Receipt epochs describe collision input, not
//! the subsequently streamed or converted cells.

use fs_exec::CancelGate;
use fs_material::state_point::ResolvedMaterialStatePoint;

use super::{
    CellMaterial, Error, IsothermalRheologyLaw, IsothermalRheologySample, LocalRheologyError,
    RheologySolve, capillary_scaling, plan_free_surface_from_material, plan_momentum_scaling,
    solve_local,
};
use crate::core2::{Cell, CollisionModel2, Grid, force_corrected_shear_rate};
use crate::freesurface::{ConversionMomentumReport2, FreeSurface, FreeSurfaceMomentumError2};
use crate::material::free_surface::MaterialFreeSurfaceScaling;
use crate::rheology::{TAU_CAP, TAU_FLOOR};

/// Reusable source/collision epochs for an existing dense native liquid.
///
/// This is fixed-temperature, shear-only generalized-Newtonian BGK flow.
/// Pressure-dependent sources and other collision models refuse. No material
/// fit, viscosity clipping, heat evolution or second integration path is added.
/// The caller retains the actual FreeSurface, including its geometry/contact
/// policy; its dimensions and source-scaled capillarity must remain unchanged.
#[derive(Clone)]
pub struct MaterialFreeSurfaceRheology2<'s, 'm> {
    law: &'s IsothermalRheologyLaw<'m>,
    scaling: MaterialFreeSurfaceScaling<'s>,
    cells: [usize; 2],
    solve: RheologySolve,
    samples: Vec<Option<CellMaterial>>,
    scratch: Vec<Option<CellMaterial>>,
}

impl<'s, 'm> MaterialFreeSurfaceRheology2<'s, 'm> {
    /// Admit units, source ownership and bounded per-cell constitutive work.
    /// Source buffers are allocated on first use and reused on later steps.
    pub fn try_new(
        law: &'s IsothermalRheologyLaw<'m>,
        cells: [usize; 2],
        cell_size_m: f64,
        timestep_s: f64,
        maximum_speed_m_s: f64,
        solve: RheologySolve,
    ) -> Result<Self, Error> {
        if cells.contains(&0) || cells[0].checked_mul(cells[1]).is_none() {
            return Err(Error::Invalid("invalid dense rheology grid extent"));
        }
        if law.pressure_axis().is_some() {
            return Err(Error::Invalid(
                "dense shear rheology has no absolute-pressure evolution owner",
            ));
        }
        solve
            .validate()
            .map_err(|_| Error::Invalid("invalid local rheology solve controls"))?;
        let mut scaling = plan_free_surface_from_material(
            &law.reference,
            cell_size_m,
            timestep_s,
            maximum_speed_m_s,
        )?;
        scaling.surface_tension_lattice = capillary_scaling(
            law.surface_tension_n_m,
            scaling.pressure_unit_pa,
            cell_size_m,
        )?;
        admit_tau(scaling.flow().tau)?;
        Ok(Self {
            law,
            scaling,
            cells,
            solve,
            samples: Vec::new(),
            scratch: Vec::new(),
        })
    }

    /// Shared units and actual sourced capillary coefficient. As in the sparse
    /// owner, the reference liquid carries density/viscosity and the complete
    /// original bundle; the law separately owns its active surface tension.
    pub const fn scaling(&self) -> &MaterialFreeSurfaceScaling<'s> {
        &self.scaling
    }

    pub const fn law(&self) -> &'s IsothermalRheologyLaw<'m> {
        self.law
    }

    /// Complete source state used by the most recent accepted collision.
    /// Streaming/conversion can subsequently change which cells are wet.
    pub fn material_at(&self, coordinates: [usize; 2]) -> Option<&ResolvedMaterialStatePoint> {
        self.cell_at(coordinates)
            .map(|cell| cell.material.material())
    }

    /// Physical rate and coefficient used by that same collision-input epoch.
    pub fn rheology_at(&self, coordinates: [usize; 2]) -> Option<&IsothermalRheologySample> {
        self.cell_at(coordinates).map(|cell| &cell.sample)
    }

    fn cell_at(&self, [x, y]: [usize; 2]) -> Option<&CellMaterial> {
        if x >= self.cells[0] || y >= self.cells[1] {
            return None;
        }
        self.samples.get(y * self.cells[0] + x)?.as_ref()
    }

    /// Resolve every incoming wet cell, then advance exactly one native step.
    ///
    /// Both Fluid and Interface collision lanes are sampled, even at zero
    /// interface mass. Unlocated carry has no invented material coordinate.
    /// Source, root, speed, native-flow or cancellation refusal preserves the
    /// entire supplied liquid and the previous accepted source epoch. Two
    /// adapter buffers retain allocation capacity; native candidate storage
    /// and original source-query allocation remain with their existing owners.
    pub fn step(
        &mut self,
        liquid: &mut FreeSurface,
        gate: &CancelGate,
    ) -> Result<ConversionMomentumReport2, Error> {
        poll(gate)?;
        let grid = &liquid.grid;
        let n = self.cells[0] * self.cells[1];
        if [grid.nx, grid.ny] != self.cells
            || [
                grid.flags.len(),
                grid.f.len(),
                grid.tau.len(),
                grid.fext.len(),
            ]
            .iter()
            .any(|&len| len != n)
            || grid.collision != CollisionModel2::Bgk
            || liquid.sigma.to_bits() != self.scaling.surface_tension_lattice().to_bits()
        {
            return Err(Error::Invalid(
                "dense source rheology requires matching BGK geometry, storage and capillarity",
            ));
        }
        if self.scratch.len() < n {
            self.scratch
                .try_reserve_exact(n - self.scratch.len())
                .map_err(|_| Error::Invalid("dense rheology source buffer allocation failed"))?;
            self.scratch.resize_with(n, || None);
        }
        for cell in 0..n {
            if cell % 256 == 0 {
                poll(gate)?;
            }
            self.scratch[cell] = None;
            if !matches!(grid.flags[cell], Cell::Fluid | Cell::Interface) {
                continue;
            }
            let moments = self.admit_cell(grid, cell)?;
            if !grid.tau[cell].is_finite() || grid.tau[cell] <= 0.5 {
                return Err(Error::Invalid("invalid incoming dense relaxation time"));
            }
            let force = std::array::from_fn(|d| grid.g[d].mul_add(moments.rho, grid.fext[cell][d]));
            let numerator =
                force_corrected_shear_rate(&grid.f[cell], moments.rho, moments.u, force, 1.0);
            let root = solve_local(
                numerator,
                self.scaling.timestep_s(),
                grid.tau[cell],
                self.solve,
                || poll(gate),
                |rate| {
                    let material = self.law.resolve(rate, None)?;
                    let (plan, _) = plan_momentum_scaling(
                        material.kinematic_viscosity_m2_s(),
                        self.scaling.cell_size_m(),
                        self.scaling.timestep_s(),
                        0.0,
                    )?;
                    admit_tau(plan.tau)?;
                    Ok(((material, plan.tau), plan.tau))
                },
            )
            .map_err(|error| match error {
                LocalRheologyError::Evaluation(error) => error,
                LocalRheologyError::Invalid(reason) => Error::Invalid(reason),
                LocalRheologyError::Convergence {
                    iterations,
                    relaxation_defect,
                } => Error::Convergence {
                    coordinates: [cell % grid.nx, cell / grid.nx, 0],
                    iterations,
                    relaxation_defect,
                },
            })?;
            self.scratch[cell] = Some(CellMaterial {
                material: root.value.0,
                sample: IsothermalRheologySample {
                    shear_rate_per_s: root.shear_rate_per_s,
                    collision_tau: root.value.1,
                    iterations: root.iterations,
                    relaxation_defect: root.relaxation_defect,
                },
            });
        }
        let (candidate, report) = liquid.momentum_candidate_with_preparation(gate, |grid| {
            for (cell, sample) in self.scratch.iter().enumerate() {
                if cell % 256 == 0 {
                    poll(gate)?;
                }
                if let Some(sample) = sample {
                    grid.tau[cell] = sample.sample.collision_tau;
                }
            }
            Ok::<(), Error>(())
        })?;
        for (cell, flag) in candidate.grid.flags.iter().enumerate() {
            if cell % 256 == 0 {
                poll(gate)?;
            }
            if matches!(flag, Cell::Fluid | Cell::Interface) {
                self.admit_cell(&candidate.grid, cell)?;
            }
        }
        poll(gate)?;
        *liquid = candidate;
        core::mem::swap(&mut self.samples, &mut self.scratch);
        Ok(report)
    }

    fn admit_cell(&self, grid: &Grid, cell: usize) -> Result<crate::core2::Moments, Error> {
        let rho = grid.f[cell].iter().sum::<f64>();
        if !rho.is_finite()
            || rho <= 0.0
            || grid.f[cell].iter().any(|f| !f.is_finite() || *f <= 0.0)
            || grid
                .g
                .iter()
                .chain(grid.fext[cell].iter())
                .any(|f| !f.is_finite())
        {
            return Err(Error::Invalid(
                "invalid dense wet-cell populations or held force",
            ));
        }
        let moments = grid.moments(cell);
        let speed = moments.u[0].hypot(moments.u[1]);
        if !speed.is_finite() || speed > self.scaling.flow().u_lattice {
            return Err(Error::Invalid(
                "dense rheology flow exceeds admitted scenario speed",
            ));
        }
        Ok(moments)
    }
}

fn admit_tau(tau: f64) -> Result<(), Error> {
    if !(TAU_FLOOR..=TAU_CAP).contains(&tau) {
        Err(Error::Invalid(
            "source viscosity exceeds the native relaxation range; rescale dx/dt",
        ))
    } else {
        Ok(())
    }
}

fn poll(gate: &CancelGate) -> Result<(), Error> {
    if gate.is_requested() {
        Err(FreeSurfaceMomentumError2::Cancelled.into())
    } else {
        Ok(())
    }
}
