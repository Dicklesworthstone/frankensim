//! Periodic dynamic SGS through the existing conservative transport owner.
//!
//! The shared Favre--Germano filter computes one same-stage viscosity field;
//! the original face stress, stress work, source-owned Cp heat flux and total
//! parabolic estimate consume it exactly as they consume a local closure.
//! This is the nonnegative deviatoric closure from `fs-lbm::sgs3::dynamic`,
//! with its explicit 3h test filter and caller-selected coefficient averaging.
//! Isotropic SGS energy, backscatter, dynamic heat-flux fitting and wall/shock
//! treatment are absent. Only fully periodic, isotropic Cartesian grids with
//! at least three cells per axis are admitted. No DNS accuracy is implied.

use super::*;
use fs_lbm::sgs3::dynamic::{DynamicSgsError, PeriodicDynamicSgs};

/// Dynamic stress controls plus explicitly authored gradient-diffusion heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DynamicSubgridConfig {
    pub model: DynamicSgsConfig,
    /// Positive finite turbulent Prandtl number; the heat model is not fitted.
    pub turbulent_prandtl: f64,
}

impl DynamicSubgridConfig {
    pub fn validate(self) -> Result<(), MolecularError> {
        self.model.validate()?;
        if !positive(self.turbulent_prandtl) {
            return Err(MolecularError::Invalid(
                "positive finite turbulent Prandtl number required",
            ));
        }
        Ok(())
    }
}

impl<L: TransportLaw> MolecularTransport3d<L> {
    /// Additional numerical payload for dynamic SGS, excluding the original
    /// molecular operator and fixed-size/allocator metadata: 176 bytes/cell.
    pub fn required_dynamic_sgs_workspace_bytes(flow: &Euler3d) -> Result<usize, MolecularError> {
        required_bytes(flow.states.len())
    }

    /// Fit periodic subgrid stress at every SSP stage and retry. The budget
    /// covers all added SGS Vec capacities, including density/velocity staging
    /// and the face viscosity cache; original molecular workspaces are outside
    /// this budget. Preparation and stepping never allocate inside the filter.
    pub fn with_dynamic_sgs(
        flow: &Euler3d,
        law: L,
        config: DynamicSubgridConfig,
        maximum_subgrid_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MolecularError> {
        poll(&mut cancelled)?;
        let sgs = SubgridTransport::new_dynamic(flow, config, maximum_subgrid_workspace_bytes)?;
        poll(&mut cancelled)?;
        let mut result = Self::new(flow, law, &mut cancelled)?;
        result.sgs = Some(sgs);
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Exact dynamic controls; `None` for molecular-only or local SGS flow.
    pub fn dynamic_subgrid_config(&self) -> Option<DynamicSubgridConfig> {
        self.sgs.as_ref().and_then(|sgs| match &sgs.closure {
            SubgridClosure::Dynamic(dynamic) => Some(dynamic.config),
            SubgridClosure::Local { .. } => None,
        })
    }
}

fn required_bytes(cells: usize) -> Result<usize, MolecularError> {
    let fields = cells
        .checked_mul(5 * core::mem::size_of::<f64>())
        .ok_or(MolecularError::Invalid("dynamic SGS workspace overflow"))?;
    PeriodicDynamicSgs::required_bytes(cells)?
        .checked_add(fields)
        .ok_or(MolecularError::Invalid("dynamic SGS workspace overflow"))
}

fn admit_bytes(required: usize, maximum: usize) -> Result<(), MolecularError> {
    if required > maximum {
        Err(DynamicSgsError::WorkspaceBytes { required, maximum }.into())
    } else {
        Ok(())
    }
}

pub(super) struct DynamicTransport {
    pub(super) config: DynamicSubgridConfig,
    filter: PeriodicDynamicSgs,
    density: Vec<f64>,
    velocity: Vec<[f64; 3]>,
}

impl SubgridTransport {
    pub(in crate::compressible::spatial::molecular) fn new_dynamic(
        flow: &Euler3d,
        config: DynamicSubgridConfig,
        maximum_bytes: usize,
    ) -> Result<Self, MolecularError> {
        config.validate()?;
        if (0..3).any(|axis| flow.boundaries.axis(axis) != AxisBoundaries::Periodic)
            || flow.grid.shape.iter().any(|&n| n < 3)
            || flow.grid.spacing[1] != flow.grid.spacing[0]
            || flow.grid.spacing[2] != flow.grid.spacing[0]
        {
            return Err(MolecularError::Invalid(
                "dynamic SGS requires a periodic isotropic grid with at least three cells per axis",
            ));
        }
        let n = flow.states.len();
        admit_bytes(required_bytes(n)?, maximum_bytes)?;
        let dynamic = DynamicTransport {
            config,
            filter: PeriodicDynamicSgs::new(flow.grid.shape, flow.grid.spacing[0], maximum_bytes)?,
            density: allocated(n, 0.0)?,
            velocity: allocated(n, [0.0; 3])?,
        };
        let dynamic_viscosity = allocated(n, 0.0)?;
        let bytes = dynamic.filter.workspace_bytes()
            + dynamic.density.capacity() * core::mem::size_of::<f64>()
            + dynamic.velocity.capacity() * core::mem::size_of::<[f64; 3]>()
            + dynamic_viscosity.capacity() * core::mem::size_of::<f64>();
        admit_bytes(bytes, maximum_bytes)?;
        Ok(Self {
            closure: SubgridClosure::Dynamic(dynamic),
            dynamic_viscosity,
        })
    }
}

impl DynamicTransport {
    pub(super) fn prepare(
        &mut self,
        primitive: &[Primitive],
        dynamic_viscosity: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), MolecularError> {
        for (cell, p) in primitive.iter().enumerate() {
            if cell % 64 == 0 {
                poll(cancelled)?;
            }
            self.density[cell] = p.density;
            self.velocity[cell] = p.velocity;
        }
        self.filter.prepare(
            self.config.model,
            &self.density,
            &self.velocity,
            &mut *cancelled,
        )?;
        for (cell, ((p, &nu), mu)) in primitive
            .iter()
            .zip(self.filter.viscosities())
            .zip(dynamic_viscosity)
            .enumerate()
        {
            if cell % 64 == 0 {
                poll(cancelled)?;
            }
            *mu = p.density * nu;
            if !mu.is_finite() {
                return Err(MolecularError::Invalid(
                    "unrepresentable dynamic subgrid viscosity",
                ));
            }
        }
        poll(cancelled)
    }
}
