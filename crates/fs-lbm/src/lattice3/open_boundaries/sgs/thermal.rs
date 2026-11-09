//! The same local turbulent viscosity in native open momentum and heat.

use super::{
    EquilibriumLattice3, GridError3, OpenBoundaryGrid3, OpenSgsWorkspace3, SgsConfig, SgsLattice3,
    SgsReport, SgsStepError,
};
use crate::lattice3::open_boundaries::thermal::{OpenThermalError3, OpenThermalStep3};
use crate::sgs3::SgsError;
use crate::thermal::open::ThermalOpenField3;
use fs_exec::Cx;

/// Explicit local closure and gradient-diffusion heat policy. The native
/// flow/thermal owners retain their independent physical-velocity admissions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OpenThermalSgsConfig3 {
    pub closure: SgsConfig,
    /// Positive finite nu_t/alpha_t, with no inferred/default value.
    pub turbulent_prandtl: f64,
    /// Nonnegative finite alpha_t/alpha_molecular admission, never a clamp.
    pub maximum_diffusivity_to_molecular: f64,
}

#[derive(Clone, Debug)]
pub enum OpenThermalSgsError3<E> {
    Sgs(SgsStepError<E>),
    Thermal(OpenThermalError3<E>),
    Invalid(&'static str),
    DiffusivityRatio {
        cell: usize,
        ratio: f64,
        maximum: f64,
    },
}
impl<E> From<SgsStepError<E>> for OpenThermalSgsError3<E> {
    fn from(error: SgsStepError<E>) -> Self {
        Self::Sgs(error)
    }
}
impl<E> From<SgsError> for OpenThermalSgsError3<E> {
    fn from(error: SgsError) -> Self {
        Self::Sgs(SgsStepError::Closure(error))
    }
}
impl<E> From<OpenThermalError3<E>> for OpenThermalSgsError3<E> {
    fn from(error: OpenThermalError3<E>) -> Self {
        Self::Thermal(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for OpenThermalSgsError3<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native open thermal SGS: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for OpenThermalSgsError3<E> {}

/// Reused momentum/macroscopic/eddy storage and one total thermal-tau field.
/// This workspace holds no independent populations, wall temperatures or time.
pub struct OpenThermalSgsWorkspace3<C> {
    flow: OpenSgsWorkspace3<C>,
    thermal: Vec<f64>,
}
impl<C> OpenThermalSgsWorkspace3<C> {
    pub fn momentum_collisions(&self) -> &[C] {
        self.flow.collisions()
    }
    pub fn thermal_relaxation(&self) -> &[f64] {
        &self.thermal
    }
    pub fn eddy_viscosities(&self) -> &[f64] {
        self.flow.eddy_viscosities()
    }
    pub fn workspace_bytes(&self) -> usize {
        self.flow.workspace_bytes() + self.thermal.capacity() * core::mem::size_of::<f64>()
    }
}

#[derive(Clone, Debug)]
pub struct OpenThermalSgsStep3 {
    /// Actual native momentum and original reference-capacity heat accounts.
    pub transport: OpenThermalStep3,
    /// Incoming resolved/modeled shear transfer; this is not deposited heat.
    pub sgs: SgsReport,
    pub maximum_eddy_diffusivity_lu: f64,
    pub maximum_diffusivity_to_molecular: f64,
}

impl<L: EquilibriumLattice3<Q> + SgsLattice3<Q>, const Q: usize> OpenBoundaryGrid3<L, Q> {
    pub fn required_thermal_sgs_workspace_bytes(&self) -> Result<usize, SgsStepError<L::Error>> {
        let heat = self
            .grid
            .f
            .len()
            .checked_mul(core::mem::size_of::<f64>())
            .ok_or(GridError3::Input("open thermal SGS workspace overflow"))?;
        self.required_sgs_workspace_bytes()?
            .checked_add(heat)
            .ok_or_else(|| GridError3::Input("open thermal SGS workspace overflow").into())
    }

    /// Reserve six f64 values plus one native collision value per grid cell.
    /// Original native and finite-wall storage retain their own budgets.
    pub fn thermal_sgs_workspace(
        &self,
        maximum_workspace_bytes: usize,
    ) -> Result<OpenThermalSgsWorkspace3<L::Collision>, SgsStepError<L::Error>> {
        if self.required_thermal_sgs_workspace_bytes()? > maximum_workspace_bytes {
            return Err(
                GridError3::Input("open thermal SGS workspace byte budget exceeded").into(),
            );
        }
        let flow = self.sgs_workspace(maximum_workspace_bytes)?;
        let mut thermal = Vec::new();
        thermal
            .try_reserve_exact(self.grid.f.len())
            .map_err(|_| GridError3::Allocation)?;
        thermal.resize(self.grid.f.len(), 1.0);
        let workspace = OpenThermalSgsWorkspace3 { flow, thermal };
        if workspace.workspace_bytes() > maximum_workspace_bytes {
            return Err(GridError3::Input("open thermal SGS capacity budget exceeded").into());
        }
        Ok(workspace)
    }

    /// Advance the original open-flow/heat/finite-wall transaction with one
    /// common incoming SGS viscosity and explicit turbulent Prandtl number.
    /// Fresh per-cell molecular collisions are optional; molecular thermal
    /// taus are an explicit x-fast field covering every cell. Never feed a
    /// previously augmented field back as molecular transport.
    ///
    /// The native reservoir gradient supplies nu_t once. D3Q7 has cs²=1/4,
    /// so alpha_t=nu_t/Pr_t adds 4*alpha_t to the original thermal tau. Every
    /// finite, representable positive increment and molecular-ratio bound is
    /// admitted before either physical field changes. Zero nu_t preserves
    /// collision and thermal tau bits. Coefficients are frozen over this step.
    ///
    /// Reservoir inflow/outflow reversal, scalar-population heat exchange,
    /// source half-kicks, advective history and finite conducting-wall graph
    /// remain with the original thermal owner. Its real fluid/solid heat
    /// balance and final cancellation gate precede joint publication. Scratch
    /// may change on refusal; physical populations, wall state, receipts and
    /// clocks do not. No adapter-owned allocation occurs during preparation.
    ///
    /// This is fixed reference-capacity, low-Mach temperature transport.
    /// Modeled mechanical transfer remains diagnostic. No heating source,
    /// wall function, pressure-work/enthalpy model, dynamic exterior filter,
    /// nonreflecting outlet or turbulent-jet qualification is supplied.
    #[allow(clippy::too_many_arguments)]
    pub fn step_serial_with_thermal_sgs(
        &mut self,
        thermal: &mut ThermalOpenField3,
        cx: &Cx<'_>,
        config: OpenThermalSgsConfig3,
        molecular_collisions: Option<&[L::Collision]>,
        heat_per_step: &[f64],
        molecular_thermal_relaxation: &[f64],
        workspace: &mut OpenThermalSgsWorkspace3<L::Collision>,
        origin: [f64; 3],
    ) -> Result<OpenThermalSgsStep3, OpenThermalSgsError3<L::Error>> {
        let (sgs, maximum_eddy_diffusivity_lu, maximum_diffusivity_to_molecular) = self
            .prepare_open_thermal_sgs(
                config,
                molecular_collisions,
                molecular_thermal_relaxation,
                workspace,
                None,
                || cx.checkpoint().map_err(|_| GridError3::Cancelled),
            )?;
        let transport = self.step_serial_with_thermal(
            thermal,
            cx,
            Some(workspace.flow.collisions()),
            heat_per_step,
            Some(&workspace.thermal),
            origin,
        )?;
        Ok(OpenThermalSgsStep3 {
            transport,
            sgs,
            maximum_eddy_diffusivity_lu,
            maximum_diffusivity_to_molecular,
        })
    }

    // One coefficient owner for passive thermal and buoyant/spatial-force
    // transactions. The explicit force slice belongs to the accepted incoming
    // populations; endpoint forces are owned by the subsequent native step.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::lattice3::open_boundaries) fn prepare_open_thermal_sgs(
        &self,
        config: OpenThermalSgsConfig3,
        molecular_collisions: Option<&[L::Collision]>,
        molecular_thermal_relaxation: &[f64],
        workspace: &mut OpenThermalSgsWorkspace3<L::Collision>,
        forces: Option<&[[f64; 3]]>,
        mut checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<(SgsReport, f64, f64), OpenThermalSgsError3<L::Error>> {
        checkpoint().map_err(SgsStepError::Grid)?;
        config.closure.validate()?;
        if !config.turbulent_prandtl.is_finite()
            || config.turbulent_prandtl <= 0.0
            || !config.maximum_diffusivity_to_molecular.is_finite()
            || config.maximum_diffusivity_to_molecular < 0.0
        {
            return Err(OpenThermalSgsError3::Invalid(
                "invalid turbulent heat controls",
            ));
        }
        let cells = self.grid.f.len();
        if molecular_thermal_relaxation.len() != cells || workspace.thermal.len() != cells {
            return Err(OpenThermalSgsError3::Invalid(
                "molecular/total thermal tau shape mismatch",
            ));
        }
        let sgs = self.prepare_open_sgs(
            config.closure,
            molecular_collisions,
            &mut workspace.flow,
            forces,
            &mut checkpoint,
        )?;
        let mut maximum_eddy_diffusivity_lu = 0.0_f64;
        let mut maximum_diffusivity_to_molecular = 0.0_f64;
        for (cell, &tau) in molecular_thermal_relaxation.iter().enumerate() {
            if cell % 64 == 0 {
                checkpoint().map_err(SgsStepError::Grid)?;
            }
            if !tau.is_finite() || tau <= 0.5 {
                return Err(OpenThermalSgsError3::Invalid(
                    "invalid molecular thermal relaxation",
                ));
            }
            workspace.thermal[cell] = tau;
            let eddy = workspace.flow.eddy[cell];
            if self.grid.solid[cell] || eddy == 0.0 {
                continue;
            }
            let alpha = eddy / config.turbulent_prandtl;
            // Form the ratio against tau-1/2 directly: dividing tiny molecular
            // diffusivity by four first could spuriously underflow that source.
            let increment = 4.0 * alpha;
            let ratio = increment / (tau - 0.5);
            if !alpha.is_finite()
                || alpha <= 0.0
                || !increment.is_finite()
                || increment <= 0.0
                || !(tau + increment).is_finite()
            {
                return Err(OpenThermalSgsError3::Invalid(
                    "unrepresentable turbulent thermal transport",
                ));
            }
            if !ratio.is_finite() || ratio > config.maximum_diffusivity_to_molecular {
                return Err(OpenThermalSgsError3::DiffusivityRatio {
                    cell,
                    ratio,
                    maximum: config.maximum_diffusivity_to_molecular,
                });
            }
            workspace.thermal[cell] = tau + increment;
            maximum_eddy_diffusivity_lu = maximum_eddy_diffusivity_lu.max(alpha);
            maximum_diffusivity_to_molecular = maximum_diffusivity_to_molecular.max(ratio);
        }
        checkpoint().map_err(SgsStepError::Grid)?;
        Ok((
            sgs,
            maximum_eddy_diffusivity_lu,
            maximum_diffusivity_to_molecular,
        ))
    }
}
