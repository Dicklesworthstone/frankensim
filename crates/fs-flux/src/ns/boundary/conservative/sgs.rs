//! Nonlinear SGS stress in the mass-consistent variable-density momentum owner.
//!
//! Density, carried inventories and the conservative advecting mass flux are
//! prepared ONCE from the authored old states. Only constitutive stress is
//! iterated: mu_eff = mu_molecular(gamma, transported cell) + rho_new * nu_sgs.
//! Pressure, inertia, SIP stress and boundary reactions retain their native
//! owner and common numerical density normalization.
//!
//! This is the same Estimated planar screen as crate::ns::sgs. A two-dimensional
//! field does not resolve a three-dimensional turbulent cascade. There is no
//! turbulent mass/species diffusion, EOS, wall model, SGS energy equation or
//! automatic thermal deposition. Isotropic modeled pressure is absorbed in the
//! incompressible mechanical pressure; only deviatoric SGS stress is added.

use super::ssprk2::DensityMomentumSspRk2Controls;
use super::{
    DensityMomentumControls, DensityMomentumError, DensityMomentumStep, PreparedDensityMomentum,
    allocate, observe_step, poll,
};
use crate::material::{MaterialFluxNumerics, MaterialNavierStokesBoundary};
use crate::ns::FluxSystem;
use crate::ns::boundary::{StokesBoundaryError, StokesRheologyControls};
use crate::ns::sgs::{PlanarSgsError, SgsModel};
use crate::ns::transport::density::DensityTransportCell;
use crate::ns::transport::density::ssprk2::DensitySspRk2Report;
use fs_lbm::sgs3::eddy_viscosity;
use std::cell::RefCell;

/// Explicit isotropic filter widths in the original triangle order.
#[derive(Clone, Copy, Debug)]
pub struct DensityMomentumSgsConfig<'a> {
    pub model: SgsModel,
    pub filter_widths_m: &'a [f64],
}

/// Original transport/flow refusal or a shared SGS-kernel refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum DensityMomentumSgsError {
    Momentum(DensityMomentumError),
    Subgrid(PlanarSgsError),
}
impl core::fmt::Display for DensityMomentumSgsError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "variable-density SGS momentum: {self:?}")
    }
}
impl std::error::Error for DensityMomentumSgsError {}
impl From<DensityMomentumError> for DensityMomentumSgsError {
    fn from(error: DensityMomentumError) -> Self {
        Self::Momentum(error)
    }
}
impl From<StokesBoundaryError> for DensityMomentumSgsError {
    fn from(error: StokesBoundaryError) -> Self {
        Self::Momentum(error.into())
    }
}
impl From<PlanarSgsError> for DensityMomentumSgsError {
    fn from(error: PlanarSgsError) -> Self {
        Self::Subgrid(error)
    }
}
type Result<T> = core::result::Result<T, DensityMomentumSgsError>;

/// Retained original cell law and shared SGS response at the accepted velocity.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DensityMomentumSgsCell {
    /// Native symmetric strain invariant sqrt(2 D:D) [1/s].
    pub strain_rate_per_s: f64,
    /// Molecular source response at this cell's centroid and actual strain [Pa s].
    pub molecular_dynamic_viscosity_pa_s: f64,
    /// Shared model's kinematic eddy viscosity [m²/s].
    pub eddy_viscosity_m2_s: f64,
    /// rho_new * nu_sgs, with the density accepted by this transport step [Pa s].
    pub modeled_dynamic_viscosity_pa_s: f64,
    /// The molecular plus modeled coefficient used by native momentum [Pa s].
    pub effective_dynamic_viscosity_pa_s: f64,
    /// Molecular-only mu * gamma² [W/m³], excluding modeled/numerical transfer.
    pub molecular_heating_w_m3: f64,
    /// Nonnegative resolved-to-modeled rho_new * nu_sgs * gamma² [W/m³].
    pub modeled_transfer_w_m3: f64,
}

/// One jointly accepted density/inventory/velocity endpoint and its SGS stress.
pub struct DensityMomentumSgsStep<const N: usize> {
    /// Original conservative momentum, mechanical pressure and numerical losses.
    pub step: DensityMomentumStep<N>,
    /// Exact endpoint source samples; reporting does not query the source again.
    pub cells: Vec<DensityMomentumSgsCell>,
    /// Volume integral of molecular heating, using the authored depth [W].
    pub molecular_heating_w: f64,
    /// Volume integral of modeled transfer, using the authored depth [W].
    pub modeled_transfer_w: f64,
    pub iterations: usize,
    pub relative_update: f64,
    /// Maximum of the original and row-scaled nonlinear endpoint residuals.
    pub nonlinear_residual: f64,
}

/// One accepted SSPRK2 density/inventory endpoint and native nonlinear SGS flow.
pub struct DensityMomentumSspRk2SgsStep<const N: usize> {
    /// Original BDF1 momentum, physical pressure and accepted constitutive
    /// diagnostics. Momentum retains first-order time integration.
    pub momentum: DensityMomentumSgsStep<N>,
    /// Actual RK work and weighted physical inventory report. `transport` is
    /// identical to `momentum.step.transport`.
    pub density_stages: DensitySspRk2Report<N>,
}

impl<'m> FluxSystem<'m> {
    /// Solve native conservative BDF1 momentum with nonlinear molecular and SGS
    /// stress at the same accepted velocity.
    ///
    /// The molecular law is a CELLWISE constitutive source: it is queried once
    /// per original triangle and nonlinear assembly, at its physical centroid,
    /// with the actual candidate shear rate and the transported density/specific
    /// inventories. Its positive finite result [Pa s] is used unchanged at all
    /// cell and face quadrature sites belonging to that cell. Authored spatial
    /// variation finer than one cell remains the original pointwise linear
    /// entry's responsibility. The law must be pure during this solve.
    ///
    /// Smagorinsky/WALE/Vreman reuse the original full-gradient kernels, with the
    /// native 2x2 gradient embedded in 3D. The transported LOCAL density converts
    /// nu_sgs to dynamic stress; the arbitrary conditioning reference does not.
    /// Force and boundary samples, density, passive inventories, old/new inertia
    /// and reconstructed mass flux are frozen once before nonlinear iteration.
    ///
    /// The native Picard owner requires an unrelaxed velocity correction and
    /// both independently reassembled endpoint residuals to converge. Retained
    /// source samples and separated molecular/SGS rates refer to that endpoint.
    /// The original temporal/upwind/SIP accounting remains separate; no modeled
    /// transfer is silently added to a passive enthalpy slot.
    ///
    /// Cancellation is polled during preparation, between cells/assemblies and
    /// before publication, not within the original synchronous sparse solve.
    /// Histories and caller inventories remain immutable on every refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn density_momentum_bdf1_sgs_step<const N: usize>(
        &self,
        numerics: MaterialFluxNumerics,
        config: DensityMomentumSgsConfig<'_>,
        molecular_dynamic_viscosity: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>, f64) -> f64,
        force_density: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        advector: &[f64],
        previous_velocity: &[f64],
        cells: &[DensityTransportCell<N>],
        inflow: &[Option<DensityTransportCell<N>>],
        timestep_s: f64,
        controls: DensityMomentumControls,
        nonlinear: StokesRheologyControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DensityMomentumSgsStep<N>> {
        self.admit_density_momentum_sgs(config, nonlinear, &mut cancelled)?;
        let prepared = self.prepare_density_momentum_bdf1(
            numerics,
            force_density,
            boundaries,
            advector,
            previous_velocity,
            cells,
            inflow,
            timestep_s,
            controls,
            &mut cancelled,
        )?;
        self.solve_prepared_density_momentum_sgs(
            prepared,
            config,
            molecular_dynamic_viscosity,
            nonlinear,
            cancelled,
        )
    }

    /// Advance density and all carried inventories with conservative SSPRK2,
    /// then solve one native BDF1 macro momentum endpoint with nonlinear
    /// molecular rheology and Smagorinsky, WALE or Vreman stress.
    ///
    /// The original RK owner dyadically subcycles when required by its stage
    /// Courant limit. It uses this system's configured density/inventory
    /// reconstruction, positive density, bounded specific inventories and fixed
    /// volume advector/exterior donors. Its actual RK-weighted mass flux, original
    /// inertia and final transported density feed every constitutive assembly;
    /// the flux is never rebuilt from the final cells or a Picard velocity.
    /// No momentum solves or SGS heat deposits occur inside density stages.
    ///
    /// The molecular source, full-gradient SGS kernels, retained endpoint
    /// samples and separated physical/numerical transfer rates have the same
    /// contract as [`Self::density_momentum_bdf1_sgs_step`]. Every molecular
    /// query sees the SAME final density and specific inventories throughout
    /// nonlinear iteration, with the current candidate strain rate. Force and
    /// boundary data are sampled once before that iteration.
    ///
    /// Only the frozen-advector density/inventory ODE is second order in time.
    /// Momentum retains its original first-order BDF1 map and the planar SGS
    /// screening limits. Failure in RK transport, source evaluation or nonlinear
    /// convergence, and cancellation, publish no partial coupled candidate.
    #[allow(clippy::too_many_arguments)]
    pub fn density_momentum_ssprk2_bdf1_sgs_step<const N: usize>(
        &self,
        numerics: MaterialFluxNumerics,
        config: DensityMomentumSgsConfig<'_>,
        molecular_dynamic_viscosity: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>, f64) -> f64,
        force_density: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        advector: &[f64],
        previous_velocity: &[f64],
        cells: &[DensityTransportCell<N>],
        inflow: &[Option<DensityTransportCell<N>>],
        timestep_s: f64,
        controls: DensityMomentumSspRk2Controls,
        nonlinear: StokesRheologyControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DensityMomentumSspRk2SgsStep<N>> {
        self.admit_density_momentum_sgs(config, nonlinear, &mut cancelled)?;
        let (prepared, density_stages) = self.prepare_density_momentum_ssprk2_bdf1(
            numerics,
            force_density,
            boundaries,
            advector,
            previous_velocity,
            cells,
            inflow,
            timestep_s,
            controls,
            &mut cancelled,
        )?;
        Ok(DensityMomentumSspRk2SgsStep {
            momentum: self.solve_prepared_density_momentum_sgs(
                prepared,
                config,
                molecular_dynamic_viscosity,
                nonlinear,
                cancelled,
            )?,
            density_stages,
        })
    }

    fn admit_density_momentum_sgs(
        &self,
        config: DensityMomentumSgsConfig<'_>,
        nonlinear: StokesRheologyControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<()> {
        poll(&mut cancelled)?;
        config.model.validate().map_err(|_| {
            PlanarSgsError::InvalidInput("finite nonnegative SGS coefficient required")
        })?;
        if config.filter_widths_m.len() != self.mesh.tris.len() || self.mesh.tris.is_empty() {
            return Err(PlanarSgsError::InvalidInput(
                "one explicit filter width per nonempty original triangle mesh required",
            )
            .into());
        }
        for (cell, &width) in config.filter_widths_m.iter().enumerate() {
            poll(&mut cancelled)?;
            if !width.is_finite()
                || width <= 0.0
                || !self.mesh.areas[cell].is_finite()
                || self.mesh.areas[cell] <= 0.0
            {
                return Err(PlanarSgsError::InvalidCell {
                    cell,
                    reason: "positive finite filter width and triangle area required",
                }
                .into());
            }
        }
        if nonlinear.maximum_iterations == 0
            || !nonlinear.relative_tolerance.is_finite()
            || nonlinear.relative_tolerance <= 0.0
            || !nonlinear.relaxation.is_finite()
            || nonlinear.relaxation <= 0.0
            || nonlinear.relaxation > 1.0
        {
            return Err(StokesBoundaryError::InvalidInput(
                "positive nonlinear budget/tolerance and relaxation in (0,1] required",
            )
            .into());
        }
        Ok(())
    }

    // One nonlinear endpoint owner for either accepted density integrator.
    // Its conservative transport data stay frozen through all Picard updates.
    fn solve_prepared_density_momentum_sgs<'a, const N: usize>(
        &'a self,
        prepared: PreparedDensityMomentum<'a, 'm, N>,
        config: DensityMomentumSgsConfig<'_>,
        molecular_dynamic_viscosity: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>, f64) -> f64,
        nonlinear: StokesRheologyControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<DensityMomentumSgsStep<N>> {
        let cancelled = RefCell::new(cancelled);
        let failure = RefCell::new(None);
        let retained = RefCell::new(allocate(
            self.mesh.tris.len(),
            DensityMomentumSgsCell::default(),
        )?);
        let assemble = |velocity: &[f64]| {
            let candidate = (|| -> Result<_> {
                poll(&mut *cancelled.borrow_mut())?;
                // The source/SGS field is one reusable native-order buffer.
                // A refused assembly ends the solve, so only successful
                // candidates restore it for the next Picard update.
                let mut response = core::mem::take(&mut *retained.borrow_mut());
                for (cell, output) in response.iter_mut().enumerate() {
                    poll(&mut *cancelled.borrow_mut())?;
                    let kinematics = self
                        .try_cell_kinematics(velocity, cell)
                        .ok_or(PlanarSgsError::UnrepresentableTransport { cell })?;
                    let rate = kinematics.shear_rate();
                    let gradient = kinematics.gradient();
                    let eddy = eddy_viscosity(
                        config.model,
                        config.filter_widths_m[cell],
                        [
                            [gradient[0][0], gradient[0][1], 0.0],
                            [gradient[1][0], gradient[1][1], 0.0],
                            [0.0; 3],
                        ],
                    )
                    .map_err(|_| PlanarSgsError::UnrepresentableTransport { cell })?;
                    let point = self.mesh.centroids[cell];
                    let molecular =
                        molecular_dynamic_viscosity(cell, point, &prepared.cells[cell], rate);
                    if !molecular.is_finite() || molecular <= 0.0 {
                        return Err(StokesBoundaryError::InvalidViscosity {
                            cell,
                            point,
                            value: molecular,
                        }
                        .into());
                    }
                    let modeled = prepared.cells[cell].density_kg_m3 * eddy;
                    let effective = molecular + modeled;
                    if !modeled.is_finite()
                        || (eddy > 0.0 && modeled == 0.0)
                        || !effective.is_finite()
                    {
                        return Err(PlanarSgsError::UnrepresentableTransport { cell }.into());
                    }
                    *output = DensityMomentumSgsCell {
                        strain_rate_per_s: rate,
                        molecular_dynamic_viscosity_pa_s: molecular,
                        eddy_viscosity_m2_s: eddy,
                        modeled_dynamic_viscosity_pa_s: modeled,
                        effective_dynamic_viscosity_pa_s: effective,
                        ..DensityMomentumSgsCell::default()
                    };
                }
                let momentum = prepared.assemble(
                    &|cell, _| response[cell].effective_dynamic_viscosity_pa_s,
                    prepared.current_weights.clone(),
                )?;
                poll(&mut *cancelled.borrow_mut())?;
                *retained.borrow_mut() = response;
                Ok(momentum)
            })();
            candidate.map_err(|error| {
                *failure.borrow_mut() = Some(error);
                StokesBoundaryError::InvalidInput(
                    "variable-density SGS constitutive assembly refused",
                )
            })
        };
        let solved = self.solve_mixed_picard(prepared.old_velocity, nonlinear, &assemble, &|| {
            failure.borrow().is_none()
        });
        if let Some(error) = failure.into_inner() {
            return Err(error);
        }
        let solved = solved?;
        let mut response = retained.into_inner();
        let mut totals = super::Sums::<2>::new();
        for (cell, output) in response.iter_mut().enumerate() {
            poll(&mut *cancelled.borrow_mut())?;
            let rate = output.strain_rate_per_s;
            output.molecular_heating_w_m3 =
                checked_power(cell, output.molecular_dynamic_viscosity_pa_s, rate)?;
            output.modeled_transfer_w_m3 =
                checked_power(cell, output.modeled_dynamic_viscosity_pa_s, rate)?;
            let powers = [output.molecular_heating_w_m3, output.modeled_transfer_w_m3];
            let mut contribution = [0.0; 2];
            for i in 0..2 {
                contribution[i] = (powers[i] * self.mesh.areas[cell]) * prepared.controls.depth_m;
                if !contribution[i].is_finite() || (powers[i] > 0.0 && contribution[i] == 0.0) {
                    return Err(PlanarSgsError::UnrepresentableDissipation { cell }.into());
                }
            }
            totals.add(contribution);
        }
        let totals = totals.values();
        if totals.iter().any(|value| !value.is_finite()) {
            return Err(
                DensityMomentumError::Invalid("unrepresentable integrated SGS power").into(),
            );
        }
        let step = observe_step(
            self,
            solved.solution,
            &prepared.cells,
            prepared.transport,
            prepared.old_cells,
            prepared.old_velocity,
            &prepared.conservative.faces,
            prepared.controls,
            &mut *cancelled.borrow_mut(),
        )?;
        poll(&mut *cancelled.borrow_mut())?;
        Ok(DensityMomentumSgsStep {
            step,
            cells: response,
            molecular_heating_w: totals[0],
            modeled_transfer_w: totals[1],
            iterations: solved.iterations,
            relative_update: solved.relative_update,
            nonlinear_residual: solved.nonlinear_residual,
        })
    }
}

fn checked_power(cell: usize, viscosity: f64, rate: f64) -> Result<f64> {
    let power = (viscosity * rate) * rate;
    if !power.is_finite() || (viscosity > 0.0 && rate > 0.0 && power == 0.0) {
        Err(PlanarSgsError::UnrepresentableDissipation { cell }.into())
    } else {
        Ok(power)
    }
}
