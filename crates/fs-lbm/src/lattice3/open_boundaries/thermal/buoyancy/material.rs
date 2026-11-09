//! Source-resolved molecular transport on the original native buoyant owner.
//!
//! Viscosity and conductivity use one original material response per fluid
//! cell and query epoch. Reference density/Cp retain the scalar's sensible
//! storage. The existing affine law remains the sole buoyancy force owner.
//! Opening coefficients are held for one step; endpoint queries admit the
//! candidate before its original physical publication.

use super::*;
use crate::lattice3::ViscousLattice3;
use crate::material::transport::boussinesq::{BoussinesqConfig, BoussinesqMaterial};
use crate::material::transport::{MaterialTransportError, TransportSample};
use crate::material::{MaterialScaling, plan_scaling_from_material};
use fs_material::{
    MaterialError,
    fluid_transport::{FluidTransport, FluidTransportState},
};

mod sgs;
pub use sgs::{BoussinesqOpenMaterialSgsError3, BoussinesqOpenMaterialSgsStep3};

/// Explicit physical coordinates, fixed-storage approximation and memory caps.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoussinesqOpenMaterialConfig3 {
    /// Kelvin at theta=0. This is an affine coordinate origin, not a reservoir.
    pub temperature_at_zero_k: f64,
    /// Signed nonzero Kelvin per theta unit; reversed coordinates are admitted.
    pub temperature_scale_k: f64,
    /// Maximum abs(rho(T)-rho_ref)/rho_ref; finite and in [0,1).
    pub maximum_relative_density_departure: f64,
    /// Maximum abs(cp(T)-cp_ref)/cp_ref; finite and in [0,1).
    pub maximum_relative_specific_heat_departure: f64,
    /// Retained payload: one native collision and two f64 values per grid cell.
    pub maximum_workspace_bytes: usize,
    /// Inline payload of BOTH returned material-sample vectors. Original source
    /// owners' internal allocations and allocator metadata are outside this cap.
    pub maximum_report_bytes: usize,
}
impl BoussinesqOpenMaterialConfig3 {
    fn projection(self) -> BoussinesqConfig<3> {
        BoussinesqConfig {
            // Only the original projection's storage admission is used here.
            // Source density is not a second gravitational force contribution.
            gravity_m_s2: [0.0; 3],
            max_relative_density_departure: self.maximum_relative_density_departure,
            max_relative_specific_heat_departure: self.maximum_relative_specific_heat_departure,
        }
    }

    fn temperature(self, theta: f64) -> Option<f64> {
        let kelvin = self.temperature_at_zero_k + theta * self.temperature_scale_k;
        (theta.is_finite() && kelvin.is_finite() && kelvin > 0.0).then_some(kelvin)
    }
}

/// Reused molecular collision/tau and temperature-query scratch.
///
/// The reference snapshot and dx/dt come from an original MaterialScaling.
/// The original source owner of each later query remains in the returned
/// BoussinesqMaterial; this workspace does not replace it with a material name.
/// Scratch is not accepted physical state and may contain a refused prefix.
pub struct BoussinesqOpenMaterialWorkspace3<C> {
    dimensions: [usize; 3],
    temperature_range: [f64; 2],
    reference: FluidTransportState,
    cell_size_m: f64,
    timestep_s: f64,
    config: BoussinesqOpenMaterialConfig3,
    collisions: Vec<C>,
    thermal_relaxation: Vec<f64>,
    temperatures_k: Vec<f64>,
}

impl<C> BoussinesqOpenMaterialWorkspace3<C> {
    pub const fn reference_transport(&self) -> &FluidTransportState {
        &self.reference
    }
    /// Opening molecular field from the most recent successful call. A later
    /// refused call can overwrite it with scratch; consume only after success.
    pub fn momentum_collisions(&self) -> &[C] {
        &self.collisions
    }
    /// Actual D3Q7 opening molecular tau, including unused solid padding.
    /// A subsequent refused call can overwrite this scratch as well.
    pub fn thermal_relaxation(&self) -> &[f64] {
        &self.thermal_relaxation
    }
}

/// Original transport publication and its source-bearing material query epochs.
#[derive(Debug)]
pub struct BoussinesqOpenMaterialStep3<F> {
    pub transport: BoussinesqOpenStep3,
    /// Opening mu/k actually used in the step, with fixed reference storage.
    pub applied: Vec<TransportSample<BoussinesqMaterial<F>>>,
    /// Actual final-temperature queries admitted before physical publication.
    /// These do not transport heat again or imply a nonlinear iteration.
    pub endpoint: Vec<TransportSample<BoussinesqMaterial<F>>>,
}

#[derive(Clone, Debug)]
pub enum BoussinesqOpenMaterialError3<E> {
    Transport(BoussinesqOpenError3<E>),
    Material(MaterialTransportError),
    Invalid(&'static str),
    Allocation,
}
impl<E> From<BoussinesqOpenError3<E>> for BoussinesqOpenMaterialError3<E> {
    fn from(error: BoussinesqOpenError3<E>) -> Self {
        Self::Transport(error)
    }
}
impl<E> From<GridError3<E>> for BoussinesqOpenMaterialError3<E> {
    fn from(error: GridError3<E>) -> Self {
        Self::Transport(error.into())
    }
}
impl<E> From<MaterialTransportError> for BoussinesqOpenMaterialError3<E> {
    fn from(error: MaterialTransportError) -> Self {
        Self::Material(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for BoussinesqOpenMaterialError3<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native buoyant material transport: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for BoussinesqOpenMaterialError3<E> {}

impl<L: EquilibriumLattice3<Q> + ViscousLattice3<Q>, const Q: usize> BoussinesqOpenGrid3<L, Q> {
    /// Bind the original reference material scaling and allocate reusable
    /// transport arrays without changing temperature, momentum or force history.
    ///
    /// The complete native theta envelope must map to finite positive Kelvin,
    /// so its original source-substep and actual boundary-trace gates remain
    /// physical under this coordinate map. The mapping does not replace those
    /// limits or infer an EOS. Supported native lattices retain cs2=1/3.
    ///
    /// Source rho/Cp departures are bounded at the opening and endpoint query
    /// epochs. Source density is NOT substituted into native populations or the
    /// independent affine buoyancy law. Native heat/source/impulse reports keep
    /// their lattice units; this adapter changes molecular transport only.
    pub fn material_transport_workspace<R: FluidTransport + ?Sized>(
        &self,
        reference: &MaterialScaling<'_, R>,
        config: BoussinesqOpenMaterialConfig3,
    ) -> Result<
        BoussinesqOpenMaterialWorkspace3<L::Collision>,
        BoussinesqOpenMaterialError3<L::Error>,
    > {
        config.projection().validate()?;
        if !config.temperature_at_zero_k.is_finite()
            || !config.temperature_scale_k.is_finite()
            || config.temperature_scale_k == 0.0
            || L::CS2 != 1.0 / 3.0
        {
            return Err(BoussinesqOpenMaterialError3::Invalid(
                "material transport requires finite nonzero Kelvin coordinates and native cs2=1/3",
            ));
        }
        let temperature_range = self.thermal.temperature_range();
        if temperature_range
            .iter()
            .any(|&theta| config.temperature(theta).is_none())
        {
            return Err(BoussinesqOpenMaterialError3::Invalid(
                "the complete native thermal envelope must map to positive finite Kelvin",
            ));
        }
        let cells = self.thermal.populations().len();
        let payload = cells
            .checked_mul(core::mem::size_of::<L::Collision>() + 2 * core::mem::size_of::<f64>())
            .ok_or(BoussinesqOpenMaterialError3::Allocation)?;
        if payload > config.maximum_workspace_bytes {
            return Err(BoussinesqOpenMaterialError3::Invalid(
                "material transport workspace exceeds its payload budget",
            ));
        }
        let mut collisions = Vec::new();
        let mut thermal_relaxation = Vec::new();
        let mut temperatures_k = Vec::new();
        collisions
            .try_reserve_exact(cells)
            .map_err(|_| BoussinesqOpenMaterialError3::Allocation)?;
        thermal_relaxation
            .try_reserve_exact(cells)
            .map_err(|_| BoussinesqOpenMaterialError3::Allocation)?;
        temperatures_k
            .try_reserve_exact(cells)
            .map_err(|_| BoussinesqOpenMaterialError3::Allocation)?;
        collisions.resize(cells, self.flow.inner.grid.collision);
        thermal_relaxation.resize(cells, 1.0);
        temperatures_k.resize(cells, 0.0);
        Ok(BoussinesqOpenMaterialWorkspace3 {
            dimensions: self.thermal.dimensions(),
            temperature_range,
            reference: *reference.transport(),
            cell_size_m: reference.cell_size_m(),
            timestep_s: reference.timestep_s(),
            config,
            collisions,
            thermal_relaxation,
            temperatures_k,
        })
    }

    /// Advance with opening-temperature molecular viscosity and conductivity.
    ///
    /// The callback returns the original source-bearing material at the supplied
    /// Kelvin temperature and retained reference state. Its transport is read
    /// once. BoussinesqMaterial admits density/Cp departures and projects mu/k
    /// onto reference storage; the existing material planner supplies native
    /// momentum tau and the distinct D3Q7 thermal tau.
    ///
    /// Optional collision templates retain their equilibrium, bulk/higher rates
    /// and TRT wall-placement product through the original ViscousLattice3
    /// update. Only molecular shear relaxation changes. No turbulence closure
    /// or temperature-dependent force law is inferred from the source.
    ///
    /// Coefficients are held for this interval. Actual endpoint samples must
    /// pass source, storage and numerical admission before the original final
    /// cancellation gate and joint flow/force/heat/wall publication. Refusal
    /// leaves accepted state intact. Callbacks must be repeatable; their own
    /// captured receipts/effects cannot be rolled back by this owner.
    ///
    /// Scratch arrays are reused. The two bounded output sample vectors are
    /// allocated per call, retaining the concrete source values without erasure.
    #[allow(clippy::too_many_arguments)]
    pub fn step_serial_with_material_transport<F: FluidTransport>(
        &mut self,
        cx: &Cx<'_>,
        collision_templates: Option<&[L::Collision]>,
        heat_per_step: &[f64],
        workspace: &mut BoussinesqOpenMaterialWorkspace3<L::Collision>,
        origin: [f64; 3],
        mut resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<BoussinesqOpenMaterialStep3<F>, BoussinesqOpenMaterialError3<L::Error>> {
        self.advance_material_with_poll(
            None,
            collision_templates,
            heat_per_step,
            workspace,
            origin,
            &mut resolve,
            &mut || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    /// The same original material transaction with a prescribed next mechanical
    /// background. Endpoint source refusal leaves the caller's background buffer
    /// unchanged; only successful joint publication exchanges it for the old one.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_serial_with_material_transport<F: FluidTransport>(
        &mut self,
        cx: &Cx<'_>,
        collision_templates: Option<&[L::Collision]>,
        background_at_end: &mut Vec<[f64; 3]>,
        heat_per_step: &[f64],
        workspace: &mut BoussinesqOpenMaterialWorkspace3<L::Collision>,
        origin: [f64; 3],
        mut resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<BoussinesqOpenMaterialStep3<F>, BoussinesqOpenMaterialError3<L::Error>> {
        self.advance_material_with_poll(
            Some(background_at_end),
            collision_templates,
            heat_per_step,
            workspace,
            origin,
            &mut resolve,
            &mut || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn advance_material_with_poll<F: FluidTransport>(
        &mut self,
        background_at_end: Option<&mut Vec<[f64; 3]>>,
        collision_templates: Option<&[L::Collision]>,
        heat_per_step: &[f64],
        workspace: &mut BoussinesqOpenMaterialWorkspace3<L::Collision>,
        origin: [f64; 3],
        resolve: &mut impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<BoussinesqOpenMaterialStep3<F>, BoussinesqOpenMaterialError3<L::Error>> {
        let applied = self.prepare_material(
            background_at_end.as_deref().map(Vec::as_slice),
            collision_templates,
            heat_per_step,
            workspace,
            resolve,
            poll,
        )?;
        let (transport, step) = self.stage_with_diffusion_poll(
            background_at_end.as_deref().map(Vec::as_slice),
            Some(&workspace.collisions),
            heat_per_step,
            Some(&workspace.thermal_relaxation),
            None,
            origin,
            poll,
        )?;
        self.finish_material_with_poll(
            background_at_end,
            collision_templates,
            workspace,
            applied,
            transport,
            step,
            resolve,
            poll,
        )
    }

    // Both molecular-only and SGS composition validate the same source payload
    // and resolve fresh opening coefficients before staging physical state.
    #[allow(clippy::too_many_arguments)]
    fn prepare_material<F: FluidTransport>(
        &self,
        background_at_end: Option<&[[f64; 3]]>,
        collision_templates: Option<&[L::Collision]>,
        heat_per_step: &[f64],
        workspace: &mut BoussinesqOpenMaterialWorkspace3<L::Collision>,
        resolve: &mut impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<Vec<TransportSample<BoussinesqMaterial<F>>>, BoussinesqOpenMaterialError3<L::Error>>
    {
        poll()?;
        let cells = self.thermal.populations().len();
        if workspace.dimensions != self.thermal.dimensions()
            || workspace.temperature_range != self.thermal.temperature_range()
            || workspace.collisions.len() != cells
            || workspace.thermal_relaxation.len() != cells
            || workspace.temperatures_k.len() != cells
            || heat_per_step.len() != cells
            || collision_templates.is_some_and(|field| field.len() != cells)
            || background_at_end.is_some_and(|field| field.len() != cells)
        {
            return Err(BoussinesqOpenMaterialError3::Invalid(
                "material workspace, geometry, temperature envelope or input field shape mismatch",
            ));
        }
        let fluid_cells = self.thermal.solid.iter().filter(|&&solid| !solid).count();
        let payload = fluid_cells
            .checked_mul(core::mem::size_of::<TransportSample<BoussinesqMaterial<F>>>())
            .and_then(|bytes| bytes.checked_mul(2))
            .ok_or(BoussinesqOpenMaterialError3::Allocation)?;
        if payload > workspace.config.maximum_report_bytes {
            return Err(BoussinesqOpenMaterialError3::Invalid(
                "material source samples exceed their explicit output payload budget",
            ));
        }
        self.sample_material(
            false,
            collision_templates,
            workspace,
            fluid_cells,
            resolve,
            poll,
        )
    }

    // This final source admission and cancellation gate is shared by every
    // material path; publication remains with the original coupled owner.
    #[allow(clippy::too_many_arguments)]
    fn finish_material_with_poll<F: FluidTransport>(
        &mut self,
        background_at_end: Option<&mut Vec<[f64; 3]>>,
        collision_templates: Option<&[L::Collision]>,
        workspace: &mut BoussinesqOpenMaterialWorkspace3<L::Collision>,
        applied: Vec<TransportSample<BoussinesqMaterial<F>>>,
        transport: BoussinesqOpenStep3,
        step: usize,
        resolve: &mut impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<BoussinesqOpenMaterialStep3<F>, BoussinesqOpenMaterialError3<L::Error>> {
        let endpoint = self.sample_material(
            true,
            collision_templates,
            workspace,
            applied.len(),
            resolve,
            poll,
        )?;
        poll()?;
        self.publish_step(step, transport, background_at_end);
        Ok(BoussinesqOpenMaterialStep3 {
            transport,
            applied,
            endpoint,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_material<F: FluidTransport>(
        &self,
        endpoint: bool,
        templates: Option<&[L::Collision]>,
        workspace: &mut BoussinesqOpenMaterialWorkspace3<L::Collision>,
        fluid_cells: usize,
        resolve: &mut impl FnMut(usize, f64, &FluidTransportState) -> Result<F, MaterialError>,
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<Vec<TransportSample<BoussinesqMaterial<F>>>, BoussinesqOpenMaterialError3<L::Error>>
    {
        // Admit all actual fluid temperatures before invoking a source. Solid
        // mask entries are neither queried nor assigned fluid storage.
        for cell in 0..self.thermal.populations().len() {
            if cell % 64 == 0 {
                poll()?;
            }
            if self.thermal.solid[cell] {
                workspace.temperatures_k[cell] = 0.0;
                if !endpoint {
                    workspace.collisions[cell] =
                        templates.map_or(self.flow.inner.grid.collision, |field| field[cell]);
                    workspace.thermal_relaxation[cell] = 1.0;
                }
                continue;
            }
            let theta = if endpoint {
                self.thermal.trial_temperature(cell)
            } else {
                self.thermal.temperature(cell)
            }
            .expect("original fluid cell");
            workspace.temperatures_k[cell] = workspace.config.temperature(theta).ok_or(
                BoussinesqOpenMaterialError3::Invalid(
                    "material query requires a positive finite Kelvin temperature",
                ),
            )?;
        }
        let mut samples = Vec::new();
        samples
            .try_reserve_exact(fluid_cells)
            .map_err(|_| BoussinesqOpenMaterialError3::Allocation)?;
        for cell in 0..self.thermal.populations().len() {
            if self.thermal.solid[cell] {
                continue;
            }
            poll()?;
            let temperature_k = workspace.temperatures_k[cell];
            let state_error = |source| MaterialTransportError::State {
                cell,
                temperature_k,
                source,
            };
            let material = resolve(cell, temperature_k, &workspace.reference)
                .map_err(|source| state_error(source.into()))?;
            let material = BoussinesqMaterial::try_new(
                material,
                &workspace.reference,
                workspace.config.projection(),
            )
            .map_err(|source| state_error(source.into()))?;
            let transport = *material.reference_storage_transport();
            // Query the frozen snapshot, not the original material a second time.
            let scaling = plan_scaling_from_material(
                &transport,
                workspace.cell_size_m,
                workspace.timestep_s,
                0.0,
            )
            .map_err(state_error)?;
            let momentum_relaxation = scaling.flow().tau;
            let thermal_relaxation = scaling.thermal_relaxation_time3().map_err(state_error)?;
            let original = templates.map_or(self.flow.inner.grid.collision, |field| field[cell]);
            let collision = L::with_shear_relaxation(original, momentum_relaxation)
                .map_err(|source| GridError3::Collision { cell, source })?;
            if !endpoint {
                workspace.collisions[cell] = collision;
                workspace.thermal_relaxation[cell] = thermal_relaxation;
            }
            samples.push(TransportSample {
                cell,
                temperature_k,
                material,
                transport,
                momentum_relaxation,
                thermal_relaxation,
            });
        }
        Ok(samples)
    }
}
