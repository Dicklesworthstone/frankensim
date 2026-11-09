//! Temperature drives native open-flow momentum in the original force frame.
//!
//! Collision uses the accepted F_n. Endpoint interpretation and reservoir
//! reconstruction use F_(n+1) computed from the transported scalar plus its
//! remaining source half-kick. The resulting body impulse is trapezoidal.
//! Source-coordinate attachment is performed only at construction, preserving
//! physical velocity/stress; it is not an extra half impulse at every step.
//!
//! The last thermal kick has a prescribed zeroth moment, so its temperature
//! can be previewed before the momentum solve. A declared force-consistency
//! gate compares the actual final scalar sum with the preview before joint
//! publication. No second heat solve, post-hoc force replacement or hidden
//! deposition of buoyancy work is introduced. Scalar transport retains its
//! opening-velocity split; no new coupled temporal-order claim is made.

use super::super::forcing::ForcedOpenBoundaryGrid3;
use super::{
    EquilibriumLattice3, GridError3, OpenBoundaryGrid3, OpenBoundaryMomentumTransport3,
    OpenThermalError3, ThermalOpenBalance3, ThermalOpenError3, ThermalOpenField3, macro_cell,
};
use crate::thermal::boundary::ThermalWall;
use crate::thermal::boundary::conduction::{WallConductionControls, WallConductionLink};
pub use crate::thermal::open::buoyancy::{BoussinesqLaw3, BoussinesqLawError3};
use crate::thermal::open::immersed::ImmersedThermalBody3;
use crate::thermal::open::immersed::conduction::ImmersedConductionLink3;
use crate::thermal::open::{ThermalOpenBoundaries3, ThermalOpenFace3};
use crate::thermal::tensor::DiffusivityField;
use fs_exec::Cx;
use fs_math::dd::Dd;

mod sgs;
pub use sgs::{BoussinesqOpenSgsError3, BoussinesqOpenSgsStep3};

/// Additional resource and constitutive-consistency limits for this owner.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoussinesqOpenLimits3 {
    pub maximum_cells: usize,
    /// Payload of the three retained force fields: accepted background, accepted
    /// total force, and candidate total force. Original grid/thermal storage,
    /// caller-supplied endpoint buffers and allocator metadata/spare capacity
    /// are outside this bound.
    pub maximum_force_storage_bytes: usize,
    /// Absolute per-component force-density difference between the preview
    /// and the actual final scalar endpoint. This resolves rounding of the
    /// final equilibrium source kick, not a lagged constitutive iteration.
    pub force_consistency_tolerance: f64,
}

#[derive(Debug, Clone)]
pub enum BoussinesqOpenError3<E> {
    Flow(GridError3<E>),
    Thermal(ThermalOpenError3),
    Law {
        cell: usize,
        source: BoussinesqLawError3,
    },
    /// Face 0..6 follows the native ordering; 6 identifies internal walls.
    BoundaryLaw {
        face: usize,
        source: BoussinesqLawError3,
    },
    ForceConsistency {
        cell: usize,
        axis: usize,
        applied: f64,
        actual: f64,
        tolerance: f64,
    },
    Invalid(&'static str),
    Allocation,
}
impl<E> From<GridError3<E>> for BoussinesqOpenError3<E> {
    fn from(error: GridError3<E>) -> Self {
        Self::Flow(error)
    }
}
impl<E> From<ThermalOpenError3> for BoussinesqOpenError3<E> {
    fn from(error: ThermalOpenError3) -> Self {
        Self::Thermal(error)
    }
}
impl<E> From<OpenThermalError3<E>> for BoussinesqOpenError3<E> {
    fn from(error: OpenThermalError3<E>) -> Self {
        match error {
            OpenThermalError3::Flow(e) => Self::Flow(e),
            OpenThermalError3::Thermal(e) => Self::Thermal(e),
        }
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for BoussinesqOpenError3<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native Boussinesq open flow: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for BoussinesqOpenError3<E> {}

/// Independently measured heat/reservoir accounts and prescribed body impulse.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoussinesqOpenStep3 {
    pub momentum: OpenBoundaryMomentumTransport3,
    pub heat: ThermalOpenBalance3,
    /// Sum of (F_n+F_(n+1))/2 over original fluid cells for one lattice step.
    /// Boundary/wall momentum exchange is separate from this body source.
    pub body_impulse: [f64; 3],
    /// Sum of (B_n+B_(n+1))/2 for the authored mechanical background.
    /// This remains separate from thermal buoyancy and wall/reservoir transport.
    pub background_impulse: [f64; 3],
    /// Body impulse after subtracting the trapezoidal mechanical background.
    pub buoyancy_impulse: [f64; 3],
    pub maximum_force_defect: f64,
}

/// One accepted flow/temperature/force epoch, with bounded retained scratch.
///
/// Reference heat capacity, reference-density buoyancy and external mechanical
/// force are explicit. The background can be held or supplied at the next
/// accepted time. Population density remains the native weakly compressible
/// pressure state, never the thermal density projection. There is no mutable
/// escape that detaches force history from accepted momentum or permits an
/// independent thermal time advance.
pub struct BoussinesqOpenGrid3<L: EquilibriumLattice3<Q>, const Q: usize> {
    flow: ForcedOpenBoundaryGrid3<L, Q>,
    thermal: ThermalOpenField3,
    law: BoussinesqLaw3,
    background: Vec<[f64; 3]>,
    pending_force: Vec<[f64; 3]>,
    limits: BoussinesqOpenLimits3,
    last_step: Option<BoussinesqOpenStep3>,
}

fn force_buffer<E>(cells: usize) -> Result<Vec<[f64; 3]>, BoussinesqOpenError3<E>> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(cells)
        .map_err(|_| BoussinesqOpenError3::Allocation)?;
    values.resize(cells, [0.0; 3]);
    Ok(values)
}

fn admit_boundary_temperature<E>(
    law: BoussinesqLaw3,
    face: usize,
    temperature: f64,
) -> Result<(), BoussinesqOpenError3<E>> {
    law.sample(temperature, [0.0; 3])
        .map(|_| ())
        .map_err(|source| BoussinesqOpenError3::BoundaryLaw { face, source })
}

fn admit_boundary_data<E>(
    law: BoussinesqLaw3,
    boundaries: ThermalOpenBoundaries3,
) -> Result<(), BoussinesqOpenError3<E>> {
    for (face, boundary) in boundaries.faces.into_iter().enumerate() {
        match boundary {
            ThermalOpenFace3::Reservoir { temperature }
            | ThermalOpenFace3::Wall(ThermalWall::Temperature(temperature)) => {
                admit_boundary_temperature(law, face, temperature)?;
            }
            _ => {}
        }
    }
    if let ThermalWall::Temperature(temperature) = boundaries.internal_wall {
        admit_boundary_temperature(law, 6, temperature)?;
    }
    Ok(())
}

fn admit_boundary_traces<E>(
    law: BoussinesqLaw3,
    balance: ThermalOpenBalance3,
) -> Result<(), BoussinesqOpenError3<E>> {
    for (face, reservoir) in balance.reservoirs.into_iter().enumerate() {
        if reservoir.links != 0 {
            admit_boundary_temperature(law, face, reservoir.min_temperature)?;
            admit_boundary_temperature(law, face, reservoir.max_temperature)?;
        }
    }
    for (face, wall) in balance
        .walls
        .into_iter()
        .chain([balance.internal_wall])
        .enumerate()
    {
        if let Some(wall) = wall {
            admit_boundary_temperature(law, face, wall.min_temperature)?;
            admit_boundary_temperature(law, face, wall.max_temperature)?;
        }
    }
    Ok(())
}

impl<L: EquilibriumLattice3<Q>, const Q: usize> BoussinesqOpenGrid3<L, Q> {
    /// Bind an existing matching native flow and thermal field without an
    /// initial velocity kick. Solid entries of the background must be zero.
    /// All finite wall/conduction controls and accepted thermal history survive.
    pub fn try_new(
        flow: OpenBoundaryGrid3<L, Q>,
        thermal: ThermalOpenField3,
        cx: &Cx<'_>,
        law: BoussinesqLaw3,
        background: Vec<[f64; 3]>,
        limits: BoussinesqOpenLimits3,
    ) -> Result<Self, BoussinesqOpenError3<L::Error>> {
        Self::try_new_with_poll(flow, thermal, law, background, limits, &mut || {
            cx.checkpoint().map_err(|_| GridError3::Cancelled)
        })
    }

    pub(super) fn try_new_with_poll(
        flow: OpenBoundaryGrid3<L, Q>,
        thermal: ThermalOpenField3,
        law: BoussinesqLaw3,
        background: Vec<[f64; 3]>,
        limits: BoussinesqOpenLimits3,
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<Self, BoussinesqOpenError3<L::Error>> {
        poll()?;
        let cells = flow.grid.f.len();
        let bytes = cells
            .checked_mul(3 * core::mem::size_of::<[f64; 3]>())
            .ok_or(BoussinesqOpenError3::Invalid("force storage overflow"))?;
        if cells > limits.maximum_cells
            || bytes > limits.maximum_force_storage_bytes
            || !limits.force_consistency_tolerance.is_finite()
            || limits.force_consistency_tolerance < 0.0
            || background.len() != cells
        {
            return Err(BoussinesqOpenError3::Invalid(
                "force resource cap, background shape or consistency tolerance",
            ));
        }
        if thermal.dims != flow.grid.dims
            || thermal.solid != flow.grid.solid
            || thermal.completed_steps() != flow.grid.completed_steps
        {
            return Err(BoussinesqOpenError3::Invalid(
                "native and thermal geometry or accepted time differ",
            ));
        }
        flow.admit_thermal_faces(thermal.boundaries())?;
        admit_boundary_data(law, thermal.boundaries())?;
        let mut initial_force = force_buffer(cells)?;
        let pending_force = force_buffer(cells)?;
        for cell in 0..cells {
            if cell % 64 == 0 {
                poll()?;
            }
            if background[cell].iter().any(|v| !v.is_finite())
                || (thermal.solid[cell] && background[cell] != [0.0; 3])
            {
                return Err(BoussinesqOpenError3::Invalid(
                    "nonfinite background or force on a solid cell",
                ));
            }
            if let Some(temperature) = thermal.temperature(cell) {
                initial_force[cell] = law
                    .sample(temperature, background[cell])
                    .map_err(|source| BoussinesqOpenError3::Law { cell, source })?
                    .total;
            }
        }
        let flow = ForcedOpenBoundaryGrid3::try_new_preserving_velocity_with_poll(
            flow,
            initial_force,
            poll,
        )?;
        poll()?;
        Ok(Self {
            flow,
            thermal,
            law,
            background,
            pending_force,
            limits,
            last_step: None,
        })
    }

    pub fn flow(&self) -> &ForcedOpenBoundaryGrid3<L, Q> {
        &self.flow
    }
    pub fn thermal(&self) -> &ThermalOpenField3 {
        &self.thermal
    }
    pub fn law(&self) -> BoussinesqLaw3 {
        self.law
    }
    /// Mechanical background in the same accepted epoch as total force and heat.
    pub fn background_forces(&self) -> &[[f64; 3]] {
        &self.background
    }
    pub fn last_step(&self) -> Option<BoussinesqOpenStep3> {
        self.last_step
    }

    /// Change held reservoir heat data without detaching accepted temperatures.
    pub fn set_reservoir_temperature(
        &mut self,
        face: usize,
        temperature: f64,
    ) -> Result<(), BoussinesqOpenError3<L::Error>> {
        admit_boundary_temperature(self.law, face, temperature)?;
        self.thermal
            .set_reservoir_temperature(face, temperature)
            .map_err(Into::into)
    }
    /// Change an existing stateless wall's held thermal input. Prescribed
    /// temperatures also obey the buoyancy envelope; film/flux laws are
    /// checked against their actual fluid traces during the original step.
    /// Finite wall storage cannot be replaced through this input setter.
    pub fn set_wall_boundary(
        &mut self,
        face: usize,
        wall: ThermalWall,
    ) -> Result<(), BoussinesqOpenError3<L::Error>> {
        if let ThermalWall::Temperature(temperature) = wall {
            admit_boundary_temperature(self.law, face, temperature)?;
        }
        self.thermal
            .set_wall_boundary(face, wall)
            .map_err(Into::into)
    }

    /// Update only the stateless law on unassigned fluid-to-solid links.
    /// Original immersed bodies retain their temperature and heat ownership.
    pub fn set_internal_wall_boundary(
        &mut self,
        wall: ThermalWall,
    ) -> Result<(), BoussinesqOpenError3<L::Error>> {
        if let ThermalWall::Temperature(temperature) = wall {
            admit_boundary_temperature(self.law, 6, temperature)?;
        }
        self.thermal
            .set_internal_wall_boundary(wall)
            .map_err(Into::into)
    }

    /// Set external heat into the retained finite wall cells, in original
    /// x-fast surface order. This changes inputs without adding energy or
    /// detaching accepted scalar, momentum or buoyancy history.
    pub fn set_finite_wall_heat(
        &mut self,
        face: usize,
        heat_per_step: &[f64],
    ) -> Result<(), ThermalOpenError3> {
        self.thermal.set_finite_wall_heat(face, heat_per_step)
    }

    pub fn set_wall_conduction(
        &mut self,
        axis: usize,
        links: &[WallConductionLink],
        controls: WallConductionControls,
    ) -> Result<(), ThermalOpenError3> {
        self.thermal.set_wall_conduction(axis, links, controls)
    }
    pub fn clear_wall_conduction(&mut self, axis: usize) -> Result<(), ThermalOpenError3> {
        self.thermal.clear_wall_conduction(axis)
    }

    /// Install shared lumped storage on original solid obstacles. The body
    /// temperature follows the thermal model's domain; its actual fluid-side
    /// traces must also satisfy the Boussinesq envelope at every step.
    pub fn set_immersed_thermal_bodies(
        &mut self,
        bodies: &[ImmersedThermalBody3<'_>],
        maximum_storage_bytes: usize,
    ) -> Result<(), ThermalOpenError3> {
        self.thermal
            .set_immersed_thermal_bodies(bodies, maximum_storage_bytes)
    }

    /// Install reciprocal solid conduction and original immersed thermal nodes
    /// together, including interior nodes without fluid interfaces. Native
    /// heat/flow publication and fluid-side buoyancy admission remain joint.
    pub fn set_immersed_thermal_network(
        &mut self,
        nodes: &[ImmersedThermalBody3<'_>],
        links: &[ImmersedConductionLink3],
        maximum_storage_bytes: usize,
        controls: WallConductionControls,
    ) -> Result<(), ThermalOpenError3> {
        self.thermal
            .set_immersed_thermal_network(nodes, links, maximum_storage_bytes, controls)
    }

    /// Prepare spatial solid conduction from the original native voxel geometry.
    /// Each selected voxel retains one capacity/source, with the same joint
    /// thermal, momentum and buoyancy publication as a manually authored graph.
    pub fn set_immersed_thermal_voxels(
        &mut self,
        voxels: &[crate::thermal::open::immersed::voxel::ImmersedThermalVoxel3],
        controls: crate::thermal::open::immersed::voxel::ImmersedVoxelControls3,
    ) -> Result<crate::thermal::open::immersed::voxel::ImmersedVoxelMap3, ThermalOpenError3> {
        self.thermal.set_immersed_thermal_voxels(voxels, controls)
    }

    /// Change the heat supplied to one whole immersed body without detaching
    /// thermal, momentum or buoyancy history.
    pub fn set_immersed_body_heat(
        &mut self,
        body: usize,
        heat_per_step: f64,
    ) -> Result<(), ThermalOpenError3> {
        self.thermal.set_immersed_body_heat(body, heat_per_step)
    }

    /// Advance heat and buoyant flow once, then publish their common epoch.
    /// Collision coefficients and heat are held over this interval. Their
    /// original numerical domains, native Mach bounds and scalar temperature,
    /// solid/conductive energy and backflow checks remain authoritative.
    pub fn step_serial(
        &mut self,
        cx: &Cx<'_>,
        collisions: Option<&[L::Collision]>,
        heat_per_step: &[f64],
        thermal_relaxation: Option<&[f64]>,
        origin: [f64; 3],
    ) -> Result<BoussinesqOpenStep3, BoussinesqOpenError3<L::Error>> {
        self.step_with_poll(
            collisions,
            heat_per_step,
            thermal_relaxation,
            origin,
            &mut || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    /// Advance heat and buoyant flow with a full SPD thermal diffusivity tensor.
    ///
    /// The borrowed field is expressed in the fixed Cartesian lattice axes and
    /// is held for this interval. Off-diagonal conduction changes the actual
    /// transported temperature used by the endpoint buoyancy force; no scalar
    /// surrogate or second heat solve is used. Uniform fields share one
    /// prepared tensor and per-cell fields cover the original grid, including
    /// solid entries. The thermal carrier's original tau still relaxes even
    /// kinetic moments, and its reference heat capacity remains fixed.
    ///
    /// This retains the existing opening-velocity thermal split, source kicks,
    /// force consistency and reservoir/wall accounts. The complete temperature,
    /// flow, finite walls and force history publish together or remain unchanged
    /// on refusal or cancellation. Tensor stability and time/spatial refinement
    /// remain caller responsibilities; this is not a variable-storage or total-
    /// energy model. A subsequent scalar step uses its own original coefficients.
    pub fn step_serial_with_tensor_diffusion(
        &mut self,
        cx: &Cx<'_>,
        collisions: Option<&[L::Collision]>,
        heat_per_step: &[f64],
        diffusion: DiffusivityField<'_, 3>,
        origin: [f64; 3],
    ) -> Result<BoussinesqOpenStep3, BoussinesqOpenError3<L::Error>> {
        self.step_with_diffusion_poll(
            collisions,
            heat_per_step,
            None,
            Some(diffusion),
            origin,
            &mut || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    /// Advance mixed buoyant/forced convection with a prescribed mechanical
    /// background at the next accepted time.
    ///
    /// Total force is B + F_b(T): collision uses its retained opening epoch,
    /// while native reconstruction and the final heat-source kick use the
    /// endpoint background plus the buoyancy preview from transported heat.
    /// The body/background/buoyancy impulses retain the trapezoidal source
    /// convention. Neither mechanical work nor buoyancy work is deposited
    /// as heat. The original opening-velocity thermal split is unchanged.
    ///
    /// The endpoint field must cover every original cell, be finite and vanish
    /// on solid cells. Successful publication exchanges it with the previous
    /// background buffer for reuse. Refusal or cancellation leaves this input
    /// buffer and every accepted flow/thermal/force/wall epoch unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_serial(
        &mut self,
        cx: &Cx<'_>,
        collisions: Option<&[L::Collision]>,
        background_at_end: &mut Vec<[f64; 3]>,
        heat_per_step: &[f64],
        thermal_relaxation: Option<&[f64]>,
        origin: [f64; 3],
    ) -> Result<BoussinesqOpenStep3, BoussinesqOpenError3<L::Error>> {
        self.advance_with_diffusion_poll(
            Some(background_at_end),
            collisions,
            heat_per_step,
            thermal_relaxation,
            None,
            origin,
            &mut || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    /// Combine a changing mechanical background with the original full-tensor
    /// heat transport and temperature-dependent buoyancy in one transaction.
    /// Tensor axes, capacity, limits and history retain their native meanings;
    /// background-buffer exchange matches [`Self::advance_serial`].
    #[allow(clippy::too_many_arguments)]
    pub fn advance_serial_with_tensor_diffusion(
        &mut self,
        cx: &Cx<'_>,
        collisions: Option<&[L::Collision]>,
        background_at_end: &mut Vec<[f64; 3]>,
        heat_per_step: &[f64],
        diffusion: DiffusivityField<'_, 3>,
        origin: [f64; 3],
    ) -> Result<BoussinesqOpenStep3, BoussinesqOpenError3<L::Error>> {
        self.advance_with_diffusion_poll(
            Some(background_at_end),
            collisions,
            heat_per_step,
            None,
            Some(diffusion),
            origin,
            &mut || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    pub(super) fn step_with_poll(
        &mut self,
        collisions: Option<&[L::Collision]>,
        heat_per_step: &[f64],
        thermal_relaxation: Option<&[f64]>,
        origin: [f64; 3],
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<BoussinesqOpenStep3, BoussinesqOpenError3<L::Error>> {
        self.step_with_diffusion_poll(
            collisions,
            heat_per_step,
            thermal_relaxation,
            None,
            origin,
            poll,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_with_diffusion_poll(
        &mut self,
        collisions: Option<&[L::Collision]>,
        heat_per_step: &[f64],
        thermal_relaxation: Option<&[f64]>,
        diffusion: Option<DiffusivityField<'_, 3>>,
        origin: [f64; 3],
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<BoussinesqOpenStep3, BoussinesqOpenError3<L::Error>> {
        self.advance_with_diffusion_poll(
            None,
            collisions,
            heat_per_step,
            thermal_relaxation,
            diffusion,
            origin,
            poll,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_with_diffusion_poll(
        &mut self,
        background_at_end: Option<&mut Vec<[f64; 3]>>,
        collisions: Option<&[L::Collision]>,
        heat_per_step: &[f64],
        thermal_relaxation: Option<&[f64]>,
        diffusion: Option<DiffusivityField<'_, 3>>,
        origin: [f64; 3],
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<BoussinesqOpenStep3, BoussinesqOpenError3<L::Error>> {
        poll()?;
        let cells = self.thermal.populations().len();
        let endpoint_background = background_at_end
            .as_deref()
            .map_or(self.background.as_slice(), Vec::as_slice);
        if endpoint_background.len() != cells {
            return Err(BoussinesqOpenError3::Invalid(
                "endpoint background field shape",
            ));
        }
        for cell in 0..cells {
            if cell % 64 == 0 {
                poll()?;
            }
            if background_at_end.is_some()
                && (endpoint_background[cell].iter().any(|v| !v.is_finite())
                    || (self.thermal.solid[cell] && endpoint_background[cell] != [0.0; 3]))
            {
                return Err(BoussinesqOpenError3::Invalid(
                    "nonfinite endpoint background or force on a solid cell",
                ));
            }
            self.thermal.velocities[cell] = if self.thermal.solid[cell] {
                [0.0; 3]
            } else {
                self.flow.macroscopic(cell)?.1
            };
        }
        let candidate = if let Some(diffusion) = diffusion {
            self.thermal
                .begin_step_tensor::<BoussinesqOpenError3<L::Error>>(
                    heat_per_step,
                    diffusion,
                    || poll().map_err(Into::into),
                )?
        } else {
            self.thermal.begin_step::<BoussinesqOpenError3<L::Error>>(
                heat_per_step,
                thermal_relaxation,
                || poll().map_err(Into::into),
            )?
        };
        // These are actual fluid traces, including outflow and implicit wall
        // conduction. A finite solid's bulk temperature is a separate model.
        admit_boundary_traces(self.law, candidate)?;
        for cell in 0..cells {
            if cell % 64 == 0 {
                poll()?;
            }
            if self.thermal.solid[cell] {
                self.pending_force[cell] = [0.0; 3];
                continue;
            }
            let opening = self.thermal.temperature(cell).expect("active thermal cell");
            let transported = self
                .thermal
                .trial_temperature(cell)
                .expect("active thermal cell");
            let half_heat = 0.5 * heat_per_step[cell];
            // The approximation envelope covers source substeps as well as
            // retained endpoints, even if wall cooling later removes a spike.
            for temperature in [opening + half_heat, transported] {
                self.law
                    .sample(temperature, self.background[cell])
                    .map_err(|source| BoussinesqOpenError3::Law { cell, source })?;
            }
            self.pending_force[cell] = self
                .law
                .sample(transported + half_heat, endpoint_background[cell])
                .map_err(|source| BoussinesqOpenError3::Law { cell, source })?
                .total;
        }
        let (momentum, step) =
            self.flow
                .stage_serial_with_poll(&self.pending_force, collisions, origin, poll)?;
        for cell in 0..cells {
            if cell % 64 == 0 {
                poll()?;
            }
            if !self.thermal.solid[cell] {
                self.thermal.velocities[cell] = macro_cell::<L, Q>(
                    &self.flow.inner.grid.next[cell],
                    self.pending_force[cell],
                    cell,
                    self.flow.inner.grid.budget.maximum_mach,
                )?
                .1;
            }
        }
        let heat = self.thermal.finish_step::<BoussinesqOpenError3<L::Error>>(
            heat_per_step,
            candidate,
            || poll().map_err(Into::into),
        )?;
        let mut body = [Dd::ZERO; 3];
        let mut background = [Dd::ZERO; 3];
        let mut maximum_force_defect = 0.0_f64;
        for cell in 0..cells {
            if cell % 64 == 0 {
                poll()?;
            }
            let Some(temperature) = self.thermal.trial_temperature(cell) else {
                continue;
            };
            let actual = self
                .law
                .sample(temperature, endpoint_background[cell])
                .map_err(|source| BoussinesqOpenError3::Law { cell, source })?
                .total;
            for axis in 0..3 {
                let applied = self.pending_force[cell][axis];
                let defect = (actual[axis] - applied).abs();
                if !defect.is_finite() || defect > self.limits.force_consistency_tolerance {
                    return Err(BoussinesqOpenError3::ForceConsistency {
                        cell,
                        axis,
                        applied,
                        actual: actual[axis],
                        tolerance: self.limits.force_consistency_tolerance,
                    });
                }
                maximum_force_defect = maximum_force_defect.max(defect);
                body[axis] = body[axis]
                    + Dd::from_f64(0.5 * self.flow.forces[cell][axis])
                    + Dd::from_f64(0.5 * applied);
                // Preserve the held-background arithmetic, including subnormal
                // values, when both endpoint components are identical.
                if endpoint_background[cell][axis] == self.background[cell][axis] {
                    background[axis] = background[axis] + Dd::from_f64(self.background[cell][axis]);
                } else {
                    background[axis] = background[axis]
                        + Dd::from_f64(0.5 * self.background[cell][axis])
                        + Dd::from_f64(0.5 * endpoint_background[cell][axis]);
                }
            }
        }
        let body_impulse = body.map(Dd::to_f64);
        let background_impulse = background.map(Dd::to_f64);
        let buoyancy_impulse =
            core::array::from_fn(|axis| (body[axis] - background[axis]).to_f64());
        if body_impulse
            .iter()
            .chain(&background_impulse)
            .chain(&buoyancy_impulse)
            .any(|v| !v.is_finite())
        {
            return Err(BoussinesqOpenError3::Invalid(
                "summed body impulse is not representable",
            ));
        }
        let report = BoussinesqOpenStep3 {
            momentum,
            heat,
            body_impulse,
            background_impulse,
            buoyancy_impulse,
            maximum_force_defect,
        };
        poll()?;
        self.flow
            .publish_with_force_swap(step, &mut self.pending_force);
        if let Some(background_at_end) = background_at_end {
            core::mem::swap(&mut self.background, background_at_end);
        }
        self.thermal.publish(step, heat);
        self.last_step = Some(report);
        Ok(report)
    }
}

#[cfg(test)]
mod tests;
