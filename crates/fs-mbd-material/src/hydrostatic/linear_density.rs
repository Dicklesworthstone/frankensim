//! Continuous stable stratification on the original clipped-volume owner.
//!
//! With downward depth d, rho(d) = rho_surface + beta*d below the free surface.
//! The same degree-two integrals supply mass, mass centroid and pressure
//! potential: U = |g| integral(rho_surface*d + beta*d²/2) dV. Consequently
//! buoyancy and torque are the negative pose derivatives of that potential.

use super::{
    DisplacementMesh, HydrostaticError, HydrostaticLoad, HydrostaticMedium, HydrostaticReservoir,
    HydrostaticStep, HydrostaticStepControls, HydrostaticStepError, finite, finite_vec, poll,
};
use crate::principal::{PrincipalBody, PrincipalFrame};
use fs_ga::Vec3;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_mbd::{Gravity, RigidBodyState};

/// Prescribed hydrostatic fluid whose density increases linearly with depth.
///
/// The surface density is resolved once from the original material source;
/// `density_gradient_kg_m4` is an explicitly authored nonnegative spatial
/// gradient, not inferred from temperature, pressure or an equation of state.
/// Density above the free surface is neglected. This frozen profile supplies
/// buoyant restoring motion through the native rigid-body integrator without
/// mixing, entrainment, diffusion, drag, added mass, waves or fluid inertia.
/// Source accuracy, profile applicability and envelope coverage remain caller
/// obligations. The construction does not certify the physical stratification.
#[derive(Debug)]
pub struct LinearDensityHydrostaticMedium<'a, F: FluidTransport + ?Sized> {
    surface: HydrostaticMedium<'a, F>,
    density_gradient_kg_m4: f64,
}

impl<'a, F: FluidTransport + ?Sized> LinearDensityHydrostaticMedium<'a, F> {
    /// Resolve the surface material and admit a finite, stably stratified law.
    /// Depth is measured parallel to gravity from `surface_point_world` [m].
    /// A zero gradient reuses the exact homogeneous implementation.
    pub fn new(
        surface_source: &'a F,
        density_gradient_kg_m4: f64,
        gravity: Gravity,
        surface_point_world: Vec3,
    ) -> Result<Self, HydrostaticError> {
        if !density_gradient_kg_m4.is_finite() || density_gradient_kg_m4 < 0.0 {
            return Err(HydrostaticError::Invalid(
                "hydrostatic density gradient must be finite and nonnegative",
            ));
        }
        Ok(Self {
            surface: HydrostaticMedium::new(surface_source, gravity, surface_point_world)?,
            density_gradient_kg_m4,
        })
    }

    /// Original owner of the density at the free surface, retaining its receipts.
    pub const fn surface_source(&self) -> &'a F {
        self.surface.source()
    }
    /// Frozen coherent transport view at the free surface, not at depth.
    pub const fn surface_transport(&self) -> FluidTransportState {
        self.surface.transport()
    }
    /// Authored stable density slope [kg/m4] in the direction of gravity.
    pub const fn density_gradient_kg_m4(&self) -> f64 {
        self.density_gradient_kg_m4
    }
    /// Gravity used exactly once for body weight by the native mechanics owner.
    pub const fn gravity(&self) -> Gravity {
        self.surface.gravity()
    }
    /// Point on the prescribed planar free surface [m].
    pub const fn surface_point_world(&self) -> Vec3 {
        self.surface.surface_point_world()
    }

    /// Positive density at a nonnegative downward depth [m]. Negative depth,
    /// nonfinite input or an unrepresentable density refuses explicitly.
    pub fn density_at_depth_kg_m3(&self, depth_m: f64) -> Result<f64, HydrostaticError> {
        if !depth_m.is_finite() || depth_m < 0.0 {
            return Err(HydrostaticError::Invalid(
                "fluid density requires finite nonnegative depth",
            ));
        }
        finite(
            self.density_gradient_kg_m4
                .mul_add(depth_m, self.surface.transport.density_kg_m3()),
        )
    }

    /// Integrate displaced fluid mass and its centroid over the actual immersed
    /// geometry. A centroid density sample alone would give the wrong torque.
    /// Work is O(cells), allocation-free, with the original per-cell cancellation.
    /// The returned wrench excludes body weight.
    pub fn load(
        &self,
        mesh: &DisplacementMesh,
        body: &PrincipalBody<'_>,
        state: RigidBodyState,
        cancelled: impl FnMut() -> bool,
    ) -> Result<HydrostaticLoad, HydrostaticError> {
        self.load_frame(mesh, body.frame(), state, cancelled)
    }

    /// Native endpoint-load Verlet motion with the original absolute energy
    /// defect gate. The density profile remains frozen across both endpoints.
    pub fn step(
        &self,
        mesh: &DisplacementMesh,
        body: &PrincipalBody<'_>,
        state: RigidBodyState,
        controls: HydrostaticStepControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<HydrostaticStep, HydrostaticStepError> {
        self.step_frame(mesh, body.frame(), state, controls, cancelled)
    }
}

impl<F: FluidTransport + ?Sized> HydrostaticReservoir for LinearDensityHydrostaticMedium<'_, F> {
    fn gravity(&self) -> Gravity {
        self.gravity()
    }

    fn load_frame(
        &self,
        mesh: &DisplacementMesh,
        body: &PrincipalFrame,
        state: RigidBodyState,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HydrostaticLoad, HydrostaticError> {
        if self.density_gradient_kg_m4 == 0.0 {
            return self.surface.load_frame(mesh, body, state, cancelled);
        }
        poll(&mut cancelled)?;
        let (displaced, depth) = mesh.submerged_depth_moments(
            body.world_from_body(state)?,
            self.surface.down_world,
            self.surface.surface_point_world,
            &mut cancelled,
        )?;
        // A tiny volume must not conceal an overflowing local density. The
        // stable affine profile reaches its maximum at a clipped vertex.
        self.density_at_depth_kg_m3(depth.maximum_depth_m)?;
        let rho = self.surface.transport.density_kg_m3();
        let uniform_mass = finite(rho * displaced.volume_m3())?;
        let gradient_mass = finite(self.density_gradient_kg_m4 * displaced.depth_integral_m4)?;
        let mass = finite(uniform_mass + gradient_mass)?;
        if displaced.volume_m3() > 0.0 && uniform_mass <= 0.0 {
            return Err(HydrostaticError::Invalid(
                "displaced fluid mass is unrepresentable",
            ));
        }
        let center = displaced
            .centroid_body_m()
            .map(|center| {
                finite_vec(center + (depth.center_body_m - center).scale(gradient_mass / mass))
            })
            .transpose()?;
        let potential = finite(
            self.surface.gravity_magnitude
                * (rho * displaced.depth_integral_m4 + 0.5 * gradient_mass * depth.mean_depth_m),
        )?;
        let load = HydrostaticLoad::from_mass_moments(
            displaced,
            mass,
            center,
            potential,
            self.gravity(),
            body,
            state,
        )?;
        poll(&mut cancelled)?;
        Ok(load)
    }
}
