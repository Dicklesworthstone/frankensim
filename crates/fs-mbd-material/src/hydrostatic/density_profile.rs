//! Continuous density profiles on the original hydrostatic geometry and motion.
//!
//! In each slab, with local downward depth s, rho = rho_top + delta_rho*s/h
//! and p/|g| = head_top + rho_top*s + delta_rho*s²/(2h). Direct two-plane
//! clipping integrates positive volume/depth moments, without subtracting
//! almost equal deep volumes or signed changes of slope. Thus density can
//! flatten or steepen through multiple thermoclines without cancellation of
//! large fictitious masses. Force, torque and potential use the same moments.

use super::{
    DisplacementMesh, HydrostaticError, HydrostaticLoad, HydrostaticMedium, HydrostaticReservoir,
    HydrostaticStep, HydrostaticStepControls, HydrostaticStepError, finite, finite_vec, poll,
};
use crate::principal::{PrincipalBody, PrincipalFrame};
use fs_ga::Vec3;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_mbd::{Gravity, RigidBodyState};

/// Maximum source samples in one frozen continuous density profile.
pub const MAX_DENSITY_SAMPLES: usize = 32;

/// One original fluid source at a prescribed downward depth from the free surface.
#[derive(Clone, Copy)]
pub struct HydrostaticDensitySample<'a> {
    /// Nonnegative depth [m], measured along the same gravity as mechanics.
    pub depth_m: f64,
    /// Original material owner evaluated once during profile preparation.
    pub source: &'a dyn FluidTransport,
}

struct ResolvedSample<'a> {
    input: HydrostaticDensitySample<'a>,
    transport: FluidTransportState,
    // Integral of density from the free surface to this sample [kg/m²].
    pressure_head: f64,
}

/// Source-bound, stable, piecewise-linear hydrostatic density stratification.
///
/// The first sample is at depth zero; density is linearly interpolated between
/// strictly ordered samples and held CONSTANT below the final sample. Density
/// above the free surface is neglected. Density must not decrease downward,
/// although its gradient may change in either direction. A single sample, or
/// identical sample densities, delegates exactly to the uniform-fluid owner.
///
/// These are prescribed spatial material snapshots, not a solved pressure,
/// composition or temperature EOS. The caller owns profile applicability and
/// the constant continuation. No mixing, fluid inertia, drag, waves or moving
/// thermocline is inferred. Body reference mass remains owned by mechanics.
pub struct DensityProfileHydrostaticMedium<'a> {
    surface: HydrostaticMedium<'a, dyn FluidTransport + 'a>,
    samples: Vec<ResolvedSample<'a>>,
    uniform: bool,
}

impl<'a> DensityProfileHydrostaticMedium<'a> {
    /// Resolve the original sources and prepare continuous hydrostatic pressure.
    /// All sample densities must be positive and nondecreasing. Preparation is
    /// O(samples); its allocation and source queries never recur during loading.
    pub fn new(
        samples: &[HydrostaticDensitySample<'a>],
        gravity: Gravity,
        surface_point_world: Vec3,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, HydrostaticError> {
        poll(&mut cancelled)?;
        if samples.is_empty() || samples.len() > MAX_DENSITY_SAMPLES || samples[0].depth_m != 0.0 {
            return Err(HydrostaticError::Invalid(
                "density profile needs 1..=32 samples starting at zero depth",
            ));
        }
        let surface = HydrostaticMedium::new(samples[0].source, gravity, surface_point_world)?;
        let mut resolved = Vec::with_capacity(samples.len());
        let mut previous_depth = 0.0;
        let mut previous_density = surface.transport.density_kg_m3();
        let mut pressure_head = 0.0;
        let mut uniform = true;
        for (index, &input) in samples.iter().enumerate() {
            poll(&mut cancelled)?;
            if !input.depth_m.is_finite()
                || input.depth_m < 0.0
                || index > 0 && input.depth_m <= previous_depth
            {
                return Err(HydrostaticError::Invalid(
                    "density profile depths must be finite and strictly increasing",
                ));
            }
            let transport = if index == 0 {
                surface.transport
            } else {
                input.source.transport_state()?
            };
            let density = transport.density_kg_m3();
            if density < previous_density {
                return Err(HydrostaticError::Invalid(
                    "hydrostatic density must not decrease downward",
                ));
            }
            let thickness = input.depth_m - previous_depth;
            pressure_head =
                finite(pressure_head + thickness * (0.5 * previous_density + 0.5 * density))?;
            uniform &= density == surface.transport.density_kg_m3();
            resolved.push(ResolvedSample {
                input,
                transport,
                pressure_head,
            });
            previous_depth = input.depth_m;
            previous_density = density;
        }
        poll(&mut cancelled)?;
        Ok(Self {
            surface,
            samples: resolved,
            uniform,
        })
    }

    /// Original sources and their consumed transport snapshots, in depth order.
    pub fn samples(
        &self,
    ) -> impl ExactSizeIterator<Item = (HydrostaticDensitySample<'a>, FluidTransportState)> + '_
    {
        self.samples
            .iter()
            .map(|sample| (sample.input, sample.transport))
    }

    /// Gravity used separately and exactly once for body weight.
    pub const fn gravity(&self) -> Gravity {
        self.surface.gravity()
    }

    /// Point on the prescribed free-surface plane [m].
    pub const fn surface_point_world(&self) -> Vec3 {
        self.surface.surface_point_world()
    }

    /// Piecewise-linear density [kg/m³], with the declared constant deep tail.
    pub fn density_at_depth_kg_m3(&self, depth_m: f64) -> Result<f64, HydrostaticError> {
        if !depth_m.is_finite() || depth_m < 0.0 {
            return Err(HydrostaticError::Invalid(
                "fluid density requires finite nonnegative depth",
            ));
        }
        let index = self
            .samples
            .partition_point(|sample| sample.input.depth_m <= depth_m)
            - 1;
        let sample = &self.samples[index];
        let density = sample.transport.density_kg_m3();
        match self.samples.get(index + 1) {
            None => Ok(density),
            Some(next) => {
                let fraction =
                    (depth_m - sample.input.depth_m) / (next.input.depth_m - sample.input.depth_m);
                finite((next.transport.density_kg_m3() - density).mul_add(fraction, density))
            }
        }
    }

    /// Integrate actual displaced mass, its center, force, torque and pressure
    /// potential. O(cells × samples), with no load-time allocations or queries.
    /// No centroid-density approximation or thin-slab volume subtraction occurs.
    pub fn load(
        &self,
        mesh: &DisplacementMesh,
        body: &PrincipalBody<'_>,
        state: RigidBodyState,
        cancelled: impl FnMut() -> bool,
    ) -> Result<HydrostaticLoad, HydrostaticError> {
        self.load_frame(mesh, body.frame(), state, cancelled)
    }

    /// Advance the original endpoint-load Verlet integrator with this profile.
    /// Both endpoint loads and the absolute energy-defect gate share its
    /// continuous pressure potential; no second dynamics state is introduced.
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

impl HydrostaticReservoir for DensityProfileHydrostaticMedium<'_> {
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
        if self.uniform {
            return self.surface.load_frame(mesh, body, state, cancelled);
        }
        poll(&mut cancelled)?;
        let pose = body.world_from_body(state)?;
        let displaced = mesh.submerged(
            pose,
            self.surface.down_world,
            self.surface.surface_point_world,
            &mut cancelled,
        )?;
        let mut mass = 0.0;
        let mut center = Vec3::new(0.0, 0.0, 0.0);
        let mut potential = 0.0;
        for (index, sample) in self.samples.iter().enumerate() {
            poll(&mut cancelled)?;
            let next = self.samples.get(index + 1);
            let thickness = next.map(|next| next.input.depth_m - sample.input.depth_m);
            let (slab, depth) = mesh.slab_depth_moments(
                pose,
                self.surface.down_world,
                self.surface.surface_point_world,
                sample.input.depth_m,
                thickness,
                &mut cancelled,
            )?;
            let Some(slab_center) = slab.centroid_body_m() else {
                continue;
            };
            let rho = sample.transport.density_kg_m3();
            let uniform_mass = finite(rho * slab.volume_m3())?;
            let gradient_mass = match next {
                None => 0.0,
                Some(next) => finite(
                    (next.transport.density_kg_m3() - rho)
                        * (slab.depth_integral_m4 / (next.input.depth_m - sample.input.depth_m)),
                )?,
            };
            let slab_mass = finite(uniform_mass + gradient_mass)?;
            if slab_mass <= 0.0 {
                return Err(HydrostaticError::Invalid(
                    "displaced fluid mass is unrepresentable",
                ));
            }
            let mass_center = finite_vec(
                slab_center + (depth.center_body_m - slab_center).scale(gradient_mass / slab_mass),
            )?;
            let next_mass = finite(mass + slab_mass)?;
            center = finite_vec(center + (mass_center - center).scale(slab_mass / next_mass))?;
            mass = next_mass;
            potential = finite(
                potential
                    + sample.pressure_head * slab.volume_m3()
                    + rho * slab.depth_integral_m4
                    + 0.5 * gradient_mass * depth.mean_depth_m,
            )?;
        }
        let load = HydrostaticLoad::from_mass_moments(
            displaced,
            mass,
            (mass > 0.0).then_some(center),
            finite(self.surface.gravity_magnitude * potential)?,
            self.gravity(),
            body,
            state,
        )?;
        poll(&mut cancelled)?;
        Ok(load)
    }
}
