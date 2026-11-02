//! A finite-volume rectangular bath on the existing clipped hydrostatic owner.
//!
//! The free surface is an unknown: A h - V_sub(h,q) = V_liquid. The original
//! tetrahedral clipping supplies both displaced volume and its first moment.
//! Relative to the empty bath at h0=V_liquid/A, the liquid gravitational
//! potential is rho*g*(integral_sub(h-z)dV - V_sub^2/(2A)). Its negative pose
//! derivative is the SAME Archimedes wrench at the solved waterline.
//!
//! This is instantaneous hydrostatic equilibration of one homogeneous liquid,
//! not a resolved sloshing/added-mass model. The stationary rectangular tank
//! has vertical walls and a horizontal floor. The complete displacement mesh
//! must stay strictly inside its lateral footprint and above its floor; the
//! body may protrude above the rim. Overflow, contact and a sealing piston
//! refuse. DisplacementMesh's non-overlapping partition contract still applies.

use super::{
    DisplacedVolume, DisplacementMesh, HydrostaticError, HydrostaticLoad, HydrostaticReservoir,
    HydrostaticStep, HydrostaticStepControls, HydrostaticStepError, finite, finite_vec,
    gravity_direction, poll,
};
use crate::principal::{PrincipalBody, PrincipalFrame};
use fs_ga::{Se3, Vec3};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_mbd::{Gravity, RigidBodyState};

mod bodies;
mod density_change;
mod solve;
pub use bodies::{FiniteTankBodiesLoad, FiniteTankBodiesStep, FiniteTankBody, MAX_TANK_BODIES};
pub use density_change::FiniteTankDensityChange;
pub use solve::FiniteTankSurface;

/// Hard bound on geometry evaluations per finite-bath load query.
pub const MAX_TANK_VOLUME_EVALUATIONS: usize = 128;

/// Rectangular tank interior. Local z=0 is the floor; local +z points up.
/// Local x/y are centered on the floor, with the stated full inside lengths.
#[derive(Clone, Copy, Debug)]
pub struct HydrostaticTankGeometry {
    /// Fixed tank frame in world metres.
    pub world_from_tank: Se3,
    /// Positive inside width, depth and floor-to-rim height [m].
    pub inside_dimensions_m: Vec3,
}

/// Root-solve work and conservation allowance, independent of dynamics dt.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteTankControls {
    /// Positive absolute liquid-volume residual allowance [m3].
    /// This numerical allowance is not a liquid source or overflow model.
    pub maximum_volume_defect_m3: f64,
    /// Includes both bracket endpoints; must lie in 2..=128.
    pub maximum_volume_evaluations: usize,
}

/// The solved bath level and original analytic displacement at one body pose.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteTankFill {
    /// Height of the horizontal free surface above the tank floor [m].
    pub surface_height_m: f64,
    /// Corresponding point in the world frame [m].
    pub surface_point_world_m: Vec3,
    /// Original body-frame displacement integral and centroid.
    pub displaced: DisplacedVolume,
    /// A*h - V_sub - V_liquid [m3], evaluated without subtracting large totals.
    pub liquid_volume_defect_m3: f64,
    /// Actual bounded geometry evaluations, including bracket endpoints.
    pub volume_evaluations: usize,
    /// Liquid gravitational energy relative to the same liquid in an empty
    /// tank [J]. Do not also add an infinite-reservoir pressure potential.
    pub liquid_gravitational_energy_change_j: f64,
}

/// One frozen source snapshot and one conserved liquid inventory.
/// `rebind_source_at_pose` carries this same mass to another supplied liquid
/// state and reports the complete gravitational-energy change at a held pose.
#[derive(Debug)]
pub struct FiniteHydrostaticTank<'a, F: FluidTransport + ?Sized> {
    source: &'a F,
    transport: FluidTransportState,
    gravity: Gravity,
    gravity_magnitude: f64,
    geometry: HydrostaticTankGeometry,
    area_m2: f64,
    liquid_volume_m3: f64,
    liquid_mass_kg: f64,
    empty_height_m: f64,
    controls: FiniteTankControls,
}

impl<'a, F: FluidTransport + ?Sized> FiniteHydrostaticTank<'a, F> {
    /// Bind a homogeneous, constant-density liquid to a fixed finite tank.
    /// The original fluid transport source is evaluated once; no body density,
    /// default water properties or buoyancy coefficient is substituted.
    pub fn new(
        source: &'a F,
        gravity: Gravity,
        geometry: HydrostaticTankGeometry,
        liquid_volume_m3: f64,
        controls: FiniteTankControls,
    ) -> Result<Self, HydrostaticError> {
        let sides = finite_vec(geometry.inside_dimensions_m)?;
        finite_vec(geometry.world_from_tank.translation())?;
        if sides.x <= 0.0
            || sides.y <= 0.0
            || sides.z <= 0.0
            || !liquid_volume_m3.is_finite()
            || liquid_volume_m3 <= 0.0
            || !controls.maximum_volume_defect_m3.is_finite()
            || controls.maximum_volume_defect_m3 <= 0.0
            || !(2..=MAX_TANK_VOLUME_EVALUATIONS).contains(&controls.maximum_volume_evaluations)
        {
            return Err(HydrostaticError::Invalid(
                "finite tank requires positive dimensions, inventory and volume tolerance, and bounded evaluations",
            ));
        }
        let (down, gravity_magnitude) = gravity_direction(gravity)?;
        let up = geometry
            .world_from_tank
            .rotation()
            .rotate(Vec3::new(0.0, 0.0, 1.0))?;
        let misalignment = finite_vec(up + down)?;
        if misalignment.dot(misalignment) > 1e-24 {
            return Err(HydrostaticError::Invalid(
                "finite tank local +z must oppose gravity; tilted floors are unsupported",
            ));
        }
        let area_m2 = finite(sides.x * sides.y)?;
        let capacity_m3 = finite(area_m2 * sides.z)?;
        if area_m2 <= 0.0 || liquid_volume_m3 > capacity_m3 {
            return Err(HydrostaticError::Invalid(
                "liquid inventory exceeds empty tank capacity",
            ));
        }
        let empty_height_m = finite(liquid_volume_m3 / area_m2)?;
        if empty_height_m <= 0.0 {
            return Err(HydrostaticError::Invalid(
                "empty tank liquid height is unrepresentable",
            ));
        }
        let transport = source.transport_state()?;
        let liquid_mass_kg = finite(transport.density_kg_m3() * liquid_volume_m3)?;
        if liquid_mass_kg <= 0.0 {
            return Err(HydrostaticError::Invalid(
                "finite tank liquid mass is unrepresentable",
            ));
        }
        Ok(Self {
            source,
            transport,
            gravity,
            gravity_magnitude,
            geometry,
            area_m2,
            liquid_volume_m3,
            liquid_mass_kg,
            empty_height_m,
            controls,
        })
    }

    /// Original immutable source, with its context and provenance intact.
    pub const fn source(&self) -> &'a F {
        self.source
    }
    /// Exact source snapshot used by every load and energy calculation.
    pub const fn transport(&self) -> FluidTransportState {
        self.transport
    }
    /// Fixed tank geometry; moving the tank is outside this model.
    pub const fn geometry(&self) -> HydrostaticTankGeometry {
        self.geometry
    }
    /// Constant liquid inventory [m3].
    pub const fn liquid_volume_m3(&self) -> f64 {
        self.liquid_volume_m3
    }
    /// Constant source-derived liquid mass [kg].
    pub const fn liquid_mass_kg(&self) -> f64 {
        self.liquid_mass_kg
    }
    /// Body gravity to pass exactly once to the original rigid-body owner.
    pub const fn gravity(&self) -> Gravity {
        self.gravity
    }

    /// Solve the free surface at the actual ORIGINAL body-frame pose.
    /// A safeguarded secant uses the exact original clipped volume, never a
    /// sampled submerged fraction. Dry/full cases reuse cached mesh moments.
    /// Every geometry pass polls cancellation once per tetrahedron.
    pub fn fill_at_pose(
        &self,
        mesh: &DisplacementMesh,
        world_from_body: Se3,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FiniteTankFill, HydrostaticError> {
        let prepared = self.prepare_displacement(mesh, world_from_body, &mut cancelled)?;
        self.check_open_area(prepared.lower, prepared.upper)?;
        let mut displaced = DisplacedVolume::dry();
        let surface = self.solve_surface(mesh.volume_m3(), &mut cancelled, |height, cancel| {
            displaced = prepared.at_height(height, cancel)?;
            Ok((displaced.volume_m3(), displaced.depth_integral_m4))
        })?;
        Ok(FiniteTankFill {
            surface_height_m: surface.surface_height_m,
            surface_point_world_m: surface.surface_point_world_m,
            displaced,
            liquid_volume_defect_m3: surface.liquid_volume_defect_m3,
            volume_evaluations: surface.volume_evaluations,
            liquid_gravitational_energy_change_j: surface.liquid_gravitational_energy_change_j,
        })
    }

    fn prepare_displacement<'m>(
        &self,
        mesh: &'m DisplacementMesh,
        world_from_body: Se3,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedDisplacement<'m>, HydrostaticError> {
        poll(cancelled)?;
        // Subtract origins before rotating to retain small tank-relative offsets.
        let inverse_rotation = self.geometry.world_from_tank.rotation().inverse();
        let tank_from_body = Se3::from_parts(
            inverse_rotation.compose(world_from_body.rotation())?,
            inverse_rotation.rotate(
                world_from_body.translation() - self.geometry.world_from_tank.translation(),
            )?,
        )?;
        let sides = self.geometry.inside_dimensions_m;
        let mut lower = Vec3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut upper = Vec3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        for cell in mesh.cells() {
            poll(cancelled)?;
            for &point in cell {
                let p = finite_vec(tank_from_body.transform_point(point)?)?;
                if p.x <= -0.5 * sides.x
                    || p.x >= 0.5 * sides.x
                    || p.y <= -0.5 * sides.y
                    || p.y >= 0.5 * sides.y
                    || p.z <= 0.0
                {
                    return Err(HydrostaticError::Invalid(
                        "displacement envelope touches or crosses the tank floor or side walls",
                    ));
                }
                lower = Vec3::new(lower.x.min(p.x), lower.y.min(p.y), lower.z.min(p.z));
                upper = Vec3::new(upper.x.max(p.x), upper.y.max(p.y), upper.z.max(p.z));
            }
        }
        Ok(PreparedDisplacement {
            mesh,
            tank_from_body,
            lower,
            upper,
        })
    }

    fn check_open_area(&self, lower: Vec3, upper: Vec3) -> Result<(), HydrostaticError> {
        // The complete envelope footprint leaves an open waterplane at every
        // height. A sealing piston requires a different pressure model.
        let footprint = finite((upper.x - lower.x) * (upper.y - lower.y))?;
        if finite(self.area_m2 - footprint)? <= 0.0 {
            return Err(HydrostaticError::Invalid(
                "finite tank has no representable open waterplane area",
            ));
        }
        Ok(())
    }

    /// Buoyancy at the conserved-inventory waterline, excluding body weight.
    pub fn load(
        &self,
        mesh: &DisplacementMesh,
        body: &PrincipalBody<'_>,
        state: RigidBodyState,
        cancelled: impl FnMut() -> bool,
    ) -> Result<HydrostaticLoad, HydrostaticError> {
        self.load_frame(mesh, body.frame(), state, cancelled)
    }

    /// Original rigid-body Verlet mechanics with both waterlines recomputed
    /// from the same liquid inventory. The existing absolute energy-defect
    /// gate now measures body energy PLUS this finite liquid's potential.
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

impl<F: FluidTransport + ?Sized> HydrostaticReservoir for FiniteHydrostaticTank<'_, F> {
    fn gravity(&self) -> Gravity {
        self.gravity
    }

    fn load_frame(
        &self,
        mesh: &DisplacementMesh,
        body: &PrincipalFrame,
        state: RigidBodyState,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HydrostaticLoad, HydrostaticError> {
        poll(&mut cancelled)?;
        let fill = self.fill_at_pose(mesh, body.world_from_body(state)?, &mut cancelled)?;
        let mass = finite(self.transport.density_kg_m3() * fill.displaced.volume_m3())?;
        if fill.displaced.volume_m3() > 0.0 && mass <= 0.0 {
            return Err(HydrostaticError::Invalid(
                "displaced finite tank mass is unrepresentable",
            ));
        }
        let load = HydrostaticLoad::from_mass_moments(
            fill.displaced,
            mass,
            fill.displaced.centroid_body_m(),
            fill.liquid_gravitational_energy_change_j,
            self.gravity,
            body,
            state,
        )?;
        poll(&mut cancelled)?;
        Ok(load)
    }
}

/// Prepared original mesh/frame, retaining the native dry/full clipping paths.
struct PreparedDisplacement<'a> {
    mesh: &'a DisplacementMesh,
    tank_from_body: Se3,
    lower: Vec3,
    upper: Vec3,
}
impl PreparedDisplacement<'_> {
    fn at_height(
        &self,
        height: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DisplacedVolume, HydrostaticError> {
        let surface = Vec3::new(0.0, 0.0, height);
        let down = Vec3::new(0.0, 0.0, -1.0);
        if height <= self.lower.z {
            Ok(DisplacedVolume::dry())
        } else if height >= self.upper.z {
            self.mesh.whole(self.tank_from_body, down, surface)
        } else {
            self.mesh
                .submerged(self.tank_from_body, down, surface, cancelled)
        }
    }
}
