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
/// Temperature/source replacement is an externally accounted model change.
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
        poll(&mut cancelled)?;
        // Subtract origins before rotating, avoiding inversion of a large
        // world translation followed by cancellation of equally large values.
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
            poll(&mut cancelled)?;
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
        // The footprint bound guarantees positive open area at EVERY height
        // for the caller's non-overlapping tetrahedral partition. A sealing
        // piston needs a moving-boundary/pressure model, not this free surface.
        let footprint = finite((upper.x - lower.x) * (upper.y - lower.y))?;
        if finite(self.area_m2 - footprint)? <= 0.0 {
            return Err(HydrostaticError::Invalid(
                "finite tank has no representable open waterplane area",
            ));
        }
        let mut evaluations = 0;
        let mut evaluate = |height: f64| -> Result<(DisplacedVolume, f64), HydrostaticError> {
            poll(&mut cancelled)?;
            evaluations += 1;
            let surface = Vec3::new(0.0, 0.0, height);
            let down = Vec3::new(0.0, 0.0, -1.0);
            let displaced = if height <= lower.z {
                DisplacedVolume::dry()
            } else if height >= upper.z {
                mesh.whole(tank_from_body, down, surface)?
            } else {
                mesh.submerged(tank_from_body, down, surface, &mut cancelled)?
            };
            let residual =
                finite(self.area_m2 * (height - self.empty_height_m) - displaced.volume_m3())?;
            Ok((displaced, residual))
        };
        let mut low = self.empty_height_m;
        let (low_displaced, mut f_low) = evaluate(low)?;
        if f_low.abs() <= self.controls.maximum_volume_defect_m3 {
            drop(evaluate);
            return self.finish_fill(low, low_displaced, f_low, evaluations, &mut cancelled);
        }
        let unconstrained_high = finite(self.empty_height_m + mesh.volume_m3() / self.area_m2)?;
        let mut high = unconstrained_high.min(sides.z);
        let (high_displaced, mut f_high) = evaluate(high)?;
        if f_high.abs() <= self.controls.maximum_volume_defect_m3 {
            drop(evaluate);
            return self.finish_fill(high, high_displaced, f_high, evaluations, &mut cancelled);
        }
        if f_high < 0.0 {
            drop(evaluate);
            if high == sides.z {
                return Err(HydrostaticError::TankOverflow {
                    excess_volume_m3: -f_high,
                });
            }
            return Err(HydrostaticError::TankVolumeNotConverged {
                defect_m3: f_high,
                maximum_m3: self.controls.maximum_volume_defect_m3,
                evaluations,
            });
        }
        let mut best_defect = if f_low.abs() < f_high.abs() {
            f_low
        } else {
            f_high
        };
        for evaluation in 3..=self.controls.maximum_volume_evaluations {
            let span = high - low;
            let fraction = -f_low / (f_high - f_low);
            // Secants solve an upright prismatic body's linear interval in one
            // evaluation. Periodic bisection also bounds progress near tiny
            // caps where an endpoint secant can otherwise stagnate.
            let trial = if evaluation % 4 == 0 || !(0.01..0.99).contains(&fraction) {
                low + 0.5 * span
            } else {
                low + fraction * span
            };
            if trial <= low || trial >= high {
                break;
            }
            let (displaced, residual) = evaluate(trial)?;
            if residual.abs() <= self.controls.maximum_volume_defect_m3 {
                drop(evaluate);
                return self.finish_fill(trial, displaced, residual, evaluations, &mut cancelled);
            }
            if residual.abs() < best_defect.abs() {
                best_defect = residual;
            }
            if residual < 0.0 {
                low = trial;
                f_low = residual;
            } else {
                high = trial;
                f_high = residual;
            }
        }
        Err(HydrostaticError::TankVolumeNotConverged {
            defect_m3: best_defect,
            maximum_m3: self.controls.maximum_volume_defect_m3,
            evaluations,
        })
    }

    fn finish_fill(
        &self,
        height: f64,
        displaced: DisplacedVolume,
        residual: f64,
        evaluations: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FiniteTankFill, HydrostaticError> {
        let volume = displaced.volume_m3();
        let geometric_potential =
            finite(displaced.depth_integral_m4 - 0.5 * volume * (volume / self.area_m2))?;
        if geometric_potential < 0.0 {
            return Err(HydrostaticError::Invalid(
                "finite tank pressure potential is unrepresentable",
            ));
        }
        let liquid_gravitational_energy_change_j =
            finite(self.transport.density_kg_m3() * self.gravity_magnitude * geometric_potential)?;
        let surface_point_world_m = finite_vec(
            self.geometry
                .world_from_tank
                .transform_point(Vec3::new(0.0, 0.0, height))?,
        )?;
        poll(cancelled)?;
        Ok(FiniteTankFill {
            surface_height_m: height,
            surface_point_world_m,
            displaced,
            liquid_volume_defect_m3: residual,
            volume_evaluations: evaluations,
            liquid_gravitational_energy_change_j,
        })
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
