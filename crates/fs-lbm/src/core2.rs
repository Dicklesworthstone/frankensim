//! The general D2Q9 core (bead tfz.19): cell flags, VECTOR gravity
//! (tilt schedules rotate it), per-cell relaxation time (the
//! non-Newtonian hook), Guo forcing, pull streaming with halfway
//! bounce-back at walls — the substrate the thermal, rheology,
//! refinement, and free-surface extensions all share. Deterministic:
//! fixed row-major cell order, no RNG.

use crate::{CS2, E, OPP, Q, W};
use std::f64::consts::{PI, SQRT_2};

const MAX_REGULARIZED_BOUNDARY_SPEED_SQ: f64 = 0.03;
const UNIT_NORMAL_TOLERANCE: f64 = 128.0 * f64::EPSILON;

/// Collision operator for the D2Q9 core (naming follows the d3q19
/// [`crate::d3q19::CollisionModel3`] precedent; tau stays per-cell,
/// so the variants carry no parameters).
///
/// Default is [`CollisionModel2::Bgk`], whose arithmetic path is
/// byte-for-byte the historical one — every existing pin and golden
/// is unchanged. [`CollisionModel2::Regularized`] is the
/// Latt–Chopard regularized collision: the non-equilibrium part is
/// projected onto the second-order Hermite subspace
/// (`f_neq_reg = w_i (9/2) Q_i : Pi_neq`, `Q_i = e_i e_i - cs^2 I`),
/// discarding the non-hydrodynamic ghost content that drives the
/// plain-BGK instability at high cell Reynolds number (executed
/// boundary on the jet-labium family: `Re/delta` in (36, 48)). The
/// shear viscosity keeps the identical `nu = (tau - 1/2)/3` law.
///
/// Scope: honored by the plain [`Grid::collide_into`] path (and
/// everything built on it, including the stepping drivers and the
/// partial-saturation proposal). Combining `Regularized` with the
/// partial-saturation solid operator or body forces is possible but
/// UNVALIDATED here: with Guo forcing the bare non-equilibrium is
/// regularized without subtracting the half-force moment — exact
/// for unforced cells, a disclosed approximation otherwise.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum CollisionModel2 {
    /// Plain BGK single-relaxation-time collision (historical path).
    #[default]
    Bgk,
    /// Latt–Chopard regularized collision (second-order Hermite
    /// projection of the non-equilibrium). EXECUTED BOUNDARIES: on
    /// the jet-labium rig it roughly doubles the stable window
    /// (Re_cell (36,48) -> (48,96)) but its plate-force spectrum is
    /// dominated by a temporal-Nyquist parity artifact with
    /// bounce-back walls, and on an UNRESOLVED (about one cell)
    /// free shear layer it destabilizes where plain BGK survives.
    Regularized,
    /// Recursive regularized collision (Coreixas-style RR3): as
    /// [`CollisionModel2::Regularized`], plus the two D2Q9-supported
    /// third-order Hermite coefficients reconstructed RECURSIVELY
    /// from the second-order ones
    /// (`a3_xxy = 2 ux a2_xy + uy a2_xx`,
    /// `a3_xyy = 2 uy a2_xy + ux a2_yy`), which repairs the
    /// projected operator's high-shear fragility.
    RecursiveRegularized,
    /// Two-relaxation-time collision (Ginzburg/d'Humieres): the
    /// even (symmetric) non-equilibrium relaxes at the viscous rate
    /// `1/tau` (so `nu = (tau - 1/2)/3` exactly as BGK), the odd
    /// (antisymmetric) part at the rate fixed by the magic
    /// combination `magic = (tau_even - 1/2)(tau_odd - 1/2)`.
    /// `magic = 3/16` places straight bounce-back walls exactly
    /// halfway; `Trt { magic: (tau - 1/2)^2 }` is analytically BGK.
    /// EXECUTED WARNING (2026-08-10 scan, this crate's example
    /// driver): the folklore "magic = 1/4 is stability-optimal" is
    /// FALSE on impulsively started shear layers at tau near 1/2 —
    /// every tested magic above the BGK point REDUCED the stable
    /// window on both the periodic box and the jet-labium rig
    /// (tau_odd = 1/2 + magic/(tau - 1/2) grows huge, freezing odd
    /// modes). Use [`CollisionModel2::CentralMoment`] when the goal
    /// is stability.
    Trt {
        /// The magic parameter `(tau_even - 1/2)(tau_odd - 1/2)`;
        /// must be finite and positive.
        magic: f64,
    },
    /// Central-moment (cascaded, Geier/De Rosis class) collision,
    /// the D2Q9 analog of the d3q19 `CollisionModel3::CentralMoment`
    /// precedent: moments are shifted into the local velocity frame,
    /// the shear pair `{k11, k20 - k02}` relaxes at the viscous rate
    /// `1/tau` (so `nu = (tau - 1/2)/3` exactly), and the bulk
    /// (`k20 + k02`), third (`k21`, `k12`), and fourth (`k22`)
    /// central moments are EQUILIBRATED each step to their
    /// Maxwellian values — the maximally dissipative ghost/bulk
    /// choice. Note the enhanced bulk dissipation damps in-domain
    /// acoustics harder than BGK; for the Curle rigs (surface-force
    /// sources, sponge-absorbed domain waves) that is a stability
    /// feature, but it is NOT a faithful acoustic-propagation
    /// operator.
    CentralMoment,
}

/// Cell classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Cell {
    /// Bulk fluid.
    Fluid,
    /// Solid wall (halfway bounce-back).
    Wall,
    /// Free-surface interface cell (carries partial mass).
    Interface,
    /// Gas cell (no populations).
    Gas,
}

/// The general D2Q9 lattice.
#[derive(Debug, Clone)]
pub struct Grid {
    /// Cells in x.
    pub nx: usize,
    /// Cells in y.
    pub ny: usize,
    /// Cell flags.
    pub flags: Vec<Cell>,
    /// Populations.
    pub f: Vec<[f64; Q]>,
    /// Per-cell relaxation time.
    pub tau: Vec<f64>,
    /// Gravity vector (lattice units).
    pub g: [f64; 2],
    /// Per-cell external force (Boussinesq buoyancy etc.), added to
    /// ρ·g.
    pub fext: Vec<[f64; 2]>,
    /// Periodic in x?
    pub periodic_x: bool,
    /// Periodic in y?
    pub periodic_y: bool,
    /// Collision operator (default [`CollisionModel2::Bgk`], the
    /// byte-identical historical path).
    pub collision: CollisionModel2,
}

/// Macroscopic moments of one cell.
#[derive(Debug, Clone, Copy)]
pub struct Moments {
    /// Density.
    pub rho: f64,
    /// Velocity (force-corrected).
    pub u: [f64; 2],
}

/// Momentum delivered to a selected set of stationary halfway-bounce-back
/// wall cells during one D2Q9 stream.
///
/// `wall_impulse` uses lattice momentum units. With the lattice time step
/// equal to one it is also the raw lattice force; callers remain responsible
/// for physical-unit and drag/lift normalization.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
#[must_use]
pub struct MomentumExchange2 {
    /// Net `(x, y)` impulse delivered by the fluid to selected walls.
    pub wall_impulse: [f64; 2],
    /// Number of selected fluid-wall links included in the sum.
    pub measured_links: usize,
}

/// Momentum, torque, and work delivered to selected moving D2Q9 walls.
///
/// `wall_impulse` uses the boundary-relative momentum-exchange expression
/// `(c_out - u_wall) f_out - (c_in - u_wall) f_in`. The separate population
/// and wall-velocity mass terms retain the discrete balance
///
/// `wall_impulse + fluid_population_impulse = wall_velocity_mass_impulse`
///
/// up to floating-point roundoff. All values are raw lattice quantities for
/// one stream step; no physical-unit or reference-area normalization is
/// applied. The relative-velocity force convention follows Wen et al.,
/// *Galilean Invariant Fluid-Solid Interfacial Dynamics in Lattice Boltzmann
/// Simulations* (2014; <https://arxiv.org/abs/1303.0625>).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
#[must_use]
pub struct MovingWallMomentumExchange2 {
    /// Boundary-relative hydrodynamic impulse delivered to selected walls.
    pub wall_impulse: [f64; 2],
    /// Change in resolved fluid-population momentum across selected links.
    pub fluid_population_impulse: [f64; 2],
    /// Net incoming-minus-outgoing population mass across selected links.
    pub fluid_mass_change: f64,
    /// Sum of `u_wall * (f_in - f_out)` across selected links.
    pub wall_velocity_mass_impulse: [f64; 2],
    /// Scalar 2-D wall angular impulse about the requested origin.
    pub wall_angular_impulse: f64,
    /// Scalar angular impulse of the resolved fluid-population change.
    pub fluid_population_angular_impulse: f64,
    /// Scalar angular impulse of `wall_velocity_mass_impulse`.
    pub wall_velocity_mass_angular_impulse: f64,
    /// Work delivered to selected walls, `sum(wall_impulse_link dot u_wall)`.
    pub wall_work: f64,
    /// Number of selected fluid-wall links included in the receipt.
    pub measured_links: usize,
}

/// Caller-supplied solid coverage and velocity for one partially saturated
/// D2Q9 cell.
///
/// `solid_fraction` is the solid volume fraction `epsilon_s` in `[0, 1]`.
/// The zero-fraction endpoint requires exactly zero solid velocity so an
/// inactive entry cannot carry latent motion into a later interpretation.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct PartialSaturationCell2 {
    solid_fraction: f64,
    solid_velocity: [f64; 2],
}

impl PartialSaturationCell2 {
    /// Construct a checked partial-saturation cell descriptor.
    #[must_use]
    pub fn new(solid_fraction: f64, solid_velocity: [f64; 2]) -> Self {
        assert!(
            solid_fraction.is_finite() && (0.0..=1.0).contains(&solid_fraction),
            "D2Q9 partial-saturation solid fraction must be finite and in [0, 1]"
        );
        assert!(
            solid_velocity.into_iter().all(f64::is_finite),
            "D2Q9 partial-saturation solid velocity must be finite"
        );
        let speed_sq =
            solid_velocity[0].mul_add(solid_velocity[0], solid_velocity[1] * solid_velocity[1]);
        assert!(
            speed_sq < MAX_REGULARIZED_BOUNDARY_SPEED_SQ,
            "D2Q9 partial-saturation solid velocity exceeds the low-Mach admission envelope"
        );
        assert!(
            solid_fraction > 0.0 || solid_velocity == [0.0; 2],
            "zero D2Q9 solid fraction requires zero solid velocity"
        );
        Self {
            solid_fraction,
            solid_velocity,
        }
    }

    /// Solid volume fraction in this lattice cell.
    #[must_use]
    pub const fn solid_fraction(self) -> f64 {
        self.solid_fraction
    }

    /// Prescribed local solid velocity in lattice units.
    #[must_use]
    pub const fn solid_velocity(self) -> [f64; 2] {
        self.solid_velocity
    }

    fn coupling_weight(self, tau: f64) -> f64 {
        let relaxation_margin = tau - 0.5;
        self.solid_fraction * relaxation_margin / ((1.0 - self.solid_fraction) + relaxation_margin)
    }
}

/// Aggregate exchange generated by one partially saturated D2Q9 collision.
///
/// The fluid terms are the population mass, momentum, and angular-momentum
/// changes contributed by the solid-coupling correction relative to ordinary
/// force-free BGK collision. Solid impulse and angular impulse are their
/// equal-and-opposite counterparts. All values are raw lattice quantities for
/// one collision; physical-unit conversion and body attribution remain with
/// the caller.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
#[must_use]
pub struct PartialSaturationExchange2 {
    /// Population mass change contributed by solid coupling (roundoff only).
    pub fluid_population_mass_change: f64,
    /// Population momentum change contributed by solid coupling.
    pub fluid_population_impulse: [f64; 2],
    /// Equal-and-opposite hydrodynamic impulse delivered to the solid field.
    pub solid_impulse: [f64; 2],
    /// Fluid angular impulse about the caller's requested origin.
    pub fluid_population_angular_impulse: f64,
    /// Equal-and-opposite solid angular impulse about that origin.
    pub solid_angular_impulse: f64,
    /// Work delivered to the solid field, `sum(solid_impulse_cell dot u_s)`.
    pub solid_work: f64,
    /// Number of cells with strictly positive solid fraction.
    pub coupled_cells: usize,
    /// Sum of the relaxation-dependent Noble-Torczynski coupling weights.
    pub coupling_weight_sum: f64,
}

/// Checked local geometry and motion for one under-resolved pair of circular
/// cylinders in a D2Q9 plane.
///
/// The correction is per unit out-of-plane depth and acts only along the
/// caller-supplied common normal. `reduced_radius` is
/// `(1 / radius_first + 1 / radius_second)^-1`; a planar second surface is
/// represented by the first cylinder's radius. Geometry discovery, pair
/// enumeration, and rigid-body integration remain caller responsibilities.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UnderResolvedCylinderPair2 {
    gap: f64,
    cutoff_gap: f64,
    minimum_gap: f64,
    reduced_radius: f64,
    dynamic_viscosity: f64,
    normal_first_to_second: [f64; 2],
    line_of_action_point: [f64; 2],
    surface_velocities: [[f64; 2]; 2],
}

impl UnderResolvedCylinderPair2 {
    /// Construct a checked pair-local lubrication descriptor.
    #[must_use]
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        gap: f64,
        cutoff_gap: f64,
        minimum_gap: f64,
        reduced_radius: f64,
        dynamic_viscosity: f64,
        normal_first_to_second: [f64; 2],
        line_of_action_point: [f64; 2],
        surface_velocities: [[f64; 2]; 2],
    ) -> Self {
        assert!(
            gap.is_finite() && gap >= 0.0,
            "D2Q9 lubrication gap must be finite and nonnegative"
        );
        assert!(
            cutoff_gap.is_finite() && cutoff_gap > 0.0,
            "D2Q9 lubrication cutoff gap must be positive and finite"
        );
        assert!(
            minimum_gap.is_finite() && minimum_gap > 0.0,
            "D2Q9 lubrication minimum gap must be positive and finite"
        );
        assert!(
            minimum_gap < cutoff_gap,
            "D2Q9 lubrication minimum gap must be below the cutoff gap"
        );
        assert!(
            reduced_radius.is_finite() && reduced_radius > 0.0,
            "D2Q9 lubrication reduced radius must be positive and finite"
        );
        assert!(
            dynamic_viscosity.is_finite() && dynamic_viscosity > 0.0,
            "D2Q9 lubrication dynamic viscosity must be positive and finite"
        );
        assert!(
            normal_first_to_second.into_iter().all(f64::is_finite),
            "D2Q9 lubrication normal must be finite"
        );
        let normal_norm = normal_first_to_second[0]
            .mul_add(
                normal_first_to_second[0],
                normal_first_to_second[1] * normal_first_to_second[1],
            )
            .sqrt();
        assert!(
            (normal_norm - 1.0).abs() <= UNIT_NORMAL_TOLERANCE,
            "D2Q9 lubrication normal must have unit length"
        );
        assert!(
            line_of_action_point.into_iter().all(f64::is_finite),
            "D2Q9 lubrication line-of-action point must be finite"
        );
        for surface_velocity in surface_velocities {
            assert!(
                surface_velocity.into_iter().all(f64::is_finite),
                "D2Q9 lubrication surface velocity must be finite"
            );
            let speed_sq = surface_velocity[0].mul_add(
                surface_velocity[0],
                surface_velocity[1] * surface_velocity[1],
            );
            assert!(
                speed_sq < MAX_REGULARIZED_BOUNDARY_SPEED_SQ,
                "D2Q9 lubrication surface velocity exceeds the low-Mach admission envelope"
            );
        }
        Self {
            gap,
            cutoff_gap,
            minimum_gap,
            reduced_radius,
            dynamic_viscosity,
            normal_first_to_second,
            line_of_action_point,
            surface_velocities,
        }
    }

    /// Surface-to-surface gap in lattice units.
    #[must_use]
    pub const fn gap(self) -> f64 {
        self.gap
    }

    /// Gap at which the explicit unresolved correction vanishes.
    #[must_use]
    pub const fn cutoff_gap(self) -> f64 {
        self.cutoff_gap
    }

    /// Positive gap floor used to regularize the singular resistance.
    #[must_use]
    pub const fn minimum_gap(self) -> f64 {
        self.minimum_gap
    }

    /// Pair reduced radius in lattice units.
    #[must_use]
    pub const fn reduced_radius(self) -> f64 {
        self.reduced_radius
    }

    /// Dynamic viscosity in raw lattice units.
    #[must_use]
    pub const fn dynamic_viscosity(self) -> f64 {
        self.dynamic_viscosity
    }

    /// Unit normal directed from the first surface to the second.
    #[must_use]
    pub const fn normal_first_to_second(self) -> [f64; 2] {
        self.normal_first_to_second
    }

    /// Point on the declared common normal used for both angular impulses.
    #[must_use]
    pub const fn line_of_action_point(self) -> [f64; 2] {
        self.line_of_action_point
    }

    /// First and second surface velocities at the line-of-action point.
    #[must_use]
    pub const fn surface_velocities(self) -> [[f64; 2]; 2] {
        self.surface_velocities
    }

    /// Evaluate the pair-local correction about `moment_origin`.
    #[must_use]
    #[allow(clippy::too_many_lines)]
    pub fn correction(self, moment_origin: [f64; 2]) -> UnderResolvedLubricationExchange2 {
        assert!(
            moment_origin.into_iter().all(f64::is_finite),
            "D2Q9 lubrication moment origin must be finite"
        );

        let effective_gap = self.gap.max(self.minimum_gap);
        let within_cutoff = self.gap < self.cutoff_gap;
        let resistance_increment = if within_cutoff {
            let radius_three_halves = self.reduced_radius * self.reduced_radius.sqrt();
            let effective_gap_inverse_three_halves = 1.0 / (effective_gap * effective_gap.sqrt());
            let cutoff_inverse_three_halves = 1.0 / (self.cutoff_gap * self.cutoff_gap.sqrt());
            3.0 * PI
                * SQRT_2
                * self.dynamic_viscosity
                * radius_three_halves
                * (effective_gap_inverse_three_halves - cutoff_inverse_three_halves)
        } else {
            0.0
        };
        assert!(
            resistance_increment.is_finite() && resistance_increment >= 0.0,
            "D2Q9 lubrication resistance increment must be finite and nonnegative"
        );

        let relative_velocity = [
            self.surface_velocities[1][0] - self.surface_velocities[0][0],
            self.surface_velocities[1][1] - self.surface_velocities[0][1],
        ];
        let normal_gap_velocity = relative_velocity[0].mul_add(
            self.normal_first_to_second[0],
            relative_velocity[1] * self.normal_first_to_second[1],
        );
        assert!(
            normal_gap_velocity.is_finite(),
            "D2Q9 lubrication normal gap velocity must be finite"
        );
        let impulse_scale = resistance_increment * normal_gap_velocity;
        assert!(
            impulse_scale.is_finite(),
            "D2Q9 lubrication impulse scale must be finite"
        );
        let first_solid_impulse = if impulse_scale == 0.0 {
            [0.0; 2]
        } else {
            [
                impulse_scale * self.normal_first_to_second[0],
                impulse_scale * self.normal_first_to_second[1],
            ]
        };
        assert!(
            first_solid_impulse.into_iter().all(f64::is_finite),
            "D2Q9 lubrication solid impulse must be finite"
        );
        let second_solid_impulse = if first_solid_impulse == [0.0; 2] {
            [0.0; 2]
        } else {
            [-first_solid_impulse[0], -first_solid_impulse[1]]
        };
        let offset = [
            self.line_of_action_point[0] - moment_origin[0],
            self.line_of_action_point[1] - moment_origin[1],
        ];
        assert!(
            offset.into_iter().all(f64::is_finite),
            "D2Q9 lubrication line-of-action offset must be finite"
        );
        let first_solid_angular_impulse = offset[0].mul_add(
            first_solid_impulse[1],
            -(offset[1] * first_solid_impulse[0]),
        );
        assert!(
            first_solid_angular_impulse.is_finite(),
            "D2Q9 lubrication angular impulse must be finite"
        );
        let second_solid_angular_impulse = if first_solid_angular_impulse == 0.0 {
            0.0
        } else {
            -first_solid_angular_impulse
        };
        let solid_pair_work = if impulse_scale == 0.0 {
            0.0
        } else {
            -(impulse_scale * normal_gap_velocity)
        };
        assert!(
            solid_pair_work.is_finite() && solid_pair_work <= 0.0,
            "D2Q9 lubrication pair work must be finite and nonpositive"
        );

        UnderResolvedLubricationExchange2 {
            first_solid_impulse,
            second_solid_impulse,
            first_solid_angular_impulse,
            second_solid_angular_impulse,
            solid_pair_work,
            normal_gap_velocity,
            resistance_increment,
            effective_gap,
            within_cutoff,
            gap_floor_applied: self.gap < self.minimum_gap,
        }
    }
}

/// Equal-and-opposite impulse receipt for one local D2Q9 lubrication pair.
///
/// All values are raw lattice quantities for one lattice time step and per
/// unit out-of-plane depth. The receipt contains no resolved-fluid population
/// term because it is an explicit correction to resistance missing below the
/// declared cutoff.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
#[must_use]
pub struct UnderResolvedLubricationExchange2 {
    /// Correction impulse delivered to the first solid.
    pub first_solid_impulse: [f64; 2],
    /// Exact equal-and-opposite correction impulse delivered to the second.
    pub second_solid_impulse: [f64; 2],
    /// First-solid angular impulse about the requested origin.
    pub first_solid_angular_impulse: f64,
    /// Exact equal-and-opposite second-solid angular impulse.
    pub second_solid_angular_impulse: f64,
    /// Pair work, `-resistance_increment * normal_gap_velocity^2`.
    pub solid_pair_work: f64,
    /// Signed gap rate `(u_second - u_first) dot normal_first_to_second`.
    pub normal_gap_velocity: f64,
    /// Unresolved normal resistance after subtracting its cutoff value.
    pub resistance_increment: f64,
    /// Gap used in the singular resistance after applying the positive floor.
    pub effective_gap: f64,
    /// Whether the supplied gap is strictly below the correction cutoff.
    pub within_cutoff: bool,
    /// Whether `minimum_gap` replaced the supplied gap.
    pub gap_floor_applied: bool,
}

/// Geometry and velocity for one off-lattice D2Q9 fluid-to-wall link.
///
/// `fraction` is the distance from the destination fluid-cell center to the
/// wall intersection divided by the full lattice-link length. A value of
/// `0.5` is the ordinary halfway wall. The velocity is evaluated at the wall
/// intersection rather than at the center of a stair-step wall cell.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurvedMovingWallLink2 {
    fraction: f64,
    wall_velocity: [f64; 2],
}

impl CurvedMovingWallLink2 {
    /// Construct a checked off-lattice wall-link description.
    #[must_use]
    pub fn new(fraction: f64, wall_velocity: [f64; 2]) -> Self {
        assert!(
            fraction.is_finite() && fraction > 0.0 && fraction <= 1.0,
            "D2Q9 curved-wall link fraction must be finite and in (0, 1]"
        );
        assert!(
            wall_velocity.into_iter().all(f64::is_finite),
            "D2Q9 curved-wall velocity must be finite"
        );
        let speed_sq =
            wall_velocity[0].mul_add(wall_velocity[0], wall_velocity[1] * wall_velocity[1]);
        assert!(
            speed_sq < MAX_REGULARIZED_BOUNDARY_SPEED_SQ,
            "D2Q9 curved-wall velocity exceeds the low-Mach admission envelope"
        );
        Self {
            fraction,
            wall_velocity,
        }
    }

    /// Fluid-center-to-wall distance as a fraction of the lattice link.
    #[must_use]
    pub const fn fraction(self) -> f64 {
        self.fraction
    }

    /// Wall velocity at the link intersection in lattice units.
    #[must_use]
    pub const fn wall_velocity(self) -> [f64; 2] {
        self.wall_velocity
    }
}

/// Atomic receipt for a D2Q9 fluid/wall cell-topology transition.
///
/// Newly uncovered fluid cells are initialized by an equal-weight average of
/// unique, surviving one-ring fluid donors from the pre-transition state.
/// Covered-cell removal and fresh-cell insertion are reported separately so a
/// caller can ledger the exact active-fluid mass and momentum delta.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
#[must_use]
pub struct WallTopologyTransition2 {
    /// Number of pre-transition fluid cells changed to wall cells.
    pub covered_fluid_cells: usize,
    /// Number of pre-transition wall cells initialized as fresh fluid.
    pub fresh_fluid_cells: usize,
    /// Total unique donor-cell samples used across all fresh cells.
    pub fresh_donor_samples: usize,
    /// Active-fluid population mass removed by newly covered cells.
    pub removed_mass: f64,
    /// Active-fluid population mass inserted into fresh cells.
    pub fresh_mass: f64,
    /// `fresh_mass - removed_mass` for the complete transition.
    pub net_mass_change: f64,
    /// Raw population momentum removed by newly covered cells.
    pub removed_momentum: [f64; 2],
    /// Raw population momentum inserted into fresh cells.
    pub fresh_momentum: [f64; 2],
    /// `fresh_momentum - removed_momentum` for the complete transition.
    pub net_momentum_change: [f64; 2],
}

/// Low-Mach regularized velocity inlet at `x = 0` and isothermal
/// pressure/density outlet at `x = nx - 1`.
///
/// The inlet copies density and non-equilibrium stress from `x = 1` while
/// prescribing velocity. The outlet copies velocity and stress from
/// `x = nx - 2` while prescribing density. This boundary pair is intended for
/// x-directed crossflow fixtures with periodic y closure.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VelocityPressureX2 {
    inlet_velocity: [f64; 2],
    outlet_density: f64,
}

impl VelocityPressureX2 {
    /// Construct a checked low-Mach boundary pair.
    #[must_use]
    pub fn new(inlet_velocity: [f64; 2], outlet_density: f64) -> Self {
        assert!(
            inlet_velocity.iter().all(|component| component.is_finite()),
            "D2Q9 inlet velocity must be finite"
        );
        let speed_sq = inlet_velocity
            .iter()
            .map(|component| component * component)
            .sum::<f64>();
        assert!(
            speed_sq < MAX_REGULARIZED_BOUNDARY_SPEED_SQ,
            "D2Q9 inlet velocity exceeds the low-Mach boundary envelope"
        );
        assert!(
            outlet_density.is_finite() && outlet_density > 0.0,
            "D2Q9 outlet density must be positive and finite"
        );
        Self {
            inlet_velocity,
            outlet_density,
        }
    }

    /// Prescribed inlet velocity in lattice units.
    #[must_use]
    pub const fn inlet_velocity(self) -> [f64; 2] {
        self.inlet_velocity
    }

    /// Prescribed outlet density in lattice units.
    #[must_use]
    pub const fn outlet_density(self) -> f64 {
        self.outlet_density
    }
}

fn raw_population_moments(populations: &[f64; Q]) -> (f64, [f64; 2]) {
    let mut mass = 0.0;
    let mut momentum = [0.0; 2];
    for (q, &population) in populations.iter().enumerate() {
        mass += population;
        momentum[0] += f64::from(E[q].0) * population;
        momentum[1] += f64::from(E[q].1) * population;
    }
    (mass, momentum)
}

impl Grid {
    /// A grid of fluid at rest (unit density), uniform `tau`.
    #[must_use]
    pub fn uniform(nx: usize, ny: usize, tau: f64) -> Grid {
        assert!(nx > 0 && ny > 0, "grid dimensions must be positive");
        assert!(
            tau.is_finite() && tau > 0.5,
            "relaxation time tau must be finite and greater than 0.5"
        );
        let f0 = crate::equilibrium(1.0, 0.0, 0.0);
        Grid {
            nx,
            ny,
            flags: vec![Cell::Fluid; nx * ny],
            f: vec![f0; nx * ny],
            tau: vec![tau; nx * ny],
            g: [0.0, 0.0],
            fext: vec![[0.0; 2]; nx * ny],
            periodic_x: true,
            periodic_y: true,
            collision: CollisionModel2::default(),
        }
    }

    /// Row-major index.
    #[must_use]
    pub fn idx(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }

    /// Moments of cell `i` (Guo half-force correction).
    #[must_use]
    pub fn moments(&self, i: usize) -> Moments {
        let f = &self.f[i];
        let rho: f64 = f.iter().sum();
        assert!(
            rho.is_finite() && rho > 0.0,
            "moments require positive finite density"
        );
        let mut m = [0.0f64; 2];
        for (q, fi) in f.iter().enumerate() {
            m[0] += f64::from(E[q].0) * fi;
            m[1] += f64::from(E[q].1) * fi;
        }
        let fx = self.g[0].mul_add(rho, self.fext[i][0]);
        let fy = self.g[1].mul_add(rho, self.fext[i][1]);
        Moments {
            rho,
            u: [(m[0] + 0.5 * fx) / rho, (m[1] + 0.5 * fy) / rho],
        }
    }

    /// Total mass over non-gas cells.
    #[must_use]
    pub fn total_mass(&self) -> f64 {
        self.f
            .iter()
            .zip(&self.flags)
            .filter(|&(_, &fl)| fl != Cell::Gas && fl != Cell::Wall)
            .map(|(c, _)| c.iter().sum::<f64>())
            .sum()
    }

    /// Atomically replace the D2Q9 wall mask and initialize fresh fluid cells.
    ///
    /// This narrow moving-topology rung accepts only a fluid/wall domain.
    /// Fresh cells average unique surviving one-ring fluid populations,
    /// relaxation times, and external forces from the pre-transition state;
    /// freshly covered cells have their populations cleared. Every donor and
    /// receipt is validated before any grid field is published.
    pub fn transition_wall_topology(&mut self, next_walls: &[bool]) -> WallTopologyTransition2 {
        let cell_count = self.nx * self.ny;
        assert_eq!(
            next_walls.len(),
            cell_count,
            "next wall mask length must match the grid"
        );
        assert_eq!(
            self.flags.len(),
            cell_count,
            "cell flags must cover the grid"
        );
        assert_eq!(self.f.len(), cell_count, "populations must cover the grid");
        assert_eq!(
            self.tau.len(),
            cell_count,
            "relaxation times must cover the grid"
        );
        assert_eq!(
            self.fext.len(),
            cell_count,
            "external forces must cover the grid"
        );
        assert!(
            self.flags
                .iter()
                .all(|flag| matches!(*flag, Cell::Fluid | Cell::Wall)),
            "wall topology transition currently requires a fluid/wall-only domain"
        );
        assert!(
            next_walls.iter().any(|&is_wall| !is_wall),
            "wall topology transition must leave at least one fluid cell"
        );

        let mut next_flags = self.flags.clone();
        let mut next_populations = self.f.clone();
        let mut next_tau = self.tau.clone();
        let mut next_external_force = self.fext.clone();
        let mut receipt = WallTopologyTransition2::default();

        for y in 0..self.ny {
            for x in 0..self.nx {
                let index = self.idx(x, y);
                match (self.flags[index], next_walls[index]) {
                    (Cell::Fluid, true) => {
                        assert!(
                            self.f[index].into_iter().all(f64::is_finite),
                            "covered fluid populations must be finite"
                        );
                        let (mass, momentum) = raw_population_moments(&self.f[index]);
                        assert!(
                            mass.is_finite()
                                && mass > 0.0
                                && momentum.into_iter().all(f64::is_finite),
                            "covered fluid cell must have positive finite mass and finite momentum"
                        );
                        receipt.covered_fluid_cells += 1;
                        receipt.removed_mass += mass;
                        receipt.removed_momentum[0] += momentum[0];
                        receipt.removed_momentum[1] += momentum[1];
                        next_flags[index] = Cell::Wall;
                        next_populations[index] = [0.0; Q];
                    }
                    (Cell::Wall, false) => {
                        let mut donors = Vec::with_capacity(Q - 1);
                        for direction in 1..Q {
                            let Some(donor) = self.source(x, y, direction) else {
                                continue;
                            };
                            if self.flags[donor] == Cell::Fluid
                                && !next_walls[donor]
                                && !donors.contains(&donor)
                            {
                                donors.push(donor);
                            }
                        }
                        assert!(
                            !donors.is_empty(),
                            "fresh fluid cell requires a surviving one-ring fluid donor"
                        );

                        let mut populations = [0.0; Q];
                        let mut relaxation_time = 0.0;
                        let mut external_force = [0.0; 2];
                        for &donor in &donors {
                            assert!(
                                self.f[donor].into_iter().all(f64::is_finite),
                                "fresh-cell donor populations must be finite"
                            );
                            let (donor_mass, donor_momentum) =
                                raw_population_moments(&self.f[donor]);
                            assert!(
                                donor_mass.is_finite()
                                    && donor_mass > 0.0
                                    && donor_momentum.into_iter().all(f64::is_finite),
                                "fresh-cell donor must have positive finite mass and finite momentum"
                            );
                            assert!(
                                self.tau[donor].is_finite() && self.tau[donor] > 0.5,
                                "fresh-cell donor relaxation time must be finite and greater than 0.5"
                            );
                            assert!(
                                self.fext[donor].into_iter().all(f64::is_finite),
                                "fresh-cell donor external force must be finite"
                            );
                            for q in 0..Q {
                                populations[q] += self.f[donor][q];
                            }
                            relaxation_time += self.tau[donor];
                            external_force[0] += self.fext[donor][0];
                            external_force[1] += self.fext[donor][1];
                        }
                        let inverse_donor_count = 1.0 / donors.len() as f64;
                        for population in &mut populations {
                            *population *= inverse_donor_count;
                        }
                        relaxation_time *= inverse_donor_count;
                        external_force[0] *= inverse_donor_count;
                        external_force[1] *= inverse_donor_count;

                        let (mass, momentum) = raw_population_moments(&populations);
                        assert!(
                            populations.into_iter().all(f64::is_finite)
                                && mass.is_finite()
                                && mass > 0.0
                                && momentum.into_iter().all(f64::is_finite)
                                && relaxation_time.is_finite()
                                && relaxation_time > 0.5
                                && external_force.into_iter().all(f64::is_finite),
                            "fresh-cell averaged state must remain physically admissible"
                        );
                        receipt.fresh_fluid_cells += 1;
                        receipt.fresh_donor_samples += donors.len();
                        receipt.fresh_mass += mass;
                        receipt.fresh_momentum[0] += momentum[0];
                        receipt.fresh_momentum[1] += momentum[1];
                        next_flags[index] = Cell::Fluid;
                        next_populations[index] = populations;
                        next_tau[index] = relaxation_time;
                        next_external_force[index] = external_force;
                    }
                    (Cell::Fluid, false) | (Cell::Wall, true) => {}
                    (Cell::Interface | Cell::Gas, _) => unreachable!(
                        "fluid/wall-only transition was validated before proposal construction"
                    ),
                }
            }
        }

        receipt.net_mass_change = receipt.fresh_mass - receipt.removed_mass;
        receipt.net_momentum_change = [
            receipt.fresh_momentum[0] - receipt.removed_momentum[0],
            receipt.fresh_momentum[1] - receipt.removed_momentum[1],
        ];
        assert!(
            [
                receipt.removed_mass,
                receipt.fresh_mass,
                receipt.net_mass_change,
                receipt.removed_momentum[0],
                receipt.removed_momentum[1],
                receipt.fresh_momentum[0],
                receipt.fresh_momentum[1],
                receipt.net_momentum_change[0],
                receipt.net_momentum_change[1],
            ]
            .into_iter()
            .all(f64::is_finite),
            "wall topology transition receipt overflowed"
        );

        self.flags = next_flags;
        self.f = next_populations;
        self.tau = next_tau;
        self.fext = next_external_force;
        receipt
    }

    /// Collide (per-cell tau, vector Guo forcing) into `post`,
    /// honoring [`Grid::collision`]. The BGK branch is byte-for-byte
    /// the historical arithmetic.
    pub fn collide_into(&self, post: &mut Vec<[f64; Q]>) {
        post.clear();
        post.resize(self.nx * self.ny, [0.0; Q]);
        let model = self.collision;
        if let CollisionModel2::Trt { magic } = model {
            assert!(
                magic.is_finite() && magic > 0.0,
                "TRT magic parameter must be finite and positive"
            );
        }
        for (i, out) in post.iter_mut().enumerate().take(self.nx * self.ny) {
            if !matches!(self.flags[i], Cell::Fluid | Cell::Interface) {
                *out = self.f[i];
                continue;
            }
            let mm = self.moments(i);
            let (rho, ux, uy) = (mm.rho, mm.u[0], mm.u[1]);
            let feq = crate::equilibrium(rho, ux, uy);
            let tau = self.tau[i];
            assert!(
                tau.is_finite() && tau > 0.5,
                "cell relaxation time tau must be finite and greater than 0.5"
            );
            let coef = 1.0 - 0.5 / tau;
            let (gx, gy) = (
                self.g[0].mul_add(rho, self.fext[i][0]),
                self.g[1].mul_add(rho, self.fext[i][1]),
            );
            if matches!(model, CollisionModel2::CentralMoment) {
                // Central moments K[m][n] = Sum f (ex-ux)^m (ey-uy)^n.
                let mut k = [[0.0f64; 3]; 3];
                for q in 0..Q {
                    let cx = f64::from(E[q].0) - ux;
                    let cy = f64::from(E[q].1) - uy;
                    let (cx2, cy2) = (cx * cx, cy * cy);
                    let fi = self.f[i][q];
                    let xs = [1.0, cx, cx2];
                    let ys = [1.0, cy, cy2];
                    for (m, xv) in xs.iter().enumerate() {
                        for (n, yv) in ys.iter().enumerate() {
                            k[m][n] += fi * xv * yv;
                        }
                    }
                }
                // Relax: conserved k00/k10/k01 untouched; shear pair
                // at 1/tau toward 0 (k11) and 0 (k20 - k02); bulk,
                // third, fourth equilibrated to Maxwellian values.
                let omega = 1.0 / tau;
                let shear_diff = (1.0 - omega) * (k[2][0] - k[0][2]);
                let bulk = 2.0 * rho * CS2; // equilibrium k20 + k02
                k[1][1] *= 1.0 - omega;
                k[2][0] = 0.5 * (bulk + shear_diff);
                k[0][2] = 0.5 * (bulk - shear_diff);
                k[2][1] = 0.0;
                k[1][2] = 0.0;
                k[2][2] = rho * CS2 * CS2;
                // Central -> raw moments (binomial shift by +u).
                let mut m_raw = [[0.0f64; 3]; 3];
                let binom = [[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [1.0, 2.0, 1.0]];
                let uxp = [1.0, ux, ux * ux];
                let uyp = [1.0, uy, uy * uy];
                for m in 0..3 {
                    for n in 0..3 {
                        let mut acc = 0.0;
                        for p in 0..=m {
                            for r in 0..=n {
                                acc +=
                                    binom[m][p] * binom[n][r] * uxp[m - p] * uyp[n - r] * k[p][r];
                            }
                        }
                        m_raw[m][n] = acc;
                    }
                }
                // Exact tensor-product reconstruction on {-1,0,1}^2:
                // 1D: g(-1) = (mu2 - mu1)/2, g(0) = mu0 - mu2,
                // g(1) = (mu2 + mu1)/2, applied per axis.
                let g = |mu: [f64; 3], a: i32| -> f64 {
                    match a {
                        -1 => 0.5 * (mu[2] - mu[1]),
                        0 => mu[0] - mu[2],
                        _ => 0.5 * (mu[2] + mu[1]),
                    }
                };
                for q in 0..Q {
                    let (ax, ay) = (E[q].0, E[q].1);
                    let partial = [
                        g([m_raw[0][0], m_raw[1][0], m_raw[2][0]], ax),
                        g([m_raw[0][1], m_raw[1][1], m_raw[2][1]], ax),
                        g([m_raw[0][2], m_raw[1][2], m_raw[2][2]], ax),
                    ];
                    let (ex, ey) = (f64::from(ax), f64::from(ay));
                    let eu = ex * ux + ey * uy;
                    let fx = (ex - ux) / CS2 + eu * ex / (CS2 * CS2);
                    let fy = (ey - uy) / CS2 + eu * ey / (CS2 * CS2);
                    let force = coef * W[q] * (fx * gx + fy * gy);
                    out[q] = g(partial, ay) + force;
                }
            } else if let CollisionModel2::Trt { magic } = model {
                // Even/odd split of the non-equilibrium; the odd
                // relaxation rate comes from the magic combination.
                let omega_even = 1.0 / tau;
                let tau_odd = 0.5 + magic / (tau - 0.5);
                let omega_odd = 1.0 / tau_odd;
                let mut fneq = [0.0f64; Q];
                for q in 0..Q {
                    fneq[q] = self.f[i][q] - feq[q];
                }
                for q in 0..Q {
                    let (ex, ey) = (f64::from(E[q].0), f64::from(E[q].1));
                    let eu = ex * ux + ey * uy;
                    let fx = (ex - ux) / CS2 + eu * ex / (CS2 * CS2);
                    let fy = (ey - uy) / CS2 + eu * ey / (CS2 * CS2);
                    let force = coef * W[q] * (fx * gx + fy * gy);
                    let even = 0.5 * (fneq[q] + fneq[OPP[q]]);
                    let odd = 0.5 * (fneq[q] - fneq[OPP[q]]);
                    out[q] = self.f[i][q] - omega_even * even - omega_odd * odd + force;
                }
            } else if matches!(
                model,
                CollisionModel2::Regularized | CollisionModel2::RecursiveRegularized
            ) {
                // Second moment of the non-equilibrium.
                let (mut pxx, mut pyy, mut pxy) = (0.0f64, 0.0f64, 0.0f64);
                for q in 0..Q {
                    let (ex, ey) = (f64::from(E[q].0), f64::from(E[q].1));
                    let fneq = self.f[i][q] - feq[q];
                    pxx += ex * ex * fneq;
                    pyy += ey * ey * fneq;
                    pxy += ex * ey * fneq;
                }
                let relax = 1.0 - 1.0 / tau;
                // Recursive third-order Hermite coefficients (RR3);
                // zero for the plain projected operator.
                let (a3_xxy, a3_xyy) = if matches!(model, CollisionModel2::RecursiveRegularized) {
                    (2.0 * ux * pxy + uy * pxx, 2.0 * uy * pxy + ux * pyy)
                } else {
                    (0.0, 0.0)
                };
                for q in 0..Q {
                    let (ex, ey) = (f64::from(E[q].0), f64::from(E[q].1));
                    let eu = ex * ux + ey * uy;
                    let fx = (ex - ux) / CS2 + eu * ex / (CS2 * CS2);
                    let fy = (ey - uy) / CS2 + eu * ey / (CS2 * CS2);
                    let force = coef * W[q] * (fx * gx + fy * gy);
                    // w_i / (2 cs^4) = 4.5 w_i; Q_i : Pi_neq with the
                    // symmetric off-diagonal counted twice.
                    let qxx = ex * ex - CS2;
                    let qyy = ey * ey - CS2;
                    // Third-order Hermite polynomials H_xxy, H_xyy;
                    // w_i / (2 cs^6) = 13.5 w_i.
                    let h_xxy = qxx * ey;
                    let h_xyy = ex * qyy;
                    let fneq_reg = 4.5 * W[q] * (qxx * pxx + qyy * pyy + 2.0 * ex * ey * pxy)
                        + 13.5 * W[q] * (h_xxy * a3_xxy + h_xyy * a3_xyy);
                    out[q] = feq[q] + relax * fneq_reg + force;
                }
            } else {
                for q in 0..Q {
                    let (ex, ey) = (f64::from(E[q].0), f64::from(E[q].1));
                    let eu = ex * ux + ey * uy;
                    // Guo forcing, vector form.
                    let fx = (ex - ux) / CS2 + eu * ex / (CS2 * CS2);
                    let fy = (ey - uy) / CS2 + eu * ey / (CS2 * CS2);
                    let force = coef * W[q] * (fx * gx + fy * gy);
                    out[q] = self.f[i][q] + (feq[q] - self.f[i][q]) / tau + force;
                }
            }
        }
    }

    /// Collide with the Noble-Torczynski partially saturated cell operator.
    ///
    /// For solid fraction `epsilon_s`, lattice time step one, and local
    /// relaxation time `tau`, the coupling weight is
    ///
    /// `B = epsilon_s (tau - 1/2) / ((1 - epsilon_s) + (tau - 1/2))`.
    ///
    /// The ordinary BGK collision `Omega_f` is blended with the published
    /// non-equilibrium-bounce-back solid collision
    ///
    /// `Omega_s[q] = f[opp(q)] - f[q] + feq[q](rho, u_s) - feq[opp(q)](rho, u)`.
    ///
    /// The output and aggregate exchange receipt are published only after the
    /// full force-free request and resulting proposal are finite. See Noble
    /// and Torczynski (1998), <https://doi.org/10.1142/S0129183198001084>.
    #[allow(clippy::too_many_lines)] // Keep the auditable collision/receipt order contiguous.
    pub fn collide_into_with_partial_saturation(
        &self,
        post: &mut Vec<[f64; Q]>,
        cells: &[PartialSaturationCell2],
        moment_origin: [f64; 2],
    ) -> PartialSaturationExchange2 {
        self.validate_partial_saturation_fields(cells, moment_origin);

        // Build in a private proposal so a later arithmetic refusal cannot
        // leave the caller's reusable scratch buffer partially overwritten.
        let mut proposal = Vec::new();
        self.collide_into(&mut proposal);
        let mut receipt = PartialSaturationExchange2::default();

        for y in 0..self.ny {
            for x in 0..self.nx {
                let index = self.idx(x, y);
                let cell = cells[index];
                if cell.solid_fraction <= 0.0 {
                    continue;
                }

                let moments = self.moments(index);
                let fluid_equilibrium = crate::equilibrium(moments.rho, moments.u[0], moments.u[1]);
                let solid_equilibrium =
                    crate::equilibrium(moments.rho, cell.solid_velocity[0], cell.solid_velocity[1]);
                let tau = self.tau[index];
                let coupling_weight = cell.coupling_weight(tau);
                assert!(
                    coupling_weight.is_finite() && (0.0..=1.0).contains(&coupling_weight),
                    "partial-saturation coupling weight at cell {index} left [0, 1]"
                );

                let mut fluid_mass_change = 0.0;
                let mut fluid_impulse = [0.0; 2];
                for q in 0..Q {
                    let fluid_collision = (fluid_equilibrium[q] - self.f[index][q]) / tau;
                    let solid_collision = self.f[index][OPP[q]] - self.f[index][q]
                        + solid_equilibrium[q]
                        - fluid_equilibrium[OPP[q]];
                    let ordinary_population = proposal[index][q];
                    let blended_population = self.f[index][q]
                        + (1.0 - coupling_weight) * fluid_collision
                        + coupling_weight * solid_collision;
                    let coupling_correction = blended_population - ordinary_population;
                    proposal[index][q] = blended_population;
                    fluid_mass_change += coupling_correction;
                    fluid_impulse[0] += f64::from(E[q].0) * coupling_correction;
                    fluid_impulse[1] += f64::from(E[q].1) * coupling_correction;
                }
                assert!(
                    proposal[index].into_iter().all(f64::is_finite)
                        && fluid_mass_change.is_finite()
                        && fluid_impulse.into_iter().all(f64::is_finite),
                    "partial-saturation collision at cell {index} produced a non-finite proposal"
                );

                let solid_impulse = [-fluid_impulse[0], -fluid_impulse[1]];
                let offset = [x as f64 - moment_origin[0], y as f64 - moment_origin[1]];
                let fluid_angular_impulse =
                    offset[0] * fluid_impulse[1] - offset[1] * fluid_impulse[0];
                let solid_angular_impulse = -fluid_angular_impulse;

                receipt.fluid_population_mass_change += fluid_mass_change;
                receipt.fluid_population_impulse[0] += fluid_impulse[0];
                receipt.fluid_population_impulse[1] += fluid_impulse[1];
                receipt.solid_impulse[0] += solid_impulse[0];
                receipt.solid_impulse[1] += solid_impulse[1];
                receipt.fluid_population_angular_impulse += fluid_angular_impulse;
                receipt.solid_angular_impulse += solid_angular_impulse;
                receipt.solid_work += solid_impulse[0].mul_add(
                    cell.solid_velocity[0],
                    solid_impulse[1] * cell.solid_velocity[1],
                );
                receipt.coupled_cells += 1;
                receipt.coupling_weight_sum += coupling_weight;
            }
        }

        assert!(
            proposal.iter().flatten().copied().all(f64::is_finite),
            "partial-saturation collision produced a non-finite full-grid proposal"
        );
        assert!(
            [
                receipt.fluid_population_mass_change,
                receipt.fluid_population_impulse[0],
                receipt.fluid_population_impulse[1],
                receipt.solid_impulse[0],
                receipt.solid_impulse[1],
                receipt.fluid_population_angular_impulse,
                receipt.solid_angular_impulse,
                receipt.solid_work,
                receipt.coupling_weight_sum,
            ]
            .into_iter()
            .all(f64::is_finite),
            "partial-saturation aggregate receipt overflowed"
        );

        *post = proposal;
        receipt
    }

    /// Source cell for pull-streaming direction `q` into (x, y);
    /// `None` when the pull crosses a non-periodic boundary (treated
    /// as wall bounce-back).
    #[must_use]
    pub fn source(&self, x: usize, y: usize, q: usize) -> Option<usize> {
        let (ex, ey) = E[q];
        let sx = match ex {
            1 => {
                if x == 0 {
                    if self.periodic_x {
                        self.nx - 1
                    } else {
                        return None;
                    }
                } else {
                    x - 1
                }
            }
            -1 => {
                if x + 1 == self.nx {
                    if self.periodic_x {
                        0
                    } else {
                        return None;
                    }
                } else {
                    x + 1
                }
            }
            _ => x,
        };
        let sy = match ey {
            1 => {
                if y == 0 {
                    if self.periodic_y {
                        self.ny - 1
                    } else {
                        return None;
                    }
                } else {
                    y - 1
                }
            }
            -1 => {
                if y + 1 == self.ny {
                    if self.periodic_y {
                        0
                    } else {
                        return None;
                    }
                } else {
                    y + 1
                }
            }
            _ => y,
        };
        Some(self.idx(sx, sy))
    }

    fn validate_stream_input(&self, post: &[[f64; Q]]) {
        assert!(
            post.len() >= self.nx * self.ny,
            "post-collision populations must cover every grid cell"
        );
    }

    fn validate_partial_saturation_fields(
        &self,
        cells: &[PartialSaturationCell2],
        moment_origin: [f64; 2],
    ) {
        assert!(
            self.nx > 0 && self.ny > 0,
            "partial-saturation grid dimensions must be positive"
        );
        let cell_count = self
            .nx
            .checked_mul(self.ny)
            .expect("D2Q9 partial-saturation grid size overflowed");
        assert_eq!(
            cells.len(),
            cell_count,
            "partial-saturation field length must match the grid"
        );
        assert_eq!(
            self.flags.len(),
            cell_count,
            "partial-saturation grid flags must cover the grid"
        );
        assert_eq!(
            self.f.len(),
            cell_count,
            "partial-saturation populations must cover the grid"
        );
        assert_eq!(
            self.tau.len(),
            cell_count,
            "partial-saturation relaxation times must cover the grid"
        );
        assert_eq!(
            self.fext.len(),
            cell_count,
            "partial-saturation external forces must cover the grid"
        );
        assert!(
            moment_origin.into_iter().all(f64::is_finite),
            "partial-saturation moment origin must be finite"
        );
        assert!(
            self.g == [0.0; 2] && self.fext.iter().all(|&force| force == [0.0; 2]),
            "partial-saturation collision currently requires a force-free grid"
        );

        for (index, cell) in cells.iter().copied().enumerate() {
            assert!(
                self.f[index].into_iter().all(f64::is_finite),
                "partial-saturation populations at cell {index} must be finite"
            );
            if matches!(self.flags[index], Cell::Fluid | Cell::Interface) {
                let (density, momentum) = raw_population_moments(&self.f[index]);
                assert!(
                    density.is_finite()
                        && density > 0.0
                        && momentum.into_iter().all(f64::is_finite),
                    "partial-saturation active cell {index} requires positive finite density and finite momentum"
                );
                assert!(
                    self.tau[index].is_finite() && self.tau[index] > 0.5,
                    "partial-saturation active cell {index} requires finite tau greater than 0.5"
                );
            }
            assert!(
                cell.solid_fraction <= 0.0 || self.flags[index] == Cell::Fluid,
                "partial saturation may be positive only on fluid cell {index}"
            );
        }
    }

    fn validate_measured_walls(&self, measured_walls: &[bool]) {
        assert_eq!(
            measured_walls.len(),
            self.nx * self.ny,
            "measured-wall mask length must match the grid"
        );
        for (index, (&measured, &flag)) in measured_walls.iter().zip(&self.flags).enumerate() {
            assert!(
                !measured || flag == Cell::Wall,
                "measured-wall mask selects non-wall cell {index}"
            );
        }
    }

    fn validate_moving_wall_fields(&self, wall_velocities: &[[f64; 2]], moment_origin: [f64; 2]) {
        assert_eq!(
            wall_velocities.len(),
            self.nx * self.ny,
            "moving-wall velocity field length must match the grid"
        );
        assert!(
            moment_origin.into_iter().all(f64::is_finite),
            "moving-wall moment origin must be finite"
        );
        for (index, (&velocity, &flag)) in wall_velocities.iter().zip(&self.flags).enumerate() {
            assert!(
                velocity.into_iter().all(f64::is_finite),
                "moving-wall velocity at cell {index} must be finite"
            );
            let speed_sq = velocity[0].mul_add(velocity[0], velocity[1] * velocity[1]);
            assert!(
                speed_sq < MAX_REGULARIZED_BOUNDARY_SPEED_SQ,
                "moving-wall velocity at cell {index} exceeds the low-Mach admission envelope"
            );
            assert!(
                flag == Cell::Wall || velocity == [0.0; 2],
                "moving-wall velocity field assigns motion to non-wall cell {index}"
            );
        }
    }

    fn validate_moving_wall_post(&self, post: &[[f64; Q]], wall_velocities: &[[f64; 2]]) {
        self.validate_stream_input(post);
        for y in 0..self.ny {
            for x in 0..self.nx {
                let i = self.idx(x, y);
                if !matches!(self.flags[i], Cell::Fluid | Cell::Interface) {
                    continue;
                }
                let mut needs_density = false;
                for q in 0..Q {
                    let Some(source) = self.source(x, y, q) else {
                        continue;
                    };
                    if self.flags[source] != Cell::Wall {
                        continue;
                    }
                    assert!(
                        post[i][OPP[q]].is_finite(),
                        "moving-wall outgoing population must be finite"
                    );
                    needs_density |= wall_velocities[source] != [0.0; 2];
                }
                if needs_density {
                    let rho_post = post[i].iter().sum::<f64>();
                    assert!(
                        rho_post.is_finite() && rho_post > 0.0,
                        "moving-wall bounce-back requires positive finite post-collision density"
                    );
                }
            }
        }
    }

    fn validate_curved_moving_wall_fields(
        &self,
        wall_links: &[[Option<CurvedMovingWallLink2>; Q]],
        moment_origin: [f64; 2],
    ) {
        assert_eq!(
            wall_links.len(),
            self.nx * self.ny,
            "curved-wall link field length must match the grid"
        );
        assert!(
            moment_origin.into_iter().all(f64::is_finite),
            "curved-wall moment origin must be finite"
        );
        for y in 0..self.ny {
            for x in 0..self.nx {
                let index = self.idx(x, y);
                let active = matches!(self.flags[index], Cell::Fluid | Cell::Interface);
                for q in 0..Q {
                    let wall_source = active
                        .then(|| self.source(x, y, q))
                        .flatten()
                        .filter(|&source| self.flags[source] == Cell::Wall);
                    match (wall_source, wall_links[index][q]) {
                        (Some(_), Some(link)) => {
                            assert!(
                                link.fraction.is_finite()
                                    && link.fraction > 0.0
                                    && link.fraction <= 1.0,
                                "curved-wall link fraction at cell {index}, direction {q} must be finite and in (0, 1]"
                            );
                            assert!(
                                link.wall_velocity.into_iter().all(f64::is_finite),
                                "curved-wall velocity at cell {index}, direction {q} must be finite"
                            );
                            let speed_sq = link.wall_velocity[0].mul_add(
                                link.wall_velocity[0],
                                link.wall_velocity[1] * link.wall_velocity[1],
                            );
                            assert!(
                                speed_sq < MAX_REGULARIZED_BOUNDARY_SPEED_SQ,
                                "curved-wall velocity at cell {index}, direction {q} exceeds the low-Mach admission envelope"
                            );
                        }
                        (Some(_), None) => panic!(
                            "curved-wall geometry omits wall link at cell {index}, direction {q}"
                        ),
                        (None, Some(_)) => panic!(
                            "curved-wall geometry assigns a link without a wall at cell {index}, direction {q}"
                        ),
                        (None, None) => {}
                    }
                }
            }
        }
    }

    fn validate_curved_moving_wall_post(
        &self,
        post: &[[f64; Q]],
        wall_links: &[[Option<CurvedMovingWallLink2>; Q]],
    ) {
        self.validate_stream_input(post);
        for y in 0..self.ny {
            for x in 0..self.nx {
                let index = self.idx(x, y);
                if !matches!(self.flags[index], Cell::Fluid | Cell::Interface) {
                    continue;
                }
                let mut needs_density = false;
                for q in 0..Q {
                    let Some(link) = wall_links[index][q] else {
                        continue;
                    };
                    assert!(
                        post[index][OPP[q]].is_finite(),
                        "curved-wall outgoing population at cell {index}, direction {q} must be finite"
                    );
                    if link.fraction < 0.5 {
                        let far = self.source(x, y, OPP[q]).unwrap_or_else(|| {
                            panic!(
                                "curved-wall short-link branch at cell {index}, direction {q} requires an in-domain far-fluid donor"
                            )
                        });
                        assert!(
                            matches!(self.flags[far], Cell::Fluid | Cell::Interface),
                            "curved-wall short-link branch at cell {index}, direction {q} requires a far-fluid donor"
                        );
                        assert!(
                            post[far][OPP[q]].is_finite(),
                            "curved-wall far-fluid population at cell {index}, direction {q} must be finite"
                        );
                    } else if link.fraction > 0.5 {
                        assert!(
                            post[index][q].is_finite(),
                            "curved-wall long-link population at cell {index}, direction {q} must be finite"
                        );
                    }
                    needs_density |= link.wall_velocity != [0.0; 2];
                }
                if needs_density {
                    let rho_post = post[index].iter().sum::<f64>();
                    assert!(
                        rho_post.is_finite() && rho_post > 0.0,
                        "curved moving-wall bounce-back requires positive finite post-collision density"
                    );
                }
            }
        }
    }

    fn validate_velocity_pressure_x(&self) {
        assert!(
            self.nx >= 3,
            "D2Q9 x-open flow requires at least three columns"
        );
        assert!(
            !self.periodic_x && self.periodic_y,
            "D2Q9 x-open flow requires non-periodic x and periodic y"
        );
        assert!(
            self.g.into_iter().all(|component| component == 0.0)
                && self
                    .fext
                    .iter()
                    .flatten()
                    .all(|component| *component == 0.0),
            "D2Q9 regularized open boundaries currently require zero body force"
        );
        for y in 0..self.ny {
            for x in [0, 1, self.nx - 2, self.nx - 1] {
                assert_eq!(
                    self.flags[self.idx(x, y)],
                    Cell::Fluid,
                    "D2Q9 open faces and first-interior columns must be fluid"
                );
            }
        }
    }

    fn regularized_source(&self, index: usize) -> (f64, [f64; 2], [[f64; 2]; 2]) {
        let moments = self.moments(index);
        let equilibrium = crate::equilibrium(moments.rho, moments.u[0], moments.u[1]);
        let mut stress = [[0.0; 2]; 2];
        for q in 0..Q {
            let nonequilibrium = self.f[index][q] - equilibrium[q];
            let e = [f64::from(E[q].0), f64::from(E[q].1)];
            for row in 0..2 {
                for column in 0..2 {
                    stress[row][column] += e[row] * e[column] * nonequilibrium;
                }
            }
        }
        (moments.rho, moments.u, stress)
    }

    fn apply_velocity_pressure_x(&mut self, boundary: VelocityPressureX2) {
        for y in 0..self.ny {
            let inlet = self.idx(0, y);
            let inlet_source = self.idx(1, y);
            let outlet_source = self.idx(self.nx - 2, y);
            let outlet = self.idx(self.nx - 1, y);
            let (inlet_density, _, inlet_stress) = self.regularized_source(inlet_source);
            let (_, outlet_velocity, outlet_stress) = self.regularized_source(outlet_source);
            self.f[inlet] =
                regularized_populations(inlet_density, boundary.inlet_velocity, inlet_stress);
            self.f[outlet] =
                regularized_populations(boundary.outlet_density, outlet_velocity, outlet_stress);
        }
    }

    fn stream_from_inner(
        &mut self,
        post: &[[f64; Q]],
        measured_walls: Option<&[bool]>,
    ) -> MomentumExchange2 {
        let mut receipt = MomentumExchange2::default();
        for y in 0..self.ny {
            for x in 0..self.nx {
                let i = self.idx(x, y);
                if !matches!(self.flags[i], Cell::Fluid | Cell::Interface) {
                    continue;
                }
                for q in 0..Q {
                    let pulled = match self.source(x, y, q) {
                        Some(s) if matches!(self.flags[s], Cell::Wall | Cell::Gas) => {
                            let reflected = post[i][OPP[q]];
                            if self.flags[s] == Cell::Wall
                                && measured_walls.is_some_and(|mask| mask[s])
                            {
                                // Pull direction q points from the wall back into
                                // the fluid. The fluid momentum change is
                                // +2 f_post c_q, so the opposite impulse delivered
                                // to the stationary wall is -2 f_post c_q.
                                receipt.wall_impulse[0] -= 2.0 * reflected * f64::from(E[q].0);
                                receipt.wall_impulse[1] -= 2.0 * reflected * f64::from(E[q].1);
                                receipt.measured_links += 1;
                            }
                            reflected
                        }
                        Some(s) => post[s][q],
                        None => post[i][OPP[q]],
                    };
                    self.f[i][q] = pulled;
                }
            }
        }
        receipt
    }

    fn stream_from_moving_wall_inner(
        &mut self,
        post: &[[f64; Q]],
        measured_walls: &[bool],
        wall_velocities: &[[f64; 2]],
        moment_origin: [f64; 2],
    ) -> MovingWallMomentumExchange2 {
        let mut receipt = MovingWallMomentumExchange2::default();
        for y in 0..self.ny {
            for x in 0..self.nx {
                let i = self.idx(x, y);
                if !matches!(self.flags[i], Cell::Fluid | Cell::Interface) {
                    continue;
                }
                let mut rho_post = None;
                for q in 0..Q {
                    let pulled = match self.source(x, y, q) {
                        Some(source) if self.flags[source] == Cell::Wall => {
                            let outgoing = post[i][OPP[q]];
                            let wall_velocity = wall_velocities[source];
                            let incoming = if wall_velocity == [0.0; 2] {
                                outgoing
                            } else {
                                let rho =
                                    *rho_post.get_or_insert_with(|| post[i].iter().sum::<f64>());
                                let e = [f64::from(E[q].0), f64::from(E[q].1)];
                                let e_dot_wall =
                                    e[0].mul_add(wall_velocity[0], e[1] * wall_velocity[1]);
                                outgoing + 2.0 * W[q] * rho * e_dot_wall / CS2
                            };

                            if measured_walls[source] {
                                let e = [f64::from(E[q].0), f64::from(E[q].1)];
                                let wall_impulse = if wall_velocity == [0.0; 2] {
                                    [-2.0 * outgoing * e[0], -2.0 * outgoing * e[1]]
                                } else {
                                    [
                                        (-e[0] - wall_velocity[0]) * outgoing
                                            - (e[0] - wall_velocity[0]) * incoming,
                                        (-e[1] - wall_velocity[1]) * outgoing
                                            - (e[1] - wall_velocity[1]) * incoming,
                                    ]
                                };
                                let fluid_population_impulse =
                                    [e[0] * (incoming + outgoing), e[1] * (incoming + outgoing)];
                                let fluid_mass_change = incoming - outgoing;
                                let wall_velocity_mass_impulse = [
                                    wall_velocity[0] * fluid_mass_change,
                                    wall_velocity[1] * fluid_mass_change,
                                ];
                                let link_offset = [
                                    (x as f64 - 0.5 * e[0]) - moment_origin[0],
                                    (y as f64 - 0.5 * e[1]) - moment_origin[1],
                                ];

                                receipt.wall_impulse[0] += wall_impulse[0];
                                receipt.wall_impulse[1] += wall_impulse[1];
                                receipt.fluid_population_impulse[0] += fluid_population_impulse[0];
                                receipt.fluid_population_impulse[1] += fluid_population_impulse[1];
                                receipt.fluid_mass_change += fluid_mass_change;
                                receipt.wall_velocity_mass_impulse[0] +=
                                    wall_velocity_mass_impulse[0];
                                receipt.wall_velocity_mass_impulse[1] +=
                                    wall_velocity_mass_impulse[1];
                                receipt.wall_angular_impulse += link_offset[0] * wall_impulse[1]
                                    - link_offset[1] * wall_impulse[0];
                                receipt.fluid_population_angular_impulse += link_offset[0]
                                    * fluid_population_impulse[1]
                                    - link_offset[1] * fluid_population_impulse[0];
                                receipt.wall_velocity_mass_angular_impulse += link_offset[0]
                                    * wall_velocity_mass_impulse[1]
                                    - link_offset[1] * wall_velocity_mass_impulse[0];
                                receipt.wall_work += wall_impulse[0]
                                    .mul_add(wall_velocity[0], wall_impulse[1] * wall_velocity[1]);
                                receipt.measured_links += 1;
                            }
                            incoming
                        }
                        Some(source) if self.flags[source] == Cell::Gas => post[i][OPP[q]],
                        Some(source) => post[source][q],
                        None => post[i][OPP[q]],
                    };
                    self.f[i][q] = pulled;
                }
            }
        }
        receipt
    }

    fn stream_from_curved_moving_wall_inner(
        &mut self,
        post: &[[f64; Q]],
        measured_walls: &[bool],
        wall_links: &[[Option<CurvedMovingWallLink2>; Q]],
        moment_origin: [f64; 2],
    ) -> MovingWallMomentumExchange2 {
        let mut receipt = MovingWallMomentumExchange2::default();
        for y in 0..self.ny {
            for x in 0..self.nx {
                let index = self.idx(x, y);
                if !matches!(self.flags[index], Cell::Fluid | Cell::Interface) {
                    continue;
                }
                let mut rho_post = None;
                for q in 0..Q {
                    let pulled = match self.source(x, y, q) {
                        Some(source) if self.flags[source] == Cell::Wall => {
                            let link = wall_links[index][q].expect(
                                "validated curved-wall geometry must cover every wall link",
                            );
                            let fraction = link.fraction;
                            let wall_velocity = link.wall_velocity;
                            let outgoing = post[index][OPP[q]];
                            let e = [f64::from(E[q].0), f64::from(E[q].1)];
                            let moving_correction = if wall_velocity == [0.0; 2] {
                                0.0
                            } else {
                                let rho = *rho_post
                                    .get_or_insert_with(|| post[index].iter().sum::<f64>());
                                let e_dot_wall =
                                    e[0].mul_add(wall_velocity[0], e[1] * wall_velocity[1]);
                                2.0 * W[q] * rho * e_dot_wall / CS2
                            };
                            let incoming = if fraction == 0.5 {
                                // Preserve the established halfway operation tree exactly.
                                outgoing + moving_correction
                            } else if fraction < 0.5 {
                                let far = self.source(x, y, OPP[q]).expect(
                                    "validated curved-wall short link must have a far-fluid donor",
                                );
                                2.0 * fraction * outgoing
                                    + (1.0 - 2.0 * fraction) * post[far][OPP[q]]
                                    + moving_correction
                            } else {
                                let denominator = 2.0 * fraction;
                                (outgoing + moving_correction) / denominator
                                    + ((denominator - 1.0) / denominator) * post[index][q]
                            };

                            if measured_walls[source] {
                                let wall_impulse = if fraction == 0.5 && wall_velocity == [0.0; 2] {
                                    // Preserve the established stationary-halfway bits.
                                    [-2.0 * outgoing * e[0], -2.0 * outgoing * e[1]]
                                } else {
                                    [
                                        (-e[0] - wall_velocity[0]) * outgoing
                                            - (e[0] - wall_velocity[0]) * incoming,
                                        (-e[1] - wall_velocity[1]) * outgoing
                                            - (e[1] - wall_velocity[1]) * incoming,
                                    ]
                                };
                                let fluid_population_impulse =
                                    [e[0] * (incoming + outgoing), e[1] * (incoming + outgoing)];
                                let fluid_mass_change = incoming - outgoing;
                                let wall_velocity_mass_impulse = [
                                    wall_velocity[0] * fluid_mass_change,
                                    wall_velocity[1] * fluid_mass_change,
                                ];
                                let link_offset = [
                                    (x as f64 - fraction * e[0]) - moment_origin[0],
                                    (y as f64 - fraction * e[1]) - moment_origin[1],
                                ];

                                receipt.wall_impulse[0] += wall_impulse[0];
                                receipt.wall_impulse[1] += wall_impulse[1];
                                receipt.fluid_population_impulse[0] += fluid_population_impulse[0];
                                receipt.fluid_population_impulse[1] += fluid_population_impulse[1];
                                receipt.fluid_mass_change += fluid_mass_change;
                                receipt.wall_velocity_mass_impulse[0] +=
                                    wall_velocity_mass_impulse[0];
                                receipt.wall_velocity_mass_impulse[1] +=
                                    wall_velocity_mass_impulse[1];
                                receipt.wall_angular_impulse += link_offset[0] * wall_impulse[1]
                                    - link_offset[1] * wall_impulse[0];
                                receipt.fluid_population_angular_impulse += link_offset[0]
                                    * fluid_population_impulse[1]
                                    - link_offset[1] * fluid_population_impulse[0];
                                receipt.wall_velocity_mass_angular_impulse += link_offset[0]
                                    * wall_velocity_mass_impulse[1]
                                    - link_offset[1] * wall_velocity_mass_impulse[0];
                                receipt.wall_work += wall_impulse[0]
                                    .mul_add(wall_velocity[0], wall_impulse[1] * wall_velocity[1]);
                                receipt.measured_links += 1;
                            }
                            incoming
                        }
                        Some(source) if self.flags[source] == Cell::Gas => post[index][OPP[q]],
                        Some(source) => post[source][q],
                        None => post[index][OPP[q]],
                    };
                    self.f[index][q] = pulled;
                }
            }
        }
        receipt
    }

    /// Stream `post` into `self.f` (fluid pull; wall and out-of-domain
    /// pulls bounce back).
    pub fn stream_from(&mut self, post: &[[f64; Q]]) {
        self.validate_stream_input(post);
        let _ = self.stream_from_inner(post, None);
    }

    /// Stream while measuring momentum delivered to selected wall cells.
    ///
    /// Only bounce-back links whose source cell is both [`Cell::Wall`] and
    /// `true` in `measured_walls` contribute. Gas-boundary and non-periodic
    /// exterior bounces are deliberately excluded.
    pub fn stream_from_with_wall_momentum(
        &mut self,
        post: &[[f64; Q]],
        measured_walls: &[bool],
    ) -> MomentumExchange2 {
        self.validate_stream_input(post);
        self.validate_measured_walls(measured_walls);
        self.stream_from_inner(post, Some(measured_walls))
    }

    /// Stream with per-wall-cell velocities and a moving-wall exchange receipt.
    ///
    /// Moving halfway bounce-back is applied to every wall cell according to
    /// `wall_velocities`; `measured_walls` only selects which links contribute
    /// to the returned receipt. Torque uses each destination-local halfway-link
    /// midpoint about `moment_origin`. All request and post-collision inputs are
    /// validated before populations are mutated.
    pub fn stream_from_with_moving_wall_momentum(
        &mut self,
        post: &[[f64; Q]],
        measured_walls: &[bool],
        wall_velocities: &[[f64; 2]],
        moment_origin: [f64; 2],
    ) -> MovingWallMomentumExchange2 {
        self.validate_measured_walls(measured_walls);
        self.validate_moving_wall_fields(wall_velocities, moment_origin);
        self.validate_moving_wall_post(post, wall_velocities);
        self.stream_from_moving_wall_inner(post, measured_walls, wall_velocities, moment_origin)
    }

    /// Stream with linear Bouzidi-Firdaouss-Lallemand off-lattice wall links.
    ///
    /// `wall_links[destination][q]` must contain exactly one description for
    /// every pull link whose source is a wall, and `None` everywhere else.
    /// Link fractions below one half interpolate the outgoing population from
    /// the destination and next fluid node; fractions above one half
    /// interpolate the two local post-collision directions. At exactly one
    /// half this is bit-compatible with
    /// [`Grid::stream_from_with_moving_wall_momentum`]. Torque is evaluated at
    /// the declared off-lattice wall intersection. The linear rule follows
    /// Bouzidi, Firdaouss, and Lallemand (2001),
    /// <https://doi.org/10.1063/1.1399290>.
    pub fn stream_from_with_curved_moving_wall_momentum(
        &mut self,
        post: &[[f64; Q]],
        measured_walls: &[bool],
        wall_links: &[[Option<CurvedMovingWallLink2>; Q]],
        moment_origin: [f64; 2],
    ) -> MovingWallMomentumExchange2 {
        self.validate_measured_walls(measured_walls);
        self.validate_curved_moving_wall_fields(wall_links, moment_origin);
        self.validate_curved_moving_wall_post(post, wall_links);
        self.stream_from_curved_moving_wall_inner(post, measured_walls, wall_links, moment_origin)
    }

    /// One plain step (no free-surface bookkeeping).
    pub fn step(&mut self, scratch: &mut Vec<[f64; Q]>) {
        self.collide_into(scratch);
        let post = std::mem::take(scratch);
        self.stream_from(&post);
        *scratch = post;
    }

    /// One force-free partially saturated collision followed by plain stream.
    ///
    /// The complete cell field and collision proposal are admitted before the
    /// grid populations advance, so a malformed request leaves both `self`
    /// and the caller's scratch buffer unchanged.
    pub fn step_with_partial_saturation(
        &mut self,
        scratch: &mut Vec<[f64; Q]>,
        cells: &[PartialSaturationCell2],
        moment_origin: [f64; 2],
    ) -> PartialSaturationExchange2 {
        let receipt = self.collide_into_with_partial_saturation(scratch, cells, moment_origin);
        let post = std::mem::take(scratch);
        self.stream_from(&post);
        *scratch = post;
        receipt
    }

    /// One plain step plus a raw stationary-wall momentum-exchange receipt.
    ///
    /// The mask is validated before collision, so a malformed measurement
    /// request cannot partially advance the grid.
    pub fn step_with_wall_momentum(
        &mut self,
        scratch: &mut Vec<[f64; Q]>,
        measured_walls: &[bool],
    ) -> MomentumExchange2 {
        self.validate_measured_walls(measured_walls);
        self.collide_into(scratch);
        let post = std::mem::take(scratch);
        let receipt = self.stream_from_inner(&post, Some(measured_walls));
        *scratch = post;
        receipt
    }

    /// One collide-stream step with moving walls and a selected-wall receipt.
    ///
    /// The velocity field, measurement mask, moment origin, and post-collision
    /// state are all admitted before streaming can mutate `self.f`.
    pub fn step_with_moving_wall_momentum(
        &mut self,
        scratch: &mut Vec<[f64; Q]>,
        measured_walls: &[bool],
        wall_velocities: &[[f64; 2]],
        moment_origin: [f64; 2],
    ) -> MovingWallMomentumExchange2 {
        self.validate_measured_walls(measured_walls);
        self.validate_moving_wall_fields(wall_velocities, moment_origin);
        self.collide_into(scratch);
        self.validate_moving_wall_post(scratch, wall_velocities);
        let post = std::mem::take(scratch);
        let receipt = self.stream_from_moving_wall_inner(
            &post,
            measured_walls,
            wall_velocities,
            moment_origin,
        );
        *scratch = post;
        receipt
    }

    /// One collide-stream step with checked linear off-lattice moving walls.
    ///
    /// The measurement mask, complete link geometry, link velocities, and
    /// moment origin are admitted before collision. The resulting
    /// post-collision populations are admitted before grid publication.
    pub fn step_with_curved_moving_wall_momentum(
        &mut self,
        scratch: &mut Vec<[f64; Q]>,
        measured_walls: &[bool],
        wall_links: &[[Option<CurvedMovingWallLink2>; Q]],
        moment_origin: [f64; 2],
    ) -> MovingWallMomentumExchange2 {
        self.validate_measured_walls(measured_walls);
        self.validate_curved_moving_wall_fields(wall_links, moment_origin);
        self.collide_into(scratch);
        self.validate_curved_moving_wall_post(scratch, wall_links);
        let post = std::mem::take(scratch);
        let receipt = self.stream_from_curved_moving_wall_inner(
            &post,
            measured_walls,
            wall_links,
            moment_origin,
        );
        *scratch = post;
        receipt
    }

    fn step_velocity_pressure_x_inner(
        &mut self,
        scratch: &mut Vec<[f64; Q]>,
        boundary: VelocityPressureX2,
        measured_walls: Option<&[bool]>,
    ) -> MomentumExchange2 {
        self.collide_into(scratch);
        let post = std::mem::take(scratch);
        let receipt = self.stream_from_inner(&post, measured_walls);
        self.apply_velocity_pressure_x(boundary);
        *scratch = post;
        receipt
    }

    /// One collide-stream step followed by regularized x-face reconstruction.
    pub fn step_velocity_pressure_x(
        &mut self,
        scratch: &mut Vec<[f64; Q]>,
        boundary: VelocityPressureX2,
    ) {
        self.validate_velocity_pressure_x();
        let _ = self.step_velocity_pressure_x_inner(scratch, boundary, None);
    }

    /// Regularized x-flow step plus a selected-wall momentum receipt.
    ///
    /// Grid topology, periodicity, forcing, and the wall mask are validated
    /// before collision, so a refused request cannot partially advance state.
    pub fn step_velocity_pressure_x_with_wall_momentum(
        &mut self,
        scratch: &mut Vec<[f64; Q]>,
        boundary: VelocityPressureX2,
        measured_walls: &[bool],
    ) -> MomentumExchange2 {
        self.validate_measured_walls(measured_walls);
        self.validate_velocity_pressure_x();
        self.step_velocity_pressure_x_inner(scratch, boundary, Some(measured_walls))
    }
}

fn regularized_populations(rho: f64, velocity: [f64; 2], stress: [[f64; 2]; 2]) -> [f64; Q] {
    let mut populations = crate::equilibrium(rho, velocity[0], velocity[1]);
    let coefficient = 1.0 / (2.0 * CS2 * CS2);
    for q in 0..Q {
        let e = [f64::from(E[q].0), f64::from(E[q].1)];
        let mut contraction = 0.0;
        for row in 0..2 {
            for column in 0..2 {
                let isotropic = if row == column { CS2 } else { 0.0 };
                contraction += (e[row] * e[column] - isotropic) * stress[row][column];
            }
        }
        populations[q] += coefficient * W[q] * contraction;
    }
    populations
}

/// Strain-rate magnitude (sqrt(2 S:S)) of one cell from its
/// non-equilibrium populations — the LOCAL quantity non-Newtonian
/// relaxation adapts to. `feq` must match the cell's moments.
#[must_use]
pub fn shear_rate(f: &[f64; Q], feq: &[f64; Q], rho: f64, tau: f64) -> f64 {
    // S_ab = −(3 / (2 ρ τ)) Σ_q e_qa e_qb (f_q − feq_q)   (c_s² = 1/3)
    let mut sxx = 0.0f64;
    let mut sxy = 0.0f64;
    let mut syy = 0.0f64;
    for q in 0..Q {
        let neq = f[q] - feq[q];
        let (ex, ey) = (f64::from(E[q].0), f64::from(E[q].1));
        sxx += ex * ex * neq;
        sxy += ex * ey * neq;
        syy += ey * ey * neq;
    }
    let c = -3.0 / (2.0 * rho * tau);
    let (sxx, sxy, syy) = (c * sxx, c * sxy, c * syy);
    let ss = 2.0f64.mul_add(sxy * sxy, sxx.mul_add(sxx, syy * syy));
    (2.0 * ss).sqrt()
}
