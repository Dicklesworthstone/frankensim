//! Viscous wall drag of the normal component in a reduced channel model.
//!
//! With stationary walls, scalar permeability `k` and dynamic shear viscosity
//! `eta`, the force per FLUID volume is `F_n = -eta*v_n/k`, while `F_s = 0`.
//! The normal inertia is `rho_n`, not total density or counterflow inertia.
//! Wall impulse is external to the two-fluid system; this is not mutual
//! friction. Its lost normal kinetic energy is exported for one thermal owner.
//!
//! This is a quasisteady, laminar, fully developed transverse-flow closure.
//! Circular/slit helpers use the exact steady Poiseuille resistance. Their
//! transient source is exact for this reduced ODE, not a resolved transverse
//! velocity profile. The caller must establish small normal Reynolds and
//! Womersley numbers, a continuum/no-slip regime, and negligible entrance and
//! profile-distortion effects. Thin oscillatory boundary layers, nonlinear
//! wall turbulence, longitudinal and second viscosities, porous tortuosity,
//! moving walls and evolving material properties require other models.
//!
//! Physics references: Kekutia and Chkhaidze,
//! <https://arxiv.org/abs/cond-mat/0604193>, section 3, derives normal-component
//! wall drag from a quadratic dissipation function and distinguishes static
//! permeability from the frequency-dependent correction. Baggaley and Laizet,
//! DOI: 10.1063/1.4828892, explicitly distinguish laminar Poiseuille normal flow
//! from turbulent normal profiles in counterflow. Neither source supplies a
//! universal permeability or validates a supplied geometry and operating point.

use super::{
    TwoFluidState, TwoFluidVelocity, finite_vector, invalid, nonnegative, positive, speed,
};
use crate::MaterialError;
use fs_math::det;

/// Admitted scalar permeability for normal-fluid drag against fixed walls.
///
/// Velocities and forces are averaged over the FLUID cross-section/volume.
/// A bulk porous-medium permeability must first be converted to this basis;
/// porosity and tortuosity are not silently set to one. For circular and slit
/// helpers, the velocity must be parallel to the channel. The scalar source
/// acts on each vector component with the same declared permeability.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NormalViscousWall {
    permeability_m2: f64,
}

impl NormalViscousWall {
    /// Authored permeability [m²], finite and strictly positive.
    pub fn try_new(permeability_m2: f64) -> Result<Self, MaterialError> {
        if !positive(permeability_m2) {
            return Err(invalid(
                "normal wall drag requires positive finite permeability per fluid cross-section",
            ));
        }
        Ok(Self { permeability_m2 })
    }

    /// Straight circular capillary with no-slip walls: `k = radius²/8`.
    /// Averaging `v_n(r)=(-dp/dx)*(R²-r²)/(4*eta)` gives this resistance.
    pub fn circular_capillary(radius_m: f64) -> Result<Self, MaterialError> {
        if !positive(radius_m) {
            return Err(invalid(
                "normal-fluid capillary radius must be positive finite",
            ));
        }
        // Divide before multiplying so a representable permeability does not
        // overflow merely because the unscaled radius square would overflow.
        Self::try_new((radius_m / 8.0) * radius_m)
    }

    /// Straight slit between no-slip parallel plates: `k = full_gap²/12`.
    /// The width is assumed much greater than the full gap supplied here.
    pub fn parallel_plates(full_gap_m: f64) -> Result<Self, MaterialError> {
        if !positive(full_gap_m) {
            return Err(invalid("normal-fluid plate gap must be positive finite"));
        }
        Self::try_new((full_gap_m / 12.0) * full_gap_m)
    }

    /// Authored permeability [m²] on the fluid-volume basis.
    pub const fn permeability_m2(&self) -> f64 {
        self.permeability_m2
    }
}

/// Exact stationary-wall source result for one local volume.
///
/// The superfluid velocity is bitwise unchanged. A spatial owner multiplies
/// impulse and heat densities by its dual volume, routes the opposite wall
/// reaction to its mechanical owner, and deposits the exported heat once.
/// No temperature or thermal state is modified here.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NormalWallStep {
    pub velocity: TwoFluidVelocity,
    /// Impulse on the normal component, hence on the whole fluid [kg/(m² s)].
    pub normal_impulse_kg_m2_s: [f64; 3],
    /// Exact opposite impulse on the stationary wall [kg/(m² s)].
    pub wall_impulse_kg_m2_s: [f64; 3],
    /// Lost normal kinetic energy [J/m³], nonnegative.
    pub dissipated_energy_j_m3: f64,
}

/// Prepared exact source factors at one frozen material state and duration.
///
/// Reuse this immutable value over a spatial batch: admission and exponential
/// evaluation happen once, and [`Self::apply`] allocates nothing. Its inertia
/// remains the normal density of the state used to prepare it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NormalWallRelaxation {
    normal_density_kg_m3: f64,
    decay_rate_s_inv: f64,
    duration_s: f64,
    attenuation: f64,
    velocity_loss: f64,
    energy_loss: f64,
}

impl NormalWallRelaxation {
    pub const fn normal_density_kg_m3(&self) -> f64 {
        self.normal_density_kg_m3
    }

    /// `eta/(rho_n*k)` [1/s]. Zero viscosity gives exactly zero rate.
    pub const fn decay_rate_s_inv(&self) -> f64 {
        self.decay_rate_s_inv
    }

    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }

    /// Advance `rho_n*dv_n/dt = -eta*v_n/k`; leave `v_s` unchanged.
    ///
    /// Uses `expm1`-prepared loss for very short intervals, preserving impulse
    /// and heat below one velocity ulp. Large decay uses the exponential tail
    /// directly rather than subtracting nearly equal velocities. The exact
    /// local ODE has no explicit stability step restriction; a spatial split
    /// composition must still control its own time-discretization error.
    pub fn apply(&self, velocity: TwoFluidVelocity) -> Result<NormalWallStep, MaterialError> {
        finite_vector(velocity.normal_m_s)?;
        finite_vector(velocity.superfluid_m_s)?;
        if self.velocity_loss == 0.0 {
            return Ok(NormalWallStep {
                velocity,
                normal_impulse_kg_m2_s: [0.0; 3],
                wall_impulse_kg_m2_s: [0.0; 3],
                dissipated_energy_j_m3: 0.0,
            });
        }
        let normal_speed = speed(velocity.normal_m_s);
        let normal_energy = 0.5 * self.normal_density_kg_m3 * normal_speed * normal_speed;
        if !nonnegative(normal_energy) {
            return Err(invalid(
                "normal-fluid kinetic energy density is unrepresentable",
            ));
        }
        let impulse = finite_vector(velocity.normal_m_s.map(|v| {
            // Form the unsuppressed momentum first when representable, so a
            // tiny loss cannot underflow before multiplication by a large v.
            let momentum = self.normal_density_kg_m3 * v;
            if momentum.is_finite() {
                -momentum * self.velocity_loss
            } else {
                -(self.normal_density_kg_m3 * self.velocity_loss) * v
            }
        }))?;
        let next = TwoFluidVelocity {
            normal_m_s: velocity.normal_m_s.map(|v| {
                if self.velocity_loss <= 0.5 {
                    v - self.velocity_loss * v
                } else {
                    self.attenuation * v
                }
            }),
            superfluid_m_s: velocity.superfluid_m_s,
        };
        finite_vector(next.normal_m_s)?;
        let dissipated_energy_j_m3 = normal_energy * self.energy_loss;
        if !nonnegative(dissipated_energy_j_m3) {
            return Err(invalid("normal wall heat density is unrepresentable"));
        }
        Ok(NormalWallStep {
            velocity: next,
            normal_impulse_kg_m2_s: impulse,
            wall_impulse_kg_m2_s: impulse.map(|x| -x),
            dissipated_energy_j_m3,
        })
    }
}

impl TwoFluidState {
    /// Normal decay rate `eta/(rho_n*k)` for the authored stationary walls.
    /// Positive viscosity may not silently underflow to a zero rate.
    pub fn normal_wall_decay_rate_s_inv(
        &self,
        wall: NormalViscousWall,
    ) -> Result<f64, MaterialError> {
        let viscosity = self.normal_viscosity_pa_s();
        if viscosity == 0.0 {
            return Ok(0.0);
        }
        let denominator = self.normal_density_kg_m3() * wall.permeability_m2;
        let rate = if denominator.is_normal() {
            viscosity / denominator
        } else {
            // A subnormal product can lose substantial relative precision;
            // try division in two stages before accepting that denominator.
            let per_density = viscosity / self.normal_density_kg_m3();
            let per_permeability = viscosity / wall.permeability_m2;
            if per_density.is_normal() {
                per_density / wall.permeability_m2
            } else if per_permeability.is_normal() {
                per_permeability / self.normal_density_kg_m3()
            } else {
                viscosity / denominator
            }
        };
        if !positive(rate) {
            return Err(invalid("normal wall decay rate is unrepresentable"));
        }
        Ok(rate)
    }

    /// Prepare one exact normal-wall source for reuse over a spatial batch.
    pub fn prepare_normal_wall_relaxation(
        &self,
        wall: NormalViscousWall,
        dt_s: f64,
    ) -> Result<NormalWallRelaxation, MaterialError> {
        if !nonnegative(dt_s) {
            return Err(invalid(
                "normal wall duration must be finite and nonnegative",
            ));
        }
        let decay_rate_s_inv = self.normal_wall_decay_rate_s_inv(wall)?;
        let exponent = decay_rate_s_inv * dt_s;
        if !nonnegative(exponent) || (decay_rate_s_inv > 0.0 && dt_s > 0.0 && exponent == 0.0) {
            return Err(invalid(
                "normal wall rate times duration is unrepresentable",
            ));
        }
        let velocity_loss = -det::expm1(-exponent);
        let attenuation = if velocity_loss <= 0.5 {
            1.0 - velocity_loss
        } else {
            det::exp(-exponent)
        };
        // 1-exp(-2*x) = (1-exp(-x))*(1+exp(-x)). This retains tiny heat
        // fractions without doubling a possibly very large finite exponent.
        let energy_loss = velocity_loss * (2.0 - velocity_loss);
        Ok(NormalWallRelaxation {
            normal_density_kg_m3: self.normal_density_kg_m3(),
            decay_rate_s_inv,
            duration_s: dt_s,
            attenuation,
            velocity_loss,
            energy_loss,
        })
    }

    /// Exact source convenience entry point; prepare once for spatial batches.
    pub fn relax_normal_wall(
        &self,
        velocity: TwoFluidVelocity,
        wall: NormalViscousWall,
        dt_s: f64,
    ) -> Result<NormalWallStep, MaterialError> {
        self.prepare_normal_wall_relaxation(wall, dt_s)?
            .apply(velocity)
    }
}
