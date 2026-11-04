//! Finite physical nozzle transport on the original variable-inertia trajectory.
//!
//! The native step owns mass, orientation, recoil, burnout and its work budget.
//! This observer integrates physical particle moments on those accepted stages.
//! For balanced ports the local-to-inertial energy correction is
//! q*(|V|²/2 + V dot R*pi + omega dot K*omega/2 + omega dot j).
//! Local relative energy already includes relative jet kinetic energy and its
//! caller-owned caloric datum. It is never replaced by effective pressure thrust
//! or supplemented with a second combustion/propellant energy model.

use super::super::super::{RigidRocketError, RigidRocketSolve, ScalarSum, VectorSum, poll};
use super::{PreparedNozzleArray, TensorRocketState, TensorRocketStep};
use crate::Vec3;
use crate::rocket::exhaust::velocity_moments;
use crate::rocket::{RocketForcing, mass_integrals};
use crate::spatial_loads::PreparedSpatialLoads;

/// Held physical nozzles and the original complete local particle-energy datum.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NozzleTransportForcing {
    /// The same prepared mass, particle momentum and co-rotation moments used
    /// by native mechanics. Pressure-resultant loads remain separate.
    pub nozzles: PreparedNozzleArray,
    /// Complete energy rate relative to EACH rotating material nozzle [W],
    /// summed over the array. Includes relative kinetic energy exactly once,
    /// plus internal/random energy in the source's original caloric datum.
    /// A signed formation datum is permitted. Zero mass rate requires zero.
    pub local_relative_energy_rate_w: f64,
    /// Uniform acceleration of retained mass, including the coast [m/s²].
    pub gravity_world_m_s2: Vec3,
}

/// Native motion and actual interval fluxes, all in one declared world frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NozzleTransportStep {
    /// Original trajectory and canonical emitted mass/burnout/substep account.
    /// Its exported spin is the same native angular outflow integral used here.
    pub mechanics: TensorRocketStep,
    /// Complete caller-supplied local relative energy during emission [J].
    pub local_relative_energy_j: f64,
    /// Co-rotation and tangential-jet correction q*(omega.K.omega/2+omega.j)
    /// integrated with the native midpoint attitude/inertia [J]. May be signed.
    pub co_rotating_transport_energy_j: f64,
    /// COM frame correction, integrating q*|V|²/2 + q*V.R*pi along the
    /// analytic changing COM velocity in every native stage [J]. May be signed.
    pub translational_transport_energy_j: f64,
    /// Local relative energy plus BOTH frame corrections [J]. The composing
    /// material owner must account for this export and its own stored energy;
    /// this is not a claim that engine, tank or thermal energy has been closed.
    pub exhaust_energy_world_j: f64,
    /// Physical emitted momentum, including mass advection [kg m/s].
    pub exhaust_linear_momentum_world_kg_m_s: Vec3,
    /// Orbital part about the supplied fixed world origin [kg m²/s].
    pub exhaust_orbital_angular_momentum_world_kg_m2_s: Vec3,
    /// Orbital plus the original native spin-flow integral [kg m²/s].
    pub exhaust_angular_momentum_world_kg_m2_s: Vec3,
    /// Held external force and changing retained-mass gravity impulse [N s].
    pub external_linear_impulse_world_n_s: Vec3,
    /// External force/gravity orbital impulse plus native applied spin impulse,
    /// about the same fixed world origin [N m s]. Includes all coasting stages.
    pub external_angular_impulse_world_n_m_s: Vec3,
    /// Actual retained momentum change + physical export - external impulse.
    /// Computed in an initial-velocity chart to retain small recoil defects.
    pub linear_momentum_residual_world_kg_m_s: Vec3,
    /// Actual spin/orbital change + physical export - external angular impulse.
    /// Includes quadrature and represented endpoint defects without repairing
    /// the export to cancel them [kg m²/s].
    pub angular_momentum_residual_world_kg_m2_s: Vec3,
}

impl TensorRocketState {
    /// Advance and integrate finite physical nozzle exchange in one transaction.
    ///
    /// Every observation comes from the ORIGINAL bounded native step. Its
    /// endpoint is identical to `advance_spatial` with `forcing.nozzles.forcing`
    /// and the same loads/controls. There is no second attitude integration,
    /// per-nozzle traversal, stored stage history, or post-burnout emission.
    ///
    /// Native analytic mass/velocity moments integrate linear transport and the
    /// translational energy correction exactly within each frozen-mid-attitude
    /// stage. Co-rotation uses the native midpoint spin and current full tensor.
    /// Two-point Gauss quadrature samples the SAME analytic COM path for orbital
    /// export and external orbital impulse, including coast. Spin export uses
    /// the original native outflow kicks. With rotation the full composition is
    /// second order; temporal refinement remains required, especially at rapid
    /// depletion. Momentum residuals retain that discretization error.
    ///
    /// Local relative energy, co-rotation and COM corrections stay separate so
    /// an arbitrary formation datum cannot erase the mechanical transport.
    /// Additional pressure loads may be supplied through `loads`, but must not
    /// be folded into particle velocity or the relative energy a second time.
    /// Internal flow/slosh, shifting COM and engine work remain caller-owned.
    ///
    /// Storage and extra work per native substep are constant. Original work
    /// limits and cancellation apply; no partial state or stream is published.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance_nozzles_transport(
        &self,
        forcing: NozzleTransportForcing,
        loads: &PreparedSpatialLoads,
        origin_world_m: Vec3,
        duration_s: f64,
        controls: RigidRocketSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NozzleTransportStep, RigidRocketError> {
        poll(&mut cancelled)?;
        let nozzles = forcing.nozzles;
        if !forcing.local_relative_energy_rate_w.is_finite()
            || !origin_world_m.is_finite()
            || (nozzles.mass_flow_kg_s == 0.0 && forcing.local_relative_energy_rate_w != 0.0)
        {
            return Err(RigidRocketError::Invalid(
                "finite nozzle energy/origin and closed-valve energy required",
            ));
        }
        let boost = self.translation().velocity_m_s;
        let origin_arm = self.translation().position_m.sub(origin_world_m);
        let force = loads.force_world();
        let gravity = forcing.gravity_world_m_s2;
        let (mut relative_energy, mut corotating_energy, mut centered_energy) = (
            ScalarSum::default(),
            ScalarSum::default(),
            ScalarSum::default(),
        );
        let (mut centered_momentum, mut orbital, mut external_linear, mut external_orbital) = (
            VectorSum::default(),
            VectorSum::default(),
            VectorSum::default(),
            VectorSum::default(),
        );
        let mechanics = self.advance_spatial_observed(
            nozzles.forcing(gravity),
            loads,
            duration_s,
            controls,
            |stage, cancelled| {
                let h = stage.duration_s;
                let dm = stage.motion.burned_propellant_mass_kg;
                let m0 = stage.before.dry_mass_kg + stage.before.propellant_mass_kg;
                let m1 = stage.motion.state_after.dry_mass_kg
                    + stage.motion.state_after.propellant_mass_kg;
                let relative_velocity = stage
                    .midpoint_orientation
                    .rotate_body_to_world(nozzles.mean_exhaust_velocity_body_m_s);
                let external = force.scale(h).add(gravity.scale(h).scale(m0 - 0.5 * dm));
                let stage_arm = stage.before.position_m.sub(self.translation().position_m);
                external_linear.add(external);
                external_orbital.add(stage_arm.cross(external));
                if dm > 0.0 {
                    let x = dm / m0;
                    let (_, a, c) = mass_integrals(x, m0, m1);
                    let d = force.scale(h / m0).sub(relative_velocity.scale(x));
                    let (mean_increment, variance, _) =
                        velocity_moments(x, a, c, m1 / m0, d, gravity.scale(h))?;
                    let mean = stage.before.velocity_m_s.sub(boost).add(mean_increment);
                    let relative_momentum = relative_velocity.scale(dm);
                    let momentum = mean.scale(dm).add(relative_momentum);
                    centered_momentum.add(momentum);
                    orbital.add(stage_arm.cross(momentum.add(boost.scale(dm))));
                    centered_energy.add(
                        (0.5 * dm) * (mean.norm_squared() + variance) + mean.dot(relative_momentum),
                    );
                    relative_energy.add(
                        forcing.local_relative_energy_rate_w * stage.motion.powered_duration_s,
                    );
                    let omega = stage
                        .midpoint_properties
                        .angular_velocity_body_checked(stage.midpoint_momentum_body)?;
                    let k = nozzles.co_rotating_specific_inertia_m2.tensor_body_m2().m;
                    let k_omega = Vec3::new(
                        k[0].mul_add(omega.x, k[1].mul_add(omega.y, k[2] * omega.z)),
                        k[3].mul_add(omega.x, k[4].mul_add(omega.y, k[5] * omega.z)),
                        k[6].mul_add(omega.x, k[7].mul_add(omega.y, k[8] * omega.z)),
                    );
                    corotating_energy.add(
                        dm * omega.dot(
                            k_omega
                                .scale(0.5)
                                .add(nozzles.relative_specific_angular_momentum_body_m2_s),
                        ),
                    );
                }
                // Two original analytic translations, never extra rigid stages.
                // Integrate only the displacement from this stage's origin;
                // the lever arm above uses its EXACT analytic linear impulse.
                // A fixed-world origin shift is applied once after reduction,
                // so it cannot amplify orbital quadrature error.
                let orbit_active = dm > 0.0 || force != Vec3::ZERO || gravity != Vec3::ZERO;
                for fraction in [0.211_324_865_405_187_13, 0.788_675_134_594_812_9]
                    .into_iter()
                    .filter(|_| orbit_active)
                {
                    let sample = stage
                        .before
                        .advance(
                            RocketForcing {
                                mass_flow_kg_s: nozzles.mass_flow_kg_s,
                                effective_exhaust_velocity_m_s: relative_velocity.scale(-1.0),
                                gravity_m_s2: gravity,
                                external_force_n: force,
                            },
                            h * fraction,
                            &mut *cancelled,
                        )?
                        .state_after;
                    let arm = sample.position_m.sub(stage.before.position_m);
                    if dm > 0.0 {
                        orbital.add(
                            arm.cross(sample.velocity_m_s.add(relative_velocity).scale(0.5 * dm)),
                        );
                    }
                    let mass = sample.dry_mass_kg + sample.propellant_mass_kg;
                    external_orbital.add(arm.cross(force.add(gravity.scale(mass)).scale(0.5 * h)));
                }
                Ok(())
            },
            &mut cancelled,
        )?;
        let dm = mechanics.burned_propellant_mass_kg;
        let after = mechanics.state_after.translation();
        let m1 = after.dry_mass_kg + after.propellant_mass_kg;
        let actual_mass_loss = self.translation().propellant_mass_kg - after.propellant_mass_kg;
        let actual_velocity_change = after.velocity_m_s.sub(boost);
        let linear = boost.scale(dm).add(centered_momentum.value);
        let translation_energy = if dm == 0.0 {
            0.0
        } else {
            centered_energy.value
                + boost.dot(centered_momentum.value)
                + (0.5 * dm) * boost.norm_squared()
        };
        let linear_residual = actual_velocity_change
            .scale(m1)
            .add(centered_momentum.value)
            .sub(external_linear.value)
            .add(boost.scale(dm - actual_mass_loss));
        let orbital_change = after
            .position_m
            .sub(self.translation().position_m)
            .cross(after.velocity_m_s.scale(m1));
        let angular_residual = mechanics
            .angular_momentum_residual_world_kg_m2_s
            .add(orbital_change)
            .add(orbital.value)
            .sub(external_orbital.value)
            .add(origin_arm.cross(linear_residual));
        let orbital_export = orbital.value.add(origin_arm.cross(linear));
        let angular = orbital_export.add(mechanics.exported_angular_momentum_world_kg_m2_s);
        let external_angular = external_orbital
            .value
            .add(origin_arm.cross(external_linear.value))
            .add(mechanics.applied_angular_impulse_world_kg_m2_s);
        let energy = relative_energy.value + corotating_energy.value + translation_energy;
        if [
            relative_energy.value,
            corotating_energy.value,
            translation_energy,
            energy,
        ]
        .into_iter()
        .any(|value| !value.is_finite())
            || [
                linear,
                orbital_export,
                angular,
                external_linear.value,
                external_angular,
                linear_residual,
                angular_residual,
            ]
            .into_iter()
            .any(|value| !value.is_finite())
        {
            return Err(RigidRocketError::Invalid(
                "finite nozzle transport overflow",
            ));
        }
        poll(&mut cancelled)?;
        Ok(NozzleTransportStep {
            mechanics,
            local_relative_energy_j: relative_energy.value,
            co_rotating_transport_energy_j: corotating_energy.value,
            translational_transport_energy_j: translation_energy,
            exhaust_energy_world_j: energy,
            exhaust_linear_momentum_world_kg_m_s: linear,
            exhaust_orbital_angular_momentum_world_kg_m2_s: orbital_export,
            exhaust_angular_momentum_world_kg_m2_s: angular,
            external_linear_impulse_world_n_s: external_linear.value,
            external_angular_impulse_world_n_m_s: external_angular,
            linear_momentum_residual_world_kg_m_s: linear_residual,
            angular_momentum_residual_world_kg_m2_s: angular_residual,
        })
    }
}

#[cfg(test)]
#[path = "transport/tests.rs"]
mod tests;
