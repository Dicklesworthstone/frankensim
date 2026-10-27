//! Source-bound inertial solids with spatial plastic and reversible heating.
//!
//! The existing thermoplastic owner retains all material claims, original
//! caloric objects, fixed reference mass, P1 displacement and J2 history. Its
//! shared thermal fixed-point now calls the original solid Newton owner with
//! consistent backward-Euler inertia. Every thermal proposal restarts from
//! the same accepted displacement, velocity, plastic history and enthalpy.
//!
//! Actual yield/overstress dissipation and reversible thermoelastic heat reach the SAME
//! thermal timeline once. Constitutive and kinetic numerical losses remain
//! separate from heat. The combined gate includes caloric storage, elastic
//! and hardening free energy, the mechanical entropy contribution, kinetic
//! energy, actual nodal/support work and both numerical losses.
//!
//! This inherits the original source/phase/small-strain admission. The model
//! uses rate-independent J2 by default; explicit material viscosity selects the
//! same source-bound linear-overstress law as the quasistatic thermal owner.
//! Both constitutive kinetics and inertia use the actual physical interval.
//! Loads and support displacements are endpoint values; resolve pulses and
//! refine physical intervals to assess first-order time accuracy. The thermal
//! owner's own substeps do not add mechanical timesteps. Geometry, transport
//! areas and reference mass remain fixed. Melt flow and finite rotations need
//! their existing distinct owners.

use super::*;
pub use fs_solid::plastic3::dynamics::InertialEnergy;

/// Solved-pressure thermal supports on the same retained inertial trajectory.
pub mod contact;
/// Circuit heating and pressure-dependent interfaces with retained inertia.
pub mod electrical;

pub(super) struct Motion {
    pub velocity_mps: Vec<[f64; 3]>,
    pub energy: InertialEnergy,
}

/// One accepted material/thermal/mechanical trajectory, with full nodal velocity.
#[derive(Clone, Debug)]
pub struct ThermoplasticDynamics<'a> {
    solid: ThermoplasticSolid<'a>,
    velocity_mps: Vec<[f64; 3]>,
}

/// Both original physical accounts for the same converged endpoint.
#[derive(Clone, Debug)]
pub struct ThermoplasticDynamicsReport<R = ThermoplasticReport> {
    /// Original thermal, constitutive, external-work and complete combined gates.
    /// Material reports keep constitutive numerical loss separate from inertia.
    pub material: R,
    /// Consistent reference-mass kinetic work and separate inertial BE loss.
    pub inertia: InertialEnergy,
}

pub struct ThermoplasticDynamicsCandidate<'m, 't, O> {
    thermal: DepositionSimulation<'t>,
    dynamics: ThermoplasticDynamics<'m>,
    observer: O,
    report: ThermoplasticDynamicsReport,
}
impl<'m, 't, O> ThermoplasticDynamicsCandidate<'m, 't, O> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn dynamics(&self) -> &ThermoplasticDynamics<'m> {
        &self.dynamics
    }
    pub const fn observer(&self) -> &O {
        &self.observer
    }
    pub const fn report(&self) -> &ThermoplasticDynamicsReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        ThermoplasticDynamics<'m>,
        O,
        ThermoplasticDynamicsReport,
    ) {
        (self.thermal, self.dynamics, self.observer, self.report)
    }
}

impl<'a> ThermoplasticDynamics<'a> {
    /// Attach actual full nodal velocity to an admitted material solid at its
    /// current clock. This establishes initial kinetic energy; it does not
    /// claim that a prior loading or acceleration process was simulated.
    pub fn try_new(
        solid: ThermoplasticSolid<'a>,
        velocity_mps: Vec<[f64; 3]>,
        cx: &Cx<'_>,
    ) -> Result<Self, ThermoplasticError> {
        poll(cx)?;
        if velocity_mps.len() != solid.nodes_m.len()
            || velocity_mps.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(ThermoplasticError::Invalid(
                "one finite velocity per native structural node required",
            ));
        }
        Ok(Self {
            solid,
            velocity_mps,
        })
    }

    pub const fn solid(&self) -> &ThermoplasticSolid<'a> {
        &self.solid
    }
    pub fn velocity_mps(&self) -> &[[f64; 3]] {
        &self.velocity_mps
    }
    pub const fn time_s(&self) -> f64 {
        self.solid.time_s()
    }

    /// Advance one physical interval through the shared Newton/thermal path.
    /// Force and support work include inertia. An unrestrained solid is valid;
    /// no artificial supports or dummy mass are inserted. All thermal/material
    /// errors, observer errors and cancellation leave every input unchanged.
    /// Accept thermal state, motion and returned observer together.
    #[allow(clippy::too_many_arguments)]
    pub fn advance<'t, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        observer: &O,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        thermal_config: StepConfig,
        config: ThermoplasticConfig,
        cx: &Cx<'_>,
    ) -> Result<ThermoplasticDynamicsCandidate<'a, 't, O>, ThermoplasticError<O::Error>> {
        let (candidate, motion) = self.solid.advance_with_motion(
            thermal,
            observer,
            target_s,
            prescribed_m,
            nodal_forces_n,
            settings,
            thermal_config,
            config,
            Some(&self.velocity_mps),
            cx,
        )?;
        let motion = motion.ok_or(ThermoplasticError::Invalid(
            "inertial thermoplastic candidate must retain its native kinetic account",
        ))?;
        let (thermal, solid, observer, material) = candidate.into_parts();
        poll(cx)?;
        Ok(ThermoplasticDynamicsCandidate {
            thermal,
            dynamics: Self {
                solid,
                velocity_mps: motion.velocity_mps,
            },
            observer,
            report: ThermoplasticDynamicsReport {
                material,
                inertia: motion.energy,
            },
        })
    }
}
