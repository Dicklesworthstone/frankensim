//! Native membrane momentum with retained finite-strain relaxation history.
//!
//! Backward-Euler inertia and the frozen material history enter the same
//! matrix-free Newton solve. Reference areal inventory is unchanged by strain
//! or relaxed thickness; both lumped and exact consistent P1 mass are supported.
//! Physical relaxation heat, signed material endpoint-work error, and numerical
//! kinetic damping are separate accounts. This is a first-order isothermal
//! interval, not a conservative time integrator or a thermal evolution law.

use super::super::{
    MembraneError, MembraneEvaluation, MembraneInertialPotential, MembranePatch, MembraneSolution,
    MembraneSolveSettings, finite_vectors, inertia::MembraneMassImage, poll,
};
use super::{
    FiniteViscoMembrane, FiniteViscoMembraneEnergy, FiniteViscoMembraneError, energy_gate,
};

pub use crate::hyper3::visco::dynamics::FiniteViscoMotionControls;
pub use crate::inertia::InertialEnergy;

/// Accepted material memory, motion and the selected original reference mass.
#[derive(Clone, Debug)]
pub struct FiniteViscoMembraneMotion<'a> {
    membrane: FiniteViscoMembrane<'a>,
    velocity_mps: Vec<[f64; 3]>,
    mass_image: MembraneMassImage,
}

/// Complete proposed interval; dropping it leaves the accepted motion intact.
#[derive(Debug)]
pub struct FiniteViscoMembraneMotionStep<'a> {
    pub motion: FiniteViscoMembraneMotion<'a>,
    /// Includes inertial forces and reactions on moving supports.
    pub equilibrium: MembraneSolution,
    /// Physical branch loss is independent of numerical kinetic damping.
    pub material_energy: FiniteViscoMembraneEnergy,
    /// Physical Maxwell loss per original triangle for this interval [J].
    /// The original increments are retained independently of cumulative history.
    pub element_dissipated_j: Vec<f64>,
    pub inertia: InertialEnergy,
    /// Momentum work residual plus the independent kinetic-work identity error.
    pub balance_error_j: f64,
}

impl<'a> FiniteViscoMembraneMotion<'a> {
    /// Retain an accepted membrane and bind its initial velocity to the original
    /// areal inventory. Consistent mass integrates the full P1 density field;
    /// lumping is an explicit alternative, never an implicit mesh-size fallback.
    pub fn try_new(
        membrane: FiniteViscoMembrane<'a>,
        velocity_mps: &[[f64; 3]],
        mass_image: MembraneMassImage,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, FiniteViscoMembraneError> {
        poll(&mut cancelled)?;
        if !finite_vectors(velocity_mps, membrane.patch.mesh().nodes.len()) {
            return Err(MembraneError::Invalid(
                "one finite initial velocity per original membrane node required",
            )
            .into());
        }
        let momentum =
            membrane
                .patch
                .reference_mass_apply(mass_image, velocity_mps, &mut cancelled)?;
        let kinetic = 0.5 * inner(velocity_mps, &momentum);
        if !kinetic.is_finite() || kinetic < 0.0 {
            return Err(MembraneError::Invalid(
                "finite nonnegative membrane kinetic energy required",
            )
            .into());
        }
        poll(&mut cancelled)?;
        Ok(Self {
            membrane,
            velocity_mps: velocity_mps.to_vec(),
            mass_image,
        })
    }

    pub fn membrane(&self) -> &FiniteViscoMembrane<'a> {
        &self.membrane
    }

    pub fn velocity_mps(&self) -> &[[f64; 3]] {
        &self.velocity_mps
    }

    pub fn mass_image(&self) -> MembraneMassImage {
        self.mass_image
    }

    /// Original reference mass applied to accepted velocity [kg m/s per node].
    /// Retains inter-node momentum when consistent inertia was selected.
    pub fn nodal_momentum_kg_m_s(
        &self,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<[f64; 3]>, FiniteViscoMembraneError> {
        Ok(self.membrane.patch.reference_mass_apply(
            self.mass_image,
            &self.velocity_mps,
            cancelled,
        )?)
    }

    /// Advance one physical interval with held nodal dead loads. Displacement,
    /// velocity, time and every branch state are published together only after
    /// equilibrium, material-work and numerical-kinetic-loss gates succeed.
    /// Prescribed values are endpoint displacements, including moving supports;
    /// their velocity is the actual displacement increment divided by `dt_s`.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        dt_s: f64,
        dead_load_n: &[[f64; 3]],
        prescribed: &[(usize, f64)],
        settings: MembraneSolveSettings,
        controls: FiniteViscoMotionControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<FiniteViscoMembraneMotionStep<'a>, FiniteViscoMembraneError> {
        self.advance_with_pressure(
            dt_s,
            dead_load_n,
            0.0,
            prescribed,
            settings,
            controls,
            cancelled,
        )
    }

    /// Advance with the original membrane's uniform current-area pressure.
    /// Nonzero pressure retains its complete prescribed-perimeter requirement.
    /// The reference mass and material histories remain fixed through all
    /// Newton and backtracking trials; no moving-edge pressure closure is added.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_pressure(
        &self,
        dt_s: f64,
        dead_load_n: &[[f64; 3]],
        pressure_pa: f64,
        prescribed: &[(usize, f64)],
        settings: MembraneSolveSettings,
        controls: FiniteViscoMotionControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FiniteViscoMembraneMotionStep<'a>, FiniteViscoMembraneError> {
        poll(&mut cancelled)?;
        if !controls.maximum_numerical_kinetic_loss_j.is_finite()
            || controls.maximum_numerical_kinetic_loss_j < 0.0
        {
            return Err(MembraneError::Invalid(
                "finite nonnegative membrane numerical kinetic-loss allowance required",
            )
            .into());
        }
        let inertial_potential = MembraneBackwardEuler::new(
            self.membrane.patch,
            self.mass_image,
            &self.membrane.displacement_m,
            &self.velocity_mps,
            dt_s,
        )?;
        let step = self.membrane.advance_loaded(
            dt_s,
            dead_load_n,
            pressure_pa,
            prescribed,
            &self.membrane.displacement_m,
            settings,
            controls.material,
            Some(&inertial_potential),
            &mut cancelled,
        )?;
        let (velocity_mps, inertia) =
            inertial_potential.endpoint(&step.membrane.displacement_m, &mut cancelled)?;
        energy_gate(
            "membrane backward-Euler numerical kinetic loss",
            inertia.numerical_loss_j,
            controls.maximum_numerical_kinetic_loss_j,
        )?;
        let balance_error_j = step.energy.equilibrium_work_residual_j + inertia.balance_error_j;
        energy_gate(
            "membrane viscoelastic inertial work",
            balance_error_j,
            controls.material.equilibrium_work_tolerance_j,
        )?;
        poll(&mut cancelled)?;
        Ok(FiniteViscoMembraneMotionStep {
            motion: Self {
                membrane: step.membrane,
                velocity_mps,
                mass_image: self.mass_image,
            },
            equilibrium: step.equilibrium,
            material_energy: step.energy,
            element_dissipated_j: step.element_dissipated_j,
            inertia,
            balance_error_j,
        })
    }
}

/// Borrow the original membrane geometry, inventory and accepted motion.
/// All Newton trials see exactly the same old displacement and velocity.
struct MembraneBackwardEuler<'a> {
    patch: &'a MembranePatch,
    mass_image: MembraneMassImage,
    old_displacement_m: &'a [[f64; 3]],
    old_velocity_mps: &'a [[f64; 3]],
    inverse_dt: f64,
}

fn inner(a: &[[f64; 3]], b: &[[f64; 3]]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(a, b)| a.iter().zip(b).map(|(a, b)| a * b).sum::<f64>())
        .sum()
}

impl<'a> MembraneBackwardEuler<'a> {
    fn new(
        patch: &'a MembranePatch,
        mass_image: MembraneMassImage,
        old_displacement_m: &'a [[f64; 3]],
        old_velocity_mps: &'a [[f64; 3]],
        dt_s: f64,
    ) -> Result<Self, MembraneError> {
        let inverse_dt = dt_s.recip();
        let coefficient = inverse_dt * inverse_dt;
        let n = patch.mesh().nodes.len();
        if !dt_s.is_finite()
            || dt_s <= 0.0
            || !coefficient.is_finite()
            || coefficient <= 0.0
            || !finite_vectors(old_displacement_m, n)
            || !finite_vectors(old_velocity_mps, n)
        {
            return Err(MembraneError::Invalid(
                "finite membrane motion and representable positive inertial interval required",
            ));
        }
        Ok(Self {
            patch,
            mass_image,
            old_displacement_m,
            old_velocity_mps,
            inverse_dt,
        })
    }

    fn velocities(
        &self,
        displacement_m: &[[f64; 3]],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<(Vec<[f64; 3]>, Vec<[f64; 3]>), MembraneError> {
        poll(&mut *cancelled)?;
        let n = self.old_displacement_m.len();
        if !finite_vectors(displacement_m, n) {
            return Err(MembraneError::Invalid(
                "one finite inertial displacement per original membrane node required",
            ));
        }
        let mut velocity = Vec::with_capacity(n);
        let mut delta = Vec::with_capacity(n);
        for (node, ((u, old_u), old_v)) in displacement_m
            .iter()
            .zip(self.old_displacement_m)
            .zip(self.old_velocity_mps)
            .enumerate()
        {
            if node % 256 == 0 {
                poll(&mut *cancelled)?;
            }
            // Subtract displacement before dividing; reference coordinates and
            // the large rigid offset never enter the physical acceleration.
            let v = core::array::from_fn(|axis| (u[axis] - old_u[axis]) * self.inverse_dt);
            delta.push(core::array::from_fn(|axis| v[axis] - old_v[axis]));
            velocity.push(v);
        }
        if !finite_vectors(&velocity, n) || !finite_vectors(&delta, n) {
            return Err(MembraneError::Invalid(
                "unrepresentable endpoint membrane velocity",
            ));
        }
        Ok((velocity, delta))
    }

    fn endpoint(
        &self,
        displacement_m: &[[f64; 3]],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<(Vec<[f64; 3]>, InertialEnergy), MembraneError> {
        let (velocity, delta) = self.velocities(displacement_m, cancelled)?;
        let old_momentum = self.patch.reference_mass_apply(
            self.mass_image,
            self.old_velocity_mps,
            &mut *cancelled,
        )?;
        let momentum =
            self.patch
                .reference_mass_apply(self.mass_image, &velocity, &mut *cancelled)?;
        let impulse = self
            .patch
            .reference_mass_apply(self.mass_image, &delta, &mut *cancelled)?;
        let kinetic_before_j = 0.5 * inner(self.old_velocity_mps, &old_momentum);
        let kinetic_after_j = 0.5 * inner(&velocity, &momentum);
        let numerical_loss_j = 0.5 * inner(&delta, &impulse);
        let mut kinetic_change_j = 0.0;
        let mut inertial_work_j = 0.0;
        for (node, (((v, old_v), dp), (u, old_u))) in velocity
            .iter()
            .zip(self.old_velocity_mps)
            .zip(&impulse)
            .zip(displacement_m.iter().zip(self.old_displacement_m))
            .enumerate()
        {
            if node % 256 == 0 {
                poll(&mut *cancelled)?;
            }
            for axis in 0..3 {
                kinetic_change_j += (0.5 * v[axis] + 0.5 * old_v[axis]) * dp[axis];
                inertial_work_j += (u[axis] - old_u[axis]) * (dp[axis] * self.inverse_dt);
            }
        }
        let energy = InertialEnergy {
            reference_mass_kg: self.patch.reference_mass_kg(),
            kinetic_before_j,
            kinetic_after_j,
            kinetic_change_j,
            numerical_loss_j,
            inertial_work_j,
            balance_error_j: inertial_work_j - kinetic_change_j - numerical_loss_j,
        };
        if [
            energy.reference_mass_kg,
            kinetic_before_j,
            kinetic_after_j,
            kinetic_change_j,
            numerical_loss_j,
            inertial_work_j,
            energy.balance_error_j,
        ]
        .iter()
        .any(|value| !value.is_finite())
            || kinetic_before_j < 0.0
            || kinetic_after_j < 0.0
            || numerical_loss_j < 0.0
        {
            return Err(MembraneError::Invalid(
                "unrepresentable membrane kinetic-work account",
            ));
        }
        poll(&mut *cancelled)?;
        Ok((velocity, energy))
    }
}

impl MembraneInertialPotential for MembraneBackwardEuler<'_> {
    fn evaluate(
        &self,
        displacement_m: &[[f64; 3]],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<MembraneEvaluation, MembraneError> {
        let (_, delta) = self.velocities(displacement_m, cancelled)?;
        let mut impulse =
            self.patch
                .reference_mass_apply(self.mass_image, &delta, &mut *cancelled)?;
        let potential = 0.5 * inner(&delta, &impulse);
        for (node, value) in impulse.iter_mut().enumerate() {
            if node % 256 == 0 {
                poll(&mut *cancelled)?;
            }
            for component in value {
                *component *= self.inverse_dt;
            }
        }
        if !potential.is_finite() || potential < 0.0 || !finite_vectors(&impulse, delta.len()) {
            return Err(MembraneError::Invalid(
                "unrepresentable membrane inertial potential or force",
            ));
        }
        Ok(MembraneEvaluation {
            strain_energy_j: potential,
            internal_gradient_n: impulse,
        })
    }

    fn tangent_apply(
        &self,
        direction_m: &[[f64; 3]],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Vec<[f64; 3]>, MembraneError> {
        let mut output =
            self.patch
                .reference_mass_apply(self.mass_image, direction_m, &mut *cancelled)?;
        let coefficient = self.inverse_dt * self.inverse_dt;
        for (node, value) in output.iter_mut().enumerate() {
            if node % 256 == 0 {
                poll(&mut *cancelled)?;
            }
            for component in value {
                *component *= coefficient;
            }
        }
        if !finite_vectors(&output, direction_m.len()) {
            return Err(MembraneError::Invalid(
                "unrepresentable membrane inertial tangent",
            ));
        }
        Ok(output)
    }
}
