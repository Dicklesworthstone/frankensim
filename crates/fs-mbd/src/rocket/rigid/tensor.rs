//! Full-tensor depletion and angular outflow in fixed material body coordinates.
//!
//! Dry inertia, propellant inertia per mass and exhaust inertia per mass need
//! not commute or share principal axes. The body frame never changes when the
//! principal axes of I_dry + m K_propellant move. Native rigid rotation owns
//! the gyroscopic term; the existing rocket owner owns burnout and translation.
//!
//! Each frozen-attitude outflow submap solves
//! (I_mid + dm K_exit/2) omega_mid = H_before,
//! H_after = H_before - dm K_exit omega_mid.
//! This is implicit midpoint in angular momentum, including noncommuting full
//! tensors. It is second order for smooth depletion, not an exact matrix power
//! or an L-stable filter. Resolve rapid exhaust damping in the declared step.
//!
//! The fixed-COM model prescribes equivalent thrust, external torque, angular
//! co-rotation flux and swirl independently. It does not derive these inputs
//! from an internal flow field, model slosh or a shifting COM, or account for
//! combustion/plume energy. See Eke, NASA/CR-1998-208246, chapters 2–3:
//! https://ntrs.nasa.gov/citations/19980210404.

use super::{
    RigidRocketError, RigidRocketForcing, RigidRocketState, RigidRocketStep, RocketInertiaModel,
    physical_diagonal, sealed,
};
use crate::{MassProperties, UnitQuaternion, Vec3, rocket::RocketState};
use fs_ga::Mat3;

/// A prepared physically admissible specific rotational inertia, in m².
///
/// Principal moments may be zero, admitting a point at the COM or an ideal line
/// distribution. Their triangle inequalities are checked before the tensor is
/// expressed in material body axes. This same type describes either the
/// retained propellant distribution or an independently declared outflow law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpecificInertia {
    tensor_body_m2: Mat3,
}

impl SpecificInertia {
    /// Zero inertia, including propellant concentrated at the common COM.
    pub const ZERO: Self = Self {
        tensor_body_m2: Mat3::zero(),
    };

    /// Prepare nonnegative principal moments in arbitrary material body axes.
    pub fn from_principal_frame(
        principal_m2: Vec3,
        principal_to_body: UnitQuaternion,
    ) -> Result<Self, RigidRocketError> {
        physical_diagonal(principal_m2, "specific inertia principal moments")?;
        let axes = [
            principal_to_body.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
            principal_to_body.rotate_body_to_world(Vec3::new(0.0, 1.0, 0.0)),
            principal_to_body.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0)),
        ];
        let values = [principal_m2.x, principal_m2.y, principal_m2.z];
        let mut tensor = Mat3::zero();
        for (axis, value) in axes.into_iter().zip(values) {
            let v = [axis.x, axis.y, axis.z];
            for row in 0..3 {
                for col in row..3 {
                    let entry = (value * v[row]).mul_add(v[col], tensor.m[3 * row + col]);
                    tensor.m[3 * row + col] = entry;
                    tensor.m[3 * col + row] = entry;
                }
            }
        }
        if tensor.m.iter().any(|value| !value.is_finite()) {
            return Err(RigidRocketError::Invalid("specific inertia overflow"));
        }
        Ok(Self {
            tensor_body_m2: tensor,
        })
    }

    /// Admit a positive-definite full tensor through the original inertia owner.
    /// Use [`Self::from_principal_frame`] for semidefinite distributions.
    pub fn from_inertia_tensor(tensor_body_m2: Mat3) -> Result<Self, RigidRocketError> {
        let properties = MassProperties::from_inertia_tensor(1.0, Vec3::ZERO, tensor_body_m2)?;
        Ok(Self {
            tensor_body_m2: properties.inertia_tensor_body(),
        })
    }

    /// Actual material-body-frame tensor, in m²; no diagonal approximation.
    #[must_use]
    pub const fn tensor_body_m2(self) -> Mat3 {
        self.tensor_body_m2
    }
}

/// A fixed common COM with independently oriented dry and propellant inertias.
///
/// I(m) = I_dry + m K_propellant. The dry tensor is positive definite and both
/// distributions satisfy the rigid-body inertia inequalities. Principal axes
/// may move as propellant is consumed; orientation and momenta stay in the same
/// material body coordinates. Only native 3×3 preparation is needed per changed
/// inventory; no body-sized dense matrix or heap workspace is introduced.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TensorRocketInertia {
    dry: MassProperties,
    propellant: SpecificInertia,
}

impl TensorRocketInertia {
    /// Prepare a dry full inertia [kg m²] and specific propellant tensor [m²].
    /// The dry mass remains solely in the original [`RocketState`] inventory.
    pub fn new(
        dry_inertia_body_kg_m2: Mat3,
        propellant: SpecificInertia,
    ) -> Result<Self, RigidRocketError> {
        Ok(Self {
            dry: MassProperties::from_inertia_tensor(1.0, Vec3::ZERO, dry_inertia_body_kg_m2)?,
            propellant,
        })
    }

    /// Prepared dry full inertia in the fixed material body axes, in kg m².
    #[must_use]
    pub fn dry_inertia_body_kg_m2(self) -> Mat3 {
        self.dry.inertia_tensor_body()
    }

    /// Prepared propellant distribution per unit remaining mass, in m².
    #[must_use]
    pub const fn propellant_specific_inertia(self) -> SpecificInertia {
        self.propellant
    }

    /// Derive native full mass properties at the original propellant inventory.
    pub fn mass_properties(
        self,
        translation: RocketState,
    ) -> Result<MassProperties, RigidRocketError> {
        if !translation.dry_mass_kg.is_finite()
            || translation.dry_mass_kg <= 0.0
            || !translation.propellant_mass_kg.is_finite()
            || translation.propellant_mass_kg < 0.0
        {
            return Err(RigidRocketError::Invalid(
                "positive dry and nonnegative propellant masses required",
            ));
        }
        let mass = translation.dry_mass_kg + translation.propellant_mass_kg;
        if translation.propellant_mass_kg == 0.0 || self.propellant == SpecificInertia::ZERO {
            return Ok(MassProperties::from_principal_frame(
                mass,
                Vec3::ZERO,
                self.dry.principal_inertia_body(),
                self.dry.principal_axes_to_body(),
            )?);
        }
        let mut inertia = self.dry.inertia_tensor_body();
        for (value, specific) in inertia.m.iter_mut().zip(self.propellant.tensor_body_m2.m) {
            *value = translation.propellant_mass_kg.mul_add(specific, *value);
        }
        Ok(MassProperties::from_inertia_tensor(
            mass,
            Vec3::ZERO,
            inertia,
        )?)
    }
}

/// The original rocket trajectory owner with full-tensor mass depletion.
pub type TensorRocketState = RigidRocketState<TensorRocketInertia>;
/// Original propulsion/force inputs with a prepared full angular outflow tensor.
pub type TensorRocketForcing = RigidRocketForcing<SpecificInertia>;
/// Original burnout, thrust and world-spin account for the tensor trajectory.
pub type TensorRocketStep = RigidRocketStep<TensorRocketInertia>;

impl sealed::Sealed for TensorRocketInertia {}
impl RocketInertiaModel for TensorRocketInertia {
    type ExhaustInertia = SpecificInertia;

    fn mass_properties(self, translation: RocketState) -> Result<MassProperties, RigidRocketError> {
        TensorRocketInertia::mass_properties(self, translation)
    }

    fn validate_exhaust(self, _exhaust: SpecificInertia) -> Result<(), RigidRocketError> {
        // All values are private and admitted by the prepared constructor.
        Ok(())
    }

    fn exhaust_spin_step(
        self,
        momentum: Vec3,
        before: MassProperties,
        after: MassProperties,
        exhaust: SpecificInertia,
        exported_mass: f64,
    ) -> Result<(Vec3, Vec3), RigidRocketError> {
        if !exported_mass.is_finite() || exported_mass < 0.0 || !momentum.is_finite() {
            return Err(RigidRocketError::Invalid(
                "finite angular outflow inputs required",
            ));
        }
        if exported_mass == 0.0 || exhaust == SpecificInertia::ZERO || momentum == Vec3::ZERO {
            return Ok((momentum, Vec3::ZERO));
        }
        let a = before.inertia_tensor_body().m;
        let b = after.inertia_tensor_body().m;
        let e = exhaust.tensor_body_m2.m;
        let mut system = [[0.0; 3]; 3];
        for row in 0..3 {
            for col in 0..3 {
                let i = 3 * row + col;
                system[row][col] = (0.5 * exported_mass).mul_add(e[i], 0.5 * a[i] + 0.5 * b[i]);
            }
        }
        // Reuse the original native rotational 3×3 solve. There is no separate
        // factorizer, nonlinear iteration, diagonal projection or explicit kick.
        let omega = crate::solve_rotation_3(system, [momentum.x, momentum.y, momentum.z])?;
        let exported = Vec3::new(
            e[0].mul_add(omega[0], e[1].mul_add(omega[1], e[2] * omega[2])),
            e[3].mul_add(omega[0], e[4].mul_add(omega[1], e[5] * omega[2])),
            e[6].mul_add(omega[0], e[7].mul_add(omega[1], e[8] * omega[2])),
        )
        .scale(exported_mass);
        let retained = momentum.sub(exported);
        if !retained.is_finite() || !exported.is_finite() {
            return Err(RigidRocketError::Invalid("angular outflow overflow"));
        }
        Ok((retained, exported))
    }
}
