//! Original-source Galilean traces and laboratory moving-face exchanges.

use super::{
    Conserved, FaceFlux, Primitive, RealGas, RealGasBoundary, RealGasFlowError, RealGasWallControl,
    kinematics, poll,
};
use fs_qty::Velocity;

impl RealGas {
    /// Original-source HLLC/HLL flux through a prescribed translating face.
    ///
    /// Returns laboratory mass/momentum/energy transport F-wU; its signal speed
    /// is RELATIVE to the face. Both traces and every numerical star state use
    /// this owner's original density-dependent calorics and one energy gauge.
    /// Relative kinetic energy is reconstructed from original rho/T, avoiding
    /// subtraction of large laboratory kinetic energies. No mesh is advanced.
    pub fn moving_hllc_flux(
        &self,
        left: Conserved,
        right: Conserved,
        face_velocity: Velocity,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, RealGasFlowError> {
        let pl = self.decode(left, &mut cancelled)?;
        let pr = self.decode(right, &mut cancelled)?;
        interface(
            self,
            (left, pl),
            (right, pr),
            face_velocity.value(),
            &mut cancelled,
        )
    }

    /// Physical impermeable moving-wall pressure and work from the original EOS.
    ///
    /// The outward normal is -1 at the left end and +1 at the right. The
    /// original shock/isentrope relation sees fluid velocity relative to this
    /// wall. The laboratory flux has exactly zero mass, pressure momentum
    /// traction and energy equal to pressure times wall velocity. Signal speed
    /// is relative to the wall. Original branch/control limits remain in force.
    pub fn moving_reflecting_wall_flux(
        &self,
        inside: Conserved,
        outward_normal: f64,
        wall_velocity: Velocity,
        control: RealGasWallControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FaceFlux, RealGasFlowError> {
        control.validate(self)?;
        let primitive = self.decode(inside, &mut cancelled)?;
        wall_flux(
            self,
            (inside, primitive),
            outward_normal,
            wall_velocity.value(),
            RealGasBoundary::PhysicalReflecting(control),
            &mut cancelled,
        )
    }
}

pub(in crate::dynamics::nonideal::evolve) fn interface(
    gas: &RealGas,
    left: (Conserved, Primitive),
    right: (Conserved, Primitive),
    velocity: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<FaceFlux, RealGasFlowError> {
    if velocity == 0.0 {
        return gas.interface(left.0, left.1, right.0, right.1, cancelled);
    }
    poll(cancelled)?;
    if !velocity.is_finite() {
        return Err(RealGasFlowError::Invalid(
            "nonfinite real-gas face velocity",
        ));
    }
    let (ql, pl) = gas.state_at(
        left.1.density,
        left.1.velocity - velocity,
        left.1.temperature,
        cancelled,
    )?;
    let (qr, pr) = gas.state_at(
        right.1.density,
        right.1.velocity - velocity,
        right.1.temperature,
        cancelled,
    )?;
    let relative = gas.interface(ql, pl, qr, pr, cancelled)?;
    laboratory(relative, velocity, cancelled)
}

pub(in crate::dynamics::nonideal::evolve) fn wall_flux(
    gas: &RealGas,
    inside: (Conserved, Primitive),
    outward_normal: f64,
    velocity: f64,
    boundary: RealGasBoundary,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<FaceFlux, RealGasFlowError> {
    poll(cancelled)?;
    if !matches!(outward_normal, -1.0 | 1.0)
        || !velocity.is_finite()
        || !(inside.1.velocity - velocity).is_finite()
    {
        return Err(RealGasFlowError::Invalid("invalid real-gas wall frame"));
    }
    let relative = Primitive {
        velocity: inside.1.velocity - velocity,
        ..inside.1
    };
    let mut flux = match boundary {
        RealGasBoundary::PhysicalReflecting(control) => {
            gas.physical_wall_flux(relative, outward_normal, control, cancelled)?
        }
        RealGasBoundary::Reflecting => {
            let (q, p) = if velocity == 0.0 {
                inside
            } else {
                gas.state_at(
                    relative.density,
                    relative.velocity,
                    relative.temperature,
                    cancelled,
                )?
            };
            let mirror_q = Conserved {
                momentum: -q.momentum,
                ..q
            };
            let mirror_p = Primitive {
                velocity: -p.velocity,
                ..p
            };
            if outward_normal < 0.0 {
                gas.interface(mirror_q, mirror_p, q, p, cancelled)?
            } else {
                gas.interface(q, p, mirror_q, mirror_p, cancelled)?
            }
        }
        _ => {
            return Err(RealGasFlowError::Invalid(
                "moving real-gas flux requires an impermeable wall",
            ));
        }
    };
    // Impose exactly the declared impermeability and wall work. Pressure
    // traction and relative wave speeds still come from the original relation.
    flux.mass = 0.0;
    flux.energy = velocity * flux.momentum;
    if flux.array().iter().any(|value| !value.is_finite()) {
        return Err(RealGasFlowError::Invalid(
            "unrepresentable real-gas wall work",
        ));
    }
    poll(cancelled)?;
    Ok(flux)
}

fn laboratory(
    relative: FaceFlux,
    velocity: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<FaceFlux, RealGasFlowError> {
    let [mass, momentum, energy] = kinematics::laboratory_flux(relative.array(), velocity);
    let flux = FaceFlux {
        mass,
        momentum,
        energy,
        ..relative
    };
    if flux.array().iter().any(|value| !value.is_finite()) {
        return Err(RealGasFlowError::Invalid(
            "unrepresentable moving real-gas flux",
        ));
    }
    poll(cancelled)?;
    Ok(flux)
}
