//! Species streams and bulk impulse on the original finite moving gas.
//!
//! Stream energy may use any explicitly declared INERTIAL frame. The original
//! gas owner admits its final mass, then the same callback removes the actual
//! kinetic change in that frame before its one caloric reconstruction. This
//! avoids both a second thermal solve and assigning kinetic energy to an
//! arbitrary species. No phase rate, mixing-heat correction, pressure load,
//! moving-volume work or force law is inferred here.

use super::{EvaporatingLiquidGasError, MovingGasParcel};
use fs_material::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::{
        GasEnergyError,
        parcel::transfer::species::boundary::{
            GasSpeciesBoundarySolve, GasSpeciesBoundaryStep, SpeciesBoundaryRelativeTransfer,
        },
    },
};
use fs_math::dd::Dd;

#[derive(Clone, Copy, Debug)]
pub struct MovingGasBoundarySolve {
    pub gas: GasSpeciesBoundarySolve,
    /// Per-component absolute allowance on ACTUAL minus requested impulse.
    pub momentum_absolute_n_s: f64,
    /// Relative allowance scales only the requested component, never common
    /// drift or a large pressure impulse on a different axis.
    pub momentum_relative_tolerance: f64,
}

#[derive(Debug)]
pub enum MovingGasBoundaryError {
    Invalid(&'static str),
    Gas(GasEnergyError),
    Motion(EvaporatingLiquidGasError),
    Momentum {
        component: usize,
        residual_n_s: f64,
        allowance_n_s: f64,
    },
    Cancelled,
}
impl core::fmt::Display for MovingGasBoundaryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "moving gas boundary refused: {self:?}")
    }
}
impl std::error::Error for MovingGasBoundaryError {}
impl From<GasEnergyError> for MovingGasBoundaryError {
    fn from(error: GasEnergyError) -> Self {
        Self::Gas(error)
    }
}

/// Accounts from actual represented gas mass and momentum endpoints.
#[derive(Clone, Debug, PartialEq)]
pub struct MovingGasBoundaryStep {
    /// Original species/internal-energy update, including the independently
    /// reported subtraction of bulk kinetic change.
    pub caloric: GasSpeciesBoundaryStep,
    /// Velocity of the stream-energy frame relative to the momentum frame.
    pub energy_frame_velocity_m_s: [f64; 3],
    pub requested_impulse_n_s: [f64; 3],
    pub actual_impulse_n_s: [f64; 3],
    pub momentum_residual_n_s: [f64; 3],
    /// Actual kinetic change in the declared stream-energy frame.
    pub kinetic_energy_change_j: f64,
    /// Same actual endpoints, in the original momentum frame.
    pub laboratory_kinetic_energy_change_j: f64,
    /// Absolute source-datum total energy, for reporting only. Small closure
    /// is audited in the original formation-independent caloric coordinate.
    pub total_energy_into_gas_j: f64,
    pub total_energy_residual_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MovingGasBoundaryCandidate<'catalog> {
    gas: MovingGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    step: MovingGasBoundaryStep,
}
impl<'catalog> MovingGasBoundaryCandidate<'catalog> {
    pub const fn gas(&self) -> &MovingGasParcel {
        &self.gas
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.catalog
    }
    pub const fn step(&self) -> &MovingGasBoundaryStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        MovingGasParcel,
        &'catalog ThermochemicalGasModel,
        MovingGasBoundaryStep,
    ) {
        (self.gas, self.catalog, self.step)
    }
}

fn finite(value: f64) -> Result<f64, GasEnergyError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(GasEnergyError::Invalid(
            "moving gas kinetic increment is unrepresentable",
        ))
    }
}

// Exact identity for changing mass:
// Delta K = (v0+v1)/2 . Delta P - (v0.v1)/2 Delta m.
// All velocities and momenta below use the declared inertial energy frame.
// Subtract frame momentum BEFORE division and retain double-double increments;
// neither large absolute kinetic total is ever formed and subtracted.
fn kinetic_change(
    m0: f64,
    m1: f64,
    p0: [f64; 3],
    p1: [f64; 3],
    frame: [f64; 3],
) -> Result<f64, GasEnergyError> {
    if !m1.is_finite() || m1 <= 0.0 {
        return Err(GasEnergyError::Invalid(
            "positive admitted moving gas mass required",
        ));
    }
    let dm = Dd::from_f64(m1) - Dd::from_f64(m0);
    let mut change = Dd::ZERO;
    for d in 0..3 {
        let reference = Dd::from_f64(frame[d]);
        let before = Dd::from_f64(p0[d]) - Dd::from_f64(m0) * reference;
        let after = Dd::from_f64(p1[d]) - Dd::from_f64(m1) * reference;
        let v0 = before / Dd::from_f64(m0);
        let v1 = after / Dd::from_f64(m1);
        let dp = (Dd::from_f64(p1[d]) - Dd::from_f64(p0[d])) - reference * dm;
        change = change + (v0 * Dd::from_f64(0.5) + v1 * Dd::from_f64(0.5)) * dp
            - (dm * Dd::from_f64(0.5) * v0) * v1;
    }
    finite(change.to_f64())
}

impl MovingGasParcel {
    /// Advance arbitrary admitted species streams and their complete impulse.
    ///
    /// `impulse_n_s` uses the same inertial axes/frame as the retained momentum.
    /// Each energy row means total molecular/stream energy in the declared
    /// energy frame, minus dm_i times ORIGINAL h_i(T_gas). For the laboratory
    /// frame supply zero frame velocity. A common translating frame can retain
    /// thermal energy beside a large shared drift: construct energy in that
    /// frame directly, never by subtracting already rounded laboratory totals.
    ///
    /// The original mass/source owner invokes kinetic accounting once with its
    /// exact final mass. It then performs ONE original gas caloric inversion.
    /// Additional receiver-mixing heat would count this energy twice. Callers
    /// must admit gross inventories and close the opposite phase/body exchanges
    /// atomically. Original fixed gas volume and auxiliary source axes remain.
    #[allow(clippy::too_many_arguments)]
    pub fn exchange_species_momentum_relative<'catalog>(
        &self,
        catalog: &'catalog ThermochemicalGasModel,
        transfers: &[SpeciesBoundaryRelativeTransfer],
        impulse_n_s: [f64; 3],
        energy_frame_velocity_m_s: [f64; 3],
        controls: MovingGasBoundarySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MovingGasBoundaryCandidate<'catalog>, MovingGasBoundaryError> {
        if cancelled() {
            return Err(MovingGasBoundaryError::Cancelled);
        }
        if !impulse_n_s
            .into_iter()
            .chain(energy_frame_velocity_m_s)
            .all(f64::is_finite)
            || !controls.momentum_absolute_n_s.is_finite()
            || controls.momentum_absolute_n_s <= 0.0
            || !controls.momentum_relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&controls.momentum_relative_tolerance)
        {
            return Err(MovingGasBoundaryError::Invalid(
                "finite impulse/frame and positive momentum bounds required",
            ));
        }
        let p0 = self.momentum_kg_m_s();
        let p1 = core::array::from_fn(|d| p0[d] + impulse_n_s[d]);
        let mut actual_impulse = [0.0; 3];
        let mut residual = [0.0; 3];
        for d in 0..3 {
            let change = Dd::from_f64(p1[d]) - Dd::from_f64(p0[d]);
            actual_impulse[d] = finite(change.to_f64())?;
            residual[d] = finite((change - Dd::from_f64(impulse_n_s[d])).to_f64())?;
            let allowance = controls.momentum_absolute_n_s
                + controls.momentum_relative_tolerance * impulse_n_s[d].abs();
            if impulse_n_s[d] != 0.0 && actual_impulse[d] == 0.0 {
                return Err(MovingGasBoundaryError::Invalid(
                    "requested molecular impulse vanished at the original gas momentum",
                ));
            }
            if !allowance.is_finite() || residual[d].abs() > allowance {
                return Err(MovingGasBoundaryError::Momentum {
                    component: d,
                    residual_n_s: residual[d],
                    allowance_n_s: allowance,
                });
            }
        }
        let m0 = self.parcel().mass_kg();
        let mut evaluated = None;
        let candidate = self
            .parcel()
            .exchange_species_boundary_relative_with_energy(
                catalog,
                transfers,
                |mass| {
                    let kinetic = kinetic_change(m0, mass, p0, p1, energy_frame_velocity_m_s)?;
                    evaluated = Some((mass, kinetic));
                    Ok(-kinetic)
                },
                controls.gas,
                &mut cancelled,
            )?;
        let (parcel, catalog, caloric) = candidate.into_parts();
        let (mass, kinetic) = evaluated.ok_or(MovingGasBoundaryError::Invalid(
            "original gas owner omitted the energy callback",
        ))?;
        if mass != parcel.mass_kg() {
            return Err(MovingGasBoundaryError::Invalid(
                "original gas changed admitted mass after kinetic accounting",
            ));
        }
        let laboratory_kinetic_energy_change_j = kinetic_change(m0, mass, p0, p1, [0.0; 3])?;
        let total_energy_into_gas_j =
            finite((Dd::from_f64(caloric.energy_into_gas_j) + Dd::from_f64(kinetic)).to_f64())?;
        let total_energy_residual_j = caloric.energy_residual_j;
        let gas =
            MovingGasParcel::from_momentum(parcel, p1).map_err(MovingGasBoundaryError::Motion)?;
        if cancelled() {
            return Err(MovingGasBoundaryError::Cancelled);
        }
        Ok(MovingGasBoundaryCandidate {
            gas,
            catalog,
            step: MovingGasBoundaryStep {
                caloric,
                energy_frame_velocity_m_s,
                requested_impulse_n_s: impulse_n_s,
                actual_impulse_n_s: actual_impulse,
                momentum_residual_n_s: residual,
                kinetic_energy_change_j: kinetic,
                laboratory_kinetic_energy_change_j,
                total_energy_into_gas_j,
                total_energy_residual_j,
            },
        })
    }
}

#[cfg(test)]
mod tests;
