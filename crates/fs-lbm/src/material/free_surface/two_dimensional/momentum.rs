//! Source-backed access to the original conservative D2Q9 conversion owner.

use super::{MaterialFreeSurface2, MaterialLbmError, finite};
use crate::core2::Cell;
use crate::freesurface::{ConversionMomentumReport2, FreeSurfaceMomentumError2};
use crate::material::free_surface::wetting::FreeSurfaceWettingError;
use fs_exec::CancelGate;

/// Accepted conversion account per unit out-of-plane depth. Collision,
/// streaming, pressure traction and body-force impulse are outside this
/// report; it is not a complete-step or angular-momentum certificate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialFreeSurfaceMomentumReport2 {
    /// Unmodified account from the original native conversion owner.
    pub native: ConversionMomentumReport2,
    /// Absolute adjacent-cell mass redistribution per depth [kg/m].
    pub redistributed_mass_kg_per_m: f64,
    /// Absolute local mass-representation remainder per depth [kg/m].
    pub absolute_mass_defect_kg_per_m: f64,
    /// Absolute physical linear-momentum remainder per depth [kg/s].
    pub absolute_momentum_defect_kg_per_s: [f64; 2],
    /// Signed numerical kinetic change per depth [J/m]. Mixing and signed debt
    /// can change resolved kinetic storage; this is NEVER deposited as heat.
    pub kinetic_energy_change_j_per_m: f64,
}

/// A refusal preserves the accepted flow and, when attached, wetting history.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialFreeSurfaceMomentumError2 {
    /// Original conversion, force-frame, numerical-domain or cancellation error.
    Native(FreeSurfaceMomentumError2),
    /// The original physical-unit conversion could not represent its result.
    Units(MaterialLbmError),
    /// The attached contact-line geometry or selected material source refused.
    Wetting(FreeSurfaceWettingError),
}
impl core::fmt::Display for MaterialFreeSurfaceMomentumError2 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Native(e) => write!(f, "material free-surface momentum: {e}"),
            Self::Units(e) => write!(f, "material free-surface momentum units: {e}"),
            Self::Wetting(e) => write!(f, "material free-surface momentum: {e}"),
        }
    }
}
impl std::error::Error for MaterialFreeSurfaceMomentumError2 {}
impl From<FreeSurfaceMomentumError2> for MaterialFreeSurfaceMomentumError2 {
    fn from(error: FreeSurfaceMomentumError2) -> Self {
        Self::Native(error)
    }
}
impl From<MaterialLbmError> for MaterialFreeSurfaceMomentumError2 {
    fn from(error: MaterialLbmError) -> Self {
        Self::Units(error)
    }
}
impl From<FreeSurfaceWettingError> for MaterialFreeSurfaceMomentumError2 {
    fn from(error: FreeSurfaceWettingError) -> Self {
        Self::Wetting(error)
    }
}

impl MaterialFreeSurface2<'_> {
    /// Resolved liquid linear momentum per unit out-of-plane depth [kg/s].
    /// Uses tracked interface mass and physical Guo-corrected velocity; gas,
    /// walls and unlocated carry have no represented liquid momentum.
    /// External forces and wall/pressure traction can change this observation;
    /// it is not a whole-step or angular-momentum conservation certificate.
    pub fn resolved_linear_momentum_kg_per_s(&self) -> Result<[f64; 2], MaterialLbmError> {
        let mut momentum = [0.0; 2];
        for (i, &cell) in self.lattice.grid.flags.iter().enumerate() {
            if !matches!(cell, Cell::Fluid | Cell::Interface) {
                continue;
            }
            let state = self.lattice.grid.moments(i);
            let mass = if cell == Cell::Fluid {
                state.rho
            } else {
                self.lattice.mass[i]
            };
            for d in 0..2 {
                momentum[d] += mass * state.u[d];
            }
        }
        let velocity_unit = self.scaling.velocity_si(1.0)?;
        let convert = |p| {
            finite(
                (p * self.mass_per_depth_unit_kg_m) * velocity_unit,
                "2-D liquid momentum is unrepresentable",
            )
        };
        Ok([convert(momentum[0])?, convert(momentum[1])?])
    }

    /// Advance the same material-bound fluid with local momentum-conserving
    /// conversion. The original density, viscosity, tension, force and wetting
    /// bindings are retained. Unsupported capacity or signed debt refuses;
    /// no epsilon mass, velocity clipping, remote redistribution or heat credit
    /// is introduced. The caller still owns spatial/time refinement and its
    /// declared physical maximum-speed bound.
    ///
    /// Native dynamics and finite SI reporting use one candidate. Cancellation,
    /// native refusal or unrepresentable physical units publish no state.
    /// Two-dimensional quantities are per depth: rho_ref*dx^2, not a 3-D cell
    /// mass, supplies the mass unit. Native collision remains a dense sweep;
    /// cancellation polls surround it and occur in streaming/conversion.
    pub fn step_with_momentum(
        &mut self,
        gate: &CancelGate,
    ) -> Result<MaterialFreeSurfaceMomentumReport2, MaterialFreeSurfaceMomentumError2> {
        let (candidate, native) = self.lattice.momentum_candidate(gate)?;
        let report = self.momentum_report_si(native)?;
        if gate.is_requested() {
            return Err(FreeSurfaceMomentumError2::Cancelled.into());
        }
        self.lattice = candidate;
        Ok(report)
    }

    pub(super) fn momentum_report_si(
        &self,
        native: ConversionMomentumReport2,
    ) -> Result<MaterialFreeSurfaceMomentumReport2, MaterialLbmError> {
        MaterialFreeSurfaceMomentumReport2::from_native(
            native,
            self.mass_per_depth_unit_kg_m,
            self.scaling.velocity_si(1.0)?,
        )
    }
}

impl MaterialFreeSurfaceMomentumReport2 {
    fn from_native(
        native: ConversionMomentumReport2,
        mass_unit: f64,
        velocity_unit: f64,
    ) -> Result<Self, MaterialLbmError> {
        // Scale the measured value first so a zero/tiny numerical defect does
        // not overflow through an unnecessary product of two large units.
        let check = |value| finite(value, "D2Q9 conversion account overflows physical units");
        Ok(Self {
            redistributed_mass_kg_per_m: check(native.redistributed_mass * mass_unit)?,
            absolute_mass_defect_kg_per_m: check(native.absolute_mass_defect * mass_unit)?,
            absolute_momentum_defect_kg_per_s: [
                check((native.absolute_momentum_defect[0] * mass_unit) * velocity_unit)?,
                check((native.absolute_momentum_defect[1] * mass_unit) * velocity_unit)?,
            ],
            kinetic_energy_change_j_per_m: check(
                ((native.kinetic_energy_change * mass_unit) * velocity_unit) * velocity_unit,
            )?,
            native,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g1_reports_use_mass_momentum_and_signed_energy_per_depth() {
        let native = ConversionMomentumReport2 {
            redistributed_mass: 0.5,
            absolute_mass_defect: 1e-12,
            absolute_momentum_defect: [2e-12, 3e-12],
            kinetic_energy_change: -0.025,
            routes: 3,
        };
        // rho=800 kg/m^3, dx=.002 m, dt=.0005 s: mass/depth=.0032
        // kg/m, speed=4 m/s, momentum/depth=.0128 kg/s, energy/depth=.0512 J/m.
        let report = MaterialFreeSurfaceMomentumReport2::from_native(native, 0.0032, 4.0).unwrap();
        assert_eq!(report.native, native);
        assert!((report.redistributed_mass_kg_per_m - 0.0016).abs() < 1e-18);
        assert!((report.absolute_mass_defect_kg_per_m - 3.2e-15).abs() < 1e-29);
        assert!((report.absolute_momentum_defect_kg_per_s[0] - 2.56e-14).abs() < 1e-28);
        assert!((report.absolute_momentum_defect_kg_per_s[1] - 3.84e-14).abs() < 1e-28);
        assert!((report.kinetic_energy_change_j_per_m + 0.00128).abs() < 1e-18);
        assert!(MaterialFreeSurfaceMomentumReport2::from_native(native, f64::MAX, 16.0).is_err());
        let zero = MaterialFreeSurfaceMomentumReport2::from_native(
            ConversionMomentumReport2::default(),
            f64::MAX,
            f64::MAX,
        )
        .unwrap();
        assert_eq!(zero.kinetic_energy_change_j_per_m, 0.0);
    }
}
