//! Source-backed access to the original conservative D2Q9 conversion owner.

use super::{
    MaterialFreeSurface2, MaterialGasPressureLink2, MaterialLbmError, finite, pressure_perturbation,
};
use crate::core2::Cell;
use crate::freesurface::{
    ConversionMomentumReport2, ConversionRedistribution2, FreeSurface, FreeSurfaceMomentumError2,
    FreeSurfacePressureError2, FreeSurfacePressureMomentumError2,
};
use crate::material::free_surface::wetting::FreeSurfaceWettingError;
use fs_exec::CancelGate;

/// Accepted conversion account per unit out-of-plane depth. Collision,
/// streaming, pressure traction and body-force impulse are outside this
/// report. Angular closure, when selected, concerns conversion alone; it is
/// not a complete-step or material first-moment conservation certificate.
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
    /// Absolute cell-centre angular defect per depth [kg m/s], summed about
    /// each route component's fixed local origin. Present only in angular mode;
    /// translating residuals also requires the separately reported linear part.
    pub absolute_angular_momentum_defect_kg_m_per_s: Option<f64>,
    /// Signed numerical kinetic change per depth [J/m]. Mixing and signed debt
    /// can change resolved kinetic storage; this is NEVER deposited as heat.
    pub kinetic_energy_change_j_per_m: f64,
}

/// A refusal preserves the accepted flow and, when attached, wetting history.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialFreeSurfaceMomentumError2 {
    /// Original conversion, force-frame, numerical-domain or cancellation error.
    Native(FreeSurfaceMomentumError2),
    /// The prescribed gas load, pressure domain or cancellation gate refused.
    Pressure(FreeSurfacePressureError2),
    /// The material scaling or physical-unit conversion refused its input/result.
    Units(MaterialLbmError),
    /// The attached contact-line geometry or selected material source refused.
    Wetting(FreeSurfaceWettingError),
}
impl core::fmt::Display for MaterialFreeSurfaceMomentumError2 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Native(e) => write!(f, "material free-surface momentum: {e}"),
            Self::Pressure(e) => write!(f, "material free-surface momentum pressure: {e}"),
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
impl From<FreeSurfacePressureMomentumError2> for MaterialFreeSurfaceMomentumError2 {
    fn from(error: FreeSurfacePressureMomentumError2) -> Self {
        match error {
            FreeSurfacePressureMomentumError2::Momentum(e) => Self::Native(e),
            FreeSurfacePressureMomentumError2::Pressure(e) => Self::Pressure(e),
        }
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
    /// Policy shared by every explicit native/material momentum step, including
    /// pressure loading and a subsequently attached source-owned wetting model.
    pub fn conversion_redistribution(&self) -> ConversionRedistribution2 {
        self.lattice.conversion_redistribution()
    }

    /// Select local linear or linear-and-angular conversion conservation.
    /// This changes no accepted physical state. Native angular admission
    /// refuses periodic geometry; unsupported local spin or population states
    /// refuse the whole later candidate. Configure this before attaching a
    /// wetting owner. Legacy native stepping keeps its historical behavior.
    pub fn set_conversion_redistribution(
        &mut self,
        policy: ConversionRedistribution2,
    ) -> Result<(), MaterialFreeSurfaceMomentumError2> {
        self.lattice.set_conversion_redistribution(policy)?;
        Ok(())
    }

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

    /// Advance with absolute gas pressure [Pa] and the original conservative
    /// mass/momentum conversion in the same step. The callback is sampled at
    /// exposed liquid cell centres [m] about the fixed pressure reference.
    /// Pressure is held for this timestep; the combined gas/capillary density
    /// excursion must satisfy the caller's native numerical-domain bound.
    ///
    /// One candidate admits traction, conversion and physical reporting before
    /// publication. The returned account covers conversion only; gas dynamics
    /// and pressure-work heating remain outside this isothermal liquid model.
    pub fn step_with_gas_pressure_and_momentum(
        &mut self,
        gate: &CancelGate,
        maximum_density_excursion: f64,
        pressure_pa: impl Fn([f64; 2]) -> f64,
    ) -> Result<MaterialFreeSurfaceMomentumReport2, MaterialFreeSurfaceMomentumError2> {
        let (candidate, report) =
            self.gas_pressure_momentum_candidate(gate, maximum_density_excursion, pressure_pa)?;
        if gate.is_requested() {
            return Err(FreeSurfaceMomentumError2::Cancelled.into());
        }
        self.lattice = candidate;
        Ok(report)
    }

    /// Apply distinct absolute pressures [Pa] to actual liquid/gas links while
    /// retaining conservative local conversion. Cell centres and link midpoint
    /// are in metres. Uniform link loads reproduce the cell-pressure entry.
    /// The same fixed reference, density bound and atomic publication apply.
    pub fn step_with_gas_link_pressure_and_momentum(
        &mut self,
        gate: &CancelGate,
        maximum_density_excursion: f64,
        pressure_pa: impl Fn(MaterialGasPressureLink2) -> f64,
    ) -> Result<MaterialFreeSurfaceMomentumReport2, MaterialFreeSurfaceMomentumError2> {
        let (candidate, report) = self.gas_link_pressure_momentum_candidate(
            gate,
            maximum_density_excursion,
            pressure_pa,
        )?;
        if gate.is_requested() {
            return Err(FreeSurfaceMomentumError2::Cancelled.into());
        }
        self.lattice = candidate;
        Ok(report)
    }

    pub(super) fn gas_pressure_momentum_candidate(
        &self,
        gate: &CancelGate,
        maximum_density_excursion: f64,
        pressure_pa: impl Fn([f64; 2]) -> f64,
    ) -> Result<(FreeSurface, MaterialFreeSurfaceMomentumReport2), MaterialFreeSurfaceMomentumError2>
    {
        let reference = self
            .scaling
            .reference_pressure_pa
            .ok_or(MaterialLbmError::Invalid(
                "absolute gas loading requires an explicit fixed pressure reference",
            ))?;
        let dx = self.scaling.cell_size_m;
        let unit = self.scaling.pressure_unit_pa;
        let (candidate, native) = self.lattice.gas_pressure_momentum_candidate(
            gate,
            maximum_density_excursion,
            |coordinates| {
                let position = coordinates.map(|c| (c as f64 + 0.5) * dx);
                pressure_perturbation(
                    reference,
                    unit,
                    pressure_pa(position),
                    coordinates,
                    &position,
                )
            },
        )?;
        Ok((candidate, self.momentum_report_si(native)?))
    }

    pub(super) fn gas_link_pressure_momentum_candidate(
        &self,
        gate: &CancelGate,
        maximum_density_excursion: f64,
        pressure_pa: impl Fn(MaterialGasPressureLink2) -> f64,
    ) -> Result<(FreeSurface, MaterialFreeSurfaceMomentumReport2), MaterialFreeSurfaceMomentumError2>
    {
        let reference = self
            .scaling
            .reference_pressure_pa
            .ok_or(MaterialLbmError::Invalid(
                "absolute gas loading requires an explicit fixed pressure reference",
            ))?;
        let dx = self.scaling.cell_size_m;
        let unit = self.scaling.pressure_unit_pa;
        let (candidate, native) = self.lattice.gas_link_pressure_momentum_candidate(
            gate,
            maximum_density_excursion,
            |link| {
                let liquid_center_m = link.liquid_cell.map(|c| (c as f64 + 0.5) * dx);
                let gas_center_m = link.gas_cell.map(|c| (c as f64 + 0.5) * dx);
                let midpoint_m =
                    core::array::from_fn(|d| f64::midpoint(liquid_center_m[d], gas_center_m[d]));
                let physical = MaterialGasPressureLink2 {
                    link,
                    liquid_center_m,
                    gas_center_m,
                    midpoint_m,
                };
                pressure_perturbation(
                    reference,
                    unit,
                    pressure_pa(physical),
                    link.liquid_cell,
                    &[
                        liquid_center_m[0],
                        liquid_center_m[1],
                        gas_center_m[0],
                        gas_center_m[1],
                        midpoint_m[0],
                        midpoint_m[1],
                    ],
                )
            },
        )?;
        Ok((candidate, self.momentum_report_si(native)?))
    }

    /// The same source-bound step with cell-centre angular conservation on
    /// each actual conversion component. Periodic coordinates, unresolved spin
    /// or inadmissible corrected populations refuse before flow publication.
    /// The selected inertial-frame and first-mass-moment limitations of
    /// `FreeSurface::step_with_angular_momentum` apply. Numerical remap energy
    /// remains separate from physical heat and wetting dissipation.
    pub fn step_with_angular_momentum(
        &mut self,
        gate: &CancelGate,
    ) -> Result<MaterialFreeSurfaceMomentumReport2, MaterialFreeSurfaceMomentumError2> {
        let (candidate, native) = self.lattice.angular_momentum_candidate(gate)?;
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
            self.scaling.cell_size_m,
        )
    }
}

impl MaterialFreeSurfaceMomentumReport2 {
    fn from_native(
        native: ConversionMomentumReport2,
        mass_unit: f64,
        velocity_unit: f64,
        length_unit: f64,
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
            absolute_angular_momentum_defect_kg_m_per_s: native
                .absolute_angular_momentum_defect
                .map(|defect| check(((defect * mass_unit) * velocity_unit) * length_unit))
                .transpose()?,
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
            absolute_angular_momentum_defect: Some(4e-12),
            kinetic_energy_change: -0.025,
            routes: 3,
        };
        // rho=800 kg/m^3, dx=.002 m, dt=.0005 s: mass/depth=.0032
        // kg/m, speed=4 m/s, momentum/depth=.0128 kg/s, energy/depth=.0512 J/m,
        // and angular momentum/depth=.0000256 kg m/s.
        let report =
            MaterialFreeSurfaceMomentumReport2::from_native(native, 0.0032, 4.0, 0.002).unwrap();
        assert_eq!(report.native, native);
        assert!((report.redistributed_mass_kg_per_m - 0.0016).abs() < 1e-18);
        assert!((report.absolute_mass_defect_kg_per_m - 3.2e-15).abs() < 1e-29);
        assert!((report.absolute_momentum_defect_kg_per_s[0] - 2.56e-14).abs() < 1e-28);
        assert!((report.absolute_momentum_defect_kg_per_s[1] - 3.84e-14).abs() < 1e-28);
        assert!(
            (report.absolute_angular_momentum_defect_kg_m_per_s.unwrap() - 1.024e-16).abs() < 1e-30
        );
        assert!((report.kinetic_energy_change_j_per_m + 0.00128).abs() < 1e-18);
        assert!(
            MaterialFreeSurfaceMomentumReport2::from_native(native, f64::MAX, 16.0, 1.0).is_err()
        );
        let zero = MaterialFreeSurfaceMomentumReport2::from_native(
            ConversionMomentumReport2::default(),
            f64::MAX,
            f64::MAX,
            f64::MAX,
        )
        .unwrap();
        assert_eq!(zero.kinetic_energy_change_j_per_m, 0.0);
        assert_eq!(zero.absolute_angular_momentum_defect_kg_m_per_s, None);
        let angular_zero = ConversionMomentumReport2 {
            absolute_angular_momentum_defect: Some(0.0),
            ..Default::default()
        };
        assert_eq!(
            MaterialFreeSurfaceMomentumReport2::from_native(
                angular_zero,
                f64::MAX,
                f64::MAX,
                f64::MAX,
            )
            .unwrap()
            .absolute_angular_momentum_defect_kg_m_per_s,
            Some(0.0)
        );
        let angular_only = ConversionMomentumReport2 {
            absolute_angular_momentum_defect: Some(1.0),
            ..Default::default()
        };
        assert!(
            MaterialFreeSurfaceMomentumReport2::from_native(angular_only, 1.0, 1.0, f64::INFINITY)
                .is_err()
        );
    }
}
