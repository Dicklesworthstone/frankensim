//! Prescribed normal gas traction on the original dense liquid VOF owner.
//!
//! The load is a perturbation about the pressure represented by liquid density
//! one. The liquid reconstruction satisfies p_liquid = p_gas + sigma*kappa;
//! numerical lattice compressibility is not a gas or liquid equation of state.

use super::*;
use fs_exec::CancelGate;

#[cfg(test)]
mod tests;
#[cfg(test)]
mod momentum_tests;

/// One actual missing incoming D2Q9 population at a liquid/gas link.
/// Coordinates are resolved native cell indices, including periodic wrapping.
/// They identify cell centres, not a reconstructed interface or surface normal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GasPressureLink2 {
    /// Liquid collision cell receiving this population.
    pub liquid_cell: [usize; 2],
    /// Resolved neighboring gas cell supplying this incoming characteristic.
    pub gas_cell: [usize; 2],
    /// Native incoming population index in `1..9`.
    pub incoming_direction: usize,
}

/// An inadmissible pressure/geometry input or cancelled unpublished step.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum FreeSurfacePressureError2 {
    /// The declared pressure, numerical domain or existing wall geometry refused.
    Invalid {
        /// Rejected liquid cell, when the error concerns one local load.
        coordinates: Option<[usize; 2]>,
        /// Concrete reason the candidate could not be admitted.
        reason: &'static str,
    },
    /// Cancellation observed before publishing the complete candidate.
    Cancelled,
}
impl core::fmt::Display for FreeSurfacePressureError2 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "2-D free-surface pressure: {self:?}")
    }
}
impl std::error::Error for FreeSurfacePressureError2 {}

/// A prescribed-traction step failed before publishing its conservative
/// momentum candidate. The original owner and its detailed error are retained.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum FreeSurfacePressureMomentumError2 {
    /// Pressure callback, combined gas/capillary bound, geometry or cancellation.
    Pressure(FreeSurfacePressureError2),
    /// Original liquid admission, signed conversion, capacity or cancellation.
    Momentum(FreeSurfaceMomentumError2),
}
impl core::fmt::Display for FreeSurfacePressureMomentumError2 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Pressure(error) => write!(f, "2-D pressure/momentum step: {error}"),
            Self::Momentum(error) => write!(f, "2-D pressure/momentum step: {error}"),
        }
    }
}
impl std::error::Error for FreeSurfacePressureMomentumError2 {}
impl From<FreeSurfacePressureError2> for FreeSurfacePressureMomentumError2 {
    fn from(error: FreeSurfacePressureError2) -> Self {
        Self::Pressure(error)
    }
}
impl From<FreeSurfaceMomentumError2> for FreeSurfacePressureMomentumError2 {
    fn from(error: FreeSurfaceMomentumError2) -> Self {
        Self::Momentum(error)
    }
}

pub(super) fn check_cancel(gate: Option<&CancelGate>) -> Result<(), FreeSurfacePressureError2> {
    if gate.is_some_and(CancelGate::is_requested) {
        Err(FreeSurfacePressureError2::Cancelled)
    } else {
        Ok(())
    }
}

#[derive(Clone, Copy)]
enum Field<'a> {
    Cell(&'a dyn Fn([usize; 2]) -> Result<f64, FreeSurfacePressureError2>),
    Link(&'a dyn Fn(GasPressureLink2) -> Result<f64, FreeSurfacePressureError2>),
}

pub(super) struct GasPressure2<'a> {
    field: Field<'a>,
    maximum_density_excursion: f64,
}
impl<'a> GasPressure2<'a> {
    fn new(
        field: Field<'a>,
        maximum_density_excursion: f64,
    ) -> Result<Self, FreeSurfacePressureError2> {
        if !maximum_density_excursion.is_finite()
            || maximum_density_excursion <= 0.0
            || maximum_density_excursion >= 1.0
        {
            return Err(FreeSurfacePressureError2::Invalid {
                coordinates: None,
                reason: "density-excursion bound must be finite and in (0,1)",
            });
        }
        Ok(Self {
            field,
            maximum_density_excursion,
        })
    }

    fn density(
        &self,
        capillary: f64,
        link: GasPressureLink2,
    ) -> Result<f64, FreeSurfacePressureError2> {
        let pressure = match self.field {
            Field::Cell(field) => field(link.liquid_cell)?,
            Field::Link(field) => field(link)?,
        };
        let density = capillary + pressure / CS2;
        if !pressure.is_finite()
            || !density.is_finite()
            || density <= 0.0
            || (density - 1.0).abs() > self.maximum_density_excursion
        {
            return Err(FreeSurfacePressureError2::Invalid {
                coordinates: Some(link.liquid_cell),
                reason: "gas plus capillary pressure exceeds the declared density-excursion bound",
            });
        }
        Ok(density)
    }
}

#[derive(Default)]
pub(super) struct GasEquilibrium2 {
    capillary: Option<f64>,
    density: Option<f64>,
    equilibrium: Option<[f64; Q]>,
}
impl GasEquilibrium2 {
    pub(super) fn at_link(
        &mut self,
        reference: impl FnOnce() -> f64,
        velocity: [f64; 2],
        force: [f64; 2],
        boundary: Option<&GasPressure2<'_>>,
        link: GasPressureLink2,
    ) -> Result<&[f64; Q], FreeSurfacePressureError2> {
        if self.density.is_none() || boundary.is_some_and(|p| matches!(p.field, Field::Link(_))) {
            let capillary = *self.capillary.get_or_insert_with(reference);
            let density = match boundary {
                Some(p) => p.density(capillary, link)?,
                None => capillary,
            };
            if self
                .density
                .is_none_or(|old| old.to_bits() != density.to_bits())
            {
                let mut eq = equilibrium(density, velocity[0], velocity[1]);
                // Stored populations are f_phys-S(u,F)/2. Using an unshifted
                // equilibrium imposes a spurious u*F normal stress whenever
                // a forced liquid moves along the loaded interface. This is
                // the same physical source frame as the original 3-D owner.
                if force != [0.0; 2] {
                    crate::thermal::forcing::reframe2(&mut eq, velocity, [0.0; 2], force);
                }
                self.equilibrium = Some(eq);
                self.density = Some(density);
            }
        }
        Ok(self.equilibrium.as_ref().expect("gas reconstruction"))
    }
}

impl FreeSurface {
    /// Advance with pressure perturbations in lattice units about density one.
    /// The callback is sampled once per exposed liquid cell, in row-major
    /// order, and describes one frozen time level. Zero reproduces `step`.
    ///
    /// The combined gas and capillary pressure must satisfy
    /// `|rho_boundary-1| <= maximum_density_excursion` in (0,1). This is an
    /// explicit numerical domain, not a physical compressibility-error bound.
    /// The original collision, mass exchange and local conversion remain the
    /// owners; gas dynamics and pressure-work storage are not introduced.
    ///
    /// All pressure, geometry and cancellation refusals preserve populations,
    /// mass, flags and conversion counters. Cancellation is checked before and
    /// after the candidate, and at each streaming row. No pressure map is kept
    /// between steps, so advancing interfaces cannot inherit stale loads.
    pub fn step_with_gas_pressure(
        &mut self,
        gate: &CancelGate,
        maximum_density_excursion: f64,
        pressure: impl Fn([usize; 2]) -> Result<f64, FreeSurfacePressureError2>,
    ) -> Result<(), FreeSurfacePressureError2> {
        self.step_with_pressure(gate, maximum_density_excursion, Field::Cell(&pressure))
    }

    /// Advance with an independently prescribed pressure on every actual gas
    /// link, in row-major cell order and then increasing incoming direction.
    /// A uniform link load exactly reproduces the cell-pressure arithmetic.
    /// Wall and absent sources remain walls and never call the pressure field.
    /// The same combined-pressure bound and atomicity as `step_with_gas_pressure`
    /// apply. Lattice directions are not geometric surface normals.
    pub fn step_with_gas_link_pressure(
        &mut self,
        gate: &CancelGate,
        maximum_density_excursion: f64,
        pressure: impl Fn(GasPressureLink2) -> Result<f64, FreeSurfacePressureError2>,
    ) -> Result<(), FreeSurfacePressureError2> {
        self.step_with_pressure(gate, maximum_density_excursion, Field::Link(&pressure))
    }

    /// Advance prescribed gas traction and capacity-bounded momentum
    /// conversion on ONE native liquid candidate. Pressure perturbations,
    /// physical Guo stress and the combined gas/capillary density bound have
    /// the same meaning as in [`Self::step_with_gas_pressure`]. The callback
    /// is sampled once per exposed liquid cell, only during native streaming.
    /// Cells exposed by subsequent conversion are first sampled next step.
    ///
    /// Admitted zero pressure reproduces [`Self::step_with_momentum`] bit for bit.
    /// Signed conversion mass carries its donor's physical velocity through
    /// the original local, capacity-bounded redistribution. Unsupported carry,
    /// capacity, populations or low-Mach states refuse instead of falling back
    /// to legacy mass-only conversion. Every refusal leaves accepted state intact.
    ///
    /// The returned account covers CONVERSION alone, excluding pressure and
    /// body-force impulses. The selected `ConservativeLocalAngular` policy
    /// additionally closes conversion angular momentum per route-connected
    /// region on nonperiodic grids; it does not remove physical pressure torque
    /// or certify angular conservation of the full step. Gas dynamics and
    /// pressure-work/thermal storage are not introduced by this composition.
    pub fn step_with_gas_pressure_and_momentum(
        &mut self,
        gate: &CancelGate,
        maximum_density_excursion: f64,
        pressure: impl Fn([usize; 2]) -> Result<f64, FreeSurfacePressureError2>,
    ) -> Result<ConversionMomentumReport2, FreeSurfacePressureMomentumError2> {
        let (candidate, report) =
            self.gas_pressure_momentum_candidate(gate, maximum_density_excursion, pressure)?;
        check_cancel(Some(gate))?;
        *self = candidate;
        Ok(report)
    }

    /// Compose distinct pressure at each actual gas link with the same native
    /// momentum conversion. Callback order remains row-major wet cell then
    /// increasing incoming direction. Walls never sample the field, and no
    /// pressure map survives the step. A uniform link load exactly reproduces
    /// [`Self::step_with_gas_pressure_and_momentum`]. Its numerical bounds,
    /// conversion-only accounting and all-or-nothing publication also apply.
    pub fn step_with_gas_link_pressure_and_momentum(
        &mut self,
        gate: &CancelGate,
        maximum_density_excursion: f64,
        pressure: impl Fn(GasPressureLink2) -> Result<f64, FreeSurfacePressureError2>,
    ) -> Result<ConversionMomentumReport2, FreeSurfacePressureMomentumError2> {
        let (candidate, report) =
            self.gas_link_pressure_momentum_candidate(gate, maximum_density_excursion, pressure)?;
        check_cancel(Some(gate))?;
        *self = candidate;
        Ok(report)
    }

    // Material units and contact-line history admit this same private proposal
    // before publishing. The callbacks and pressure-owner internals do not
    // escape or become persistent state on the accepted free surface.
    pub(crate) fn gas_pressure_momentum_candidate(
        &self,
        gate: &CancelGate,
        maximum_density_excursion: f64,
        pressure: impl Fn([usize; 2]) -> Result<f64, FreeSurfacePressureError2>,
    ) -> Result<(Self, ConversionMomentumReport2), FreeSurfacePressureMomentumError2> {
        self.pressure_momentum_candidate(gate, maximum_density_excursion, Field::Cell(&pressure))
    }

    pub(crate) fn gas_link_pressure_momentum_candidate(
        &self,
        gate: &CancelGate,
        maximum_density_excursion: f64,
        pressure: impl Fn(GasPressureLink2) -> Result<f64, FreeSurfacePressureError2>,
    ) -> Result<(Self, ConversionMomentumReport2), FreeSurfacePressureMomentumError2> {
        self.pressure_momentum_candidate(gate, maximum_density_excursion, Field::Link(&pressure))
    }

    fn pressure_momentum_candidate(
        &self,
        gate: &CancelGate,
        maximum_density_excursion: f64,
        field: Field<'_>,
    ) -> Result<(Self, ConversionMomentumReport2), FreeSurfacePressureMomentumError2> {
        check_cancel(Some(gate))?;
        let boundary = GasPressure2::new(field, maximum_density_excursion)?;
        self.momentum_candidate_with_stream(gate, |candidate| {
            candidate
                .stream_step_with_pressure(Some(gate), Some(&boundary))
                .map_err(FreeSurfacePressureMomentumError2::Pressure)
        })
    }

    fn step_with_pressure(
        &mut self,
        gate: &CancelGate,
        maximum_density_excursion: f64,
        field: Field<'_>,
    ) -> Result<(), FreeSurfacePressureError2> {
        check_cancel(Some(gate))?;
        let boundary = GasPressure2::new(field, maximum_density_excursion)?;
        let mut candidate = self.clone();
        candidate.stream_step_with_pressure(Some(gate), Some(&boundary))?;
        candidate.apply_conversions();
        check_cancel(Some(gate))?;
        *self = candidate;
        Ok(())
    }
}
