//! Borrowed sources for one original thermal wall transaction.

use super::*;
use crate::material::free_surface::{MaterialFreeSurfaceScaling, MaterialGasPressureLink3};

/// Prescribed positive absolute gas pressure [Pa] on actual exposed links.
///
/// The field is held over one physical tick and is never retained by the
/// liquid. It must return the same value for the same coordinate or link
/// throughout any finite-body pilot and final replay. Interior mutation of a
/// callback must not change that physical input. Cell fields are sampled once
/// per exposed liquid cell; link fields may distinguish its gas neighbors.
///
/// The construction-time scaling must carry an explicit pressure reference.
/// Gas plus capillary loading obeys the original numerical density-excursion
/// bound. This supplies mechanical traction, without evolving a gas reservoir,
/// compressive heating or pressure-work energy closure; applicability of the
/// fixed source-pressure caloric law remains with the caller.
#[derive(Clone, Copy)]
pub enum FreeSurfaceWallGasPressure<'a> {
    /// Field at exposed liquid cell centers in metres.
    Cell {
        /// Maximum `|rho_boundary - 1|`, finite and strictly between zero and one.
        maximum_density_excursion: f64,
        pressure_pa: &'a dyn Fn([f64; 3]) -> f64,
    },
    /// Field at each exact native gas link, with physical cell and midpoint coordinates.
    Link {
        /// Maximum `|rho_boundary - 1|`, finite and strictly between zero and one.
        maximum_density_excursion: f64,
        pressure_pa: &'a dyn Fn(MaterialGasPressureLink3) -> f64,
    },
}

/// All explicitly selected sources for one moving-wall thermal candidate.
///
/// Authored and optional molecular heat are applied once at the accepted
/// initial state. The same collision then receives the refreshed local
/// material/force fields, gas loading and every wall reflection. Its sole mass
/// trace transports extensive enthalpy. Wall and pressure work are not added
/// to caloric heat. The default selects no heat or prescribed gas load.
///
/// This borrowed, copyable input can be reused by a finite-body pilot and
/// replay. Each trial begins from the same accepted liquid; only the final
/// candidate may publish heat, momentum and the thermal clock together.
#[derive(Clone, Copy, Default)]
pub struct FreeSurfaceWallSources<'a> {
    /// Finite authored cell heat [J], with unique positive-mass coordinates.
    pub heat: &'a [FreeSurfaceCellHeat],
    /// Select the original held initial-state molecular stress-power source.
    pub molecular_heating: bool,
    /// Optional frozen mechanical gas loading on the actual native interfaces.
    pub gas_pressure: Option<FreeSurfaceWallGasPressure<'a>>,
}

impl FreeSurfaceWallSources<'_> {
    // Keep the original physical-coordinate conversion and native validation.
    // A stack-local conversion closure reaches the same original collision;
    // no sampled map, allocation, retained pressure or extra flow clone exists.
    pub(super) fn with_lattice_gas_pressure<R>(
        self,
        scaling: MaterialFreeSurfaceScaling<'_>,
        advance: impl FnOnce(Option<(f64, GasPressureField<'_>)>) -> Result<R, Error>,
    ) -> Result<R, Error> {
        let Some(pressure) = self.gas_pressure else {
            return advance(None);
        };
        let reference = scaling.reference_pressure_pa.ok_or(Error::Invalid(
            "absolute gas loading requires an explicit fixed pressure reference",
        ))?;
        match pressure {
            FreeSurfaceWallGasPressure::Cell {
                maximum_density_excursion,
                pressure_pa,
            } => {
                let field = crate::material::free_surface::pressure::lattice_pressure_field(
                    reference,
                    scaling.cell_size_m,
                    scaling.pressure_unit_pa,
                    pressure_pa,
                );
                advance(Some((
                    maximum_density_excursion,
                    GasPressureField::Cell(&field),
                )))
            }
            FreeSurfaceWallGasPressure::Link {
                maximum_density_excursion,
                pressure_pa,
            } => {
                let field = crate::material::free_surface::pressure::lattice_link_pressure_field(
                    reference,
                    scaling.cell_size_m,
                    scaling.pressure_unit_pa,
                    pressure_pa,
                );
                advance(Some((
                    maximum_density_excursion,
                    GasPressureField::Link(&field),
                )))
            }
        }
    }
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    pub(super) fn prepare_wall_sources(
        &self,
        sources: FreeSurfaceWallSources<'_>,
        gate: &CancelGate,
    ) -> Result<viscous::PreparedViscousHeat, Error> {
        if sources.molecular_heating {
            self.prepare_viscous_heat(sources.heat, gate)
        } else {
            let (heated, external, absolute_heat) = self.prepare_heat(sources.heat, gate)?;
            Ok(viscous::PreparedViscousHeat {
                heated,
                external,
                absolute_heat,
                viscous: 0.0,
                applied: Vec::new(),
            })
        }
    }
}
