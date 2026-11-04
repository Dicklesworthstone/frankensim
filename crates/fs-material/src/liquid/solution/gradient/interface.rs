//! The same liquid-liquid interface energy and original coexisting sources.
//!
//! The finite inventory determines whether two phases actually coexist. The
//! planar interfacial calculation supplies an area energy for that same pair,
//! temperature and square-gradient law. It does not supply a wall contact angle
//! or claim that a thick interface is equivalent to a sharp curved meniscus.

use super::super::{BinaryLiquidError, ResolvedBinaryLiquidEquilibrium, poll};
use super::BinaryLiquidSquareGradientLaw;
use fs_thermochem::solution::gradient::interface::{
    BinaryPlanarInterfaceOptionsV1, BinaryPlanarInterfaceV1,
};

/// An actual two-phase inventory, original material receipts and same-pair
/// planar interface calculation. Private construction prevents unrelated
/// phases and surface energies from being silently joined by a consumer.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedBinaryLiquidPlanarInterface {
    equilibrium: ResolvedBinaryLiquidEquilibrium,
    interface: BinaryPlanarInterfaceV1,
}
impl ResolvedBinaryLiquidPlanarInterface {
    pub const fn equilibrium(&self) -> &ResolvedBinaryLiquidEquilibrium {
        &self.equilibrium
    }
    pub const fn interface(&self) -> &BinaryPlanarInterfaceV1 {
        &self.interface
    }
}

impl BinaryLiquidSquareGradientLaw<'_, '_> {
    /// Resolve the actual coexisting phase sources and their planar interface.
    /// A one-phase inventory returns None, even if another overall composition
    /// of this pair could split. Every actual phase must satisfy the same
    /// source-density/fixed-site-volume requirement as a spatial cell.
    ///
    /// Work is bounded by two original coexistence solves, two source refreshes
    /// and the explicit quadrature budget. Sources are resolved once per phase;
    /// the volume check does not repeat or replace their property evaluation.
    pub fn equilibrium_interface(
        &self,
        component_moles: [f64; 2],
        options: BinaryPlanarInterfaceOptionsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Option<ResolvedBinaryLiquidPlanarInterface>, BinaryLiquidError> {
        poll(&mut cancelled)?;
        let equilibrium = self.liquid.equilibrate(
            self.gradient.temperature(),
            component_moles,
            options.coexistence,
            &mut cancelled,
        )?;
        if equilibrium.phases().len() != 2 {
            poll(&mut cancelled)?;
            return Ok(None);
        }
        for phase in equilibrium.phases() {
            poll(&mut cancelled)?;
            self.validate_phase_volume(phase)?;
        }
        let interface = self
            .gradient
            .equilibrium_interface(options, &mut cancelled)?
            .ok_or(BinaryLiquidError::Arithmetic(
                "two coexisting liquid phases require the same original planar interface",
            ))?;
        poll(&mut cancelled)?;
        Ok(Some(ResolvedBinaryLiquidPlanarInterface {
            equilibrium,
            interface,
        }))
    }
}
