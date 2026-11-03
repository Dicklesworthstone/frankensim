//! Source-resolved liquid coexistence supplies the original tank inventories.

use super::{FiniteHydrostaticLayer, FiniteLayeredHydrostaticTank};
use crate::hydrostatic::{FiniteTankControls, HydrostaticError, HydrostaticTankGeometry, poll};
use fs_material::liquid::solution::ResolvedBinaryLiquidEquilibrium;
use fs_mbd::Gravity;

impl<'a> FiniteLayeredHydrostaticTank<'a> {
    /// Bind the actual one- or two-phase inventory of a binary liquid equilibrium.
    ///
    /// A homogeneous solution supplies one layer. A miscibility gap supplies
    /// two layers at the coexistence compositions, ordered by their resolved
    /// densities rather than species name or composition. Each layer borrows
    /// its original phase owner and uses that owner's volume and exact mass;
    /// mole-based lever fractions are never interpreted as volume fractions.
    /// Equal densities retain the equilibrium owner's deterministic order.
    ///
    /// This is the existing fixed-temperature, incompressible, quasistatic tank
    /// model. Equilibrium and source states remain frozen during mechanics;
    /// hydrostatic pressure variations do not re-evaluate the phase diagram.
    /// Changing temperature/composition requires a new equilibrium and tank.
    /// No mixing kinetics, entrainment, capillarity or thermal work is inferred.
    ///
    /// # Errors
    /// Propagates the original source, inventory, geometry and cancellation
    /// refusals. A failed construction publishes no partial tank.
    pub fn from_binary_equilibrium(
        equilibrium: &'a ResolvedBinaryLiquidEquilibrium,
        gravity: Gravity,
        geometry: HydrostaticTankGeometry,
        controls: FiniteTankControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, HydrostaticError> {
        poll(&mut cancelled)?;
        let (ordered, count) = match equilibrium.phases() {
            [phase] => ([phase, phase], 1),
            [first, second] => {
                if first.liquid().density_kg_m3() < second.liquid().density_kg_m3() {
                    ([second, first], 2)
                } else {
                    ([first, second], 2)
                }
            }
            _ => {
                return Err(HydrostaticError::Invalid(
                    "binary liquid equilibrium requires one or two phases",
                ));
            }
        };
        let inputs = ordered.map(|phase| FiniteHydrostaticLayer {
            source: phase,
            liquid_volume_m3: phase.volume_m3(),
        });
        let mut tank = Self::new(
            &inputs[..count],
            gravity,
            geometry,
            controls,
            &mut cancelled,
        )?;
        // The generic constructor derives mass from rho*V. Preserve the exact
        // conserved inventory supplied by the composition owner instead of a
        // floating-point roundtrip through volume = mass/rho.
        for (layer, phase) in tank.layers.iter_mut().zip(ordered) {
            layer.liquid_mass_kg = phase.mass_kg();
        }
        tank.liquid_mass_kg = equilibrium.total_mass_kg();
        poll(&mut cancelled)?;
        Ok(tank)
    }
}
