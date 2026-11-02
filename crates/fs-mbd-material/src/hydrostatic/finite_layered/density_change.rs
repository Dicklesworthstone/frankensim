//! Source updates of closed immiscible liquids and their coupled interfaces.

use super::{
    FiniteHydrostaticLayer, FiniteLayeredHydrostaticTank, FiniteTankBodiesInterface,
    FiniteTankInterface, ResolvedFiniteLayer, set_density_jumps,
};
use crate::hydrostatic::{
    DisplacementMesh, FiniteTankBody, HydrostaticError, finite, gravity_direction, poll,
};
use fs_ga::Se3;
use fs_material::fluid_transport::FluidTransport;
use fs_mbd::RigidBodyState;

/// One closed liquid's original mass and changed source-defined volume.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteLayerDensityChange {
    /// Conserved original inventory, copied bit for bit [kg].
    pub liquid_mass_kg: f64,
    /// Original and proposed source densities [kg/m3].
    pub density_kg_m3: [f64; 2],
    /// Original and proposed physical liquid volumes [m3].
    pub liquid_volume_m3: [f64; 2],
}

/// Complete gravitational-storage change at held body geometry and a fixed tank.
/// The original and returned tanks retain their full sources and snapshots.
/// Defaults retain one `Se3` and `FiniteTankInterface`; a batch carries ordered
/// `Vec<Se3>` transforms and original-body `FiniteTankBodiesInterface` records.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteLayeredTankDensityChange<Pose = Se3, Interface = FiniteTankInterface> {
    /// Original total liquid mass, copied without another sum [kg].
    pub liquid_mass_kg: f64,
    /// Per-liquid changes in the unchanged bottom-to-top order.
    pub layers: Vec<FiniteLayerDensityChange>,
    /// Same original-body transform(s) used by both interface queries.
    /// Batch transforms follow the original body binding order.
    pub world_from_body: Pose,
    /// Original interface locations and cumulative volume residuals.
    pub interfaces_before: Vec<Interface>,
    /// Proposed interface locations and cumulative volume residuals.
    pub interfaces_after: Vec<Interface>,
    /// Change of sum_i M_i*|g|*(h_empty,i-1+h_empty,i)/2 [J].
    /// Expansion of a lower liquid also lifts the liquids above it.
    pub empty_tank_gravitational_energy_change_j: f64,
    /// Change of the native layered potential relative to its empty bath [J].
    pub displacement_potential_change_j: f64,
    /// Complete liquid gravitational change, the sum of the preceding terms [J].
    /// External source-state evolution must account for this energy; it is
    /// neither dissipation nor automatically supplied caloric heat.
    pub liquid_gravitational_energy_change_j: f64,
}

impl FiniteLayeredHydrostaticTank<'_> {
    /// Supply one new state of EACH SAME immiscible liquid in retained order.
    ///
    /// Each source is queried once; every original layer mass and the total
    /// mass are copied exactly. V_i=M_i/rho_i sets each new physical volume,
    /// and cumulative volumes drive the original finite-tank interface solves.
    /// An unchanged density retains that layer's volume bits. If all densities
    /// are unchanged, interface records are reused and gravitational change is zero.
    ///
    /// The full gravitational audit includes the changing empty-stack datum.
    /// A lower liquid's expansion can lift every upper liquid even when the
    /// displacement envelope is dry. World elevation cancels because the tank,
    /// gravity and EACH material mass are unchanged.
    ///
    /// `FluidTransport` cannot establish composition identity or immiscibility.
    /// Binding the same physical liquids is a caller obligation. Source owners
    /// supply thermal/EOS closure and account for caloric/ambient-pressure work;
    /// no heat, phase change, mixing or fluid momentum transfer is inferred.
    /// The returned tank continues through the existing native load and step.
    ///
    /// # Errors
    /// A changed layer count, source refusal, density-order inversion, lost
    /// volume, overflow, invalid geometry, unrepresentable energy or cancellation
    /// returns neither a new tank nor a partial audit. Liquids are never reordered.
    pub fn rebind_sources_at_pose<'b>(
        &self,
        sources: &[&'b dyn FluidTransport],
        mesh: &DisplacementMesh,
        world_from_body: Se3,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<
        (
            FiniteLayeredHydrostaticTank<'b>,
            FiniteLayeredTankDensityChange,
        ),
        HydrostaticError,
    > {
        let (next, changes, all_same_density, empty_change) =
            self.rebound_sources(sources, &mut cancelled)?;
        let (interfaces_before, potential_before) =
            self.interfaces_and_potential_at_pose(mesh, world_from_body, &mut cancelled)?;
        let (interfaces_after, potential_after) = if all_same_density {
            (interfaces_before.clone(), potential_before)
        } else {
            next.interfaces_and_potential_at_pose(mesh, world_from_body, &mut cancelled)?
        };
        let displacement_change = finite(potential_after - potential_before)?;
        let gravitational_change = finite(empty_change + displacement_change)?;
        poll(&mut cancelled)?;
        Ok((
            next,
            FiniteLayeredTankDensityChange {
                liquid_mass_kg: self.liquid_mass_kg,
                layers: changes,
                world_from_body,
                interfaces_before,
                interfaces_after,
                empty_tank_gravitational_energy_change_j: empty_change,
                displacement_potential_change_j: displacement_change,
                liquid_gravitational_energy_change_j: gravitational_change,
            },
        ))
    }

    /// Change every closed liquid source while ALL original bodies stay held.
    ///
    /// The same source/inventory preparation as [`Self::rebind_sources_at_pose`]
    /// retains each exact material mass and sets `V_i=M_i/rho_i` once. Both
    /// complete stacks solve `A*h_i - sum(V_sub,b(h_i)) = sum(V_liquid,<=i)`.
    /// Shared interfaces and their potential are computed together using one
    /// prepared geometry/cuts buffer per stack; no force assembly or extra
    /// per-body tank solve is required. Equal-density refreshes retain the
    /// original interface records and zero complete gravitational change.
    ///
    /// The audit includes changing empty-stack storage and displacement cross
    /// terms. Thus a lower liquid's expansion also lifts every upper inventory.
    /// The external source owner supplies/accounts for that work alongside its
    /// caloric and ambient-pressure terms. Bodies retain their poses and momenta
    /// across the update; [`Self::step_bodies`] continues with the proposed tank.
    /// The original same-liquid, immiscibility, nonoverlap and source-domain
    /// obligations still apply. No thermal step, mixing, fluid inertia or
    /// expansion-driven momentum impulse is inferred by this quasistatic model.
    /// Source/order/geometry/overflow/cancellation refusal returns no candidate.
    pub fn rebind_sources_at_bodies<'b>(
        &self,
        sources: &[&'b dyn FluidTransport],
        bodies: &[FiniteTankBody<'_>],
        states: &[RigidBodyState],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<
        (
            FiniteLayeredHydrostaticTank<'b>,
            FiniteLayeredTankDensityChange<Vec<Se3>, FiniteTankBodiesInterface>,
        ),
        HydrostaticError,
    > {
        let (next, changes, all_same_density, empty_change) =
            self.rebound_sources(sources, &mut cancelled)?;
        let (interfaces_before, potential_before) =
            self.interfaces_and_potential_at_bodies(bodies, states, &mut cancelled)?;
        let mut world_from_body = Vec::with_capacity(bodies.len());
        for (body, &state) in bodies.iter().zip(states) {
            poll(&mut cancelled)?;
            world_from_body.push(body.frame.world_from_body(state)?);
        }
        let (interfaces_after, potential_after) = if all_same_density {
            (interfaces_before.clone(), potential_before)
        } else {
            next.interfaces_and_potential_at_bodies(bodies, states, &mut cancelled)?
        };
        let displacement_change = finite(potential_after - potential_before)?;
        let gravitational_change = finite(empty_change + displacement_change)?;
        poll(&mut cancelled)?;
        Ok((
            next,
            FiniteLayeredTankDensityChange {
                liquid_mass_kg: self.liquid_mass_kg,
                layers: changes,
                world_from_body,
                interfaces_before,
                interfaces_after,
                empty_tank_gravitational_energy_change_j: empty_change,
                displacement_potential_change_j: displacement_change,
                liquid_gravitational_energy_change_j: gravitational_change,
            },
        ))
    }

    /// Shared material inventory/source and empty-stack gravitational owner.
    fn rebound_sources<'b>(
        &self,
        sources: &[&'b dyn FluidTransport],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<
        (
            FiniteLayeredHydrostaticTank<'b>,
            Vec<FiniteLayerDensityChange>,
            bool,
            f64,
        ),
        HydrostaticError,
    > {
        poll(cancelled)?;
        if sources.len() != self.layers.len() {
            return Err(HydrostaticError::Invalid(
                "finite liquid source update must retain every layer in order",
            ));
        }
        let mut resolved = Vec::with_capacity(self.layers.len());
        let mut changes = Vec::with_capacity(self.layers.len());
        let mut volume = 0.0;
        let mut previous_density = f64::INFINITY;
        let mut all_same_density = true;
        for (old, &source) in self.layers.iter().zip(sources) {
            poll(cancelled)?;
            let transport = source.transport_state()?;
            let density = transport.density_kg_m3();
            if density > previous_density {
                return Err(HydrostaticError::Invalid(
                    "updated finite liquid density must not increase upward",
                ));
            }
            let old_density = old.cumulative.transport().density_kg_m3();
            let same_density = density.to_bits() == old_density.to_bits();
            all_same_density &= same_density;
            let liquid_volume_m3 = if same_density {
                old.input.liquid_volume_m3
            } else {
                finite(old.liquid_mass_kg / density)?
            };
            let cumulative_volume = finite(volume + liquid_volume_m3)?;
            if liquid_volume_m3 <= 0.0 || cumulative_volume <= volume {
                return Err(HydrostaticError::Invalid(
                    "updated finite liquid layer volume is unrepresentable",
                ));
            }
            let cumulative = old.cumulative.with_snapshot_inventory(
                source,
                transport,
                cumulative_volume,
                finite(density * cumulative_volume)?,
            )?;
            resolved.push(ResolvedFiniteLayer {
                input: FiniteHydrostaticLayer {
                    source,
                    liquid_volume_m3,
                },
                cumulative,
                liquid_mass_kg: old.liquid_mass_kg,
                density_jump: 0.0,
            });
            changes.push(FiniteLayerDensityChange {
                liquid_mass_kg: old.liquid_mass_kg,
                density_kg_m3: [old_density, density],
                liquid_volume_m3: [old.input.liquid_volume_m3, liquid_volume_m3],
            });
            volume = cumulative_volume;
            previous_density = density;
        }
        set_density_jumps(&mut resolved);
        let next = FiniteLayeredHydrostaticTank {
            layers: resolved,
            liquid_mass_kg: self.liquid_mass_kg,
        };
        let mut empty_change = 0.0;
        if !all_same_density {
            let sides = self.geometry().inside_dimensions_m;
            let area = sides.x * sides.y;
            let (_, gravity_magnitude) = gravity_direction(self.gravity())?;
            let mut lower_change = 0.0;
            for (old, new) in self.layers.iter().zip(&next.layers) {
                poll(cancelled)?;
                let upper_change = new.cumulative.liquid_volume_m3() / area
                    - old.cumulative.liquid_volume_m3() / area;
                let twice_center_change = lower_change + upper_change;
                if twice_center_change != 0.0 {
                    empty_change = finite(
                        empty_change
                            + (0.5 * old.liquid_mass_kg * gravity_magnitude) * twice_center_change,
                    )?;
                }
                lower_change = upper_change;
            }
        }
        Ok((next, changes, all_same_density, empty_change))
    }
}
