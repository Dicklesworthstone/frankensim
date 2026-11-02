//! Conserved liquid mass when an external thermal/EOS owner changes density.

use super::{FiniteHydrostaticTank, FiniteTankFill, HydrostaticError, finite, poll};
use crate::hydrostatic::{DisplacementMesh, FiniteTankBody, FiniteTankSurface};
use fs_ga::Se3;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_mbd::RigidBodyState;

/// Proposed density change at fixed displacement-envelope geometry.
///
/// The source thermal/EOS owner supplies the new state. This accounts for the
/// complete liquid gravitational storage at the fixed tank floor datum, not
/// caloric heat, ambient-pressure work or fluid kinetic energy. The omitted
/// world-elevation constant cancels because tank, gravity and mass are fixed.
///
/// The default record contains one `Se3` and `FiniteTankFill`; a body batch
/// contains ordered `Vec<Se3>` transforms and one shared `FiniteTankSurface`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteTankDensityChange<Pose = Se3, Fill = FiniteTankFill> {
    /// Original liquid inventory, retained bit for bit [kg].
    pub liquid_mass_kg: f64,
    /// Original and proposed source densities [kg/m3].
    pub density_kg_m3: [f64; 2],
    /// Original and proposed liquid volumes, excluding the solid envelope [m3].
    pub liquid_volume_m3: [f64; 2],
    /// The same held body transform(s) used for both native waterline solves.
    /// Batch transforms follow the original body binding order.
    pub world_from_body: Pose,
    /// Original waterline, displacement and relative pressure potential.
    pub fill_before: Fill,
    /// Proposed waterline, displacement and relative pressure potential.
    pub fill_after: Fill,
    /// Change of M*|g|*h_empty/2, required even when the body is dry [J].
    pub empty_tank_gravitational_energy_change_j: f64,
    /// Change of the native finite-bath potential relative to its empty bath [J].
    pub displacement_potential_change_j: f64,
    /// Sum of the two preceding changes, positive for increasing storage [J].
    /// An external source update must supply/account for this energy; it is
    /// not numerical dissipation and is not automatically deposited as heat.
    pub liquid_gravitational_energy_change_j: f64,
}

impl<'a, F: FluidTransport + ?Sized> FiniteHydrostaticTank<'a, F> {
    /// Retain validated tank geometry while consuming one actual source snapshot.
    /// The caller owns the inventory relation: either M=rho*V for a new volume
    /// or V=M/rho with the original mass copied exactly. No source is re-queried.
    pub(in crate::hydrostatic) fn with_snapshot_inventory<'b, G: FluidTransport + ?Sized>(
        &self,
        source: &'b G,
        transport: FluidTransportState,
        liquid_volume_m3: f64,
        liquid_mass_kg: f64,
    ) -> Result<FiniteHydrostaticTank<'b, G>, HydrostaticError> {
        if finite(liquid_volume_m3)? <= 0.0 || finite(liquid_mass_kg)? <= 0.0 {
            return Err(HydrostaticError::Invalid(
                "updated finite tank inventory is unrepresentable",
            ));
        }
        let capacity = finite(self.area_m2 * self.geometry.inside_dimensions_m.z)?;
        if liquid_volume_m3 > capacity {
            return Err(HydrostaticError::TankOverflow {
                excess_volume_m3: liquid_volume_m3 - capacity,
            });
        }
        let empty_height_m = if liquid_volume_m3.to_bits() == self.liquid_volume_m3.to_bits() {
            self.empty_height_m
        } else {
            finite(liquid_volume_m3 / self.area_m2)?
        };
        if empty_height_m <= 0.0 {
            return Err(HydrostaticError::Invalid(
                "updated empty-tank height is unrepresentable",
            ));
        }
        Ok(FiniteHydrostaticTank {
            source,
            transport,
            gravity: self.gravity,
            gravity_magnitude: self.gravity_magnitude,
            geometry: self.geometry,
            area_m2: self.area_m2,
            liquid_volume_m3,
            liquid_mass_kg,
            empty_height_m,
            controls: self.controls,
        })
    }

    /// Propose another source state for the SAME homogeneous liquid inventory.
    ///
    /// The immutable original tank remains available. The new source is
    /// evaluated once through its existing domain/coherence checks and retained
    /// by the returned tank. The tank geometry, gravity, volume-solve controls
    /// and exact liquid mass remain fixed; V_new=M/rho_new sets the next volume.
    /// Equal-density updates preserve the original volume and fill bits.
    ///
    /// `FluidTransport` exposes no material/composition identity. The caller
    /// must bind another admissible state of the same single-phase liquid;
    /// equal Rust source types or numerical tuples cannot prove that identity.
    /// No phase selection, thermal expansion coefficient or EOS is inferred.
    ///
    /// Both waterlines use the existing analytic clipping and bounded root
    /// solve at `world_from_body`. The audit adds the changing empty-bath
    /// reference energy to the displacement potential: subtracting only the
    /// latter would lose the gravitational work of liquid expansion. This is
    /// a quasistatic parameter update, not an energy-free density refresh, a
    /// heat integrator, sloshing, or momentum transfer during fluid motion.
    /// Continue native rigid-body steps using the returned tank.
    ///
    /// # Errors
    /// Source refusal, unrepresentable volume/energy, overflow, invalid body
    /// geometry or cancellation publishes neither a new tank nor a partial audit.
    pub fn rebind_source_at_pose<'b>(
        &self,
        source: &'b F,
        mesh: &DisplacementMesh,
        world_from_body: Se3,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(FiniteHydrostaticTank<'b, F>, FiniteTankDensityChange), HydrostaticError> {
        poll(&mut cancelled)?;
        let (next, same_density, empty_change) = self.rebound_source(source)?;
        let fill_before = self.fill_at_pose(mesh, world_from_body, &mut cancelled)?;
        let fill_after = if same_density {
            fill_before
        } else {
            next.fill_at_pose(mesh, world_from_body, &mut cancelled)?
        };
        let displacement_change = finite(
            fill_after.liquid_gravitational_energy_change_j
                - fill_before.liquid_gravitational_energy_change_j,
        )?;
        let gravitational_change = finite(empty_change + displacement_change)?;
        let audit = FiniteTankDensityChange {
            liquid_mass_kg: self.liquid_mass_kg,
            density_kg_m3: [
                self.transport.density_kg_m3(),
                next.transport.density_kg_m3(),
            ],
            liquid_volume_m3: [self.liquid_volume_m3, next.liquid_volume_m3],
            world_from_body,
            fill_before,
            fill_after,
            empty_tank_gravitational_energy_change_j: empty_change,
            displacement_potential_change_j: displacement_change,
            liquid_gravitational_energy_change_j: gravitational_change,
        };
        poll(&mut cancelled)?;
        Ok((next, audit))
    }

    /// Change the source of one conserved liquid shared by ALL held bodies.
    ///
    /// The same source/inventory owner as [`Self::rebind_source_at_pose`]
    /// retains exact mass and queries the proposed source once. Each waterline
    /// solves `A*h - sum(V_sub,b) = M/rho`; the full fluid potential therefore
    /// retains the cross-body displacement terms. Body transforms, bounds and
    /// one cuts buffer are prepared once and reused for both density states.
    /// No unused individual wrenches or independently filled tanks are built.
    /// Equal-density updates reuse the complete original surface record.
    ///
    /// Bodies and momenta remain held during this quasistatic parameter change.
    /// Continue native [`Self::step_bodies`] motion with the returned tank;
    /// its external thermal/EOS owner accounts for the complete gravitational
    /// change together with its caloric and ambient-pressure work. The original
    /// same-liquid identity, nonoverlap and fixed-tank obligations still apply.
    /// Failure returns neither a new tank nor a partial geometry/energy audit.
    pub fn rebind_source_at_bodies<'b>(
        &self,
        source: &'b F,
        bodies: &[FiniteTankBody<'_>],
        states: &[RigidBodyState],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<
        (
            FiniteHydrostaticTank<'b, F>,
            FiniteTankDensityChange<Vec<Se3>, FiniteTankSurface>,
        ),
        HydrostaticError,
    > {
        poll(&mut cancelled)?;
        let (next, same_density, empty_change) = self.rebound_source(source)?;
        let mut prepared = self.prepare_bodies(bodies, states, &mut cancelled)?;
        let mut world_from_body = Vec::with_capacity(bodies.len());
        for (body, &state) in bodies.iter().zip(states) {
            poll(&mut cancelled)?;
            world_from_body.push(body.frame.world_from_body(state)?);
        }
        let fill_before = self.fill_prepared_bodies(&mut prepared, &mut cancelled)?;
        let fill_after = if same_density {
            fill_before
        } else {
            next.fill_prepared_bodies(&mut prepared, &mut cancelled)?
        };
        let displacement_change = finite(
            fill_after.liquid_gravitational_energy_change_j
                - fill_before.liquid_gravitational_energy_change_j,
        )?;
        let gravitational_change = finite(empty_change + displacement_change)?;
        let audit = FiniteTankDensityChange {
            liquid_mass_kg: self.liquid_mass_kg,
            density_kg_m3: [
                self.transport.density_kg_m3(),
                next.transport.density_kg_m3(),
            ],
            liquid_volume_m3: [self.liquid_volume_m3, next.liquid_volume_m3],
            world_from_body,
            fill_before,
            fill_after,
            empty_tank_gravitational_energy_change_j: empty_change,
            displacement_potential_change_j: displacement_change,
            liquid_gravitational_energy_change_j: gravitational_change,
        };
        poll(&mut cancelled)?;
        Ok((next, audit))
    }

    /// One source snapshot, exact inventory and empty-bath datum for either
    /// original-pose or shared-body geometry. Caller polls before the query.
    fn rebound_source<'b>(
        &self,
        source: &'b F,
    ) -> Result<(FiniteHydrostaticTank<'b, F>, bool, f64), HydrostaticError> {
        let transport = source.transport_state()?;
        let old_density = self.transport.density_kg_m3();
        let density = transport.density_kg_m3();
        let same_density = density.to_bits() == old_density.to_bits();
        let liquid_volume_m3 = if same_density {
            self.liquid_volume_m3
        } else {
            finite(self.liquid_mass_kg / density)?
        };
        let next = self.with_snapshot_inventory(
            source,
            transport,
            liquid_volume_m3,
            // Never round-trip the inventory through rho_new*V_new.
            self.liquid_mass_kg,
        )?;
        let empty_change = if same_density {
            // A zero change must not first form an overflowing absolute scale.
            0.0
        } else {
            finite(
                (0.5 * self.liquid_mass_kg * self.gravity_magnitude)
                    * (next.empty_height_m - self.empty_height_m),
            )?
        };
        Ok((next, same_density, empty_change))
    }
}
