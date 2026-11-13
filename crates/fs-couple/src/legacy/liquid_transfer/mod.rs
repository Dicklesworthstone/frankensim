//! Conservative matter transfer between two original moving-liquid owners.
//!
//! A donor parcel leaves at its actual Guo-corrected velocity. The original
//! receiver mixes that momentum into its populations and deposits the lost
//! relative kinetic energy as heat exactly once. Removing mass at the donor's
//! own velocity requires no recoil actuator or mechanical work.
//!
//! This is an instantaneous transfer between coincident interface cells in
//! overlapping representations, not a resolved pipe, jet or phase transition.
//! Each grid has an explicit origin and orientation in one inertial frame;
//! their mapped cell centers must coincide. Thus the stream has one application
//! point, with no missing connector torque or gravitational transport work.
//! Coincidence uses compensated origin differences and rotated local offsets,
//! with roundoff bounded by the local geometry, never absolute world position.
//! A large common translation therefore cannot hide a real separation.
//! Owners retain the same original caloric object, complete material claims
//! and external context. Their temperature references, physical densities and
//! lattice scalings may differ. The accepted receiver still has to satisfy its
//! own storage-density law; no thermodynamic density is imposed as pressure.
//! Ordinary native steps advance their synchronized physical clocks.
//!
//! The audited energy is original thermal enthalpy plus resolved kinetic
//! energy. No compressible pressure-work, capillary-energy, remeshing accuracy
//! or universal total-energy claim is added to the original native model.

/// Remove a no-storage rigid constraint through the original thermal VOF handoff.
pub mod rigid;

use fs_exec::CancelGate;
use fs_lbm::d3q19::freesurface3::FreeSurfaceError3;
use fs_lbm::material::free_surface::enthalpy::{
    FreeSurfaceEnthalpyError as Error, FreeSurfaceMassExchange, FreeSurfaceMassExchangeControls,
    FreeSurfaceMassInflow, FreeSurfaceMassInflowReport, FreeSurfaceMassOutflow,
    FreeSurfaceMassOutflowReport, MaterialFreeSurfaceEnthalpy3,
};
use fs_math::dd::Dd;
use fs_mbd::{UnitQuaternion, Vec3};

/// One finite packet moved from donor to receiver, without advancing time.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LiquidMassTransfer {
    /// Already occupied native donor interface; complete exhaustion is refused.
    pub donor: [usize; 3],
    /// Already occupied receiver interface; overflow/nucleation is refused.
    pub receiver: [usize; 3],
    /// World location of the donor grid's lower corner [m].
    pub donor_origin_m: [f64; 3],
    /// World location of the receiver grid's lower corner [m], in the same frame.
    pub receiver_origin_m: [f64; 3],
    /// Fixed donor-grid axes mapped into the common inertial frame. Native
    /// velocities and forces remain expressed in their original grid axes.
    pub donor_axes_to_world: UnitQuaternion,
    /// Fixed receiver-grid axes mapped into the same inertial frame. This is
    /// a coordinate rotation, not a translating or rotating reference frame.
    pub receiver_axes_to_world: UnitQuaternion,
    /// Positive requested packet mass [kg]. Actual native represented amounts
    /// remain explicit in both reports and enter all joint balance gates.
    pub mass_kg: f64,
}

/// Actual native state changes after a jointly accepted transfer.
#[derive(Clone, Debug, PartialEq)]
pub struct LiquidMassTransferReport {
    /// Common accepted clock [s]; neither original timeline advances here.
    pub physical_time_s: f64,
    /// Common mapped application point [m].
    pub point_m: [f64; 3],
    /// Native withdrawal, including carried momentum in DONOR grid axes and
    /// kinetic energy.
    pub donor: FreeSurfaceMassOutflowReport,
    /// Native momentum mixing in RECEIVER grid axes and its single physical
    /// heat deposition.
    pub receiver: FreeSurfaceMassInflowReport,
    /// Signed sum of the two actual represented mass changes [kg].
    pub mass_defect_kg: f64,
    /// Signed sum of actual native endpoint momentum changes, both rotated
    /// into the common WORLD axes [kg m/s].
    pub momentum_defect_kg_m_s: [f64; 3],
    /// Signed sum of actual enthalpy and kinetic changes [J], not a pressure
    /// work/internal-energy or capillary-energy theorem.
    pub enthalpy_kinetic_defect_j: f64,
    pub reference_volumes: LiquidMassTransferVolumes,
}

/// Reference-volume inventory carried between the original physical scalings.
/// These use actual represented kg divided by each owner's reference kg/m3.
/// They are not pressure-density-derived VOF geometry, an inferred pressure
/// impulse, or a replacement for an unmodeled thermal expansion solve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LiquidMassTransferVolumes {
    pub withdrawn_m3: f64,
    pub received_m3: f64,
    pub net_change_m3: f64,
}

/// Several finite original donors mixed at ONE coincident receiver interface.
#[derive(Clone, Debug, PartialEq)]
pub struct LiquidMassTransferBatchReport {
    pub physical_time_s: f64,
    pub point_m: [f64; 3],
    /// One actual withdrawal per authored donor/transfer pair, in that order.
    pub donors: Vec<FreeSurfaceMassOutflowReport>,
    /// One native multistream mixture, including stream velocity dispersion.
    pub receiver: FreeSurfaceMassInflowReport,
    pub mass_defect_kg: f64,
    pub momentum_defect_kg_m_s: [f64; 3],
    pub enthalpy_kinetic_defect_j: f64,
    pub reference_volumes: LiquidMassTransferVolumes,
}

/// Transfer actual donor matter into the original receiver's momentum state.
///
/// `controls` bounds each original native transaction and the combined actual
/// mass, componentwise-L1 momentum and enthalpy-plus-kinetic residuals. Receiver
/// credit uses the donor's represented withdrawal, not its requested amount;
/// any remaining receiver representation error enters the joint balance.
/// Original source, capacity, speed, force-frame and thermal-history admission
/// remains active. Any refusal or cancellation preserves BOTH original owners.
/// The physical donor velocity is rotated into the receiver's axes before
/// mixing; local momentum increments are rotated into world axes before their
/// joint check. No grid resampling, moving-frame work or extra mixing heat is
/// introduced. Original native per-owner reports retain their local axes.
/// Different thermal references/densities are admitted only when the complete
/// original material claims and caloric object agree and the receiver's mixed
/// endpoint still fits its own reference-density model. Its population density
/// is never reset to make an otherwise unsupported mixture pass.
/// Each liquid builds only its original native proposal; neither complete
/// trajectory nor its unchanged material/history bindings are cloned.
pub fn transfer_interface_mass(
    donor: &mut MaterialFreeSurfaceEnthalpy3<'_, '_>,
    receiver: &mut MaterialFreeSurfaceEnthalpy3<'_, '_>,
    transfer: LiquidMassTransfer,
    controls: FreeSurfaceMassExchangeControls,
    gate: &CancelGate,
) -> Result<LiquidMassTransferReport, Error> {
    poll(gate)?;
    if !transfer.mass_kg.is_finite()
        || transfer.mass_kg <= 0.0
        || [
            controls.mass_tolerance_kg,
            controls.energy_tolerance_j,
            controls.momentum_tolerance_kg_m_s,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
    {
        return Err(Error::Invalid(
            "liquid transfer needs positive mass and finite nonnegative budgets",
        ));
    }
    if !donor.shares_mass_transfer_material(receiver) {
        return Err(Error::Invalid(
            "liquid transfer must retain the original caloric object, complete source claims and nonthermal context",
        ));
    }
    let physical_time_s = donor.physical_time_s();
    if !physical_time_s.is_finite() || physical_time_s != receiver.physical_time_s() {
        return Err(Error::Invalid(
            "liquid transfer requires synchronized accepted physical clocks",
        ));
    }
    // Project before preparing either trajectory: reject exhaustion, overflow,
    // zero-mass ghosts and unrepresentable transfer increments at the owners.
    let outgoing = -donor
        .project_interface_mass_exchange(transfer.donor, -transfer.mass_kg)?
        .mass_change_kg;
    let incoming = receiver
        .project_interface_mass_exchange(transfer.receiver, outgoing)?
        .mass_change_kg;
    if (incoming - outgoing).abs() > controls.mass_tolerance_kg {
        return Err(Error::Invalid(
            "native liquid packet masses disagree beyond the declared budget",
        ));
    }
    let donor_center = center(
        donor,
        transfer.donor,
        transfer.donor_origin_m,
        transfer.donor_axes_to_world,
    )?;
    let receiver_center = center(
        receiver,
        transfer.receiver,
        transfer.receiver_origin_m,
        transfer.receiver_axes_to_world,
    )?;
    if !donor_center.coincides(&receiver_center) {
        return Err(Error::Invalid(
            "liquid transfer requires coincident mapped cell centers in one Cartesian frame",
        ));
    }
    let point_m = donor_center.point_m;
    let [x, y, z] = transfer.donor;
    let velocity_m_s = donor
        .flow()
        .sample_at_si(x, y, z)?
        .ok_or(Error::Invalid("liquid transfer requires an occupied donor"))?
        .velocity_m_s;
    // Identical axes need no round-trip rotation. In particular the existing
    // parallel-grid path keeps its original population and mixing-heat bits.
    let receiver_velocity_m_s = if transfer.donor_axes_to_world == transfer.receiver_axes_to_world {
        velocity_m_s
    } else {
        vector(
            transfer.receiver_axes_to_world.rotate_world_to_body(
                transfer
                    .donor_axes_to_world
                    .rotate_body_to_world(vec3(velocity_m_s)),
            ),
        )
    };
    let specific_enthalpy_j_kg =
        donor
            .specific_enthalpy_at(transfer.donor)
            .ok_or(Error::Invalid(
                "liquid transfer requires retained donor enthalpy",
            ))?;
    let donor_density = donor.flow().scaling().material().density_kg_m3();
    let receiver_density = receiver.flow().scaling().material().density_kg_m3();
    let donor_proposal = donor.prepare_interface_outflow(
        &[FreeSurfaceMassOutflow {
            exchange: FreeSurfaceMassExchange {
                coordinates: transfer.donor,
                mass_kg: -transfer.mass_kg,
                specific_enthalpy_j_kg,
                heat_j: 0.0,
            },
            velocity_m_s,
        }],
        // The native owner samples the identical velocity_si operation; the
        // relative exhaust velocity and authored work are exactly zero.
        0.0,
        controls,
        gate,
    )?;
    let receiver_proposal = receiver.prepare_interface_inflow(
        &[FreeSurfaceMassInflow {
            exchange: FreeSurfaceMassExchange {
                coordinates: transfer.receiver,
                mass_kg: -donor_proposal.report().exchange.mass_change_kg,
                specific_enthalpy_j_kg,
                heat_j: 0.0,
            },
            velocity_m_s: receiver_velocity_m_s,
        }],
        controls,
        gate,
    )?;
    let withdrawal = donor_proposal.report();
    let mixing = receiver_proposal.report();
    let mass_defect_kg = withdrawal.exchange.mass_change_kg + mixing.exchange.mass_change_kg;
    let donor_change = transfer
        .donor_axes_to_world
        .rotate_body_to_world(vec3(withdrawal.exchange.momentum_change_kg_m_s));
    let receiver_change = transfer
        .receiver_axes_to_world
        .rotate_body_to_world(vec3(mixing.exchange.momentum_change_kg_m_s));
    let momentum_defect_kg_m_s = vector(donor_change.add(receiver_change));
    // Sum the opposite transported enthalpies before the smaller kinetic
    // changes, retaining their cancellation without a large absolute store.
    let enthalpy_kinetic_defect_j = sum([
        withdrawal.exchange.enthalpy_change_j,
        mixing.exchange.enthalpy_change_j,
        withdrawal.kinetic_energy_change_j,
        mixing.kinetic_energy_change_j,
    ]);
    if !mass_defect_kg.is_finite()
        || mass_defect_kg.abs() > controls.mass_tolerance_kg
        || momentum_defect_kg_m_s.iter().any(|x| !x.is_finite())
        || momentum_defect_kg_m_s.iter().map(|x| x.abs()).sum::<f64>()
            > controls.momentum_tolerance_kg_m_s
        || !enthalpy_kinetic_defect_j.is_finite()
        || enthalpy_kinetic_defect_j.abs() > controls.energy_tolerance_j
    {
        return Err(Error::Invalid(
            "joint native liquid mass, momentum or enthalpy-plus-kinetic balance exceeded budget",
        ));
    }
    let reference_volumes = reference_volumes(
        -withdrawal.exchange.mass_change_kg / donor_density,
        mixing.exchange.mass_change_kg / receiver_density,
    )?;
    poll(gate)?;
    // Both proposals exclusively borrow their original owners. No state can
    // become stale, and publication after the joint check cannot fail halfway.
    let withdrawal = donor_proposal.publish();
    let mixing = receiver_proposal.publish();
    Ok(LiquidMassTransferReport {
        physical_time_s,
        point_m,
        donor: withdrawal,
        receiver: mixing,
        mass_defect_kg,
        momentum_defect_kg_m_s,
        enthalpy_kinetic_defect_j,
        reference_volumes,
    })
}

/// Mix finite streams from distinct original donors in one atomic transaction.
///
/// Each donor has one authored packet; all packets meet the same receiver cell,
/// pose and world point. The existing multistream receiver combines them ONCE,
/// retaining both their mean momentum and relative kinetic variance as physical
/// mixing heat. Sequential pairwise publication would reject a temporary
/// unsupported temperature or change the inelastic mixing result.
///
/// Temperature/reference density and dx/dt may differ across original owners.
/// The exact original caloric object, every optional source claim, typed axes,
/// nonthermal conditions, synchronized clock and each local native admission
/// remain required. In particular, an accepted mixture needing an unmodeled
/// receiver density/volume change refuses the whole batch; balanced hot/cold
/// feeds or an admitted constant-density part of the source law can pass.
///
/// `maximum_streams` bounds authored participants and proposal storage. Controls
/// bound every native transaction and the joint actual mass, world momentum and
/// enthalpy-plus-kinetic defects. The source owners debit their actual native
/// withdrawals; receiver rounding remains an explicit joint defect. No caller
/// trajectory is cloned, no clock advances, and no original participant changes
/// until every candidate and the final cancellation poll have passed.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn transfer_interface_mass_streams<'s, 'm>(
    donors: &mut [&mut MaterialFreeSurfaceEnthalpy3<'s, 'm>],
    receiver: &mut MaterialFreeSurfaceEnthalpy3<'_, '_>,
    transfers: &[LiquidMassTransfer],
    maximum_streams: usize,
    controls: FreeSurfaceMassExchangeControls,
    gate: &CancelGate,
) -> Result<LiquidMassTransferBatchReport, Error> {
    poll(gate)?;
    if donors.is_empty()
        || donors.len() != transfers.len()
        || donors.len() > maximum_streams
        || [
            controls.mass_tolerance_kg,
            controls.energy_tolerance_j,
            controls.momentum_tolerance_kg_m_s,
        ]
        .into_iter()
        .any(|x| !x.is_finite() || x < 0.0)
    {
        return Err(Error::Invalid(
            "liquid stream batch requires bounded matching donors and finite budgets",
        ));
    }
    let receiver_pose = transfers[0];
    let physical_time_s = receiver.physical_time_s();
    if !physical_time_s.is_finite() {
        return Err(Error::Invalid(
            "liquid stream batch needs a finite accepted clock",
        ));
    }
    let receiver_center = center(
        receiver,
        receiver_pose.receiver,
        receiver_pose.receiver_origin_m,
        receiver_pose.receiver_axes_to_world,
    )?;
    let point_m = receiver_center.point_m;
    let receiver_density = receiver.flow().scaling().material().density_kg_m3();
    let mut outgoing = Vec::new();
    let mut incoming = Vec::new();
    let mut donor_densities = Vec::new();
    outgoing
        .try_reserve_exact(donors.len())
        .map_err(|_| Error::Invalid("liquid stream withdrawal allocation refused"))?;
    incoming
        .try_reserve_exact(donors.len())
        .map_err(|_| Error::Invalid("liquid stream inflow allocation refused"))?;
    donor_densities
        .try_reserve_exact(donors.len())
        .map_err(|_| Error::Invalid("liquid stream source allocation refused"))?;
    let mut total_incoming = 0.0;
    for (donor, transfer) in donors.iter().zip(transfers) {
        poll(gate)?;
        if !transfer.mass_kg.is_finite()
            || transfer.mass_kg <= 0.0
            || transfer.receiver != receiver_pose.receiver
            || transfer.receiver_origin_m != receiver_pose.receiver_origin_m
            || transfer.receiver_axes_to_world != receiver_pose.receiver_axes_to_world
            || !donor.shares_mass_transfer_material(receiver)
            || donor.physical_time_s() != physical_time_s
        {
            return Err(Error::Invalid(
                "liquid streams changed original source, clock or common receiver",
            ));
        }
        let donor_center = center(
            donor,
            transfer.donor,
            transfer.donor_origin_m,
            transfer.donor_axes_to_world,
        )?;
        if !donor_center.coincides(&receiver_center) {
            return Err(Error::Invalid(
                "liquid stream centers must coincide with the receiver in one world frame",
            ));
        }
        let mass = -donor
            .project_interface_mass_exchange(transfer.donor, -transfer.mass_kg)?
            .mass_change_kg;
        total_incoming += mass;
        let p = transfer.donor;
        let velocity_m_s = donor
            .flow()
            .sample_at_si(p[0], p[1], p[2])?
            .ok_or(Error::Invalid("liquid stream requires an occupied donor"))?
            .velocity_m_s;
        let receiver_velocity = if transfer.donor_axes_to_world == transfer.receiver_axes_to_world {
            velocity_m_s
        } else {
            vector(
                transfer.receiver_axes_to_world.rotate_world_to_body(
                    transfer
                        .donor_axes_to_world
                        .rotate_body_to_world(vec3(velocity_m_s)),
                ),
            )
        };
        let h = donor.specific_enthalpy_at(p).ok_or(Error::Invalid(
            "liquid stream lost its original donor enthalpy",
        ))?;
        outgoing.push(FreeSurfaceMassOutflow {
            exchange: FreeSurfaceMassExchange {
                coordinates: p,
                mass_kg: -transfer.mass_kg,
                specific_enthalpy_j_kg: h,
                heat_j: 0.0,
            },
            velocity_m_s,
        });
        incoming.push(FreeSurfaceMassInflow {
            exchange: FreeSurfaceMassExchange {
                coordinates: transfer.receiver,
                mass_kg: mass,
                specific_enthalpy_j_kg: h,
                heat_j: 0.0,
            },
            velocity_m_s: receiver_velocity,
        });
        donor_densities.push(donor.flow().scaling().material().density_kg_m3());
    }
    // Collective capacity, before allocating a native proposal for any donor.
    receiver.project_interface_mass_exchange(receiver_pose.receiver, total_incoming)?;
    let mut prepared = Vec::new();
    prepared
        .try_reserve_exact(donors.len())
        .map_err(|_| Error::Invalid("liquid stream proposal allocation refused"))?;
    for ((donor, outflow), inflow) in donors.iter_mut().zip(&outgoing).zip(&mut incoming) {
        poll(gate)?;
        let candidate =
            donor.prepare_interface_outflow(core::slice::from_ref(outflow), 0.0, controls, gate)?;
        inflow.exchange.mass_kg = -candidate.report().exchange.mass_change_kg;
        prepared.push(candidate);
    }
    let receiver_proposal = receiver.prepare_interface_inflow(&incoming, controls, gate)?;
    let mixing = receiver_proposal.report();
    let mass_defect_kg = sum(prepared
        .iter()
        .map(|p| p.report().exchange.mass_change_kg)
        .chain(core::iter::once(mixing.exchange.mass_change_kg)));
    let receiver_change = vector(
        receiver_pose
            .receiver_axes_to_world
            .rotate_body_to_world(vec3(mixing.exchange.momentum_change_kg_m_s)),
    );
    let momentum_defect_kg_m_s = core::array::from_fn(|axis| {
        sum(prepared
            .iter()
            .zip(transfers)
            .map(|(p, transfer)| {
                vector(
                    transfer
                        .donor_axes_to_world
                        .rotate_body_to_world(vec3(p.report().exchange.momentum_change_kg_m_s)),
                )[axis]
            })
            .chain(core::iter::once(receiver_change[axis])))
    });
    // Cancel the transported caloric datums before the smaller kinetic terms.
    let enthalpy_kinetic_defect_j = sum(prepared
        .iter()
        .map(|p| p.report().exchange.enthalpy_change_j)
        .chain(core::iter::once(mixing.exchange.enthalpy_change_j))
        .chain(prepared.iter().map(|p| p.report().kinetic_energy_change_j))
        .chain(core::iter::once(mixing.kinetic_energy_change_j)));
    if !mass_defect_kg.is_finite()
        || mass_defect_kg.abs() > controls.mass_tolerance_kg
        || momentum_defect_kg_m_s.iter().any(|x| !x.is_finite())
        || momentum_defect_kg_m_s.iter().map(|x| x.abs()).sum::<f64>()
            > controls.momentum_tolerance_kg_m_s
        || !enthalpy_kinetic_defect_j.is_finite()
        || enthalpy_kinetic_defect_j.abs() > controls.energy_tolerance_j
    {
        return Err(Error::Invalid(
            "joint liquid stream mass, momentum or enthalpy-plus-kinetic balance exceeded budget",
        ));
    }
    let reference_volumes = reference_volumes(
        sum(prepared
            .iter()
            .zip(&donor_densities)
            .map(|(p, density)| -p.report().exchange.mass_change_kg / density)),
        mixing.exchange.mass_change_kg / receiver_density,
    )?;
    // Reserve receipt storage before the first infallible native publication.
    let mut reports = Vec::new();
    reports
        .try_reserve_exact(prepared.len())
        .map_err(|_| Error::Invalid("liquid stream receipt allocation refused"))?;
    poll(gate)?;
    for candidate in prepared {
        reports.push(candidate.publish());
    }
    let mixing = receiver_proposal.publish();
    Ok(LiquidMassTransferBatchReport {
        physical_time_s,
        point_m,
        donors: reports,
        receiver: mixing,
        mass_defect_kg,
        momentum_defect_kg_m_s,
        enthalpy_kinetic_defect_j,
        reference_volumes,
    })
}

fn reference_volumes(
    withdrawn_m3: f64,
    received_m3: f64,
) -> Result<LiquidMassTransferVolumes, Error> {
    let net_change_m3 = received_m3 - withdrawn_m3;
    if [withdrawn_m3, received_m3, net_change_m3]
        .into_iter()
        .any(|v| !v.is_finite())
        || withdrawn_m3 <= 0.0
        || received_m3 <= 0.0
    {
        return Err(Error::Invalid(
            "liquid transfer reference volume is unrepresentable",
        ));
    }
    Ok(LiquidMassTransferVolumes {
        withdrawn_m3,
        received_m3,
        net_change_m3,
    })
}

struct MappedCenter {
    origin_m: [f64; 3],
    offset_m: [f64; 3],
    point_m: [f64; 3],
    roundoff_m: f64,
}

impl MappedCenter {
    fn coincides(&self, other: &Self) -> bool {
        let tolerance_m = self.roundoff_m + other.roundoff_m;
        (0..3).all(|d| {
            // Subtract origins BEFORE adding the small grid-local offsets.
            // Equal rounded world points alone can hide entire cell spacings
            // and silently lose the transported parcel's angular momentum.
            let separation = (Dd::from_f64(self.origin_m[d]) - Dd::from_f64(other.origin_m[d]))
                + (Dd::from_f64(self.offset_m[d]) - Dd::from_f64(other.offset_m[d]));
            separation.to_f64().abs() <= tolerance_m
        })
    }
}

fn center(
    liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>,
    coordinates: [usize; 3],
    origin_m: [f64; 3],
    axes_to_world: UnitQuaternion,
) -> Result<MappedCenter, Error> {
    let dx = liquid.flow().scaling().cell_size_m();
    let local = core::array::from_fn(|d| (coordinates[d] as f64 + 0.5) * dx);
    let offset_m = if axes_to_world == UnitQuaternion::IDENTITY {
        local
    } else {
        vector(axes_to_world.rotate_body_to_world(vec3(local)))
    };
    // Preserve the existing rounded world-space receipt. It is deliberately
    // separate from the relative geometry used to authorize the transfer.
    let point_m = if axes_to_world == UnitQuaternion::IDENTITY {
        core::array::from_fn(|d| (coordinates[d] as f64 + 0.5).mul_add(dx, origin_m[d]))
    } else {
        core::array::from_fn(|d| offset_m[d] + origin_m[d])
    };
    // A conservative allowance for coordinate scaling and the two quaternion
    // products: 64 binary64 epsilons times the grid-local L1 extent per owner.
    // The absolute world origin must NEVER enlarge this geometric allowance.
    let roundoff_m = local
        .into_iter()
        .map(|x| x.abs() * (64.0 * f64::EPSILON))
        .sum::<f64>();
    if origin_m
        .iter()
        .chain(&offset_m)
        .chain(&point_m)
        .any(|x| !x.is_finite())
        || !roundoff_m.is_finite()
    {
        return Err(Error::Invalid(
            "liquid transfer Cartesian origin or center is unrepresentable",
        ));
    }
    Ok(MappedCenter {
        origin_m,
        offset_m,
        point_m,
        roundoff_m,
    })
}

fn vec3(value: [f64; 3]) -> Vec3 {
    Vec3::new(value[0], value[1], value[2])
}

fn vector(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}

fn poll(gate: &CancelGate) -> Result<(), Error> {
    if gate.is_requested() {
        Err(FreeSurfaceError3::Cancelled.into())
    } else {
        Ok(())
    }
}

fn sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let mut result = 0.0;
    let mut correction = 0.0;
    for value in values {
        let next = result + value;
        correction += if result.abs() >= value.abs() {
            (result - next) + value
        } else {
            (value - next) + result
        };
        result = next;
    }
    result + correction
}
