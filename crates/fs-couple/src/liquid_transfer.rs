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
//! Both owners retain the same original material/temperature/caloric objects
//! and accepted physical time. Ordinary native steps advance their clocks.
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
    if !donor.shares_mass_transfer_source(receiver) {
        return Err(Error::Invalid(
            "liquid transfer must retain the same original surface, temperature and caloric sources",
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
    let point_m = center(
        donor,
        transfer.donor,
        transfer.donor_origin_m,
        transfer.donor_axes_to_world,
    )?;
    if point_m
        != center(
            receiver,
            transfer.receiver,
            transfer.receiver_origin_m,
            transfer.receiver_axes_to_world,
        )?
    {
        return Err(Error::Invalid(
            "liquid transfer requires coincident mapped cell centers in one Cartesian frame",
        ));
    }
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
    })
}

fn center(
    liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>,
    coordinates: [usize; 3],
    origin_m: [f64; 3],
    axes_to_world: UnitQuaternion,
) -> Result<[f64; 3], Error> {
    let dx = liquid.flow().scaling().cell_size_m();
    let point = if axes_to_world == UnitQuaternion::IDENTITY {
        core::array::from_fn(|d| (coordinates[d] as f64 + 0.5).mul_add(dx, origin_m[d]))
    } else {
        let local = core::array::from_fn(|d| (coordinates[d] as f64 + 0.5) * dx);
        let world = vector(axes_to_world.rotate_body_to_world(vec3(local)));
        core::array::from_fn(|d| world[d] + origin_m[d])
    };
    if origin_m.iter().chain(&point).any(|x| !x.is_finite()) {
        return Err(Error::Invalid(
            "liquid transfer Cartesian origin or center is unrepresentable",
        ));
    }
    Ok(point)
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

fn sum(values: [f64; 4]) -> f64 {
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
