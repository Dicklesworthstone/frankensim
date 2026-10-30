//! Original material sources drive the native restrained-member/bond solve.
//! Bulk temperatures may evolve along each source's admitted temperature line.
//! Interface conditions remain explicit and fixed: changing them requires a
//! thermodynamic cohesive-history transfer, not replacement of coefficients.

use std::collections::BTreeMap;

use fs_exec::Cx;
use fs_matdb::QueryPoint;
use fs_material::cohesive::source::ResolvedBilinearCohesive;
use fs_material::state_point::MaterialStatePointError;
use fs_material::thermovisco::source::ResolvedThermoMaxwell;

use super::{Budget, Config, Error as JointError, Joint, Member, Report};

/// One original member geometry and its complete source-backed thermal law.
#[derive(Clone, Copy, Debug)]
pub struct SourceMember<'a> {
    pub id: u64,
    pub reference_length_m: f64,
    pub reference_area_m2: f64,
    pub material: &'a ResolvedThermoMaxwell,
    pub birth_temperature_k: f64,
}

/// Source admission and native numerical failures retain their original causes.
#[derive(Debug)]
pub enum SourceJointError {
    Invalid(&'static str),
    Interface(MaterialStatePointError),
    Joint(JointError),
    Cancelled,
}
impl From<JointError> for SourceJointError {
    fn from(error: JointError) -> Self {
        match error {
            JointError::Cancelled => Self::Cancelled,
            other => Self::Joint(other),
        }
    }
}
impl core::fmt::Display for SourceJointError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-bound thermal joint refused: {self:?}")
    }
}
impl std::error::Error for SourceJointError {}

/// Native history with the exact bulk and ordered-interface source owners.
/// Cloning checkpoints the mechanical state; immutable sources stay borrowed.
#[derive(Clone, Debug)]
pub struct SourceJoint<'a> {
    joint: Joint<'a>,
    bulk_sources: BTreeMap<u64, &'a ResolvedThermoMaxwell>,
    interface: &'a ResolvedBilinearCohesive,
    surface_members: [u64; 2],
    interface_query: QueryPoint,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), SourceJointError> {
    if cancelled() {
        Err(SourceJointError::Cancelled)
    } else {
        Ok(())
    }
}

impl<'a> SourceJoint<'a> {
    /// Build the existing series-joint owner from source laws without copying
    /// or selecting coefficients. `surface_members` names the actual receiving
    /// member on each ordered interface side; both IDs must exist. Geometry and
    /// physical bonding are declared inputs, never inferred from material names.
    ///
    /// The cohesive source's complete query is the fixed external interface
    /// state. It may differ from the bulk temperatures. This constructor does
    /// not impose thermal equilibrium or add interface heat transport.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        time_s: f64,
        members: Vec<SourceMember<'a>>,
        interface_id: u64,
        interface_area_m2: f64,
        interface: &'a ResolvedBilinearCohesive,
        surface_members: [u64; 2],
        budget: Budget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, SourceJointError> {
        poll(&mut cancelled)?;
        if members.is_empty() || members.len() > budget.max_members {
            return Err(JointError::Budget.into());
        }
        let mut bulk_sources = BTreeMap::new();
        let mut native = Vec::with_capacity(members.len());
        for member in members {
            poll(&mut cancelled)?;
            if bulk_sources.insert(member.id, member.material).is_some() {
                return Err(SourceJointError::Invalid("duplicate source member ID"));
            }
            // Birth and every later native update enforce this exact source's
            // temperature, strain and rate envelope; no endpoint extrapolation.
            native.push(Member::new(
                member.id,
                member.reference_length_m,
                member.reference_area_m2,
                member.material.law(),
                member.birth_temperature_k,
            )?);
        }
        let receiving = surface_members.map(|id| {
            bulk_sources
                .get(&id)
                .and_then(|source| source.samples().first())
                .map(|sample| sample.material())
                .ok_or(SourceJointError::Invalid(
                    "ordered interface side must name an original source member",
                ))
        });
        let [first, second] = receiving;
        interface
            .require_surface_materials(first?, second?)
            .map_err(SourceJointError::Interface)?;
        let joint = Joint::new(
            time_s,
            native,
            interface_id,
            interface_area_m2,
            interface.law(),
            budget,
        )?;
        let candidate = Self {
            joint,
            bulk_sources,
            interface,
            surface_members,
            interface_query: interface.resolved().query(),
        };
        poll(&mut cancelled)?;
        Ok(candidate)
    }

    /// Original native displacement, member memory, damage and stored energy.
    pub const fn joint(&self) -> &Joint<'a> {
        &self.joint
    }
    /// Original card identities, pinned claims and full temperature-line sources.
    pub fn member_sources(&self) -> &BTreeMap<u64, &'a ResolvedThermoMaxwell> {
        &self.bulk_sources
    }
    pub const fn interface_source(&self) -> &'a ResolvedBilinearCohesive {
        self.interface
    }
    pub const fn surface_members(&self) -> [u64; 2] {
        self.surface_members
    }
    pub const fn interface_query(&self) -> &QueryPoint {
        &self.interface_query
    }

    /// Advance the original joint, retaining every branch strain and cohesive
    /// damage variable. Bulk temperatures vary within their original sources;
    /// all other bulk coordinates stay fixed. The caller supplies the actual
    /// interface query explicitly: changing any value or quantity convention
    /// refuses before evolution. Bulk heating therefore cannot silently stand
    /// in for an unsupported temperature-dependent cohesive-history update.
    ///
    /// The native force, compatibility, storage/work and cancellation checks
    /// remain authoritative. Refusal publishes neither material history.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_controlled(
        &self,
        target_s: f64,
        extension_m: f64,
        temperatures_k: &BTreeMap<u64, f64>,
        interface_query: &QueryPoint,
        config: Config,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, Report), SourceJointError> {
        poll(&mut cancelled)?;
        if interface_query != &self.interface_query {
            return Err(SourceJointError::Invalid(
                "changing cohesive external state requires a thermodynamic history transfer",
            ));
        }
        let (joint, report) = self.joint.advance_controlled(
            target_s,
            extension_m,
            temperatures_k,
            config,
            &mut cancelled,
        )?;
        let mut candidate = self.clone();
        candidate.joint = joint;
        poll(&mut cancelled)?;
        Ok((candidate, report))
    }

    /// Same operation in the existing executor cancellation scope.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to(
        &self,
        cx: &Cx<'_>,
        target_s: f64,
        extension_m: f64,
        temperatures_k: &BTreeMap<u64, f64>,
        interface_query: &QueryPoint,
        config: Config,
    ) -> Result<(Self, Report), SourceJointError> {
        self.advance_controlled(
            target_s,
            extension_m,
            temperatures_k,
            interface_query,
            config,
            || cx.checkpoint().is_err(),
        )
    }
}
