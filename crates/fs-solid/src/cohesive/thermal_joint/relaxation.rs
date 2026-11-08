//! Intact relaxing adhesive layers between the original thermal series members.
//!
//! The existing bulk owner supplies each endpoint affine stress and retained
//! Maxwell memory; the existing adhesive owner supplies its exact straight-ramp
//! Prony response. Common force and total extension determine one scalar affine
//! compatibility equation. No substitute elastic bond or fitted damping enters.
//!
//! Bulk temperatures may change within their original source domains. The
//! interface's complete source query remains fixed, including its temperature.
//! This is a massless, small-strain axial joint with endpoint load/extension
//! control; refine paths for transient accuracy. There is no damage, healing,
//! thermal transport, automatic heat deposition or finite-body peel model.

use std::collections::BTreeMap;

use fs_exec::Cx;
use fs_matdb::QueryPoint;
use fs_material::cohesive::{
    CohesiveError,
    relaxation::{ViscoelasticBondPoint, source::ResolvedViscoelasticBond},
};
use fs_material::thermovisco::{Energy, source::ResolvedThermoMaxwell};
use fs_math::dd::Dd;

use super::source::SourceMember;
use super::{Error as JointError, Member, accumulate, finite, poll, positive};

/// Combined original bulk and adhesive memory admission.
#[derive(Clone, Copy, Debug)]
pub struct Budget {
    pub maximum_members: usize,
    /// Includes normal and both shear branch banks of the original bond point.
    pub maximum_branch_states: usize,
}

/// Independent mechanical and constitutive arithmetic gates.
#[derive(Clone, Copy, Debug)]
pub struct Config {
    pub force_tolerance_n: f64,
    pub extension_tolerance_m: f64,
    pub local_bond_energy_tolerance_j_m2: f64,
    pub energy_tolerance_j: f64,
}

/// Absolute endpoint control; neither variant changes the accepted old state.
#[derive(Clone, Copy, Debug)]
pub enum Control {
    ExtensionM(f64),
    ForceN(f64),
}

#[derive(Debug)]
pub enum Error {
    Joint(JointError),
    Bond {
        interface: u64,
        error: CohesiveError,
    },
    ControlForce {
        residual_n: f64,
    },
    Source(&'static str),
}
impl From<JointError> for Error {
    fn from(error: JointError) -> Self {
        Self::Joint(error)
    }
}
impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-bound relaxing thermal joint: {self:?}")
    }
}
impl std::error::Error for Error {}

/// Immutable accepted geometry, original sources, member and bond memories.
#[derive(Clone, Debug)]
pub struct SourceRelaxingJoint<'a> {
    time_s: f64,
    extension_m: f64,
    interface_id: u64,
    area_m2: f64,
    members: BTreeMap<u64, Member<'a>>,
    bulk_sources: BTreeMap<u64, &'a ResolvedThermoMaxwell>,
    interface: &'a ResolvedViscoelasticBond,
    surface_members: [u64; 2],
    interface_query: QueryPoint,
    bond: ViscoelasticBondPoint<'a>,
}

/// Actual extensive increments; physical adhesive heat is distinct from both
/// bulk numerical loss and the signed endpoint-work quadrature correction.
#[derive(Clone, Debug)]
pub struct Report {
    pub force_n: f64,
    pub opening_m: f64,
    pub extension_m: f64,
    pub extension_residual_m: f64,
    pub maximum_force_residual_n: f64,
    pub bulk_energy_j: Energy,
    pub bond_work_j: f64,
    pub bond_stored_change_j: f64,
    pub bond_relaxation_heat_j: f64,
    pub bond_energy_residual_j: f64,
    /// F_end*(D_end-D_old); a time-discrete quadrature, not exact actuator work.
    pub endpoint_work_j: f64,
    /// F_end*delta_gap minus exact Prony ramp work. Signed and never heat.
    pub bond_work_correction_j: f64,
    /// Common endpoint force times member extension changes minus original
    /// bulk mechanical work. Includes exact-ramp quadrature, never heat.
    pub bulk_work_correction_j: f64,
    pub energy_residual_j: f64,
}

impl<'a> SourceRelaxingJoint<'a> {
    /// Initialize a stress-free assembly with original source-backed members.
    /// Ordered interface sides name the actual two receiving member IDs; their
    /// material identities must match the original interface card exactly.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        time_s: f64,
        members: Vec<SourceMember<'a>>,
        interface_id: u64,
        area_m2: f64,
        interface: &'a ResolvedViscoelasticBond,
        surface_members: [u64; 2],
        budget: Budget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        if !time_s.is_finite() || time_s < 0.0 {
            return Err(JointError::Invalid("invalid start time").into());
        }
        positive(area_m2)?;
        if members.is_empty() || members.len() > budget.maximum_members {
            return Err(JointError::Budget.into());
        }
        if surface_members[0] == surface_members[1] {
            return Err(Error::Source(
                "two distinct original receiving members required",
            ));
        }
        let law = interface.law();
        let mut branches = law
            .normal_spectrum()
            .terms
            .len()
            .checked_add(
                law.shear_spectrum()
                    .terms
                    .len()
                    .checked_mul(2)
                    .ok_or(JointError::Budget)?,
            )
            .ok_or(JointError::Budget)?;
        if branches > budget.maximum_branch_states {
            return Err(JointError::Budget.into());
        }
        let mut native = BTreeMap::new();
        let mut bulk_sources = BTreeMap::new();
        for member in members {
            poll(&mut cancelled)?;
            branches = branches
                .checked_add(member.material.law().branch_count())
                .ok_or(JointError::Budget)?;
            if branches > budget.maximum_branch_states {
                return Err(JointError::Budget.into());
            }
            if bulk_sources.insert(member.id, member.material).is_some() {
                return Err(Error::Source("duplicate original member ID"));
            }
            native.insert(
                member.id,
                Member::new(
                    member.id,
                    member.reference_length_m,
                    member.reference_area_m2,
                    member.material.law(),
                    member.birth_temperature_k,
                )?,
            );
        }
        for (side, id) in surface_members.into_iter().enumerate() {
            let material = bulk_sources
                .get(&id)
                .and_then(|source| source.samples().first())
                .ok_or(Error::Source(
                    "interface side must name an original source member",
                ))?
                .material();
            if interface.surfaces()[side].material != *material {
                return Err(Error::Source(
                    "ordered interface material differs from receiving member",
                ));
            }
        }
        poll(&mut cancelled)?;
        let candidate = Self {
            time_s,
            extension_m: 0.0,
            interface_id,
            area_m2,
            members: native,
            bulk_sources,
            interface,
            surface_members,
            interface_query: interface.resolved().query(),
            bond: law.bonded_point(),
        };
        poll(&mut cancelled)?;
        Ok(candidate)
    }

    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn extension_m(&self) -> f64 {
        self.extension_m
    }
    pub const fn interface_id(&self) -> u64 {
        self.interface_id
    }
    pub const fn original_bonded_area_m2(&self) -> f64 {
        self.area_m2
    }
    pub fn members(&self) -> &BTreeMap<u64, Member<'a>> {
        &self.members
    }
    pub fn member_sources(&self) -> &BTreeMap<u64, &'a ResolvedThermoMaxwell> {
        &self.bulk_sources
    }
    pub const fn interface_source(&self) -> &'a ResolvedViscoelasticBond {
        self.interface
    }
    pub const fn surface_members(&self) -> [u64; 2] {
        self.surface_members
    }
    pub const fn interface_query(&self) -> &QueryPoint {
        &self.interface_query
    }
    pub const fn bond(&self) -> &ViscoelasticBondPoint<'a> {
        &self.bond
    }
    pub fn stored_energy_j(&self) -> f64 {
        let mut total = Dd::from_f64(self.area_m2 * self.bond.stored_energy_j_m2());
        for member in self.members.values() {
            total = total
                + Dd::from_f64(member.length_m * member.area_m2)
                    * Dd::from_f64(member.point.stored_energy_j_m3());
        }
        total.to_f64()
    }

    /// Advance all original memories together at the requested absolute time.
    /// Bulk temperature is evaluated by its original implicit material owner;
    /// the adhesive integrates its straight jump ramp exactly. Its endpoint
    /// stiffness enters compatibility, while its work-mean traction owns work.
    /// The complete interface query must remain fixed. Refusal, including late
    /// cancellation, publishes no clock, displacement, source or branch update.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_controlled(
        &self,
        target_s: f64,
        control: Control,
        temperatures_k: &BTreeMap<u64, f64>,
        interface_query: &QueryPoint,
        config: Config,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, Report), Error> {
        poll(&mut cancelled)?;
        if interface_query != &self.interface_query {
            return Err(Error::Source(
                "interface conditions remain fixed; memory transfer is not inferred",
            ));
        }
        let dt = finite(target_s - self.time_s)?;
        positive(dt)?;
        for value in [
            config.force_tolerance_n,
            config.extension_tolerance_m,
            config.local_bond_energy_tolerance_j_m2,
            config.energy_tolerance_j,
        ] {
            positive(value)?;
        }
        if temperatures_k.len() != self.members.len() {
            return Err(JointError::Invalid("one temperature per original member required").into());
        }
        let mut affine = Vec::with_capacity(self.members.len());
        let (mut compliance, mut free) = (Dd::ZERO, Dd::ZERO);
        for (&id, member) in &self.members {
            poll(&mut cancelled)?;
            let temperature = *temperatures_k
                .get(&id)
                .ok_or(JointError::Invalid("missing member temperature"))?;
            let law = member
                .point
                .linearize(temperature, dt)
                .map_err(|error| JointError::Material { member: id, error })?;
            let s = finite(member.length_m / member.area_m2 / law.tangent_pa)?;
            positive(s)?;
            compliance = compliance + Dd::from_f64(s);
            free =
                free - Dd::from_f64(member.length_m) * Dd::from_f64(law.offset_pa / law.tangent_pa);
            affine.push((id, temperature, law));
        }
        let s = finite(compliance.to_f64())?;
        positive(s)?;
        finite(free.to_f64())?;
        let bond_error = |error| Error::Bond {
            interface: self.interface_id,
            error,
        };
        let prepared = self.interface.law().prepare(dt).map_err(bond_error)?;
        let old_gap = self.bond.jump_m()[0];
        let held = self
            .bond
            .trial(
                &prepared,
                [old_gap, 0.0, 0.0],
                config.local_bond_energy_tolerance_j_m2,
            )
            .map_err(bond_error)?;
        poll(&mut cancelled)?;
        let held_force = finite(self.area_m2 * held.traction_pa[0])?;
        let stiffness = finite(self.area_m2 * prepared.endpoint_stiffness_pa_m()[0])?;
        positive(stiffness)?;
        let delta = match control {
            Control::ExtensionM(extension) => {
                finite(extension)?;
                let numerator = Dd::from_f64(extension)
                    - free
                    - Dd::from_f64(old_gap)
                    - compliance * Dd::from_f64(held_force);
                let denominator = Dd::ONE + compliance * Dd::from_f64(stiffness);
                finite((numerator / denominator).to_f64())?
            }
            Control::ForceN(force) => {
                finite(force)?;
                finite(
                    ((Dd::from_f64(force) - Dd::from_f64(held_force)) / Dd::from_f64(stiffness))
                        .to_f64(),
                )?
            }
        };
        let gap = finite(old_gap + delta)?;
        let bond = self
            .bond
            .trial(
                &prepared,
                [gap, 0.0, 0.0],
                config.local_bond_energy_tolerance_j_m2,
            )
            .map_err(bond_error)?;
        let force = finite(self.area_m2 * bond.traction_pa[0])?;
        let mut maximum_force_residual_n = if let Control::ForceN(target) = control {
            finite(force - target)?.abs()
        } else {
            0.0
        };
        if maximum_force_residual_n > config.force_tolerance_n {
            return Err(Error::ControlForce {
                residual_n: maximum_force_residual_n,
            });
        }
        let mut members = BTreeMap::new();
        let mut energy = Energy::default();
        let mut extension = Dd::from_f64(gap);
        let mut bulk_correction = Dd::ZERO;
        for (id, temperature, law) in affine {
            poll(&mut cancelled)?;
            let member = &self.members[&id];
            let strain = finite(
                ((Dd::from_f64(force / member.area_m2) - Dd::from_f64(law.offset_pa))
                    / Dd::from_f64(law.tangent_pa))
                .to_f64(),
            )?;
            let volume = member.length_m * member.area_m2;
            let local_gate = config.energy_tolerance_j / self.members.len() as f64 / volume;
            let step = member
                .point
                .advance(temperature, strain, dt, local_gate)
                .map_err(|error| JointError::Material { member: id, error })?;
            let residual = finite(member.area_m2 * step.stress_pa - force)?;
            maximum_force_residual_n = maximum_force_residual_n.max(residual.abs());
            if residual.abs() > config.force_tolerance_n {
                return Err(JointError::MemberForce {
                    member: id,
                    residual_n: residual,
                }
                .into());
            }
            extension = extension + Dd::from_f64(member.length_m) * Dd::from_f64(strain);
            bulk_correction = bulk_correction
                + Dd::from_f64(force)
                    * Dd::from_f64(member.length_m)
                    * (Dd::from_f64(strain) - Dd::from_f64(member.point.total_strain()))
                - Dd::from_f64(volume) * Dd::from_f64(step.energy.mechanical_work);
            accumulate(&mut energy, step.energy, volume)?;
            members.insert(
                id,
                Member {
                    id,
                    length_m: member.length_m,
                    area_m2: member.area_m2,
                    point: step.point,
                },
            );
        }
        let actual_extension = finite(extension.to_f64())?;
        let extension_m = match control {
            Control::ExtensionM(target) => target,
            Control::ForceN(_) => actual_extension,
        };
        let extension_residual_m = finite((extension - Dd::from_f64(extension_m)).to_f64())?;
        if extension_residual_m.abs() > config.extension_tolerance_m {
            return Err(JointError::Compatibility(extension_residual_m).into());
        }
        let work = finite(self.area_m2 * bond.work_j_m2)?;
        let storage = finite(self.area_m2 * bond.stored_change_j_m2)?;
        let heat = finite(self.area_m2 * bond.viscous_dissipated_j_m2)?;
        let bond_residual = finite(self.area_m2 * bond.defect_j_m2)?;
        let endpoint_work = finite(
            (Dd::from_f64(force) * (Dd::from_f64(extension_m) - Dd::from_f64(self.extension_m)))
                .to_f64(),
        )?;
        let correction = finite(
            (Dd::from_f64(force) * (Dd::from_f64(gap) - Dd::from_f64(old_gap))
                - Dd::from_f64(work))
            .to_f64(),
        )?;
        let bulk_correction = finite(bulk_correction.to_f64())?;
        let residual = finite(
            (Dd::from_f64(energy.stored_change)
                + Dd::from_f64(storage)
                + Dd::from_f64(energy.dissipated)
                + Dd::from_f64(energy.numerical_dissipation)
                + Dd::from_f64(heat)
                + Dd::from_f64(correction)
                + Dd::from_f64(bulk_correction)
                - Dd::from_f64(endpoint_work)
                - Dd::from_f64(energy.thermal_strain_work)
                - Dd::from_f64(energy.modulus_work))
            .to_f64(),
        )?;
        let conjugacy = finite(
            (Dd::from_f64(endpoint_work)
                - Dd::from_f64(correction)
                - Dd::from_f64(bulk_correction)
                - Dd::from_f64(energy.mechanical_work)
                - Dd::from_f64(work))
            .to_f64(),
        )?;
        for defect in [residual, conjugacy, bond_residual] {
            if defect.abs() > config.energy_tolerance_j {
                return Err(JointError::Energy(defect).into());
            }
        }
        poll(&mut cancelled)?;
        let next = Self {
            time_s: target_s,
            extension_m,
            interface_id: self.interface_id,
            area_m2: self.area_m2,
            members,
            bulk_sources: self.bulk_sources.clone(),
            interface: self.interface,
            surface_members: self.surface_members,
            interface_query: self.interface_query.clone(),
            bond: bond.point,
        };
        finite(next.stored_energy_j())?;
        poll(&mut cancelled)?;
        Ok((
            next,
            Report {
                force_n: force,
                opening_m: gap,
                extension_m,
                extension_residual_m,
                maximum_force_residual_n,
                bulk_energy_j: energy,
                bond_work_j: work,
                bond_stored_change_j: storage,
                bond_relaxation_heat_j: heat,
                bond_energy_residual_j: bond_residual,
                endpoint_work_j: endpoint_work,
                bond_work_correction_j: correction,
                bulk_work_correction_j: bulk_correction,
                energy_residual_j: residual,
            },
        ))
    }

    /// Same atomic operation within the existing executor cancellation scope.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to(
        &self,
        cx: &Cx<'_>,
        target_s: f64,
        control: Control,
        temperatures_k: &BTreeMap<u64, f64>,
        interface_query: &QueryPoint,
        config: Config,
    ) -> Result<(Self, Report), Error> {
        self.advance_controlled(
            target_s,
            control,
            temperatures_k,
            interface_query,
            config,
            || cx.checkpoint().is_err(),
        )
    }
}
