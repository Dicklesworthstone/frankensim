//! Simultaneous Coulomb impacts through complete frozen modal contact frames.
//!
//! All normal and tangent rows enter ONE original fs-contact solve with the
//! full 3N-by-3N Gram mobility. The same rows transpose its impulses into native
//! momenta. Normal/shear and cross-contact response are never split into clips.
//! No fictitious inertia repairs locked tangents or a redundant contact set.
//!
//! Frames, friction coefficients and normal endpoint laws are explicit inputs.
//! This owner supplies no collision search, frame transport, friction history,
//! material inference, heat allocation, shock model or persistent time stepper.
//! Spatial momentum conservation is limited to symmetries in the original
//! retained basis; a constrained mode does not become a free rigid body.

use super::{
    CoupledModalSystem, ModalAcousticTimeModel, ModalCouplingError, audit_modal_momentum, checked,
    finite, invalid, poll, relative_column,
};
use fs_contact::frictional_multicontact::{
    FrictionalContact, FrictionalContactConfig, FrictionalContactError, FrictionalContactResiduals,
    FrictionalContactResult, NormalVelocityLaw, solve_frictional_contacts,
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;

/// One component's displacement/velocity traces in a common physical contact
/// frame, ordered [normal, tangent 1, tangent 2]. Each vector covers ALL of
/// that component's original modes. The supplied physical frame is orthonormal;
/// its projected modal rows need not be orthogonal or independent. All-zero
/// tangent traces explicitly retain a locked direction.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalContactFrameTrace {
    pub component: usize,
    pub shapes: [Vec<f64>; 3],
}

/// One identified contact with a single endpoint pair for all three axes.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalFrictionalImpactContact {
    pub id: u64,
    /// Receives +p in the supplied normal-first frame.
    pub first: ModalContactFrameTrace,
    /// Receives -p in that SAME frame.
    pub second: ModalContactFrameTrace,
    /// Actual separation is this reference plus the normal relative trace of q.
    pub reference_gap_m: f64,
    pub activation_tolerance_m: f64,
    /// Explicit coefficient for this selected event; no material law is inferred.
    pub friction_coefficient: f64,
    /// Original Newton restitution or an explicit nonnegative endpoint target.
    /// Either choice must pass the joint passive-event work admission.
    pub normal_law: NormalVelocityLaw,
}

/// Fixed frames, complete mobility and reusable native projection scratch.
/// Reprepare when the physical frame or modal basis changes, retaining the
/// continuing network state. The original bounded contact solver allocates its
/// candidate and receipt; this is not an allocation-free solver API.
pub struct PreparedModalFrictionalImpacts {
    columns: Vec<Vec<f64>>,
    delassus: Vec<f64>,
    layout: Vec<usize>,
    components: Vec<[usize; 2]>,
    reference_gaps: Vec<f64>,
    sources: Vec<FrictionalContact>,
    order: Vec<usize>,
    config: FrictionalContactConfig,
    before: Vec<[Dd; 3]>,
    requested: Vec<Dd>,
    impulses: Vec<f64>,
}

/// Original Coulomb/determinacy result and separate actual native-state audits.
/// A velocity-only result retains its unresolved local load allocation. Total
/// impact dissipation never becomes an automatic per-contact thermal credit.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalFrictionalImpactReceipt {
    pub sample: u64,
    pub contact_ids: Vec<u64>,
    pub contacts: FrictionalContactResult,
    pub gaps_m: Vec<f64>,
    pub relative_velocities_before_m_s: Vec<[f64; 3]>,
    pub relative_velocities_after_m_s: Vec<[f64; 3]>,
    pub kinetic_energy_change_j: f64,
    /// Original contact impulses paired with actual mean endpoint flows [J].
    pub contact_work_j: f64,
    /// Actual kinetic increment minus actual contact work [J].
    pub energy_residual_j: f64,
    /// Actual contact work minus the original solver's predicted work [J].
    pub contact_work_residual_j: f64,
    pub maximum_velocity_residual_m_s: f64,
    pub actual_residuals: FrictionalContactResiduals,
    pub maximum_relative_momentum_residual: f64,
}

#[derive(Debug)]
pub enum ModalFrictionalImpactError {
    Coupling(ModalCouplingError),
    Contact(FrictionalContactError),
    Frame {
        index: usize,
        axis: usize,
        source: ModalCouplingError,
    },
    Velocity {
        index: usize,
        axis: usize,
        residual_m_s: f64,
        tolerance_m_s: f64,
    },
}
impl From<ModalCouplingError> for ModalFrictionalImpactError {
    fn from(value: ModalCouplingError) -> Self {
        Self::Coupling(value)
    }
}
impl From<FrictionalContactError> for ModalFrictionalImpactError {
    fn from(value: FrictionalContactError) -> Self {
        Self::Contact(value)
    }
}
impl core::fmt::Display for ModalFrictionalImpactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Coupling(error) => core::fmt::Display::fmt(error, f),
            Self::Contact(error) => core::fmt::Display::fmt(error, f),
            Self::Frame {
                index,
                axis,
                source,
            } => write!(f, "modal contact {index}, axis {axis}: {source}"),
            Self::Velocity {
                index,
                axis,
                residual_m_s,
                tolerance_m_s,
            } => write!(
                f,
                "modal contact {index}, axis {axis}, velocity residual {residual_m_s:e} m/s exceeds {tolerance_m_s:e} m/s"
            ),
        }
    }
}
impl std::error::Error for ModalFrictionalImpactError {}

impl PreparedModalFrictionalImpacts {
    /// Prepare at most 21 normal-first triples within the original solver's
    /// limits. The network's existing max_setup_terms bounds (3N)^2*n_modes
    /// before assembly. Every normal must have positive finite retained
    /// mobility. Tangential nullspaces remain exactly present for native PSD
    /// admission/classification; no diagonal regularization changes the model.
    pub fn new(
        system: &CoupledModalSystem,
        contacts: Vec<ModalFrictionalImpactContact>,
        config: FrictionalContactConfig,
        gate: &CancelGate,
    ) -> Result<Self, ModalFrictionalImpactError> {
        poll(Some(gate))?;
        let count = contacts.len();
        if count > 21 || count > config.maximum_contacts {
            return Err(FrictionalContactError::Invalid(
                "modal contact count exceeds Coulomb capacity",
            )
            .into());
        }
        let dimension = 3 * count;
        let work = dimension
            .checked_mul(dimension)
            .and_then(|square| square.checked_mul(system.mode_count()))
            .ok_or_else(|| invalid("modal Coulomb setup extent overflow"))?;
        if work > system.config.max_setup_terms {
            return Err(invalid("modal Coulomb Gram work exceeds max_setup_terms").into());
        }
        let _ = solve_frictional_contacts(&[], &[], config, || gate.is_requested())?;
        let mut columns = Vec::with_capacity(dimension);
        let mut sources = Vec::with_capacity(count);
        let mut components = Vec::with_capacity(count);
        let mut reference_gaps = Vec::with_capacity(count);
        for (index, contact) in contacts.into_iter().enumerate() {
            poll(Some(gate))?;
            if !contact.reference_gap_m.is_finite()
                || !contact.activation_tolerance_m.is_finite()
                || contact.activation_tolerance_m < 0.0
                || !contact.friction_coefficient.is_finite()
                || contact.friction_coefficient < 0.0
                || !match contact.normal_law {
                    NormalVelocityLaw::NewtonRestitution { coefficient } => {
                        coefficient.is_finite() && (0.0..=1.0).contains(&coefficient)
                    }
                    NormalVelocityLaw::Target { speed } => speed.is_finite() && speed >= 0.0,
                }
            {
                return Err(FrictionalContactError::Invalid(
                    "modal frame geometry, friction or normal law",
                )
                .into());
            }
            for axis in 0..3 {
                columns.push(
                    relative_column(
                        system,
                        (contact.first.component, &contact.first.shapes[axis]),
                        (contact.second.component, &contact.second.shapes[axis]),
                        gate,
                    )
                    .map_err(|source| ModalFrictionalImpactError::Frame {
                        index,
                        axis,
                        source,
                    })?,
                );
            }
            let mut mobility = Dd::ZERO;
            for (mode, &shape) in columns[3 * index].iter().enumerate() {
                if mode % 64 == 0 {
                    poll(Some(gate))?;
                }
                mobility = mobility + Dd::from_f64(shape) * Dd::from_f64(shape);
            }
            if checked(mobility)? <= 0.0 {
                return Err(invalid(
                    "modal Coulomb normal requires positive retained inertia mobility",
                )
                .into());
            }
            sources.push(FrictionalContact {
                id: contact.id,
                gap: 0.0,
                activation_tolerance: contact.activation_tolerance_m,
                relative_velocity: [0.0; 3],
                friction_coefficient: contact.friction_coefficient,
                normal_law: contact.normal_law,
            });
            components.push([contact.first.component, contact.second.component]);
            reference_gaps.push(contact.reference_gap_m);
        }
        let mut order: Vec<_> = (0..count).collect();
        order.sort_unstable_by_key(|&index| sources[index].id);
        if order
            .windows(2)
            .any(|pair| sources[pair[0]].id == sources[pair[1]].id)
        {
            return Err(FrictionalContactError::Invalid("duplicate contact ID").into());
        }
        let mut delassus = vec![0.0; dimension * dimension];
        for i in 0..dimension {
            for j in 0..=i {
                poll(Some(gate))?;
                let mut entry = Dd::ZERO;
                for (mode, (&a, &b)) in columns[i].iter().zip(&columns[j]).enumerate() {
                    if mode % 64 == 0 {
                        poll(Some(gate))?;
                    }
                    entry = entry + Dd::from_f64(a) * Dd::from_f64(b);
                }
                let entry = checked(entry)?;
                delassus[i * dimension + j] = entry;
                delassus[j * dimension + i] = entry;
            }
        }
        poll(Some(gate))?;
        Ok(Self {
            columns,
            delassus,
            components,
            reference_gaps,
            sources,
            order,
            config,
            layout: system
                .models
                .iter()
                .map(|model| model.modes().len())
                .collect(),
            before: vec![[Dd::ZERO; 3]; count],
            requested: vec![Dd::ZERO; system.mode_count()],
            impulses: vec![0.0; system.mode_count()],
        })
    }

    pub fn contact_count(&self) -> usize {
        self.sources.len()
    }
    /// Complete row-major 3N-by-3N mobility, in caller/contact/axis order [1/kg].
    pub fn delassus_kg_inv(&self) -> &[f64] {
        &self.delassus
    }
    pub fn component_pairs(&self) -> impl ExactSizeIterator<Item = [usize; 2]> + '_ {
        self.components.iter().copied()
    }

    fn sample(
        &self,
        index: usize,
        models: &[ModalAcousticTimeModel],
        gate: &CancelGate,
    ) -> Result<(f64, [Dd; 3]), ModalCouplingError> {
        let mut gap = Dd::from_f64(self.reference_gaps[index]);
        let mut velocity = [Dd::ZERO; 3];
        for (mode, state) in models.iter().flat_map(|model| model.states()).enumerate() {
            if mode % 64 == 0 {
                poll(Some(gate))?;
            }
            gap = gap
                + Dd::from_f64(self.columns[3 * index][mode])
                    * Dd::from_f64(state.displacement_m_sqrt_kg);
            for axis in 0..3 {
                velocity[axis] = velocity[axis]
                    + Dd::from_f64(self.columns[3 * index + axis][mode])
                        * Dd::from_f64(state.velocity_m_sqrt_kg_per_s);
            }
        }
        for value in velocity {
            checked(value)?;
        }
        Ok((checked(gap)?, velocity))
    }
}

impl CoupledModalSystem {
    /// Apply one simultaneous original Signorini--Coulomb event. All three
    /// axes, both finite endpoints and all contacts share the same solve and
    /// native reset. Actual final-velocity Coulomb residuals, global work and
    /// per-coordinate momentum pass independently before state publication.
    ///
    /// Frozen q, material/radiation histories, held sources, preceding frame
    /// and clock remain unchanged. Cancellation, ambiguous/unresolved native
    /// determinacy, passivity or any original state limit publishes nothing.
    pub fn apply_frictional_impacts(
        &mut self,
        impacts: &mut PreparedModalFrictionalImpacts,
        gate: &CancelGate,
    ) -> Result<ModalFrictionalImpactReceipt, ModalFrictionalImpactError> {
        self.apply_frictional_impacts_with_admission(impacts, gate, |_| Ok(()))
    }

    pub(crate) fn apply_frictional_impacts_with_admission(
        &mut self,
        impacts: &mut PreparedModalFrictionalImpacts,
        gate: &CancelGate,
        admit: impl FnOnce(&[ModalAcousticTimeModel]) -> Result<(), ModalCouplingError>,
    ) -> Result<ModalFrictionalImpactReceipt, ModalFrictionalImpactError> {
        poll(Some(gate))?;
        self.require_time_domain()?;
        if impacts.impulses.len() != self.mode_count()
            || self.models.len() != impacts.layout.len()
            || self
                .models
                .iter()
                .zip(&impacts.layout)
                .any(|(model, &count)| model.modes().len() != count)
        {
            return Err(
                invalid("prepared Coulomb frames no longer match the native modal layout").into(),
            );
        }
        for &index in &impacts.order {
            let (gap, before) = impacts.sample(index, &self.models, gate)?;
            impacts.sources[index].gap = gap;
            for axis in 0..3 {
                impacts.sources[index].relative_velocity[axis] = checked(before[axis])?;
            }
            impacts.before[index] = before;
        }
        let contacts =
            solve_frictional_contacts(&impacts.sources, &impacts.delassus, impacts.config, || {
                gate.is_requested()
            })?;
        for mode in 0..self.mode_count() {
            if mode % 64 == 0 {
                poll(Some(gate))?;
            }
            let mut requested = Dd::ZERO;
            // Match the original law's stable-ID order, including cancellation
            // tails, so contact enumeration does not change native recoil.
            for &index in &impacts.order {
                for axis in 0..3 {
                    let shape = impacts.columns[3 * index + axis][mode];
                    let impulse = contacts.impulses[index][axis];
                    let term = Dd::from_f64(shape) * Dd::from_f64(impulse);
                    if shape != 0.0 && impulse != 0.0 && term.hi == 0.0 && term.lo == 0.0 {
                        return Err(
                            invalid("Coulomb projection loses nonzero modal momentum").into()
                        );
                    }
                    requested = requested + term;
                }
            }
            impacts.requested[mode] = requested;
            impacts.impulses[mode] = checked(requested)?;
        }
        let changes_state = contacts
            .impulses
            .iter()
            .flatten()
            .any(|&impulse| impulse != 0.0);
        let jump = if changes_state {
            Some(self.stage_modal_impulse(&impacts.impulses, Some(gate))?)
        } else {
            None
        };
        let after_models = if changes_state {
            &self.candidates
        } else {
            &self.models
        };
        let maximum_relative_momentum_residual =
            audit_modal_momentum(self, after_models, &impacts.requested, gate)?
                .max(jump.map_or(0.0, |jump| jump.maximum_relative_momentum_residual));
        let mut velocities_after = vec![[0.0; 3]; impacts.sources.len()];
        let mut actual_residuals = FrictionalContactResiduals::default();
        let mut maximum_velocity_residual_m_s = 0.0_f64;
        let mut work = Dd::ZERO;
        for &index in &impacts.order {
            let (_, after) = impacts.sample(index, after_models, gate)?;
            for axis in 0..3 {
                let actual = checked(after[axis])?;
                velocities_after[index][axis] = actual;
                let predicted = contacts.relative_velocities_after[index][axis];
                let residual_m_s = checked(after[axis] - Dd::from_f64(predicted))?;
                let tolerance_m_s = finite(
                    impacts.config.velocity_tolerance
                        + self.config.solve_relative_tolerance
                            * actual
                                .abs()
                                .max(predicted.abs())
                                .max(impacts.sources[index].relative_velocity[axis].abs()),
                )?;
                if residual_m_s.abs() > tolerance_m_s {
                    return Err(ModalFrictionalImpactError::Velocity {
                        index,
                        axis,
                        residual_m_s,
                        tolerance_m_s,
                    });
                }
                maximum_velocity_residual_m_s =
                    maximum_velocity_residual_m_s.max(residual_m_s.abs());
                let mean = impacts.before[index][axis] * Dd::from_f64(0.5)
                    + after[axis] * Dd::from_f64(0.5);
                work = work + Dd::from_f64(contacts.impulses[index][axis]) * mean;
            }
            if impacts.sources[index].gap <= impacts.sources[index].activation_tolerance {
                audit_contact(
                    impacts.sources[index].friction_coefficient,
                    contacts.normal_targets[index],
                    contacts.impulses[index],
                    after,
                    impacts.config,
                    &mut actual_residuals,
                )?;
            }
        }
        let r = actual_residuals;
        if r.normal_impulse_violation > impacts.config.impulse_tolerance
            || r.normal_velocity_violation > impacts.config.velocity_tolerance
            || r.normal_complementarity > impacts.config.complementarity_tolerance
            || r.cone_violation > impacts.config.impulse_tolerance
            || r.sliding_impulse_residual > impacts.config.impulse_tolerance
            || r.tangential_duality_gap > impacts.config.complementarity_tolerance
        {
            return Err(FrictionalContactError::ResidualNotSatisfied(r).into());
        }
        let contact_work_j = checked(work)?;
        if contact_work_j > impacts.config.energy_tolerance {
            return Err(FrictionalContactError::EnergyIncrease {
                work: contact_work_j,
                allowed: impacts.config.energy_tolerance,
            }
            .into());
        }
        let kinetic_energy_change_j = jump.map_or(0.0, |jump| jump.kinetic_energy_change_j);
        let energy_residual_j = checked(Dd::from_f64(kinetic_energy_change_j) - work)?;
        let contact_work_residual_j = checked(work - Dd::from_f64(contacts.relative_work))?;
        for (residual_j, scale) in [
            (
                energy_residual_j,
                kinetic_energy_change_j.abs().max(contact_work_j.abs()),
            ),
            (
                contact_work_residual_j,
                contact_work_j.abs().max(contacts.relative_work.abs()),
            ),
        ] {
            let tolerance_j = finite(
                self.config.energy_absolute_tolerance_j
                    + self.config.energy_relative_tolerance * scale,
            )?;
            if residual_j.abs() > tolerance_j {
                return Err(ModalCouplingError::EnergyBalance {
                    residual_j,
                    tolerance_j,
                }
                .into());
            }
        }
        let receipt = ModalFrictionalImpactReceipt {
            sample: self.samples_rendered(),
            contact_ids: impacts.sources.iter().map(|source| source.id).collect(),
            contacts,
            gaps_m: impacts.sources.iter().map(|source| source.gap).collect(),
            relative_velocities_before_m_s: impacts
                .sources
                .iter()
                .map(|source| source.relative_velocity)
                .collect(),
            relative_velocities_after_m_s: velocities_after,
            kinetic_energy_change_j,
            contact_work_j,
            energy_residual_j,
            contact_work_residual_j,
            maximum_velocity_residual_m_s,
            actual_residuals,
            maximum_relative_momentum_residual,
        };
        admit(after_models)?;
        poll(Some(gate))?;
        if changes_state {
            self.publish_modal_impulse();
        }
        Ok(receipt)
    }
}

// Recheck the original physical Coulomb law at the ACTUAL staged velocity,
// independently of the native solve's natural map and predicted W*p response.
fn audit_contact(
    mu: f64,
    target: f64,
    p: [f64; 3],
    velocity: [Dd; 3],
    config: FrictionalContactConfig,
    r: &mut FrictionalContactResiduals,
) -> Result<(), ModalCouplingError> {
    let g = velocity[0] - Dd::from_f64(target);
    let radius = Dd::from_f64(mu) * Dd::from_f64(p[0]);
    let v = [
        checked(velocity[0])?,
        checked(velocity[1])?,
        checked(velocity[2])?,
    ];
    let speed = finite(v[1].hypot(v[2]))?;
    r.normal_impulse_violation = r.normal_impulse_violation.max((-p[0]).max(0.0));
    r.normal_velocity_violation = r.normal_velocity_violation.max((-checked(g)?).max(0.0));
    r.normal_complementarity = r
        .normal_complementarity
        .max(checked(Dd::from_f64(p[0]) * g)?.abs());
    r.cone_violation = r
        .cone_violation
        .max(checked(Dd::from_f64(p[1].hypot(p[2])) - radius)?.max(0.0));
    let duality_gap = radius * Dd::from_f64(speed)
        + Dd::from_f64(p[1]) * velocity[1]
        + Dd::from_f64(p[2]) * velocity[2];
    r.tangential_duality_gap = r.tangential_duality_gap.max(checked(duality_gap)?.abs());
    if speed > config.velocity_tolerance {
        let a = checked(Dd::from_f64(p[1]) + radius * Dd::from_f64(v[1] / speed))?;
        let b = checked(Dd::from_f64(p[2]) + radius * Dd::from_f64(v[2] / speed))?;
        r.sliding_impulse_residual = r.sliding_impulse_residual.max(finite(a.hypot(b))?);
    }
    Ok(())
}
