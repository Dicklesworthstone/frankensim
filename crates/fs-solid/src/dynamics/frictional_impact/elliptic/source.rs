//! Original material-friction sources sampled at native incoming contact speed.
//!
//! Source candidates and the native impulse share one acceptance boundary.
//! Instantaneous sampling updates original claim receipts and coefficients;
//! it does not advance slip, travel, temperature history or a thermal store.

use super::*;
use fs_contact::interface_binding::thermal_friction::elliptic::impulse::{
    SourceFrictionConditions, SourceImpulseCoefficientRole, SourceImpulseError,
};
use fs_contact::interface_binding::thermal_friction::elliptic::{
    SourceThermalEllipticCoulomb, SourceThermalEllipticState,
};

/// Original contact geometry, before its material source chooses coefficients.
pub type SolidImpactGeometry<'a> = SolidFrictionalContact<'a, ()>;

/// Actual preimpact kinematics in the geometry-owned normal-first frame.
/// `original_index` indexes both the input geometry and accepted history slices.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolidImpactKinematics {
    pub original_index: usize,
    pub contact_id: u64,
    pub gap_m: f64,
    pub relative_velocity_mps: [f64; 3],
}

/// One explicit ordered-interface source, matched to the original contact ID.
/// Its retained texture must already use that geometry's tangent coordinates.
#[derive(Clone, Copy, Debug)]
pub struct SolidEllipticSourceBinding<'a> {
    pub contact_id: u64,
    pub law: &'a SourceThermalEllipticCoulomb,
    /// Complete caller-owned temperature, nominal pressure and external query.
    /// Actual incoming slip speed replaces declared speed placeholders; each
    /// accumulated-distance coordinate stays at its own accepted value.
    pub conditions: &'a SourceFrictionConditions,
    pub coefficient_role: SourceImpulseCoefficientRole,
}

/// Source identity/index survive a failure; no native or history candidate is
/// published. A cancellation in a source may retain that source's own error.
#[derive(Debug)]
pub enum SolidEllipticSourceError<E = SourceImpulseError> {
    Impact(SolidFrictionalImpactError),
    Source {
        contact_id: u64,
        original_index: usize,
        source: E,
    },
}
impl<E> From<SolidFrictionalImpactError> for SolidEllipticSourceError<E> {
    fn from(error: SolidFrictionalImpactError) -> Self {
        Self::Impact(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for SolidEllipticSourceError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-bound solid directional impact: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for SolidEllipticSourceError<E> {}

type Result<T, E = SourceImpulseError> = core::result::Result<T, SolidEllipticSourceError<E>>;

struct SourceProposal<'a, S> {
    impulse: FrictionalImpulseProposal<EllipticFrictionalContactResult>,
    contacts: Vec<SolidEllipticContact<'a>>,
    candidates: Vec<(usize, S)>,
}

fn prepare<'a, S, E>(
    velocity_mps: &[f64],
    contacts: &[SolidImpactGeometry<'a>],
    histories: &[S],
    controls: SolidFrictionalImpactControls,
    inverse_mass_action: impl FnMut(
        &[f64],
        &Cx<'_>,
    ) -> core::result::Result<Vec<f64>, LinearSolidImpulseError>,
    mut resolve: impl FnMut(
        &S,
        &SolidImpactKinematics,
    ) -> core::result::Result<(EllipticFriction, S), E>,
    cx: &Cx<'_>,
) -> Result<SourceProposal<'a, S>, E> {
    poll(cx)?;
    if contacts.len() != histories.len()
        || !(1..=21).contains(&contacts.len())
        || contacts.len() > controls.contact.maximum_contacts
    {
        return Err(SolidFrictionalImpactError::Invalid(
            "one accepted history per bounded original contact required",
        )
        .into());
    }
    // Validate both quantities before any inactive-source shortcut: NaN,
    // infinity or negative geometry must never evade admission by comparison.
    for contact in contacts {
        poll(cx)?;
        if !contact.gap_m.is_finite()
            || contact.gap_m < 0.0
            || !contact.activation_tolerance_m.is_finite()
            || contact.activation_tolerance_m < 0.0
        {
            return Err(SolidFrictionalImpactError::Invalid(
                "finite nonnegative gap and activation tolerance required",
            )
            .into());
        }
    }
    let mut selected: Vec<_> = contacts
        .iter()
        .map(|c| SolidEllipticContact {
            id: c.id,
            velocity_weights: c.velocity_weights,
            gap_m: c.gap_m,
            activation_tolerance_m: c.activation_tolerance_m,
            prescribed_surface_velocity_mps: c.prescribed_surface_velocity_mps,
            friction_coefficient: EllipticFriction {
                coefficients_mu: [0.0; 2],
                principal_axis: [1.0, 0.0],
            },
            normal_law: c.normal_law,
        })
        .collect();
    let mut candidates = Vec::with_capacity(contacts.len());
    let impulse = propose_with(
        velocity_mps,
        contacts,
        controls,
        inverse_mass_action,
        |velocities, order, delassus, cx| -> Result<_, E> {
            for &i in order {
                poll(cx)?;
                let contact = contacts[i];
                if contact.gap_m > contact.activation_tolerance_m {
                    continue;
                }
                let query = SolidImpactKinematics {
                    original_index: i,
                    contact_id: contact.id,
                    gap_m: contact.gap_m,
                    relative_velocity_mps: velocities[i],
                };
                let (friction, next) = resolve(&histories[i], &query).map_err(|source| {
                    SolidEllipticSourceError::Source {
                        contact_id: contact.id,
                        original_index: i,
                        source,
                    }
                })?;
                selected[i].friction_coefficient = friction;
                candidates.push((i, next));
            }
            let native: Vec<_> = selected
                .iter()
                .zip(velocities)
                .map(|(contact, &velocity)| EllipticFriction::contact(contact, velocity))
                .collect();
            Ok(EllipticFriction::solve(
                &native,
                delassus,
                controls.contact,
                cx,
            )?)
        },
        cx,
    )?;
    Ok(SourceProposal {
        impulse,
        contacts: selected,
        candidates,
    })
}

fn publish<S>(histories: &mut [S], candidates: &mut [(usize, S)]) {
    // Hold old histories until ALL swaps complete; a caller-owned destructor
    // cannot run between publishing individual accepted source states.
    for (index, next) in candidates {
        core::mem::swap(&mut histories[*index], next);
    }
}

fn admit_bindings(
    contacts: &[SolidImpactGeometry<'_>],
    sources: &[SolidEllipticSourceBinding<'_>],
    cx: &Cx<'_>,
) -> Result<()> {
    poll(cx)?;
    if contacts.len() != sources.len() || contacts.len() > 21 {
        return Err(SolidFrictionalImpactError::Invalid(
            "one original source binding per bounded contact required",
        )
        .into());
    }
    for (contact, source) in contacts.iter().zip(sources) {
        poll(cx)?;
        if contact.id != source.contact_id {
            return Err(SolidFrictionalImpactError::Invalid(
                "source binding IDs must match original contact input order",
            )
            .into());
        }
    }
    Ok(())
}

fn sample(
    binding: SolidEllipticSourceBinding<'_>,
    history: &SourceThermalEllipticState,
    query: &SolidImpactKinematics,
    cx: &Cx<'_>,
) -> core::result::Result<(EllipticFriction, SourceThermalEllipticState), SourceImpulseError> {
    let sample = binding.law.sample_preimpact_coefficients(
        history,
        [
            query.relative_velocity_mps[1],
            query.relative_velocity_mps[2],
        ],
        binding.conditions,
        binding.coefficient_role,
        &mut || cx.checkpoint().is_err(),
    )?;
    let friction = EllipticFriction {
        coefficients_mu: sample.coefficients_mu,
        principal_axis: sample.state.constitutive_state().principal_axis(),
    };
    Ok((friction, sample.state))
}

impl ExactLinearSolid {
    /// Query original sources in canonical contact-ID order using actual
    /// incoming native contact velocity, then resolve one simultaneous event.
    ///
    /// Geometry and histories share input order. Inactive contacts are retained
    /// mechanically but never call the source or change their history. Each
    /// callback returns a private coefficient/history candidate; the complete
    /// native impulse, endpoint law and work must pass before either body or
    /// source state publishes. Histories need not implement Clone.
    ///
    /// Callbacks own their source-evaluation/history budgets and must treat
    /// accepted history as read-only. Interior mutation, external side effects
    /// and panics are outside this candidate transaction. No time, pressure,
    /// travel or thermal dissipation is inferred from an impulse.
    pub fn apply_elliptic_impacts_with<S, E>(
        &mut self,
        contacts: &[SolidImpactGeometry<'_>],
        histories: &mut [S],
        controls: SolidFrictionalImpactControls,
        resolve: impl FnMut(
            &S,
            &SolidImpactKinematics,
        ) -> core::result::Result<(EllipticFriction, S), E>,
        cx: &Cx<'_>,
    ) -> Result<SolidEllipticImpactReport, E> {
        let SourceProposal {
            impulse,
            contacts,
            mut candidates,
        } = prepare(
            self.velocity_mps(),
            contacts,
            histories,
            controls,
            |row, cx| impulse::action(&self.inverse_mass, row, cx),
            resolve,
            cx,
        )?;
        let candidate = self
            .trial_momentum_impulse(
                &impulse.impulse_kg_m_s,
                controls.maximum_impulse_error_kg_m_s,
                cx,
            )
            .map_err(SolidFrictionalImpactError::from)?;
        let contact = audit(
            impulse,
            &contacts,
            &candidate.velocity_mps,
            candidate.report.energy.supplied_work_j,
            controls,
            cx,
        )?;
        let impulse = self.accept_momentum_impulse(candidate);
        publish(histories, &mut candidates);
        Ok(SolidEllipticImpactReport { contact, impulse })
    }

    /// Bind the original material-interface sampler to the native event.
    ///
    /// Both original principal claim sets and their authority remain pinned.
    /// Explicit temperature/pressure and fixed external conditions must agree
    /// with each source. Only declared speed coordinates use actual preimpact
    /// relative slip; accepted path, native slip, transported texture and last
    /// finite-step temperature survive unchanged. The role fixes one impact
    /// ellipse and is not an inferred intrastep static/kinetic transition.
    pub fn apply_sourced_elliptic_impacts(
        &mut self,
        contacts: &[SolidImpactGeometry<'_>],
        histories: &mut [SourceThermalEllipticState],
        sources: &[SolidEllipticSourceBinding<'_>],
        controls: SolidFrictionalImpactControls,
        cx: &Cx<'_>,
    ) -> Result<SolidEllipticImpactReport> {
        admit_bindings(contacts, sources, cx)?;
        self.apply_elliptic_impacts_with(
            contacts,
            histories,
            controls,
            |history, query| sample(sources[query.original_index], history, query, cx),
            cx,
        )
    }
}

impl CompleteSpectralSolid {
    /// The same private source/native transaction on every original rigid and
    /// elastic mode. See the exact owner's method for callback ownership and
    /// inactive-contact semantics. No preparation or inertia is cloned.
    pub fn apply_elliptic_impacts_with<S, E>(
        &mut self,
        contacts: &[SolidImpactGeometry<'_>],
        histories: &mut [S],
        controls: SolidFrictionalImpactControls,
        resolve: impl FnMut(
            &S,
            &SolidImpactKinematics,
        ) -> core::result::Result<(EllipticFriction, S), E>,
        cx: &Cx<'_>,
    ) -> Result<SolidEllipticImpactReport, E> {
        let SourceProposal {
            impulse,
            contacts,
            mut candidates,
        } = prepare(
            self.velocity_mps(),
            contacts,
            histories,
            controls,
            |row, cx| self.momentum_response(row, cx),
            resolve,
            cx,
        )?;
        let candidate = self
            .trial_momentum_impulse(
                &impulse.impulse_kg_m_s,
                controls.maximum_impulse_error_kg_m_s,
                cx,
            )
            .map_err(SolidFrictionalImpactError::from)?;
        let contact = audit(
            impulse,
            &contacts,
            &candidate.velocity_mps,
            candidate.report.energy.supplied_work_j,
            controls,
            cx,
        )?;
        let impulse = self.accept_momentum_impulse(candidate);
        publish(histories, &mut candidates);
        Ok(SolidEllipticImpactReport { contact, impulse })
    }

    /// Refresh the original material sources at the complete native incoming
    /// contact speeds and explicit physical conditions. Both source candidates
    /// and all native state publish together after the original physical audit.
    /// See the exact owner's method for role and retained-history boundaries.
    pub fn apply_sourced_elliptic_impacts(
        &mut self,
        contacts: &[SolidImpactGeometry<'_>],
        histories: &mut [SourceThermalEllipticState],
        sources: &[SolidEllipticSourceBinding<'_>],
        controls: SolidFrictionalImpactControls,
        cx: &Cx<'_>,
    ) -> Result<SolidEllipticImpactReport> {
        admit_bindings(contacts, sources, cx)?;
        self.apply_elliptic_impacts_with(
            contacts,
            histories,
            controls,
            |history, query| sample(sources[query.original_index], history, query, cx),
            cx,
        )
    }
}
