//! Directional friction on the original complete solid impulse response.
//!
//! The shared contact owner maps the physical tangent ellipse to its existing
//! disk solver. Only its physical impulse returns through the original solid
//! transpose. Every normal/shear and intercontact mobility block is retained;
//! the native witness remains labeled in the contact owner's mapped basis.

/// Original material source queries and native impulses accepted together.
pub mod source;

use super::*;
use fs_contact::frictional_multicontact::elliptic::{
    EllipticContactResiduals, EllipticFriction, EllipticFrictionalContact,
    EllipticFrictionalContactError, EllipticFrictionalContactResult,
    solve_elliptic_frictional_contacts,
};

/// Geometry and covectors retain the original contact meanings. Principal
/// coefficients and their first axis are expressed in that contact's two
/// tangent coordinates, not in a separate world or numerical basis.
pub type SolidEllipticContact<'a> = SolidFrictionalContact<'a, EllipticFriction>;

/// Original-coordinate contact response with the explicitly mapped native law.
pub type SolidEllipticContactImpulseReport =
    SolidFrictionalContactImpulseReport<EllipticFrictionalContactResult, EllipticContactResiduals>;

/// One accepted directional event, including the original solid impulse audit.
pub type SolidEllipticImpactReport =
    SolidFrictionalImpactReport<EllipticFrictionalContactResult, EllipticContactResiduals>;

impl From<EllipticFrictionalContactError> for SolidFrictionalImpactError {
    fn from(error: EllipticFrictionalContactError) -> Self {
        if matches!(error, EllipticFrictionalContactError::Cancelled) {
            Self::Cancelled
        } else {
            Self::EllipticContact(error)
        }
    }
}

impl ImpactResponse for EllipticFrictionalContactResult {
    fn physical(&self) -> PhysicalImpulseView<'_> {
        PhysicalImpulseView {
            impulses: &self.impulses,
            relative_velocities_after: &self.relative_velocities_after,
            relative_work: self.relative_work,
        }
    }
}

impl ImpactFriction for EllipticFriction {
    type Contact = EllipticFrictionalContact;
    type Response = EllipticFrictionalContactResult;
    type Residuals = EllipticContactResiduals;

    fn contact(c: &SolidEllipticContact<'_>, velocity: [f64; 3]) -> Self::Contact {
        EllipticFrictionalContact {
            id: c.id,
            gap: c.gap_m,
            activation_tolerance: c.activation_tolerance_m,
            relative_velocity: velocity,
            normal_law: c.normal_law,
            friction: c.friction_coefficient,
        }
    }

    fn solve(
        contacts: &[Self::Contact],
        delassus: &[f64],
        config: FrictionalContactConfig,
        cx: &Cx<'_>,
    ) -> Result<Self::Response> {
        Ok(solve_elliptic_frictional_contacts(
            contacts,
            delassus,
            config,
            || cx.checkpoint().is_err(),
        )?)
    }

    fn measure(
        proposal: &FrictionalImpulseProposal<Self::Response>,
        contacts: &[SolidEllipticContact<'_>],
        velocities: &[[f64; 3]],
        config: FrictionalContactConfig,
        cx: &Cx<'_>,
    ) -> Result<Self::Residuals> {
        measured_law(proposal, contacts, velocities, config, cx)
    }
}

/// Reconstruct the physical ellipse at the ACTUAL native endpoint. A rounded
/// modal basis must not turn a solved mapped slip into an unadmitted material
/// traction. The ellipse support function supplies maximum-dissipation work:
/// p_t dot v_t + p_n sqrt((mu_1 v_1)^2 + (mu_2 v_2)^2) = 0.
fn measured_law(
    proposal: &FrictionalImpulseProposal<EllipticFrictionalContactResult>,
    contacts: &[SolidEllipticContact<'_>],
    velocities: &[[f64; 3]],
    config: FrictionalContactConfig,
    cx: &Cx<'_>,
) -> Result<EllipticContactResiduals> {
    let mut r = EllipticContactResiduals::default();
    let (mut work, mut work_scale) = (Dd::ZERO, Dd::ZERO);
    for &i in &proposal.order {
        poll(cx)?;
        let delta: [f64; 3] = core::array::from_fn(|a| {
            velocities[i][a] - proposal.law.relative_velocities_after[i][a]
        });
        r.maximum_velocity_residual_m_s = r.maximum_velocity_residual_m_s.max(finite(det::hypot(
            delta[0],
            det::hypot(delta[1], delta[2]),
        ))?);
        for a in 0..3 {
            let term = Dd::from_f64(proposal.law.impulses[i][a])
                * (Dd::from_f64(0.5 * proposal.relative_velocities_before_mps[i][a])
                    + Dd::from_f64(0.5 * velocities[i][a]));
            work = work + term;
            work_scale = work_scale + Dd::from_f64(rounded(term)?.abs());
        }
    }
    for &i in &proposal.law.active_indices {
        poll(cx)?;
        let (p, v) = (proposal.law.impulses[i], velocities[i]);
        let g = finite(v[0] - proposal.law.normal_targets[i])?;
        let friction = contacts[i].friction_coefficient;
        // The source coefficients/axis were already admitted by fs-contact.
        // Normalize directly from the original material axis, independently of
        // the solver's congruence map, without squaring a large input axis.
        let scale = friction.principal_axis[0]
            .abs()
            .max(friction.principal_axis[1].abs());
        let axis = friction.principal_axis.map(|x| x / scale);
        let length = det::hypot(axis[0], axis[1]);
        let [c, s] = axis.map(|x| x / length);
        let material = |vector: [f64; 3]| -> Result<[f64; 2]> {
            Ok([
                rounded(
                    Dd::from_f64(c) * Dd::from_f64(vector[1])
                        + Dd::from_f64(s) * Dd::from_f64(vector[2]),
                )?,
                rounded(
                    -Dd::from_f64(s) * Dd::from_f64(vector[1])
                        + Dd::from_f64(c) * Dd::from_f64(vector[2]),
                )?,
            ])
        };
        let (pt, vt) = (material(p)?, material(v)?);
        let mu = friction.coefficients_mu;
        let reference = mu[0].max(mu[1]);
        let ratios = if reference == 0.0 {
            [0.0; 2]
        } else {
            mu.map(|x| x / reference)
        };
        let weighted = [finite(ratios[0] * vt[0])?, finite(ratios[1] * vt[1])?];
        let weighted_speed = finite(det::hypot(weighted[0], weighted[1]))?;
        let mut gauge = [0.0; 2];
        let mut forbidden = 0.0_f64;
        for a in 0..2 {
            if mu[a] == 0.0 {
                forbidden = forbidden.max(pt[a].abs());
            } else {
                gauge[a] = finite(pt[a] / mu[a])?;
            }
        }
        let violation = if reference == 0.0 {
            finite(det::hypot(pt[0], pt[1]))?
        } else {
            finite(reference * (det::hypot(gauge[0], gauge[1]) - p[0]).max(0.0))?.max(forbidden)
        };
        let support_work = finite(finite(reference * p[0])? * weighted_speed)?;
        let duality = rounded(
            Dd::from_f64(support_work)
                + Dd::from_f64(p[1]) * Dd::from_f64(v[1])
                + Dd::from_f64(p[2]) * Dd::from_f64(v[2]),
        )?;
        r.normal_impulse_violation = r.normal_impulse_violation.max((-p[0]).max(0.0));
        r.normal_velocity_violation = r.normal_velocity_violation.max((-g).max(0.0));
        r.normal_complementarity = r.normal_complementarity.max(finite(p[0] * g)?.abs());
        r.ellipse_impulse_violation = r.ellipse_impulse_violation.max(violation);
        r.tangential_duality_gap = r.tangential_duality_gap.max(duality.abs());
        // Free slip along a zero-mu axis must not saturate the other axis's
        // sublimit static traction. For two positive coefficients, test the
        // unscaled physical speed so a weak principal axis cannot hide motion.
        let resisted_speed = if mu.contains(&0.0) {
            weighted_speed
        } else {
            finite(det::hypot(v[1], v[2]))?
        };
        if resisted_speed > config.velocity_tolerance {
            let expected = if weighted_speed == 0.0 {
                [0.0; 2]
            } else {
                [
                    finite(-mu[0] * p[0] * (weighted[0] / weighted_speed))?,
                    finite(-mu[1] * p[0] * (weighted[1] / weighted_speed))?,
                ]
            };
            r.sliding_impulse_residual = r.sliding_impulse_residual.max(finite(det::hypot(
                finite(pt[0] - expected[0])?,
                finite(pt[1] - expected[1])?,
            ))?);
        }
    }
    r.coordinate_work_residual_j =
        rounded(work - Dd::from_f64(proposal.law.transformed.relative_work))?;
    let work_allowance =
        finite(config.energy_tolerance + 128.0 * f64::EPSILON * rounded(work_scale)?)?;
    if r.normal_impulse_violation > config.impulse_tolerance
        || r.normal_velocity_violation > config.velocity_tolerance
        || r.normal_complementarity > config.complementarity_tolerance
        || r.ellipse_impulse_violation > config.impulse_tolerance
        || r.sliding_impulse_residual > config.impulse_tolerance
        || r.tangential_duality_gap > config.complementarity_tolerance
        || r.maximum_velocity_residual_m_s > config.velocity_tolerance
        || r.coordinate_work_residual_j.abs() > work_allowance
    {
        return Err(EllipticFrictionalContactError::ResidualNotSatisfied(r).into());
    }
    Ok(r)
}

impl ExactLinearSolid {
    /// Resolve simultaneous directional friction through the original complete
    /// solid mass, native momentum trial and physical endpoint/work audits.
    ///
    /// Equal principal coefficients recover the original circular law. A zero
    /// coefficient admits exact free slip along that material axis. Original
    /// geometry, positions, clock, material preparation and viscosity accounts
    /// stay with the trajectory; impact loss is not assigned to physical heat.
    /// The law's determinacy witness retains its explicitly mapped coordinates.
    /// Cancellation or any late native refusal leaves all accepted state intact.
    pub fn apply_elliptic_impacts(
        &mut self,
        contacts: &[SolidEllipticContact<'_>],
        controls: SolidFrictionalImpactControls,
        cx: &Cx<'_>,
    ) -> Result<SolidEllipticImpactReport> {
        let proposal = propose(
            self.velocity_mps(),
            contacts,
            controls,
            |row, cx| impulse::action(&self.inverse_mass, row, cx),
            cx,
        )?;
        let candidate = self.trial_momentum_impulse(
            &proposal.impulse_kg_m_s,
            controls.maximum_impulse_error_kg_m_s,
            cx,
        )?;
        let contact = audit(
            proposal,
            contacts,
            &candidate.velocity_mps,
            candidate.report.energy.supplied_work_j,
            controls,
            cx,
        )?;
        Ok(SolidEllipticImpactReport {
            contact,
            impulse: self.accept_momentum_impulse(candidate),
        })
    }
}

impl CompleteSpectralSolid {
    /// Apply the same directional event to every admitted original rigid and
    /// elastic mode, without truncation or re-preparation. Native readback gates
    /// the physical ellipse even when an admitted basis has rounding error.
    /// See the dense owner's method for work and no-claim boundaries.
    pub fn apply_elliptic_impacts(
        &mut self,
        contacts: &[SolidEllipticContact<'_>],
        controls: SolidFrictionalImpactControls,
        cx: &Cx<'_>,
    ) -> Result<SolidEllipticImpactReport> {
        let proposal = propose(
            self.velocity_mps(),
            contacts,
            controls,
            |row, cx| self.momentum_response(row, cx),
            cx,
        )?;
        let candidate = self.trial_momentum_impulse(
            &proposal.impulse_kg_m_s,
            controls.maximum_impulse_error_kg_m_s,
            cx,
        )?;
        let contact = audit(
            proposal,
            contacts,
            &candidate.velocity_mps,
            candidate.report.energy.supplied_work_j,
            controls,
            cx,
        )?;
        Ok(SolidEllipticImpactReport {
            contact,
            impulse: self.accept_momentum_impulse(candidate),
        })
    }
}
