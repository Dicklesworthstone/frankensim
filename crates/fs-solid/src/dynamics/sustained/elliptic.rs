//! Sustained directional friction on the original complete solid trajectory.
//!
//! The original point trace is transformed by `A^T`, with physical force
//! `f = A q`. The existing simultaneous Coulomb solve then sees a disk while
//! the original solid receives the complete transpose load. No body, inertia,
//! penalty spring or dissipative integrator is substituted. Physical witnesses
//! and transformation/work checks precede the native atomic acceptance.

use super::{
    ExactLinearSolid, Native, SolidSustainedContact, SolidSustainedControls, SolidSustainedError,
    SolidSustainedStep, bound, dot, finite, gate, poll, rounding,
};
use fs_contact::frictional_multicontact::elliptic::EllipticFriction;
use fs_contact::sustained::elliptic::EllipticForceCoordinates;
use fs_contact::sustained::{ContactLawError, NaturalMapScales};
use fs_exec::Cx;
use fs_math::dd::Dd;

mod static_kinetic;
pub(crate) use static_kinetic::advance as advance_static_kinetic;

/// Original orthonormal point trace with an oriented Coulomb friction ellipse.
///
/// Coefficients and the principal axis are expressed in this contact's two
/// physical tangent directions. The law is held throughout the interval.
#[derive(Clone, Copy, Debug)]
pub struct SolidSustainedEllipticContact<'a> {
    pub id: u64,
    pub velocity_weights: [&'a [f64]; 3],
    pub gap_m: f64,
    pub prescribed_surface_velocity_mps: [f64; 3],
    pub friction: EllipticFriction,
    /// Numerical N/m scales in solver coordinates, never physical stiffness.
    pub scales: NaturalMapScales,
}

/// Physical witnesses for mandatory complete-interval geometry/mode admission.
///
/// All contacts, forces and paths use the original physical contact frame and
/// stable ID order. The caller must admit fixed traces/principal directions,
/// whole-path nonpenetration, loaded dwell and sticking/sliding direction;
/// endpoints alone do not establish them.
pub struct SolidSustainedEllipticPath<'a, 'contact, B = ExactLinearSolid> {
    pub before: &'a B,
    pub displacement_after_m: &'a [f64],
    pub velocity_after_mps: &'a [f64],
    pub displacement_change_m: &'a [f64],
    pub contacts: &'a [SolidSustainedEllipticContact<'contact>],
    pub forces_n: &'a [[f64; 3]],
    /// Actual endpoint gap and accumulated two-tangent slip in metres.
    pub paths_m: &'a [[f64; 3]],
    pub external_force_n: &'a [f64],
    pub total_force_n: &'a [f64],
    pub duration_s: f64,
}

/// Admitted original motion plus explicitly separated physical/solver frames.
#[derive(Debug)]
pub struct SolidSustainedEllipticStep<P, S = super::ExactLinearSolidStep> {
    /// Original native motion and work accounts. Only `solve`'s force/path
    /// coordinates, residual and nonlinear diagnostics use the transformed
    /// tangent basis; its normal coordinates and contact IDs are unchanged.
    pub transformed: SolidSustainedStep<P, S>,
    /// Physical normal/two-tangent forces, in `transformed.contact_ids` order.
    pub forces_n: Vec<[f64; 3]>,
    /// Actual original physical gap/slip, independently reconstructed before
    /// acceptance, in the same contact order.
    pub paths_m: Vec<[f64; 3]>,
    /// Actual physical ellipse selected for each contact, in the same order.
    /// A static/kinetic solve may select different pairs at different contacts.
    pub friction: Vec<EllipticFriction>,
    pub coordinates: Vec<EllipticForceCoordinates>,
}

impl ExactLinearSolid {
    /// Resolve static capture and kinetic sliding with distinct principal
    /// coefficient pairs on each contact's original material axes.
    ///
    /// `contacts` carries the static ellipse. Each kinetic coefficient must
    /// lie between zero and its corresponding static coefficient; the two
    /// principal directions may weaken by different proportions. Incoming
    /// relative motion immediately uses kinetic strength. A resting contact
    /// first attempts static capture, then every coupled force is re-solved
    /// from the same native predecessor if capture fails. The original
    /// held-trajectory admission still excludes interior stops or reversals.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_sustained_static_kinetic_elliptic_contacts<P, E>(
        &mut self,
        contacts: &[SolidSustainedEllipticContact<'_>],
        kinetic_coefficients_mu: &[[f64; 2]],
        initial_forces_n: &[[f64; 3]],
        external_force_n: &[f64],
        controls: SolidSustainedControls,
        cx: &Cx<'_>,
        admit_path: impl FnOnce(&SolidSustainedEllipticPath<'_, '_>) -> Result<P, E>,
    ) -> Result<SolidSustainedEllipticStep<P>, SolidSustainedError<E>> {
        advance_static_kinetic(
            self,
            contacts,
            kinetic_coefficients_mu,
            initial_forces_n,
            external_force_n,
            controls,
            cx,
            admit_path,
        )
    }

    /// Advance held directional contact through the original exact solid flow.
    ///
    /// Two positive principal coefficients, or equal coefficients including
    /// zero, are supported. A single free-slip axis is refused explicitly;
    /// no artificial coefficient regularizes a singular ellipse. Sticking
    /// tolerances are conservatively tightened by the smallest tangent scale.
    /// Every physical/solver witness and caller admission is completed before
    /// the original state, clock and work accounts are published together.
    pub fn advance_sustained_elliptic_contacts<P, E>(
        &mut self,
        contacts: &[SolidSustainedEllipticContact<'_>],
        initial_forces_n: &[[f64; 3]],
        external_force_n: &[f64],
        controls: SolidSustainedControls,
        cx: &Cx<'_>,
        admit_path: impl FnOnce(&SolidSustainedEllipticPath<'_, '_>) -> Result<P, E>,
    ) -> Result<SolidSustainedEllipticStep<P>, SolidSustainedError<E>> {
        advance(
            self,
            contacts,
            initial_forces_n,
            external_force_n,
            controls,
            cx,
            admit_path,
        )
    }
}

fn law_error<E, N>(error: ContactLawError) -> SolidSustainedError<E, N> {
    match error {
        ContactLawError::Invalid(reason) | ContactLawError::NonFinite(reason) => {
            SolidSustainedError::Invalid(reason)
        }
    }
}

// A wrapper admission may refuse its own physical coordinate/work checks as
// well as forwarding a caller refusal. Preserve either cause without nesting
// the public error type or performing fallible work after native acceptance.
fn flatten<E, N>(
    error: SolidSustainedError<SolidSustainedError<E, N>, N>,
) -> SolidSustainedError<E, N> {
    match error {
        SolidSustainedError::Invalid(reason) => SolidSustainedError::Invalid(reason),
        SolidSustainedError::Native(error) => SolidSustainedError::Native(error),
        SolidSustainedError::Solve(error) => SolidSustainedError::Solve(error),
        SolidSustainedError::Budget {
            resource,
            required,
            maximum,
        } => SolidSustainedError::Budget {
            resource,
            required,
            maximum,
        },
        SolidSustainedError::ImpactRequired {
            contact_id,
            normal_speed_mps,
        } => SolidSustainedError::ImpactRequired {
            contact_id,
            normal_speed_mps,
        },
        SolidSustainedError::ContactMode { contact_id, reason } => {
            SolidSustainedError::ContactMode { contact_id, reason }
        }
        SolidSustainedError::Balance {
            quantity,
            defect,
            allowance,
        } => SolidSustainedError::Balance {
            quantity,
            defect,
            allowance,
        },
        SolidSustainedError::Path(error) => error,
        SolidSustainedError::Cancelled => SolidSustainedError::Cancelled,
    }
}

fn staging_requirements(n: usize, m: usize) -> (u128, u128) {
    let (nn, mm) = (n as u128, m as u128);
    let bytes_per_contact = core::mem::size_of::<SolidSustainedEllipticContact<'_>>()
        + core::mem::size_of::<SolidSustainedContact<'_>>()
        + core::mem::size_of::<EllipticForceCoordinates>()
        + core::mem::size_of::<EllipticFriction>()
        + core::mem::size_of::<[Vec<f64>; 2]>()
        + core::mem::size_of::<usize>()
        + 3 * core::mem::size_of::<[f64; 3]>();
    (
        16 * nn * mm + mm * bytes_per_contact as u128 + 512,
        16 * nn * mm + 256 * mm,
    )
}

pub(crate) fn advance<B: Native, P, E>(
    body: &mut B,
    contacts: &[SolidSustainedEllipticContact<'_>],
    initial_forces_n: &[[f64; 3]],
    external_force_n: &[f64],
    controls: SolidSustainedControls,
    cx: &Cx<'_>,
    admit_path: impl FnOnce(&SolidSustainedEllipticPath<'_, '_, B::Owner>) -> Result<P, E>,
) -> Result<SolidSustainedEllipticStep<P, B::Step>, SolidSustainedError<E, B::Error>> {
    advance_selected(
        body,
        contacts,
        initial_forces_n,
        external_force_n,
        controls,
        cx,
        |_, _, _| Ok(()),
        admit_path,
    )
}

#[allow(clippy::too_many_arguments)]
fn advance_selected<B: Native, P, E>(
    body: &mut B,
    contacts: &[SolidSustainedEllipticContact<'_>],
    initial_forces_n: &[[f64; 3]],
    external_force_n: &[f64],
    controls: SolidSustainedControls,
    cx: &Cx<'_>,
    select_trial: impl FnOnce(
        &super::SolidSustainedPath<'_, '_, B::Owner>,
        &super::SolveReport,
        &[EllipticForceCoordinates],
    ) -> Result<(), E>,
    admit_path: impl FnOnce(&SolidSustainedEllipticPath<'_, '_, B::Owner>) -> Result<P, E>,
) -> Result<SolidSustainedEllipticStep<P, B::Step>, SolidSustainedError<E, B::Error>> {
    poll(cx)?;
    let n = body.state().displacement.len();
    let m = contacts.len();
    if m == 0
        || m > super::MAX_CONTACTS
        || m > controls.solve.maximum_contacts
        || initial_forces_n.len() != m
        || external_force_n.len() != n
    {
        return Err(SolidSustainedError::Invalid(
            "complete bounded elliptic contact inputs required",
        ));
    }
    let (staging_bytes, staging_terms) = staging_requirements(n, m);
    bound(
        "elliptic contact staging bytes",
        staging_bytes,
        controls.maximum_workspace_bytes as u128,
    )?;
    bound(
        "elliptic contact transformation work",
        staging_terms,
        controls.maximum_response_multiply_terms as u128,
    )?;
    let mut inner_controls = controls;
    inner_controls.maximum_workspace_bytes -= staging_bytes as usize;
    inner_controls.maximum_response_multiply_terms -= staging_terms as u64;

    let mut order: Vec<_> = (0..m).collect();
    order.sort_unstable_by_key(|&i| contacts[i].id);
    let ordered: Vec<_> = order.iter().map(|&i| contacts[i]).collect();
    let mut coordinates = Vec::with_capacity(m);
    let mut tangent_rows = Vec::with_capacity(m);
    let mut initial = Vec::with_capacity(m);
    for (&i, contact) in order.iter().zip(&ordered) {
        poll(cx)?;
        if contact
            .velocity_weights
            .iter()
            .any(|row| row.len() != n || row.iter().any(|v| !v.is_finite()))
        {
            return Err(SolidSustainedError::Invalid(
                "finite complete original elliptic point traces required",
            ));
        }
        let coordinate = EllipticForceCoordinates::prepare(contact.friction).map_err(law_error)?;
        inner_controls.displacement_tolerance_m = inner_controls
            .displacement_tolerance_m
            .min(controls.displacement_tolerance_m * coordinate.minimum_scale());
        inner_controls.velocity_tolerance_mps = inner_controls
            .velocity_tolerance_mps
            .min(controls.velocity_tolerance_mps * coordinate.minimum_scale());
        initial.push(
            coordinate
                .force_to_solver(initial_forces_n[i])
                .map_err(law_error)?,
        );
        let mut rows = [Vec::with_capacity(n), Vec::with_capacity(n)];
        for k in 0..n {
            let row = coordinate
                .path_to_solver([
                    0.0,
                    contact.velocity_weights[1][k],
                    contact.velocity_weights[2][k],
                ])
                .map_err(law_error)?;
            rows[0].push(row[1]);
            rows[1].push(row[2]);
        }
        tangent_rows.push(rows);
        coordinates.push(coordinate);
    }
    let mut transformed_contacts = Vec::with_capacity(m);
    for ((contact, coordinate), rows) in ordered.iter().zip(&coordinates).zip(&tangent_rows) {
        transformed_contacts.push(SolidSustainedContact {
            id: contact.id,
            velocity_weights: [contact.velocity_weights[0], &rows[0], &rows[1]],
            gap_m: contact.gap_m,
            prescribed_surface_velocity_mps: coordinate
                .path_to_solver(contact.prescribed_surface_velocity_mps)
                .map_err(law_error)?,
            friction_coefficient: coordinate.coefficient(),
            scales: contact.scales,
        });
    }
    let mut forces_n = Vec::with_capacity(m);
    let mut paths_m = Vec::with_capacity(m);
    let friction = ordered.iter().map(|contact| contact.friction).collect();
    let transformed = super::advance_with_selection(
        body,
        &transformed_contacts,
        None,
        &initial,
        external_force_n,
        inner_controls,
        cx,
        |path, solve| select_trial(path, solve, &coordinates).map_err(SolidSustainedError::Path),
        |path| {
            poll(cx)?;
            let h = Dd::from_f64(path.duration_s);
            let mut physical_work = Dd::ZERO;
            let mut solver_work = Dd::ZERO;
            let mut physical_source = Dd::ZERO;
            let mut solver_source = Dd::ZERO;
            let mut round = Dd::ZERO;
            for (index, ((contact, coordinate), q)) in ordered
                .iter()
                .zip(&coordinates)
                .zip(path.forces_n)
                .enumerate()
            {
                let force = coordinate.force_to_physical(*q).map_err(law_error)?;
                let mut physical_path = [0.0; 3];
                for a in 0..3 {
                    let relative = dot(contact.velocity_weights[a], path.displacement_change_m)
                        - h * Dd::from_f64(contact.prescribed_surface_velocity_mps[a]);
                    let transformed_relative = dot(
                        path.contacts[index].velocity_weights[a],
                        path.displacement_change_m,
                    ) - h * Dd::from_f64(
                        path.contacts[index].prescribed_surface_velocity_mps[a],
                    );
                    physical_path[a] =
                        finite(relative + Dd::from_f64(if a == 0 { contact.gap_m } else { 0.0 }))?;
                    let pw = Dd::from_f64(force[a]) * relative;
                    let qw = Dd::from_f64(q[a]) * transformed_relative;
                    physical_work = physical_work + pw;
                    solver_work = solver_work + qw;
                    round = round + rounding(&[pw, qw]);
                }
                let ps = h * dot(&force, &contact.prescribed_surface_velocity_mps);
                let qs = h * dot(q, &path.contacts[index].prescribed_surface_velocity_mps);
                physical_source = physical_source + ps;
                solver_source = solver_source + qs;
                round = round + rounding(&[ps, qs]);
                forces_n.push(force);
                paths_m.push(physical_path);
            }
            let allowance = Dd::from_f64(controls.maximum_work_error_j) + round;
            gate(
                "elliptic physical relative work",
                physical_work - solver_work,
                allowance,
            )?;
            gate(
                "elliptic physical moving-surface work",
                physical_source - solver_source,
                allowance,
            )?;
            // Independently reconstruct the actual original transpose load before
            // exposing physical force witnesses or allowing the native commit.
            for k in 0..n {
                let mut original = Dd::from_f64(external_force_n[k]);
                let mut force_round = rounding(&[original]);
                for (contact, force) in ordered.iter().zip(&forces_n) {
                    for a in 0..3 {
                        let term =
                            Dd::from_f64(contact.velocity_weights[a][k]) * Dd::from_f64(force[a]);
                        original = original + term;
                        force_round = force_round + rounding(&[term]);
                    }
                }
                let actual = Dd::from_f64(path.total_force_n[k]);
                gate(
                    "elliptic original nodal force",
                    original - actual,
                    Dd::from_f64(controls.force_tolerance_n) + force_round + rounding(&[actual]),
                )?;
            }
            poll(cx)?;
            admit_path(&SolidSustainedEllipticPath {
                before: path.before,
                displacement_after_m: path.displacement_after_m,
                velocity_after_mps: path.velocity_after_mps,
                displacement_change_m: path.displacement_change_m,
                contacts: &ordered,
                forces_n: &forces_n,
                paths_m: &paths_m,
                external_force_n: path.external_force_n,
                total_force_n: path.total_force_n,
                duration_s: path.duration_s,
            })
            .map_err(SolidSustainedError::Path)
        },
    )
    .map_err(flatten)?;
    Ok(SolidSustainedEllipticStep {
        transformed,
        forces_n,
        paths_m,
        friction,
        coordinates,
    })
}
