//! Sustained normal and Coulomb forces on the complete exact solid trajectory.
//!
//! Each contact trial evaluates the original prepared, held-force flow from the
//! same accepted physical state. All original consistent inertia, stiffness and
//! optional physical viscosity remain present. Contact forces are not impulses
//! or artificial springs. A final native replay, physical mode/work checks and
//! caller-owned whole-path admission precede the original atomic acceptance.

use super::{ExactLinearSolid, ExactLinearSolidError, ExactLinearSolidStep};
use fs_contact::sustained::NaturalMapScales;
use fs_contact::sustained::solve::{MAX_CONTACTS, SolveControls, SolveError, SolveReport};
use fs_exec::Cx;
use fs_math::dd::Dd;

mod exact;

/// A fixed physical point trace in one orthonormal normal/two-tangent frame.
///
/// Covectors are dimensionless. Free owners use `free_dofs()` order; the supported
/// owner uses every `represented_dofs()` coordinate, including supports. Their
/// transpose maps force [N] to nodal force, including both signed traces for
/// self contact. Geometry owns interpolation, frame orthonormality, feature
/// coverage and the validity of holding the trace over the complete interval.
#[derive(Clone, Copy, Debug)]
pub struct SolidSustainedContact<'a> {
    pub id: u64,
    pub velocity_weights: [&'a [f64]; 3],
    /// Actual initial nonnegative separation [m], never repaired here.
    pub gap_m: f64,
    /// Constant opposing point velocity in the same contact frame [m/s].
    pub prescribed_surface_velocity_mps: [f64; 3],
    pub friction_coefficient: f64,
    /// Numerical conditioning only; these are not physical stiffnesses.
    pub scales: NaturalMapScales,
}

/// Original solver budgets plus separate physical acceptance and staging limits.
#[derive(Clone, Copy, Debug)]
pub struct SolidSustainedControls {
    pub solve: SolveControls,
    /// Adapter/numerical response scratch, excluding immutable prepared data
    /// and the shared solver's independently bounded workspace.
    pub maximum_workspace_bytes: usize,
    /// Conservative accumulated scalar-product work for every native response
    /// and one final physical audit. Preparation retains its original budget.
    pub maximum_response_multiply_terms: u64,
    pub displacement_tolerance_m: f64,
    pub velocity_tolerance_mps: f64,
    pub force_tolerance_n: f64,
    pub complementarity_tolerance_j: f64,
    pub maximum_work_error_j: f64,
}

/// Original physical witnesses supplied to the mandatory whole-path owner.
///
/// Endpoints cannot prove absence of interior penetration or a hidden slip
/// reversal. The callback must admit the complete trajectory's geometry,
/// fixed point traces, loaded dwell, sticking/sliding mode, and any needed
/// branch/determinacy evidence. Caller state must remain staged until success.
/// Exactness is for the original spatial pencil and held force, not continuum
/// mechanics or an arbitrary force history.
pub struct SolidSustainedPath<'a, 'contact, B = ExactLinearSolid> {
    pub before: &'a B,
    pub displacement_after_m: &'a [f64],
    pub velocity_after_mps: &'a [f64],
    pub displacement_change_m: &'a [f64],
    pub contacts: &'a [SolidSustainedContact<'contact>],
    pub forces_n: &'a [[f64; 3]],
    pub external_force_n: &'a [f64],
    /// Actual summed nodal force used by the original native response.
    pub total_force_n: &'a [f64],
    pub duration_s: f64,
}

/// One physically admitted interval on the original trajectory and clock.
#[derive(Debug)]
pub struct SolidSustainedStep<P, S = ExactLinearSolidStep> {
    pub native: S,
    /// Stable contact-ID order, also used by `solve` and the path callback.
    pub contact_ids: Vec<u64>,
    pub solve: SolveReport,
    pub path_admission: P,
    pub displacement_change_m: Vec<f64>,
    /// Original consistent mass times actual represented velocity change.
    pub momentum_change_kg_m_s: Vec<f64>,
    pub external_work_j: f64,
    pub prescribed_surface_work_j: f64,
    /// Independently integrated actuator work at prescribed native supports.
    /// Zero for owners with homogeneous stationary supports.
    pub support_work_j: f64,
    pub normal_contact_work_j: f64,
    pub frictional_work_j: f64,
    /// Relative sliding loss, separate from native physical viscosity. No heat
    /// destination or material loss partition is inferred from this value.
    pub frictional_dissipation_j: f64,
    pub passivity_residual_j: f64,
    pub contact_work_balance_error_j: f64,
    /// Original energy change + viscosity - external - surface - support - relative work.
    pub energy_balance_error_j: f64,
}

/// Every refusal preserves displacement, velocity, clock and cumulative accounts.
#[derive(Debug)]
pub enum SolidSustainedError<E, N = ExactLinearSolidError> {
    Invalid(&'static str),
    Native(N),
    Solve(Box<SolveError<N>>),
    Budget {
        resource: &'static str,
        required: u128,
        maximum: u128,
    },
    ImpactRequired {
        contact_id: u64,
        normal_speed_mps: f64,
    },
    ContactMode {
        contact_id: u64,
        reason: &'static str,
    },
    Balance {
        quantity: &'static str,
        defect: f64,
        allowance: f64,
    },
    Path(E),
    Cancelled,
}
impl<E: core::fmt::Debug, N: core::fmt::Debug> core::fmt::Display for SolidSustainedError<E, N> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "solid sustained contact refused: {self:?}")
    }
}
impl<E: core::fmt::Debug + 'static, N: core::fmt::Debug + 'static> core::error::Error
    for SolidSustainedError<E, N>
{
}
impl<E> From<ExactLinearSolidError> for SolidSustainedError<E> {
    fn from(error: ExactLinearSolidError) -> Self {
        if matches!(error, ExactLinearSolidError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Native(error)
        }
    }
}

fn finite(value: Dd) -> Result<f64, &'static str> {
    let rounded = value.to_f64();
    if value.hi.is_finite() && value.lo.is_finite() && rounded.is_finite() {
        Ok(rounded)
    } else {
        Err("unrepresentable contact response or work")
    }
}
fn dot(a: &[f64], b: &[f64]) -> Dd {
    a.iter().zip(b).fold(Dd::ZERO, |sum, (&a, &b)| {
        sum + Dd::from_f64(a) * Dd::from_f64(b)
    })
}
fn negative(value: Dd) -> bool {
    value.hi < 0.0 || (value.hi == 0.0 && value.lo < 0.0)
}
fn rounding(terms: &[Dd]) -> Dd {
    terms.iter().fold(Dd::ZERO, |sum, term| {
        sum + term.abs() * Dd::from_f64(128.0 * f64::EPSILON)
    })
}
fn gate<E, N>(
    quantity: &'static str,
    defect: Dd,
    allowance: Dd,
) -> Result<(), SolidSustainedError<E, N>> {
    let defect = finite(defect)?;
    let allowance = finite(allowance)?;
    if allowance < 0.0 || defect.abs() > allowance {
        Err(SolidSustainedError::Balance {
            quantity,
            defect,
            allowance,
        })
    } else {
        Ok(())
    }
}
fn bound<E, N>(
    resource: &'static str,
    required: u128,
    maximum: u128,
) -> Result<(), SolidSustainedError<E, N>> {
    if required > maximum {
        Err(SolidSustainedError::Budget {
            resource,
            required,
            maximum,
        })
    } else {
        Ok(())
    }
}
fn mode<E, N>(contact_id: u64, reason: &'static str) -> SolidSustainedError<E, N> {
    SolidSustainedError::ContactMode { contact_id, reason }
}

impl ExactLinearSolid {
    /// Resolve simultaneous sustained forces through one original held interval.
    ///
    /// Every response shares the retained prepared matrices and reuses O(n+m)
    /// adapter scratch; there is no new exponential, factorization, truncated
    /// model, velocity reset or contact compliance. The shared solver owns its
    /// Newton/Krylov and response limits. Each bounded native matrix contraction
    /// is indivisible; cancellation is polled before and after, and by the solver.
    ///
    /// Only held-force preparations are admitted. Exact physical viscosity,
    /// when present, remains independently integrated by its original owner.
    /// Contact loss is accounted through force work and is never added again to
    /// native viscosity. Initial impacts and unresolved changes of contact mode
    /// require event resolution outside this local sustained-contact interval.
    pub fn advance_sustained_contacts<P, E>(
        &mut self,
        contacts: &[SolidSustainedContact<'_>],
        initial_forces_n: &[[f64; 3]],
        external_force_n: &[f64],
        controls: SolidSustainedControls,
        cx: &Cx<'_>,
        admit_path: impl FnOnce(&SolidSustainedPath<'_, '_>) -> Result<P, E>,
    ) -> Result<SolidSustainedStep<P>, SolidSustainedError<E>> {
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

// Sealed adapter to original complete-coordinate dynamics. Every response is
// read-only; only the final native acceptance may publish a state or clock.
pub(crate) struct State<'a> {
    pub displacement: &'a [f64],
    pub velocity: &'a [f64],
    pub mass: &'a [f64],
    pub stiffness: &'a [f64],
    pub dt_s: f64,
}
pub(crate) struct Work {
    pub applied_j: f64,
    pub viscosity_j: f64,
    pub support_j: f64,
}
pub(crate) trait Native {
    type Owner;
    type Step;
    type Error;
    type Workspace;
    fn owner(&self) -> &Self::Owner;
    fn state(&self) -> State<'_>;
    fn admit(&self, cx: &Cx<'_>) -> Result<(), Self::Error>;
    fn cancelled() -> Self::Error;
    fn is_cancelled(error: &Self::Error) -> bool;
    fn invalid(what: &'static str) -> Self::Error;
    fn workspace_bytes(&self) -> u128;
    fn response_multiply_terms(&self) -> u128;
    fn workspace(&self) -> Result<Self::Workspace, Self::Error>;
    fn evaluate(
        &self,
        force: &[f64],
        scratch: &mut Self::Workspace,
        cx: &Cx<'_>,
    ) -> Result<(), Self::Error>;
    fn displacement<'a>(&self, scratch: &'a Self::Workspace) -> &'a [f64];
    fn velocity<'a>(&self, scratch: &'a Self::Workspace) -> &'a [f64];
    fn work(
        &self,
        force: &[f64],
        scratch: &Self::Workspace,
        cx: &Cx<'_>,
    ) -> Result<Work, Self::Error>;
    fn accept(
        &mut self,
        scratch: Self::Workspace,
        force: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Self::Step, Self::Error>;
}
impl<E, N> From<&'static str> for SolidSustainedError<E, N> {
    fn from(what: &'static str) -> Self {
        Self::Invalid(what)
    }
}
fn poll<E, N>(cx: &Cx<'_>) -> Result<(), SolidSustainedError<E, N>> {
    cx.checkpoint().map_err(|_| SolidSustainedError::Cancelled)
}
fn native_error<E, B: Native>(error: B::Error) -> SolidSustainedError<E, B::Error> {
    if B::is_cancelled(&error) {
        SolidSustainedError::Cancelled
    } else {
        SolidSustainedError::Native(error)
    }
}

pub(crate) fn advance<B: Native, P, E>(
    body: &mut B,
    contacts: &[SolidSustainedContact<'_>],
    initial_forces_n: &[[f64; 3]],
    external_force_n: &[f64],
    controls: SolidSustainedControls,
    cx: &Cx<'_>,
    admit_path: impl FnOnce(&SolidSustainedPath<'_, '_, B::Owner>) -> Result<P, E>,
) -> Result<SolidSustainedStep<P, B::Step>, SolidSustainedError<E, B::Error>> {
    poll(cx)?;
    body.admit(cx).map_err(native_error::<E, B>)?;
    let before = body.state();
    let n = before.displacement.len();
    let m = contacts.len();
    let h = before.dt_s;
    if m == 0
        || m > MAX_CONTACTS
        || m > controls.solve.maximum_contacts
        || initial_forces_n.len() != m
        || external_force_n.len() != n
        || external_force_n.iter().any(|f| !f.is_finite())
        || [
            controls.displacement_tolerance_m,
            controls.velocity_tolerance_mps,
            controls.force_tolerance_n,
            controls.complementarity_tolerance_j,
            controls.maximum_work_error_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
    {
        return Err(SolidSustainedError::Invalid(
            "held preparation, complete finite loads, contact count and tolerances required",
        ));
    }
    let (nn, mm) = (n as u128, m as u128);
    bound(
        "native response scratch bytes",
        body.workspace_bytes() + 8 * (8 * nn + 48 * mm + 4),
        controls.maximum_workspace_bytes as u128,
    )?;
    // Includes source transforms, transition, the optional full loss form,
    // contact transpose/path reductions and one complete physical audit.
    let terms_per_response = body.response_multiply_terms() + 12 * nn * mm + 64 * nn + 128 * mm;
    bound(
        "native response multiply terms",
        2 * terms_per_response,
        controls.maximum_response_multiply_terms as u128,
    )?;
    let mut order: Vec<_> = (0..m).collect();
    order.sort_unstable_by_key(|&i| contacts[i].id);
    if order
        .windows(2)
        .any(|w| contacts[w[0]].id == contacts[w[1]].id)
    {
        return Err(SolidSustainedError::Invalid("duplicate contact identity"));
    }
    let ordered: Vec<_> = order.iter().map(|&i| contacts[i]).collect();
    let initial: Vec<_> = order.iter().map(|&i| initial_forces_n[i]).collect();
    for c in &ordered {
        if !c.gap_m.is_finite()
            || c.gap_m < 0.0
            || !c.friction_coefficient.is_finite()
            || c.friction_coefficient < 0.0
            || c.prescribed_surface_velocity_mps
                .iter()
                .any(|v| !v.is_finite())
            || c.velocity_weights
                .iter()
                .any(|w| w.len() != n || w.iter().any(|v| !v.is_finite()))
        {
            return Err(SolidSustainedError::Invalid(
                "finite complete point traces, nonnegative gap and friction required",
            ));
        }
        let vn = finite(
            dot(c.velocity_weights[0], before.velocity)
                - Dd::from_f64(c.prescribed_surface_velocity_mps[0]),
        )?;
        if (c.gap_m == 0.0 && vn < 0.0)
            || (c.gap_m <= controls.displacement_tolerance_m
                && vn < -controls.velocity_tolerance_mps)
        {
            return Err(SolidSustainedError::ImpactRequired {
                contact_id: c.id,
                normal_speed_mps: vn,
            });
        }
    }
    let friction: Vec<_> = ordered.iter().map(|c| c.friction_coefficient).collect();
    let scales: Vec<_> = ordered.iter().map(|c| c.scales).collect();
    let mut workspace = body.workspace().map_err(native_error::<E, B>)?;
    let mut total_force = vec![0.0; n];
    let mut displacement = vec![0.0; n];
    let mut final_forces = vec![[0.0; 3]; m];
    let mut evaluations = 0usize;
    let mut exhausted = None;
    let mut cancelled = || cx.checkpoint().is_err();
    let result = fs_contact::sustained::solve::solve(
        &friction,
        &scales,
        &initial,
        controls.solve,
        &mut |forces, stop| {
            if stop() {
                return Err(B::cancelled());
            }
            let required = (evaluations as u128 + 2) * terms_per_response;
            if required > controls.maximum_response_multiply_terms as u128 {
                exhausted = Some(required);
                return Err(B::invalid(
                    "sustained native response work budget exhausted",
                ));
            }
            for (j, load) in total_force.iter_mut().enumerate() {
                let mut sum = Dd::from_f64(external_force_n[j]);
                for (c, force) in ordered.iter().zip(forces) {
                    for a in 0..3 {
                        sum = sum + Dd::from_f64(c.velocity_weights[a][j]) * Dd::from_f64(force[a]);
                    }
                }
                *load = finite(sum).map_err(B::invalid)?;
            }
            evaluations += 1;
            body.evaluate(&total_force, &mut workspace, cx)?;
            if stop() {
                return Err(B::cancelled());
            }
            for (i, du) in displacement.iter_mut().enumerate() {
                *du = body.displacement(&workspace)[i] - before.displacement[i];
            }
            let mut paths = Vec::with_capacity(m);
            for c in &ordered {
                let mut path = [0.0; 3];
                for a in 0..3 {
                    let relative = dot(c.velocity_weights[a], &displacement)
                        - Dd::from_f64(h) * Dd::from_f64(c.prescribed_surface_velocity_mps[a]);
                    path[a] = finite(relative + Dd::from_f64(if a == 0 { c.gap_m } else { 0.0 }))
                        .map_err(B::invalid)?;
                }
                paths.push(path);
            }
            final_forces.copy_from_slice(forces);
            Ok(paths)
        },
        &mut cancelled,
    );
    if let Some(required) = exhausted {
        return Err(SolidSustainedError::Budget {
            resource: "native response multiply terms",
            required,
            maximum: controls.maximum_response_multiply_terms as u128,
        });
    }
    let solve = result.map_err(|error| match error {
        SolveError::Cancelled => SolidSustainedError::Cancelled,
        SolveError::Response(error) if B::is_cancelled(&error) => SolidSustainedError::Cancelled,
        error => SolidSustainedError::Solve(Box::new(error)),
    })?;
    if evaluations != solve.response_evaluations || final_forces != solve.forces_n {
        return Err(SolidSustainedError::Invalid(
            "missing original final force replay",
        ));
    }
    poll(cx)?;
    let work = body
        .work(&total_force, &workspace, cx)
        .map_err(native_error::<E, B>)?;
    let measured = audit(
        before,
        body.displacement(&workspace),
        body.velocity(&workspace),
        &displacement,
        &ordered,
        &solve,
        external_force_n,
        work.applied_j,
        work.viscosity_j,
        work.support_j,
        controls,
        cx,
    )?;
    let path_admission = admit_path(&SolidSustainedPath {
        before: body.owner(),
        displacement_after_m: body.displacement(&workspace),
        velocity_after_mps: body.velocity(&workspace),
        displacement_change_m: &displacement,
        contacts: &ordered,
        forces_n: &solve.forces_n,
        external_force_n,
        total_force_n: &total_force,
        duration_s: h,
    })
    .map_err(SolidSustainedError::Path)?;
    poll(cx)?;
    let contact_ids = ordered.iter().map(|c| c.id).collect();
    // Original physical, cumulative-work and viscosity gates precede its
    // sole publication; no fallible contact operation follows acceptance.
    let native = body
        .accept(workspace, &total_force, cx)
        .map_err(native_error::<E, B>)?;
    Ok(SolidSustainedStep {
        native,
        contact_ids,
        solve,
        path_admission,
        displacement_change_m: displacement,
        momentum_change_kg_m_s: measured.momentum,
        external_work_j: measured.external,
        prescribed_surface_work_j: measured.surface,
        support_work_j: work.support_j,
        normal_contact_work_j: measured.normal,
        frictional_work_j: measured.friction,
        frictional_dissipation_j: (-measured.friction).max(0.0),
        passivity_residual_j: measured.friction.max(0.0),
        contact_work_balance_error_j: measured.contact_error,
        energy_balance_error_j: measured.energy_error,
    })
}

struct Measured {
    momentum: Vec<f64>,
    external: f64,
    surface: f64,
    normal: f64,
    friction: f64,
    contact_error: f64,
    energy_error: f64,
}

#[allow(clippy::too_many_arguments)]
fn audit<E, N>(
    before: State<'_>,
    displacement_after: &[f64],
    velocity_after: &[f64],
    du: &[f64],
    contacts: &[SolidSustainedContact<'_>],
    solve: &SolveReport,
    external: &[f64],
    supplied: f64,
    viscosity_loss: f64,
    support_work: f64,
    controls: SolidSustainedControls,
    cx: &Cx<'_>,
) -> Result<Measured, SolidSustainedError<E, N>> {
    let n = before.displacement.len();
    let h = Dd::from_f64(before.dt_s);
    let mut surface = Dd::ZERO;
    let mut normal = Dd::ZERO;
    let mut friction = Dd::ZERO;
    let mut friction_round = Dd::ZERO;
    let mut contact_world = Dd::ZERO;
    let mut round = Dd::ZERO;
    let complementarity = Dd::from_f64(controls.complementarity_tolerance_j);
    for (index, (c, force)) in contacts.iter().zip(&solve.forces_n).enumerate() {
        poll(cx)?;
        let mut relative = [Dd::ZERO; 3];
        let mut velocity0 = [0.0; 3];
        let mut velocity1 = [0.0; 3];
        for a in 0..3 {
            relative[a] = dot(c.velocity_weights[a], du)
                - h * Dd::from_f64(c.prescribed_surface_velocity_mps[a]);
            velocity0[a] = finite(
                dot(c.velocity_weights[a], before.velocity)
                    - Dd::from_f64(c.prescribed_surface_velocity_mps[a]),
            )?;
            velocity1[a] = finite(
                dot(c.velocity_weights[a], velocity_after)
                    - Dd::from_f64(c.prescribed_surface_velocity_mps[a]),
            )?;
            let path = finite(relative[a] + Dd::from_f64(if a == 0 { c.gap_m } else { 0.0 }))?;
            if path.to_bits() != solve.paths_m[index][a].to_bits() {
                return Err(SolidSustainedError::Invalid(
                    "final native path differs from solver replay",
                ));
            }
        }
        let gap = relative[0] + Dd::from_f64(c.gap_m);
        if negative(gap) || force[0] < 0.0 {
            return Err(mode(
                c.id,
                "negative actual gap or normal force; no projection is allowed",
            ));
        }
        let pn = Dd::from_f64(force[0]);
        let cap = Dd::from_f64(c.friction_coefficient) * pn;
        let capacity = finite(cap)?;
        let ft = force[1].hypot(force[2]);
        let slip = finite(relative[1])?.hypot(finite(relative[2])?);
        if !ft.is_finite() || !slip.is_finite() || ft > capacity + controls.force_tolerance_n {
            return Err(mode(
                c.id,
                "tangential force exceeds the original Coulomb disk",
            ));
        }
        gate(
            "normal force-gap complementarity",
            pn * gap,
            complementarity + rounding(&[pn * gap]),
        )?;
        let tangent_work =
            Dd::from_f64(force[1]) * relative[1] + Dd::from_f64(force[2]) * relative[2];
        gate(
            "tangent displacement duality",
            cap * Dd::from_f64(slip) + tangent_work,
            complementarity + rounding(&[cap * Dd::from_f64(slip), tangent_work]),
        )?;
        if force[0] > controls.force_tolerance_n {
            if c.gap_m > controls.displacement_tolerance_m
                || finite(gap)? > controls.displacement_tolerance_m
                || velocity0[0].abs() > controls.velocity_tolerance_mps
                || velocity1[0].abs() > controls.velocity_tolerance_mps
            {
                return Err(mode(
                    c.id,
                    "loaded contact must dwell; arrival or detachment needs event resolution",
                ));
            }
            gate(
                "initial normal complementarity",
                pn * Dd::from_f64(c.gap_m),
                complementarity + rounding(&[pn * Dd::from_f64(c.gap_m)]),
            )?;
            if capacity > controls.force_tolerance_n {
                let inside = ft < capacity - controls.force_tolerance_n;
                for velocity in [velocity0, velocity1] {
                    let speed = velocity[1].hypot(velocity[2]);
                    if !speed.is_finite()
                        || (inside
                            && (slip > controls.displacement_tolerance_m
                                || speed > controls.velocity_tolerance_mps))
                    {
                        return Err(mode(
                            c.id,
                            "interior-cone force requires sustained sticking",
                        ));
                    }
                    let power = Dd::from_f64(force[1]) * Dd::from_f64(velocity[1])
                        + Dd::from_f64(force[2]) * Dd::from_f64(velocity[2]);
                    let capacity_power = cap * Dd::from_f64(speed);
                    gate(
                        "endpoint sliding direction",
                        power + capacity_power,
                        cap * Dd::from_f64(controls.velocity_tolerance_mps)
                            + Dd::from_f64(controls.force_tolerance_n) * Dd::from_f64(speed)
                            + rounding(&[power, capacity_power]),
                    )?;
                }
            }
        }
        let source = h * dot(force, &c.prescribed_surface_velocity_mps);
        let normal_work = pn * relative[0];
        let world = (0..3).fold(Dd::ZERO, |sum, a| {
            sum + Dd::from_f64(force[a]) * dot(c.velocity_weights[a], du)
        });
        surface = surface + source;
        normal = normal + normal_work;
        friction = friction + tangent_work;
        friction_round = friction_round + rounding(&[tangent_work]);
        contact_world = contact_world + world;
        round = round + rounding(&[source, normal_work, tangent_work, world]);
    }
    let external_work = dot(external, du);
    let work_allowance = Dd::from_f64(controls.maximum_work_error_j)
        + round
        + rounding(&[external_work, Dd::from_f64(supplied)]);
    let friction_value = finite(friction)?;
    let friction_allowance = Dd::from_f64(controls.maximum_work_error_j) + friction_round;
    if friction_value > finite(friction_allowance)? {
        return Err(SolidSustainedError::Balance {
            quantity: "positive relative friction work",
            defect: friction_value,
            allowance: finite(friction_allowance)?,
        });
    }
    let contact_error = Dd::from_f64(supplied) - external_work - contact_world;
    gate(
        "actual applied force versus contact work",
        contact_error,
        work_allowance,
    )?;
    gate(
        "contact source and relative work",
        contact_world - surface - normal - friction,
        work_allowance,
    )?;
    let dv: Vec<_> = velocity_after
        .iter()
        .zip(before.velocity)
        .map(|(a, b)| a - b)
        .collect();
    let mut momentum = Vec::with_capacity(n);
    let mut change = Dd::ZERO;
    for i in 0..n {
        poll(cx)?;
        let dp = dot(&before.mass[i * n..(i + 1) * n], &dv);
        momentum.push(finite(dp)?);
        change = change
            + (Dd::from_f64(0.5 * velocity_after[i]) + Dd::from_f64(0.5 * before.velocity[i])) * dp
            + (Dd::from_f64(0.5 * displacement_after[i])
                + Dd::from_f64(0.5 * before.displacement[i]))
                * dot(&before.stiffness[i * n..(i + 1) * n], du);
    }
    let energy_error = change + Dd::from_f64(viscosity_loss)
        - external_work
        - surface
        - normal
        - friction
        - Dd::from_f64(support_work);
    gate(
        "complete original energy and contact work",
        energy_error,
        work_allowance
            + rounding(&[
                change,
                Dd::from_f64(viscosity_loss),
                Dd::from_f64(support_work),
            ]),
    )?;
    Ok(Measured {
        momentum,
        external: finite(external_work)?,
        surface: finite(surface)?,
        normal: finite(normal)?,
        friction: friction_value,
        contact_error: finite(contact_error)?,
        energy_error: finite(energy_error)?,
    })
}
