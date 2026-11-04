//! Sustained planar contact through original Newmark motion and force laws.
//!
//! This is a bounded, caller-admitted local held-force step. The shared
//! fs-contact owner solves force/displacement complementarity; the original
//! frame owns every mechanical trial, constitutive update and finite rate.

use super::{DynamicHistory, Frame2d, Geometry, Rayleigh, SectionModel, SolidError};
use fs_contact::sustained::NaturalMapScales;
use fs_contact::sustained::solve::{SolveControls, SolveError, SolveReport};
use fs_math::dd::Dd;
use std::cell::Cell;

mod audit;

/// Hard bounds on each indivisible dense native Newmark interval.
const MAX_DOFS: usize = 256;
const MAX_MEMBERS: usize = 256;

/// A fixed planar point trace and translating obstacle for one interval.
///
/// Both original node-major covectors describe the SAME physical contact in
/// an orthonormal normal/tangent frame. Translation weights are dimensionless,
/// rotation weights have metres. Their transpose maps held forces [N] to
/// generalized nodal loads and moments. Internal contacts include both signed
/// traces. Geometry, unit directions and complete feature coverage belong to
/// the mandatory path admission, not to the algebraic rows alone.
#[derive(Clone, Copy, Debug)]
pub struct FrameSustainedContact<'a> {
    pub id: u64,
    pub normal_covector: &'a [f64],
    pub tangent_covector: &'a [f64],
    /// Actual nonnegative initial signed gap [m].
    pub gap_m: f64,
    /// Constant prescribed [normal, tangent] obstacle velocity [m/s].
    pub obstacle_velocity: [f64; 2],
    pub friction_coefficient: f64,
    /// Numerical N/m conditioning only; no physical spring or compliance.
    pub scales: NaturalMapScales,
}

/// Separate original-response budgets and physical acceptance tolerances.
#[derive(Clone, Copy, Debug)]
pub struct FrameSustainedControls {
    pub duration_s: f64,
    pub damping: Rayleigh,
    pub solve: SolveControls,
    /// At most 256 original coordinates and 256 members per native response.
    pub maximum_dofs: usize,
    pub maximum_members: usize,
    /// Conservative f64-equivalent bound for native/adapter numeric staging.
    /// Source metadata, caller inputs and allocator overhead are excluded.
    /// The shared contact solver separately owns its workspace budget.
    pub maximum_workspace_entries: usize,
    /// Actual translational increment, loaded gap and sticking tolerance [m].
    /// Actual negative endpoint gaps are ALWAYS refused.
    pub displacement_tolerance_m: f64,
    /// Actual rotational Newmark increment tolerance [rad], independent of
    /// contact gap/slip tolerances, which remain in metres.
    pub rotation_tolerance_rad: f64,
    /// Loaded normal dwell and sticking/sliding endpoint speed tolerance [m/s].
    pub velocity_tolerance_m_s: f64,
    /// Contact cone/unloaded-force and native translational equilibrium [N].
    pub force_tolerance_n: f64,
    /// Native rotational equilibrium [N m].
    pub moment_tolerance_nm: f64,
    /// Normal force-gap and tangent duality-gap allowance [J].
    pub complementarity_tolerance_j: f64,
    /// Work/energy and numerical positive friction-work allowance [J].
    pub energy_tolerance_j: f64,
}
impl Default for FrameSustainedControls {
    fn default() -> Self {
        Self {
            duration_s: 0.01,
            damping: Rayleigh { a0: 0.0, a1: 0.0 },
            solve: SolveControls::default(),
            maximum_dofs: MAX_DOFS,
            maximum_members: MAX_MEMBERS,
            maximum_workspace_entries: 4_194_304,
            displacement_tolerance_m: 1e-9,
            rotation_tolerance_rad: 1e-9,
            velocity_tolerance_m_s: 1e-8,
            force_tolerance_n: 1e-8,
            moment_tolerance_nm: 1e-8,
            complementarity_tolerance_j: 1e-8,
            energy_tolerance_j: 1e-8,
        }
    }
}

/// Actual native witnesses supplied to the mandatory physical path owner.
///
/// Contacts and held forces use canonical ID order. The owner must establish
/// complete-interval nonpenetration/feature coverage, validity of the fixed
/// point traces, and contact mode: loaded points dwell; sticking has zero slip;
/// sliding has no interior velocity reversal. Newmark endpoint samples alone
/// do NOT establish these claims. The owner also supplies any branch/determinacy
/// evidence needed to adopt this LOCAL numerical candidate.
///
/// Callback-owned state must stay staged until the complete call succeeds.
pub struct FrameSustainedPath<'a, 'contact> {
    pub before: &'a Frame2d,
    pub after: &'a Frame2d,
    pub native_history: &'a DynamicHistory,
    pub contacts: &'a [FrameSustainedContact<'contact>],
    pub forces_n: &'a [[f64; 3]],
    pub external_load: &'a [f64],
    pub actual_displacement: &'a [f64],
    pub duration_s: f64,
    pub damping: Rayleigh,
}

/// One locally solved, physically admitted and atomically published interval.
/// No global existence/uniqueness or automatic contact discovery is claimed.
#[derive(Debug)]
pub struct FrameSustainedStep<P> {
    pub duration_s: f64,
    /// Original Newmark history; time is local to this one interval.
    pub native: DynamicHistory,
    pub contact_ids: Vec<u64>,
    /// Unmodified shared solve report in canonical contact-ID order.
    /// The third tangent force/path coordinate is exactly zero.
    pub solve: SolveReport,
    /// ALL original Newmark responses, including probes and final replay.
    pub response_evaluations: usize,
    pub path_admission: P,
    pub displacement_change: Vec<f64>,
    pub velocity_change: Vec<f64>,
    pub generalized_momentum_change: Vec<f64>,
    /// Integrated support reaction ON the frame at fixed DOFs; zero otherwise.
    /// Includes original consistent-mass cross terms. Supports do zero work.
    pub support_impulse: Vec<f64>,
    pub kinetic_energy_change_j: f64,
    /// Original elastic internal-force work at the actual stored increment.
    pub strain_energy_change_j: f64,
    pub material_damping_energy_j: f64,
    pub rayleigh_damping_energy_j: f64,
    pub external_work_j: f64,
    pub obstacle_work_j: f64,
    /// Signed normal force times relative normal displacement, kept separate
    /// from friction and never relabeled as heat.
    pub normal_contact_work_j: f64,
    pub frictional_work_j: f64,
    pub frictional_dissipation_j: f64,
    pub passivity_residual_j: f64,
    /// Delta kinetic + elastic + damping - external - obstacle - normal - friction.
    pub energy_balance_residual_j: f64,
    /// Stable actual kinetic change against the original native work ledger.
    pub native_energy_balance_residual_j: f64,
}

/// No response, geometry, mode, work or late callback refusal changes the frame.
#[derive(Debug)]
pub enum FrameSustainedError<E> {
    Invalid(&'static str),
    Frame(SolidError),
    Solve(Box<SolveError<SolidError>>),
    Budget {
        quantity: &'static str,
        requested: usize,
        cap: usize,
    },
    ImpactRequired {
        contact_id: u64,
        normal_speed_m_s: f64,
    },
    ContactMode {
        contact_id: u64,
        reason: &'static str,
    },
    Balance {
        quantity: &'static str,
        residual: f64,
        allowed: f64,
    },
    Path(E),
    Cancelled,
}
impl<E: core::fmt::Debug> core::fmt::Display for FrameSustainedError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "frame sustained contact: {self:?}")
    }
}
impl<E: core::fmt::Debug + 'static> core::error::Error for FrameSustainedError<E> {}
impl<E> From<SolidError> for FrameSustainedError<E> {
    fn from(error: SolidError) -> Self {
        if matches!(error, SolidError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Frame(error)
        }
    }
}
fn invalid(what: &str) -> SolidError {
    SolidError::InvalidInput { what: what.into() }
}
fn poll(cancelled: &mut dyn FnMut() -> bool) -> Result<(), SolidError> {
    if cancelled() {
        Err(SolidError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(value: Dd) -> Result<f64, SolidError> {
    let rounded = value.to_f64();
    if value.hi.is_finite() && value.lo.is_finite() && rounded.is_finite() {
        Ok(rounded)
    } else {
        Err(invalid("unrepresentable sustained frame response or work"))
    }
}
fn dot(a: &[f64], b: &[f64]) -> Dd {
    a.iter().zip(b).fold(Dd::ZERO, |sum, (&a, &b)| {
        sum + Dd::from_f64(a) * Dd::from_f64(b)
    })
}
fn actual_paths(
    before: &Frame2d,
    after: &Frame2d,
    contacts: &[FrameSustainedContact<'_>],
    duration: f64,
) -> Result<Vec<[f64; 3]>, SolidError> {
    // Physical geometry uses actual represented positions, never an internal
    // Newmark increment that might disappear on a large absolute displacement.
    let du: Vec<_> = after.u.iter().zip(&before.u).map(|(a, b)| a - b).collect();
    contacts
        .iter()
        .map(|c| {
            Ok([
                finite(
                    Dd::from_f64(c.gap_m) + dot(c.normal_covector, &du)
                        - Dd::from_f64(duration) * Dd::from_f64(c.obstacle_velocity[0]),
                )?,
                finite(
                    dot(c.tangent_covector, &du)
                        - Dd::from_f64(duration) * Dd::from_f64(c.obstacle_velocity[1]),
                )?,
                0.0,
            ])
        })
        .collect()
}
struct Trial {
    frame: Frame2d,
    history: DynamicHistory,
    forces: Vec<[f64; 3]>,
    load: Vec<f64>,
    paths: Vec<[f64; 3]>,
}

impl Frame2d {
    /// Advance ONE sustained-contact interval using original Newmark motion.
    ///
    /// Held contact forces and the complete held external load act together at
    /// both native force-balance endpoints. Resting support therefore balances
    /// gravity without a velocity reset, impact impulse or artificial kick loss.
    /// Every shared force-solver response starts from this same accepted frame.
    /// The final freshly replayed native frame is retained, avoiding an extra,
    /// uncounted motion solve after force convergence.
    ///
    /// Fixed linear geometry, elastic sections, positive free inertia and
    /// stationary supports are required. Source viscosity, retained Rayleigh
    /// stiffness, native material admission and finite endpoint acceleration
    /// keep their original meanings. Generalized loads use N and N m.
    ///
    /// The original Newmark interval is an INDIVISIBLE bounded response:
    /// cancellation is checked before/after it and within the shared contact
    /// solver, not inside native Newmark's Newton iteration. Coordinate/member,
    /// numeric-workspace and original-response limits bound that work.
    ///
    /// Initial closing impacts, inconsistent loaded contact modes, resolved
    /// endpoint slip reversals and negative actual final gaps refuse. Complete
    /// path/feature/mode and any additional determinacy admission are mandatory.
    /// Endpoint force convergence alone is never authority to adopt a trajectory.
    pub fn advance_sustained_contacts<P, E>(
        &mut self,
        contacts: &[FrameSustainedContact<'_>],
        initial_forces_n: &[[f64; 2]],
        external_load: &[f64],
        controls: FrameSustainedControls,
        mut cancelled: impl FnMut() -> bool,
        admit_path: impl FnOnce(&FrameSustainedPath<'_, '_>) -> Result<P, E>,
    ) -> Result<FrameSustainedStep<P>, FrameSustainedError<E>> {
        poll(&mut cancelled)?;
        let h = controls.duration_s;
        if !h.is_finite()
            || h <= 0.0
            || !(0.25 * h * h).is_finite()
            || 0.25 * h * h == 0.0
            || [
                controls.displacement_tolerance_m,
                controls.rotation_tolerance_rad,
                controls.velocity_tolerance_m_s,
                controls.force_tolerance_n,
                controls.moment_tolerance_nm,
                controls.complementarity_tolerance_j,
                controls.energy_tolerance_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
            || contacts.len() != initial_forces_n.len()
            || contacts.is_empty()
            || contacts.len() > 21
            || contacts.len() > controls.solve.maximum_contacts
        {
            return Err(FrameSustainedError::Invalid(
                "duration, tolerances, contact count or warm-start shape",
            ));
        }
        let n = self.ndof();
        for (quantity, requested, cap) in [
            ("native coordinates", n, controls.maximum_dofs.min(MAX_DOFS)),
            (
                "native members",
                self.elements.len(),
                controls.maximum_members.min(MAX_MEMBERS),
            ),
        ] {
            if requested > cap {
                return Err(FrameSustainedError::Budget {
                    quantity,
                    requested,
                    cap,
                });
            }
        }
        let entries = n
            .checked_mul(n)
            .and_then(|v| v.checked_mul(32))
            .and_then(|v| n.checked_mul(128).and_then(|x| v.checked_add(x)))
            .and_then(|v| {
                self.elements
                    .len()
                    .checked_mul(512)
                    .and_then(|x| v.checked_add(x))
            })
            .ok_or(FrameSustainedError::Invalid("native workspace overflow"))?;
        if entries > controls.maximum_workspace_entries {
            return Err(FrameSustainedError::Budget {
                quantity: "native numeric workspace entries",
                requested: entries,
                cap: controls.maximum_workspace_entries,
            });
        }
        if self.geometry!=Geometry::Linear || !self.fixed.iter().any(|fixed|!fixed)
            || self.elements.iter().any(|e|e.sections.iter().any(|s|
                !matches!(s,SectionModel::Elastic{ea,ei} if ea.is_finite()&&*ea>0.0&&ei.is_finite()&&*ei>0.0)))
        {
            return Err(FrameSustainedError::Invalid("sustained frame requires fixed linear geometry, elastic sections and free inertial coordinates"));
        }
        let mass = self.impulse_mass(external_load)?;
        let damping = self.prepare_damping(&mass, controls.damping, h)?;
        let mut order: Vec<_> = (0..contacts.len()).collect();
        order.sort_unstable_by_key(|&i| contacts[i].id);
        if order
            .windows(2)
            .any(|pair| contacts[pair[0]].id == contacts[pair[1]].id)
        {
            return Err(FrameSustainedError::Invalid("duplicate contact ID"));
        }
        let ordered: Vec<_> = order.iter().map(|&i| contacts[i]).collect();
        let initial: Vec<_> = order
            .iter()
            .map(|&i| [initial_forces_n[i][0], initial_forces_n[i][1], 0.0])
            .collect();
        for c in &ordered {
            if !c.gap_m.is_finite()
                || c.gap_m < 0.0
                || [c.normal_covector, c.tangent_covector]
                    .iter()
                    .any(|row| row.len() != n || row.iter().any(|v| !v.is_finite()))
                || c.obstacle_velocity.iter().any(|v| !v.is_finite())
                || !c.friction_coefficient.is_finite()
                || c.friction_coefficient < 0.0
            {
                return Err(FrameSustainedError::Invalid(
                    "finite physical traces, nonnegative gap/friction and finite obstacle motion required",
                ));
            }
            let vn =
                finite(dot(c.normal_covector, &self.vel) - Dd::from_f64(c.obstacle_velocity[0]))?;
            if (c.gap_m == 0.0 && vn < 0.0)
                || (c.gap_m <= controls.displacement_tolerance_m
                    && vn < -controls.velocity_tolerance_m_s)
            {
                return Err(FrameSustainedError::ImpactRequired {
                    contact_id: c.id,
                    normal_speed_m_s: vn,
                });
            }
        }
        let mu: Vec<_> = ordered.iter().map(|c| c.friction_coefficient).collect();
        let scales: Vec<_> = ordered.iter().map(|c| c.scales).collect();
        let count = Cell::new(0usize);
        let mut last = None;
        let solve = fs_contact::sustained::solve::solve(
            &mu,
            &scales,
            &initial,
            controls.solve,
            &mut |forces, stop| {
                poll(stop)?;
                let mut load = external_load.to_vec();
                for (j, value) in load.iter_mut().enumerate() {
                    let mut sum = Dd::from_f64(*value);
                    for (c, force) in ordered.iter().zip(forces) {
                        sum = sum
                            + Dd::from_f64(c.normal_covector[j]) * Dd::from_f64(force[0])
                            + Dd::from_f64(c.tangent_covector[j]) * Dd::from_f64(force[1]);
                    }
                    *value = finite(sum)?;
                }
                let mut frame = self.clone();
                // Charge every original interval, even one returning a native
                // source/solve refusal. No native trial is silently cached away.
                count.set(count.get() + 1);
                let history =
                    frame.newmark_with_load(&[0.0, 0.0], h, (0.0, 0.0), controls.damping, &load)?;
                poll(stop)?;
                let paths = actual_paths(self, &frame, &ordered, h)?;
                last = Some(Trial {
                    frame,
                    history,
                    forces: forces.to_vec(),
                    load,
                    paths: paths.clone(),
                });
                Ok(paths)
            },
            &mut cancelled,
        )
        .map_err(|error| match error {
            SolveError::Cancelled | SolveError::Response(SolidError::Cancelled) => {
                FrameSustainedError::Cancelled
            }
            error => FrameSustainedError::Solve(Box::new(error)),
        })?;
        let trial = last.ok_or(FrameSustainedError::Invalid(
            "missing final native force replay",
        ))?;
        if count.get() != solve.response_evaluations
            || trial.forces != solve.forces_n
            || trial.paths != solve.paths_m
            || solve.forces_n.iter().any(|force| force[2] != 0.0)
            || solve.paths_m.iter().any(|path| path[2] != 0.0)
        {
            return Err(FrameSustainedError::Invalid(
                "final native response or planar auxiliary coordinate mismatch",
            ));
        }
        poll(&mut cancelled)?;
        let measured = audit::check(
            self,
            &trial,
            &ordered,
            external_load,
            controls,
            &mass,
            &damping,
            &mut cancelled,
        )?;
        let path_admission = admit_path(&FrameSustainedPath {
            before: self,
            after: &trial.frame,
            native_history: &trial.history,
            contacts: &ordered,
            forces_n: &solve.forces_n,
            external_load,
            actual_displacement: &measured.displacement_change,
            duration_s: h,
            damping: controls.damping,
        })
        .map_err(FrameSustainedError::Path)?;
        poll(&mut cancelled)?;
        let result = FrameSustainedStep {
            duration_s: h,
            native: trial.history,
            contact_ids: ordered.iter().map(|c| c.id).collect(),
            solve,
            response_evaluations: count.get(),
            path_admission,
            displacement_change: measured.displacement_change,
            velocity_change: measured.velocity_change,
            generalized_momentum_change: measured.generalized_momentum_change,
            support_impulse: measured.support_impulse,
            kinetic_energy_change_j: measured.kinetic_energy_change_j,
            strain_energy_change_j: measured.strain_energy_change_j,
            material_damping_energy_j: measured.material_damping_energy_j,
            rayleigh_damping_energy_j: measured.rayleigh_damping_energy_j,
            external_work_j: measured.external_work_j,
            obstacle_work_j: measured.obstacle_work_j,
            normal_contact_work_j: measured.normal_contact_work_j,
            frictional_work_j: measured.frictional_work_j,
            frictional_dissipation_j: (-measured.frictional_work_j).max(0.0),
            passivity_residual_j: measured.frictional_work_j.max(0.0),
            energy_balance_residual_j: measured.energy_balance_residual_j,
            native_energy_balance_residual_j: measured.native_energy_balance_residual_j,
        };
        // Sole publication: includes the original material/element trial,
        // finite rates, supports and retained/captured Rayleigh reference.
        *self = trial.frame;
        Ok(result)
    }
}
