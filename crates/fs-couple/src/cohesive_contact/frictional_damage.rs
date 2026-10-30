//! Progressive bond damage with compression-only friction on the damaged area.
//!
//! Explicit model: the ORIGINAL cohesive free energy, damage law and surviving
//! bond viscosity are retained. A damaged fraction `d` supports an additional
//! rigid Coulomb traction of capacity `mu d Kc max(-g_n,0)`. Compression storage
//! is counted once. This friction has no elastic spring or extra stored energy;
//! its nonnegative sliding work is reported separately, without a heat sink.
//! The intact/damaged-area decomposition follows the model family of Alfano &
//! Sacco (2006), DOI 10.1002/nme.1728. This declared rigid-friction specialization
//! does NOT reproduce that paper's full elastic-slip internal-variable model.
//!
//! Mean friction uses the material owner's analytic integral of `d(s) p_c(s)`
//! on its complete straight path, including damage growth and opening. Neither
//! endpoint damage nor total normal traction is used as a pressure shortcut.
//! Coefficients are the independently sourced residual-surface coefficients;
//! selecting this owner explicitly extends them to each damaged area fraction.
//!
//! All entries of reciprocal SPD normal/two-tangent compliance are retained.
//! A bounded scalar friction multiplier calls the ORIGINAL cohesive load-line
//! solve through a Schur elimination. Softening can admit multiple equilibria:
//! this is a local discrete-path solve, not a uniqueness or stability theorem.
//! Refusals publish no state. Resolve intrastep reversals with timestep refinement.

use super::reclosure::{ReclosureActivation, ResidualCoulombSource};
use super::{CohesiveContactConfig, CohesiveContactError, CohesiveLoadLine};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::MaterialStateId;
use fs_material::cohesive::source::ResolvedBilinearCohesive;
use fs_material::cohesive::{CohesiveError, CohesivePoint, CohesiveRateStep};
use fs_material::state_point::MaterialStatePointError;
use fs_tribo::TriboError;

pub(crate) mod local;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DamageFrictionRegime {
    /// No damaged-area compression anywhere on the accepted path. Surviving
    /// cohesive bonds and their physical viscosity may still carry traction.
    Inactive,
    Sticking,
    Sliding,
}

/// One total material-evaluation budget includes every inner solve, path
/// measure, friction probe and final audit. No hidden per-probe reset.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DamageFrictionConfig {
    pub contact: CohesiveContactConfig,
    pub max_friction_iterations: usize,
    pub force_absolute_tolerance_n: f64,
    pub force_relative_tolerance: f64,
    pub work_absolute_tolerance_j: f64,
    pub work_relative_tolerance: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum DamageFrictionError {
    Invalid(&'static str),
    Source(MaterialStatePointError),
    Material(CohesiveError),
    Contact(CohesiveContactError),
    Friction(TriboError),
    StateMismatch,
    Cancelled,
    Budget,
    Range(&'static str),
    NotConverged,
    ForceMismatch {
        component: usize,
        residual_n: f64,
        tolerance_n: f64,
    },
    WorkMismatch {
        residual_j: f64,
        tolerance_j: f64,
    },
}
impl core::fmt::Display for DamageFrictionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "damaged-area cohesive friction refused: {self:?}")
    }
}
impl std::error::Error for DamageFrictionError {}
impl From<MaterialStatePointError> for DamageFrictionError {
    fn from(e: MaterialStatePointError) -> Self {
        Self::Source(e)
    }
}
impl From<CohesiveError> for DamageFrictionError {
    fn from(e: CohesiveError) -> Self {
        Self::Material(e)
    }
}
impl From<CohesiveContactError> for DamageFrictionError {
    fn from(e: CohesiveContactError) -> Self {
        match e {
            CohesiveContactError::Cancelled => Self::Cancelled,
            CohesiveContactError::Budget(_) => Self::Budget,
            e => Self::Contact(e),
        }
    }
}
impl From<TriboError> for DamageFrictionError {
    fn from(e: TriboError) -> Self {
        Self::Friction(e)
    }
}

/// One canonical original material point plus friction branch occupancy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DamageFrictionState<'a> {
    identity: ContentHash,
    point: CohesivePoint<'a>,
    occupancy: DamageFrictionRegime,
}
impl<'a> DamageFrictionState<'a> {
    #[must_use]
    pub const fn point(&self) -> CohesivePoint<'a> {
        self.point
    }
    #[must_use]
    pub const fn occupancy(&self) -> DamageFrictionRegime {
        self.occupancy
    }

    /// A zero net step displacement does not imply rest at its endpoint.
    /// Preserve kinetic occupancy from the receiving bodies' actual velocity.
    pub fn with_endpoint_tangential_velocity(
        mut self,
        velocity_m_s: [f64; 2],
    ) -> Result<Self, DamageFrictionError> {
        for v in velocity_m_s {
            finite(v)?;
        }
        if self.point.jump_m()[0] < 0.0
            && self.point.damage() > 0.0
            && velocity_m_s.iter().any(|v| *v != 0.0)
        {
            self.occupancy = DamageFrictionRegime::Sliding;
        }
        Ok(self)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct DamageFrictionAudit<'a> {
    pub state: DamageFrictionState<'a>,
    pub material: CohesiveRateStep<'a>,
    pub regime: DamageFrictionRegime,
    /// Original bonded area times the exact mean of damaged-area compression.
    pub damaged_compression_force_n: f64,
    pub friction_force_n: [f64; 2],
    pub expected_force_n: [f64; 3],
    pub force_residual_n: [f64; 3],
    pub force_tolerance_n: [f64; 3],
    pub friction_dissipation_j: f64,
    pub work_residual_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct DamageFrictionStep<'a> {
    pub audit: DamageFrictionAudit<'a>,
    pub mean_force_n: [f64; 3],
    pub compatibility_residual_m: [f64; 3],
    pub compatibility_tolerance_m: f64,
    pub iterations: usize,
    pub material_evaluations: usize,
}

/// Reusable source-bound solve. Only the original small LU workspace allocates
/// at construction; trials retain scratch and never clone either source owner.
pub struct DamageFrictionLoadLine<'a> {
    cohesive: &'a ResolvedBilinearCohesive,
    friction: &'a ResidualCoulombSource,
    area: f64,
    compliance: [[f64; 3]; 3],
    config: DamageFrictionConfig,
    identity: ContentHash,
    static_mu: f64,
    kinetic_mu: f64,
    tangent_scale: f64,
    inner: CohesiveLoadLine,
}
impl<'a> DamageFrictionLoadLine<'a> {
    pub fn new(
        cohesive: &'a ResolvedBilinearCohesive,
        friction: &'a ResidualCoulombSource,
        area_m2: f64,
        compliance_m_per_n: [[f64; 3]; 3],
        config: DamageFrictionConfig,
    ) -> Result<Self, DamageFrictionError> {
        if !(1..=256).contains(&config.max_friction_iterations)
            || config.contact.max_material_evaluations < 3
        {
            return Err(DamageFrictionError::Invalid(
                "bounded friction iterations and total material budget required",
            ));
        }
        for x in [
            config.force_absolute_tolerance_n,
            config.work_absolute_tolerance_j,
        ] {
            positive(x)?;
        }
        for x in [
            config.force_relative_tolerance,
            config.work_relative_tolerance,
        ] {
            if !x.is_finite() || x <= 0.0 || x >= 1.0 {
                return Err(DamageFrictionError::Invalid(
                    "relative gates must be in (0,1)",
                ));
            }
        }
        let inner = CohesiveLoadLine::new(area_m2, compliance_m_per_n, config.contact)?;
        friction.require_cohesive_source(cohesive)?;
        let static_mu = friction.law().static_coefficient()?;
        let kinetic_mu = friction.law().zero_speed_kinetic_coefficient()?;
        let s = compliance_m_per_n;
        let scale = s.iter().flatten().map(|x| x.abs()).fold(0.0, f64::max);
        positive(scale)?;
        for i in 0..3 {
            for j in 0..3 {
                finite(s[i][j])?;
                if (s[i][j] - s[j][i]).abs() > 64.0 * f64::EPSILON * scale {
                    return Err(DamageFrictionError::Invalid(
                        "reciprocal mechanical compliance required",
                    ));
                }
            }
        }
        positive(s[1][1])?;
        positive(s[2][2])?;
        let inv_t = inverse2([[s[1][1], s[1][2]], [s[2][1], s[2][2]]])?;
        positive(s[0][0] - dot2([s[0][1], s[0][2]], mul2(inv_t, [s[1][0], s[2][0]])))?;
        let tangent_scale = positive(s[1][1].sqrt() * s[2][2].sqrt())?;
        let mut hash = DomainHasher::new("org.frankensim.fs-couple.damaged-area-rigid-friction.v1");
        hash.update(cohesive.resolved().identity().as_bytes());
        hash.update(friction.resolved().identity().as_bytes());
        for x in [
            area_m2,
            config.contact.jump_absolute_tolerance_m,
            config.contact.jump_relative_tolerance,
            config.contact.maximum_jump_increment_m,
            config.contact.maximum_force_n,
            config.contact.material_energy_tolerance_j_m2,
            config.force_absolute_tolerance_n,
            config.force_relative_tolerance,
            config.work_absolute_tolerance_j,
            config.work_relative_tolerance,
        ] {
            hash.update(&x.to_bits().to_le_bytes());
        }
        for x in s.iter().flatten() {
            hash.update(&x.to_bits().to_le_bytes());
        }
        for n in [
            config.contact.max_iterations,
            config.contact.max_material_evaluations,
            config.contact.max_backtracks,
            config.max_friction_iterations,
        ] {
            hash.update(&(n as u64).to_le_bytes());
        }
        Ok(Self {
            cohesive,
            friction,
            area: area_m2,
            compliance: s,
            config,
            identity: hash.finalize(),
            static_mu,
            kinetic_mu,
            tangent_scale,
            inner,
        })
    }

    #[must_use]
    pub const fn friction_source(&self) -> &'a ResidualCoulombSource {
        self.friction
    }
    #[must_use]
    pub const fn config(&self) -> DamageFrictionConfig {
        self.config
    }

    /// Explicit model selection works at any original damage, including zero.
    /// No bond is created and no prior dissipation is replayed on adoption.
    pub fn adopt(
        &self,
        point: CohesivePoint<'a>,
        receiving: [&MaterialStateId; 2],
        activation: ReclosureActivation,
    ) -> Result<DamageFrictionState<'a>, DamageFrictionError> {
        self.local_law().adopt(point, receiving, activation)
    }

    pub fn trial(
        &mut self,
        previous: DamageFrictionState<'a>,
        free: [f64; 3],
        dt_s: f64,
        mut keep_going: impl FnMut() -> bool,
    ) -> Result<DamageFrictionStep<'a>, DamageFrictionError> {
        self.require_state(previous)?;
        positive(dt_s)?;
        for x in free {
            finite(x)?;
        }
        let mut evaluations = 0;
        let unresisted = self.candidate(
            previous.point,
            free,
            dt_s,
            0.0,
            &mut evaluations,
            &mut keep_going,
        )?;
        let (chosen, regime, iterations) = if unresisted.pressure == 0.0 {
            (unresisted, DamageFrictionRegime::Inactive, 0)
        } else {
            let fixed = self.candidate(
                previous.point,
                free,
                dt_s,
                1.0,
                &mut evaluations,
                &mut keep_going,
            )?;
            let magnitude = norm2(fixed.friction);
            let cap = finite(self.capture_mu(previous) * fixed.pressure)?;
            if magnitude <= cap + self.force_tolerance(magnitude, cap)? {
                (fixed, DamageFrictionRegime::Sticking, 0)
            } else if self.kinetic_mu == 0.0 {
                (unresisted, DamageFrictionRegime::Sliding, 0)
            } else {
                // Endpoint signs bracket a local Coulomb root. Unlike the
                // fully failed compression-only law, cohesive softening does
                // not imply global monotonicity or uniqueness of this branch.
                let mut lo = 0.0;
                let mut hi = 1.0;
                let mut found = None;
                for iteration in 1..=self.config.max_friction_iterations {
                    let z = f64::midpoint(lo, hi);
                    let next = self.candidate(
                        previous.point,
                        free,
                        dt_s,
                        z,
                        &mut evaluations,
                        &mut keep_going,
                    )?;
                    let magnitude = norm2(next.friction);
                    let capacity = finite(self.kinetic_mu * next.pressure)?;
                    let residual = finite(magnitude - capacity)?;
                    if residual.abs() <= self.force_tolerance(magnitude, capacity)? {
                        found = Some((next, DamageFrictionRegime::Sliding, iteration));
                        break;
                    }
                    if residual > 0.0 {
                        hi = z;
                    } else {
                        lo = z;
                    }
                    if lo.next_up() >= hi {
                        break;
                    }
                }
                found.ok_or(DamageFrictionError::NotConverged)?
            }
        };
        self.spend(&mut evaluations, &mut keep_going)?;
        let audit = self.audit(previous, chosen.jump, dt_s, chosen.force, regime)?;
        let residual = core::array::from_fn(|i| {
            chosen.jump[i] - free[i] + dot3(self.compliance[i], chosen.force)
        });
        let tolerance = self.jump_tolerance(previous.point.jump_m(), free, chosen.jump)?;
        for r in residual {
            if finite(r)?.abs() > tolerance {
                return Err(DamageFrictionError::NotConverged);
            }
        }
        poll(&mut keep_going)?;
        Ok(DamageFrictionStep {
            audit,
            mean_force_n: chosen.force,
            compatibility_residual_m: residual,
            compatibility_tolerance_m: tolerance,
            iterations,
            material_evaluations: evaluations,
        })
    }

    /// Independent original-material and friction check on actual body motion.
    /// The accepted successor is solely the original timed material trial.
    pub fn audit(
        &self,
        previous: DamageFrictionState<'a>,
        jump: [f64; 3],
        dt_s: f64,
        applied: [f64; 3],
        regime: DamageFrictionRegime,
    ) -> Result<DamageFrictionAudit<'a>, DamageFrictionError> {
        self.local_law()
            .audit(previous, jump, dt_s, applied, regime)
    }

    fn candidate(
        &mut self,
        point: CohesivePoint<'a>,
        free: [f64; 3],
        dt: f64,
        z: f64,
        evaluations: &mut usize,
        keep_going: &mut impl FnMut() -> bool,
    ) -> Result<Candidate, DamageFrictionError> {
        poll(keep_going)?;
        // lambda = z/((1-z)*ell). Eliminate rigid friction f=lambda*delta_t
        // without forming infinity at stick (z=1), or inverting the full S.
        let s = self.compliance;
        let old = point.jump_m();
        let t = [[s[1][1], s[1][2]], [s[2][1], s[2][2]]];
        let row = [s[0][1], s[0][2]];
        let col = [s[1][0], s[2][0]];
        let v = [free[1] - old[1], free[2] - old[2]];
        let h = (1.0 - z) * self.tangent_scale;
        let inv = inverse2([
            [h + z * t[0][0], z * t[0][1]],
            [z * t[1][0], h + z * t[1][1]],
        ])?;
        let b = inv.map(|r| r.map(|x| z * x));
        let bcol = mul2(b, col);
        let bv = mul2(b, v);
        let mut effective = [[0.0; 3]; 3];
        effective[0][0] = positive(s[0][0] - dot2(row, bcol))?;
        for i in 0..2 {
            effective[0][i + 1] = h * (row[0] * inv[0][i] + row[1] * inv[1][i]);
            effective[i + 1][0] = h * dot2(inv[i], col);
            for j in 0..2 {
                effective[i + 1][j + 1] = h * (inv[i][0] * t[0][j] + inv[i][1] * t[1][j]);
            }
        }
        let hv = mul2(inv, v).map(|x| h * x);
        let transformed = [
            finite(free[0] - dot2(row, bv))?,
            finite(old[1] + hv[0])?,
            finite(old[2] + hv[1])?,
        ];
        let remaining = self
            .config
            .contact
            .max_material_evaluations
            .saturating_sub(*evaluations);
        if remaining <= 2 {
            return Err(DamageFrictionError::Budget);
        }
        self.inner.compliance_m_n = effective;
        self.inner.config = self.config.contact;
        self.inner.config.max_material_evaluations = remaining - 2;
        // Leave room for outer compatibility/force gates, especially near stick.
        self.inner.config.jump_absolute_tolerance_m *= 0.0625;
        self.inner.config.jump_relative_tolerance *= 0.0625;
        let solved = self.inner.trial(point, transformed, dt, &mut *keep_going)?;
        *evaluations += solved.material_evaluations;
        self.spend(evaluations, keep_going)?;
        // Immutable reevaluation shares the original exact path. This is a
        // constitutive probe, not a second damage commit or auxiliary history.
        let jump = solved.material.equilibrium.point.jump_m();
        let (_, mean, pressure) = point.trial_over_time_with_damage_weighted_compression(
            jump,
            dt,
            self.config.contact.material_energy_tolerance_j_m2,
        )?;
        let material_force = mean.map(|x| self.area * x);
        let rhs = [
            v[0] - col[0] * material_force[0] - dot2(t[0], [material_force[1], material_force[2]]),
            v[1] - col[1] * material_force[0] - dot2(t[1], [material_force[1], material_force[2]]),
        ];
        let friction = mul2(b, rhs);
        let force = [
            material_force[0],
            finite(material_force[1] + friction[0])?,
            finite(material_force[2] + friction[1])?,
        ];
        Ok(Candidate {
            jump,
            force,
            friction,
            pressure: finite(self.area * pressure)?,
        })
    }
    fn local_law(&self) -> local::LocalDamageFriction<'a> {
        local::LocalDamageFriction {
            cohesive: self.cohesive,
            friction: self.friction,
            area: self.area,
            config: self.config,
            identity: self.identity,
            static_mu: self.static_mu,
            kinetic_mu: self.kinetic_mu,
        }
    }
    fn require_state(&self, s: DamageFrictionState<'a>) -> Result<(), DamageFrictionError> {
        self.local_law().require_state(s)
    }
    fn capture_mu(&self, s: DamageFrictionState<'a>) -> f64 {
        self.local_law().capture_mu(s)
    }
    fn force_tolerance(&self, a: f64, b: f64) -> Result<f64, DamageFrictionError> {
        self.local_law().force_tolerance(a, b)
    }
    fn jump_tolerance(
        &self,
        old: [f64; 3],
        free: [f64; 3],
        jump: [f64; 3],
    ) -> Result<f64, DamageFrictionError> {
        self.local_law().jump_tolerance(old, free, jump)
    }

    fn spend(
        &self,
        n: &mut usize,
        keep_going: &mut impl FnMut() -> bool,
    ) -> Result<(), DamageFrictionError> {
        poll(keep_going)?;
        if *n >= self.config.contact.max_material_evaluations {
            return Err(DamageFrictionError::Budget);
        }
        *n += 1;
        Ok(())
    }
}

#[derive(Clone, Copy)]
struct Candidate {
    jump: [f64; 3],
    force: [f64; 3],
    friction: [f64; 2],
    pressure: f64,
}
fn inverse2(m: [[f64; 2]; 2]) -> Result<[[f64; 2]; 2], DamageFrictionError> {
    let scale = positive(m.iter().flatten().map(|x| x.abs()).fold(0.0, f64::max))?;
    let a = m[0][0] / scale;
    let b = m[0][1] / scale;
    let c = m[1][0] / scale;
    let d = m[1][1] / scale;
    let factor = positive((1.0 / scale) / positive(a * d - b * c)?)?;
    let result = [[d * factor, -b * factor], [-c * factor, a * factor]];
    for x in result.iter().flatten() {
        finite(*x)?;
    }
    Ok(result)
}
fn mul2(m: [[f64; 2]; 2], x: [f64; 2]) -> [f64; 2] {
    [dot2(m[0], x), dot2(m[1], x)]
}
fn dot2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn norm2(x: [f64; 2]) -> f64 {
    x[0].hypot(x[1])
}
fn finite(x: f64) -> Result<f64, DamageFrictionError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(DamageFrictionError::Invalid(
            "nonfinite input or arithmetic",
        ))
    }
}
fn positive(x: f64) -> Result<f64, DamageFrictionError> {
    if x.is_finite() && x > 0.0 {
        Ok(x)
    } else {
        Err(DamageFrictionError::Invalid(
            "positive finite physical input required",
        ))
    }
}
fn poll(keep_going: &mut impl FnMut() -> bool) -> Result<(), DamageFrictionError> {
    if keep_going() {
        Ok(())
    } else {
        Err(DamageFrictionError::Cancelled)
    }
}
