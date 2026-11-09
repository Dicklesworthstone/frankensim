//! Shared intrinsic Kelvin creep for simultaneous contact consumers.
//! The consumer supplies the nonlinear spring and the geometric A/L conversion.
//! No spring envelope or material parameters are invented to realize a spectrum.
use super::{CreepBranch, MAX_FELT_RATE_BRANCHES, bad, dd, finite};
use crate::{MaterialError, visco::GeneralizedMaxwell};
use fs_math::{dd::Dd, det};

/// Positive series compliance in strain/Pa units; at most eight material arms.
#[derive(Clone, Debug)]
pub struct CreepSpectrum {
    branches: Vec<CreepBranch>,
}

/// A held nonnegative stress trial. Energies are densities [J/m^3].
/// Work and stored change use the accepted floating-point checkpoint. Heat
/// integrates the original Kelvin loss law independently. The residual and
/// roundoff estimate remain visible when a tiny update lies below state resolution.
#[derive(Clone, Debug, PartialEq)]
pub struct CreepAdvance {
    pub strain: [f64; MAX_FELT_RATE_BRANCHES],
    pub work_density_j_m3: f64,
    pub stored_energy_change_j_m3: f64,
    pub dissipated_energy_j_m3: f64,
    pub balance_residual_j_m3: f64,
    /// Numerical estimate for accepted-state rounding and scalar arithmetic;
    /// this is not an interval certificate or a model/discretization allowance.
    pub roundoff_bound_j_m3: f64,
}

#[derive(Clone, Copy, Debug, Default)]
pub(super) struct PreparedBranch {
    pub(super) modulus: f64,
    pub(super) fraction: f64,
    pub(super) decay: f64,
    pub(super) loss_fraction: f64,
}

pub(super) struct BranchAdvance {
    pub(super) after: f64,
    pub(super) change: Dd,
    pub(super) actual: Dd,
    pub(super) work: Dd,
    pub(super) stored_change: Dd,
    pub(super) stored: Dd,
    pub(super) heat: Dd,
}

/// Immutable fixed-clock creep response. Successful warm calls use only bounded
/// scalar arithmetic and fixed-size checkpoint storage; no allocation or exp.
#[derive(Clone, Debug)]
pub struct PreparedCreep {
    branches: [PreparedBranch; MAX_FELT_RATE_BRANCHES],
    count: usize,
    duration_s: f64,
    compliance: f64,
}

pub(super) fn validate_branches(branches: &[CreepBranch]) -> Result<(), MaterialError> {
    if branches.len() > MAX_FELT_RATE_BRANCHES
        || branches.iter().any(|b| {
            !b.modulus_pa.is_finite()
                || b.modulus_pa <= 0.0
                || !b.viscosity_pa_s.is_finite()
                || b.viscosity_pa_s <= 0.0
        })
    {
        return Err(MaterialError::Parameters {
            what: "felt creep needs at most eight positive finite modulus/viscosity pairs".into(),
        });
    }
    Ok(())
}

pub(super) fn prepare_branches(
    source: &[CreepBranch],
    duration_s: f64,
) -> Result<([PreparedBranch; MAX_FELT_RATE_BRANCHES], f64), MaterialError> {
    validate_branches(source)?;
    if !duration_s.is_finite() || duration_s < 0.0 {
        return Err(bad("finite nonnegative duration required"));
    }
    let mut branches = [PreparedBranch::default(); MAX_FELT_RATE_BRANCHES];
    let mut compliance = Dd::ZERO;
    for (i, branch) in source.iter().enumerate() {
        let tau = branch.viscosity_pa_s / branch.modulus_pa;
        if !tau.is_finite() || tau <= 0.0 {
            return Err(bad("unrepresentable retardation time"));
        }
        let exposure = duration_s / tau;
        if duration_s > 0.0 && exposure == 0.0 {
            return Err(bad("unresolved positive creep exposure"));
        }
        let fraction = -det::expm1(-exposure);
        let decay = det::exp(-exposure);
        let loss_fraction = -det::expm1(-2.0 * exposure);
        let response = fraction / branch.modulus_pa;
        if !response.is_finite() || (fraction > 0.0 && response == 0.0) {
            return Err(bad("unrepresentable positive creep response"));
        }
        branches[i] = PreparedBranch {
            modulus: branch.modulus_pa,
            fraction,
            decay,
            loss_fraction,
        };
        compliance = compliance + dd(response);
    }
    Ok((
        branches,
        finite(compliance.to_f64(), "creep compliance overflow")?,
    ))
}

pub(super) fn free_strain(
    branches: &[PreparedBranch],
    state: &[f64; MAX_FELT_RATE_BRANCHES],
) -> Result<f64, MaterialError> {
    finite(
        branches
            .iter()
            .zip(state)
            .fold(Dd::ZERO, |sum, (b, z)| {
                if b.fraction <= 0.5 {
                    sum + dd(*z) - dd(b.fraction) * dd(*z)
                } else {
                    sum + dd(b.decay) * dd(*z)
                }
            })
            .to_f64(),
        "free creep strain overflow",
    )
}

impl PreparedBranch {
    pub(super) fn advance(&self, strain: f64, stress: f64) -> Result<BranchAdvance, MaterialError> {
        let before = dd(strain);
        let target = dd(stress) / dd(self.modulus);
        let difference = target - before;
        // Preserve both a small increment and a small surviving recovery tail.
        let (change, endpoint) = if self.fraction <= 0.5 {
            let change = dd(self.fraction) * difference;
            (change, before + change)
        } else {
            let remainder = dd(self.decay) * difference;
            (difference - remainder, target - remainder)
        };
        let after = finite(endpoint.to_f64(), "creep strain overflow")?;
        if after < 0.0 {
            return Err(bad("negative physical creep strain"));
        }
        let actual = dd(after) - before;
        let heat = dd(0.5)
            * dd(self.modulus)
            * dd(self.loss_fraction)
            * (before - target)
            * (before - target);
        let work = dd(stress) * actual;
        let stored_change = dd(0.5) * dd(self.modulus) * (dd(after) + before) * actual;
        let stored = dd(0.5) * dd(self.modulus) * dd(after) * dd(after);
        for (value, what) in [
            (change, "creep increment overflow"),
            (heat, "creep heat overflow"),
            (work, "creep work overflow"),
            (stored_change, "creep stored change overflow"),
            (stored, "creep storage overflow"),
        ] {
            finite(value.to_f64(), what)?;
        }
        if heat.to_f64() < 0.0 {
            return Err(bad("negative physical Kelvin loss"));
        }
        Ok(BranchAdvance {
            after,
            change,
            actual,
            work,
            stored_change,
            stored,
            heat,
        })
    }
}

impl CreepSpectrum {
    pub fn new(branches: Vec<CreepBranch>) -> Result<Self, MaterialError> {
        validate_branches(&branches)?;
        Ok(Self { branches })
    }
    /// Original source compliance, excluding its instantaneous spring.
    /// Exact equal source times coalesce; unresolved weak modes refuse.
    pub fn from_prony(card: &GeneralizedMaxwell) -> Result<Self, MaterialError> {
        Self::new(super::prony::realize(card)?)
    }
    pub fn branches(&self) -> &[CreepBranch] {
        &self.branches
    }
    pub fn prepare(&self, duration_s: f64) -> Result<PreparedCreep, MaterialError> {
        let (branches, compliance) = prepare_branches(&self.branches, duration_s)?;
        Ok(PreparedCreep {
            branches,
            count: self.branches.len(),
            duration_s,
            compliance,
        })
    }
}

impl PreparedCreep {
    pub fn duration_s(&self) -> f64 {
        self.duration_s
    }
    pub fn branch_count(&self) -> usize {
        self.count
    }
    pub fn compliance_pa_inverse(&self) -> f64 {
        self.compliance
    }

    fn admit(&self, state: &[f64; MAX_FELT_RATE_BRANCHES]) -> Result<(), MaterialError> {
        if state.iter().any(|z| !z.is_finite() || *z < 0.0)
            || state[self.count..].iter().any(|z| *z != 0.0)
        {
            return Err(bad(
                "creep history needs finite nonnegative active strains and zero inactive entries",
            ));
        }
        Ok(())
    }
    pub fn strain(&self, state: &[f64; MAX_FELT_RATE_BRANCHES]) -> Result<f64, MaterialError> {
        self.admit(state)?;
        finite(
            state[..self.count]
                .iter()
                .fold(Dd::ZERO, |sum, z| sum + dd(*z))
                .to_f64(),
            "total creep strain overflow",
        )
    }
    pub fn free_strain(&self, state: &[f64; MAX_FELT_RATE_BRANCHES]) -> Result<f64, MaterialError> {
        self.admit(state)?;
        free_strain(&self.branches[..self.count], state)
    }
    pub fn stored_energy_density(
        &self,
        state: &[f64; MAX_FELT_RATE_BRANCHES],
    ) -> Result<f64, MaterialError> {
        self.admit(state)?;
        finite(
            self.branches[..self.count]
                .iter()
                .zip(state)
                .fold(Dd::ZERO, |sum, (b, z)| {
                    sum + dd(0.5) * dd(b.modulus) * dd(*z) * dd(*z)
                })
                .to_f64(),
            "creep storage overflow",
        )
    }
    pub fn advance(
        &self,
        state: &[f64; MAX_FELT_RATE_BRANCHES],
        stress_pa: f64,
    ) -> Result<CreepAdvance, MaterialError> {
        self.admit(state)?;
        if !stress_pa.is_finite() || stress_pa < 0.0 {
            return Err(bad("finite compressive creep stress required"));
        }
        let mut strain = [0.0; MAX_FELT_RATE_BRANCHES];
        let mut work = Dd::ZERO;
        let mut stored_change = Dd::ZERO;
        let mut heat = Dd::ZERO;
        let mut roundoff = Dd::ZERO;
        for (i, b) in self.branches[..self.count].iter().enumerate() {
            let step = b.advance(state[i], stress_pa)?;
            strain[i] = step.after;
            work = work + step.work;
            stored_change = stored_change + step.stored_change;
            heat = heat + step.heat;
            // Difference between the analytic increment and the accepted
            // checkpoint. Its effect on work-minus-storage is bounded by
            // (|stress|+E|z1|)|delta| + E delta^2/2.
            let error = dd((step.actual - step.change).to_f64().abs());
            let state_rounding = dd(stress_pa) * error
                + dd(b.modulus) * dd(step.after.abs()) * error
                + dd(0.5) * dd(b.modulus) * error * error;
            let scale = step
                .work
                .to_f64()
                .abs()
                .max(step.stored_change.to_f64().abs())
                .max(step.heat.to_f64().abs());
            roundoff = roundoff + state_rounding + dd(128.0 * f64::EPSILON * scale);
        }
        let residual = work - stored_change - heat;
        let bound = finite(roundoff.to_f64(), "creep roundoff estimate overflow")?;
        let defect = finite(residual.to_f64(), "creep work balance overflow")?;
        if defect.abs() > bound {
            return Err(bad(
                "physical Kelvin work and heat exceed the state-roundoff estimate",
            ));
        }
        Ok(CreepAdvance {
            strain,
            work_density_j_m3: finite(work.to_f64(), "creep work overflow")?,
            stored_energy_change_j_m3: finite(
                stored_change.to_f64(),
                "creep stored change overflow",
            )?,
            dissipated_energy_j_m3: finite(heat.to_f64(), "creep heat overflow")?,
            balance_residual_j_m3: defect,
            roundoff_bound_j_m3: bound,
        })
    }
}
