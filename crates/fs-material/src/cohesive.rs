//! Bilinear cohesive failure in a fixed local interface frame.
//!
//! The jump is [normal opening, shear 1, shear 2] in metres. Positive normal
//! opening and both shear components drive a common irreversible damage variable.
//! By default all three modes share the declared stiffness, peak traction and
//! fracture energy. `new_mixed_mode` explicitly selects separate opening/shear
//! data, quadratic initiation and frozen-onset BK toughness through `mixed_mode`.
//! Both paths retain the same immutable point, work and structural solver APIs.
//! `with_viscosity` explicitly adds physical bond dashpots and requires timed steps.
//!
//! Unloading/reloading is elastic at the retained damaged stiffness. Compression
//! has a separate, undamaged penalty; after complete failure only that normal
//! penalty remains (frictionless reclosure, not exact nonpenetration). No healing,
//! fatigue, automatic bond creation, temperature dependence, calibrated
//! adhesion, evolving contact geometry or finite-rotation objectivity is claimed.
//!
//! Reference formulation: Camanho and Davila (2002), as described by the Abaqus
//! cohesive-contact damage-evolution documentation. Fracture energy is the area
//! under the tensile/shear traction-separation envelope. No foreign code is used.

use fs_math::det;

/// Ordered interface data bound to the existing cohesive constitutive law.
pub mod source;

/// Anisotropic initiation and source-defined mode-dependent fracture energy.
pub mod mixed_mode;

/// Explicit physical bond viscosity and timed constitutive steps.
pub mod rate;
pub use rate::{CohesiveRateStep, CohesiveViscosity};

/// A refusal never mutates an accepted point.
#[derive(Debug, Clone, PartialEq)]
pub enum CohesiveError {
    /// Invalid data or nonrepresentable arithmetic.
    Invalid(&'static str),
    /// Separation exceeds the caller's declared local-model range.
    SeparationLimit,
    /// Actual work/storage/fracture mismatch exceeds the supplied J/m2 gate.
    EnergyDefect(f64),
}
impl core::fmt::Display for CohesiveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "cohesive interface refused: {self:?}")
    }
}
impl std::error::Error for CohesiveError {}

/// An immutable, explicitly supplied cohesive law. There are no material defaults.
#[derive(Debug, Clone, PartialEq)]
pub struct BilinearCohesive {
    stiffness: f64,
    strength: f64,
    fracture_energy: f64,
    compression_stiffness: f64,
    maximum_separation: f64,
    onset: f64,
    failure: f64,
    mixed_mode: Option<mixed_mode::MixedModeParameters>,
    viscosity: Option<CohesiveViscosity>,
}
impl BilinearCohesive {
    /// Coherent SI inputs: Pa/m, Pa, J/m2, Pa/m and m, respectively.
    /// Requires positive finite data and `2 Gc / strength > strength / stiffness`.
    /// The separation limit must include complete failure. Supplying parameters
    /// is not evidence that they are calibrated for a particular printed bond.
    pub fn new(stiffness_pa_m: f64, strength_pa: f64, fracture_energy_j_m2: f64,
        compression_stiffness_pa_m: f64, maximum_separation_m: f64,
    ) -> Result<Self, CohesiveError> {
        for value in [stiffness_pa_m, strength_pa, fracture_energy_j_m2,
            compression_stiffness_pa_m, maximum_separation_m] { positive(value)?; }
        let onset = positive(strength_pa / stiffness_pa_m)?;
        let failure = positive(2.0 * (fracture_energy_j_m2 / strength_pa))?;
        if failure <= onset || failure > maximum_separation_m {
            return Err(CohesiveError::Invalid("require onset < failure <= separation limit"));
        }
        // Refuse data whose bounded storage or tractions cannot be represented.
        positive((stiffness_pa_m * maximum_separation_m) * maximum_separation_m)?;
        positive((compression_stiffness_pa_m * maximum_separation_m) * maximum_separation_m)?;
        positive(strength_pa / (failure - onset))?;
        Ok(Self { stiffness: stiffness_pa_m, strength: strength_pa,
            fracture_energy: fracture_energy_j_m2, compression_stiffness: compression_stiffness_pa_m,
            maximum_separation: maximum_separation_m, onset, failure, mixed_mode: None, viscosity: None })
    }

    #[must_use]
    pub fn onset_separation_m(&self) -> f64 { self.onset }
    #[must_use]
    pub fn failure_separation_m(&self) -> f64 { self.failure }
    #[must_use]
    pub fn maximum_separation_m(&self) -> f64 { self.maximum_separation }
    #[must_use]
    pub fn strength_pa(&self) -> f64 { self.strength }
    #[must_use]
    pub fn fracture_energy_j_m2(&self) -> f64 { self.fracture_energy }
    /// Magnitude of the pure-opening negative softening tangent, Pa/m.
    /// For a mixed-mode law this is not a bound on every coupled tangent.
    #[must_use]
    pub fn softening_stiffness_pa_m(&self) -> f64 { self.strength / (self.failure-self.onset) }

    /// Explicitly create an initially bonded, stress-free point at zero jump.
    /// Geometry/contact/healing owners, not this constructor, establish whether
    /// an actual interface is bonded and which physical area it represents.
    #[must_use]
    pub fn bonded_point(&self) -> CohesivePoint<'_> {
        CohesivePoint { law: self, jump: [0.0; 3], maximum: 0.0, mixed_history: None }
    }

    fn intact(&self, maximum: f64) -> f64 {
        if maximum <= self.onset { 1.0 }
        else if maximum >= self.failure { 0.0 }
        else { (self.onset/maximum) * ((self.failure-maximum)/(self.failure-self.onset)) }
    }
    fn fracture(&self, maximum: f64) -> f64 {
        if maximum <= self.onset { 0.0 }
        else if maximum >= self.failure { self.fracture_energy }
        else { self.fracture_energy * ((maximum-self.onset)/(self.failure-self.onset)) }
    }
    // Exact radial work on the loading envelope; stable difference of integrals.
    fn envelope_work(&self, start: f64, end: f64) -> f64 {
        let a = start.min(self.onset);
        let b = end.min(self.onset);
        let elastic = 0.5 * self.stiffness * (b-a) * (b+a);
        let a = start.max(self.onset).min(self.failure);
        let b = end.max(self.onset).min(self.failure);
        let softening = self.strength * ((b-a)/(self.failure-self.onset))
            * (0.5*(self.failure-a) + 0.5*(self.failure-b));
        elastic + softening
    }
}

/// Immutable accepted history, bound to its original law and fixed local frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CohesivePoint<'a> {
    law: &'a BilinearCohesive,
    jump: [f64; 3],
    maximum: f64,
    mixed_history: Option<mixed_mode::MixedModeHistory>,
}

/// Incremental energy per original bonded area, J/m2.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CohesiveEnergy {
    /// Analytically integrated traction work on the declared straight jump path.
    pub work: f64,
    /// Change in recoverable cohesive plus compression-penalty energy.
    pub stored_change: f64,
    /// Irreversible fracture dissipation, not numerical damping.
    pub dissipated: f64,
    /// Work minus storage change minus fracture dissipation.
    pub defect: f64,
}

/// A trial can be discarded without committing any damage.
#[derive(Debug, Clone, PartialEq)]
pub struct CohesiveStep<'a> {
    pub point: CohesivePoint<'a>,
    /// Work-conjugate traction, Pa. Resisting force on the positive side is -A*t.
    pub traction_pa: [f64; 3],
    /// Exact local d(traction)/d(jump), Pa/m, at smooth points. None at active
    /// loading/unloading, failure-initiation or tension/compression switches.
    /// Softening tangents can be indefinite: do not pass them to an SPD solver.
    pub tangent_pa_m: Option<[[f64; 3]; 3]>,
    /// Selected piecewise tangent for a semismooth structural iteration. At a
    /// switch it selects tension and retained-history unloading; it is NOT a
    /// classical derivative there. Trial iterates never commit damage.
    pub branch_tangent_pa_m: [[f64; 3]; 3],
    pub energy: CohesiveEnergy,
}

impl<'a> CohesivePoint<'a> {
    #[must_use]
    pub fn law(&self) -> &'a BilinearCohesive { self.law }
    #[must_use]
    pub fn jump_m(&self) -> [f64; 3] { self.jump }
    #[must_use]
    pub fn damage(&self) -> f64 {
        if self.law.mixed_mode.is_some() { return 1.0-self.mixed_integrity(); }
        1.0-self.law.intact(self.maximum)
    }
    #[must_use]
    pub fn failed(&self) -> bool {
        if self.law.mixed_mode.is_some() {
            return self.mixed_history.is_some_and(|h| self.maximum >= h.failure_sqrt_energy);
        }
        self.maximum >= self.law.failure
    }
    #[must_use]
    pub fn dissipated_j_m2(&self) -> f64 {
        if self.law.mixed_mode.is_some() {
            return self.mixed_history.map_or(0.0, |h| h.fracture(self.maximum));
        }
        self.law.fracture(self.maximum)
    }
    #[must_use]
    pub fn stored_energy_j_m2(&self) -> f64 {
        if self.law.mixed_mode.is_some() { return self.mixed_stored_energy(); }
        let active = active(self.jump);
        let stiffness = self.law.stiffness*self.law.intact(self.maximum);
        let cohesive = active.iter().map(|&x| (0.5*stiffness*x)*x).sum::<f64>();
        let compression = self.jump[0].min(0.0);
        cohesive + (0.5*self.law.compression_stiffness*compression)*compression
    }
    #[must_use]
    pub fn traction_pa(&self) -> [f64; 3] {
        if self.law.mixed_mode.is_some() { return self.mixed_traction(); }
        let stiffness = self.law.stiffness*self.law.intact(self.maximum);
        [if self.jump[0] < 0.0 { self.law.compression_stiffness*self.jump[0] }
            else { stiffness*self.jump[0] }, stiffness*self.jump[1], stiffness*self.jump[2]]
    }

    /// Advance along a straight local separation path, with an absolute J/m2
    /// numerical energy gate. Work is integrated independently from the change
    /// in fracture dissipation. Endpoint-only data cannot encode unreported
    /// excursions: split a physical path at each reversal/contact-frame change.
    /// This fixed-size operation is bounded O(1); batch owners poll cancellation.
    pub fn trial(&self, jump_m: [f64; 3], energy_tolerance_j_m2: f64)
        -> Result<CohesiveStep<'a>, CohesiveError> {
        if self.law.viscosity.is_some() {
            return Err(CohesiveError::Invalid("physical bond viscosity requires trial_over_time and an explicit time interval"));
        }
        self.trial_equilibrium(jump_m, energy_tolerance_j_m2)
    }

    // Shared rate-independent owner. Timed updates add physical dashpot traction
    // without replacing its irreversible history, storage or fracture work.
    fn trial_equilibrium(&self, jump_m: [f64; 3], energy_tolerance_j_m2: f64)
        -> Result<CohesiveStep<'a>, CohesiveError> {
        positive(energy_tolerance_j_m2)?;
        if norm(jump_m)? > self.law.maximum_separation { return Err(CohesiveError::SeparationLimit); }
        if self.law.mixed_mode.is_some() { return self.trial_mixed(jump_m, energy_tolerance_j_m2); }
        let radius = norm(active(jump_m))?;
        let old_radius = norm(active(self.jump))?;
        let maximum = self.maximum.max(radius).min(self.law.failure);
        let point = Self { law: self.law, jump: jump_m, maximum, mixed_history: None };
        let stiffness_old = self.law.stiffness*self.law.intact(self.maximum);
        // The effective radius is convex on a straight jump segment, so it
        // cannot exceed both endpoints. Below the old maximum, constant-damage
        // work telescopes even if the path passes through compression or turns.
        let work = if radius > self.maximum {
            0.5*stiffness_old*(self.maximum-old_radius)*(self.maximum+old_radius)
                + self.law.envelope_work(self.maximum, radius)
        } else {
            0.5*stiffness_old*(radius-old_radius)*(radius+old_radius)
        };
        let old_compression = self.jump[0].min(0.0);
        let compression = jump_m[0].min(0.0);
        let work = finite(work + 0.5*self.law.compression_stiffness
            *(compression-old_compression)*(compression+old_compression))?;
        let stored_change = finite(point.stored_energy_j_m2()-self.stored_energy_j_m2())?;
        let dissipated = finite(point.dissipated_j_m2()-self.dissipated_j_m2())?;
        let defect = finite(work-stored_change-dissipated)?;
        if defect.abs() > energy_tolerance_j_m2 { return Err(CohesiveError::EnergyDefect(defect)); }
        let traction_pa = point.traction_pa();
        for value in traction_pa { finite(value)?; }
        let (branch_tangent_pa_m, classical) = self.tangent(&point, radius)?;
        let tangent_pa_m = classical.then_some(branch_tangent_pa_m);
        Ok(CohesiveStep { point, traction_pa, tangent_pa_m, branch_tangent_pa_m,
            energy: CohesiveEnergy { work, stored_change, dissipated, defect } })
    }

    fn tangent(&self, next: &Self, radius: f64) -> Result<([[f64; 3]; 3], bool), CohesiveError> {
        let law = self.law;
        let classical = !(next.jump[0] == 0.0
            || (self.maximum < law.failure && radius >= self.maximum
                && (radius == law.onset || radius == law.failure
                    || (radius == self.maximum && radius > law.onset))));
        let stiffness = law.stiffness*law.intact(next.maximum);
        let mut tangent = [[0.0; 3]; 3];
        tangent[0][0] = if next.jump[0] < 0.0 { law.compression_stiffness } else { stiffness };
        tangent[1][1] = stiffness;
        tangent[2][2] = stiffness;
        if radius > self.maximum && radius > law.onset && radius < law.failure {
            let unit = active(next.jump).map(|value| value/radius);
            let change = -law.softening_stiffness_pa_m()-stiffness;
            for i in 0..3 { for j in 0..3 { tangent[i][j] += change*unit[i]*unit[j]; } }
        }
        for row in tangent { for value in row { finite(value)?; } }
        Ok((tangent, classical))
    }
}

fn active(jump: [f64; 3]) -> [f64; 3] { [jump[0].max(0.0), jump[1], jump[2]] }
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 { a[0]*b[0]+a[1]*b[1]+a[2]*b[2] }
fn norm(value: [f64; 3]) -> Result<f64, CohesiveError> {
    for x in value { finite(x)?; }
    let scale = value.iter().fold(0.0_f64, |a, b| a.max(b.abs()));
    if scale == 0.0 { return Ok(0.0); }
    let v = value.map(|x| x/scale);
    finite(scale*det::sqrt(dot(v,v)))
}
fn finite(value: f64) -> Result<f64, CohesiveError> {
    if value.is_finite() { Ok(value) } else { Err(CohesiveError::Invalid("nonfinite arithmetic/input")) }
}
fn positive(value: f64) -> Result<f64, CohesiveError> {
    if value.is_finite() && value > 0.0 { Ok(value) }
    else { Err(CohesiveError::Invalid("requires positive finite data")) }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn law() -> BilinearCohesive { BilinearCohesive::new(1e9, 1e6, 2e3, 2e9, 0.01).unwrap() }
    fn close(a: f64, b: f64) { assert!((a-b).abs() < 1e-9*(1.0+a.abs()+b.abs()), "{a} != {b}"); }

    #[test]
    fn pure_and_mixed_proportional_failure_spends_exact_declared_fracture_energy() {
        let law = law();
        for direction in [[1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,-1.0], [0.6,0.0,0.8]] {
            let mut point = law.bonded_point();
            let mut work = 0.0;
            for i in 1..=100 {
                let r = law.failure_separation_m()*f64::from(i)/80.0;
                let step = point.trial(direction.map(|v| v*r), 1e-8).unwrap();
                assert!(step.point.damage() >= point.damage());
                assert!(step.energy.dissipated >= 0.0);
                work += step.energy.work;
                point = step.point;
            }
            assert!(point.failed());
            close(work, 2000.0);
            close(point.stored_energy_j_m2(), 0.0);
            close(point.dissipated_j_m2(), 2000.0);
        }
    }

    #[test]
    fn unloading_does_not_heal_and_complete_failure_retains_compression_contact() {
        let law = law();
        let partial = law.bonded_point().trial([0.002,0.0,0.0], 1e-8).unwrap().point;
        let unload = partial.trial([0.0005,0.0,0.0], 1e-8).unwrap();
        close(unload.point.damage(), partial.damage());
        close(unload.energy.dissipated, 0.0);
        close(unload.traction_pa[0], partial.traction_pa()[0]/4.0);
        let reload = unload.point.trial([0.002,0.0,0.0], 1e-8).unwrap();
        close(unload.energy.work+reload.energy.work, 0.0);
        assert!(reload.tangent_pa_m.is_none());
        let broken = partial.trial([0.006,0.0,0.0], 1e-8).unwrap().point;
        let closed = broken.trial([-0.001,0.002,-0.001], 1e-8).unwrap();
        assert_eq!(closed.traction_pa, [-2e6,0.0,0.0]);
        close(closed.energy.dissipated, 0.0);
        close(closed.point.damage(), 1.0);
    }

    #[test]
    fn compression_alone_never_damages() {
        let law = law();
        let compression = law.bonded_point().trial([-0.006,0.0,0.0], 1e-8).unwrap();
        assert_eq!(compression.point.damage(), 0.0);
        assert_eq!(compression.energy.dissipated, 0.0);
        close(compression.energy.work, 36000.0);
    }

    #[test]
    fn consistent_tangent_matches_independent_endpoint_perturbations() {
        let law = law();
        let virgin = law.bonded_point();
        let damaged = virgin.trial([0.0025,0.0,0.0], 1e-8).unwrap().point;
        for (point, jump) in [(virgin,[0.0003,0.0002,-0.0001]),
            (virgin,[0.0016,0.0008,0.0004]), (virgin,[-0.001,0.0018,0.0005]),
            (damaged,[0.001,0.0004,-0.0002])] {
            let tangent = point.trial(jump,1e-8).unwrap().tangent_pa_m.unwrap();
            for j in 0..3 {
                let mut plus = jump; let mut minus = jump;
                plus[j] += 1e-8; minus[j] -= 1e-8;
                let a = point.trial(plus,1e-8).unwrap().traction_pa;
                let b = point.trial(minus,1e-8).unwrap().traction_pa;
                for i in 0..3 {
                    let fd = (a[i]-b[i])/2e-8;
                    assert!((fd-tangent[i][j]).abs() < 1e-7*1e9);
                    close(tangent[i][j],tangent[j][i]);
                }
            }
        }
    }

    #[test]
    fn nonproportional_straight_path_work_matches_traction_quadrature() {
        let law = law();
        let point = law.bonded_point().trial([0.0012,0.0009,0.0],1e-8).unwrap().point;
        let end = [-0.0004,0.0028,0.0006];
        let delta = std::array::from_fn(|i| end[i]-point.jump[i]);
        let mut work = 0.0;
        let count = 20_000;
        for i in 0..count {
            let s = (f64::from(i)+0.5)/f64::from(count);
            let jump = std::array::from_fn(|j| point.jump[j]+s*delta[j]);
            let traction = point.trial(jump,1e-8).unwrap().traction_pa;
            work += dot(traction,delta)/f64::from(count);
        }
        let result = point.trial(end,1e-8).unwrap();
        assert!((work-result.energy.work).abs() < 1e-4);
    }

    #[test]
    fn rejected_trials_and_restart_preserve_accepted_damage() {
        let law = law();
        let point = law.bonded_point().trial([0.002,0.0001,0.0],1e-8).unwrap().point;
        let saved = point;
        assert!(point.trial([0.02,0.0,0.0],1e-8).is_err());
        assert!(point.trial([f64::NAN,0.0,0.0],1e-8).is_err());
        assert!(point.trial([0.003,0.0,0.0],0.0).is_err());
        assert_eq!(point,saved);
        assert_eq!(point.trial([0.003,0.0002,0.0],1e-8), saved.trial([0.003,0.0002,0.0],1e-8));
        assert!(BilinearCohesive::new(1e9,1e6,100.0,1e9,0.01).is_err());
        assert!(BilinearCohesive::new(1e9,1e6,2000.0,1e9,0.001).is_err());
    }
}
