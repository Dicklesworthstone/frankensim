//! Fixed-plane brittle crack bands in the existing small-strain material path.
//!
//! One declared material plane and width `h` split isotropic elastic energy
//! exactly into a positive plane-stress remainder and three transverse modes.
//! With plane normal `n`, tangents `s,t`, `r = nu/(1-nu)`, and `M = lambda+2mu`,
//! the work-conjugate mode tensors are
//!
//! ```text
//! A0 = n*n + r*(s*s+t*t), A1 = n*s+s*n, A2 = n*t+t*n;
//! q_i = A_i : strain, jump_i = h*q_i;
//! Kn = M/h, Kt = mu/h;
//! stress = C_parallel*strain + sum_i A_i*traction_i.
//! ```
//!
//! The ORIGINAL mixed-mode cohesive owner supplies onset, irreversible damage,
//! BK fracture work, unloading, compression and exact straight-path work. There
//! is no second fracture evolution law. Its area energy divided by `h` is the
//! volumetric energy of this band. Intact response is the original isotropic
//! elasticity; failed opening/shear modes lose transfer while in-plane elastic
//! response and undamaged normal compression remain. Free-sided compression
//! does not become tensile damage merely through Poisson lateral expansion.
//!
//! `h` is an explicitly declared physical/numerical crack-band width, not an
//! automatically inferred mesh size. The normal band coordinate includes the
//! Poisson correction: it is NOT a measured crack-face gap. Constitutive
//! reclosure is not a geometric contact or nonpenetration algorithm. The plane
//! does not rotate or nucleate automatically; no multiple-crack interaction,
//! arbitrary-mesh objectivity, fragment topology, remeshing/state transfer,
//! friction, crushing, rate evolution or experimental qualification is claimed.
//! Ordinary strain updates hold coefficients fixed; [`temperature`] transfers
//! damage and energy through explicit thermal coefficient changes at fixed
//! Poisson ratio and geometry. The caller resolves load reversals/crack events.
//!
//! Crack-band energy scaling follows Bazant and Oh (1983),
//! doi:10.1007/BF02486267. The fixed-plane elastic projection here is an explicit
//! composition with the native Camanho--Davila/BK cohesive owner, not a claim to
//! reproduce a particular calibrated concrete or published constitutive law.
//!
//! This fixed-size constitutive update allocates no memory. Batch/element owners
//! remain responsible for cancellation and accepted-state publication.

use crate::cohesive::mixed_mode::MixedModeParameters;
use crate::cohesive::{BilinearCohesive, CohesiveError, CohesivePoint, CohesiveStep};
use crate::elastic::IsotropicElastic;
use crate::tensor::contract;
use crate::{MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::{Ambition, ModelCard, ValidityDomain};
use fs_math::det;

/// Original bulk-material coefficients bound to this native constitutive law.
pub mod source;

/// History-preserving changes of source-defined thermal coefficients.
pub mod temperature;

/// Bounded tensile onset discovery and elastic-to-crack-band handoff.
pub mod initiation;

/// Exact path-mean stress and its native derivative for work-conjugate dynamics.
pub mod path;

/// Explicit pure-mode strengths and fracture work. No values are inferred from
/// a material name, hardness, flexural modulus of rupture or yield stress.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleCrackBandParameters {
    /// Pure normal peak traction on the declared plane, Pa.
    pub tensile_strength_pa: f64,
    /// Pure shear peak traction, Pa; isotropic within the declared plane.
    pub shear_strength_pa: f64,
    /// Complete mode-I work per fracture area, J/m2.
    pub mode_i_fracture_energy_j_m2: f64,
    /// Complete mode-II work per fracture area, J/m2.
    pub mode_ii_fracture_energy_j_m2: f64,
    /// Benzeggagh--Kenane exponent, at least one. Mode mix freezes at onset.
    pub bk_exponent: f64,
}

/// Explicit receiving geometry and numerical/model domain.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleCrackBandSettings {
    /// Positive band width, metres. Energy per volume is fracture work / width.
    pub band_width_m: f64,
    /// Finite nonzero plane direction, normalized once during construction.
    pub normal: [f64; 3],
    /// Maximum Frobenius norm of total tensor strain, in (0, 0.05], also
    /// restricted by the supplied elastic law's declared strain domain.
    pub maximum_strain: f64,
    /// Absolute gate on incremental work minus storage and loss, J/m3.
    pub energy_tolerance_j_m3: f64,
}

/// Immutable prepared material and geometric projection. Histories borrow this
/// particular owner, so coefficients or band width cannot change under damage.
#[derive(Clone, Debug, PartialEq)]
pub struct BrittleCrackBand {
    elastic: IsotropicElastic,
    parameters: BrittleCrackBandParameters,
    settings: BrittleCrackBandSettings,
    axes: [[f64; 3]; 3],
    modes: [Voigt; 3],
    parallel_modes: [Voigt; 3],
    parallel_moduli: [f64; 3],
    parallel_tangent: Tangent6,
    cohesive: BilinearCohesive,
    area_energy_tolerance: f64,
}

/// A checked refusal never modifies accepted material history.
#[derive(Clone, Debug, PartialEq)]
pub enum BrittleError {
    Invalid(&'static str),
    /// Another prepared owner, including an identical clone, owns this history.
    StateMismatch,
    /// Frobenius tensor-strain norm exceeds the explicitly declared domain.
    StrainLimit { norm: f64, maximum: f64 },
    Cohesive(CohesiveError),
    /// Reported volumetric work/storage/loss do not meet the absolute gate.
    EnergyDefect(f64),
}
impl core::fmt::Display for BrittleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "brittle crack band refused: {self:?}")
    }
}
impl std::error::Error for BrittleError {}
impl From<CohesiveError> for BrittleError {
    fn from(error: CohesiveError) -> Self { Self::Cohesive(error) }
}

/// Accepted strain and the ORIGINAL cohesive history, borrowing its owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleState<'a> {
    owner: &'a BrittleCrackBand,
    strain: Voigt,
    point: CohesivePoint<'a>,
}

/// Incremental energies per band volume, J/m3. Fracture loss is not implicitly
/// converted into heat, numerical damping or kinetic energy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleEnergy {
    pub work: f64,
    pub stored_change: f64,
    pub dissipated: f64,
    pub defect: f64,
}

/// A disposable trial with no mutation of the supplied accepted history.
#[derive(Clone, Debug, PartialEq)]
pub struct BrittleStep<'a> {
    pub stress_pa: Voigt,
    /// Classical derivative at smooth points; absent at damage/contact switches.
    pub tangent_pa: Option<Tangent6>,
    /// Original native piecewise selection for semismooth Newton iterations.
    /// Mixed-mode onset can make this tangent nonsymmetric or indefinite.
    pub branch_tangent_pa: Tangent6,
    pub state: BrittleState<'a>,
    pub energy: BrittleEnergy,
}

/// Borrowed law for existing tetrahedron, hexahedron and reduced solid owners.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleSmallStrain<'a> {
    owner: &'a BrittleCrackBand,
}

impl BrittleCrackBand {
    /// Prepare one fixed plane and width. A width whose elastic energy at
    /// initiation exhausts the supplied fracture work refuses: the native law
    /// requires a resolved positive softening interval in EVERY mode mix.
    pub fn new(
        elastic: IsotropicElastic,
        parameters: BrittleCrackBandParameters,
        settings: BrittleCrackBandSettings,
    ) -> Result<Self, BrittleError> {
        positive(elastic.youngs)?;
        positive(elastic.strain_limit)?;
        if !(elastic.poisson > -1.0 && elastic.poisson < 0.5) {
            return Err(BrittleError::Invalid("Poisson ratio must lie in (-1, 0.5)"));
        }
        positive(settings.band_width_m)?;
        positive(settings.maximum_strain)?;
        positive(settings.energy_tolerance_j_m3)?;
        if settings.maximum_strain > 0.05 || settings.maximum_strain > elastic.strain_limit {
            return Err(BrittleError::Invalid(
                "strain domain must not exceed 0.05 or the elastic source domain",
            ));
        }
        let normal = normalize(settings.normal)?;
        let mut least = 0;
        for i in 1..3 {
            if normal[i].abs() < normal[least].abs() { least = i; }
        }
        let tangent = normalize(core::array::from_fn(|i| {
            if i == least { 1.0 - normal[least]*normal[i] }
            else { -normal[least]*normal[i] }
        }))?;
        let second = normalize(cross(normal, tangent))?;
        let axes = [normal, tangent, second];
        let nn = dyad(normal, normal);
        let ss = dyad(tangent, tangent);
        let tt = dyad(second, second);
        let parallel_modes = [
            core::array::from_fn(|i| ss[i] + tt[i]),
            core::array::from_fn(|i| ss[i] - tt[i]),
            symmetric_dyad(tangent, second),
        ];
        let e = elastic.youngs;
        let nu = elastic.poisson;
        let mu = positive((0.5*e)/(1.0+nu))?;
        let normal_modulus = positive((e/(1.0+nu))*((1.0-nu)/(1.0-2.0*nu)))?;
        let r = nu/(1.0-nu);
        let modes = [
            core::array::from_fn(|i| nn[i] + r*parallel_modes[0][i]),
            symmetric_dyad(normal, tangent),
            symmetric_dyad(normal, second),
        ];
        // Positive rank-one plane-stress decomposition avoids subtracting two
        // large bulk tensors, especially near incompressibility or auxeticity.
        let parallel_moduli = [positive((0.5*e)/(1.0-nu))?, mu, mu];
        let mut parallel_tangent = [[0.0; 6]; 6];
        for k in 0..3 {
            add_outer(&mut parallel_tangent, &parallel_modes[k], &parallel_modes[k], parallel_moduli[k]);
        }
        finite_tangent(&parallel_tangent)?;
        let width = settings.band_width_m;
        let normal_stiffness = positive(normal_modulus/width)?;
        let shear_stiffness = positive(mu/width)?;
        // The exact operator norm bounds every admitted total strain. The
        // native envelope additionally requires room for its pure-mode failure
        // distances; the ACTUAL total strain is still checked on every trial.
        let mode_norm = det::sqrt(1.0+2.0*r*r).max(core::f64::consts::SQRT_2);
        let maximum = positive(width*(mode_norm*settings.maximum_strain))?
            .max(positive(2.0*(parameters.mode_i_fracture_energy_j_m2/parameters.tensile_strength_pa))?)
            .max(positive(2.0*(parameters.mode_ii_fracture_energy_j_m2/parameters.shear_strength_pa))?);
        let maximum = positive(maximum*(1.0+64.0*f64::EPSILON))?;
        let cohesive = BilinearCohesive::new_mixed_mode(
            MixedModeParameters {
                normal_stiffness_pa_m: normal_stiffness,
                tangential_stiffness_pa_m: shear_stiffness,
                normal_strength_pa: parameters.tensile_strength_pa,
                tangential_strength_pa: parameters.shear_strength_pa,
                mode_i_fracture_energy_j_m2: parameters.mode_i_fracture_energy_j_m2,
                mode_ii_fracture_energy_j_m2: parameters.mode_ii_fracture_energy_j_m2,
                bk_exponent: parameters.bk_exponent,
            },
            normal_stiffness,
            maximum,
        )?;
        positive(parameters.mode_i_fracture_energy_j_m2/width)?;
        positive(parameters.mode_ii_fracture_energy_j_m2/width)?;
        let area_energy_tolerance = positive(settings.energy_tolerance_j_m3*width)?;
        Ok(Self {
            elastic, parameters, settings: BrittleCrackBandSettings { normal, ..settings },
            axes, modes, parallel_modes, parallel_moduli, parallel_tangent,
            cohesive, area_energy_tolerance,
        })
    }

    #[must_use]
    pub const fn elastic(&self) -> &IsotropicElastic { &self.elastic }
    #[must_use]
    pub const fn parameters(&self) -> BrittleCrackBandParameters { self.parameters }
    /// Settings with the admitted unit normal; no coefficients change in history.
    #[must_use]
    pub const fn settings(&self) -> BrittleCrackBandSettings { self.settings }
    #[must_use]
    pub const fn band_width_m(&self) -> f64 { self.settings.band_width_m }
    #[must_use]
    pub const fn normal(&self) -> [f64; 3] { self.settings.normal }
    /// Rows are the fixed normal and two in-plane tangents.
    #[must_use]
    pub const fn local_axes(&self) -> [[f64; 3]; 3] { self.axes }
    #[must_use]
    pub const fn cohesive_law(&self) -> &BilinearCohesive { &self.cohesive }
    #[must_use]
    pub fn small_strain_law(&self) -> BrittleSmallStrain<'_> {
        BrittleSmallStrain { owner: self }
    }
    #[must_use]
    pub fn initial_state(&self) -> BrittleState<'_> {
        BrittleState { owner: self, strain: [0.0; 6], point: self.cohesive.bonded_point() }
    }

    /// Evolve along a straight TOTAL strain path from the accepted state.
    /// Histories cannot be silently moved between widths, planes or materials.
    pub fn trial<'a>(
        &'a self,
        strain: &Voigt,
        state: &BrittleState<'a>,
    ) -> Result<BrittleStep<'a>, BrittleError> {
        let jump = self.jump_for_trial(strain, state)?;
        let step = state.point.trial(jump, self.area_energy_tolerance)?;
        self.finish_trial(strain, state, step)
    }

    // Shared admission keeps endpoint and path-mean operations on exactly the
    // same original owner, geometry and Frobenius strain domain.
    fn jump_for_trial(
        &self,
        strain: &Voigt,
        state: &BrittleState<'_>,
    ) -> Result<[f64; 3], BrittleError> {
        if !core::ptr::eq(state.owner, self) || !core::ptr::eq(state.point.law(), &self.cohesive) {
            return Err(BrittleError::StateMismatch);
        }
        self.check_strain(&state.strain)?;
        self.check_strain(strain)?;
        Ok(self.modes.map(|mode| self.settings.band_width_m*contract(&mode, strain)))
    }

    // Project one already accepted native trial. Path operations reuse this
    // endpoint assembly instead of repeating constitutive evolution.
    fn finish_trial<'a>(
        &'a self,
        strain: &Voigt,
        state: &BrittleState<'a>,
        step: CohesiveStep<'a>,
    ) -> Result<BrittleStep<'a>, BrittleError> {
        let width = self.settings.band_width_m;
        let mut stress = self.parallel_stress(strain);
        for (mode, traction) in self.modes.iter().zip(step.traction_pa) {
            for i in 0..6 { stress[i] += mode[i]*traction; }
        }
        finite_voigt(&stress)?;
        let branch_tangent = self.tangent_from_modes(step.branch_tangent_pa_m)?;
        let classical = step.tangent_pa_m.is_some();
        let parallel_work = self.parallel_modes.iter().zip(self.parallel_moduli)
            .map(|(mode, modulus)| {
                let a = contract(mode, &state.strain);
                let b = contract(mode, strain);
                (0.5*modulus*(b-a))*(b+a)
            }).sum::<f64>();
        let work = finite(parallel_work + step.energy.work/width)?;
        let stored_change = finite(parallel_work + step.energy.stored_change/width)?;
        let dissipated = finite(step.energy.dissipated/width)?;
        let defect = finite((work-stored_change)-dissipated)?;
        if defect.abs() > self.settings.energy_tolerance_j_m3 {
            return Err(BrittleError::EnergyDefect(defect));
        }
        let next = BrittleState { owner: self, strain: *strain, point: step.point };
        finite(next.stored_energy_j_m3())?;
        finite(next.dissipated_j_m3())?;
        Ok(BrittleStep {
            stress_pa: stress,
            tangent_pa: classical.then_some(branch_tangent),
            branch_tangent_pa: branch_tangent,
            state: next,
            energy: BrittleEnergy { work, stored_change, dissipated, defect },
        })
    }

    fn check_strain(&self, strain: &Voigt) -> Result<(), BrittleError> {
        finite_voigt(strain)?;
        let scale = strain.iter().fold(0.0_f64, |a, b| a.max(b.abs()));
        if scale == 0.0 { return Ok(()); }
        let scaled = strain.map(|v| v/scale);
        let norm = finite(scale*det::sqrt(contract(&scaled, &scaled)))?;
        if norm > self.settings.maximum_strain {
            return Err(BrittleError::StrainLimit { norm, maximum: self.settings.maximum_strain });
        }
        Ok(())
    }

    fn parallel_stress(&self, strain: &Voigt) -> Voigt {
        let mut stress = [0.0; 6];
        for (mode, modulus) in self.parallel_modes.iter().zip(self.parallel_moduli) {
            let value = modulus*contract(mode, strain);
            for i in 0..6 { stress[i] += mode[i]*value; }
        }
        stress
    }

    fn parallel_energy(&self, strain: &Voigt) -> f64 {
        self.parallel_modes.iter().zip(self.parallel_moduli).map(|(mode, modulus)| {
            let value = contract(mode, strain);
            (0.5*modulus*value)*value
        }).sum()
    }

    fn tangent_from_modes(&self, local: [[f64; 3]; 3]) -> Result<Tangent6, BrittleError> {
        let mut tangent = self.parallel_tangent;
        for (i, row) in local.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                add_outer(&mut tangent, &self.modes[i], &self.modes[j],
                    finite(self.settings.band_width_m*value)?);
            }
        }
        finite_tangent(&tangent)?;
        Ok(tangent)
    }
}

impl BrittleState<'_> {
    /// Infallible trait failures carry a nonfinite strain marker; checked trials
    /// never return such a state. Reusing a marker is always a checked refusal.
    #[must_use]
    pub fn is_valid(&self) -> bool { self.strain.iter().all(|v| v.is_finite()) }
    #[must_use]
    pub const fn strain(&self) -> &Voigt { &self.strain }
    #[must_use]
    pub fn damage(&self) -> f64 {
        if self.is_valid() { self.point.damage() } else { f64::NAN }
    }
    #[must_use]
    pub fn failed(&self) -> bool { self.is_valid() && self.point.failed() }
    #[must_use]
    pub fn dissipated_j_m3(&self) -> f64 {
        if self.is_valid() { self.point.dissipated_j_m2()/self.owner.settings.band_width_m }
        else { f64::NAN }
    }
    #[must_use]
    pub fn stored_energy_j_m3(&self) -> f64 {
        if self.is_valid() {
            self.owner.parallel_energy(&self.strain)
                + self.point.stored_energy_j_m2()/self.owner.settings.band_width_m
        } else { f64::NAN }
    }
    /// Original immutable cohesive observation, including its frozen onset mix.
    #[must_use]
    pub const fn cohesive_point(&self) -> &CohesivePoint<'_> { &self.point }
}

impl<'a> SmallStrainLaw for BrittleSmallStrain<'a> {
    type State = BrittleState<'a>;

    fn initial_state(&self) -> Self::State { self.owner.initial_state() }
    fn stress(&self, strain: &Voigt, state: &Self::State) -> Voigt {
        self.owner.trial(strain, state).map(|step| step.stress_pa).unwrap_or([f64::NAN; 6])
    }
    fn tangent(&self, strain: &Voigt, state: &Self::State) -> Tangent6 {
        // The trait uses the declared semismooth branch at switches. Callers
        // requiring a classical derivative use the checked Option above.
        self.owner.trial(strain, state).map(|step| step.branch_tangent_pa)
            .unwrap_or([[f64::NAN; 6]; 6])
    }
    fn update_state(&self, strain: &Voigt, state: &Self::State) -> Self::State {
        self.owner.trial(strain, state).map(|step| step.state).unwrap_or(BrittleState {
            owner: self.owner, strain: [f64::NAN; 6], point: self.owner.cohesive.bonded_point(),
        })
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            has_stored_energy: true,
            dissipation_nonnegative: true,
            polyconvex: None,
            tangent_symmetric: false,
            failure_envelope: "fixed-plane mixed-mode brittle crack band; no crushing, crack nucleation or fragment topology",
        }
    }
    fn card(&self) -> ModelCard {
        ModelCard::new(
            "material.brittle-fixed-plane-crack-band", "0.1.0", Ambition::Solid,
            vec![
                "small strain/rotation and one fixed material plane with explicit width".into(),
                "isotropic elastic backbone; native irreversible mixed-mode cohesive evolution".into(),
                "quadratic strength initiation; BK mode mix frozen at first onset".into(),
                "constitutive reclosure with undamaged compression; in-plane response stays elastic".into(),
            ],
            ValidityDomain::unconstrained().with("strain-magnitude", 0.0, self.owner.settings.maximum_strain),
            vec![
                "no automatic nucleation, evolving crack plane, fragments or geometric contact".into(),
                "no arbitrary-mesh objectivity, crushing, friction or rate dependence; thermal transfer is an explicit fixed-projection step".into(),
                "fracture-energy matching is not experimental material qualification".into(),
            ],
            // No calibration establishes a finite model-form discrepancy band.
            f64::INFINITY,
        )
    }
}

fn add_outer(tangent: &mut Tangent6, a: &Voigt, b: &Voigt, modulus: f64) {
    for i in 0..6 {
        for j in 0..6 {
            let weight = if j < 3 { 1.0 } else { 2.0 };
            tangent[i][j] += a[i]*modulus*b[j]*weight;
        }
    }
}
fn dyad(a: [f64; 3], b: [f64; 3]) -> Voigt {
    [a[0]*b[0], a[1]*b[1], a[2]*b[2], a[0]*b[1], a[1]*b[2], a[2]*b[0]]
}
fn symmetric_dyad(a: [f64; 3], b: [f64; 3]) -> Voigt {
    let ab = dyad(a, b);
    let ba = dyad(b, a);
    core::array::from_fn(|i| ab[i]+ba[i])
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
}
fn normalize(v: [f64; 3]) -> Result<[f64; 3], BrittleError> {
    for x in v { finite(x)?; }
    let scale = positive(v.iter().fold(0.0_f64, |a, b| a.max(b.abs())))?;
    let scaled = v.map(|x| x/scale);
    let norm = positive(det::sqrt(scaled.iter().map(|x| x*x).sum()))?;
    Ok(scaled.map(|x| x/norm))
}
fn finite(v: f64) -> Result<f64, BrittleError> {
    if v.is_finite() { Ok(v) } else { Err(BrittleError::Invalid("nonfinite input/arithmetic")) }
}
fn positive(v: f64) -> Result<f64, BrittleError> {
    if v.is_finite() && v > 0.0 { Ok(v) }
    else { Err(BrittleError::Invalid("positive finite parameter required")) }
}
fn finite_voigt(v: &Voigt) -> Result<(), BrittleError> {
    for &value in v { finite(value)?; }
    Ok(())
}
fn finite_tangent(v: &Tangent6) -> Result<(), BrittleError> {
    for row in v { finite_voigt(row)?; }
    Ok(())
}

#[cfg(test)]
mod tests;
