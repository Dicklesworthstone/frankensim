//! Progressive ductile degradation on the ORIGINAL energetic native J2 law.
//!
//! After initiation, an explicitly selected profile g(x) multiplies the
//! effective deviatoric and tensile volumetric response, where
//! `x = h/Gp * integral (sigma_y + H_iso*alpha) d alpha`. The compressive bulk
//! response remains intact. Existing constructors retain `g=exp(-x)`; an
//! authored finite-support option `g=max(1-x/3,0)^2` reaches exact zero at x=3.
//! Both integrate to unit degraded equivalent work in x. Thus Gp remains
//! POST-INITIATION EQUIVALENT PLASTIC-WORK CAPACITY per area, not total
//! thermodynamic fracture dissipation. Its primitive equals damage only for
//! the exponential profile. Elastic and energetic hardening release remain
//! separate; no quantity silently becomes heat, crack area or fragment energy.
//!
//! The default exponential energy form and nondeleted-3D compression split are in
//! SIMULIA, *Damage Evolution and Element Removal for Ductile Metals*:
//! <https://docs.software.vt.edu/abaqusv2025/English/SIMACAEMATRefMap/simamat-c-damageevolductile.htm>.
//! The native J2 return map, Prager history and algorithmic tangent are reused.
//! Damage supplies their analytic chain rule; no return map is duplicated.
//!
//! Initiation retains the parent's endpoint quadrature on the undegraded trial
//! stress, which is the physical stress before damage. A first crossing spends
//! only the estimated post-onset part of that plastic increment. Its moving
//! onset is differentiated, including triaxiality, Lode and rate dependence. Once
//! reached, the initiation observation freezes; ongoing plastic strain belongs
//! to `native_state()`. Subsequent evolution does not reapply an initiation
//! rate/triaxiality envelope to an already initiated point. Original temperature
//! support remains mandatory, even during elastic unloading.
//!
//! Endpoint work closes into storage, native plastic dissipation, damage
//! release, and backward-Euler numerical loss. This is a discrete work identity,
//! not exact integration of an unresolved strain path. An explicit damage-step
//! allowance refuses large trials rather than clipping damage or its tangent.
//! Native continuation must resolve onset, reversals and softening with physical
//! step refinement. The tangent can be nonsymmetric and indefinite.
//!
//! Characteristic-length scaling alone does not establish mesh-objective crack
//! paths. No element deletion, residual-stiffness floor, finite-strain law,
//! crack-face contact, heat partition or measured validation
//! is supplied. Exponential integrity stays positive; an unrepresentable value
//! refuses instead of creating a hidden erosion threshold. Finite-support
//! failure retains the original compressive response and all native histories.
//! A reduced plane-stress owner must refuse a singular transverse constraint;
//! the scalar profile does not select a terminal thickness or crack opening.

use core::fmt::Debug;
use fs_evidence::{Ambition, ModelCard, ValidityDomain};
use fs_math::det;

use super::{
    DuctileError, DuctileInitiationLaw, DuctileState, DuctileStressState, JohnsonCookInitiation,
    TabulatedDuctileLocus,
};
use crate::plastic::energy::{EnergyError, EnergyTolerances};
use crate::plastic::{J2Plasticity, J2State};
use crate::tensor::{contract, trace};
use crate::{MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};

/// Explicit profiles in the normalized original equivalent-work coordinate.
mod work;
pub use work::DuctileWorkEvolution;

/// Original source-bound post-initiation work, not a toughness alias.
pub mod source;

/// Native plane-stress section work with separate plastic and damage losses.
pub mod section;

/// Endpoint-locus derivatives needed only while a trial crosses initiation.
/// Derivatives use tensile-positive triaxiality, the normalized third
/// deviatoric invariant xi and physical inverse seconds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuctileLocusGradient {
    pub triaxiality: f64,
    pub lode_parameter: f64,
    pub strain_rate_s_inv: f64,
}

/// A differentiable initiation owner for the exact constitutive chain rule.
/// At a tabulation knot, the right branch is selected (left at the last knot).
/// That selection is a semismooth branch, not a classical derivative there.
pub trait DifferentiableDuctileInitiationLaw: DuctileInitiationLaw {
    fn fracture_strain_gradient(
        &self,
        eta: f64,
        rate: f64,
        temperature_k: f64,
    ) -> Result<DuctileLocusGradient, DuctileError>;

    /// Full-invariant derivative on the same branch as `fracture_strain_at`.
    /// Lode-dependent laws override this entry point and refuse the scalar
    /// derivative when xi is absent. Reduced laws retain a zero xi derivative.
    fn fracture_strain_gradient_at(
        &self,
        stress: DuctileStressState,
        rate: f64,
        temperature_k: f64,
    ) -> Result<DuctileLocusGradient, DuctileError> {
        stress.validate()?;
        self.fracture_strain_gradient(stress.triaxiality, rate, temperature_k)
    }
}

impl<L: DifferentiableDuctileInitiationLaw + ?Sized> DifferentiableDuctileInitiationLaw for &L {
    fn fracture_strain_gradient(
        &self,
        eta: f64,
        rate: f64,
        temperature: f64,
    ) -> Result<DuctileLocusGradient, DuctileError> {
        (**self).fracture_strain_gradient(eta, rate, temperature)
    }

    fn fracture_strain_gradient_at(
        &self,
        stress: DuctileStressState,
        rate: f64,
        temperature: f64,
    ) -> Result<DuctileLocusGradient, DuctileError> {
        (**self).fracture_strain_gradient_at(stress, rate, temperature)
    }
}

impl DifferentiableDuctileInitiationLaw for TabulatedDuctileLocus {
    fn fracture_strain_gradient(
        &self,
        eta: f64,
        rate: f64,
        temperature: f64,
    ) -> Result<DuctileLocusGradient, DuctileError> {
        self.fracture_strain(eta, rate, temperature)?;
        let knots = self.knots();
        let upper = knots
            .partition_point(|knot| knot.0 <= eta)
            .min(knots.len() - 1)
            .max(1);
        let (a, fa) = knots[upper - 1];
        let (b, fb) = knots[upper];
        locus_gradient((fb - fa) / (b - a), 0.0)
    }
}

impl DifferentiableDuctileInitiationLaw for JohnsonCookInitiation {
    fn fracture_strain_gradient(
        &self,
        eta: f64,
        rate: f64,
        temperature: f64,
    ) -> Result<DuctileLocusGradient, DuctileError> {
        let strain = self.fracture_strain(eta, rate, temperature)?;
        let p = self.parameters();
        let d_eta = if p.d2 == 0.0 || p.d3 == 0.0 {
            0.0
        } else {
            let exponential = p.d2 * det::exp(-p.d3 * eta);
            strain * ((-p.d3 * exponential) / (p.d1 + exponential))
        };
        let d_rate = if p.d4 == 0.0 {
            0.0
        } else {
            let ratio = rate / p.reference_strain_rate_s_inv;
            let logarithm = if ratio.is_finite() && ratio > 0.0 {
                det::ln(ratio)
            } else {
                det::ln(rate) - det::ln(p.reference_strain_rate_s_inv)
            };
            (strain / rate) * (p.d4 / (1.0 + p.d4 * logarithm))
        };
        locus_gradient(d_eta, d_rate)
    }
}

fn locus_gradient(
    triaxiality: f64,
    strain_rate_s_inv: f64,
) -> Result<DuctileLocusGradient, DuctileError> {
    if triaxiality.is_finite() && strain_rate_s_inv.is_finite() {
        Ok(DuctileLocusGradient {
            triaxiality,
            lode_parameter: 0.0,
            strain_rate_s_inv,
        })
    } else {
        Err(DuctileError::Numerical(
            "initiation-locus derivative is not representable",
        ))
    }
}

/// Explicit local regularization and numerical gates; no material defaults.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuctileJ2Settings {
    /// Fixed reference band length, m. The geometry owner supplies its meaning.
    pub characteristic_length_m: f64,
    /// Post-initiation equivalent plastic work per area, J/m2. This governs
    /// `integral g*(sigma_y+H_iso*alpha)*h*dalpha`, not total fracture loss.
    pub post_initiation_plastic_work_j_m2: f64,
    /// Maximum accepted increase in `1-g`, in (0,1]. A larger trial refuses;
    /// native continuation must cut back. No damage or tangent is clipped.
    pub maximum_damage_increment: f64,
    pub energy_tolerances: EnergyTolerances,
}

/// Actual interval and fixed endpoint temperature of one native solve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuctileJ2Conditions {
    pub duration_s: f64,
    pub temperature_k: f64,
}

/// Native energetic J2 material, initiation law and immutable receiving length.
#[derive(Clone, Debug, PartialEq)]
pub struct DuctileJ2<L> {
    native: J2Plasticity,
    initiation: L,
    settings: DuctileJ2Settings,
    bulk_modulus: f64,
    work_scale: f64,
    evolution: DuctileWorkEvolution,
}

/// Accepted history bound to this exact owner, not an interchangeable clone.
#[derive(Clone, Debug, PartialEq)]
pub struct DuctileJ2State<'a, L> {
    owner: &'a DuctileJ2<L>,
    strain: Voigt,
    native: J2State,
    initiation: DuctileState,
    // Original normalized work x; an exponent only for the default profile.
    exponent: f64,
    stored: f64,
    plastic_dissipated: f64,
    damage_dissipated: f64,
    numerical_dissipated: f64,
}

/// Incremental work per reference volume, J/m3, except the named area accounts.
/// The equivalent-work area accounts are diagnostics of the attenuation law,
/// not extra terms to add to the mechanical work identity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuctileJ2StepEnergy {
    pub initial_stored_j_m3: f64,
    pub final_stored_j_m3: f64,
    pub mechanical_work_j_m3: f64,
    pub plastic_dissipation_j_m3: f64,
    pub damage_release_j_m3: f64,
    pub numerical_dissipation_j_m3: f64,
    pub residual_j_m3: f64,
    pub equivalent_plastic_work_increment_j_m2: f64,
    pub degraded_equivalent_plastic_work_increment_j_m2: f64,
}

/// Disposable constitutive proposal; the old state is unchanged on every path.
#[derive(Clone, Debug, PartialEq)]
pub struct DuctileJ2Step<'a, L> {
    pub stress_pa: Voigt,
    /// Exact algorithmic derivative on smooth branches; semismooth selection
    /// at yield, onset, table knots and zero volumetric strain. May be indefinite
    /// and nonsymmetric, so an SPD-only equilibrium solver is unsuitable.
    pub branch_tangent_pa: Tangent6,
    pub state: DuctileJ2State<'a, L>,
    pub energy: DuctileJ2StepEnergy,
}

/// Borrowed adapter for the existing tet, hex and reduced small-strain owners.
/// The duration and temperature belong to this solve, not to material history.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuctileJ2Law<'a, L> {
    owner: &'a DuctileJ2<L>,
    conditions: DuctileJ2Conditions,
}

#[derive(Clone, Debug, PartialEq)]
pub enum DuctileEvolutionError {
    Invalid(&'static str),
    StateMismatch,
    Initiation(DuctileError),
    Native(EnergyError),
    DamageIncrement { proposed: f64, maximum: f64 },
    EnergyDefect(f64),
}
impl core::fmt::Display for DuctileEvolutionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "ductile J2 evolution refused: {self:?}")
    }
}
impl std::error::Error for DuctileEvolutionError {}
impl From<DuctileError> for DuctileEvolutionError {
    fn from(value: DuctileError) -> Self {
        Self::Initiation(value)
    }
}
impl From<EnergyError> for DuctileEvolutionError {
    fn from(value: EnergyError) -> Self {
        Self::Native(value)
    }
}

impl<L: DifferentiableDuctileInitiationLaw> DuctileJ2<L> {
    pub fn new(
        native: J2Plasticity,
        initiation: L,
        settings: DuctileJ2Settings,
    ) -> Result<Self, DuctileEvolutionError> {
        native.validate_energy_law()?;
        initiation.domain().validate()?;
        positive(settings.characteristic_length_m)?;
        positive(settings.post_initiation_plastic_work_j_m2)?;
        positive(settings.maximum_damage_increment)?;
        if settings.maximum_damage_increment > 1.0 {
            return Err(DuctileEvolutionError::Invalid(
                "maximum damage increment must be in (0,1]",
            ));
        }
        positive(settings.energy_tolerances.energy_j_m3)?;
        positive(settings.energy_tolerances.yield_pa)?;
        let (lambda, mu) = native.elastic.lame();
        let bulk_modulus = positive(lambda + 2.0 * mu / 3.0)?;
        let work_scale = positive(
            settings.characteristic_length_m / settings.post_initiation_plastic_work_j_m2,
        )?;
        Ok(Self {
            native,
            initiation,
            settings,
            bulk_modulus,
            work_scale,
            evolution: DuctileWorkEvolution::Exponential,
        })
    }

    /// Explicitly author the post-initiation shape on a new immutable owner.
    /// The original Gp remains degraded equivalent-work capacity. Selecting a
    /// shape neither changes its units nor qualifies that shape experimentally.
    pub fn with_evolution(mut self, evolution: DuctileWorkEvolution) -> Self {
        self.evolution = evolution;
        self
    }

    pub const fn evolution(&self) -> DuctileWorkEvolution {
        self.evolution
    }

    pub const fn native(&self) -> &J2Plasticity {
        &self.native
    }
    pub const fn initiation_law(&self) -> &L {
        &self.initiation
    }
    pub const fn settings(&self) -> DuctileJ2Settings {
        self.settings
    }
    /// Cheap ownership admission with no constitutive trial or plastic flow.
    /// Trait-failure markers still require the separate `state.is_valid()` gate.
    pub fn owns_state(&self, state: &DuctileJ2State<'_, L>) -> bool {
        core::ptr::eq(state.owner, self)
    }

    pub fn initial_state(&self) -> DuctileJ2State<'_, L> {
        DuctileJ2State {
            owner: self,
            strain: [0.0; 6],
            native: self.native.initial_state(),
            initiation: DuctileState::virgin(),
            exponent: 0.0,
            stored: 0.0,
            plastic_dissipated: 0.0,
            damage_dissipated: 0.0,
            numerical_dissipated: 0.0,
        }
    }

    pub fn small_strain_law(
        &self,
        conditions: DuctileJ2Conditions,
    ) -> Result<DuctileJ2Law<'_, L>, DuctileEvolutionError> {
        self.check_conditions(conditions)?;
        Ok(DuctileJ2Law {
            owner: self,
            conditions,
        })
    }

    fn check_conditions(
        &self,
        conditions: DuctileJ2Conditions,
    ) -> Result<(), DuctileEvolutionError> {
        positive(conditions.duration_s)?;
        positive(conditions.temperature_k)?;
        self.initiation
            .domain()
            .check_temperature(conditions.temperature_k)?;
        Ok(())
    }

    /// Reuse one native update and its energy audit. This fixed-size operation
    /// allocates no memory. A batch owner supplies cancellation and commit gates.
    pub fn trial<'a>(
        &'a self,
        strain: &Voigt,
        old: &DuctileJ2State<'a, L>,
        conditions: DuctileJ2Conditions,
    ) -> Result<DuctileJ2Step<'a, L>, DuctileEvolutionError> {
        if !core::ptr::eq(old.owner, self) {
            return Err(DuctileEvolutionError::StateMismatch);
        }
        self.check_conditions(conditions)?;
        if !old.is_valid() {
            return Err(DuctileEvolutionError::Invalid("invalid accepted history"));
        }
        // Native equilibrium evaluates the exact accepted displacement before
        // imposing its next load. A second radial return there can invent a
        // positive plastic rate from yield roundoff. This rate-independent law
        // must replay an unchanged strain exactly. No tolerance is used: every
        // nonzero strain change still takes the complete checked update below.
        if strain == &old.strain {
            return self.replay_accepted(old);
        }
        let native = self.native.advance_with_energy(
            &old.strain,
            &old.native,
            strain,
            self.settings.energy_tolerances,
        )?;
        let native_tangent = self.native.tangent(strain, &old.native);
        let alpha = native.state.alpha;
        let delta_alpha = nonnegative(alpha - old.native.alpha)?;
        let mut initiation = old.initiation;
        let mut onset_alpha = alpha;
        let mut onset_gradient = [0.0; 6];
        let alpha_gradient = self
            .native
            .accumulated_plastic_strain_gradient(strain, &old.native);
        if old.initiation.initiated() {
            onset_alpha = old.native.alpha;
        } else {
            let step = old.initiation.trial(
                &self.initiation,
                alpha,
                conditions.duration_s,
                conditions.temperature_k,
                &native.stress_pa,
            )?;
            initiation = step.state;
            if initiation.utilization() > 1.0 {
                let fracture_strain =
                    step.fracture_strain.ok_or(DuctileEvolutionError::Invalid(
                        "a new initiation crossing needs positive native plastic flow",
                    ))?;
                let remaining = 1.0 - old.initiation.utilization();
                onset_alpha = finite(old.native.alpha + remaining * fracture_strain)?;
                // Differentiate the same endpoint-held onset estimate, rather
                // than silently freezing it during the first Newton increment.
                let stress_state = DuctileStressState {
                    triaxiality: step
                        .triaxiality
                        .ok_or(DuctileEvolutionError::Invalid("missing onset triaxiality"))?,
                    lode_parameter: step.lode_parameter.ok_or(DuctileEvolutionError::Invalid(
                        "missing onset Lode parameter",
                    ))?,
                };
                let gradient = self.initiation.fracture_strain_gradient_at(
                    stress_state,
                    step.strain_rate_s_inv,
                    conditions.temperature_k,
                )?;
                let stress_gradient =
                    if gradient.triaxiality == 0.0 && gradient.lode_parameter == 0.0 {
                        [0.0; 6]
                    } else {
                        locus_stress_gradient(
                            &native.stress_pa,
                            &native_tangent,
                            stress_state,
                            gradient,
                        )?
                    };
                for j in 0..6 {
                    onset_gradient[j] = finite(
                        remaining
                            * (stress_gradient[j]
                                + gradient.strain_rate_s_inv
                                    * (alpha_gradient[j] / conditions.duration_s)),
                    )?;
                }
            }
        }
        let post_alpha = if onset_alpha < alpha {
            alpha - onset_alpha
        } else {
            0.0
        };
        if post_alpha > delta_alpha * (1.0 + 64.0 * f64::EPSILON) {
            return Err(DuctileEvolutionError::Invalid(
                "onset estimate exceeds the native plastic increment",
            ));
        }
        let resistance = |a: f64| self.native.yield_stress + self.native.hardening * a;
        let equivalent_work = nonnegative(
            post_alpha
                * (self.native.yield_stress
                    + self.native.hardening * (0.5 * onset_alpha + 0.5 * alpha)),
        )?;
        let exponent_increment = nonnegative(self.work_scale * equivalent_work)?;
        let exponent = nonnegative(old.exponent + exponent_increment)?;
        finite(self.settings.post_initiation_plastic_work_j_m2 * exponent)?;
        let integrity = self.evolution.integrity(exponent);
        let integrity = match self.evolution {
            DuctileWorkEvolution::Exponential => positive(integrity)?,
            DuctileWorkEvolution::FiniteSupportQuadratic => nonnegative(integrity)?,
        };
        let old_integrity = old.integrity();
        let profile_increment = self.evolution.increment(old.exponent, exponent_increment);
        let damage_increment = nonnegative(profile_increment.damage)?;
        let degradation_slope = self.evolution.degradation_slope(exponent);
        if damage_increment > self.settings.maximum_damage_increment {
            return Err(DuctileEvolutionError::DamageIncrement {
                proposed: damage_increment,
                maximum: self.settings.maximum_damage_increment,
            });
        }
        let old_compression = self.compression(&old.strain)?;
        let compression = self.compression(strain)?;
        let old_native_storage = native.energy.initial.total_j_m3();
        let native_storage = native.energy.final_state.total_j_m3();
        let old_storage =
            finite(old_integrity * old_native_storage + (1.0 - old_integrity) * old_compression.1)?;
        let stored = finite(integrity * native_storage + (1.0 - integrity) * compression.1)?;
        let stress_pa = core::array::from_fn(|i| {
            integrity * native.stress_pa[i]
                + if i < 3 {
                    (1.0 - integrity) * compression.0
                } else {
                    0.0
                }
        });
        finite_voigt(&stress_pa)?;
        let mut branch_tangent_pa = native_tangent;
        for i in 0..6 {
            let effective_active_stress =
                native.stress_pa[i] - if i < 3 { compression.0 } else { 0.0 };
            for j in 0..6 {
                let compression_tangent = if i < 3 && j < 3 && trace(strain) < 0.0 {
                    self.bulk_modulus
                } else {
                    0.0
                };
                let exponent_gradient = if post_alpha > 0.0 {
                    self.work_scale
                        * (resistance(alpha) * alpha_gradient[j]
                            - resistance(onset_alpha) * onset_gradient[j])
                } else {
                    0.0
                };
                branch_tangent_pa[i][j] = finite(
                    integrity * native_tangent[i][j] + (1.0 - integrity) * compression_tangent
                        - effective_active_stress * (degradation_slope * exponent_gradient),
                )?;
            }
        }
        let delta_strain = core::array::from_fn(|i| strain[i] - old.strain[i]);
        let work = finite(contract(&stress_pa, &delta_strain))?;
        let plastic_dissipation = nonnegative(integrity * native.energy.physical_dissipation_j_m3)?;
        let available_release = old_native_storage - old_compression.1;
        // Both are independently evaluated convex energies. Their difference
        // is nonnegative mathematically; permit only the declared floating-point
        // work allowance when hydrostatic cancellation rounds it below zero.
        if available_release < -self.settings.energy_tolerances.energy_j_m3 {
            return Err(DuctileEvolutionError::EnergyDefect(available_release));
        }
        let damage_release = nonnegative(damage_increment * available_release.max(0.0))?;
        let compression_loss = self.compression_loss(trace(&old.strain), trace(strain))?;
        let numerical = nonnegative(
            integrity * native.energy.numerical_dissipation_j_m3
                + (1.0 - integrity) * compression_loss,
        )?;
        let residual = finite(
            work - (stored - old_storage) - plastic_dissipation - damage_release - numerical,
        )?;
        if residual.abs() > self.settings.energy_tolerances.energy_j_m3 {
            return Err(DuctileEvolutionError::EnergyDefect(residual));
        }
        let state = DuctileJ2State {
            owner: self,
            strain: *strain,
            native: native.state,
            initiation,
            exponent,
            stored,
            plastic_dissipated: nonnegative(old.plastic_dissipated + plastic_dissipation)?,
            damage_dissipated: nonnegative(old.damage_dissipated + damage_release)?,
            numerical_dissipated: nonnegative(old.numerical_dissipated + numerical)?,
        };
        Ok(DuctileJ2Step {
            stress_pa,
            branch_tangent_pa,
            state,
            energy: DuctileJ2StepEnergy {
                initial_stored_j_m3: old_storage,
                final_stored_j_m3: stored,
                mechanical_work_j_m3: work,
                plastic_dissipation_j_m3: plastic_dissipation,
                damage_release_j_m3: damage_release,
                numerical_dissipation_j_m3: numerical,
                residual_j_m3: residual,
                equivalent_plastic_work_increment_j_m2: nonnegative(
                    self.settings.characteristic_length_m * equivalent_work,
                )?,
                degraded_equivalent_plastic_work_increment_j_m2: nonnegative(
                    self.settings.post_initiation_plastic_work_j_m2
                        * profile_increment.degraded_work,
                )?,
            },
        })
    }

    // The accepted point lies on the loading/unloading switch when plastic.
    // Select its retained-history elastic branch, as at a cohesive unloading
    // switch, rather than differentiating fictitious roundoff plastic flow.
    fn replay_accepted<'a>(
        &'a self,
        old: &DuctileJ2State<'a, L>,
    ) -> Result<DuctileJ2Step<'a, L>, DuctileEvolutionError> {
        let elastic_strain = core::array::from_fn(|i| old.strain[i] - old.native.plastic_strain[i]);
        let effective_stress = self.native.elastic.stress(&elastic_strain, &());
        let elastic_tangent = self.native.elastic.stiffness();
        let compression = self.compression(&old.strain)?;
        let compressing = trace(&old.strain) < 0.0;
        let integrity = old.integrity();
        let stress_pa = core::array::from_fn(|i| {
            integrity * effective_stress[i]
                + if i < 3 {
                    (1.0 - integrity) * compression.0
                } else {
                    0.0
                }
        });
        finite_voigt(&stress_pa)?;
        let mut branch_tangent_pa = elastic_tangent;
        for i in 0..6 {
            for j in 0..6 {
                branch_tangent_pa[i][j] = finite(
                    integrity * elastic_tangent[i][j]
                        + if i < 3 && j < 3 && compressing {
                            (1.0 - integrity) * self.bulk_modulus
                        } else {
                            0.0
                        },
                )?;
            }
        }
        Ok(DuctileJ2Step {
            stress_pa,
            branch_tangent_pa,
            state: DuctileJ2State {
                owner: self,
                strain: old.strain,
                native: old.native.clone(),
                initiation: old.initiation,
                exponent: old.exponent,
                stored: old.stored,
                plastic_dissipated: old.plastic_dissipated,
                damage_dissipated: old.damage_dissipated,
                numerical_dissipated: old.numerical_dissipated,
            },
            energy: DuctileJ2StepEnergy {
                initial_stored_j_m3: old.stored,
                final_stored_j_m3: old.stored,
                mechanical_work_j_m3: 0.0,
                plastic_dissipation_j_m3: 0.0,
                damage_release_j_m3: 0.0,
                numerical_dissipation_j_m3: 0.0,
                residual_j_m3: 0.0,
                equivalent_plastic_work_increment_j_m2: 0.0,
                degraded_equivalent_plastic_work_increment_j_m2: 0.0,
            },
        })
    }

    // Pressure and stored energy of the undegraded compressive bulk potential.
    fn compression(&self, strain: &Voigt) -> Result<(f64, f64), DuctileEvolutionError> {
        let volume = finite(trace(strain))?.min(0.0);
        let pressure = finite(self.bulk_modulus * volume)?;
        Ok((pressure, nonnegative((0.5 * pressure) * volume)?))
    }

    // Convex backward-Euler loss of K/2*min(tr(epsilon),0)^2. Evaluating each
    // sign branch avoids subtracting almost equal large compression energies.
    fn compression_loss(&self, old: f64, new: f64) -> Result<f64, DuctileEvolutionError> {
        let loss = if old < 0.0 && new < 0.0 {
            (0.5 * self.bulk_modulus * (new - old)) * (new - old)
        } else if old < 0.0 {
            (0.5 * self.bulk_modulus * old) * old
        } else if new < 0.0 {
            (self.bulk_modulus * new) * (0.5 * new - old)
        } else {
            0.0
        };
        nonnegative(loss)
    }
}

impl<L> DuctileJ2State<'_, L> {
    pub fn is_valid(&self) -> bool {
        self.strain.iter().all(|x| x.is_finite())
            && self.stored.is_finite()
            && self.exponent.is_finite()
            && self.exponent >= 0.0
    }
    pub const fn strain(&self) -> &Voigt {
        &self.strain
    }
    pub const fn native_state(&self) -> &J2State {
        &self.native
    }
    /// Initiation observation, frozen after the first accepted omega>=1. Its
    /// alpha is the onset-step endpoint, not the later native plastic clock.
    pub const fn initiation(&self) -> &DuctileState {
        &self.initiation
    }
    pub fn integrity(&self) -> f64 {
        if self.is_valid() {
            self.owner.evolution.integrity(self.exponent)
        } else {
            f64::NAN
        }
    }
    /// Stable `1-g`; rounding to one does not mean an element was deleted.
    pub fn damage(&self) -> f64 {
        if self.is_valid() {
            self.owner.evolution.damage(self.exponent)
        } else {
            f64::NAN
        }
    }
    /// Exact zero of the explicitly selected finite-support profile. Rounded
    /// damage=1 in the positive exponential tail is never complete failure.
    /// This predicate does not erase history or authorize material removal.
    pub fn is_fully_failed(&self) -> bool {
        self.is_valid() && self.owner.evolution.is_fully_failed(self.exponent)
    }
    pub const fn stored_energy_j_m3(&self) -> f64 {
        self.stored
    }
    pub const fn dissipated_plastic_j_m3(&self) -> f64 {
        self.plastic_dissipated
    }
    pub const fn dissipated_damage_j_m3(&self) -> f64 {
        self.damage_dissipated
    }
    pub const fn numerical_dissipation_j_m3(&self) -> f64 {
        self.numerical_dissipated
    }
    /// Undegraded resistance integral times the original characteristic length.
    pub fn equivalent_plastic_work_j_m2(&self) -> f64 {
        self.owner.settings.post_initiation_plastic_work_j_m2 * self.exponent
    }
    /// Exact scalar attenuation-coordinate work, not additional physical loss.
    pub fn degraded_equivalent_plastic_work_j_m2(&self) -> f64 {
        if self.is_valid() {
            self.owner.settings.post_initiation_plastic_work_j_m2
                * self.owner.evolution.degraded_work(self.exponent)
        } else {
            f64::NAN
        }
    }
}

impl<'a, L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq> SmallStrainLaw
    for DuctileJ2Law<'a, L>
{
    type State = DuctileJ2State<'a, L>;
    fn initial_state(&self) -> Self::State {
        self.owner.initial_state()
    }
    fn stress(&self, strain: &Voigt, state: &Self::State) -> Voigt {
        self.owner
            .trial(strain, state, self.conditions)
            .map(|step| step.stress_pa)
            .unwrap_or([f64::NAN; 6])
    }
    fn tangent(&self, strain: &Voigt, state: &Self::State) -> Tangent6 {
        self.owner
            .trial(strain, state, self.conditions)
            .map(|step| step.branch_tangent_pa)
            .unwrap_or([[f64::NAN; 6]; 6])
    }
    fn update_state(&self, strain: &Voigt, state: &Self::State) -> Self::State {
        self.owner
            .trial(strain, state, self.conditions)
            .map(|step| step.state)
            .unwrap_or_else(|_| {
                let mut invalid = self.owner.initial_state();
                invalid.strain = [f64::NAN; 6];
                invalid.stored = f64::NAN;
                invalid
            })
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            has_stored_energy: true,
            dissipation_nonnegative: true,
            polyconvex: None,
            tangent_symmetric: false,
            failure_envelope: "progressive length-scaled ductile degradation; compressive bulk retained; no erosion or crack topology",
        }
    }
    fn card(&self) -> ModelCard {
        let (id, attenuation) = match self.owner.evolution {
            DuctileWorkEvolution::Exponential => (
                "material.ductile-j2-exponential-work",
                "exponential equivalent-plastic-work attenuation and intact compressive bulk",
            ),
            DuctileWorkEvolution::FiniteSupportQuadratic => (
                "material.ductile-j2-finite-support-work",
                "authored quadratic finite-support equivalent-work attenuation and intact compressive bulk",
            ),
        };
        ModelCard::new(
            id,
            "0.1.0",
            Ambition::Solid,
            vec![
                "native small-strain J2 with energetic isotropic/Prager hardening".into(),
                "endpoint initiation quadrature and fixed characteristic length".into(),
                attenuation.into(),
            ],
            ValidityDomain::unconstrained().with(
                "strain-magnitude",
                0.0,
                self.owner.native.elastic.strain_limit.min(0.05),
            ),
            vec![
                "no mesh-objective crack path, geometric contact, erosion or fragments".into(),
                "post-initiation equivalent plastic work is not total fracture dissipation".into(),
                "no finite deformation, heat partition or experimental qualification".into(),
            ],
            f64::INFINITY,
        )
    }
}

// Chain rule of the initiation locus through the native endpoint stress. With
// U = s/q, d(xi) = [27/2 dev(U^2) - 9/2 xi U] : d(sigma) / q. This polynomial
// derivative remains regular at xi=+/-1; no acos/sin Lode-angle singularity is
// introduced. Voigt shears receive their contraction weight exactly once.
pub(crate) fn locus_stress_gradient<const N: usize>(
    stress: &Voigt,
    tangent: &[[f64; N]; 6],
    state: DuctileStressState,
    locus: DuctileLocusGradient,
) -> Result<[f64; N], DuctileEvolutionError> {
    let scale = positive(stress.iter().fold(0.0_f64, |a, b| a.max(b.abs())))?;
    let s = stress.map(|x| x / scale);
    let differences = [s[0] - s[1], s[1] - s[2], s[2] - s[0]];
    let dev_scale = positive(
        differences
            .iter()
            .chain(&s[3..])
            .fold(0.0_f64, |a, b| a.max(b.abs())),
    )?;
    let diagonal = differences.map(|x| x / dev_scale);
    let dev = [
        (diagonal[0] - diagonal[2]) / 3.0,
        (diagonal[1] - diagonal[0]) / 3.0,
        (diagonal[2] - diagonal[1]) / 3.0,
        s[3] / dev_scale,
        s[4] / dev_scale,
        s[5] / dev_scale,
    ];
    let q = positive(det::sqrt(1.5 * contract(&dev, &dev)))?;
    let unit = dev.map(|value| value / q);
    let xi_derivative = if locus.lode_parameter == 0.0 {
        [0.0; 6]
    } else {
        let square = [
            unit[0] * unit[0] + unit[3] * unit[3] + unit[5] * unit[5],
            unit[1] * unit[1] + unit[3] * unit[3] + unit[4] * unit[4],
            unit[2] * unit[2] + unit[4] * unit[4] + unit[5] * unit[5],
            unit[3] * (unit[0] + unit[1]) + unit[5] * unit[4],
            unit[4] * (unit[1] + unit[2]) + unit[3] * unit[5],
            unit[5] * (unit[0] + unit[2]) + unit[3] * unit[4],
        ];
        let mean_square = trace(&square) / 3.0;
        core::array::from_fn(|i| {
            13.5 * (square[i] - if i < 3 { mean_square } else { 0.0 })
                - 4.5 * state.lode_parameter * unit[i]
        })
    };
    let mut gradient = [0.0; N];
    for j in 0..N {
        let mut derivative = 0.0;
        if locus.triaxiality != 0.0 {
            let mean_gradient = tangent[0][j] / 3.0 + tangent[1][j] / 3.0 + tangent[2][j] / 3.0;
            let mut q_gradient = 0.0;
            for i in 0..6 {
                q_gradient += (1.5 * dev[i] / q) * tangent[i][j] * if i < 3 { 1.0 } else { 2.0 };
            }
            derivative += locus.triaxiality
                * finite(
                    (((mean_gradient - state.triaxiality * q_gradient) / scale) / dev_scale) / q,
                )?;
        }
        if locus.lode_parameter != 0.0 {
            let mut xi_gradient = 0.0;
            for i in 0..6 {
                xi_gradient += xi_derivative[i] * tangent[i][j] * if i < 3 { 1.0 } else { 2.0 };
            }
            derivative += locus.lode_parameter * finite(((xi_gradient / scale) / dev_scale) / q)?;
        }
        gradient[j] = finite(derivative)?;
    }
    Ok(gradient)
}

fn finite(value: f64) -> Result<f64, DuctileEvolutionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DuctileEvolutionError::Invalid(
            "nonfinite evolution input or arithmetic",
        ))
    }
}
fn positive(value: f64) -> Result<f64, DuctileEvolutionError> {
    finite(value)?;
    if value > 0.0 {
        Ok(value)
    } else {
        Err(DuctileEvolutionError::Invalid(
            "positive evolution parameter or integrity required",
        ))
    }
}
fn nonnegative(value: f64) -> Result<f64, DuctileEvolutionError> {
    finite(value)?;
    if value >= 0.0 {
        Ok(value)
    } else {
        Err(DuctileEvolutionError::Invalid(
            "negative evolution energy or increment",
        ))
    }
}
fn finite_voigt(value: &Voigt) -> Result<(), DuctileEvolutionError> {
    for &x in value {
        finite(x)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests;
