//! Progressive ductile attenuation of the original multiplicative finite J2 law.
//!
//! The original plastic metric, Kirchhoff yield radius, spectral return and
//! analytic derivatives remain with [`FiniteJ2`]. After initiation, integrity
//! is `g=exp[-h/Gp integral R(alpha) d alpha]`, with the original linear/Voce/tabulated
//! energetic resistance R. h is a fixed REFERENCE band width and Gp is
//! post-initiation equivalent Kirchhoff plastic work per reference area. Neither
//! is current crack area or total fracture dissipation. The shared initiation
//! owner retains triaxiality, Lode and actual physical-rate dependence.
//!
//! Deviatoric/tensile storage degrades; compressive bulk storage remains intact:
//! `psi=g*psi_native+(1-g)*K*min(log(J),0)^2/2`. Physical stress is its derivative
//! at fixed internal history. The complete algorithmic tangent also includes
//! moving onset and damage. It can be nonsymmetric and is NOT the Hessian of a
//! reduced scalar potential. This owner deliberately implements no
//! `FiniteStrainPotential` interface.
//!
//! All storage and work refer to the original volume. Endpoint `P:delta F`
//! leaves a SIGNED integration defect, including during finite rotation of a
//! prestressed body. Physical plastic loss and released damage storage stay
//! nonnegative and separate. Callers bound/refine the absolute defect; it is
//! never plastic heat. Damage-step limits likewise refuse without clipping.
//! Fixed-temperature Perzyna overstress can be selected explicitly. The same
//! physical interval then drives native flow and damage initiation, including
//! stress relaxation and continued damage at held deformation. The attenuation
//! coordinate retains the EQUILIBRIUM resistance R; viscous overstress contributes
//! to native physical plastic dissipation, not an inferred work calibration.
//! Creep, caloric feedback, crack topology, contact and fragment transfer are
//! not supplied. Temperature dependence belongs only to the initiation locus.

use super::{
    FiniteJ2, FiniteJ2Error, FiniteJ2State, RadialReturn, ReturnMap, ReturnSensitivities, dot, pack,
};
use super::rate::{FiniteOverstressJ2, FiniteOverstressParameters, validate_viscous_committed};
use crate::ductile::evolution::{
    DifferentiableDuctileInitiationLaw, DuctileEvolutionError, DuctileJ2Conditions,
    DuctileJ2Settings, locus_stress_gradient,
};
use crate::ductile::{DuctileError, DuctileState, DuctileStressState};
use fs_math::det;

/// Original finite constitutive card, reference-work calibration and initiation.
pub mod source;

const IDENTITY: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

/// Immutable native finite material, original initiation law and reference width.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteDuctileJ2<L> {
    native: FiniteJ2,
    initiation: L,
    settings: DuctileJ2Settings,
    work_scale: f64,
    overstress: Option<FiniteOverstressParameters>,
}

/// Complete accepted history bound to this exact immutable owner.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteDuctileState<'a, L> {
    owner: &'a FiniteDuctileJ2<L>,
    deformation: [f64; 9],
    native: FiniteJ2State,
    initiation: DuctileState,
    exponent: f64,
    stored: f64,
    native_stored: f64,
    active_stored: f64,
    compression_stored: f64,
    first_piola: [f64; 9],
    kirchhoff: [f64; 9],
    plastic_dissipated: f64,
    damage_dissipated: f64,
    endpoint_defect: f64,
}

/// Per-step reference-volume work. The two area quantities are diagnostics,
/// not additional terms in the mechanical balance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteDuctileEnergy {
    pub initial_stored_j_m3: f64,
    pub final_stored_j_m3: f64,
    pub mechanical_work_j_m3: f64,
    pub plastic_dissipation_j_m3: f64,
    pub damage_release_j_m3: f64,
    /// Signed first-order Piola work error; never physical dissipation or heat.
    pub endpoint_work_defect_j_m3: f64,
    /// Independently assembled split-work identity, not the integration defect.
    pub residual_j_m3: f64,
    pub equivalent_plastic_work_increment_j_m2: f64,
    pub degraded_equivalent_plastic_work_increment_j_m2: f64,
}

/// Disposable full proposal; caller-owned old history never changes.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteDuctileStep<'a, L> {
    pub first_piola_pa: [f64; 9],
    pub kirchhoff_pa: [f64; 9],
    /// Exact branch dP/dF, including the moving initiation event.
    pub tangent_pa: [[f64; 9]; 9],
    pub state: FiniteDuctileState<'a, L>,
    pub energy: FiniteDuctileEnergy,
}

#[derive(Clone, Debug, PartialEq)]
pub enum FiniteDuctileError {
    Invalid(&'static str),
    StateMismatch,
    Native(FiniteJ2Error),
    Initiation(DuctileError),
    DamageIncrement { proposed: f64, maximum: f64 },
    EnergyDefect(f64),
}
impl core::fmt::Display for FiniteDuctileError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite ductile response refused: {self:?}")
    }
}
impl std::error::Error for FiniteDuctileError {}
impl From<FiniteJ2Error> for FiniteDuctileError {
    fn from(error: FiniteJ2Error) -> Self {
        Self::Native(error)
    }
}
impl From<DuctileError> for FiniteDuctileError {
    fn from(error: DuctileError) -> Self {
        Self::Initiation(error)
    }
}

struct Returned {
    native: ReturnMap,
    /// Derivative of the selected inverse-rate law at this returned increment.
    viscous_modulus_pa: f64,
    initiation: DuctileState,
    onset_alpha: f64,
    onset: Option<(DuctileStressState, f64, f64)>,
    post_alpha: f64,
    equivalent_work: f64,
    exponent: f64,
    integrity: f64,
    damage_increment: f64,
    compression_pressure: f64,
    compression_stored: f64,
    active_stored: f64,
    stored: f64,
    first_piola: [f64; 9],
    kirchhoff: [f64; 9],
}

impl<L: DifferentiableDuctileInitiationLaw + Clone> FiniteDuctileJ2<L> {
    pub fn new(
        native: FiniteJ2,
        initiation: L,
        settings: DuctileJ2Settings,
    ) -> Result<Self, FiniteDuctileError> {
        if native.parameters.yield_stress_pa <= 0.0 || native.maximum_creep_alpha.is_some() {
            return Err(FiniteDuctileError::Invalid(
                "finite ductile plasticity requires the original positive-yield law without creep history",
            ));
        }
        initiation.domain().validate()?;
        positive(settings.characteristic_length_m)?;
        positive(settings.post_initiation_plastic_work_j_m2)?;
        positive(settings.maximum_damage_increment)?;
        positive(settings.energy_tolerances.energy_j_m3)?;
        positive(settings.energy_tolerances.yield_pa)?;
        if settings.maximum_damage_increment > 1.0 {
            return Err(FiniteDuctileError::Invalid(
                "maximum damage increment must be in (0,1]",
            ));
        }
        let work_scale = positive(
            settings.characteristic_length_m / settings.post_initiation_plastic_work_j_m2,
        )?;
        Ok(Self {
            native,
            initiation,
            settings,
            work_scale,
            overstress: None,
        })
    }

    /// Select fixed-temperature Perzyna flow on the same native backbone.
    ///
    /// The physical duration is supplied by each trial, never inferred from
    /// solver iterations. This creates a distinct immutable owner: an accepted
    /// history from another owner cannot be relabelled as a viscous specimen.
    /// Viscous work does not change the supplied equilibrium-resistance damage
    /// calibration. No rate coefficients or temperature law are inferred.
    pub fn with_overstress(
        mut self,
        parameters: FiniteOverstressParameters,
    ) -> Result<Self, FiniteDuctileError> {
        parameters.validate()?;
        self.overstress = Some(parameters);
        Ok(self)
    }

    /// Explicit physical kinetics, or the original rate-independent selection.
    pub const fn overstress_parameters(&self) -> Option<&FiniteOverstressParameters> {
        self.overstress.as_ref()
    }

    /// Native storage and equilibrium resistance. If kinetics is selected this
    /// backbone alone does not reproduce the complete constitutive response.
    pub const fn native(&self) -> &FiniteJ2 {
        &self.native
    }
    pub const fn initiation_law(&self) -> &L {
        &self.initiation
    }
    pub const fn settings(&self) -> DuctileJ2Settings {
        self.settings
    }
    pub fn owns_state(&self, state: &FiniteDuctileState<'_, L>) -> bool {
        core::ptr::eq(state.owner, self)
    }

    /// Inspect retained history without advancing the physical clock. Selected
    /// kinetics permits legitimate overstress; only the rate-independent owner
    /// requires the accepted endpoint to lie inside the equilibrium yield surface.
    pub fn validate_committed(
        &self,
        state: &FiniteDuctileState<'_, L>,
        yield_tolerance_pa: f64,
    ) -> Result<(), FiniteDuctileError> {
        if !self.owns_state(state) {
            return Err(FiniteDuctileError::StateMismatch);
        }
        if !state.is_valid() {
            return Err(FiniteDuctileError::Invalid("invalid accepted history"));
        }
        if self.overstress.is_some() {
            validate_viscous_committed(
                &self.native,
                &state.deformation,
                &state.native,
                yield_tolerance_pa,
            )?;
        } else {
            self.native.validate_committed(
                &state.deformation,
                &state.native,
                yield_tolerance_pa,
            )?;
        }
        Ok(())
    }

    /// Admit physical time and temperature without proposing a relaxation step.
    /// Rate/stress-state support is checked on actual plastic flow in `trial`.
    pub fn validate_conditions(
        &self,
        conditions: DuctileJ2Conditions,
    ) -> Result<(), FiniteDuctileError> {
        positive(conditions.duration_s)?;
        positive(conditions.temperature_k)?;
        self.initiation
            .domain()
            .check_temperature(conditions.temperature_k)?;
        Ok(())
    }

    pub fn initial_state(&self) -> FiniteDuctileState<'_, L> {
        FiniteDuctileState {
            owner: self,
            deformation: IDENTITY,
            native: self.native.initial_state(),
            initiation: DuctileState::virgin(),
            exponent: 0.0,
            stored: 0.0,
            native_stored: 0.0,
            active_stored: 0.0,
            compression_stored: 0.0,
            first_piola: [0.0; 9],
            kirchhoff: [0.0; 9],
            plastic_dissipated: 0.0,
            damage_dissipated: 0.0,
            endpoint_defect: 0.0,
        }
    }

    fn admit(
        &self,
        old: &FiniteDuctileState<'_, L>,
        conditions: DuctileJ2Conditions,
    ) -> Result<(), FiniteDuctileError> {
        if !self.owns_state(old) {
            return Err(FiniteDuctileError::StateMismatch);
        }
        if !old.is_valid() {
            return Err(FiniteDuctileError::Invalid("invalid accepted history"));
        }
        self.validate_conditions(conditions)
    }

    /// Residual/backtracking path: same return, source support and damage gates,
    /// without assembling either 9-by-9 derivative or a second spectral solve.
    pub fn stress(
        &self,
        f: &[f64; 9],
        old: &FiniteDuctileState<'_, L>,
        conditions: DuctileJ2Conditions,
    ) -> Result<[f64; 9], FiniteDuctileError> {
        self.admit(old, conditions)?;
        if self.overstress.is_none() && f == &old.deformation {
            return Ok(old.first_piola);
        }
        Ok(self.return_map(f, old, conditions)?.first_piola)
    }

    /// Complete native finite update and analytical damage chain. No accepted
    /// state, calibration or source is changed by success, failure or retry.
    pub fn trial<'a>(
        &'a self,
        f: &[f64; 9],
        old: &FiniteDuctileState<'a, L>,
        conditions: DuctileJ2Conditions,
    ) -> Result<FiniteDuctileStep<'a, L>, FiniteDuctileError> {
        self.admit(old, conditions)?;
        if self.overstress.is_none() && f == &old.deformation {
            return self.replay(old);
        }
        let returned = self.return_map(f, old, conditions)?;
        let mut sensitivities = ReturnSensitivities::default();
        let native_tangent = self.native.return_derivatives(
            f,
            &old.native,
            &returned.native,
            Some(&mut sensitivities),
        )?;
        let alpha = returned.native.state.alpha;
        let mut alpha_gradient = [0.0; 9];
        if returned.native.dgamma > 0.0 {
            let modulus = positive(
                3.0 * self.native.parameters.shear_modulus_pa
                    + self.native.hardening_tangent(alpha)
                    + returned.viscous_modulus_pa,
            )?;
            alpha_gradient = sensitivities.trial_equivalent.map(|dq| dq / modulus);
        }
        let mut onset_gradient = [0.0; 9];
        if let Some((stress_state, rate, remaining)) = returned.onset {
            let gradient = self.initiation.fracture_strain_gradient_at(
                stress_state,
                rate,
                conditions.temperature_k,
            )?;
            let mut spatial = [[0.0; 9]; 6];
            for j in 0..9 {
                let column = pack(&core::array::from_fn(|k| sensitivities.kirchhoff[k][j]));
                for i in 0..6 {
                    spatial[i][j] = column[i];
                }
            }
            let stress_gradient = if gradient.triaxiality == 0.0 && gradient.lode_parameter == 0.0 {
                [0.0; 9]
            } else {
                locus_stress_gradient(
                    &pack(&returned.native.tau),
                    &spatial,
                    stress_state,
                    gradient,
                )
                .map_err(|error| match error {
                    DuctileEvolutionError::Invalid(message) => FiniteDuctileError::Invalid(message),
                    _ => FiniteDuctileError::Invalid(
                        "stress-invariant derivative is not representable",
                    ),
                })?
            };
            for j in 0..9 {
                onset_gradient[j] = finite(
                    remaining
                        * (stress_gradient[j]
                            + gradient.strain_rate_s_inv
                                * (alpha_gradient[j] / conditions.duration_s)),
                )?;
            }
        }
        let mut exponent_gradient = [0.0; 9];
        if returned.post_alpha > 0.0 {
            for j in 0..9 {
                exponent_gradient[j] = finite(
                    self.work_scale
                        * (self.native.yield_radius(alpha) * alpha_gradient[j]
                            - self.native.yield_radius(returned.onset_alpha) * onset_gradient[j]),
                )?;
            }
        }
        let tangent_pa = self.degraded_tangent(
            &returned.native,
            &native_tangent,
            returned.integrity,
            &exponent_gradient,
        )?;
        let delta_f = core::array::from_fn(|i| f[i] - old.deformation[i]);
        let work = finite(dot(&returned.first_piola, &delta_f))?;
        let plastic = nonnegative(returned.integrity * returned.native.dissipation)?;
        let damage = nonnegative(returned.damage_increment * old.active_stored)?;
        let compression_piola = returned
            .native
            .inverse_transpose
            .map(|x| returned.compression_pressure * x);
        let native_defect = finite(
            dot(&returned.native.first_piola, &delta_f)
                - (returned.native.stored - old.native_stored)
                - returned.native.dissipation,
        )?;
        let compression_defect = finite(
            dot(&compression_piola, &delta_f)
                - (returned.compression_stored - old.compression_stored),
        )?;
        let defect = finite(
            returned.integrity * native_defect + (1.0 - returned.integrity) * compression_defect,
        )?;
        let residual = finite(work - (returned.stored - old.stored) - plastic - damage - defect)?;
        if residual.abs() > self.settings.energy_tolerances.energy_j_m3 {
            return Err(FiniteDuctileError::EnergyDefect(residual));
        }
        let energy = FiniteDuctileEnergy {
            initial_stored_j_m3: old.stored,
            final_stored_j_m3: returned.stored,
            mechanical_work_j_m3: work,
            plastic_dissipation_j_m3: plastic,
            damage_release_j_m3: damage,
            endpoint_work_defect_j_m3: defect,
            residual_j_m3: residual,
            equivalent_plastic_work_increment_j_m2: nonnegative(
                self.settings.characteristic_length_m * returned.equivalent_work,
            )?,
            degraded_equivalent_plastic_work_increment_j_m2: nonnegative(
                self.settings.post_initiation_plastic_work_j_m2 * returned.damage_increment,
            )?,
        };
        let state = FiniteDuctileState {
            owner: self,
            deformation: *f,
            native: returned.native.state,
            initiation: returned.initiation,
            exponent: returned.exponent,
            stored: returned.stored,
            native_stored: returned.native.stored,
            active_stored: returned.active_stored,
            compression_stored: returned.compression_stored,
            first_piola: returned.first_piola,
            kirchhoff: returned.kirchhoff,
            plastic_dissipated: nonnegative(old.plastic_dissipated + plastic)?,
            damage_dissipated: nonnegative(old.damage_dissipated + damage)?,
            endpoint_defect: finite(old.endpoint_defect + defect)?,
        };
        Ok(FiniteDuctileStep {
            first_piola_pa: returned.first_piola,
            kirchhoff_pa: returned.kirchhoff,
            tangent_pa,
            state,
            energy,
        })
    }

    fn return_map(
        &self,
        f: &[f64; 9],
        old: &FiniteDuctileState<'_, L>,
        conditions: DuctileJ2Conditions,
    ) -> Result<Returned, FiniteDuctileError> {
        // Native hardening tables are immutable Arc-backed arrays: preparing
        // this interval neither allocates nor copies a measured curve. The
        // existing rate owner retains its scalar root, metric and work semantics.
        let kinetics = self.overstress.map(|parameters| {
            FiniteOverstressJ2::new(self.native.clone(), parameters, conditions.duration_s)
        }).transpose()?;
        let native = self.native.return_map_with_rate(f, &old.native, kinetics.as_ref())?;
        let viscous_modulus_pa = if native.dgamma > 0.0 {
            kinetics.as_ref().map(|law| {
                law.endpoint_overstress_pa(native.dgamma)
                    .map(|stress| (stress / native.dgamma) * (1.0 / law.parameters().exponent))
            }).transpose()?.unwrap_or(0.0)
        } else {
            0.0
        };
        nonnegative(viscous_modulus_pa)?;
        let alpha = native.state.alpha;
        let delta_alpha = nonnegative(alpha - old.native.alpha)?;
        let mut initiation = old.initiation;
        let mut onset_alpha = alpha;
        let mut onset = None;
        if old.initiation.initiated() {
            onset_alpha = old.native.alpha;
        } else {
            let step = old.initiation.trial(
                &self.initiation,
                alpha,
                conditions.duration_s,
                conditions.temperature_k,
                &pack(&native.tau),
            )?;
            initiation = step.state;
            if initiation.utilization() > 1.0 {
                let remaining = 1.0 - old.initiation.utilization();
                onset_alpha = finite(
                    old.native.alpha
                        + remaining
                            * step.fracture_strain.ok_or(FiniteDuctileError::Invalid(
                                "initiation crossing needs plastic flow",
                            ))?,
                )?;
                onset = Some((
                    DuctileStressState {
                        triaxiality: step
                            .triaxiality
                            .ok_or(FiniteDuctileError::Invalid("missing onset triaxiality"))?,
                        lode_parameter: step
                            .lode_parameter
                            .ok_or(FiniteDuctileError::Invalid("missing onset Lode parameter"))?,
                    },
                    step.strain_rate_s_inv,
                    remaining,
                ));
            }
        }
        let post_alpha = (alpha - onset_alpha).max(0.0);
        if post_alpha > delta_alpha * (1.0 + 64.0 * f64::EPSILON) {
            return Err(FiniteDuctileError::Invalid(
                "onset exceeds the native plastic increment",
            ));
        }
        let p = self.native.parameters;
        let mut equivalent_work = post_alpha
            * (p.yield_stress_pa + p.hardening_modulus_pa * (0.5 * onset_alpha + 0.5 * alpha));
        if let Some(table) = self.native.tabulated_hardening() {
            equivalent_work += table.hardening_energy_increment_j_m3(onset_alpha, post_alpha)
                .map_err(FiniteJ2Error::Tabulated)?;
        }
        if let Some(voce) = self.native.voce {
            // Exact translated primitive. No subtraction of large hardening
            // endpoint energies; the native small-increment primitive is reused.
            equivalent_work += post_alpha * voce.strength(onset_alpha)
                + det::exp(-voce.rate_per_strain() * onset_alpha) * voce.stored_energy(post_alpha);
        }
        nonnegative(equivalent_work)?;
        let increment = nonnegative(self.work_scale * equivalent_work)?;
        let exponent = nonnegative(old.exponent + increment)?;
        if equivalent_work > 0.0 && (increment == 0.0 || exponent == old.exponent) {
            return Err(FiniteDuctileError::Invalid(
                "damage increment is below retained-history resolution",
            ));
        }
        finite(self.settings.post_initiation_plastic_work_j_m2 * exponent)?;
        let integrity = positive(det::exp(-exponent))?;
        let damage_increment = nonnegative(old.integrity() * -det::expm1(-increment))?;
        if damage_increment > self.settings.maximum_damage_increment {
            return Err(FiniteDuctileError::DamageIncrement {
                proposed: damage_increment,
                maximum: self.settings.maximum_damage_increment,
            });
        }
        // Native log(J) preserves sub-ulp volume strain even if reported J=1.
        let compression_log = native.trial.log_j.min(0.0);
        let compression_pressure = finite(p.bulk_modulus_pa * compression_log)?;
        let compression_stored = nonnegative(0.5 * compression_pressure * compression_log)?;
        let returned_deviator = native.trial.deviator.map(|x| native.scale * x);
        let tension_log = native.trial.log_j.max(0.0);
        let active_stored = nonnegative(
            p.shear_modulus_pa * dot(&returned_deviator, &returned_deviator)
                + self.native.hardening_energy(alpha)
                + 0.5 * p.bulk_modulus_pa * tension_log * tension_log,
        )?;
        let stored =
            nonnegative(integrity * native.stored + (1.0 - integrity) * compression_stored)?;
        let first_piola = core::array::from_fn(|i| {
            integrity * native.first_piola[i]
                + (1.0 - integrity) * compression_pressure * native.inverse_transpose[i]
        });
        let kirchhoff = core::array::from_fn(|i| {
            integrity * native.tau[i]
                + if i == 0 || i == 4 || i == 8 {
                    (1.0 - integrity) * compression_pressure
                } else {
                    0.0
                }
        });
        finite_tensor(&first_piola)?;
        finite_tensor(&kirchhoff)?;
        Ok(Returned {
            native,
            viscous_modulus_pa,
            initiation,
            onset_alpha,
            onset,
            post_alpha,
            equivalent_work,
            exponent,
            integrity,
            damage_increment,
            compression_pressure,
            compression_stored,
            active_stored,
            stored,
            first_piola,
            kirchhoff,
        })
    }

    fn degraded_tangent(
        &self,
        returned: &ReturnMap,
        native: &[[f64; 9]; 9],
        g: f64,
        exponent_gradient: &[f64; 9],
    ) -> Result<[[f64; 9]; 9], FiniteDuctileError> {
        let inverse = &returned.inverse_transpose;
        let bulk = self.native.parameters.bulk_modulus_pa;
        let pressure = bulk * returned.trial.log_j.min(0.0);
        let active_bulk = if returned.trial.log_j < 0.0 {
            bulk
        } else {
            0.0
        };
        let mut tangent = [[0.0; 9]; 9];
        for ij in 0..9 {
            let (i, j) = (ij / 3, ij % 3);
            let active_stress = returned.first_piola[ij] - pressure * inverse[ij];
            for kl in 0..9 {
                let (k, l) = (kl / 3, kl % 3);
                let compression = active_bulk * inverse[ij] * inverse[kl]
                    - pressure * inverse[3 * i + l] * inverse[3 * k + j];
                tangent[ij][kl] = finite(
                    g * native[ij][kl] + (1.0 - g) * compression
                        - active_stress * (g * exponent_gradient[kl]),
                )?;
            }
        }
        Ok(tangent)
    }

    fn replay<'a>(
        &'a self,
        old: &FiniteDuctileState<'a, L>,
    ) -> Result<FiniteDuctileStep<'a, L>, FiniteDuctileError> {
        // Exact holds take the elastic frozen-metric branch. Reconstructing a
        // metric must not generate a new plastic clock or change retained work.
        let trial = self
            .native
            .elastic_endpoint(&old.deformation, &old.native)?;
        let q =
            det::sqrt(6.0) * self.native.parameters.shear_modulus_pa * super::norm(&trial.deviator);
        let returned = self.native.finish_return(
            &old.deformation,
            &old.native,
            trial,
            q,
            RadialReturn {
                plastic_increment: 0.0,
                creep_increment: 0.0,
                scale: 1.0,
                tangent_correction: 0.0,
                plastic_dissipation: 0.0,
                creep_dissipation: 0.0,
                primary_potential: 0.0,
                rate_potential: 0.0,
            },
        )?;
        let native = self
            .native
            .return_tangent(&old.deformation, &old.native, &returned)?;
        let tangent_pa = self.degraded_tangent(&returned, &native, old.integrity(), &[0.0; 9])?;
        Ok(FiniteDuctileStep {
            first_piola_pa: old.first_piola,
            kirchhoff_pa: old.kirchhoff,
            tangent_pa,
            state: old.clone(),
            energy: FiniteDuctileEnergy {
                initial_stored_j_m3: old.stored,
                final_stored_j_m3: old.stored,
                mechanical_work_j_m3: 0.0,
                plastic_dissipation_j_m3: 0.0,
                damage_release_j_m3: 0.0,
                endpoint_work_defect_j_m3: 0.0,
                residual_j_m3: 0.0,
                equivalent_plastic_work_increment_j_m2: 0.0,
                degraded_equivalent_plastic_work_increment_j_m2: 0.0,
            },
        })
    }
}

impl<L> FiniteDuctileState<'_, L> {
    pub fn is_valid(&self) -> bool {
        [
            self.exponent,
            self.stored,
            self.native_stored,
            self.active_stored,
            self.compression_stored,
            self.plastic_dissipated,
            self.damage_dissipated,
            self.native.alpha,
        ]
        .iter()
        .all(|v| v.is_finite() && *v >= 0.0)
            && self.endpoint_defect.is_finite()
            && self.native.creep_alpha == 0.0
            && self
                .deformation
                .iter()
                .chain(&self.first_piola)
                .chain(&self.kirchhoff)
                .chain(&self.native.inverse_plastic_cauchy_green)
                .all(|v| v.is_finite())
    }
    pub const fn deformation(&self) -> &[f64; 9] {
        &self.deformation
    }
    pub const fn native_state(&self) -> &FiniteJ2State {
        &self.native
    }
    /// Physical accepted Piola stress. Reading it advances neither rate flow
    /// nor damage, even when the retained native stress is overstressed.
    pub const fn first_piola_pa(&self) -> &[f64; 9] {
        &self.first_piola
    }
    /// Physical accepted Kirchhoff stress with the same no-advance semantics.
    pub const fn kirchhoff_pa(&self) -> &[f64; 9] {
        &self.kirchhoff
    }
    pub const fn initiation(&self) -> &DuctileState {
        &self.initiation
    }
    pub fn integrity(&self) -> f64 {
        det::exp(-self.exponent)
    }
    pub fn damage(&self) -> f64 {
        -det::expm1(-self.exponent)
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
    pub const fn endpoint_work_defect_j_m3(&self) -> f64 {
        self.endpoint_defect
    }
    pub fn equivalent_plastic_work_j_m2(&self) -> f64 {
        self.owner.settings.post_initiation_plastic_work_j_m2 * self.exponent
    }
    pub fn degraded_equivalent_plastic_work_j_m2(&self) -> f64 {
        self.owner.settings.post_initiation_plastic_work_j_m2 * self.damage()
    }
}

fn finite(value: f64) -> Result<f64, FiniteDuctileError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(FiniteDuctileError::Invalid(
            "finite constitutive result required",
        ))
    }
}
fn positive(value: f64) -> Result<f64, FiniteDuctileError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(FiniteDuctileError::Invalid(
            "positive finite constitutive value required",
        ))
    }
}
fn nonnegative(value: f64) -> Result<f64, FiniteDuctileError> {
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(FiniteDuctileError::Invalid(
            "nonnegative finite constitutive value required",
        ))
    }
}
fn finite_tensor(value: &[f64; 9]) -> Result<(), FiniteDuctileError> {
    if value.iter().all(|x| x.is_finite()) {
        Ok(())
    } else {
        Err(FiniteDuctileError::Invalid("finite stress tensor required"))
    }
}
