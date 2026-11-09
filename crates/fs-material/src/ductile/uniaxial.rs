//! Progressive ductility on the original free-sided axial J2 return.
//!
//! This explicitly one-dimensional model degrades ALL native axial storage and
//! stress by the selected equivalent-work profile `g(x)`, in tension and
//! compression. It is not a reduction of the three-dimensional intact-bulk
//! split in [`super::evolution`]. It supplies neither crack closure/contact nor
//! a separated member. Native Prager plasticity and its accumulated path length
//! survive reversal; the original energetic owner supplies plastic and numerical
//! losses. The declared `Gp` is integrated degraded equivalent plastic work per
//! area, not total fracture dissipation or a heat source.
//!
//! The exact discrete identity is `W = delta(g*psi) + g*D +
//! (g_old-g)*psi_old + g*N`. The third term releases previously stored energy;
//! it is kept separate from plastic loss and backward-Euler numerical loss.
//! Characteristic length belongs to the original receiving band. Scaling it
//! alone does not prove mesh-objective localization or an experimentally
//! calibrated failure path.
//!
//! Immutable models and histories share their original binding. Cloning a
//! model or preparing a new physical interval preserves history ownership;
//! selecting another evolution profile creates a new binding. Trials allocate
//! no memory and never mutate accepted state. A whole-section zero tangent is
//! a structural mechanism for the equilibrium owner, not permission to insert
//! a residual stiffness floor or erase material inventory.

use std::sync::Arc;

use fs_evidence::{Ambition, ModelCard, ValidityDomain};

use super::evolution::{
    DifferentiableDuctileInitiationLaw, DuctileEvolutionError, DuctileJ2Conditions,
    DuctileJ2Settings, DuctileJ2StepEnergy, DuctileWorkEvolution,
};
use super::{DuctileState, DuctileStressState};
use crate::MaterialAdmissibility;
use crate::fiber::Uniaxial;
use crate::plastic::reduction::{UniaxialJ2, UniaxialJ2State};

#[derive(Debug)]
struct Binding {
    native: UniaxialJ2,
    settings: DuctileJ2Settings,
    evolution: DuctileWorkEvolution,
    work_scale: f64,
}

// Equal coefficients do not authorize transplanting a committed damage history
// into an independently constructed material. Clones share this exact binding.
impl PartialEq for Binding {
    fn eq(&self, other: &Self) -> bool {
        core::ptr::eq(self, other)
    }
}

#[derive(Debug, PartialEq)]
struct Model<L> {
    binding: Arc<Binding>,
    initiation: L,
}

/// Immutable native axial plasticity, initiation calibration and band length.
#[derive(Debug, PartialEq)]
pub struct UniaxialDuctileJ2<L> {
    model: Arc<Model<L>>,
}

impl<L> Clone for UniaxialDuctileJ2<L> {
    fn clone(&self) -> Self {
        Self {
            model: Arc::clone(&self.model),
        }
    }
}

/// Owned accepted axial history. Its original model binding cannot be replaced.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialDuctileJ2State {
    binding: Arc<Binding>,
    strain: f64,
    stress_pa: f64,
    native: UniaxialJ2State,
    initiation: DuctileState,
    work_coordinate: f64,
    stored: f64,
    plastic_dissipated: f64,
    damage_dissipated: f64,
    numerical_dissipated: f64,
    last_step_energy: DuctileJ2StepEnergy,
}

/// Disposable complete proposal with an analytic, possibly negative tangent.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialDuctileJ2Step {
    /// Physical axial Cauchy stress in the declared small-strain model [Pa].
    pub stress_pa: f64,
    /// Exact smooth-branch derivative at the same incoming accepted history.
    pub tangent_pa: f64,
    pub state: UniaxialDuctileJ2State,
    pub energy: DuctileJ2StepEnergy,
}

/// Owned adapter for the existing native fiber interface and one physical step.
#[derive(Debug, PartialEq)]
pub struct UniaxialDuctileJ2Law<L> {
    model: UniaxialDuctileJ2<L>,
    conditions: DuctileJ2Conditions,
}

impl<L> Clone for UniaxialDuctileJ2Law<L> {
    fn clone(&self) -> Self {
        Self {
            model: self.model.clone(),
            conditions: self.conditions,
        }
    }
}

impl<L: DifferentiableDuctileInitiationLaw> UniaxialDuctileJ2<L> {
    pub fn new(
        native: UniaxialJ2,
        initiation: L,
        settings: DuctileJ2Settings,
    ) -> Result<Self, DuctileEvolutionError> {
        native.backbone().validate_energy_law()?;
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
        let work_scale = positive(
            settings.characteristic_length_m / settings.post_initiation_plastic_work_j_m2,
        )?;
        Ok(Self {
            model: Arc::new(Model {
                binding: Arc::new(Binding {
                    native,
                    settings,
                    evolution: DuctileWorkEvolution::Exponential,
                    work_scale,
                }),
                initiation,
            }),
        })
    }

    /// Select an authored shape on a new binding. Previously lent histories
    /// retain their original profile and cannot enter the new material.
    pub fn with_evolution(self, evolution: DuctileWorkEvolution) -> Self
    where
        L: Clone,
    {
        Self {
            model: Arc::new(Model {
                binding: Arc::new(Binding {
                    native: self.native().clone(),
                    settings: self.settings(),
                    evolution,
                    work_scale: self.model.binding.work_scale,
                }),
                initiation: self.initiation_law().clone(),
            }),
        }
    }

    pub fn native(&self) -> &UniaxialJ2 {
        &self.model.binding.native
    }
    pub fn initiation_law(&self) -> &L {
        &self.model.initiation
    }
    pub fn settings(&self) -> DuctileJ2Settings {
        self.model.binding.settings
    }
    pub fn evolution(&self) -> DuctileWorkEvolution {
        self.model.binding.evolution
    }
    pub fn same_model(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.model, &other.model)
    }
    pub fn owns_state(&self, state: &UniaxialDuctileJ2State) -> bool {
        Arc::ptr_eq(&self.model.binding, &state.binding)
    }

    pub fn initial_state(&self) -> UniaxialDuctileJ2State {
        UniaxialDuctileJ2State {
            binding: Arc::clone(&self.model.binding),
            strain: 0.0,
            stress_pa: 0.0,
            native: self.native().initial_state(),
            initiation: DuctileState::virgin(),
            work_coordinate: 0.0,
            stored: 0.0,
            plastic_dissipated: 0.0,
            damage_dissipated: 0.0,
            numerical_dissipated: 0.0,
            last_step_energy: zero_energy(0.0),
        }
    }

    pub fn prepared(
        &self,
        conditions: DuctileJ2Conditions,
    ) -> Result<UniaxialDuctileJ2Law<L>, DuctileEvolutionError> {
        self.check_conditions(conditions)?;
        Ok(UniaxialDuctileJ2Law {
            model: self.clone(),
            conditions,
        })
    }

    fn check_conditions(
        &self,
        conditions: DuctileJ2Conditions,
    ) -> Result<(), DuctileEvolutionError> {
        positive(conditions.duration_s)?;
        positive(conditions.temperature_k)?;
        self.initiation_law()
            .domain()
            .check_temperature(conditions.temperature_k)?;
        Ok(())
    }

    /// One complete native axial trial and its original energetic audit.
    /// Temperature and duration stay fixed through all equilibrium iterations.
    pub fn trial(
        &self,
        strain: f64,
        old: &UniaxialDuctileJ2State,
        conditions: DuctileJ2Conditions,
    ) -> Result<UniaxialDuctileJ2Step, DuctileEvolutionError> {
        if !self.owns_state(old) {
            return Err(DuctileEvolutionError::StateMismatch);
        }
        self.check_conditions(conditions)?;
        if !old.is_valid() || !strain.is_finite() {
            return Err(DuctileEvolutionError::Invalid(
                "invalid axial strain or accepted history",
            ));
        }
        if strain == old.strain {
            return self.replay(old);
        }
        let response = self.native().response(strain, &old.native);
        let native = self.native().advance_with_energy(
            old.strain,
            &old.native,
            strain,
            self.settings().energy_tolerances,
        )?;
        let alpha = response.state.equivalent_plastic_strain();
        let old_alpha = old.native.equivalent_plastic_strain();
        let delta_alpha = nonnegative(alpha - old_alpha)?;
        let backbone = self.native().backbone();
        let e = backbone.elastic.youngs;
        let trial = e * (strain - old.native.plastic_strain())
            - backbone.kinematic_hardening * old.native.plastic_strain();
        let alpha_gradient = if delta_alpha > 0.0 {
            trial.signum() * (e / (e + (backbone.hardening + backbone.kinematic_hardening)))
        } else {
            0.0
        };
        let mut initiation = old.initiation;
        let mut onset_alpha = alpha;
        let mut onset_gradient = 0.0;
        if old.initiation.initiated() {
            onset_alpha = old_alpha;
        } else {
            let step = old.initiation.trial(
                self.initiation_law(),
                alpha,
                conditions.duration_s,
                conditions.temperature_k,
                &[response.stress_pa, 0.0, 0.0, 0.0, 0.0, 0.0],
            )?;
            initiation = step.state;
            if initiation.utilization() > 1.0 {
                let fracture = step.fracture_strain.ok_or(DuctileEvolutionError::Invalid(
                    "axial onset requires positive native plastic flow",
                ))?;
                let remaining = 1.0 - old.initiation.utilization();
                onset_alpha = finite(old_alpha + remaining * fracture)?;
                let stress = DuctileStressState {
                    triaxiality: step.triaxiality.ok_or(DuctileEvolutionError::Invalid(
                        "missing physical axial triaxiality",
                    ))?,
                    lode_parameter: step.lode_parameter.ok_or(DuctileEvolutionError::Invalid(
                        "missing physical axial Lode parameter",
                    ))?,
                };
                let gradient = self.initiation_law().fracture_strain_gradient_at(
                    stress,
                    step.strain_rate_s_inv,
                    conditions.temperature_k,
                )?;
                // eta=+/-1/3 and xi=+/-1 are constant on each nonzero
                // physical axial-stress branch; only rate moves this onset.
                onset_gradient = finite(
                    remaining * gradient.strain_rate_s_inv * alpha_gradient / conditions.duration_s,
                )?;
            }
        }
        let post_alpha = (alpha - onset_alpha).max(0.0);
        if post_alpha > delta_alpha * (1.0 + 64.0 * f64::EPSILON) {
            return Err(DuctileEvolutionError::Invalid(
                "axial onset exceeds the native plastic increment",
            ));
        }
        let equivalent_work = nonnegative(
            post_alpha
                * (backbone.yield_stress + backbone.hardening * (0.5 * onset_alpha + 0.5 * alpha)),
        )?;
        let increment = nonnegative(self.model.binding.work_scale * equivalent_work)?;
        let coordinate = nonnegative(old.work_coordinate + increment)?;
        finite(self.settings().post_initiation_plastic_work_j_m2 * coordinate)?;
        let profile = self.evolution();
        let integrity = match profile {
            DuctileWorkEvolution::Exponential => positive(profile.integrity(coordinate))?,
            DuctileWorkEvolution::FiniteSupportQuadratic => {
                nonnegative(profile.integrity(coordinate))?
            }
        };
        let profile_increment = profile.increment(old.work_coordinate, increment);
        let damage_increment = nonnegative(profile_increment.damage)?;
        if damage_increment > self.settings().maximum_damage_increment {
            return Err(DuctileEvolutionError::DamageIncrement {
                proposed: damage_increment,
                maximum: self.settings().maximum_damage_increment,
            });
        }
        let degradation_slope = profile.degradation_slope(coordinate);
        // Exact finite-support failure has no active derivative. Do not form
        // a potentially unrepresentable derivative only to multiply it by zero.
        let coordinate_gradient = if post_alpha > 0.0 && degradation_slope > 0.0 {
            let resistance = |a| backbone.yield_stress + backbone.hardening * a;
            self.model.binding.work_scale
                * (resistance(alpha) * alpha_gradient - resistance(onset_alpha) * onset_gradient)
        } else {
            0.0
        };
        let stress_pa = finite(integrity * response.stress_pa)?;
        let tangent_pa = finite(
            integrity * response.tangent_pa
                - response.stress_pa * degradation_slope * coordinate_gradient,
        )?;
        let stored = nonnegative(integrity * native.energy.final_state.total_j_m3())?;
        let work = finite(stress_pa * (strain - old.strain))?;
        let plastic = nonnegative(integrity * native.energy.physical_dissipation_j_m3)?;
        let damage = nonnegative(damage_increment * native.energy.initial.total_j_m3())?;
        let numerical = nonnegative(integrity * native.energy.numerical_dissipation_j_m3)?;
        let residual = finite(work - (stored - old.stored) - plastic - damage - numerical)?;
        if residual.abs() > self.settings().energy_tolerances.energy_j_m3 {
            return Err(DuctileEvolutionError::EnergyDefect(residual));
        }
        let energy = DuctileJ2StepEnergy {
            initial_stored_j_m3: old.stored,
            final_stored_j_m3: stored,
            mechanical_work_j_m3: work,
            plastic_dissipation_j_m3: plastic,
            damage_release_j_m3: damage,
            numerical_dissipation_j_m3: numerical,
            residual_j_m3: residual,
            equivalent_plastic_work_increment_j_m2: nonnegative(
                self.settings().characteristic_length_m * equivalent_work,
            )?,
            degraded_equivalent_plastic_work_increment_j_m2: nonnegative(
                self.settings().post_initiation_plastic_work_j_m2 * profile_increment.degraded_work,
            )?,
        };
        Ok(UniaxialDuctileJ2Step {
            stress_pa,
            tangent_pa,
            state: UniaxialDuctileJ2State {
                binding: Arc::clone(&self.model.binding),
                strain,
                stress_pa,
                native: response.state,
                initiation,
                work_coordinate: coordinate,
                stored,
                plastic_dissipated: nonnegative(old.plastic_dissipated + plastic)?,
                damage_dissipated: nonnegative(old.damage_dissipated + damage)?,
                numerical_dissipated: nonnegative(old.numerical_dissipated + numerical)?,
                last_step_energy: energy,
            },
            energy,
        })
    }

    fn replay(
        &self,
        old: &UniaxialDuctileJ2State,
    ) -> Result<UniaxialDuctileJ2Step, DuctileEvolutionError> {
        let e = self.native().backbone().elastic.youngs;
        let energy = zero_energy(old.stored);
        let mut state = old.clone();
        state.last_step_energy = energy;
        Ok(UniaxialDuctileJ2Step {
            stress_pa: old.stress_pa,
            tangent_pa: finite(old.integrity() * e)?,
            state,
            energy,
        })
    }
}

impl UniaxialDuctileJ2State {
    pub fn is_valid(&self) -> bool {
        [
            self.strain,
            self.stress_pa,
            self.work_coordinate,
            self.stored,
            self.plastic_dissipated,
            self.damage_dissipated,
            self.numerical_dissipated,
        ]
        .iter()
        .all(|x| x.is_finite())
            && self.work_coordinate >= 0.0
            && self.stored >= 0.0
    }
    pub const fn strain(&self) -> f64 {
        self.strain
    }
    pub const fn stress_pa(&self) -> f64 {
        self.stress_pa
    }
    pub const fn native_state(&self) -> &UniaxialJ2State {
        &self.native
    }
    pub const fn equivalent_plastic_strain(&self) -> f64 {
        self.native.equivalent_plastic_strain()
    }
    /// Original accepted increment; an unchanged-strain hold contains zero work.
    pub const fn last_step_energy(&self) -> &DuctileJ2StepEnergy {
        &self.last_step_energy
    }
    /// Frozen after the accepted initiation crossing; native alpha keeps growing.
    pub const fn initiation(&self) -> &DuctileState {
        &self.initiation
    }
    pub fn integrity(&self) -> f64 {
        self.binding.evolution.integrity(self.work_coordinate)
    }
    pub fn damage(&self) -> f64 {
        self.binding.evolution.damage(self.work_coordinate)
    }
    pub fn is_fully_failed(&self) -> bool {
        self.is_valid() && self.binding.evolution.is_fully_failed(self.work_coordinate)
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
    pub fn equivalent_plastic_work_j_m2(&self) -> f64 {
        self.binding.settings.post_initiation_plastic_work_j_m2 * self.work_coordinate
    }
    pub fn degraded_equivalent_plastic_work_j_m2(&self) -> f64 {
        self.binding.settings.post_initiation_plastic_work_j_m2
            * self.binding.evolution.degraded_work(self.work_coordinate)
    }
}

impl<L: DifferentiableDuctileInitiationLaw> UniaxialDuctileJ2Law<L> {
    pub const fn model(&self) -> &UniaxialDuctileJ2<L> {
        &self.model
    }
    pub const fn conditions(&self) -> DuctileJ2Conditions {
        self.conditions
    }
    pub fn with_conditions(
        &self,
        conditions: DuctileJ2Conditions,
    ) -> Result<Self, DuctileEvolutionError> {
        self.model.prepared(conditions)
    }
    pub fn initial_state(&self) -> UniaxialDuctileJ2State {
        self.model.initial_state()
    }
    pub fn response(
        &self,
        strain: f64,
        old: &UniaxialDuctileJ2State,
    ) -> Result<UniaxialDuctileJ2Step, DuctileEvolutionError> {
        self.model.trial(strain, old, self.conditions)
    }
}

impl<L: DifferentiableDuctileInitiationLaw> Uniaxial for UniaxialDuctileJ2Law<L> {
    type State = UniaxialDuctileJ2State;
    fn initial_state(&self) -> Self::State {
        self.model.initial_state()
    }
    fn stress(&self, strain: f64, state: &Self::State) -> f64 {
        self.response(strain, state)
            .map_or(f64::NAN, |v| v.stress_pa)
    }
    fn tangent(&self, strain: f64, state: &Self::State) -> f64 {
        self.response(strain, state)
            .map_or(f64::NAN, |v| v.tangent_pa)
    }
    fn update_state(&self, strain: f64, state: &Self::State) -> Self::State {
        self.response(strain, state)
            .map(|v| v.state)
            .unwrap_or_else(|_| {
                let mut invalid = self.model.initial_state();
                invalid.strain = f64::NAN;
                invalid
            })
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            has_stored_energy: true,
            dissipation_nonnegative: true,
            polyconvex: None,
            tangent_symmetric: true,
            failure_envelope: "native axial J2 with explicit length-scaled work degradation; both stress signs degrade; no crack contact or deletion",
        }
    }
    fn card(&self) -> ModelCard {
        ModelCard::new(
            "material.ductile-uniaxial-j2-work", "0.1.0", Ambition::Solid,
            vec![
                "original small-strain axial J2 with energetic isotropic and Prager hardening".into(),
                "endpoint initiation and explicit characteristic-length equivalent-work attenuation".into(),
                "uniform axial energy degradation in tension and compression".into(),
            ],
            ValidityDomain::unconstrained().with(
                "strain-magnitude", 0.0,
                self.model.native().backbone().elastic.strain_limit.min(0.05),
            ),
            vec![
                "no necking, finite strain, crack closure, separated members or experimental qualification".into(),
                "Gp is degraded equivalent plastic work per area, not total fracture loss or heat".into(),
                "negative and zero tangent require an appropriate native structural continuation owner".into(),
            ],
            f64::INFINITY,
        )
    }
}

fn finite(value: f64) -> Result<f64, DuctileEvolutionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DuctileEvolutionError::Invalid(
            "axial ductile arithmetic is not representable",
        ))
    }
}
fn zero_energy(stored: f64) -> DuctileJ2StepEnergy {
    DuctileJ2StepEnergy {
        initial_stored_j_m3: stored,
        final_stored_j_m3: stored,
        mechanical_work_j_m3: 0.0,
        plastic_dissipation_j_m3: 0.0,
        damage_release_j_m3: 0.0,
        numerical_dissipation_j_m3: 0.0,
        residual_j_m3: 0.0,
        equivalent_plastic_work_increment_j_m2: 0.0,
        degraded_equivalent_plastic_work_increment_j_m2: 0.0,
    }
}
fn nonnegative(value: f64) -> Result<f64, DuctileEvolutionError> {
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(DuctileEvolutionError::Invalid(
            "axial ductile quantity must be finite and nonnegative",
        ))
    }
}
fn positive(value: f64) -> Result<f64, DuctileEvolutionError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(DuctileEvolutionError::Invalid(
            "axial ductile quantity must be finite and positive",
        ))
    }
}
