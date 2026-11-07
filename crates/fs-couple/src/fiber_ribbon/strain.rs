//! The existing extensional strain owner in the original fiber modal basis.
//!
//! Normal deflection changes `u_x + y_x²/2`, hence both local tension and axial
//! motion. Each width cell keeps its original share of mass, EA and preload.
//! This is the moderate-slope, fixed-end rod limit: omitted axial harmonics
//! equilibrate statically, and the retained axial/normal modes must be refined
//! together. The original axial viscosity acts on this same complete strain
//! rate, including normal/axial cross terms. Its linear axial damping is removed
//! once; bending viscosity remains independent. Thermal evolution is absent.

use super::{FiberRibbonError, TensionedFiberRibbon};
use crate::render::schedule::force::coupled::ModalCouplingError;
use crate::render::schedule::force::coupled::contact::friction::nonlinear::{
    MappedModalPotentialBinding, ModalCoordinate, ModalPotential, NonlinearModalConfig,
};
use fs_math::det;
use fs_nlmodal::{KcStringParams, longitudinal::LongitudinalString};
use fs_phs::Storage;

pub(super) mod relaxation;
/// Source-bound elastic tensile admission on the original fiber continuum.
pub mod strength;

/// Continuous-span constitutive bounds, independent of numerical resolution.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiberStrainLimits {
    /// Conservative absolute normal-slope ceiling in (0, 0.3].
    pub maximum_slope: f64,
    /// Conservative additional material-strain ceiling in (0, 0.1].
    pub maximum_strain: f64,
}

impl FiberStrainLimits {
    /// Admit explicit motion bounds before reading material sources or preparing modes.
    pub fn validate(self) -> Result<(), FiberRibbonError> {
        if !self.maximum_slope.is_finite()
            || !(0.0 < self.maximum_slope && self.maximum_slope <= 0.3)
            || !self.maximum_strain.is_finite()
            || !(0.0 < self.maximum_strain && self.maximum_strain <= 0.1)
        {
            return Err(FiberRibbonError::Invalid(
                "fiber strain requires slope in (0, 0.3] and strain in (0, 0.1]",
            ));
        }
        Ok(())
    }
}

/// Physical strain and equilibrium-elastic tension of one original strip.
/// Additional Maxwell stress/storage requires the joint owner's accepted memory
/// and is not included in this displacement-only observation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiberStrainObservation {
    /// Additional mean axial material strain, dimensionless.
    pub mean_dynamic_strain: f64,
    /// Conservative absolute strain bound over the complete fiber span.
    pub maximum_abs_dynamic_strain: f64,
    /// Lower bound of equilibrium-elastic tensile force [N], including preload.
    pub minimum_tension_n: f64,
    /// Upper bound of equilibrium-elastic tensile force [N], including preload.
    pub maximum_tension_n: f64,
    /// Equilibrium extensional energy, including linear axial storage [J].
    pub strain_energy_j: f64,
}

/// Source-derived extension in axial-first, normal-second mass coordinates.
/// The coupled owner subtracts the existing axial quadratic storage exactly
/// once. No second mass, longitudinal oscillator, or nonlinear solver is added.
#[derive(Clone, Debug)]
pub struct FiberStrain {
    model: LongitudinalString,
    slope_scales: Vec<f64>,
    tension_n: f64,
    axial_rigidity_n: f64,
    /// C_a L for this original strip [J s], from the same axial source bundle.
    viscous_coefficient_j_s: f64,
    limits: FiberStrainLimits,
    config: NonlinearModalConfig,
    tensile: Option<std::sync::Arc<strength::FiberTensileGuard>>,
}

impl TensionedFiberRibbon {
    /// Couple the already retained axial and normal waves through material strain.
    ///
    /// The unchanged material section supplies each strip's EA and mass. Normal
    /// deformation can now drive longitudinal waves and change the same contact
    /// reaction and pressure history. When attached to a string bank, this law
    /// joins its original exact-linear/discrete-gradient contact solve. Existing
    /// nonlinear constituents must use this SAME explicit solve configuration.
    /// The selected axial viscosity follows the complete material strain rate;
    /// its original linear modal damping is replaced, never counted twice.
    ///
    /// No modes are introduced or silently discarded. This fixed-end,
    /// moderate-slope rung refuses loss of local tension, excess strain, and
    /// nonconvergence; it does not model stick camber or migrating contact.
    pub fn with_axial_transverse_strain(
        mut self,
        limits: FiberStrainLimits,
        config: NonlinearModalConfig,
    ) -> Result<Self, FiberRibbonError> {
        limits.validate()?;
        if self.strain.is_some() || config.validate().is_err() {
            return Err(FiberRibbonError::Invalid(
                "fiber strain needs one cold binding and finite slope, strain and solve bounds",
            ));
        }
        let share = self.geometry.fiber_count as f64 / self.resolution.cells as f64;
        let tension_n = self.geometry.tension_n / self.resolution.cells as f64;
        let axial_rigidity_n = self.axial.elastic().axial_stiffness_n() * share;
        let model = LongitudinalString::new(
            KcStringParams {
                length: self.geometry.length_m,
                lin_density: self.cell_line_density_kg_m,
                tension: tension_n,
                ea: axial_rigidity_n,
            },
            self.resolution.transverse_modes,
            1,
            &(1..=self.resolution.axial_modes).collect::<Vec<_>>(),
        )
        .map_err(FiberRibbonError::Strain)?;
        if model
            .axial_omegas()
            .iter()
            .zip(self.cell.modes())
            .any(|(a, b)| {
                (a - b.angular_frequency_rad_s).abs()
                    > 32.0 * f64::EPSILON * a.abs().max(b.angular_frequency_rad_s.abs())
            })
        {
            return Err(FiberRibbonError::Invalid(
                "fiber strain must retain the original material and modal mass normalization",
            ));
        }
        let normalization = det::sqrt(2.0 / (self.cell_line_density_kg_m * self.geometry.length_m));
        let slope_scales: Vec<_> = (1..=self.resolution.transverse_modes)
            .map(|n| normalization * n as f64 * core::f64::consts::PI / self.geometry.length_m)
            .collect();
        if slope_scales.iter().any(|x| !x.is_finite() || *x <= 0.0) {
            return Err(FiberRibbonError::Invalid(
                "fiber strain map is unrepresentable",
            ));
        }
        let viscosity = self.axial.viscous_rigidity_n_s();
        let viscous_coefficient_j_s = viscosity * share * self.geometry.length_m;
        if !viscous_coefficient_j_s.is_finite()
            || (viscosity > 0.0 && viscous_coefficient_j_s == 0.0)
        {
            return Err(FiberRibbonError::Invalid(
                "fiber strain viscosity is unrepresentable",
            ));
        }
        if viscosity > 0.0 {
            let mut damping = self.cell.damping_rates_per_s().to_vec();
            damping[..self.resolution.axial_modes].fill(0.0);
            // The nonlinear strain owner now supplies axial viscosity; all
            // other original rates, states, budgets and the exact clock remain.
            self.cell = self.cell.with_damping_rates(&damping)
                .map_err(FiberRibbonError::Modal)?;
        }
        let strain = FiberStrain {
            model,
            slope_scales,
            tension_n,
            axial_rigidity_n,
            viscous_coefficient_j_s,
            limits,
            config,
            tensile: None,
        };
        if let Some(memory) = &mut self.relaxation {
            memory.with_finite_strain(
                &strain,
                self.geometry.length_m,
                share,
                self.cell.modes(),
                self.cell.sample_period_s(),
            )?;
        }
        self.strain = Some(strain);
        Ok(self)
    }

    /// Original conservative strain source; absent is the explicitly linear rung.
    pub const fn strain(&self) -> Option<&FiberStrain> {
        self.strain.as_ref()
    }

    pub(crate) fn strain_bindings(
        &self,
        first_component: usize,
    ) -> Option<(Vec<MappedModalPotentialBinding>, NonlinearModalConfig)> {
        self.strain.as_ref().map(|strain| {
            // Atomic ribbon assembly binds the complete guarded history next.
            // A stand-alone public strain clone retains its refusal when that
            // original memory owner is missing; only this private path hands
            // total-stress admission to the actual hereditary law.
            let mut potential = strain.clone();
            if potential.tensile.as_ref().is_some_and(|guard| guard.hereditary) {
                potential.tensile = None;
            }
            let bindings = (0..self.resolution.cells)
                .map(|cell| MappedModalPotentialBinding {
                    coordinates: (0..self.cell.modes().len())
                        .map(|mode| ModalCoordinate {
                            component: first_component + cell,
                            mode,
                        })
                        .collect(),
                    replaces_linear: (0..self.resolution.axial_modes).collect(),
                    potential: Box::new(potential.clone()),
                })
                .collect();
            (bindings, strain.config)
        })
    }
}

fn invalid() -> ModalCouplingError {
    ModalCouplingError::Invalid {
        what: "fiber motion exceeds its finite slope/strain domain or loses local tension",
    }
}

impl FiberStrain {
    pub(in crate::fiber_ribbon) fn require_tensile_history(&mut self, required: bool) {
        if let Some(guard) = &mut self.tensile {
            std::sync::Arc::make_mut(guard).hereditary = required;
        }
    }

    fn validate_kinematics(&self, q: &[f64]) -> Result<(f64, f64), ModalCouplingError> {
        if q.len() != self.coordinate_count() || q.iter().any(|x| !x.is_finite()) {
            return Err(invalid());
        }
        let slopes = self.slopes(q);
        let slope_bound: f64 = slopes.iter().map(|x| x.abs()).sum();
        let (mean, variation) = self
            .model
            .strain_bounds(&slopes, &q[..self.model.axial_modes().len()])
            .map_err(|_| invalid())?;
        self.validate_kinematic_bounds(slope_bound, mean, variation)?;
        Ok((mean, variation))
    }

    fn validate_kinematic_bounds(&self, slope_bound: f64, mean: f64, variation: f64)
        -> Result<(), ModalCouplingError>
    {
        if !slope_bound.is_finite()
            || !mean.is_finite()
            || !variation.is_finite()
            || slope_bound > self.limits.maximum_slope
            || mean.abs() + variation > self.limits.maximum_strain
            || self.tension_n + self.axial_rigidity_n * (mean - variation) <= 0.0
        {
            return Err(invalid());
        }
        Ok(())
    }

    /// Original continuous-span motion bounds.
    pub const fn limits(&self) -> FiberStrainLimits {
        self.limits
    }
    /// Joint constitutive solve limits retained by every participating strip.
    pub const fn config(&self) -> NonlinearModalConfig {
        self.config
    }

    fn slopes(&self, q: &[f64]) -> Vec<f64> {
        let mut slopes = vec![0.0; self.slope_scales.len() + 1];
        for ((value, scale), q) in slopes[1..]
            .iter_mut()
            .zip(&self.slope_scales)
            .zip(&q[self.model.axial_modes().len()..])
        {
            *value = scale * q;
        }
        slopes
    }

    /// Read the supplied original strip coordinates without advancing history.
    /// `q` has all axial coordinates followed by all normal coordinates.
    pub fn observe(&self, q: &[f64]) -> Result<FiberStrainObservation, ModalCouplingError> {
        self.validate_coordinates(q)?;
        let slopes = self.slopes(q);
        let axial = &q[..self.model.axial_modes().len()];
        let (mean, variation) = self
            .model
            .strain_bounds(&slopes, axial)
            .map_err(|_| invalid())?;
        Ok(FiberStrainObservation {
            mean_dynamic_strain: mean,
            maximum_abs_dynamic_strain: mean.abs() + variation,
            minimum_tension_n: self.tension_n + self.axial_rigidity_n * (mean - variation),
            maximum_tension_n: self.tension_n + self.axial_rigidity_n * (mean + variation),
            strain_energy_j: self.model.energy(&slopes, axial).map_err(|_| invalid())?,
        })
    }
}

impl Storage for FiberStrain {
    fn hamiltonian(&self, q: &[f64]) -> f64 {
        if q.len() != self.coordinate_count() {
            return f64::NAN;
        }
        self.model
            .energy(&self.slopes(q), &q[..self.model.axial_modes().len()])
            .unwrap_or(f64::NAN)
    }
    fn gradient(&self, q: &[f64], out: &mut [f64]) {
        if q.len() != self.coordinate_count() || out.len() != q.len() {
            out.fill(f64::NAN);
            return;
        }
        let slopes = self.slopes(q);
        let axial = &q[..self.model.axial_modes().len()];
        let mut gradient = vec![0.0; slopes.len()];
        let (axial_out, normal_out) = out.split_at_mut(axial.len());
        if self
            .model
            .discrete_gradient(&slopes, &slopes, axial, axial, &mut gradient, axial_out)
            .is_err()
        {
            out.fill(f64::NAN);
            return;
        }
        for ((out, scale), derivative) in normal_out
            .iter_mut()
            .zip(&self.slope_scales)
            .zip(&gradient[1..])
        {
            *out = scale * derivative;
        }
    }
}

impl ModalPotential for FiberStrain {
    fn coordinate_count(&self) -> usize {
        self.model.axial_modes().len() + self.slope_scales.len()
    }
    fn validate_coordinates(&self, q: &[f64]) -> Result<(), ModalCouplingError> {
        let (mean, variation) = self.validate_kinematics(q)?;
        if let Some(guard) = &self.tensile {
            // Rate stress can oppose elastic stress. The native endpoint gate
            // owns that combined decision; storage/loss quadrature sees only q.
            if guard.has_viscosity() || guard.hereditary { return Ok(()); }
            let axial_strain = self.tension_n / self.axial_rigidity_n + mean + variation;
            if guard
                .assess(&q[self.model.axial_modes().len()..], axial_strain)?
                .exceeds_limit()
            {
                return Err(ModalCouplingError::Invalid {
                    what: "fiber tensile envelope exceeds its original material limit",
                });
            }
        }
        Ok(())
    }

    fn requires_velocity_validation(&self) -> bool {
        self.tensile.as_ref().is_some_and(|guard| guard.has_viscosity() || guard.hereditary)
    }

    fn validate_state(&self, q: &[f64], velocity: &[f64]) -> Result<(), ModalCouplingError> {
        if self.tensile_state_assessment(q, velocity)?.is_some_and(|assessment| assessment.exceeds_limit()) {
            return Err(ModalCouplingError::Invalid {
                what: "fiber tensile envelope exceeds its original material limit",
            });
        }
        Ok(())
    }

    fn has_constitutive_loss(&self) -> bool {
        self.viscous_coefficient_j_s > 0.0
    }

    fn dissipative_reaction(
        &self,
        q: &[f64],
        velocity: &[f64],
        out: &mut [f64],
    ) -> Result<(), ModalCouplingError> {
        self.validate_coordinates(q)?;
        if velocity.len() != q.len()
            || out.len() != q.len()
            || velocity.iter().any(|v| !v.is_finite())
        {
            return Err(ModalCouplingError::Invalid {
                what: "fiber strain loss needs complete original modal velocities",
            });
        }
        let axial_count = self.model.axial_modes().len();
        let slopes = self.slopes(q);
        let mut reaction = vec![0.0; slopes.len()];
        let (axial_out, normal_out) = out.split_at_mut(axial_count);
        self.model
            .viscous_strain_reaction(
                &slopes,
                &self.slopes(velocity),
                &velocity[..axial_count],
                self.viscous_coefficient_j_s,
                &mut reaction,
                axial_out,
            )
            .map_err(|_| invalid())?;
        for ((force, &scale), &reaction) in normal_out
            .iter_mut()
            .zip(&self.slope_scales)
            .zip(&reaction[1..])
        {
            *force = scale * reaction;
        }
        if out.iter().any(|force| !force.is_finite()) {
            return Err(invalid());
        }
        Ok(())
    }
}
