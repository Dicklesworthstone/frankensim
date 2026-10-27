//! The existing extensional strain owner in the original fiber modal basis.
//!
//! Normal deflection changes `u_x + y_x²/2`, hence both local tension and axial
//! motion. Each width cell keeps its original share of mass, EA and preload.
//! This is the moderate-slope, fixed-end rod limit: omitted axial harmonics
//! equilibrate statically, and the retained axial/normal modes must be refined
//! together. Linear material viscosities remain the original small-strain
//! operators; finite-slope viscous cross terms and thermal evolution are absent.

use super::{FiberRibbonError, TensionedFiberRibbon};
use crate::render::schedule::force::coupled::ModalCouplingError;
use crate::render::schedule::force::coupled::contact::friction::nonlinear::{
    MappedModalPotentialBinding, ModalCoordinate, ModalPotential, NonlinearModalConfig,
};
use fs_math::det;
use fs_nlmodal::{KcStringParams, longitudinal::LongitudinalString};
use fs_phs::Storage;

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

/// Physical strain and total-tension bounds of one original quadrature strip.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiberStrainObservation {
    /// Additional mean axial material strain, dimensionless.
    pub mean_dynamic_strain: f64,
    /// Conservative absolute strain bound over the complete fiber span.
    pub maximum_abs_dynamic_strain: f64,
    /// Lower bound of strip tensile force [N], including installed preload.
    pub minimum_tension_n: f64,
    /// Upper bound of strip tensile force [N], including installed preload.
    pub maximum_tension_n: f64,
    /// Full extensional energy, including retained linear axial storage [J].
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
    limits: FiberStrainLimits,
    config: NonlinearModalConfig,
}

impl TensionedFiberRibbon {
    /// Couple the already retained axial and normal waves through material strain.
    ///
    /// The unchanged material section supplies each strip's EA and mass. Normal
    /// deformation can now drive longitudinal waves and change the same contact
    /// reaction and pressure history. When attached to a string bank, this law
    /// joins its original exact-linear/discrete-gradient contact solve. Existing
    /// nonlinear constituents must use this SAME explicit solve configuration.
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
        self.strain = Some(FiberStrain {
            model,
            slope_scales,
            tension_n,
            axial_rigidity_n,
            limits,
            config,
        });
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
            let bindings = (0..self.resolution.cells)
                .map(|cell| MappedModalPotentialBinding {
                    coordinates: (0..self.cell.modes().len())
                        .map(|mode| ModalCoordinate {
                            component: first_component + cell,
                            mode,
                        })
                        .collect(),
                    replaces_linear: (0..self.resolution.axial_modes).collect(),
                    potential: Box::new(strain.clone()),
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
        if q.len() != self.coordinate_count() || q.iter().any(|x| !x.is_finite()) {
            return Err(invalid());
        }
        let slopes = self.slopes(q);
        let slope_bound: f64 = slopes.iter().map(|x| x.abs()).sum();
        let (mean, variation) = self
            .model
            .strain_bounds(&slopes, &q[..self.model.axial_modes().len()])
            .map_err(|_| invalid())?;
        if !slope_bound.is_finite()
            || slope_bound > self.limits.maximum_slope
            || mean.abs() + variation > self.limits.maximum_strain
            || self.tension_n + self.axial_rigidity_n * (mean - variation) <= 0.0
        {
            return Err(invalid());
        }
        Ok(())
    }
}
