//! Tensile envelope of the original fibers, independent of width cells.
//!
//! Each quadrature strip represents parallel fibers: its force and EA scale
//! together, while radius stays that of ONE physical fiber. The native strain
//! owner supplies the continuous-span axial upper bound; sine curvature obeys
//! `|y''| <= sum k_n^2 sqrt(2/(mu_strip L)) |q_n|`. Perfect bonding then gives
//! the elastic per-layer bound `E_i (epsilon_max + r_outer_i curvature_bound)`.
//! With sourced viscosity, signed elastic and rate stress are combined in each
//! retained harmonic BEFORE taking the continuous-span magnitude envelope.

use super::*;
use fs_material::section::strength::{AxialTensileCapacity, KelvinVoigtTensileCapacity};
use fs_nlmodal::longitudinal::{MAX_LONGITUDINAL_MODES, MAX_TRANSVERSE_MODES};
use std::sync::Arc;

/// Conservative normal-stress envelope on one original fiber.
/// This is a spatial modal envelope, not an interval certificate or rupture law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiberTensileAssessment {
    /// Axial force upper bound in one physical fiber [N], including preload
    /// and the supplied state's Kelvin--Voigt stress when velocities are given.
    pub per_fiber_tension_bound_n: f64,
    /// Continuous-span curvature upper bound [1/m].
    pub curvature_bound_per_m: f64,
    /// Continuous-span curvature-rate upper bound [1/(m s)]; zero for an
    /// explicitly displacement-only elastic assessment.
    pub curvature_rate_bound_per_m_s: f64,
    /// Original layer with greatest utilization; first wins an exact tie.
    pub governing_layer: usize,
    /// Governing tensile normal-stress bound [Pa].
    pub governing_tensile_stress_bound_pa: f64,
    /// Greatest layer stress bound divided by that layer's sourced limit.
    pub utilization: f64,
}

impl FiberTensileAssessment {
    /// Strict comparison to the selected source value, without a safety factor.
    #[must_use]
    pub fn exceeds_limit(self) -> bool {
        self.utilization > 1.0
    }
}

#[derive(Clone, Debug)]
pub(super) struct FiberTensileGuard {
    capacity: KelvinVoigtTensileCapacity,
    curvature_scales: Vec<f64>,
}

fn invalid_strength() -> ModalCouplingError {
    ModalCouplingError::Invalid {
        what: "fiber tensile stress or curvature is unrepresentable",
    }
}

impl FiberTensileGuard {
    pub(super) fn has_viscosity(&self) -> bool {
        self.capacity.axial_viscous_rigidity_n_s() != 0.0
            || self.capacity.bending_viscous_stiffness_n_m2_s() != 0.0
    }

    pub(super) fn assess(
        &self,
        normal_q: &[f64],
        axial_strain: f64,
    ) -> Result<FiberTensileAssessment, ModalCouplingError> {
        let mut curvature = 0.0;
        for (&scale, &q) in self.curvature_scales.iter().zip(normal_q) {
            let contribution = scale * q.abs();
            if q != 0.0 && contribution == 0.0 {
                return Err(invalid_strength());
            }
            curvature += contribution;
        }
        let capacity = self.capacity.elastic_capacity();
        let tension = capacity.section().axial_stiffness_n() * axial_strain;
        if !curvature.is_finite()
            || !tension.is_finite()
            || !axial_strain.is_finite()
            || (axial_strain != 0.0 && tension == 0.0)
        {
            return Err(invalid_strength());
        }
        let mut result = FiberTensileAssessment {
            per_fiber_tension_bound_n: tension,
            curvature_bound_per_m: curvature,
            curvature_rate_bound_per_m_s: 0.0,
            governing_layer: 0,
            governing_tensile_stress_bound_pa: 0.0,
            utilization: 0.0,
        };
        for (index, ((layer, integrals), &limit)) in capacity
            .section()
            .layers()
            .iter()
            .zip(capacity.section().layer_integrals())
            .zip(capacity.layer_limits_pa())
            .enumerate()
        {
            let bending = layer.outer_radius_m * curvature;
            let strain = (axial_strain + bending).max(0.0);
            let stress = integrals.young_modulus_pa() * strain;
            let utilization = stress / limit;
            if !stress.is_finite()
                || !utilization.is_finite()
                || (curvature > 0.0 && bending == 0.0)
                || (strain > 0.0 && stress == 0.0)
                || (stress > 0.0 && utilization == 0.0)
            {
                return Err(invalid_strength());
            }
            if utilization > result.utilization {
                result.governing_layer = index;
                result.governing_tensile_stress_bound_pa = stress;
                result.utilization = utilization;
            }
        }
        Ok(result)
    }

    fn assess_state(
        &self,
        normal_q: &[f64],
        normal_velocity: &[f64],
        mean_strain: f64,
        mean_strain_rate: f64,
        strain_harmonics: &[f64],
        rate_harmonics: &[f64],
    ) -> Result<FiberTensileAssessment, ModalCouplingError> {
        let capacity = self.capacity.elastic_capacity();
        let ea = capacity.section().axial_stiffness_n();
        let ca = self.capacity.axial_viscous_rigidity_n_s();
        let mut tension = response(ea, mean_strain, ca, mean_strain_rate)?;
        for (&strain, &rate) in strain_harmonics.iter().zip(rate_harmonics) {
            tension = sum(tension, response(ea, strain, ca, rate)?.abs())?;
        }
        let mut curvature = 0.0;
        let mut curvature_rate = 0.0;
        for ((&scale, &q), &v) in self
            .curvature_scales
            .iter()
            .zip(normal_q)
            .zip(normal_velocity)
        {
            curvature = sum(curvature, product(scale, q)?.abs())?;
            curvature_rate = sum(curvature_rate, product(scale, v)?.abs())?;
        }
        let mut result = FiberTensileAssessment {
            per_fiber_tension_bound_n: tension,
            curvature_bound_per_m: curvature,
            curvature_rate_bound_per_m_s: curvature_rate,
            governing_layer: 0,
            governing_tensile_stress_bound_pa: 0.0,
            utilization: 0.0,
        };
        for (index, ((layer, geometry), &limit)) in self
            .capacity
            .layers()
            .iter()
            .zip(capacity.section().layers())
            .zip(capacity.layer_limits_pa())
            .enumerate()
        {
            let young = layer.young_modulus_pa();
            let axial_eta = layer.axial_viscosity_pa_s();
            let bending_eta = layer.bending_viscosity_pa_s();
            let mut stress = response(young, mean_strain, axial_eta, mean_strain_rate)?;
            for (&strain, &rate) in strain_harmonics.iter().zip(rate_harmonics) {
                stress = sum(stress, response(young, strain, axial_eta, rate)?.abs())?;
            }
            let mut bending = 0.0;
            for ((&scale, &q), &v) in self
                .curvature_scales
                .iter()
                .zip(normal_q)
                .zip(normal_velocity)
            {
                let mode_stress = response(young, q, bending_eta, v)?;
                bending = sum(bending, product(scale, mode_stress)?.abs())?;
            }
            let stress = sum(stress, product(geometry.outer_radius_m, bending)?)?.max(0.0);
            let utilization = stress / limit;
            if !utilization.is_finite() || (stress > 0.0 && utilization == 0.0) {
                return Err(invalid_strength());
            }
            if utilization > result.utilization {
                result.governing_layer = index;
                result.governing_tensile_stress_bound_pa = stress;
                result.utilization = utilization;
            }
        }
        Ok(result)
    }
}

fn product(a: f64, b: f64) -> Result<f64, ModalCouplingError> {
    let value = a * b;
    if !value.is_finite() || (a != 0.0 && b != 0.0 && value == 0.0) {
        Err(invalid_strength())
    } else {
        Ok(value)
    }
}

fn sum(a: f64, b: f64) -> Result<f64, ModalCouplingError> {
    let value = a + b;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid_strength())
    }
}

fn response(
    modulus: f64,
    strain: f64,
    viscosity: f64,
    rate: f64,
) -> Result<f64, ModalCouplingError> {
    sum(product(modulus, strain)?, product(viscosity, rate)?)
}

impl TensionedFiberRibbon {
    /// Select a sourced tensile limit on every original fiber layer.
    ///
    /// Call after `with_axial_transverse_strain`. Preload and every candidate
    /// accepted by the existing coupled contact solver use the SAME retained
    /// axial waves, geometric extension and per-fiber bending curvature. A
    /// crossed envelope refuses transactionally; forces, energy, and coordinates
    /// are never clipped. A cancelling shape can cross this conservative bound.
    ///
    /// This compares total elastic and Kelvin--Voigt normal stress at fixed
    /// source state. Every original layer keeps its own axial and bending
    /// viscosities. Rate stress uses actual native endpoint velocities, and
    /// signed elastic/rate modes are combined before spatial bounding. The
    /// source strength itself has no inferred strain-rate dependence. Yielding,
    /// fatigue and fracture are absent.
    /// Hereditary material combinations refuse in either construction order:
    /// equilibrium-only stress must not stand in for missing layer memories.
    /// Source values are retained unchanged, never converted to design allowables.
    pub fn with_tensile_limit(mut self, property: &str) -> Result<Self, FiberRibbonError> {
        if self.relaxation.is_some() {
            return Err(FiberRibbonError::Invalid(
                "fiber tensile admission needs actual hereditary layer stresses",
            ));
        }
        let strain = self.strain.as_mut().ok_or(FiberRibbonError::Invalid(
            "fiber tensile admission requires its original nonlinear strain owner",
        ))?;
        if strain.tensile.is_some() {
            return Err(FiberRibbonError::Invalid(
                "fiber tensile limit is already selected",
            ));
        }
        let capacity = self
            .axial
            .elastic()
            .axial_tensile_capacity(property)
            .and_then(|capacity| capacity.with_kelvin_voigt(&self.axial, &self.bending))
            .map_err(FiberRibbonError::Material)?;
        let curvature_scales: Vec<_> = strain
            .slope_scales
            .iter()
            .enumerate()
            .map(|(i, scale)| {
                scale * (i + 1) as f64 * core::f64::consts::PI / self.geometry.length_m
            })
            .collect();
        if curvature_scales
            .iter()
            .any(|scale| !scale.is_finite() || *scale <= 0.0)
        {
            return Err(FiberRibbonError::Invalid(
                "fiber tensile curvature map is unrepresentable",
            ));
        }
        let guard = FiberTensileGuard {
            capacity,
            curvature_scales,
        };
        let preload_strain = strain.tension_n / strain.axial_rigidity_n;
        if !preload_strain.is_finite() || preload_strain <= 0.0 {
            return Err(FiberRibbonError::Invalid(
                "fiber tensile preload strain is unrepresentable",
            ));
        }
        let preload = guard
            .assess(&vec![0.0; self.resolution.transverse_modes], preload_strain)
            .map_err(|_| FiberRibbonError::Invalid("fiber tensile preload is unrepresentable"))?;
        if preload.exceeds_limit() {
            return Err(FiberRibbonError::Invalid(
                "fiber tensile preload exceeds its material limit",
            ));
        }
        strain.tensile = Some(Arc::new(guard));
        Ok(self)
    }
}

impl FiberStrain {
    /// Complete source section and named per-layer tensile limits, if selected.
    #[must_use]
    pub fn tensile_capacity(&self) -> Option<&AxialTensileCapacity> {
        self.tensile
            .as_ref()
            .map(|guard| guard.capacity.elastic_capacity())
    }

    /// Diagnose elastic stress alone without enforcing its strength limit.
    /// Native slope, strain and positive-tension bounds still apply. Coordinates
    /// retain axial-first, normal-second order and strip inertia. For viscous
    /// material, use `tensile_state_assessment` to include total rate stress.
    pub fn tensile_assessment(
        &self,
        q: &[f64],
    ) -> Result<Option<FiberTensileAssessment>, ModalCouplingError> {
        let (mean, variation) = self.validate_kinematics(q)?;
        self.tensile
            .as_ref()
            .map(|guard| {
                guard.assess(
                    &q[self.model.axial_modes().len()..],
                    self.tension_n / self.axial_rigidity_n + mean + variation,
                )
            })
            .transpose()
    }

    /// Diagnose total sourced elastic/Kelvin--Voigt stress at an actual state.
    /// Both slices use the original axial-first, normal-second mass basis.
    /// The complete retained strain rate includes `y_x*y_xt` and its mean.
    /// Elastic and viscous coefficients combine with their signs BEFORE taking
    /// span bounds, so unloading may cancel tensile stress. No force, storage,
    /// rate-dependent strength or history is inferred. The viscous assessment
    /// reuses the original longitudinal kinematics in bounded stack workspaces
    /// without heap allocation.
    pub fn tensile_state_assessment(
        &self,
        q: &[f64],
        velocity: &[f64],
    ) -> Result<Option<FiberTensileAssessment>, ModalCouplingError> {
        if velocity.len() != self.coordinate_count() || velocity.iter().any(|v| !v.is_finite()) {
            return Err(ModalCouplingError::Invalid {
                what: "fiber tensile stress needs complete original modal velocities",
            });
        }
        let Some(guard) = &self.tensile else {
            self.validate_kinematics(q)?;
            return Ok(None);
        };
        if !guard.has_viscosity() {
            let mut assessment = self.tensile_assessment(q)?;
            if let Some(assessment) = &mut assessment {
                for (&scale, &v) in guard
                    .curvature_scales
                    .iter()
                    .zip(&velocity[self.model.axial_modes().len()..])
                {
                    assessment.curvature_rate_bound_per_m_s = sum(
                        assessment.curvature_rate_bound_per_m_s,
                        product(scale, v)?.abs(),
                    )?;
                }
            }
            return Ok(assessment);
        }
        if q.len() != self.coordinate_count() || q.iter().any(|q| !q.is_finite()) {
            return Err(invalid());
        }
        let axial_count = self.model.axial_modes().len();
        let normal_count = self.slope_scales.len();
        let mut slopes = [0.0; MAX_TRANSVERSE_MODES + 1];
        let mut slope_velocity = [0.0; MAX_TRANSVERSE_MODES + 1];
        let mut slope_bound = 0.0;
        for (index, &scale) in self.slope_scales.iter().enumerate() {
            slopes[index + 1] = product(scale, q[axial_count + index])?;
            slope_velocity[index + 1] = product(scale, velocity[axial_count + index])?;
            slope_bound = sum(slope_bound, slopes[index + 1].abs())?;
        }
        let slopes = &slopes[..=normal_count];
        let slope_velocity = &slope_velocity[..=normal_count];
        let mut strain = [0.0; MAX_LONGITUDINAL_MODES];
        let strain = &mut strain[..axial_count];
        let mean = self
            .model
            .geometric_strain(slopes, strain)
            .map_err(|_| invalid_strength())?;
        let mut variation = 0.0;
        for ((harmonic, &scale), &q) in strain
            .iter_mut()
            .zip(self.model.axial_slope_scales())
            .zip(&q[..axial_count])
        {
            *harmonic = sum(*harmonic, product(scale, q)?)?;
            variation = sum(variation, harmonic.abs())?;
        }
        self.validate_kinematic_bounds(slope_bound, mean, variation)?;
        let mut rates = [0.0; MAX_LONGITUDINAL_MODES];
        let rates = &mut rates[..axial_count];
        let mean_rate = self
            .model
            .strain_rate(slopes, slope_velocity, &velocity[..axial_count], rates)
            .map_err(|_| invalid_strength())?;
        guard
            .assess_state(
                &q[axial_count..],
                &velocity[axial_count..],
                sum(self.tension_n / self.axial_rigidity_n, mean)?,
                mean_rate,
                strain,
                rates,
            )
            .map(Some)
    }
}
