//! Elastic tensile envelope of the original fibers, independent of width cells.
//!
//! Each quadrature strip represents parallel fibers: its force and EA scale
//! together, while radius stays that of ONE physical fiber. The native strain
//! owner supplies the continuous-span axial upper bound; sine curvature obeys
//! `|y''| <= sum k_n^2 sqrt(2/(mu_strip L)) |q_n|`. Perfect bonding then gives
//! the per-layer bound `E_i (epsilon_max + r_outer_i curvature_bound)`.

use super::*;
use fs_material::section::strength::AxialTensileCapacity;
use std::sync::Arc;

/// Conservative elastic normal-stress envelope on one original fiber.
/// This is a spatial modal envelope, not an interval certificate or rupture law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiberTensileAssessment {
    /// Axial force upper bound in one physical fiber [N], including preload.
    pub per_fiber_tension_bound_n: f64,
    /// Continuous-span curvature upper bound [1/m].
    pub curvature_bound_per_m: f64,
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
    capacity: AxialTensileCapacity,
    curvature_scales: Vec<f64>,
}

fn invalid_strength() -> ModalCouplingError {
    ModalCouplingError::Invalid {
        what: "fiber tensile stress or curvature is unrepresentable",
    }
}

impl FiberTensileGuard {
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
        let tension = self.capacity.section().axial_stiffness_n() * axial_strain;
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
            governing_layer: 0,
            governing_tensile_stress_bound_pa: 0.0,
            utilization: 0.0,
        };
        for (index, ((layer, integrals), &limit)) in self
            .capacity
            .section()
            .layers()
            .iter()
            .zip(self.capacity.section().layer_integrals())
            .zip(self.capacity.layer_limits_pa())
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
    /// This explicitly compares ELASTIC normal stress at fixed source state.
    /// Both sourced axial and bending Kelvin--Voigt viscosities must be zero:
    /// total rate stress needs a velocity-aware candidate owner. Rate effects,
    /// yielding, fatigue and fracture are absent.
    /// Hereditary material combinations refuse in either construction order:
    /// equilibrium-only stress must not stand in for missing layer memories.
    /// Source values are retained unchanged, never converted to design allowables.
    pub fn with_tensile_limit(mut self, property: &str) -> Result<Self, FiberRibbonError> {
        if self.relaxation.is_some() {
            return Err(FiberRibbonError::Invalid(
                "fiber tensile admission needs actual hereditary layer stresses",
            ));
        }
        if self.axial.viscous_rigidity_n_s() != 0.0
            || self.bending.viscous_stiffness_n_m2_s() != 0.0
        {
            return Err(FiberRibbonError::Invalid(
                "fiber tensile admission needs actual Kelvin-Voigt strain-rate stresses",
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
        self.tensile.as_ref().map(|guard| &guard.capacity)
    }

    /// Diagnose the supplied strip coordinates without enforcing its strength
    /// limit. Native slope, strain and positive-tension bounds still apply.
    /// Coordinates retain axial-first, normal-second order and strip inertia.
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
}
