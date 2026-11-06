//! Elastic normal-stress envelope on the original string's continuous span.
//!
//! For mass-normalized sine coordinates, the vector curvature satisfies
//! `|w''(x)| <= sum_n k_n^2 sqrt(2/(mu L)) ||q_n||`. Pairing the two
//! polarizations before taking their norm keeps this bound independent of the
//! chosen transverse frame. An affine endpoint lift has zero curvature.
//! The original KC owner supplies the shared tension T. Perfect bonding then
//! gives the layer envelope `E_i (T/EA + r_outer_i curvature_bound)`.
//! This is a conservative spatial envelope, not a sampled or interval-certified
//! maximum and not a constitutive fracture or necking law.

use super::{ImpactError, StringPotential, det, invalid};
use fs_material::section::strength::AxialTensileCapacity;
use std::sync::Arc;

/// Elastic tensile demand on a homogeneous or perfectly bonded string section.
/// A crossed conservative envelope can refuse a cancelling modal shape even
/// when its actual maximum is lower. No material state is altered.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StringTensileAssessment {
    /// Axial tension used for the elastic envelope [N]. This is the original
    /// uniform KC tension, or the continuous-span upper tension bound when
    /// the original support retains longitudinal waves.
    pub tension_n: f64,
    /// Continuous-span upper envelope of the transverse curvature norm [1/m].
    pub curvature_bound_per_m: f64,
    /// Original layer index with largest utilization; first wins an exact tie.
    pub governing_layer: usize,
    /// Tensile normal-stress envelope in the governing layer [Pa].
    pub governing_tensile_stress_bound_pa: f64,
    /// Governing elastic normal-stress envelope divided by its sourced limit.
    pub utilization: f64,
}

impl StringTensileAssessment {
    /// Strict comparison with the declared source value; no safety factor.
    #[must_use]
    pub fn exceeds_limit(self) -> bool {
        self.utilization > 1.0
    }
}

#[derive(Clone, Debug)]
pub(in crate::render::plate::impact) struct StringTensileGuard {
    pub(in crate::render::plate::impact) capacity: AxialTensileCapacity,
    curvature_scales: Vec<f64>,
}

impl StringTensileGuard {
    pub(in crate::render::plate::impact) fn new(
        capacity: AxialTensileCapacity,
        axial_rigidity_n: f64,
        curvature_scales: Vec<f64>,
    ) -> Result<Self, ImpactError> {
        let ea = capacity.section().axial_stiffness_n();
        if !axial_rigidity_n.is_finite()
            || (ea - axial_rigidity_n).abs() > 32.0 * f64::EPSILON * ea
            || curvature_scales.is_empty()
            || curvature_scales.iter().any(|s| !s.is_finite() || *s <= 0.0)
        {
            return Err(invalid(
                "string tensile capacity requires its original section EA and finite modal curvature scales",
            ));
        }
        Ok(Self {
            capacity,
            curvature_scales,
        })
    }

    pub(in crate::render::plate::impact) fn admit(
        &self,
        q: &[f64],
        planes: usize,
        tension_n: f64,
    ) -> Result<(), ImpactError> {
        let assessment = self.assess(q, planes, tension_n)?;
        if assessment.exceeds_limit() {
            return Err(ImpactError::Owner(format!(
                "string tensile envelope exceeds '{}' in layer {}: {} Pa, utilization {}",
                self.capacity.limit_property(),
                assessment.governing_layer,
                assessment.governing_tensile_stress_bound_pa,
                assessment.utilization,
            )));
        }
        Ok(())
    }

    pub(in crate::render::plate::impact) fn assess(
        &self,
        q: &[f64],
        planes: usize,
        tension_n: f64,
    ) -> Result<StringTensileAssessment, ImpactError> {
        let modes = self.curvature_scales.len();
        let mut curvature = 0.0;
        for (mode, &scale) in self.curvature_scales.iter().enumerate() {
            let amplitude = if planes == 2 {
                det::hypot(q[mode], q[modes + mode])
            } else {
                q[mode].abs()
            };
            let contribution = scale * amplitude;
            if amplitude > 0.0 && contribution == 0.0 {
                return Err(invalid("string tensile curvature is unrepresentable"));
            }
            curvature += contribution;
        }
        let section = self.capacity.section();
        let axial_strain = tension_n / section.axial_stiffness_n();
        if !curvature.is_finite()
            || !axial_strain.is_finite()
            || (tension_n != 0.0 && axial_strain == 0.0)
        {
            return Err(invalid("string tensile curvature or axial strain overflow"));
        }
        let mut result = StringTensileAssessment {
            tension_n,
            curvature_bound_per_m: curvature,
            governing_layer: 0,
            governing_tensile_stress_bound_pa: 0.0,
            utilization: 0.0,
        };
        for (index, ((layer, integrals), &limit)) in section
            .layers()
            .iter()
            .zip(section.layer_integrals())
            .zip(self.capacity.layer_limits_pa())
            .enumerate()
        {
            let bending_strain = layer.outer_radius_m * curvature;
            if curvature > 0.0 && bending_strain == 0.0 {
                return Err(invalid("string tensile bending strain is unrepresentable"));
            }
            let strain = axial_strain + bending_strain;
            let stress = integrals.young_modulus_pa() * strain.max(0.0);
            let utilization = stress / limit;
            if !stress.is_finite()
                || !utilization.is_finite()
                || (strain > 0.0 && stress == 0.0)
                || (stress > 0.0 && utilization == 0.0)
            {
                return Err(invalid(
                    "string tensile stress or utilization is unrepresentable",
                ));
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

impl StringPotential {
    /// Original per-layer strength claims, geometry and elastic load sharing.
    /// `None` means no tensile criterion was selected.
    #[must_use]
    pub fn tensile_capacity(&self) -> Option<&AxialTensileCapacity> {
        self.tensile.as_ref().map(|guard| &guard.capacity)
    }

    /// Assess the original coordinates without enforcing the selected tensile
    /// limit, so an overlimit trial can be diagnosed. The original finite-state
    /// and moderate-slope domain still apply. This calculation allocates nothing.
    ///
    /// Only elastic axial and bending normal stress are assessed. Viscous
    /// stress, torsional/shear failure, fatigue and fracture evolution are absent.
    pub fn tensile_assessment(
        &self,
        q: &[f64],
    ) -> Result<Option<StringTensileAssessment>, ImpactError> {
        let observation = self.observe_kinematics(q)?;
        self.tensile
            .as_ref()
            .map(|guard| guard.assess(q, self.transverse_planes, observation.tension_n))
            .transpose()
    }

    // Only source-retaining material adapters may attach this guard. It shares
    // the immutable storage and never changes force, energy or modal inertia.
    pub(crate) fn with_tensile_capacity(
        mut self,
        capacity: AxialTensileCapacity,
        length_m: f64,
    ) -> Result<Self, ImpactError> {
        if self.tensile.is_some() || !length_m.is_finite() || length_m <= 0.0 {
            return Err(invalid(
                "string tensile capacity requires a finite span and exclusive limit",
            ));
        }
        let total = self.mode_count();
        let count = (total - usize::from(self.affine_endpoint)) / self.transverse_planes;
        let channel = &self.storage.channels[0];
        let mut curvature_scales = Vec::with_capacity(count);
        for mode in 0..count {
            // Reuse the actual owner's modal normalization; do not rebuild it
            // from a material density that may differ from transverse inertia.
            let wave = (mode + 1) as f64 * core::f64::consts::PI / length_m;
            let scale = wave * det::sqrt(channel.coupling[mode * total + mode]);
            if !scale.is_finite() || scale <= 0.0 {
                return Err(invalid("string tensile curvature scale is unrepresentable"));
            }
            curvature_scales.push(scale);
        }
        self.tensile = Some(Arc::new(StringTensileGuard::new(
            capacity,
            self.stretching.axial_rigidity_n,
            curvature_scales,
        )?));
        Ok(self)
    }
}
