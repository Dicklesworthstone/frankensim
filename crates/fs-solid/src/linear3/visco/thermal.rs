//! Imposed eigenstrain on the original, fixed-coefficient material histories.
//!
//! The material advances epsilon_mechanical = B u - epsilon_T. Its exact
//! strain-ramp work is nodal work plus -integral sigma_mean:delta epsilon_T.
//! The latter is prescribed eigenstrain work, not reciprocal caloric heat.
//! Reference geometry, density, tensors, relaxation times and internal variables
//! remain with the same owner. Temperature is retained source metadata here.

use std::sync::Arc;

use super::*;
use crate::linear3::{QuadraticThermalStrainField, TetThermalStrainField, TetThermalStrainState};

// Immutable source samples are shared by frozen-history trials. Only B-bar
// quadrature input needs an additional scalar per element; deviatoric strain
// and all six original source components remain at their own material points.
#[derive(Clone)]
pub(super) struct ThermalSnapshot {
    states: Vec<TetThermalStrainState>,
    points_per_element: Option<usize>,
    mean_normal: Option<Vec<f64>>,
}

fn normal_mean(free: [f64; 6]) -> f64 {
    free[0] / 3.0 + free[1] / 3.0 + free[2] / 3.0
}

impl ThermalSnapshot {
    pub(super) fn prepare(
        system: &ViscoElements,
        field: QuadraticThermalStrainField<'_>,
        previous: Option<&Self>,
        cx: &Cx<'_>,
    ) -> Result<Arc<Self>, ViscoElementError> {
        checkpoint(cx)?;
        if matches!(
            field,
            QuadraticThermalStrainField::ThermomechanicalAtQuadrature { .. }
        ) {
            return Err(TetElasticError::InvalidThermalStrain {
                what: "thermal eigenstrain input cannot replace the original Prony coefficients",
            }
            .into());
        }
        let points = system.elements.first().map_or(0, |e| e.points.len());
        field.validate(system.elements.len(), points)?;
        let points_per_element = match field {
            QuadraticThermalStrainField::Elementwise(_) => None,
            QuadraticThermalStrainField::AtQuadrature {
                points_per_element, ..
            } => Some(points_per_element),
            QuadraticThermalStrainField::ThermomechanicalAtQuadrature { .. } => unreachable!(),
        };
        let count = system
            .elements
            .len()
            .checked_mul(points_per_element.unwrap_or(1))
            .ok_or(TetElasticError::SizeOverflow)?;
        let mut states = Vec::with_capacity(count);
        for (element, material) in system.elements.iter().enumerate() {
            for point in 0..material.points.len() {
                checkpoint(cx)?;
                let state = field.at(element, point);
                if state.elastic_material_state_identity() != material.material_state_identity {
                    return Err(TetElasticError::ThermalElasticStateMismatch { element }.into());
                }
                if previous.is_some_and(|old| {
                    old.at(element, point).reference_temperature_k().to_bits()
                        != state.reference_temperature_k().to_bits()
                }) {
                    return Err(TetElasticError::InvalidThermalStrain {
                        what: "viscoelastic thermal history must retain each original reference temperature",
                    }.into());
                }
                if points_per_element.is_some() || point == 0 {
                    states.push(*state);
                }
            }
        }
        let mut snapshot = Self {
            states,
            points_per_element,
            mean_normal: None,
        };
        // Each prepared system has one element family/formulation. Use the
        // actual positive physical volumes already retained by its B-bar map.
        let bbar = system
            .elements
            .first()
            .and_then(|e| e.points.first())
            .is_some_and(|p| p.compatible_gradients.is_some());
        if bbar && points_per_element.is_some() {
            let mut means = Vec::with_capacity(system.elements.len());
            for (element, material) in system.elements.iter().enumerate() {
                checkpoint(cx)?;
                let volume = material.points.iter().map(|p| p.volume_m3).sum::<f64>();
                if !volume.is_finite() || volume <= 0.0 {
                    return Err(ViscoElementError::NonFiniteResult);
                }
                let anchor = normal_mean(snapshot.at(element, 0).free_strain_mandel());
                let mut mean = anchor;
                for (point, geometry) in material.points.iter().enumerate() {
                    checkpoint(cx)?;
                    let local = normal_mean(snapshot.at(element, point).free_strain_mandel());
                    mean = (geometry.volume_m3 / volume).mul_add(local - anchor, mean);
                }
                if !mean.is_finite() {
                    return Err(ViscoElementError::NonFiniteResult);
                }
                means.push(mean);
            }
            snapshot.mean_normal = Some(means);
            for (element, material) in system.elements.iter().enumerate() {
                for point in 0..material.points.len() {
                    checkpoint(cx)?;
                    if snapshot
                        .free_strain(element, point)
                        .iter()
                        .any(|v| !v.is_finite())
                    {
                        return Err(ViscoElementError::NonFiniteResult);
                    }
                }
            }
        }
        checkpoint(cx)?;
        Ok(Arc::new(snapshot))
    }

    // A coefficient refresh may replace only the elastic identity. Every
    // physical source sample and its assumed strain must survive unchanged.
    pub(super) fn rebound(
        &self,
        next: &ViscoElements,
        field: QuadraticThermalStrainField<'_>,
        cx: &Cx<'_>,
    ) -> Result<Arc<Self>, ViscoElementError> {
        checkpoint(cx)?;
        let points = next.elements.first().map_or(0, |e| e.points.len());
        self.field().validate(next.elements.len(), points)?;
        let rebound = Self::prepare(next, field, Some(self), cx)?;
        for (element, material) in next.elements.iter().enumerate() {
            for point in 0..material.points.len() {
                checkpoint(cx)?;
                let old = self.at(element, point);
                let new = rebound.at(element, point);
                if old.free_strain_mandel().map(f64::to_bits)
                    != new.free_strain_mandel().map(f64::to_bits)
                    || old.temperature_k().to_bits() != new.temperature_k().to_bits()
                    || old.reference_temperature_k().to_bits()
                        != new.reference_temperature_k().to_bits()
                    || old.thermal_law_identity() != new.thermal_law_identity()
                    || old.maximum_absolute_mandel_strain().to_bits()
                        != new.maximum_absolute_mandel_strain().to_bits()
                    || self.free_strain(element, point).map(f64::to_bits)
                        != rebound.free_strain(element, point).map(f64::to_bits)
                {
                    return Err(TetElasticError::InvalidThermalStrain {
                        what: "coefficient refresh must retain each original thermal sample and projected free strain",
                    }
                    .into());
                }
            }
        }
        checkpoint(cx)?;
        Ok(rebound)
    }

    pub(super) fn field(&self) -> QuadraticThermalStrainField<'_> {
        match self.points_per_element {
            Some(points_per_element) => QuadraticThermalStrainField::AtQuadrature {
                states: &self.states,
                points_per_element,
            },
            None => QuadraticThermalStrainField::Elementwise(TetThermalStrainField::PerElement(
                &self.states,
            )),
        }
    }

    pub(super) fn at(&self, element: usize, point: usize) -> &TetThermalStrainState {
        &self.states[self
            .points_per_element
            .map_or(element, |n| element * n + point)]
    }

    pub(super) fn free_strain(&self, element: usize, point: usize) -> [f64; 6] {
        let mut free = self.at(element, point).free_strain_mandel();
        if let Some(means) = &self.mean_normal {
            let correction = means[element] - normal_mean(free);
            if correction != 0.0 {
                for value in &mut free[..3] {
                    *value += correction;
                }
            }
        }
        free
    }
}

impl ViscoElementState<'_> {
    /// Accepted raw thermal states, retaining their complete identities.
    /// `None` means zero eigenstrain before the first thermal input. Thereafter
    /// reference temperature is fixed at each original material point. Samples
    /// use the layout returned by `thermal_field`; B-bar never replaces them by
    /// projected values or an element average.
    pub fn thermal_strains(&self) -> Option<&[TetThermalStrainState]> {
        self.thermal_strains.as_ref().map(|s| s.states.as_slice())
    }

    /// Accepted source samples with their original element/quadrature layout.
    pub fn thermal_field(&self) -> Option<QuadraticThermalStrainField<'_>> {
        self.thermal_strains.as_ref().map(|s| s.field())
    }

    /// Cumulative signed imposed-eigenstrain work since the original relaxed
    /// state [J]. It is already included in `energy().work_j`, not added to heat.
    pub const fn eigenstrain_work_j(&self) -> f64 {
        self.eigenstrain_work_j
    }

    /// Advance displacement and source-integrated free strain over one interval.
    ///
    /// Both strain fields interpolate affinely between their accepted endpoints;
    /// nonlinear expansion versus temperature requires interval refinement.
    /// Each input must name the ORIGINAL frozen coefficient material identity,
    /// as returned by `TetViscoMaterial::equilibrium_material`. A newly resolved
    /// target-temperature stiffness identity is refused. The integrated-path
    /// identity and strain ceiling are retained, and reference temperature may
    /// not change after the first accepted thermal step. Source selection,
    /// complete temperature-path support and mesh-frame orientation remain the
    /// responsibility of the original expansion resolver.
    ///
    /// Uniform, elementwise and original-quadrature inputs share this path.
    /// Quadratic source-resolved samples convert directly into it. B-bar applies
    /// the same physical-volume-mean volumetric projection to total and free
    /// strain, retaining pointwise deviatoric expansion. Both projected and raw
    /// compatible total/mechanical strains independently meet source limits.
    /// Local replacement elastic tensors, evolving Prony coefficients and
    /// reciprocal caloric heat/capacity are outside this eigenstrain input.
    /// Failure leaves displacement, thermal state, all memories and work intact.
    pub fn step_internal_force_with_thermal_strain<'t>(
        &mut self,
        u_new: &[f64],
        dt_s: f64,
        thermal: impl Into<QuadraticThermalStrainField<'t>>,
        cx: &Cx<'_>,
    ) -> Result<ViscoForceStep, ViscoElementError> {
        let thermal = self.prepare_thermal_field(thermal.into(), cx)?;
        self.step_internal_force_impl(u_new, dt_s, Some(&thermal), cx)
    }

    pub(super) fn prepare_thermal_field(
        &self,
        thermal: QuadraticThermalStrainField<'_>,
        cx: &Cx<'_>,
    ) -> Result<Arc<ThermalSnapshot>, ViscoElementError> {
        ThermalSnapshot::prepare(self.system, thermal, self.thermal_strains.as_deref(), cx)
    }

    pub(super) fn thermal_increment(
        &self,
        element: usize,
        point: usize,
        thermal: Option<&ThermalSnapshot>,
    ) -> [f64; 6] {
        // No thermal input holds the accepted field exactly, without subtracting
        // or reconstructing its source-integrated endpoint.
        let Some(field) = thermal else {
            return [0.0; 6];
        };
        let next = field.free_strain(element, point);
        let old = self
            .thermal_strains
            .as_ref()
            .map_or([0.0; 6], |states| states.free_strain(element, point));
        core::array::from_fn(|i| next[i] - old[i])
    }

    // Exact linear offset in the original mean-force operator. This avoids
    // evaluating a fictitious clamped thermal state, whose strain may exceed
    // the source domain even when the actual free expansion is admissible.
    pub(super) fn add_thermal_mean_force(
        &self,
        thermal: Option<&ThermalSnapshot>,
        tangents: &[Mandel6],
        reduced: &[Option<usize>],
        force: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ViscoElementError> {
        let Some(thermal) = thermal else {
            return Ok(());
        };
        let constant = thermal.points_per_element.is_none()
            && self
                .thermal_strains
                .as_ref()
                .is_none_or(|s| s.points_per_element.is_none());
        for (index, (element, tangent)) in self.system.elements.iter().zip(tangents).enumerate() {
            checkpoint(cx)?;
            let stress_at = |point| {
                let delta = self.thermal_increment(index, point, Some(thermal));
                core::array::from_fn::<f64, 6, _>(|i| {
                    -(0..6).map(|j| tangent[i][j] * delta[j]).sum::<f64>()
                })
            };
            let constant_stress = constant.then(|| stress_at(0));
            for (q, point) in element.points.iter().enumerate() {
                checkpoint(cx)?;
                let stress = constant_stress.unwrap_or_else(|| stress_at(q));
                for (b, &dof) in point.b.iter().zip(&element.dofs) {
                    if let Some(row) = reduced[dof] {
                        force[row] +=
                            point.volume_m3 * b.iter().zip(stress).map(|(b, s)| b * s).sum::<f64>();
                    }
                }
            }
        }
        Ok(())
    }
}

impl IntegrationPoint {
    pub(super) fn thermal_mechanical_strain(
        &self,
        element: &Element,
        displacement: &[f64],
        total: [f64; 6],
        thermal: TetThermalStrainState,
        assumed_free: [f64; 6],
    ) -> Result<[f64; 6], ViscoElementError> {
        let free = thermal.free_strain_mandel();
        let mechanical = core::array::from_fn(|i| total[i] - assumed_free[i]);
        let admit = |total: &[f64; 6], mechanical: &[f64; 6]| -> Result<(), ViscoElementError> {
            if total
                .iter()
                .chain(mechanical)
                .any(|v| !v.is_finite() || v.abs() > thermal.maximum_absolute_mandel_strain())
            {
                return Err(TetElasticError::InvalidThermalStrain {
                    what: "total or mechanical strain exceeds the original thermal small-strain ceiling",
                }.into());
            }
            element.law.validate_mandel_strain(total)?;
            element.law.validate_mandel_strain(mechanical)?;
            Ok(())
        };
        admit(&total, &mechanical)?;
        if self.compatible_gradients.is_some() {
            let raw = self.compatible_strain(&element.dofs, displacement, total);
            let raw_mechanical = core::array::from_fn(|i| raw[i] - free[i]);
            admit(&raw, &raw_mechanical)?;
        }
        Ok(mechanical)
    }
}
