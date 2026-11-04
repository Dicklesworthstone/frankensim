//! Imposed eigenstrain on the original, fixed-coefficient material histories.
//!
//! The material advances epsilon_mechanical = B u - epsilon_T. Its exact
//! strain-ramp work is nodal work plus -integral sigma_mean:delta epsilon_T.
//! The latter is prescribed eigenstrain work, not reciprocal caloric heat.
//! Reference geometry, density, tensors, relaxation times and internal variables
//! remain with the same owner. Temperature is retained source metadata here.

use super::*;
use crate::linear3::{TetThermalStrainField, TetThermalStrainState};

impl ViscoElementState<'_> {
    /// Accepted elementwise thermal states, retaining their complete identities.
    /// `None` means zero eigenstrain before the first thermal input. Thereafter
    /// the original reference temperature is fixed independently for each element.
    pub fn thermal_strains(&self) -> Option<&[TetThermalStrainState]> {
        self.thermal_strains.as_deref()
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
    /// Input is uniform or elementwise, including for B-bar. Pointwise thermal
    /// pressure treatment, evolving Prony coefficients and reciprocal heat or
    /// heat capacity are outside this imposed-expansion owner. Both assumed and
    /// raw compatible TOTAL and mechanical strains must meet source limits.
    /// Failure leaves displacement, thermal state, all memories and work intact.
    pub fn step_internal_force_with_thermal_strain(
        &mut self,
        u_new: &[f64],
        dt_s: f64,
        thermal: TetThermalStrainField<'_>,
        cx: &Cx<'_>,
    ) -> Result<ViscoForceStep, ViscoElementError> {
        self.step_internal_force_impl(u_new, dt_s, Some(thermal), cx)
    }

    pub(super) fn validate_thermal_field(
        &self,
        thermal: Option<TetThermalStrainField<'_>>,
        cx: &Cx<'_>,
    ) -> Result<(), ViscoElementError> {
        let Some(field) = thermal else { return Ok(()) };
        field.validate(self.system.elements.len())?;
        for (element, material) in self.system.elements.iter().enumerate() {
            checkpoint(cx)?;
            let state = field.at(element);
            if state.elastic_material_state_identity() != material.material_state_identity {
                return Err(TetElasticError::ThermalElasticStateMismatch { element }.into());
            }
            if let Some(old) = &self.thermal_strains {
                if old[element].reference_temperature_k().to_bits()
                    != state.reference_temperature_k().to_bits()
                {
                    return Err(TetElasticError::InvalidThermalStrain {
                        what: "viscoelastic thermal history must retain its original reference temperature",
                    }.into());
                }
            }
        }
        Ok(())
    }

    pub(super) fn target_thermal_strain(
        &self,
        element: usize,
        thermal: Option<TetThermalStrainField<'_>>,
    ) -> Option<TetThermalStrainState> {
        thermal.map(|field| *field.at(element))
            .or_else(|| self.thermal_strains.as_ref().map(|states| states[element]))
    }

    pub(super) fn thermal_increment(
        &self,
        element: usize,
        thermal: Option<TetThermalStrainField<'_>>,
    ) -> [f64; 6] {
        // No thermal input holds the accepted field exactly, without subtracting
        // or reconstructing its source-integrated endpoint.
        let Some(field) = thermal else { return [0.0; 6] };
        let next = field.at(element).free_strain_mandel();
        let old = self.thermal_strains.as_ref()
            .map_or([0.0; 6], |states| states[element].free_strain_mandel());
        core::array::from_fn(|i| next[i] - old[i])
    }

    pub(super) fn copy_thermal_field(
        &self,
        thermal: Option<TetThermalStrainField<'_>>,
        cx: &Cx<'_>,
    ) -> Result<Option<Vec<TetThermalStrainState>>, ViscoElementError> {
        let Some(field) = thermal else { return Ok(None) };
        let mut states = Vec::with_capacity(self.system.elements.len());
        for element in 0..self.system.elements.len() {
            checkpoint(cx)?;
            states.push(*field.at(element));
        }
        Ok(Some(states))
    }

    // Exact linear offset in the original mean-force operator. This avoids
    // evaluating a fictitious clamped thermal state, whose strain may exceed
    // the source domain even when the actual free expansion is admissible.
    pub(super) fn add_thermal_mean_force(
        &self,
        thermal: Option<TetThermalStrainField<'_>>,
        tangents: &[Mandel6],
        reduced: &[Option<usize>],
        force: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ViscoElementError> {
        if thermal.is_none() {
            return Ok(());
        }
        for (index, (element, tangent)) in self.system.elements.iter().zip(tangents).enumerate() {
            checkpoint(cx)?;
            let delta = self.thermal_increment(index, thermal);
            let stress: [f64; 6] = core::array::from_fn(|i| {
                -(0..6).map(|j| tangent[i][j] * delta[j]).sum::<f64>()
            });
            for point in &element.points {
                checkpoint(cx)?;
                for (b, &dof) in point.b.iter().zip(&element.dofs) {
                    if let Some(row) = reduced[dof] {
                        force[row] += point.volume_m3
                            * b.iter().zip(stress).map(|(b, s)| b * s).sum::<f64>();
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
    ) -> Result<[f64; 6], ViscoElementError> {
        let free = thermal.free_strain_mandel();
        let mechanical = core::array::from_fn(|i| total[i] - free[i]);
        let admit = |total: &[f64; 6], mechanical: &[f64; 6]| -> Result<(), ViscoElementError> {
            if total.iter().chain(mechanical).any(|v| {
                !v.is_finite() || v.abs() > thermal.maximum_absolute_mandel_strain()
            }) {
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
