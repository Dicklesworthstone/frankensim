//! Original-source binding for the alternative molecular-kinetic wetting law.
//!
//! The named temperature coordinate belongs to the actual liquid source. Wet
//! thermal cells use their accepted local temperature and native capillary
//! coefficient; dry numerical ghosts retain the reference query. Jump kinetics
//! remain explicit model calibration, separate from static angle receipts.

use super::*;
use fs_qty::Temperature;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

#[derive(Clone)]
pub(super) struct MolecularKineticWetting {
    pub(super) law: MolecularKineticLaw,
    temperature_axis: String,
    reference_temperature_k: f64,
    surface_tension_n_m: f64,
}

impl MolecularKineticWetting {
    pub(super) fn sample(
        &self,
        endpoint_angle_radians: f64,
        signed_speed_m_s: f64,
        model_speed_m_s: f64,
        temperature: Option<(&str, WettingThermalPoint)>,
    ) -> Result<MolecularKineticSample, Error> {
        // A numerical deadband cannot extend the calibrated operating-speed
        // domain. Its modeled speed is zero, while the outer sample keeps the
        // actual geometric speed and does not claim exact immobilization.
        if !signed_speed_m_s.is_finite()
            || signed_speed_m_s.abs() > self.law.controls().maximum_speed_m_s
        {
            return Err(MolecularKineticError::OutsideDomain(
                "observed contact speed exceeds the molecular calibration",
            )
            .into());
        }
        let (temperature_k, surface_tension_n_m) = if let Some((axis, point)) = temperature {
            if axis != self.temperature_axis {
                return Err(Error::Invalid(
                    "molecular wetting requires its admitted liquid-temperature axis",
                ));
            }
            (
                point.temperature_k,
                point.surface_tension_n_m.ok_or(Error::Invalid(
                    "molecular wetting requires the same sourced tension as native pressure",
                ))?,
            )
        } else {
            (self.reference_temperature_k, self.surface_tension_n_m)
        };
        Ok(self.law.apparent_angle(
            endpoint_angle_radians,
            model_speed_m_s,
            temperature_k,
            surface_tension_n_m,
        )?)
    }
}

impl WettingPhysics<'_> {
    pub(in crate::material::free_surface) fn molecular_kinetic_law(
        &self,
    ) -> Option<&MolecularKineticLaw> {
        self.molecular_kinetic.as_ref().map(|binding| &binding.law)
    }

    pub(in crate::material::free_surface) fn with_molecular_kinetic(
        mut self,
        law: MolecularKineticLaw,
        temperature_axis: &str,
        scaling: &MaterialFreeSurfaceScaling<'_>,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        checkpoint(gate)?;
        if self.hydrodynamic.is_some() || self.molecular_kinetic.is_some() {
            return Err(Error::Invalid(
                "contact line already owns a dynamic model; molecular kinetics and Cox cannot be stacked",
            ));
        }
        if self
            .equilibrium
            .as_ref()
            .is_some_and(|binding| !binding.dynamic)
        {
            return Err(Error::Invalid(
                "equilibrium wetting cannot acquire an implicit molecular or hysteresis rate law",
            ));
        }
        let point = scaling.material().material().query();
        let quantity = point.axis_quantities().get(temperature_axis).copied();
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        // Match the material owner's established legacy T convention. Other
        // names need an explicit absolute-temperature role; a temperature
        // difference or an unrelated numeric coordinate is not k_B*T.
        let temperature_ok = quantity == Some(absolute)
            || temperature_axis == "T"
                && (quantity.is_none()
                    || quantity == Some(QuantitySpec::dimensional(Temperature::DIMS)));
        let reference_temperature_k = point
            .axes()
            .get(temperature_axis)
            .copied()
            .filter(|t| temperature_ok && t.is_finite() && *t > 0.0)
            .ok_or(Error::Invalid(
                "molecular wetting requires a named positive absolute-temperature source coordinate",
            ))?;
        let surface_tension_n_m =
            scaling
                .material()
                .surface_tension_n_m()
                .ok_or(Error::Invalid(
                    "molecular wetting requires sourced positive surface tension",
                ))?;
        for property in self.source.properties() {
            if property
                .source_axes()
                .iter()
                .any(|(axis, _)| axis == DYNAMIC_WETTING_SPEED_AXIS)
            {
                return Err(Error::Invalid(
                    "molecular wetting cannot add rate dependence to a measured dynamic-angle curve",
                ));
            }
            if self.equilibrium.is_none() {
                law.apparent_angle(
                    property.value_si(),
                    0.0,
                    reference_temperature_k,
                    surface_tension_n_m,
                )?;
            }
        }
        if let Some(binding) = &self.equilibrium_hysteresis {
            let (receding, advancing) = binding.angles()?;
            law.admit_angle(receding)?;
            law.admit_angle(advancing)?;
        }
        if let Some(equilibrium) = &self.equilibrium {
            let sample = self.equilibrium_sample(equilibrium, 0.0, None, gate)?;
            law.apparent_angle(
                sample.applied_angle_radians,
                0.0,
                reference_temperature_k,
                surface_tension_n_m,
            )?;
        }
        self.molecular_kinetic = Some(MolecularKineticWetting {
            law,
            temperature_axis: temperature_axis.to_owned(),
            reference_temperature_k,
            surface_tension_n_m,
        });
        checkpoint(gate)?;
        Ok(self)
    }
}

impl WettingBinding<'_> {
    pub(in crate::material::free_surface) fn molecular_kinetic_law(
        &self,
    ) -> Option<&MolecularKineticLaw> {
        self.physics.molecular_kinetic_law()
    }

    pub(in crate::material::free_surface) fn with_molecular_kinetic(
        mut self,
        law: MolecularKineticLaw,
        temperature_axis: &str,
        scaling: &MaterialFreeSurfaceScaling<'_>,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        for physics in self.physics_mut() {
            *physics = physics.clone().with_molecular_kinetic(law, temperature_axis, scaling, gate)?;
        }
        let mut samples = BTreeMap::new();
        for (&position, current) in &self.geometry {
            checkpoint(gate)?;
            let speed = self.samples[&position].signed_speed_m_s;
            samples.insert(
                position,
                self.physics_at(position)
                    .sample(current.angle_radians, current.tangent_gradient > 0.0, speed, None, gate)?,
            );
        }
        self.samples = samples;
        Ok(self)
    }
}

impl MaterialFreeSurfaceWetting3<'_, '_> {
    /// Apply calibrated molecular jump kinetics to the STATIC endpoint claims.
    ///
    /// The original source supplies the named absolute temperature and surface
    /// tension. Jump length/frequency are separate model calibration, held only
    /// within the law's declared temperature interval. No Arrhenius dependence,
    /// measured dynamic-angle receipt or bulk-viscosity substitute is inferred.
    /// Existing Cox or measured rate-dependent endpoints refuse. Local line
    /// dissipation is reported per length, never injected as unaccounted heat.
    pub fn with_molecular_kinetic(
        mut self,
        law: MolecularKineticLaw,
        temperature_axis: &str,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        let binding =
            self.binding
                .with_molecular_kinetic(law, temperature_axis, &self.flow.scaling, gate)?;
        binding.install(&mut self.flow.lattice, gate)?;
        checkpoint(gate)?;
        self.binding = binding;
        Ok(self)
    }

    /// Calibration and domain of the optional molecular rate law.
    #[must_use]
    pub fn molecular_kinetic_law(&self) -> Option<&MolecularKineticLaw> {
        self.binding.molecular_kinetic_law()
    }
}
