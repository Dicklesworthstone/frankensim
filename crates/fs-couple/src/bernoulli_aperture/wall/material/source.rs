//! Shared source-property selection for material walls in wind performances.
//! This binds original material data to the existing wall law; it never fits a
//! spectrum to sound or adds a second damping/evolution path.

use super::{AxialConstraint, LocalMaterialWall, WallMaterialLoss, WallMaterialModel,
    WallRelaxation, YOUNG_VISCOSITY_PROPERTY};
use crate::acoustic_realize::AcousticRealizeError;
use fs_material::state_point::{
    MaterialStatePointError, ResolvedMaterialStatePoint, ScalarAdmissibility as Admit,
    ScalarPropertyRequirement as Requirement, MAX_MATERIAL_PROPERTY_NAME_BYTES,
};
use fs_qty::{Density, Dims, DynViscosity, Pressure};
use fs_vfit::relaxation::MAX_RELAXATION_BRANCHES;
use std::str::SplitAsciiWhitespace;

/// One source selection shared by material assemblies, DG, graph reed and jet
/// ingress. These are original property selectors, not fitted wall coefficients.
#[derive(Clone, Debug, PartialEq)]
pub enum WallConstitutiveSource {
    /// Isotropic Young/Poisson membrane with the declared internal loss source.
    Isotropic(WallMaterialLossSource),
    /// Principal material directions aligned with local hoop and axial motion.
    Orthotropic {
        /// Distinct zero-based hoop and axial directions in `0..3`.
        axes: [usize; 2],
        /// Exact time property for viscosity proportional to elastic stiffness.
        /// `None` selects an elastic wall.
        proportional_viscosity_time_property: Option<String>,
    },
    /// Independent normal-plane Kelvin viscosity with its exact reduced memory.
    OrthotropicNormal {
        /// Distinct zero-based hoop and axial principal directions.
        axes: [usize; 2],
        /// Exact `[D_hh, D_ha, D_aa]` viscosity properties in Pa s.
        properties: [String; 3],
        /// Declared time-domain resolution limit, in `(0, 2]`.
        max_dt_over_tau: f64,
        /// Declared angular-step limit, in `(0, 1]`.
        max_angular_step: f64,
    },
}

impl WallConstitutiveSource {
    pub(crate) fn read(orthotropic: bool, fields: &mut SplitAsciiWhitespace<'_>)
        -> Result<Self, &'static str>
    {
        if !orthotropic { return WallMaterialLossSource::read(fields).map(Self::Isotropic); }
        let mut axes = [0; 2];
        for axis in &mut axes {
            *axis = fields.next().and_then(|s| s.parse::<usize>().ok())
                .filter(|n| (1..=3).contains(n))
                .ok_or("orthotropic wall axes must be distinct integers in 1..=3")? - 1;
        }
        if axes[0] == axes[1] { return Err("orthotropic wall hoop and axial directions must differ"); }
        let proportional_viscosity_time_property = match fields.next() {
            Some("normal-kelvin-voigt") => {
                let max_dt_over_tau = fields.next().and_then(|s| s.parse::<f64>().ok())
                    .filter(|n| n.is_finite() && *n > 0.0 && *n <= 2.0)
                    .ok_or("normal wall viscosity requires finite 0 < dt/tau <= 2")?;
                let max_angular_step = fields.next().and_then(|s| s.parse::<f64>().ok())
                    .filter(|n| n.is_finite() && *n > 0.0 && *n <= 1.0)
                    .ok_or("normal wall viscosity requires finite 0 < angular step <= 1")?;
                let mut name = || fields.next()
                    .filter(|name| !name.is_empty() && name.len() <= MAX_MATERIAL_PROPERTY_NAME_BYTES)
                    .map(str::to_owned).ok_or("normal wall viscosity needs three exact source properties");
                let properties = [name()?, name()?, name()?];
                LocalMaterialWall::orthotropic_normal_viscosity_requirements(
                    axes, properties.each_ref().map(String::as_str))
                    .map_err(|_| "normal wall viscosity properties conflict with the elastic quantities")?;
                return Ok(Self::OrthotropicNormal { axes, properties, max_dt_over_tau, max_angular_step });
            }
            Some("elastic") => None,
            Some("proportional-kelvin-voigt") => Some(fields.next()
                .filter(|name| !name.is_empty() && name.len() <= MAX_MATERIAL_PROPERTY_NAME_BYTES)
                .ok_or("orthotropic wall needs an exact proportional-viscosity time property")?.to_owned()),
            _ => return Err("orthotropic wall loss requires elastic, proportional-kelvin-voigt or normal-kelvin-voigt with explicit source properties"),
        };
        LocalMaterialWall::orthotropic_requirements(axes, proportional_viscosity_time_property.as_deref())
            .map_err(|_| "orthotropic wall has an invalid proportional-viscosity source property")?;
        Ok(Self::Orthotropic { axes, proportional_viscosity_time_property })
    }

    pub(crate) fn requirements(&self) -> Result<Vec<Requirement>, AcousticRealizeError> {
        match self {
            Self::Isotropic(loss) => loss.requirements()
                .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string())),
            Self::Orthotropic { axes, proportional_viscosity_time_property } =>
                LocalMaterialWall::orthotropic_requirements(*axes, proportional_viscosity_time_property.as_deref()),
            Self::OrthotropicNormal { axes, properties, .. } =>
                LocalMaterialWall::orthotropic_normal_viscosity_requirements(*axes,
                    properties.each_ref().map(String::as_str)),
        }
    }

    pub(crate) fn selection(&self, material: &ResolvedMaterialStatePoint, thickness_m: f64,
        axial: AxialConstraint, band_hz: [f64; 2])
        -> Result<(WallMaterialModel, WallMaterialLoss), AcousticRealizeError>
    {
        match self {
            Self::Isotropic(loss) => Ok((WallMaterialModel::Isotropic, loss.resolve(material, band_hz)?)),
            Self::Orthotropic { axes, proportional_viscosity_time_property } => Ok((
                WallMaterialModel::Orthotropic { axes: *axes,
                    proportional_viscosity_time_property: proportional_viscosity_time_property.clone() },
                if proportional_viscosity_time_property.is_some() { WallMaterialLoss::KelvinVoigt }
                else { WallMaterialLoss::Elastic },
            )),
            Self::OrthotropicNormal { axes, properties, max_dt_over_tau, max_angular_step } => {
                let wall = LocalMaterialWall::new_orthotropic_normal_kelvin_voigt(material.clone(),
                    thickness_m, axial, *axes, properties.each_ref().map(String::as_str),
                    band_hz, *max_dt_over_tau, *max_angular_step)?;
                Ok((WallMaterialModel::OrthotropicNormal { axes:*axes, properties:properties.clone(),
                    max_dt_over_tau:*max_dt_over_tau, max_angular_step:*max_angular_step }, wall.loss().clone()))
            }
        }
    }

    pub(crate) fn resolve(&self, material: ResolvedMaterialStatePoint, thickness_m: f64,
        axial: AxialConstraint, band_hz: [f64; 2]) -> Result<LocalMaterialWall, AcousticRealizeError>
    {
        let (model, loss) = self.selection(&material, thickness_m, axial, band_hz)?;
        model.resolve(material, thickness_m, axial, loss, band_hz)
    }

    pub(crate) fn normal_viscosity_properties(&self) -> Option<[String; 3]> {
        match self { Self::OrthotropicNormal { properties, .. } => Some(properties.clone()), _ => None }
    }
}

/// Original material properties supplying isotropic membrane loss and memory.
#[derive(Clone, Debug, PartialEq)]
pub enum WallMaterialLossSource {
    /// Explicitly lossless elastic membrane.
    Elastic,
    /// Resolve `kelvin_voigt_young_viscosity` at the same material state.
    KelvinVoigt,
    /// Resolve each generalized-Maxwell arm from the same original card.
    Maxwell {
        /// Exact relaxing-modulus and relaxation-time property pairs.
        properties: Vec<(String, String)>,
        /// Declared time-domain resolution limit, in `(0, 2]`.
        max_dt_over_tau: f64,
        /// Declared angular-step limit, in `(0, 1]`.
        max_angular_step: f64,
    },
}

impl WallMaterialLossSource {
    // The enclosing reader retains line numbers, checks trailing fields, and
    // owns artifact I/O. Both wind families use exactly this source grammar.
    pub(crate) fn read(fields: &mut SplitAsciiWhitespace<'_>) -> Result<Self, &'static str> {
        match fields.next() {
            Some("elastic") => Ok(Self::Elastic),
            Some("kelvin-voigt") => Ok(Self::KelvinVoigt),
            Some("maxwell") => {
                let count = fields.next().and_then(|s| s.parse::<usize>().ok())
                    .filter(|n| (1..=MAX_RELAXATION_BRANCHES).contains(n))
                    .ok_or("Maxwell wall needs 1..=8 explicit material-property pairs")?;
                let max_dt_over_tau = fields.next().and_then(|s| s.parse::<f64>().ok())
                    .filter(|x| x.is_finite() && *x > 0.0 && *x <= 2.0)
                    .ok_or("Maxwell wall needs finite 0 < max_dt_over_tau <= 2")?;
                let max_angular_step = fields.next().and_then(|s| s.parse::<f64>().ok())
                    .filter(|x| x.is_finite() && *x > 0.0 && *x <= 1.0)
                    .ok_or("Maxwell wall needs finite 0 < max_angular_step <= 1")?;
                let mut properties = Vec::with_capacity(count);
                for _ in 0..count {
                    let modulus = fields.next().filter(|s| s.len() <= MAX_MATERIAL_PROPERTY_NAME_BYTES)
                        .ok_or("Maxwell wall needs a modulus property name of at most 128 bytes")?;
                    let time = fields.next().filter(|s| s.len() <= MAX_MATERIAL_PROPERTY_NAME_BYTES)
                        .ok_or("Maxwell wall needs a relaxation-time property name of at most 128 bytes")?;
                    if properties.iter().any(|(e, tau)| e == modulus && tau == time) {
                        return Err("Maxwell wall must not duplicate a sourced relaxation branch");
                    }
                    properties.push((modulus.to_owned(), time.to_owned()));
                }
                Ok(Self::Maxwell { properties, max_dt_over_tau, max_angular_step })
            }
            _ => Err("wall loss requires elastic, kelvin-voigt, or maxwell COUNT MAX_DT_OVER_TAU MAX_ANGULAR_STEP MODULUS TIME ..."),
        }
    }

    pub(crate) fn requirements(&self) -> Result<Vec<Requirement>, MaterialStatePointError> {
        if let Self::Maxwell { properties, .. } = self
            && !(1..=MAX_RELAXATION_BRANCHES).contains(&properties.len())
        {
            return Err(MaterialStatePointError::InvalidRequirement {
                property: "wall Maxwell spectrum".into(),
                reason: "wall Maxwell source requires 1..=8 explicit property pairs",
            });
        }
        let mut required = vec![
            Requirement::try_new("density", Density::DIMS, Admit::StrictlyPositive)?,
            Requirement::try_new("young_modulus", Pressure::DIMS, Admit::StrictlyPositive)?,
            Requirement::try_new("poisson_ratio", Dims::NONE,
                Admit::OpenInterval { lower: -1.0, upper: 0.5 })?,
        ];
        match self {
            Self::Elastic => {}
            Self::KelvinVoigt => required.push(Requirement::try_new(
                YOUNG_VISCOSITY_PROPERTY, DynViscosity::DIMS, Admit::NonNegative)?),
            Self::Maxwell { properties, .. } => {
                for (modulus, time) in properties {
                    // One source query per distinct property; every declared
                    // branch still enters the constitutive spectrum separately.
                    if !required.iter().any(|p| p.name() == modulus) {
                        required.push(Requirement::try_new(modulus, Pressure::DIMS, Admit::NonNegative)?);
                    }
                    if !required.iter().any(|p| p.name() == time) {
                        required.push(Requirement::try_new(time, Dims([0, 0, 1, 0, 0, 0]),
                            Admit::StrictlyPositive)?);
                    }
                }
            }
        }
        Ok(required)
    }

    pub(crate) fn resolve(&self, material: &ResolvedMaterialStatePoint, band_hz: [f64; 2])
        -> Result<WallMaterialLoss, AcousticRealizeError>
    {
        match self {
            Self::Elastic => Ok(WallMaterialLoss::Elastic),
            Self::KelvinVoigt => Ok(WallMaterialLoss::KelvinVoigt),
            Self::Maxwell { properties, max_dt_over_tau, max_angular_step } => {
                WallRelaxation::from_state_point(material, properties, band_hz,
                    *max_dt_over_tau, *max_angular_step).map(WallMaterialLoss::Maxwell)
            }
        }
    }
}
