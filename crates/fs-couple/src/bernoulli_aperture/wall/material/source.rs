//! Shared source-property selection for material walls in wind performances.
//! This binds original material data to the existing wall law; it never fits a
//! spectrum to sound or adds a second damping/evolution path.

use super::{WallMaterialLoss, WallRelaxation, YOUNG_VISCOSITY_PROPERTY};
use crate::acoustic_realize::AcousticRealizeError;
use fs_material::state_point::{
    MaterialStatePointError, ResolvedMaterialStatePoint, ScalarAdmissibility as Admit,
    ScalarPropertyRequirement as Requirement, MAX_MATERIAL_PROPERTY_NAME_BYTES,
};
use fs_qty::{Density, Dims, DynViscosity, Pressure};
use fs_vfit::relaxation::MAX_RELAXATION_BRANCHES;
use std::str::SplitAsciiWhitespace;

pub(crate) enum WallMaterialLossSource {
    Elastic,
    KelvinVoigt,
    Maxwell {
        properties: Vec<(String, String)>,
        max_dt_over_tau: f64,
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
