//! Source-bound mixed mechanics with the same held thermal reference as the
//! displacement material. No native constitutive law or pressure inverse is
//! replaced, and this material adapter does not advance temperature.
//!
//! With F=Fe Fth and Jth=det(Fth), the original-volume saddle potential is
//! Psi0=Jth[Wrem(Fe)+p(Je-1)-U*(p)]. Thus the pressure equation is
//! J-Jth J(p), its tangent is -Jth/U''(J(p)), and both cross blocks are cof(F).
//! Exact elastic incompressibility imposes J=Jth. It does not prevent the
//! source-declared thermal expansion. Pressure remains tensile-positive in Pa.
//! The public original-reference pressure conjugate is
//! U0*(p)=Jth Ue*(p)+p(Jth-1), with derivative Jth Je(p)-1. This linear
//! offset changes neither the full saddle potential nor the pressure compliance.
//!
//! Every original energy and volume diagnostic is pulled back along with the
//! Piola stress and full consistent tangent. Native P2/P1 mixed equilibrium,
//! motion and vibration consumers may borrow this exact material directly;
//! they retain their own geometry, quadrature, mass, stability and work gates.

use super::FiniteSolidMaterial;
use super::thermal::ThermalReference;
use fs_material::MaterialError;
use fs_material::hyper::Hyperelastic;
use fs_material::hyper::mixed::{
    MixedFiniteStrainPotential, MixedPointResponse, MixedPointState, PressurePotential,
    RemainderResponse,
    VolumetricTreatment,
};

fn finite(value: f64) -> Result<f64, MaterialError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(MaterialError::State {
            what: "nonfinite thermal-reference mixed response".into(),
        })
    }
}

fn conjugate(volume: f64, pressure: f64, elastic: f64) -> Result<f64, MaterialError> {
    finite(pressure.mul_add(volume - 1.0, finite(volume * elastic)?))
}

fn compliance(
    volume: f64,
    elastic: f64,
    treatment: VolumetricTreatment,
) -> Result<f64, MaterialError> {
    let value = finite(volume * elastic)?;
    if treatment == VolumetricTreatment::OriginalMaterial && value <= 0.0 {
        return Err(MaterialError::State {
            what: "thermal-reference pressure compliance underflows".into(),
        });
    }
    Ok(value)
}

fn pullback_state(
    thermal: &ThermalReference,
    pressure: f64,
    mut state: MixedPointState,
) -> Result<MixedPointState, MaterialError> {
    let volume = thermal.volume_ratio();
    state.potential_pa = finite(volume * state.potential_pa)?;
    state.remainder_potential_pa = finite(volume * state.remainder_potential_pa)?;
    state.bulk_stored_pa = finite(volume * state.bulk_stored_pa)?;
    state.bulk_complementary_pa = conjugate(volume, pressure, state.bulk_complementary_pa)?;
    state.pressure_residual = finite(volume * state.pressure_residual)?;
    state.volume_ratio = finite(volume * state.volume_ratio)?;
    state.bulk_volume_ratio = finite(volume * state.bulk_volume_ratio)?;
    if state.volume_ratio <= 0.0 || state.bulk_volume_ratio <= 0.0 {
        return Err(MaterialError::State {
            what: "thermal-reference mixed volumes must remain representably positive".into(),
        });
    }
    state.first_piola_pa = thermal.pullback_piola(state.first_piola_pa)?;
    Ok(state)
}

impl MixedFiniteStrainPotential for FiniteSolidMaterial {
    fn equilibrium_law(&self) -> &Hyperelastic {
        self.law.law()
    }

    fn remainder_stress_and_potential(
        &self,
        f: &[f64; 9],
    ) -> Result<([f64; 9], f64), MaterialError> {
        match &self.thermal {
            None => self.law.remainder_stress_and_potential(f),
            Some(thermal) => {
                let (p, w) = self
                    .law
                    .remainder_stress_and_potential(&thermal.elastic_deformation(f))?;
                Ok((
                    thermal.pullback_piola(p)?,
                    finite(thermal.volume_ratio() * w)?,
                ))
            }
        }
    }

    fn remainder_response(&self, f: &[f64; 9]) -> Result<RemainderResponse, MaterialError> {
        match &self.thermal {
            None => self.law.remainder_response(f),
            Some(thermal) => {
                let r = self
                    .law
                    .remainder_response(&thermal.elastic_deformation(f))?;
                Ok(RemainderResponse {
                    potential_pa: finite(thermal.volume_ratio() * r.potential_pa)?,
                    first_piola_pa: thermal.pullback_piola(r.first_piola_pa)?,
                    tangent_pa: thermal.pullback_tangent(r.tangent_pa)?,
                })
            }
        }
    }

    fn mixed_pressure_potential(
        &self,
        f: &[f64; 9],
        pressure_pa: f64,
        treatment: VolumetricTreatment,
    ) -> Result<PressurePotential, MaterialError> {
        match &self.thermal {
            None => self.law.mixed_pressure_potential(f, pressure_pa, treatment),
            Some(thermal) => {
                let elastic = self.law.mixed_pressure_potential(
                    &thermal.elastic_deformation(f), pressure_pa, treatment,
                )?;
                let volume = thermal.volume_ratio();
                let target = finite(volume * elastic.volume_ratio)?;
                if target <= 0.0 {
                    return Err(MaterialError::State {
                        what: "thermal-reference pressure volume underflows".into(),
                    });
                }
                Ok(PressurePotential {
                    complementary_pa: conjugate(volume, pressure_pa, elastic.complementary_pa)?,
                    bulk_stored_pa: finite(volume * elastic.bulk_stored_pa)?,
                    volume_ratio: target,
                    volume_increment: finite(
                        volume.mul_add(elastic.volume_increment, volume - 1.0),
                    )?,
                    compliance_per_pa: compliance(volume, elastic.compliance_per_pa, treatment)?,
                })
            }
        }
    }

    fn mixed_state(
        &self,
        f: &[f64; 9],
        pressure_pa: f64,
        treatment: VolumetricTreatment,
    ) -> Result<MixedPointState, MaterialError> {
        match &self.thermal {
            None => self.law.mixed_state(f, pressure_pa, treatment),
            Some(thermal) => pullback_state(
                thermal,
                pressure_pa,
                self.law.mixed_state(
                    &thermal.elastic_deformation(f),
                    pressure_pa,
                    treatment,
                )?,
            ),
        }
    }

    fn mixed_response(
        &self,
        f: &[f64; 9],
        pressure_pa: f64,
        treatment: VolumetricTreatment,
    ) -> Result<MixedPointResponse, MaterialError> {
        match &self.thermal {
            None => self.law.mixed_response(f, pressure_pa, treatment),
            Some(thermal) => {
                let r = self.law.mixed_response(
                    &thermal.elastic_deformation(f),
                    pressure_pa,
                    treatment,
                )?;
                let pressure_tangent = -compliance(
                    thermal.volume_ratio(), -r.pressure_tangent_per_pa, treatment,
                )?;
                Ok(MixedPointResponse {
                    state: pullback_state(thermal, pressure_pa, r.state)?,
                    tangent_pa: thermal.pullback_tangent(r.tangent_pa)?,
                    pressure_coupling: thermal.pullback_piola(r.pressure_coupling)?,
                    pressure_tangent_per_pa: pressure_tangent,
                })
            }
        }
    }
}
