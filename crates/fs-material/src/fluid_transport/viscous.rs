//! Three-dimensional viscous response on the original transport/rheology laws.
//!
//! With D = sym(grad u), d = dev(D) and theta = tr(D),
//! tau = 2 mu d + zeta theta I and tau:D = 2 mu d:d + zeta theta^2.
//! Bulk viscosity is an independent nonnegative input; choosing zero is explicit.
//! Thermodynamic pressure is not included in this dissipative stress.
//! The generalized-Newtonian extension evaluates the existing scalar law at
//! g = sqrt(2 d:d), agreeing with the incompressible shear convention without
//! treating volume change as shear thinning.
//!
//! This isotropic split follows the Newtonian viscous law described in
//! COMSOL's "Viscous Damping: Space-Fixed Frame":
//! <https://doc.comsol.com/6.4/doc/com.comsol.help.rotor/rotor_ug_theory.4.13.html>.
//! Applying a shear-calibrated law to general flow is an explicit constitutive
//! assumption, not an extensional-flow calibration or a viscoelastic model.
//! All responses use physical SI units and tensor-Voigt xy/yz/zx ordering.
//! Local evaluation and tangent action need no allocation, finite differences,
//! source refit or gradient stencil.

use super::FluidTransport;
use super::rheology::{ThermalCarreauYasuda, ThermalRheologyState};
use crate::tensor::contract;
use crate::{MaterialError, Tangent6, Voigt};
use fs_math::det;

#[cfg(test)]
mod tests;

fn invalid(what: &'static str) -> MaterialError {
    MaterialError::State { what: what.into() }
}

/// Symmetric velocity gradient and its physical shear/volume decomposition.
/// Components have units 1/s. Rigid spin is discarded on entry.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DeformationRate {
    tensor: Voigt,
    deviator: Voigt,
    dilatation: f64,
    shear_rate: f64,
}

impl DeformationRate {
    /// Admit gradient[i][j] = du_i/dx_j in 1/s.
    pub fn from_velocity_gradient(gradient: [[f64; 3]; 3]) -> Result<Self, MaterialError> {
        if gradient.iter().flatten().any(|v| !v.is_finite()) {
            return Err(invalid("finite physical velocity gradient required"));
        }
        Self::from_tensor([
            gradient[0][0], gradient[1][1], gradient[2][2],
            0.5 * gradient[0][1] + 0.5 * gradient[1][0],
            0.5 * gradient[1][2] + 0.5 * gradient[2][1],
            0.5 * gradient[2][0] + 0.5 * gradient[0][2],
        ])
    }

    /// Admit a physical tensor in xx/yy/zz/xy/yz/zx order.
    /// Shear entries are D_ij, not the engineering rate 2 D_ij.
    pub fn from_tensor(tensor: Voigt) -> Result<Self, MaterialError> {
        if tensor.iter().any(|v| !v.is_finite()) {
            return Err(invalid("finite tensor deformation rate required"));
        }
        let mean = tensor[0] / 3.0 + tensor[1] / 3.0 + tensor[2] / 3.0;
        let dilatation = 3.0 * mean;
        // Pairwise differences retain a small deviator on a large spherical
        // rate, and give exactly zero for isotropic expansion/compression.
        let deviator = [
            (tensor[0] - tensor[1]) / 3.0 + (tensor[0] - tensor[2]) / 3.0,
            (tensor[1] - tensor[0]) / 3.0 + (tensor[1] - tensor[2]) / 3.0,
            (tensor[2] - tensor[0]) / 3.0 + (tensor[2] - tensor[1]) / 3.0,
            tensor[3], tensor[4], tensor[5],
        ];
        let scale = deviator.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
        let shear_rate = if scale == 0.0 { 0.0 } else {
            let unit = deviator.map(|x| x / scale);
            scale * det::sqrt(2.0 * contract(&unit, &unit))
        };
        if !dilatation.is_finite() || !shear_rate.is_finite()
            || deviator.iter().any(|v| !v.is_finite())
        {
            return Err(invalid("unrepresentable deformation-rate invariants"));
        }
        Ok(Self { tensor, deviator, dilatation, shear_rate })
    }

    /// Complete symmetric physical strain-rate tensor [1/s].
    pub const fn tensor_per_s(&self) -> &Voigt { &self.tensor }
    /// Trace-free symmetric physical strain-rate tensor [1/s].
    pub const fn deviator_per_s(&self) -> &Voigt { &self.deviator }
    /// Physical velocity divergence [1/s].
    pub const fn dilatation_per_s(&self) -> f64 { self.dilatation }
    /// sqrt(2 dev(D):dev(D)); equals |du/dy| for simple shear.
    pub const fn shear_rate_per_s(&self) -> f64 { self.shear_rate }
}

/// Checked Newtonian coefficients. A material source owns their provenance;
/// this small numerical kernel does not create experimental authority.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NewtonianViscosity {
    shear_pa_s: f64,
    bulk_pa_s: f64,
}

impl NewtonianViscosity {
    /// Admit positive shear and independently nonnegative bulk viscosity [Pa s].
    pub fn new(shear_pa_s: f64, bulk_pa_s: f64) -> Result<Self, MaterialError> {
        if !shear_pa_s.is_finite() || shear_pa_s <= 0.0
            || !bulk_pa_s.is_finite() || bulk_pa_s < 0.0
        {
            return Err(invalid("positive shear viscosity and nonnegative bulk viscosity required"));
        }
        Ok(Self { shear_pa_s, bulk_pa_s })
    }

    /// Consume the original evaluated gas/liquid transport, preserving all of
    /// its admission failures. Bulk viscosity is never inferred from mu or cp.
    pub fn from_transport<F: FluidTransport + ?Sized>(
        source: &F,
        bulk_pa_s: f64,
    ) -> Result<Self, MaterialError> {
        Self::new(source.transport_state()?.dynamic_viscosity_pa_s(), bulk_pa_s)
    }

    /// Original dynamic shear viscosity [Pa s].
    pub const fn shear_viscosity_pa_s(&self) -> f64 { self.shear_pa_s }
    /// Independently authored bulk viscosity [Pa s].
    pub const fn bulk_viscosity_pa_s(&self) -> f64 { self.bulk_pa_s }

    /// Evaluate local stress and irreversible power at the supplied physical rate.
    pub fn evaluate(&self, rate: DeformationRate) -> Result<ViscousResponse, MaterialError> {
        ViscousResponse::evaluate(rate, self.shear_pa_s, self.bulk_pa_s, self.shear_pa_s)
    }
}

/// Local physical viscous stress, heating and its exact material derivative.
/// A Newtonian source is frozen at its evaluated temperature/composition.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ViscousResponse {
    rate: DeformationRate,
    shear_pa_s: f64,
    bulk_pa_s: f64,
    differential_shear_pa_s: f64,
    deviatoric_stress: Voigt,
    isotropic_stress: f64,
    stress: Voigt,
    dissipation: f64,
}

impl ViscousResponse {
    fn evaluate(rate: DeformationRate, mu: f64, bulk: f64, differential: f64)
        -> Result<Self, MaterialError>
    {
        let deviatoric_stress = rate.deviator.map(|d| (mu * d) * 2.0);
        let isotropic_stress = bulk * rate.dilatation;
        let stress: Voigt = core::array::from_fn(|i| {
            deviatoric_stress[i] + if i < 3 { isotropic_stress } else { 0.0 }
        });
        // Keep shear and bulk work separate: no cancellation against pressure
        // or a large spherical rate can hide irreversible shear heating.
        let dissipation = contract(&deviatoric_stress, &rate.deviator)
            + isotropic_stress * rate.dilatation;
        if stress.iter().any(|v| !v.is_finite()) || !dissipation.is_finite()
            || dissipation < 0.0
        {
            return Err(invalid("unrepresentable physical viscous stress or heating"));
        }
        Ok(Self { rate, shear_pa_s: mu, bulk_pa_s: bulk,
            differential_shear_pa_s: differential, deviatoric_stress,
            isotropic_stress, stress, dissipation })
    }

    /// Original admitted rate.
    pub const fn rate(&self) -> &DeformationRate { &self.rate }
    /// Viscous Cauchy stress, excluding the separate thermodynamic -p I.
    pub const fn stress_pa(&self) -> &Voigt { &self.stress }
    /// Physical deviatoric stress [Pa].
    pub const fn deviatoric_stress_pa(&self) -> &Voigt { &self.deviatoric_stress }
    /// Signed spherical viscous stress [Pa].
    pub const fn isotropic_stress_pa(&self) -> f64 { self.isotropic_stress }
    /// Physical irreversible power per volume [W/m^3], never numerical loss.
    pub const fn dissipation_w_m3(&self) -> f64 { self.dissipation }
    /// Evaluated dynamic shear viscosity [Pa s].
    pub const fn shear_viscosity_pa_s(&self) -> f64 { self.shear_pa_s }
    /// Authored bulk viscosity [Pa s].
    pub const fn bulk_viscosity_pa_s(&self) -> f64 { self.bulk_pa_s }

    /// Matrix-free derivative of stress with respect to symmetric D. The
    /// input uses tensor shear, so its contraction carries the factor two.
    /// Decomposing along/perpendicular to dev(D) uses the original positive
    /// d(mu*g)/dg rather than subtracting almost equal thinning coefficients.
    pub fn apply_tangent(&self, increment: &Voigt) -> Result<Voigt, MaterialError> {
        let increment = DeformationRate::from_tensor(*increment)?;
        let direction = if self.rate.shear_rate == 0.0
            || self.differential_shear_pa_s == self.shear_pa_s
        { [0.0; 6] } else {
            self.rate.deviator.map(|x| (x / self.rate.shear_rate) * core::f64::consts::SQRT_2)
        };
        let along = contract(&direction, &increment.deviator);
        let spherical = self.bulk_pa_s * increment.dilatation;
        let result: Voigt = core::array::from_fn(|i| {
            let parallel = along * direction[i];
            2.0 * (self.shear_pa_s * (increment.deviator[i] - parallel)
                + self.differential_shear_pa_s * parallel)
                + if i < 3 { spherical } else { 0.0 }
        });
        if result.iter().any(|v| !v.is_finite()) {
            return Err(invalid("unrepresentable viscous tangent action"));
        }
        Ok(result)
    }

    /// Optional 6x6 tensor-Voigt derivative [Pa s]. Normal/shear coupling is
    /// self-adjoint under tensor contraction, not ordinary matrix transpose.
    pub fn tangent_pa_s(&self) -> Result<Tangent6, MaterialError> {
        let mut tangent = [[0.0; 6]; 6];
        for column in 0..6 {
            let mut unit = [0.0; 6];
            unit[column] = 1.0;
            let values = self.apply_tangent(&unit)?;
            for row in 0..6 { tangent[row][column] = values[row]; }
        }
        Ok(tangent)
    }

    /// Resolve traction on an explicitly supplied unit normal. Tangential
    /// traction uses the deviator alone, so large bulk stress cannot contaminate
    /// physical shear through pressure cancellation.
    pub fn traction(&self, normal: [f64; 3]) -> Result<ViscousTraction, MaterialError> {
        let norm2 = normal.iter().map(|x| x * x).sum::<f64>();
        if !norm2.is_finite() || (norm2 - 1.0).abs() > 1e-12 {
            return Err(invalid("viscous traction requires a finite unit normal"));
        }
        let s = &self.deviatoric_stress;
        let deviatoric = [
            s[0] * normal[0] + s[3] * normal[1] + s[5] * normal[2],
            s[3] * normal[0] + s[1] * normal[1] + s[4] * normal[2],
            s[5] * normal[0] + s[4] * normal[1] + s[2] * normal[2],
        ];
        let normal_deviatoric = (0..3).map(|i| normal[i] * deviatoric[i]).sum::<f64>();
        let result = ViscousTraction {
            vector_pa: core::array::from_fn(|i| deviatoric[i] + self.isotropic_stress * normal[i]),
            normal_pa: normal_deviatoric + self.isotropic_stress,
            shear_pa: core::array::from_fn(|i| deviatoric[i] - normal_deviatoric * normal[i]),
        };
        if !result.normal_pa.is_finite()
            || result.vector_pa.iter().chain(&result.shear_pa).any(|v| !v.is_finite())
        {
            return Err(invalid("unrepresentable viscous traction"));
        }
        Ok(result)
    }
}

/// Oriented viscous traction [Pa]. Reversing the normal reverses both vectors;
/// normal_pa is the signed tensile normal stress and retains its sign.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ViscousTraction {
    /// Complete viscous traction tau*n.
    pub vector_pa: [f64; 3],
    /// Signed normal component n*tau*n.
    pub normal_pa: f64,
    /// Tangential viscous traction.
    pub shear_pa: [f64; 3],
}

/// The original temperature/rate evaluation accompanies the complete tensor
/// response, retaining support, temperature-shift choice and derivatives.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalViscousResponse {
    /// Original scalar material evaluation at the actual deviatoric rate.
    pub rheology: ThermalRheologyState,
    /// Complete local tensor response.
    pub viscous: ViscousResponse,
    /// d(tau)/dT at fixed physical D and fixed authored bulk viscosity [Pa/K].
    pub temperature_tangent_pa_per_k: Voigt,
}

impl ThermalCarreauYasuda {
    /// Evaluate the original bounded shear law and its exact tensor tangent on
    /// dev(D). Bulk viscosity remains an independent fixed input. No previous
    /// scalar rheology sample is reused at a different shear rate.
    pub fn viscous_response(
        &self,
        temperature_k: f64,
        rate: DeformationRate,
        bulk_viscosity_pa_s: f64,
    ) -> Result<ThermalViscousResponse, MaterialError> {
        let rheology = self.evaluate(temperature_k, rate.shear_rate)?;
        let coefficients = NewtonianViscosity::new(
            rheology.dynamic_viscosity_pa_s(), bulk_viscosity_pa_s,
        )?;
        let viscous = ViscousResponse::evaluate(rate, coefficients.shear_pa_s,
            coefficients.bulk_pa_s, rheology.stress_tangent_pa_s())?;
        let temperature_tangent_pa_per_k = rate.deviator.map(|d| {
            (rheology.temperature_derivative_pa_s_per_k() * d) * 2.0
        });
        if temperature_tangent_pa_per_k.iter().any(|v| !v.is_finite()) {
            return Err(invalid("unrepresentable viscous temperature tangent"));
        }
        Ok(ThermalViscousResponse { rheology, viscous, temperature_tangent_pa_per_k })
    }
}
