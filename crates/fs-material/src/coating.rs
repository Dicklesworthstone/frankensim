//! Source-bound, isotropic viscoelastic films on layered structures.
//!
//! A measured uniaxial film spectrum is not a measured three-dimensional
//! tensor. The caller explicitly supplies a real, frequency-independent
//! Poisson ratio and a small-strain envelope. These assumptions stay attached
//! to the film; no missing mechanical, thermal, curing or moisture data is
//! inferred from a material name. The source's conditioning and frequency
//! support remain authoritative.

use crate::elastic::OrthotropicElastic;
use crate::laminate::{LaminaLayer, LayerLaw};
use crate::state_point::{ResolvedMaterialStatePoint, ScalarAdmissibility, ScalarPropertyRequirement};
use crate::visco::modified_zener::ModifiedFractionalZener;
use crate::visco::{GeneralizedMaxwell, LoweredModel, TensorProny};
use crate::MaterialError;
use fs_qty::{Density, Dims, Pressure, QuantitySpec};

/// Density and the five published modified-fractional-Zener fit coefficients.
pub const COATING_PROPERTIES: [(&str, Dims); 6] = [
    ("density", Density::DIMS),
    ("varnish_e0", Pressure::DIMS),
    ("varnish_einf", Pressure::DIMS),
    ("varnish_alpha", Dims([0; 6])),
    ("varnish_beta", Dims([0; 6])),
    ("varnish_tau", Dims([0, 0, 1, 0, 0, 0])),
];

fn invalid(message: impl Into<String>) -> MaterialError {
    MaterialError::Parameters { what: message.into() }
}

/// A frozen material source and an explicitly authored isotropic tensor lift.
#[derive(Clone, Debug)]
pub struct ResolvedCoating {
    source: ResolvedMaterialStatePoint,
    law: ModifiedFractionalZener,
    density_kg_m3: f64,
    poisson_ratio: f64,
    strain_limit: f64,
}

/// One frequency-domain slice of a film. Its viscosity is E″/ω, so its
/// complex constitutive response agrees at the declared reference frequency.
/// Using this slice over a broad band is a separate approximation, not the
/// original hereditary law.
#[derive(Clone, Debug)]
pub struct CoatingFrequencySlice {
    /// Retained source, density, thickness and storage-modulus tensor.
    pub layer: LaminaLayer,
    /// Principal reduced plane-stress viscosity [Pa s], engineering shear.
    pub plane_viscosity_pa_s: [f64; 9],
    /// Principal G13 and G23 viscosities [Pa s].
    pub transverse_shear_viscosity_pa_s: [f64; 2],
    /// The frequency at which this slice matches the source spectrum [rad/s].
    pub reference_omega_rad_s: f64,
}

impl ResolvedCoating {
    /// Requirements for the existing common material resolver. All six
    /// coefficients are dimensioned source data, never fallback defaults.
    pub fn requirements() -> Result<Vec<ScalarPropertyRequirement>, MaterialError> {
        COATING_PROPERTIES.iter().map(|&(name, dims)| {
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
                .map_err(|e| invalid(e.to_string()))
        }).collect()
    }

    /// Bind an existing resolved source. The fitted coefficient claims must
    /// declare one intersecting frequency band; density need not depend on
    /// frequency. Typed fixed-environment curves may supply the fit parameters;
    /// frequency-dependent coefficient samples cannot replace the fitted law.
    /// The source's original query includes its actual conditioning coordinates,
    /// which are retained unchanged. Existing fixed-temperature film sources
    /// keep their original applicability-only conditioning and provenance flags.
    pub fn from_material(source: ResolvedMaterialStatePoint, poisson_ratio: f64,
        strain_limit: f64) -> Result<Self, MaterialError>
    {
        if !poisson_ratio.is_finite() || poisson_ratio <= -1.0 || poisson_ratio >= 0.5
            || !strain_limit.is_finite() || strain_limit <= 0.0 {
            return Err(invalid("coating tensor lift requires explicit -1 < nu < 0.5 and positive strain limit"));
        }
        let mut values = [0.0; 6];
        let mut band = (0.0_f64, f64::INFINITY);
        for (i, &(name, dims)) in COATING_PROPERTIES.iter().enumerate() {
            let p = source.property(name).ok_or_else(|| invalid(format!("missing coating property {name}")))?;
            // Legacy measured films retain applicability-only source flags;
            // additional response axes need the shared typed environment policy.
            let held = p.is_constant_at_fixed_temperature() || p.is_constant_at_fixed_environment();
            if p.requirement().quantity() != QuantitySpec::dimensional(dims)
                || !held || !p.value_si().is_finite() || p.value_si() <= 0.0 {
                return Err(invalid("coating coefficients must be positive fixed-state source scalars with original units"));
            }
            let frequency = p.angular_frequency_band_rad_s().map_err(|e| invalid(e.to_string()))?;
            if let Some((lo, hi)) = frequency {
                band.0 = band.0.max(lo); band.1 = band.1.min(hi);
            } else if i != 0 {
                return Err(invalid("every coating fit coefficient requires explicit frequency support"));
            }
            values[i] = p.value_si();
        }
        if !band.0.is_finite() || !band.1.is_finite() || band.0 <= 0.0 || band.1 <= band.0 {
            return Err(invalid("coating source frequency supports do not intersect"));
        }
        let [density, e0, efit, alpha, beta, tau] = values;
        let law = ModifiedFractionalZener::new(e0, efit, alpha, beta, tau,
            (band.0/core::f64::consts::TAU, band.1/core::f64::consts::TAU))?;
        Ok(Self { source, law, density_kg_m3: density, poisson_ratio, strain_limit })
    }

    /// Original source object, including all individual usage receipts.
    pub const fn source(&self) -> &ResolvedMaterialStatePoint { &self.source }
    /// Full measured-fit frequency law, with its original band.
    pub const fn law(&self) -> &ModifiedFractionalZener { &self.law }
    /// Film-group density at the retained conditioning state [kg/m³].
    pub const fn density_kg_m3(&self) -> f64 { self.density_kg_m3 }
    /// Authored real Poisson ratio; this is not a measured source property.
    pub const fn poisson_ratio(&self) -> f64 { self.poisson_ratio }
    /// Authored small-strain envelope.
    pub const fn strain_limit(&self) -> f64 { self.strain_limit }

    /// Equilibrium layer of the declared fitted constitutive model. This uses
    /// its explicit E0 coefficient, not an out-of-band evaluation of E*(omega).
    /// Transient use is an authored causal continuation; the original measured
    /// fit band and conditioning remain unchanged on `source` and `law`.
    pub fn equilibrium_layer(&self, thickness_m: f64) -> Result<LaminaLayer, MaterialError> {
        if !thickness_m.is_finite() || thickness_m <= 0.0 {
            return Err(invalid("coating layer thickness must be finite and positive"));
        }
        let e = self.law.equilibrium_modulus_pa();
        let g = e / (2.0 * (1.0 + self.poisson_ratio));
        let elastic = OrthotropicElastic::new([e; 3], [self.poisson_ratio; 3],
            [g; 3], self.strain_limit)?;
        Ok(LaminaLayer { thickness_m, orientation_rad: 0.0,
            density_kg_m3: self.density_kg_m3, material: Some(self.source.clone()),
            law: LayerLaw::Elastic(elastic) })
    }

    /// Form a passive equivalent elastic/viscous section at one explicitly
    /// selected frequency. Layer geometry remains an independent input.
    pub fn at_frequency(&self, thickness_m: f64, omega: f64)
        -> Result<CoatingFrequencySlice, MaterialError>
    {
        if !thickness_m.is_finite() || thickness_m <= 0.0 {
            return Err(invalid("coating layer thickness must be finite and positive"));
        }
        let (storage, loss) = self.law.modulus(omega)?;
        let nu = self.poisson_ratio;
        let shear = storage / (2.0*(1.0+nu));
        let elastic = OrthotropicElastic::new([storage;3], [nu;3], [shear;3], self.strain_limit)?;
        let viscosity = loss/omega;
        let q = viscosity/(1.0-nu*nu);
        let g = viscosity/(2.0*(1.0+nu));
        if !q.is_finite() || !g.is_finite() || q < 0.0 || g < 0.0 {
            return Err(invalid("unrepresentable coating loss tensor"));
        }
        Ok(CoatingFrequencySlice {
            layer: LaminaLayer { thickness_m, orientation_rad:0.0,
                density_kg_m3:self.density_kg_m3, material:Some(self.source.clone()),
                law:LayerLaw::Elastic(elastic) },
            plane_viscosity_pa_s:[q,nu*q,0.0,nu*q,q,0.0,0.0,0.0,g],
            transverse_shear_viscosity_pa_s:[g;2], reference_omega_rad_s:omega,
        })
    }

    /// Lower the original scalar spectrum into the existing passive tensor
    /// Prony owner. Bulk and shear share the declared constant Poisson ratio.
    /// Both storage and loss errors must meet their own authored tolerances.
    /// This returns material memory; it does not silently substitute a frozen
    /// frequency slice or claim that a shell runtime already consumes it.
    pub fn lower_to_tensor_prony(&self, terms: usize, storage_relative_tolerance: f64,
        loss_relative_tolerance: f64) -> Result<(TensorProny, LoweredModel), MaterialError>
    {
        if terms == 0 || terms > 32 {
            return Err(invalid("isotropic coating lowering admits 1..=32 terms (at most 64 tensor branches)"));
        }
        let lowered = self.law.lower_to_prony(terms, storage_relative_tolerance, loss_relative_tolerance)?;
        let scale = |factor:f64| GeneralizedMaxwell::new(lowered.model.e_inf*factor,
            lowered.model.terms.iter().map(|&(e,t)| (e*factor,t)).collect());
        let bulk = scale(1.0/(3.0*(1.0-2.0*self.poisson_ratio)))?;
        let shear = scale(1.0/(2.0*(1.0+self.poisson_ratio)))?;
        let tensor = TensorProny::isotropic(bulk,shear)?.with_band(lowered.band)?;
        Ok((tensor,lowered))
    }
}
