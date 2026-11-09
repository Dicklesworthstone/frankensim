//! Nonisothermal plane-parallel absorption and emission with angular memory.
//!
//! Each layer has a uniform temperature and sourced absorption in inverse
//! metres. The default constructor uses refractive index one. The explicit
//! refractive constructor supplies one sourced real index to the stack.
//! Individual layers may instead declare different real indices, retaining
//! every smooth Fresnel interface and all internal reflections. There is no
//! scattering or coherent interference.
//! Thermal incident radiance is diffuse; `SlabBeamResponse` additionally
//! resolves collimated, polarized illumination without angular averaging.
//! Emergent hemispherical flux is resolved, but
//! is generally NOT Lambertian; a scalar view-factor consumer must explicitly
//! own its diffuse-port closure.
//!
//! At each positive direction cosine mu the exact ray-layer factors are
//! t_i = exp(-kappa_i*d_i/mu) and a_i = -expm1(-kappa_i*d_i/mu).
//! The SAME conserved transverse optical direction traverses every interface
//! before angular integration, with Snell-adjusted mu in each material.
//! Dividing a material layer does not introduce artificial angular mixing.
//! Positive hemispherical weights 2*mu*dmu retain the grazing-ray limit and
//! the optically thin diffuse absorptance 2*kappa*d.
//! Positive GL8/GL16 rules prepare boundary escape and symmetric interlayer
//! exchange coefficients once; subsequent heat applications need no quadrature.
//!
//! For a homogeneous slab the diffuse transmittance is 2*E_3(kappa*d), not
//! exp(-kappa*d). E_n is defined in NIST DLMF 8.19.3:
//! https://dlmf.nist.gov/8.19.E3
//! The embedded-rule discrepancy is an estimate of coefficient integration
//! error, not a certified bound on material, angular, spectral or spatial error.

use std::sync::Arc;

use fs_math::decay::{decayed_fraction, remaining_fraction};

use super::radiation_error;
use super::spectral::planck::{GL8, GL16};
use crate::ConductionError;

mod beam;
mod depth;
mod diffuse;
mod layered;
mod refraction;
mod thermal;
pub use beam::{SlabBeamExchange, SlabBeamResponse};
pub use depth::{SlabDepthReport, SlabDepthResponse};
pub use diffuse::{SlabDiffuseDepthReport, SlabDiffuseDepthResponse};
pub use thermal::{GraySlabThermal, MAX_SLAB_SPECTRAL_BANDS, SlabThermalBand, SpectralSlabThermal};

/// Maximum physical layers in one prepared slab.
pub const MAX_SLAB_LAYERS: usize = 64;
const MAX_ANGULAR_EVALUATIONS: usize = 65_536;

/// One non-scattering, uniform-absorption material layer in front-to-back order.
#[derive(Clone, Debug, PartialEq)]
pub struct AbsorbingLayer {
    thickness_m: f64,
    absorption_m_inv: f64,
    optical_depth: f64,
    temperature_support_k: [f64; 2],
    source: Arc<str>,
    refractive_index: f64,
    refractive_index_source: Option<Arc<str>>,
}

impl AbsorbingLayer {
    /// Coefficients are explicitly gray or constant over a caller-owned band;
    /// a total-gray material claim must not silently become spectral data.
    pub fn new(
        thickness_m: f64,
        absorption_m_inv: f64,
        temperature_support_k: [f64; 2],
        source: impl Into<String>,
    ) -> Result<Self, ConductionError> {
        let source = source.into();
        let [low, high] = temperature_support_k;
        if !positive(thickness_m)
            || !nonnegative(absorption_m_inv)
            || !positive(low)
            || !high.is_finite()
            || high <= low
            || source.trim().is_empty()
        {
            return Err(invalid(
                "need positive finite thickness, nonnegative finite absorption, a positive temperature interval and a source",
            ));
        }
        Ok(Self {
            thickness_m,
            absorption_m_inv,
            // Overflow of finite positive factors has the optically thick
            // limit. No infinity is accepted as an authored material property.
            optical_depth: thickness_m * absorption_m_inv,
            temperature_support_k,
            source: source.into(),
            refractive_index: 1.0,
            refractive_index_source: None,
        })
    }
    /// Declare this layer's real index for smooth incoherent Fresnel transport.
    /// The source must support this layer's temperature and gray/band scope.
    /// The default index-one layer has no inferred optical-material source.
    pub fn with_refractive_index(
        mut self,
        refractive_index: f64,
        source: impl Into<String>,
    ) -> Result<Self, ConductionError> {
        let source = source.into();
        if !refractive_index.is_finite()
            || refractive_index < 1.0
            || !(refractive_index * refractive_index).is_finite()
            || source.trim().is_empty()
        {
            return Err(invalid(
                "need a sourced real refractive index >= 1 with representable n²",
            ));
        }
        self.refractive_index = refractive_index;
        self.refractive_index_source = Some(source.into());
        Ok(self)
    }
    pub const fn refractive_index(&self) -> f64 {
        self.refractive_index
    }
    pub fn refractive_index_source(&self) -> Option<&str> {
        self.refractive_index_source.as_deref()
    }
    pub const fn thickness_m(&self) -> f64 {
        self.thickness_m
    }
    pub const fn absorption_m_inv(&self) -> f64 {
        self.absorption_m_inv
    }
    /// Nonnegative normal optical depth; infinity denotes the opaque limit.
    pub const fn optical_depth(&self) -> f64 {
        self.optical_depth
    }
    pub const fn temperature_support_k(&self) -> [f64; 2] {
        self.temperature_support_k
    }
    pub fn source(&self) -> &str {
        &self.source
    }
}

/// Preparation accuracy and work limits, independent of thermal solve gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlabQuadrature {
    /// Absolute embedded-rule discrepancy allowance per dimensionless coefficient.
    pub absolute_coefficient_error: f64,
    /// Relative allowance per coefficient; all retained integrands are positive.
    pub relative_coefficient_error: f64,
    /// Counts every ordinate, including rejected panels; at most 65,536.
    /// A completely transparent, index-matched slab needs zero evaluations.
    pub maximum_evaluations: usize,
    /// Maximum dyadic angular subdivision depth, from 1 through 48.
    pub maximum_depth: usize,
}
impl Default for SlabQuadrature {
    fn default() -> Self {
        Self {
            absolute_coefficient_error: 1e-10,
            relative_coefficient_error: 1e-8,
            maximum_evaluations: 16_384,
            maximum_depth: 40,
        }
    }
}

impl SlabQuadrature {
    fn is_valid(self) -> bool {
        nonnegative(self.absolute_coefficient_error)
            && nonnegative(self.relative_coefficient_error)
            && (self.absolute_coefficient_error > 0.0 || self.relative_coefficient_error > 0.0)
            && self.maximum_evaluations <= MAX_ANGULAR_EVALUATIONS
            && (1..=48).contains(&self.maximum_depth)
    }
}

/// Actual preparation work and numerical diagnostics, never a continuum certificate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlabQuadratureReport {
    pub controls: SlabQuadrature,
    pub angular_evaluations: usize,
    pub accepted_panels: usize,
    pub maximum_estimated_coefficient_error: f64,
    /// Largest absolute r + t + sum(layer escape fractions) - 1 over both sides.
    pub maximum_partition_residual: f64,
}

/// Immutable, source-bearing diffuse-incident slab response.
///
/// The default matched-index material has zero reflection. An explicitly
/// refractive stack includes all declared Fresnel interfaces and trapped rays.
/// Internal layer-to-layer heat still exists and is retained independently of
/// the two outgoing port fields. Clones share source and coefficient storage.
#[derive(Clone, Debug, PartialEq)]
pub struct AbsorbingSlab {
    layers: Arc<[AbsorbingLayer]>,
    exterior_refractive_index: f64,
    exterior_index_source: Option<Arc<str>>,
    reflectance: [f64; 2],
    transmittance: f64,
    face_coupling: Arc<[[f64; 2]]>,
    layer_coupling: Arc<[f64]>,
    quadrature: SlabQuadratureReport,
}

impl AbsorbingSlab {
    /// Prepare the declared layers, retaining every per-layer refractive index.
    /// Unqualified layers use the original matched-index model; use
    /// `AbsorbingLayer::with_refractive_index` to declare optical material data.
    pub fn new(
        layers: Vec<AbsorbingLayer>,
        controls: SlabQuadrature,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ConductionError> {
        Self::prepare(layers, 1.0, None, controls, cancelled)
    }

    /// A plane-parallel dielectric stack in an index-one environment.
    ///
    /// One real index n >= 1 applies to every layer and both polarizations.
    /// Its source must support the layers' declared temperature intervals and
    /// the caller's gray or finite-band scope. Absorption is the geometric
    /// Beer--Lambert bulk law; real-index Fresnel interfaces require negligible
    /// interfacial extinction. Conductors, within-layer index gradients, rough surfaces,
    /// coherent films and evanescent tunnelling are outside this model.
    ///
    /// Each polarization retains its own reflection series before averaging.
    /// Internally trapped thermal rays contribute to layer-to-layer exchange,
    /// with the n² radiance/etendue factor, even though they cannot escape.
    pub fn new_refractive(
        mut layers: Vec<AbsorbingLayer>,
        refractive_index: f64,
        index_source: impl Into<String>,
        controls: SlabQuadrature,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ConductionError> {
        let source = index_source.into();
        if !refractive_index.is_finite()
            || refractive_index < 1.0
            || !(refractive_index * refractive_index).is_finite()
            || source.trim().is_empty()
        {
            return Err(invalid(
                "need a sourced real refractive index >= 1 with representable n²",
            ));
        }
        let source: Arc<str> = source.into();
        for layer in &mut layers {
            layer.refractive_index = refractive_index;
            layer.refractive_index_source = Some(source.clone());
        }
        Self::prepare(layers, 1.0, None, controls, cancelled)
    }

    /// Prepare diffuse illumination and thermal emission in one lossless
    /// dielectric shared by both exterior faces. Layer indices stay absolute.
    ///
    /// Irradiance is physical flux in that medium. Layer blackbody exitances
    /// remain vacuum Planck/Stefan--Boltzmann values: escape emission includes
    /// the exterior n² radiance factor, while interlayer exchange retains each
    /// material's available optical modes. Equal-temperature equilibrium uses
    /// exterior irradiance n² times that common blackbody exitance.
    ///
    /// The exterior source must cover the layers' temperature and gray/band
    /// domain. Smooth incoherent interfaces, bulk Beer absorption and trapped
    /// propagating rays are supported; evanescent coupling is not represented.
    pub fn new_in_medium(
        layers: Vec<AbsorbingLayer>,
        exterior_refractive_index: f64,
        exterior_index_source: impl Into<String>,
        controls: SlabQuadrature,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ConductionError> {
        let source = exterior_index_source.into();
        if !exterior_refractive_index.is_finite()
            || exterior_refractive_index < 1.0
            || !(exterior_refractive_index * exterior_refractive_index).is_finite()
            || source.trim().is_empty()
        {
            return Err(invalid(
                "need a sourced real exterior index >= 1 with representable n²",
            ));
        }
        Self::prepare(
            layers,
            exterior_refractive_index,
            Some(source.into()),
            controls,
            cancelled,
        )
    }

    fn prepare(
        layers: Vec<AbsorbingLayer>,
        exterior_refractive_index: f64,
        exterior_index_source: Option<Arc<str>>,
        controls: SlabQuadrature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ConductionError> {
        checkpoint(&mut cancelled, 0)?;
        let n = layers.len();
        if !(1..=MAX_SLAB_LAYERS).contains(&n) || !controls.is_valid() {
            return Err(invalid(
                "invalid layer count or bounded angular quadrature controls",
            ));
        }
        let mut report = SlabQuadratureReport {
            controls,
            angular_evaluations: 0,
            accepted_panels: 0,
            maximum_estimated_coefficient_error: 0.0,
            maximum_partition_residual: 0.0,
        };
        let refractive_index = layers[0].refractive_index;
        let uniform_index = layers
            .iter()
            .all(|layer| layer.refractive_index == refractive_index);
        let general_index = !uniform_index || exterior_refractive_index != 1.0;
        if uniform_index
            && refractive_index == exterior_refractive_index
            && layers.iter().all(|layer| layer.optical_depth == 0.0)
        {
            checkpoint(&mut cancelled, 0)?;
            return Ok(Self {
                layers: layers.into(),
                exterior_refractive_index,
                exterior_index_source,
                reflectance: [0.0; 2],
                transmittance: 1.0,
                face_coupling: vec![[0.0; 2]; n].into(),
                layer_coupling: vec![0.0; n * n].into(),
                quadrature: report,
            });
        }
        // Only the upper triangular exchange entries are accumulated. Mirror
        // once at publication so reciprocity is bit-identical by construction.
        let channels = 3 + 2 * n + n * n;
        let mut total = vec![0.0; channels];
        let mut correction = vec![0.0; channels];
        let mut errors = vec![0.0; channels];
        let mut coarse = vec![0.0; channels];
        let mut fine = vec![0.0; channels];
        let mut rule_correction = vec![0.0; channels];
        let mut pending = Vec::new();
        let mut bands = Vec::new();
        if general_index {
            let mut edges = vec![0.0, exterior_refractive_index];
            edges.extend(layers.iter().map(AbsorbingLayer::refractive_index));
            edges.sort_by(f64::total_cmp);
            edges.dedup();
            bands.extend(edges.windows(2).map(|edge| [edge[0], edge[1]]));
            for band in (0..bands.len()).rev() {
                pending.push((0.0, 1.0, 0_usize, false, band));
            }
        } else if refractive_index > 1.0
            && n > 1
            && layers.iter().any(|layer| layer.optical_depth > 0.0)
        {
            let inverse = refractive_index.recip();
            let critical_mu = ((1.0 - inverse) * (1.0 + inverse)).sqrt();
            pending.push((0.0, critical_mu, 0_usize, true, 0));
        }
        // External incidence is integrated in external mu, avoiding a square
        // root singularity at the escape-cone edge. Trapped rays use internal mu.
        if !general_index {
            pending.push((0.0, 1.0, 0_usize, false, 0));
        }
        let integration_width: f64 = pending.iter().map(|(lo, hi, _, _, _)| hi - lo).sum();
        while let Some((lo, hi, depth, trapped, band)) = pending.pop() {
            checkpoint(&mut cancelled, report.angular_evaluations)?;
            if report.angular_evaluations + 24 > controls.maximum_evaluations {
                return Err(invalid(format!(
                    "angular evaluation budget exhausted: {} of {}",
                    report.angular_evaluations, controls.maximum_evaluations
                )));
            }
            for (nodes, output) in [(&GL8[..], &mut coarse), (&GL16[..], &mut fine)] {
                if general_index {
                    layered::rule(
                        &layers,
                        exterior_refractive_index,
                        bands[band],
                        lo,
                        hi,
                        nodes,
                        output,
                        &mut rule_correction,
                        &mut cancelled,
                        &mut report.angular_evaluations,
                        None,
                    )?;
                } else if refractive_index == 1.0 {
                    rule(
                        &layers,
                        lo,
                        hi,
                        nodes,
                        output,
                        &mut rule_correction,
                        &mut cancelled,
                        &mut report.angular_evaluations,
                    )?;
                } else {
                    refraction::rule(
                        &layers,
                        refractive_index,
                        lo,
                        hi,
                        trapped,
                        nodes,
                        output,
                        &mut rule_correction,
                        &mut cancelled,
                        &mut report.angular_evaluations,
                    )?;
                }
            }
            let admitted = coarse.iter().zip(&fine).all(|(&a, &b)| {
                (b - a).abs()
                    <= controls.absolute_coefficient_error * (hi - lo) / integration_width
                        + controls.relative_coefficient_error * b.abs()
            });
            if admitted {
                for i in 0..channels {
                    add(&mut total[i], &mut correction[i], fine[i]);
                    errors[i] += (fine[i] - coarse[i]).abs();
                }
                report.accepted_panels += 1;
            } else {
                if depth >= controls.maximum_depth {
                    return Err(invalid("angular subdivision depth exhausted"));
                }
                let mid = 0.5 * (lo + hi);
                if mid <= lo || mid >= hi {
                    return Err(invalid("angular panel cannot be subdivided representably"));
                }
                // Fixed left-first traversal, independent of any task scheduling.
                pending.push((mid, hi, depth + 1, trapped, band));
                pending.push((lo, mid, depth + 1, trapped, band));
            }
        }
        if total.iter().any(|&x| !nonnegative(x)) {
            return Err(invalid("nonfinite angular response"));
        }
        let transmittance = total[0];
        let reflectance = [total[channels - 2], total[channels - 1]];
        let face: Vec<[f64; 2]> = (0..n)
            .map(|i| [total[1 + 2 * i], total[2 + 2 * i]])
            .collect();
        let mut exchange = vec![0.0; n * n];
        for i in 0..n {
            for j in i + 1..n {
                let coefficient = total[1 + 2 * n + i * n + j];
                exchange[i * n + j] = coefficient;
                exchange[j * n + i] = coefficient;
            }
        }
        for side in 0..2 {
            let mut partition = reflectance[side] + transmittance;
            let mut c = 0.0;
            for a in &face {
                add(&mut partition, &mut c, a[side]);
            }
            report.maximum_partition_residual = report
                .maximum_partition_residual
                .max((partition - 1.0).abs());
        }
        if report.maximum_partition_residual > 128.0 * f64::EPSILON * (n + 1) as f64 {
            return Err(invalid(
                "angular response violates incident energy partition",
            ));
        }
        report.maximum_estimated_coefficient_error = errors.into_iter().fold(0.0_f64, f64::max);
        checkpoint(&mut cancelled, report.angular_evaluations)?;
        Ok(Self {
            layers: layers.into(),
            exterior_refractive_index,
            exterior_index_source,
            reflectance,
            transmittance,
            face_coupling: face.into(),
            layer_coupling: exchange.into(),
            quadrature: report,
        })
    }

    pub fn layer_count(&self) -> usize {
        self.layers.len()
    }
    pub fn layers(&self) -> &[AbsorbingLayer] {
        &self.layers
    }
    /// Absolute real index of the common exterior dielectric.
    pub const fn exterior_refractive_index(&self) -> f64 {
        self.exterior_refractive_index
    }
    /// Source of the exterior index; absent for the original index-one model.
    pub fn exterior_index_source(&self) -> Option<&str> {
        self.exterior_index_source.as_deref()
    }
    /// Common layer index, or None when the physical layers differ.
    /// Original per-layer indices and sources remain on `layers()`.
    pub fn refractive_index(&self) -> Option<f64> {
        let index = self.layers[0].refractive_index;
        self.layers
            .iter()
            .all(|layer| layer.refractive_index == index)
            .then_some(index)
    }
    /// Common source only when every layer declares the same index and source.
    /// Per-layer source access remains available for heterogeneous stacks.
    pub fn refractive_index_source(&self) -> Option<&str> {
        let source = self.layers[0].refractive_index_source.as_deref()?;
        (self.refractive_index().is_some()
            && self
                .layers
                .iter()
                .all(|layer| layer.refractive_index_source.as_deref() == Some(source)))
        .then_some(source)
    }
    /// Front diffuse incident reflectance. Use `reflectances()` for both sides:
    /// an absorbing heterogeneous stack can have unequal front/back reflection.
    pub const fn reflectance(&self) -> f64 {
        self.reflectance[0]
    }
    pub const fn reflectances(&self) -> [f64; 2] {
        self.reflectance
    }
    pub const fn transmittance(&self) -> f64 {
        self.transmittance
    }
    /// Layer-major [front, back] fractions of physical incident irradiance.
    /// Multiply by the exterior index squared for reciprocal emission
    /// coefficients relative to vacuum blackbody exitance.
    pub fn face_coupling(&self) -> &[[f64; 2]] {
        &self.face_coupling
    }
    /// Symmetric row-major interlayer exchange; its diagonal is zero.
    pub fn layer_coupling(&self) -> &[f64] {
        &self.layer_coupling
    }
    pub const fn quadrature(&self) -> &SlabQuadratureReport {
        &self.quadrature
    }

    /// Each physical layer keeps its own temperature and optical source support.
    pub fn validate_temperatures(&self, temperatures_k: &[f64]) -> Result<(), ConductionError> {
        if temperatures_k.len() != self.layer_count() {
            return Err(invalid(
                "slab temperature axis differs from original layers",
            ));
        }
        for (layer, &temperature) in self.layers.iter().zip(temperatures_k) {
            let [lo, hi] = layer.temperature_support_k;
            if !temperature.is_finite() || !(lo..=hi).contains(&temperature) {
                return Err(invalid(
                    "layer temperature is outside its optical source support",
                ));
            }
        }
        Ok(())
    }

    /// Layer emission escaping each face [W/m²]. Inputs are VACUUM blackbody
    /// exitances [W/m²]; the exterior n² radiance factor is included here.
    /// Finite signed deviations from ONE common background are supported for
    /// near-equilibrium solvers. This does not add transmitted external light.
    /// The consumer owns source-temperature and spectral-window admission.
    pub fn emission_to_faces(
        &self,
        blackbody_exitance: &[f64],
    ) -> Result<[f64; 2], ConductionError> {
        self.admit_fields(blackbody_exitance)?;
        let mut result = [0.0; 2];
        let mut correction = [0.0; 2];
        for (a, &b) in self.face_coupling.iter().zip(blackbody_exitance) {
            for side in 0..2 {
                add(
                    &mut result[side],
                    &mut correction[side],
                    finite((a[side] * self.exterior_refractive_index.powi(2)) * b)?,
                );
            }
        }
        for value in result {
            finite(value)?;
        }
        Ok(result)
    }

    /// Net heat INTO each layer per physical port area. Both external fields
    /// and layer exitances must use the SAME gray/band/per-log-wavelength scope
    /// and the SAME optional common background. No angular integration occurs.
    /// Output is written only after the complete finite result is admitted.
    pub fn layer_net_fluxes(
        &self,
        irradiation: [f64; 2],
        blackbody_exitance: &[f64],
        out: &mut [f64],
    ) -> Result<(), ConductionError> {
        self.admit_fields(blackbody_exitance)?;
        let n = self.layer_count();
        if out.len() != n || irradiation.iter().any(|x| !x.is_finite()) {
            return Err(invalid("invalid slab irradiation or layer heat axis"));
        }
        let mut heat = [0.0; MAX_SLAB_LAYERS];
        let mut correction = [0.0; MAX_SLAB_LAYERS];
        for i in 0..n {
            for side in 0..2 {
                let value = exterior_difference_product(
                    self.face_coupling[i][side],
                    irradiation[side],
                    blackbody_exitance[i],
                    self.exterior_refractive_index.powi(2),
                )?;
                add(&mut heat[i], &mut correction[i], value);
            }
            for j in i + 1..n {
                let into_i = difference_product(
                    self.layer_coupling[i * n + j],
                    blackbody_exitance[j],
                    blackbody_exitance[i],
                )?;
                add(&mut heat[i], &mut correction[i], into_i);
                add(&mut heat[j], &mut correction[j], -into_i);
            }
        }
        for &value in &heat[..n] {
            finite(value)?;
        }
        out.copy_from_slice(&heat[..n]);
        Ok(())
    }

    fn admit_fields(&self, values: &[f64]) -> Result<(), ConductionError> {
        if values.len() != self.layer_count() || values.iter().any(|x| !x.is_finite()) {
            Err(invalid(
                "need one finite blackbody exitance per original layer",
            ))
        } else {
            Ok(())
        }
    }
}

// q = 2*mu*dmu. No hemispherical averaging occurs inside this traversal.
#[allow(clippy::too_many_arguments)]
fn rule(
    layers: &[AbsorbingLayer],
    lo: f64,
    hi: f64,
    nodes: &[(f64, f64)],
    output: &mut [f64],
    correction: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
    evaluations: &mut usize,
) -> Result<(), ConductionError> {
    output.fill(0.0);
    correction.fill(0.0);
    let n = layers.len();
    let midpoint = 0.5 * (lo + hi);
    let half = 0.5 * (hi - lo);
    for &(node, weight) in nodes {
        for sign in [-1.0, 1.0] {
            checkpoint(cancelled, *evaluations)?;
            let mu = midpoint + sign * half * node;
            let q = 2.0 * mu * half * weight;
            let mut t = [0.0; MAX_SLAB_LAYERS];
            let mut a = [0.0; MAX_SLAB_LAYERS];
            let mut through = 1.0;
            for (i, layer) in layers.iter().enumerate() {
                let depth = layer.optical_depth / mu;
                t[i] = remaining_fraction(depth);
                a[i] = decayed_fraction(depth);
                let index = 1 + 2 * i;
                add(
                    &mut output[index],
                    &mut correction[index],
                    q * a[i] * through,
                );
                through *= t[i];
            }
            add(&mut output[0], &mut correction[0], q * through);
            through = 1.0;
            for i in (0..n).rev() {
                let index = 2 + 2 * i;
                add(
                    &mut output[index],
                    &mut correction[index],
                    q * a[i] * through,
                );
                through *= t[i];
            }
            for i in 0..n {
                if a[i] == 0.0 {
                    continue;
                }
                let mut between = 1.0;
                for j in i + 1..n {
                    let index = 1 + 2 * n + i * n + j;
                    add(
                        &mut output[index],
                        &mut correction[index],
                        q * a[i] * a[j] * between,
                    );
                    between *= t[j];
                }
            }
            *evaluations += 1;
        }
    }
    Ok(())
}

fn add(sum: &mut f64, correction: &mut f64, value: f64) {
    let y = value - *correction;
    let next = *sum + y;
    *correction = (next - *sum) - y;
    *sum = next;
}

// Keep exact equilibrium cancellation before scaling, with a finite-product
// fallback when the exterior blackbody flux overflows but its escape is tiny.
fn exterior_difference_product(
    coefficient: f64,
    incoming: f64,
    source: f64,
    index_squared: f64,
) -> Result<f64, ConductionError> {
    let exterior_source = index_squared * source;
    if coefficient == 0.0 {
        Ok(0.0)
    } else if exterior_source.is_finite() {
        difference_product(coefficient, incoming, exterior_source)
    } else {
        finite(coefficient * incoming - (coefficient * index_squared) * source)
    }
}

fn difference_product(
    coefficient: f64,
    incoming: f64,
    source: f64,
) -> Result<f64, ConductionError> {
    if coefficient == 0.0 {
        return Ok(0.0);
    }
    let difference = incoming - source;
    finite(if difference.is_finite() {
        coefficient * difference
    } else {
        coefficient * incoming - coefficient * source
    })
}

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn nonnegative(x: f64) -> bool {
    x.is_finite() && x >= 0.0
}
fn finite(x: f64) -> Result<f64, ConductionError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(invalid("unrepresentable slab energy arithmetic"))
    }
}
fn invalid(reason: impl Into<String>) -> ConductionError {
    radiation_error(
        "layered-slab-radiation",
        reason,
        "supply sourced nonscattering layers, an admitted refractive model, diffuse incident fields and sufficient angular work",
    )
}
fn checkpoint(cancelled: &mut impl FnMut() -> bool, at: usize) -> Result<(), ConductionError> {
    if cancelled() {
        Err(ConductionError::Cancelled {
            stage: "slab-angular-prepare",
            at,
        })
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
