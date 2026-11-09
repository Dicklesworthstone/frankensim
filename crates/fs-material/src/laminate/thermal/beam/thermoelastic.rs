//! Reversible thermal potential of the original free-sided narrow beam.
//!
//! The original affine-source lamina potential supplies all material derivatives.
//! Its transverse strains relax with temperature, so its fixed-full-strain heat
//! capacity is not the beam capacity. Writing q=[epsilon_y,gamma_xy] and
//! b=d(sigma)/dT at fixed full strain, local stationarity gives
//! `Psi_bar_TT = Psi_TT - V*b_q^T*Q_qq^-1*b_q` and
//! `d(sigma_x)/dT = b_x - Q_xq*Q_qq^-1*b_q`. Entropy and internal energy retain
//! the original potential; the capacity receives the positive Schur correction.
//!
//! Each original layer has its own temperature. All outputs are per ORIGINAL
//! beam length; the consumer supplies original length quadrature, the caloric
//! potential and heat transport. These signed elastic contributions do not add
//! a second sensible-heat inventory or implement thermal time integration.
//! The complete smooth source interval must be admitted by the existing affine
//! owner. Principal axes must align with the beam (within rotation roundoff),
//! or both elastic and expansion source symmetries must be explicitly isotropic.
//! Accidental off-axis alignment at one temperature cannot admit a thermal path.

use super::PreparedLaminateBeam;
use crate::elastic::thermal::ThermoelasticError;
use crate::laminate::thermal::LaminateThermalError;
use crate::laminate::thermal::source::{
    LaminaElasticSource, LaminaExpansionSource, ResolvedLaminateThermal,
};
use crate::laminate::thermal::thermoelastic::{
    AffineLaminaThermoelastic, ResolvedLaminaThermoelastic,
};
use fs_math::det;
use std::sync::Arc;

/// Refusal before a complete potential or thermal response is exposed.
#[derive(Debug)]
pub enum LaminateBeamThermoelasticError {
    /// The original layer's smooth source potential refused admission.
    Source {
        layer: usize,
        source: ThermoelasticError,
    },
    /// The original narrow-beam reduction or recovery refused admission.
    Beam(LaminateThermalError),
    /// Source identity, planar family, strain or arithmetic is inadmissible.
    Invalid(&'static str),
    /// Cancellation occurred before a complete result was available.
    Cancelled,
}
impl core::fmt::Display for LaminateBeamThermoelasticError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-MATERIAL-LAMINATE-BEAM-THERMOELASTIC: {self:?}")
    }
}
impl std::error::Error for LaminateBeamThermoelasticError {}

#[derive(Clone, Debug)]
struct Layer {
    resolved: ResolvedLaminaThermoelastic,
    sin: f64,
    cos: f64,
    limit: f64,
    minimum_capacity: f64,
    // Scaled Cholesky metric for a nonnegative b_q^T Q_qq^-1 b_q.
    lateral_scale: [f64; 2],
    lateral_ratio: f64,
    lateral_root: f64,
}

/// Current source-pinned potential of the same admitted beam closure.
/// Preparation is O(layers); evaluation allocates nothing and performs neither
/// source queries nor trigonometry. Clones retain the original prepared Arc.
#[derive(Clone, Debug)]
pub struct ResolvedLaminateBeamThermoelastic {
    current: ResolvedLaminateThermal,
    beam: PreparedLaminateBeam,
    layers: Vec<Layer>,
}

/// One layer, or the sum of all original layers, per reference beam length.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct LaminateBeamThermoelasticState {
    pub axial_force_n: f64,
    pub bending_moment_n_m: f64,
    pub free_energy_j_m: f64,
    pub entropy_j_m_k: f64,
    /// U = Psi + sum_i T_i S_i; different layer temperatures remain distinct.
    pub internal_energy_j_m: f64,
    /// Signed elastic capacity at fixed eps0/kappa, with transverse relaxation.
    /// For a whole-section state this is sum_i C_i, not a single-temperature law.
    pub heat_capacity_j_m_k: f64,
    /// d[N,M]/dT_i in [N/K,N m/K] for one layer. The whole-section sum is the
    /// response when every original layer temperature increases by one kelvin.
    /// The corresponding entropy derivatives in eps0/kappa are its NEGATIVE.
    pub generalized_force_temperature_derivative: [f64; 2],
}
impl LaminateBeamThermoelasticState {
    fn add(&mut self, other: Self) {
        self.axial_force_n += other.axial_force_n;
        self.bending_moment_n_m += other.bending_moment_n_m;
        self.free_energy_j_m += other.free_energy_j_m;
        self.entropy_j_m_k += other.entropy_j_m_k;
        self.internal_energy_j_m += other.internal_energy_j_m;
        self.heat_capacity_j_m_k += other.heat_capacity_j_m_k;
        for (a, b) in self
            .generalized_force_temperature_derivative
            .iter_mut()
            .zip(other.generalized_force_temperature_derivative)
        {
            *a += b;
        }
    }
    fn checked(self) -> Result<Self, LaminateBeamThermoelasticError> {
        if [
            self.axial_force_n,
            self.bending_moment_n_m,
            self.free_energy_j_m,
            self.entropy_j_m_k,
            self.internal_energy_j_m,
            self.heat_capacity_j_m_k,
        ]
        .into_iter()
        .chain(self.generalized_force_temperature_derivative)
        .all(f64::is_finite)
            && self.free_energy_j_m >= 0.0
        {
            Ok(self)
        } else {
            Err(invalid("beam thermodynamic arithmetic is unrepresentable"))
        }
    }
}

impl ResolvedLaminateBeamThermoelastic {
    /// Adopt the exact current beam snapshot and each original layer potential.
    /// Equal-valued independent thermal preparations are refused: the current
    /// source receipts and prepared mechanical allocation must be identical.
    pub fn try_new(
        laws: &[AffineLaminaThermoelastic],
        current: &ResolvedLaminateThermal,
        beam: &PreparedLaminateBeam,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LaminateBeamThermoelasticError> {
        check(&mut cancelled)?;
        let original = current.source().layers();
        if !core::ptr::eq(current.prepared().as_ref(), beam.source())
            || laws.len() != original.len()
            || current.layers().len() != laws.len()
            || beam.layers.len() != laws.len()
        {
            return Err(invalid(
                "beam potential requires its exact thermal snapshot and original layers",
            ));
        }
        let mut layers = Vec::with_capacity(laws.len());
        for (layer, ((law, source), material)) in
            laws.iter().zip(original).zip(current.layers()).enumerate()
        {
            check(&mut cancelled)?;
            let input = source.input();
            let bound = law.source();
            if !Arc::ptr_eq(&bound.input().card, &input.card)
                || bound.input().thickness_m.to_bits() != input.thickness_m.to_bits()
                || bound.input().orientation_rad.to_bits() != input.orientation_rad.to_bits()
                || bound.input().elasticity != input.elasticity
                || bound.input().expansion != input.expansion
                || bound.reference_material() != source.reference_material()
                || bound.reference_expansion().reference()
                    != source.reference_expansion().reference()
            {
                return Err(invalid(
                    "beam potential must retain each original card, reference and layer geometry",
                ));
            }
            let (sin, cos) = det::sin_cos(input.orientation_rad);
            let aligned = sin.abs().min(cos.abs()) <= 16.0 * f64::EPSILON;
            let isotropic = matches!(&input.elasticity, LaminaElasticSource::Isotropic(_))
                && input.expansion == LaminaExpansionSource::Isotropic;
            if !aligned && !isotropic {
                return Err(invalid(
                    "off-axis beam thermal paths require explicitly isotropic elastic and expansion sources",
                ));
            }
            let source_error = |source| LaminateBeamThermoelasticError::Source { layer, source };
            let resolved = law.at_materials(material).map_err(source_error)?;
            let minimum_capacity = law
                .minimum_heat_capacity_j_k(beam.width_m * input.thickness_m)
                .map_err(source_error)?;
            let q = beam.source.samples[layer].q;
            let lateral_scale = [det::sqrt(q[4]), det::sqrt(q[8])];
            let lateral_ratio = (q[5] / lateral_scale[0]) / lateral_scale[1];
            let lateral_root = det::sqrt((1.0 - lateral_ratio) * (1.0 + lateral_ratio));
            if lateral_scale
                .into_iter()
                .chain([lateral_root])
                .any(|v| !v.is_finite() || v <= 0.0)
                || !lateral_ratio.is_finite()
            {
                return Err(invalid("beam lateral relaxation metric is unrepresentable"));
            }
            layers.push(Layer {
                resolved,
                sin,
                cos,
                limit: law.maximum_absolute_engineering_strain(),
                minimum_capacity,
                lateral_scale,
                lateral_ratio,
                lateral_root,
            });
        }
        check(&mut cancelled)?;
        Ok(Self {
            current: current.clone(),
            beam: beam.clone(),
            layers,
        })
    }

    /// Exact actual source snapshot, including every layer temperature and receipt.
    pub const fn source(&self) -> &ResolvedLaminateThermal {
        &self.current
    }
    /// Original mechanical closure, axis, width, fixed mass and prepared Arc.
    pub const fn beam(&self) -> &PreparedLaminateBeam {
        &self.beam
    }
    pub fn layer_count(&self) -> usize {
        self.layers.len()
    }
    pub fn layer_temperature_k(&self, layer: usize) -> Option<f64> {
        self.layers.get(layer).map(|p| p.resolved.temperature_k())
    }
    /// Conservative elastic capacity floor over the original affine-source
    /// interval and full total/elastic principal-strain envelope, per original
    /// beam length. The original fixed-strain bound remains valid because
    /// lateral relaxation adds a nonnegative capacity correction. Add the
    /// original caloric capacity before claiming positive total capacity.
    pub fn layer_minimum_heat_capacity_j_m_k(
        &self,
        layer: usize,
    ) -> Result<f64, LaminateBeamThermoelasticError> {
        Ok(self.layer(layer)?.minimum_capacity)
    }

    /// One original layer's contribution. The full affine strain envelope is
    /// checked at both original faces before the two exact Gauss evaluations.
    /// Temperature is held uniform within that declared source layer/cell.
    pub fn layer_state(
        &self,
        layer: usize,
        eps0: f64,
        kappa: f64,
    ) -> Result<LaminateBeamThermoelasticState, LaminateBeamThermoelasticError> {
        let prepared = self.layer(layer)?;
        let faces = self
            .beam
            .layer_strain_faces(layer, eps0, kappa)
            .map_err(LaminateBeamThermoelasticError::Beam)?;
        let free_principal = prepared.resolved.free_strain();
        for face in faces {
            for (total, free) in principal_strain(face, prepared.sin, prepared.cos)
                .into_iter()
                .zip(free_principal)
            {
                if !total.is_finite()
                    || total.abs() > prepared.limit
                    || !(total - free).is_finite()
                    || (total - free).abs() > prepared.limit
                {
                    return Err(invalid(
                        "whole beam layer leaves the original total/elastic strain envelope",
                    ));
                }
            }
        }
        let sample = &self.beam.source.samples[layer];
        let reduction = self.beam.layers[layer];
        let volume = self.beam.width_m * sample.weight_m;
        let mut response = LaminateBeamThermoelasticState::default();
        for p in 0..2 {
            let y = sample.z[p] - self.beam.axis_offset_m;
            let total_x = (-y).mul_add(kappa, eps0);
            let elastic_x = total_x - sample.free[p][0];
            let total = [
                total_x,
                (-reduction.relaxation[0]).mul_add(elastic_x, sample.free[p][1]),
                (-reduction.relaxation[1]).mul_add(elastic_x, sample.free[p][2]),
            ];
            let at = prepared
                .resolved
                .in_plane_state(principal_strain(total, prepared.sin, prepared.cos), volume)
                .map_err(|source| LaminateBeamThermoelasticError::Source { layer, source })?;
            let b = section_stress(
                at.stress_temperature_derivative_pa_k,
                prepared.sin,
                prepared.cos,
            );
            let axial_t = (-reduction.relaxation[0])
                .mul_add(b[1], (-reduction.relaxation[1]).mul_add(b[2], b[0]));
            let first = b[1] / prepared.lateral_scale[0];
            let second = (-prepared.lateral_ratio).mul_add(first, b[2] / prepared.lateral_scale[1])
                / prepared.lateral_root;
            let correction =
                prepared.resolved.temperature_k() * volume * (first * first + second * second);
            let stress = reduction.axial_pa * elastic_x;
            let force = volume * stress;
            let force_t = volume * axial_t;
            response.add(LaminateBeamThermoelasticState {
                axial_force_n: force,
                bending_moment_n_m: -y * force,
                free_energy_j_m: at.free_energy_j,
                entropy_j_m_k: at.entropy_j_k,
                internal_energy_j_m: at.internal_energy_j,
                heat_capacity_j_m_k: at.heat_capacity_j_k + correction,
                generalized_force_temperature_derivative: [force_t, -y * force_t],
            });
        }
        response.checked()
    }

    /// Sum all original layers, retaining per-layer temperature weights in U.
    /// For independent layer-temperature heat balances use `layer_state`.
    pub fn state(
        &self,
        eps0: f64,
        kappa: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LaminateBeamThermoelasticState, LaminateBeamThermoelasticError> {
        let mut total = LaminateBeamThermoelasticState::default();
        for layer in 0..self.layers.len() {
            check(&mut cancelled)?;
            total.add(self.layer_state(layer, eps0, kappa)?);
        }
        check(&mut cancelled)?;
        total.checked()
    }
    fn layer(&self, layer: usize) -> Result<&Layer, LaminateBeamThermoelasticError> {
        self.layers
            .get(layer)
            .ok_or(LaminateBeamThermoelasticError::Beam(
                LaminateThermalError::LayerIndex {
                    index: layer,
                    count: self.layers.len(),
                },
            ))
    }
}

fn principal_strain(e: [f64; 3], s: f64, c: f64) -> [f64; 3] {
    let (ss, cc, sc) = (s * s, c * c, s * c);
    [
        cc * e[0] + ss * e[1] + sc * e[2],
        ss * e[0] + cc * e[1] - sc * e[2],
        2.0 * sc * (e[1] - e[0]) + (cc - ss) * e[2],
    ]
}
fn section_stress(p: [f64; 5], s: f64, c: f64) -> [f64; 3] {
    let (ss, cc, sc) = (s * s, c * c, s * c);
    [
        cc * p[0] + ss * p[1] - 2.0 * sc * p[2],
        ss * p[0] + cc * p[1] + 2.0 * sc * p[2],
        sc * (p[0] - p[1]) + (cc - ss) * p[2],
    ]
}
fn invalid(message: &'static str) -> LaminateBeamThermoelasticError {
    LaminateBeamThermoelasticError::Invalid(message)
}
fn check(cancelled: &mut impl FnMut() -> bool) -> Result<(), LaminateBeamThermoelasticError> {
    if cancelled() {
        Err(LaminateBeamThermoelasticError::Cancelled)
    } else {
        Ok(())
    }
}
