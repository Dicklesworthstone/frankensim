//! Coherent passive multilayers on the existing opaque conductor substrate.
//!
//! Backwards propagation of the two boundary fields composes all reflected
//! amplitudes without multiplying ill-conditioned full transfer matrices.
//! Each layer uses the single-film analytic critical-angle limit and passive
//! exponential scaling; a common field/flux normalization after every layer
//! prevents overflow without changing its complex reflection coefficient.
//! Conventions follow Byrnes, <https://arxiv.org/abs/1603.02720>:
//! exp(-i omega t), n = eta + i k, and incident-side-first layer order.

use super::{
    CoherentAbsorbingCoating, CoherentDielectricCoating, Complex, ComplexIor, ConductorDataStatus,
    ConductorError, ConductorOptics, ConductorSource, MAX_COHERENT_COATING_THICKNESS_M,
    MAX_IOR_COMPONENT, canonical_zero,
};
use crate::dielectric::CauchyIor;
use crate::polarization::Mueller;
use crate::spectral::{LAMBDA_MAX, LAMBDA_MIN};
use fs_blake3::{ContentHash, DomainHasher};

/// Maximum physical layers, after expanding all repeats. No amplitude-based
/// truncation is used: even weak returning waves can control a coherent null.
pub const MAX_COHERENT_STACK_LAYERS: usize = 64;
/// Compact scene materials store at most four separately authored layer laws.
pub const MAX_COHERENT_UNIT_CELL_LAYERS: usize = 4;

/// One layer in a compact scene coating: Cauchy-dispersive phase index and an
/// explicitly wavelength-independent extinction coefficient over 380--780 nm.
/// This declared model never substitutes for an arbitrary measured spectrum.
/// Full spectral layers are admitted by [`CoherentCoatingStack`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CauchyCoatingLayer {
    ior: CauchyIor,
    extinction_k: f64,
    thickness_m: f64,
}

impl CauchyCoatingLayer {
    /// Physical thickness is normal to the layer in metres, from zero to
    /// 100 micrometres. Extinction is dimensionless and in `[0, 100]`.
    pub fn try_new(
        ior: CauchyIor,
        extinction_k: f64,
        thickness_m: f64,
    ) -> Result<Self, ConductorError> {
        if !extinction_k.is_finite() || !(0.0..=MAX_IOR_COMPONENT).contains(&extinction_k) {
            return Err(ConductorError::InvalidSample);
        }
        if !thickness_m.is_finite()
            || !(0.0..=MAX_COHERENT_COATING_THICKNESS_M).contains(&thickness_m)
        {
            return Err(ConductorError::InvalidCoating);
        }
        Ok(Self {
            ior,
            extinction_k: canonical_zero(extinction_k),
            thickness_m: canonical_zero(thickness_m),
        })
    }

    /// Declared visible-band phase-index law.
    #[must_use]
    pub const fn ior(self) -> CauchyIor {
        self.ior
    }

    /// Declared wavelength-independent imaginary index.
    #[must_use]
    pub const fn extinction_k(self) -> f64 {
        self.extinction_k
    }

    /// Uniform normal thickness in metres.
    #[must_use]
    pub const fn thickness_m(self) -> f64 {
        self.thickness_m
    }
}

/// Scene-integrated coherent stack made by repeating one ordered unit cell.
/// Layers run from the incident medium towards the opaque substrate: `[H, L]`
/// repeated twice means `[H, L, H, L]`, never a reordered effective medium.
/// One repeat admits any ordering of up to four compact layer laws.
///
/// The source binds the entire authored coating, including the dispersion and
/// constant-k assumption. Boundaries are locally parallel and laterally uniform;
/// the caller declares sufficient coherence over the whole stack. No finite
/// bandwidth, thickness distribution, scattering or inter-path coherence is
/// inferred. Evaluation allocates nothing, prepares each cell layer once per
/// wavelength/angle, and performs at most 64 field compositions. Its storage
/// fits inside the existing single absorbing-film material payload.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentPeriodicCoating {
    cell: [CauchyCoatingLayer; MAX_COHERENT_UNIT_CELL_LAYERS],
    cell_len: u8,
    repetitions: u8,
    source: ConductorSource,
}

impl CoherentPeriodicCoating {
    /// Admit one to four explicitly ordered layer laws and a positive repeat
    /// count. The expanded stack may contain at most 64 physical layers.
    pub fn try_new(
        cell: &[CauchyCoatingLayer],
        repetitions: usize,
        source: ConductorSource,
    ) -> Result<Self, ConductorError> {
        if cell.is_empty()
            || cell.len() > MAX_COHERENT_UNIT_CELL_LAYERS
            || repetitions == 0
            || repetitions > MAX_COHERENT_STACK_LAYERS / cell.len()
        {
            return Err(ConductorError::InvalidCoatingStack);
        }
        let mut storage = [cell[0]; MAX_COHERENT_UNIT_CELL_LAYERS];
        storage[..cell.len()].copy_from_slice(cell);
        Ok(Self {
            cell: storage,
            cell_len: cell.len() as u8,
            repetitions: repetitions as u8,
            source,
        })
    }

    /// Canonical unit cell, incident side first.
    #[must_use]
    pub fn unit_cell(&self) -> &[CauchyCoatingLayer] {
        &self.cell[..self.cell_len as usize]
    }

    /// Positive number of complete ordered repeats.
    #[must_use]
    pub const fn repetitions(self) -> usize {
        self.repetitions as usize
    }

    /// Total physical layers, bounded by [`MAX_COHERENT_STACK_LAYERS`].
    #[must_use]
    pub const fn layer_count(self) -> usize {
        self.cell_len as usize * self.repetitions as usize
    }

    /// Caller-declared source and status for the complete authored coating.
    #[must_use]
    pub const fn source(self) -> ConductorSource {
        self.source
    }

    pub(super) fn update_identity(&self, hasher: &mut DomainHasher) {
        hasher.update(b"coherent-periodic-coating.v1");
        hasher.update(&[self.cell_len, self.repetitions]);
        for layer in self.unit_cell() {
            for coefficient in layer.ior.coefficients() {
                hasher.update(&coefficient.to_bits().to_le_bytes());
            }
            hasher.update(&layer.extinction_k.to_bits().to_le_bytes());
            hasher.update(&layer.thickness_m.to_bits().to_le_bytes());
        }
        hasher.update(self.source.identity().as_bytes());
        hasher.update(&[match self.source.status() {
            ConductorDataStatus::Representative => 0,
            ConductorDataStatus::CallerAssertedMeasured => 1,
            ConductorDataStatus::MaterialStateClaims => 2,
        }]);
    }

    pub(super) fn fresnel(
        &self,
        wavelength_nm: f64,
        incident_eta: f64,
        substrate: ComplexIor,
        cosine: f64,
    ) -> Result<f64, ConductorError> {
        let geometry = Geometry::new(wavelength_nm, incident_eta, substrate, cosine)?;
        self.reflection(geometry)?.reflectance(geometry)
    }

    pub(super) fn amplitudes(
        &self,
        wavelength_nm: f64,
        incident_eta: f64,
        substrate: ComplexIor,
        cosine: f64,
    ) -> Result<[Complex; 2], ConductorError> {
        Ok(self
            .reflection(Geometry::new(
                wavelength_nm,
                incident_eta,
                substrate,
                cosine,
            )?)?
            .amplitudes)
    }

    fn reflection(&self, geometry: Geometry) -> Result<Reflection, ConductorError> {
        let mut prepared = [LayerTransfer::IDENTITY; MAX_COHERENT_UNIT_CELL_LAYERS];
        for (slot, layer) in prepared.iter_mut().zip(self.unit_cell()) {
            let ior = ComplexIor {
                eta: layer.ior.eval_unchecked(geometry.wavelength_nm),
                k: layer.extinction_k,
            };
            *slot = LayerTransfer::new(ior, layer.thickness_m, geometry);
        }
        reflection(geometry, self.layer_count(), |index| {
            Ok(prepared[index % self.cell_len as usize])
        })
    }
}

/// Borrow one fully declared transparent or absorbing spectral layer. Its
/// original optical law, physical thickness and source identity are retained.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CoatingLayer<'a> {
    /// Shared Cauchy, Sellmeier or explicit tabulated real phase index.
    Dielectric(&'a CoherentDielectricCoating),
    /// Nine-knot eta and k, including exact zero extinction.
    Absorbing(&'a CoherentAbsorbingCoating),
}

/// Borrowed general stack with arbitrary layer ordering and full spectral
/// laws, evaluated with the same kernel as compact scene coatings. This API
/// does not attach heap-backed stacks to Copy scene materials. It is useful
/// for spectral design/evaluation without approximating measured eta or k.
/// Layers are incident-side first; an empty stack is exactly bare substrate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentCoatingStack<'a> {
    layers: &'a [CoatingLayer<'a>],
}

impl<'a> CoherentCoatingStack<'a> {
    /// Admit up to 64 fully authored layers. Constructors of the borrowed
    /// films have already validated spectral support and each thickness.
    pub fn try_new(layers: &'a [CoatingLayer<'a>]) -> Result<Self, ConductorError> {
        if layers.len() > MAX_COHERENT_STACK_LAYERS {
            return Err(ConductorError::InvalidCoatingStack);
        }
        Ok(Self { layers })
    }

    /// Borrowed layers in physical incident-to-substrate order.
    #[must_use]
    pub const fn layers(self) -> &'a [CoatingLayer<'a>] {
        self.layers
    }

    /// Unpolarized reflected power. Only the supplied substrate's complex-IOR
    /// table is used: this explicit stack replaces any coating attached to it.
    /// The opaque substrate makes total absorption `1-R`; this does not resolve
    /// absorption or heat deposition separately in individual layers.
    pub fn fresnel(
        self,
        substrate: &ConductorOptics,
        wavelength_nm: f64,
        incident_eta: f64,
        cosine: f64,
    ) -> Result<f64, ConductorError> {
        let geometry = Geometry::new(
            wavelength_nm,
            incident_eta,
            substrate.complex_ior(wavelength_nm)?,
            cosine,
        )?;
        self.reflection(geometry)?.reflectance(geometry)
    }

    /// Reflected Mueller matrix with the same s/p frames, phase convention
    /// and opaque-substrate semantics as [`ConductorOptics::fresnel_mueller`].
    /// The supplied explicit stack replaces any coating on the substrate.
    pub fn fresnel_mueller(
        self,
        substrate: &ConductorOptics,
        wavelength_nm: f64,
        incident_eta: f64,
        cosine: f64,
    ) -> Result<Mueller, ConductorError> {
        let geometry = Geometry::new(
            wavelength_nm,
            incident_eta,
            substrate.complex_ior(wavelength_nm)?,
            cosine,
        )?;
        super::thin_film::polarization::mueller_from_amplitudes(
            self.reflection(geometry)?.amplitudes,
            cosine,
        )
    }

    /// Identity of ordered optical laws, physical thicknesses and their
    /// original sources. The substrate retains its separate content identity.
    #[must_use]
    pub fn content_identity(self) -> ContentHash {
        let mut hasher = DomainHasher::new("org.frankensim.coherent-coating-stack.v1");
        hasher.update(&(self.layers.len() as u64).to_le_bytes());
        for layer in self.layers {
            match layer {
                CoatingLayer::Dielectric(layer) => layer.update_identity(&mut hasher),
                CoatingLayer::Absorbing(layer) => layer.update_identity(&mut hasher),
            }
        }
        hasher.finalize()
    }

    fn reflection(self, geometry: Geometry) -> Result<Reflection, ConductorError> {
        reflection(geometry, self.layers.len(), |index| {
            let (ior, thickness) = match self.layers[index] {
                CoatingLayer::Dielectric(layer) => (
                    ComplexIor {
                        eta: layer
                            .ior()
                            .eval(geometry.wavelength_nm)
                            .map_err(|_| ConductorError::NumericalFailure)?,
                        k: 0.0,
                    },
                    layer.thickness_m(),
                ),
                CoatingLayer::Absorbing(layer) => (
                    layer.complex_ior(geometry.wavelength_nm)?,
                    layer.thickness_m(),
                ),
            };
            Ok(LayerTransfer::new(ior, thickness, geometry))
        })
    }
}

#[derive(Clone, Copy)]
struct Geometry {
    wavelength_nm: f64,
    incident_eta: f64,
    substrate: ComplexIor,
    cosine: f64,
    transverse_squared: f64,
}

impl Geometry {
    fn new(
        wavelength_nm: f64,
        incident_eta: f64,
        substrate: ComplexIor,
        cosine: f64,
    ) -> Result<Self, ConductorError> {
        if !wavelength_nm.is_finite() || !(LAMBDA_MIN..=LAMBDA_MAX).contains(&wavelength_nm) {
            return Err(ConductorError::InvalidWavelength);
        }
        if !incident_eta.is_finite() || incident_eta <= 0.0 || incident_eta > MAX_IOR_COMPONENT {
            return Err(ConductorError::InvalidIncidentMedium);
        }
        if !cosine.is_finite() || !(0.0..=1.0).contains(&cosine) {
            return Err(ConductorError::InvalidIncidence);
        }
        Ok(Self {
            wavelength_nm,
            incident_eta,
            substrate,
            cosine,
            transverse_squared: incident_eta * incident_eta * (1.0 - cosine * cosine),
        })
    }
}

#[derive(Clone, Copy)]
struct LayerTransfer {
    cosine: Complex,
    upper: [Complex; 2],
    lower: [Complex; 2],
    inactive: bool,
    incident_match: bool,
    substrate_match: bool,
}

impl LayerTransfer {
    const IDENTITY: Self = Self {
        cosine: Complex::new(1.0, 0.0),
        upper: [Complex::new(0.0, 0.0); 2],
        lower: [Complex::new(0.0, 0.0); 2],
        inactive: true,
        incident_match: true,
        substrate_match: true,
    };

    fn new(ior: ComplexIor, thickness_m: f64, geometry: Geometry) -> Self {
        // A zero-thickness layer has no physical interfaces. Removing it
        // exactly also avoids catastrophic cancellation between its interfaces.
        if thickness_m == 0.0 {
            return Self::IDENTITY;
        }
        let n_squared = Complex::new(ior.eta, ior.k).square();
        let q_squared = n_squared - Complex::new(geometry.transverse_squared, 0.0);
        let phase_scale = core::f64::consts::TAU * (thickness_m * 1.0e9) / geometry.wavelength_nm;
        let (cosine, sine_over_q) = super::absorbing_film::propagation(q_squared, phase_scale);
        let minus_i = Complex::new(0.0, -1.0);
        Self {
            cosine,
            upper: [minus_i * sine_over_q, minus_i * n_squared * sine_over_q],
            lower: [
                minus_i * q_squared * sine_over_q,
                minus_i * q_squared * sine_over_q / n_squared,
            ],
            inactive: false,
            incident_match: ior.k == 0.0 && ior.eta == geometry.incident_eta,
            substrate_match: ior == geometry.substrate,
        }
    }
}

struct Reflection {
    amplitudes: [Complex; 2],
    bare_equivalent: bool,
}

impl Reflection {
    fn reflectance(self, geometry: Geometry) -> Result<f64, ConductorError> {
        if self.bare_equivalent {
            return super::fresnel_conductor(
                geometry.incident_eta,
                geometry.substrate.eta,
                geometry.substrate.k,
                geometry.cosine,
            );
        }
        let value = 0.5 * (self.amplitudes[0].norm_squared() + self.amplitudes[1].norm_squared());
        if !value.is_finite() || !(-1.0e-12..=1.0 + 1.0e-12).contains(&value) {
            return Err(ConductorError::NumericalFailure);
        }
        Ok(value.clamp(0.0, 1.0))
    }
}

fn reflection(
    geometry: Geometry,
    layer_count: usize,
    mut layer: impl FnMut(usize) -> Result<LayerTransfer, ConductorError>,
) -> Result<Reflection, ConductorError> {
    if geometry.cosine == 0.0 {
        return Ok(Reflection {
            amplitudes: [Complex::new(-1.0, 0.0); 2],
            bare_equivalent: true,
        });
    }
    let n_squared = Complex::new(geometry.substrate.eta, geometry.substrate.k).square();
    let q =
        super::thin_film::forward_sqrt(n_squared - Complex::new(geometry.transverse_squared, 0.0));
    let mut field = [Complex::new(1.0, 0.0); 2];
    let mut flux = [q, q / n_squared];
    let mut bare_equivalent = true;
    let mut inside_incident = false;
    for index in (0..layer_count).rev() {
        let transfer = layer(index)?;
        if transfer.inactive {
            continue;
        }
        // Incident-matched prefixes and substrate-matched suffixes alter only
        // a common phase or the arbitrary substrate reference plane. Preserve
        // the established scalar and Mueller bits when these are all layers.
        if transfer.incident_match {
            inside_incident = true;
        } else if !transfer.substrate_match || inside_incident {
            bare_equivalent = false;
        }
        for pol in 0..2 {
            let next_field = transfer.cosine * field[pol] + transfer.upper[pol] * flux[pol];
            let next_flux = transfer.lower[pol] * field[pol] + transfer.cosine * flux[pol];
            let scale = next_field
                .re
                .abs()
                .max(next_field.im.abs())
                .max(next_flux.re.abs())
                .max(next_flux.im.abs());
            if !scale.is_finite() || scale == 0.0 {
                return Err(ConductorError::NumericalFailure);
            }
            // Divide components separately: a tiny nonzero scale must not be
            // inverted into infinity before normalizing the projective fields.
            field[pol] = Complex::new(next_field.re / scale, next_field.im / scale);
            flux[pol] = Complex::new(next_flux.re / scale, next_flux.im / scale);
        }
    }
    if bare_equivalent {
        return Ok(Reflection {
            amplitudes: super::thin_film::polarization::bare_amplitudes(
                geometry.incident_eta,
                geometry.substrate,
                geometry.cosine,
            ),
            bare_equivalent,
        });
    }
    let incident = [
        geometry.incident_eta * geometry.cosine,
        geometry.cosine / geometry.incident_eta,
    ];
    let mut amplitudes = core::array::from_fn(|pol| {
        (field[pol] * incident[pol] - flux[pol]) / (field[pol] * incident[pol] + flux[pol])
    });
    if geometry.cosine == 1.0 {
        amplitudes[1] = Complex::new(-amplitudes[0].re, -amplitudes[0].im);
    }
    Ok(Reflection {
        amplitudes,
        bare_equivalent,
    })
}

#[cfg(test)]
mod tests;
