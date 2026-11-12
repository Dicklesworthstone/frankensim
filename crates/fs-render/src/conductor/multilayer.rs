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
    ConductorError, ConductorOptics, ConductorPowerBalance, ConductorSource, MAX_COHERENT_COATING_THICKNESS_M,
    MAX_IOR_COMPONENT, canonical_zero,
};
use crate::dielectric::CauchyIor;
use crate::polarization::Mueller;
use fs_blake3::{ContentHash, DomainHasher};
use fs_math::det;

mod absorption;

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

    pub(crate) fn update_identity(&self, hasher: &mut DomainHasher) {
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
        self.response::<false>(geometry, false)
    }

    pub(crate) fn transmission(
        &self,
        wavelength_nm: f64,
        incident_eta: f64,
        exit_eta: f64,
        cosine: f64,
        reverse: bool,
        exit_normal: [f64; 2],
    ) -> Result<CoherentTransmission, ConductorError> {
        // Both exterior half spaces carry propagating real-index ray power.
        // The bounded original films may absorb; their complex propagation
        // must remain in the transmitted field and must never be normalized
        // into a lossless R+T partition.
        if !crate::dielectric::admitted_dielectric_index(exit_eta) {
            return Err(ConductorError::InvalidSample);
        }
        let mut geometry = Geometry::new(
            wavelength_nm,
            incident_eta,
            ComplexIor {
                eta: exit_eta,
                k: 0.0,
            },
            cosine,
        )?;
        geometry.exit_normal = Some(Complex::new(exit_normal[0], exit_normal[1]));
        Ok(self
            .response::<true>(geometry, reverse)?
            .transmission_result(self.unit_cell().iter().all(|layer| {
                layer.thickness_m == 0.0 || layer.extinction_k == 0.0
            })))
    }

    fn response<const TRANSMISSION: bool>(
        &self,
        geometry: Geometry,
        reverse: bool,
    ) -> Result<Reflection, ConductorError> {
        self.response_as::<TRANSMISSION, ReflectionOnly>(geometry, reverse)
    }

    pub(super) fn power_balance(
        &self,
        wavelength_nm: f64,
        incident_eta: f64,
        substrate: ComplexIor,
        cosine: f64,
    ) -> Result<ConductorPowerBalance, ConductorError> {
        self.response_as::<false, AbsorbedPower>(
            Geometry::new(wavelength_nm, incident_eta, substrate, cosine)?, false,
        )
    }

    pub(super) fn layer_power_balance(
        &self,
        wavelength_nm: f64,
        incident_eta: f64,
        substrate: ComplexIor,
        cosine: f64,
    ) -> Result<super::ConductorLayerPowerBalance, ConductorError> {
        let geometry = Geometry::new(wavelength_nm, incident_eta, substrate, cosine)?;
        let mut prepared = [absorption::Layer::IDENTITY; MAX_COHERENT_UNIT_CELL_LAYERS];
        for (slot, layer) in prepared.iter_mut().zip(self.unit_cell()) {
            *slot = absorption::Layer::new::<false>(
                ComplexIor {
                    eta: layer.ior.eval(wavelength_nm).map_err(|_| ConductorError::InvalidWavelength)?,
                    k: layer.extinction_k,
                },
                layer.thickness_m,
                geometry,
            )?;
        }
        absorption::partition(geometry, self.layer_count(), |index| {
            Ok(prepared[index % self.cell_len as usize])
        })
    }

    fn response_as<const TRANSMISSION: bool, Kind: ResponseKind>(
        &self,
        geometry: Geometry,
        reverse: bool,
    ) -> Result<Kind::Output, ConductorError> {
        let mut prepared = [LayerTransfer::IDENTITY; MAX_COHERENT_UNIT_CELL_LAYERS];
        for (slot, layer) in prepared.iter_mut().zip(self.unit_cell()) {
            let ior = ComplexIor {
                eta: layer.ior.eval(geometry.wavelength_nm).map_err(|_| ConductorError::InvalidWavelength)?,
                k: layer.extinction_k,
            };
            *slot = LayerTransfer::new::<TRANSMISSION>(ior, layer.thickness_m, geometry);
        }
        response_as::<TRANSMISSION, Kind>(geometry, self.layer_count(), |index| {
            let index = if reverse {
                self.layer_count() - 1 - index
            } else {
                index
            };
            Ok(prepared[index % self.cell_len as usize])
        })
    }
}

/// Internal phase-preserving result shared with the transmitting-dielectric
/// owner. Transmission amplitudes are normalized to normal Poynting flux,
/// including the p-wave magnetic/electric field conversion.
pub(crate) struct CoherentTransmission {
    pub(crate) reflected: [[f64; 2]; 2],
    pub(crate) transmitted: [[f64; 2]; 2],
    pub(crate) bare_equivalent: bool,
    pub(crate) lossless: bool,
    /// First-order rounding bound on |R + T - 1| of this solve, per
    /// polarization (`ClosureBound`). `None` where its derivation does not
    /// apply: an evanescent exit, a step that is not a propagating lossless
    /// layer with the exact conjugate matrix structure, or a state past the
    /// first-order gate u kappa <= 2^-10.
    pub(crate) closure_allowance: [Option<f64>; 2],
}

/// Native transmitted/reflected fields plus positive work in finite films.
/// The real-index exit half-space receives transmitted radiation, never heat.
pub(crate) struct CoherentLayerTransmission {
    pub(crate) response: CoherentTransmission,
    pub(crate) absorbed: [f64; 2],
    pub(crate) layers: [[f64; 2]; MAX_COHERENT_STACK_LAYERS],
}

pub(crate) fn lossless_stack_transmission(
    layers: &[CoherentDielectricCoating],
    wavelength_nm: f64,
    incident_eta: f64,
    exit_eta: f64,
    cosine: f64,
    reverse: bool,
    exit_normal: [f64; 2],
) -> Result<CoherentTransmission, ConductorError> {
    if layers.len() > MAX_COHERENT_STACK_LAYERS {
        return Err(ConductorError::InvalidCoatingStack);
    }
    for layer in layers {
        super::validate_supported_wavelength(wavelength_nm, layer.wavelength_bounds_nm())?;
    }
    let mut geometry = Geometry::new(
        wavelength_nm,
        incident_eta,
        ComplexIor {
            eta: exit_eta,
            k: 0.0,
        },
        cosine,
    )?;
    geometry.exit_normal = Some(Complex::new(exit_normal[0], exit_normal[1]));
    Ok(response::<true>(geometry, layers.len(), |index| {
        let index = if reverse {
            layers.len() - 1 - index
        } else {
            index
        };
        let layer = layers[index];
        let eta = layer
            .ior()
            .eval(wavelength_nm)
            .map_err(|_| ConductorError::InvalidWavelength)?;
        Ok(LayerTransfer::new::<true>(
            ComplexIor { eta, k: 0.0 },
            layer.thickness_m(),
            geometry,
        ))
    })?
    .transmission_result(true))
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

impl CoatingLayer<'_> {
    /// Closed support of this original physical layer's optical law, in nm.
    /// Even a zero-thickness layer retains its source's declared domain.
    #[must_use]
    pub const fn wavelength_bounds_nm(self) -> (f64, f64) {
        match self {
            Self::Dielectric(film) => film.wavelength_bounds_nm(),
            Self::Absorbing(film) => film.wavelength_bounds_nm(),
        }
    }
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

    fn validate_wavelength(self, wavelength_nm: f64) -> Result<(), ConductorError> {
        for layer in self.layers {
            super::validate_supported_wavelength(wavelength_nm, layer.wavelength_bounds_nm())?;
        }
        Ok(())
    }

    /// Original complex films between two real-index surrounding media. The
    /// same scaled field kernel owns reflection and transmitted Poynting flux.
    /// Film order reverses together with the caller's incident/exit media.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn transmission(
        self,
        wavelength_nm: f64,
        incident_eta: f64,
        exit_eta: f64,
        cosine: f64,
        reverse: bool,
        exit_normal: [f64; 2],
    ) -> Result<CoherentTransmission, ConductorError> {
        self.validate_wavelength(wavelength_nm)?;
        if !crate::dielectric::admitted_dielectric_index(exit_eta) {
            return Err(ConductorError::InvalidSample);
        }
        let mut geometry = Geometry::new(
            wavelength_nm, incident_eta, ComplexIor { eta: exit_eta, k: 0.0 }, cosine,
        )?;
        geometry.exit_normal = Some(Complex::new(exit_normal[0], exit_normal[1]));
        let mut lossless = true;
        let result = response::<true>(geometry, self.layers.len(), |index| {
            let index = if reverse { self.layers.len() - 1 - index } else { index };
            let (ior, thickness) = self.layer(index, wavelength_nm)?;
            lossless &= thickness == 0.0 || ior.k == 0.0;
            Ok(LayerTransfer::new::<true>(ior, thickness, geometry))
        })?;
        Ok(result.transmission_result(lossless))
    }

    /// Retain the original power-normalized transmission solve while resolving
    /// each film's work from its incident port waves, rather than an R+T deficit.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn transmission_with_absorption(
        self,
        wavelength_nm: f64,
        incident_eta: f64,
        exit_eta: f64,
        cosine: f64,
        reverse: bool,
        exit_normal: [f64; 2],
    ) -> Result<CoherentLayerTransmission, ConductorError> {
        self.validate_wavelength(wavelength_nm)?;
        if !crate::dielectric::admitted_dielectric_index(exit_eta) {
            return Err(ConductorError::InvalidSample);
        }
        let mut geometry = Geometry::new(
            wavelength_nm, incident_eta, ComplexIor { eta: exit_eta, k: 0.0 }, cosine,
        )?;
        geometry.exit_normal = Some(Complex::new(exit_normal[0], exit_normal[1]));
        let mut lossless = true;
        let (response, mut work) = absorption::transmitting_partition(
            geometry, self.layers.len(), |index| {
                let index = if reverse { self.layers.len() - 1 - index } else { index };
                let (ior, thickness) = self.layer(index, wavelength_nm)?;
                lossless &= thickness == 0.0 || ior.k == 0.0;
                absorption::Layer::new::<true>(ior, thickness, geometry)
            },
        )?;
        if reverse {
            work.layers[..self.layers.len()].reverse();
        }
        Ok(CoherentLayerTransmission {
            response: response.transmission_result(lossless),
            absorbed: work.absorbed,
            layers: work.layers,
        })
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
        self.validate_wavelength(wavelength_nm)?;
        let geometry = Geometry::new(
            wavelength_nm,
            incident_eta,
            substrate.complex_ior(wavelength_nm)?,
            cosine,
        )?;
        self.reflection(geometry)?.reflectance(geometry)
    }

    /// Reflected and absorbed fractions from the same original boundary-field
    /// solve. This explicit ordered stack replaces any attached substrate film.
    /// Absorption is aggregate coating/opaque-substrate power, not per-layer heat.
    pub fn power_balance(
        self,
        substrate: &ConductorOptics,
        wavelength_nm: f64,
        incident_eta: f64,
        cosine: f64,
    ) -> Result<ConductorPowerBalance, ConductorError> {
        self.validate_wavelength(wavelength_nm)?;
        let geometry = Geometry::new(
            wavelength_nm, incident_eta, substrate.complex_ior(wavelength_nm)?, cosine,
        )?;
        response_as::<false, AbsorbedPower>(geometry, self.layers.len(), |index| {
            let (ior, thickness) = self.layer(index, wavelength_nm)?;
            Ok(LayerTransfer::new::<false>(ior, thickness, geometry))
        })
    }

    /// Resolve positive absorption in each authored film and in the opaque
    /// substrate from the original coherent boundary fields. This explicit
    /// ordered stack replaces any coating attached to the supplied substrate.
    pub fn layer_power_balance(
        self,
        substrate: &ConductorOptics,
        wavelength_nm: f64,
        incident_eta: f64,
        cosine: f64,
    ) -> Result<super::ConductorLayerPowerBalance, ConductorError> {
        self.validate_wavelength(wavelength_nm)?;
        let geometry = Geometry::new(
            wavelength_nm,
            incident_eta,
            substrate.complex_ior(wavelength_nm)?,
            cosine,
        )?;
        absorption::partition(geometry, self.layers.len(), |index| {
            let (ior, thickness) = self.layer(index, wavelength_nm)?;
            absorption::Layer::new::<false>(ior, thickness, geometry)
        })
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
        self.validate_wavelength(wavelength_nm)?;
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
            let (ior, thickness) = self.layer(index, geometry.wavelength_nm)?;
            Ok(LayerTransfer::new::<false>(ior, thickness, geometry))
        })
    }

    fn layer(self, index: usize, wavelength_nm: f64) -> Result<(ComplexIor, f64), ConductorError> {
        Ok(match self.layers[index] {
            CoatingLayer::Dielectric(layer) => (
                ComplexIor {
                    eta: layer.ior().eval(wavelength_nm)
                        .map_err(|_| ConductorError::InvalidWavelength)?,
                    k: 0.0,
                },
                layer.thickness_m(),
            ),
            CoatingLayer::Absorbing(layer) => (
                layer.complex_ior(wavelength_nm)?, layer.thickness_m(),
            ),
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
    // Transmitting boundaries retain the native Snell/TIR decision: its
    // ratio-based rounding can differ from n_t²-k_parallel² at criticality.
    exit_normal: Option<Complex>,
}

impl Geometry {
    fn new(
        wavelength_nm: f64,
        incident_eta: f64,
        substrate: ComplexIor,
        cosine: f64,
    ) -> Result<Self, ConductorError> {
        if !wavelength_nm.is_finite() || wavelength_nm <= 0.0 {
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
            exit_normal: None,
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
    propagation_log_scale: f64,
    // Off-diagonal coefficient in the fixed incident-admittance wave basis.
    // It is factored by n_layer²-n_incident², retaining weak matched-film
    // reflection without subtracting almost equal total fields.
    wave_mixing: [Complex; 2],
}

impl LayerTransfer {
    const IDENTITY: Self = Self {
        cosine: Complex::new(1.0, 0.0),
        upper: [Complex::new(0.0, 0.0); 2],
        lower: [Complex::new(0.0, 0.0); 2],
        inactive: true,
        incident_match: true,
        substrate_match: true,
        propagation_log_scale: 0.0,
        wave_mixing: [Complex::new(0.0, 0.0); 2],
    };

    fn new<const TRANSMISSION: bool>(ior: ComplexIor, thickness_m: f64, geometry: Geometry) -> Self {
        // A zero-thickness layer has no physical interfaces. Removing it
        // exactly also avoids catastrophic cancellation between its interfaces.
        if thickness_m == 0.0 {
            return Self::IDENTITY;
        }
        let n_squared = Complex::new(ior.eta, ior.k).square();
        let q_squared = super::longitudinal_wave_squared(
            geometry.incident_eta, ior, geometry.cosine,
        );
        let phase_scale = core::f64::consts::TAU * (thickness_m * 1.0e9) / geometry.wavelength_nm;
        let (cosine, sine_over_q, propagation_log_scale) =
            super::absorbing_film::propagation_with_scale(q_squared, phase_scale);
        let minus_i = Complex::new(0.0, -1.0);
        // Opaque reflection retains its field/flux evaluation and does not
        // prepare coefficients used only by the transmitting wave basis.
        let wave_mixing = if TRANSMISSION && geometry.cosine > 0.0 {
            let incident_squared = geometry.incident_eta * geometry.incident_eta;
            let contrast = Complex::new(
                (ior.eta - geometry.incident_eta) * (ior.eta + geometry.incident_eta)
                    - ior.k * ior.k,
                2.0 * ior.eta * ior.k,
            );
            let p_contrast = contrast
                * ((n_squared + Complex::new(incident_squared, 0.0))
                    * (geometry.cosine * geometry.cosine / incident_squared)
                    - Complex::new(1.0, 0.0));
            [
                minus_i * sine_over_q * contrast
                    * (-0.5 / (geometry.incident_eta * geometry.cosine)),
                minus_i * sine_over_q * p_contrast / n_squared
                    * (0.5 * geometry.incident_eta / geometry.cosine),
            ]
        } else { [Complex::new(0.0, 0.0); 2] };
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
            propagation_log_scale,
            wave_mixing,
        }
    }
}

struct Reflection {
    amplitudes: [Complex; 2],
    bare_equivalent: bool,
    transmission: [Complex; 2],
    closure_allowance: [Option<f64>; 2],
}

/// The lossless closure allowance, accumulated along the transmitting solve
/// of one polarization.
///
/// It applies where every active step is a propagating lossless layer
/// (`propagation_log_scale == 0`) whose computed matrix keeps the exact
/// conjugate structure [[A, B], [conj B, conj A]], and the exit wave
/// propagates. Then J = |field|^2 - |flux|^2 is conserved in exact
/// arithmetic, and to first order in u = eps / 2
///   |R + T - 1| <= |t|^2 u [18 S(|A|^2 + |B|^2)
///       + 20 S(|A| + |B|) sqrt(kappa_{k-1} kappa_k) + 8 / |t_0|^2
///       + 2 (Lambda + N max_j |L_j|)] + 16 u,
/// with kappa = (|field|^2 + |flux|^2) / J of the stored state.
/// - 18: each computed matrix's defect |A|^2 - |B|^2 - 1. D and m each
///   carry at most 7u, and c and the sine at most 4u on terms <= 1.
/// - 20: one step's arithmetic. A complex product with its add is within
///   3.83u (|A| + |B|) |v| per component. Cauchy-Schwarz against the
///   computed output norm gives 2 sqrt 2 x 3.83 = 10.83u (|A| + |B|)
///   sqrt(kappa_{k-1} kappa_k) J. The rescale division adds at most
///   4u kappa_k <= 4u (|A| + |B|) sqrt(kappa_{k-1} kappa_k). That makes
///   14.83, kept at 20.
/// - 8 / |t_0|^2: the seed's cancellation in 1 - |r_0|^2.
/// - 2 (Lambda + N max |L_j|): the logarithmic scale. Lambda is the sum of
///   |ln s_k|, and L_j are its N partial sums.
/// - 16u: the final divisions and the four-square closure sum.
///
/// The bound is first order. It is derived only while
/// u (|A_k| + |B_k|)^2 kappa_{k-1} <= 2^-10 at every step, the a-priori bound
/// on kappa_k from the previous computed state, and u kappa_0 <= 2^-10 at
/// the seed. Inside that gate the dropped terms are O(u kappa) relative to
/// the kept ones, and the 20 kept against 14.83 absorbs a second-order
/// constant up to 356. Rounding in kappa itself is O(u kappa) relative,
/// also second order. Past the gate the stored J has no resolution left (a
/// deep reflector reaches u kappa ~ 1e14), the allowance is None and the
/// caller keeps the flat check. Inside the gate a stored J <= 0 can only
/// come from a defect, and the solve refuses. A seed with J <= 0 has lost J
/// to rounding (|r_0| rounds to 1 only for an admittance ratio beyond 1 / u)
/// and gives None.
#[derive(Clone, Copy)]
struct ClosureBound {
    derived: bool,
    nonpositive: bool,
    kappa: f64,
    entries: f64,
    pending: f64,
    states: f64,
    logs: f64,
    partial: f64,
    steps: f64,
}

impl ClosureBound {
    /// The first-order regime: u kappa <= 2^-10.
    const REGIME: f64 = 1.0 / 1024.0;

    fn resolved(kappa: f64) -> bool {
        0.5 * f64::EPSILON * kappa <= Self::REGIME
    }

    fn seed(field: Complex, flux: Complex, propagating: bool) -> Self {
        let mut bound = Self {
            derived: false,
            nonpositive: false,
            kappa: 0.0,
            entries: 0.0,
            pending: 0.0,
            states: 0.0,
            logs: 0.0,
            partial: 0.0,
            steps: 0.0,
        };
        let (field_squared, flux_squared) = (field.norm_squared(), flux.norm_squared());
        let conserved = field_squared - flux_squared;
        if propagating && conserved > 0.0 {
            bound.kappa = (field_squared + flux_squared) / conserved;
            bound.derived = Self::resolved(bound.kappa);
        }
        bound
    }

    fn condition(&mut self, field: Complex, flux: Complex) -> f64 {
        let (field_squared, flux_squared) = (field.norm_squared(), flux.norm_squared());
        let conserved = field_squared - flux_squared;
        if conserved.is_nan() || conserved <= 0.0 {
            self.nonpositive = true;
            return f64::INFINITY;
        }
        (field_squared + flux_squared) / conserved
    }

    fn matrix(
        &mut self,
        upper: Complex,
        lower: Complex,
        mixing: Complex,
        propagation_log_scale: f64,
    ) {
        self.derived &= propagation_log_scale == 0.0
            && lower.re == upper.re
            && lower.im == -upper.im
            && mixing.re == 0.0;
        if self.derived {
            let (upper_squared, mixing_squared) = (upper.norm_squared(), mixing.norm_squared());
            self.entries += upper_squared + mixing_squared;
            self.pending = det::sqrt(upper_squared) + det::sqrt(mixing_squared);
            self.derived = Self::resolved(self.pending * self.pending * self.kappa);
        }
    }

    fn state(
        &mut self,
        field: Complex,
        flux: Complex,
        boundary_scale: f64,
        logarithmic_scale: f64,
    ) {
        if !self.derived {
            return;
        }
        let kappa = self.condition(field, flux);
        self.states += self.pending * det::sqrt(self.kappa * kappa);
        self.kappa = kappa;
        self.logs += boundary_scale.abs();
        self.partial = self.partial.max(logarithmic_scale.abs());
        self.steps += 1.0;
    }

    fn finish(
        self,
        transmitted: f64,
        seed_transmitted: f64,
    ) -> Result<Option<f64>, ConductorError> {
        // A J <= 0 found inside the gate refuses even though the next step's
        // gate then fails on its infinite kappa.
        if self.nonpositive {
            return Err(ConductorError::NumericalFailure);
        }
        if !self.derived {
            return Ok(None);
        }
        let u = 0.5 * f64::EPSILON;
        let allowance = transmitted
            * u
            * (18.0 * self.entries
                + 20.0 * self.states
                + 8.0 / seed_transmitted
                + 2.0 * (self.logs + self.steps * self.partial))
            + 16.0 * u;
        if !allowance.is_finite() {
            return Err(ConductorError::NumericalFailure);
        }
        Ok(Some(allowance))
    }
}

// Static output selection leaves the existing reflected/transmitted paths
// without absorption arithmetic or extra response storage. New heat consumers
// retain normal flux from the SAME scaled field solve.
trait ResponseKind {
    type Output;
    type History;
    const BOUNDARY_SCALES: bool = false;
    fn history(layer_count: usize) -> Self::History;
    fn record(_: &mut Self::History, _: usize, _: [Complex; 2], _: [Complex; 2], _: [f64; 2]) {}
    fn finish(
        geometry: Geometry,
        reflected: Reflection,
        field: [Complex; 2],
        flux: [Complex; 2],
        history: Self::History,
    ) -> Result<Self::Output, ConductorError>;
}
struct ReflectionOnly;
impl ResponseKind for ReflectionOnly {
    type Output = Reflection;
    type History = ();
    fn history(_: usize) {}
    fn finish(
        _: Geometry,
        reflected: Reflection,
        _: [Complex; 2],
        _: [Complex; 2],
        _: (),
    ) -> Result<Self::Output, ConductorError> {
        Ok(reflected)
    }
}
struct AbsorbedPower;
impl ResponseKind for AbsorbedPower {
    type Output = ConductorPowerBalance;
    type History = ();
    fn history(_: usize) {}
    fn finish(
        geometry: Geometry,
        reflected: Reflection,
        field: [Complex; 2],
        flux: [Complex; 2],
        _: (),
    ) -> Result<Self::Output, ConductorError> {
        ConductorPowerBalance::from_boundary(
            geometry.cosine,
            [
                geometry.incident_eta * geometry.cosine,
                geometry.cosine / geometry.incident_eta,
            ],
            field,
            flux,
            reflected.amplitudes,
        )
    }
}

impl Reflection {
    fn transmission_result(self, lossless: bool) -> CoherentTransmission {
        CoherentTransmission {
            reflected: self.amplitudes.map(|value| [value.re, value.im]),
            transmitted: self.transmission.map(|value| [value.re, value.im]),
            bare_equivalent: self.bare_equivalent,
            lossless,
            closure_allowance: self.closure_allowance,
        }
    }

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
    layer: impl FnMut(usize) -> Result<LayerTransfer, ConductorError>,
) -> Result<Reflection, ConductorError> {
    response::<false>(geometry, layer_count, layer)
}

fn response<const TRANSMISSION: bool>(
    geometry: Geometry,
    layer_count: usize,
    layer: impl FnMut(usize) -> Result<LayerTransfer, ConductorError>,
) -> Result<Reflection, ConductorError> {
    response_as::<TRANSMISSION, ReflectionOnly>(geometry, layer_count, layer)
}

fn response_as<const TRANSMISSION: bool, Kind: ResponseKind>(
    geometry: Geometry,
    layer_count: usize,
    mut layer: impl FnMut(usize) -> Result<LayerTransfer, ConductorError>,
) -> Result<Kind::Output, ConductorError> {
    let mut history = Kind::history(layer_count);
    if geometry.cosine == 0.0 {
        return Kind::finish(
            geometry,
            Reflection {
                amplitudes: [Complex::new(-1.0, 0.0); 2],
                bare_equivalent: true,
                transmission: [Complex::new(0.0, 0.0); 2],
                closure_allowance: [None; 2],
            },
            [Complex::new(1.0, 0.0); 2],
            [Complex::new(0.0, 0.0); 2],
            history,
        );
    }
    let n_squared = Complex::new(geometry.substrate.eta, geometry.substrate.k).square();
    let q = geometry.exit_normal.unwrap_or_else(|| {
        super::thin_film::forward_sqrt(super::longitudinal_wave_squared(
            geometry.incident_eta, geometry.substrate, geometry.cosine,
        ))
    });
    let incident = [
        geometry.incident_eta * geometry.cosine,
        geometry.cosine / geometry.incident_eta,
    ];
    let exit = [q, q / n_squared];
    let mut initial_transmission = [Complex::new(0.0, 0.0); 2];
    let mut field = [Complex::new(1.0, 0.0); 2];
    let mut flux = exit;
    if TRANSMISSION {
        // Here field/flux hold incoming/reflected wave coordinates, both
        // relative to the ORIGINAL incident admittance, rather than F/G.
        // Factor the bare contrast too: adjacent real exterior indices must
        // not seed a spurious or vanished weak reflection.
        let ni = geometry.incident_eta;
        let nt = geometry.substrate.eta;
        let ni2 = ni * ni;
        let nt2 = nt * nt;
        let contrast = (nt - ni) * (nt + ni);
        let s_denominator = Complex::new(incident[0], 0.0) + q;
        let p_denominator = Complex::new(incident[0] * nt2, 0.0) + q * ni2;
        flux = [
            Complex::new(-contrast, 0.0) / s_denominator.square(),
            Complex::new(ni2 * contrast
                * (geometry.cosine * geometry.cosine * (nt2 + ni2) - ni2), 0.0)
                / p_denominator.square(),
        ];
        for pol in 0..2 {
            if exit[pol].re > 0.0 {
                initial_transmission[pol] = Complex::new(
                    2.0 * det::sqrt(incident[pol]) * det::sqrt(exit[pol].re), 0.0,
                ) / (Complex::new(incident[pol], 0.0) + exit[pol]);
            }
        }
    }
    let mut closure: [ClosureBound; 2] = core::array::from_fn(|pol| {
        ClosureBound::seed(field[pol], flux[pol], TRANSMISSION && exit[pol].re > 0.0)
    });
    let mut logarithmic_scale = [0.0; 2];
    let mut bare_equivalent = true;
    let mut inside_incident = false;
    Kind::record(&mut history, layer_count, field, flux, [0.0; 2]);
    for index in (0..layer_count).rev() {
        let transfer = layer(index)?;
        if transfer.inactive {
            Kind::record(&mut history, index, field, flux, [0.0; 2]);
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
        let mut boundary_scale = [0.0; 2];
        for pol in 0..2 {
            let (next_field, next_flux) = if TRANSMISSION {
                let diagonal = (transfer.upper[pol] * incident[pol]
                    + transfer.lower[pol] * (1.0 / incident[pol])) * 0.5;
                let mixing = transfer.wave_mixing[pol];
                let upper_entry = transfer.cosine + diagonal;
                let lower_entry = transfer.cosine - diagonal;
                closure[pol].matrix(
                    upper_entry, lower_entry, mixing, transfer.propagation_log_scale,
                );
                (upper_entry * field[pol] - mixing * flux[pol],
                    mixing * field[pol] + lower_entry * flux[pol])
            } else {
                (transfer.cosine * field[pol] + transfer.upper[pol] * flux[pol],
                    transfer.lower[pol] * field[pol] + transfer.cosine * flux[pol])
            };
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
            if TRANSMISSION || Kind::BOUNDARY_SCALES {
                // The reflected ratio forgets this common scale; transmitted
                // flux must retain it. Log scaling also covers thick evanescent
                // barriers without constructing a growing exponential.
                boundary_scale[pol] = det::ln(scale) + transfer.propagation_log_scale;
                if TRANSMISSION {
                    logarithmic_scale[pol] += boundary_scale[pol];
                }
            }
            if TRANSMISSION {
                closure[pol].state(
                    field[pol], flux[pol], boundary_scale[pol], logarithmic_scale[pol],
                );
            }
        }
        Kind::record(&mut history, index, field, flux, boundary_scale);
    }
    if bare_equivalent && !TRANSMISSION {
        return Kind::finish(
            geometry,
            Reflection {
                amplitudes: super::thin_film::polarization::bare_amplitudes(
                    geometry.incident_eta,
                    geometry.substrate,
                    geometry.cosine,
                ),
                bare_equivalent,
                transmission: [Complex::new(0.0, 0.0); 2],
                closure_allowance: [None; 2],
            },
            [Complex::new(1.0, 0.0); 2],
            exit,
            history,
        );
    }
    let mut amplitudes = core::array::from_fn(|pol| {
        if TRANSMISSION {
            flux[pol] / field[pol]
        } else {
            (field[pol] * incident[pol] - flux[pol]) / (field[pol] * incident[pol] + flux[pol])
        }
    });
    if geometry.cosine == 1.0 {
        amplitudes[1] = Complex::new(-amplitudes[0].re, -amplitudes[0].im);
    }
    let mut transmission = [Complex::new(0.0, 0.0); 2];
    let mut closure_allowance = [None; 2];
    if TRANSMISSION {
        for pol in 0..2 {
            transmission[pol] = initial_transmission[pol] * det::exp(-logarithmic_scale[pol])
                / field[pol];
            closure_allowance[pol] = closure[pol].finish(
                transmission[pol].norm_squared(),
                initial_transmission[pol].norm_squared(),
            )?;
        }
        if geometry.cosine == 1.0 {
            transmission[1] = transmission[0];
            closure_allowance[1] = closure_allowance[0];
        }
    }
    Kind::finish(
        geometry,
        Reflection {
            amplitudes,
            bare_equivalent,
            transmission,
            closure_allowance,
        },
        field,
        flux,
        history,
    )
}

#[cfg(test)]
mod tests;
