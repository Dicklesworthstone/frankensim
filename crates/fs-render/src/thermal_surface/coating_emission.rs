//! Reciprocal far-field emission from a nonisothermal passive planar stack.
//!
//! For the outgoing vacuum mode, the native coherent owner computes positive
//! per-film absorption of the time-reversed incoming plane wave. Each film is
//! internally in local thermal equilibrium at its own uniform temperature.
//! Its emitted s/p radiance is `A_[s,p] * B_lambda(T) / 2`. The opaque substrate,
//! when present, contributes its own absorbed fraction at its own temperature.
//! This is the planar-layer integral of the local Kirchhoff law in Greffet,
//! Bouchon, Brucoli and Marquier, *Light Emission by Nonequilibrium Bodies:
//! Local Kirchhoff Law*, Phys. Rev. X 8, 021008 (2018), equations (5)--(6):
//! <https://doi.org/10.1103/PhysRevX.8.021008>.
//!
//! The optical fields remain coherent; statistically independent thermal
//! sources add radiances. A single averaged temperature generally gives the
//! wrong spectrum. Both direction and polarization must match the reciprocal
//! absorption problem. Isotropic reciprocal local constitutive laws, planar
//! parallel boundaries and the original finite-band optical assumptions apply.
//! This is an outgoing radiance observation, not internal radiative exchange,
//! near-field transport, a rough-surface angular average or a thermal debit.

use super::{BOLTZMANN_J_K, LIGHT_M_S, ThermalSurfaceError, planck_radiance_w_m2_sr_nm};
use crate::conductor::{
    CoatingLayer, CoherentAbsorbingCoating, CoherentCoatingStack, CoherentDielectricCoating,
    ConductorError, ConductorOptics, MAX_COHERENT_STACK_LAYERS,
};
use crate::dielectric::{DielectricCoatingStack, DielectricError, DielectricIor};
use crate::polarization::{Mueller, Stokes};
use crate::spectral::{LAMBDA_MAX, LAMBDA_MIN};
use fs_blake3::{ContentHash, DomainHasher};
use std::sync::Arc;

/// Owned finite-band emitter with one temperature per original physical film.
///
/// The selected outgoing half space is exactly vacuum (`eta = 1`), consistent
/// with the existing thermal renderer's exterior-air approximation. Radiance
/// is W/(m² sr nm) per vacuum wavelength, without an exterior-index-squared
/// factor. Emission into a material medium is not admitted by these constructors.
/// The other transmitting half space may retain its original real dispersive
/// index, but receives no invented temperature or thermal source.
///
/// Preparation copies at most 64 full optical laws once. Clones share immutable
/// storage; wavelength/angle queries allocate nothing and borrow the original
/// numerical owners. The source identity binds the caller's physical state and
/// geometry declaration without upgrading the optical data's authority.
#[derive(Clone, Debug)]
pub struct CoatingThermalEmission {
    films: Arc<[Film]>,
    exit: Exit,
    source_identity: ContentHash,
    identity: ContentHash,
    optical_identity: ContentHash,
}

#[derive(Clone, Copy, Debug)]
struct Film {
    optics: OwnedLayer,
    temperature_k: f64,
}

#[derive(Clone, Copy, Debug)]
enum OwnedLayer {
    Absorbing(CoherentAbsorbingCoating),
    Lossless(CoherentDielectricCoating),
}

impl OwnedLayer {
    fn borrowed(&self) -> CoatingLayer<'_> {
        match self {
            Self::Absorbing(film) => CoatingLayer::Absorbing(film),
            Self::Lossless(film) => CoatingLayer::Dielectric(film),
        }
    }
}

#[derive(Clone, Copy, Debug)]
enum Exit {
    Opaque {
        substrate: ConductorOptics,
        temperature_k: f64,
    },
    Transmitting {
        opposite_ior: DielectricIor,
        reverse: bool,
    },
}

impl CoatingThermalEmission {
    /// Prepare emission through films authored from vacuum to an opaque
    /// substrate. Every film and the entire emitting substrate half space has
    /// its own uniform temperature. As in [`CoherentCoatingStack`], this explicit
    /// stack replaces any coating attached to `substrate`; only the substrate's
    /// original bulk complex-index law enters the solve.
    pub fn try_opaque(
        layers: &[CoatingLayer<'_>],
        film_temperatures_k: &[f64],
        substrate: ConductorOptics,
        substrate_temperature_k: f64,
        source_identity: ContentHash,
    ) -> Result<Self, CoatingEmissionError> {
        validate_temperature(substrate_temperature_k)?;
        Self::prepare(
            layers,
            film_temperatures_k,
            Exit::Opaque {
                substrate,
                temperature_k: substrate_temperature_k,
            },
            source_identity,
        )
    }

    /// Prepare finite-film emission into a selected vacuum half space.
    ///
    /// `opposite_ior` is the real-index medium on the OTHER side of the stack.
    /// `reverse = false` observes the authored first-film side; `true` observes
    /// the authored last-film side. The selected observed side must physically
    /// be vacuum in either case. Switching this flag alone does not turn a
    /// glass-side observation into a vacuum observation. Film temperatures and
    /// returned contributions always stay in authored order.
    ///
    /// The reciprocal optical query uses `eta_i = 1` and the original opposite
    /// phase-index law at the requested wavelength. Transmitted radiation is
    /// never emissivity: an empty or entirely lossless stack emits exactly zero.
    pub fn try_transmitting(
        layers: &[CoatingLayer<'_>],
        film_temperatures_k: &[f64],
        opposite_ior: DielectricIor,
        reverse: bool,
        source_identity: ContentHash,
    ) -> Result<Self, CoatingEmissionError> {
        Self::prepare(
            layers,
            film_temperatures_k,
            Exit::Transmitting {
                opposite_ior,
                reverse,
            },
            source_identity,
        )
    }

    fn prepare(
        layers: &[CoatingLayer<'_>],
        temperatures: &[f64],
        exit: Exit,
        source_identity: ContentHash,
    ) -> Result<Self, CoatingEmissionError> {
        let stack = CoherentCoatingStack::try_new(layers)?;
        if layers.len() != temperatures.len() {
            return Err(CoatingEmissionError::LayerTemperatures);
        }
        if source_identity == ContentHash([0; 32]) {
            return Err(ThermalSurfaceError::SourceIdentity.into());
        }
        for &temperature in temperatures {
            validate_temperature(temperature)?;
        }
        // A temperature change with fixed optical laws does not invalidate an
        // angular/spectral absorption response. Keep the existing complete
        // emission identity's byte sequence unchanged.
        let mut optical_hash =
            DomainHasher::new("org.frankensim.render.coating-emission-optics.v1");
        optical_hash.update(source_identity.as_bytes());
        optical_hash.update(stack.content_identity().as_bytes());
        match exit {
            Exit::Opaque { substrate, .. } => {
                optical_hash.update(&[0]);
                optical_hash.update(substrate.content_identity().as_bytes());
            }
            Exit::Transmitting {
                opposite_ior,
                reverse,
            } => {
                optical_hash.update(&[1, u8::from(reverse)]);
                hash_ior(&mut optical_hash, opposite_ior);
            }
        }
        let mut hash = DomainHasher::new("org.frankensim.render.coating-thermal-emission.v1");
        hash.update(source_identity.as_bytes());
        hash.update(stack.content_identity().as_bytes());
        for &temperature in temperatures {
            hash.update(&temperature.to_bits().to_le_bytes());
        }
        match exit {
            Exit::Opaque {
                substrate,
                temperature_k,
            } => {
                hash.update(&[0]);
                hash.update(substrate.content_identity().as_bytes());
                hash.update(&temperature_k.to_bits().to_le_bytes());
            }
            Exit::Transmitting {
                opposite_ior,
                reverse,
            } => {
                hash.update(&[1, u8::from(reverse)]);
                hash_ior(&mut hash, opposite_ior);
            }
        }
        let films = layers
            .iter()
            .zip(temperatures)
            .map(|(layer, &temperature_k)| Film {
                optics: match layer {
                    CoatingLayer::Absorbing(film) => OwnedLayer::Absorbing(**film),
                    CoatingLayer::Dielectric(film) => OwnedLayer::Lossless(**film),
                },
                temperature_k,
            })
            .collect();
        Ok(Self {
            films,
            exit,
            source_identity,
            identity: hash.finalize(),
            optical_identity: optical_hash.finalize(),
        })
    }

    /// Evaluate the outgoing vacuum mode at `0 < outgoing_cosine <= 1`.
    ///
    /// The cosine is with the outward normal of the SELECTED observed side.
    /// Its reciprocal incoming mode has the same cosine and the opposite ray
    /// direction. No absolute value, side swap, area or cosine projection is
    /// applied. Exact tangent directions have no positive normal incident flux
    /// and are refused, consistently for opaque and transmitting stacks.
    /// Wavelength support is exactly the original closed 380--780 nm interval.
    pub fn evaluate(
        &self,
        wavelength_nm: f64,
        outgoing_cosine: f64,
    ) -> Result<CoatingEmissionSample, CoatingEmissionError> {
        if !wavelength_nm.is_finite() || !(LAMBDA_MIN..=LAMBDA_MAX).contains(&wavelength_nm) {
            return Err(ThermalSurfaceError::Wavelength.into());
        }
        if !outgoing_cosine.is_finite() || outgoing_cosine <= 0.0 || outgoing_cosine > 1.0 {
            return Err(CoatingEmissionError::EmissionDirection);
        }
        self.with_layers(|layers| self.evaluate_layers(layers, wavelength_nm, outgoing_cosine))
    }

    fn with_layers<T>(
        &self,
        evaluate: impl FnOnce(&[CoatingLayer<'_>]) -> Result<T, CoatingEmissionError>,
    ) -> Result<T, CoatingEmissionError> {
        // Store references only: full spectral laws stay in their shared owned
        // allocation. The empty stack needs no placeholder physical layer.
        if self.films.is_empty() {
            return evaluate(&[]);
        }
        let mut layers = [self.films[0].optics.borrowed(); MAX_COHERENT_STACK_LAYERS];
        for (slot, film) in layers.iter_mut().zip(self.films.iter()) {
            *slot = film.optics.borrowed();
        }
        evaluate(&layers[..self.films.len()])
    }

    /// Positive reciprocal work for preparing exchange; no Planck source is
    /// evaluated and no temperature is substituted into the original optics.
    pub(super) fn reciprocal_absorbed_fraction(
        &self,
        wavelength_nm: f64,
        cosine: f64,
    ) -> Result<[f64; 2], CoatingEmissionError> {
        self.with_layers(|layers| match self.exit {
            Exit::Opaque { substrate, .. } => {
                let work = CoherentCoatingStack::try_new(layers)?.layer_power_balance(
                    &substrate,
                    wavelength_nm,
                    1.0,
                    cosine,
                )?;
                let mut absorbed = [PositiveSum::default(); 2];
                for fraction in work
                    .layer_absorbed()
                    .iter()
                    .chain(core::iter::once(&work.substrate_absorbed))
                {
                    for pol in 0..2 {
                        absorbed[pol].add(fraction[pol]);
                    }
                }
                Ok(absorbed.map(PositiveSum::value))
            }
            Exit::Transmitting {
                opposite_ior,
                reverse,
            } => Ok(DielectricCoatingStack::try_mixed(layers)?
                .layer_power_balance(
                    wavelength_nm,
                    reverse,
                    cosine,
                    1.0,
                    opposite_ior.eval(wavelength_nm)?,
                )?
                .absorbed),
        })
    }

    fn evaluate_layers(
        &self,
        layers: &[CoatingLayer<'_>],
        wavelength_nm: f64,
        cosine: f64,
    ) -> Result<CoatingEmissionSample, CoatingEmissionError> {
        let mut sample = CoatingEmissionSample {
            stokes: Stokes::unpolarized(0.0),
            films: [[0.0; 2]; MAX_COHERENT_STACK_LAYERS],
            film_count: layers.len(),
            substrate: None,
            reflected: [0.0; 2],
            transmitted: [0.0; 2],
            absorbed: [0.0; 2],
        };
        match self.exit {
            Exit::Opaque {
                substrate,
                temperature_k,
            } => {
                let work = CoherentCoatingStack::try_new(layers)?.layer_power_balance(
                    &substrate,
                    wavelength_nm,
                    1.0,
                    cosine,
                )?;
                sample.reflected = work.total.reflected;
                sample.absorbed = work.total.absorbed;
                self.emit_films(&mut sample, work.layer_absorbed(), wavelength_nm)?;
                sample.substrate = Some(emission(
                    temperature_k,
                    wavelength_nm,
                    work.substrate_absorbed,
                )?);
            }
            Exit::Transmitting {
                opposite_ior,
                reverse,
            } => {
                let work = DielectricCoatingStack::try_mixed(layers)?.layer_power_balance(
                    wavelength_nm,
                    reverse,
                    cosine,
                    1.0,
                    opposite_ior.eval(wavelength_nm)?,
                )?;
                sample.reflected = polarized_powers(work.fresnel.reflection);
                sample.transmitted = polarized_powers(work.fresnel.transmission);
                sample.absorbed = work.absorbed;
                self.emit_films(&mut sample, work.layer_absorbed(), wavelength_nm)?;
            }
        }
        // Positive compensated sums retain smaller emitting layers without
        // borrowing a single averaged temperature or subtracting outgoing flux.
        let mut sums = [PositiveSum::default(); 2];
        for radiance in sample
            .film_radiance_w_m2_sr_nm()
            .iter()
            .chain(sample.substrate.iter())
        {
            for pol in 0..2 {
                sums[pol].add(radiance[pol]);
            }
        }
        let [s, p] = sums.map(PositiveSum::value);
        sample.stokes =
            Stokes::try_new(s + p, s - p, 0.0, 0.0).ok_or(ThermalSurfaceError::Unrepresentable)?;
        Ok(sample)
    }

    fn emit_films(
        &self,
        sample: &mut CoatingEmissionSample,
        absorbed: &[[f64; 2]],
        wavelength_nm: f64,
    ) -> Result<(), CoatingEmissionError> {
        for ((slot, film), &fraction) in
            sample.films.iter_mut().zip(self.films.iter()).zip(absorbed)
        {
            *slot = emission(film.temperature_k, wavelength_nm, fraction)?;
        }
        Ok(())
    }

    /// Total emitted radiance, `I = L_s + L_p`, in W/(m² sr nm).
    /// Use [`Self::evaluate`] when polarization or individual film contributions
    /// are required; this scalar projection does not propagate Q/U/V.
    pub fn radiance_w_m2_sr_nm(
        &self,
        wavelength_nm: f64,
        outgoing_cosine: f64,
    ) -> Result<f64, CoatingEmissionError> {
        Ok(self.evaluate(wavelength_nm, outgoing_cosine)?.stokes.i)
    }

    /// Number of original films, including explicit lossless/zero-thickness ones.
    #[must_use]
    pub fn film_count(&self) -> usize {
        self.films.len()
    }

    /// Uniform temperature of one physical film in authored order, K.
    #[must_use]
    pub fn film_temperature_k(&self, index: usize) -> Option<f64> {
        self.films.get(index).map(|film| film.temperature_k)
    }

    /// Uniform opaque half-space temperature; absent for transmitting stacks.
    #[must_use]
    pub fn substrate_temperature_k(&self) -> Option<f64> {
        match self.exit {
            Exit::Opaque { temperature_k, .. } => Some(temperature_k),
            Exit::Transmitting { .. } => None,
        }
    }

    /// Binding supplied by the caller, without an inferred measurement claim.
    #[must_use]
    pub const fn source_identity(&self) -> ContentHash {
        self.source_identity
    }

    /// Original optical laws, film order, temperatures, selected side and source.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Original optical sources, films, geometry and selected vacuum side,
    /// excluding temperatures. A native-state consumer may reuse an absorption
    /// response while this identity is unchanged; temperature-dependent optics
    /// must first be resolved at the actual accepted state.
    #[must_use]
    pub const fn optical_identity(&self) -> ContentHash {
        self.optical_identity
    }
}

/// One outgoing spectral mode, including the independently emitting layers.
/// All radiances are W/(m² sr nm); reciprocal power fractions are dimensionless.
#[derive(Clone, Debug, PartialEq)]
pub struct CoatingEmissionSample {
    stokes: Stokes,
    films: [[f64; 2]; MAX_COHERENT_STACK_LAYERS],
    film_count: usize,
    substrate: Option<[f64; 2]>,
    reflected: [f64; 2],
    transmitted: [f64; 2],
    absorbed: [f64; 2],
}

impl CoatingEmissionSample {
    /// Physical incoherent Stokes radiance in the outgoing interface frame:
    /// +Q is s (perpendicular to the plane containing the normal and ray),
    /// `I=L_s+L_p`, `Q=L_s-L_p`, and U=V=0. Rotate that transverse frame before
    /// combining with a different ray frame or downstream Mueller element.
    /// At normal incidence s/p are degenerate and Q is zero.
    #[must_use]
    pub const fn stokes(&self) -> Stokes {
        self.stokes
    }

    /// Per-film `[L_s, L_p]` in authored order, even under reverse observation.
    #[must_use]
    pub fn film_radiance_w_m2_sr_nm(&self) -> &[[f64; 2]] {
        &self.films[..self.film_count]
    }

    /// Opaque-substrate `[L_s, L_p]`; no source exists in a real transmitting exit.
    #[must_use]
    pub const fn substrate_radiance_w_m2_sr_nm(&self) -> Option<[f64; 2]> {
        self.substrate
    }

    /// Native reflected `[R_s, R_p]` of the reciprocal incident mode.
    #[must_use]
    pub const fn reciprocal_reflected_fraction(&self) -> [f64; 2] {
        self.reflected
    }

    /// Native transmitted `[T_s, T_p]` of the reciprocal incident mode.
    /// These are ray-power fractions, without a camera-radiance index factor.
    #[must_use]
    pub const fn reciprocal_transmitted_fraction(&self) -> [f64; 2] {
        self.transmitted
    }

    /// Positive native absorbed `[A_s, A_p]`, not the rounded residual `1-R-T`.
    /// Only an isothermal stack has total emission `A_[s,p] * B_lambda(T)/2`.
    #[must_use]
    pub const fn reciprocal_absorbed_fraction(&self) -> [f64; 2] {
        self.absorbed
    }
}

/// Refusal before an unsupported mode or state becomes emitted radiance.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CoatingEmissionError {
    /// Every original film needs exactly one uniform temperature.
    LayerTemperatures,
    /// Outgoing cosine must be finite and in `(0,1]` for the selected vacuum side.
    EmissionDirection,
    /// Native Planck source, wavelength support or source-identity admission.
    Thermal(ThermalSurfaceError),
    /// Original conductor/coherent-film optical refusal.
    Conductor(ConductorError),
    /// Original transmitting-stack or opposite-index refusal.
    Dielectric(DielectricError),
}

impl From<ThermalSurfaceError> for CoatingEmissionError {
    fn from(error: ThermalSurfaceError) -> Self {
        Self::Thermal(error)
    }
}
impl From<ConductorError> for CoatingEmissionError {
    fn from(error: ConductorError) -> Self {
        Self::Conductor(error)
    }
}
impl From<DielectricError> for CoatingEmissionError {
    fn from(error: DielectricError) -> Self {
        Self::Dielectric(error)
    }
}
impl core::fmt::Display for CoatingEmissionError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            formatter,
            "invalid reciprocal coating thermal emission: {self:?}"
        )
    }
}
impl core::error::Error for CoatingEmissionError {}

fn validate_temperature(temperature_k: f64) -> Result<(), ThermalSurfaceError> {
    if !temperature_k.is_finite() || temperature_k <= 0.0 {
        return Err(ThermalSurfaceError::Temperature);
    }
    // Same finite-band Rayleigh--Jeans envelope as the original spectral
    // thermal owner. No T^4 or unobserved infrared continuation is assumed.
    let lambda2 = LAMBDA_MIN * LAMBDA_MIN;
    let upper_per_kelvin = 2.0 * LIGHT_M_S * BOLTZMANN_J_K * 1.0e27 / (lambda2 * lambda2);
    if temperature_k > f64::MAX / upper_per_kelvin {
        return Err(ThermalSurfaceError::Unrepresentable);
    }
    planck_radiance_w_m2_sr_nm(temperature_k, LAMBDA_MIN)?;
    Ok(())
}

fn emission(
    temperature_k: f64,
    wavelength_nm: f64,
    absorbed: [f64; 2],
) -> Result<[f64; 2], ThermalSurfaceError> {
    if absorbed == [0.0; 2] {
        return Ok([0.0; 2]);
    }
    let blackbody = planck_radiance_w_m2_sr_nm(temperature_k, wavelength_nm)?;
    // B already includes both independent polarizations. Do not halve a tiny
    // absorptance before multiplying by a bright Planck source.
    let radiance = absorbed.map(|fraction| 0.5 * (fraction * blackbody));
    if radiance
        .iter()
        .all(|value| value.is_finite() && *value >= 0.0)
    {
        Ok(radiance)
    } else {
        Err(ThermalSurfaceError::Unrepresentable)
    }
}

fn polarized_powers(matrix: Mueller) -> [f64; 2] {
    [
        matrix.0[0][0] + matrix.0[0][1],
        matrix.0[0][0] - matrix.0[0][1],
    ]
}

#[derive(Clone, Copy, Default)]
struct PositiveSum {
    sum: f64,
    correction: f64,
}
impl PositiveSum {
    fn add(&mut self, value: f64) {
        let next = self.sum + value;
        self.correction += if self.sum >= value {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
    }
    fn value(self) -> f64 {
        self.sum + self.correction
    }
}

fn hash_ior(hash: &mut DomainHasher, ior: DielectricIor) {
    match ior {
        DielectricIor::Cauchy(law) => {
            hash.update(&[0]);
            for value in law.coefficients() {
                hash.update(&value.to_bits().to_le_bytes());
            }
        }
        DielectricIor::Sellmeier(law) => {
            hash.update(&[1]);
            let (strengths, resonances_um2) = law.coefficients();
            for value in strengths.into_iter().chain(resonances_um2) {
                hash.update(&value.to_bits().to_le_bytes());
            }
        }
        DielectricIor::Tabulated(law) => {
            hash.update(&[2]);
            for sample in law.samples() {
                hash.update(&sample.wavelength_nm().to_bits().to_le_bytes());
                hash.update(&sample.index().to_bits().to_le_bytes());
            }
        }
    }
}

#[cfg(test)]
mod tests;
