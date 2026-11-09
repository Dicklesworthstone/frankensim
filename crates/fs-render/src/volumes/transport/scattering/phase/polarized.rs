//! Full-Stokes dispatch using the original scalar spectral mixture weights.

use super::PhaseFunction;
use crate::polarization::dielectric::transport::{BeamError, physical_roundoff};
use crate::polarization::{Mueller, Stokes};
use crate::volumes::transport::VolumeTransportError;
use crate::volumes::transport::scattering::polarized::PolarizedScatteringError;

impl PhaseFunction {
    /// Admit a retained molecular or spherical full-Stokes response. Spectral
    /// endpoints and active independent populations must all carry their own
    /// physical matrices; intensity-only tables/HG/isotropic laws still refuse.
    /// The native table/mixture limits bound recursion to three levels.
    pub fn admit_polarized(&self) -> Result<(), PolarizedScatteringError> {
        self.validate()?;
        match self {
            Self::Rayleigh | Self::RayleighDepolarized { .. } | Self::SphericalTable(_) => Ok(()),
            Self::SpectralTable(table) => {
                for (_, phase) in table.knots() {
                    phase.admit_polarized()?;
                }
                Ok(())
            }
            Self::SpectralMixture(mixture) => {
                for component in mixture.components() {
                    if component.mass_weight > 0.0 && !component.scattering.is_zero() {
                        component.phase.admit_polarized()?;
                    }
                }
                Ok(())
            }
            _ => Err(PolarizedScatteringError::UnsupportedPhase),
        }
    }

    /// Differential Mueller response per steradian at the admitted wavelength.
    /// Independent species use the SAME scattering weights as their phase PDF.
    /// This matrix is for inspection/composition; `apply_polarized_at` retains
    /// weak positive channels that a rounded M00-M01 cannot reconstruct.
    pub fn mueller_at(
        &self,
        wavelength_nm: f64,
        cosine: f64,
    ) -> Result<Mueller, PolarizedScatteringError> {
        self.validate_wavelength(wavelength_nm)?;
        match self {
            Self::Rayleigh | Self::RayleighDepolarized { .. } => self.rayleigh_mueller(cosine),
            Self::SphericalTable(table) => Ok(table.evaluate(cosine)?.mueller()),
            Self::SpectralTable(table) => {
                let (a, b, t) = table.interval(wavelength_nm);
                if t == 0.0 {
                    return a.mueller_at(wavelength_nm, cosine);
                }
                if t == 1.0 {
                    return b.mueller_at(wavelength_nm, cosine);
                }
                let x = a.mueller_at(wavelength_nm, cosine)?;
                let y = b.mueller_at(wavelength_nm, cosine)?;
                Ok(Mueller(core::array::from_fn(|r| {
                    core::array::from_fn(|c| (1.0 - t) * x.0[r][c] + t * y.0[r][c])
                })))
            }
            Self::SpectralMixture(mixture) => {
                let weights = mixture.weights_at(wavelength_nm)?;
                if weights.active == 0 {
                    return Err(PolarizedScatteringError::UnsupportedPhase);
                }
                let mut matrix = [[0.0; 4]; 4];
                for (component, &weight) in mixture.components().iter().zip(&weights.probability) {
                    if weight == 0.0 {
                        continue;
                    }
                    let value = component.phase.mueller_at(wavelength_nm, cosine)?;
                    for r in 0..4 {
                        for c in 0..4 {
                            matrix[r][c] += weight * value.0[r][c];
                        }
                    }
                }
                Ok(Mueller(matrix))
            }
            _ => Err(PolarizedScatteringError::UnsupportedPhase),
        }
    }

    /// Apply the physical differential response in a common s/p scattering
    /// frame. Pure weak channels and exact extinction are retained explicitly.
    /// Spectral interpolation and particle mixtures sum Stokes responses, never
    /// complex amplitudes. No PDF, concentration or extinction is applied here.
    pub fn apply_polarized_at(
        &self,
        wavelength_nm: f64,
        cosine: f64,
        input: Stokes,
    ) -> Result<Stokes, PolarizedScatteringError> {
        self.validate_wavelength(wavelength_nm)?;
        if !cosine.is_finite() || cosine.abs() > 1.0 + 2e-10 {
            return Err(VolumeTransportError::Phase.into());
        }
        if Stokes::try_new(input.i, input.q, input.u, input.v).is_none() {
            return Err(BeamError::InvalidStokes.into());
        }
        match self {
            Self::Rayleigh | Self::RayleighDepolarized { .. } => {
                self.validate()?;
                let delta = match self {
                    Self::RayleighDepolarized { depolarization } => *depolarization,
                    _ => 0.0,
                };
                let c = cosine.clamp(-1.0, 1.0);
                let scale = 3.0 / (4.0 * core::f64::consts::PI * (2.0 + delta));
                // Half-sums keep valid near-maximum finite Stokes powers from
                // overflowing before the bounded differential response acts.
                let perpendicular = 0.5 * input.i + 0.5 * input.q;
                let parallel = (0.5 * input.i - 0.5 * input.q) * c * c;
                let out = Stokes {
                    i: scale * ((1.0 - delta) * (perpendicular + parallel) + delta * input.i),
                    q: scale * (1.0 - delta) * (perpendicular - parallel),
                    u: scale * (1.0 - delta) * c * input.u,
                    v: scale * (1.0 - 2.0 * delta) * c * input.v,
                };
                if out.i == 0.0
                    && input.i > 0.0
                    && (perpendicular != 0.0 || c != 0.0 || delta != 0.0)
                {
                    return Err(BeamError::NumericalRange.into());
                }
                // Near a dipole null, frame-rotation roundoff in I+Q is
                // magnified relative to the weak outgoing intensity. Bound
                // the coherency correction by the original differential
                // scale, preserving I and the explicit positive-channel
                // underflow guard above. Spectral/mixture endpoints must not
                // reject a physical ray before their responses are combined.
                admit(physical_roundoff(out, scale * input.i)?)
            }
            Self::SphericalTable(table) => Ok(table.evaluate(cosine)?.apply(input)?),
            Self::SpectralTable(table) => {
                let (a, b, t) = table.interval(wavelength_nm);
                if t == 0.0 {
                    return a.apply_polarized_at(wavelength_nm, cosine, input);
                }
                if t == 1.0 {
                    return b.apply_polarized_at(wavelength_nm, cosine, input);
                }
                weighted_sum([
                    (1.0 - t, a.apply_polarized_at(wavelength_nm, cosine, input)?),
                    (t, b.apply_polarized_at(wavelength_nm, cosine, input)?),
                ])
            }
            Self::SpectralMixture(mixture) => {
                let weights = mixture.weights_at(wavelength_nm)?;
                if weights.active == 0 {
                    return Err(PolarizedScatteringError::UnsupportedPhase);
                }
                let mut sum = [0.0; 4];
                let mut correction = [0.0; 4];
                for (component, &weight) in mixture.components().iter().zip(&weights.probability) {
                    if weight == 0.0 {
                        continue;
                    }
                    accumulate(
                        &mut sum,
                        &mut correction,
                        weight,
                        component
                            .phase
                            .apply_polarized_at(wavelength_nm, cosine, input)?,
                    )?;
                }
                finish(sum, correction)
            }
            _ => Err(PolarizedScatteringError::UnsupportedPhase),
        }
    }
}

fn weighted_sum<const N: usize>(
    values: [(f64, Stokes); N],
) -> Result<Stokes, PolarizedScatteringError> {
    let mut sum = [0.0; 4];
    let mut correction = [0.0; 4];
    for (weight, value) in values {
        accumulate(&mut sum, &mut correction, weight, value)?;
    }
    finish(sum, correction)
}
fn accumulate(
    sum: &mut [f64; 4],
    correction: &mut [f64; 4],
    weight: f64,
    value: Stokes,
) -> Result<(), PolarizedScatteringError> {
    let input = [value.i, value.q, value.u, value.v];
    for k in 0..4 {
        let v = weight * input[k];
        if !v.is_finite() || (k == 0 && input[k] > 0.0 && weight > 0.0 && v == 0.0) {
            return Err(BeamError::NumericalRange.into());
        }
        let next = sum[k] + v;
        correction[k] += if sum[k].abs() >= v.abs() {
            (sum[k] - next) + v
        } else {
            (v - next) + sum[k]
        };
        sum[k] = next;
    }
    Ok(())
}
fn finish(sum: [f64; 4], correction: [f64; 4]) -> Result<Stokes, PolarizedScatteringError> {
    let s = core::array::from_fn::<_, 4, _>(|k| sum[k] + correction[k]);
    admit(Stokes {
        i: s[0],
        q: s[1],
        u: s[2],
        v: s[3],
    })
}
fn admit(out: Stokes) -> Result<Stokes, PolarizedScatteringError> {
    if Stokes::try_new(out.i, out.q, out.u, out.v).is_none() {
        return Err(BeamError::NumericalRange.into());
    }
    Ok(out)
}
