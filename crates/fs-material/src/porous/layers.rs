//! Coherent rigid-backed porous/fluid multilayers and complex reflection.
//!
//! Ordered homogeneous isotropic equivalent-fluid layers share ONE saturating
//! fluid snapshot. Pore drag, density, bulk modulus and source-frequency support
//! remain owned by PorousFluidLaw. Gaps use the same fluid's lossless acoustics;
//! no elastic frame, lateral finite-size effect or gap wall loss is invented.
//!
//! exp(+i omega t), forward exp(-i k_n z). Tangential wavenumber is conserved,
//! k_n^2=k^2-k_0^2 sin^2(theta), Z_n=Z_c k/k_n. Pressure and normal velocity
//! are continuous between layers, terminating in v_n=0 at the rigid backing.
//! A bounded attenuating-wave transfer replaces growing cos/sin matrices:
//! E=exp(-2i k_n d); [p; Z_ref v] propagates through
//! [[1+E, (Z_n/Z_ref)(1-E)],[(Z_ref/Z_n)(1-E),1+E]], up to a common scalar.
//! Each step is rescaled, retaining impedance poles without dividing by p or v.
//!
//! Reference: COMSOL Porous Absorber, Eq. (5), rigid-frame oblique response:
//! <https://doc.comsol.com/6.1/doc/com.comsol.help.models.aco.porous_absorber/models.aco.porous_absorber.pdf>.
//! Plane reflection here is extended reaction. A local BEM boundary must
//! explicitly select the NORMAL-incidence admittance instead.

use super::{PorousFluidLaw, PorousFluidResponse};
use crate::MaterialError;
use crate::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_math::{c64::C64, det};

/// Bounded preparation storage and per-angle work.
pub const MAX_LAYERS: usize = 64;

/// Ordered from the exposed face toward the rigid backing. Zero thickness is
/// an exact identity, but an included porous source still must admit frequency.
#[derive(Clone, Copy, Debug)]
pub enum Layer<'a> {
    Porous { law: &'a PorousFluidLaw, thickness_m: f64 },
    /// Lossless gap filled by the original coherent saturating/exterior fluid.
    FluidGap { thickness_m: f64 },
}
impl Layer<'_> {
    fn thickness_m(self) -> f64 {
        match self {
            Self::Porous { thickness_m, .. } | Self::FluidGap { thickness_m } => thickness_m,
        }
    }
}

#[derive(Clone, Copy)]
struct PreparedLayer<'a> {
    source: Layer<'a>,
    porous: Option<PorousFluidResponse>,
}

/// Prepared material responses at ONE frequency and ONE coherent fluid state.
/// The fluid source and all porous source laws remain borrowed, not replaced by
/// equivalent fitted coefficients. Each angle query is allocation-free.
pub struct RigidBackedLayers<'a> {
    source: &'a dyn ThermoviscousFluid,
    fluid: ThermoviscousFluidState,
    omega: f64,
    k0: f64,
    z0: f64,
    layers: Vec<PreparedLayer<'a>>,
}

/// Complex pressure reflection and normal absorbed/incident power fraction.
/// Values are not clipped; zero absorption may have floating-point residue
/// within the declared arithmetic balance check.
#[derive(Clone, Copy, Debug)]
pub struct LayeredReflection {
    pub pressure_reflection: C64,
    pub absorbed_fraction: f64,
    pressure: C64,
    normalized_velocity: C64,
    cosine: f64,
}
impl LayeredReflection {
    /// beta=Z0*v_n/p in exp(+i omega t). For a BEM owner with exp(-i omega t),
    /// conjugate ONCE. A pressure node has no representable finite admittance;
    /// reflection itself remains valid there. There is no panel-area factor.
    pub fn normalized_surface_admittance(&self) -> Result<C64, MaterialError> {
        if self.pressure == C64::ZERO {
            return Err(bad("layered boundary has a pressure node, not a finite admittance"));
        }
        let beta = self.normalized_velocity.scale(self.cosine) / self.pressure;
        if !finite(beta) {
            return Err(bad("layered surface admittance is unrepresentable"));
        }
        Ok(beta)
    }
}

fn bad(message: &'static str) -> MaterialError {
    MaterialError::State { what: message.into() }
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn magnitude_scale(z: C64) -> f64 {
    z.re.abs().max(z.im.abs())
}
fn divided(z: C64, scale: f64) -> C64 {
    C64::new(z.re / scale, z.im / scale)
}

impl<'a> RigidBackedLayers<'a> {
    pub fn new(
        fluid: &'a dyn ThermoviscousFluid,
        omega: f64,
        layers: &[Layer<'a>],
    ) -> Result<Self, MaterialError> {
        if layers.len() > MAX_LAYERS || !omega.is_finite() || omega <= 0.0
            || layers.iter().any(|l| !l.thickness_m().is_finite() || l.thickness_m() < 0.0)
        {
            return Err(bad("layered reflection needs <=64 layers, positive frequency and finite nonnegative thickness"));
        }
        // The actual source is evaluated exactly once. Re-querying it per
        // layer could accidentally compose different evolving fluid states.
        let snapshot = fluid.thermoviscous_state()?;
        let acoustic = snapshot.acoustic();
        let k0 = omega / acoustic.sound_speed_m_s();
        let z0 = acoustic.characteristic_impedance_pa_s_m();
        if !k0.is_finite() || k0 <= 0.0 || !z0.is_finite() || z0 <= 0.0 {
            return Err(bad("layered exterior wave data are unrepresentable"));
        }
        let mut prepared = Vec::with_capacity(layers.len());
        for &source in layers {
            let porous = match source {
                Layer::Porous { law, .. } => Some(law.effective_fluid(&snapshot, omega)?),
                Layer::FluidGap { .. } => None,
            };
            prepared.push(PreparedLayer { source, porous });
        }
        Ok(Self { source: fluid, fluid: snapshot, omega, k0, z0, layers: prepared })
    }

    pub const fn fluid_source(&self) -> &'a dyn ThermoviscousFluid {
        self.source
    }
    pub const fn fluid(&self) -> &ThermoviscousFluidState {
        &self.fluid
    }
    pub const fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    pub const fn exterior_wavenumber_rad_m(&self) -> f64 {
        self.k0
    }
    pub fn layer_count(&self) -> usize {
        self.layers.len()
    }
    pub fn layer(&self, index: usize) -> Option<Layer<'a>> {
        self.layers.get(index).map(|l| l.source)
    }
    pub fn material_response(&self, index: usize) -> Option<&PorousFluidResponse> {
        self.layers.get(index).and_then(|l| l.porous.as_ref())
    }

    /// Exact planar extended-reaction reflection within the admitted material
    /// model. Positive incidence cosine is measured from the interface normal.
    /// Empty layers give the rigid wall. Source laws are never re-evaluated here.
    pub fn reflection(&self, cosine: f64) -> Result<LayeredReflection, MaterialError> {
        if !cosine.is_finite() || cosine <= 0.0 || cosine > 1.0 {
            return Err(bad("layered reflection needs incidence cosine in (0,1]"));
        }
        // Reference is Z0/cos(theta), the incident normal wave impedance.
        // Represent its ratio algebraically; forming Z0/cos could overflow.
        let (mut pressure, mut velocity) = (C64::ONE, C64::ZERO);
        for layer in self.layers.iter().rev() {
            let d = layer.source.thickness_m();
            if d == 0.0 {
                continue;
            }
            let (kn, ratio) = match layer.porous {
                None => (C64::from_re(self.k0 * cosine), C64::ONE),
                Some(response) => {
                    let (zc, k) = response.characteristic();
                    let kn = if cosine == 1.0 {
                        k
                    } else {
                        let tangential = self.k0 * det::sqrt((1.0-cosine)*(1.0+cosine));
                        let mut value = (k*k-C64::from_re(tangential*tangential)).sqrt();
                        if value.im > 0.0 { value = -value; }
                        value
                    };
                    let zn = zc * (k / kn);
                    (kn, divided(zn, self.z0).scale(cosine))
                }
            };
            if !finite(kn) || kn.re < 0.0 || kn.im > 0.0 || kn == C64::ZERO
                || !finite(ratio) || ratio.re <= 0.0
            {
                return Err(bad("layer normal wave impedance is unrepresentable or nonpassive"));
            }
            let (sum, difference) = bounded_propagator(kn, d)?;
            let inverse = ratio.recip();
            if !finite(inverse) {
                return Err(bad("layer inverse wave impedance is unrepresentable"));
            }
            // A fixed common power-of-two factor leaves headroom for complex
            // products. Contrast-dependent scaling would erase O(cosine^2)
            // entries before state normalization at extreme grazing angles.
            let diagonal = sum.scale(1.0/16.0);
            let off_diagonal = difference.scale(1.0/16.0);
            let upper = ratio * off_diagonal;
            let lower = inverse * off_diagonal;
            let next_p = diagonal*pressure + upper*velocity;
            let next_v = lower*pressure + diagonal*velocity;
            let normalization = magnitude_scale(next_p).max(magnitude_scale(next_v));
            if !finite(next_p) || !finite(next_v) || normalization == 0.0 {
                return Err(bad("layer transfer state is unrepresentable"));
            }
            pressure = divided(next_p, normalization);
            velocity = divided(next_v, normalization);
        }
        let denominator = pressure + velocity;
        let norm = denominator.norm_sq();
        if !norm.is_finite() || norm <= 0.0 {
            return Err(bad("layered incident pressure is unrepresentable"));
        }
        // Preserve the small complex correction when R approaches +/-1;
        // direct subtraction and division would lose grazing reflection phase.
        let reflection = if magnitude_scale(pressure) < 0.25*magnitude_scale(velocity) {
            -C64::ONE + (pressure/denominator).scale(2.0)
        } else if magnitude_scale(velocity) < 0.25*magnitude_scale(pressure) {
            C64::ONE - (velocity/denominator).scale(2.0)
        } else {
            (pressure-velocity)/denominator
        };
        // Evaluate actual surface work. 1-|R|^2 would erase weak absorption
        // whenever the reflected magnitude rounds to unity.
        let absorbed = 4.0*(pressure.re*velocity.re + pressure.im*velocity.im)/norm;
        let tolerance = 128.0*f64::EPSILON;
        if !finite(reflection) || !absorbed.is_finite() || absorbed < -tolerance
            || (reflection.norm_sq()+absorbed-1.0).abs() > tolerance
        {
            return Err(bad("layered reflection failed passive surface-work balance"));
        }
        Ok(LayeredReflection {
            pressure_reflection: reflection,
            absorbed_fraction: absorbed,
            pressure,
            normalized_velocity: velocity,
            cosine,
        })
    }
}

/// Bounded common-scaled cos/sin pair. The lossless branch is purely real /
/// imaginary, so a stack of fluid gaps retains exactly zero dissipative work.
fn bounded_propagator(kn: C64, thickness: f64) -> Result<(C64, C64), MaterialError> {
    let decay = 2.0*(kn.im*thickness);
    let attenuation = if decay == f64::NEG_INFINITY { 0.0 } else { det::exp(decay) };
    if attenuation == 0.0 {
        return Ok((C64::ONE, C64::ONE));
    }
    let phase = kn.re*thickness;
    if !phase.is_finite() {
        return Err(bad("layer phase is unrepresentable"));
    }
    let (sine, cosine) = det::sin_cos(phase);
    if kn.im == 0.0 {
        return Ok((C64::from_re(cosine), C64::new(0.0, sine)));
    }
    let loss = -det::expm1(decay);
    let cross = 2.0*attenuation*sine*cosine;
    // Stable 1-E and 1+E even at thin layers and near quarter-wave poles.
    Ok((
        C64::new(loss+2.0*attenuation*cosine*cosine, -cross),
        C64::new(loss+2.0*attenuation*sine*sine, cross),
    ))
}
