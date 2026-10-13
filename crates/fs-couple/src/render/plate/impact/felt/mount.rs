//! A shared finite-mass carrier for preloaded, compression-only felt supports.
//!
//! Each pad keeps its own geometric row, material and history. All pad backings
//! move with ONE carrier, rather than independent squeeze jaws or a fixed wall.
//! The original body coordinates remain unchanged; the carrier is appended.
//! This is cold composition of existing linear inertia, WoolFelt, Kelvin and
//! pHS owners, not another constitutive law or a partitioned force-feedback loop.
use super::FeltPad;
use super::super::{BodyPotential, ImpactBody, ImpactError, MAX_IMPACT_MODES, invalid};

/// One fixed-axis translating carrier, with explicit passive ground attachment.
/// Mass excludes the instrument, whose original inertia remains in its body.
/// Stiffness/rest position describe a physical support, not an audio parameter.
#[derive(Debug, Clone, Copy)]
pub struct Carrier {
    pub mass_kg: f64,
    pub stiffness_n_m: f64,
    pub resistance_n_s_m: f64,
    /// Initial displacement from the support spring's rest position [m].
    pub initial_position_m: f64,
    pub initial_velocity_m_s: f64,
}

/// One felt patch and its backing's signed fixed-axis compression derivative.
#[derive(Debug, Clone)]
pub struct MountSite {
    /// Weights address the ORIGINAL structural basis, before the carrier.
    /// Precompression is defined at zero original coordinates AND zero carrier.
    /// Creep coefficients belong to this patch, not an entire washer.
    pub pad: FeltPad,
    /// d(compression)/d(carrier displacement), the signed projection of the
    /// carrier translation onto this pad's compression direction, in [-1,1].
    /// Opposite backings normally have opposite signs. A value of zero denotes
    /// a pad orthogonal to carrier travel; at least one site must couple.
    pub compression_per_carrier_m: f64,
}

/// A physical force/motion port, with no independent dynamic state.
#[derive(Debug, Clone, Copy)]
pub struct CarrierPort {
    pub coordinate: usize,
    pub inverse_sqrt_mass: f64,
}
impl CarrierPort {
    /// Observe physical displacement/velocity from an accepted interleaved q,p
    /// state. Kelvin and other memory coordinates may follow the body prefix.
    pub fn motion(&self, state: &[f64]) -> Result<[f64; 2], ImpactError> {
        let at=self.coordinate.checked_mul(2)
            .ok_or_else(||invalid("carrier state address overflow"))?;
        let end=at.checked_add(2).ok_or_else(||invalid("carrier state address overflow"))?;
        let values=state.get(at..end).ok_or_else(||invalid("carrier is outside the state"))?;
        let result=[values[0]*self.inverse_sqrt_mass, values[1]*self.inverse_sqrt_mass];
        if !self.inverse_sqrt_mass.is_finite() || self.inverse_sqrt_mass<=0.0
            || result.iter().any(|v|!v.is_finite()) {
            return Err(invalid("carrier motion is not representable"));
        }
        Ok(result)
    }
}

/// Append `body` after the original bodies and use these extended pad rows.
/// Extend EVERY other contact, volume, force, damper and observation row with
/// zero before constructing ImpactSystem. Do not leave the old fixed pads in
/// the assembly: these pads REPLACE them, without adding their inertia twice.
#[derive(Debug)]
pub struct MountedFelt {
    pub body: ImpactBody,
    pub pads: Vec<FeltPad>,
    pub port: CarrierPort,
}
impl MountedFelt {
    /// Compile a complete reciprocal support. There is no per-step adapter work:
    /// carrier and pad histories participate in the same existing implicit solve.
    /// Nonnegative ground stiffness/resistance may both be zero for a free
    /// carrier. Gravity, finite rotations, stand flexural modes and a solved
    /// preload equilibrium are NOT inferred from these inputs.
    ///
    /// # Errors
    /// Invalid geometry rows, passive coefficients, arithmetic or owner budgets.
    pub fn new(structural_modes: usize, carrier: Carrier, sites: &[MountSite])
        -> Result<Self, ImpactError>
    {
        if structural_modes==0 || structural_modes>=MAX_IMPACT_MODES
            || sites.is_empty() || sites.len()>16
            || !carrier.stiffness_n_m.is_finite() || carrier.stiffness_n_m<0.0
            || !carrier.resistance_n_s_m.is_finite() || carrier.resistance_n_s_m<0.0 {
            return Err(invalid("felt carrier needs a bounded basis, pads and passive ground coefficients"));
        }
        let (mut body,weight)=ImpactBody::free_mass(carrier.mass_kg,
            carrier.initial_position_m,carrier.initial_velocity_m_s)?;
        let frequency_squared=carrier.stiffness_n_m/carrier.mass_kg;
        let drag=carrier.resistance_n_s_m/carrier.mass_kg;
        if !weight.is_finite() || weight<=0.0 || !frequency_squared.is_finite() || !drag.is_finite()
            || (carrier.stiffness_n_m>0.0 && frequency_squared==0.0)
            || (carrier.resistance_n_s_m>0.0 && drag==0.0)
            || !body.initial[0].displacement_m_sqrt_kg.is_finite()
            || !body.initial[0].velocity_m_sqrt_kg_per_s.is_finite()
            || (carrier.initial_position_m!=0.0 && body.initial[0].displacement_m_sqrt_kg==0.0)
            || (carrier.initial_velocity_m_s!=0.0 && body.initial[0].velocity_m_sqrt_kg_per_s==0.0) {
            return Err(invalid("felt carrier inertia or support cannot be represented"));
        }
        body.potential=BodyPotential::Linear(vec![frequency_squared.sqrt()]);
        body.damping_per_s[0]=drag;
        let mut pads=Vec::with_capacity(sites.len());let mut coupled=false;
        for site in sites {
            site.pad.validate(structural_modes)?;
            let projection=site.compression_per_carrier_m;
            if !projection.is_finite() || projection.abs()>1.0
                || site.pad.weights.iter().all(|v|*v==0.0) {
                return Err(invalid("felt mount needs an actual structural row and a unit-axis compression projection"));
            }
            let coefficient=projection*weight;
            if !coefficient.is_finite() || (projection!=0.0 && coefficient==0.0) {
                return Err(invalid("felt carrier contact projection cannot be represented"));
            }
            coupled|=projection!=0.0;
            let mut pad=site.pad.clone();pad.weights.push(coefficient);
            pad.validate(structural_modes+1)?;pads.push(pad);
        }
        if !coupled {return Err(invalid("all felt pads are decoupled from the carrier"));}
        Ok(Self {body,pads,port:CarrierPort {coordinate:structural_modes,inverse_sqrt_mass:weight}})
    }
}

#[cfg(test)]
mod tests;
