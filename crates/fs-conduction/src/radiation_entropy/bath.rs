//! Opaque gray/diffuse surface in a black, isotropic, view-factor-one bath.
//!
//! Outgoing photon occupations are epsilon*n_surface+(1-epsilon)*n_bath.
//! Entropy is evaluated AFTER this addition, never by adding entropy spectra.
//! Surface and remote bath irreversibility are separate nonnegative Bose
//! relative-entropy sums. This observation neither advances heat nor adds a
//! second cooling debit. Unpolarized chaotic modes and vacuum propagation are
//! explicit assumptions; directional illumination and coherent light are absent.

use super::{C1_NM, C2_NM_K, C3_NM, GrayPhotonEmission, PI, RadiationEntropyError,
    det, finite, log_mode_entropy, log_occupation, represented_exp};

/// Immutable surface/bath temperatures and reciprocal gray reflection.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GrayBlackBathExchange {
    surface: GrayPhotonEmission,
    bath_temperature_k: f64,
    log_bath_temperature: f64,
    log_reflectance: f64,
}

/// Hemispherical spectral fluxes per area of one surface side and per nm.
/// Positive heat leaves the surface. Photon entropy flux is NOT heat/T.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhotonBathSpectrum {
    /// Incoming blackbody energy, W/(m² nm).
    pub incoming_energy_w_m2_nm: f64,
    /// Combined reflected/emitted energy, W/(m² nm).
    pub outgoing_energy_w_m2_nm: f64,
    /// Original incoming photon entropy, W/(m² K nm).
    pub incoming_entropy_w_m2_k_nm: f64,
    /// Entropy of the combined outgoing occupation, W/(m² K nm).
    pub outgoing_entropy_w_m2_k_nm: f64,
    /// Net outward energy, evaluated before subtracting rounded spectra.
    pub heat_out_w_m2_nm: f64,
    /// Material surface entropy storage rate -heat/T_surface.
    pub surface_storage_w_m2_k_nm: f64,
    /// Remote black bath entropy storage rate +heat/T_bath.
    pub bath_storage_w_m2_k_nm: f64,
    /// Irreversibility at the emitting/absorbing/reflecting surface.
    pub surface_generation_w_m2_k_nm: f64,
    /// Irreversibility where the outgoing radiation thermalizes in the bath.
    pub bath_generation_w_m2_k_nm: f64,
    /// Independent stable heat*(1/T_bath-1/T_surface).
    pub total_generation_w_m2_k_nm: f64,
    /// Sum of local productions minus total; arithmetic, not physical heat.
    pub generation_residual_w_m2_k_nm: f64,
}

impl GrayPhotonEmission {
    /// Add an isotropic black reservoir without changing this surface model.
    pub fn in_black_bath(self, bath_temperature_k: f64)
        -> Result<GrayBlackBathExchange, RadiationEntropyError>
    {
        if !bath_temperature_k.is_finite() || bath_temperature_k <= 0.0 {
            return Err(RadiationEntropyError::Temperature);
        }
        Ok(GrayBlackBathExchange {
            surface: self, bath_temperature_k,
            log_bath_temperature: det::ln(bath_temperature_k),
            log_reflectance: if self.emissivity == 1.0 { f64::NEG_INFINITY }
                else { signed_log1p(-self.emissivity) },
        })
    }
}

impl GrayBlackBathExchange {
    /// Original material temperature in kelvin.
    pub const fn surface_temperature_k(self) -> f64 { self.surface.temperature_k }
    /// Original black reservoir temperature in kelvin.
    pub const fn bath_temperature_k(self) -> f64 { self.bath_temperature_k }
    /// Unchanged reciprocal gray emissivity.
    pub const fn emissivity(self) -> f64 { self.surface.emissivity }

    /// Allocation-free spectrum; small physical production is never recovered
    /// by subtracting two large entropy fluxes. Equal temperatures or an exact
    /// nonabsorbing reflector give exact zero production. A deep-tail zero can
    /// also mean floating-point underflow, as in the original emission closure.
    pub fn spectrum_nm(self, wavelength_nm: f64)
        -> Result<PhotonBathSpectrum, RadiationEntropyError>
    {
        if !wavelength_nm.is_finite() || wavelength_nm <= 0.0 {
            return Err(RadiationEntropyError::Wavelength);
        }
        let log_lambda = det::ln(wavelength_nm);
        let hot = self.surface.temperature_k.max(self.bath_temperature_k);
        let log_x = det::ln(C2_NM_K)-log_lambda-det::ln(hot);
        let direct = (C2_NM_K/wavelength_nm)/hot;
        let x = if direct.is_finite() && direct > 0.0 { direct } else { det::exp(log_x) };
        let modes = self.modes(x, log_x)?;
        let le = det::ln(PI*C1_NM)-5.0*log_lambda;
        let ls = det::ln(PI*C3_NM)-4.0*log_lambda;
        let heat = self.direction()*represented_exp(le+modes.log_net_occupation)?;
        let surface_generation = represented_exp(ls+modes.log_surface_generation)?;
        let bath_generation = represented_exp(ls+modes.log_bath_generation)?;
        let total = represented_exp(ls+modes.log_total_generation)?;
        Ok(PhotonBathSpectrum {
            incoming_energy_w_m2_nm: represented_exp(le+modes.log_incoming)?,
            outgoing_energy_w_m2_nm: represented_exp(le+modes.log_outgoing)?,
            incoming_entropy_w_m2_k_nm: represented_exp(ls+log_mode_entropy(modes.log_incoming))?,
            outgoing_entropy_w_m2_k_nm: represented_exp(ls+log_mode_entropy(modes.log_outgoing))?,
            heat_out_w_m2_nm: heat,
            surface_storage_w_m2_k_nm: -self.direction()*represented_exp(le+modes.log_net_occupation-self.surface.log_temperature)?,
            bath_storage_w_m2_k_nm: self.direction()*represented_exp(le+modes.log_net_occupation-self.log_bath_temperature)?,
            surface_generation_w_m2_k_nm: surface_generation,
            bath_generation_w_m2_k_nm: bath_generation,
            total_generation_w_m2_k_nm: total,
            generation_residual_w_m2_k_nm: finite((surface_generation+bath_generation)-total)?,
        })
    }

    fn direction(self) -> f64 {
        if self.surface.temperature_k > self.bath_temperature_k { 1.0 }
        else if self.surface.temperature_k < self.bath_temperature_k { -1.0 }
        else { 0.0 }
    }

    fn modes(self, x_hot: f64, log_x_hot: f64) -> Result<ModeExchange, RadiationEntropyError> {
        let ts = self.surface.temperature_k;
        let tb = self.bath_temperature_k;
        let hot = ts.max(tb);
        let cold = ts.min(tb);
        let log_hot = self.surface.log_temperature.max(self.log_bath_temperature);
        let log_cold = self.surface.log_temperature.min(self.log_bath_temperature);
        let log_x_cold = log_x_hot+(log_hot-log_cold);
        let raw_cold = x_hot*(hot/cold);
        let x_cold = if raw_cold.is_finite() && raw_cold > 0.0 { raw_cold }
            else { det::exp(log_x_cold) };
        if x_hot.is_infinite() { return Ok(ModeExchange::zero()); }
        if !x_hot.is_finite() || !x_cold.is_finite() { return Err(RadiationEntropyError::Unrepresentable); }
        let occupation = |x, lx| if x == 0.0 { -lx } else { log_occupation(x, 0.0) };
        let lh = occupation(x_hot, log_x_hot);
        let lc = occupation(x_cold, log_x_cold);
        let (ns, nb, xb) = if ts >= tb { (lh,lc,x_cold) } else { (lc,lh,x_hot) };
        if ts == tb || self.surface.emissivity == 0.0 {
            return Ok(ModeExchange { log_incoming: nb, log_outgoing: nb, ..ModeExchange::zero() });
        }
        // Preserve the temperature departure BEFORE forming rounded Planck logs.
        let delta_x = x_hot*((hot-cold)/cold);
        let dx = -self.direction()*delta_x;
        let z = if dx.abs() < 0.5 && xb > 0.0 {
            det::expm1(dx)/(-det::expm1(-xb))
        } else { f64::NAN };
        let d = if z.is_finite() && z.abs() < 0.5 { -signed_log1p(z) }
            else if x_hot == 0.0 && x_cold == 0.0 {
                let relative = (ts-tb)/tb;
                if relative.abs() < 0.5 { signed_log1p(relative) }
                else { self.surface.log_temperature-self.log_bath_temperature }
            } else { ns-nb };
        if !d.is_finite() || d == 0.0 { return Err(RadiationEntropyError::Unrepresentable); }
        let e = self.surface.emissivity;
        let out_bath = mixed_log_ratio(d, e)?;
        let out_surface = mixed_log_ratio(-d, 1.0-e)?;
        let no = if d > 0.0 { ns+out_surface } else { nb+out_bath };
        let mut generation = self.surface.log_emissivity+bose_divergence(ns, -d)?;
        if out_surface != 0.0 {
            generation = log_add(generation, self.surface.log_emissivity+bose_divergence(no, -out_surface)?);
        }
        if e < 1.0 {
            generation = log_add(generation, self.log_reflectance+bose_divergence(no, -out_bath)?);
        }
        let net = self.surface.log_emissivity+if d > 0.0 {
            ns+det::ln(-det::expm1(-d))
        } else { nb+det::ln(-det::expm1(d)) };
        let log_delta_x = log_x_hot+det::ln(hot-cold)-log_cold;
        Ok(ModeExchange {
            log_incoming: nb, log_outgoing: no,
            log_net_occupation: net,
            log_surface_generation: generation,
            log_bath_generation: bose_divergence(nb, out_bath)?,
            log_total_generation: net+log_delta_x,
        })
    }
}

#[derive(Clone, Copy)]
struct ModeExchange {
    log_incoming: f64, log_outgoing: f64, log_net_occupation: f64,
    log_surface_generation: f64, log_bath_generation: f64, log_total_generation: f64,
}
impl ModeExchange {
    fn zero() -> Self {
        Self { log_incoming: f64::NEG_INFINITY, log_outgoing: f64::NEG_INFINITY,
            log_net_occupation: f64::NEG_INFINITY, log_surface_generation: f64::NEG_INFINITY,
            log_bath_generation: f64::NEG_INFINITY, log_total_generation: f64::NEG_INFINITY }
    }
}

fn signed_log1p(x: f64) -> f64 {
    let y = 1.0+x;
    if y == 1.0 { x } else { det::ln(y)*(x/(y-1.0)) }
}
fn log_add(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY { return b; }
    if b == f64::NEG_INFINITY { return a; }
    a.max(b)+signed_log1p(det::exp(a.min(b)-a.max(b)))
}
fn softplus(x: f64) -> f64 { x.max(0.0)+signed_log1p(det::exp(-x.abs())) }
fn log_softplus(x: f64) -> f64 { if x < -36.0 { x } else { det::ln(softplus(x)) } }
fn log_difference(a: f64, b: f64) -> Result<f64, RadiationEntropyError> {
    if b == f64::NEG_INFINITY { return Ok(a); }
    if a <= b { return Err(RadiationEntropyError::Unrepresentable); }
    finite(a+det::ln(-det::expm1(b-a)))
}
fn mixed_log_ratio(d: f64, weight: f64) -> Result<f64, RadiationEntropyError> {
    if weight == 0.0 { return Ok(0.0); }
    if weight == 1.0 { return Ok(d); }
    let result = if d.abs() < 0.5 { signed_log1p(weight*det::expm1(d)) }
        else { log_add(det::ln(weight)+d, signed_log1p(-weight)) };
    if result == 0.0 && d != 0.0 { return Err(RadiationEntropyError::Unrepresentable); }
    finite(result)
}

// log D(a||b), D=a ln(a/b)-(1+a)ln((1+a)/(1+b)). The separately retained
// d=ln(a/b) protects a near-equilibrium difference smaller than an ulp of ln b.
// The local series factors b/(1+b) and delta² before summation; unlike two
// convex remainders it remains conditioned for both dilute and occupied modes.
fn bose_divergence(log_b: f64, d: f64) -> Result<f64, RadiationEntropyError> {
    if d == 0.0 { return Ok(f64::NEG_INFINITY); }
    if d == f64::NEG_INFINITY { return Ok(log_softplus(log_b)); }
    let log_a = finite(log_b+d)?;
    if d.abs() < 0.1 {
        let delta = det::expm1(d);
        let log_r = -softplus(-log_b);
        let r = det::exp(log_r);
        let (mut power, mut geometric) = (1.0,1.0);
        let mut sum = super::Sum::default();
        for k in 2..22 {
            sum.add(power*geometric/f64::from(k*(k-1)))?;
            power *= -delta;
            geometric = 1.0+r*geometric;
        }
        return finite(log_r+2.0*det::ln(delta.abs())+det::ln(sum.total()?));
    }
    if d > 0.0 {
        let s = det::ln(-det::expm1(-d));
        let c = log_softplus(s-log_b-softplus(-log_a));
        let l = log_softplus(log_a+s-softplus(log_b));
        log_difference(log_a+c,l)
    } else {
        let s = det::ln(-det::expm1(d));
        let l = log_softplus(log_b+s-softplus(log_a));
        let c = log_softplus(s-log_a-softplus(-log_b));
        log_difference(l,log_a+c)
    }
}

#[cfg(test)]
mod tests {
use super::*;
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance*a.abs().max(b.abs()).max(1e-300), "{a:e} != {b:e}");
}
fn bath(ts: f64, tb: f64, e: f64) -> GrayBlackBathExchange {
    GrayPhotonEmission::new(ts,e).unwrap().in_black_bath(tb).unwrap()
}
#[test]
fn black_bath_reflection_and_emission_share_outgoing_modes() {
    let (ts,tb,e,lambda) = (900.0,300.0,0.3,10_000.0);
    let s = bath(ts,tb,e).spectrum_nm(lambda).unwrap();
    let emitted = GrayPhotonEmission::new(ts,e).unwrap().spectrum_nm(lambda).unwrap();
    let reflected = GrayPhotonEmission::new(tb,1.0-e).unwrap().spectrum_nm(lambda).unwrap();
    close(s.outgoing_energy_w_m2_nm,PI*(emitted.energy_w_m2_sr_nm+reflected.energy_w_m2_sr_nm),2e-12);
    assert!(s.outgoing_entropy_w_m2_k_nm < PI*(emitted.entropy_w_m2_sr_k_nm+reflected.entropy_w_m2_sr_k_nm));
    let ns = 1.0/(C2_NM_K/(lambda*ts)).exp_m1();
    let nb = 1.0/(C2_NM_K/(lambda*tb)).exp_m1();
    let no = e*ns+(1.0-e)*nb;
    let g = |n:f64| (1.0+n)*n.ln_1p()-n*n.ln();
    let scale = PI*C3_NM/lambda.powi(4);
    close(s.outgoing_entropy_w_m2_k_nm,scale*g(no),2e-12);
    close(s.surface_generation_w_m2_k_nm,scale*(g(no)-g(nb)-C2_NM_K/(lambda*ts)*(no-nb)),2e-11);
    close(s.bath_generation_w_m2_k_nm,scale*(g(nb)-g(no)+C2_NM_K/(lambda*tb)*(no-nb)),2e-11);
}
#[test]
fn both_directions_close_material_and_photon_accounts_without_clamping() {
    for (ts,tb) in [(900.0,300.0),(300.0,900.0)] {
        for e in [1e-8,0.01,0.5,0.99,1.0] {
            for lambda in [500.0,5_000.0,50_000.0,500_000.0] {
                let s = bath(ts,tb,e).spectrum_nm(lambda).unwrap();
                assert!(s.surface_generation_w_m2_k_nm>0.0 && s.bath_generation_w_m2_k_nm>0.0);
                assert!(s.heat_out_w_m2_nm*(ts-tb)>0.0);
                close(s.total_generation_w_m2_k_nm,s.heat_out_w_m2_nm*(1.0/tb-1.0/ts),2e-11);
                close(s.surface_generation_w_m2_k_nm+s.bath_generation_w_m2_k_nm,s.total_generation_w_m2_k_nm,2e-10);
                close(s.surface_storage_w_m2_k_nm,-s.heat_out_w_m2_nm/ts,2e-12);
                close(s.bath_storage_w_m2_k_nm,s.heat_out_w_m2_nm/tb,2e-12);
            }
        }
    }
}
#[test]
fn equilibrium_and_nonabsorbing_limits_are_exact() {
    for e in [0.0,1e-200,0.5,1.0] {
        for model in [bath(400.0,400.0,e),bath(900.0,400.0,0.0)] {
            let s = model.spectrum_nm(10_000.0).unwrap();
            assert_eq!(s.outgoing_energy_w_m2_nm,s.incoming_energy_w_m2_nm);
            assert_eq!(s.outgoing_entropy_w_m2_k_nm,s.incoming_entropy_w_m2_k_nm);
            assert_eq!(s.heat_out_w_m2_nm,0.0);
            assert_eq!(s.surface_generation_w_m2_k_nm,0.0);
            assert_eq!(s.bath_generation_w_m2_k_nm,0.0);
            assert_eq!(s.total_generation_w_m2_k_nm,0.0);
        }
    }
}
#[test]
fn adjacent_temperatures_keep_second_order_production() {
    let tb = 300.0_f64;
    for bits in [tb.to_bits()-1,tb.to_bits()+1] {
        let ts = f64::from_bits(bits);
        let e = 0.4;
        let s = bath(ts,tb,e).spectrum_nm(C2_NM_K/(2.0*tb)).unwrap();
        assert!(s.surface_generation_w_m2_k_nm>0.0 && s.bath_generation_w_m2_k_nm>0.0);
        close(s.bath_generation_w_m2_k_nm/s.total_generation_w_m2_k_nm,e/2.0,2e-11);
        close(s.surface_generation_w_m2_k_nm/s.total_generation_w_m2_k_nm,1.0-e/2.0,2e-11);
    }
}
#[test]
fn relative_entropy_retains_dilute_occupied_and_nearby_modes() {
    let d = 2.0_f64.ln();
    close(bose_divergence(0.0,d).unwrap().exp(),2.0*d-3.0*1.5_f64.ln(),2e-12);
    close(bose_divergence(700.0,d).unwrap().exp(),1.0-d,5e-12);
    close(bose_divergence(-1000.0,d).unwrap()+1000.0,(2.0*d-1.0).ln(),5e-12);
    for lb in [-100.0,0.0,100.0] {
        let d = 1e-16_f64;
        let log_reference = -softplus(-lb)-2.0_f64.ln()+2.0*d.ln();
        close(bose_divergence(lb,d).unwrap(),log_reference,2e-14);
    }
}
#[test]
fn invalid_coordinates_and_unrepresentable_departures_refuse() {
    let emitter = GrayPhotonEmission::new(300.0,0.5).unwrap();
    for t in [0.0,-1.0,f64::NAN,f64::INFINITY] {
        assert_eq!(emitter.in_black_bath(t),Err(RadiationEntropyError::Temperature));
    }
    for lambda in [0.0,-1.0,f64::NAN,f64::INFINITY] {
        assert_eq!(bath(400.0,300.0,0.5).spectrum_nm(lambda),Err(RadiationEntropyError::Wavelength));
    }
    assert!(matches!(mixed_log_ratio(1e-100,f64::from_bits(1)),Err(RadiationEntropyError::Unrepresentable)));
}
}
