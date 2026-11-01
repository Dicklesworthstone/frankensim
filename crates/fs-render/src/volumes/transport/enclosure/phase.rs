//! Numerical enclosures of the ORIGINAL propagation-angle phase laws.
//! Spectral interpolation and independent-population weights stay in their
//! native owners. These are retained-model bounds, not phase-data calibration.

use super::{Interval, VolumeTransportError, checkpoint, finite, nonnegative};
use crate::volumes::transport::scattering::PhaseFunction;
use crate::volumes::transport::MAX_OPTICAL_MIXTURE_COMPONENTS;
use fs_exec::Cx;

/// Version of angular, spectral and mixture numerical enclosure semantics.
pub const PHASE_ENCLOSURE_SEMANTICS_VERSION: u32 = 1;

impl PhaseFunction {
    /// Enclose a wavelength-independent law over a propagation-cosine interval.
    /// Spectral laws refuse `Wavelength`; use `enclosure_at` for those laws.
    /// No allocation, fitted angular law, or measured-error assumption is used.
    pub fn enclosure_cosine(&self, cx: &Cx<'_>, cosine: Interval)
        -> Result<Interval, VolumeTransportError> {
        checkpoint(cx)?;
        validate_cosine(cosine)?;
        gray(self, cosine, cx)
    }

    /// Enclose the retained phase interpolant in sr^-1 at a vacuum wavelength.
    /// Tables use outward spectral interpolation; mixtures use the original
    /// loading times scattering cross section, never equal component weights.
    /// Wavelength extrapolation refuses. Angular tables have their native knot
    /// cap, and mixtures cannot recurse into mixtures, so work is bounded.
    pub fn enclosure_at(&self, cx: &Cx<'_>, wavelength_nm: f64, cosine: Interval)
        -> Result<Interval, VolumeTransportError> {
        checkpoint(cx)?;
        validate_cosine(cosine)?;
        self.validate_wavelength(wavelength_nm)?;
        let result = match self {
            Self::SpectralTable(table) => {
                let knots = table.knots();
                let upper = knots.partition_point(|(lambda, _)| *lambda < wavelength_nm);
                if upper < knots.len() && knots[upper].0 == wavelength_nm {
                    return gray(&knots[upper].1, cosine, cx);
                }
                let (left, a) = &knots[upper - 1];
                let (right, b) = &knots[upper];
                let denominator = p(*right) - p(*left);
                if denominator.lo() <= 0.0 { return Err(VolumeTransportError::Unrepresentable); }
                let t = (p(wavelength_nm) - p(*left)) / denominator;
                let t = Interval::new(t.lo().max(0.0), t.hi().min(1.0));
                add(mul(nonnegative(p(1.0)-t), gray(a, cosine, cx)?),
                    mul(t, gray(b, cosine, cx)?))
            }
            Self::SpectralMixture(mixture) => {
                let mut weights = [p(0.0); MAX_OPTICAL_MIXTURE_COMPONENTS];
                let mut scale = 0.0_f64;
                for (i, component) in mixture.components().iter().enumerate() {
                    checkpoint(cx)?;
                    if component.mass_weight == 0.0 { continue; }
                    weights[i] = mul(p(component.mass_weight),
                        component.scattering.coefficient_enclosure(wavelength_nm)?);
                    if !finite(weights[i]) { return Err(VolumeTransportError::Unrepresentable); }
                    scale = scale.max(weights[i].hi());
                }
                if scale == 0.0 { return gray(&Self::Isotropic, cosine, cx); }
                let mut numerator = p(0.0);
                let mut denominator = p(0.0);
                let mut envelope: Option<Interval> = None;
                for (component, weight) in mixture.components().iter().zip(weights) {
                    checkpoint(cx)?;
                    if weight.hi() == 0.0 { continue; }
                    let phase = component.phase.enclosure_at(cx, wavelength_nm, cosine)?;
                    envelope = Some(envelope.map_or(phase, |old| old.hull(phase)));
                    let scaled = nonnegative(weight / p(scale));
                    numerator = add(numerator, mul(scaled, phase));
                    denominator = add(denominator, scaled);
                }
                let envelope = envelope.ok_or(VolumeTransportError::Unrepresentable)?;
                if denominator.lo() <= 0.0 {
                    // Positive weights may be below the subnormal resolution.
                    // A zero-scattering lane has the native isotropic convention.
                    envelope.hull(gray(&Self::Isotropic, cosine, cx)?)
                } else {
                    (numerator / denominator).intersect(envelope)
                        .ok_or(VolumeTransportError::Unrepresentable)?
                }
            }
            _ => return gray(self, cosine, cx),
        };
        if !finite(result) { return Err(VolumeTransportError::Unrepresentable); }
        Ok(nonnegative(result))
    }
}

fn validate_cosine(mu: Interval) -> Result<(), VolumeTransportError> {
    if !finite(mu) || mu.lo() < -1.0 || mu.hi() > 1.0 { Err(VolumeTransportError::Phase) }
    else { Ok(()) }
}

// Preserve the former homogeneous enclosure's operation order verbatim.
fn gray(phase: &PhaseFunction, mu: Interval, cx: &Cx<'_>) -> Result<Interval, VolumeTransportError> {
    checkpoint(cx)?;
    let pi = Interval::new(fs_math::next_down(core::f64::consts::PI), fs_math::next_up(core::f64::consts::PI));
    let inverse_sphere = p(1.0)/(p(4.0)*pi);
    let result = match phase {
        PhaseFunction::Isotropic => inverse_sphere,
        PhaseFunction::HenyeyGreenstein { g } => {
            if !g.is_finite() || g.abs() >= 1.0 { return Err(VolumeTransportError::Phase); }
            let a = p(g.abs());
            let gap = p(1.0)-a;
            let angular = if *g >= 0.0 { p(1.0)-mu } else { p(1.0)+mu };
            let denominator = add(mul(gap,gap), mul(p(2.0)*a, nonnegative(angular)));
            if denominator.lo() <= 0.0 { return Err(VolumeTransportError::Unrepresentable); }
            mul(mul(gap,p(1.0)+a),inverse_sphere)/denominator/denominator.sqrt()
        }
        PhaseFunction::Rayleigh => mul(p(0.75)*inverse_sphere, p(1.0)+square(mu)),
        PhaseFunction::RayleighDepolarized { depolarization: d } => {
            if !d.is_finite() || !(0.0..6.0/7.0).contains(d) { return Err(VolumeTransportError::Phase); }
            mul(p(1.5)*inverse_sphere,
                add(p(1.0)+p(*d), mul(p(1.0)-p(*d),square(mu))))/(p(2.0)+p(*d))
        }
        PhaseFunction::Tabulated(table) => {
            let knots = table.knots_per_sr();
            let first = knots.partition_point(|&(x,_)| x < mu.lo()).saturating_sub(1).min(knots.len()-2);
            let last = knots.partition_point(|&(x,_)| x <= mu.hi()).min(knots.len()-1);
            let mut range: Option<Interval> = None;
            for i in first..last.max(first+1) {
                checkpoint(cx)?;
                let (a,fa)=knots[i]; let (b,fb)=knots[i+1];
                let lo=mu.lo().max(a); let hi=mu.hi().min(b);
                if lo>hi { continue; }
                let denominator=p(b)-p(a);
                if denominator.lo()<=0.0 { return Err(VolumeTransportError::Unrepresentable); }
                let value = if fa==fb { p(fa) } else {
                    let t=(Interval::new(lo,hi)-p(a))/denominator;
                    let t=Interval::new(t.lo().max(0.0),t.hi().min(1.0));
                    add(mul(nonnegative(p(1.0)-t),p(fa)),mul(t,p(fb)))
                };
                range=Some(range.map_or(value,|old|old.hull(value)));
            }
            range.ok_or(VolumeTransportError::Unrepresentable)?
        }
        PhaseFunction::SpectralTable(_) | PhaseFunction::SpectralMixture(_) => return Err(VolumeTransportError::Wavelength),
    };
    if !finite(result) { return Err(VolumeTransportError::Unrepresentable); }
    Ok(nonnegative(result))
}
fn square(x: Interval) -> Interval {
    let lo=if x.contains_zero() {0.0} else {x.lo().abs().min(x.hi().abs())};
    let hi=x.lo().abs().max(x.hi().abs());
    Interval::new(nonnegative(p(lo)*p(lo)).lo(),(p(hi)*p(hi)).hi())
}
fn p(v: f64) -> Interval { Interval::point(v) }
fn add(a: Interval, b: Interval) -> Interval {
    if a.hi()==0.0 { b } else if b.hi()==0.0 { a } else { nonnegative(a+b) }
}
fn mul(a: Interval, b: Interval) -> Interval {
    if a.hi()==0.0 || b.hi()==0.0 { p(0.0) } else { nonnegative(a*b) }
}

#[cfg(test)]
mod tests;
