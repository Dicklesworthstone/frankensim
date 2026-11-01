//! Outward evaluation of the original gray propagation-angle laws.
use super::*;

pub(super) fn enclose(phase: &PhaseFunction, mu: Interval, cx: &Cx<'_>) -> Result<Interval, ScatterEnclosureError> {
    cx.checkpoint()?;
    let pi = Interval::new(fs_math::next_down(core::f64::consts::PI), fs_math::next_up(core::f64::consts::PI));
    let inverse_sphere = point(1.0)/(point(4.0)*pi);
    let result = match phase {
        PhaseFunction::Isotropic => inverse_sphere,
        PhaseFunction::HenyeyGreenstein { g } => {
            let a = point(g.abs());
            let gap = point(1.0)-a;
            let angular = if *g >= 0.0 { point(1.0)-mu } else { point(1.0)+mu };
            let denominator = add(mul(gap,gap), mul(point(2.0)*a, nonnegative(angular)));
            if denominator.lo() <= 0.0 { return Err(ScatterEnclosureError::NumericalRange); }
            mul(mul(gap,point(1.0)+a),inverse_sphere)/denominator/denominator.sqrt()
        }
        PhaseFunction::Rayleigh => mul(point(0.75)*inverse_sphere, point(1.0)+square(mu)),
        PhaseFunction::RayleighDepolarized { depolarization: d } => {
            mul(point(1.5)*inverse_sphere,
                add(point(1.0)+point(*d), mul(point(1.0)-point(*d),square(mu))))/(point(2.0)+point(*d))
        }
        PhaseFunction::Tabulated(table) => {
            let knots = table.knots_per_sr();
            let first = knots.partition_point(|&(x,_)| x < mu.lo()).saturating_sub(1).min(knots.len()-2);
            let last = knots.partition_point(|&(x,_)| x <= mu.hi()).min(knots.len()-1);
            let mut range: Option<Interval> = None;
            for i in first..last.max(first+1) {
                cx.checkpoint()?;
                let (a,fa)=knots[i]; let (b,fb)=knots[i+1];
                let lo=mu.lo().max(a); let hi=mu.hi().min(b);
                if lo>hi { continue; }
                let denominator=point(b)-point(a);
                if denominator.lo()<=0.0 { return Err(ScatterEnclosureError::NumericalRange); }
                let value = if fa==fb { point(fa) } else {
                    let t=(Interval::new(lo,hi)-point(a))/denominator;
                    let t=Interval::new(t.lo().max(0.0),t.hi().min(1.0));
                    add(mul(nonnegative(point(1.0)-t),point(fa)),mul(t,point(fb)))
                };
                range=Some(range.map_or(value,|old|old.hull(value)));
            }
            range.ok_or(ScatterEnclosureError::NumericalRange)?
        }
        PhaseFunction::SpectralTable(_) | PhaseFunction::SpectralMixture(_) => return Err(ScatterEnclosureError::UnsupportedPhase),
    };
    if !finite(result) { return Err(ScatterEnclosureError::NumericalRange); }
    Ok(nonnegative(result))
}
fn square(x: Interval) -> Interval {
    let lo=if x.contains_zero() {0.0} else {x.lo().abs().min(x.hi().abs())};
    let hi=x.lo().abs().max(x.hi().abs());
    Interval::new(nonnegative(point(lo)*point(lo)).lo(),(point(hi)*point(hi)).hi())
}
