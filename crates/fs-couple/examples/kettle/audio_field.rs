//! Dispatch only: each acoustic owner keeps its original pressure representation.
//! Area-mean traces are for P0 surface WORK, never a substitute receiver field.
use super::{C64, Error, Medium, SpherePanels, exterior_pressure_at_points_integrated};
use fs_bem::{helmholtz::RadiationSolution, reciprocal_radiation::{self, Prepared}};
use fs_exec::CancelGate;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Image { Collocation, Reciprocal }

pub(super) enum Field<'p,'s> {
    Collocation(RadiationSolution),
    Reciprocal(reciprocal_radiation::Solution<'p,'s>),
}
impl Field<'_, '_> {
    pub(super) fn trace(&self) -> &[C64] {
        match self { Self::Collocation(s)=>&s.pressure, Self::Reciprocal(s)=>&s.pressure }
    }
    /// Preserve each owner's diagnostics; no synthetic collocation solution is
    /// made from the P1 means or substituted after a failed power check.
    pub(super) fn admit(&self, panels:usize) -> Result<f64,Error> {
        let ppw=match self {
            Self::Collocation(s)=>{
                if !s.radiated_power_roundoff_interval.1.is_finite() || s.radiated_power_roundoff_interval.1<0.0 {
                    return Err("collocation source power refused".into());
                }
                s.panels_per_wavelength
            },
            Self::Reciprocal(s)=>{
                if !s.surface_power.re.is_finite() || !s.surface_power.im.is_finite()
                    || !s.outgoing_power_w.is_finite() || s.outgoing_power_w<0.0
                    || !s.backward_residual.is_finite() || !(0.0..=1e-10).contains(&s.backward_residual)
                    || !s.flux_balance.is_finite() || s.flux_balance.abs()>1e-8 {
                    return Err("reciprocal source flux or backward residual refused".into());
                }
                s.panels_per_wavelength
            },
        };
        if !ppw.is_finite() || ppw<6.0 || self.trace().len()!=panels
            || self.trace().iter().any(|v|!v.re.is_finite() || !v.im.is_finite()) {
            return Err("exterior resolution or pressure trace refused".into());
        }
        Ok(ppw)
    }
    pub(super) fn pressure_at(&self,surface:&SpherePanels,medium:Medium,point:[f64;3]) -> Result<C64,Error> {
        Ok(match self {
            Self::Collocation(s)=>exterior_pressure_at_points_integrated(surface,s,medium,&[point])?[0],
            // This method integrates the original P1 double layer and P0
            // auxiliary density. The receiver gets no extra distance gain.
            Self::Reciprocal(s)=>s.pressure_at(point)?,
        })
    }
}

/// Bound the complete sweep; the existing owner enforces its per-frequency
/// budget, P1/P0 dimension, cache, wavelength and interior-resonance chart.
/// Rejection propagates without falling back to a different equation.
pub(super) fn prepare<'s>(surface:&'s SpherePanels,k_max:f64,sources:usize,samples:usize,
    maximum_work:u64,gate:&CancelGate) -> Result<Prepared<'s>,Error> {
    if gate.is_requested() {return Err("reciprocal kettle preparation cancelled".into());}
    if samples==0 || !(1..=32).contains(&sources) || maximum_work==0 || maximum_work>100_000_000_000_000 {
        return Err("invalid reciprocal kettle sweep budget".into());
    }
    let triangles=surface.triangles().ok_or("reciprocal kettle needs explicit exterior triangles")?;
    let mut low=[f64::INFINITY;3];let mut high=[f64::NEG_INFINITY;3];
    for p in triangles.iter().flatten() {for d in 0..3 {low[d]=low[d].min(p[d]);high[d]=high[d].max(p[d]);}}
    let width:[f64;3]=core::array::from_fn(|d|(high[d]-low[d]).next_up());
    let diameter=width[0].hypot(width[1]).hypot(width[2]);
    let angular_refinement=(0..=2_u32).find(|r|k_max*diameter<=4.0*f64::from(1_u32<<*r))
        .ok_or("reciprocal kettle band exceeds angular quadrature budget")?;
    let per_frequency=maximum_work/samples as u64;
    let source=Prepared::new(surface,k_max,reciprocal_radiation::Budget {
        maximum_unknowns:4096,maximum_cache_bytes:256*1024*1024,
        maximum_work:per_frequency,angular_refinement },||gate.is_requested())?;
    let work=source.frequency_work(sources)?;
    if work>per_frequency {return Err("complete reciprocal kettle batch exceeds the declared sweep work".into());}
    eprintln!("reciprocal kettle: unknowns={}, cache_bytes={}, angular_refinement={angular_refinement}, sweep_work_units={}; no acoustic operator fallback",
        source.unknowns(),source.cache_bytes(),work*samples as u64);
    Ok(source)
}
