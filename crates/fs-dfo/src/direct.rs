//! Deterministic DIRECT box search (Jones, Perttunen and Stuckman, 1993).
//! Select every rectangle admitting a positive Lipschitz slope in the original
//! potentially-optimal inequalities, then trisect all its longest sides in
//! order of their best sampled values. No Lipschitz constant is supplied or
//! estimated. See Jones and Martins (2021), doi:10.1007/s10898-020-00952-6.
//!
//! A finite run returns observations, NOT a certified global lower bound or
//! an optimum. Selection is ordinary binary64 with explicit bounded pair work.
//! Successful callbacks survive pauses/cancellation; a failed callback is a
//! typed error, never a fabricated penalty. The accepted partition changes only
//! after all samples for a split exist. Keep the objective fixed on continuation.

use fs_exec::Cx;
use std::{collections::VecDeque, fmt};

#[derive(Debug, Clone, PartialEq)]
pub struct DirectConfig {
    /// Lifetime callback attempts, including failures. Zero is allowed.
    pub max_evaluations: usize,
    pub max_rectangles: usize,
    /// Rectangle count times dimension; not a byte or total-heap bound.
    pub max_rectangle_entries: usize,
    /// Maximum rectangle pairs inspected in ONE selection round.
    pub max_selection_comparisons: usize,
    /// Relative improvement parameter in (0,1], commonly 1e-4.
    pub epsilon: f64,
    /// Caller-supplied finite goal, not an inferred optimum.
    pub target: Option<f64>,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DirectError { Invalid(&'static str), NonFinite, Allocation, Cancelled }
impl fmt::Display for DirectError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(s) => f.write_str(s),
            Self::NonFinite => f.write_str("DIRECT arithmetic or objective is not finite"),
            Self::Allocation => f.write_str("DIRECT allocation refused"),
            Self::Cancelled => f.write_str("DIRECT cancelled; completed samples retained"),
        }
    }
}
impl std::error::Error for DirectError {}
#[derive(Debug)]
pub enum DirectRunError<E> { Solver(DirectError), Evaluation(E) }
impl<E> From<DirectError> for DirectRunError<E> {
    fn from(e: DirectError) -> Self { Self::Solver(e) }
}
impl<E: fmt::Display> fmt::Display for DirectRunError<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self { Self::Solver(e) => write!(f,"{e}"), Self::Evaluation(e) => write!(f,"DIRECT objective: {e}") }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for DirectRunError<E> {}
#[derive(Debug, Clone, PartialEq)]
pub struct DirectPoint { pub coordinates: Vec<f64>, pub value: f64 }
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DirectStop {
    TargetReached, EvaluationLimit, Paused, RectangleLimit, SelectionLimit,
    /// A child center cannot be distinguished in normalized or physical coordinates.
    ResolutionLimit,
}
#[derive(Debug, Clone, PartialEq)]
pub struct DirectReport {
    pub stop: DirectStop,
    /// Best consumed sample, including a partially sampled split.
    pub best: Option<DirectPoint>,
    pub evaluations: usize,
    pub iterations: usize,
    pub rectangles: usize,
    pub pending_samples: usize,
}
#[derive(Debug)]
struct Rectangle { center: Vec<f64>, half: Vec<f64>, radius: f64, value: f64 }
#[derive(Debug)]
struct Sample { center: Vec<f64>, physical: Vec<f64>, value: Option<f64> }
#[derive(Debug)]
struct Split { parent: usize, axes: Vec<usize>, samples: Vec<Sample>, consumed: usize }

/// Serial, resumable in-process DIRECT. There is no stochastic seed, derivative
/// tape, nonlinear constraint adapter, checkpoint codec or parallel callback.
/// Parents and split queues are private so callers cannot refund spent work.
#[derive(Debug)]
pub struct DirectState {
    bounds: Vec<[f64;2]>, config: DirectConfig,
    root: Option<Sample>, rectangles: Vec<Rectangle>, selected: VecDeque<usize>,
    split: Option<Split>, best: Option<DirectPoint>, evaluations: usize, iterations: usize,
}
impl DirectState {
    /// Finite, strictly ordered bounds with finite widths; fixed variables must
    /// be removed upstream. Constructor invokes no callback.
    pub fn new(cx: &Cx<'_>, bounds: &[[f64;2]], config: DirectConfig) -> Result<Self, DirectError> {
        Self::new_with(bounds, config, &mut || cx.checkpoint().is_ok())
    }
    fn new_with(bounds: &[[f64;2]], config: DirectConfig, poll: &mut impl FnMut()->bool)
        -> Result<Self, DirectError>
    {
        check(poll)?;
        if bounds.is_empty() || config.max_rectangles == 0 || bounds.len() > config.max_rectangle_entries
            || config.max_selection_comparisons == 0 || !config.epsilon.is_finite()
            || config.epsilon <= 0.0 || config.epsilon > 1.0
            || config.target.is_some_and(|v| !v.is_finite())
        { return Err(DirectError::Invalid("DIRECT needs a nonempty bounded shape, positive selection cap and epsilon in (0,1]")); }
        for (i,&[lo,hi]) in bounds.iter().enumerate() {
            if i%128==0 {check(poll)?;}
            if !lo.is_finite() || !hi.is_finite() || lo >= hi || !(hi-lo).is_finite() {
                return Err(DirectError::Invalid("DIRECT needs finite ordered bounds with representable widths"));
            }
        }
        let mut stored=reserve(bounds.len())?; stored.extend_from_slice(bounds);
        let mut center=reserve(bounds.len())?; center.resize(bounds.len(),0.5);
        let physical=map(&stored,&center,poll)?;
        Ok(Self { bounds:stored, config, root:Some(Sample {center,physical,value:None}),
            rectangles:Vec::new(), selected:VecDeque::new(), split:None, best:None, evaluations:0, iterations:0 })
    }
    pub const fn evaluations(&self)->usize {self.evaluations}
    pub const fn evaluation_limit(&self)->usize {self.config.max_evaluations}
    pub fn best(&self)->Option<&DirectPoint> {self.best.as_ref()}
    pub fn extend_evaluation_budget(&mut self,total:usize)->Result<(),DirectError> {
        if total<=self.config.max_evaluations {return Err(DirectError::Invalid("DIRECT evaluation allowance must increase"));}
        self.config.max_evaluations=total;Ok(())
    }
    /// Increase partition/selection caps without discarding completed samples.
    pub fn extend_partition_budget(&mut self,rectangles:usize,entries:usize,pairs:usize)->Result<(),DirectError> {
        if rectangles<self.config.max_rectangles || entries<self.config.max_rectangle_entries
            || pairs<self.config.max_selection_comparisons
            || (rectangles==self.config.max_rectangles && entries==self.config.max_rectangle_entries
                && pairs==self.config.max_selection_comparisons)
        {return Err(DirectError::Invalid("DIRECT partition allowances may only increase"));}
        self.config.max_rectangles=rectangles;self.config.max_rectangle_entries=entries;
        self.config.max_selection_comparisons=pairs;Ok(())
    }
    pub fn advance<E>(&mut self,cx:&Cx<'_>,additional_evaluations:usize,
        evaluate:&mut impl FnMut(&[f64])->Result<f64,E>)->Result<DirectReport,DirectRunError<E>>
    {self.advance_with(additional_evaluations,evaluate,&mut ||cx.checkpoint().is_ok())}

    fn advance_with<E>(&mut self,additional:usize,evaluate:&mut impl FnMut(&[f64])->Result<f64,E>,
        poll:&mut impl FnMut()->bool)->Result<DirectReport,DirectRunError<E>>
    {
        let mut spent=0usize;
        loop {
            check(poll)?;
            if self.root.as_ref().is_some_and(|s|s.value.is_some()) {self.consume_root(poll)?;continue;}
            if let Some(split)=&self.split {
                if split.consumed==split.samples.len() {self.publish_split(poll)?;continue;}
                if split.samples[split.consumed].value.is_some() {self.consume_sample(poll)?;continue;}
            }
            let stop=if self.best.as_ref().is_some_and(|p|self.config.target.is_some_and(|t|p.value<=t)) {
                Some(DirectStop::TargetReached)
            } else if self.evaluations>=self.config.max_evaluations {Some(DirectStop::EvaluationLimit)}
            else if spent>=additional {Some(DirectStop::Paused)} else {None};
            if let Some(stop)=stop {return Ok(self.report(stop,poll)?);}
            if self.root.is_none() && self.split.is_none() {
                if self.selected.is_empty() {
                    if self.rectangles.len().checked_mul(self.rectangles.len())
                        .is_none_or(|n|n>self.config.max_selection_comparisons)
                    {return Ok(self.report(DirectStop::SelectionLimit,poll)?);}
                    self.selected=select(&self.rectangles,self.best.as_ref().expect("scored root").value,
                        self.config.epsilon,poll)?;
                    if self.selected.is_empty() {return Err(DirectError::Invalid("DIRECT selection found no representable positive slope").into());}
                }
                if let Some(stop)=self.prepare_split(poll)? {return Ok(self.report(stop,poll)?);}
            }
            check(poll)?;
            let sample=match &mut self.root {
                Some(root)=>root,
                None=>{let split=self.split.as_mut().expect("prepared split");&mut split.samples[split.consumed]},
            };
            self.evaluations+=1;spent+=1;
            let value=evaluate(&sample.physical).map_err(DirectRunError::Evaluation)?;
            // Cache before cancellation can be observed; retries never repeat a
            // finite successful callback, even in the middle of a split.
            sample.value=Some(finite(value)?);
        }
    }
    fn consume_root(&mut self,poll:&mut impl FnMut()->bool)->Result<(),DirectError> {
        let s=self.root.as_ref().expect("root");
        let value=s.value.expect("scored root");
        let half=copy(&s.center,poll)?;
        let root=Rectangle {center:copy(&s.center,poll)?,radius:radius(&half,poll)?,half,value};
        let best=DirectPoint {coordinates:copy(&s.physical,poll)?,value};
        self.rectangles.try_reserve(1).map_err(|_|DirectError::Allocation)?;
        check(poll)?;
        self.rectangles.push(root);self.best=Some(best);self.root=None;Ok(())
    }
    fn consume_sample(&mut self,poll:&mut impl FnMut()->bool)->Result<(),DirectError> {
        let split=self.split.as_ref().expect("split");let s=&split.samples[split.consumed];
        let value=s.value.expect("scored sample");
        let best=if self.best.as_ref().is_none_or(|p|value.total_cmp(&p.value).is_lt()) {
            Some(DirectPoint {coordinates:copy(&s.physical,poll)?,value})
        } else {None};
        check(poll)?;
        if let Some(best)=best {self.best=Some(best);}
        self.split.as_mut().expect("split").consumed+=1;Ok(())
    }
    fn prepare_split(&mut self,poll:&mut impl FnMut()->bool)->Result<Option<DirectStop>,DirectError> {
        let parent=*self.selected.front().expect("selected rectangle");let r=&self.rectangles[parent];
        let largest=r.half.iter().copied().fold(0.0,f64::max);
        let mut axes=reserve(self.bounds.len())?;
        for (axis,&half) in r.half.iter().enumerate() {if axis%128==0{check(poll)?;} if half==largest {axes.push(axis);}}
        let added=axes.len().checked_mul(2).ok_or(DirectError::Allocation)?;
        if self.rectangles.len().checked_add(added).is_none_or(|n|n>self.config.max_rectangles
            || n.checked_mul(self.bounds.len()).is_none_or(|m|m>self.config.max_rectangle_entries))
        {return Ok(Some(DirectStop::RectangleLimit));}
        let delta=(2.0/3.0)*largest;
        let base=map(&self.bounds,&r.center,poll)?;
        let mut samples=reserve(added)?;
        for &axis in &axes {
            check(poll)?;
            for sign in [-1.0,1.0] {
                let mut center=copy(&r.center,poll)?;center[axis]+=sign*delta;
                let physical=map(&self.bounds,&center,poll)?;
                if !(0.0..=1.0).contains(&center[axis]) || center[axis]==r.center[axis]
                    || physical[axis]==base[axis]
                {return Ok(Some(DirectStop::ResolutionLimit));}
                samples.push(Sample {center,physical,value:None});
            }
        }
        check(poll)?;self.split=Some(Split {parent,axes,samples,consumed:0});Ok(None)
    }
    fn publish_split(&mut self,poll:&mut impl FnMut()->bool)->Result<(),DirectError> {
        let split=self.split.as_ref().expect("complete split");let parent=&self.rectangles[split.parent];
        let mut order=reserve(split.axes.len())?;order.extend(0..split.axes.len());
        order.sort_by(|&a,&b| {
            let score=|i:usize|split.samples[2*i].value.unwrap().min(split.samples[2*i+1].value.unwrap());
            score(a).total_cmp(&score(b)).then_with(||split.axes[a].cmp(&split.axes[b]))
        });
        let mut half=copy(&parent.half,poll)?;let mut children=reserve(split.samples.len())?;
        for i in order {
            check(poll)?;half[split.axes[i]]/=3.0;
            for s in &split.samples[2*i..2*i+2] {
                children.push(Rectangle {center:copy(&s.center,poll)?,half:copy(&half,poll)?,
                    radius:radius(&half,poll)?,value:s.value.expect("complete split")});
            }
        }
        let middle=Rectangle {center:copy(&parent.center,poll)?,radius:radius(&half,poll)?,half,value:parent.value};
        let iteration=if self.selected.len()==1 {self.iterations.checked_add(1)
            .ok_or(DirectError::Invalid("DIRECT iteration counter overflow"))?} else {self.iterations};
        self.rectangles.try_reserve(children.len()).map_err(|_|DirectError::Allocation)?;
        check(poll)?;
        let parent=split.parent;
        self.rectangles[parent]=middle;self.rectangles.extend(children);
        self.selected.pop_front();self.split=None;self.iterations=iteration;Ok(())
    }
    fn report(&self,stop:DirectStop,poll:&mut impl FnMut()->bool)->Result<DirectReport,DirectError> {
        let best=self.best.as_ref().map(|p|Ok::<_,DirectError>(DirectPoint {
            coordinates:copy(&p.coordinates,poll)?,value:p.value
        })).transpose()?;
        check(poll)?;
        Ok(DirectReport {stop,best,evaluations:self.evaluations,iterations:self.iterations,
            rectangles:self.rectangles.len(),pending_samples:self.root.as_ref().map_or_else(
                ||self.split.as_ref().map_or(0,|s|s.samples.len()-s.consumed),|_|1)})
    }
}

// Definition 1 in Jones/Martins: intersect admissible K intervals, including
// the epsilon-improvement inequality. Equal-size ties all remain candidates.
fn select(rectangles:&[Rectangle],best:f64,epsilon:f64,poll:&mut impl FnMut()->bool)
    ->Result<VecDeque<usize>,DirectError>
{
    let scale=rectangles.iter().map(|r|r.value.abs()).fold(1.0,f64::max);
    let goal=best/scale-epsilon*(best/scale).abs();
    let mut selected=VecDeque::new();selected.try_reserve(rectangles.len()).map_err(|_|DirectError::Allocation)?;
    for (j,r) in rectangles.iter().enumerate() {
        check(poll)?;let fj=r.value/scale;
        let mut lower=finite((fj-goal)/r.radius)?.max(0.0);let mut upper=f64::INFINITY;let mut possible=true;
        for (i,other) in rectangles.iter().enumerate() {
            if i%128==0{check(poll)?;}
            if other.radius==r.radius {if other.value<r.value {possible=false;break;}continue;}
            let slope=finite((fj-other.value/scale)/(r.radius-other.radius))?;
            if r.radius>other.radius {lower=lower.max(slope);} else {upper=upper.min(slope);}
            if lower>upper {possible=false;break;}
        }
        if possible && upper>0.0 && lower<=upper {selected.push_back(j);}
    }
    check(poll)?;Ok(selected)
}
fn radius(half:&[f64],poll:&mut impl FnMut()->bool)->Result<f64,DirectError> {
    // Permutations of the same side lengths need identical size classes.
    let mut sides=copy(half,poll)?;sides.sort_by(f64::total_cmp);
    let mut sum=0.0;for (i,x) in sides.into_iter().enumerate(){if i%128==0{check(poll)?;}sum+=x*x;}
    let result=finite(fs_math::det::sqrt(sum))?;
    if result<=0.0 {return Err(DirectError::Invalid("DIRECT rectangle radius underflow"));}Ok(result)
}
fn map(bounds:&[[f64;2]],unit:&[f64],poll:&mut impl FnMut()->bool)->Result<Vec<f64>,DirectError> {
    let mut out=reserve(bounds.len())?;
    for (i,(&[lo,hi],&u)) in bounds.iter().zip(unit).enumerate(){
        if i%128==0{check(poll)?;}
        let x=finite(u.mul_add(hi-lo,lo))?;
        if x<lo || x>hi{return Err(DirectError::Invalid("DIRECT center left its declared box"));}out.push(x);
    }Ok(out)
}
fn reserve<T>(n:usize)->Result<Vec<T>,DirectError>{let mut v=Vec::new();v.try_reserve_exact(n).map_err(|_|DirectError::Allocation)?;Ok(v)}
fn copy(v:&[f64],poll:&mut impl FnMut()->bool)->Result<Vec<f64>,DirectError>{
    let mut out=reserve(v.len())?;for (i,&x) in v.iter().enumerate(){if i%128==0{check(poll)?;}out.push(x);}Ok(out)
}
fn finite(x:f64)->Result<f64,DirectError>{if x.is_finite(){Ok(x)}else{Err(DirectError::NonFinite)}}
fn check(poll:&mut impl FnMut()->bool)->Result<(),DirectError>{if poll(){Ok(())}else{Err(DirectError::Cancelled)}}

#[cfg(test)]
mod tests;
