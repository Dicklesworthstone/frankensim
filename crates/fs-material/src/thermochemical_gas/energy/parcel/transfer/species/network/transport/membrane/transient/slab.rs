//! Finite sorbed inventory in an isothermal, planar, layered dense membrane.
//!
//! One permeant obeys c = S q, P = D S and S dq/dt = div(P grad q).
//! q is the gas-equivalent Henry pressure, not polymer stress or pore gas.
//! Finite-volume capacities are S dx; face conductances include both half-cell
//! resistances and any authored interlayer permeance. Backward Euler gives a
//! positive M-matrix for every representable positive time step. The original
//! tridiagonal solver owns elimination; this module checks original cell and
//! whole-membrane amount balances before returning a candidate.
//!
//! Material coefficients and geometry stay fixed. Reservoir pressures are
//! prescribed over each step; the reservoirs and an isothermal thermostat are
//! external. No finite gas energy update, sorption heat, cross-diffusion,
//! nonlinear sorption, swelling or concentration-dependent D is inferred.
//! P = D S and the initially degassed homogeneous-slab lag L^2/(6D):
//! <https://doi.org/10.3390/membranes8030073>.

use std::sync::Arc;
use fs_la::tridiagonal::{TridiagonalError, TridiagonalWorkspace};
use super::super::coefficients::{MembraneCoefficientError, MembraneDiffusivity,
    MolarHenrySolubility, MolarPermeability, MolarPermeance};

pub const MAX_TRANSIENT_MEMBRANE_LAYERS: usize = 64;
pub const MAX_TRANSIENT_MEMBRANE_CELLS: usize = 4096;

#[derive(Clone, Debug, PartialEq)]
pub enum MembraneDiffusionError {
    Invalid(&'static str),
    Coefficient(MembraneCoefficientError),
    Linear(TridiagonalError),
    AmountBalance { residual_mol: f64, tolerance_mol: f64 },
    Cancelled,
}
impl core::fmt::Display for MembraneDiffusionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "transient dense membrane: {self:?}")
    }
}
impl std::error::Error for MembraneDiffusionError {}
impl From<MembraneCoefficientError> for MembraneDiffusionError {
    fn from(error: MembraneCoefficientError) -> Self { Self::Coefficient(error) }
}
impl From<TridiagonalError> for MembraneDiffusionError {
    fn from(error: TridiagonalError) -> Self {
        if error == TridiagonalError::Cancelled { Self::Cancelled } else { Self::Linear(error) }
    }
}
use MembraneDiffusionError as Error;
fn positive(x: f64) -> bool { x.is_finite() && x > 0.0 }
fn nonnegative(x: f64) -> bool { x.is_finite() && x >= 0.0 }
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() { Err(Error::Cancelled) } else { Ok(()) }
}
fn finite(value: f64) -> Result<f64, Error> {
    if value.is_finite() { Ok(value) } else { Err(Error::Invalid("membrane arithmetic is unrepresentable")) }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MembraneDiffusionLayer {
    pub thickness_m: f64,
    pub cells: usize,
    pub diffusivity: MembraneDiffusivity,
    /// Finite-storage cells require positive S. A zero D still stores its
    /// authored initial amount and blocks transport exactly.
    pub solubility: MolarHenrySolubility,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MembraneDiffusionSolve {
    pub absolute_amount_tolerance_mol: f64,
    pub relative_amount_tolerance: f64,
}
impl MembraneDiffusionSolve {
    fn tolerance(self, scale: f64) -> Result<f64, Error> {
        finite(self.absolute_amount_tolerance_mol + self.relative_amount_tolerance * scale)
    }
    fn validate(self) -> Result<(), Error> {
        if !nonnegative(self.absolute_amount_tolerance_mol)
            || !nonnegative(self.relative_amount_tolerance) || self.relative_amount_tolerance >= 1.0
            || (self.absolute_amount_tolerance_mol == 0.0 && self.relative_amount_tolerance == 0.0) {
            return Err(Error::Invalid("a finite nonzero absolute/relative amount budget is required"));
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct Cell {
    layer: usize,
    width_m: f64,
    capacity_mol_m2_pa: f64,
    half_permeance: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct Face {
    permeance: f64,
    left_half: f64,
    right_half: f64,
    interface: Option<f64>,
}
fn series(values: &[f64]) -> Result<f64, Error> {
    if values.contains(&0.0) { return Ok(0.0); }
    let scale = values.iter().copied().fold(f64::INFINITY, f64::min);
    let mut sum = 0.0;
    for value in values {
        let ratio = scale / value;
        if !positive(ratio) { return Err(Error::Invalid("face resistance ratio is unrepresentable")); }
        sum += ratio;
    }
    let value = scale / sum;
    if !positive(value) { return Err(Error::Invalid("face permeance is unrepresentable")); }
    Ok(value)
}
impl Face {
    fn traces(self, left: f64, right: f64) -> [Option<f64>; 2] {
        if self.permeance > 0.0 {
            let a = self.permeance / self.left_half;
            let b = self.permeance / self.right_half;
            return [Some((1.0-a)*left+a*right), Some((1.0-b)*right+b*left)];
        }
        // A single seal leaves determinate traces on its open sides. Coincident
        // seals may leave a zero-volume trace undetermined; do not invent it.
        let open_interface = self.interface.is_none_or(|g| g > 0.0);
        let a = if self.left_half > 0.0 { Some(left) }
            else if open_interface && self.right_half > 0.0 { Some(right) } else { None };
        let b = if self.right_half > 0.0 { Some(right) }
            else if open_interface && self.left_half > 0.0 { Some(left) } else { None };
        [a, b]
    }
}

/// Immutable 1-D geometry and coefficients; one state cannot be applied to a
/// different grid because it retains this exact model through an Arc.
#[derive(Clone, Debug, PartialEq)]
pub struct MembraneDiffusion {
    area_m2: f64,
    layers: Vec<MembraneDiffusionLayer>,
    layer_ranges: Vec<std::ops::Range<usize>>,
    cells: Vec<Cell>,
    faces: Vec<Face>,
}
impl MembraneDiffusion {
    /// None between two layers means equilibrium contact. Some(0) is an exact
    /// interlayer seal, with finite inventories retained on both sides.
    pub fn new(area_m2: f64, layers: &[MembraneDiffusionLayer],
        interfaces: &[Option<MolarPermeance>], mut cancelled: impl FnMut() -> bool)
        -> Result<Self, Error>
    {
        poll(&mut cancelled)?;
        if !positive(area_m2) || layers.is_empty() || layers.len() > MAX_TRANSIENT_MEMBRANE_LAYERS
            || interfaces.len()+1 != layers.len() {
            return Err(Error::Invalid("positive area, bounded nonempty layers and one interface slot per layer pair required"));
        }
        let count = layers.iter().try_fold(0usize, |count, layer| count.checked_add(layer.cells));
        let Some(count) = count.filter(|&n| n > 0 && n <= MAX_TRANSIENT_MEMBRANE_CELLS) else {
            return Err(Error::Invalid("transient membrane cell count exceeds its bound"));
        };
        let mut cells = Vec::with_capacity(count);
        let mut layer_ranges = Vec::with_capacity(layers.len());
        for (index, layer) in layers.iter().enumerate() {
            poll(&mut cancelled)?;
            if !positive(layer.thickness_m) || layer.cells == 0 || !positive(layer.solubility.value_si()) {
                return Err(Error::Invalid("finite-storage layers need positive thickness, cell count and Henry solubility"));
            }
            let width_m = layer.thickness_m / layer.cells as f64;
            let capacity = layer.solubility.value_si() * width_m;
            let permeability = MolarPermeability::from_diffusion(layer.diffusivity,layer.solubility)?;
            let half = permeability.permeance(0.5*width_m)?.value_si();
            if !positive(capacity) || !positive(capacity*area_m2) {
                return Err(Error::Invalid("membrane storage capacity is unrepresentable"));
            }
            let start = cells.len();
            for local in 0..layer.cells {
                if local % 128 == 0 { poll(&mut cancelled)?; }
                cells.push(Cell { layer:index, width_m, capacity_mol_m2_pa:capacity, half_permeance:half });
            }
            layer_ranges.push(start..cells.len());
        }
        let mut faces = Vec::with_capacity(count+1);
        faces.push(Face { permeance:cells[0].half_permeance, left_half:0.0,
            right_half:cells[0].half_permeance, interface:None });
        for i in 1..count {
            if i % 128 == 0 { poll(&mut cancelled)?; }
            let (a,b) = (cells[i-1],cells[i]);
            let interface = if a.layer == b.layer { None }
                else { interfaces[a.layer].map(MolarPermeance::value_si) };
            let permeance = match interface {
                Some(g) => series(&[a.half_permeance,g,b.half_permeance])?,
                None => series(&[a.half_permeance,b.half_permeance])?,
            };
            faces.push(Face {permeance,left_half:a.half_permeance,right_half:b.half_permeance,interface});
        }
        let half = cells[count-1].half_permeance;
        faces.push(Face {permeance:half,left_half:half,right_half:0.0,interface:None});
        poll(&mut cancelled)?;
        Ok(Self {area_m2,layers:layers.to_vec(),layer_ranges,cells,faces})
    }
    pub const fn area_m2(&self) -> f64 { self.area_m2 }
    pub fn layers(&self) -> &[MembraneDiffusionLayer] { &self.layers }
    pub fn layer_cell_ranges(&self) -> &[std::ops::Range<usize>] { &self.layer_ranges }
    pub fn cell_count(&self) -> usize { self.cells.len() }
    pub fn face_permeances(&self) -> impl ExactSizeIterator<Item=f64> + '_ {
        self.faces.iter().map(|face|face.permeance)
    }
    /// Initial loading is explicit per-cell Henry pressure q=c/S. Zero is an
    /// initially degassed specimen, not an instantaneous equilibrium state.
    pub fn initial_state(self, partial_pressures_pa: &[f64],
        mut cancelled: impl FnMut() -> bool) -> Result<MembraneDiffusionState,Error>
    {
        poll(&mut cancelled)?;
        if partial_pressures_pa.len()!=self.cells.len() {
            return Err(Error::Invalid("initial membrane loading must cover every cell"));
        }
        for (i,&q) in partial_pressures_pa.iter().enumerate() {
            if i%128==0 {poll(&mut cancelled)?;}
            if !nonnegative(q) {
                return Err(Error::Invalid("initial Henry pressures must be nonnegative and representable"));
            }
            finite(self.layers[self.cells[i].layer].solubility.value_si()*q)?;
            finite(self.cells[i].capacity_mol_m2_pa*self.area_m2*q)?;
        }
        let state=MembraneDiffusionState {model:Arc::new(self),partial_pressures_pa:partial_pressures_pa.to_vec(),elapsed_s:0.0};
        state.stored_amount_mol()?;
        poll(&mut cancelled)?;
        Ok(state)
    }
}

/// Retained linear factors and scratch. Repeated equal-dt steps on the same
/// immutable grid reuse the factorization; failed proposals change no state.
#[derive(Debug)]
pub struct MembraneDiffusionWorkspace {
    linear: TridiagonalWorkspace,
    mass: Vec<f64>,
    left: Vec<f64>,
    right: Vec<f64>,
    diagonal: Vec<f64>,
    lower: Vec<f64>,
    upper: Vec<f64>,
    rhs: Vec<f64>,
    solution: Vec<f64>,
    cached: Option<(Arc<MembraneDiffusion>,u64)>,
}
impl MembraneDiffusionWorkspace {
    pub fn new(cells:usize)->Result<Self,Error> {
        if cells==0 || cells>MAX_TRANSIENT_MEMBRANE_CELLS {return Err(Error::Invalid("bounded nonempty membrane workspace required"));}
        Ok(Self {linear:TridiagonalWorkspace::new(cells)?,mass:vec![0.0;cells],left:vec![0.0;cells],
            right:vec![0.0;cells],diagonal:vec![0.0;cells],lower:vec![0.0;cells-1],upper:vec![0.0;cells-1],
            rhs:vec![0.0;cells],solution:vec![0.0;cells],cached:None})
    }
    fn factor(&mut self,model:&Arc<MembraneDiffusion>,dt:f64,cancelled:&mut impl FnMut()->bool)->Result<(),Error> {
        if self.mass.len()!=model.cells.len() {return Err(Error::Invalid("membrane workspace/grid dimensions differ"));}
        if self.cached.as_ref().is_some_and(|(old,bits)|Arc::ptr_eq(old,model)&&*bits==dt.to_bits()) {return Ok(());}
        self.cached=None;
        for (i,cell) in model.cells.iter().enumerate() {
            if i%128==0 {poll(cancelled)?;}
            let l=finite(dt*model.faces[i].permeance)?;
            let r=finite(dt*model.faces[i+1].permeance)?;
            if (model.faces[i].permeance>0.0&&l==0.0)||(model.faces[i+1].permeance>0.0&&r==0.0) {
                return Err(Error::Invalid("integrated membrane face conductance underflowed"));
            }
            let scale=cell.capacity_mol_m2_pa.max(l).max(r);
            let (mass,left,right)=(cell.capacity_mol_m2_pa/scale,l/scale,r/scale);
            let diagonal=mass+left+right;
            if !positive(mass)||diagonal<=left+right||(l>0.0&&left==0.0)||(r>0.0&&right==0.0) {
                return Err(Error::Invalid("time step cannot retain positive membrane storage in floating-point arithmetic"));
            }
            self.mass[i]=mass;self.left[i]=left;self.right[i]=right;self.diagonal[i]=diagonal;
            if i>0 {self.lower[i-1]=-left;}
            if i+1<model.cells.len() {self.upper[i]=-right;}
        }
        self.linear.factor(&self.lower,&self.diagonal,&self.upper,&mut *cancelled)?;
        poll(cancelled)?;
        self.cached=Some((model.clone(),dt.to_bits()));
        Ok(())
    }
}

#[derive(Clone,Debug,PartialEq)]
pub struct MembraneDiffusionState {
    model:Arc<MembraneDiffusion>,
    partial_pressures_pa:Vec<f64>,
    elapsed_s:f64,
}
#[derive(Clone,Debug,PartialEq)]
pub struct MembraneDiffusionFlux {
    pub reservoir_partial_pressures_pa:[f64;2],
    pub face_molar_fluxes_mol_m2_s:Vec<f64>,
    pub face_partial_pressures_pa:Vec<[Option<f64>;2]>,
}
#[derive(Clone,Debug,PartialEq)]
pub struct MembraneDiffusionStep {
    pub duration_s:f64,
    pub reservoir_partial_pressures_pa:[f64;2],
    /// Each face has one signed left-to-right flux, shared by neighboring cells.
    pub face_molar_fluxes_mol_m2_s:Vec<f64>,
    /// The two sides of each internal face; None is an undetermined inactive
    /// zero-volume trace between coincident seals, never a guessed pressure.
    pub face_partial_pressures_pa:Vec<[Option<f64>;2]>,
    pub stored_before_mol:f64,
    pub stored_after_mol:f64,
    /// Sum of cell changes avoids subtracting two large stored totals.
    pub stored_change_mol:f64,
    pub left_reservoir_to_membrane_mol:f64,
    pub membrane_to_right_reservoir_mol:f64,
    pub amount_residual_mol:f64,
    pub maximum_cell_amount_residual_mol:f64,
}
#[derive(Clone,Debug,PartialEq)]
pub struct MembraneDiffusionCandidate {
    state:MembraneDiffusionState,
    step:MembraneDiffusionStep,
}
impl MembraneDiffusionCandidate {
    pub const fn state(&self)->&MembraneDiffusionState {&self.state}
    pub const fn step(&self)->&MembraneDiffusionStep {&self.step}
    pub fn into_parts(self)->(MembraneDiffusionState,MembraneDiffusionStep) {(self.state,self.step)}
}
#[derive(Default)]
struct Sum {value:f64,correction:f64}
impl Sum {
    fn add(&mut self,x:f64)->Result<(),Error> {
        let y=finite(x-self.correction)?;let sum=finite(self.value+y)?;
        self.correction=finite((sum-self.value)-y)?;self.value=sum;Ok(())
    }
}
impl MembraneDiffusionState {
    pub fn model(&self)->&MembraneDiffusion {&self.model}
    pub fn partial_pressures_pa(&self)->&[f64] {&self.partial_pressures_pa}
    pub const fn elapsed_s(&self)->f64 {self.elapsed_s}
    pub fn concentrations_mol_m3(&self)->impl ExactSizeIterator<Item=f64>+'_ {
        self.partial_pressures_pa.iter().zip(&self.model.cells)
            .map(|(&q,cell)|q*self.model.layers[cell.layer].solubility.value_si())
    }
    pub fn stored_amount_mol(&self)->Result<f64,Error> {
        let mut amount=Sum::default();
        for (&q,cell) in self.partial_pressures_pa.iter().zip(&self.model.cells) {
            amount.add((cell.capacity_mol_m2_pa*self.model.area_m2)*q)?;
        }
        Ok(amount.value)
    }
    /// Instantaneous resolved-grid fluxes and exact interlayer Henry traces.
    /// Sampling does not equilibrate or change the sorbed inventory.
    pub fn sample_fluxes(&self,reservoirs:[f64;2],mut cancelled:impl FnMut()->bool)
        ->Result<MembraneDiffusionFlux,Error>
    {
        poll(&mut cancelled)?;
        if !reservoirs.into_iter().all(nonnegative) {
            return Err(Error::Invalid("nonnegative finite reservoir pressures required"));
        }
        let n=self.partial_pressures_pa.len();
        let mut fluxes=Vec::with_capacity(n+1);let mut traces=Vec::with_capacity(n+1);
        for (i,face) in self.model.faces.iter().enumerate() {
            if i%128==0 {poll(&mut cancelled)?;}
            let left=if i==0 {reservoirs[0]}else{self.partial_pressures_pa[i-1]};
            let right=if i==n {reservoirs[1]}else{self.partial_pressures_pa[i]};
            let flux=finite(face.permeance*(left-right))?;
            if face.permeance>0.0&&left!=right&&flux==0.0 {
                return Err(Error::Invalid("membrane molar flux underflowed"));
            }
            fluxes.push(flux);
            traces.push(if i==0 {[Some(reservoirs[0]);2]}else if i==n {[Some(reservoirs[1]);2]}
                else {face.traces(left,right)});
        }
        poll(&mut cancelled)?;
        Ok(MembraneDiffusionFlux {reservoir_partial_pressures_pa:reservoirs,
            face_molar_fluxes_mol_m2_s:fluxes,face_partial_pressures_pa:traces})
    }
    /// First-order backward Euler with fixed reservoir pressures throughout
    /// the step. Stability does not imply resolved breakthrough: refine time
    /// and cells for accuracy. Source-domain admission belongs to the binder.
    pub fn step_candidate(&self,reservoirs:[f64;2],duration_s:f64,controls:MembraneDiffusionSolve,
        work:&mut MembraneDiffusionWorkspace,mut cancelled:impl FnMut()->bool)
        ->Result<MembraneDiffusionCandidate,Error>
    {
        poll(&mut cancelled)?;controls.validate()?;
        let elapsed_s=self.elapsed_s+duration_s;
        if !positive(duration_s)||!reservoirs.into_iter().all(nonnegative)||!elapsed_s.is_finite()||elapsed_s<=self.elapsed_s {
            return Err(Error::Invalid("positive representable duration and nonnegative reservoir pressures required"));
        }
        work.factor(&self.model,duration_s,&mut cancelled)?;
        let n=self.partial_pressures_pa.len();
        for i in 0..n {
            if i%128==0 {poll(&mut cancelled)?;}
            let mut rhs=work.mass[i]*self.partial_pressures_pa[i];
            if i==0 {rhs=work.left[i].mul_add(reservoirs[0],rhs);}
            if i+1==n {rhs=work.right[i].mul_add(reservoirs[1],rhs);}
            work.rhs[i]=finite(rhs)?;
        }
        work.linear.solve_factored_into(&work.rhs,&mut work.solution,&mut cancelled)?;
        for (i,&q) in work.solution.iter().enumerate() {
            if i%128==0 {poll(&mut cancelled)?;}
            if !nonnegative(q) {
                return Err(Error::Invalid("implicit membrane solution lost nonnegative sorbed amount"));
            }
            finite(self.model.layers[self.model.cells[i].layer].solubility.value_si()*q)?;
        }
        let state=Self {model:self.model.clone(),partial_pressures_pa:work.solution.clone(),elapsed_s};
        let sampled=state.sample_fluxes(reservoirs,&mut cancelled)?;
        let mut face_amount=Vec::with_capacity(n+1);
        for (i,&flux) in sampled.face_molar_fluxes_mol_m2_s.iter().enumerate() {
            if i%128==0 {poll(&mut cancelled)?;}
            let amount=finite((flux*self.model.area_m2)*duration_s)?;
            if flux!=0.0&&amount==0.0 {
                return Err(Error::Invalid("membrane transfer amount underflowed"));
            }
            face_amount.push(amount);
        }
        let mut change=Sum::default();let mut maximum:f64=0.0;
        for (i,cell) in self.model.cells.iter().enumerate() {
            if i%128==0 {poll(&mut cancelled)?;}
            let delta=finite((cell.capacity_mol_m2_pa*self.model.area_m2)
                *(work.solution[i]-self.partial_pressures_pa[i]))?;
            change.add(delta)?;
            let residual=finite(delta-(face_amount[i]-face_amount[i+1]))?;
            let tolerance=controls.tolerance(delta.abs()+face_amount[i].abs()+face_amount[i+1].abs())?;
            if residual.abs()>tolerance {return Err(Error::AmountBalance {residual_mol:residual,tolerance_mol:tolerance});}
            maximum=maximum.max(residual.abs());
        }
        let residual=finite(change.value-(face_amount[0]-face_amount[n]))?;
        let tolerance=controls.tolerance(change.value.abs()+face_amount[0].abs()+face_amount[n].abs())?;
        if residual.abs()>tolerance {return Err(Error::AmountBalance {residual_mol:residual,tolerance_mol:tolerance});}
        let step=MembraneDiffusionStep {duration_s,reservoir_partial_pressures_pa:reservoirs,
            face_molar_fluxes_mol_m2_s:sampled.face_molar_fluxes_mol_m2_s,
            face_partial_pressures_pa:sampled.face_partial_pressures_pa,
            stored_before_mol:self.stored_amount_mol()?,stored_after_mol:state.stored_amount_mol()?,stored_change_mol:change.value,
            left_reservoir_to_membrane_mol:face_amount[0],membrane_to_right_reservoir_mol:face_amount[n],
            amount_residual_mol:residual,maximum_cell_amount_residual_mol:maximum};
        poll(&mut cancelled)?;
        Ok(MembraneDiffusionCandidate {state,step})
    }
}

#[cfg(test)]
#[path="slab/tests.rs"]
mod tests;
