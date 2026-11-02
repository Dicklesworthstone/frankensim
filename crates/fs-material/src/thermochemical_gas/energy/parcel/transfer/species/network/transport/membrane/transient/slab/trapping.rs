//! Mobile Henry diffusion coupled to finite, immobile Langmuir site kinetics.
//!
//! Only the mobile concentration `c_H=S*q` diffuses with the original `D*S`
//! face conductances. Each cell additionally stores `n_T`, with
//! `dn_T/dt = k_d*b*q*(C_T-n_T)-k_d*n_T`. The original shared site owner
//! eliminates its backward-Euler endpoint analytically. Its positive analytic
//! derivative leaves one nonlinear tridiagonal pressure system, solved by the
//! existing `fs_la` factorization. There is no splitting of diffusion and uptake.
//!
//! Saturation, finite uptake/desorption and retained loading history therefore
//! change breakthrough without inventing a pressure-dependent diffusivity.
//! Temperature and coefficients remain fixed; traps do not move, compete across
//! species, swell the polymer or supply an inferred caloric/heat closure.

use super::*;
use super::finite::MolarPressureCapacity;
use crate::adsorption::competitive::sites::scalar::{
    ScalarSiteError, ScalarSiteStep, backward_euler_site,
};

/// Trap capacity is per whole membrane volume, distinct from mobile S*q.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LangmuirMembraneTrap {
    capacity_mol_m3: f64,
    affinity_pa_inverse: f64,
    desorption_rate_s_inverse: f64,
    adsorption_rate_pa_inverse_s_inverse: f64,
}
impl LangmuirMembraneTrap {
    pub fn try_new(capacity_mol_m3: f64, affinity_pa_inverse: f64,
        desorption_rate_s_inverse: f64) -> Result<Self, TrappingMembraneError> {
        let rate = desorption_rate_s_inverse * affinity_pa_inverse;
        if !positive(capacity_mol_m3) || !positive(affinity_pa_inverse)
            || !nonnegative(desorption_rate_s_inverse) || !nonnegative(rate)
            || (desorption_rate_s_inverse > 0.0 && rate == 0.0) {
            return Err(TrappingMembraneError::Invalid("positive trap capacity/affinity and representable nonnegative kinetics required"));
        }
        Ok(Self { capacity_mol_m3, affinity_pa_inverse, desorption_rate_s_inverse,
            adsorption_rate_pa_inverse_s_inverse: rate })
    }
    pub const fn capacity_mol_m3(self) -> f64 { self.capacity_mol_m3 }
    pub const fn affinity_pa_inverse(self) -> f64 { self.affinity_pa_inverse }
    pub const fn desorption_rate_s_inverse(self) -> f64 { self.desorption_rate_s_inverse }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TrappingMembraneSolve {
    /// Original molar balance budget, used in every nonlinear row and globally.
    pub balance: MembraneDiffusionSolve,
    pub maximum_iterations: usize,
    pub maximum_backtracks: usize,
}
impl Default for TrappingMembraneSolve {
    fn default() -> Self {
        Self { balance: MembraneDiffusionSolve { absolute_amount_tolerance_mol: 1e-14,
            relative_amount_tolerance: 1e-10 }, maximum_iterations: 64, maximum_backtracks: 32 }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum TrappingMembraneError {
    Invalid(&'static str),
    Diffusion(MembraneDiffusionError),
    Sites(ScalarSiteError),
    NotConverged { iterations: usize, residual_ratio: f64 },
    Cancelled,
}
impl core::fmt::Display for TrappingMembraneError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "kinetic membrane traps: {self:?}")
    }
}
impl std::error::Error for TrappingMembraneError {}
impl From<MembraneDiffusionError> for TrappingMembraneError {
    fn from(error: MembraneDiffusionError) -> Self {
        if error == MembraneDiffusionError::Cancelled { Self::Cancelled }
        else { Self::Diffusion(error) }
    }
}
impl From<TridiagonalError> for TrappingMembraneError {
    fn from(error: TridiagonalError) -> Self { MembraneDiffusionError::from(error).into() }
}
impl From<ScalarSiteError> for TrappingMembraneError {
    fn from(error: ScalarSiteError) -> Self { Self::Sites(error) }
}
type Result<T, E = TrappingMembraneError> = core::result::Result<T, E>;

#[derive(Clone, Debug, PartialEq)]
pub struct TrappingMembraneDiffusionState {
    mobile: MembraneDiffusionState,
    traps: Arc<[Option<LangmuirMembraneTrap>]>,
    trapped_moles: Vec<f64>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct TrappingReservoirBalance {
    pub pressures_before_pa: [f64; 2],
    pub pressures_after_pa: [f64; 2],
    pub amount_changes_mol: [f64; 2],
    pub amount_residuals_mol: [f64; 2],
    pub closed_amount_residual_mol: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct TrappingMembraneStep {
    pub duration_s: f64,
    /// Original mobile-phase face law, including ordered interface traces.
    pub transport: MembraneDiffusionFlux,
    pub mobile_before_mol: f64,
    pub mobile_after_mol: f64,
    pub mobile_change_mol: f64,
    pub trapped_before_mol: f64,
    pub trapped_after_mol: f64,
    pub trapped_change_mol: f64,
    pub gross_adsorption_mol: Vec<f64>,
    pub gross_desorption_mol: Vec<f64>,
    /// [Left reservoir into membrane, membrane into right reservoir].
    pub boundary_transfers_mol: [f64; 2],
    /// Mobile plus trapped change minus the net boundary transfer.
    pub membrane_amount_residual_mol: f64,
    /// Maximum stored-amount or independently evaluated site-kinetic residual.
    pub maximum_cell_amount_residual_mol: f64,
    pub finite_reservoirs: Option<TrappingReservoirBalance>,
    pub nonlinear_iterations: usize,
    pub residual_evaluations: usize,
}
#[derive(Clone, Debug, PartialEq)]
pub struct TrappingMembraneCandidate<S = TrappingMembraneDiffusionState> {
    state: S,
    step: TrappingMembraneStep,
}
impl<S> TrappingMembraneCandidate<S> {
    pub const fn state(&self) -> &S { &self.state }
    pub const fn step(&self) -> &TrappingMembraneStep { &self.step }
    pub fn into_parts(self) -> (S, TrappingMembraneStep) { (self.state, self.step) }
}

/// Reusable O(cells) nonlinear scratch; the original LU owns all elimination.
#[derive(Debug)]
pub struct TrappingMembraneWorkspace {
    linear: TridiagonalWorkspace,
    diagonal: Vec<f64>, lower: Vec<f64>, upper: Vec<f64>, rhs: Vec<f64>,
    current: Vec<f64>, trial: Vec<f64>, correction: Vec<f64>, sites: Vec<ScalarSiteStep>,
}
impl TrappingMembraneWorkspace {
    fn with_rows(rows: usize) -> Result<Self> {
        Ok(Self { linear: TridiagonalWorkspace::new(rows)?, diagonal: vec![0.0; rows],
            lower: vec![0.0; rows-1], upper: vec![0.0; rows-1], rhs: vec![0.0; rows],
            current: vec![0.0; rows], trial: vec![0.0; rows], correction: vec![0.0; rows],
            sites: vec![ScalarSiteStep::default(); rows] })
    }
}

#[derive(Clone, Copy)]
enum Boundary {
    Prescribed([f64; 2]),
    Finite { pressures: [f64; 2], capacities: [f64; 2] },
}
impl Boundary {
    fn pressures(self) -> [f64; 2] {
        match self { Self::Prescribed(p) | Self::Finite { pressures: p, .. } => p }
    }
    fn offset(self) -> usize { usize::from(matches!(self, Self::Finite { .. })) }
}

impl TrappingMembraneDiffusionState {
    /// Each layer explicitly supplies one immobile site population or None.
    /// Initial trapped amounts are independent of the mobile Henry loading.
    pub fn new(mobile: MembraneDiffusionState, traps: &[Option<LangmuirMembraneTrap>],
        initial_trapped_moles: &[f64], mut cancelled: impl FnMut() -> bool) -> Result<Self> {
        poll(&mut cancelled)?;
        if traps.len() != mobile.model.layers.len() || initial_trapped_moles.len() != mobile.model.cells.len() {
            return Err(TrappingMembraneError::Invalid("one trap binding per layer and one loading per cell required"));
        }
        for (i, cell) in mobile.model.cells.iter().enumerate() {
            if i % 128 == 0 { poll(&mut cancelled)?; }
            let capacity = match traps[cell.layer] {
                Some(trap) => finite(trap.capacity_mol_m3 * (cell.width_m * mobile.model.area_m2))?,
                None => 0.0,
            };
            if !nonnegative(initial_trapped_moles[i]) || initial_trapped_moles[i] > capacity
                || (traps[cell.layer].is_some() && !positive(capacity)) {
                return Err(TrappingMembraneError::Invalid("initial trapped amount must fit its finite authored cell capacity"));
            }
        }
        let state = Self { mobile, traps: Arc::from(traps), trapped_moles: initial_trapped_moles.to_vec() };
        state.stored_amount_mol()?;
        poll(&mut cancelled)?;
        Ok(state)
    }
    pub const fn mobile(&self) -> &MembraneDiffusionState { &self.mobile }
    pub fn trapped_moles(&self) -> &[f64] { &self.trapped_moles }
    pub fn traps(&self) -> &[Option<LangmuirMembraneTrap>] { &self.traps }
    pub fn trapped_amount_mol(&self) -> Result<f64> {
        let mut sum = Sum::default();
        for &amount in &self.trapped_moles { sum.add(amount)?; }
        Ok(sum.value)
    }
    pub fn stored_amount_mol(&self) -> Result<f64> {
        Ok(finite(self.mobile.stored_amount_mol()? + self.trapped_amount_mol()?)?)
    }
    pub fn workspace(&self) -> Result<TrappingMembraneWorkspace> {
        TrappingMembraneWorkspace::with_rows(self.mobile.model.cells.len())
    }
    pub fn step_candidate(&self, reservoirs_pa: [f64; 2], duration_s: f64,
        controls: TrappingMembraneSolve, work: &mut TrappingMembraneWorkspace,
        cancelled: impl FnMut() -> bool) -> Result<TrappingMembraneCandidate> {
        self.solve(Boundary::Prescribed(reservoirs_pa), duration_s, controls, work, cancelled)
    }

    fn row(&self, row: usize, boundary: Boundary) -> (f64, f64, f64, f64, Option<usize>) {
        let model = &self.mobile.model;
        let n = model.cells.len();
        if let Boundary::Finite { pressures, capacities } = boundary {
            if row == 0 { return (capacities[0], pressures[0], 0.0, model.faces[0].permeance, None); }
            if row == n+1 { return (capacities[1], pressures[1], model.faces[n].permeance, 0.0, None); }
        }
        let i = row-boundary.offset();
        (model.cells[i].capacity_mol_m2_pa * model.area_m2, self.mobile.partial_pressures_pa[i],
            model.faces[i].permeance, model.faces[i+1].permeance, Some(i))
    }

    fn evaluate(&self, boundary: Boundary, duration: f64, controls: MembraneDiffusionSolve,
        work: &mut TrappingMembraneWorkspace, trial: bool,
        cancelled: &mut impl FnMut() -> bool) -> Result<(f64, f64)> {
        let q = if trial { &work.trial } else { &work.current };
        let rows = q.len(); let model = &self.mobile.model;
        let mut ratio = 0.0_f64; let mut maximum_residual = 0.0_f64;
        let mut changes = Sum::default(); let mut scales = Sum::default();
        for row in 0..rows {
            if row % 128 == 0 { poll(cancelled)?; }
            let (capacity, old, gl, gr, cell) = self.row(row, boundary);
            let l = finite((duration * gl) * model.area_m2)?;
            let r = finite((duration * gr) * model.area_m2)?;
            if (gl > 0.0 && l == 0.0) || (gr > 0.0 && r == 0.0) {
                return Err(TrappingMembraneError::Invalid("integrated trap-film conductance underflowed"));
            }
            let left = if row > 0 { q[row-1] } else { boundary.pressures()[0] };
            let right = if row+1 < rows { q[row+1] } else { boundary.pressures()[1] };
            let mut sites = ScalarSiteStep::default(); let mut derivative = 0.0;
            if let Some(i) = cell {
                if let Some(trap) = self.traps[model.cells[i].layer] {
                    let site_capacity = trap.capacity_mol_m3 * (model.cells[i].width_m * model.area_m2);
                    let rate = finite(trap.adsorption_rate_pa_inverse_s_inverse * q[row])?;
                    if trap.adsorption_rate_pa_inverse_s_inverse > 0.0 && q[row] > 0.0 && rate == 0.0 {
                        return Err(TrappingMembraneError::Invalid("nonzero membrane trapping rate underflowed"));
                    }
                    sites = backward_euler_site(self.trapped_moles[i], site_capacity, rate,
                        trap.desorption_rate_s_inverse, duration)?;
                    derivative = finite(sites.adsorption_rate_derivative_mol_s
                        * trap.adsorption_rate_pa_inverse_s_inverse)?;
                }
            }
            let mobile_change = finite(capacity * (q[row]-old))?;
            let change = finite(mobile_change + sites.occupied_change_mol)?;
            let incoming = finite(l * (left-q[row]))?;
            let outgoing = finite(r * (q[row]-right))?;
            let residual = finite(change-(incoming-outgoing))?;
            let scale = finite(mobile_change.abs()+sites.occupied_change_mol.abs()+incoming.abs()+outgoing.abs())?;
            ratio = ratio.max(balance_ratio(residual, controls.tolerance(scale)?));
            maximum_residual = maximum_residual.max(residual.abs());
            changes.add(change)?; scales.add(mobile_change.abs()+sites.occupied_change_mol.abs())?;
            let mass = finite(capacity+derivative)?;
            let scale = mass.max(l).max(r);
            let (mass, left, right) = (mass/scale, l/scale, r/scale);
            let diagonal = mass+left+right;
            if !positive(mass) || diagonal <= left+right
                || (l > 0.0 && left == 0.0) || (r > 0.0 && right == 0.0) {
                return Err(TrappingMembraneError::Invalid("Newton matrix cannot retain positive storage at this time scale"));
            }
            work.diagonal[row] = diagonal; work.rhs[row] = finite(-residual/scale)?;
            if row > 0 { work.lower[row-1] = -left; }
            if row+1 < rows { work.upper[row] = -right; }
            work.sites[row] = sites;
        }
        let (incoming, outgoing) = match boundary {
            Boundary::Prescribed(p) => (
                finite(((model.faces[0].permeance * model.area_m2) * duration) * (p[0]-q[0]))?,
                finite(((model.faces[rows].permeance * model.area_m2) * duration) * (q[rows-1]-p[1]))?),
            Boundary::Finite { .. } => (0.0, 0.0),
        };
        let residual = finite(changes.value-(incoming-outgoing))?;
        ratio = ratio.max(balance_ratio(residual, controls.tolerance(scales.value+incoming.abs()+outgoing.abs())?));
        Ok((ratio, maximum_residual.max(residual.abs())))
    }

    fn solve(&self, boundary: Boundary, duration: f64, controls: TrappingMembraneSolve,
        work: &mut TrappingMembraneWorkspace, mut cancelled: impl FnMut() -> bool)
        -> Result<TrappingMembraneCandidate> {
        poll(&mut cancelled)?; controls.balance.validate()?;
        let n = self.mobile.model.cells.len(); let offset = boundary.offset(); let rows = n+2*offset;
        let elapsed = self.mobile.elapsed_s+duration;
        if !positive(duration) || !elapsed.is_finite() || elapsed <= self.mobile.elapsed_s
            || !boundary.pressures().into_iter().all(nonnegative) || work.current.len() != rows
            || !(1..=256).contains(&controls.maximum_iterations)
            || !(1..=64).contains(&controls.maximum_backtracks) {
            return Err(TrappingMembraneError::Invalid("positive representable step, matched workspace and bounded nonlinear work required"));
        }
        for row in 0..rows { work.current[row] = self.row(row,boundary).1; }
        let mut evaluations = 1;
        let mut residual = self.evaluate(boundary,duration,controls.balance,work,false,&mut cancelled)?;
        // A zero-pressure state is a subsolution of the concave, monotone
        // eliminated site law. Newton then approaches the physical root from
        // below. Keep already exact disconnected equilibria at their original
        // datum, so sealed/uniform components do not acquire roundoff flux.
        if residual.0 > 1.0 {
            let mut start = 0;
            while start < rows {
                let mut end = start+1;
                while end < rows && work.upper[end-1] != 0.0 { end += 1; }
                if work.rhs[start..end].iter().any(|&r| r != 0.0) { work.current[start..end].fill(0.0); }
                start = end;
            }
            residual = self.evaluate(boundary,duration,controls.balance,work,false,&mut cancelled)?;
            evaluations += 1;
        }
        let mut iterations = 0;
        while residual.0 > 1.0 {
            poll(&mut cancelled)?;
            if iterations >= controls.maximum_iterations {
                return Err(TrappingMembraneError::NotConverged { iterations, residual_ratio: residual.0 });
            }
            work.linear.solve_into(&work.lower,&work.diagonal,&work.upper,&work.rhs,
                &mut work.correction,&mut cancelled)?;
            let mut fraction = 1.0;
            let mut accepted = None;
            for _ in 0..controls.maximum_backtracks {
                poll(&mut cancelled)?;
                let mut positive_trial = true;
                for row in 0..rows {
                    if row % 128 == 0 { poll(&mut cancelled)?; }
                    work.trial[row] = work.current[row]+fraction*work.correction[row];
                    positive_trial &= nonnegative(work.trial[row]);
                }
                if positive_trial {
                    let next = self.evaluate(boundary,duration,controls.balance,work,true,&mut cancelled)?;
                    evaluations += 1;
                    if next.0 <= 1.0 || next.1 < residual.1 {
                        accepted = Some(next); break;
                    }
                }
                fraction *= 0.5;
            }
            let Some(next) = accepted else {
                return Err(TrappingMembraneError::NotConverged { iterations, residual_ratio: residual.0 });
            };
            core::mem::swap(&mut work.current,&mut work.trial);
            residual = next; iterations += 1;
        }
        let pressures = if offset == 0 { boundary.pressures() } else { [work.current[0],work.current[rows-1]] };
        let mobile = MembraneDiffusionState { model: self.mobile.model.clone(),
            partial_pressures_pa: work.current[offset..offset+n].to_vec(), elapsed_s: elapsed };
        let state = Self { mobile, traps: self.traps.clone(),
            trapped_moles: work.sites[offset..offset+n].iter().map(|s| s.occupied_moles).collect() };
        let transport = state.mobile.sample_fluxes(pressures,&mut cancelled)?;
        let area = self.mobile.model.area_m2;
        let boundary_transfers = [finite((transport.face_molar_fluxes_mol_m2_s[0]*area)*duration)?,
            finite((transport.face_molar_fluxes_mol_m2_s[n]*area)*duration)?];
        let mut mobile_change = Sum::default(); let mut trapped_change = Sum::default();
        let mut maximum_cell: f64 = 0.0;
        for i in 0..n {
            if i % 128 == 0 { poll(&mut cancelled)?; }
            let mobile = finite(self.mobile.model.cells[i].capacity_mol_m2_pa*area
                * (state.mobile.partial_pressures_pa[i]-self.mobile.partial_pressures_pa[i]))?;
            let trapped = finite(state.trapped_moles[i]-self.trapped_moles[i])?;
            mobile_change.add(mobile)?; trapped_change.add(trapped)?;
            let incoming = finite((transport.face_molar_fluxes_mol_m2_s[i]*area)*duration)?;
            let outgoing = finite((transport.face_molar_fluxes_mol_m2_s[i+1]*area)*duration)?;
            let residual = finite((mobile+trapped)-(incoming-outgoing))?;
            let tolerance=controls.balance.tolerance(mobile.abs()+trapped.abs()+incoming.abs()+outgoing.abs())?;
            check_balance(residual,tolerance)?;
            // Stored-amount conservation alone cannot detect a kinetic change
            // rounded away beside a large preloaded site inventory. Audit the
            // original site's independent backward-Euler balance at this final
            // candidate, using the caller's increment-based amount allowance.
            let site_residual=work.sites[offset+i].kinetic_balance_residual_mol.abs();
            check_balance(site_residual,tolerance)?;
            maximum_cell = maximum_cell.max(residual.abs()).max(site_residual);
        }
        let membrane_residual = finite((mobile_change.value+trapped_change.value)
            -(boundary_transfers[0]-boundary_transfers[1]))?;
        check_balance(membrane_residual,controls.balance.tolerance(mobile_change.value.abs()
            +trapped_change.value.abs()+boundary_transfers[0].abs()+boundary_transfers[1].abs())?)?;
        let finite_reservoirs = if let Boundary::Finite { pressures: before, capacities } = boundary {
            let changes = [finite(capacities[0]*(pressures[0]-before[0]))?,
                finite(capacities[1]*(pressures[1]-before[1]))?];
            let residuals = [finite(changes[0]+boundary_transfers[0])?,finite(changes[1]-boundary_transfers[1])?];
            for side in 0..2 { check_balance(residuals[side],controls.balance.tolerance(
                changes[side].abs()+boundary_transfers[side].abs())?)?; }
            let mut total = Sum::default();
            for value in [changes[0],mobile_change.value,trapped_change.value,changes[1]] { total.add(value)?; }
            check_balance(total.value,controls.balance.tolerance(changes[0].abs()+changes[1].abs()
                +mobile_change.value.abs()+trapped_change.value.abs())?)?;
            Some(TrappingReservoirBalance { pressures_before_pa: before, pressures_after_pa: pressures,
                amount_changes_mol: changes, amount_residuals_mol: residuals,
                closed_amount_residual_mol: total.value })
        } else { None };
        let step = TrappingMembraneStep { duration_s: duration, transport,
            mobile_before_mol: self.mobile.stored_amount_mol()?, mobile_after_mol: state.mobile.stored_amount_mol()?,
            mobile_change_mol: mobile_change.value, trapped_before_mol: self.trapped_amount_mol()?,
            trapped_after_mol: state.trapped_amount_mol()?, trapped_change_mol: trapped_change.value,
            gross_adsorption_mol: work.sites[offset..offset+n].iter().map(|s|s.gross_adsorption_mol).collect(),
            gross_desorption_mol: work.sites[offset..offset+n].iter().map(|s|s.gross_desorption_mol).collect(),
            boundary_transfers_mol: boundary_transfers, membrane_amount_residual_mol: membrane_residual,
            maximum_cell_amount_residual_mol: maximum_cell, finite_reservoirs,
            nonlinear_iterations: iterations, residual_evaluations: evaluations };
        poll(&mut cancelled)?;
        Ok(TrappingMembraneCandidate { state, step })
    }
}

fn balance_ratio(residual: f64, tolerance: f64) -> f64 {
    if residual == 0.0 { 0.0 } else { residual.abs()/tolerance }
}
fn check_balance(residual: f64, tolerance: f64) -> Result<()> {
    if residual.abs() > tolerance {
        Err(MembraneDiffusionError::AmountBalance { residual_mol: residual, tolerance_mol: tolerance }.into())
    } else { Ok(()) }
}

/// Both reservoir inventories participate in the SAME nonlinear pressure solve.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteTrappingMembraneDiffusionState {
    membrane: TrappingMembraneDiffusionState,
    capacities: [MolarPressureCapacity; 2],
    pressures: [f64; 2],
}
impl FiniteTrappingMembraneDiffusionState {
    pub fn new(membrane: TrappingMembraneDiffusionState, capacities: [MolarPressureCapacity; 2],
        pressures: [f64; 2]) -> Result<Self> {
        if !pressures.into_iter().all(nonnegative) {
            return Err(TrappingMembraneError::Invalid("finite trapped-membrane reservoir pressures must be nonnegative"));
        }
        let state = Self { membrane, capacities, pressures };
        state.total_amount_mol()?;
        Ok(state)
    }
    pub const fn membrane(&self) -> &TrappingMembraneDiffusionState { &self.membrane }
    /// Retain the mobile and trapped inventories while reseeding actual bath pressures.
    pub fn into_membrane(self) -> TrappingMembraneDiffusionState { self.membrane }
    pub const fn reservoir_pressures_pa(&self) -> [f64; 2] { self.pressures }
    pub const fn capacities(&self) -> [MolarPressureCapacity; 2] { self.capacities }
    pub fn total_amount_mol(&self) -> Result<f64> {
        let mut sum = Sum::default(); sum.add(self.membrane.stored_amount_mol()?)?;
        for side in 0..2 { sum.add(finite(self.capacities[side].value_si()*self.pressures[side])?)?; }
        Ok(sum.value)
    }
    pub fn workspace(&self) -> Result<TrappingMembraneWorkspace> {
        TrappingMembraneWorkspace::with_rows(self.membrane.mobile.model.cells.len()+2)
    }
    pub fn step_candidate(&self, duration_s: f64, controls: TrappingMembraneSolve,
        work: &mut TrappingMembraneWorkspace, cancelled: impl FnMut() -> bool)
        -> Result<TrappingMembraneCandidate<Self>> {
        let (membrane, step) = self.membrane.solve(Boundary::Finite { pressures: self.pressures,
            capacities: self.capacities.map(MolarPressureCapacity::value_si) },duration_s,controls,work,cancelled)?.into_parts();
        let pressures = step.finite_reservoirs.as_ref().expect("finite reservoir report").pressures_after_pa;
        let state = Self::new(membrane,self.capacities,pressures)?;
        Ok(TrappingMembraneCandidate { state, step })
    }
}

pub mod adaptive;

#[cfg(test)]
mod tests;
