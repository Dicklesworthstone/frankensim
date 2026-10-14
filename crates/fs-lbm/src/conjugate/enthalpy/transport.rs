//! Shared enthalpy-state evaluation for current generation and thermal transport.
//! The existing network owns storage/CFL/heat exchange; the current owner owns
//! bulk and electrical-contact heating. No mass or apparent Cp is re-assigned.
use super::*;
use std::sync::Arc;
use fs_conduction::transient::deposition::ContactPatch;

/// One constitutive evaluation with its original concrete material payload.
/// Unlike fixed-Cp fluid transport, this carries no replacement storage law:
/// the original enthalpy chart continues to own sensible and latent energy.
#[derive(Debug, PartialEq)]
pub struct EnthalpyTransport<F> {
    material: F,
    resistivity_ohm_m: Option<f64>,
    conductivity_w_m_k: f64,
}
impl<F> EnthalpyTransport<F> {
    /// Explicit insulation is None electrically or zero thermally. Missing
    /// source support must be an error, never interpreted as either value.
    pub fn try_new(material: F, resistivity_ohm_m: Option<f64>, conductivity_w_m_k: f64)
        -> Result<Self, JouleEnthalpyError> {
        if resistivity_ohm_m.is_some_and(|r| !r.is_finite() || r <= 0.0)
            || !conductivity_w_m_k.is_finite() || conductivity_w_m_k < 0.0 {
            return Err(invalid("invalid electrical resistivity or thermal conductivity"));
        }
        Ok(Self { material, resistivity_ohm_m, conductivity_w_m_k })
    }
    /// Original material/phase/source evaluation, without type erasure.
    pub const fn material(&self) -> &F { &self.material }
    /// Bulk Ohmic response [ohm m]; unused inside ideal electrode cells.
    pub const fn resistivity_ohm_m(&self) -> Option<f64> { self.resistivity_ohm_m }
    /// Isotropic thermal conductivity [W/(m K)] supplied at the same state.
    pub const fn conductivity_w_m_k(&self) -> f64 { self.conductivity_w_m_k }
}

/// A mapped cell's immutable evaluation point and material response.
#[derive(Debug, PartialEq)]
pub struct EnthalpyTransportSample<F> {
    /// Original thermal cell, including geometry, mass, material ID and h.
    pub cell: DepositedCell,
    /// Temperature queried from that cell's existing enthalpy law [K].
    pub temperature_k: f64,
    /// Response used by both the electrical solve and thermal-link resolver.
    pub response: EnthalpyTransport<F>,
}

/// A thermal link's applied response, retaining its original interface payload.
#[derive(Debug, PartialEq)]
pub struct EnthalpyLinkSample<I> {
    /// Original authored link and both surface coverages, not a rewritten input.
    pub patch: ContactPatch,
    /// Total conductance [W/K] actually used by the thermal kernel.
    pub conductance_w_k: f64,
    /// Original concrete link law/source evaluation.
    pub material: I,
}

/// One shared source-evaluation epoch. Objects are immutable; no Clone bound on
/// material payloads. An Arc retains them without cloning source-owned state.
#[derive(Debug, PartialEq)]
pub struct EnthalpyTransportEpoch<F, I> {
    /// Physical source evaluation time [s].
    pub time_s: f64,
    /// Electrical voxel order, including thermally evolving ideal electrodes.
    pub cells: Vec<EnthalpyTransportSample<F>>,
    /// Explicit requested thermal-link order.
    pub links: Vec<EnthalpyLinkSample<I>>,
}

/// Original electrical/thermal accounts and matching material evaluations.
#[derive(Debug)]
pub struct JouleTransportReport<F, I> {
    /// Same energy/work/source report as the original fixed-link entry point.
    pub joule: JouleEnthalpyReport,
    /// Last actually applied epoch; absent on a zero-duration advance.
    pub last_applied: Option<Arc<EnthalpyTransportEpoch<F, I>>>,
    /// Fresh endpoint support evaluation, not applied in the preceding step.
    pub endpoint: Arc<EnthalpyTransportEpoch<F, I>>,
}

type BulkLaw<'a, F> = dyn Fn(usize, &DepositedCell, f64)
    -> Result<EnthalpyTransport<F>, JouleEnthalpyError> + 'a;
type LinkLaw<'a, F, I> = dyn Fn(&ContactPatch, [&EnthalpyTransportSample<F>; 2])
    -> Result<(f64, I), JouleEnthalpyError> + 'a;

struct TransportSource<'a, F, I> {
    account: Account<'a>,
    selected: &'a [u64],
    indices: BTreeMap<CellId, usize>,
    bulk: &'a BulkLaw<'a, F>,
    link: &'a LinkLaw<'a, F, I>,
    last: Option<Arc<EnthalpyTransportEpoch<F, I>>>,
}
impl<F, I> Clone for TransportSource<'_, F, I> {
    fn clone(&self) -> Self {
        Self { account:self.account.clone(), selected:self.selected,
            indices:self.indices.clone(), bulk:self.bulk, link:self.link, last:self.last.clone() }
    }
}
impl<F, I> TransportSource<'_, F, I> {
    fn prepare(&self, thermal: &DepositionSimulation<'_>)
        -> Result<EnthalpyTransportEpoch<F, I>, JouleEnthalpyError> {
        checkpoint(self.account.gate)?;
        // Bound constitutive work before sampling, not just before the solve.
        if self.account.solves >= self.account.input.max_current_solves {
            return Err(JouleEnthalpyError::CurrentBudget);
        }
        admit_map(thermal,self.account.input,self.account.gate)?;
        if self.selected.len() > thermal.links().len() {
            return Err(invalid("selected thermal-link count exceeds active links"));
        }
        let mut patches=Vec::with_capacity(self.selected.len());
        let mut seen=BTreeSet::new();
        for &id in self.selected {
            checkpoint(self.account.gate)?;
            let patch=thermal.contact_patch(id).ok_or(invalid("selected thermal link is inactive"))?;
            if !seen.insert(id) || !self.indices.contains_key(&patch.link.a)
                || !self.indices.contains_key(&patch.link.b) {
                return Err(invalid("selected thermal links must be unique with two mapped receiving cells"));
            }
            patches.push(patch);
        }
        let mut cells=Vec::with_capacity(self.account.input.cell_ids.len());
        for (c,&id) in self.account.input.cell_ids.iter().enumerate() {
            checkpoint(self.account.gate)?;
            let cell=&thermal.cells()[&id];
            let temperature_k=thermal.temperature_k(id)?;
            let response=(self.bulk)(c,cell,temperature_k)?;
            cells.push(EnthalpyTransportSample { cell:cell.clone(), temperature_k, response });
        }
        let mut links=Vec::with_capacity(patches.len());
        for patch in patches {
            checkpoint(self.account.gate)?;
            let sides=[&cells[self.indices[&patch.link.a]], &cells[self.indices[&patch.link.b]]];
            let (conductance_w_k,material)=(self.link)(&patch,sides)?;
            if !conductance_w_k.is_finite() || conductance_w_k < 0.0 {
                return Err(invalid("thermal link law returned an invalid conductance"));
            }
            links.push(EnthalpyLinkSample {patch,conductance_w_k,material});
        }
        checkpoint(self.account.gate)?;
        Ok(EnthalpyTransportEpoch {time_s:thermal.time_s(),cells,links})
    }
    fn evaluate(&mut self, thermal:&DepositionSimulation<'_>)
        -> Result<(EnthalpyCurrentState, Arc<EnthalpyTransportEpoch<F,I>>), JouleEnthalpyError> {
        let epoch=Arc::new(self.prepare(thermal)?);
        let electrical=|c:usize,_:&DepositedCell,_:f64| Ok(epoch.cells[c].response.resistivity_ohm_m());
        let state=self.account.evaluate(thermal,&electrical)?;
        Ok((state,epoch))
    }
}
impl<F, I> DepositionObserver for TransportSource<'_, F, I> {
    type Error=JouleEnthalpyError;
    fn event(&mut self,_:&Event,_:&DepositionSimulation<'_>)->Result<(),Self::Error>{Ok(())}
    fn sources(&mut self,thermal:&DepositionSimulation<'_>,check:&mut dyn FnMut()->Result<(),DepositionError>)
        ->Result<BTreeMap<CellId,f64>,Self::Error> {
        check()?;
        let (state,epoch)=self.evaluate(thermal)?;
        let power=self.account.apply(state)?;
        self.last=Some(epoch);
        Ok(power)
    }
    fn conductances(&mut self,thermal:&DepositionSimulation<'_>,check:&mut dyn FnMut()->Result<(),DepositionError>)
        ->Result<BTreeMap<u64,f64>,Self::Error> {
        check()?;
        let epoch=self.last.as_ref().ok_or(invalid("thermal transport has no prepared material epoch"))?;
        if epoch.time_s != thermal.time_s() || self.account.last_dt.is_some() {
            return Err(invalid("material transport epoch does not match the thermal substep"));
        }
        Ok(epoch.links.iter().map(|s|(s.patch.link.id,s.conductance_w_k)).collect())
    }
    fn step(&mut self,start:f64,dt:f64,_:&BTreeMap<CellId,f64>,_:&DepositionSimulation<'_>,
        check:&mut dyn FnMut()->Result<(),DepositionError>)->Result<(),Self::Error> {
        check()?;
        self.account.record(start,dt)
    }
}

/// Advance generated current heating AND local thermal transport on one enthalpy
/// timeline. `resolve` is called once per mapped cell/epoch, including electrodes,
/// and can use h to distinguish equal-temperature phase states. `resolve_link`
/// receives those SAME ordered cell evaluations and the original contact patch.
/// It explicitly supplies total W/K, retaining any material/interface payload.
/// No heat-transfer geometry, contact law or mixing rule is silently inferred.
///
/// Only `link_ids` are rebound. Each must stay active with two mapped endpoints;
/// unselected links keep authored conductance. Original per-side coverage and
/// templates survive stepping, including zero conductance. Selected-link changes
/// at events are admitted on the next evaluation; missing links refuse.
///
/// Endpoint current, bulk and selected-link support are checked before returning
/// a candidate. Current budgets cover every substep plus the endpoint, and limit
/// sampling epochs too. No callback may debit mutable history during trials.
/// Original typed failures and accepted-state rollback are retained. Extra memory
/// is O(mapped cells + selected links), plus existing thermal/current scratch.
///
/// Fixed masses, geometry and caloric laws remain with the original thermal owner.
/// The existing CFL uses the actual applied conductances but does not bound source
/// stiffness or arbitrary nonlinear constitutive feedback; refine max_dt_s. This
/// is explicit stationary caloric coupling, not melt motion, implicit feedback,
/// pressure/volume work, induction, total thermomechanics or empirical validation.
#[allow(clippy::too_many_arguments)]
pub fn advance_joule_enthalpy_with_material_transport<'t, F, I>(
    thermal:&DepositionSimulation<'t>, target_s:f64, config:StepConfig,
    input:JouleEnthalpyInput<'_>, link_ids:&[u64], gate:&CancelGate,
    resolve:impl Fn(usize,&DepositedCell,f64)->Result<EnthalpyTransport<F>,JouleEnthalpyError>,
    resolve_link:impl Fn(&ContactPatch,[&EnthalpyTransportSample<F>;2])->Result<(f64,I),JouleEnthalpyError>,
) -> Result<(DepositionSimulation<'t>,JouleTransportReport<F,I>),JouleEnthalpyError> {
    admit_map(thermal,input,gate)?;
    let source=TransportSource {account:Account::new(input,gate),selected:link_ids,
        indices:input.cell_ids.iter().copied().enumerate().map(|(c,id)|(id,c)).collect(),
        bulk:&resolve,link:&resolve_link,last:None};
    let (next,mut source,audit)=thermal.advance_observed(&source,target_s,config,||checkpoint(gate))?;
    let (endpoint_current,endpoint)=source.evaluate(&next)?;
    let joule=source.account.finish(audit,endpoint_current)?;
    checkpoint(gate)?;
    Ok((next,JouleTransportReport {joule,last_applied:source.last,endpoint}))
}
