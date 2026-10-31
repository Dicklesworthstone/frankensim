//! Reciprocal motion of windings sharing one uniform magnetic gap.
//!
//! With signed turns `n`, fixed leakage `D`, gap `g` and `k = mu A`, the complete
//! inductance without a core is `L(g) = D + (k/g) n n^T`. An optional original
//! reversible H(B) source adds a finite, saturating return-path reluctance.
//! Every retained winding has explicit, strictly positive leakage;
//! ideal rank-deficient flux coordinates are refused.
//! One magnetic energy supplies all terminal currents and the gap force. The
//! original descriptor advances electrical linkage, gap and source-derived body
//! momentum together, so changes in mutual inductance generate reciprocal emf.

use fs_material::field_storage::ResolvedInductor;
use fs_mbd_material::ResolvedBody;
use fs_phs::circuit::{Branch, CircuitDae, CircuitGraph, MutualInductance};
use fs_phs::{DescriptorPortHamiltonian, QuadraticStorage, Storage, step_descriptor};

use super::super::contacts::{self, SeriesElectricalContact};
use super::super::{
    MAX_MATERIAL_CIRCUIT_BRANCHES, MAX_PORTS, MAX_STATES, account_step, admit_graph,
    bind_with_mutual, extended, poll, positive,
};
use super::{
    AxialSupport, FieldMotionInput, FieldMotionStep, MaterialCircuit, MaterialCircuitError,
    MaterialCircuitSolve, MaterialResistor, MaterialStorageElement, coupling_scale, invalid, phs,
    within,
};

mod core;
pub use core::ReversibleMagneticCore;

/// A winding in original graph branch order, with a physical leakage store.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GapWinding {
    /// Original inductor branch index, before material/contact lowering.
    pub branch: usize,
    /// Nonzero integer turns. The sign fixes the dot convention relative to
    /// positive branch current from the first terminal to the second.
    pub signed_turns: i32,
    /// Constant uncoupled leakage inductance [H], strictly positive.
    /// This is an explicit model input, never a numerical regularization.
    pub leakage_h: f64,
}

/// Source-resolved gap medium and a set of windings sharing its common flux.
///
/// The one-turn `medium` retains the actual permeability receipt, constant area
/// and present gap. The return path has negligible reluctance unless an explicit
/// reversible core is bound. Leakage is fixed and does not act on this gap.
/// Fringing, hysteresis, skin loss, changing leakage and thermal coefficient
/// evolution need additional owners.
#[derive(Clone, Debug, PartialEq)]
pub struct CommonMagneticGap {
    medium: ResolvedInductor,
    windings: Vec<GapWinding>,
    limits_m: [f64; 2],
    scale: f64,
    inverse_leakage_turns: f64,
    core: Option<ReversibleMagneticCore>,
}

impl CommonMagneticGap {
    /// Bind at least two windings to the same source-resolved one-turn gap.
    /// The caller supplies positive OPEN geometric limits and all physical
    /// leakage inductances. A validity-wide constant gap permeability is required
    /// so that magnetic energy does not omit a temperature-dependent entropy.
    pub fn new(
        medium: ResolvedInductor,
        mut windings: Vec<GapWinding>,
        limits_m: [f64; 2],
    ) -> Result<Self, MaterialCircuitError> {
        if medium.turns() != 1
            || windings.len() < 2
            || windings.len() > MAX_MATERIAL_CIRCUIT_BRANCHES
            || limits_m.iter().any(|&g| !positive(g))
            || !within(medium.path_length_m(), limits_m)
        {
            return Err(invalid(
                "shared magnetic gap requires one-turn geometry, bounded windings and positive open gap limits",
            ));
        }
        let scale = coupling_scale(&MaterialStorageElement::Inductor(medium.clone()))?;
        windings.sort_by_key(|w| w.branch);
        let mut inverse_leakage_turns = 0.0;
        for (index, winding) in windings.iter().enumerate() {
            if winding.signed_turns == 0
                || !positive(winding.leakage_h)
                || !positive(1.0 / winding.leakage_h)
                || (index > 0 && windings[index - 1].branch == winding.branch)
            {
                return Err(invalid(
                    "shared-gap windings require unique branches, nonzero signed turns and positive explicit leakage",
                ));
            }
            let turns = f64::from(winding.signed_turns);
            inverse_leakage_turns += turns * (turns / winding.leakage_h);
        }
        if !positive(inverse_leakage_turns)
            || !positive(1.0 / inverse_leakage_turns)
            || limits_m.iter().any(|&g| {
                !positive(scale / g)
                    || !positive(g / scale)
                    || !positive(inverse_leakage_turns + g / scale)
            })
        {
            return Err(invalid(
                "shared-gap reluctance or leakage is unrepresentable",
            ));
        }
        Ok(Self {
            medium,
            windings,
            limits_m,
            scale,
            inverse_leakage_turns,
            core: None,
        })
    }

    /// Add one source-resolved reversible core in series with this gap.
    /// Its complete H(B) source range remains authoritative; neither a magnetic
    /// trial nor an accepted endpoint may extrapolate beyond that range.
    pub fn with_reversible_core(
        mut self,
        core: ReversibleMagneticCore,
    ) -> Result<Self, MaterialCircuitError> {
        if self.core.is_some() {
            return Err(invalid("shared gap already has a reversible core"));
        }
        for bound in self.limits_m {
            core.admit_reluctance(self.inverse_leakage_turns + bound / self.scale)?;
        }
        self.core = Some(core);
        Ok(self)
    }

    /// Original reversible core source and physical geometry, if present.
    pub const fn core(&self) -> Option<&ReversibleMagneticCore> {
        self.core.as_ref()
    }

    /// Original material receipt, constant area and current one-turn gap geometry.
    pub const fn medium(&self) -> &ResolvedInductor {
        &self.medium
    }
    /// Winding identities, signed turns and leakage, sorted by original branch.
    pub fn windings(&self) -> &[GapWinding] {
        &self.windings
    }
    /// Accepted open interval for pole separation [m].
    pub const fn limits_m(&self) -> [f64; 2] {
        self.limits_m
    }
    /// Current shared magnetic gap [m].
    pub fn gap_m(&self) -> f64 {
        self.medium.path_length_m()
    }

    fn flux(&self, gap: f64, linkage: impl Fn(usize) -> f64) -> f64 {
        let numerator = self
            .windings
            .iter()
            .enumerate()
            .map(|(i, w)| f64::from(w.signed_turns) * (linkage(i) / w.leakage_h))
            .sum::<f64>();
        let reluctance = self.inverse_leakage_turns + gap / self.scale;
        self.core.as_ref().map_or_else(
            || numerator / reluctance,
            |core| core.flux(reluctance, numerator),
        )
    }

    fn current(&self, winding: usize, linkage: f64, common_flux: f64) -> f64 {
        let w = &self.windings[winding];
        (linkage - f64::from(w.signed_turns) * common_flux) / w.leakage_h
    }

    fn energy(&self, gap: f64, linkage: impl Fn(usize) -> f64) -> f64 {
        let common_flux = self.flux(gap, &linkage);
        // Positive leakage + gap energies avoid subtracting nearly equal
        // quadratic forms when the common magnetic coupling is strong.
        let leakage = self
            .windings
            .iter()
            .enumerate()
            .map(|(i, w)| {
                let current = self.current(i, linkage(i), common_flux);
                0.5 * current * (w.leakage_h * current)
            })
            .sum::<f64>();
        leakage
            + 0.5 * common_flux * ((gap / self.scale) * common_flux)
            + self
                .core
                .as_ref()
                .map_or(0.0, |core| core.energy(common_flux))
    }

    // The original CircuitDae needs the exact accepted-state winding currents
    // to admit ideal constraints. A secant stamp supplies those currents only;
    // its quadratic energy is never used by the nonlinear coupled Storage.
    fn secant_permeance(
        &self,
        gap: f64,
        linkage: impl Fn(usize) -> f64,
    ) -> Result<f64, MaterialCircuitError> {
        let permeance = match &self.core {
            None => self.scale / gap,
            Some(core) => {
                let flux = self.flux(gap, linkage);
                1.0 / (gap / self.scale + core.secant_reluctance(flux))
            }
        };
        if !positive(permeance) {
            return Err(invalid(
                "shared-gap operating point leaves the core source range or has unrepresentable permeance",
            ));
        }
        Ok(permeance)
    }

    fn linkage_slots(&self, graph: &CircuitGraph) -> Result<Vec<usize>, MaterialCircuitError> {
        self.windings
            .iter()
            .map(|winding| {
                if graph
                    .branches
                    .get(winding.branch)
                    .is_none_or(|(_, _, branch)| !matches!(branch, Branch::Inductor { .. }))
                {
                    return Err(invalid(
                        "shared-gap winding must name an original inductor branch",
                    ));
                }
                Ok(graph.branches[..winding.branch]
                    .iter()
                    .filter(|(_, _, branch)| matches!(branch, Branch::Inductor { .. }))
                    .count())
            })
            .collect()
    }

    fn at_gap(&self, gap: f64) -> Result<Self, MaterialCircuitError> {
        if !within(gap, self.limits_m) {
            return Err(invalid("shared magnetic gap left its admitted open domain"));
        }
        let mut next = self.clone();
        next.medium = ResolvedInductor::try_from_material(
            self.medium.material().clone(),
            self.medium.area_m2(),
            gap,
            1,
        )?;
        Ok(next)
    }

    fn lower_at_permeance(
        &self,
        graph: &CircuitGraph,
        permeance: f64,
    ) -> Result<(CircuitGraph, Vec<MutualInductance>), MaterialCircuitError> {
        let mut lowered = graph.clone();
        let mut mutual = Vec::with_capacity(self.windings.len() * (self.windings.len() - 1) / 2);
        for (index, winding) in self.windings.iter().enumerate() {
            let branch = lowered
                .branches
                .get_mut(winding.branch)
                .ok_or_else(|| invalid("shared-gap winding branch is out of range"))?;
            if !matches!(branch.2, Branch::Inductor { .. }) {
                return Err(invalid(
                    "shared-gap winding must name an original inductor branch",
                ));
            }
            let turns = f64::from(winding.signed_turns);
            let self_l = winding.leakage_h + permeance * turns * turns;
            if !positive(self_l) || !positive(1.0 / self_l) {
                return Err(invalid("shared-gap self inductance is unrepresentable"));
            }
            branch.2 = Branch::Inductor { henries: self_l };
            for other in &self.windings[..index] {
                let henries = permeance * turns * f64::from(other.signed_turns);
                if !henries.is_finite() || henries == 0.0 {
                    return Err(invalid("shared-gap mutual inductance is unrepresentable"));
                }
                mutual.push(MutualInductance {
                    branches: (other.branch, winding.branch),
                    henries,
                });
            }
        }
        Ok((lowered, mutual))
    }
}

/// Complete conserved electrical coordinates and guided-body velocity at assembly.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledGapInitial {
    /// Total winding flux linkages in original inductor branch order [Wb].
    pub fluxes_wb: Vec<f64>,
    /// Charges in original capacitor branch order [C].
    pub charges_c: Vec<f64>,
    /// Original electrical ports used to check ideal algebraic constraints.
    pub electrical_ports: Vec<f64>,
    /// Initial velocity, positive toward increasing separation [m/s].
    pub velocity_m_s: f64,
}

/// One accepted circuit, shared gap, material body and optional frozen support.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledGapMotion {
    circuit: MaterialCircuit,
    gap: CommonMagneticGap,
    body: ResolvedBody,
    support: Option<AxialSupport>,
    flux_slots: Vec<usize>,
    momentum_kg_m_s: f64,
}

/// An atomic electrical/mechanical proposal and the original physical work ledger.
pub struct CoupledGapCandidate {
    model: CoupledGapMotion,
    step: FieldMotionStep,
}
impl CoupledGapCandidate {
    /// Joint proposed endpoint, retaining all field, body and contact sources.
    pub const fn model(&self) -> &CoupledGapMotion {
        &self.model
    }
    /// Branch heat and both external supplies from the same descriptor solve.
    pub const fn step(&self) -> &FieldMotionStep {
        &self.step
    }
    /// Accept the complete model and its work/heat record together.
    pub fn into_parts(self) -> (CoupledGapMotion, FieldMotionStep) {
        (self.model, self.step)
    }
}

impl CoupledGapMotion {
    /// Replace the selected graph inductors with the complete shared-gap matrix.
    /// Other C/L branches and transformer constraints remain in the original
    /// circuit. Every resistor still requires its actual material source.
    ///
    /// All common and leakage magnetic energy belongs to `gap`; no independent
    /// per-winding material storage is invented. The graph's selected numerical
    /// L values are placeholders, overwritten before native matrix admission.
    pub fn new(
        graph: CircuitGraph,
        sources: Vec<MaterialResistor>,
        gap: CommonMagneticGap,
        body: ResolvedBody,
        initial: CoupledGapInitial,
    ) -> Result<Self, MaterialCircuitError> {
        Self::new_with_contacts_and_support(graph, sources, Vec::new(), gap, body, None, initial)
    }

    /// Assemble the same gap with original series contacts and an optional
    /// source-bound axial return member. The member contributes both elastic
    /// energy and consistent modal mass. Thermally live supports refuse because
    /// this owner does not evolve their enthalpy or reciprocal reversible heat.
    pub fn new_with_contacts_and_support(
        graph: CircuitGraph,
        sources: Vec<MaterialResistor>,
        contacts: Vec<SeriesElectricalContact>,
        gap: CommonMagneticGap,
        body: ResolvedBody,
        support: Option<AxialSupport>,
        initial: CoupledGapInitial,
    ) -> Result<Self, MaterialCircuitError> {
        admit_graph(&graph)?;
        if !initial.velocity_m_s.is_finite()
            || support
                .as_ref()
                .is_some_and(AxialSupport::is_thermally_live)
        {
            return Err(invalid(
                "shared-gap motion requires finite velocity and a frozen support",
            ));
        }
        if let Some(support) = &support {
            support.force_and_energy(gap.gap_m())?;
        }
        let mass =
            body.inertia().mass() + support.as_ref().map_or(0.0, AxialSupport::modal_mass_kg);
        let momentum = mass * initial.velocity_m_s;
        let kinetic = 0.5 * momentum * (momentum / mass);
        if !positive(mass)
            || !positive(1.0 / mass)
            || !momentum.is_finite()
            || (initial.velocity_m_s != 0.0 && momentum == 0.0)
            || !kinetic.is_finite()
            || (momentum != 0.0 && kinetic == 0.0)
        {
            return Err(invalid(
                "source-derived shared-gap body mass or momentum is unrepresentable",
            ));
        }
        let flux_slots = gap.linkage_slots(&graph)?;
        if initial.fluxes_wb.len()
            != graph
                .branches
                .iter()
                .filter(|(_, _, branch)| matches!(branch, Branch::Inductor { .. }))
                .count()
        {
            return Err(invalid(
                "shared-gap initial linkages must cover all original inductors",
            ));
        }
        let permeance = gap.secant_permeance(gap.gap_m(), |i| initial.fluxes_wb[flux_slots[i]])?;
        let (graph, mutual) = gap.lower_at_permeance(&graph, permeance)?;
        let circuit = MaterialCircuit::new_with_contacts_and_mutual_inductance(
            graph,
            sources,
            contacts,
            mutual,
            initial.fluxes_wb,
            initial.charges_c,
            &initial.electrical_ports,
        )?;
        // Admit the entire geometric domain through the native full-matrix
        // checker, including resolved leakage at the most strongly coupled end.
        // Positive core reluctance only reduces coupling, so the no-core bound
        // remains conservative for every supported saturation segment. No initial
        // flux is projected onto an endpoint current clamp here.
        for bound in gap.limits_m {
            let (bounded, links) = gap.lower_at_permeance(&circuit.graph, gap.scale / bound)?;
            bind_with_mutual(
                &bounded,
                circuit.resistors.clone(),
                &circuit.contacts,
                &links,
            )?;
        }
        let model = Self {
            circuit,
            gap,
            body,
            support,
            flux_slots,
            momentum_kg_m_s: momentum,
        };
        model.magnetic_energy_j()?;
        model.field_force_n()?;
        Ok(model)
    }

    /// Original circuit topology with coefficients reproducing the accepted
    /// endpoint currents. With a nonlinear core these are SECANT coefficients
    /// for native constraint admission, not its physical magnetic energy or
    /// differential inductance. Use this owner's corresponding accessors.
    pub const fn circuit(&self) -> &MaterialCircuit {
        &self.circuit
    }
    /// Source-bound common gap and every winding's physical leakage and polarity.
    pub const fn gap(&self) -> &CommonMagneticGap {
        &self.gap
    }
    /// Original regional material/geometry receipts and spatial inertia.
    pub const fn body(&self) -> &ResolvedBody {
        &self.body
    }
    /// Original layered support, if present.
    pub const fn support(&self) -> Option<&AxialSupport> {
        self.support.as_ref()
    }
    /// Actual separation shared by all windings [m].
    pub fn gap_m(&self) -> f64 {
        self.gap.gap_m()
    }
    /// Body plus consistent support mass [kg], never an authored mass override.
    pub fn generalized_mass_kg(&self) -> f64 {
        self.body.inertia().mass()
            + self
                .support
                .as_ref()
                .map_or(0.0, AxialSupport::modal_mass_kg)
    }
    /// Generalized translational momentum [kg m/s].
    pub const fn momentum_kg_m_s(&self) -> f64 {
        self.momentum_kg_m_s
    }
    /// Endpoint velocity, positive toward opening the gap [m/s].
    pub fn velocity_m_s(&self) -> f64 {
        self.momentum_kg_m_s / self.generalized_mass_kg()
    }
    /// Physical flux common to all the windings [Wb], distinct from linkage.
    pub fn common_flux_wb(&self) -> Result<f64, MaterialCircuitError> {
        let flux = self
            .gap
            .flux(self.gap_m(), |i| self.circuit.fluxes_wb[self.flux_slots[i]]);
        if !flux.is_finite() {
            return Err(invalid("common magnetic flux is unrepresentable"));
        }
        Ok(flux)
    }
    /// Endpoint branch currents in `gap().windings()` order [A]. Each depends
    /// on every linkage; dividing linkage by the branch's self L is incorrect.
    pub fn winding_currents_a(&self) -> Result<Vec<f64>, MaterialCircuitError> {
        let flux = self.common_flux_wb()?;
        let currents: Vec<_> = self
            .flux_slots
            .iter()
            .enumerate()
            .map(|(i, &slot)| self.gap.current(i, self.circuit.fluxes_wb[slot], flux))
            .collect();
        if currents.iter().any(|i| !i.is_finite()) {
            return Err(invalid("shared-gap winding current is unrepresentable"));
        }
        Ok(currents)
    }
    /// Common-gap, reversible core and all declared leakage energy [J], counted
    /// exactly once. A nonlinear core uses its integrated source potential.
    pub fn magnetic_energy_j(&self) -> Result<f64, MaterialCircuitError> {
        let energy = self
            .gap
            .energy(self.gap_m(), |i| self.circuit.fluxes_wb[self.flux_slots[i]]);
        if !energy.is_finite() || energy < 0.0 {
            return Err(invalid("shared-gap magnetic energy is unrepresentable"));
        }
        Ok(energy)
    }

    /// Bias differential inductance `dLambda / di` [H], a symmetric row-major
    /// matrix in `gap().windings()` order. Core saturation uses source `dH/dB`,
    /// never secant `H/B`; each diagonal includes its explicit physical leakage.
    /// At a source knee the material's right-hand magnitude slope is retained.
    pub fn differential_inductance_h(&self) -> Result<Vec<f64>, MaterialCircuitError> {
        let flux = self.common_flux_wb()?;
        let core_reluctance = self
            .gap
            .core
            .as_ref()
            .map_or(0.0, |core| core.differential_reluctance(flux));
        let permeance = 1.0 / (self.gap_m() / self.gap.scale + core_reluctance);
        if !positive(permeance) {
            return Err(invalid(
                "shared-gap differential permeance is unrepresentable",
            ));
        }
        let windings = &self.gap.windings;
        let mut matrix = vec![0.0; windings.len() * windings.len()];
        for (i, a) in windings.iter().enumerate() {
            for (j, b) in windings[..=i].iter().enumerate() {
                let value = permeance * f64::from(a.signed_turns) * f64::from(b.signed_turns)
                    + if i == j { a.leakage_h } else { 0.0 };
                if !value.is_finite() || value == 0.0 {
                    return Err(invalid(
                        "shared-gap differential inductance is unrepresentable",
                    ));
                }
                matrix[i * windings.len() + j] = value;
                matrix[j * windings.len() + i] = value;
            }
        }
        Ok(matrix)
    }
    /// Fixed-linkage force `-Phi^2/(2 mu A)` [N]. Leakage produces no force on
    /// this coordinate because its declared inductances are gap-independent.
    pub fn field_force_n(&self) -> Result<f64, MaterialCircuitError> {
        let flux = self.common_flux_wb()?;
        let force = -0.5 * flux * (flux / self.gap.scale);
        if !force.is_finite() || (flux != 0.0 && force == 0.0) {
            return Err(invalid("shared-gap magnetic force is unrepresentable"));
        }
        Ok(force)
    }

    /// Advance all winding linkages, gap and momentum in ONE native descriptor
    /// step. Material losses, electrical work and mechanical displacement work
    /// use that same solve; no separate motion update or force gain is added.
    pub fn step_candidate(
        &self,
        sources: Vec<MaterialResistor>,
        input: FieldMotionInput<'_>,
        controls: MaterialCircuitSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CoupledGapCandidate, MaterialCircuitError> {
        poll(&mut cancelled)?;
        if !self.circuit.contacts.is_empty() {
            return Err(invalid(
                "contact-bound shared-gap motion requires refreshed contact sources",
            ));
        }
        self.step_with_contacts_candidate(sources, Vec::new(), input, controls, cancelled)
    }

    /// Refresh the existing bulk/contact sources while retaining every interface
    /// identity and branch. Failed solves, gap limits, incompatible endpoint
    /// clamps, imbalance and cancellation leave the complete model unchanged.
    pub fn step_with_contacts_candidate(
        &self,
        sources: Vec<MaterialResistor>,
        mut contacts: Vec<SeriesElectricalContact>,
        input: FieldMotionInput<'_>,
        controls: MaterialCircuitSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CoupledGapCandidate, MaterialCircuitError> {
        poll(&mut cancelled)?;
        contacts::admit(&mut contacts)?;
        if contacts.len() != self.circuit.contacts.len()
            || contacts
                .iter()
                .zip(&self.circuit.contacts)
                .any(|(a, b)| (a.id(), a.branch()) != (b.id(), b.branch()))
        {
            return Err(invalid(
                "shared-gap contact identities changed inside a physical step",
            ));
        }
        let dt = input.duration_s;
        if !positive(dt)
            || !input.external_force_n.is_finite()
            || !positive(controls.energy_tolerance_j)
            || !positive(controls.power_tolerance_w)
        {
            return Err(invalid(
                "shared-gap step requires finite force and positive duration/work/power tolerances",
            ));
        }
        let (graph, resistors, dae) = bind_with_mutual(
            &self.circuit.graph,
            sources,
            &contacts,
            &self.circuit.mutual,
        )?;
        let old_x = extended(
            &dae,
            &self.circuit.fluxes_wb,
            &self.circuit.charges_c,
            input.electrical_ports,
        )?;
        let (system, law) = self.assemble(&dae)?;
        let nd = dae.system.differential_dim();
        let map = |i: usize| if i < nd { i } else { i + 2 };
        let mut x0 = vec![0.0; system.state_dim()];
        for (i, &value) in old_x.iter().enumerate() {
            x0[map(i)] = value;
        }
        x0[nd] = self.gap_m();
        x0[nd + 1] = self.momentum_kg_m_s;
        let mut ports = input.electrical_ports.to_vec();
        ports.push(input.external_force_n);
        poll(&mut cancelled)?;
        let record = step_descriptor(&system, &x0, &ports, dt).map_err(phs)?;
        let gap = self.gap.at_gap(record.x[nd])?;
        let electrical = [law.electrical_energy(&x0), law.electrical_energy(&record.x)];
        let kinetic = [law.kinetic_energy(&x0), law.kinetic_energy(&record.x)];
        let elastic = [law.elastic_energy(&x0), law.elastic_energy(&record.x)];
        let nodes: Vec<_> = dae.node_potential_index.iter().map(|&i| map(i)).collect();
        let (heat, contact_heat, heat_j, residual) = account_step(
            &graph,
            &resistors,
            &contacts,
            &nodes,
            &x0,
            &record,
            electrical[0] + kinetic[0] + elastic[0],
            electrical[1] + kinetic[1] + elastic[1],
            dt,
            controls,
            &mut cancelled,
        )?;
        let electrical_supply = dt
            * input
                .electrical_ports
                .iter()
                .zip(&record.y)
                .map(|(&u, &y)| u * y)
                .sum::<f64>();
        let mechanical_supply = input.external_force_n * (gap.gap_m() - self.gap_m());
        let native_mechanical_supply = dt * input.external_force_n * record.y[dae.ports];
        let field_work = kinetic[1] - kinetic[0] + elastic[1] - elastic[0] - mechanical_supply;
        let electrical_work = electrical_supply - heat_j - (electrical[1] - electrical[0]);
        let total_residual = electrical[1] - electrical[0] + kinetic[1] - kinetic[0] + elastic[1]
            - elastic[0]
            + heat_j
            - electrical_supply
            - mechanical_supply;
        if [
            electrical_supply,
            mechanical_supply,
            native_mechanical_supply,
            field_work,
            electrical_work,
            total_residual,
        ]
        .iter()
        .any(|value| !value.is_finite())
            || (mechanical_supply - native_mechanical_supply).abs() > controls.energy_tolerance_j
            || (field_work - electrical_work).abs() > controls.energy_tolerance_j
            || (total_residual - residual).abs() > controls.energy_tolerance_j
            || total_residual.abs() > controls.energy_tolerance_j
        {
            return Err(MaterialCircuitError::Balance(
                "shared-gap electrical, magnetic, mechanical work and branch heat disagree",
            ));
        }
        let permeance = gap.secant_permeance(gap.gap_m(), |i| {
            record.x[dae.flux_index[self.flux_slots[i]]]
        })?;
        let (graph, mutual) = gap.lower_at_permeance(&graph, permeance)?;
        poll(&mut cancelled)?;
        let circuit = MaterialCircuit::new_with_contacts_and_mutual_inductance(
            graph,
            resistors,
            contacts,
            mutual,
            dae.flux_index.iter().map(|&i| record.x[i]).collect(),
            dae.charge_index.iter().map(|&i| record.x[i]).collect(),
            input.electrical_ports,
        )?;
        let model = Self {
            circuit,
            gap,
            body: self.body.clone(),
            support: self.support.clone(),
            flux_slots: self.flux_slots.clone(),
            momentum_kg_m_s: record.x[nd + 1],
        };
        model.magnetic_energy_j()?;
        model.field_force_n()?;
        poll(&mut cancelled)?;
        Ok(CoupledGapCandidate {
            model,
            step: FieldMotionStep {
                native: record,
                electrical_energy_j: electrical,
                kinetic_energy_j: kinetic,
                elastic_energy_j: elastic,
                electrical_supply_j: electrical_supply,
                mechanical_supply_j: mechanical_supply,
                field_work_on_body_j: field_work,
                resistors: heat,
                contacts: contact_heat,
                heat_j,
                energy_residual_j: total_residual,
            },
        })
    }

    fn assemble(
        &self,
        dae: &CircuitDae,
    ) -> Result<(DescriptorPortHamiltonian, CoupledGapStorage), MaterialCircuitError> {
        let old_n = dae.system.state_dim();
        let old_d = dae.system.differential_dim();
        let (n, nd, m) = (old_n + 2, old_d + 2, dae.ports + 1);
        if n > MAX_STATES || m > MAX_PORTS {
            return Err(invalid(
                "shared-gap motion exceeds bounded descriptor capacity",
            ));
        }
        let indices: Vec<_> = self
            .flux_slots
            .iter()
            .map(|&slot| dae.flux_index[slot])
            .collect();
        let mut common = vec![false; old_d];
        for &index in &indices {
            common[index] = true;
        }
        let old_q = dae
            .system
            .storage()
            .exact_quadratic_matrix()
            .ok_or_else(|| {
                invalid("shared-gap coupling requires original quadratic circuit storage")
            })?;
        if old_q.len() != old_n * old_n {
            return Err(invalid(
                "shared-gap original storage has the wrong dimension",
            ));
        }
        let mut q = vec![0.0; nd * nd];
        for i in 0..old_d {
            for k in 0..old_d {
                // Replace the ENTIRE common winding block, including mutual
                // cross terms, exactly once. Retain every unrelated C/L store.
                if !common[i] && !common[k] {
                    q[i * nd + k] = old_q[i * old_n + k];
                }
            }
        }
        let law = CoupledGapStorage {
            other: QuadraticStorage::new(q, nd).map_err(phs)?,
            gap: self.gap.clone(),
            indices,
            gap_index: old_d,
            momentum_index: old_d + 1,
            mass: self.generalized_mass_kg(),
            support: self.support.clone(),
        };
        let (old_j, old_r, old_g) = dae.system.structure();
        let (mut j, mut r, mut g) = (vec![0.0; n * n], vec![0.0; n * n], vec![0.0; n * m]);
        let map = |i: usize| if i < old_d { i } else { i + 2 };
        for i in 0..old_n {
            for k in 0..old_n {
                j[map(i) * n + map(k)] = old_j[i * old_n + k];
                r[map(i) * n + map(k)] = old_r[i * old_n + k];
            }
            for port in 0..dae.ports {
                g[map(i) * m + port] = old_g[i * dae.ports + port];
            }
        }
        j[old_d * n + old_d + 1] = 1.0;
        j[(old_d + 1) * n + old_d] = -1.0;
        g[(old_d + 1) * m + dae.ports] = 1.0;
        let system = DescriptorPortHamiltonian::new(n, nd, m, j, r, g, Box::new(law.clone()))
            .map_err(phs)?;
        Ok((system, law))
    }
}

#[derive(Clone)]
struct CoupledGapStorage {
    other: QuadraticStorage,
    gap: CommonMagneticGap,
    indices: Vec<usize>,
    gap_index: usize,
    momentum_index: usize,
    mass: f64,
    support: Option<AxialSupport>,
}
impl CoupledGapStorage {
    fn flux(&self, x: &[f64]) -> f64 {
        self.gap.flux(x[self.gap_index], |i| x[self.indices[i]])
    }
    fn electrical_energy(&self, x: &[f64]) -> f64 {
        self.other.hamiltonian(x) + self.gap.energy(x[self.gap_index], |i| x[self.indices[i]])
    }
    fn kinetic_energy(&self, x: &[f64]) -> f64 {
        0.5 * x[self.momentum_index] * (x[self.momentum_index] / self.mass)
    }
    fn elastic_energy(&self, x: &[f64]) -> f64 {
        self.support.as_ref().map_or(0.0, |support| {
            support
                .force_and_energy(x[self.gap_index])
                .map_or(f64::NAN, |(_, energy)| energy)
        })
    }
}
impl Storage for CoupledGapStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        if !within(x[self.gap_index], self.gap.limits_m) {
            return f64::NAN;
        }
        self.electrical_energy(x) + self.kinetic_energy(x) + self.elastic_energy(x)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        if !within(x[self.gap_index], self.gap.limits_m) {
            out.fill(f64::NAN);
            return;
        }
        self.other.gradient(x, out);
        let flux = self.flux(x);
        for (winding, &index) in self.indices.iter().enumerate() {
            out[index] += self.gap.current(winding, x[index], flux);
        }
        out[self.gap_index] += 0.5 * flux * (flux / self.gap.scale);
        out[self.momentum_index] += x[self.momentum_index] / self.mass;
        if let Some(support) = &self.support {
            match support.force_and_energy(x[self.gap_index]) {
                Ok((force, _)) => out[self.gap_index] -= force,
                Err(_) => out.fill(f64::NAN),
            }
        }
    }
    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        if [a, b, mid]
            .iter()
            .any(|x| !within(x[self.gap_index], self.gap.limits_m))
        {
            return Some(f64::NAN);
        }
        let (fa, fb, fm) = (self.flux(a), self.flux(b), self.flux(mid));
        let delta_gap = b[self.gap_index] - a[self.gap_index];
        // Exact two-endpoint identity:
        // dH = avg(i_a,i_b).dLambda + Phi_a Phi_b dg/(2 mu A)
        //      + dW_core - avg(core_mmf_a,core_mmf_b)*dPhi.
        // Subtract midpoint efforts in factored form; never difference large H.
        let mut remainder = 0.5 * ((fa - fm) * fb + fm * (fb - fm)) * (delta_gap / self.gap.scale);
        for (winding, &index) in self.indices.iter().enumerate() {
            let ia = self.gap.current(winding, a[index], fa);
            let ib = self.gap.current(winding, b[index], fb);
            let im = self.gap.current(winding, mid[index], fm);
            remainder += 0.5 * ((ia - im) + (ib - im)) * (b[index] - a[index]);
        }
        if let Some(core) = &self.gap.core {
            remainder += core.trapezoid_remainder(fa, fb);
        }
        remainder += self.other.midpoint_energy_remainder(a, b, mid)?;
        let p = self.momentum_index;
        remainder += 0.5 * ((a[p] - mid[p]) + (b[p] - mid[p])) * ((b[p] - a[p]) / self.mass);
        if let Some(support) = &self.support {
            match (
                support.force_and_energy(a[self.gap_index]),
                support.force_and_energy(b[self.gap_index]),
                support.force_and_energy(mid[self.gap_index]),
            ) {
                (Ok((fa, _)), Ok((fb, _)), Ok((fm, _))) => {
                    remainder -= 0.5 * ((fa - fm) + (fb - fm)) * delta_gap;
                }
                _ => return Some(f64::NAN),
            }
        }
        Some(remainder)
    }
}
