//! One reversible dielectric potential drives charge, gap force and motion.
//!
//! Fixed area A and the source's w(D)=integral E dD give H=A*d*w(q/A),
//! V=d*E(q/A), and force -A*w(q/A) at fixed charge. The native descriptor
//! advances all electrical coordinates, gap and source-derived body momentum.
//! This uniform gap-medium model supplies no electrostriction, material flow,
//! fringing, contact, breakdown, memory or thermal field-entropy evolution.

use fs_material::field_storage::dielectric::reversible::{
    ELECTRIC_FIELD_DIMS, REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY, ReversibleDielectricMaterial,
};
use fs_mbd_material::ResolvedBody;
use fs_phs::circuit::capacitance::{CapacitorChargeCurve, CapacitorCurve};
use fs_phs::circuit::{Branch, CircuitDae, CircuitGraph, MutualInductance};
use fs_phs::{DescriptorPortHamiltonian, QuadraticStorage, StepRecord, Storage, step_descriptor};

use super::super::contacts::{self, ContactHeat, SeriesElectricalContact};
use super::super::{
    MAX_STATES, MaterialCircuit, MaterialCircuitError, MaterialCircuitSolve, MaterialResistor,
    ResistorHeat, account_step_with_storage_change, bind_with_storage_laws, extended, poll,
};
use super::support::AxialSupport;
use super::{FieldMotionInput, invalid, phs, positive, within};

/// Original constitutive source and fixed-area geometry of one moving branch.
#[derive(Clone, Debug, PartialEq)]
pub struct ReversibleDielectricGap {
    material: ReversibleDielectricMaterial,
    branch: usize,
    area_m2: f64,
    gap_limits_m: [f64; 2],
    maximum_charge_c: f64,
    maximum_energy_j: f64,
}

impl ReversibleDielectricGap {
    /// Bind a source-independent reversible field law to a uniform gap.
    /// Actual dependence on T or another external coordinate refuses: freezing
    /// a free-energy snapshot does not supply its missing reversible heat.
    pub fn new(
        material: ReversibleDielectricMaterial,
        branch: usize,
        area_m2: f64,
        gap_limits_m: [f64; 2],
    ) -> Result<Self, MaterialCircuitError> {
        let axes = material
            .material()
            .property(REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY)
            .ok_or_else(|| invalid("moving dielectric lost its constitutive source"))?
            .source_axes();
        if axes.len() != 1
            || axes[0].0 != material.source_path().axis
            || axes[0].1 != ELECTRIC_FIELD_DIMS
        {
            return Err(invalid(
                "closed dielectric motion requires a field-only source potential",
            ));
        }
        let maximum_charge_c = area_m2 * material.maximum_displacement_c_m2();
        let response = material.evaluate_displacement(material.maximum_displacement_c_m2())?;
        let maximum_force_n = area_m2 * response.energy_density_j_m3;
        let maximum_energy_j = gap_limits_m[1] * maximum_force_n;
        if branch >= super::super::MAX_MATERIAL_CIRCUIT_BRANCHES
            || [
                area_m2,
                gap_limits_m[0],
                gap_limits_m[1],
                maximum_charge_c,
                maximum_force_n,
                maximum_energy_j,
                gap_limits_m[1] * response.field_strength_v_m,
            ]
            .iter()
            .any(|&v| !positive(v) || !positive(v.recip()))
            || gap_limits_m[0] >= gap_limits_m[1]
        {
            return Err(invalid(
                "bounded dielectric branch, area, charge and open gap domain required",
            ));
        }
        Ok(Self {
            material,
            branch,
            area_m2,
            gap_limits_m,
            maximum_charge_c,
            maximum_energy_j,
        })
    }
    pub fn material(&self) -> &ReversibleDielectricMaterial {
        &self.material
    }
    pub const fn branch(&self) -> usize {
        self.branch
    }
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    pub const fn gap_limits_m(&self) -> [f64; 2] {
        self.gap_limits_m
    }
    pub const fn maximum_charge_c(&self) -> f64 {
        self.maximum_charge_c
    }

    fn displacement(&self, charge: f64) -> Result<f64, MaterialCircuitError> {
        if !charge.is_finite() || charge.abs() > self.maximum_charge_c {
            return Err(invalid(
                "dielectric charge is outside the original source support",
            ));
        }
        // Admission was in the original q=A*Dmax coordinate. Division at that
        // supported boundary can round upward; no outside charge is projected.
        let displacement = (charge.abs() / self.area_m2)
            .min(self.material.maximum_displacement_c_m2())
            .copysign(charge);
        if charge != 0.0 && displacement == 0.0 {
            return Err(invalid(
                "nonzero dielectric displacement is unrepresentable",
            ));
        }
        Ok(displacement)
    }
    fn density(&self, charge: f64) -> Result<(f64, f64), MaterialCircuitError> {
        let displacement = self.displacement(charge)?;
        let response = self.material.evaluate_displacement(displacement.abs())?;
        Ok((
            response.field_strength_v_m.copysign(displacement),
            response.energy_density_j_m3,
        ))
    }
}

/// Explicit initial physical coordinates; no charge, flux or body speed is inferred.
pub struct ReversibleDielectricMotionInitial {
    pub fluxes_wb: Vec<f64>,
    pub charges_c: Vec<f64>,
    pub gap_m: f64,
    pub velocity_m_s: f64,
    pub electrical_inputs: Vec<f64>,
}

/// Accepted joint charge/flux, material, gap and momentum history.
#[derive(Clone, Debug, PartialEq)]
pub struct ReversibleDielectricGapMotion {
    routing: MaterialCircuit,
    dielectric: ReversibleDielectricGap,
    body: ResolvedBody,
    support: Option<AxialSupport>,
    charge_slot: usize,
    gap_m: f64,
    momentum_kg_m_s: f64,
}

/// Independently audited physical energy changes; endpoint arrays are diagnostics.
pub struct ReversibleDielectricMotionStep {
    pub native: StepRecord,
    pub electrical_energy_j: [f64; 2],
    pub kinetic_energy_j: [f64; 2],
    pub elastic_energy_j: [f64; 2],
    pub electrical_energy_change_j: f64,
    pub kinetic_energy_change_j: f64,
    pub elastic_energy_change_j: f64,
    pub electrical_supply_j: f64,
    pub mechanical_supply_j: f64,
    pub field_work_on_body_j: f64,
    pub resistors: Vec<ResistorHeat>,
    pub contacts: Vec<ContactHeat>,
    pub heat_j: f64,
    pub energy_residual_j: f64,
}

/// A proposal publishes electrical and mechanical coordinates together, or neither.
pub struct ReversibleDielectricMotionCandidate {
    model: ReversibleDielectricGapMotion,
    step: ReversibleDielectricMotionStep,
}
impl ReversibleDielectricMotionCandidate {
    pub const fn model(&self) -> &ReversibleDielectricGapMotion {
        &self.model
    }
    pub const fn step(&self) -> &ReversibleDielectricMotionStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        ReversibleDielectricGapMotion,
        ReversibleDielectricMotionStep,
    ) {
        (self.model, self.step)
    }
}

fn finite(value: f64) -> Result<f64, MaterialCircuitError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("dielectric motion arithmetic is unrepresentable"))
    }
}
fn sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let (mut total, mut tail) = (0.0_f64, 0.0_f64);
    for value in values {
        let next = total + value;
        tail += if total.abs() >= value.abs() {
            (total - next) + value
        } else {
            (value - next) + total
        };
        total = next;
    }
    total + tail
}
fn offset(a: f64, b: f64, midpoint: f64) -> f64 {
    0.5 * (a - midpoint) + 0.5 * (b - midpoint)
}

impl ReversibleDielectricGapMotion {
    /// Assemble a new physical system. The active graph capacitance is only a
    /// positive routing placeholder; the retained material supplies ALL its
    /// energy and voltage. This is not a zero-time replacement of live storage.
    /// Other authored curves, mutual inductance, transformer/source topology,
    /// bulk resistors and contacts retain their native owners and ordering.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        graph: CircuitGraph,
        sources: Vec<MaterialResistor>,
        mut contacts: Vec<SeriesElectricalContact>,
        mutual: Vec<MutualInductance>,
        mut curves: Vec<CapacitorCurve>,
        dielectric: ReversibleDielectricGap,
        body: ResolvedBody,
        support: Option<AxialSupport>,
        initial: ReversibleDielectricMotionInitial,
    ) -> Result<Self, MaterialCircuitError> {
        if !graph
            .branches
            .get(dielectric.branch)
            .is_some_and(|b| matches!(b.2, Branch::Capacitor { .. }))
            || curves.len() > super::super::MAX_MATERIAL_CIRCUIT_BRANCHES
            || curves.iter().any(|curve| curve.branch == dielectric.branch)
            || !within(initial.gap_m, dielectric.gap_limits_m)
            || !initial.velocity_m_s.is_finite()
            || support
                .as_ref()
                .is_some_and(AxialSupport::is_thermally_live)
        {
            return Err(invalid(
                "one source-owned capacitor, supported gap and frozen support required",
            ));
        }
        curves.sort_by_key(|curve| curve.branch);
        contacts::admit(&mut contacts)?;
        let (graph, resistors, dae) =
            bind_with_storage_laws(&graph, sources, &contacts, &mutual, &curves)?;
        let charge_slot = graph.branches[..dielectric.branch]
            .iter()
            .filter(|b| matches!(b.2, Branch::Capacitor { .. }))
            .count();
        let mass =
            body.inertia().mass() + support.as_ref().map_or(0.0, AxialSupport::modal_mass_kg);
        let momentum = mass * initial.velocity_m_s;
        if !positive(mass)
            || !positive(mass.recip())
            || !momentum.is_finite()
            || initial.velocity_m_s != 0.0 && momentum == 0.0
        {
            return Err(invalid(
                "source-derived dielectric actuator mass/momentum is unrepresentable",
            ));
        }
        let routing = MaterialCircuit {
            graph,
            mutual,
            magnetic_couplings: Vec::new(),
            capacitor_curves: curves.into(),
            dielectric_capacitors: Vec::new().into(),
            resistors,
            contacts,
            fluxes_wb: initial.fluxes_wb,
            charges_c: initial.charges_c,
            storage: Vec::new(),
        };
        let model = Self {
            routing,
            dielectric,
            body,
            support,
            charge_slot,
            gap_m: initial.gap_m,
            momentum_kg_m_s: momentum,
        };
        let (coupled, _) = model.assemble(&dae)?;
        let state = model.extended(&dae, &coupled, &initial.electrical_inputs)?;
        let mut inputs = initial.electrical_inputs;
        inputs.push(0.0);
        admit(&coupled, &state, &inputs)?;
        Ok(model)
    }
    /// Routing/source view only. Its selected C is a placeholder; advance this
    /// joint owner, never the isolated MaterialCircuit, to obtain actual physics.
    pub const fn routing_circuit(&self) -> &MaterialCircuit {
        &self.routing
    }
    pub const fn dielectric(&self) -> &ReversibleDielectricGap {
        &self.dielectric
    }
    pub const fn body(&self) -> &ResolvedBody {
        &self.body
    }
    pub const fn support(&self) -> Option<&AxialSupport> {
        self.support.as_ref()
    }
    pub const fn gap_m(&self) -> f64 {
        self.gap_m
    }
    pub const fn momentum_kg_m_s(&self) -> f64 {
        self.momentum_kg_m_s
    }
    pub fn generalized_mass_kg(&self) -> f64 {
        self.body.inertia().mass()
            + self
                .support
                .as_ref()
                .map_or(0.0, AxialSupport::modal_mass_kg)
    }
    pub fn velocity_m_s(&self) -> f64 {
        self.momentum_kg_m_s / self.generalized_mass_kg()
    }
    pub fn charge_c(&self) -> f64 {
        self.routing.charges_c[self.charge_slot]
    }
    pub fn voltage_v(&self) -> Result<f64, MaterialCircuitError> {
        finite(self.gap_m * self.dielectric.density(self.charge_c())?.0)
    }
    /// Energy of the source-backed moving dielectric alone [J].
    pub fn dielectric_energy_j(&self) -> Result<f64, MaterialCircuitError> {
        finite(self.dielectric.area_m2 * self.gap_m * self.dielectric.density(self.charge_c())?.1)
    }
    /// Instantaneous virtual-work force at fixed charge; the native discrete
    /// gradient, not a separately sampled force, drives the finite step.
    pub fn field_force_n(&self) -> Result<f64, MaterialCircuitError> {
        finite(-self.dielectric.area_m2 * self.dielectric.density(self.charge_c())?.1)
    }

    pub fn step_candidate(
        &self,
        sources: Vec<MaterialResistor>,
        input: FieldMotionInput<'_>,
        controls: MaterialCircuitSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReversibleDielectricMotionCandidate, MaterialCircuitError> {
        if !self.routing.contacts.is_empty() {
            return Err(invalid(
                "contact-bound dielectric motion needs explicit refreshed contacts",
            ));
        }
        self.step_with_contacts_candidate(sources, Vec::new(), input, controls, cancelled)
    }
    pub fn step_with_contacts_candidate(
        &self,
        sources: Vec<MaterialResistor>,
        mut contacts: Vec<SeriesElectricalContact>,
        input: FieldMotionInput<'_>,
        controls: MaterialCircuitSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ReversibleDielectricMotionCandidate, MaterialCircuitError> {
        poll(&mut cancelled)?;
        contacts::admit(&mut contacts)?;
        if contacts.len() != self.routing.contacts.len()
            || contacts
                .iter()
                .zip(&self.routing.contacts)
                .any(|(a, b)| (a.id(), a.branch()) != (b.id(), b.branch()))
            || !positive(input.duration_s)
            || !input.external_force_n.is_finite()
            || !positive(controls.energy_tolerance_j)
            || !positive(controls.power_tolerance_w)
        {
            return Err(invalid(
                "fixed contact identities and finite positive motion controls required",
            ));
        }
        let (graph, resistors, dae) = bind_with_storage_laws(
            &self.routing.graph,
            sources,
            &contacts,
            &self.routing.mutual,
            &self.routing.capacitor_curves,
        )?;
        let (coupled, law) = self.assemble(&dae)?;
        let x0 = self.extended(&dae, &coupled, input.electrical_ports)?;
        let mut ports = input.electrical_ports.to_vec();
        ports.push(input.external_force_n);
        admit(&coupled, &x0, &ports)?;
        poll(&mut cancelled)?;
        let native =
            step_descriptor(&coupled.system, &x0, &ports, input.duration_s).map_err(phs)?;
        admit(&coupled, &native.x, &ports)?;
        let before = law.energies(&x0);
        let after = law.energies(&native.x);
        let changes = law.changes(&x0, &native.x)?;
        let total_change = finite(sum(changes))?;
        let (resistors_heat, contact_heat, heat_j, residual) = account_step_with_storage_change(
            &graph,
            &resistors,
            &contacts,
            &coupled.node_potential_index,
            &x0,
            &native,
            sum(before),
            sum(after),
            total_change,
            input.duration_s,
            controls,
            &mut cancelled,
        )?;
        let electrical_supply_j = finite(
            input.duration_s
                * sum(input
                    .electrical_ports
                    .iter()
                    .zip(&native.y)
                    .map(|(&u, &y)| u * y)),
        )?;
        let gap = native.x[law.gap];
        let mechanical_supply_j = finite(input.external_force_n * (gap - self.gap_m))?;
        let field_work_on_body_j = finite(sum([changes[1], changes[2], -mechanical_supply_j]))?;
        let electrical_transfer = finite(sum([electrical_supply_j, -heat_j, -changes[0]]))?;
        let total_residual = finite(sum([
            total_change,
            heat_j,
            -electrical_supply_j,
            -mechanical_supply_j,
        ]))?;
        if (mechanical_supply_j - input.duration_s * input.external_force_n * native.y[dae.ports])
            .abs()
            > controls.energy_tolerance_j
            || (field_work_on_body_j - electrical_transfer).abs() > controls.energy_tolerance_j
            || (total_residual - residual).abs() > controls.energy_tolerance_j
            || total_residual.abs() > controls.energy_tolerance_j
        {
            return Err(MaterialCircuitError::Balance(
                "nonlinear dielectric field/body/source work does not close",
            ));
        }
        let mut model = self.clone();
        model.routing.graph = graph;
        model.routing.resistors = resistors;
        model.routing.contacts = contacts;
        model.routing.fluxes_wb = dae.flux_index.iter().map(|&i| native.x[i]).collect();
        model.routing.charges_c = dae.charge_index.iter().map(|&i| native.x[i]).collect();
        model.gap_m = gap;
        model.momentum_kg_m_s = native.x[law.momentum];
        poll(&mut cancelled)?;
        Ok(ReversibleDielectricMotionCandidate {
            model,
            step: ReversibleDielectricMotionStep {
                native,
                electrical_energy_j: [before[0], after[0]],
                kinetic_energy_j: [before[1], after[1]],
                elastic_energy_j: [before[2], after[2]],
                electrical_energy_change_j: changes[0],
                kinetic_energy_change_j: changes[1],
                elastic_energy_change_j: changes[2],
                electrical_supply_j,
                mechanical_supply_j,
                field_work_on_body_j,
                resistors: resistors_heat,
                contacts: contact_heat,
                heat_j,
                energy_residual_j: total_residual,
            },
        })
    }

    fn extended(
        &self,
        original: &CircuitDae,
        coupled: &CircuitDae,
        ports: &[f64],
    ) -> Result<Vec<f64>, MaterialCircuitError> {
        let old = extended(
            original,
            &self.routing.fluxes_wb,
            &self.routing.charges_c,
            ports,
        )?;
        let nd = original.system.differential_dim();
        let mut state = vec![0.0; coupled.system.state_dim()];
        for (i, value) in old.into_iter().enumerate() {
            state[if i < nd { i } else { i + 2 }] = value;
        }
        state[nd] = self.gap_m;
        state[nd + 1] = self.momentum_kg_m_s;
        Ok(state)
    }
    fn assemble(
        &self,
        original: &CircuitDae,
    ) -> Result<(CircuitDae, MotionStorage), MaterialCircuitError> {
        let old_n = original.system.state_dim();
        let old_d = original.system.differential_dim();
        let (n, nd, m) = (old_n + 2, old_d + 2, original.ports + 1);
        if n > MAX_STATES || m > super::super::MAX_PORTS {
            return Err(invalid("dielectric motion exceeds native capacity"));
        }
        let charge = *original
            .charge_index
            .get(self.charge_slot)
            .ok_or_else(|| invalid("dielectric charge absent"))?;
        let part = original
            .system
            .storage()
            .quadratic_part()
            .ok_or_else(|| invalid("native circuit storage decomposition absent"))?;
        if part.matrix.len() != old_n * old_n {
            return Err(invalid("native quadratic storage dimensions changed"));
        }
        let mut q = vec![0.0; nd * nd];
        for i in 0..old_d {
            for k in 0..old_d {
                q[i * nd + k] = part.matrix[i * old_n + k];
            }
        }
        q[charge * nd + charge] = 0.0;
        let mut curves = Vec::with_capacity(self.routing.capacitor_curves.len());
        for binding in self.routing.capacitor_curves.iter() {
            let slot = self.routing.graph.branches[..binding.branch]
                .iter()
                .filter(|b| matches!(b.2, Branch::Capacitor { .. }))
                .count();
            curves.push((original.charge_index[slot], binding.curve.clone()));
        }
        let mass = self.generalized_mass_kg();
        let momentum_scale = fs_math::det::sqrt(mass * self.dielectric.maximum_energy_j);
        if !positive(momentum_scale) {
            return Err(invalid("dielectric mechanical probe scale unrepresentable"));
        }
        let law = MotionStorage {
            other: QuadraticStorage::new(q, nd).map_err(phs)?,
            curves,
            dielectric: self.dielectric.clone(),
            charge,
            gap: old_d,
            momentum: old_d + 1,
            mass,
            momentum_scale,
            support: self.support.clone(),
        };
        let (old_j, old_r, old_g) = original.system.structure();
        let (mut j, mut r, mut g) = (vec![0.0; n * n], vec![0.0; n * n], vec![0.0; n * m]);
        let map = |i| if i < old_d { i } else { i + 2 };
        for i in 0..old_n {
            for k in 0..old_n {
                j[map(i) * n + map(k)] = old_j[i * old_n + k];
                r[map(i) * n + map(k)] = old_r[i * old_n + k];
            }
            for p in 0..original.ports {
                g[map(i) * m + p] = old_g[i * original.ports + p];
            }
        }
        j[old_d * n + old_d + 1] = 1.0;
        j[(old_d + 1) * n + old_d] = -1.0;
        g[(old_d + 1) * m + original.ports] = 1.0;
        let system = DescriptorPortHamiltonian::new(n, nd, m, j, r, g, Box::new(law.clone()))
            .map_err(phs)?;
        Ok((
            CircuitDae {
                system,
                flux_index: original.flux_index.clone(),
                charge_index: original.charge_index.clone(),
                node_potential_index: original
                    .node_potential_index
                    .iter()
                    .copied()
                    .map(map)
                    .collect(),
                ports: m,
            },
            law,
        ))
    }
}

fn admit(dae: &CircuitDae, state: &[f64], inputs: &[f64]) -> Result<(), MaterialCircuitError> {
    if !dae.system.hamiltonian(state).is_finite() {
        return Err(invalid("dielectric motion left its source/geometry domain"));
    }
    let consistent = dae.consistent_initial_state(state, inputs)?;
    for i in 0..dae.system.differential_dim() {
        if consistent[i].to_bits() != state[i].to_bits() {
            return Err(invalid(
                "nonlinear dielectric constraints would change a conserved coordinate",
            ));
        }
    }
    Ok(())
}

#[derive(Clone)]
struct MotionStorage {
    other: QuadraticStorage,
    curves: Vec<(usize, CapacitorChargeCurve)>,
    dielectric: ReversibleDielectricGap,
    charge: usize,
    gap: usize,
    momentum: usize,
    mass: f64,
    momentum_scale: f64,
    support: Option<AxialSupport>,
}
impl MotionStorage {
    fn energies(&self, x: &[f64]) -> [f64; 3] {
        let electric = self
            .dielectric
            .density(x[self.charge])
            .map_or(f64::NAN, |(_, w)| self.dielectric.area_m2 * x[self.gap] * w);
        [
            sum(std::iter::once(self.other.hamiltonian(x))
                .chain(std::iter::once(electric))
                .chain(
                    self.curves
                        .iter()
                        .map(|(i, c)| c.stored_energy_j(x[*i]).unwrap_or(f64::NAN)),
                )),
            0.5 * x[self.momentum] * (x[self.momentum] / self.mass),
            self.support.as_ref().map_or(0.0, |s| {
                s.force_and_energy(x[self.gap]).map_or(f64::NAN, |(_, e)| e)
            }),
        ]
    }
    fn changes(&self, a: &[f64], b: &[f64]) -> Result<[f64; 3], MaterialCircuitError> {
        let n = self.other.matrix().len().isqrt();
        let mut mid = [0.0; MAX_STATES];
        let mut gradient = [0.0; MAX_STATES];
        for i in 0..n {
            mid[i] = f64::midpoint(a[i], b[i]);
        }
        self.other.gradient(&mid[..n], &mut gradient[..n]);
        let quadratic = sum(std::iter::once(
            self.other
                .midpoint_energy_remainder(&a[..n], &b[..n], &mid[..n])
                .unwrap(),
        )
        .chain((0..n).map(|i| gradient[i] * (b[i] - a[i]))));
        let da = self.dielectric.displacement(a[self.charge])?;
        let db = self.dielectric.displacement(b[self.charge])?;
        let wa = self
            .dielectric
            .material
            .evaluate_displacement(da.abs())?
            .energy_density_j_m3;
        let wb = self
            .dielectric
            .material
            .evaluate_displacement(db.abs())?
            .energy_density_j_m3;
        let dw = self
            .dielectric
            .material
            .energy_change_density_j_m3(da, db)?;
        let dg = b[self.gap] - a[self.gap];
        let gm = mid[self.gap];
        let wm = f64::midpoint(wa, wb);
        let active = self.dielectric.area_m2
            * sum([
                gm * dw,
                offset(a[self.gap], b[self.gap], gm) * dw,
                wm * dg,
                offset(wa, wb, wm) * dg,
            ]);
        let other_curves = sum(self
            .curves
            .iter()
            .map(|(i, c)| c.energy_change_j(a[*i], b[*i]).unwrap_or(f64::NAN)));
        let dp = b[self.momentum] - a[self.momentum];
        let pm = mid[self.momentum];
        let kinetic = sum([
            dp * (pm / self.mass),
            dp * (offset(a[self.momentum], b[self.momentum], pm) / self.mass),
        ]);
        let elastic = if let Some(s) = &self.support {
            let force = s.force_and_energy(gm)?.0;
            sum([
                -force * dg,
                s.stiffness_n_m() * offset(a[self.gap], b[self.gap], gm) * dg,
            ])
        } else {
            0.0
        };
        Ok([
            finite(sum([quadratic, other_curves, active]))?,
            finite(kinetic)?,
            finite(elastic)?,
        ])
    }
    fn remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Result<f64, MaterialCircuitError> {
        let da = self.dielectric.displacement(a[self.charge])?;
        let db = self.dielectric.displacement(b[self.charge])?;
        let dm = self.dielectric.displacement(mid[self.charge])?;
        let law = &self.dielectric.material;
        let area = self.dielectric.area_m2;
        let em = law
            .evaluate_displacement(dm.abs())?
            .field_strength_v_m
            .copysign(dm);
        let dw = law.energy_change_density_j_m3(da, db)?;
        let rq = sum([
            area * law.energy_change_remainder_density_j_m3(da, db, dm)?,
            em * area.mul_add(db - da, -(b[self.charge] - a[self.charge])),
        ]);
        let around = sum([
            0.5 * law.energy_change_remainder_density_j_m3(dm, da, dm)?,
            0.5 * law.energy_change_remainder_density_j_m3(dm, db, dm)?,
            em * offset(da, db, dm),
        ]);
        let dg = b[self.gap] - a[self.gap];
        let gap_offset = offset(a[self.gap], b[self.gap], mid[self.gap]);
        let field = sum([
            mid[self.gap] * rq,
            area * gap_offset * dw,
            area * dg * around,
        ]);
        let kinetic = (b[self.momentum] - a[self.momentum])
            * (offset(a[self.momentum], b[self.momentum], mid[self.momentum]) / self.mass);
        let elastic = self
            .support
            .as_ref()
            .map_or(0.0, |s| s.stiffness_n_m() * gap_offset * dg);
        finite(sum(std::iter::once(
            self.other.midpoint_energy_remainder(a, b, mid).unwrap(),
        )
        .chain([field, kinetic, elastic])
        .chain(self.curves.iter().map(|(i, c)| {
            c.energy_change_remainder_j(a[*i], b[*i], mid[*i])
                .unwrap_or(f64::NAN)
        }))))
    }
}
impl Storage for MotionStorage {
    fn coordinate_scale(&self, i: usize) -> Option<f64> {
        if i == self.charge {
            Some(self.dielectric.maximum_charge_c)
        } else if i == self.gap {
            Some(self.dielectric.gap_limits_m[1] - self.dielectric.gap_limits_m[0])
        } else if i == self.momentum {
            Some(self.momentum_scale)
        } else {
            self.curves.iter().find(|(k, _)| *k == i).map(|(_, c)| {
                let (a, b) = c.charge_range_c();
                a.abs().max(b.abs())
            })
        }
    }
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        if !within(x[self.gap], self.dielectric.gap_limits_m) {
            return f64::NAN;
        }
        sum(self.energies(x))
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        if !within(x[self.gap], self.dielectric.gap_limits_m) {
            out.fill(f64::NAN);
            return;
        }
        self.other.gradient(x, out);
        let Ok((e, w)) = self.dielectric.density(x[self.charge]) else {
            out.fill(f64::NAN);
            return;
        };
        out[self.charge] = x[self.gap] * e;
        out[self.gap] = self.dielectric.area_m2 * w;
        out[self.momentum] = x[self.momentum] / self.mass;
        for (i, c) in &self.curves {
            out[*i] = c.voltage_v(x[*i]).unwrap_or(f64::NAN);
        }
        if let Some(s) = &self.support {
            out[self.gap] -= s.force_and_energy(x[self.gap]).map_or(f64::NAN, |(f, _)| f);
        }
    }
    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        Some(self.remainder(a, b, mid).unwrap_or(f64::NAN))
    }
}
