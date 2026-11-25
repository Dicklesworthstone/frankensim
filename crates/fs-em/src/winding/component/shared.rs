//! The owner-independent core of fs-em's winding Components: the ports of a
//! native graph and its heat receivers, the admission of a step, the heat of
//! each receiver, the graph resistors' laws, the magnetic lane's terms, the
//! receipt's closure and the identity words. A mesh element is an index the
//! Component gives meaning to. A refusal here is a [`Refusal`], which each
//! Component turns into its own error, in its own words for its elements.

use core::fmt;
use core::num::NonZeroUsize;
use std::collections::BTreeMap;
use std::sync::Arc;

use fs_material::conductor::ResolvedConductorState;
use fs_material::field_storage::linear_magnetic_permeability;
use fs_material::field_storage::thermal::FieldThermalState;
use fs_material::field_storage::thermal::magnetic::{
    ResolvedThermalMagneticMaterial, ThermalMagneticMaterial,
};
#[cfg(feature = "magnetic-hysteresis")]
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_math::eft::two_sum;
use fs_phs::arithmetic_rounding_floor;
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_ports::component::{ContractViolation, PortLayout, PortRole, PortSpec, Receipt};
use fs_ports::{CoupleError, PortKind, StableId};
#[cfg(feature = "magnetic-hysteresis")]
use fs_qty::{
    QuantitySpec, Temperature,
    semantic::{QuantityKind, ValueForm},
};
use fs_solver::{NewtonKrylovConfig, NewtonReport};

use super::ResistorLaw;
use crate::winding::nonlinear::resistors::GraphResistorMaterial;

/// A refusal of the shared core.
#[derive(Debug)]
pub(crate) enum Refusal {
    /// A shape, a value or a construction that is refused, in words that
    /// name no mesh element.
    Invalid(&'static str),
    /// A refusal that names a mesh element, which the Component words.
    Element(ElementRefusal),
    /// Invalid or repeated port metadata.
    Port(ContractViolation),
    /// An invalid port identifier.
    Identifier(CoupleError),
    /// Mechanical work, or magnetic parameter work outside the magnetic lane
    /// or in an element without a law.
    UnaccountedWork {
        /// Signed magnetic parameter work \[J\].
        parameter_j: f64,
        /// Signed mechanical work \[J\].
        mechanical_j: f64,
    },
    /// The energy account does not close.
    Balance {
        /// Supplied work minus stored change minus heat out minus the
        /// owner's numerical dissipation \[J\].
        defect_j: f64,
        /// The allowance \[J\].
        tolerance_j: f64,
    },
    /// A law had no answer at the temperature its port received, or a
    /// magnetic law's free energy is not the field's.
    Law {
        /// The heat receiver whose port gave the temperature.
        receiver: usize,
        /// The element asked for.
        element: Element,
        /// The temperature \[K\].
        kelvin: f64,
        /// The law's reason.
        message: String,
    },
}

/// The refusals whose words name a mesh element.
#[derive(Clone, Copy, Debug)]
pub(crate) enum ElementRefusal {
    /// A receiver names an element that carries no conductor.
    NoConductor,
    /// The owner reported heat in an element no receiver maps.
    HeatUnmapped,
    /// The owner's magnetic records omit a law element.
    LawRecordMissing,
    /// The state's field has no such law element.
    NoLawElement,
}

/// The element whose law had no answer.
#[derive(Clone, Copy, Debug)]
pub(crate) enum Element {
    /// A native graph resistor's branch.
    Resistor(usize),
    /// A magnetic law's mesh element.
    Permeability(usize),
}

/// A Component's ports and the receivers' maps.
pub(crate) struct Layout {
    /// The external inputs of the native graph, which are the first ports.
    pub(crate) sources: usize,
    /// Each dissipating mesh element's receiver.
    pub(crate) element_port: BTreeMap<usize, usize>,
    /// Each graph resistor's receiver, by branch.
    pub(crate) resistor_port: BTreeMap<usize, usize>,
    /// `source<p>` for each external input `p`, then `heat<k>` for each
    /// receiver.
    pub(crate) ports: Vec<PortSpec>,
}

/// The ports of `graph`'s external inputs (`EffortIn` for a voltage source,
/// `FlowIn` for a current source), then one heat port per receiver, given as
/// its mesh elements and its graph resistors' branches. Each of `elements`
/// (the elements that dissipate) and each graph resistor maps to exactly one
/// receiver. `resistor_laws` are the laws by branch when laws supply the
/// resistors, of which every resistor or none has one.
pub(crate) fn layout(
    graph: &CircuitGraph,
    elements: &[usize],
    receivers: &[(&[usize], &[usize])],
    resistor_laws: Option<&BTreeMap<usize, Arc<dyn ResistorLaw>>>,
    energy_tolerance_j: f64,
) -> Result<Layout, Refusal> {
    if !(energy_tolerance_j.is_finite() && energy_tolerance_j > 0.0) {
        return Err(Refusal::Invalid(
            "energy tolerance must be positive and finite",
        ));
    }
    let mut roles: Vec<Option<PortRole>> = Vec::new();
    let mut resistors = Vec::new();
    for (branch, &(_, _, element)) in graph.branches.iter().enumerate() {
        let (port, role) = match element {
            Branch::Resistor { .. } => {
                resistors.push(branch);
                continue;
            }
            Branch::VoltageSource { port } => (port, PortRole::EffortIn),
            Branch::CurrentSource { port } => (port, PortRole::FlowIn),
            Branch::Inductor { .. } | Branch::Capacitor { .. } => continue,
        };
        if roles.len() <= port {
            roles.resize(port + 1, None);
        }
        if roles[port].replace(role).is_some() {
            return Err(Refusal::Invalid(
                "an external port drives more than one source",
            ));
        }
    }
    let mut element_port = BTreeMap::new();
    let mut resistor_port = BTreeMap::new();
    for (port, &(cells, branches)) in receivers.iter().enumerate() {
        if cells.is_empty() && branches.is_empty() {
            return Err(Refusal::Invalid("a heat receiver maps no element"));
        }
        for &element in cells {
            if !elements.contains(&element) {
                return Err(Refusal::Element(ElementRefusal::NoConductor));
            }
            place(&mut element_port, element, port)?;
        }
        for &branch in branches {
            if !resistors.contains(&branch) {
                return Err(Refusal::Invalid(
                    "a receiver names a branch that is not a resistor",
                ));
            }
            place(&mut resistor_port, branch, port)?;
        }
    }
    if elements
        .iter()
        .any(|element| !element_port.contains_key(element))
        || resistors
            .iter()
            .any(|branch| !resistor_port.contains_key(branch))
    {
        return Err(Refusal::Invalid(
            "a dissipating element is mapped to no receiver",
        ));
    }
    if let Some(laws) = resistor_laws {
        if laws.keys().any(|branch| !resistors.contains(branch)) {
            return Err(Refusal::Invalid(
                "a resistor law names a branch that is not a resistor",
            ));
        }
        // The owner rebinds every graph resistor or none.
        if !laws.is_empty() && resistors.iter().any(|branch| !laws.contains_key(branch)) {
            return Err(Refusal::Invalid("every graph resistor or none has a law"));
        }
    }
    let mut ports = Vec::with_capacity(roles.len() + receivers.len());
    for (port, role) in roles.iter().enumerate() {
        let role = role.ok_or(Refusal::Invalid("an external port drives no source"))?;
        let id = StableId::new(format!("source{port}")).map_err(Refusal::Identifier)?;
        ports.push(
            PortSpec::power(
                id,
                PortKind::ElectricalVoltageCurrent,
                role,
                NonZeroUsize::MIN,
            )
            .map_err(Refusal::Port)?,
        );
    }
    for k in 0..receivers.len() {
        let id = StableId::new(format!("heat{k}")).map_err(Refusal::Identifier)?;
        ports.push(PortSpec::heat(id, PortRole::EffortIn).map_err(Refusal::Port)?);
    }
    PortLayout::new(&ports).map_err(Refusal::Port)?;
    Ok(Layout {
        sources: roles.len(),
        element_port,
        resistor_port,
        ports,
    })
}

fn place(map: &mut BTreeMap<usize, usize>, element: usize, port: usize) -> Result<(), Refusal> {
    if map.insert(element, port).is_some() {
        return Err(Refusal::Invalid(
            "a dissipating element is mapped to two receivers",
        ));
    }
    Ok(())
}

/// A step of `dt_s` with `inputs` and `outputs` values on `sources` source
/// ports and `receivers` heat ports: a positive finite duration, one finite
/// value per port each way, and a positive temperature at each heat port.
pub(crate) fn admit_ports(
    dt_s: f64,
    inputs: &[f64],
    outputs: usize,
    sources: usize,
    receivers: usize,
) -> Result<(), Refusal> {
    let ports = sources + receivers;
    if !(dt_s.is_finite() && dt_s > 0.0)
        || inputs.len() != ports
        || outputs != ports
        || inputs.iter().any(|value| !value.is_finite())
        || inputs[sources..].iter().any(|&kelvin| kelvin <= 0.0)
    {
        return Err(Refusal::Invalid(
            "positive duration, complete finite ports and positive temperatures required",
        ));
    }
    Ok(())
}

/// Whether a state's graph is the Component's: the same nodes, branch ends
/// and transformer terminals, and every value bit for bit (a winding's
/// inductor placeholder may be NaN, which equals nothing by value); when
/// resistor laws set the resistors at every step, all but the resistors'
/// values.
pub(crate) fn same_graph(old: &CircuitGraph, own: &CircuitGraph, resistor_laws: bool) -> bool {
    old.node_count == own.node_count
        && old.transformers.len() == own.transformers.len()
        && old
            .transformers
            .iter()
            .zip(&own.transformers)
            .all(|(x, y)| {
                x.primary == y.primary
                    && x.secondary == y.secondary
                    && x.ratio.to_bits() == y.ratio.to_bits()
            })
        && old.branches.len() == own.branches.len()
        && old.branches.iter().zip(&own.branches).all(|(x, y)| {
            x.0 == y.0
                && x.1 == y.1
                && ((resistor_laws
                    && matches!(
                        (x.2, y.2),
                        (Branch::Resistor { .. }, Branch::Resistor { .. })
                    ))
                    || same_branch(x.2, y.2))
        })
}

/// The same element with the same value, bit for bit.
fn same_branch(a: Branch, b: Branch) -> bool {
    let (mut x, mut y) = (Vec::with_capacity(2), Vec::with_capacity(2));
    branch_words(a, &mut x);
    branch_words(b, &mut y);
    x == y
}

/// A state's temperatures: in the magnetic lane (`lane`, the number of
/// receivers) one positive finite temperature per receiver, and none outside
/// it.
pub(crate) fn admit_temperatures(
    temperatures: Option<&[f64]>,
    lane: Option<usize>,
) -> Result<(), Refusal> {
    if temperatures.map(<[f64]>::len) != lane
        || temperatures
            .into_iter()
            .flatten()
            .any(|&kelvin| !(kelvin.is_finite() && kelvin > 0.0))
    {
        return Err(Refusal::Invalid(
            "a state carries temperatures in the magnetic lane only, one positive \
             temperature per receiver",
        ));
    }
    Ok(())
}

/// The owner's clock at the end of a step of `dt_s` from `time_s`.
pub(crate) fn target_time(time_s: f64, dt_s: f64) -> Result<f64, Refusal> {
    let target_s = time_s + dt_s;
    if !(target_s.is_finite() && target_s > time_s) {
        return Err(Refusal::Invalid(
            "the step does not advance the owner's clock",
        ));
    }
    Ok(target_s)
}

/// Refuses magnetic parameter work outside the magnetic lane (`lane`), and
/// any mechanical work: the receipt has no term for either.
pub(crate) fn admit_work(parameter_j: f64, mechanical_j: f64, lane: bool) -> Result<(), Refusal> {
    if (parameter_j != 0.0 && !lane) || mechanical_j != 0.0 {
        return Err(Refusal::UnaccountedWork {
            parameter_j,
            mechanical_j,
        });
    }
    Ok(())
}

/// The graph resistors at the temperatures (`kelvin`, one per receiver) of
/// the ports their heat leaves through, from their laws.
pub(crate) fn graph_resistors(
    laws: &BTreeMap<usize, Arc<dyn ResistorLaw>>,
    resistor_port: &BTreeMap<usize, usize>,
    kelvin: &[f64],
) -> Result<Vec<GraphResistorMaterial>, Refusal> {
    let mut sources = Vec::with_capacity(laws.len());
    for (&branch, law) in laws {
        let receiver = resistor_port[&branch];
        let kelvin = kelvin[receiver];
        let conductor = law.conductor(kelvin).map_err(|message| Refusal::Law {
            receiver,
            element: Element::Resistor(branch),
            kelvin,
            message,
        })?;
        sources.push(GraphResistorMaterial { branch, conductor });
    }
    Ok(sources)
}

/// The heat \[J\] of each of `receivers` receivers over a step, from the
/// owner's `(element, heat)` records of its stranded and field cells, then
/// its `(branch, heat)` records of the graph resistors, in that order.
pub(crate) fn heat_by_receiver(
    receivers: usize,
    element_port: &BTreeMap<usize, usize>,
    resistor_port: &BTreeMap<usize, usize>,
    cells: impl IntoIterator<Item = (usize, f64)>,
    resistors: impl IntoIterator<Item = (usize, f64)>,
) -> Result<Vec<f64>, Refusal> {
    let mut heat = vec![0.0; receivers];
    for (element, heat_j) in cells {
        let port = element_port
            .get(&element)
            .copied()
            .ok_or(Refusal::Element(ElementRefusal::HeatUnmapped))?;
        heat[port] += heat_j;
    }
    for (branch, heat_j) in resistors {
        let port = resistor_port.get(&branch).copied().ok_or(Refusal::Invalid(
            "the owner reported heat in a resistor no receiver maps",
        ))?;
        heat[port] += heat_j;
    }
    if heat.iter().any(|q| !(q.is_finite() && *q >= 0.0)) {
        return Err(Refusal::Invalid(
            "a receiver's heat is not finite and nonnegative",
        ));
    }
    Ok(heat)
}

/// One mesh element's magnetic record of a step, as the magnetic lane reads
/// it.
pub(crate) struct LawCell {
    /// The mesh element.
    pub(crate) element: usize,
    /// The magnitude of the flux density \[T\] at the start and the end of
    /// the step.
    pub(crate) flux_density_t: [f64; 2],
    /// The field's free energy \[J\] at the start and the end of the step.
    pub(crate) free_energy_j: [f64; 2],
    /// The element's volume \[m³\].
    pub(crate) volume_m3: f64,
    /// The owner's parameter work `psi_new(B_old) - psi_old(B_old)` \[J\].
    pub(crate) parameter_work_j: f64,
}

/// The magnetic lane's terms of one step.
pub(crate) struct Reversible {
    /// The entropy change of each receiver's law elements \[J/K\].
    entropy_j_per_k: Vec<f64>,
    /// The change of `T s` over the law elements \[J\].
    stored_change_j: f64,
    /// The sum of `-(P + dT s0)` over the law elements \[J\].
    numerical_loss_j: f64,
    /// The size of the terms above, for the rounding floor \[J\].
    magnitude_j: f64,
}

/// The magnetic lane's terms of a step from the temperatures `before` to the
/// ports' (`kelvin`), one per receiver, read from the owner's magnetic
/// records `cells`, in mesh order. Each law's free energy must be the
/// field's at both ends of the step, within `energy_tolerance_j`; an element
/// without a law must do no parameter work.
pub(crate) fn reversible(
    laws: &BTreeMap<usize, ThermalMagneticMaterial>,
    element_port: &BTreeMap<usize, usize>,
    receivers: usize,
    before: &[f64],
    kelvin: &[f64],
    cells: &[LawCell],
    energy_tolerance_j: f64,
) -> Result<Reversible, Refusal> {
    if laws.keys().any(|&element| element >= cells.len()) {
        return Err(Refusal::Element(ElementRefusal::LawRecordMissing));
    }
    let mut terms = Reversible {
        entropy_j_per_k: vec![0.0; receivers],
        stored_change_j: 0.0,
        numerical_loss_j: 0.0,
        magnitude_j: 0.0,
    };
    for (element, cell) in cells.iter().enumerate() {
        if cell.element != element {
            return Err(Refusal::Invalid(
                "the owner's magnetic records are not in mesh order",
            ));
        }
        let Some(law) = laws.get(&element) else {
            if cell.parameter_work_j != 0.0 {
                return Err(Refusal::UnaccountedWork {
                    parameter_j: cell.parameter_work_j,
                    mechanical_j: 0.0,
                });
            }
            continue;
        };
        let receiver = element_port[&element];
        let end = |end: usize, kelvin: f64| -> Result<FieldThermalState, Refusal> {
            let refuse = |message| Refusal::Law {
                receiver,
                element: Element::Permeability(element),
                kelvin,
                message,
            };
            let state =
                law_state(law, kelvin, cell.flux_density_t[end], cell.volume_m3).map_err(refuse)?;
            let field_j = cell.free_energy_j[end];
            let gap_j = (state.free_energy_j - field_j).abs();
            if gap_j.is_nan() || gap_j > energy_tolerance_j {
                return Err(refuse(format!(
                    "its free energy {} J is not the field's {field_j} J",
                    state.free_energy_j
                )));
            }
            Ok(state)
        };
        let (t0, t1) = (before[receiver], kelvin[receiver]);
        let (s0, s1) = (end(0, t0)?.entropy_j_k, end(1, t1)?.entropy_j_k);
        terms.entropy_j_per_k[receiver] += s1 - s0;
        let (ts0, ts1) = (t0 * s0, t1 * s1);
        terms.stored_change_j += ts1 - ts0;
        // -(P + dT s0), with P = psi(T1, B0) - psi(T0, B0) = -s0 dT + O(dT^2).
        let loss = -(cell.parameter_work_j + (t1 - t0) * s0);
        terms.numerical_loss_j += loss;
        terms.magnitude_j += ts0.abs() + ts1.abs() + cell.parameter_work_j.abs() + loss.abs();
    }
    Ok(terms)
}

/// A state's `(sum of T s [J], sum of s [J/K])` over the magnetic lane's law
/// elements, measured from the state alone: each law at its receiver's
/// temperature (`temperatures`, one per receiver), at the flux density
/// magnitude and volume `cells` gives each mesh element. Both sums are
/// compensated (`two_sum`).
pub(crate) fn stored_reversible(
    laws: &BTreeMap<usize, ThermalMagneticMaterial>,
    element_port: &BTreeMap<usize, usize>,
    temperatures: &[f64],
    cells: &[(f64, f64)],
) -> Result<(f64, f64), Refusal> {
    let (mut energy, mut entropy) = ([0.0; 2], [0.0; 2]);
    for (&element, law) in laws {
        let receiver = element_port[&element];
        let kelvin = temperatures[receiver];
        let &(b, volume_m3) = cells
            .get(element)
            .ok_or(Refusal::Element(ElementRefusal::NoLawElement))?;
        let s = law_state(law, kelvin, b, volume_m3)
            .map_err(|message| Refusal::Law {
                receiver,
                element: Element::Permeability(element),
                kelvin,
                message,
            })?
            .entropy_j_k;
        compensated(&mut energy, kelvin * s);
        compensated(&mut entropy, s);
    }
    let (energy, entropy) = (energy[0] + energy[1], entropy[0] + entropy[1]);
    if !(energy.is_finite() && entropy.is_finite()) {
        return Err(Refusal::Invalid(
            "the magnetic lane's stored entropy is not finite",
        ));
    }
    Ok((energy, entropy))
}

/// A magnetic law's field state at `kelvin`, for a flux density of magnitude
/// `b` \[T\] in `volume_m3`, or why it has none.
fn law_state(
    law: &ThermalMagneticMaterial,
    kelvin: f64,
    b: f64,
    volume_m3: f64,
) -> Result<FieldThermalState, String> {
    law.at_temperature(kelvin)
        .and_then(|at: ResolvedThermalMagneticMaterial| at.state(b, volume_m3))
        .map_err(|error| error.to_string())
}

/// Adds `value` to the compensated sum `(sum, error)`.
fn compensated(total: &mut [f64; 2], value: f64) {
    let (sum, error) = two_sum(total[0], value);
    total[0] = sum;
    total[1] += error;
}

/// A step's account as its owner reports it.
pub(crate) struct Account {
    /// Stored energy \[J\] at the start of the step.
    pub(crate) stored_before_j: f64,
    /// Stored energy \[J\] at the end of the step.
    pub(crate) stored_after_j: f64,
    /// The owner's stored change \[J\].
    pub(crate) stored_change_j: f64,
    /// Work supplied through the source ports \[J\].
    pub(crate) supplied_work_j: f64,
    /// The owner's numerical (backward Euler) dissipation \[J\].
    pub(crate) numerical_dissipation_j: f64,
    /// The records of the step that carry heat: stranded and field cells
    /// and graph resistors.
    pub(crate) elements: usize,
    /// The lengths of the owner's circuit state and field potential.
    pub(crate) coordinates: [usize; 2],
}

/// A step's closed account, before its ports are written.
pub(crate) struct Closed {
    /// The heat \[J\] that leaves through each heat port: the receiver's
    /// heat, less `T dS` of its law elements in the magnetic lane.
    out: Vec<f64>,
    stored_change_j: f64,
    dissipated_j: f64,
    numerical_loss_j: f64,
    entropy_j_per_k: f64,
    rounding_j: f64,
}

/// The account of a step whose receivers took `heat` \[J\] at the ports'
/// temperatures `kelvin`: the stored change, with the magnetic lane's terms
/// and its number of laws (`lane`); the heat out; the numerical loss that
/// closes the balance, which must be the owner's numerical dissipation (with
/// the lane's terms) within the rounding floor and `energy_tolerance_j`; and
/// the entropy generated, `sum q / T` over the receivers' gross heat.
pub(crate) fn close(
    account: &Account,
    heat: &[f64],
    kelvin: &[f64],
    lane: Option<(&Reversible, usize)>,
    sources: usize,
    energy_tolerance_j: f64,
) -> Result<Closed, Refusal> {
    let heat_j: f64 = heat.iter().sum();
    let entropy: f64 = heat.iter().zip(kelvin).map(|(q, t)| q / t).sum();
    // In the magnetic lane the heat that leaves through each port is less
    // `T dS` of its law elements, and may be negative.
    let net: Option<Vec<f64>> = lane.map(|(terms, _)| {
        heat.iter()
            .zip(&terms.entropy_j_per_k)
            .zip(kelvin)
            .map(|((q, ds), t)| q - t * ds)
            .collect()
    });
    let out_j = net.as_ref().map_or(heat_j, |net| net.iter().sum());
    let stored_change_j = lane.map_or(account.stored_change_j, |(terms, _)| {
        account.stored_change_j + terms.stored_change_j
    });
    let numerical_loss_j = account.supplied_work_j - stored_change_j - out_j;
    let mut magnitude = account.stored_before_j.abs()
        + account.stored_after_j.abs()
        + account.supplied_work_j.abs()
        + account.numerical_dissipation_j.abs()
        + heat_j;
    let mut elements = account.elements;
    if let (Some((terms, laws)), Some(net)) = (lane, &net) {
        magnitude += terms.magnitude_j + net.iter().map(|q| q.abs()).sum::<f64>();
        elements += 2 * laws;
    }
    let [state, potential] = account.coordinates;
    let rounding = arithmetic_rounding_floor(
        64_usize
            .saturating_add(8_usize.saturating_mul(state))
            .saturating_add(8_usize.saturating_mul(potential))
            .saturating_add(8_usize.saturating_mul(elements + sources + heat.len())),
        magnitude,
    );
    let allowance = rounding + energy_tolerance_j;
    let expected_j = lane.map_or(account.numerical_dissipation_j, |(terms, _)| {
        account.numerical_dissipation_j + terms.numerical_loss_j
    });
    let defect_j = numerical_loss_j - expected_j;
    if ![rounding, numerical_loss_j, entropy, defect_j]
        .iter()
        .all(|value| value.is_finite())
    {
        return Err(Refusal::Invalid(
            "the winding's energy account is not finite",
        ));
    }
    if defect_j.abs() > allowance {
        return Err(Refusal::Balance {
            defect_j,
            tolerance_j: allowance,
        });
    }
    Ok(Closed {
        out: net.unwrap_or_else(|| heat.to_vec()),
        stored_change_j,
        dissipated_j: heat_j,
        numerical_loss_j,
        entropy_j_per_k: entropy,
        rounding_j: rounding,
    })
}

/// Writes a closed step's ports: the native step-average conjugates `y` of
/// the `sources` source ports, then minus each heat port's mean heat rate
/// over `dt_s`; and its receipt, with the owner's residual and the Newton
/// tolerance it met (the absolute one, or the relative one times the first
/// residual).
pub(crate) fn publish(
    closed: &Closed,
    report: &NewtonReport,
    newton: &NewtonKrylovConfig,
    y: &[f64],
    sources: usize,
    dt_s: f64,
    outputs: &mut [f64],
) -> Receipt {
    let initial = report
        .history
        .first()
        .map_or(report.residual_norm, |iteration| iteration.residual_before);
    let tolerance = newton
        .absolute_tolerance
        .max(newton.relative_tolerance * initial);
    outputs[..sources].copy_from_slice(&y[..sources]);
    for (rate, q) in outputs[sources..].iter_mut().zip(&closed.out) {
        *rate = -q / dt_s;
    }
    Receipt {
        stored_change_j: closed.stored_change_j,
        dissipated_j: closed.dissipated_j,
        numerical_loss_j: closed.numerical_loss_j,
        entropy_generated_j_per_k: closed.entropy_j_per_k,
        energy_tolerance_j: closed.rounding_j,
        residual: report.residual_norm,
        tolerance,
    }
}

/// Equal slices, bit for bit.
pub(crate) fn same_bits(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits())
}

/// Equal lists of vectors, bit for bit.
#[cfg(feature = "magnetic-hysteresis")]
pub(crate) fn same_vectors<const N: usize>(a: &[[f64; N]], b: &[[f64; N]]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| same_bits(x, y))
}

/// The fit's one query coordinate that is a static absolute temperature:
/// typed so, or the untyped `T` (as fs-em's silica optics reads it).
#[cfg(feature = "magnetic-hysteresis")]
pub(crate) fn fitted_temperature(material: &ResolvedMaterialStatePoint) -> Option<(&str, f64)> {
    let mut found = None;
    for (axis, value) in material.query_point() {
        let quantity = material.axis_quantities().get(axis).copied();
        let absolute = quantity
            .and_then(QuantitySpec::semantic_type)
            .is_some_and(|ty| {
                ty.kind() == QuantityKind::AbsoluteTemperature && ty.form() == ValueForm::Static
            })
            || axis == "T"
                && quantity.is_none_or(|q| q == QuantitySpec::dimensional(Temperature::DIMS));
        if absolute && found.replace((axis.as_str(), *value)).is_some() {
            return None;
        }
    }
    found
}

/// The graph's nodes, branches (ends and element) and transformers, each
/// list preceded by its length.
pub(crate) fn graph_words(graph: &CircuitGraph, words: &mut Vec<u64>) {
    words.extend([graph.node_count as u64, graph.branches.len() as u64]);
    for &(a, b, element) in &graph.branches {
        words.extend([a as u64, b as u64]);
        branch_words(element, words);
    }
    words.push(graph.transformers.len() as u64);
    for link in &graph.transformers {
        words.extend([
            link.primary.0 as u64,
            link.primary.1 as u64,
            link.secondary.0 as u64,
            link.secondary.1 as u64,
            link.ratio.to_bits(),
        ]);
    }
}

fn branch_words(element: Branch, words: &mut Vec<u64>) {
    match element {
        Branch::Resistor { ohms } => words.extend([0, ohms.to_bits()]),
        Branch::Inductor { henries } => words.extend([1, henries.to_bits()]),
        Branch::Capacitor { farads } => words.extend([2, farads.to_bits()]),
        Branch::VoltageSource { port } => words.extend([3, port as u64]),
        Branch::CurrentSource { port } => words.extend([4, port as u64]),
    }
}

/// A scalar conductor's resistivity bits and card identity.
pub(crate) fn conductor_words(conductor: &ResolvedConductorState, words: &mut Vec<u64>) {
    words.push(conductor.resistivity_ohm_m().to_bits());
    for chunk in conductor.material().card_identity().0.chunks(8) {
        let mut bytes = [0_u8; 8];
        bytes.copy_from_slice(chunk);
        words.push(u64::from_le_bytes(bytes));
    }
}

/// The graph resistors' laws: their count, then each branch and the law's
/// words.
pub(crate) fn resistor_law_words(
    laws: &BTreeMap<usize, Arc<dyn ResistorLaw>>,
    words: &mut Vec<u64>,
) {
    words.push(laws.len() as u64);
    for (&branch, law) in laws {
        words.push(branch as u64);
        law_words(words, |words| law.identity_words(words));
    }
}

/// The receivers: their count, then each one's mesh elements and graph
/// resistor branches, each list preceded by its length.
pub(crate) fn receiver_words(receivers: &[(&[usize], &[usize])], words: &mut Vec<u64>) {
    words.push(receivers.len() as u64);
    for &(cells, branches) in receivers {
        words.push(cells.len() as u64);
        words.extend(cells.iter().map(|&element| element as u64));
        words.push(branches.len() as u64);
        words.extend(branches.iter().map(|&branch| branch as u64));
    }
}

/// `value`'s `Debug` text (which prints each `f64` so that it parses back to
/// the same bits), preceded by its length in bytes.
pub(crate) fn debug_words(value: &impl fmt::Debug, words: &mut Vec<u64>) {
    let text = format!("{value:?}");
    words.push(text.len() as u64);
    byte_words(text.as_bytes(), words);
}

/// The magnetic lane: `b"magtemp1"` read little-endian, the number of laws,
/// then each law's mesh element and words (its card identity, temperature
/// axis, range and slope, and its source's permeability).
pub(crate) fn lane_words(laws: &BTreeMap<usize, ThermalMagneticMaterial>, words: &mut Vec<u64>) {
    words.extend([u64::from_le_bytes(*b"magtemp1"), laws.len() as u64]);
    for (&element, law) in laws {
        words.push(element as u64);
        law_words(words, |words| magnetic_words(law, words));
    }
}

/// Bytes as little-endian words, the last one padded with zeros.
fn byte_words(bytes: &[u8], words: &mut Vec<u64>) {
    for chunk in bytes.chunks(8) {
        let mut word = [0_u8; 8];
        word[..chunk.len()].copy_from_slice(chunk);
        words.push(u64::from_le_bytes(word));
    }
}

/// A magnetic law's card identity, temperature axis (preceded by its
/// length), range and slope, and its source's permeability.
fn magnetic_words(law: &ThermalMagneticMaterial, words: &mut Vec<u64>) {
    byte_words(&law.material().card_identity().0, words);
    words.push(law.temperature_axis().len() as u64);
    byte_words(law.temperature_axis().as_bytes(), words);
    let [low, high] = law.temperature_range_k();
    words.extend([
        low.to_bits(),
        high.to_bits(),
        law.permeability_slope_h_m_k().to_bits(),
        linear_magnetic_permeability(law.material()).map_or(u64::MAX, f64::to_bits),
    ]);
}

/// A law's words, preceded by their count.
pub(crate) fn law_words(words: &mut Vec<u64>, law: impl FnOnce(&mut Vec<u64>)) {
    let at = words.len();
    words.push(0);
    law(words);
    words[at] = (words.len() - at - 1) as u64;
}
