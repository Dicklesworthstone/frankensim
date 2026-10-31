//! Frequency-domain pressure/volume-flow graphs. Instruments are graph data;
//! the solver has no tone-hole, nasal, or organ-specific topology.
//!
//! Edge convention: [p_from,U_from] = M [p_to,U_to], with U directed from->to.
//! Both nodal currents point INTO the edge. The time convention is exp(-i wt).

use fs_la::eigen_complex::lu_complex;
use fs_material::gas::GasState;
use fs_math::c64::C64;
use std::collections::BTreeMap;

use crate::{Duct, DuctError, HoleState, LossModel, Segment, TabulatedLoad, Termination};

mod media;
pub use media::NetworkMedia;
mod loaded;
pub use loaded::{
    network_drive, network_drive_coherent, network_drive_coherent_regional, network_drive_coupled,
    network_drive_coupled_regional, network_drive_regional, network_drive_source_bank,
    network_drive_source_bank_regional, network_drive_sources, network_drive_sources_regional,
};

#[cfg(test)]
mod tests;

/// Dense graph node cap, checked before nodal allocation.
pub const MAX_NETWORK_NODES: usize = 128;
/// Edge cap, including parallel edges and cycles.
pub const MAX_NETWORK_EDGES: usize = 256;
/// Maximum labelled external ports.
pub const MAX_NETWORK_PORTS: usize = 32;
/// Aggregate segment cap over every edge.
pub const MAX_NETWORK_SEGMENTS: usize = 4096;
/// Cap after supplied junction multiports expand their incident-edge pressures.
pub const MAX_NETWORK_UNKNOWNS: usize = 256;
/// Pressure unknowns plus two inward-flow unknowns per retained scattering edge.
pub const MAX_NETWORK_WAVE_UNKNOWNS: usize = MAX_NETWORK_UNKNOWNS + 2 * MAX_NETWORK_EDGES;
/// Retained wave equations plus one current for every external port. Finite
/// source impedances require source currents as well as passive-load currents.
pub const MAX_NETWORK_LOADED_UNKNOWNS: usize = MAX_NETWORK_WAVE_UNKNOWNS + MAX_NETWORK_PORTS;
/// Maximum infinity-norm condition number for normalized edge and graph solves.
pub const MAX_NETWORK_CONDITION: f64 = 1.0e10;

/// Explicit treatment of the chain-local compact tone-hole corrections.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum JunctionCorrections {
    /// Preserve the existing chain's inner, series and mutual corrections.
    ChainLocal,
    /// Ideal junctions: retain the complete chimney line and mouth load, but
    /// omit inner/series/mutual compact-junction corrections.
    Disabled,
}

/// A passive boundary impedance in Pa s/m³. Radiation resolution is explicit.
#[derive(Debug, Clone)]
pub enum NetworkLoad {
    /// Zero outward flow.
    Closed,
    /// Zero acoustic pressure.
    PressureRelease,
    /// A finite nonnegative-real impedance at this solve's frequency.
    Impedance(C64),
    /// Radiation at a declared mouth radius; includes the exact LS owner.
    Radiation {
        /// Geometric mouth radius [m].
        radius: f64,
        /// Explicit radiation model and, for exact LS, quadrature resolution.
        termination: Termination,
    },
    /// A measured/baked load; its own support and passivity admission apply.
    Tabulated(TabulatedLoad),
    /// Fluid-filled rigid-frame absorber over a uniform-pressure boundary.
    /// The solve's actual fluid refreshes drag, inertia and thermal relaxation.
    /// Gas networks use their exterior gas; the shared-fluid chain uses its
    /// identical admitted fluid for both propagation and this boundary.
    /// This is a normal-incidence load, not a distributed duct-wall model.
    PorousRigid {
        /// Original porous skeleton and characterization-fluid viscosity.
        material: fs_material::porous::PorousFluidLaw,
        /// Porous layer thickness against its rigid backing [m].
        thickness_m: f64,
        /// Exposed area carrying the boundary's volume flow [m2].
        area_m2: f64,
    },
}

/// Nodal boundary condition. A bare junction has common pressure and zero net
/// injected volume flow. Compact junction physics is never added implicitly.
#[derive(Debug, Clone)]
pub enum NodeBoundary {
    /// Ideal massless junction.
    Junction,
    /// A supplied frequency-specific plane-wave multiport replacing the ideal
    /// common-pressure junction. No feature geometry is inferred here.
    JunctionMultiport(JunctionMultiport),
    /// A shunt load or pressure-release boundary.
    Termination(NetworkLoad),
    /// Labelled external volume-flow port. Other ports are not auto-terminated.
    Port {
        /// Unique external port label.
        label: String,
    },
}

/// End of an oriented edge, with current positive into that edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum EdgeEnd {
    /// Chain input plane.
    From,
    /// Chain output plane.
    To,
}

/// Stable binding of a feature port to an incident edge end.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct NetworkTerminal {
    /// Exact NetworkEdge::id.
    pub edge_id: String,
    /// Which edge end touches this junction.
    pub end: EdgeEnd,
}

/// Supplied junction admittance p -> U INTO the feature, in m³/(Pa s).
/// A 3-D/modal feature owner must perform its own valid plane-port projection;
/// this type neither invents evanescent closures nor certifies that projection.
#[derive(Debug, Clone)]
pub struct JunctionMultiport {
    /// Exact angular frequency of the supplied operator [rad/s].
    pub omega: f64,
    /// Terminal order of the row-major matrix, covering every incident end.
    pub terminals: Vec<NetworkTerminal>,
    /// Row-major plane-port admittance.
    pub y: Vec<C64>,
    /// Explicit feature/source identity.
    pub source_id: String,
}

fn admit_multiport(feature: &JunctionMultiport) -> Result<(), DuctError> {
    let n = feature.terminals.len();
    if n == 0
        || n > MAX_NETWORK_UNKNOWNS
        || feature.y.len() != n * n
        || feature.source_id.trim().is_empty()
        || !feature.omega.is_finite()
        || feature.omega <= 0.0
        || feature.y.iter().any(|&v| !finite(v))
    {
        return Err(bad(
            "network junction multiport needs finite frequency, source and square terminal matrix",
        ));
    }
    let scale = feature.y.iter().map(|z| z.abs()).fold(0.0, f64::max);
    let tol = 1.0e-12 * scale;
    let mut real = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            if (feature.y[i * n + j] - feature.y[j * n + i]).abs() > tol {
                return Err(bad("network junction multiport is not reciprocal (1e-12)"));
            }
            real[i * n + j] = if scale == 0.0 {
                0.0
            } else {
                0.5 * (feature.y[i * n + j].re / scale + feature.y[j * n + i].re / scale)
            };
        }
    }
    // Diagonal-pivoted LDL of the real Hermitian part. At a zero pivot the
    // remaining block must also vanish; positive diagonals alone are insufficient.
    for k in 0..n {
        let mut p = k;
        for i in k + 1..n {
            if real[i * n + i] > real[p * n + p] {
                p = i;
            }
        }
        if p != k {
            for j in 0..n {
                real.swap(k * n + j, p * n + j);
            }
            for i in 0..n {
                real.swap(i * n + k, i * n + p);
            }
        }
        let pivot = real[k * n + k];
        if !pivot.is_finite() || pivot < -1.0e-12 {
            return Err(bad("network junction multiport has active Hermitian part"));
        }
        if pivot <= 1.0e-12 {
            for i in k..n {
                for j in k..n {
                    if !real[i * n + j].is_finite() || real[i * n + j].abs() > 1.0e-12 {
                        return Err(bad(
                            "network junction multiport has indefinite zero-pivot block",
                        ));
                    }
                }
            }
            break;
        }
        for i in k + 1..n {
            for j in k + 1..n {
                let update = real[i * n + k] * real[k * n + j] / pivot;
                real[i * n + j] -= update;
            }
        }
    }
    Ok(())
}

/// A stable identity allows node-array reindexing without changing arithmetic.
#[derive(Debug, Clone)]
pub struct NetworkNode {
    /// Unique nonempty identity, independent of the array index.
    pub id: String,
    /// Junction, termination or external port.
    pub boundary: NodeBoundary,
}

/// An oriented edge carrying the existing duct-chain owner.
#[derive(Debug, Clone)]
pub struct NetworkEdge {
    /// Unique stable identity; fixes summation order even for parallel paths.
    pub id: String,
    /// Node index at the chain's input plane.
    pub from: usize,
    /// Node index at the chain's output plane.
    pub to: usize,
    /// Complete chain geometry, in SI.
    pub duct: Duct,
    /// Explicit compact correction policy for holes inside this edge.
    pub corrections: JunctionCorrections,
}

/// A connected graph; parallel paths and loops are admitted.
#[derive(Debug, Clone)]
pub struct DuctNetwork {
    /// Node indices are local references, not physical identities.
    pub nodes: Vec<NetworkNode>,
    /// Oriented duct chains.
    pub edges: Vec<NetworkEdge>,
}

/// A finite chain 2-port, normalized only for condition reporting.
#[derive(Debug, Clone, Copy)]
pub struct NetworkTwoPort {
    /// Physical ABCD matrix, row-major: [p_in,U_in] = M [p_out,U_out].
    pub abcd: [C64; 4],
    /// Infinity-norm condition of [A,B/Zref;C*Zref,D].
    pub condition_estimate: f64,
    /// Minimum shear number reported by the chain owner.
    pub min_shear_number: f64,
}

/// Solved open-circuit N-port impedance; all currents point into the network.
#[derive(Debug, Clone)]
pub struct NetworkResponse {
    gas: GasState,
    media: Option<NetworkMedia>,
    port_gases: Vec<Option<GasState>>,
    /// Frequency [rad/s].
    pub omega: f64,
    /// Port labels in the caller's node-array order.
    pub ports: Vec<String>,
    /// Row-major p=Z*U in Pa s/m³; unexcited ports have exactly zero injection.
    pub z: Vec<C64>,
    /// Infinity-norm condition of the grounded nodal system, or of the
    /// dimensionless row-equilibrated scattering constraints when retained.
    pub condition_estimate: f64,
    /// Maximum composition/admittance condition per edge, in stable id order.
    /// With retained wave constraints these are composition conditions only:
    /// a singular pressure-to-flow chart is not a singular physical edge.
    pub edge_conditions: Vec<(String, f64)>,
    /// Number of edges whose scattering equations were retained in the solve.
    /// Zero selects the original, smaller nodal-admittance solve.
    pub wave_constraint_edges: usize,
    /// Minimum shear number over all edges.
    pub min_shear_number: f64,
}

/// Port pressures and currents for unit flow injected at one active port.
#[derive(Debug, Clone)]
pub struct NetworkDriveResponse {
    /// Active port's pressure / injected volume flow [Pa s/m³].
    pub impedance: C64,
    /// Pressures / unit injected flow, in NetworkResponse::ports order.
    pub port_pressures: Vec<C64>,
    /// Currents into the graph / unit injected flow, in the same order.
    pub port_inflows: Vec<C64>,
}

/// One prescribed harmonic volume-flow source into a labelled external port.
/// All sources in a solve share its frequency and exp(-i omega t) clock.
#[derive(Debug, Clone)]
pub struct PortFlowSource {
    /// Exact external port label. Each driven port is specified once.
    pub label: String,
    /// Complex peak volume flow [m³/s], including the source's phase.
    /// Zero prescribes zero flow; it does not remove this boundary condition.
    pub volume_flow_m3_s: C64,
}

/// Harmonic source boundary on the solve's common exp(-i omega t) clock.
#[derive(Debug, Clone, Copy)]
pub enum PortExcitation {
    /// Ideal imposed complex peak flow into the network [m³/s].
    VolumeFlow(C64),
    /// Thevenin pressure generator: `p + Zs U = p_open`, where U points into
    /// the network. The finite passive source impedance includes the driver's
    /// own linear reaction, independently of the network and exterior loads.
    Pressure {
        /// Complex peak open-circuit pressure [Pa]. Zero leaves the source's
        /// impedance connected, so other sources can drive and load it.
        open_circuit_pressure_pa: C64,
        /// Source impedance [Pa s/m³], with finite real and imaginary parts
        /// and Re Zs >= 0. Zero is an ideal pressure source, not a refusal.
        impedance_pa_s_m3: C64,
    },
}

/// One explicit source, permitting pressure and prescribed-flow generators
/// in the same coherent solve. Every source port occurs exactly once.
#[derive(Debug, Clone)]
pub struct PortSource {
    /// Exact external port label.
    pub label: String,
    /// Physical boundary law and complex peak excitation.
    pub excitation: PortExcitation,
}

/// Reciprocal pressure generators with a complete passive internal impedance:
/// `p_i + sum_j Zs_ij U_j = p_open_i`, with every U directed into the graph.
/// The diagonal and mutual terms describe the source bank itself, separately
/// from the network and exterior load. Independent banks can be block diagonal.
#[derive(Debug, Clone)]
pub struct CoupledSourceBank {
    /// Frequency, exact port labels, source identity and full row-major source
    /// impedance [Pa s/m³]. The shared passive multiport carrier is reused;
    /// here its operator belongs to the generators, not an exterior load.
    /// Zs must be reciprocal and its Hermitian part positive semidefinite,
    /// within a relative 1e-12 numerical admission tolerance. No entry is clipped.
    pub impedance: CoupledTermination,
    /// Complex peak open-circuit pressure [Pa], in impedance-label order on
    /// the common exp(-i omega t) clock. Zero retains that generator's complete
    /// impedance and mutual interaction with every other bank port.
    pub open_circuit_pressures_pa: Vec<C64>,
}

struct ResolvedSources<'a> {
    conditions: Vec<Option<PortExcitation>>,
    bank: Option<&'a CoupledSourceBank>,
    bank_rows: Vec<Option<usize>>,
}

impl ResolvedSources<'_> {
    fn impedance(&self, port: usize, other: usize) -> Option<C64> {
        if let (Some(bank), Some(row), Some(column)) =
            (self.bank, self.bank_rows[port], self.bank_rows[other])
        {
            Some(bank.impedance.z[row * bank.impedance.labels.len() + column])
        } else if port == other {
            self.conditions[port]
                .and_then(PortExcitation::pressure)
                .map(|(_, impedance)| impedance)
        } else {
            None
        }
    }
}

impl PortExcitation {
    fn flow(self) -> Option<C64> {
        match self {
            Self::VolumeFlow(flow) => Some(flow),
            Self::Pressure { .. } => None,
        }
    }

    fn pressure(self) -> Option<(C64, C64)> {
        match self {
            Self::VolumeFlow(_) => None,
            Self::Pressure {
                open_circuit_pressure_pa,
                impedance_pa_s_m3,
            } => Some((open_circuit_pressure_pa, impedance_pa_s_m3)),
        }
    }
}

fn flow_sources(sources: &[PortFlowSource]) -> Result<Vec<PortSource>, DuctError> {
    if sources.is_empty() || sources.len() > MAX_NETWORK_PORTS {
        return Err(bad("network needs 1..=32 coherent sources"));
    }
    Ok(sources
        .iter()
        .map(|source| PortSource {
            label: source.label.clone(),
            excitation: PortExcitation::VolumeFlow(source.volume_flow_m3_s),
        })
        .collect())
}

/// Actual simultaneous harmonic fields, with interference retained before
/// calculating power. These are peak phasors, not unit-flow transfer ratios.
#[derive(Debug, Clone)]
pub struct NetworkCoherentResponse {
    /// Common angular frequency [rad/s].
    pub omega: f64,
    /// Port labels in the original graph's node-array order.
    pub ports: Vec<String>,
    /// Complex peak pressure [Pa] in port order.
    pub port_pressures: Vec<C64>,
    /// Complex peak volume flow into the graph [m³/s] in port order.
    pub port_inflows: Vec<C64>,
    /// Signed mean power into the graph, 0.5 Re(p conj(U)) [W], per port.
    /// A source may absorb power, and mutually coupled loads may exchange it.
    pub port_powers_w: Vec<f64>,
    /// Signed mean work supplied at all source/network interfaces [W]. This
    /// excludes dissipation inside any finite source impedance.
    pub source_power_w: f64,
    /// Signed mean work delivered to the external scalar/coupled loads [W].
    pub external_load_power_w: f64,
    /// Sum of all signed port powers [W], without clipping roundoff.
    /// This includes internal nodal loads, wall work and gas losses; it is not
    /// an attribution of heat to any particular physical mechanism.
    pub network_power_w: f64,
}

/// Coherent fields including finite-source back reaction and its work balance.
#[derive(Debug, Clone)]
pub struct NetworkSourceResponse {
    /// Actual network fields and signed source-interface / load / graph work.
    pub fields: NetworkCoherentResponse,
    /// Signed mean generator work per port [W]; zero at passive-load ports.
    /// Pressure sources use 0.5 Re(p_open conj(U)); ideal flow sources use
    /// their actual interface work. A generator may absorb power.
    pub generator_powers_w: Vec<f64>,
    /// Signed mean internal-source work per port [W]:
    /// 0.5 Re((Zs U)_i conj(U_i)). With scalar impedances this is the local
    /// nonnegative loss 0.5 Re(Zs) |U|². With mutual impedances an individual
    /// entry can be negative: it includes exchange through the source bank,
    /// including reactive coupling, and is not a local heat attribution.
    /// Zero for ideal flow sources and passive-load ports.
    pub source_dissipations_w: Vec<f64>,
    /// Signed total generator work [W]. It equals source-interface work plus
    /// source dissipation, up to the admitted linear solve's numerical error.
    pub generator_power_w: f64,
    /// Total internal source dissipation [W], nonnegative up to admission and
    /// solve roundoff for a passive bank. This is separate from exterior load
    /// power and the gas/wall/internal-load work inside the graph. Signed
    /// mutual terms are summed before this quantity is reported, never clipped.
    pub source_dissipation_w: f64,
}

fn source_conditions<'a>(
    labels: &[&str],
    sources: &[PortSource],
    bank: Option<&'a CoupledSourceBank>,
    omega: f64,
) -> Result<ResolvedSources<'a>, DuctError> {
    if labels.is_empty() || labels.len() > MAX_NETWORK_PORTS {
        return Err(bad("network needs 1..=32 external ports"));
    }
    let bank_count = bank.map_or(0, |bank| bank.impedance.labels.len());
    if sources.len() > MAX_NETWORK_PORTS
        || bank_count > MAX_NETWORK_PORTS
        || sources.len() + bank_count == 0
        || sources.len() + bank_count > MAX_NETWORK_PORTS
    {
        return Err(bad("network needs 1..=32 coherent sources"));
    }
    let mut conditions = vec![None; labels.len()];
    for source in sources {
        let port = labels
            .iter()
            .position(|&label| label == source.label)
            .ok_or_else(|| bad("unknown network source port label"))?;
        if conditions[port].is_some() {
            return Err(bad("network sources need unique ports"));
        }
        match source.excitation {
            PortExcitation::VolumeFlow(flow) if !finite(flow) => {
                return Err(bad("network source peak flow must be finite"));
            }
            PortExcitation::Pressure {
                open_circuit_pressure_pa,
                impedance_pa_s_m3,
            } => {
                if !finite(open_circuit_pressure_pa) {
                    return Err(bad("network source peak pressure must be finite"));
                }
                finite_impedance(impedance_pa_s_m3)?;
            }
            PortExcitation::VolumeFlow(_) => {}
        }
        conditions[port] = Some(source.excitation);
    }
    let mut bank_rows = vec![None; labels.len()];
    if let Some(bank) = bank {
        admit_coupled(&bank.impedance, omega)?;
        if bank.open_circuit_pressures_pa.len() != bank_count
            || bank.open_circuit_pressures_pa.iter().any(|&p| !finite(p))
        {
            return Err(bad("source bank needs one finite peak pressure per port"));
        }
        let scale = bank
            .impedance
            .z
            .iter()
            .map(|z| z.re.abs().max(z.im.abs()))
            .fold(0.0, f64::max);
        for row in 0..bank_count {
            for column in 0..row {
                let a = bank.impedance.z[row * bank_count + column];
                let b = bank.impedance.z[column * bank_count + row];
                if scale > 0.0
                    && (C64::new(a.re / scale, a.im / scale) - C64::new(b.re / scale, b.im / scale))
                        .abs()
                        > 1.0e-12
                {
                    return Err(bad("source bank impedance must be reciprocal"));
                }
            }
            let label = &bank.impedance.labels[row];
            let port = labels
                .iter()
                .position(|&candidate| candidate == label)
                .ok_or_else(|| bad("unknown source-bank port label"))?;
            if conditions[port].is_some() {
                return Err(bad("source-bank ports must be unique and disjoint"));
            }
            bank_rows[port] = Some(row);
            conditions[port] = Some(PortExcitation::Pressure {
                open_circuit_pressure_pa: bank.open_circuit_pressures_pa[row],
                impedance_pa_s_m3: bank.impedance.z[row * bank_count + row],
            });
        }
    }
    Ok(ResolvedSources {
        conditions,
        bank,
        bank_rows,
    })
}

// Neumaier accumulation preserves weak residuals of strongly cancelling
// sources. Callers supply canonical port order, independent of input order.
fn coherent_sum(values: impl Iterator<Item = C64>) -> C64 {
    let mut sum = C64::ZERO;
    let mut correction = C64::ZERO;
    for value in values {
        let next = sum + value;
        correction.re += if sum.re.abs() >= value.re.abs() {
            (sum.re - next.re) + value.re
        } else {
            (value.re - next.re) + sum.re
        };
        correction.im += if sum.im.abs() >= value.im.abs() {
            (sum.im - next.im) + value.im
        } else {
            (value.im - next.im) + sum.im
        };
        sum = next;
    }
    sum + correction
}

impl NetworkCoherentResponse {
    fn from_fields(
        omega: f64,
        labels: &[&str],
        sources: &[Option<PortExcitation>],
        port_pressures: Vec<C64>,
        port_inflows: Vec<C64>,
    ) -> Result<Self, DuctError> {
        let port_powers_w: Vec<f64> = port_pressures
            .iter()
            .zip(&port_inflows)
            .map(|(&p, &u)| 0.5 * (p * u.conj()).re)
            .collect();
        let mut canonical: Vec<usize> = (0..labels.len()).collect();
        canonical.sort_unstable_by(|&a, &b| labels[a].cmp(labels[b]));
        let work = |driven: Option<bool>| {
            coherent_sum(canonical.iter().filter_map(|&port| {
                driven
                    .is_none_or(|value| sources[port].is_some() == value)
                    .then(|| C64::from_re(port_powers_w[port]))
            }))
            .re
        };
        let source_power_w = work(Some(true));
        let external_load_power_w = -work(Some(false));
        let network_power_w = work(None);
        if port_pressures
            .iter()
            .chain(&port_inflows)
            .any(|&v| !finite(v))
            || port_powers_w.iter().any(|p| !p.is_finite())
            || !source_power_w.is_finite()
            || !external_load_power_w.is_finite()
            || !network_power_w.is_finite()
        {
            return Err(DuctError::Singular);
        }
        Ok(Self {
            omega,
            ports: labels.iter().map(|label| (*label).to_owned()).collect(),
            port_pressures,
            port_inflows,
            port_powers_w,
            source_power_w,
            external_load_power_w,
            network_power_w,
        })
    }
}

impl NetworkSourceResponse {
    fn from_fields(
        omega: f64,
        labels: &[&str],
        sources: &ResolvedSources<'_>,
        pressures: Vec<C64>,
        inflows: Vec<C64>,
    ) -> Result<Self, DuctError> {
        let fields = NetworkCoherentResponse::from_fields(
            omega,
            labels,
            &sources.conditions,
            pressures,
            inflows,
        )?;
        let mut generator_powers_w = vec![0.0; labels.len()];
        let mut source_dissipations_w = vec![0.0; labels.len()];
        let mut canonical: Vec<usize> = (0..labels.len()).collect();
        canonical.sort_unstable_by(|&a, &b| labels[a].cmp(labels[b]));
        for (port, source) in sources.conditions.iter().enumerate() {
            if let Some(source) = source {
                if let Some((pressure, impedance)) = source.pressure() {
                    let flow = fields.port_inflows[port];
                    generator_powers_w[port] = 0.5 * (pressure * flow.conj()).re;
                    source_dissipations_w[port] = if sources.bank_rows[port].is_some() {
                        let internal_pressure =
                            coherent_sum(canonical.iter().filter_map(|&other| {
                                sources
                                    .impedance(port, other)
                                    .map(|z| z * fields.port_inflows[other])
                            }));
                        0.5 * (internal_pressure * flow.conj()).re
                    } else {
                        (0.5 * impedance.re * flow.re) * flow.re
                            + (0.5 * impedance.re * flow.im) * flow.im
                    };
                } else {
                    generator_powers_w[port] = fields.port_powers_w[port];
                }
            }
        }
        let total = |powers: &[f64]| {
            coherent_sum(canonical.iter().map(|&port| C64::from_re(powers[port]))).re
        };
        let generator_power_w = total(&generator_powers_w);
        // Contract the Hermitian part for the total. Large reactive exchange
        // can otherwise bury a small real loss when signed port work cancels.
        // This is the same complete operator; no mutual term is removed from
        // the solve or the reported signed work at an individual source port.
        let source_dissipation_w = coherent_sum(canonical.iter().map(|&port| {
            if sources.bank_rows[port].is_none() {
                return C64::from_re(source_dissipations_w[port]);
            }
            let dissipative_pressure = coherent_sum(canonical.iter().filter_map(|&other| {
                let forward = sources.impedance(port, other)?;
                let reverse = sources
                    .impedance(other, port)
                    .expect("same complete source bank");
                let hermitian = forward.scale(0.5) + reverse.conj().scale(0.5);
                Some(hermitian * fields.port_inflows[other])
            }));
            C64::from_re(0.5 * (dissipative_pressure * fields.port_inflows[port].conj()).re)
        }))
        .re;
        if generator_powers_w
            .iter()
            .chain(&source_dissipations_w)
            .any(|power| !power.is_finite())
            || !generator_power_w.is_finite()
            || !source_dissipation_w.is_finite()
        {
            return Err(DuctError::Singular);
        }
        Ok(Self {
            fields,
            generator_powers_w,
            source_dissipations_w,
            generator_power_w,
            source_dissipation_w,
        })
    }
}

/// Explicit termination of one external port; no omitted-port default exists.
#[derive(Debug, Clone)]
pub struct PortTermination {
    /// Label in NetworkResponse::ports.
    pub label: String,
    /// Boundary physics, resolved at the graph frequency in its bound exterior.
    /// On the uniform API this is the graph's original gas.
    pub load: NetworkLoad,
}

/// One passive boundary shared by several labelled ports, in Pa s/m³:
/// `p = -Z U` with `U` into the network (so `-U` flows into the load). Mutual
/// entries carry, for example, the external radiation interaction of open tone
/// holes and a bell, taken from one multi-aperture radiation matrix.
#[derive(Debug, Clone)]
pub struct CoupledTermination {
    /// Exact angular frequency of the supplied matrix [rad/s].
    pub omega: f64,
    /// Port labels, in matrix order.
    pub labels: Vec<String>,
    /// Row-major acoustic impedance [Pa s/m³].
    pub z: Vec<C64>,
    /// Explicit source identity, for example a bake declaration hash.
    pub source_id: String,
}

/// Finite, nonempty, square, at the response frequency, with a positive
/// semidefinite Hermitian part (an active coupled load refuses). Exact
/// reciprocity is not required: a collocation radiation matrix is reciprocal
/// only to its discretization error.
pub(crate) fn admit_coupled(coupled: &CoupledTermination, omega: f64) -> Result<(), DuctError> {
    let n = coupled.labels.len();
    if n == 0
        || n > MAX_NETWORK_PORTS
        || coupled.z.len() != n * n
        || coupled.source_id.trim().is_empty()
        || coupled.omega.to_bits() != omega.to_bits()
        || coupled.z.iter().any(|&v| !finite(v))
    {
        return Err(bad(
            "network coupled termination needs the response frequency, a source and a finite square matrix",
        ));
    }
    // Normalize before forming H: finite near-maximum entries must neither
    // overflow into NaN and bypass admission nor underflow its passivity gate.
    let scale = coupled
        .z
        .iter()
        .map(|z| z.re.abs().max(z.im.abs()))
        .fold(0.0, f64::max);
    if scale == 0.0 {
        return Ok(());
    }
    let mut embedding = vec![0.0; 4 * n * n];
    for a in 0..n {
        for b in 0..n {
            let forward = coupled.z[a * n + b];
            let reverse = coupled.z[b * n + a];
            let h = C64::new(
                0.5 * (forward.re / scale + reverse.re / scale),
                0.5 * (forward.im / scale - reverse.im / scale),
            );
            embedding[a * 2 * n + b] = h.re;
            embedding[(a + n) * 2 * n + b + n] = h.re;
            embedding[a * 2 * n + b + n] = -h.im;
            embedding[(a + n) * 2 * n + b] = h.im;
        }
    }
    let (values, _) = fs_la::eigen::jacobi_eigh(&embedding, 2 * n);
    if values.iter().any(|value| !value.is_finite()) || values[0] < -1.0e-12 {
        return Err(bad(
            "network coupled termination has an active Hermitian part",
        ));
    }
    Ok(())
}

fn bad(what: &'static str) -> DuctError {
    DuctError::BadParameter { what }
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

fn validate_segment(segment: &Segment) -> Result<(), DuctError> {
    segment.validate()?;
    if let Segment::ToneHole {
        state: HoleState::Vent(sigma),
        ..
    } = *segment
    {
        if !sigma.is_finite() || !(0.0..=1.0).contains(&sigma) {
            return Err(bad("network hole fraction must be finite and in [0,1]"));
        }
    }
    Ok(())
}

fn finite_impedance(z: C64) -> Result<C64, DuctError> {
    if !finite(z) || z.re < 0.0 {
        return Err(bad("network load must be finite and passive (Re Z >= 0)"));
    }
    Ok(z)
}

fn load_impedance(
    load: &NetworkLoad,
    interior: Option<&GasState>,
    exterior: &GasState,
    omega: f64,
) -> Result<Option<C64>, DuctError> {
    match load {
        NetworkLoad::Closed => Ok(None),
        NetworkLoad::PressureRelease => Ok(Some(C64::ZERO)),
        NetworkLoad::Impedance(z) => Ok(Some(finite_impedance(*z)?)),
        NetworkLoad::Radiation {
            radius,
            termination,
        } => {
            if !radius.is_finite() || *radius <= 0.0 {
                return Err(bad("network mouth radius must be positive and finite"));
            }
            if matches!(termination, Termination::Closed | Termination::IdealOpen) {
                return Err(bad("network radiation load needs a radiating termination"));
            }
            let gas = match interior {
                Some(gas) => gas,
                None if matches!(
                    termination,
                    Termination::FlangedOpen | Termination::UnflangedThinWall { .. }
                ) =>
                {
                    return Err(bad(
                        "network radiation mouth needs one unambiguous incident gas at its reference plane",
                    ));
                }
                None => exterior,
            };
            if matches!(termination, Termination::UnflangedThinWall { .. }) && gas != exterior {
                return Err(bad(
                    "exact thin-wall radiation requires the same gas inside and outside the mouth",
                ));
            }
            let (z, _) = crate::exterior::termination_impedance(
                *termination,
                gas,
                exterior,
                *radius,
                omega,
            )?;
            Ok(Some(finite_impedance(z.ok_or_else(|| {
                bad("network radiation load has no impedance")
            })?)?))
        }
        NetworkLoad::Tabulated(table) => Ok(Some(finite_impedance(table.z_at(omega)?)?)),
        NetworkLoad::PorousRigid {
            material,
            thickness_m,
            area_m2,
        } => {
            if !area_m2.is_finite() || *area_m2 <= 0.0 {
                return Err(bad("porous boundary area must be positive and finite"));
            }
            let surface = material
                .surface_impedance_rigid(exterior, omega, *thickness_m)
                .map_err(|_| bad("porous boundary gas, geometry or material support refused"))?;
            Ok(Some(porous_volume_impedance(surface, *area_m2)?))
        }
    }
}

/// The material owner uses exp(+i wt); all duct consumers use exp(-i wt).
/// U = A v, so the pressure/volume-flow load is conj(Z_surface)/A.
/// Gas graphs and generic-fluid chains share this one conversion/admission.
pub(crate) fn porous_volume_impedance(surface: C64, area_m2: f64) -> Result<C64, DuctError> {
    if !area_m2.is_finite() || area_m2 <= 0.0 {
        return Err(bad("porous boundary area must be positive and finite"));
    }
    let load = C64::new(surface.re / area_m2, -surface.im / area_m2);
    if (surface.re != 0.0 && load.re == 0.0) || (surface.im != 0.0 && load.im == 0.0) {
        return Err(bad("porous boundary impedance conversion underflows"));
    }
    finite_impedance(load)
}

fn norm_inf(a: &[C64], n: usize) -> f64 {
    a.chunks_exact(n)
        .map(|row| row.iter().map(|z| z.abs()).sum::<f64>())
        .fold(0.0, f64::max)
}

// Full inverse columns are already needed for the N-port solve and allow an
// actual infinity-norm condition computation, rather than a diagonal heuristic.
fn inverse(a: &[C64], n: usize) -> Result<(Vec<C64>, f64), DuctError> {
    if n == 0 || a.iter().any(|&z| !finite(z)) {
        return Err(DuctError::Singular);
    }
    let norm = norm_inf(a, n);
    if !norm.is_finite() || norm == 0.0 {
        return Err(DuctError::Singular);
    }
    let scaled: Vec<C64> = a.iter().map(|z| z.scale(1.0 / norm)).collect();
    let lu = lu_complex(&scaled, n).map_err(|_| DuctError::Singular)?;
    let mut inv = vec![C64::ZERO; n * n];
    for col in 0..n {
        let mut b = vec![C64::ZERO; n];
        b[col] = C64::ONE;
        lu.solve(&mut b);
        for row in 0..n {
            inv[row * n + col] = b[row];
        }
    }
    let condition = norm_inf(&scaled, n) * norm_inf(&inv, n);
    if !condition.is_finite() || condition > MAX_NETWORK_CONDITION {
        return Err(DuctError::Singular);
    }
    for v in &mut inv {
        *v = v.scale(1.0 / norm);
    }
    if inv.iter().any(|&z| !finite(z)) {
        return Err(DuctError::Singular);
    }
    Ok((inv, condition))
}

// Consume the canonical local-wave composition. Opposite signed load probes
// cannot recover an opaque chain: their tiny through response is unobservable.
fn chain_owned(
    duct: &Duct,
    gas: &GasState,
    omega: f64,
    loss: LossModel,
    zref: f64,
) -> Result<([C64; 4], f64), DuctError> {
    let _reference = zref;
    let options = crate::compose::ComposeOptions::default();
    let chain = crate::compose::duct_chain(duct, gas, omega, loss, None, options, true, || false)?;
    let (s, _) = crate::compose::compose_scattering(&chain.blocks, options, &mut || false)?;
    Ok((s.abcd()?, chain.min_shear_number))
}

/// Stable pressure-to-inward-flow stamp, including opaque/strongly attenuating
/// edges for which a finite well-conditioned ABCD export does not exist.
pub fn network_chain_admittance(
    duct: &Duct,
    gas: &GasState,
    omega: f64,
    loss: LossModel,
    corrections: JunctionCorrections,
) -> Result<([C64; 4], f64, f64), DuctError> {
    network_chain_admittance_in_exterior(duct, gas, gas, omega, loss, corrections)
}

/// The original stable chain stamp with independent exterior gas at embedded
/// open chimneys. Propagation and closed chimneys retain the interior gas.
/// The two media must have the same static pressure; neither is averaged.
pub fn network_chain_admittance_in_exterior(
    duct: &Duct,
    interior: &GasState,
    exterior: &GasState,
    omega: f64,
    loss: LossModel,
    corrections: JunctionCorrections,
) -> Result<([C64; 4], f64, f64), DuctError> {
    let (s, condition, shear) =
        network_chain_scattering(duct, interior, exterior, omega, loss, corrections)?;
    let (y, c) = s.admittance(MAX_NETWORK_CONDITION)?;
    let (_, yc) = inverse(&y, 2)?;
    Ok((y, condition.max(c).max(yc), shear))
}

fn network_chain_scattering(
    duct: &Duct,
    interior: &GasState,
    exterior: &GasState,
    omega: f64,
    loss: LossModel,
    corrections: JunctionCorrections,
) -> Result<(crate::compose::TwoPort, f64, f64), DuctError> {
    let gas = interior;
    validate_medium(gas, omega)?;
    crate::exterior::admit(interior, exterior)?;
    if duct.segments.len() > MAX_NETWORK_SEGMENTS {
        return Err(bad("network segment cap exceeded"));
    }
    let options = crate::compose::ComposeOptions {
        maximum_condition: MAX_NETWORK_CONDITION,
        ..crate::compose::ComposeOptions::default()
    };
    let (chain, _) = crate::compose::duct_chain_in_exterior(
        duct,
        gas,
        exterior,
        omega,
        loss,
        crate::WallModel::Rigid,
        options,
        corrections == JunctionCorrections::ChainLocal,
        || false,
    )?;
    let (s, condition) = crate::compose::compose_scattering(&chain.blocks, options, &mut || false)?;
    Ok((s, condition, chain.min_shear_number))
}

struct WaveEdge {
    ends: [Option<usize>; 2],
    two_port: crate::compose::TwoPort,
    condition: f64,
}

// A half-wave edge has singular I+S although its loaded graph can be regular.
// Keep b=S*a instead of eliminating the currents through (I+S)^-1:
// (I-S)p - (I+S)R*U = 0. KCL closes the graph, including pressure-release ends.
// Each current unknown is R_edge*U; multiplying KCL by one graph reference
// puts every unknown/equation in pressure units before row equilibration.
fn wave_equations(
    boundary_y: &[C64],
    pressures: usize,
    edges: &[WaveEdge],
    external_currents: usize,
) -> Result<(Vec<C64>, usize, f64), DuctError> {
    let n = pressures + 2 * edges.len() + external_currents;
    if external_currents > MAX_NETWORK_PORTS || n > MAX_NETWORK_LOADED_UNKNOWNS {
        return Err(bad(
            "network loaded scattering-constraint cap exceeded (800)",
        ));
    }
    let reference = edges
        .first()
        .ok_or(DuctError::Singular)?
        .two_port
        .reference_impedance;
    let mut a = vec![C64::ZERO; n * n];
    for i in 0..pressures {
        for j in 0..pressures {
            a[i * n + j] = boundary_y[i * pressures + j].scale(reference);
        }
    }
    for (e, edge) in edges.iter().enumerate() {
        let current = pressures + 2 * e;
        let flow_scale = reference / edge.two_port.reference_impedance;
        for row in 0..2 {
            if let Some(node) = edge.ends[row] {
                a[node * n + current + row] =
                    a[node * n + current + row] + C64::from_re(flow_scale);
            }
            for col in 0..2 {
                let identity = if row == col { C64::ONE } else { C64::ZERO };
                let s = edge.two_port.scattering[2 * row + col];
                if let Some(node) = edge.ends[col] {
                    a[(current + row) * n + node] = a[(current + row) * n + node] + identity - s;
                }
                a[(current + row) * n + current + col] = -(identity + s);
            }
        }
    }
    Ok((a, n, reference))
}

fn wave_pressure_inverse(
    boundary_y: &[C64],
    pressures: usize,
    edges: &[WaveEdge],
) -> Result<(Vec<C64>, f64), DuctError> {
    let (mut a, n, reference) = wave_equations(boundary_y, pressures, edges, 0)?;
    let mut source_scale = vec![0.0; pressures];
    for (row, values) in a.chunks_exact_mut(n).enumerate() {
        let scale = values.iter().map(|z| z.abs()).fold(0.0, f64::max);
        if !scale.is_finite() || scale == 0.0 || values.iter().any(|&z| !finite(z)) {
            return Err(DuctError::Singular);
        }
        for value in values {
            *value = C64::new(value.re / scale, value.im / scale);
        }
        if row < pressures {
            source_scale[row] = reference / scale;
        }
    }
    let (inv, condition) = inverse(&a, n)?;
    let mut pressure_inverse = vec![C64::ZERO; pressures * pressures];
    for i in 0..pressures {
        for j in 0..pressures {
            let value = inv[i * n + j].scale(source_scale[j]);
            if !finite(value) {
                return Err(DuctError::Singular);
            }
            pressure_inverse[i * pressures + j] = value;
        }
    }
    Ok((pressure_inverse, condition))
}

/// Compute the complete chain's physical ABCD block and reject an unqualified
/// normalized condition number. Lower-level propagation refusals are preserved.
pub fn network_chain_two_port(
    duct: &Duct,
    gas: &GasState,
    omega: f64,
    loss: LossModel,
    corrections: JunctionCorrections,
) -> Result<NetworkTwoPort, DuctError> {
    validate_medium(gas, omega)?;
    let first = duct.segments.first().ok_or(DuctError::EmptyDuct)?;
    if duct.segments.len() > MAX_NETWORK_SEGMENTS {
        return Err(bad("network segment cap exceeded"));
    }
    for s in &duct.segments {
        validate_segment(s)?;
    }
    let zref = gas.density * gas.sound_speed
        / (core::f64::consts::PI * first.inlet_radius() * first.inlet_radius());
    let (m, shear) = if corrections == JunctionCorrections::ChainLocal
        || !duct
            .segments
            .iter()
            .any(|s| matches!(s, Segment::ToneHole { .. }))
    {
        chain_owned(duct, gas, omega, loss, zref)?
    } else {
        let mut m = [C64::ONE, C64::ZERO, C64::ZERO, C64::ONE];
        let mut shear = f64::INFINITY;
        for s in &duct.segments {
            let (block, rv) = if let Segment::ToneHole {
                hole_radius,
                chimney_height,
                state,
                ..
            } = *s
            {
                let chimney = Duct {
                    segments: vec![Segment::Cylinder {
                        radius: hole_radius,
                        length: chimney_height,
                    }],
                };
                let chimney_loss = if loss == LossModel::WideTube {
                    LossModel::Bessel
                } else {
                    loss
                };
                let zi = |termination| {
                    crate::input_impedance(&chimney, gas, omega, chimney_loss, termination)
                };
                let admittance = match state {
                    HoleState::Open => zi(Termination::FlangedOpen)?.impedance.recip(),
                    HoleState::Closed => zi(Termination::Closed)?.impedance.recip(),
                    HoleState::Vent(sigma) if sigma.is_finite() && (0.0..=1.0).contains(&sigma) => {
                        zi(Termination::FlangedOpen)?.impedance.recip().scale(sigma)
                            + zi(Termination::Closed)?
                                .impedance
                                .recip()
                                .scale(1.0 - sigma)
                    }
                    HoleState::Vent(_) => {
                        return Err(bad("network hole fraction must be finite and in [0,1]"));
                    }
                };
                let rv = crate::segment_wave(gas, hole_radius, omega, chimney_loss)?.shear_number;
                ([C64::ONE, C64::ZERO, admittance, C64::ONE], rv)
            } else {
                chain_owned(
                    &Duct {
                        segments: vec![s.clone()],
                    },
                    gas,
                    omega,
                    loss,
                    zref,
                )?
            };
            m = crate::mul2(m, block);
            shear = shear.min(rv);
        }
        (m, shear)
    };
    let normalized = [m[0], m[1].scale(1.0 / zref), m[2].scale(zref), m[3]];
    let (_, condition) = inverse(&normalized, 2)?;
    let det = m[0] * m[3] - m[1] * m[2];
    if (det - C64::ONE).abs() > 1.0e-10 {
        return Err(DuctError::Singular);
    }
    Ok(NetworkTwoPort {
        abcd: m,
        condition_estimate: condition,
        min_shear_number: shear,
    })
}

fn validate_medium(gas: &GasState, omega: f64) -> Result<(), DuctError> {
    crate::validate_gas(gas)?;
    if !omega.is_finite()
        || omega <= 0.0
        || !gas.density.is_finite()
        || gas.density <= 0.0
        || !gas.sound_speed.is_finite()
        || gas.sound_speed <= 0.0
    {
        return Err(bad(
            "network needs positive finite angular frequency, gas density and sound speed",
        ));
    }
    Ok(())
}

impl DuctNetwork {
    /// Validate topology and all structural caps before dense allocation.
    pub fn validate(&self) -> Result<(), DuctError> {
        let n = self.nodes.len();
        if n == 0 || n > MAX_NETWORK_NODES {
            return Err(bad("network node count must be in 1..=128"));
        }
        if self.edges.len() > MAX_NETWORK_EDGES {
            return Err(bad("network edge cap exceeded (256)"));
        }
        let ports: Vec<&str> = self
            .nodes
            .iter()
            .filter_map(|node| match &node.boundary {
                NodeBoundary::Port { label } => Some(label.as_str()),
                _ => None,
            })
            .collect();
        if ports.is_empty() || ports.len() > MAX_NETWORK_PORTS {
            return Err(bad("network port count must be in 1..=32"));
        }
        let mut ids: Vec<&str> = self.nodes.iter().map(|node| node.id.as_str()).collect();
        ids.sort_unstable();
        let mut labels = ports;
        labels.sort_unstable();
        let mut edge_ids: Vec<&str> = self.edges.iter().map(|e| e.id.as_str()).collect();
        edge_ids.sort_unstable();
        for names in [&ids, &labels, &edge_ids] {
            if names.iter().any(|s| s.trim().is_empty() || s.len() > 256)
                || names.windows(2).any(|s| s[0] == s[1])
            {
                return Err(bad(
                    "network node/edge identities and port labels must be unique, nonempty and <=256 bytes",
                ));
            }
        }
        let mut segments = 0usize;
        for e in &self.edges {
            if e.from >= n || e.to >= n {
                return Err(bad("network edge node index out of range"));
            }
            if e.duct.segments.is_empty() {
                return Err(DuctError::EmptyDuct);
            }
            segments = segments.saturating_add(e.duct.segments.len());
            if segments > MAX_NETWORK_SEGMENTS {
                return Err(bad("network aggregate segment cap exceeded (4096)"));
            }
            let mut length = 0.0;
            for s in &e.duct.segments {
                validate_segment(s)?;
                length += s.axial_length();
            }
            if !length.is_finite() || length <= 0.0 {
                return Err(bad("network edge needs positive finite axial length"));
            }
        }
        let mut seen = vec![false; n];
        seen[0] = true;
        loop {
            let mut changed = false;
            for e in &self.edges {
                if seen[e.from] != seen[e.to] {
                    seen[e.from] = true;
                    seen[e.to] = true;
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
        if seen.iter().any(|&v| !v) {
            return Err(bad("network graph is disconnected"));
        }
        let mut unknowns = 0usize;
        for (node_index, node) in self.nodes.iter().enumerate() {
            if let NodeBoundary::JunctionMultiport(feature) = &node.boundary {
                unknowns = unknowns.saturating_add(feature.terminals.len());
                if unknowns > MAX_NETWORK_UNKNOWNS {
                    return Err(bad("network expanded-pressure cap exceeded (256)"));
                }
                admit_multiport(feature)?;
                let mut expected = Vec::new();
                for e in &self.edges {
                    if e.from == node_index {
                        expected.push(NetworkTerminal {
                            edge_id: e.id.clone(),
                            end: EdgeEnd::From,
                        });
                    }
                    if e.to == node_index {
                        expected.push(NetworkTerminal {
                            edge_id: e.id.clone(),
                            end: EdgeEnd::To,
                        });
                    }
                }
                let mut actual = feature.terminals.clone();
                actual.sort();
                expected.sort();
                if actual != expected {
                    return Err(bad(
                        "junction multiport terminals must bind each incident edge end exactly once",
                    ));
                }
            } else {
                unknowns = unknowns.saturating_add(1);
            }
            if unknowns > MAX_NETWORK_UNKNOWNS {
                return Err(bad("network expanded-pressure cap exceeded (256)"));
            }
        }
        Ok(())
    }
}

/// Solve a connected frequency-domain graph. Internal nodes are eliminated by
/// dense LU; labelled ports retain their full mutual impedance matrix.
pub fn network_port_impedance(
    net: &DuctNetwork,
    gas: &GasState,
    omega: f64,
    loss: LossModel,
) -> Result<NetworkResponse, DuctError> {
    network_port_impedance_bound(net, gas, None, omega, loss)
}

/// Solve the same graph with an exact edge-id gas binding and an independent
/// homogeneous exterior. All branch propagation, viscothermal loss, embedded
/// chimneys and boundary loads consume their actual media in the original
/// nodal solve. Junction pressure/volume flow remain continuous.
///
/// Coupling between external mouths on different edges still needs an explicit
/// [`CoupledTermination`]; a gas binding does not invent a radiation matrix.
pub fn network_port_impedance_regional(
    net: &DuctNetwork,
    media: &NetworkMedia,
    omega: f64,
    loss: LossModel,
) -> Result<NetworkResponse, DuctError> {
    media.validate_for(net)?;
    network_port_impedance_bound(net, media.exterior_gas(), Some(media), omega, loss)
}

fn node_gas(
    net: &DuctNetwork,
    exterior: &GasState,
    media: Option<&NetworkMedia>,
    node: usize,
) -> Option<GasState> {
    let Some(media) = media else {
        return Some(*exterior);
    };
    let mut gas = None;
    for edge in &net.edges {
        if edge.from == node || edge.to == node {
            let current = *media
                .edge_gas(&edge.id)
                .expect("complete edge identities admitted");
            if gas.is_some_and(|previous| previous != current) {
                return None;
            }
            gas = Some(current);
        }
    }
    gas
}

fn network_port_impedance_bound(
    net: &DuctNetwork,
    exterior: &GasState,
    media: Option<&NetworkMedia>,
    omega: f64,
    loss: LossModel,
) -> Result<NetworkResponse, DuctError> {
    net.validate()?;
    validate_medium(exterior, omega)?;
    // A pressure node may join unlike gases, but that does not define the
    // single interior modal basis required by a flanged mouth's Schur load.
    // Retain absence explicitly; finite/coupled loads remain admissible there.
    let local_gas = |node| node_gas(net, exterior, media, node);
    let mut order: Vec<usize> = (0..net.nodes.len()).collect();
    order.sort_unstable_by(|&a, &b| net.nodes[a].id.cmp(&net.nodes[b].id));
    let mut ground = vec![false; net.nodes.len()];
    let mut loads = vec![None; net.nodes.len()];
    for &node in &order {
        if let NodeBoundary::Termination(load) = &net.nodes[node].boundary {
            loads[node] = load_impedance(load, local_gas(node).as_ref(), exterior, omega)?;
            ground[node] = loads[node].is_some_and(|z| z == C64::ZERO);
        }
    }
    let mut slots = vec![None; net.nodes.len()];
    let mut terminal_slots = BTreeMap::new();
    let mut n = 0usize;
    for &node in &order {
        if let NodeBoundary::JunctionMultiport(feature) = &net.nodes[node].boundary {
            if feature.omega.to_bits() != omega.to_bits() {
                return Err(bad(
                    "junction multiport frequency differs from network query",
                ));
            }
            let mut terminals = feature.terminals.clone();
            terminals.sort();
            for terminal in terminals {
                terminal_slots.insert(terminal, n);
                n += 1;
            }
        } else if !ground[node] {
            slots[node] = Some(n);
            n += 1;
        }
    }
    let mut y = vec![C64::ZERO; n * n];
    for &node in &order {
        if let (Some(i), Some(z)) = (slots[node], loads[node]) {
            y[i * n + i] = z.recip();
        }
        if let NodeBoundary::JunctionMultiport(feature) = &net.nodes[node].boundary {
            let d = feature.terminals.len();
            for (i, ti) in feature.terminals.iter().enumerate() {
                for (j, tj) in feature.terminals.iter().enumerate() {
                    let a = terminal_slots[ti];
                    let b = terminal_slots[tj];
                    y[a * n + b] = y[a * n + b] + feature.y[i * d + j];
                }
            }
        }
    }
    let mut edges: Vec<&NetworkEdge> = net.edges.iter().collect();
    edges.sort_unstable_by(|a, b| a.id.cmp(&b.id));
    let boundary_y = y.clone();
    let mut wave_edges = Vec::with_capacity(edges.len());
    let mut keep_waves = false;
    let mut edge_conditions = Vec::with_capacity(edges.len());
    let mut shear = f64::INFINITY;
    for e in edges {
        let interior = media.map_or(exterior, |bound| {
            bound
                .edge_gas(&e.id)
                .expect("complete edge identities admitted")
        });
        let (two_port, composition_condition, edge_shear) =
            network_chain_scattering(&e.duct, interior, exterior, omega, loss, e.corrections)?;
        let endpoint = |node: usize, end: EdgeEnd| {
            if matches!(net.nodes[node].boundary, NodeBoundary::JunctionMultiport(_)) {
                terminal_slots
                    .get(&NetworkTerminal {
                        edge_id: e.id.clone(),
                        end,
                    })
                    .copied()
            } else {
                slots[node]
            }
        };
        let ends = [endpoint(e.from, EdgeEnd::From), endpoint(e.to, EdgeEnd::To)];
        wave_edges.push(WaveEdge {
            ends,
            two_port,
            condition: composition_condition,
        });
        let admittance = two_port
            .admittance(MAX_NETWORK_CONDITION)
            .and_then(|(stamp, c)| {
                let (_, yc) = inverse(&stamp, 2)?;
                Ok((stamp, composition_condition.max(c).max(yc)))
            });
        let condition = match admittance {
            Ok((stamp, condition)) => {
                for row in 0..2 {
                    for col in 0..2 {
                        if let (Some(i), Some(j)) = (ends[row], ends[col]) {
                            y[i * n + j] = y[i * n + j] + stamp[2 * row + col];
                        }
                    }
                }
                condition
            }
            Err(DuctError::Singular) => {
                keep_waves = true;
                composition_condition
            }
            Err(error) => return Err(error),
        };
        edge_conditions.push((e.id.clone(), condition));
        shear = shear.min(edge_shear);
    }
    let nodal = if keep_waves {
        Err(DuctError::Singular)
    } else {
        inverse(&y, n)
    };
    let (inv, condition) = match nodal {
        Ok(result) => result,
        Err(DuctError::Singular) => {
            keep_waves = true;
            let result = wave_pressure_inverse(&boundary_y, n, &wave_edges)?;
            for ((_, condition), edge) in edge_conditions.iter_mut().zip(&wave_edges) {
                *condition = edge.condition;
            }
            result
        }
        Err(error) => return Err(error),
    };
    let port_nodes: Vec<usize> = net
        .nodes
        .iter()
        .enumerate()
        .filter_map(|(i, node)| matches!(node.boundary, NodeBoundary::Port { .. }).then_some(i))
        .collect();
    let mut ports = Vec::new();
    let mut z = Vec::new();
    for &node in &port_nodes {
        if let NodeBoundary::Port { label } = &net.nodes[node].boundary {
            ports.push(label.clone());
        }
        let i = slots[node].ok_or(DuctError::Singular)?;
        for &other in &port_nodes {
            let j = slots[other].ok_or(DuctError::Singular)?;
            z.push(inv[i * n + j]);
        }
    }
    Ok(NetworkResponse {
        gas: *exterior,
        media: media.cloned(),
        port_gases: port_nodes.iter().map(|&node| local_gas(node)).collect(),
        omega,
        ports,
        z,
        condition_estimate: condition,
        edge_conditions,
        wave_constraint_edges: if keep_waves { wave_edges.len() } else { 0 },
        min_shear_number: shear,
    })
}

impl NetworkResponse {
    /// Complete original regional records, retained with the solved operator.
    /// The uniform API returns `None` and retains its original numerical gas.
    #[must_use]
    pub const fn media(&self) -> Option<&NetworkMedia> {
        self.media.as_ref()
    }

    /// Terminate the response in its already-bound exterior. Gas-dependent
    /// radiation and porous loads cannot acquire a different atmosphere after
    /// the interior and embedded chimneys have been assembled.
    pub fn drive_bound(
        &self,
        active: &str,
        loads: &[PortTermination],
    ) -> Result<NetworkDriveResponse, DuctError> {
        self.drive_with(active, loads, None, &self.gas)
    }

    /// [`Self::drive_bound`] with a supplied full exterior impedance matrix.
    /// The matrix producer remains responsible for its source medium and
    /// geometry; the original passivity/frequency/port checks remain mandatory.
    pub fn drive_coupled_bound(
        &self,
        active: &str,
        loads: &[PortTermination],
        coupled: &CoupledTermination,
    ) -> Result<NetworkDriveResponse, DuctError> {
        self.drive_with(active, loads, Some(coupled), &self.gas)
    }

    /// Drive all prescribed sources simultaneously, retaining their relative
    /// phases. Every remaining port needs exactly one scalar or coupled load.
    /// The already-bound exterior resolves all gas-dependent load physics.
    /// One loaded factorization serves the complete source vector; powers are
    /// computed from the combined fields, never by summing isolated powers.
    pub fn drive_coherent(
        &self,
        sources: &[PortFlowSource],
        loads: &[PortTermination],
        coupled: Option<&CoupledTermination>,
    ) -> Result<NetworkCoherentResponse, DuctError> {
        Ok(self
            .drive_sources(&flow_sources(sources)?, loads, coupled)?
            .fields)
    }

    /// Simultaneously close ideal flow and finite-impedance pressure sources.
    /// Pressure-source flows are solved with every mutual network and exterior
    /// interaction; open-circuit pressures are not imposed at the graph ports.
    /// The sources share this operator's frequency and peak-phasor clock.
    pub fn drive_sources(
        &self,
        sources: &[PortSource],
        loads: &[PortTermination],
        coupled: Option<&CoupledTermination>,
    ) -> Result<NetworkSourceResponse, DuctError> {
        self.drive_sources_and_bank(sources, None, loads, coupled)
    }

    /// Close a complete mutually coupled generator bank, optionally alongside
    /// independent scalar pressure or prescribed-flow sources. The bank's
    /// full operator enters the same solve as the graph and exterior matrix;
    /// every port must belong to exactly one source or passive load.
    pub fn drive_source_bank(
        &self,
        sources: &[PortSource],
        bank: &CoupledSourceBank,
        loads: &[PortTermination],
        coupled: Option<&CoupledTermination>,
    ) -> Result<NetworkSourceResponse, DuctError> {
        self.drive_sources_and_bank(sources, Some(bank), loads, coupled)
    }

    fn drive_sources_and_bank(
        &self,
        sources: &[PortSource],
        bank: Option<&CoupledSourceBank>,
        loads: &[PortTermination],
        coupled: Option<&CoupledTermination>,
    ) -> Result<NetworkSourceResponse, DuctError> {
        let labels: Vec<&str> = self.ports.iter().map(String::as_str).collect();
        let sources = source_conditions(&labels, sources, bank, self.omega)?;
        let (pressures, inflows) = self.drive_sources_with(&sources, loads, coupled, &self.gas)?;
        NetworkSourceResponse::from_fields(self.omega, &labels, &sources, pressures, inflows)
    }

    /// Terminate EVERY other labelled port explicitly, then return the driving
    /// point p/U at `active`. Missing/duplicate/unknown terminations refuse.
    pub fn driving_point(
        &self,
        active: &str,
        loads: &[PortTermination],
        gas: &GasState,
    ) -> Result<C64, DuctError> {
        Ok(self.drive(active, loads, gas)?.impedance)
    }

    /// Apply explicit loads and return all port pressures and volume flows.
    pub fn drive(
        &self,
        active: &str,
        loads: &[PortTermination],
        gas: &GasState,
    ) -> Result<NetworkDriveResponse, DuctError> {
        self.drive_with(active, loads, None, gas)
    }

    /// [`Self::drive`] with one coupled boundary: the ports in
    /// `coupled.labels` are terminated together through its full matrix, and
    /// every other port through `loads`. Each non-active port is terminated
    /// exactly once.
    pub fn drive_coupled(
        &self,
        active: &str,
        loads: &[PortTermination],
        coupled: &CoupledTermination,
        gas: &GasState,
    ) -> Result<NetworkDriveResponse, DuctError> {
        self.drive_with(active, loads, Some(coupled), gas)
    }

    fn drive_with(
        &self,
        active: &str,
        loads: &[PortTermination],
        coupled: Option<&CoupledTermination>,
        gas: &GasState,
    ) -> Result<NetworkDriveResponse, DuctError> {
        let sources = [PortSource {
            label: active.to_owned(),
            excitation: PortExcitation::VolumeFlow(C64::ONE),
        }];
        let labels: Vec<&str> = self.ports.iter().map(String::as_str).collect();
        let sources = source_conditions(&labels, &sources, None, self.omega)?;
        let (port_pressures, port_inflows) =
            self.drive_sources_with(&sources, loads, coupled, gas)?;
        let active = self
            .ports
            .iter()
            .position(|label| label == active)
            .expect("source port was admitted");
        Ok(NetworkDriveResponse {
            impedance: port_pressures[active],
            port_pressures,
            port_inflows,
        })
    }

    fn drive_sources_with(
        &self,
        sources: &ResolvedSources<'_>,
        loads: &[PortTermination],
        coupled: Option<&CoupledTermination>,
        gas: &GasState,
    ) -> Result<(Vec<C64>, Vec<C64>), DuctError> {
        validate_medium(gas, self.omega)?;
        if *gas != self.gas {
            return Err(bad("network response and termination gas states differ"));
        }
        let n = self.ports.len();
        if self.z.len() != n * n
            || self.port_gases.len() != n
            || n == 0
            || n > MAX_NETWORK_PORTS
            || self.z.iter().any(|&z| !finite(z))
        {
            return Err(bad("invalid network response matrix"));
        }
        let conditions = &sources.conditions;
        let source_flows: Vec<Option<C64>> = conditions
            .iter()
            .map(|condition| condition.and_then(PortExcitation::flow))
            .collect();
        let coupled_count = coupled.map_or(0, |c| c.labels.len());
        let source_count = conditions.iter().filter(|source| source.is_some()).count();
        if loads.len() + coupled_count != n - source_count {
            return Err(bad(
                "every other network port needs an explicit termination",
            ));
        }
        let mut resolved: Vec<Option<C64>> = conditions
            .iter()
            .map(|condition| condition.and_then(PortExcitation::pressure).map(|(_, z)| z))
            .collect();
        let mut seen: Vec<bool> = conditions.iter().map(Option::is_some).collect();
        // Port index of each coupled-matrix row, and each port's coupled row.
        let mut coupled_rows = vec![None; n];
        if let Some(c) = coupled {
            admit_coupled(c, self.omega)?;
            for (row, label) in c.labels.iter().enumerate() {
                let index = self
                    .ports
                    .iter()
                    .position(|s| s == label)
                    .ok_or(bad("unknown network coupled-termination label"))?;
                if seen[index] {
                    return Err(bad("duplicate or active-port coupled termination"));
                }
                seen[index] = true;
                coupled_rows[index] = Some(row);
            }
        }
        for load in loads {
            let index = self
                .ports
                .iter()
                .position(|s| s == &load.label)
                .ok_or(bad("unknown network termination label"))?;
            if seen[index] {
                return Err(bad("duplicate or active-port termination"));
            }
            seen[index] = true;
            resolved[index] =
                load_impedance(&load.load, self.port_gases[index].as_ref(), gas, self.omega)?;
        }
        if seen.iter().any(|&v| !v) {
            return Err(bad("network port termination omitted"));
        }
        // A rigid external termination imposes U=0, so it needs no unknown.
        // Coupled ports are always unknowns.
        let mut live: Vec<usize> = (0..n)
            .filter(|&j| {
                source_flows[j].is_none() && (resolved[j].is_some() || coupled_rows[j].is_some())
            })
            .collect();
        live.sort_unstable_by(|&a, &b| self.ports[a].cmp(&self.ports[b]));
        let mut currents: Vec<C64> = source_flows
            .iter()
            .map(|value| value.unwrap_or(C64::ZERO))
            .collect();
        let mut canonical: Vec<usize> = (0..n).collect();
        canonical.sort_unstable_by(|&a, &b| self.ports[a].cmp(&self.ports[b]));
        // Keep the original unit-source arithmetic unchanged. Other source
        // vectors retain their phase, with compensation for deep cancellation.
        let unit_source = (source_count == 1)
            .then(|| source_flows.iter().position(|&flow| flow == Some(C64::ONE)))
            .flatten();
        let incident: Vec<C64> = (0..n)
            .map(|row| {
                unit_source.map_or_else(
                    || {
                        coherent_sum(canonical.iter().filter_map(|&port| {
                            source_flows[port].map(|u| self.z[row * n + port] * u)
                        }))
                    },
                    |active| self.z[row * n + active],
                )
            })
            .collect();
        let m = live.len();
        let mut a = vec![C64::ZERO; m * m];
        for (i, &p) in live.iter().enumerate() {
            for (j, &q) in live.iter().enumerate() {
                a[i * m + j] = self.z[p * n + q];
                if let Some(z) = sources.impedance(p, q) {
                    a[i * m + j] = a[i * m + j] + z;
                } else if let (Some(c), Some(r), Some(s)) =
                    (coupled, coupled_rows[p], coupled_rows[q])
                {
                    a[i * m + j] = a[i * m + j] + c.z[r * c.labels.len() + s];
                } else if i == j {
                    a[i * m + j] = a[i * m + j] + resolved[p].expect("finite load");
                }
            }
        }
        if m > 0 {
            let (inv, _) = inverse(&a, m)?;
            for (i, &p) in live.iter().enumerate() {
                if unit_source.is_some() {
                    for (j, &q) in live.iter().enumerate() {
                        currents[p] = currents[p] - inv[i * m + j] * incident[q];
                    }
                } else {
                    currents[p] = coherent_sum(live.iter().enumerate().map(|(j, &q)| {
                        let pressure = conditions[q]
                            .and_then(PortExcitation::pressure)
                            .map_or(C64::ZERO, |(pressure, _)| pressure);
                        inv[i * m + j] * (pressure - incident[q])
                    }));
                }
            }
        }
        let mut pressures = vec![C64::ZERO; n];
        for i in 0..n {
            if unit_source.is_some() {
                for &j in &canonical {
                    pressures[i] = pressures[i] + self.z[i * n + j] * currents[j];
                }
            } else {
                pressures[i] =
                    coherent_sum(canonical.iter().map(|&j| self.z[i * n + j] * currents[j]));
            }
        }
        if pressures.iter().chain(&currents).any(|&v| !finite(v)) {
            return Err(DuctError::Singular);
        }
        Ok((pressures, currents))
    }
}
