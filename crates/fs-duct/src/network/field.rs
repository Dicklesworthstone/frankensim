//! Actual internal pressure/flow fields of a coherently driven duct network.
//!
//! One retained loaded solve provides every source column, including at an
//! open-circuit pole. Interior fields use the original local-wave blocks from
//! both ends, never an exponentially growing reverse transfer matrix. The
//! directed volume flow is positive from the edge's `from` node to its `to`
//! node. Boundary work and finite-window flux are not local dissipated heat.

use super::*;
use crate::compose::{ComposeOptions, TwoPort};
use fs_fft::harmonic_window::{HarmonicProduct, HarmonicWindow, ProperHarmonicProduct};
use fs_math::dd::Dd;

mod welch;
pub use welch::NetworkFieldBinStatistics;

mod spectrum;
pub use spectrum::{
    MAX_FIELD_FREQUENCIES, MAX_FIELD_SPECTRAL_LINES, NetworkFieldExposureLimits,
    NetworkFieldFactorStatistics, NetworkFieldFrequency, NetworkFieldSpectralFactor,
    NetworkFieldSpectralLine, NetworkFieldSpectrumBudget, NetworkFieldSpectrumRequest,
    PreparedNetworkFieldExposure, PreparedNetworkFieldSpectrum, WindowedNetworkFieldExposure,
};

/// Maximum locations in one prepared observation bank.
pub const MAX_FIELD_PROBES: usize = 1024;
/// Maximum independent factors in one finite-window observation.
pub const MAX_FIELD_FACTORS: usize = 256;

/// A section in the ORIGINAL authored segment, without remeshing the network.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum NetworkFieldSection {
    /// Immediately inside this segment's inlet. At a radius step, its area
    /// determines particle velocity; pressure and volume flow remain continuous.
    Inlet { segment: usize },
    /// Immediately inside this segment's outlet. A compact tone hole has
    /// distinct inlet/outlet states at the same axial coordinate.
    Outlet { segment: usize },
    /// A point within a cylindrical segment, including either endpoint.
    /// Cones/profiles retain their native endpoint observation until their
    /// own spatial propagator exports a resolved interior basis.
    Cylinder { segment: usize, offset_m: f64 },
}

/// Stable physical edge identity and unambiguous section.
#[derive(Clone, Debug, PartialEq)]
pub struct NetworkFieldProbe {
    pub edge_id: String,
    pub section: NetworkFieldSection,
}

/// Resolved location and the original edge medium at that location.
#[derive(Clone, Debug, PartialEq)]
pub struct NetworkFieldLocation {
    pub probe: NetworkFieldProbe,
    pub axial_position_m: f64,
    pub area_m2: f64,
    pub density_kg_m3: f64,
    pub sound_speed_m_s: f64,
}

/// Additional spatial work bounds. The original graph caps and 1e10 loaded
/// condition bound still apply, including its at-most-800-unknown dense solve.
#[derive(Clone, Copy, Debug)]
pub struct NetworkFieldBudget {
    pub maximum_probes: usize,
    /// Total prepared pressure/flow entries plus temporary four-component
    /// edge-end entries, per source column. Original port-bank storage is
    /// separately bounded by the existing 32-port/source cap.
    pub maximum_transfer_entries: usize,
    /// Sum of native propagation blocks over selected edges; prefix and suffix
    /// workspace is at most twice this many TwoPorts.
    pub maximum_chain_blocks: usize,
}

impl Default for NetworkFieldBudget {
    fn default() -> Self {
        Self {
            maximum_probes: MAX_FIELD_PROBES,
            maximum_transfer_entries: 98_304,
            maximum_chain_blocks: crate::compose::MAX_COMPOSITION_BLOCKS,
        }
    }
}

/// All source boundaries and loads remain fixed during prepared application.
pub struct NetworkFieldRequest<'a> {
    pub omega: f64,
    pub loss: LossModel,
    pub source_clock: &'a str,
    pub sources: &'a [PortSource],
    pub source_bank: Option<&'a CoupledSourceBank>,
    pub loads: &'a [PortTermination],
    pub coupled_load: Option<&'a CoupledTermination>,
    pub probes: &'a [NetworkFieldProbe],
    pub budget: NetworkFieldBudget,
}

/// Actual coherent peak phasors at one location, using exp(-i omega t).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NetworkFieldSample {
    pub pressure_pa: C64,
    pub volume_flow_m3_s: C64,
    /// Axial cross-section-mean particle velocity U/A; not a resolved radial
    /// viscothermal velocity profile.
    pub particle_velocity_m_s: C64,
    /// Cycle mean signed power along the edge, 1/2 Re(p conj(U)).
    pub mean_power_w: f64,
}

/// Uncentered second moments of physical real fields over one finite window.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NetworkFieldMoments {
    pub mean_pressure_pa: f64,
    pub mean_volume_flow_m3_s: f64,
    pub mean_particle_velocity_m_s: f64,
    pub mean_square_pressure_pa2: f64,
    pub mean_square_volume_flow_m6_s2: f64,
    pub mean_square_particle_velocity_m2_s2: f64,
    pub mean_power_w: f64,
    pub mean_intensity_w_m2: f64,
}

/// Immutable, original-owner source-to-internal-field transfer bank.
/// Repeated source/factor observations perform contractions only.
#[derive(Clone, Debug)]
pub struct PreparedNetworkFields {
    ports: PreparedNetworkSources,
    source_clock: String,
    locations: Vec<NetworkFieldLocation>,
    pressure: Vec<C64>,
    flow: Vec<C64>,
}

impl PreparedNetworkFields {
    /// Uniform gas, with an explicit clock, point/work budget and cancellation.
    /// The existing synchronous dense graph factorization is bounded but does
    /// not poll internally; this API polls before/after it and between native
    /// propagation blocks, probes and source columns.
    pub fn new(
        net: &DuctNetwork,
        gas: &GasState,
        request: &NetworkFieldRequest<'_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DuctError> {
        Self::prepare(net, gas, None, request, cancelled)
    }

    /// Preserve the original per-edge records and independent exterior gas.
    pub fn new_regional(
        net: &DuctNetwork,
        media: &NetworkMedia,
        request: &NetworkFieldRequest<'_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DuctError> {
        media.validate_for(net)?;
        Self::prepare(net, media.exterior_gas(), Some(media), request, cancelled)
    }

    fn prepare(
        net: &DuctNetwork,
        exterior: &GasState,
        media: Option<&NetworkMedia>,
        request: &NetworkFieldRequest<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DuctError> {
        poll(&mut cancelled, 0)?;
        net.validate()?;
        validate_medium(exterior, request.omega)?;
        let budget = request.budget;
        if request.source_clock.trim().is_empty()
            || request.probes.is_empty()
            || request.probes.len() > MAX_FIELD_PROBES
            || request.probes.len() > budget.maximum_probes
            || budget.maximum_probes > MAX_FIELD_PROBES
            || budget.maximum_chain_blocks == 0
            || budget.maximum_chain_blocks > crate::compose::MAX_COMPOSITION_BLOCKS
        {
            return Err(bad(
                "spatial observation needs a clock and bounded nonempty probes",
            ));
        }
        let count = request
            .sources
            .len()
            .checked_add(
                request
                    .source_bank
                    .map_or(0, |bank| bank.impedance.labels.len()),
            )
            .ok_or_else(|| bad("spatial source count overflows"))?;
        if count == 0 || count > MAX_NETWORK_PORTS {
            return Err(bad("spatial observation needs 1..=32 sources"));
        }
        let mut edges: Vec<&NetworkEdge> = net.edges.iter().collect();
        edges.sort_unstable_by(|a, b| a.id.cmp(&b.id));
        let mut selected = Vec::new();
        let mut locations = Vec::with_capacity(request.probes.len());
        let mut probe_edges = Vec::with_capacity(request.probes.len());
        for (index, probe) in request.probes.iter().enumerate() {
            poll(&mut cancelled, index)?;
            let edge_index = edges
                .binary_search_by(|edge| edge.id.cmp(&probe.edge_id))
                .map_err(|_| bad("spatial probe names an unknown edge"))?;
            let edge = edges[edge_index];
            let gas = media.map_or(exterior, |binding| {
                binding
                    .edge_gas(&edge.id)
                    .expect("regional edge identities admitted")
            });
            locations.push(resolve_location(probe, edge, gas)?);
            probe_edges.push(edge_index);
            selected.push(edge_index);
        }
        selected.sort_unstable();
        selected.dedup();
        let entries = (2 * request.probes.len() + 4 * selected.len())
            .checked_mul(count)
            .ok_or_else(|| bad("spatial transfer count overflows"))?;
        if entries > budget.maximum_transfer_entries {
            return Err(bad("spatial field transfer budget exceeded"));
        }
        let (ports, endpoints) = super::loaded::prepare_spatial_transfers(
            net,
            exterior,
            media,
            request.omega,
            request.loss,
            request.sources,
            request.source_bank,
            request.loads,
            request.coupled_load,
            &selected,
        )?;
        poll(&mut cancelled, request.probes.len())?;
        let mut pressure = vec![C64::ZERO; request.probes.len() * count];
        let mut flow = pressure.clone();
        let mut used_blocks = 0usize;
        for (selected_index, &edge_index) in selected.iter().enumerate() {
            let edge = edges[edge_index];
            let gas = media.map_or(exterior, |binding| {
                binding
                    .edge_gas(&edge.id)
                    .expect("regional edge identities admitted")
            });
            if used_blocks == budget.maximum_chain_blocks {
                return Err(bad("spatial propagation block budget exceeded"));
            }
            let options = ComposeOptions {
                maximum_condition: MAX_NETWORK_CONDITION,
                maximum_blocks: budget.maximum_chain_blocks - used_blocks,
                ..ComposeOptions::default()
            };
            let mut boundaries = Vec::with_capacity(edge.duct.segments.len() + 1);
            let (chain, _) = crate::compose::duct_chain_with_segment_boundaries(
                &edge.duct,
                gas,
                exterior,
                request.omega,
                request.loss,
                crate::WallModel::Rigid,
                options,
                edge.corrections == JunctionCorrections::ChainLocal,
                Some(&mut boundaries),
                &mut cancelled,
            )?;
            used_blocks += chain.blocks.len();
            let reference = chain.reference_impedance;
            let mut prefix = Vec::with_capacity(chain.blocks.len() + 1);
            prefix.push(TwoPort::through(reference)?);
            for (index, &block) in chain.blocks.iter().enumerate() {
                poll(&mut cancelled, index)?;
                prefix.push(prefix[index].star(block, MAX_NETWORK_CONDITION)?.0);
            }
            let mut suffix = vec![TwoPort::through(reference)?; chain.blocks.len() + 1];
            for index in (0..chain.blocks.len()).rev() {
                poll(&mut cancelled, index)?;
                suffix[index] = chain.blocks[index]
                    .star(suffix[index + 1], MAX_NETWORK_CONDITION)?
                    .0;
            }
            for (probe_index, location) in locations.iter().enumerate() {
                if probe_edges[probe_index] != edge_index {
                    continue;
                }
                poll(&mut cancelled, probe_index)?;
                let (left, right) = match location.probe.section {
                    NetworkFieldSection::Inlet { segment } => {
                        let cut = boundaries[segment];
                        (prefix[cut], suffix[cut])
                    }
                    NetworkFieldSection::Outlet { segment } => {
                        let cut = boundaries[segment + 1];
                        (prefix[cut], suffix[cut])
                    }
                    NetworkFieldSection::Cylinder { segment, offset_m } => {
                        let Segment::Cylinder { radius, length } = edge.duct.segments[segment]
                        else {
                            unreachable!("cylinder probe admitted")
                        };
                        let before = boundaries[segment];
                        let after = boundaries[segment + 1];
                        if offset_m == 0.0 {
                            (prefix[before], suffix[before])
                        } else if offset_m == length {
                            (prefix[after], suffix[after])
                        } else {
                            let wave =
                                crate::segment_wave(gas, radius, request.omega, request.loss)?;
                            let a = TwoPort::from_wave(
                                &Segment::Cylinder {
                                    radius,
                                    length: offset_m,
                                },
                                &wave,
                                reference,
                            )?;
                            let b = TwoPort::from_wave(
                                &Segment::Cylinder {
                                    radius,
                                    length: length - offset_m,
                                },
                                &wave,
                                reference,
                            )?;
                            (
                                prefix[before].star(a, MAX_NETWORK_CONDITION)?.0,
                                b.star(suffix[after], MAX_NETWORK_CONDITION)?.0,
                            )
                        }
                    }
                };
                let weights = section_weights(left, right)?;
                for column in 0..count {
                    poll(&mut cancelled, column)?;
                    let [p0, q0, p1, q1] = endpoints[selected_index * count + column];
                    // Pressure-wave coordinates. The right endpoint's inward
                    // flow is -q1. No characteristic-impedance approximation.
                    let incoming = [
                        (p0.scale(0.5) + q0.scale(0.5 * reference)),
                        (p1.scale(0.5) - q1.scale(0.5 * reference)),
                    ];
                    let p = weighted(weights[0], incoming);
                    let q = weighted(weights[1], incoming);
                    if !finite(p) || !finite(q) {
                        return Err(DuctError::Singular);
                    }
                    pressure[probe_index * count + column] = p;
                    flow[probe_index * count + column] = q;
                }
            }
        }
        poll(&mut cancelled, used_blocks)?;
        Ok(Self {
            ports,
            source_clock: request.source_clock.to_owned(),
            locations,
            pressure,
            flow,
        })
    }

    /// Original exterior source/work bank from the SAME loaded factorization.
    pub fn port_response(&self) -> &PreparedNetworkSources {
        &self.ports
    }
    pub fn locations(&self) -> &[NetworkFieldLocation] {
        &self.locations
    }
    pub fn source_clock(&self) -> &str {
        &self.source_clock
    }
    /// Row-major pressure per canonical source-coordinate SI unit.
    pub fn pressure_transfer(&self) -> &[C64] {
        &self.pressure
    }
    /// Row-major directed volume flow per source-coordinate SI unit.
    pub fn flow_transfer(&self) -> &[C64] {
        &self.flow
    }

    /// Combine a complete source vector before taking any quadratic quantity.
    /// Cancellation/nonfinite input returns no partially accepted output.
    pub fn sample(
        &self,
        excitation: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<NetworkFieldSample>, DuctError> {
        self.validate_excitation(excitation)?;
        let mut result = Vec::with_capacity(self.locations.len());
        for (index, location) in self.locations.iter().enumerate() {
            poll(&mut cancelled, index)?;
            let [p, q] = self.parts(index, excitation)?;
            let pressure_pa = rounded(p);
            let volume_flow_m3_s = rounded(q);
            let particle_velocity_m_s = C64::new(
                volume_flow_m3_s.re / location.area_m2,
                volume_flow_m3_s.im / location.area_m2,
            );
            let mean_power_w = ((p[0] * q[0] + p[1] * q[1]) * Dd::from_f64(0.5)).to_f64();
            if !finite(pressure_pa)
                || !finite(volume_flow_m3_s)
                || !finite(particle_velocity_m_s)
                || !mean_power_w.is_finite()
            {
                return Err(DuctError::Singular);
            }
            result.push(NetworkFieldSample {
                pressure_pa,
                volume_flow_m3_s,
                particle_velocity_m_s,
                mean_power_w,
            });
        }
        poll(&mut cancelled, result.len())?;
        Ok(result)
    }

    /// Finite-window real-field moments at this prepared frequency. Mean and
    /// proper factors use the same source coordinate order and physical clock.
    /// Every factor has E[z]=E[z²]=0, E[|z|²]=1; distinct factors are uncorrelated.
    /// No source phase, bandwidth, Gaussian law or stationarity is inferred.
    /// This is a single-frequency source bank; cross-frequency coherence must
    /// be composed by a spectral observation owner before temporal averaging.
    pub fn observe_window(
        &self,
        window: HarmonicWindow,
        mean: Option<&[C64]>,
        proper_factors: &[&[C64]],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<NetworkFieldMoments>, DuctError> {
        poll(&mut cancelled, 0)?;
        if proper_factors.len() > MAX_FIELD_FACTORS {
            return Err(bad("spatial observation factor cap exceeded"));
        }
        if let Some(values) = mean {
            self.validate_excitation(values)?;
        }
        for values in proper_factors {
            self.validate_excitation(values)?;
        }
        let omega = self.ports.omega();
        let kernel = window.product_kernel(omega, omega).map_err(window_error)?;
        let dc = window.product_kernel(omega, 0.0).map_err(window_error)?;
        let proper = window
            .proper_product_kernel(omega, omega)
            .map_err(window_error)?;
        let one = [Dd::ONE, Dd::ZERO];
        let mut result = Vec::with_capacity(self.locations.len());
        for (index, location) in self.locations.iter().enumerate() {
            poll(&mut cancelled, index)?;
            let mut moments = [Dd::ZERO; 5];
            if let Some(values) = mean {
                let [p, q] = self.parts(index, values)?;
                moments[0] = dc
                    .mean_real_product_parts_dd(p, one)
                    .map_err(window_error)?;
                moments[1] = dc
                    .mean_real_product_parts_dd(q, one)
                    .map_err(window_error)?;
                accumulate(&mut moments, p, q, Kernel::Real(kernel))?;
            }
            for (factor, values) in proper_factors.iter().enumerate() {
                poll(&mut cancelled, factor)?;
                let [p, q] = self.parts(index, values)?;
                accumulate(&mut moments, p, q, Kernel::Proper(proper))?;
            }
            let [mp, mq, pp, qq, pq] = moments.map(Dd::to_f64);
            if pp < 0.0 || qq < 0.0 {
                return Err(bad("spatial self moment became negative"));
            }
            let area = location.area_m2;
            let output = NetworkFieldMoments {
                mean_pressure_pa: mp,
                mean_volume_flow_m3_s: mq,
                mean_particle_velocity_m_s: mq / area,
                mean_square_pressure_pa2: pp,
                mean_square_volume_flow_m6_s2: qq,
                mean_square_particle_velocity_m2_s2: (qq / area) / area,
                mean_power_w: pq,
                mean_intensity_w_m2: pq / area,
            };
            if [
                mp,
                mq,
                pp,
                qq,
                pq,
                output.mean_particle_velocity_m_s,
                output.mean_square_particle_velocity_m2_s2,
                output.mean_intensity_w_m2,
            ]
            .iter()
            .any(|v| !v.is_finite())
            {
                return Err(DuctError::Singular);
            }
            result.push(output);
        }
        poll(&mut cancelled, result.len())?;
        Ok(result)
    }

    fn validate_excitation(&self, excitation: &[C64]) -> Result<(), DuctError> {
        if excitation.len() != self.ports.source_coordinates().len()
            || excitation.iter().any(|&z| !finite(z))
        {
            return Err(bad(
                "spatial excitation must match the finite canonical source coordinates",
            ));
        }
        Ok(())
    }

    fn parts(&self, row: usize, excitation: &[C64]) -> Result<[[Dd; 2]; 2], DuctError> {
        let count = excitation.len();
        let p = dot_parts(&self.pressure[row * count..(row + 1) * count], excitation);
        let q = dot_parts(&self.flow[row * count..(row + 1) * count], excitation);
        if p.into_iter()
            .chain(q)
            .any(|v| !v.hi.is_finite() || !v.lo.is_finite())
        {
            return Err(DuctError::Singular);
        }
        Ok([p, q])
    }
}

fn resolve_location(
    probe: &NetworkFieldProbe,
    edge: &NetworkEdge,
    gas: &GasState,
) -> Result<NetworkFieldLocation, DuctError> {
    let (segment, offset, radius) = match probe.section {
        NetworkFieldSection::Inlet { segment } => {
            let s = edge
                .duct
                .segments
                .get(segment)
                .ok_or_else(|| bad("unknown probe segment"))?;
            (segment, 0.0, s.inlet_radius())
        }
        NetworkFieldSection::Outlet { segment } => {
            let s = edge
                .duct
                .segments
                .get(segment)
                .ok_or_else(|| bad("unknown probe segment"))?;
            (segment, s.axial_length(), s.outlet_radius())
        }
        NetworkFieldSection::Cylinder { segment, offset_m } => {
            let Some(&Segment::Cylinder { radius, length }) = edge.duct.segments.get(segment)
            else {
                return Err(bad("interior field probe requires an authored cylinder"));
            };
            if !offset_m.is_finite() || offset_m < 0.0 || offset_m > length {
                return Err(bad("cylindrical field probe lies outside its segment"));
            }
            (segment, offset_m, radius)
        }
    };
    let axial_position_m = edge.duct.segments[..segment]
        .iter()
        .fold(offset, |x, s| x + s.axial_length());
    let area_m2 = core::f64::consts::PI * radius * radius;
    if !axial_position_m.is_finite() || !area_m2.is_finite() || area_m2 <= 0.0 {
        return Err(bad(
            "field location has unrepresentable axial position or area",
        ));
    }
    Ok(NetworkFieldLocation {
        probe: probe.clone(),
        axial_position_m,
        area_m2,
        density_kg_m3: gas.density,
        sound_speed_m_s: gas.sound_speed,
    })
}

/// Incoming pressure waves at the cut obey
/// a = L21 a_left + L22 b, b = R11 a + R12 a_right.
/// Solving the two equations retains both sides of a standing wave. Pressure
/// is a+b and directed flow is (a-b)/R, even at either physical field node.
fn section_weights(left: TwoPort, right: TwoPort) -> Result<[[C64; 2]; 2], DuctError> {
    let [_, _, l21, l22] = left.scattering;
    let [r11, r12, _, _] = right.scattering;
    let reference = left.reference_impedance;
    if reference != right.reference_impedance {
        return Err(bad("spatial scattering references differ"));
    }
    let (denominator, condition) = crate::compose::feedback(l22, r11);
    if !condition.is_finite() || condition > MAX_NETWORK_CONDITION {
        return Err(DuctError::Singular);
    }
    let a = [l21 / denominator, l22 * (r12 / denominator)];
    let b = [r11 * a[0], r12 / denominator];
    let p = [a[0] + b[0], a[1] + b[1]];
    let q = [
        (a[0] - b[0]).scale(1.0 / reference),
        (a[1] - b[1]).scale(1.0 / reference),
    ];
    if p.into_iter().chain(q).any(|z| !finite(z)) {
        return Err(DuctError::Singular);
    }
    Ok([p, q])
}

fn dot_parts(a: &[C64], b: &[C64]) -> [Dd; 2] {
    let mut sum = [Dd::ZERO; 2];
    for (&a, &b) in a.iter().zip(b) {
        let [ar, ai, br, bi] = [a.re, a.im, b.re, b.im].map(Dd::from_f64);
        sum[0] = sum[0] + ar * br - ai * bi;
        sum[1] = sum[1] + ar * bi + ai * br;
    }
    sum
}
fn rounded(z: [Dd; 2]) -> C64 {
    C64::new(z[0].to_f64(), z[1].to_f64())
}
fn weighted(weights: [C64; 2], inputs: [C64; 2]) -> C64 {
    rounded(dot_parts(&weights, &inputs))
}
#[derive(Clone, Copy)]
enum Kernel {
    Real(HarmonicProduct),
    Proper(ProperHarmonicProduct),
}
fn accumulate(m: &mut [Dd; 5], p: [Dd; 2], q: [Dd; 2], kernel: Kernel) -> Result<(), DuctError> {
    let contract = |a, b| {
        match kernel {
            Kernel::Real(k) => k.mean_real_product_parts_dd(a, b),
            Kernel::Proper(k) => k.mean_real_product_parts_dd(a, b),
        }
        .map_err(window_error)
    };
    m[2] = m[2] + contract(p, p)?;
    m[3] = m[3] + contract(q, q)?;
    m[4] = m[4] + contract(p, q)?;
    Ok(())
}
fn window_error(_: fs_fft::harmonic_window::WindowError) -> DuctError {
    bad("spatial finite-window moment is invalid or unrepresentable")
}
fn poll(cancelled: &mut impl FnMut() -> bool, completed_steps: usize) -> Result<(), DuctError> {
    if cancelled() {
        Err(DuctError::Cancelled { completed_steps })
    } else {
        Ok(())
    }
}
