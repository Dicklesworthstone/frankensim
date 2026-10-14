//! Evolving ordered electrical interfaces in the existing material midpoint step.
//!
//! Bulk properties and both surface temperatures belong to one evaluation epoch.
//! The original current kernel owns resistance assembly and contact heat; the
//! original midpoint driver owns iteration and atomic publication. Interfaces
//! have fixed geometry, no storage, and an explicitly supplied heat partition.

use super::*;
use crate::conjugate::{ElectricSetup, ElectricalContact};
use std::collections::BTreeSet;

/// Select an existing electrical topology owner, independently of material names.
#[derive(Debug, Clone, Copy)]
pub enum InterfaceCircuit {
    /// Externally wired conductors; existing voltage-reference rules apply.
    Wired,
    /// Electrode-free conductors; preserve numerical gauges and check every
    /// cell's current continuity against this finite nonnegative tolerance [A].
    Isolated {
        /// Maximum local current-continuity discrepancy [A].
        current_tolerance_a: f64,
    },
}

/// Coupled bulk/interface criteria. No field is an experimental error estimate.
#[derive(Debug, Clone, Copy)]
pub struct InterfaceMidpointOptions {
    /// Existing coordinate, source, work, transport and iteration controls.
    pub material: MaterialMidpointOptions,
    /// Relative contact-resistance mismatch |a-b|/max(a,b), in [0,1]. Two zeros
    /// match exactly; a zero versus a positive resistance has discrepancy one.
    pub resistance_relative_tolerance: f64,
    /// Maximum absolute heat-partition-fraction mismatch, in [0,1].
    pub partition_absolute_tolerance: f64,
    /// Wired or electrode-free geometry. Does not infer electrical conductivity.
    pub circuit: InterfaceCircuit,
}

/// Constitutive response at the two supplied surface states. Numerical values
/// are admitted before use. `material` retains the original concrete source
/// object, without requiring Clone or erasing source receipts.
#[derive(Debug, Clone, PartialEq)]
pub struct InterfaceResponse<I> {
    /// Original interface-law evaluation, or an explicitly synthetic payload.
    pub material: I,
    /// Area-specific electrical resistance [ohm m2], finite and nonnegative.
    /// Neither bulk resistivity nor thermal contact resistance.
    pub resistance_ohm_m2: f64,
    /// Fraction of contact-generated heat sent to the first ordered cell.
    /// Must be in [0,1]; no equal split or effusivity law is inferred.
    pub heat_fraction_to_first: f64,
}

/// An actual ordered interface evaluation, not a copy of the template values.
#[derive(Debug, Clone, PartialEq)]
pub struct InterfaceSample<I> {
    /// Index in the caller's original contact array.
    pub index: usize,
    /// Surface temperatures [K] in the original cell order.
    pub temperature_k: [f64; 2],
    /// Exact geometry/resistance/partition supplied to the existing current owner.
    pub contact: ElectricalContact,
    /// Original concrete constitutive payload returned for this epoch.
    pub material: I,
}

/// Maximum local differences; different interfaces cannot cancel each other.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InterfaceResidual {
    /// Relative resistance discrepancy, with an explicit zero convention.
    pub resistance_relative: f64,
    /// Absolute discrepancy in the fraction of heat into the first surface.
    pub partition_absolute: f64,
}

/// Accepted source-producing bulk and interface states, plus endpoint checks.
#[derive(Debug, Clone, PartialEq)]
pub struct InterfaceMidpointStep<F, I> {
    /// The original material midpoint fields, balances and bulk evaluations.
    pub bulk: MaterialMidpointStep<F>,
    /// Interface evaluations used to generate the applied currents and heat.
    pub applied_interfaces: Vec<InterfaceSample<I>>,
    /// Endpoint support checks, never substituted for the applied interfaces.
    pub endpoint_interfaces: Vec<InterfaceSample<I>>,
    /// Fresh constitutive checks at the actual old/new midpoint.
    pub interface_residual: InterfaceResidual,
}

/// Original physical state survives every error. External callback side effects
/// are caller-owned. A contact law's original error type is preserved as E.
#[derive(Debug, Clone, PartialEq)]
pub enum InterfaceMidpointError<E> {
    /// Invalid topology, coefficients or numerical controls.
    Invalid(&'static str),
    /// Original bulk material, current, heat, velocity or cancellation refusal.
    Bulk(MaterialMidpointError),
    /// The ordered interface law refused this exact pair of surface states.
    Interface {
        /// Original ordered contact index.
        index: usize,
        /// Physical surface cells in caller order.
        cells: [usize; 2],
        /// Exact temperatures supplied to the failed interface query [K].
        temperature_k: [f64; 2],
        /// Unmodified constitutive failure.
        source: E,
    },
    /// The shared nonlinear budget ended before all criteria passed.
    NonConvergence {
        /// Original iteration count and coordinate/source residuals.
        iteration: MagneticMidpointError,
        /// Most recent local transport recheck, when reached.
        last_transport: Option<TransportResidual>,
        /// Most recent interface recheck, when reached.
        last_interface: Option<InterfaceResidual>,
    },
}
impl<E: core::fmt::Display> core::fmt::Display for InterfaceMidpointError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(s) => write!(f, "material electrical interface: {s}"),
            Self::Bulk(e) => write!(f, "{e}"),
            Self::Interface { index, cells, temperature_k, source } =>
                write!(f, "electrical interface {index} {cells:?} at {temperature_k:?} K: {source}"),
            Self::NonConvergence { iteration, last_transport, last_interface } =>
                write!(f, "{iteration}; transport {last_transport:?}; interfaces {last_interface:?}"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for InterfaceMidpointError<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Bulk(e) => Some(e),
            Self::Interface { source, .. } => Some(source),
            Self::NonConvergence { iteration, .. } => Some(iteration),
            Self::Invalid(_) => None,
        }
    }
}
impl<E> From<MaterialMidpointError> for InterfaceMidpointError<E> {
    fn from(e: MaterialMidpointError) -> Self { Self::Bulk(e) }
}
impl<E> From<MagneticMidpointError> for InterfaceMidpointError<E> {
    fn from(e: MagneticMidpointError) -> Self { Self::Bulk(e.into()) }
}
impl<E> From<ElectricHeatError> for InterfaceMidpointError<E> {
    fn from(e: ElectricHeatError) -> Self { Self::Bulk(e.into()) }
}
impl<E> From<MaterialTransportError> for InterfaceMidpointError<E> {
    fn from(e: MaterialTransportError) -> Self { Self::Bulk(e.into()) }
}
impl<E> From<ChtError> for InterfaceMidpointError<E> {
    fn from(e: ChtError) -> Self { Self::Bulk(e.into()) }
}

fn admit_interfaces<E>(
    domain: &VoxelDomain, plane: usize, input: ElectricHeatInput<'_>,
    options: InterfaceMidpointOptions, gate: &CancelGate,
) -> Result<Option<f64>, InterfaceMidpointError<E>> {
    poll(gate)?;
    options.material.admit()?;
    if ![options.resistance_relative_tolerance, options.partition_absolute_tolerance]
        .iter().all(|v| v.is_finite() && (0.0..=1.0).contains(v))
    { return Err(InterfaceMidpointError::Invalid("interface tolerances must be in [0,1]")); }
    let isolated = match options.circuit {
        InterfaceCircuit::Wired => None,
        InterfaceCircuit::Isolated { current_tolerance_a } => {
            if !input.setup.electrodes.is_empty() || !current_tolerance_a.is_finite() || current_tolerance_a < 0.0 {
                return Err(InterfaceMidpointError::Invalid("isolated interfaces require no electrodes and nonnegative finite current tolerance"));
            }
            Some(current_tolerance_a)
        }
    };
    let cells = domain.cell_count();
    if input.contacts.len() > cells.saturating_mul(3) {
        return Err(InterfaceMidpointError::Invalid("contact count exceeds the voxel face count"));
    }
    let mut seen = BTreeSet::new();
    for (index, site) in input.contacts.iter().enumerate() {
        if index % 256 == 0 { poll(gate)?; }
        let [a, b] = site.cells;
        if a < plane || b < plane || a >= cells-plane || b >= cells-plane
            || !(0..6).any(|face| domain.neighbor(a, face) == Some(b))
            || !seen.insert((a.min(b), a.max(b)))
        {
            return Err(InterfaceMidpointError::Invalid("contacts must be unique adjacent interior faces, not thermal wall reservoirs"));
        }
    }
    Ok(isolated)
}

fn electrode_cell(domain: &VoxelDomain, setup: &ElectricSetup, cell: usize) -> bool {
    let [nx, ny, _] = domain.dims();
    let xyz = [cell % nx, (cell / nx) % ny, cell / nx / ny];
    setup.electrodes.iter().any(|e| (0..3).all(|a| e.lo[a] <= xyz[a] && xyz[a] < e.hi[a]))
}

struct Prepared<F, I> {
    bulk: Samples<F>,
    interfaces: Vec<InterfaceSample<I>>,
}

#[allow(clippy::too_many_arguments)]
fn advance<S: MovingSlab<D, Q>, F: FluidTransport, I, E, const D: usize, const Q: usize>(
    slab: &mut MaterialSlab<S>, domain: VoxelDomain, plane: usize,
    input: ElectricHeatInput<'_>, field: &[[f64; 3]], external: &[[f64; D]], depth_m: f64,
    options: InterfaceMidpointOptions, gate: &CancelGate,
    mut bulk: impl FnMut(usize, f64, &FluidTransportState) -> Result<ElectricalTransport<F>, MaterialError>,
    mut interface: impl FnMut(usize, &ElectricalContact, [&TransportSample<ElectricalTransport<F>>; 2]) -> Result<InterfaceResponse<I>, E>,
) -> Result<InterfaceMidpointStep<F, I>, InterfaceMidpointError<E>> {
    let isolated = admit_interfaces::<E>(&domain, plane, input, options, gate)?;
    let cells = domain.cell_count();
    let transport = slab.force_units.transport;
    let dt = slab.force_units.timestep_s;
    let mut last_transport = None;
    let mut last_interface = None;
    let result = driver::solve_prepared::<S, Prepared<F, I>, (TransportResidual, InterfaceResidual), InterfaceMidpointError<E>, D, Q>(
        slab, &domain, plane, input, field, external, options.material.coupling, gate,
        |mapping, temperature| {
            let bulk = sample_temperatures::<S, ElectricalTransport<F>, D, Q>(
                mapping, temperature, cells, plane, input.limits.thermal.temperature_range_k, gate, &mut bulk,
            )?;
            let mut interfaces = Vec::with_capacity(input.contacts.len());
            for (index, site) in input.contacts.iter().enumerate() {
                poll(gate)?;
                let pair = [&bulk[site.cells[0]-plane], &bulk[site.cells[1]-plane]];
                // Check every epoch, INCLUDING the endpoint: a contact cannot
                // become supported by insulation after the last current solve.
                let electrodes = site.cells.map(|c| electrode_cell(&domain, input.setup, c));
                if electrodes.iter().all(|e| *e)
                    || pair.iter().zip(electrodes).any(|(p, e)| !e && p.material.resistivity_ohm_m().is_none())
                { return Err(InterfaceMidpointError::Invalid("electrical contact requires bulk conduction on at least one side and no insulated side")); }
                let temperature_k = pair.map(|p| p.temperature_k);
                let response = interface(index, site, pair).map_err(|source| InterfaceMidpointError::Interface {
                    index, cells: site.cells, temperature_k, source,
                })?;
                let r = response.resistance_ohm_m2;
                let fraction = response.heat_fraction_to_first;
                if !r.is_finite() || r < 0.0 || !fraction.is_finite() || !(0.0..=1.0).contains(&fraction) {
                    return Err(InterfaceMidpointError::Invalid("interface response requires nonnegative finite resistance and heat partition in [0,1]"));
                }
                interfaces.push(InterfaceSample { index, temperature_k, material: response.material,
                    contact: ElectricalContact { cells: site.cells, resistance_ohm_m2: r, heat_fraction_to_first: fraction } });
            }
            poll(gate)?;
            Ok(Prepared { bulk, interfaces })
        },
        |temperature, velocity, applied| {
            let contacts: Vec<_> = applied.interfaces.iter().map(|s| s.contact).collect();
            let mut electrical = |cell: usize, _: f64, _: &FluidTransportState| {
                Ok(applied.bulk[cell-plane].material.resistivity_ohm_m())
            };
            Ok(source_at(&domain, plane, ElectricHeatInput { contacts: &contacts, ..input }, field,
                temperature, velocity, &transport, dt, isolated, gate, &mut electrical)?)
        },
        |old, power, force, applied| {
            let (lattice, heat) = candidate_with_samples::<S, ElectricalTransport<F>, D, Q>(
                old, power, force, depth_m, input.limits.thermal, &applied.bulk, gate,
            )?;
            Ok((MaterialSlab { lattice, top_temperature_k: old.top_temperature_k,
                temperature_difference_k: old.temperature_difference_k, force_units: old.force_units }, heat))
        },
        |used, checked| {
            let t = transport_residual(&used.bulk, &checked.bulk, gate)?;
            let mut r = InterfaceResidual { resistance_relative: 0.0, partition_absolute: 0.0 };
            for (index, (a, b)) in used.interfaces.iter().zip(&checked.interfaces).enumerate() {
                if index % 256 == 0 { poll(gate)?; }
                let (x, y) = (a.contact.resistance_ohm_m2, b.contact.resistance_ohm_m2);
                let relative = if x == 0.0 && y == 0.0 { 0.0 } else { (x-y).abs()/x.max(y) };
                r.resistance_relative = r.resistance_relative.max(relative);
                r.partition_absolute = r.partition_absolute.max((a.contact.heat_fraction_to_first-b.contact.heat_fraction_to_first).abs());
            }
            last_transport = Some(t); last_interface = Some(r);
            Ok((t.conductivity_relative <= options.material.transport_relative_tolerance
                && t.viscosity_relative <= options.material.transport_relative_tolerance
                && r.resistance_relative <= options.resistance_relative_tolerance
                && r.partition_absolute <= options.partition_absolute_tolerance, (t, r)))
        },
    );
    match result {
        Ok(r) => Ok(InterfaceMidpointStep {
            bulk: MaterialMidpointStep { step: r.step, applied: r.applied.bulk, endpoint: r.endpoint.bulk,
                transport_residual: r.material_residual.0 },
            applied_interfaces: r.applied.interfaces, endpoint_interfaces: r.endpoint.interfaces,
            interface_residual: r.material_residual.1,
        }),
        Err(InterfaceMidpointError::Bulk(MaterialMidpointError::Iteration(iteration @ MagneticMidpointError::NonConvergence { .. }))) =>
            Err(InterfaceMidpointError::NonConvergence { iteration, last_transport, last_interface }),
        Err(e) => Err(e),
    }
}

impl MaterialSlab<ThermalLbm> {
    /// Couple bulk and ordered electrical-interface laws at one 2D midpoint.
    /// `input.contacts` supplies fixed sites and callback templates; ONLY the
    /// returned interface resistance/partition is used in the current solve.
    /// The callback receives the exact bulk samples for both ordered surfaces,
    /// including temperatures and concrete material objects from the same epoch.
    /// It must not advance history on a speculative query. Missing source data
    /// returns its original typed error; it never becomes perfect contact.
    ///
    /// Zero B gives electrothermal feedback. Otherwise only Bz is supported.
    /// Depth [m] sets extensive energy, not 2D dynamics; forces are N/m3. Select
    /// wired/isolated behavior explicitly in options. Endpoint queries precede
    /// publication, and interface criteria apply even when current is zero.
    ///
    /// Existing fixed-density/Cp Boussinesq, imposed-field, split-integrator and
    /// cancellation boundaries remain. No evolving contact geometry, contact
    /// storage, thermal interface law, induction or total-energy claim is added.
    #[allow(clippy::too_many_arguments)]
    pub fn step_magnetothermal_midpoint_with_interfaces<F: FluidTransport, I, E>(
        &mut self, input: ElectricHeatInput<'_>, field: &[[f64; 3]], forces: &[[f64; 2]], depth_m: f64,
        options: InterfaceMidpointOptions, gate: &CancelGate,
        bulk: impl FnMut(usize, f64, &FluidTransportState) -> Result<ElectricalTransport<F>, MaterialError>,
        interface: impl FnMut(usize, &ElectricalContact, [&TransportSample<ElectricalTransport<F>>; 2]) -> Result<InterfaceResponse<I>, E>,
    ) -> Result<InterfaceMidpointStep<F, I>, InterfaceMidpointError<E>> {
        let (_, plane) = admit::<ThermalLbm, 2, 5>(self, input, forces, 1, self.lattice.grid.ny, gate)?;
        let dx = self.force_units.cell_size_m;
        if !positive(depth_m) || !positive((dx*dx)*depth_m) {
            return Err(InterfaceMidpointError::Invalid("invalid represented 2D volume"));
        }
        let domain = VoxelDomain::graded([vec![dx; self.lattice.grid.nx], vec![dx; self.lattice.grid.ny], vec![depth_m]])?;
        advance::<ThermalLbm, F, I, E, 2, 5>(self, domain, plane, input, field, forces, depth_m, options, gate, bulk, interface)
    }
}
impl MaterialSlab<ThermalLbm3> {
    /// 3D counterpart, using the same prepared-state midpoint, electrical
    /// contact kernel, original error types and atomic bulk/interface publication.
    #[allow(clippy::too_many_arguments)]
    pub fn step_magnetothermal_midpoint_with_interfaces<F: FluidTransport, I, E>(
        &mut self, input: ElectricHeatInput<'_>, field: &[[f64; 3]], forces: &[[f64; 3]],
        options: InterfaceMidpointOptions, gate: &CancelGate,
        bulk: impl FnMut(usize, f64, &FluidTransportState) -> Result<ElectricalTransport<F>, MaterialError>,
        interface: impl FnMut(usize, &ElectricalContact, [&TransportSample<ElectricalTransport<F>>; 2]) -> Result<InterfaceResponse<I>, E>,
    ) -> Result<InterfaceMidpointStep<F, I>, InterfaceMidpointError<E>> {
        let (_, plane) = admit::<ThermalLbm3, 3, 7>(self, input, forces, 2, self.lattice.grid.nz, gate)?;
        let dx = self.force_units.cell_size_m;
        if !positive((dx*dx)*dx) { return Err(InterfaceMidpointError::Invalid("invalid represented 3D volume")); }
        let g = &self.lattice.grid;
        let domain = VoxelDomain::new(g.nx, g.ny, g.nz, dx)?;
        advance::<ThermalLbm3, F, I, E, 3, 7>(self, domain, plane, input, field, forces, dx, options, gate, bulk, interface)
    }
}
