//! Conservative molecular transport through branches and compact gas volumes.
//!
//! Each connected face uses the original branch molecular constitutive kernel.
//! A caller-specified distance from the physical port to the stationary,
//! perfectly mixed plenum bulk combines with the branch half-cell distance.
//! This is an explicit compact mixing-length closure, not a resolved junction.
//! Species diffusion, Fourier heat and partial enthalpy/viscous energy flux
//! enter the two owners with opposite signs. The plenum's fixed support receives
//! axial momentum; different branch axes are never added into scalar momentum.
//!
//! The inherited constant-area, dilute-gas and estimated-parabolic restrictions
//! remain in force. No Soret/Dufour model or entropy certificate is implied.

use super::{MixtureNetwork1d, NetworkTotals, PlenumConnection, PortExchange, empty_exchange};
use crate::multicomponent::dissipation::{
    DissipationControl, DissipationError, DissipationReport, TransportLaw,
};
use crate::multicomponent::{
    Boundaries, Boundary, MAX_COMPONENTS, Mesh1d, MixtureGas, Primitive, Totals, allocated,
    positive,
};
use fs_qty::{Length, Time};
use fs_thermochem::transport::maxwell_stefan::{MaxwellStefanV1, MaxwellStefanWorkspaceV1};

mod evolve;

/// Complete molecular interval; the accepted network clock is unchanged.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NetworkDissipationReport {
    /// Complete requested molecular duration [s].
    pub elapsed: Time,
    /// Accepted globally synchronized SSPRK2 substeps.
    pub accepted_steps: usize,
    /// Rejected and halved complete network substeps.
    pub rejected_steps: usize,
    /// Maximum branch or compact-volume estimated parabolic number.
    pub maximum_parabolic_number: f64,
    /// Maximum original Maxwell--Stefan gradient residual [1/m].
    pub maximum_gradient_residual_per_m: f64,
    /// Maximum absolute residual in the mass-averaged diffusion frame [kg/(m² s)].
    pub maximum_mass_frame_residual_kg_m2_s: f64,
    /// Original combined branch and plenum storage.
    pub before: NetworkTotals,
    /// Accepted combined branch and plenum storage.
    pub after: NetworkTotals,
    /// Inflow through unconnected branch endpoints only.
    pub boundary_change: NetworkTotals,
    /// `after - before - boundary_change`, retaining each species separately.
    pub conservation_residual: NetworkTotals,
    /// Physical energy residual including original species references [J].
    pub absolute_energy_residual_j: f64,
}

/// Both exact-duration molecular halves of a complete coupled network step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SplitNetworkDissipationReport {
    /// Molecular half before advection or the chemistry/advection composition.
    pub first: NetworkDissipationReport,
    /// Molecular half after advection or the chemistry/advection composition.
    pub second: NetworkDissipationReport,
}

struct BranchShape {
    mesh: Mesh1d,
    boundaries: Boundaries,
    molecular_boundaries: Boundaries,
    widths: Vec<f64>,
    cell_start: usize,
    face_start: usize,
}

struct Workspace {
    current: Vec<f64>,
    stage: Vec<f64>,
    candidate: Vec<f64>,
    primitive: Vec<Primitive>,
    heat_capacity_v: Vec<f64>,
    faces: Vec<f64>,
    face_rates: Vec<f64>,
    binary_diffusivities: Vec<f64>,
    diffusion: MaxwellStefanWorkspaceV1,
    plenum_rhs: Vec<f64>,
    plenum_corrections: Vec<f64>,
    plenum_rates: Vec<f64>,
    endpoints_a: Vec<f64>,
    endpoint_changes: Vec<f64>,
    endpoint_corrections: Vec<f64>,
    branch_before: Vec<Totals>,
    branch_report_candidate: Vec<DissipationReport>,
    port_exchange_candidate: Vec<PortExchange>,
}

/// Prepared original-source molecular operator for a fixed network topology.
///
/// Every branch must have constant cross section. Each connection retains its
/// own strictly positive plenum-side mixing distance. Scratch fields, one
/// constitutive owner, and reduction/report storage are allocated once and
/// reused. A requested interval publishes all owners and reports together.
pub struct NetworkDissipationOperator<L: TransportLaw> {
    gas: MixtureGas,
    branches: Vec<BranchShape>,
    connections: Vec<PlenumConnection>,
    endpoint_connections: Vec<[Option<usize>; 2]>,
    plenum_volumes: Vec<f64>,
    plenum_lengths: Vec<Length>,
    plenum_cell_start: usize,
    law: L,
    maxwell_stefan: MaxwellStefanV1,
    mass_ratio: f64,
    work: Workspace,
    last_report: Option<NetworkDissipationReport>,
    trial_report: Option<NetworkDissipationReport>,
    port_exchanges: Vec<PortExchange>,
    branch_reports: Vec<DissipationReport>,
}
impl<L: TransportLaw> NetworkDissipationOperator<L> {
    /// Bind original geometry/source and a bulk-to-port distance per connection.
    ///
    /// Distances are physical lengths in retained connection order, not one
    /// value per plenum. They specify the unresolved compact-volume gradient
    /// closure; no neck geometry or empirical mixing loss is invented.
    pub fn new(
        network: &MixtureNetwork1d,
        law: L,
        plenum_lengths: Vec<Length>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        poll(&mut cancelled)?;
        let gas = network.gas();
        if law.source() != gas.source() {
            return Err(DissipationError::Invalid(
                "network and molecular law have different original species banks",
            ));
        }
        if plenum_lengths.len() != network.connections.len()
            || plenum_lengths
                .iter()
                .any(|length| !positive(2.0 * length.value()))
        {
            return Err(DissipationError::Invalid(
                "one positive plenum mixing length per molecular connection required",
            ));
        }
        let n = gas.species_count();
        let m = n + 2;
        let mut branches = Vec::new();
        branches
            .try_reserve_exact(network.branches.len())
            .map_err(|_| DissipationError::Invalid("molecular branch allocation"))?;
        let mut cells = 0;
        let mut faces = 0;
        for (index, branch) in network.branches.iter().enumerate() {
            poll(&mut cancelled)?;
            let area = branch
                .mesh
                .uniform_cross_section()
                .ok_or(DissipationError::Invalid(
                    "network molecular transport requires constant-area branches",
                ))?;
            let mut widths = allocated(branch.cell_count(), 0.0)?;
            for (i, width) in widths.iter_mut().enumerate() {
                poll(&mut cancelled)?;
                *width = branch.mesh.uniform_cell_width().map_or_else(
                    || branch.mesh.positions()[i + 1] - branch.mesh.positions()[i],
                    |w| w.value(),
                );
                let expected = area * *width;
                if !positive(expected)
                    || (branch.mesh.volumes()[i] - expected).abs()
                        > 64.0 * f64::EPSILON * expected.max(branch.mesh.volumes()[i])
                {
                    return Err(DissipationError::Invalid(
                        "molecular branch volume differs from area times cell width",
                    ));
                }
            }
            let molecular_boundaries = match branch.boundaries {
                Boundaries::Periodic => Boundaries::Periodic,
                Boundaries::Ends { left, right } => {
                    let left = if network.endpoint_connections[index][0].is_some() {
                        Boundary::Transmissive
                    } else {
                        left
                    };
                    let right = if network.endpoint_connections[index][1].is_some() {
                        Boundary::Transmissive
                    } else {
                        right
                    };
                    if matches!(left, Boundary::Reservoir(_))
                        || matches!(right, Boundary::Reservoir(_))
                    {
                        return Err(DissipationError::Invalid(
                            "external molecular reservoir needs a species-bearing owner",
                        ));
                    }
                    Boundaries::Ends { left, right }
                }
            };
            branches.push(BranchShape {
                mesh: branch.mesh.clone(),
                boundaries: branch.boundaries,
                molecular_boundaries,
                widths,
                cell_start: cells,
                face_start: faces,
            });
            cells += branch.cell_count();
            faces += branch.cell_count() + 1;
        }
        let plenum_cell_start = cells;
        cells += network.plenums.len();
        let species = gas
            .source()
            .models()
            .iter()
            .map(|s| (s.species().clone(), s.molar_mass()))
            .collect();
        let maxwell_stefan = MaxwellStefanV1::new(species)?;
        let minimum = gas
            .source()
            .models()
            .iter()
            .map(|s| s.molar_mass().value())
            .fold(f64::INFINITY, f64::min);
        let maximum = gas
            .source()
            .models()
            .iter()
            .map(|s| s.molar_mass().value())
            .fold(0.0, f64::max);
        let mass_ratio = maximum / minimum;
        if !positive(mass_ratio) {
            return Err(DissipationError::Invalid(
                "molecular mass ratio is unrepresentable",
            ));
        }
        let branch_count = branches.len();
        let port_count = network.connections.len();
        let plenum_count = network.plenums.len();
        let work = Workspace {
            current: allocated(cells * m, 0.0)?,
            stage: allocated(cells * m, 0.0)?,
            candidate: allocated(cells * m, 0.0)?,
            primitive: allocated(cells, Primitive::default())?,
            heat_capacity_v: allocated(cells, 0.0)?,
            faces: allocated(faces * m, 0.0)?,
            face_rates: allocated(faces, 0.0)?,
            binary_diffusivities: allocated(n * n, 0.0)?,
            diffusion: maxwell_stefan.workspace()?,
            plenum_rhs: allocated(plenum_count * m, 0.0)?,
            plenum_corrections: allocated(plenum_count * m, 0.0)?,
            plenum_rates: allocated(plenum_count, 0.0)?,
            endpoints_a: allocated(branch_count * 2 * m, 0.0)?,
            endpoint_changes: allocated(branch_count * 2 * m, 0.0)?,
            endpoint_corrections: allocated(branch_count * 2 * m, 0.0)?,
            branch_before: allocated(branch_count, super::empty_totals(n))?,
            branch_report_candidate: allocated(branch_count, empty_report(n))?,
            port_exchange_candidate: allocated(port_count, empty_exchange(n))?,
        };
        let owner = Self {
            gas: gas.clone(),
            branches,
            connections: network.connections.clone(),
            endpoint_connections: network.endpoint_connections.clone(),
            plenum_volumes: network.plenums.iter().map(|p| p.volume_m3).collect(),
            plenum_lengths,
            plenum_cell_start,
            law,
            maxwell_stefan,
            mass_ratio,
            work,
            last_report: None,
            trial_report: None,
            port_exchanges: allocated(port_count, empty_exchange(n))?,
            branch_reports: allocated(branch_count, empty_report(n))?,
        };
        poll(&mut cancelled)?;
        Ok(owner)
    }

    /// The exact original caloric/source bank.
    pub const fn gas(&self) -> &MixtureGas {
        &self.gas
    }
    /// Original authored transport law and its admitted range.
    pub const fn law(&self) -> &L {
        &self.law
    }
    /// Admitted compact-side gradient distances, in connection order.
    pub fn plenum_lengths(&self) -> &[Length] {
        &self.plenum_lengths
    }
    /// Most recent complete accepted molecular interval.
    pub const fn last_report(&self) -> Option<NetworkDissipationReport> {
        self.last_report
    }
    /// Accepted molecular exchanges into each plenum in connection order.
    pub fn accepted_port_exchanges(&self) -> Option<&[PortExchange]> {
        self.last_report.map(|_| self.port_exchanges.as_slice())
    }
    /// Accepted molecular branch balances, each in its own axial coordinate.
    pub fn accepted_branch_reports(&self) -> Option<&[DissipationReport]> {
        self.last_report.map(|_| self.branch_reports.as_slice())
    }

    /// Advance all branches and plenums for exactly the requested duration.
    ///
    /// The network clock is unchanged. Its previous advection-only report is
    /// cleared on a positive completed interval; this operator owns the new
    /// molecular reports. A refusal preserves all fields and accepted reports.
    pub fn advance_for(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        controls: DissipationControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NetworkDissipationReport, DissipationError> {
        let report = self.advance_for_trial(network, duration, controls, &mut cancelled)?;
        self.publish_trial_report();
        Ok(report)
    }

    pub(super) fn advance_for_trial(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        controls: DissipationControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NetworkDissipationReport, DissipationError> {
        self.trial_report = None;
        let report = evolve::advance(self, network, duration, controls, cancelled)?;
        self.trial_report = Some(report);
        Ok(report)
    }
    pub(super) fn publish_trial_report(&mut self) {
        self.last_report = Some(self.trial_report.take().expect("completed molecular trial"));
        self.port_exchanges
            .copy_from_slice(&self.work.port_exchange_candidate);
        self.branch_reports
            .copy_from_slice(&self.work.branch_report_candidate);
    }
    pub(super) fn trial_port_exchanges(&self) -> &[PortExchange] {
        &self.work.port_exchange_candidate
    }
}

fn empty_report(n: usize) -> DissipationReport {
    let zero = super::empty_totals(n);
    DissipationReport {
        elapsed: Time::new(0.0),
        accepted_steps: 0,
        rejected_steps: 0,
        maximum_parabolic_number: 0.0,
        maximum_gradient_residual_per_m: 0.0,
        maximum_mass_frame_residual_kg_m2_s: 0.0,
        before: zero,
        after: zero,
        boundary_change: zero,
        conservation_residual: zero,
        absolute_energy_residual_j: 0.0,
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DissipationError> {
    crate::multicomponent::poll(cancelled).map_err(DissipationError::from)
}

#[cfg(test)]
pub(super) mod tests;
