//! One original-source chemistry update for axial, planar and spatial cells.

use super::*;
use crate::multicomponent::{
    Cell, MAX_SPATIAL_COMPONENTS, MixtureGas, NormalCell, Primitive, SpatialCell,
};
use fs_qty::{Energy, Temperature, Volume};

/// Prepared chemical source for repeated conservative control-volume updates.
///
/// The immutable gas and mechanism must have the identical original species
/// bank. The same conversion, caloric reconstruction and energy-reference
/// accounting serve axial flows, networks and planar/volumetric fields. Every momentum
/// component and physical total energy are held fixed during chemistry.
/// Successful advances allocate no heap storage. A refused cell leaves both
/// its conserved row and the caller's accumulated report unchanged.
pub struct CellChemistry {
    gas: MixtureGas,
    reactor: IsochoricReactor,
    amounts: Vec<Amount>,
}

impl CellChemistry {
    /// Bind the complete mechanism and its original admitted caloric path.
    pub fn new(
        gas: MixtureGas,
        mechanism: IdealGasKineticsV1,
        controls: IsochoricKineticsSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ReactiveError> {
        poll(&mut cancelled)?;
        if gas.source() != mechanism.thermodynamics() {
            return Err(ReactiveError::Invalid(
                "flow and mechanism have different original species banks",
            ));
        }
        let domain = gas.domain();
        if controls.minimum_temperature_k < domain.minimum_temperature
            || controls.maximum_temperature_k > domain.maximum_temperature
        {
            return Err(ReactiveError::Invalid(
                "chemistry interval exceeds the flow caloric path",
            ));
        }
        let amounts = allocated(gas.species_count(), Amount::new(0.0))?;
        let reactor = IsochoricReactor::new(mechanism, controls, &mut cancelled)
            .map_err(|source| chemistry_error(0, source))?;
        poll(&mut cancelled)?;
        Ok(Self {
            gas,
            reactor,
            amounts,
        })
    }

    /// Exact immutable caloric owner whose conservative rows this source accepts.
    pub const fn gas(&self) -> &MixtureGas {
        &self.gas
    }

    /// Original reaction laws and canonical species; never a rebuilt surrogate.
    pub fn mechanism(&self) -> &IdealGasKineticsV1 {
        self.reactor.mechanism()
    }

    /// Original reactor's time method, source interval, tolerances and work bounds.
    pub fn controls(&self) -> IsochoricKineticsSolve {
        self.reactor.controls()
    }

    /// Numerical heap capacity, including the full original reactor and the
    /// conversion buffer. Immutable gas/source models and prepared caloric
    /// coefficients, fixed-size owners and bounded stack scratch are excluded.
    pub fn workspace_bytes(&self) -> usize {
        self.reactor.workspace_bytes() + self.amounts.capacity() * core::mem::size_of::<Amount>()
    }

    /// Advance `[rho_i..., px, E_ref]` at fixed volume, momentum and total energy.
    /// `cell` identifies source refusals; it does not select physics.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_axial(
        &mut self,
        q: &mut [f64],
        volume: Volume,
        duration: Time,
        cell: usize,
        report: &mut ChemistryReport,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), ReactiveError> {
        react_components::<1>(
            &self.gas,
            q,
            volume.value(),
            cell,
            &mut self.reactor,
            &mut self.amounts,
            duration.value(),
            report,
            &mut cancelled,
        )
    }

    /// Advance `[rho_i..., px, py, E_ref]`; both kinetic components participate
    /// in the original internal-energy target and endpoint admission.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_planar(
        &mut self,
        q: &mut [f64],
        volume: Volume,
        duration: Time,
        cell: usize,
        report: &mut ChemistryReport,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), ReactiveError> {
        react_components::<2>(
            &self.gas,
            q,
            volume.value(),
            cell,
            &mut self.reactor,
            &mut self.amounts,
            duration.value(),
            report,
            &mut cancelled,
        )
    }

    /// Advance `[rho_i..., px, py, pz, E_ref]` at fixed physical volume and all
    /// three momenta. The source target subtracts all three kinetic components;
    /// endpoint admission independently checks sensible energy and formation
    /// references before publishing the cell and caller's report together.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_spatial(
        &mut self,
        q: &mut [f64],
        volume: Volume,
        duration: Time,
        cell: usize,
        report: &mut ChemistryReport,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), ReactiveError> {
        react_components::<3>(
            &self.gas,
            q,
            volume.value(),
            cell,
            &mut self.reactor,
            &mut self.amounts,
            duration.value(),
            report,
            &mut cancelled,
        )
    }
}

/// Reduce only transverse kinetic energy, then reuse the exact original
/// caloric inversion (including its endpoint temperature hint).
fn caloric<const D: usize>(
    gas: &MixtureGas,
    q: &[f64],
    hint: Option<f64>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Primitive, f64), FlowError> {
    let n = gas.species_count();
    let reduced = match D {
        1 => Cell {
            partial_densities: &q[..n],
            momentum: q[n],
            energy: q[n + 1],
        },
        2 => gas.normal_reduction(NormalCell {
            partial_densities: &q[..n],
            normal_momentum: q[n],
            transverse_momentum: q[n + 1],
            energy: q[n + 2],
        })?,
        3 => gas.spatial_reduction(SpatialCell {
            partial_densities: &q[..n],
            normal_momentum: q[n],
            transverse_momentum: [q[n + 1], q[n + 2]],
            energy: q[n + 3],
        })?,
        _ => return Err(FlowError::Invalid("cell chemistry momentum dimension")),
    };
    let primitive = gas.decode_cell(reduced, hint, cancelled)?;
    let thermal = reduced.energy - 0.5 * reduced.momentum * primitive.velocity;
    Ok((primitive, thermal))
}

#[allow(clippy::too_many_arguments)]
pub(super) fn react_components<const D: usize>(
    gas: &MixtureGas,
    q: &mut [f64],
    volume: f64,
    cell: usize,
    reactor: &mut IsochoricReactor,
    amounts: &mut [Amount],
    duration: f64,
    report: &mut ChemistryReport,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ReactiveError> {
    poll(cancelled)?;
    let n = gas.species_count();
    let energy_index = n + D;
    let width = energy_index + 1;
    if q.len() != width
        || amounts.len() != n
        || !positive(volume)
        || !duration.is_finite()
        || duration < 0.0
    {
        return Err(ReactiveError::Invalid(
            "cell chemistry dimensions, volume or duration",
        ));
    }
    let controls = reactor.controls();
    let domain = gas.domain();
    let (primitive, thermal) = caloric::<D>(gas, q, None, cancelled)?;
    if !(controls.minimum_temperature_k..=controls.maximum_temperature_k)
        .contains(&primitive.temperature)
    {
        return Err(ReactiveError::Invalid(
            "cell temperature is outside chemistry interval",
        ));
    }
    let reference = sum(q[..n]
        .iter()
        .zip(gas.energy_references())
        .map(|(rho, e)| rho * e));
    let internal = (thermal + reference) * volume;
    let tolerance = controls.energy_tolerance_j
        + volume
            * (primitive.density * domain.absolute_energy_tolerance
                + domain.relative_energy_tolerance * thermal.abs());
    let initial_roundtrip = (internal / volume - reference - thermal) * volume;
    if !internal.is_finite() || !tolerance.is_finite() || initial_roundtrip.abs() > tolerance {
        return Err(ReactiveError::EnergyConsistency {
            cell,
            residual_j: initial_roundtrip,
            tolerance_j: tolerance,
        });
    }
    for k in 0..n {
        amounts[k] = Amount::new((q[k] * volume) / gas.source().models()[k].molar_mass().value());
    }
    let endpoint = reactor
        .advance_adiabatic_energy(
            amounts,
            Energy::new(internal),
            Volume::new(volume),
            Temperature::new(primitive.temperature),
            Time::new(duration),
            &mut *cancelled,
        )
        .map_err(|source| chemistry_error(cell, source))?;
    let mut next = [0.0; MAX_SPATIAL_COMPONENTS];
    for k in 0..n {
        next[k] =
            endpoint.amounts()[k].value() * gas.source().models()[k].molar_mass().value() / volume;
    }
    next[n..energy_index].copy_from_slice(&q[n..energy_index]);
    // This linear conservative gauge includes heat release exactly once.
    let reference_change = sum((0..n).map(|k| (next[k] - q[k]) * gas.energy_references()[k]));
    next[energy_index] = q[energy_index] - reference_change;
    let (_, next_thermal) = caloric::<D>(
        gas,
        &next[..width],
        Some(endpoint.temperature().value()),
        cancelled,
    )?;
    let balance = endpoint.balance();
    let mut source_thermal = [0.0; MAX_SPECIES];
    for k in 0..n {
        source_thermal[k] = (next[k] / gas.source().models()[k].molar_mass().value())
            * gas.source().models()[k]
                .molar_internal_energy_change(
                    Temperature::new(domain.minimum_temperature),
                    endpoint.temperature(),
                )
                .map_err(FlowError::from)?
                .value();
    }
    // Independent reconstruction at the reactor's temperature; the conserved
    // gauge identity alone cannot establish a thermodynamically valid endpoint.
    let reconstruction_residual =
        (next_thermal - sum(source_thermal[..n].iter().copied())) * volume;
    let reconstruction_tolerance = tolerance + balance.energy_residual_j.abs();
    if !reconstruction_residual.is_finite()
        || reconstruction_residual.abs() > reconstruction_tolerance
    {
        return Err(ReactiveError::EnergyConsistency {
            cell,
            residual_j: reconstruction_residual,
            tolerance_j: reconstruction_tolerance,
        });
    }
    // Small mass error can change every kinetic term at fixed momenta. Admit
    // the actual physical discrepancy under the same declared energy budget.
    let physical_residual = ((next_thermal - thermal) + reference_change) * volume;
    if !physical_residual.is_finite() || physical_residual.abs() > tolerance {
        return Err(ReactiveError::EnergyConsistency {
            cell,
            residual_j: physical_residual,
            tolerance_j: tolerance,
        });
    }
    let mut accepted_report = *report;
    fn count(a: usize, b: usize) -> Result<usize, ReactiveError> {
        a.checked_add(b)
            .ok_or(ReactiveError::Invalid("chemistry work counter overflow"))
    }
    accepted_report.cell_solves = count(report.cell_solves, 1)?;
    accepted_report.accepted_steps = count(report.accepted_steps, balance.accepted_steps)?;
    accepted_report.rejected_steps = count(report.rejected_steps, balance.rejected_steps)?;
    accepted_report.nonlinear_iterations =
        count(report.nonlinear_iterations, balance.nonlinear_iterations)?;
    accepted_report.maximum_energy_residual_j = report
        .maximum_energy_residual_j
        .max(balance.energy_residual_j.abs());
    accepted_report.maximum_reconstruction_residual_j = report
        .maximum_reconstruction_residual_j
        .max(reconstruction_residual.abs());
    accepted_report.entropy_change_j_k += balance.entropy_change_j_k;
    for k in 0..n {
        accepted_report.species_change_kg[k] += (next[k] - q[k]) * volume;
    }
    accepted_report.referenced_energy_change_j += (next[energy_index] - q[energy_index]) * volume;
    if ![
        report.maximum_energy_residual_j,
        report.maximum_reconstruction_residual_j,
        accepted_report.entropy_change_j_k,
        accepted_report.referenced_energy_change_j,
    ]
    .iter()
    .all(|x| x.is_finite())
        || !accepted_report
            .species_change_kg
            .iter()
            .all(|x| x.is_finite())
    {
        return Err(ReactiveError::Invalid(
            "non-finite accumulated chemistry balance",
        ));
    }
    poll(cancelled)?;
    q.copy_from_slice(&next[..width]);
    *report = accepted_report;
    Ok(())
}
