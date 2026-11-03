//! Simultaneous finite-volume absorption and exsolution of dilute gas mixtures.
//!
//! At fixed T and gas volume, each partial pressure is n_i RT/V_g. Independent
//! dilute Henry laws therefore retain the original scalar affine solution even
//! while other species move in the opposite direction. All species sample the
//! SAME initial parcel, and one native boundary reconstruction updates its
//! composition, pressure, calorics, transport and acoustics. Sequentially
//! committing scalar exchanges would expose unphysical intermediate gas states.
//!
//! Liquid inventories share one physical volume and clock. Their volumes are
//! never added. Source laws must explicitly describe each solute in that same
//! solvent; this owner cannot infer mixture applicability from species names.
//! Cross-solute interactions, reactions, moving volumes and adiabatic heating
//! are outside the independent dilute model. Both thermostats remain explicit.

use super::*;
use fs_math::dd::Dd;

/// Bounds preparation, source resolution and native species reconstruction.
pub const MAX_DISSOLVING_SPECIES: usize = 64;

/// One accepted interval with a single shared gas storage/energy account.
#[derive(Clone, Debug, PartialEq)]
pub struct DissolutionMixtureStep {
    /// Original law responses in canonical species order, independent of input order.
    pub exchanges: Vec<DissolutionExchange>,
    pub gas: GasSpeciesBoundaryStep,
    /// Conservative total-pressure envelope used for original source admission.
    /// Opposite transfers can make total pressure nonmonotone within an interval.
    pub pressure_envelope_pa: [f64; 2],
    pub dissolved_internal_energy_change_j: f64,
    pub thermostat_heat_into_gas_j: f64,
    pub thermostat_heat_into_liquid_j: f64,
    pub thermostat_heat_into_system_j: f64,
    /// Gas plus all dissolved storage changes minus both thermostat inputs.
    pub energy_residual_j: f64,
}

fn total(values: impl Iterator<Item = f64>) -> Result<f64> {
    finite(
        values
            .fold(Dd::ZERO, |sum, value| sum + Dd::from_f64(value))
            .to_f64(),
    )
}

/// Advance 1..=64 independently sourced species in one liquid and finite gas.
///
/// The three slices use the same species order. Every inventory must hold the
/// exact same liquid volume and accepted time, and every law must retain the
/// same original gas catalog and temperature. Each species appears once.
/// Different diffusivities and physical films are allowed, including zero area.
///
/// `maximum_concentration_mol_m3` bounds the SUM of exchanging solutes across
/// the liquid bulk, liquid-side interface and admitted quasi-steady films.
/// Each species contributes its maximum across both locations and times, so
/// unequal film thicknesses and opposing transfer cannot hide an interior
/// concentration peak. The per-species source domains also apply.
/// A conservative envelope of each monotone partial-pressure change
/// admits the complete interval's total pressure; if that envelope is too wide
/// for a source, shorten the interval. No trajectory is clipped or rescaled.
/// All histories remain unchanged on any refusal or cancellation.
pub fn advance_mixture(
    laws: &[HenryDissolution],
    gas: &mut ClosedGasParcel,
    liquids: &mut [DissolvedGasInventory],
    films: &[DissolutionFilm],
    end_time_s: f64,
    controls: DissolvedGasControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<DissolutionMixtureStep> {
    poll(&mut cancelled)?;
    if laws.is_empty()
        || laws.len() > MAX_DISSOLVING_SPECIES
        || laws.len() != liquids.len()
        || laws.len() != films.len()
    {
        return Err(invalid(
            "matching nonempty finite species, liquid and film layouts required",
        ));
    }
    let mut order: Vec<_> = (0..laws.len()).collect();
    order.sort_by(|&a, &b| laws[a].species.cmp(&laws[b].species));
    let first = order[0];
    for (position, &index) in order.iter().enumerate() {
        if (position > 0 && laws[order[position - 1]].species == laws[index].species)
            || liquids[index].liquid_volume_m3 != liquids[first].liquid_volume_m3
            || liquids[index].time_s != liquids[first].time_s
            || laws[index].temperature_k != laws[first].temperature_k
            || laws[index].catalog != laws[first].catalog
        {
            return Err(invalid(
                "unique species must share one liquid volume, clock, temperature and original gas catalog",
            ));
        }
    }
    let mut prepared = Vec::with_capacity(laws.len());
    let mut transfers = Vec::with_capacity(laws.len());
    for &index in &order {
        let proposal = laws[index].prepare(
            gas,
            liquids[index],
            films[index],
            end_time_s,
            controls,
            &mut cancelled,
        )?;
        transfers.push(proposal.transfer(&laws[index]));
        prepared.push(proposal);
    }
    let candidate = gas.exchange_species_boundary_relative(
        &laws[first].catalog,
        &transfers,
        controls.gas,
        &mut cancelled,
    )?;
    let mut exchanges = Vec::with_capacity(laws.len());
    for (&index, proposal) in order.iter().zip(prepared) {
        poll(&mut cancelled)?;
        exchanges.push(laws[index].finish(proposal, candidate.gas(), controls)?);
    }

    // Each species' partial pressure is monotone, but their sum need not be.
    // Summing separate extrema is conservative for unlike relaxation rates.
    let initial_pressure = gas.state().state().pressure;
    let final_pressure = candidate.gas().state().state().pressure;
    let changes = || {
        exchanges.iter().map(|exchange| {
            exchange.sources[1].partial_pressure_pa - exchange.sources[0].partial_pressure_pa
        })
    };
    let lower = total(core::iter::once(initial_pressure).chain(changes().filter(|v| *v < 0.0)))?
        .min(initial_pressure)
        .min(final_pressure);
    let upper = total(core::iter::once(initial_pressure).chain(changes().filter(|v| *v > 0.0)))?
        .max(initial_pressure)
        .max(final_pressure);
    if !positive(lower) || !positive(upper) {
        return Err(invalid(
            "mixture source-pressure envelope must remain positive and finite",
        ));
    }
    // Original gas transport sources also apply throughout the interval, not
    // just at the accepted endpoint. The full catalog includes species that
    // may first appear through exsolution during this very interval.
    for pressure in [lower, upper] {
        poll(&mut cancelled)?;
        let point = gas
            .point()
            .clone()
            .with_quantity("pressure", typed(QuantityKind::Pressure), pressure)
            .map_err(|_| invalid("unrepresentable gas source-pressure envelope"))?;
        laws[first]
            .catalog
            .evaluate(&point)
            .map_err(GasEnergyError::from)?;
    }
    for (&index, exchange) in order.iter().zip(&exchanges) {
        let law = &laws[index];
        for pressure in [lower, upper] {
            poll(&mut cancelled)?;
            for sample in &exchange.sources {
                for concentration in [
                    sample.concentration_mol_m3,
                    sample.equilibrium_concentration_mol_m3,
                ] {
                    let point = law
                        .source
                        .point
                        .clone()
                        .with_quantity("pressure", typed(QuantityKind::Pressure), pressure)
                        .and_then(|p| {
                            p.with_quantity(
                                "dissolved_concentration",
                                typed(QuantityKind::AmountConcentration),
                                concentration,
                            )
                        })
                        .map_err(|_| invalid("unrepresentable mixture source envelope"))?;
                    law.resolved.at_query_point(&law.source.material, &point)?;
                }
            }
        }
    }
    // A species' linear quasi-steady film is bounded by its own bulk and
    // interface concentrations. With unequal film thicknesses, opposing
    // species gradients can create a mixture peak inside those films even
    // when the total bulk and interface concentrations are both admissible.
    // Sum maxima across BOTH locations and times for each species; checking
    // the two mixture totals separately does not bound that interior peak.
    let concentration = total(exchanges.iter().map(|exchange| {
        exchange.sources.iter().fold(0.0_f64, |maximum, sample| {
            maximum
                .max(sample.concentration_mol_m3)
                .max(sample.equilibrium_concentration_mol_m3)
        })
    }))?;
    if concentration > controls.maximum_concentration_mol_m3 {
        return Err(DissolvedGasError::Concentration {
            actual_mol_m3: concentration,
            maximum_mol_m3: controls.maximum_concentration_mol_m3,
        });
    }
    let dissolved_change = total(
        exchanges
            .iter()
            .map(|e| e.dissolved_internal_energy_change_j),
    )?;
    let gas_heat = total(exchanges.iter().map(|e| e.thermostat_heat_into_gas_j))?;
    let liquid_heat = total(exchanges.iter().map(|e| e.thermostat_heat_into_liquid_j))?;
    let bath_heat = finite(gas_heat + liquid_heat)?;
    let energy_residual =
        finite((candidate.step().stored_energy_change_j + dissolved_change) - bath_heat)?;
    budget(
        "combined mixture internal energy",
        energy_residual,
        controls.gas.energy_tolerance_j,
    )?;
    let (next_gas, _, gas_step) = candidate.into_parts();
    let report = DissolutionMixtureStep {
        exchanges,
        gas: gas_step,
        pressure_envelope_pa: [lower, upper],
        dissolved_internal_energy_change_j: dissolved_change,
        thermostat_heat_into_gas_j: gas_heat,
        thermostat_heat_into_liquid_j: liquid_heat,
        thermostat_heat_into_system_j: bath_heat,
        energy_residual_j: energy_residual,
    };
    poll(&mut cancelled)?;
    for (&index, exchange) in order.iter().zip(&report.exchanges) {
        liquids[index] = exchange.state;
    }
    *gas = next_gas;
    Ok(report)
}
