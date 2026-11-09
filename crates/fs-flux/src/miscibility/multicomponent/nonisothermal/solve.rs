//! Original coupled face rates and synchronous conservative thermal updates.

use super::*;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::solution::multicomponent::MAX_REGULAR_SOLUTION_SPECIES_V1;

/// One immutable right-hand side, shared by full/first-half trials and by
/// retries at the same accepted state. The friction workspace is reused across
/// every face; no independent binary approximation replaces the coupled solve.
pub(in crate::miscibility::multicomponent) struct Rates {
    fluxes: Vec<f64>,
    energy_fluxes: Vec<f64>,
    evaluations: usize,
    maximum_gradient_residual_per_m: f64,
    maximum_molar_frame_residual_mol_m2_s: f64,
    entropy_production_rate_w_per_k: f64,
}
impl Rates {
    pub(in crate::miscibility::multicomponent) const fn evaluations(&self) -> usize {
        self.evaluations
    }
}

pub(super) fn step(
    initial: &MulticomponentNonisothermalField,
    duration_s: f64,
    controls: MulticomponentNonisothermalSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MulticomponentNonisothermalCandidate> {
    poll(&mut cancelled)?;
    admit_controls(controls)?;
    admit_duration(initial, duration_s)?;
    let rates = rates(initial, controls, &mut cancelled)?;
    candidate_from_rates(initial, duration_s, &rates, controls, &mut cancelled)
}

fn admit_duration(initial: &MulticomponentNonisothermalField, duration_s: f64) -> Result<f64> {
    let elapsed = finite(initial.elapsed_s() + duration_s, "thermal elapsed time")?;
    if !positive(duration_s) || elapsed <= initial.elapsed_s() {
        return Err(MulticomponentNonisothermalError::Invalid(
            "positive representable physical interval required",
        ));
    }
    Ok(elapsed)
}

pub(super) fn rates(
    initial: &MulticomponentNonisothermalField,
    controls: MulticomponentNonisothermalSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Rates> {
    rates_with_gradient(initial, None, controls, cancelled)
}

/// Extra cell-major entropy-conjugate chemical potentials divided by R. Only
/// the interfacial owner supplies these; source efforts, friction, reciprocal
/// heat carriage and entropy diagnostics remain in this one original kernel.
pub(in crate::miscibility::multicomponent) fn rates_with_gradient(
    initial: &MulticomponentNonisothermalField,
    extra_reduced_potentials: Option<&[f64]>,
    controls: MulticomponentNonisothermalSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Rates> {
    poll(cancelled)?;
    admit_controls(controls)?;
    if let Some(potentials) = extra_reduced_potentials {
        if potentials.len() != initial.component_moles().len() {
            return Err(MulticomponentNonisothermalError::Invalid(
                "one extra thermal chemical potential per cell and original species required",
            ));
        }
        for row in potentials.chunks(initial.network().species_count()) {
            poll(cancelled)?;
            if row.iter().any(|value| !value.is_finite()) {
                return Err(MulticomponentNonisothermalError::Invalid(
                    "finite interfacial thermal chemical potentials required",
                ));
            }
        }
    }
    let network = initial.network();
    let n = network.species_count();
    let mut fluxes = vec![0.0; network.faces().len() * n];
    let mut energy_fluxes = vec![0.0; network.faces().len()];
    let mut workspace = network
        .transport()
        .workspace()
        .map_err(MulticomponentMiscibilityError::from)?;
    let mut evaluations = 0usize;
    let mut maximum_gradient_residual_per_m = 0.0_f64;
    let mut maximum_molar_frame_residual_mol_m2_s = 0.0_f64;
    let mut entropy = Sum::default();
    for (face_index, face) in network.faces().iter().enumerate() {
        poll(cancelled)?;
        if face.area_m2 == 0.0 {
            continue;
        }
        let [i, j] = face.cells;
        let coefficients = &initial.transport_faces()[face_index];
        let efforts = initial
            .caloric_model()
            .transport_efforts([initial.point(i), initial.point(j)], &mut *cancelled)?;
        let [tl, tr] = [
            initial.temperatures[i].value(),
            initial.temperatures[j].value(),
        ];
        // Difference first preserves adjacent temperatures; dividing in this
        // order avoids forming T_left*T_right, which can overflow or underflow.
        let delta_beta = finite(
            ((tl - tr) / tl) / tr,
            "thermal reciprocal-temperature difference",
        )?;
        let mut reduced = [0.0; MAX_REGULAR_SOLUTION_SPECIES_V1];
        let mut driving = [Sum::default(); MAX_REGULAR_SOLUTION_SPECIES_V1];
        let mut enthalpy_contrasts = [0.0; MAX_REGULAR_SOLUTION_SPECIES_V1];
        let q0 = coefficients.heat_of_transport_j_per_mol[0];
        let a0 = efforts.chemical_efforts_j_per_mol_k[0];
        for s in 0..n {
            poll(cancelled)?;
            // Common species gauges are removed analytically, not by a
            // dependent-species correction or an invented thermodiffusivity.
            let q = finite(
                coefficients.heat_of_transport_j_per_mol[s] - q0,
                "heat-of-transport contrast",
            )?;
            reduced[s] = finite(
                (efforts.chemical_efforts_j_per_mol_k[s] - a0) / R - (q / R) * delta_beta,
                "source-conjugate thermal chemical effort",
            )?;
            if let Some(potentials) = extra_reduced_potentials {
                let mut total = Sum::default();
                total.add(reduced[s])?;
                // Remove each cell's common species gauge before differencing.
                total.add(potentials[j * n + s] - potentials[j * n])?;
                total.add(-(potentials[i * n + s] - potentials[i * n]))?;
                reduced[s] = total.value()?;
            }
            enthalpy_contrasts[s] = finite(
                efforts.partial_enthalpy_contrasts_j_per_mol[s] + q,
                "original divided enthalpy and reciprocal heat contrast",
            )?;
        }
        for s in 0..n {
            poll(cancelled)?;
            for t in s + 1..n {
                let pair = super::super::solve::pair_force(
                    reduced[s] - reduced[t],
                    efforts.mean_mole_fractions[s],
                    efforts.mean_mole_fractions[t],
                    face.distance_m,
                )?;
                driving[s].add(pair)?;
                driving[t].add(-pair)?;
            }
        }
        let mut forces = [0.0; MAX_REGULAR_SOLUTION_SPECIES_V1];
        for s in 0..n {
            forces[s] = driving[s].value()?;
        }
        let result = network
            .transport()
            .evaluate_molar_driving_forces(
                &efforts.mean_mole_fractions,
                &forces[..n],
                network.concentration(),
                &network.owner.diffusivity_matrices[face.diffusivity_matrix],
                controls.maxwell_stefan,
                &mut workspace,
                &mut *cancelled,
            )
            .map_err(MulticomponentMiscibilityError::from)?;
        evaluations += 1;
        maximum_gradient_residual_per_m =
            maximum_gradient_residual_per_m.max(result.maximum_gradient_residual_per_m());
        maximum_molar_frame_residual_mol_m2_s =
            maximum_molar_frame_residual_mol_m2_s.max(result.molar_frame_residual_mol_m2_s().abs());
        let fourier = if coefficients.thermal_conductivity_w_m_k == 0.0 {
            0.0
        } else {
            super::super::solve::pair_force(
                tl - tr,
                coefficients.thermal_conductivity_w_m_k,
                1.0,
                face.distance_m,
            )?
        };
        let mut mean = Sum::default();
        for s in 0..n {
            mean.add(efforts.mean_mole_fractions[s] * reduced[s])?;
        }
        let mean = mean.value()?;
        let mut molecular_entropy = Sum::default();
        let mut energy = Sum::default();
        energy.add(fourier)?;
        for s in 0..n {
            let flux = result.molar_fluxes()[s];
            fluxes[face_index * n + s] = flux;
            energy.add(finite(
                flux * enthalpy_contrasts[s],
                "molecular enthalpy and Dufour energy flux",
            )?)?;
            molecular_entropy.add(finite(
                -flux * (reduced[s] - mean),
                "thermal chemical entropy production",
            )?)?;
        }
        energy_fluxes[face_index] = energy.value()?;
        let mut production = Sum::default();
        production.add(finite(
            R * molecular_entropy.value()?,
            "molecular face entropy production",
        )?)?;
        production.add(finite(
            fourier * delta_beta,
            "Fourier face entropy production",
        )?)?;
        let production = finite(
            production.value()? * face.area_m2,
            "thermal face entropy production rate",
        )?;
        if production < 0.0 {
            return Err(MulticomponentNonisothermalError::Invalid(
                "negative source-conjugate constitutive entropy production",
            ));
        }
        entropy.add(production)?;
    }
    poll(cancelled)?;
    Ok(Rates {
        fluxes,
        energy_fluxes,
        evaluations,
        maximum_gradient_residual_per_m,
        maximum_molar_frame_residual_mol_m2_s,
        entropy_production_rate_w_per_k: entropy.value()?,
    })
}

pub(super) fn candidate_from_rates(
    initial: &MulticomponentNonisothermalField,
    duration_s: f64,
    rates: &Rates,
    controls: MulticomponentNonisothermalSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<MulticomponentNonisothermalCandidate> {
    let trial = conservative_trial_from_rates(initial, duration_s, rates, controls, cancelled)?;
    let candidate = recover_candidate(initial, trial, rates, controls, cancelled)?;
    admit_energy_entropy(
        candidate.step.energy_residual_j,
        candidate.step.entropy_change_j_per_k,
        controls,
    )?;
    poll(cancelled)?;
    Ok(candidate)
}

/// Private synchronous inventory/energy stage. The interfacial caller may
/// deduct its half-edge storage changes from the bulk enthalpy changes before
/// the SAME original caloric inverse runs. No trial is a public accepted state.
pub(in crate::miscibility::multicomponent) struct Trial {
    pub(in crate::miscibility::multicomponent) field: MulticomponentNonisothermalField,
    pub(in crate::miscibility::multicomponent) energy_changes: Vec<Sum>,
    duration_s: f64,
    transferred_mol: Vec<Vec<f64>>,
    transferred_energy_j: Vec<f64>,
}

pub(in crate::miscibility::multicomponent) fn conservative_trial_from_rates(
    initial: &MulticomponentNonisothermalField,
    duration_s: f64,
    rates: &Rates,
    controls: MulticomponentNonisothermalSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Trial> {
    poll(cancelled)?;
    admit_controls(controls)?;
    let elapsed = admit_duration(initial, duration_s)?;
    let network = initial.network();
    let n = network.species_count();
    let mut field = initial.clone();
    let mut energy_changes = vec![Sum::default(); network.cells().len()];
    let mut transferred_mol = Vec::with_capacity(network.faces().len());
    let mut transferred_energy_j = Vec::with_capacity(network.faces().len());
    for (face_index, face) in network.faces().iter().enumerate() {
        poll(cancelled)?;
        let [i, j] = face.cells;
        let mut transfers = Vec::with_capacity(n);
        for s in 0..n {
            let transfer = transfer(rates.fluxes[face_index * n + s], face.area_m2, duration_s)?;
            for (cell, change) in [(i, -transfer), (j, transfer)] {
                let index = cell * n + s;
                (
                    field.state.amounts[index],
                    field.state.amount_remainders[index],
                ) = crate::miscibility::solve::add_amount(
                    field.state.amounts[index],
                    field.state.amount_remainders[index],
                    change,
                )
                .map_err(|_| {
                    MulticomponentMiscibilityError::Arithmetic("paired thermal species update")
                })?;
            }
            transfers.push(transfer);
        }
        transferred_mol.push(transfers);
        let energy = transfer(rates.energy_fluxes[face_index], face.area_m2, duration_s)?;
        energy_changes[i].add(-energy)?;
        energy_changes[j].add(energy)?;
        transferred_energy_j.push(energy);
    }
    field.state.elapsed_s = elapsed;
    poll(cancelled)?;
    Ok(Trial {
        field,
        energy_changes,
        duration_s,
        transferred_mol,
        transferred_energy_j,
    })
}

/// Only an interface owner needs fractions before caloric recovery. Ordinary
/// homogeneous transport avoids this additional original-law evaluation.
pub(in crate::miscibility::multicomponent) fn prepare_trial_composition(
    trial: &mut Trial,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<()> {
    let network = trial.field.state.network.clone();
    let n = network.species_count();
    for cell in 0..network.cells().len() {
        poll(cancelled)?;
        let amounts = &trial.field.state.amounts[cell * n..(cell + 1) * n];
        for (species, &trial_mol) in amounts.iter().enumerate() {
            if !positive(trial_mol) {
                return Err(MulticomponentMiscibilityError::Positivity {
                    cell,
                    species,
                    trial_mol,
                }
                .into());
            }
        }
        let fractions = network.fractions(amounts, cancelled)?;
        trial.field.state.fractions[cell * n..(cell + 1) * n]
            .copy_from_slice(&fractions[..n]);
    }
    poll(cancelled)?;
    Ok(())
}

/// Recover all local bulk enthalpies and retain the original precision account.
/// The caller must audit its complete stored energy and entropy before exposing
/// the result: homogeneous H here, H plus interface energy for phase fields.
pub(in crate::miscibility::multicomponent) fn recover_candidate(
    initial: &MulticomponentNonisothermalField,
    trial: Trial,
    rates: &Rates,
    controls: MulticomponentNonisothermalSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<MulticomponentNonisothermalCandidate> {
    poll(cancelled)?;
    admit_controls(controls)?;
    let network = initial.network();
    let n = network.species_count();
    let Trial {
        mut field,
        energy_changes,
        duration_s,
        transferred_mol,
        transferred_energy_j,
    } = trial;
    let mut caloric_iterations = 0_u64;
    // Recover each cell only after all synchronous debits/credits are present.
    // A cell is an open subsystem: only the closed graph owns the entropy gate.
    for (cell, energy) in energy_changes.iter().enumerate() {
        poll(cancelled)?;
        let amounts = &field.state.amounts[cell * n..(cell + 1) * n];
        for (species, &trial_mol) in amounts.iter().enumerate() {
            if !positive(trial_mol) {
                return Err(MulticomponentMiscibilityError::Positivity {
                    cell,
                    species,
                    trial_mol,
                }
                .into());
            }
        }
        let energy = energy.value()?;
        if energy == 0.0 && amounts == &initial.state.amounts[cell * n..(cell + 1) * n] {
            continue;
        }
        let result = initial.caloric_model().homogeneous_transport_candidate(
            &initial.caloric_state(cell, cancelled)?,
            amounts,
            energy,
            MulticomponentLiquidTransportSolveV1 {
                temperature_bracket: initial.temperature_bracket(),
                temperature_tolerance_k: controls.temperature_tolerance_k,
                // The whole-field check remains authoritative; splitting the
                // scalar allowance keeps many cells from consuming it each.
                energy_tolerance_j: controls.energy_tolerance_j / network.cells().len() as f64,
                maximum_iterations: controls.maximum_caloric_iterations,
            },
            &mut *cancelled,
        )?;
        caloric_iterations = caloric_iterations
            .checked_add(result.iterations() as u64)
            .ok_or(MulticomponentMiscibilityError::Arithmetic(
                "thermal caloric-work count",
            ))?;
        field.temperatures[cell] = result.final_state().temperature();
        let fractions = result.final_state().mixing().mole_fractions();
        if fractions.iter().any(|&x| !positive(x)) {
            return Err(MulticomponentNonisothermalError::Invalid(
                "positive thermal component fraction is unrepresentable",
            ));
        }
        field.state.fractions[cell * n..(cell + 1) * n].copy_from_slice(fractions);
    }
    let (
        component_residual_mol,
        maximum_capacity_residual_mol,
        energy_residual_j,
        entropy_change_j_per_k,
    ) = initial.caloric_balances(&field, controls, cancelled)?;
    poll(cancelled)?;
    Ok(MulticomponentNonisothermalCandidate {
        field,
        step: MulticomponentNonisothermalStep {
            duration_s,
            transferred_mol,
            transferred_energy_j,
            maxwell_stefan_evaluations: rates.evaluations,
            maximum_gradient_residual_per_m: rates.maximum_gradient_residual_per_m,
            maximum_molar_frame_residual_mol_m2_s: rates.maximum_molar_frame_residual_mol_m2_s,
            caloric_iterations,
            component_residual_mol,
            maximum_capacity_residual_mol,
            energy_residual_j,
            entropy_change_j_per_k,
            entropy_production_rate_w_per_k: rates.entropy_production_rate_w_per_k,
        },
    })
}

fn transfer(flux: f64, area: f64, duration: f64) -> Result<f64> {
    let value = super::super::solve::product3(flux, area, duration)?;
    if flux != 0.0 && area != 0.0 && value == 0.0 {
        return Err(MulticomponentMiscibilityError::Arithmetic(
            "nonzero thermal transfer underflows",
        )
        .into());
    }
    Ok(value)
}
