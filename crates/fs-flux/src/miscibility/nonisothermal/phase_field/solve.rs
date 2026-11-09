//! Frozen-temperature convex splitting and original local-energy recovery.

use super::*;
use fs_qty::AmountConcentration;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

struct ChemicalRoot {
    transfer: f64,
    amounts: [[f64; 2]; 2],
    remainders: [[f64; 2]; 2],
    fractions: [[f64; 2]; 2],
    iterations: usize,
}

/// Scratch is O(cells), reused across sparse incident-edge patches. Only the
/// touched entries are reset, so a face never scans the complete mesh.
struct Patch {
    energy: Vec<Sum>,
    marked: Vec<bool>,
    cells: Vec<usize>,
}
impl Patch {
    fn new(count: usize) -> Self {
        Self {
            energy: vec![Sum::default(); count],
            marked: vec![false; count],
            cells: Vec::with_capacity(count),
        }
    }
    fn add(&mut self, cell: usize, energy: f64) -> Result<()> {
        if !self.marked[cell] {
            self.marked[cell] = true;
            self.cells.push(cell);
        }
        self.energy[cell].add(energy)?;
        Ok(())
    }
    fn clear(&mut self) {
        for &cell in &self.cells {
            self.energy[cell] = Sum::default();
            self.marked[cell] = false;
        }
        self.cells.clear();
    }
}

pub(super) fn step(
    initial: &BinaryNonisothermalPhaseField,
    duration_s: f64,
    controls: BinaryNonisothermalSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<BinaryNonisothermalPhaseFieldCandidate> {
    poll(&mut cancelled)?;
    admit_controls(controls)?;
    let elapsed_s = finite(initial.elapsed_s() + duration_s, "phase-field elapsed time")?;
    if !positive(duration_s) || elapsed_s <= initial.elapsed_s() {
        return Err(BinaryNonisothermalError::Invalid(
            "positive representable phase-field interval required",
        ));
    }
    let count = initial.transport_faces().len();
    let mut inverse_bound = 0usize;
    for face in initial.network().owner.faces.iter() {
        poll(&mut cancelled)?;
        let [i, j] = face.physical.cells;
        let offsets = &initial.network().owner.adjacent_offsets;
        // Two exchanged cells, all their neighbors, and two Fourier inverses.
        inverse_bound = inverse_bound
            .checked_add(offsets[i + 1] - offsets[i])
            .and_then(|n| n.checked_add(offsets[j + 1] - offsets[j]))
            .and_then(|n| n.checked_add(4))
            .ok_or(BinaryNonisothermalError::Invalid(
                "phase-field inverse work overflow",
            ))?;
    }
    let thermal = BinaryLiquidTransportSolveV1 {
        temperature_bracket: initial.temperature_bracket(),
        temperature_tolerance_k: controls.temperature_tolerance_k,
        energy_tolerance_j: controls.energy_tolerance_j / (4.0 * inverse_bound.max(1) as f64),
        maximum_iterations: controls.maximum_root_iterations,
    };
    let mut field = initial.clone();
    let mut patch = Patch::new(field.temperatures().len());
    let mut transferred_a_mol = vec![0.0; count];
    let mut transferred_energy_j = vec![0.0; count];
    let mut root_iterations = 0u64;
    let mut caloric_iterations = 0u64;
    for face in 0..count {
        poll(&mut cancelled)?;
        let physical = field.network().owner.faces[face].physical;
        let authored = field.transport_faces()[face];
        let [i, j] = physical.cells;
        let diffusivity = field.bulk.face_diffusivity(face, &mut cancelled)?;
        if physical.area_m2 == 0.0 {
            continue;
        }
        let geometry = finite(
            physical.area_m2 / physical.distance_m,
            "phase-field face geometry",
        )?;
        if !positive(geometry) {
            return Err(BinaryNonisothermalError::Invalid(
                "active phase-field geometry underflows",
            ));
        }
        let original = field.caloric_model().transport_face(
            [field.bulk.point(i), field.bulk.point(j)],
            AmountConcentration::new(field.network().model().molar_concentration_mol_m3()),
            physical.distance_m,
            BinaryLiquidTransportCoefficientsV1 {
                diffusivity_m2_s: diffusivity,
                thermal_conductivity_w_m_k: authored.thermal_conductivity_w_m_k,
                heat_of_transport_j_per_mol: authored.heat_of_transport_j_per_mol,
            },
            &mut cancelled,
        )?;
        let mobility = finite(
            geometry * original.mobility_mol2_k_j_m_s,
            "phase-field face mobility",
        )?;
        if original.mobility_mol2_k_j_m_s > 0.0 && mobility == 0.0 {
            return Err(BinaryNonisothermalError::Invalid(
                "active phase-field mobility underflows",
            ));
        }
        if mobility > 0.0 {
            let beta = reciprocal_jump(field.temperatures()[i], field.temperatures()[j])?;
            let initial_bulk_force = finite(
                -original.chemical_effort_j_per_mol_k + authored.heat_of_transport_j_per_mol * beta,
                "original phase-field bulk force",
            )?;
            let h = finite(
                original.exchange_enthalpy_j_per_mol + authored.heat_of_transport_j_per_mol,
                "phase-field transported enthalpy",
            )?;
            if let Some(root) = chemical_root(
                &field,
                face,
                duration_s,
                mobility,
                initial_bulk_force,
                controls,
                &mut cancelled,
            )? {
                let transfer = root.transfer;
                root_iterations += root.iterations as u64;
                let heat = finite(h * transfer, "phase-field paired chemical energy")?;
                caloric_iterations += recover_patch(
                    &mut field,
                    face,
                    root,
                    heat,
                    thermal,
                    &mut patch,
                    &mut cancelled,
                )? as u64;
                transferred_a_mol[face] = transfer;
                transferred_energy_j[face] = heat;
                field.bulk.face_diffusivity(face, &mut cancelled)?;
            }
        }
        let [tl, tr] = [
            field.temperatures()[i].value(),
            field.temperatures()[j].value(),
        ];
        let heat_mobility = finite(
            (geometry * authored.thermal_conductivity_w_m_k * tl) * tr,
            "phase-field Fourier entropy mobility",
        )?;
        if authored.thermal_conductivity_w_m_k > 0.0 && heat_mobility == 0.0 {
            return Err(BinaryNonisothermalError::Invalid(
                "phase-field heat mobility underflows",
            ));
        }
        if heat_mobility > 0.0 {
            let heat = super::super::solve::transfer_heat(
                &mut field.bulk,
                [i, j],
                duration_s,
                heat_mobility,
                controls,
                thermal,
                &mut cancelled,
            )?;
            if let Some((transfer, iterations, thermal_iterations)) = heat {
                let mut total = Sum::default();
                total.add(transferred_energy_j[face])?;
                total.add(transfer)?;
                transferred_energy_j[face] = total.value()?;
                root_iterations += iterations as u64;
                caloric_iterations += thermal_iterations as u64;
            }
        }
    }
    field.bulk.state.elapsed_s = elapsed_s;
    field.bulk.qualify(&mut cancelled)?;
    let (
        component_residual_mol,
        maximum_capacity_residual_mol,
        energy_residual_j,
        entropy_change_j_per_k,
        homogeneous_enthalpy_change_j,
        interfacial_energy_change_j,
    ) = initial.balances(&field, controls, &mut cancelled)?;
    poll(&mut cancelled)?;
    Ok(BinaryNonisothermalPhaseFieldCandidate {
        field,
        step: BinaryNonisothermalPhaseFieldStep {
            duration_s,
            transferred_a_mol,
            transferred_energy_j,
            root_iterations,
            caloric_iterations,
            component_residual_mol,
            maximum_capacity_residual_mol,
            energy_residual_j,
            entropy_change_j_per_k,
            homogeneous_enthalpy_change_j,
            interfacial_energy_change_j,
        },
    })
}

fn reciprocal_jump(left: Temperature, right: Temperature) -> Result<f64> {
    finite(
        ((left.value() - right.value()) / left.value()) / right.value(),
        "reciprocal temperature jump",
    )
    .map_err(Into::into)
}

/// Exact frozen-temperature graph stencil; every neighboring beta participates.
fn weighted_stencil(
    field: &BinaryNonisothermalPhaseField,
    cell: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, f64)> {
    let graph = &field.network().owner;
    let mut gradient = Sum::default();
    let mut diagonal = Sum::default();
    for &edge in
        &graph.adjacent_faces[graph.adjacent_offsets[cell]..graph.adjacent_offsets[cell + 1]]
    {
        poll(cancelled)?;
        let face = &graph.faces[edge];
        let neighbor = if face.physical.cells[0] == cell {
            face.physical.cells[1]
        } else {
            face.physical.cells[0]
        };
        let weight = finite(
            face.gradient_j
                * (0.5 / field.temperatures()[cell].value()
                    + 0.5 / field.temperatures()[neighbor].value()),
            "entropy-weighted graph edge",
        )?;
        diagonal.add(weight)?;
        gradient.add(finite(
            weight * (field.mole_fractions()[cell][0] - field.mole_fractions()[neighbor][0]),
            "entropy gradient",
        )?)?;
    }
    Ok((gradient.value()?, diagonal.value()?))
}

fn log_change(new: f64, old: f64) -> f64 {
    if new == old {
        return 0.0;
    }
    let (lo, hi, sign) = if new > old {
        (old, new, 1.0)
    } else {
        (new, old, -1.0)
    };
    let delta = (hi - lo) / lo;
    let rounded = 1.0 + delta;
    // Match the original caloric law's deterministic logarithm, preserving
    // nearby fractions with Kahan's log1p correction and exact antisymmetry.
    let value = if rounded == 1.0 {
        delta
    } else if rounded.is_finite() {
        fs_math::det::ln(rounded) * (delta / (rounded - 1.0))
    } else {
        fs_math::det::ln(hi) - fs_math::det::ln(lo)
    };
    sign * value
}

fn chemical_root(
    field: &BinaryNonisothermalPhaseField,
    face_index: usize,
    duration: f64,
    mobility: f64,
    initial_bulk_force: f64,
    controls: BinaryNonisothermalSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<ChemicalRoot>> {
    let graph = &field.network().owner;
    let face = &graph.faces[face_index];
    let cells = face.physical.cells;
    let old = cells.map(|c| field.mole_fractions()[c]);
    let amounts = cells.map(|c| field.bulk.state.amounts[c]);
    let remainders = cells.map(|c| field.bulk.state.amount_remainders[c]);
    let capacity = cells.map(|c| graph.capacities_mol[c]);
    let mut chi = [0.0; 2];
    let mut gradient = [0.0; 2];
    let mut diagonal = [0.0; 2];
    for k in 0..2 {
        poll(cancelled)?;
        chi[k] = field
            .caloric_model()
            .thermodynamics()
            .chi(field.temperatures()[cells[k]])
            .map_err(BinaryCaloricErrorV1::from)?
            .min(0.0);
        (gradient[k], diagonal[k]) = weighted_stencil(field, cells[k], cancelled)?;
    }
    let cross = finite(
        face.gradient_j
            * (0.5 / field.temperatures()[cells[0]].value()
                + 0.5 / field.temperatures()[cells[1]].value()),
        "implicit entropy gradient cross coefficient",
    )?;
    let mut force = Sum::default();
    force.add(initial_bulk_force)?;
    force.add(gradient[0] / capacity[0])?;
    force.add(-gradient[1] / capacity[1])?;
    let force = force.value()?;
    if force == 0.0 {
        return Ok(None);
    }
    let scale = finite(duration * mobility, "phase-field implicit scale")?;
    if !positive(scale) {
        return Err(BinaryNonisothermalError::Invalid(
            "phase-field integrated mobility underflows",
        ));
    }
    if scale * force == 0.0 {
        return Err(BinaryNonisothermalError::Invalid(
            "nonzero phase-field chemical transfer underflows",
        ));
    }
    let mut lower = -amounts[0][1].min(amounts[1][0]);
    let mut upper = amounts[0][0].min(amounts[1][1]);
    let mut transfer = 0.0;
    for iteration in 1..=controls.maximum_root_iterations {
        poll(cancelled)?;
        if transfer <= lower || transfer >= upper {
            break;
        }
        let mut next = amounts;
        let mut lows = remainders;
        for k in 0..2 {
            for s in 0..2 {
                (next[k][s], lows[k][s]) = super::super::super::solve::add_amount(
                    amounts[k][s],
                    remainders[k][s],
                    if k == s { -transfer } else { transfer },
                )?;
            }
        }
        if next.iter().flatten().any(|&n| !positive(n)) {
            break;
        }
        let fractions = [
            field.network().fractions(next[0], capacity[0])?,
            field.network().fractions(next[1], capacity[1])?,
        ];
        let delta = [fractions[0][0] - old[0][0], fractions[1][0] - old[1][0]];
        let mut value = Sum::default();
        value.add(force)?;
        for k in 0..2 {
            let sign = if k == 0 { 1.0 } else { -1.0 };
            value.add(finite(
                sign * R
                    * (log_change(fractions[k][0], old[k][0])
                        - log_change(fractions[k][1], old[k][1])
                        - 2.0 * chi[k] * delta[k]),
                "convex-split chemical increment",
            )?)?;
            value.add(finite(
                sign * (diagonal[k] * delta[k] - cross * delta[1 - k]) / capacity[k],
                "implicit full-graph increment",
            )?)?;
        }
        let value = value.value()?;
        let residual = finite(transfer - scale * value, "phase-field implicit residual")?;
        let tolerance = finite(
            controls.absolute_transfer_tolerance_mol
                + controls.relative_transfer_tolerance * transfer.abs().max((scale * value).abs()),
            "phase-field root tolerance",
        )?;
        // The convex split has r'(q)>=1: a residual bound is already a molar
        // error bound. A nonzero force is not silently discarded at q=0.
        if residual == 0.0 || (transfer != 0.0 && residual.abs() <= tolerance) {
            return Ok(Some(ChemicalRoot {
                transfer,
                amounts: next,
                remainders: lows,
                fractions,
                iterations: iteration,
            }));
        }
        if residual < 0.0 {
            lower = transfer;
        } else {
            upper = transfer;
        }
        let midpoint = 0.5 * lower + 0.5 * upper;
        if midpoint == lower || midpoint == upper {
            break;
        }
        let curvature = (R * (1.0 / fractions[0][0] + 1.0 / fractions[0][1] - 2.0 * chi[0]))
            / capacity[0]
            + (R * (1.0 / fractions[1][0] + 1.0 / fractions[1][1] - 2.0 * chi[1])) / capacity[1]
            + diagonal[0] / capacity[0] / capacity[0]
            + diagonal[1] / capacity[1] / capacity[1]
            + 2.0 * cross / capacity[0] / capacity[1];
        let derivative = 1.0 + scale * curvature;
        let proposal = transfer - residual / derivative;
        transfer = if positive(derivative)
            && proposal.is_finite()
            && proposal > lower
            && proposal < upper
            && proposal != transfer
        {
            proposal
        } else {
            midpoint
        };
    }
    Err(BinaryNonisothermalError::PairIterationLimit)
}

fn recover_patch(
    field: &mut BinaryNonisothermalPhaseField,
    face_index: usize,
    root: ChemicalRoot,
    transported_energy: f64,
    thermal: BinaryLiquidTransportSolveV1,
    patch: &mut Patch,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<usize> {
    patch.clear();
    let cells = field.network().owner.faces[face_index].physical.cells;
    patch.add(cells[0], -transported_energy)?;
    patch.add(cells[1], transported_energy)?;
    let new_x = |cell: usize| {
        if cell == cells[0] {
            root.fractions[0][0]
        } else if cell == cells[1] {
            root.fractions[1][0]
        } else {
            field.mole_fractions()[cell][0]
        }
    };
    let graph = &field.network().owner;
    for k in 0..2 {
        for &edge in &graph.adjacent_faces
            [graph.adjacent_offsets[cells[k]]..graph.adjacent_offsets[cells[k] + 1]]
        {
            poll(cancelled)?;
            // The pair edge is the only duplicate in the two incidence lists.
            if k == 1 && edge == face_index {
                continue;
            }
            let face = &graph.faces[edge];
            let [i, j] = face.physical.cells;
            let old = field.mole_fractions()[i][0] - field.mole_fractions()[j][0];
            let new = new_x(i) - new_x(j);
            let half_change = finite(
                (face.gradient_j * (new - old)) * (0.25 * new + 0.25 * old),
                "local half-edge energy change",
            )?;
            patch.add(i, -half_change)?;
            patch.add(j, -half_change)?;
        }
    }
    let mut iterations = 0usize;
    for &cell in &patch.cells {
        poll(cancelled)?;
        let next = if cell == cells[0] {
            root.amounts[0]
        } else if cell == cells[1] {
            root.amounts[1]
        } else {
            field.component_moles()[cell]
        };
        let energy = patch.energy[cell].value()?;
        if energy == 0.0 && next == field.component_moles()[cell] {
            continue;
        }
        let balance = field
            .caloric_model()
            .homogeneous_transport_candidate(
                &field.bulk.caloric_state(cell)?,
                next,
                energy,
                thermal,
                &mut *cancelled,
            )
            .map_err(temperature_error)?;
        iterations += balance.iterations();
        field.bulk.temperatures[cell] = balance.final_state().temperature();
    }
    for k in 0..2 {
        field.bulk.state.amounts[cells[k]] = root.amounts[k];
        field.bulk.state.amount_remainders[cells[k]] = root.remainders[k];
        field.bulk.state.fractions[cells[k]] = root.fractions[k];
    }
    Ok(iterations)
}

fn temperature_error(error: BinaryCaloricErrorV1) -> BinaryNonisothermalError {
    match error {
        BinaryCaloricErrorV1::HeatOutsideBracket { .. } => {
            BinaryNonisothermalError::PairOutsideTemperatureDomain
        }
        error => error.into(),
    }
}
