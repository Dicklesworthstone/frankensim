//! Implicit scalar updates of the two positive Onsager rank-one blocks.

use super::*;
use fs_qty::AmountConcentration;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

#[derive(Clone, Copy)]
enum Mode {
    Chemical { transported_enthalpy: f64 },
    Heat,
}

#[derive(Clone)]
struct Pair {
    amounts: [[f64; 2]; 2],
    remainders: [[f64; 2]; 2],
    fractions: [[f64; 2]; 2],
    states: [BinaryCaloricStateV1; 2],
    force: f64,
    derivative: f64,
    caloric_iterations: usize,
}
struct Root {
    pair: Pair,
    transfer: f64,
    iterations: usize,
    caloric_iterations: usize,
}

pub(super) fn step(
    initial: &BinaryNonisothermalField,
    duration_s: f64,
    controls: BinaryNonisothermalSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<BinaryNonisothermalCandidate> {
    poll(&mut cancelled)?;
    admit_controls(controls)?;
    let elapsed_s = finite(
        initial.elapsed_s() + duration_s,
        "nonisothermal elapsed time",
    )?;
    if !positive(duration_s) || elapsed_s <= initial.elapsed_s() {
        return Err(BinaryNonisothermalError::Invalid(
            "positive representable physical interval required",
        ));
    }
    let mut field = initial.clone();
    let count = field.owner.faces.len();
    let mut transferred_a_mol = vec![0.0; count];
    let mut transferred_energy_j = vec![0.0; count];
    let mut root_iterations = 0u64;
    let mut caloric_iterations = 0u64;
    // A global energy budget is not multiplied by faces or temperature solves.
    let local_energy_tolerance = controls.energy_tolerance_j / (4.0 * count.max(1) as f64);
    let local = BinaryLiquidTransportSolveV1 {
        temperature_bracket: field.owner.temperature_bracket,
        temperature_tolerance_k: controls.temperature_tolerance_k,
        energy_tolerance_j: local_energy_tolerance,
        maximum_iterations: controls.maximum_root_iterations,
    };
    for face in 0..count {
        poll(&mut cancelled)?;
        let physical = field.network().owner.faces[face].physical;
        let transport = field.owner.faces[face];
        let [i, j] = physical.cells;
        let d = field.face_diffusivity(face, &mut cancelled)?;
        if physical.area_m2 == 0.0 {
            continue;
        }
        let geometry = finite(
            physical.area_m2 / physical.distance_m,
            "nonisothermal face geometry",
        )?;
        if geometry == 0.0 {
            return Err(BinaryNonisothermalError::Invalid(
                "positive nonisothermal face geometry underflows",
            ));
        }
        let flux = field.caloric_model().transport_face(
            [field.point(i), field.point(j)],
            concentration(&field),
            physical.distance_m,
            BinaryLiquidTransportCoefficientsV1 {
                diffusivity_m2_s: d,
                thermal_conductivity_w_m_k: transport.thermal_conductivity_w_m_k,
                heat_of_transport_j_per_mol: transport.heat_of_transport_j_per_mol,
            },
            &mut cancelled,
        )?;
        let mobility = finite(
            geometry * flux.mobility_mol2_k_j_m_s,
            "chemical face mobility",
        )?;
        if flux.mobility_mol2_k_j_m_s > 0.0 && mobility == 0.0 {
            return Err(BinaryNonisothermalError::Invalid(
                "positive chemical face mobility underflows",
            ));
        }
        if mobility > 0.0 {
            let h = finite(
                flux.exchange_enthalpy_j_per_mol + transport.heat_of_transport_j_per_mol,
                "transported partial enthalpy",
            )?;
            if let Some(root) = pair_root(
                &field,
                [i, j],
                Mode::Chemical {
                    transported_enthalpy: h,
                },
                duration_s,
                mobility,
                controls,
                local,
                &mut cancelled,
            )? {
                transferred_a_mol[face] = root.transfer;
                transferred_energy_j[face] = finite(h * root.transfer, "paired chemical energy")?;
                root_iterations += root.iterations as u64;
                caloric_iterations += root.caloric_iterations as u64;
                apply_pair(&mut field, [i, j], root.pair);
                // The completed chemical substage remains inside its original
                // molecular source support before another physical substage.
                field.face_diffusivity(face, &mut cancelled)?;
            }
        }
        let [tl, tr] = [field.temperatures[i].value(), field.temperatures[j].value()];
        let heat_mobility = finite(
            (geometry * transport.thermal_conductivity_w_m_k * tl) * tr,
            "Fourier entropy mobility",
        )?;
        if transport.thermal_conductivity_w_m_k > 0.0 && heat_mobility == 0.0 {
            return Err(BinaryNonisothermalError::Invalid(
                "positive Fourier entropy mobility underflows",
            ));
        }
        if heat_mobility > 0.0 {
            if let Some((transfer, iterations, thermal_iterations)) = transfer_heat(
                &mut field,
                [i, j],
                duration_s,
                heat_mobility,
                controls,
                local,
                &mut cancelled,
            )? {
                let mut energy = Sum::default();
                energy.add(transferred_energy_j[face])?;
                energy.add(transfer)?;
                transferred_energy_j[face] = energy.value()?;
                root_iterations += iterations as u64;
                caloric_iterations += thermal_iterations as u64;
            }
        }
    }
    field.state.elapsed_s = elapsed_s;
    field.qualify(&mut cancelled)?;
    let (
        component_residual_mol,
        maximum_capacity_residual_mol,
        energy_residual_j,
        entropy_change_j_per_k,
    ) = initial.balances(&field, controls, &mut cancelled)?;
    poll(&mut cancelled)?;
    Ok(BinaryNonisothermalCandidate {
        field,
        step: BinaryNonisothermalStep {
            duration_s,
            transferred_a_mol,
            transferred_energy_j,
            root_iterations,
            caloric_iterations,
            component_residual_mol,
            maximum_capacity_residual_mol,
            energy_residual_j,
            entropy_change_j_per_k,
        },
    })
}

fn concentration(field: &BinaryNonisothermalField) -> AmountConcentration {
    AmountConcentration::new(field.network().model().molar_concentration_mol_m3())
}
fn zero_coefficients() -> BinaryLiquidTransportCoefficientsV1 {
    BinaryLiquidTransportCoefficientsV1 {
        diffusivity_m2_s: 0.0,
        thermal_conductivity_w_m_k: 0.0,
        heat_of_transport_j_per_mol: 0.0,
    }
}
fn apply_pair(field: &mut BinaryNonisothermalField, cells: [usize; 2], pair: Pair) {
    for k in 0..2 {
        field.state.amounts[cells[k]] = pair.amounts[k];
        field.state.amount_remainders[cells[k]] = pair.remainders[k];
        field.state.fractions[cells[k]] = pair.fractions[k];
        field.temperatures[cells[k]] = pair.states[k].temperature();
    }
}

/// Shared fixed-composition Fourier kernel. Interface energy is unchanged, so
/// the local-temperature phase owner uses this same original caloric solve.
pub(super) fn transfer_heat(
    field: &mut BinaryNonisothermalField,
    cells: [usize; 2],
    duration: f64,
    mobility: f64,
    controls: BinaryNonisothermalSolve,
    thermal: BinaryLiquidTransportSolveV1,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<(f64, usize, usize)>> {
    let Some(root) = pair_root(
        field,
        cells,
        Mode::Heat,
        duration,
        mobility,
        controls,
        thermal,
        cancelled,
    )?
    else {
        return Ok(None);
    };
    let receipt = (root.transfer, root.iterations, root.caloric_iterations);
    apply_pair(field, cells, root.pair);
    Ok(Some(receipt))
}

fn pair_root(
    field: &BinaryNonisothermalField,
    cells: [usize; 2],
    mode: Mode,
    duration: f64,
    mobility: f64,
    controls: BinaryNonisothermalSolve,
    thermal: BinaryLiquidTransportSolveV1,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<Root>> {
    let mut best = sample(field, cells, mode, 0.0, thermal, cancelled)?
        .ok_or(BinaryNonisothermalError::PairOutsideTemperatureDomain)?;
    if best.force == 0.0 {
        return Ok(None);
    }
    let sign = best.force.signum();
    let scale = finite(duration * mobility, "implicit Onsager step scale")?;
    if !positive(scale) {
        return Err(BinaryNonisothermalError::Invalid(
            "nonzero pair transfer scale underflows",
        ));
    }
    let mut hi = match mode {
        Mode::Chemical { .. } => {
            let limit = if sign > 0.0 {
                best.amounts[0][0].min(best.amounts[1][1])
            } else {
                best.amounts[0][1].min(best.amounts[1][0])
            };
            limit.next_down()
        }
        Mode::Heat => {
            let [cold, hot] = thermal.temperature_bracket;
            let end = if sign > 0.0 { [cold, hot] } else { [hot, cold] };
            let mut limit = f64::INFINITY;
            for k in 0..2 {
                let target = field.caloric_model().homogeneous(end[k], best.amounts[k])?;
                let change = field
                    .caloric_model()
                    .homogeneous_transport_change(&best.states[k], &target)?;
                limit = limit.min(change.enthalpy_change_j.abs());
            }
            limit
        }
    };
    if !positive(hi) {
        return Err(BinaryNonisothermalError::PairOutsideTemperatureDomain);
    }
    let absolute = match mode {
        Mode::Chemical { .. } => controls.absolute_transfer_tolerance_mol,
        Mode::Heat => thermal.energy_tolerance_j,
    };
    let mut lo = 0.0;
    let mut best_x = 0.0_f64;
    let mut best_residual = finite(-scale * sign * best.force, "initial pair residual")?;
    if best_residual == 0.0 {
        return Err(BinaryNonisothermalError::Invalid(
            "nonzero integrated pair force underflows",
        ));
    }
    let mut lower = best.clone();
    let mut lower_residual = best_residual;
    let mut caloric_iterations = 0usize;
    let mut consecutive_newton = 2;
    let mut has_upper_root_bound = false;
    for iteration in 1..=controls.maximum_root_iterations {
        poll(cancelled)?;
        let tolerance = absolute + controls.relative_transfer_tolerance * best_x.abs();
        if best_residual == 0.0
            || (best_x > 0.0
                && has_upper_root_bound
                && hi - lo <= tolerance
                && best_residual.abs() <= tolerance)
        {
            return Ok(Some(Root {
                pair: best,
                transfer: sign * best_x,
                iterations: iteration - 1,
                caloric_iterations,
            }));
        }
        let midpoint = 0.5 * lo + 0.5 * hi;
        if midpoint == lo || midpoint == hi {
            break;
        }
        let derivative = 1.0 - scale * best.derivative;
        let proposal = best_x - best_residual / derivative;
        let x = if consecutive_newton < 2
            && derivative.is_finite()
            && derivative > 0.0
            && proposal.is_finite()
            && proposal > lo
            && proposal < hi
        {
            proposal
        } else {
            midpoint
        };
        let Some(next) = sample(field, cells, mode, sign * x, thermal, cancelled)? else {
            hi = x;
            has_upper_root_bound = false;
            if best_x > hi {
                best = lower.clone();
                best_x = lo;
                best_residual = lower_residual;
            }
            consecutive_newton = 2;
            continue;
        };
        caloric_iterations += next.caloric_iterations;
        let residual = finite(x - scale * sign * next.force, "implicit pair residual")?;
        if residual >= 0.0 {
            hi = x;
            has_upper_root_bound = true;
        } else {
            lo = x;
            lower = next.clone();
            lower_residual = residual;
        }
        if residual.abs() < best_residual.abs() || best_x < lo || best_x > hi {
            best_x = x;
            best_residual = residual;
            best = next;
        }
        consecutive_newton = if x == midpoint {
            0
        } else {
            consecutive_newton + 1
        };
    }
    let tolerance = absolute + controls.relative_transfer_tolerance * best_x.abs();
    if best_residual == 0.0
        || (best_x > 0.0
            && has_upper_root_bound
            && hi - lo <= tolerance
            && best_residual.abs() <= tolerance)
    {
        Ok(Some(Root {
            pair: best,
            transfer: sign * best_x,
            iterations: controls.maximum_root_iterations,
            caloric_iterations,
        }))
    } else if !has_upper_root_bound {
        Err(BinaryNonisothermalError::PairOutsideTemperatureDomain)
    } else {
        Err(BinaryNonisothermalError::PairIterationLimit)
    }
}

fn sample(
    field: &BinaryNonisothermalField,
    cells: [usize; 2],
    mode: Mode,
    transfer: f64,
    thermal: BinaryLiquidTransportSolveV1,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<Pair>> {
    poll(cancelled)?;
    let mut amounts = cells.map(|c| field.state.amounts[c]);
    let mut remainders = cells.map(|c| field.state.amount_remainders[c]);
    let initial = [
        field.caloric_state(cells[0])?,
        field.caloric_state(cells[1])?,
    ];
    let heat = match mode {
        Mode::Chemical {
            transported_enthalpy,
        } => {
            for k in 0..2 {
                let q = if k == 0 { -transfer } else { transfer };
                for s in 0..2 {
                    (amounts[k][s], remainders[k][s]) = super::super::solve::add_amount(
                        amounts[k][s],
                        remainders[k][s],
                        if s == 0 { q } else { -q },
                    )?;
                    if !positive(amounts[k][s]) {
                        return Ok(None);
                    }
                }
            }
            finite(
                transported_enthalpy * transfer,
                "chemical enthalpy transfer",
            )?
        }
        Mode::Heat => transfer,
    };
    let fractions = [
        field
            .network()
            .fractions(amounts[0], field.network().owner.capacities_mol[cells[0]])?,
        field
            .network()
            .fractions(amounts[1], field.network().owner.capacities_mol[cells[1]])?,
    ];
    let mut states = initial.clone();
    let mut caloric_iterations = 0usize;
    if transfer != 0.0 {
        for k in 0..2 {
            let balance = field.caloric_model().homogeneous_transport_candidate(
                &initial[k],
                amounts[k],
                if k == 0 { -heat } else { heat },
                thermal,
                &mut *cancelled,
            );
            match balance {
                Ok(balance) => {
                    caloric_iterations += balance.iterations();
                    states[k] = balance.final_state().clone();
                }
                Err(BinaryCaloricErrorV1::HeatOutsideBracket { .. }) => return Ok(None),
                Err(error) => return Err(error.into()),
            }
        }
    }
    let points = [0, 1].map(|k| BinaryLiquidTransportPointV1 {
        temperature: states[k].temperature(),
        mole_fractions: fractions[k],
    });
    let [tl, tr] = points.map(|p| p.temperature.value());
    let beta = finite(((tl - tr) / tl) / tr, "reciprocal temperature jump")?;
    let mut derivative = 0.0;
    let force = match mode {
        Mode::Chemical {
            transported_enthalpy: h,
        } => {
            let evaluated = field.caloric_model().transport_face(
                points,
                concentration(field),
                1.0,
                zero_coefficients(),
                &mut *cancelled,
            )?;
            for k in 0..2 {
                let local = field.caloric_model().transport_face(
                    [points[k]; 2],
                    concentration(field),
                    1.0,
                    zero_coefficients(),
                    &mut *cancelled,
                )?;
                let x = fractions[k];
                let chi = field
                    .caloric_model()
                    .thermodynamics()
                    .chi(points[k].temperature)
                    .map_err(BinaryCaloricErrorV1::from)?;
                let gamma =
                    ((x[0] - x[1]) * (x[0] - x[1]) + 2.0 * (2.0 - chi) * x[0] * x[1]) / x[0] / x[1];
                let n = amounts[k][0] + amounts[k][1];
                let cp = states[k].isobaric_heat_capacity_j_per_k()?;
                let thermal_effort =
                    (h - local.exchange_enthalpy_j_per_mol) / points[k].temperature.value();
                derivative -= R * gamma / n + thermal_effort * thermal_effort / cp;
            }
            finite(
                -evaluated.chemical_effort_j_per_mol_k
                    + (h - evaluated.exchange_enthalpy_j_per_mol) * beta,
                "conjugate chemical force",
            )?
        }
        Mode::Heat => {
            for state in &states {
                let t = state.temperature().value();
                derivative -= (1.0 / state.isobaric_heat_capacity_j_per_k()? / t) / t;
            }
            beta
        }
    };
    Ok(Some(Pair {
        amounts,
        remainders,
        fractions,
        states,
        force,
        derivative,
        caloric_iterations,
    }))
}
